#!/usr/bin/env python
"""
Created on Mon April 11 03:52:02 2016

This script is used to call do_shuffle.py for a single target.
It allows you to specify the target RA, Dec, the desired
IFU slot for it to land on, a small offset from the center of that
IFU if needed, and the track of the night.

Using the default shuffle.cfg, a visualization tool pops up for
the new shuffled coordinates and allows an inspection of the field.

Also produced are ACAM images along with the (RA, Dec) of the GP(s),
(RA, Dec) of the IHMP center, and the (x, y) pixel of the star(s) in the
ACAM image for the Telescope operator.

@author: gregz
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

try:  # python 3
    import builtins
except ImportError:
    import __builtin__ as builtins

import argparse as ap
import inspect
import logging
import os
# from pprint import pformat
import sys
import textwrap as tw
import warnings
import datetime
from pyhetdex.coordinates.tangent_projection import TangentPlane as TP
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

import numpy
from six.moves import configparser

from . import shuffle
from . import parang
from . import findStars
from . import __full_version__


# replace the standard print with a version that warns to use logging instead
print_ = print


def _print(*args, **kwargs):
    # caller = inspect.getframeinfo(inspect.stack()[1][0])

    # with warnings.catch_warnings():
    #    warnings.simplefilter('always')
    #     warnings.warn('\n\t{t.filename}:{t.lineno}: please replace the print'
    #                   ' with a log of the appropriate type'.format(t=caller))
    print_(*args, **kwargs)


setattr(builtins, 'print', _print)


class IFUException(ValueError):
    '''Exception raised when the IFU does not exist'''
    pass


def load_config(configfile):
    """Read configuration file

    Parameters
    ----------
    configfile : string
        name of the configuration file

    Returns
    -------
    config : :class:`~ConfigParser.ConfigParser` instance
    """
    config = configparser.ConfigParser(defaults={'xpa_method': 'local'})
    if not config.read(configfile):
        msg = ("File '{}' does not exist. Using the default one."
               " You can get the full set of configuration files with the"
               " command: ``shuffle_config``".format(configfile))
        warnings.warn(msg)
        config_def = os.path.join(os.path.dirname(__file__), 'configs',
                                  'shuffle.cfg')
        if not config.read(config_def):
            msg = ("ERROR: Failed to load the default configuration file.")
            raise ap.ArgumentTypeError(msg.format(configfile))

    return config


class RaToFloat(ap.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        '''Convert ra to floats.

        Parameters
        ----------
        parser : current parser
        namespace : name space
        value : string
            value from the command line
        option_string : string
        '''
        ra = values.split(':')
        try:
            if len(ra) == 1:
                ra_deg = float(ra[0])
            else:
                rah, ram, ras = ra
                ra_deg = float(rah) + float(ram)/60. + float(ras)/3600.
        except ValueError:
            raise parser.error('Ra can be provided as number or as string of'
                               ' colon separated numbers, like "hh:mm:ss".'
                               ' "{}" is not allowed'.format(values))
        setattr(namespace, self.dest, ra_deg)


class DecToFloat(ap.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        '''Convert the dec to float

        Parameters
        ----------
        parser : current parser
        namespace : name space
        value : string
            value from the command line
        option_string : string
        '''
        dec = values.strip().split(':')
        try:
            if len(dec) == 1:
                dec_deg = float(dec[0])
            else:
                decd, decm, decs = dec
                if decd.startswith('-'):
                    dec_deg = float(decd) - float(decm)/60. - float(decs)/3600.
                else:
                    dec_deg = float(decd) + float(decm)/60. + float(decs)/3600.
        except ValueError:
            raise parser.error('Dec can be provided as number or as string of'
                               ' colon separated numbers, like "dd:mm:ss".'
                               ' "{}" is not allowed'.format(values))
        setattr(namespace, self.dest, dec_deg)


class PositionalParser(ap.ArgumentParser):
    def error(self, message):
        msg = tw.dedent('''
                        Are you by any chance trying to provide a negative
                        ``dec`` in the form ``dd:mm:ss`` to ``%s``? If the
                        answer is yes, use the float representation of dec or
                        add ``--`` before the positional arguments. E.g.:

                            do_shuffle [options] 23:40:05.43 -1.32195
                                       radius {0,1,2} ifuslot [ra_offset]
                                       [dec_offset] [target]

                        or

                            do_shuffle [options] -- 23:40:05.43 -1:19:19.02
                                       radius {0,1,2} ifuslot [ra_offset]
                                       [dec_offset] [target]
                        ''' % (self.prog))
        print_(msg)

        super(PositionalParser, self).error(message)


def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "Shuffle the HETDEX shots"

    parser = PositionalParser(description=description,
                              formatter_class=ap.ArgumentDefaultsHelpFormatter,
                              )

    parser.add_argument("--verbose", '-v', action="count", default=0,
                        help="""Increase verbosity, can be called multiple
                        times""")
    parser.add_argument('--version', '-V', action='version',
                        version=__full_version__)

    parser.add_argument("ra", action=RaToFloat, help='''ra of the target in
                        hours.''')

    parser.add_argument("dec", action=DecToFloat,
                        help='dec of the target in degrees')

    parser.add_argument("radius", type=float,
                        help='''Search radius (in degrees) for proper guide
                        stars and wfs stars.''')

    parser.add_argument("track", type=int, choices=[0, 1],
                        help='Type of track: 0: East 1: West')

    parser.add_argument("ifuslot", type=int,
                        help='''ifuslot of the target. For example: 075.
                        Note: the boresight is 000.''')

    parser.add_argument("x_offset", nargs='?', type=float,
                        help='''Small X offset from center of ifuslot in
                        arcsecs''', default=0.0)

    parser.add_argument("y_offset", nargs='?', type=float,
                        help='''Small Y offset from center of ifuslot in
                        arcsecs''', default=0.0)

    parser.add_argument("target", nargs='?', type=str,
                        help='''Name of the target observation''',
                        default='Test')

    parser.add_argument("--epoch", type=float,
                        help='''Epoch of target coordinates''',
                        default=2000.0)

    parser.add_argument("--pmRA", type=float,
                        help='''Proper motion in right ascension (mas/yr)''',
                        default=None)

    parser.add_argument("--pmDEC", type=float,
                        help='''Proper motion in delination (mas/yr)''',
                        default=None)

    parser.add_argument("--gaia2radius", type=float,
                        help='''Search radius for target match to gaia2 in "''',
                        default=30.)

    parser.add_argument("-c", "--config", help="""Name of the configuration
                        file. When parsing the command line, the file is loaded
                        into a configuration object""",
                        default="./shuffle.cfg", type=load_config)

    parser.add_argument("--limited_probe_annulus",
                        help='''Use 0.31 deg for inner probe radius''',
                        action="count", default=0)

    # ARGUMENTS IN CONFIG FILE

    for i in ['1', '2']:
        parser.add_argument("--force-gp" + i, metavar='id' + i, help="""Force
                            the guide probe {} to use the star with the given
                            id""".format(i))

    for i in ['1', '2']:
        parser.add_argument("--force-wfs" + i, metavar='id' + i, help="""Force
                            the wavefront sensor probe {} to use the star with
                            the given id""".format(i))

    parser.add_argument("--localcat", type=str, help="""Use a local catalog
                        for star searches""", default=None)

    parser.add_argument("--interactive", type=str, help="""Use interactive
                        mode for guide star search.""", default=None)

    parser.add_argument("--hpf", type=str,
                        help="""Use hpf additional selection""", default=None)

    parser.add_argument("--use_brightness", type=str,
                        help="""Use the brightness
                        of the probe stars as a selection criteria.  If False,
                        then the distance from the target angle is
                        minimizes.""", default=None)

    parser.add_argument("--fplane_file", type=str, help="""Fplane file to be
                        used.""", default=None)

    parser.add_argument("--az", nargs='?', type=float,
                        help='''Manual Non-optimal AZ input''', default=None)

    for i in ['1', '2']:
        parser.add_argument("--dpatrol_g" + i + 'min',
                            help="""Guide probe minimum patrol angle""",
                            type=float, default=None)
        parser.add_argument("--dpatrol_g" + i + 'max',
                            help="""Guide probe maximum patrol angle""",
                            type=float, default=None)
        parser.add_argument("--dpatrol_g" + i + 'targ',
                            help="""Guide probe target patrol angle""",
                            type=float, default=None)

    for i in ['1', '2']:
        parser.add_argument("--dpatrol_w" + i + 'min',
                            help="""Wavefront probe minimum patrol angle""",
                            type=float, default=None)
        parser.add_argument("--dpatrol_w" + i + 'max',
                            help="""Wavefront probe maximum patrol angle""",
                            type=float, default=None)
        parser.add_argument("--dpatrol_w" + i + 'targ',
                            help="""Wavefront probe target patrol angle""",
                            type=float, default=None)

    parser.add_argument("--catalog", type=str, help="""Star catalog to use,
                        for example: GAIA""", default=None)

    for st in ['gp', 'wfs']:
        parser.add_argument("--" + st + "pickcol", type=int,
                            help="""Probe filter to be used.
                            For SDSS DR9 1 = u, 2 = g, 3 = r, 4 = i , 5 = z
                            For USNO A2  2 = B, 3 = R""", default=None)

    parser.add_argument("--visualize", type=str, help="""Make visualization
                        for the shot.""", default=None)

    parser.add_argument("--visualize_ACAM", type=str, help="""Make visualization
                        for the ACAM.""", default=None)

    parser.add_argument("--visualize_probestars", type=str,
                        help="""Make visualization for the probestars.""",
                        default=None)

    parser.add_argument("--acam", type=float, help="""Offset angle for ACAM""",
                        default=None)

    parser.add_argument("--fplane", type=float,
                        help="""Offset angle for fplane""",
                        default=None)

    parser.add_argument("--probes", type=float,
                        help="""Offset angle for probes""",
                        default=None)

    for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'ifu', 'acam', 'fplane',
              'allprobes', 'allcams']:
        for p in ['magadd', 'nstars', 'minsep', 'magmin', 'magmax']:
            parser.add_argument('--' + o + '_' + p,
                                help="""{} {}""".format(o, p), default=None)

    args = parser.parse_args(args=argv)

    if args.allprobes_magadd is not None:
        for o in ['guide1', 'guide2', 'wfs1', 'wfs2']:
            if getattr(args, o+'_'+'magadd') is None:
                setattr(args, o+'_'+'magadd', args.allprobes_magadd)

    if args.allcams_magadd is not None:
        for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'acam']:
            if getattr(args, o+'_'+'magadd') is None:
                setattr(args, o+'_'+'magadd', args.allcams_magadd)

    if args.fplane is not None:
        args.probes = args.fplane

    args.ifuslot = "{:03d}".format(args.ifuslot)

    if args.limited_probe_annulus:
        args.config.set('General', 'dpatrol_min', '0.31')
    return args


def cl_to_config(args):
    """Copy the values of some command line options into the configuration.

    Values copied: force_* into the corresponding entries in the ``General``
    section if the value is not ``None``

    Parameters
    ----------
    args : Namespace
        parsed command line arguments

    Returns
    -------
    args : Namespace
        parsed command line arguments
    """
    # move not None values into the General, visualization, offsets section
    # of the configuration
    general_args, viz_args, offset_args, mag_args = [], [], [], []

    list_type = ['General', 'visualisation', 'offsets', 'MagLimits']
    for o in ['gp1', 'gp2', 'wfs1', 'wfs2']:
        general_args.append('force_'+o)
    general_args.extend(['localcat', 'interactive', 'use_brightness',
                         'fplane_file', 'gppickcol', 'wfspickcol', 'visualize',
                         'visualize_ACAM', 'catalog', 'hpf'])
    for o in ['g1', 'g2', 'w1', 'w2']:
        for p in ['min', 'max', 'targ']:
            general_args.append('dpatrol_' + o + p)

    viz_args.append('visualize_probestars')

    offset_args.extend(['acam', 'fplane', 'probes'])

    for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'ifu', 'acam', 'fplane']:
        for p in ['magadd', 'nstars', 'minsep', 'magmin', 'magmax']:
            mag_args.append(o+'_'+p)

    big_list = [general_args, viz_args, offset_args, mag_args]
    for list_name, p in zip(big_list, list_type):
        for o in list_name:
            value = getattr(args, o)
            if value is not None:
                args.config.set(p, o, value)

    return args


def get_ifuslot_list(args):
    """Get the ifu centers and ids of the VIRUS+LRS IFUs for the adjusted
    fplane file.
    """
    fplane_adj_file = args.config.get("General", "fplane_file")
    fplane = shuffle.fplane_parser.FPlane(fplane_adj_file)
    ifu_centers = numpy.array([[ifu.x, ifu.y] for ifu in fplane.ifus])
    ifu_centers = ifu_centers / 3600.
    ifu_ids = numpy.array([[ifu.ifuslot] for ifu in fplane.ifus])
    return ifu_centers, ifu_ids


def check_fplane_file(args):
    '''If the fplane file does not exist, in the config file point to the one
    shipped with shuffle.

    Parameters
    ----------
    args : Namespace
        parsed configuration options

    Returns
    -------
    args : Namespace
        parsed configuration options
    '''
    fplane_file = args.config.get("General", "fplane_file")
    if not os.path.exists(fplane_file):
        msg = ("File '{}' does not exist. Using the default one."
               " You can get the full set of configuration files with the"
               " command: ``shuffle_config``".format(fplane_file))
        warnings.warn(msg)
        fplane_def = os.path.join(os.path.dirname(__file__), 'configs',
                                  'fplane.txt')
        args.config.set('General', 'fplane_file', fplane_def)
    return args


def make_dirs(args):
    '''To run shuffle, some directories are needed.
    The name of the directories is given in the configuration file.
    This function tries to make make them, and ignore any error

    Parameters
    ----------
    args : Namespace
        parsed configuration options
    '''
    for i in ['cache', 'images', 'acam_images', 'ds9_regions', 'gp_images']:
        dir_name = args.config.get("directories", i)
        try:
            os.mkdir(dir_name)
        except OSError:
            pass


def setup_logging(args):
    '''Set up a logger for shuffle with a name ``shuffle``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    fmt = '[%(levelname)s - %(asctime)s] %(message)s'
    if args.verbose == 0:
        level = logging.WARNING
    elif args.verbose == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
        fmt = '[%(levelname)s - %(filename)s - %(asctime)s] %(message)s'
    fmt = logging.Formatter(fmt)

    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    handler.setLevel(level)

    log = logging.getLogger('shuffle')
    log.setLevel(logging.DEBUG)
    log.addHandler(handler)


def set_xpa_method(config):
    '''Push the xpa method to the environment before loading ds9. By default is
    set to 'local'. It can be customised with the [ds9][xpa_method] option in
    the configuration file.

    Parameters
    ----------
    config : :class:`configparser.ConfigParser` instance
        configuration
    '''

    # set the xpa method; defaults to local to avoid locks when loading pyds9
    os.environ['XPA_METHOD'] = config.get('ds9', 'xpa_method')


def ifu_index(ifu_ids, ifuslot):
    '''Return the index of ifuslot in the list of ifu_ids.

    Parameters
    ----------
    ifu_ids : nd-array
        array containing the list of ifuslots
    ifuslot : string
        IFUSLOT to find

    Returns
    -------
    ifu_ind : nd-array
        array of one element containing the index of ifuslot
    '''
    ifu_ind = numpy.where(ifu_ids == ifuslot)[0]

    if not list(ifu_ind):
        msg = 'The required IFUSLOT "{id_}" is not in the fplane file.\n'
        msg += 'Check the fplane file for a list of available IFUSLOTS'
        raise IFUException(msg.format(id_=ifuslot))

    return ifu_ind


def update_coords_for_proper_motion(ra, dec, epoch, log, args):
    current_epoch = Time(datetime.datetime.now()).decimalyear
    if (args.pmRA is not None) and (args.pmDEC is not None):
        deltaRA = ((current_epoch - epoch) * args.pmRA /
                   1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
        deltaDE = ((current_epoch - epoch) * args.pmDEC /
                   1e3 / 3600.)
        ra += deltaRA
        dec += deltaDE
        log.info('Change in RA from proper motion: %0.2f"' %
                 (deltaRA * 3600.))
        log.info('Change in Dec from proper motion: %0.2f"' %
                 (deltaDE * 3600.))
        return ra, dec
    try:
        pmcat = findStars.queryGAIA2(ra, dec, args.gaia2radius / 3600.)
        c1 = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
        deltaRA = ((epoch - 2000.) * pmcat[:, -2] /
                   1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
        deltaDE = ((epoch - 2000.) * pmcat[:, -1] /
                   1e3 / 3600.)
        c2 = SkyCoord((pmcat[:, 2]+deltaRA)*u.deg, (pmcat[:, 3]+deltaDE)*u.deg,
                      frame='icrs')
        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        if d2d.arcsec < 2.:
            deltaRA = ((current_epoch - epoch) * pmcat[idx, -2] /
                       1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
            deltaDE = ((current_epoch - epoch) * pmcat[idx, -1] /
                       1e3 / 3600.)

            idx, d2d, d3d = c1.match_to_catalog_sky(c2)
            ra += deltaRA
            dec += deltaDE
            log.info('Change in RA from proper motion: %0.2f"' %
                     (deltaRA * 3600.))
            log.info('Change in Dec from proper motion: %0.2f"' %
                     (deltaDE * 3600.))
            return ra, dec
        else:
            return ra, dec
    except Exception:
        log.warning("Could NOT update astrometry with proper motion")
        return ra, dec


def main():

    args = parse_args()
    args = cl_to_config(args)

    args = check_fplane_file(args)
    setup_logging(args)
    set_xpa_method(args.config)
    make_dirs(args)

    log = logging.getLogger('shuffle')

    # read the parameters from command line
    # assume input ra in hms, and convert to degrees by *15
    ra = args.ra * 15
    dec = args.dec
    ra, dec = update_coords_for_proper_motion(ra, dec, args.epoch, log, args)
    radius = args.radius
    track = args.track
    obstime = 0
    ill = 0
    trans = 0
    seeing = 0
    ifu_centers, ifu_ids = get_ifuslot_list(args)
    ifu_ind = ifu_index(ifu_ids, args.ifuslot)
    dec_center = dec*1.
    for i in numpy.arange(4):
        pas = parang.parang(dec_center, az_input=args.az, verbose=False)
        if pas[0] == 0:
            log.error("No track is available for DEC=%f.", dec)
            sys.exit(0)
        if pas[0] == 1:
            log.warning("Warning: Only one track is available for DEC= %f, two"
                        " tracks are identical.", dec)
        pa = pas[int(track) + 1] + args.config.getfloat('offsets', 'fplane')
        tan_plane = TP(ra, dec, pa, 1., 1.)
        ra_center, dec_center = tan_plane.xy2raDec(-3600. *
                                                   (ifu_centers[ifu_ind, 0][0])
                                                   - args.x_offset,
                                                   -3600. *
                                                   (ifu_centers[ifu_ind, 1][0])
                                                   - args.y_offset)
        test = TP(ra_center, dec_center, pa, 1., 1.)
        r, d = test.xy2raDec(+3600. * (ifu_centers[ifu_ind, 0][0])
                             + args.x_offset,
                             +3600. * (ifu_centers[ifu_ind, 1][0])
                             + args.y_offset)
        ra_center += ra - r
        dec_center += dec - d
    ra_center = numpy.max([ra_center, 0.])
    if track == 2:
        # create the data array
        data = numpy.array([[1, ra_center, dec_center, obstime, ill, trans,
                             seeing, 0],
                            [2, ra_center, dec_center, obstime, ill, trans,
                             seeing, 1]])
    else:
        # create the data array
        data = numpy.array([1, ra_center, dec_center, obstime, ill, trans,
                            seeing, track])
        # reshape to 2-d array (1,data.size)
        data = data.reshape(1, data.size)
    shuffle.do(args, data, radius=radius, start=1, orig_loc=[ra, dec],
               targetID=args.target)


if __name__ == '__main__':
    main()
