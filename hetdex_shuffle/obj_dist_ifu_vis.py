#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse as ap
import itertools as it
import logging
import os

from astropy.io import fits
import numpy

from . import visualize
from . import shuffle
from . import __full_version__
from .do_shuffle_target import load_config, setup_logging


failure_code_dict = {-1: "No IFUs within close proximity."}


def writeHeader(f):
    """Write the header to file ``f``

    Parameters
    ----------
    f : file-like object
        where to write to; must have a ``write`` method
    """
    s = []
    s.append("#RA       \t shot RA[deg]")
    s.append("#DEC      \t shot DEC[deg]")
    s.append("#ID       \t shot id for closest VIRUS IFU")
    s.append("#VID      \t virus IFU id")
    s.append("#DRA_VIRUS  \t delta RA for closest VIRUS IFU [arcsec]")
    s.append("#DDEC_VIRUS \t delta Dec for closest VIRUS IFU [arcsec]")
    s.append("#ID       \t shot id for closest LRS IFU")
    s.append("#LID        \t lrs IFU id")
    s.append("#DRA_LRS    \t delta RA for closest LRS IFU [arcsec]")
    s.append("#DDEC_LRS   \t delta Dec for closest LRS IFU [arcsec]")
    f.write('\n'.join(s) + "\n")


def writeShuffle(f, ra, dec, sid, vid, drav, ddecv, sidl, lid, dral,
                 ddecl, failure_code):
    """Write something to file ``f``

    Parameters
    ----------
    f : file-like object
        where to write to; must have a ``write`` method
    """
    s = (" %12.6f %12.6f %5d %3d %12.2f %12.2f %5d %3d %12.2f %12.2f %6d" %
         (ra, dec, sid, int(vid[0]), drav, ddecv, sidl, int(lid[0]), dral,
          ddecl, failure_code))
    f.write(s)
    f.write(" \"%s\"" % "" + "\n")
    f.flush()


def writeFailed(f, ra, dec, failure_code):
    """Write failures to file ``f``

    Parameters
    ----------
    f : file-like object
        where to write to; must have a ``write`` method
    """
    s = (" %12.6f %12.6f %5s %3s %12s %12s %5s %3s %12s %12s %6d" %
         (ra, dec, "-", "-", "-", "-", "-", "-", "-", "-", failure_code))
    f.write(s)
    f.write(" \"%s\"" % failure_code_dict[failure_code] + "\n")

    f.flush()


def get_ifu_list_with_lrs(config):
    """
    Get the ifu centers and ids of the VIRUS+LRS IFUs for the adjusted fplane
    file.
    """
    fplane_adj_file = config.get("General", "fplane_file")
    fplane = shuffle.fplane_parser.FPlane(fplane_adj_file)
    ifu_centers_lrs = numpy.array([[ifu.x, ifu.y] for ifu in fplane.ifus])
    ifu_centers_lrs = ifu_centers_lrs / 3600.
    ifu_id_lrs = numpy.array([[ifu.ifuid] for ifu in fplane.ifus])
    return ifu_centers_lrs, ifu_id_lrs


def parse(argv=None):
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
    parser = ap.ArgumentParser(description=description,
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--verbose", '-v', action="count", default=0,
                        help="""Increase verbosity, can be called multiple
                        times""")
    parser.add_argument('--version', '-V', action='version',
                        version=__full_version__)

    parser.add_argument("-c", "--config", help="""Name of the configuration
                        file. When parsing the command line, the file is loaded
                        into a configuration object""",
                        default="./shuffle.cfg", type=load_config)

    parser.add_argument('fitsfile', help='''name of the file for objects of
                        interest''')

    return parser.parse_args(args=argv)


def main(argv=None):
    '''This code is meant to use tools pre-built and produce images of a list
    of objects with respect to the closest shot, as well as the distance of the
    object to the closest IFU and LRS. LRS is now put in the fplane_w_lrs.txt
    file
    '''
    args = parse(argv=argv)
    setup_logging(args)
    config = args.config
    log = logging.getLogger('shuffle')

    outfile = config.get("General", "outfile")

    # load the data from the outfile with specific interest in RA, Dec, and PA
    data = numpy.loadtxt(outfile, comments='#')

    # load the data from the file of objs (first two columns should be RA and
    # DEC)
    objfile = fits.open(args.fitsfile)
    objdata = numpy.zeros((len(objfile[1].data['ra']), 2))
    objdata[:, 0] = objfile[1].data['ra']
    objdata[:, 1] = objfile[1].data['dec']
    # load ifu centers for VIRUS only
    # ifu_centers = shuffle.load_fplane_file(config)

    # load ifu centers and ids for VIRUS+lRS
    ifu_centers_lrs, ifu_id_lrs = get_ifu_list_with_lrs(config)
    # data=numpy.matrix(data)
    # objdata=numpy.matrix(objdata)
    # spherical distance requires radians for calculation; defining shuffle
    # shot position angles outside of loop
    phi1 = data[:, 0] * numpy.pi / 180
    delta1 = data[:, 1] * numpy.pi / 180

    outcatalog = config.get("General", "closestcatalog")
    log.info("Catalog of closest VIRUS and LRS IFUs: %s", outcatalog)
    with open(outcatalog, 'a') as f_closest:
        writeHeader(f_closest)
        # running through objects in the input file to find closest shuffle
        # position and later IFU
        for j in range(len(objdata[:, 0])):
            failure_code = 0
            ra = objdata[j, 0]
            dec = objdata[j, 1]
            phi2 = ra * numpy.pi / 180.
            delta2 = dec * numpy.pi / 180.

            # calculated distance in radians (spherical distance formula)
            dist_rad = (numpy.arccos(numpy.sin(phi1) * numpy.sin(phi2) +
                        numpy.cos(phi1) * numpy.cos(phi2) *
                        numpy.cos(delta1-delta2)))  # in Radians
            dist_deg = dist_rad * 180. / numpy.pi  # convert to degrees

            # minimum distance required for a possible match is 1.5x diameter
            # of focal plane
            mindist = config.getfloat("General", "dfplane") * 0.55
            # find the shuffles that are with mdist of our object
            loc = [i for (i, val) in enumerate(dist_deg) if val < mindist]

            if len(loc) > 0:
                # if there is a match, then we will make an image and find
                # closest IFU
                # every IFU's center RA value in radians near object
                # every IFU's center Dec value in radians near object
                ifu_ra, ifu_dec = [], []
                for k, (r, d) in it.product(loc, ifu_centers_lrs[:, 0:2]):
                    ifu_ra.append((data[k, 0] +
                                   (numpy.cos(-1.*numpy.pi/180.*data[k, 4]) *
                                    r -
                                    numpy.sin(-1.*numpy.pi/180.*data[k, 4]) *
                                    d) /
                                   numpy.cos(dec*numpy.pi/180.))*numpy.pi/180.)

                ifu_dec.append((data[k, 1] +
                                numpy.sin(-1.*numpy.pi/180. *
                                          data[k, 4])*r +
                                numpy.cos(-1.*numpy.pi/180. *
                                          data[k, 4])*d) *
                               numpy.pi/180.)

                ifu_ra = numpy.array(ifu_ra)
                ifu_dec = numpy.array(ifu_dec)

                # track shot id
                shot_id = numpy.array([k for k in loc
                                       for d in ifu_centers_lrs[:, 1]])
                # track IFU id
                ifu_ids = numpy.array([i for k in loc for i in ifu_id_lrs])

                # hardcoded for LRS to be 079 and 080 (Note: this is not the a
                # flexible approach and could lead to problems as files change)
                sel_lrs = numpy.logical_or(ifu_ids == '079', ifu_ids == '080')
                sel_virus = ~sel_lrs
                sel_lrs.shape = (len(sel_lrs),)
                sel_virus.shape = (len(sel_virus),)
                # distance in arcsecs
                dist_lrs = (180. / numpy.pi * 3600. *
                            numpy.arccos(numpy.sin(ifu_ra[sel_lrs]) *
                                         numpy.sin(phi2) +
                                         numpy.cos(ifu_ra[sel_lrs]) *
                                         numpy.cos(phi2) *
                                         numpy.cos(numpy.abs(ifu_dec[sel_lrs]
                                                             - delta2))))
                # distance in arcsecs
                dist_virus = (180. / numpy.pi * 3600. *
                              numpy.arccos(numpy.sin(ifu_ra[sel_virus]) *
                                           numpy.sin(phi2) +
                                           numpy.cos(ifu_ra[sel_virus]) *
                                           numpy.cos(phi2) *
                                           numpy.cos(numpy.abs(ifu_dec[sel_virus]
                                                               - delta2))))
                # get closest VIRUS IFU center
                virus_closest = numpy.argmin(dist_virus)
                # get closest LRS IFU center
                lrs_closest = numpy.argmin(dist_lrs)

                ra_vir = ifu_ra[sel_virus]
                ra_lrs = ifu_ra[sel_lrs]
                dec_vir = ifu_dec[sel_virus]
                dec_lrs = ifu_dec[sel_lrs]
                # RA offset from closest VIRUS IFU (")
                dra_vir = (-1. * numpy.cos(delta2) *
                           (ra_vir[virus_closest]*180./numpy.pi-ra) * 3600.)
                # RA offset from closest LRS IFU (")
                dra_lrs = (-1. * numpy.cos(delta2) *
                           (ra_lrs[lrs_closest]*180./numpy.pi-ra) * 3600.)
                # Dec offset from closest VIRUS IFU (")
                ddec_vir = (-1. * (dec_vir[virus_closest]*180./numpy.pi-dec) *
                            3600.)
                # Dec offset from closest LRS IFU (")
                ddec_lrs = (-1. * (dec_lrs[lrs_closest]*180./numpy.pi - dec) *
                            3600.)

                sid = shot_id[sel_virus]
                sidl = shot_id[sel_lrs]
                vid = ifu_ids[sel_virus]
                lid = ifu_ids[sel_lrs]

                writeShuffle(f_closest, ra, dec, sid[virus_closest],
                             vid[virus_closest], dra_vir, ddec_vir,
                             sidl[lrs_closest], lid[lrs_closest], dra_lrs,
                             ddec_lrs, failure_code)

                # Get Visualization for Object in the focal plane
                loc1 = sid[virus_closest]
                image_name = os.path.join(config.get('directories', 'images'),
                                          'objimage%04d_%.3f%+.3f.png' %
                                          (j, ra, dec))
                # Shuffle RA,  Dec,  PA
                s = [data[loc1, 0], data[loc1, 1], data[loc1, 4]]
                # Object's RA and Dec
                o = [ra, dec]
                visualize.visualize_nearby_gal(s, o, ifu_centers_lrs,
                                               ifu_id_lrs, config, image_name)
            else:
                failure_code = -1
                writeFailed(f_closest, ra, dec, failure_code)


if __name__ == '__main__':
    main()
