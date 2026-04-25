from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import os
import math
try:
    from StringIO import StringIO  # for Python2
except ImportError:
    from io import BytesIO as StringIO  # python3
import json

from astropy.io import fits
import astropy.units as aunit
from astroquery.skyview import SkyView
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import (RegularPolygon, Circle, Rectangle, FancyArrow, Polygon,
                                Arc, Arrow)
from matplotlib.lines import Line2D
from reproject import reproject_interp
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
from astropy.wcs import WCS
from astropy.stats import biweight_location, biweight_midvariance, mad_std
from matplotlib.collections import PatchCollection
import numpy
import six
from six.moves import configparser
import re
from distutils.dir_util import mkpath
from .astrometry import TangentPlane as TP
import os.path as op

try:
    from urllib import urlencode
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen
    from urllib.parse import urlencode

from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO

from . import findStars

matplotlib.rcParams['text.usetex'] = False


class DS9region(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []

        s.append('global color=green dashlist=8 3 width=1 '
                 'font="helvetica 10 normal roman" select=1 '
                 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 '
                 'include=1 source=1')
        s.append('physical')
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeRegion(cls, f, x, y, radius):
        s = ("circle(%0.2f,%0.2f,%0.2f)" % (x, y, radius))
        f.write(s+"\n")
        f.flush()


try:
    # Python 2 compatibility (safe in Py3 as NameError is caught)
    input = raw_input  # noqa: F821
except NameError:
    pass


def deg2pix(degree, scale=1.698):
    return degree*3600./scale


def wcs2pix(ra, dec, ra0, dec0, scale=1.698, im_size=848, CD=None):
    if CD is not None:
        pixvec = CD.I*numpy.matrix([[(ra-ra0)*numpy.cos(numpy.deg2rad(dec))],
                                    [dec-dec0]])
        x = pixvec[0, 0] + im_size/2.
        y = pixvec[1, 0] + im_size/2.
    else:
        x = -deg2pix(ra-ra0, scale)*numpy.cos(numpy.deg2rad(dec)) + im_size/2.
        y = deg2pix(dec-dec0, scale) + im_size/2.

    return x, y


def deltadeg_to_pix(dra, ddec, dec, scale, size):
    '''Convert the input delta ra/dec around a certain dec into pixels (e.g. in
    an image) given a scale and a size

    Parameters
    ----------
    dra, ddec : floats
        delta ra and dec
    dec : float
        dec at which the deltas are
    scale : float
        pixel scale
    size : float
        size of, e.g., the image

    Returns
    -------
    xc, yc : floats
        pixels positions of the input ``dra/ddec`` in the image frame
    '''
    xc = size/2. - deg2pix(dra, scale)*numpy.cos(numpy.deg2rad(dec))
    yc = size/2. + deg2pix(ddec, scale)
    return xc, yc


def label(xy, text):
    plt.text(xy[0], xy[1], text, ha="center", va="center", family='sans-serif',
             size=4, color='magenta')


def SDSS_coverage(ra, dec):

    log = logging.getLogger('shuffle')

    url_sdssCoverage = 'http://www.sdss3.org/dr9/index.php'
    request_sdssCoverage = urlencode({'coverageRA': ra,
                                      'coverageDec': dec})

    try:
        for line in urlopen(url_sdssCoverage, request_sdssCoverage.encode()):
            if b'overlaps with the SDSS DR9 survey area.' in line:
                return True
    except IOError:
        log.error("Could not connect to %s", url_sdssCoverage)
        return False

    return False


def retrieve_image(ra, dec, size, config, yflip, scale=1.5):
    """Wrapper function for retrieving image from SDSS. If region outside SDSS
    coverage, it uses DSS image instead.
    (ra, dec, size) in degree
    """
    return retrieve_image_hips(ra, dec, size, yflip, scale)


def return_blank_image(size_pix, scale):
    imarray = numpy.zeros((size_pix, size_pix))
    CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
    return imarray, CD, '', ''


def retrieve_image_PS1(ra, dec, size, yflip, scale):
    """
    Load an image from the PS1 image cut out service. Choose
    a scale that is a multiple of 0.25. Bin image up to match
    the input scale.
    """

    log = logging.getLogger("shuffle")
    url_ps1_fn = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url_ps1_fitscut = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    request = urlencode({'ra': ra, 'dec': dec, 'filters': 'gri'})

    # image params
    scale_factor = int(scale/0.25)
    size_pix = int(size*3600./scale)*scale_factor

    if scale < 0.25:
        raise ValueError("Please choose a scale thats a multiple of 0.25")

    # first get the names of the images we need
    try:
        imlist = urlopen(url_ps1_fn, request.encode())
    except IOError:
        log.warning("HTTP ERROR request of: %s %s", url_ps1_fn, request)
        return return_blank_image(size_pix/scale_factor, scale)

    # produce a dictionary of filters and filenames
    imdict = {}
    next(imlist)
    for line in imlist:
        els = line.split()
        if len(els) < 7:
            log.warning("Query failed!")
            return return_blank_image(size_pix/scale_factor, scale)
        imdict[els[4]] = els[7]

    # make the image cutout
    try:
        # just one band as otherwise too much noise to see overplot
        opts = {'Green': imdict['r'], 'Red': imdict['i'], 'Blue': imdict['g'],
                'Size': size_pix, 'Ra': ra, 'Dec': dec, 'AutoscaleMax': 99.5,
                'AutoscaleMin': 50.0}
    except KeyError:
        log.warning("Required PS1 images not found!")
        return return_blank_image(size_pix/scale_factor, scale)

    request_cutout = urlencode(opts)

    # get image
    try:
        imfile = urlopen(url_ps1_fitscut, request_cutout.encode())
    except IOError:
        log.warning("HTTP ERROR request of: %s %s", url_ps1_fitscut,
                    request_cutout)
        return return_blank_image(size_pix/scale_factor, scale)

    try:
        imarray = plt.imread(StringIO(imfile.read()), format='jpeg')
    except IOError:
        log.warning("Image aquisition failed!")
        return return_blank_image(size_pix/scale_factor, scale)

    if yflip:
        imarray = imarray[::-1, :]

    # bin it up
    rgb_arr = []
    nx, ny, nc = imarray.shape
    for c in range(nc):
        channel = imarray[:, :, c]
        nx, ny = channel.shape
        # nsize = int(nx/scale_factor)
        im_view = channel.reshape(int(nx/scale_factor), scale_factor,
                                  int(ny/scale_factor), scale_factor)
        binned = im_view.mean(axis=3).mean(axis=1)
        rgb_arr.append(numpy.array(255*(1.0 -
                                        binned/numpy.max(numpy.array(binned))),
                       dtype=int))

    # stick binned up channels into single RGB image
    rgb_binned = numpy.zeros((rgb_arr[0].shape[0], rgb_arr[0].shape[1], 3),
                             dtype=int)
    rgb_binned[..., 0] = rgb_arr[0]
    rgb_binned[..., 1] = rgb_arr[1]
    rgb_binned[..., 2] = rgb_arr[2]

    # Get WCS
    opts['getWCS'] = 'true'
    request_wcs = urlencode(opts)
    try:
        wcsr = urlopen(url_ps1_fitscut, request_wcs.encode())
    except IOError:
        log.warning("HTTP ERROR request of: %s %s", url_ps1_fitscut,
                    request_wcs)
        return return_blank_image(size_pix/scale_factor, scale)

    try:
        wcs_dict = json.load(wcsr)
        CD = numpy.matrix(scale_factor *
                          numpy.array(wcs_dict['cdmatrix']).reshape(2, 2))
    except KeyError:
        log.warning("Could not get WCS info!")
        return return_blank_image(size_pix/scale_factor, scale)

    return rgb_binned, CD, url_ps1_fitscut+'?'+request_cutout, 'PS1'


def retrieve_image_PANSTARRS(ra, dec, size, yflip, scale):

    def getimages(ra, dec, size=240, filters="grizy"):

        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
               "&filters={filters}").format(**locals())
        table = Table.read(url, format='ascii')
        return table

    def geturl(ra, dec, size=240, output_size=None, filters="grizy",
               format="jpg", color=False):

        """Get URL for images in the table

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        color = if True, creates a color image (only for jpg or png format).
                Default is return a list of URLs for single-filter
                grayscale images.
        Returns a string with the URL
        """

        if color and format == "fits":
            raise ValueError("color images are available only "
                             "for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = getimages(ra, dec, size=size, filters=filters)
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[numpy.argsort(flist)]
        if color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0, len(table)//2, len(table)-1]]
            for i, param in enumerate(["red", "green", "blue"]):
                url = url + "&{}={}".format(param, table['filename'][i])
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase+filename)
        return url

    def getgrayim(ra, dec, size=240, output_size=None, filter="g",
                  format="jpg"):

        """Get grayscale image at a sky position

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filter = string with filter to extract (one of grizy)
        format = data format (options are "jpg", "png")
        Returns the image
        """

        if format not in ("jpg", "png"):
            raise ValueError("format must be jpg or png")
        if filter not in list("grizy"):
            raise ValueError("filter must be one of grizy")
        url = geturl(ra, dec, size=size, filters=filter,
                     output_size=output_size, format=format)
        r = requests.get(url[0])
        im = Image.open(BytesIO(r.content))
        return im
    size_pix = int(size*3600./scale)
    gim = getgrayim(ra, dec, size=size_pix, filter="g")
    imarray = numpy.array(gim)
    if yflip:
        imarray = imarray[::-1, :]
    CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
    return imarray, CD, '', 'PS1'


def retrieve_image_SDSS(ra, dec, size, config, yflip, scale=1.5):
    """
    Load image from sdss-dr9 or dss server(jpeg) and return the image array
    and the url. Note that the transformation from world coordinate(ra,dec) to
    pixel position(x,y) is simple projection without rotation, i.e.
    x=-scale*(ra-ra0)*cos(dec)+x0; y=scale*(dec-dec0)+y0
    """
    log = logging.getLogger('shuffle')

    url_sdssJPEG = 'http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx'
    size_pix = int(size*3600./scale)
    # options for the finding chart (see
    # http://skyserver.sdss3.org/dr9/en/tools/chart/chart.asp )
    opt = config.get("General", "SDSS_FC_options")

    if size_pix > 2048:
        size_pix = 2048
        scale = size*3600./size_pix

    request_sdss = urlencode({'ra': ra, 'dec': dec, 'scale': scale,
                              'height': size_pix, 'width': size_pix,
                              'opt': opt})
    # url =
    # "http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx?"+request_sdss

    try:
        imfile = urlopen(url_sdssJPEG, request_sdss.encode())
    except IOError:
        log.warning("HTTP ERROR request of: %s %s", url_sdssJPEG, request_sdss)
        raise IOError
        #return -99, -99, 'bogus', 'bogus'
    imarray = plt.imread(StringIO(imfile.read()), format='jpeg')

    if yflip:
        imarray = imarray[::-1, :]
    CD = numpy.matrix([[-1.*scale/3600., 0], [0, 1.*scale/3600.]])
    return imarray, CD, url_sdssJPEG+'?'+request_sdss, 'SDSS'


def retrieve_image_SkyView(ra, dec, size, yflip, scale):
    size_pix = int(size*3600./scale)

    paths = SkyView.get_images(position='%3.5f %2.5f' % (ra, dec),
                               survey='DSS2 Blue', coordinates='J2000',
                               height=size*aunit.deg, width=size*aunit.deg,
                               pixels=str(size_pix), cache=False)
    print(paths)
    if paths:
        hdu = paths[0][0].header
        imarray = paths[0][0].data
        CD = numpy.matrix([[hdu['CDELT1'], 0], [0, hdu['CDELT2']]])
    else:
        # Default to PS1
        log = logging.getLogger('shuffle')
        log.warning("Could not contact SDSS or DSS. Falling back to PS1!")
        return retrieve_image_PS1(ra, dec, size, yflip, scale)

    return imarray, CD, '', 'DSS2'


def retrieve_image_DSS(ra, dec, size, yflip):
    """Load image from dss server(fits) and return the image array, the url,
    and the CD matrix (there is rotation in DSS images).
    CD matrix transforms the pixel position(x,y) to world coordinate(ra,dec).
    """
    url_dss = 'http://archive.eso.org/dss/dss/image'
    # request_dss = urlencode({'ra':ra, 'dec':dec, 'x':size*60,
    # 'y':size*60, 'mime-type':'download-gif'})
    request_dss = urlencode({'ra': ra, 'dec': dec, 'x': size*60,
                             'y': size*60,
                             'mime-type': 'download-fits'})
    # url = "http://archive.eso.org/dss/dss/image?"+request_dss
    # imfile = urlopen(url_dss, request_dss.encode())
    # imarray = plt.imread(StringIO(imfile.read()), format='gif')
    # print(url_dss+'?'+request_dss_fits)
    hdulist = fits.open(url_dss+'?'+request_dss)
    hdu = hdulist[0].header
    imarray = hdulist[0].data
    if yflip:
            imarray = imarray[::-1, :]
    CD = numpy.matrix([[hdu['CD1_1'], hdu['CD1_2']], [hdu['CD2_1'],
                       hdu['CD2_2']]])
    return imarray, CD, url_dss+'?'+request_dss, 'DSS'


def plotFocalPlane(dra, ddec, pa, scale, ifu_centers, ifu_ids, dec, im_size,
                   config, color='green', linewidth=0.2):
    """Plot the region of IFUs and patrol circle and return as a
    PatchCollection, which can be added to axes by axes.add_collection.
    """

    ifu_size = config.getfloat("General", "ifu_size")
    ifu_edge = config.getfloat("General", "ifu_edge")

    patches = []
    xc, yc = deltadeg_to_pix(dra, ddec, dec, scale, im_size)
    rpa = numpy.deg2rad(pa + config.getfloat('offsets', 'fplane'))

    # plot all IFU regions
    for id_, c in zip(ifu_ids, ifu_centers):
        x, y = c*3600./scale
        x *= -1.
        xr = numpy.cos(rpa)*x - numpy.sin(rpa)*y + xc
        yr = numpy.sin(rpa)*x + numpy.cos(rpa)*y + yc
        # still need to correct the xr?
        patches.append(RegularPolygon((xr, yr), 4,
                                      radius=deg2pix(ifu_size,
                                                     scale)/numpy.sqrt(2.),
                                      orientation=rpa-numpy.pi/4.))
        patches.append(RegularPolygon((xr, yr), 4,
                                      radius=deg2pix(ifu_size-ifu_edge,
                                                     scale)/numpy.sqrt(2.),
                                      orientation=rpa-numpy.pi/4.))

    # patrol field and field center
    patches.extend(patrol_and_zenit(config, dra, ddec, dec, pa, scale, im_size,
                                    color=color, linewidth=linewidth,
                                    return_patches=True))

    return PatchCollection(patches, edgecolor=color, facecolor='none',
                           linewidth=linewidth)


def plotFocalPlaneLRS(dra, ddec, pa, scale, ifu_centers, ifu_ids, dec, im_size,
                      config, color='green', linewidth=1.0, text_ax=None,
                      text_rot=None):
    """Plot the region of IFUs and patrol circle and return as a
    PatchCollection, which can be added to axes by axes.add_collection.

    Parameters
    ----------
    dra, ddec : floats
        delta ra and dec
    pa : float
        PA in degrees
    scale : float
        pixel scale
    ifu_centers, ifu_ids : lists
        list of centers and id (which one?) for the IFUS from the fplane file
    dec : float
        dec at which the deltas are
    im_size : float
        size of the image
    config : :class:`configparser.ConfigParser` instance
        configuration object
    color : any matplotlib color, optional
        color to use for the patches
    linewidth : float, optional
        with of the lines used to plot the parts
    text_ax : None, "none" or a :class:`matplotlib.axes` instance
        if None: add the text with the id of the IFU using ``plt.text``
        if 'none': no text added
        if a :class:`matplotlib.axes` instance: adds the text using
        ``text_ax.text``
    text_rot : float
        rotate the text by the give amount
    """
    # log = logging.getLogger('shuffle')

    lrs_sizex = config.getfloat("General", "lrs_sizex")
    lrs_sizey = config.getfloat("General", "lrs_sizey")
    lrs_sizex = deg2pix(lrs_sizex, scale)
    lrs_sizey = deg2pix(lrs_sizey, scale)

    patches = []
    xc, yc = deltadeg_to_pix(dra, ddec, dec, scale, im_size)
    rpa = numpy.deg2rad(pa + config.getfloat('offsets', 'fplane'))

    if text_ax is None:
        text = plt.text
    elif text_ax == 'none':  # do nothing
        text = None
    else:
        text = text_ax.text

    lrs_ifus = [(id_, cen) for id_, cen in zip(ifu_ids, ifu_centers)
                if id_ in ['056', '066', '600', '603']]
    labels = ['LRS2B', 'LRS2R', 'HPFACAM', 'HPFSCI']
    for thing, lab in zip(lrs_ifus, labels):
        id_, cen = thing
        x, y = cen * 3600. / scale
        # x1 = x + deg2pix(lrs_sizex, scale)/2.
        # y1 = y - deg2pix(lrs_sizey, scale)/2.
        x = -x  # RA and x run in opposite directions
        # get the LRS lower left vertex
        if 'LRS2' in lab:
            x_vertex, y_vertex = x - lrs_sizex/2., y - lrs_sizey/2.
        elif 'HPFACAM' in lab:
            x_vertex, y_vertex = x - 15./2., y - 15//2.
        else:
             x_vertex, y_vertex = x, y
        # get the coordinates of the vertex of the LRS in the rotated reference
        # system
        x_rect = numpy.cos(rpa)*x_vertex - numpy.sin(rpa)*y_vertex + xc
        y_rect = numpy.sin(rpa)*x_vertex + numpy.cos(rpa)*y_vertex + yc
        # still need to correct the xr?
        dpalrs2 = config.getfloat('offsets', 'fplane')
        if 'LRS2' in lab:
            patches.append(Rectangle((x_rect, y_rect),
                           width=lrs_sizex, height=lrs_sizey,
                           angle=pa+dpalrs2))
        elif 'HPFACAM' in lab:
            patches.append(Rectangle((x_rect, y_rect),
                           width=15., height=15.,
                           angle=pa+dpalrs2))
        else: 
            patches.append(Circle((x_rect, y_rect),
                           5.))
            
        if text:
            # x and y coordinates of the text
            xl = x_rect - 0.*lrs_sizex/2.
            yl = y_rect - 5.*lrs_sizey/8.
            text(xl, yl, lab, ha="center", va="center", family='sans-serif',
                 size=5, color='magenta', rotation=text_rot)

        # get the coordinates of the center of the LRS in the rotated reference
        # system
        xr = numpy.cos(rpa)*x - numpy.sin(rpa)*y + xc
        yr = numpy.sin(rpa)*x + numpy.cos(rpa)*y + yc
        rpaflip = rpa + numpy.pi / 2.
        xln = numpy.cos(rpaflip)*(599.-xr) - numpy.sin(rpaflip)*(yr-599.) + 599
        yln = numpy.sin(rpaflip)*(599.-xr) + numpy.cos(rpaflip)*(yr-599.) + 599
        # xn = (xln-599.)*0.25/0.271+388.
        # yn = (yln-599.)*0.25/0.271+386.
        # LRS2B should be at 710,183; LRS2R at 729,550
        # log.info("LRS2 centered at: %s, %s, %s", id_, xn, yn)

    return PatchCollection(patches, edgecolor=color, facecolor='none',
                           linewidth=linewidth)


def patrol_and_zenit(config, dra, ddec, dec, pa, scale, size, color='green',
                     linewidth=1., return_patches=False):
    '''Create patches with the patrol regions, the center of the field and the
    zenit

    Parameters
    ----------
    config : :class:`configparser.ConfigParser` instance
        configuration object
    dra, ddec : floats
        delta ra and dec
    dec : float
        dec at which the deltas are
    pa : float
        PA in degrees
    scale : float
        pixel scale
    size : float
        size of, e.g., the image
    color : any matplotlib color, optional
        color to use for the patches
    linewidth : float, optional
        with of the lines used to plot the parts
    return_patches : bool, optional
        if True, returns the list of patches, otherwise returns a
        :class:`matplotlib.collections.PatchCollection` instance

    Returns
    -------
    :class:`matplotlib.collections.PatchCollection` instance or list of patches
        collection of circles
    '''
    # patrol field and field center
    xc, yc = deltadeg_to_pix(dra, ddec, dec, scale, size)
    rpa = numpy.deg2rad(pa)

    patches = []
    rpatrol_min = config.getfloat("General", "dpatrol_min")/2.
    rpatrol_max = config.getfloat("General", "dpatrol_max")/2.
    patches.append(Circle((xc, yc), deg2pix(rpatrol_min, scale)))
    patches.append(Circle((xc, yc), deg2pix(rpatrol_max, scale)))
    patches.append(Circle((xc, yc), deg2pix(0.001, scale)))

    # indication of zenith
    x, y = 0., deg2pix(0.17,  scale)
    xr = numpy.cos(rpa)*x - numpy.sin(rpa)*y + xc
    yr = numpy.sin(rpa)*x + numpy.cos(rpa)*y + yc
    patches.append(Circle((xr, yr), deg2pix(0.01, scale)))

    if return_patches:
        return patches
    else:
        return PatchCollection(patches, edgecolor=color, facecolor='none',
                               linewidth=linewidth)


def pa_arrow(ax, x, y, pa, rarrow, scale, va='bottom', text_rot=0):
    '''Add to the give axes the arrow and the text for the PA arrow

    Parameters
    ----------
    ax : :class:`matplotlib.axes` instance
        axes to which add the arrow and the text
    x, y : floats
        starting point of the arrow
    pa : float
        PA angle
    rarrow : float
        length of the arrow
    scale : float
        pixel scale of the image
    va : string, optional
        vertical alignment of the text added near the arrow
    text_rot : float, optional
        rotation of the text; if None
    '''
    rpa = numpy.deg2rad(pa + 90)

    x_arrow = x+deg2pix(6./3600., scale)*(numpy.cos(rpa))
    y_arrow = y+deg2pix(6./3600., scale)*(numpy.sin(rpa))
    dx_arrow, dy_arrow = rarrow*numpy.cos(rpa), rarrow*numpy.sin(rpa)
    arrow = FancyArrow(x_arrow, y_arrow, dx_arrow, dy_arrow, width=0.5,
                       linewidth=0.5, edgecolor=[0.8, 0.8, 0.8],
                       facecolor=[1.0, 1.0, 1.0])
    ax.add_patch(arrow)
    # write the PA angle
    ax.text(x_arrow + dx_arrow / 2., y_arrow + dy_arrow / 2.,
            'PA:%.1f\u00b0' % (pa), ha="center", va=va,
            family='sans-serif', size=3, color='white', rotation=text_rot)


def plotFocalPlaneDS9(ds9, ra, dec, pa, ifu_centers, config, color='green',
                      ifu_ids=None, usename=False):

    def circle(ds9, ra, dec, r=0.183, color='red', nsub=30, start=0,
               stop=2*numpy.pi, width=1, name=None):

        decrad = dec * numpy.pi/180.0
        if stop < start:
            start -= 2*numpy.pi
        log = logging.getLogger('shuffle')
        log.debug(start, stop)
        angle = numpy.linspace(start, stop, nsub)
        x = -(r*numpy.cos(angle)/numpy.cos(decrad)) + ra
        y = r*numpy.sin(angle) + dec
        # Do this once here to add a title to the circle
        if name:
            ds9.set('region', "wcs; line %f %f %f %f # color='%s' "
                    "select=0 width=%i text=%s"
                    % (x[0], y[0], x[1], y[1], color, width, name))
        else:
            ds9.set('region', "wcs; line %f %f %f %f # color='%s' "
                    "select=0 width=%i"
                    % (x[0], y[0], x[1], y[1], color, width))

        for x1, x2, y1, y2 in zip(x[1:-1], x[2:], y[1:-1], y[2:]):
            ds9.set('region', "wcs; line %f %f %f %f # color='%s' "
                    "select=0 width=%i"
                    % (x1, y1, x2, y2, color, width))

    def line(ds9, ra, dec, r1, r2, angle, color='red', width=1):
        decrad = dec * numpy.pi/180.0

        x1 = -(r1*numpy.cos(angle)/numpy.cos(decrad)) + ra
        y1 = r1*numpy.sin(angle) + dec
        x2 = -(r2*numpy.cos(angle)/numpy.cos(decrad)) + ra
        y2 = r2*numpy.sin(angle) + dec
        ds9.set('region', "wcs; line %f %f %f %f # color='%s' "
                          "select=0 width=%i" % (x1, y1, x2, y2, color, width))

    def ifu(ds9, x, y, pa, size=0.014, color='green', ifu_imhpid=None):
        if ifu_imhpid is not None:
            ds9.set('region',
                    "wcs; box %f %f %f %f %f # color='%s'"
                    " text={%s}" % (x, y, size, size, pa, color,
                                    "{:03d}".format(int(ifu_imhpid))))
        else:
            ds9.set('region',
                    "wcs; box %f %f %f %f %f #"
                    " color='%s'" % (x, y, size, size, pa, color))

    ifu_size = config.getfloat("General", "ifu_size")
    ifu_edge = config.getfloat("General", "ifu_edge")

    rpa = numpy.deg2rad(pa)
    if usename:
        i = 0
        for c in ifu_centers:
            x, y = c
            xr = numpy.cos(rpa) * x + numpy.sin(rpa) * y
            yr = - numpy.sin(rpa) * x + numpy.cos(rpa) * y
            xl, yl = xr/numpy.cos(numpy.deg2rad(yr+dec)), yr
            if "{:03d}".format(int(ifu_ids[i])) != "000":
                ifu(ds9, xl+ra, yl+dec, pa, size=ifu_size, color=color,
                    ifu_imhpid=ifu_ids[i])
                ifu(ds9, xl+ra, yl+dec, pa, size=ifu_size+ifu_edge,
                    color=color)
            i += 1
    else:
        i = 0
        for c in ifu_centers:
            x, y = c
            xr = numpy.cos(rpa) * x + numpy.sin(rpa) * y
            yr = - numpy.sin(rpa) * x + numpy.cos(rpa) * y
            xl, yl = xr/numpy.cos(numpy.deg2rad(yr+dec)), yr
            if "{:03d}".format(int(ifu_ids[i])) != "000":
                ifu(ds9, xl+ra, yl+dec, pa, size=ifu_size, color=color)
                ifu(ds9, xl+ra, yl+dec, pa, size=ifu_size+ifu_edge,
                    color=color)
            i += 1

    rpatrol_min = config.getfloat("General", "dpatrol_min")/2.
    rpatrol_max = config.getfloat("General", "dpatrol_max")/2.

    # guide patrol regions
    g1_min = config.getfloat("General", "dpatrol_g1min")
    g1_max = config.getfloat("General", "dpatrol_g1max")
    g2_min = config.getfloat("General", "dpatrol_g2min")
    g2_max = config.getfloat("General", "dpatrol_g2max")
    w1_min = config.getfloat("General", "dpatrol_w1min")
    w1_max = config.getfloat("General", "dpatrol_w1max")
    w2_min = config.getfloat("General", "dpatrol_w2min")
    w2_max = config.getfloat("General", "dpatrol_w2max")

    arc_width = [[1.00, 1.00], [1.00, 1.00], [1.005, 0.995], [1.005, 0.995]]
    colors = ['yellow', 'yellow', 'white', 'white']
    names = ['gd1', 'gd2', 'wfs1', 'wfs2']
    thetas = [[g1_min, g1_max], [g2_min, g2_max],
              [w1_min, w1_max], [w2_min, w2_max]]
    pa_rotation = pa + 90 + config.getfloat('offsets', 'probes')
    for width, theta, col, name in zip(arc_width, thetas, colors, names):
        circle(ds9, ra, dec, r=rpatrol_min*width[0], color=col, width=3,
               start=(theta[0]+pa_rotation)*numpy.pi/180.,
               stop=(theta[1]+pa_rotation)*numpy.pi/180.)
        circle(ds9, ra, dec, r=rpatrol_max*width[1], color=col, width=3,
               start=(theta[0]+pa_rotation)*numpy.pi/180.,
               stop=(theta[1]+pa_rotation)*numpy.pi/180.)
        line(ds9, ra, dec, r1=rpatrol_min*width[0], r2=rpatrol_max*width[1],
             angle=(theta[0]+pa_rotation)*numpy.pi/180., color=col,
             width=2)
        line(ds9, ra, dec, r1=rpatrol_min*width[0], r2=rpatrol_max*width[1],
             angle=(theta[1]+pa_rotation)*numpy.pi/180., color=col,
             width=2)

    # patrol field
    # circle(ds9, ra, dec, r=rpatrol_min, color=color)
    # circle(ds9, ra, dec, r=rpatrol_max, color=color)
    # field center
    circle(ds9, ra, dec, r=0.001, color=color)

    # indication of zenith
    x, y = 0., 0.17
    xr = numpy.cos(rpa) * x + numpy.sin(rpa) * y
    yr = - numpy.sin(rpa) * x + numpy.cos(rpa) * y
    xl, yl = xr/numpy.cos(numpy.deg2rad(yr+dec)), yr
    circle(ds9, xl+ra, yl+dec, r=0.01, color=color)


def visualize(s, r, ifu_centers, ifu_ids, guideWFSSol, cal_star_cand, config,
              filename, ifu_imhpid):
    """
    Plot a suffle solution with matplotlib and save image to [filename].
    r is [ra,dec,pa] before shuffling.
    s is [ra',dec',pa'] after shuffling.
    """
    log = logging.getLogger('shuffle')

    ra, dec, pa = s[0], s[1], s[2]
    log.debug("PA: %f", pa)
    # size of the image in degree
    size = config.getfloat("General", "dfplane")
    size += 1.2*config.getfloat("General", "radius")
    dpi = 300.

    imarray, CD, url, img_src = retrieve_image(ra, dec, size, config,
                                               config.getboolean("General",
                                                                 "yflip"), 1.5)
    if (len(imarray) == 1) and imarray == -99:
        return
    size_pix = len(imarray)
    scale = size*3600./size_pix

    fg = plt.figure(figsize=(size_pix/dpi, size_pix/dpi), dpi=dpi)
    ax = fg.add_axes([0, 0, 1, 1], frameon=False)

    # plot the focal plane of IFUs
    try:
        no_plot_original = config.getboolean('visualisation',
                                             'no_original_ifu_pos')
    except (configparser.NoSectionError, configparser.NoOptionError):
        no_plot_original = False

    if not no_plot_original:
        ax.add_collection(plotFocalPlane(0, 0, pa, scale, ifu_centers, ifu_ids,
                                         dec, size_pix, config, color='grey'))
    ax.add_collection(plotFocalPlane(r[0]-ra, r[1]-dec, pa, scale, ifu_centers,
                                     ifu_ids, dec, size_pix, config,
                                     color='green'))

    # put the central RA and Dec in the center. \u00b0 : degree symbol
    ax.text(size_pix/2., size_pix/2.,
            u'RA:%.3f\u00b0\n\nDec:%.3f\u00b0' % (ra, dec), color='grey',
            size=5/scale, ha='center', va='center')
    # plot the scale and direction indicators for DSS image
    if img_src == 'DSS':
        ax.text(size_pix*1./40.+60./scale, size_pix*9./10.+50./scale,
                'DSS image',  color='white', size=5./scale, ha='center')
        # plot the scale
        ax.add_line(Line2D([size_pix*1./40., size_pix*1./40.+120./scale],
                    [size_pix*9./10., size_pix*9./10.], linewidth=0.3,
                    color='white', marker='|', markersize=2))
        ax.text(size_pix*1./40.+60./scale, size_pix*9./10.+15./scale, "2'",
                color='white', size=5./scale, ha='center')
        # add the direction indicator
        ax.text(size_pix/2., 40./scale, 'S', color='white', size=5./scale,
                ha='center', va='center')
        ax.text(size_pix/2., size_pix-40./scale, 'N', color='white',
                size=5./scale, ha='center', va='center')
        ax.text(40./scale, size_pix/2., 'E', color='white', size=5./scale,
                ha='center', va='center')
        ax.text(size_pix-40./scale, size_pix/2., 'W', color='white',
                size=5./scale, ha='center', va='center')

    # print("Number of suitable IFU calibration in the FOV:",
    # len(cal_star_cand))
    ifu_size = config.getfloat("General", "ifu_size")
    ifu_edge = config.getfloat("General", "ifu_edge")
    insiders = findStars.inside(r[0], r[1], pa, ifu_centers, cal_star_cand,
                                ifu_size=ifu_size, edge=ifu_edge, PLOT=False)

    for o in cal_star_cand[~insiders]:
        ax.add_patch(Circle(wcs2pix(o[2], o[3], ra, dec, scale=scale,
                                    im_size=size_pix, CD=CD),
                            deg2pix(0.0015,  scale), edgecolor='blue',
                            facecolor='none', linewidth=0.2))

    for o in cal_star_cand[insiders]:
        x, y = wcs2pix(o[2], o[3], ra, dec, scale=scale, im_size=size_pix,
                       CD=CD)
        ax.add_patch(Circle((x, y), deg2pix(0.0015,  scale),
                            edgecolor='magenta', facecolor='none',
                            linewidth=0.2))
        ax.text(x, y+10/scale, '%d:%d' % (o[0], o[1]), color='magenta',
                size=5/scale, ha='center')

    # plot the guide star and WFS solution
    n_guide = 0
    colordict = {0: 'yellow', 1: 'yellow', 2: 'white', 3: 'white'}
    textprefix = {0: 'gd1', 1: 'gd2', 2: 'wfs1', 3: 'wfs2'}
    for s in guideWFSSol:
        x, y = wcs2pix(s[2], s[3], ra, dec, scale=scale, im_size=size_pix,
                       CD=CD)
        ax.add_patch(Circle((x, y), deg2pix(0.0017,  scale),
                            edgecolor=colordict[n_guide], facecolor='none',
                            linewidth=0.2))
        ax.text(x, y+10/scale, textprefix[n_guide]+':%d:%d' % (s[0], s[1]),
                color=colordict[n_guide], size=5/scale, ha='center')
        n_guide += 1

    # add the arrow of the PA. This arrow should *not* be corrected for the
    # small offset
    x, y = wcs2pix(r[0], r[1], ra, dec, scale=scale, im_size=size_pix, CD=CD)

    # add the arrow showing the PA
    pa_arrow(ax, x, y, pa, size_pix / 10, scale)

    # plot the arcs of the patrol regions
    rpatrol_min = config.getfloat("General", "dpatrol_min")
    rpatrol_max = config.getfloat("General", "dpatrol_max")

    g1_min = config.getfloat("General", "dpatrol_g1min")
    g1_max = config.getfloat("General", "dpatrol_g1max")
    g2_min = config.getfloat("General", "dpatrol_g2min")
    g2_max = config.getfloat("General", "dpatrol_g2max")
    w1_min = config.getfloat("General", "dpatrol_w1min")
    w1_max = config.getfloat("General", "dpatrol_w1max")
    w2_min = config.getfloat("General", "dpatrol_w2min")
    w2_max = config.getfloat("General", "dpatrol_w2max")

    arc_width = [[1.01, 0.99], [1.01, 0.99], [1.02, 0.98], [1.02, 0.98]]
    colors = ['y', 'y', 'w', 'w']
    thetas = [[g1_min, g1_max], [g2_min, g2_max],
              [w1_min, w1_max], [w2_min, w2_max]]
    pa_rotation = pa + 90 + config.getfloat('offsets', 'probes')

    for width, theta, col in zip(arc_width, thetas, colors):
        ax.add_patch(Arc((x, y), deg2pix(rpatrol_min*width[0], scale),
                     deg2pix(rpatrol_min*width[0],  scale), angle=pa_rotation,
                     theta1=theta[0], theta2=theta[1], color=col))
        ax.add_patch(Arc((x, y), deg2pix(rpatrol_max*width[1], scale),
                     deg2pix(rpatrol_max*width[1],  scale), angle=pa_rotation,
                     theta1=theta[0], theta2=theta[1], color=col))
        x1 = (numpy.cos(numpy.deg2rad(theta[0] + pa_rotation)) *
              deg2pix(rpatrol_min*width[0], scale)/2. + x)
        x2 = (numpy.cos(numpy.deg2rad(theta[0] + pa_rotation)) *
              deg2pix(rpatrol_max*width[1], scale)/2. + x)
        x3 = (numpy.cos(numpy.deg2rad(theta[1] + pa_rotation)) *
              deg2pix(rpatrol_min*width[0], scale)/2. + x)
        x4 = (numpy.cos(numpy.deg2rad(theta[1] + pa_rotation)) *
              deg2pix(rpatrol_max*width[1], scale)/2. + x)
        y1 = (numpy.sin(numpy.deg2rad(theta[0] + pa_rotation)) *
              deg2pix(rpatrol_min*width[0], scale)/2. + y)
        y2 = (numpy.sin(numpy.deg2rad(theta[0] + pa_rotation)) *
              deg2pix(rpatrol_max*width[1], scale)/2. + y)
        y3 = (numpy.sin(numpy.deg2rad(theta[1] + pa_rotation)) *
              deg2pix(rpatrol_min*width[0], scale)/2. + y)
        y4 = (numpy.sin(numpy.deg2rad(theta[1] + pa_rotation)) *
              deg2pix(rpatrol_max*width[1], scale)/2. + y)
        ax.plot([x1, x2], [y1, y2], '-', color=col)
        ax.plot([x3, x4], [y3, y4], '-', color=col)

    ax.imshow(imarray, origin='lower', cmap='gray')
    log.debug('Image retrieved from "%s"', url)

    fg.savefig(filename, dpi=dpi)
    log.debug('Visualization image is saved to %s', filename)

    plt.close(fg)


def visualize_nearby_gal(s, objloc, ifu_centers, ifu_ids, config, filename,
                         ifu_ihmpid):
    """Plot a shuffle solution with matplotlib along side object of interest
    and save image to [filename].
    objloc is [ra,dec]  of the object.
    s is [ra',dec',pa'] of the shuffe.
    """
    ra, dec, pa = s[0], s[1], s[2]
    ra_obj, dec_obj = objloc[0], objloc[1]
    print("PA:", pa)
    # print(ra,dec,ra_obj,dec_obj)
    # The size of the image is the max of 2x size of IFU or 4x Petrosian Radius
    size = 6*config.getfloat("General", "ifu_size")
    dpi = 300.

    imarray, CD, url, img_src = retrieve_image(ra_obj, dec_obj, size, config,
                                               config.getboolean("General",
                                                                 "yflip"), 0.8)
    size_pix = len(imarray)
    scale = size*3600./size_pix

    fg = plt.figure(figsize=(size_pix/dpi, size_pix/dpi), dpi=dpi)
    ax = fg.add_axes([0, 0, 1, 1], frameon=False)

    # plot the focal plane of IFUs
    ax.add_collection(plotFocalPlaneLRS(ra-ra_obj, dec-dec_obj, pa, scale,
                                        ifu_centers, ifu_ids, dec_obj,
                                        size_pix, config, color='green'))
    ax.add_collection(patrol_and_zenit(config, ra-ra_obj, dec-dec_obj, dec_obj,
                                       pa, scale, size_pix))
    ax.imshow(imarray, origin='lower', cmap='gray')

    # ax.text(size_pix*0.05,size_pix*0.95,
    # u'RA:%.3f\u00b0\nDec:%.3f\u00b0'%(ra,dec), color='grey', size=1./scale,
    # ha='left', va='top')

    # ax.add_line(Line2D([size_pix*39./40.,size_pix*39./40.-60./scale],
    # [size_pix*1./10.,size_pix*1./10.], linewidth=0.3, color='white',
    # marker='|', markersize=2))
    # ax.text(size_pix*39./40.-30./scale, size_pix*1./10.+3./scale, '1\'',
    # color='white', size=1./scale, ha='center')
    # add the direction indicator
    # ax.text(size_pix/2.,		4./scale, 		'S',
    # color='white', size=1./scale, ha='center', va='center')
    # ax.text(size_pix/2.,		size_pix-4./scale, 	'N',
    # color='white', size=1./scale, ha='center', va='center')
    # ax.text(4./scale,		size_pix/2., 		'E', color='white',
    # size=1./scale, ha='center', va='center')
    # ax.text(size_pix-4./scale, 	size_pix/2., 		'W',
    # color='white', size=1./scale, ha='center', va='center')

    n_center = 0
    colordict = {0: 'red', 1: 'yellow', 2: 'white', 3: 'white'}
    x, y = wcs2pix(ra, dec, ra, dec, scale=scale, im_size=size_pix, CD=CD)
    ax.add_patch(Circle((x, y), deg2pix(5./3600.,  scale),
                        edgecolor=colordict[n_center], facecolor='none',
                        linewidth=0.2))

    # add the arrow showing the PA
    pa_arrow(ax, x, y, pa, size_pix / 10, scale)

    # idn = [i for (i,val) in enumerate(ifu_ids) if val=='006']
    # xoff,yoff=ifu_centers[idn][0]
    # rpa=-1.*pa*numpy.pi/180.
    # xnoff = (1.*xoff*numpy.cos(rpa)-numpy.sin(rpa)*yoff)
    # ynoff = 1.*xoff*numpy.sin(rpa)+numpy.cos(rpa)*yoff
    # print(xoff/numpy.cos(dec/180.*numpy.pi),yoff,xnoff/numpy.cos(dec*numpy.pi/180.),ynoff)
    # x,y =
    # wcs2pix(ra+xoff/numpy.cos(dec/180.*numpy.pi),dec+yoff,ra_obj,dec_obj,scale=scale,im_size=size_pix,CD=CD)
    # ax.add_patch(Circle((x,y),deg2pix(10./3600.,
    # scale),edgecolor=colordict[n_center],facecolor='none',linewidth=0.2))
    # x,y =
    # wcs2pix(ra+xnoff/numpy.cos(dec*numpy.pi/180.),dec+ynoff,ra_obj,dec_obj,scale=scale,im_size=size_pix,CD=CD)
    # ax.add_patch(Circle((x,y),deg2pix(10./3600.,
    # scale),edgecolor=colordict[n_center],facecolor='none',linewidth=0.2))

    print('Image retrieved from ' + url)
    fg.savefig(filename, dpi=dpi)

    print('Visualization image is saved to %s' % filename)
    plt.close(fg)


def rotate(angle, size, points, resample=0, expand=0, center=None,
           translate=None):

    """
    mF: essential rotation code stolen from PIL.Image to 
    compute x,y coordinates after rotation.
    """
    angle = angle % 360.0

    w, h = size

    post_trans = (0, 0)
    rotn_center = (w / 2.0, h / 2.0)  

    angle = - math.radians(angle)
    matrix = [
        round(math.cos(angle), 15), round(math.sin(angle), 15), 0.0,
        round(-math.sin(angle), 15), round(math.cos(angle), 15), 0.0
    ]

    def transform(x, y, matrix):
        (a, b, c, d, e, f) = matrix
        return a*x + b*y + c, d*x + e*y + f

    matrix[2], matrix[5] = transform(-rotn_center[0] - post_trans[0],
                                     -rotn_center[1] - post_trans[1], matrix)
    matrix[2] += rotn_center[0]
    matrix[5] += rotn_center[1]


    augm = numpy.array( [[matrix[0],matrix[1],matrix[2]],[matrix[3],matrix[4],matrix[5]],[0.,0.,1.]] )
    result = []
    for p in points:
        result.append( augm.dot( numpy.array([p[0],p[1],1.]))[0:2] )

    return result


def visualize_acam(s, acamloc, ifu_centers, ifu_ids, guideWFSSol,
                   cal_star_cand, config, filename, targetID, fplane_provider=None):
    """Plot a shuffle solution with matplotlib along side object of interest
    and save image.

    Parameters
    ----------
    s : list
        ra, dec of the center of the field input into shuffle and pa
    acamloc : list
        ra and dec of the center of the field after running shuffle, pa, ?, ?
    ifu_centers : ndarray
        list of [x, y] positions of the IFUs in (units?)
    ifu_ids : ndarray
        list of slot ids for each of the ifu_centers
    guideWFSSol : ndarray
        ?
    cal_star_cand : ndarray
        ?
    config : :class:`configparser.ConfigParser`
        configuration object
    filename : string
        name of the file where to save the plot
    """
    log = logging.getLogger('shuffle')
    acam_offset = config.getfloat('offsets', 'acam')
    acam_size = 209.0   # arcseconds
    # image constants
    img_size = acam_size * 2. / 3600.
    dpi = 300.
    # acam  constants
    acam_scale = config.getfloat('General', 'acam_pix_scale')  # arsec * pix
    acam_corner_center_x = (config.getfloat('General', 'acam_x_length')+1)/2.
    acam_corner_center_y = (config.getfloat('General', 'acam_y_length')+1)/2.
    # x = -222, y = 24 is the pixel position in the ACAM of the center of the
    # focal plane
    acam_ihmp_x = config.getfloat('General', 'acam_x_origin')
    acam_ihmp_y = config.getfloat('General', 'acam_y_origin')
    dx_fplane_acam = config.getfloat('General', 'acam_x_origin') \
        - acam_corner_center_x
    dy_fplane_acam = config.getfloat('General', 'acam_y_origin') \
        - acam_corner_center_y

    ra_shuffle, dec_shuffle = acamloc[:2]
    ra_input, dec_input, pa = s[:3]
    # rotate the pa by 90 degrees plus the offset
    pa_offset = pa + 90 + acam_offset
    pa_offset_rad = numpy.deg2rad(pa_offset)

    # get all the labels in the same place
    vertical_alignment = 'top' if (pa_offset % 360) < 180 else 'bottom'

    # Get image at centered at the center of the ACAM
    imarray, CD, url, img_src = retrieve_image(ra_shuffle, dec_shuffle,
                                               img_size, config,
                                               config.getboolean("General",
                                                                 "yflip"),
                                               0.25)

    # NEW: delegate ACAM visualization to the cleaner implementation
    try:
        # Build an input WCS for the retrieved image using the authoritative
        # sky center and the returned CD/pixel scale. Define WCS entirely from
        # CRVAL=(ra_shuffle,dec_shuffle), CRPIX=(nx+1)/2,(ny+1)/2 and CD.
        ny, nx = imarray.shape[:2]
        w_in = WCS(naxis=2)
        w_in.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w_in.wcs.crval = [ra_shuffle, dec_shuffle]
        w_in.wcs.crpix = [(nx + 1.0) / 2.0, (ny + 1.0) / 2.0]
        # CD is a numpy.matrix; convert to ndarray
        w_in.wcs.cd = numpy.array([[float(CD[0, 0]), float(CD[0, 1])],
                                   [float(CD[1, 0]), float(CD[1, 1])]])
        hdu = fits.PrimaryHDU(data=imarray, header=w_in.to_header())
        hdul = fits.HDUList([hdu])

        # Build candidate star tuples: (ra, dec, label, g, r, i)
        cand = []
        for o in cal_star_cand:
            try:
                label = f"{int(o[0])}:{int(o[1])}"
            except Exception:
                label = "star"
            # Try to extract g, r, i if available (legacy indices: -6,-5,-4)
            try:
                gmag = float(o[-6])
                rmag = float(o[-5])
                imag = float(o[-4])
            except Exception:
                gmag = numpy.nan
                rmag = numpy.nan
                imag = numpy.nan
            cand.append((float(o[2]), float(o[3]), label, gmag, rmag, imag))

        visualize_acam_clean(
            ra_shuffle,
            dec_shuffle,
            pa,
            hdul,
            cand,
            ifu_centers,
            ifu_ids,
            config,
            filename,
            target_ra=ra_input,
            target_dec=dec_input,
            fplane_provider=fplane_provider,
        )
        return None
    except Exception:
        # Fall back to legacy plotting below if anything goes wrong
        pass

    # size, center and pixel scale for the image

    y_pix_img, x_pix_img = imarray.shape[:2]
    img_center_x, img_center_y = x_pix_img / 2., y_pix_img / 2.

    # pixel scale from CD, in arcsec / pix
    img_scale_x = 3600.0 * numpy.sqrt(CD[0, 0] ** 2 + CD[1, 0] ** 2)
    img_scale_y = 3600.0 * numpy.sqrt(CD[0, 1] ** 2 + CD[1, 1] ** 2)
    # create figure
    fg = plt.figure(figsize=(x_pix_img/dpi, y_pix_img/dpi), dpi=dpi,
                    frameon=False)
    ax = fg.add_axes([0, 0, 1, 1], frameon=False)
    ax.axis('off')

    # add the image
    # Apply robust asinh stretch with black background and percentile upper bound
    ax.imshow(imarray, origin='lower', cmap='gray')
    # add a circle with the center of the focal plane
    ax.add_patch(Circle((img_center_x, img_center_y),
                        deg2pix(5./3600., img_scale_x),
                        edgecolor='red', facecolor='none',
                        linewidth=0.2))

    ax.plot([img_center_x, img_center_x],
            [img_center_y, img_center_y+30], 'r-')
    ax.text(img_center_x, img_center_y+40, 'N', color='r',
            rotation=-pa_offset, ha='center')
    ax.plot([img_center_x, img_center_x-30],
            [img_center_y, img_center_y], 'r-')
    ax.text(img_center_x-30, img_center_y, 'E', color='r',
            rotation=-pa_offset, va='center')
    # add the arrow showing the PA
    pa_arrow(ax, img_center_y, img_center_y, pa, x_pix_img / 13, img_scale_x,
             va=vertical_alignment, text_rot=-pa-(90 + acam_offset))

    # get the center of the acam w.r.t the 0,0 coordinate of the figure
    # convert the dx and dx into image pixel
    # diff_fplane_acam = (dx_fplane_acam**2 + dy_fplane_acam**2) ** 0.5
    # x_center_acam = diff_fplane_acam * numpy.cos(pa_offset_rad - numpy.pi)
    # y_center_acam = diff_fplane_acam * numpy.sin(pa_offset_rad - numpy.pi)

    acam_center_angle = pa_offset_rad  # - numpy.pi
    x_center_acam = (dx_fplane_acam * numpy.cos(acam_center_angle) +
                     dy_fplane_acam * numpy.sin(acam_center_angle))
    y_center_acam = (dx_fplane_acam * numpy.sin(acam_center_angle) -
                     dy_fplane_acam * numpy.cos(acam_center_angle))
    # convert to the plotted image scale
    x_center_acam *= acam_scale / img_scale_x
    y_center_acam *= acam_scale / img_scale_y
    # the offset is from the center of the image
    x_center_acam = img_center_x + x_center_acam
    y_center_acam = img_center_y + y_center_acam
    # half diagonal of the acam square
    radius = acam_size * numpy.sqrt(2) / (2 * img_scale_x)

    # add the rectangle of the acam camera
    acam_rect = RegularPolygon((x_center_acam, y_center_acam), 4,
                               radius=radius,
                               orientation=pa_offset_rad+numpy.pi/4., ec='g',
                               color='none')
    ax.add_patch(acam_rect)
    ds9_regions = config.get('directories', 'ds9_regions')
    mkpath(ds9_regions)
    region_file = os.path.join(ds9_regions, '{:s}_acam.reg'.format(targetID))
    f_ds9region = open(region_file, 'w')
    DS9region.writeHeader(f_ds9region)
    # enhance and label the stars
    tan_plane = TP(ra_shuffle, dec_shuffle, pa_offset, acam_scale,
                   acam_scale)
    positions, mags = ([], [])
    for o in cal_star_cand:
        x, y = wcs2pix(o[2], o[3], ra_shuffle, dec_shuffle,
                       scale=img_scale_x, im_size=x_pix_img, CD=CD)
        if not acam_rect.contains_point((x, y)):
            continue  # only plot stars within acam
        try:
            if numpy.abs(o[-1]) > 0.0:
                x2, y2 = wcs2pix(o[2]-o[-2], o[3]-o[-1], ra_shuffle,
                                 dec_shuffle, scale=img_scale_x,
                                 im_size=x_pix_img, CD=CD)
                ax.add_patch(Arrow(x2, y2, (x-x2), (y-y2), color='magenta',
                                   alpha=0.5))
        except Exception:
            dummy = 0.0
        ax.add_patch(Circle((x, y), deg2pix(0.0015, img_scale_x),
                            edgecolor='magenta', facecolor='none',
                            linewidth=0.2))
        # distance from the star to the center in image coordinates
        xn, yn = tan_plane.raDec2xy(o[2], o[3])
        xn += acam_ihmp_x
        yn += acam_ihmp_y
        positions.append([xn, yn])
        mags.append(o[5])
        DS9region.writeRegion(f_ds9region, xn, yn, 5.)
        ax.text(x, y, '{0:.1f},{1:.1f}'.format(float(xn), float(yn)),
            color='magenta', size=1.5/img_scale_x, ha='center',
                va=vertical_alignment, rotation=-pa_offset, zorder=10,
                clip_on=False)
    positions = numpy.array(positions)
    mags = numpy.array(mags)
    inds = numpy.argsort(mags)
    allind = numpy.arange(len(positions), dtype=int)
    if len(inds):
        goodchoice = positions[inds[0]]
    else:
        goodchoice = None
    for i, ind in enumerate(inds):
        rest = numpy.delete(allind, ind)
        sel1 = numpy.abs(positions[ind, 0] - positions[rest, 0]) < 20
        sel2 = numpy.abs(positions[ind, 1] - positions[rest, 1]) < 20
        sel3 = numpy.abs(mags[ind] - mags[rest]) < 3.
        if not numpy.any(sel1 * sel2 * sel3):
            goodchoice = positions[ind]
            break
    x_input, y_input = tan_plane.raDec2xy(ra_input, dec_input)
    x_input += acam_ihmp_x
    y_input += acam_ihmp_y
    log.info("Source centered at: %0.1f, %0.1f" % (x_input, y_input))
    ax.add_collection(plotFocalPlaneLRS(0, 0, pa, img_scale_x, ifu_centers,
                                        ifu_ids, dec_shuffle, x_pix_img,
                                        config, color='green', text_ax=ax,
                                        text_rot=-pa_offset))

    # Compute the x/y limits and flip the image along x
    verts = acam_rect.get_verts()
    padd = numpy.array([acam_corner_center_x / 5., acam_corner_center_y / 5.])
    x_min_v, y_min_v = verts.min(axis=0) - padd
    x_max_v, y_max_v = verts.max(axis=0) + padd
    x_min_ax, x_max_ax = ax.get_xlim()
    y_min_ax, y_may_ax = ax.get_ylim()
    x_min = max(x_min_v, x_min_ax)
    x_max = min(x_max_v, x_max_ax)
    y_min = max(y_min_v, y_min_ax)
    y_max = min(y_max_v, y_may_ax)
    ax.set_xlim(left=x_max, right=x_min)
    ax.set_ylim(bottom=y_min, top=y_max)

    # Add a small buffer so labels near the ACAM edge are not cut off
    try:
        bx0, bx1 = ax.get_xlim()
        by0, by1 = ax.get_ylim()
        buf = 20
        # Expand regardless of axis direction (x may be inverted)
        ax.set_xlim(bx0 - buf, bx1 + buf)
        ax.set_ylim(by0 - buf, by1 + buf)
    except Exception:
        pass

    # mF
    # we maintain a list of guide star coordinates in the image
    # such that we can later calculate their rotated coordinates
    cal_star_cand_radec = []
    cal_star_cand_acamxy = []
    cal_star_cand_imxy = []
    cal_star_cand_gri = []
    a = ax.transAxes
    cal_star = Table(dtype=[float,float,float,float,float,float], names=["RA","Dec","acam_x","acam_y","img_x","img_y" ])
    for o in cal_star_cand:
            x, y = wcs2pix(o[2], o[3], ra_shuffle, dec_shuffle,
                           scale=img_scale_x, im_size=x_pix_img, CD=CD)
            if not acam_rect.contains_point((x, y)):
                continue  # only plot stars within acam
            xn, yn = tan_plane.raDec2xy(o[2], o[3])
            xn += acam_ihmp_x
            yn += acam_ihmp_y
            cal_star_cand_radec.append([o[2], o[3]])
            cal_star_cand_gri.append([o[-6], o[-5], o[-4]])
            cal_star_cand_acamxy.append([xn,yn])
            imxy = ax.transData.transform([x,y])
            cal_star_cand_imxy.append([imxy[0],imxy[1]])

    fg.savefig(filename, dpi=dpi)

    log.debug('Visualization image is saved to %s', filename)
    plt.close(fg)

    # rotate the image to have the ACAM rectangle straight
    from PIL import Image
    im1 = Image.open(filename)

    rot_cal_star_cand_imxy = rotate(pa_offset, im1.size, cal_star_cand_imxy)
    with open(filename.replace(".jpg",".csv"),"w") as f: 
        s = "#RA, Dec, acam_x, acam_y, img_x, img_y, g, r, i\n"
        f.write(s)
        for (ra,dec),(x,y),(ix,iy), (g, r, i) in zip(cal_star_cand_radec, 
                                                     cal_star_cand_acamxy, 
                                                     rot_cal_star_cand_imxy,
                                                     cal_star_cand_gri):
            s = "{:10.6f}, {:11.6f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}\n".format( ra,dec,x,y,ix,iy,g,r,i )
            f.write(s)
    
    im1.rotate(pa_offset).save(filename)
    log.debug('Visualization rotated image is saved to %s', filename)
    return goodchoice




def visualize_probestars(guideWFSSol, config, basename, targetID, tr):
    """Plot a shuffle solution with matplotlib along side object of interest
    and save image to ``filename``.

    Parameters
    ----------
    guideWFSSol : nd-array
        list of 4 stars, with ID, ra and dec in the columns [1, 2, 3]
    config : :class:`configparser.ConfigParser`
        configuration option
    filename : string
        name of the file where to save the plot
    """
    log = logging.getLogger('shuffle')

    # The size of the image is the max of 2x size of IFU or 4x Petrosian Radius
    size = 20./3600.
    dpi = 300.
    colors = ['yellow', 'yellow', 'white', 'white']
    labels = ['gd1', 'gd2', 'wfs1', 'wfs2']
    if not targetID:
        targetID = 'blank'

    for (ID, ra, dec), color, label in zip(guideWFSSol[:, 1:4],
                                           colors, labels):
        fig, ax = plt.subplots(1, figsize=(4, 4))

        imarray, CD, url, img_src = retrieve_image(ra, dec, size, config,
                                                   config.getboolean("General",
                                                                     "yflip"),
                                                   0.05)
        log.debug('Guide/Wf star image retrieved from %s', url)
        size_pix = len(imarray)
        scale = size*3600./size_pix

        ax.imshow(imarray, origin='lower', cmap='gray')

        # plot the guide star and WFS solution
        x, y = wcs2pix(ra, dec, ra, dec, scale=scale, im_size=size_pix, CD=CD)
        ax.add_patch(Circle((x, y), deg2pix(0.00055,  scale),
                            edgecolor=color, facecolor='none',
                            linewidth=0.2))
        ax.set_title('{l}, id: {id_}'.format(l=label, id_=ID))
        ax.set_axis_off()
        fn = op.join(basename, '%s_%s_%s.png' % (targetID, label, tr))
        fig.savefig(fn, dpi=dpi)
        log.debug('Visualization image is saved to %s', fn)
    plt.close(fig)


def visualizeDS9(s, r, ifu_centers, ifu_ids, guideWFSSol, cal_star_cand,
                 config, ifu_ihmpid, orig_loc, targetID):
    """Plot a shuffle solution in ds9"""
    # NOTE: leave the pyds9 import in the functions using it, this way the
    # setting of the xpa_method in the code has effect
    log = logging.getLogger('shuffle')
    import pyds9
    ds9 = pyds9.DS9()
    ds9.set('height 1000')
    ds9.set('width 1200')

    ra, dec, pa = s[0], s[1], s[2]

    # get image from sdss
    size = config.getfloat('General', 'dfplane')
    size += 2. * config.getfloat('General', 'radius')
    ds9.set('dsseso size %f %f degrees' % (size, size))
    ds9.set('frame delete all')
    ds9.set('frame new')
    ds9.set('dsseso coord %f %f degrees' % (ra, dec))
    ds9.set('dsseso close')

    # plot the original and shuffled focal plane in ds9
    plotFocalPlaneDS9(ds9, r[0], r[1], pa, ifu_centers, config, color='green',
                      ifu_ids=ifu_ihmpid, usename=True)

    # compass and ruler
    width_img = float(ds9.get('fits width'))
    height_img = float(ds9.get('fits height'))
    length = (width_img + height_img) / 2. / size / 60.*2.
    # TODO: Set compass to wcs coordinate system, such that it points to
    # correct direction
    # ds9.set('region', 'wcs; compass %fi %fi %f # compass=wcs N E color=white'
    # % (width_img/17., height_img/15., 2./60.))
    ds9.set('region', 'image; ruler %f %f %f %f # ruler=arcmin color=white' %
            (width_img/15., height_img/30., width_img/15.+length,
             height_img/30.))

    # RA and Dec at center
    ds9.set('region', 'image; text %i %i {%s} # color=grey' %
            (width_img/2, height_img/2+20, 'RA:%.3f deg' % ra))
    ds9.set('region', 'image; text %i %i {%s} # color=grey' %
            (width_img/2, height_img/2-20, 'Dec:%.3f deg' % dec))

    # find those stars which are suitable candidates for IFU calibration stars
    insiders = findStars.inside(r[0], r[1], pa, ifu_centers, cal_star_cand,
                                ifu_size=0.012, edge=0.000833, PLOT=False)
    for o in cal_star_cand[~insiders]:
        # ds9.set('region', 'wcs; circle %f %f .0015 # color=blue width=2
        # text={%d:%d} tag={%s}' % (o[2],o[3],o[0],o[1],'DR9'))
        ds9.set('region',
                'wcs; circle %f %f .0015 # color=blue width=2' %
                (o[2], o[3]))
    for o in cal_star_cand[insiders]:
        ds9.set('region', 'wcs; circle %f %f .0015 # color=magenta width=2'
                ' text={%d:%d}' % (o[2], o[3], o[0], o[1]))

    # if statement if guidewfssol doesn't exist

    colors = ['yellow', 'yellow', 'white', 'white']
    names = ['gd1', 'gd2', 'wfs1', 'wfs2']
    for s, t, c in zip(guideWFSSol, names, colors):
        ds9.set('region', 'wcs; circle %f %f .005 # color=%s width=2'
                ' text={%s:%d:%d}' % (s[2], s[3], c, t, s[0], s[1]))

    if orig_loc is not None:
        ds9.set('region', 'wcs; point %f %f # point=x color=red size=8'
                ' ' % (orig_loc[0], orig_loc[1]))

    ds9.set('scale asinh')
    ds9.set('scale mode zmax')
    ds9.set('regions system wcs')
    ds9.set('regions sky fk5')
    ds9.set('wcs align yes')
    ds9.set('zoom 0.75 0.75')

    # wait for the user to do stuff
    input('Hit ENTER...')
    # save the region file
    ds9_regions = config.get('directories', 'ds9_regions')
    mkpath(ds9_regions)
    region_file = os.path.join(ds9_regions, '{:s}'.format(targetID))
    try:
        ds9.set('regions save {:s}.reg'.format(region_file))
    except Exception:
        log.warning("Can't save regions file: %s" % region_file)
    ds9.set('regions skyformat degrees')
    regions = ds9.get('regions list')
    str_fib = regions.split('\n')
    for i, strf in enumerate(str_fib):
        if strf.split('(')[0] == 'circle':
            if "color=white" in strf.split(' '):
                log.info("Found WFS star.")
                s = re.split('[\( \) \,]', strf)
                loc = strf.find('wfs')
                log.info("RA: %0.5f, Dec: %0.5f" % (float(s[1]), float(s[2])))
                log.info("Probe: %s" % strf[loc:loc+4])
                for j, name in enumerate(names):
                    if name == strf[loc:loc+4]:
                        guideWFSSol[j][2] = float(s[1])
                        guideWFSSol[j][3] = float(s[2])
            if "color=yellow" in strf.split(' '):
                log.info("Found Guide star.")
                s = re.split('[\( \) \,]', strf)
                loc = strf.find('gd')
                log.info("RA: %0.5f, Dec: %0.5f" % (float(s[1]), float(s[2])))
                log.info("Probe: %s" % strf[loc:loc+3])
                for j, name in enumerate(names):
                    if name == strf[loc:loc+3]:
                        guideWFSSol[j][2] = float(s[1])
                        guideWFSSol[j][3] = float(s[2])
    return guideWFSSol

def retrieve_image_hips(ra, dec, fov_deg, yflip=True, scale=1.5):
    surveys = ["CDS/P/PanSTARRS/DR1/g", "CDS/P/DSS2/blue", "CDS/P/SDSS9/g"]
    url = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
    width = int(fov_deg * 3600. / scale)
    height = int(fov_deg * 3600. / scale)
    # hips2fits defines fov as the size of the largest image dimension
    pixscale = float(fov_deg) / float(max(width, height))   # deg / pixel
    cd11, cd22 = -pixscale, pixscale
    if yflip:
        cd22 *= -1.0

    for survey in surveys:
        params = {
            "hips": survey,
            "ra": ra,
            "dec": dec,
            "fov": fov_deg,
            "width": width,
            "height": height,
            "format": "fits",
            "projection": "TAN",
        }
        r = requests.get(url, params=params, timeout=20)
        r.raise_for_status()

        hdul = fits.open(BytesIO(r.content))
        imarray = hdul[0].data

        if imarray is not None and min(imarray.shape[-2:]) > 32:
            if yflip:
                imarray = imarray[::-1]

            CD = numpy.matrix([[cd11, 0.0],
                           [0.0,  cd22]])

            return imarray, CD, r.url, survey
    return return_blank_image(width, scale)



def _build_acam_wcs(ra0, dec0, pa_deg, acam_offset_deg,
                    acam_x_origin, acam_y_origin,
                    acam_pix_scale, acam_x_length, acam_y_length):
    """
    Build the output WCS in ACAM detector coordinates.

    The pixel (acam_x_origin, acam_y_origin) is the shuffle center (ra0, dec0).
    ACAM pixel scale is in arcsec/pixel.
    """
    theta = numpy.deg2rad(pa_deg - 90.0 + acam_offset_deg)
    s = acam_pix_scale / 3600.0  # deg / pixel
    c, sn = numpy.cos(theta), numpy.sin(theta)

    # Detector x/y to sky east/north, with east left on the sky image
    # Apply a full-frame Y-axis inversion so increasing pixel Y points toward
    # detector-bottom (south on the sky). This matches the historical ACAM
    # orientation where the scene appears vertically flipped (LRS2-R above LRS2-B).
    # Implemented by negating the second column of the CD matrix.
    cd = numpy.array([
        [-s * c, -s * sn],
        [ s * sn, -s * c],
    ])

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [acam_x_origin + 1.0, acam_y_origin + 1.0]  # FITS is 1-based
    w.wcs.cd = cd
    w.array_shape = (int(round(acam_y_length)), int(round(acam_x_length)))
    return w


def _sky_to_acam(ra, dec, ra0, dec0, pa_deg, acam_offset_deg,
                 acam_pix_scale, acam_x_origin, acam_y_origin, wcs=None):
    """
    Convert sky coordinates to ACAM x,y in pixels.

    If an ACAM WCS (output grid) is provided via `wcs`, use it directly
    for the most accurate sky->pixel transform. Otherwise, fall back to
    the analytic small-angle mapping consistent with _build_acam_wcs.
    """
    # Prefer authoritative WCS if available
    if wcs is not None:
        try:
            x, y = wcs.wcs_world2pix(float(ra), float(dec), 0)  # origin=0 (pixel coords match imshow)
            return float(x), float(y)
        except Exception:
            # Fall through to analytic mapping on any failure
            pass

    east = (ra - ra0) * numpy.cos(numpy.deg2rad(dec0)) * 3600.0
    north = (dec - dec0) * 3600.0

    theta = numpy.deg2rad(pa_deg - 90.0 + acam_offset_deg)
    c, s = numpy.cos(theta), numpy.sin(theta)

    dx_arc = -c * east + s * north
    dy_arc =  s * east + c * north

    x = acam_x_origin + dx_arc / acam_pix_scale
    # Inverted detector Y: increasing pixel Y goes toward south, so subtract
    y = acam_y_origin - dy_arc / acam_pix_scale
    return x, y


def _draw_cardinal_directions(ax, x0, y0, pa_deg, acam_offset_deg, length=60):
    """
    Draw N and E arrows directly in ACAM coordinates.
    """
    theta = numpy.deg2rad(pa_deg - 90.0 + acam_offset_deg)
    c, s = numpy.cos(theta), numpy.sin(theta)

    # north in ACAM pixels (Y is inverted: increasing pixel Y is south)
    nx =  s * length
    ny = -c * length

    # east in ACAM pixels (Y is inverted)
    ex = -c * length
    ey = -s * length

    ax.plot([x0, x0 + nx], [y0, y0 + ny], "r-", lw=0.6)
    ax.plot([x0, x0 + ex], [y0, y0 + ey], "r-", lw=0.6)
    ax.text(x0 + nx * 1.15, y0 + ny * 1.15, "N", color="r", ha="center", va="center", fontsize=8)
    ax.text(x0 + ex * 1.15, y0 + ey * 1.15, "E", color="r", ha="center", va="center", fontsize=8)


def _lrs2_box(center_arcsec, size_arcsec=(12.0, 6.0)):
    """
    Simple rectangle around an IFU center in ACAM-native tangent-plane units.
    Adjust size if you want a more exact footprint.
    """
    cx, cy = center_arcsec
    hx, hy = 0.5 * size_arcsec[0], 0.5 * size_arcsec[1]
    return numpy.array([
        [cx - hx, cy - hy],
        [cx + hx, cy - hy],
        [cx + hx, cy + hy],
        [cx - hx, cy + hy],
    ])


def _arcsec_to_acam_poly(poly_arcsec, acam_pix_scale, acam_x_origin, acam_y_origin, acam_offset_deg):
    """
    Convert tangent-plane offsets (arcsec) to ACAM pixel coordinates.

    Preferred path: if `wcs` and `center_ra/center_dec` are provided, convert
    each EN offset (arcsec) to absolute RA/Dec around (center_ra, center_dec),
    transform with the authoritative ACAM WCS via world2pix, and return the
    pixel coordinates. When the caller intends to add these to the center
    pixel separately, this function will instead return RELATIVE pixel offsets
    (vertex_pix - center_pix) so existing call sites that add (cx, cy) still
    work unchanged.

    Fallback path: if WCS is not provided, apply a rotation by acam_offset
    (degrees) only — no additional 90° term — then scale by the ACAM pixel
    scale and translate by the provided ACAM origin.

    Parameters
    ----------
    poly_arcsec : (N, 2) array-like
        Offsets in arcsec, columns are (east, north) relative to shuffle center
        or IFU center in the EN plane.
    acam_pix_scale : float
        ACAM pixel scale in arcsec/pixel.
    acam_x_origin, acam_y_origin : float
        ACAM detector pixel coordinates of the shuffle center (used only in
        the fallback, non-WCS path).
    acam_offset_deg : float
        Rotation to apply (degrees), consistent with ACAM offset definition
        (used only in the fallback, non-WCS path).
    wcs : astropy.wcs.WCS, optional
        ACAM output WCS. If provided with center_ra/center_dec, use for the
        preferred world2pix mapping.
    center_ra, center_dec : float, optional
        Sky coordinates of the polygon center, in degrees. Required when wcs
        is provided.

    Returns
    -------
    ndarray of shape (N, 2)
        If WCS path is used: relative pixel offsets from the center pixel
        (so you should add the center pixel coordinates externally).
        If fallback path is used: absolute ACAM pixel coordinates.
    """


    # Fallback analytic path (no WCS)
    theta = numpy.deg2rad(acam_offset_deg)
    c, s = numpy.cos(theta), numpy.sin(theta)
    dx = numpy.asarray(poly_arcsec)[:, 0]
    dy = numpy.asarray(poly_arcsec)[:, 1]
    dx_r = c * dx - s * dy
    dy_r = s * dx + c * dy
    x = acam_x_origin + dx_r / acam_pix_scale
    # Inverted detector Y in fallback as well
    y = acam_y_origin - dy_r / acam_pix_scale
    return numpy.column_stack([x, y])



def plotFocalPlaneLRS_acam(pa, ra0, dec0, ifu_centers_deg, ifu_ids, config, color='lime', 
                           linewidth=0.8, text_ax=None, text_rot=None, wcs=None, fplane_provider=None):
    """
    Plot LRS2/HPF overlays using ACAM pixel coordinates.

    - ifu_centers_deg are in degrees (IFU-frame x,y relative to shuffle center).
    - Convert IFU centers to RA/Dec using astrometry TangentPlane, then map to ACAM using the ACAM WCS/orientation.

    Returns a PatchCollection that can be added to axes.
    """
    log = logging.getLogger('shuffle')

    # Geometry sizes
    lrs_sizex_deg = config.getfloat("General", "lrs_sizex")
    lrs_sizey_deg = config.getfloat("General", "lrs_sizey")
    lrs_sizex_arc = lrs_sizex_deg * 3600.0
    lrs_sizey_arc = lrs_sizey_deg * 3600.0

    # Offsets and ACAM mapping params
    fplane_off = config.getfloat('offsets', 'fplane')
    acam_offset = config.getfloat('offsets', 'acam')
    acam_x_origin = config.getfloat('General', 'acam_x_origin')
    acam_y_origin = config.getfloat('General', 'acam_y_origin')
    acam_pix_scale = config.getfloat('General', 'acam_pix_scale')

    # Precompute rotations
    rpa_fplane = numpy.deg2rad(pa + fplane_off)
    cf, sf = numpy.cos(rpa_fplane), numpy.sin(rpa_fplane)

    if text_ax is None:
        text = plt.text
    elif text_ax == 'none':
        text = None
    else:
        text = text_ax.text

    patches = []

    # Normalize IFU IDs to strings and extract a 3-digit core for matching
    raw_ids = []
    core_ids = []
    for id_ in ifu_ids:
        try:
            s = str(id_, 'utf-8') if isinstance(id_, (bytes, bytearray)) else str(id_)
        except Exception:
            s = str(id_)
        s = s.strip()
        raw_ids.append(s)
        # Extract the last 3 consecutive digits (handles formats like "[056]", "IFU056", etc.)
        m = re.findall(r"(\d{3})", s)
        core = m[-1] if m else s
        core_ids.append(core)

    wanted = ['056', '066', '600', '603']
    # Build list: (raw_id, core_id, center)
    lrs_ifus = [(rid, cid, cen) for rid, cid, cen in zip(raw_ids, core_ids, ifu_centers_deg)
                if cid in wanted]

    # If none matched (unexpected formatting), fall back to plotting all IFUs
    plot_all = False
    if len(lrs_ifus) == 0:
        plot_all = True
        lrs_ifus = list(zip(raw_ids, core_ids, ifu_centers_deg))

    label_map = {'056': 'LRS2B', '066': 'LRS2R', '600': 'HPFACAM', '603': 'HPFSCI'}

    # Tangent plane at shuffle center for IFU->sky conversion (fallback when no provider)
    tp = TP(ra0, dec0, pa + fplane_off, 1.0, 1.0)

    for (ifu_raw, ifu_core, cen_deg) in lrs_ifus:
        lab = label_map.get(ifu_core, ifu_core)
        # IFU center in degrees relative to shuffle center, convert via TP to RA/Dec
        x_ifu_deg = float(cen_deg[0])
        y_ifu_deg = float(cen_deg[1])
        if fplane_provider is not None:
            # Prefer lookup by core id, then try raw id
            sky = fplane_provider.sky_centroid(ifu_core)
            if sky is None:
                sky = fplane_provider.sky_centroid(ifu_raw)
            if sky is None:
                # Fallback to local TP if provider lacks this ID
                ra_c, dec_c = tp.xy2raDec(3600.0 * x_ifu_deg, 3600.0 * y_ifu_deg)
            else:
                ra_c, dec_c = sky
        else:
            ra_c, dec_c = tp.xy2raDec(3600.0 * x_ifu_deg, 3600.0 * y_ifu_deg)

        # Map center sky position to ACAM pixels
        cx_pix, cy_pix = _sky_to_acam(ra_c, dec_c, ra0, dec0, pa, acam_offset,
                                      acam_pix_scale, acam_x_origin, acam_y_origin, wcs=wcs)

        # Build rectangle around center in IFU frame, then rotate to EN
        if 'LRS2' in lab:
            hy, hx = 0.5 * lrs_sizex_arc, 0.5 * lrs_sizey_arc
            rect_offsets = numpy.array([[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]])
        elif 'HPFACAM' in lab:
            # 15" square
            hx = hy = 7.5
            rect_offsets = numpy.array([[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]])
        else:
            rect_offsets = None

        if rect_offsets is not None:
            # Keep rectangles axis-aligned in DETECTOR pixel coordinates.
            # Build vertices directly in pixels around the center, avoiding any WCS/rotation.
            hx_pix = (abs(rect_offsets[1,0]) if rect_offsets[1,0] != 0 else abs(rect_offsets[2,0])) / acam_pix_scale
            hy_pix = (abs(rect_offsets[0,1]) if rect_offsets[0,1] != 0 else abs(rect_offsets[2,1])) / acam_pix_scale
            poly_pix = numpy.array([
                [cx_pix - hx_pix, cy_pix - hy_pix],
                [cx_pix + hx_pix, cy_pix - hy_pix],
                [cx_pix + hx_pix, cy_pix + hy_pix],
                [cx_pix - hx_pix, cy_pix + hy_pix],
            ])
            patches.append(Polygon(poly_pix, closed=True))
        # Always add a small center marker so IFU centers are visible even if polygon fails
        patches.append(Circle((cx_pix, cy_pix), 3.0))

        # Label at center
        log.debug("%s center RA/Dec=(%.6f, %.6f) -> ACAM(%.1f, %.1f)", lab, ra_c, dec_c, cx_pix, cy_pix)
        if text:
            text(cx_pix, cy_pix+30., lab, ha="center", va="center", family='sans-serif', size=7, color=color, rotation=text_rot)

    pc = PatchCollection(patches, facecolor='none', linewidth=linewidth)
    try:
        pc.set_edgecolor(color)
    except Exception:
        pass
    pc.set_zorder(20)
    pc.set_alpha(1.0)
    return pc


def visualize_acam_clean(
    ra_shuffle,
    dec_shuffle,
    pa,
    im_hdul,
    candidate_stars,
    ifu_centers,
    ifu_ids,
    config,
    outfile,
    target_ra=None,
    target_dec=None,
    fplane_provider=None,
):
    """
    Parameters
    ----------
    im_hdul : FITS HDUList
        Downloaded HiPS image with valid celestial WCS.
    candidate_stars : iterable
        Each entry should provide at least (ra, dec, label).
    ifu_centers : iterable
        Tangent-plane IFU centers in degrees (east, north) relative to shuffle center.
    ifu_ids : iterable
        IFU IDs matching ifu_centers.
    """

    acam_offset = config.getfloat("offsets", "acam")
    acam_x_origin = config.getfloat("General", "acam_x_origin")
    acam_y_origin = config.getfloat("General", "acam_y_origin")
    acam_pix_scale = config.getfloat("General", "acam_pix_scale")
    acam_x_length = config.getint("General", "acam_x_length")
    acam_y_length = config.getint("General", "acam_y_length")

    # 1) Build ACAM-native output WCS
    w_out = _build_acam_wcs(
        ra_shuffle, dec_shuffle, pa, acam_offset,
        acam_x_origin, acam_y_origin,
        acam_pix_scale, acam_x_length, acam_y_length
    )

    # 2) Reproject downloaded image onto ACAM grid
    acam_image, footprint = reproject_interp(
        im_hdul,
        w_out,
        shape_out=(acam_y_length, acam_x_length),
        return_footprint=True,
        order="bilinear",
    )

    # Optional: mask places with no input coverage
    acam_image = numpy.array(acam_image, dtype=float)
    acam_image[footprint <= 0] = numpy.nan

    # Work with finite values only for stats
    finite = numpy.isfinite(acam_image)
    vals = acam_image[finite]
    bkg = biweight_location(vals)
    try:
        var = biweight_midvariance(vals)
        sigma = numpy.sqrt(var) if numpy.isfinite(var) else mad_std(vals)
    except Exception:
        sigma = mad_std(vals)
    vmin = bkg - 1.0 * sigma
    vmax = bkg + 25.0 * sigma
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(a=0.2))

    # 4) Plot directly in ACAM coordinates
    dpi = 150
    fig = plt.figure(figsize=(acam_x_length / dpi, acam_y_length / dpi), dpi=dpi, frameon=False)
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.axis("off")
    ax.imshow(acam_image, origin="lower", cmap="gray", norm=norm)

    # Add a small pixel buffer around the ACAM to avoid cutting off labels near edges
    pad = max(10, int(0.02 * min(acam_x_length, acam_y_length)))
    ax.set_xlim(-pad, acam_x_length - 1 + pad)
    ax.set_ylim(-pad, acam_y_length - 1 + pad)

    # Draw a border showing the true ACAM footprint for visual reference
    ax.add_patch(Rectangle((0, 0), acam_x_length, acam_y_length,
                           ec="white", fc="none", lw=0.8, zorder=5))

    # shuffle center
    ax.add_patch(Circle((acam_x_origin, acam_y_origin), 8, ec="red", fc="none", lw=0.7))

    # cardinal directions
    _draw_cardinal_directions(ax, acam_x_origin, acam_y_origin, pa, acam_offset, length=70)

    # If target RA/Dec provided, mark where it lands on the ACAM image with a star
    if (target_ra is not None) and (target_dec is not None):
        tx, ty = _sky_to_acam(
            float(target_ra), float(target_dec),
            ra_shuffle, dec_shuffle,
            pa, acam_offset,
            acam_pix_scale, acam_x_origin, acam_y_origin,
            wcs=w_out,
        )
        ax.plot(tx, ty, marker='*', color='yellow', markersize=4, markeredgewidth=1.5, zorder=16)
        ax.text(tx + 8, ty + 8, f"Target ({tx:.1f}, {ty:.1f})", color='yellow', fontsize=8, ha='right', va='bottom', zorder=16)

    # candidate stars
    # Accumulate CSV rows: RA, Dec, acam_x, acam_y, img_x, img_y, g, r, i
    csv_rows = []
    for star in candidate_stars:
        ra_star, dec_star, label = star[:3]
        # Optional g, r, i provided as 4th-6th elements
        gmag = rmag = imag = numpy.nan
        if len(star) >= 6:
            try:
                gmag = float(star[3])
                rmag = float(star[4])
                imag = float(star[5])
            except Exception:
                gmag = rmag = imag = numpy.nan
        x_acam, y_acam = _sky_to_acam(
            ra_star, dec_star,
            ra_shuffle, dec_shuffle,
            pa, acam_offset,
            acam_pix_scale, acam_x_origin, acam_y_origin,
            wcs=w_out,
        )

        if (0 <= x_acam < acam_x_length) and (0 <= y_acam < acam_y_length):
            ax.add_patch(Circle((x_acam, y_acam), 10, ec="magenta", fc="none", lw=0.8))
            ax.text(
                x_acam, y_acam + 12,
                f"({x_acam:.1f}, {y_acam:.1f})",
                color="magenta", fontsize=7, ha="center", va="bottom",
                clip_on=False, zorder=10
            )
            # For the clean schema, img_x/img_y are the same as acam_x/acam_y
            csv_rows.append((ra_star, dec_star, x_acam, y_acam, x_acam, y_acam, gmag, rmag, imag))

    # LRS2 / HPF overlays in ACAM coordinates (ifu_centers are in degrees)
    ax.add_collection(
        plotFocalPlaneLRS_acam(
            pa,
            ra_shuffle,
            dec_shuffle,
            ifu_centers,
            ifu_ids,
            config,
            color='lime',
            linewidth=0.8,
            text_ax=ax,
            text_rot=None,
            wcs=w_out,
            fplane_provider=fplane_provider,
        )
    )

    # Write CSV next to the image output, matching legacy header but with img_x/img_y = acam_x/acam_y
    try:
        csv_path = op.splitext(outfile)[0] + ".csv"
        with open(csv_path, "w") as f:
            f.write("#RA, Dec, acam_x, acam_y, img_x, img_y, g, r, i\n")
            for ra_s, dec_s, ax_s, ay_s, ix_s, iy_s, g_s, r_s, i_s in csv_rows:
                f.write(
                    "{:10.6f}, {:11.6f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}\n".format(
                        ra_s, dec_s, ax_s, ay_s, ix_s, iy_s, g_s, r_s, i_s
                    )
                )
    except Exception:
        # Be robust: CSV creation should not break visualization
        pass

    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)