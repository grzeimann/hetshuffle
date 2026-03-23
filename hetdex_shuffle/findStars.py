from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import os
import time
import warnings
import sys
import datetime
import requests
import uuid
import psycopg2
import psycopg2.extras

from astropy.io.votable import parse_single_table
from astroquery.mast import Catalogs

import numpy
import scipy.spatial.distance as cdist
from astroquery.vo_conesearch import vos_catalog
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.io import ascii

import json
try:
    from urllib import pathname2url as urlencode
except ImportError:
    from urllib.request import pathname2url as urlencode
try:
    import httplib
except ImportError:
    import http.client as httplib
from astropy.table import Table
from astropy.time import Time

from . import sqlcl
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from astroquery.sdss import SDSS


def force(config, probe, candidates):
    forcestr = config.get("General", probe)
    # cut to the last ten digits since that is how the candidate catalog does
    # it
    shortforce10 = float(forcestr[-10:])
    # cut to the last nine digits since that is how the candidate catalog does
    # it
    shortforce9 = float(forcestr[-9:])
    sel10 = shortforce10 == candidates[:, 1]
    sel9 = shortforce9 == candidates[:, 1]
    sel1 = sel10 + sel9
    if sum(sel1) > 0:
        return candidates[sel1]
    else:
        # adapt the message
        raise ValueError("The selected " + probe + " id is not valid")


def aatof(tt):
    '''Convert the input string or list of stings to a list of floats

    Parameters
    ----------
    tt : string or list of strings
        strings to convert

    Returns
    -------
    list of floats
    '''
    return [float(t) for t in tt]


def toLocal(ra, dec, pa, sources):
    # translate source positions into focal plane coordianate frame
    lsources = sources.copy()
    pa_rad = pa/180.*numpy.pi
    # xl = (sources[:,2] - ra)
    xl = aDist2(sources[:, 2],  ra)
    yl = (sources[:, 3] - dec)
    xl = xl*numpy.cos(dec/180.*numpy.pi)
    # rotate source position into folcal place frame
    xr = xl*numpy.cos(pa_rad) - yl*numpy.sin(pa_rad)
    yr = xl*numpy.sin(pa_rad) + yl*numpy.cos(pa_rad)
    lsources[:, 2] = xr
    lsources[:, 3] = yr
    return lsources


def toGlobal(ra, dec, pa, ls):
    # translate source positions into focal plane coordinate frame
    xr, yr = ls[:, 2], ls[:, 3]
    s = ls.copy()
    pa_rad = pa/180.*numpy.pi
    # rotate source position into global place frame
    xl = xr*numpy.cos(-pa_rad) - yr*numpy.sin(-pa_rad)
    yl = xr*numpy.sin(-pa_rad) + yr*numpy.cos(-pa_rad)
    x = xl/numpy.cos((dec)/180.*numpy.pi)

    x = x + ra

    x[x < 0] += 360.

    y = yl + dec

    s[:, 2] = x
    s[:, 3] = y
    return s


def xy2rv(xo, yo):
    r = numpy.sqrt(xo**2. + yo**2.)
    t = numpy.zeros(xo.shape)

    ii = xo > 0.
    t[ii] = numpy.arctan(yo[ii]/xo[ii])
    ii = (xo < 0.) * (yo >= 0.)
    t[ii] = (numpy.arctan(yo[ii]/xo[ii]) + numpy.pi)
    ii = (xo < 0.) * (yo < 0.)
    t[ii] = (numpy.arctan(yo[ii]/xo[ii]) - numpy.pi)
    ii = (xo == 0.) * (yo > 0.)
    t[ii] = numpy.pi/2.
    ii = (xo == 0.) * (yo < 0.)
    t[ii] = -numpy.pi/2.
    ii = (xo == 0.) * (yo == 0.)
    t[ii] = 0.
    return r, (-(t-numpy.pi/2.)/numpy.pi*180.) % 360.


def rv2xy(r, v):
    x = - r*numpy.sin(-v/180.*numpy.pi)
    y = r*numpy.cos(-v/180.*numpy.pi)
    return x, y


def toLocalRV(ra, dec, pa, sources):
    # translate source positions into focal plane polar coordinate frame
    lsources = toLocal(ra, dec, pa, sources)
    rvsources = lsources.copy()
    r, v = xy2rv(lsources[:, 2], lsources[:, 3])
    rvsources[:, 2] = r
    rvsources[:, 3] = v
    return rvsources


def rvToGlobal(ra, dec, pa, rvsources):
    # translate local focal plane polar coordinates to global
    r, v = rvsources[:, 2], rvsources[:, 3]
    x, y = rv2xy(r, v)
    lsources = rvsources.copy()
    lsources[:, 2] = x
    lsources[:, 3] = y
    return toGlobal(ra, dec, pa, lsources)


def dist(s, ss):
    dsq = (ss[:, 1]-s[1])**2. + (ss[:, 2]-s[2])**2.
    t = numpy.sqrt(dsq)
    return t


def aDist(s, ss):
    # calculates angular separation
    ad = ss[:, 3] - s[3]
    ad[ad >= 180.] -= 360.
    ad[ad < -180.] += 360.
    return ad


def aDist2(aa, ra):
    # calculates angular separation
    ad = aa - ra

    ad[ad >= 180.] -= 360.
    ad[ad < -180.] += 360.
    return ad


def stars_inside_acam(ra, dec, oo, pa, config):
    '''
    Identify stars that fall within the ACAM
    ra : float, list, tuple, or array
        RA for candidate stars
    dec : float, list, tuple, or array
        Dec for candidate stars
    poly : RA and Dec corners of
    '''

    log = logging.getLogger('shuffle')

    xc = config.getfloat("General", "acam_x_origin")
    yc = config.getfloat("General", "acam_y_origin")
    scale = config.getfloat("General", "acam_pix_scale") / 3600.
    xlen = config.getfloat("General", "acam_x_length")
    ylen = config.getfloat("General", "acam_y_length")
    acam_rot = config.getfloat("offsets", "acam")
    pa_offset = pa + 90. + acam_rot
    pa_offset_rad = numpy.deg2rad(pa_offset)

    loo = toLocal(ra, dec, 0., oo)

    xd = [0., xlen]
    yd = [0., ylen]
    x_acam_corners_unrotated = [(x-xc)*scale for x in xd]
    y_acam_corners_unrotated = [(y-yc)*scale for y in yd]
    x_acam_corners_rotated = [x * numpy.cos(pa_offset_rad)
                              + y * numpy.sin(pa_offset_rad)
                              for x in x_acam_corners_unrotated
                              for y in y_acam_corners_unrotated]
    y_acam_corners_rotated = [x * -1.*numpy.sin(pa_offset_rad)
                              + y * numpy.cos(pa_offset_rad)
                              for x in x_acam_corners_unrotated
                              for y in y_acam_corners_unrotated]
    poly = []
    order = [0, 1, 3, 2]
    for i in range(len(order)):
        poly.append([x_acam_corners_rotated[order[i]],
                     y_acam_corners_rotated[order[i]]])

    n = len(poly)
    # log.info("The Four ACAM corners are: (%0.4f,%0.4f), (%0.4f,%0.4f), "
    #          "(%0.4f,%0.4f), (%0.4f,%0.4f)", *tuple([i for j in poly
    #                                                  for i in j]))
    # import matplotlib.pyplot as plt
    # fig = plt.figure(figsize=(6,6))
    # plt.plot([poly[0][0],poly[1][0]],[poly[0][1],poly[1][1]],'k-')
    # plt.plot([poly[1][0],poly[2][0]],[poly[1][1],poly[2][1]],'k-')
    # plt.plot([poly[2][0],poly[3][0]],[poly[2][1],poly[3][1]],'k-')
    # plt.plot([poly[3][0],poly[0][0]],[poly[3][1],poly[0][1]],'k-')
    # plt.scatter(poly[0][0],poly[0][1],marker='x',color='k',s=50,linewidth=2)
    # plt.scatter(loo[:,2],loo[:,3],marker='*',color=[0.95,0.34,0.45],
    #             edgecolor=[0.75,0.34,0.45],s=35)
    # plt.axis('equal')
    # plt.xlim([0.2,-0.2])
    # plt.ylim([-0.2,0.2])
    # fig.savefig('acam_quick_visual.pdf',dpi=150)
    # plt.close(fig)

    inside = numpy.array(numpy.zeros((len(loo),)), dtype=bool)
    for j in range(len(loo)):
        x = loo[j, 2]
        y = loo[j, 3]
        p1x, p1y = poly[0]
        for i in range(n+1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            inside[j] = not inside[j]
            p1x, p1y = p2x, p2y

    return inside


def queryGAIA2(ra, dec, boxsize, maxsources=10000):
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = logging.getLogger('shuffle')
    log.debug("queryGAIA: ra      = %f ", ra)
    log.debug("queryGAIA: dec     = %f ", dec)
    log.debug("queryGAIA: boxsize = %f ", boxsize)

    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS',
                             'phot_g_mean_mag', '_RAJ2000', '_DEJ2000',
                             'pmRA', 'pmDE'],
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='fk5')
    D = vquery.query_region(field, width=("%fd" % boxsize),
                            catalog="I/345/gaia2")
    try:
        Data = D[0]
    except Exception:
        return numpy.array([])
    oo = []
    g_i = 1.  # g-i color
    for i, obj in enumerate(Data['Source']):
        oid_a = Data['Source'][i]
        oid_b = 0
        ra = Data['_RAJ2000'][i]
        dec = Data['_DEJ2000'][i]
        pmra = Data['pmRA'][i]
        pmde = Data['pmDE'][i]
        G = Data['Gmag'][i]
        g = G + 0.0939 + 0.6758 * g_i + 0.04 * g_i**2 - 0.003 * g_i**3
        if numpy.any([j for j in Data.mask[i]]):
            continue
        oo.append([oid_a, oid_b, ra, dec, 0., g, 0., 0., 0., pmra, pmde])
    return numpy.array(oo)


def queryGAIA(ra, dec, boxsize, maxsources=10000):
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = logging.getLogger('shuffle')
    log.debug("queryGAIA: ra      = %f ", ra)
    log.debug("queryGAIA: dec     = %f ", dec)
    log.debug("queryGAIA: boxsize = %f ", boxsize)

    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS',
                             'phot_g_mean_mag', '_RAJ2000', '_DEJ2000'],
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra, dec=dec,
                           unit=(u.deg, u.deg),
                           frame='fk5')
    Data = vquery.query_region(field, width=("%fd" % boxsize),
                               catalog="I/337/gaia")[0]
    oo = []
    g_i = 1.  # g-i color
    for i, obj in enumerate(Data['Source']):
        oid_a = Data['Source'][i]
        oid_b = 0
        ra = Data['_RAJ2000'][i]
        dec = Data['_DEJ2000'][i]
        G = Data['__Gmag_'][i]
        g = G + 0.0939 + 0.6758 * g_i + 0.04 * g_i**2 - 0.003 * g_i**3
        if numpy.any([j for j in Data.mask[i]]):
            continue
        oo.append([oid_a, oid_b, ra, dec, 0., g, 0., 0., 0.])
    return numpy.array(oo)

def queryUSNO_B1(ra, dec, boxsize, stellarity_sel=2, maxsources=10000):
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = logging.getLogger('shuffle')
    log.debug("queryUSNO_B1: ra      = %f ", ra)
    log.debug("queryUSNO_B1: dec     = %f ", dec)
    log.debug("queryUSNO_B1: boxsize = %f ", boxsize)

    vquery = Vizier(columns=['USNO-B1.0', 'RAJ2000', 'DEJ2000',
                             'B1s/g', 'B1mag', 'R1mag'],
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra, dec=dec,
                           unit=(u.deg, u.deg),
                           frame='fk5')
    Data = vquery.query_region(field, width=("%fd" % boxsize),
                               catalog="I/284/out")[0]
    oo = []
    for i, obj in enumerate(Data['USNO-B1.0']):
        oid_a, oid_b = aatof(obj.split(b'-'))
        ra = Data['RAJ2000'][i]
        dec = Data['DEJ2000'][i]
        B = Data['B1mag'][i]
        R = Data['R1mag'][i]
        stellarity = Data['B1s_g'][i]
        if numpy.any([j for j in Data.mask[i]]):
            continue
        if stellarity >= stellarity_sel:
            oo.append([oid_a, oid_b, ra, dec, 0., B, R, 0., 0.])
    return numpy.array(oo)


def mastQuery(request):
    server = 'mast.stsci.edu'

    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent": "python-requests/" + version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')
    # Close the https connection
    conn.close()

    return head, content


def mastJson2Table(jsonObj):

    dataTable = Table()
    for col, atype in [(x['name'], x['type']) for x in jsonObj['fields']]:
        if atype == 'string':
            atype = 'str'
        if atype == 'boolean':
            atype = 'bool'
        tcol = []
        for x in jsonObj['data']:
            col_val = x.get(col, None)
            if col_val is None:
                tcol.append(-999.)
            else:
                tcol.append(col_val)
        dataTable[col] = numpy.array(tcol, dtype=atype)
    return dataTable


def queryTESS_IC(ra, dec, config, radius):
    log = logging.getLogger('shuffle')
    log.info('Looking for stars at %0.5f, %0.5f with radius %0.2f' % (ra, dec,
                                                                      radius))
    mashupRequest = {'service': 'Mast.Catalogs.Tic.Cone',
                     'params': {'ra': ra,
                                'dec': dec,
                                'radius': radius},
                     'format': 'json',
                     'pagesize': 10000,
                     'page': 1}

    headers, outString = mastQuery(mashupRequest)

    outData = json.loads(outString)

    table = mastJson2Table(outData)
    oo = []
    sel = numpy.where(table['objType'] == 'STAR')[0]
    # a = datetime.datetime.now()
    # b = datetime.datetime(2000, 1, 1)
    # tdiff = (a-b).total_seconds() / (365.24*3600.*24.)
    # nsel = obj>-999.
    # dra = ss[nsel, 12] / 1000. * tdiff
    # ddec = ss[nsel, 13] / 1000. * tdiff
    # ss[nsel, 2] = ss[nsel, 2] + dra
    # ss[nsel, 3] = ss[nsel, 3] + ddec
    for i, obj in enumerate(table[sel]):
        f = '%018d' % int(obj['ID'])
        oid_a = f[:9]
        oid_b = f[9:]
        ra = obj['ra']
        dec = obj['dec']

        if obj['gmag'] != -999.:
            u = obj['umag']
            g = obj['gmag']
            r = obj['rmag']
            i = obj['imag']
            z = obj['zmag']
        else:
            u, r, i, z = (0., 0., 0., 0.)
            if obj['Bmag'] > -999. and obj['Vmag'] > -999.:
                B = obj['Bmag']
                V = obj['Vmag']
                g = V + 0.06 * (B-V) - 0.12
            else:
                T = obj['Tmag']
                g = -1.054 + T * 1.19
        oo.append([oid_a, oid_b, ra, dec, u, g, r, i, z, obj['Teff'],
                   obj['mass'], obj['rad'], obj['pmRA'], obj['pmDEC']])
    return numpy.array(oo, dtype=float)

def ps1_cone_table(
    ra, dec,
    radius_deg=9.0/60.0,
    table="stack",
    release="dr2",
    columns=["objID", "raMean", "decMean",
            "gApMag", "rApMag", "iApMag", "zApMag", "yApMag",],
    baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs",
    **kw,
):
    """
    Cone search PS1 and return an astropy Table.
    Uses CSV output (per STScI example) and is quite robust.
    """
    params = dict(ra=ra, dec=dec, radius=radius_deg, **kw)

    if columns:
        # PS1 API expects '[col1,col2,...]'
        params["columns"] = "[" + ",".join(columns) + "]"

    url = f"{baseurl}/{release}/{table}.csv"
    r = requests.get(url, params=params)
    r.raise_for_status()
    text = r.text

    # No results
    if not text.strip():
        return Table(names=columns or [])

    tab = ascii.read(text)

    # Optional: turn PS1 sentinel -999.0 into NaN for magnitudes
    for col in tab.colnames:
        if tab[col].dtype.kind in "fi":
            mask = tab[col] <= -999.0
            if numpy.any(mask):
                tab[col][mask] = numpy.nan

    return tab

def panstarrs_query1(ra_deg, dec_deg, rad_deg, mindet=1,
                    maxsources=30000,
                    server=('https://catalogs.mast.stsci.edu/api/v0.1/'
                            'panstarrs/dr2/mean.votable')):
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                mindet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
    """
    r = requests.get(server, params={'ra': ra_deg, 'dec': dec_deg,
                                     'radius': rad_deg, 'pagesize': maxsources,
                                     'nDetections.gte': ('%d' % mindet),
                                     'verify': False})

    # write query data into local file
    name = str(uuid.uuid4()) + '.xml'
    outf = open(name, 'w')
    outf.write(r.text)
    outf.close()
    # parse local file into astropy.table object
    data = parse_single_table(name)
    os.remove(name)
    return data.to_table(use_names_over_ids=True)

def panstarrs_query(ra_deg, dec_deg, rad_deg, mindet=1,
                    maxsources=30000):
    Pan = Catalogs.query_region("%s %s" % (ra_deg, dec_deg), radius=rad_deg, 
                               catalog="Panstarrs") 
    return Pan

def get_panstarrs_data(RA, DEC, SR, DB_USER, DB_PWD, DB_HOST, DB_DATABASE):
    try:
        dbLink = psycopg2.connect(user=DB_USER, password=DB_PWD, host=DB_HOST, database=DB_DATABASE)
    except (Exception, psycopg2.Error) as err:
        print(err, file=sys.stderr)
        exit(1)
    cur = dbLink.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    query = (
        'select "objID", "raMean", "decMean", "gPSFMag" as "gMeanPSFMag", '
        '"rPSFMag" as "rMeanPSFMag", "iPSFMag" as "iMeanPSFMag", "zPSFMag" as "zMeanPSFMag" '
        'from pshet '
        'where '
            'q3c_radial_query("raMean", "decMean", {}, {}, {}) '
            'and "primaryDetection">0 '
            'and ("rPSFMag"-"rKronMag")<0.05 '
            'and "gpsfQfPerfect">0.85 and "rpsfQfPerfect">0.85 and "ipsfQfPerfect">0.85 '
            'and "gPSFMag">0 and "rPSFMag">0 and "iPSFMag">0 '
    ).format(RA, DEC, SR)
    cur.execute(query)
    data = cur.fetchall()
    if (dbLink):
        cur.close()
        dbLink.close()
    return data


def het_panstarrs_data(ra_deg, dec_deg, rad_deg, config):
    """
    Query HET local stacked version of Pan-STARRS DR1
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
    returns: astropy.table object
    """
    DB_USER = config.get('database','db_user')
    DB_PWD = config.get('database','db_pwd')
    DB_HOST = config.get('database','db_host')
    DB_DATABASE = config.get('database','db_database')

    data = get_panstarrs_data(ra_deg, dec_deg, rad_deg,
                              DB_USER, DB_PWD, DB_HOST, DB_DATABASE)
    if len(data):
        T = Table(rows=data)
        return T
    else:
        return []


def queryPANSTARRS_new(ra, dec, radius, debug=False):
    """
    Queries PANSTARRS
    """
    T = ps1_cone_table(ra, dec, radius/2.)
    S = numpy.zeros((len(T), 9))
    for i, objid in enumerate(T['objID']):
        objid = str(objid)
        S[i, 0] = objid[:9]
        S[i, 1] = objid[9:]
    S[:, 2] = T['raMean']
    S[:, 3] = T['decMean']
    S[:, 4] = 0.0
    S[:, 5] = numpy.where(T['gApMag'] > 0.0, T['gApMag'], 0.0)
    S[:, 6] = numpy.where(T['rApMag'] > 0.0, T['rApMag'], 0.0)
    S[:, 7] = numpy.where(T['iApMag'] > 0.0, T['iApMag'], 0.0)
    S[:, 8] = numpy.where(T['zApMag'] > 0.0, T['zApMag'], 0.0)
    return S
    

def queryPANSTARRS(ra, dec, radius, debug=False):
    """
    Queries PANSTARRS
    """
    T = panstarrs_query(ra, dec, radius/2.)
    sel1 = (~T['gMeanPSFMag'].mask) * (~T['rMeanPSFMag'].mask) \
        * (~T['iMeanPSFMag'].mask)
    T = T[sel1]
    S = numpy.zeros((len(T), 9))
    for i, objid in enumerate(T['objID']):
        objid = str(objid)
        S[i, 0] = objid[:9]
        S[i, 1] = objid[9:]
    S[:, 2] = T['raMean']
    S[:, 3] = T['decMean']
    S[:, 4] = 0.0
    S[:, 5] = numpy.where(T['gMeanPSFMag'] > 0.0, T['gMeanPSFMag'], 0.0)
    S[:, 6] = numpy.where(T['rMeanPSFMag'] > 0.0, T['rMeanPSFMag'], 0.0)
    S[:, 7] = numpy.where(T['iMeanPSFMag'] > 0.0, T['iMeanPSFMag'], 0.0)
    S[:, 8] = numpy.where(T['zMeanPSFMag'] > 0.0, T['zMeanPSFMag'], 0.0)
    return S

def queryPANSTARRS_LOCAL(ra, dec, radius, config, debug=False):
    """
    Queries PANSTARRS
    """
    T = het_panstarrs_data(ra, dec, radius/2., config)
    if len(T) == 0:
        S = numpy.zeros((len(T), 9))
        return S
    S = numpy.zeros((len(T), 9))
    for i, objid in enumerate(T['objID']):
        objid = str(objid)
        S[i, 0] = objid[:9]
        S[i, 1] = objid[9:]
    S[:, 2] = T['raMean']
    S[:, 3] = T['decMean']
    S[:, 4] = 0.0
    S[:, 5] = numpy.where(T['gMeanPSFMag'] > 0.0, T['gMeanPSFMag'], 0.0)
    S[:, 6] = numpy.where(T['rMeanPSFMag'] > 0.0, T['rMeanPSFMag'], 0.0)
    S[:, 7] = numpy.where(T['iMeanPSFMag'] > 0.0, T['iMeanPSFMag'], 0.0)
    S[:, 8] = numpy.where(T['zMeanPSFMag'] > 0.0, T['zMeanPSFMag'], 0.0)
    return S


def queryUSNO_A2(ra, dec, radius, debug=False):
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = logging.getLogger('shuffle')
    log.debug("queryUSNO_A2: ra      = %f ", ra)
    log.debug("queryUSNO_A2: dec     = %f ", dec)
    log.debug("queryUSNO_A2: radius = %f ", radius)

    usno_a2_name = 'The USNO-A2.0 Catalogue (Monet+ 1998) 1'
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = vos_catalog.call_vo_service('conesearch_good', verbose=False,
                                             kwargs={'RA': ra, 'DEC': dec,
                                                     'SR': radius, 'VERB': 1,
                                                     'cftype': 'ASCII'},
                                             catalog_db=usno_a2_name)
    Data = result.array.data
    oo = []
    for i, obj in enumerate(Data['USNO-A2.0']):
        oid_a, oid_b = aatof(obj.split(b'-'))
        ra = Data['RAJ2000'][i]
        dec = Data['DEJ2000'][i]
        B = Data['Bmag'][i]
        V = Data['Rmag'][i]
        g = V + 0.06 * (B-V) - 0.12
        oo.append([oid_a, oid_b, ra, dec, 0., g, 0., 0., 0.])

    return numpy.array(oo)


def query2MASS(ra, dec, radius, debug=False):
    """
    Queries 2MASS.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. J H K 0.
    """
    log = logging.getLogger('shuffle')
    log.debug("query2MASS: ra      = %f ", ra)
    log.debug("query2MASS: dec     = %f ", dec)
    log.debug("query2MASS: radius = %f ", radius)

    twomass_aspsc_name = '2MASS All-Sky Point Source Catalog 2'
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = vos_catalog.call_vo_service('conesearch_good', verbose=False,
                                             kwargs={'RA': ra, 'DEC': dec,
                                                     'SR': radius, 'VERB': 1,
                                                     'cftype': 'ASCII'},
                                             catalog_db=twomass_aspsc_name)

    Data = result.array.data
    oo = []
    for i, obj in enumerate(Data['designation']):
        oid_a, oid_b = aatof(obj.split(b'-'))
        ra = Data['ra'][i]
        dec = Data['dec'][i]
        J = Data['j_m'][i]
        H = Data['h_m'][i]
        K = Data['k_m'][i]
        oo.append([oid_a, oid_b, ra, dec, 0., J, H, K, 0.])

    return numpy.array(oo)


def query_cwise(ra, dec, boxsize, debug=True, cachedir='cache'):
    '''Query the curewise database and fallback to the SDSS DR9 one if
    ``cx_Oracle`` package cannot be found

    Parameters
    ----------
    ra, dec : floats
        center of the field
    boxsize : float
        size of the box around ra/dec in which to look for sources
    debug : bool, optional
        enable extra debut prints
    cachedir : string, optional
        name of the directory where to cache sdss dr9 results

    Returns
    -------
    :class:`numpy.ndarray`
        return array of stars with format
        IDa IDb RA DEC u g r i z
    '''
    try:
        import cx_Oracle as cxo
    except ImportError:
        print('Failed to import cx_Oracle, falling back to SDSS.')
        return querySDSSDR9(ra, dec, boxsize, debug=debug, cachedir=cachedir)

    conn = cxo.connect('AWANONYMOUS', 'anonymous', 'hetdex')

    query = 'SELECT "SDSS_objID", ra, dec, "SDSS_modelMag_u",'
    query += ' "SDSS_modelMag_g", "SDSS_modelMag_r", "SDSS_modelMag_i",'
    query += ' "SDSS_modelMag_z"'
    query += ' FROM AWOPER."SOURCELIST*SOURCES**01"'
    query += ' WHERE "SDSS_type"=6'
    query += ' AND BITAND( "SDSS_flags", 262144 ) = 0'
    query += ' AND ra BETWEEN {0} AND {1}'
    query += ' AND dec BETWEEN {2} AND {3}'
    query += ' ORDER BY "SDSS_objID"'

    cur = conn.cursor()

    rasize = boxsize/numpy.cos(dec/180.*numpy.pi)/2.
    decsize = boxsize/2.

    cur.execute(query.format(ra-rasize, ra+rasize, dec-decsize, dec+decsize))

    res = cur.fetchall()

    oo = []

    if res:
        for l in res:
            oid_a = float(str(l[0])[:9])
            oid_b = float(str(l[0])[9:])
            oo.append([oid_a, oid_b] + list(l[1:]))

    return numpy.array(oo)


def querySDSS(ra, dec, radius):
    ''' using astroquery sdss system '''
    pos = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
    table = SDSS.query_region(pos, radius=radius*u.deg
                              / numpy.cos(dec*numpy.pi/180.),
                              photoobj_fields=['ra', 'dec', 'objid', 'type',
                                               'u', 'g', 'r', 'i', 'z'])
    table = table[numpy.where(table['type'] == 6)[0]]
    oo = numpy.zeros((len(table), 9))
    for i, t in enumerate(table):
        oo[i, 0] = str(t['objid'])[:9]
        oo[i, 1] = str(t['objid'])[9:]
        oo[i, 2] = t['ra']
        oo[i, 3] = t['dec']
        oo[i, 4] = t['u']
        oo[i, 5] = t['g']
        oo[i, 6] = t['r']
        oo[i, 7] = t['i']
        oo[i, 8] = t['z']

    return oo


def querySDSS_all(ra, dec, radius):
    ''' using astroquery sdss system '''
    pos = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
    table = SDSS.query_region(pos, radius=radius*u.deg
                              / numpy.cos(dec*numpy.pi/180.),
                              photoobj_fields=['ra', 'dec', 'objid', 'type',
                                               'u', 'g', 'r', 'i', 'z'])
    oo = numpy.zeros((len(table), 9))
    for i, t in enumerate(table):
        oo[i, 0] = str(t['objid'])[:9]
        oo[i, 1] = str(t['objid'])[9:]
        oo[i, 2] = t['ra']
        oo[i, 3] = t['dec']
        oo[i, 4] = t['u']
        oo[i, 5] = t['g']
        oo[i, 6] = t['r']
        oo[i, 7] = t['i']
        oo[i, 8] = t['z']

    return oo


def querySDSSDR9(ra, dec, boxsize, cachedir="cache"):
    '''Query the SDSS DR9 database for stars.

    Parameters
    ----------
    ra, dec : floats
        center of the field
    boxsize : float
        size of the box around ra/dec in which to look for sources
    cachedir : string, optional
        name of the directory where to cache sdss dr9 results

    Returns
    -------
    :class:`numpy.ndarray`
        return array of stars with format
        IDa IDb RA DEC u g r i z
    '''
    log = logging.getLogger('shuffle')
    USE_PHOTOOBJ = True
    flag = 0

    if USE_PHOTOOBJ:
        sql = ""
        sql += "select p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z "
        sql += "FROM PhotoPrimary AS p "
        sql += "WHERE type=6 "
        sql += "AND ((flags_r & 0x10000000) != 0)"
        sql += "AND ((flags_r & 0x8100000c00a4) = 0)"
        sql += "AND (((flags_r & 0x400000000000) = 0) or"
        sql += " (psfmagerr_g <= 0.03))"
        sql += "AND (((flags_r & 0x100000000000) = 0) or"
        sql += " (flags_r & 0x1000) = 0)"
        sql += "AND p.ra BETWEEN (%f-%f) AND (%f+%f) "
        sql += "AND p.dec BETWEEN (%f-%f) AND (%f+%f) "
        sql += "ORDER BY p.objid"
        sql = sql % (ra, boxsize/numpy.cos(dec/180.*numpy.pi)/2., ra,
                     boxsize/numpy.cos(dec/180.*numpy.pi)/2., dec, boxsize/2.,
                     dec, boxsize/2.)
    else:
        sql = ""
        sql += "select p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z "
        sql += "FROM PhotoPrimary AS p "
        sql += "WHERE type=6 "
        sql += "AND p.ra BETWEEN (%f-%f) AND (%f+%f) "
        sql += "AND p.dec BETWEEN (%f-%f) AND (%f+%f) "
        sql += "ORDER BY ra,dec"
        sql = sql % (ra, boxsize/numpy.cos(dec/180.*numpy.pi)/2., ra,
                     boxsize/numpy.cos(dec/180.*numpy.pi)/2., dec, boxsize/2.,
                     dec, boxsize/2.)

    log.debug(sql)

    ll = None

    cache_file = os.path.join(cachedir,
                              "{}_{}_{}.dat".format(ra, dec, boxsize))
    if os.path.exists(cache_file):
        log.debug("Found data in cache. File: %s", cache_file)
        with open(cache_file) as f:
            ll = f.readlines()
    else:
        url = sqlcl.public_url
        fmt = 'csv'

        file_ = sqlcl.query(sql, url, fmt)
        ll = file_.readlines()
        file_.close()

        # cache for future use
        with open(cache_file, 'w') as f:
            f.writelines([str(l.decode()) for l in ll])

    oo = []
    for l in ll[2:]:
        tt = l.split(b',')
        soid = tt[0]
        # the object ID has 18 digits, this is too long for a float, we
        # therefore split it into 2
        try:
            oid_a = float(soid[:9])
        except ValueError:
            flag = -1
            return numpy.array([]), flag
        oid_b = float(soid[9:])
        a = [oid_a, oid_b] + aatof(tt[1:])
        # print "%d" % a[0]
        oo.append(a)
    return numpy.array(oo), flag


def update_coords(cat1, cat2, current_epoch, log, epoch_cat1=2000.0):
    try:
        c1 = SkyCoord(cat1[:, 2]*u.deg, cat1[:, 3]*u.deg, frame='icrs')
        c2 = SkyCoord(cat2[:, 2]*u.deg, cat2[:, 3]*u.deg, frame='icrs')
        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        sel1 = d2d.arcsec < 5.
        sel2 = numpy.abs(cat1[:, 5] - cat2[idx, 5]) < 1.
        log.info('Number of sources with GAIA2 within 5": %i / %i' %
                 (sel1.sum(), len(cat1)))
        log.info('And with |g_cat - g_gaia2|<1: %i / %i' %
                 ((sel1*sel2).sum(), len(cat1)))
        sel = numpy.where(sel1*sel2)[0]
        deltaRA = ((current_epoch - epoch_cat1) * cat2[idx[sel], -2] /
                   1e3 / 3600. / numpy.cos(cat1[sel, 3] * numpy.pi / 180.))
        deltaDE = ((current_epoch - epoch_cat1) * cat2[idx[sel], -1] /
                   1e3 / 3600.)
        cat1[sel, 2] += deltaRA
        cat1[sel, 3] += deltaDE
        A = numpy.zeros((cat1.shape[0], cat1.shape[1]+2))
        A[:, :-2] = cat1
        A[sel, -2] = deltaRA
        A[sel, -1] = deltaDE
        return A
    except Exception:
        log.warning("Could NOT update astrometry with proper motion")
        return cat1


def getStars(sid, ra, dec,  pa, config, boxsize, catalog="SDSSDR9",
             local_catalogue=None, cachedir='cache'):
    '''Get stars from the desired catalogue.

    Parameters
    ----------
    sid : ?
        ???
    ra, dec : floats
        center of the field
    boxsize : float
        size of the box around ra/dec in which to look for sources
    catalog : string, optional
        name of the catalogue
    local_catalogue : string, optional
        path to a local catalogue
    cachedir : string, optional
        name of the directory where to cache sdss dr9 results

    Returns
    -------
    :class:`numpy.ndarray`
        return array of stars with format
        IDa IDb RA DEC u g r i z
    '''
    log = logging.getLogger('shuffle')
    log.debug("getStars: ra      = %f ", ra)
    log.debug("getStars: dec     = %f ", dec)
    log.debug("getStars: boxsize = %f ", boxsize)
    try:
        pmcat = queryGAIA2(ra, dec, boxsize/2.)
    except Exception:
        pmcat = []
    current_time = Time(datetime.datetime.now()).decimalyear
    # Try the local catalogue first as it gives SDSS survers a break
    if local_catalogue:
        oo = local_catalogue.query_preloaded_catalogue(ra, dec, boxsize)
        if len(oo) != 0:
            oo = update_coords(oo, pmcat, current_time, log)
            return oo, "LOCAL"
        else:
            log.warning("Using local catalogue failed! Defaulting to online"
                        " servers.")

    if catalog == "CWISE":
        oo = query_cwise(ra, dec, boxsize=boxsize, cachedir=cachedir)
        if len(oo) != 0:
            oo = update_coords(oo, pmcat, current_time, log)
            return oo, "CWISE"
        else:
            log.warning("No star is found in the local SDSS catalog, use"
                        " USNOA2 catalog instead.")
            return queryUSNO_A2(ra, dec, radius=boxsize/2.), "USNOA2"

    elif catalog == "SDSSDR9":
        # continue through loop if flag < 0 which means network failure
        flag = -1
        cnt = 0  # keep cnt
        maxcnt = 5  # go through loop at most 5 times
        while flag < 0:
            t1 = time.time()
            oo, flag = querySDSSDR9(ra, dec, boxsize=boxsize,
                                    cachedir=cachedir)
            if cnt >= maxcnt:
                log.warning("Max tries reached for SDSS catalog, use USNOA2"
                            " catalog instead.")
                return queryUSNO_A2(ra, dec, radius=boxsize/2.), "USNOA2"
            if len(oo) != 0:
                co = fastfindStarCand(oo, config, "acam", forcemagmax=0.0)
                inside_acam = stars_inside_acam(ra, dec, co, pa, config)
                log.info("%i ACAM stars found using SDSS" % (sum(inside_acam)))
                if sum(inside_acam) > 0:
                    oo = update_coords(oo, pmcat, current_time, log)
                    return oo, "SDSS"
                else:
                    log.warning("No SDSS found inside ACAM, using USNOA2"
                                " instead")
                    oo = queryUSNO_A2(ra, dec, radius=boxsize/2.)
                    co = fastfindStarCand(oo, config, "acam", forcemagmax=0.0)
                    inside_acam = stars_inside_acam(ra, dec, co, pa, config)
                    log.info("%i ACAM stars found using USNO_A2"
                             % (sum(inside_acam)))
                    oo = update_coords(oo, pmcat, current_time, log)
                    return oo, "USNOA2"

            else:
                if flag == 0:
                    log.warning("No star is found from SDSS catalog, use"
                                " USNOA2 catalog instead.")
                    return queryUSNO_A2(ra, dec, radius=boxsize/2.), "USNOA2"
                else:
                    cache_file = "%s/%f_%f_%f.dat" % ("cache", ra, dec,
                                                      boxsize)
                    if os.path.exists(cache_file):
                        log.debug("Removing cache. File: ", cache_file)
                        os.remove(cache_file)
                    cnt = cnt+1  # increase cnt
                    t2 = time.time()
                    dt = t2 - t1
                    # limit to one query per s (SDSS complains otherwise)
                    while (dt < 1.):
                        t2 = time.time()
                        dt = t2 - t1

    elif catalog == "USNOA2":
        oo = queryUSNO_A2(ra, dec, radius=boxsize/2.)
        if not len(oo):
            log.warning("No star found in primary search of USNOA2 catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "USNOA2"

    elif catalog == "SDSS":
        oo = querySDSS(ra, dec, radius=boxsize/2.)
        if not len(oo):
            log.warning("No star found in primary search of SDSS catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "SDSS"
    elif catalog == "SDSS_all":
        oo = querySDSS_all(ra, dec, radius=boxsize/2.)
        if not len(oo):
            log.warning("No star found in primary search of SDSS catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "SDSS_all"
    elif catalog == "TIC":
        oo = queryTESS_IC(ra, dec, config, radius=boxsize/2.)
        if not len(oo):
            log.warning("No star found in primary search of TESS catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "TIC"

    elif catalog == "USNOB1":
        oo = queryUSNO_B1(ra, dec, boxsize)
        if not len(oo):
            log.warning("No star found in primary search of USNOB1 catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "USNOB1"

    elif catalog == "GAIA":
        oo = queryGAIA(ra, dec, boxsize)
        if not len(oo):
            log.warning("No star found in primary search of GAIA catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "GAIA"
    elif catalog == 'PANSTARRS':
        oo = queryPANSTARRS_new(ra, dec, boxsize)
        if not len(oo):
            log.warning("No star found in primary search of "
                        "PANSTARRS catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "PANSTARRS"
    elif catalog == 'PANSTARRS_LOCAL':
        oo = queryPANSTARRS_LOCAL(ra, dec, boxsize, config)
        if not len(oo):
            log.warning("No star found in primary search of "
                        "PANSTARRS_LOCAL catalog.")
            log.error("Shuffle failed: Insufficient number of guider or WFS stars.")
            sys.exit(0)
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "PANSTARRS_LOCAL"
    elif catalog == "2MASS":
        oo = query2MASS(ra, dec, radius=boxsize/2.)
        if not len(oo):
            log.warning("No star found in primary search of 2MASS catalog.")
        oo = update_coords(oo, pmcat, current_time, log)
        return oo, "2MASS"
    else:
        raise ValueError("ERROR: Catalog %s unknown." % (catalog))


def removeDuplicates(oo, d=10./3600.):
    """
    Eliminates close pairs.
    """
    v = numpy.copy(oo[:, 2:4])
    v[:, 0] = numpy.cos(v[:, 1]*numpy.pi/180.)*v[:, 0]
    a = cdist.cdist(v, v)
    b = numpy.tril(a)
    b[b == 0] = 99
    flag = numpy.min(b, axis=1) > d
    result = oo[flag, :]
    return numpy.array(result)


def eliminateClose(oo, d=5./3600.):
    """
    Eliminates close pairs.
    """
    result = []
    for o in oo:
        ddsq = (numpy.cos(o[3]*numpy.pi/180.)*(oo[:, 2]-o[2]))**2.
        ddsq += (oo[:, 3]-o[3])**2.
        ii = (ddsq < d**2.)
        if sum(ii) == 1:
            result.append(o)
    return numpy.array(result)


def fastfindStarCand(oo, config, kind, forcemagmax=0.0):
    """
    Takes a list of stars and returns those which
    are suitable for IFU calibration stars.
    I.e. close pairs are eliminated,  and
    only those are returned which fit the magnitude limits.
    """
    dd = config.getfloat("MagLimits", kind+"_minsep")
    if dd > 0:
        ss = removeDuplicates(oo, d=dd)
    else:
        ss = oo
    # apply magnitude limits,  for IFU calibration stars
    if config.getboolean("General", "hpf") and (kind in ['guide1', 'guide2']):
        try:
            sel = numpy.where((ss[:, 9] > -999.) * (ss[:, 9] < 4000.) *
                              (ss[:, 10] < 0.5) * (ss[:, 11] < 0.5))[0]
            ss = ss[sel, :]
        except Exception:
            dummy = False
    magmin = aatof(config.get("MagLimits", kind+"_magmin").split(","))
    magmax = aatof(config.get("MagLimits", kind+"_magmax").split(","))
    for i, mag in enumerate(magmax):
        magmax[i] = magmax[i] + config.getfloat("MagLimits", kind+"_magadd")

    # Sungryong Added
    # forcefully overwrite the magmax to forcemagmax
    if forcemagmax > 0.0:
        magmax[1] = forcemagmax  # g band max 21
        isguide = kind.find("guide")
        if isguide > -1:
            magmax[2] = forcemagmax  # g band max 21

    # SDSSDR9 has ugriz 5 magnitudes, while USNOA2 have only B,R magnitudes
    ii = numpy.ones((len(ss),), dtype=bool)
    for i in numpy.arange(4, 9):
        if ss[:, i].all() != 0:
            ii *= (ss[:, i] >= magmin[i-4]) * (ss[:, i] <= magmax[i-4])
    return ss[ii]


def findStarCand(oo, config, kind):
    """
    Takes a list of stars and returns those which
    are suitable for IFU calibration stars.
    I.e. close pairs are eliminated, and
    only those are returned which fit the magnitude limits.
    """
    # discard close pairs (disable check by setting minsep to zero)
    dd = config.getfloat("MagLimits", kind+"_minsep")
    if dd > 0:
        ss = removeDuplicates(oo, d=dd)
    else:
        ss = oo
    # apply magnitude limits, for IFU calibration stars
    magmin = aatof(config.get("MagLimits", kind+"_magmin").split(","))
    magmax = aatof(config.get("MagLimits", kind+"_magmax").split(","))
    # SDSSDR9 has ugriz 5 magnitudes, while USNOA2 have only B,R magnitudes
    if ss[:, [5, 6, 7]].all() != 0:
        ii = (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        ii *= (ss[:, 7] >= magmin[3]) * (ss[:, 7] <= magmax[3])
        return ss[ii]
    else:
        ii = (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        return ss[ii]


def findStarCand20(oo, config, kind):
    """
    Takes a list of stars and returns those which
    are suitable for IFU calibration stars.
    I.e. close pairs are eliminated, and
    only those are returned which fit the magnitude limits.
    """
    # discard close pairs (disable check by setting minsep to zero)
    dd = config.getfloat("MagLimits", kind+"_minsep")
    if dd > 0:
        ss = removeDuplicates(oo, d=dd)
    else:
        ss = oo
    # apply magnitude limits, for IFU calibration stars
    magmin = aatof(config.get("MagLimits", kind+"_magmin").split(","))
    magmax = aatof(config.get("MagLimits", kind+"_magmax").split(","))

    # Sungryong Added
    # forcefully overwrite the magmax to 20
    magmax[1] = 20.0  # g band max 20
    isguide = kind.find("guide")
    if isguide > -1:
        magmax[2] = 20.0  # g band max 20
    # print "findStarCand20 : kind = ", kind
    # print "findStarCand20 : MagMaxLimits = ", magmax

    # SDSSDR9 has ugriz 5 magnitudes, while USNOA2 have only B,R magnitudes
    if ss[:, [4, 7, 8]].all() != 0:
        ii = (ss[:, 4] >= magmin[0]) * (ss[:, 4] <= magmax[0])
        ii *= (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        ii *= (ss[:, 7] >= magmin[3]) * (ss[:, 7] <= magmax[3])
        ii *= (ss[:, 8] >= magmin[4]) * (ss[:, 8] <= magmax[4])
        return ss[ii]
    else:
        ii = (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        return ss[ii]


def findStarCand21(oo, config, kind):
    """
    Takes a list of stars and returns those which
    are suitable for IFU calibration stars.
    I.e. close pairs are eliminated,  and
    only those are returned which fit the magnitude limits.
    """
    # discard close pairs (disable check by setting minsep to zero)
    dd = config.getfloat("MagLimits", kind+"_minsep")
    if dd > 0:
        ss = removeDuplicates(oo, d=dd)
    else:
        ss = oo
    # apply magnitude limits, for IFU calibration stars
    magmin = aatof(config.get("MagLimits", kind+"_magmin").split(","))
    magmax = aatof(config.get("MagLimits", kind+"_magmax").split(","))

    # Sungryong Added
    # forcefully overwrite the magmax to 21
    magmax[1] = 21.0  # g band max 21
    isguide = kind.find("guide")
    if isguide > -1:
        magmax[2] = 21.0  # g band max 21
    # print "findStarCand21 : kind = ", kind
    # print "findStarCand21 : MagMaxLimits = ", magmax

    # SDSSDR9 has ugriz 5 magnitudes, while USNOA2 have only B,R magnitudes
    if ss[:, [4, 7, 8]].all() != 0:
        ii = (ss[:, 4] >= magmin[0]) * (ss[:, 4] <= magmax[0])
        ii *= (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        ii *= (ss[:, 7] >= magmin[3]) * (ss[:, 7] <= magmax[3])
        ii *= (ss[:, 8] >= magmin[4]) * (ss[:, 8] <= magmax[4])
        return ss[ii]
    else:
        ii = (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        return ss[ii]


def findStarCand22(oo, config, kind):
    """
    Takes a list of stars and returns those which
    are suitable for IFU calibration stars.
    I.e. close pairs are eliminated,  and
    only those are returned which fit the magnitude limits.
    """
    # discard close pairs (disable check by setting minsep to zero)
    dd = config.getfloat("MagLimits", kind+"_minsep")
    if dd > 0:
        ss = removeDuplicates(oo, d=dd)
    else:
        ss = oo
    # apply magnitude limits, for IFU calibration stars
    magmin = aatof(config.get("MagLimits", kind+"_magmin").split(","))
    magmax = aatof(config.get("MagLimits", kind+"_magmax").split(","))

    # Sungryong Added
    # forcefully overwrite the magmax to 21
    magmax[1] = 22.0  # g band max 21
    isguide = kind.find("guide")
    if isguide > -1:
        magmax[2] = 22.0  # g band max 21
    # print "findStarCand21 : kind = ", kind
    # print "findStarCand21 : MagMaxLimits = ", magmax

    # SDSSDR9 has ugriz 5 magnitudes, while USNOA2 have only B,R magnitudes
    if ss[:, [4, 7, 8]].all() != 0:
        ii = (ss[:, 4] >= magmin[0]) * (ss[:, 4] <= magmax[0])
        ii *= (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        ii *= (ss[:, 7] >= magmin[3]) * (ss[:, 7] <= magmax[3])
        ii *= (ss[:, 8] >= magmin[4]) * (ss[:, 8] <= magmax[4])
        return ss[ii]
    else:
        ii = (ss[:, 5] >= magmin[1]) * (ss[:, 5] <= magmax[1])
        ii *= (ss[:, 6] >= magmin[2]) * (ss[:, 6] <= magmax[2])
        return ss[ii]


def inside(ra, dec, pa, ifu_centers, sources, ifu_size=0.012, edge=0.000833,
           PLOT=False):
    """
    Takes a list of stars and checks whether they fall into IFUs.
    ra = shot ra
    dec = shot dec
    ifu_centers = list with IFU center positions
    sources = list of sources that should be checked
    ifu_size = IFU size
    edge = width of edge. Stars will not be considered if they fall within this
    distannce of the edge of the IFU negative values are allowed here.
    returns list of booleans
    """
    xl = aDist2(sources[:, 2],  ra)
    yl = (sources[:, 3] - dec)
    xl = xl*numpy.cos(sources[:, 3]/180.*numpy.pi)
    yl = yl
    pa_rad = pa/180.*numpy.pi
    w = (ifu_size - edge)/2.
    # rotate source positions into focal plane frame
    xr = xl*numpy.cos(pa_rad) - yl*numpy.sin(pa_rad)
    yr = xl*numpy.sin(pa_rad) + yl*numpy.cos(pa_rad)
    result = []
    # for s in zip(xr,yr):
    #   ii =   (s[0] >= (ifu_centers[:,0]-w) ) * (s[0] <= (ifu_centers[:,0]+w))
    #   ii *=  (s[1] >= (ifu_centers[:,1]-w) ) * (s[1] <= (ifu_centers[:,1]+w))
    #     result.append(any(ii))
    xr.shape = (len(xr), 1)
    yr.shape = (len(xr), 1)

    a = numpy.abs(cdist.cdist(xr, ifu_centers[:, 0:1])) < w
    b = numpy.abs(cdist.cdist(yr, ifu_centers[:, 1:2])) < w
    c = a*b
    result = c.any(axis=1)
    # result = numpy.array(result)
    # return array([True])

    if PLOT:
        import pylab

        def ifu(s, x, y, pa, size):
            xx = numpy.array([-size/2., -size/2., size/2., size/2.])
            yy = numpy.array([-size/2., size/2., size/2., -size/2.])
            xxr = numpy.cos(pa) * xx + numpy.sin(pa) * yy
            yyr = - numpy.sin(pa) * xx + numpy.cos(pa) * yy
            # print xx,yy
            s.fill(xxr+x, yyr+y, facecolor='None', linewidth=1.0)

        s = pylab.subplot(111, aspect='equal')
        # s = pylab.subplot(111)
        for ic in ifu_centers:
            ifu(s, ic[0], ic[1], 0., ifu_size)
            ifu(s, ic[0], ic[1], 0., ifu_size-edge)

        pylab.plot(xr[result], yr[result], 'g.')
        pylab.plot(xr[~result], yr[~result], 'r.')
        xlim = pylab.xlim()
        pylab.xlim((xlim[1], xlim[0]))
        pylab.show()

    return result


def getGuideWFScand(ra, dec, pa, oo, config, inmagmax=0, debug=False):
    """
    Implements the logic of picking guide probe and wavefront sensor stars for
    one shot position.

    9/2015, implement new guide finding method by Karl, Sungryong, and Greg.

    In short:
        In "toLocalRV" angles, the allowable angles for each G1, G2, W1, W2, in
        the searching annulus,
        G1      15-190
        G2     200-10
        W1     290-100
        W2     110-280

        1) find all G1 and G2 stars, make pairs between G1 and G2,
            then find whether any good W1 W2 pair to support the chosen
            (trying) g1-g2 pair

    """
    # log = logging.getLogger('shuffle')

    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")
    dpatrol_g1min = config.getfloat("General", "dpatrol_g1min")
    dpatrol_g1max = config.getfloat("General", "dpatrol_g1max")
    dpatrol_g2min = config.getfloat("General", "dpatrol_g2min")
    dpatrol_g2max = config.getfloat("General", "dpatrol_g2max")
    dpatrol_w1min = config.getfloat("General", "dpatrol_w1min")
    dpatrol_w1max = config.getfloat("General", "dpatrol_w1max")
    dpatrol_w2min = config.getfloat("General", "dpatrol_w2min")
    dpatrol_w2max = config.getfloat("General", "dpatrol_w2max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = fastfindStarCand(oo, config, "guide1",
                                        forcemagmax=inmagmax)
    # print "Number of suitable guidestars for probe 1:",
    # len(guide_star1_cand) find those stars which are suitable candidates for
    # guiding in probe 2
    guide_star2_cand = fastfindStarCand(oo, config, "guide2",
                                        forcemagmax=inmagmax)
    # print "Number of suitable guidestars for probe 2:",
    # len(guide_star2_cand) find those stars which are suitable candidates for
    # WFS in probe 1
    wfs_star1_cand = fastfindStarCand(oo, config, "wfs1", forcemagmax=inmagmax)
    # print "Number of suitable wfsstars for probe 1:",  len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = fastfindStarCand(oo, config, "wfs2", forcemagmax=inmagmax)
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec,  pa,  guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec,  pa,  guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec,  pa,  wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec,  pa,  wfs_star2_cand)
    # pre-select all stars that fall into probe patrol region
    # TODO: use or/and
    cand_list = [guide_star1_cand, guide_star2_cand, wfs_star1_cand,
                 wfs_star2_cand]
    lcand_list = [lguide_star1_cand, lguide_star2_cand, lwfs_star1_cand,
                  lwfs_star2_cand]
    patrol_list = [[dpatrol_g1min, dpatrol_g1max],
                   [dpatrol_g2min, dpatrol_g2max],
                   [dpatrol_w1min, dpatrol_w1max],
                   [dpatrol_w2min, dpatrol_w2max]]
    empty_list = []
    for patrol, cand, lcand in zip(patrol_list, cand_list, lcand_list):
        if patrol[0] < patrol[1]:
            ii = ((lcand[:, 2] >= dpatrol_min/2.) *
                  (lcand[:, 2] <= dpatrol_max/2.) *
                  (lcand[:, 3] >= patrol[0]) *
                  (lcand[:, 3] <= patrol[1]))
        else:
            ii = ((lcand[:, 2] >= dpatrol_min/2.) *
                  (lcand[:, 2] <= dpatrol_max/2.) *
                  ((lcand[:, 3] >= patrol[0]) +
                  (lcand[:, 3] <= patrol[1])))
        empty_list.append(cand[ii])

    return empty_list


def fastfindGuideWFSSol(ra, dec, pa, oo, config, inmagmax=0, debug=False):
    """
    Implements the logic of picking guide probe and wavefront sensor stars for
    one shot position.

    9/2015, implement new guide finding method by Karl, Sungryong, and Greg.

    In short:
        In "toLocalRV" angles, the allowable angles for each G1, G2, W1, W2, in
        the searching annulus,
        G1      15-190
        G2     200-10
        W1     290-100
        W2     110-280

        1) find all G1 and G2 stars, make pairs between G1 and G2,
            then find whether any good W1 W2 pair to support the chosen
            (trying) g1-g2 pair

    """
    log = logging.getLogger('shuffle')

    # pair search version /Karl/Sungryong/Greg
    log.debug(">>>> Greg/Sungryong findGuideWFSSol ... ")

    probemind0 = config.getfloat("General", "probemind0")
    probemind1 = config.getfloat("General", "probemind1")
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")
    dpatrol_g1min = config.getfloat("General", "dpatrol_g1min")
    dpatrol_g1max = config.getfloat("General", "dpatrol_g1max")
    dpatrol_g2min = config.getfloat("General", "dpatrol_g2min")
    dpatrol_g2max = config.getfloat("General", "dpatrol_g2max")
    dpatrol_w1min = config.getfloat("General", "dpatrol_w1min")
    dpatrol_w1max = config.getfloat("General", "dpatrol_w1max")
    dpatrol_w2min = config.getfloat("General", "dpatrol_w2min")
    dpatrol_w2max = config.getfloat("General", "dpatrol_w2max")
    dpatrol_g1targ = config.getfloat("General", "dpatrol_g1targ")
    dpatrol_w1targ = config.getfloat("General", "dpatrol_w1targ")
    dpatrol_g2targ = config.getfloat("General", "dpatrol_g2targ")
    dpatrol_w2targ = config.getfloat("General", "dpatrol_w2targ")
    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = fastfindStarCand(oo, config, "guide1",
                                        forcemagmax=inmagmax)
    # print "Number of suitable guidestars for probe 1:",
    # len(guide_star1_cand) find those stars which are suitable candidates for
    # guiding in probe 2
    guide_star2_cand = fastfindStarCand(oo, config, "guide2",
                                        forcemagmax=inmagmax)
    # print "Number of suitable guidestars for probe 2:",
    # len(guide_star2_cand) find those stars which are suitable candidates for
    # WFS in probe 1
    wfs_star1_cand = fastfindStarCand(oo, config, "wfs1", forcemagmax=inmagmax)
    # print "Number of suitable wfsstars for probe 1:",  len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = fastfindStarCand(oo, config, "wfs2", forcemagmax=inmagmax)
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec,  pa,  guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec,  pa,  guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec,  pa,  wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec,  pa,  wfs_star2_cand)
    # pre-select all stars that fall into probe patrol region
    # TODO: use or/and
    if dpatrol_g1min < dpatrol_g1max:
        ii = ((lguide_star1_cand[:, 2] >= dpatrol_min/2.) *
              (lguide_star1_cand[:, 2] <= dpatrol_max/2.) *
              (lguide_star1_cand[:, 3] >= dpatrol_g1min) *
              (lguide_star1_cand[:, 3] <= dpatrol_g1max))
    else:
        ii = ((lguide_star1_cand[:, 2] >= dpatrol_min/2.) *
              (lguide_star1_cand[:, 2] <= dpatrol_max/2.) *
              ((lguide_star1_cand[:, 3] >= dpatrol_g1min) +
               (lguide_star1_cand[:, 3] <= dpatrol_g1max)))
    log.debug("Number of guide_star1 stars in patrol region: %s, %s, %s, %s",
              len(ii[ii]), lguide_star1_cand[ii, 2], dpatrol_min/2.,
              dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]
    else:
        log.info("Insufficients stars in guide_star1 patrol region.")
        return None

    if dpatrol_g2min < dpatrol_g2max:
        ii = ((lguide_star2_cand[:, 2] >= dpatrol_min/2.) *
              (lguide_star2_cand[:, 2] <= dpatrol_max/2.) *
              (lguide_star2_cand[:, 3] >= dpatrol_g2min) *
              (lguide_star2_cand[:, 3] <= dpatrol_g2max))
    else:
        ii = ((lguide_star2_cand[:, 2] >= dpatrol_min/2.) *
              (lguide_star2_cand[:, 2] <= dpatrol_max/2.) *
              ((lguide_star2_cand[:, 3] >= dpatrol_g2min) +
               (lguide_star2_cand[:, 3] <= dpatrol_g2max)))
    log.debug("Number of guide_star2 stars in patrol region: %s, %s, %s, %s",
              len(ii[ii]), lguide_star2_cand[ii, 2], dpatrol_min/2.,
              dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]
    else:
        log.info("Insufficients stars in guide_star2 patrol region.")
        return None

    if dpatrol_w1min < dpatrol_w1max:
        ii = ((lwfs_star1_cand[:, 2] >= dpatrol_min/2.) *
              (lwfs_star1_cand[:, 2] <= dpatrol_max/2.) *
              (lwfs_star1_cand[:, 3] >= dpatrol_w1min) *
              (lwfs_star1_cand[:, 3] <= dpatrol_w1max))
    else:
        ii = ((lwfs_star1_cand[:, 2] >= dpatrol_min/2.) *
              (lwfs_star1_cand[:, 2] <= dpatrol_max/2.) *
              ((lwfs_star1_cand[:, 3] >= dpatrol_w1min) +
               (lwfs_star1_cand[:, 3] <= dpatrol_w1max)))
    log.debug("Number of wfs1 stars in patrol region: %s, %s, %s, %s",
              len(ii[ii]), lwfs_star1_cand[ii, 2], dpatrol_min/2.,
              dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]
    else:
        log.info("Insufficients stars in wfs_star1 patrol region.")
        return None

    if dpatrol_w2min < dpatrol_w2max:
        ii = ((lwfs_star2_cand[:, 2] >= dpatrol_min/2.) *
              (lwfs_star2_cand[:, 2] <= dpatrol_max/2.) *
              (lwfs_star2_cand[:, 3] >= dpatrol_w2min) *
              (lwfs_star2_cand[:, 3] <= dpatrol_w2max))
    else:
        ii = ((lwfs_star2_cand[:, 2] >= dpatrol_min/2.) *
              (lwfs_star2_cand[:, 2] <= dpatrol_max/2.) *
              ((lwfs_star2_cand[:, 3] >= dpatrol_w2min) +
               (lwfs_star2_cand[:, 3] <= dpatrol_w2max)))
    log.debug("Number of wfs2 stars in patrol region: %s, %s, %s, %s",
              len(ii[ii]), lwfs_star2_cand[ii, 2], dpatrol_min/2.,
              dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]
    else:
        log.info("Insufficients stars in wfs_star2 patrol region.")
        return None

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    #######################
    # Force gp1, gp2, wfs1, and/or wfs2 if SDSS ID is given
    #######################

    if config.getint("General", "force_gp1") > 0:
        lguide_star1_cand = force(config, "force_gp1", lguide_star1_cand)
    if config.getint("General", "force_gp2") > 0:
        lguide_star2_cand = force(config, "force_gp2", lguide_star2_cand)
    if config.getint("General", "force_wfs1") > 0:
        lwfs_star1_cand = force(config, "force_wfs1", lwfs_star1_cand)
    if config.getint("General", "force_wfs2") > 0:
        lwfs_star2_cand = force(config, "force_wfs2", lwfs_star2_cand)
    gl1 = len(lguide_star1_cand)
    gl2 = len(lguide_star2_cand)
    gp_pair = numpy.zeros((gl1*gl2, 7))
    igp = 0
    for i in range(gl1):
        for j in range(gl2):
            gp_pair[igp, 0] = i
            gp_pair[igp, 1] = j
            gp_pair[igp, 2] = lguide_star1_cand[i, 3]
            gp_pair[igp, 3] = lguide_star2_cand[j, 3]
            gp_pair[igp, 4] = (lguide_star1_cand[i, gppickcol] +
                               lguide_star2_cand[j, gppickcol])/2.
            ab = numpy.abs(lguide_star1_cand[i, 3] - lguide_star2_cand[j, 3])
            if ab >= 180.:
                ab = 360. - ab
            gp_pair[igp, 5] = ab
            gp_pair[igp, 6] = (numpy.abs(gp_pair[igp, 2] - dpatrol_g1targ) +
                               numpy.abs(gp_pair[igp, 3] - dpatrol_g2targ))
            igp = igp + 1
    if config.getboolean("General", "use_brightness"):
        use_col = 4
    else:
        use_col = 6
    sortind = numpy.argsort(gp_pair[:, use_col])
    gp_pair = gp_pair[sortind]
    gp_pair = gp_pair[gp_pair[:, 5] > probemind1]

    wl1 = len(lwfs_star1_cand)
    wl2 = len(lwfs_star2_cand)
    wp_pair = numpy.zeros((wl1*wl2, 7))
    iwp = 0
    for i in range(wl1):
        for j in range(wl2):
            wp_pair[iwp, 0] = i
            wp_pair[iwp, 1] = j
            wp_pair[iwp, 2] = lwfs_star1_cand[i, 3]
            wp_pair[iwp, 3] = lwfs_star2_cand[j, 3]
            wp_pair[iwp, 4] = (lwfs_star1_cand[i, wfspickcol] +
                               lwfs_star2_cand[j, wfspickcol])/2.
            ab = numpy.abs(lwfs_star1_cand[i, 3] - lwfs_star2_cand[j, 3])
            if ab >= 180.:
                ab = 360. - ab
            wp_pair[iwp, 5] = ab
            wp_pair[iwp, 6] = (numpy.abs(wp_pair[iwp, 2] - dpatrol_w1targ) +
                               numpy.abs(wp_pair[iwp, 3] - dpatrol_w2targ))
            iwp = iwp + 1
    sortind = numpy.argsort(wp_pair[:, use_col])
    wp_pair = wp_pair[sortind]
    wp_pair = wp_pair[wp_pair[:, 5] > probemind1]
    for i in range(len(gp_pair)):
        for j in range(len(wp_pair)):
            dup = ((gp_pair[i, 2]-wp_pair[j, 2]) *
                   (gp_pair[i, 2]-wp_pair[j, 3]) *
                   (gp_pair[i, 3]-wp_pair[j, 2]) *
                   (gp_pair[i, 3]-wp_pair[j, 3]))
            if dup == 0:
                continue
            temp = numpy.zeros((4, ))
            temp[0:2] = gp_pair[i, 2:4]
            temp[2:4] = wp_pair[j, 2:4]
            temp_order = numpy.zeros((4,))
            temp_order[2] = 1
            temp_order[3] = 1
            args = numpy.argsort(temp)
            check = numpy.abs(numpy.diff(temp_order[args]))

            if numpy.sum(check) < 3:
                log.debug("Probe Angles for GP1, GP2, WF1, WF2: %.1f, %.1f,"
                          " %.1f, %.1f", gp_pair[i, 2], gp_pair[i, 3],
                          wp_pair[j, 2], wp_pair[j, 3])
                continue
            # CHECK THAT THE GP1 IS CCW OF WF1
            ang_sep = gp_pair[i, 2] - wp_pair[j, 2]
            if numpy.abs(ang_sep) >= 180.:
                a1 = gp_pair[i, 2]
                a2 = wp_pair[j, 2]
                if a1 > 180.:
                    a1 -= 360.
                if a2 > 180.:
                    a2 -= 360.
                adiff = a1 - a2
                if adiff < 0.:
                    continue
            else:
                if ang_sep < 0.:
                    continue
            # CHECK THAT THE GP2 IS CCW OF WF2
            ang_sep = gp_pair[i, 3] - wp_pair[j, 3]
            if numpy.abs(ang_sep) >= 180.:
                a1 = gp_pair[i, 2]
                a2 = wp_pair[j, 2]
                if a1 > 180.:
                    a1 -= 360.
                if a2 > 180.:
                    a2 -= 360.
                adiff = a1 - a2
                if adiff < 0:
                    continue
            else:
                if ang_sep < 0.:
                    continue
            ab = numpy.abs(gp_pair[i, 2]-wp_pair[j, 2])
            if ab >= 180.:
                ab = 360.-ab
            if ab < probemind0:
                log.debug("Probe Angles between GP1 and WF1 < 22.1 degs: %.1f",
                          gp_pair[i, 2]-wp_pair[j, 2])
                continue
            ab = numpy.abs(gp_pair[i, 3]-wp_pair[j, 2])
            if ab >= 180.:
                ab = 360.-ab
            if ab < probemind0:
                log.debug("Probe Angles between GP2 and WF1 < 22.1 degs: %.1f",
                          gp_pair[i, 3]-wp_pair[j, 2])
                continue
            ab = numpy.abs(gp_pair[i, 2]-wp_pair[j, 3])
            if ab >= 180.:
                ab = 360.-ab
            if ab < probemind0:
                log.debug("Probe Angles between GP1 and WF2 < 22.1 degs: %.1f",
                          gp_pair[i, 2]-wp_pair[j, 3])
                continue
            ab = numpy.abs(gp_pair[i, 3]-wp_pair[j, 3])
            if ab >= 180.:
                ab = 360.-ab
            if ab < probemind0:
                log.debug("Probe Angles between GP2 and WF2 < 22.1 degs: %.1f",
                          gp_pair[i, 3]-wp_pair[j, 3])
                continue

            gs1 = lguide_star1_cand[int(gp_pair[i, 0])]
            gs2 = lguide_star2_cand[int(gp_pair[i, 1])]
            wfs1 = lwfs_star1_cand[int(wp_pair[j, 0])]
            wfs2 = lwfs_star2_cand[int(wp_pair[j, 1])]

            sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
            log.debug("Guideprobe/WFS solution:")
            log.debug("gp1  %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f", *tuple(sol[0][0:9]))
            log.debug("gp2  %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f", *tuple(sol[1][0:9]))
            log.debug("wfs1 %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f", *tuple(sol[2][0:9]))
            log.debug("wfs2 %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f", *tuple(sol[3][0:9]))

            log.info("Probe Angles for GP1, GP2, WF1, WF2: %.1f, %.1f, %.1f,"
                     " %.1f", gp_pair[i, 2], gp_pair[i, 3], wp_pair[j, 2],
                     wp_pair[j, 3])

            return sol


def findGuideWFSSol(ra, dec, pa, oo, config, debug=True):
    """
    Implements the logic of picking guide probe and wavefront sensor stars for
    one shot position.

    9/2015, implement new guide finding method by Karl, Sungryong, and Greg.

    In short:
        In "toLocalRV" angles, the allowable angles for each G1, G2, W1, W2,
        in the searching annulus,
        G1      15-190
        G2     200-10
        W1     290-100
        W2     110-280

        1) find all G1 and G2 stars, make pairs between G1 and G2,
            then find whether any good W1 W2 pair to support the chosen
            (trying) g1-g2 pair

    """

    # pair search version /Karl/Sungryong/Greg
    if debug:
        print(">>>> Greg/Sungryong findGuideWFSSol ...")

    probemind0 = config.getfloat("General", "probemind0")
    # probemind1 = config.getfloat("General", "probemind1")
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")
    dpatrol_g1min = config.getfloat("General", "dpatrol_g1min")
    dpatrol_g1max = config.getfloat("General", "dpatrol_g1max")
    dpatrol_g2min = config.getfloat("General", "dpatrol_g2min")
    dpatrol_g2max = config.getfloat("General", "dpatrol_g2max")
    dpatrol_w1min = config.getfloat("General", "dpatrol_w1min")
    dpatrol_w1max = config.getfloat("General", "dpatrol_w1max")
    dpatrol_w2min = config.getfloat("General", "dpatrol_w2min")
    dpatrol_w2max = config.getfloat("General", "dpatrol_w2max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = fastfindStarCand(oo, config, "guide1")
    # print "Number of suitable guidestars for probe 1:", len(guide_star1_cand)
    # find those stars which are suitable candidates for guiding in probe 2
    guide_star2_cand = fastfindStarCand(oo, config, "guide2")
    # print "Number of suitable guidestars for probe 2:", len(guide_star2_cand)
    # find those stars which are suitable candidates for WFS in probe 1
    wfs_star1_cand = fastfindStarCand(oo, config, "wfs1")
    # print "Number of suitable wfsstars for probe 1:", len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = fastfindStarCand(oo, config, "wfs2")
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec, pa, guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec, pa, guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec, pa, wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec, pa, wfs_star2_cand)
    # pre-select all stars that fall into probe patrol region
    if dpatrol_g1min < dpatrol_g1max:
        ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
            * (lguide_star1_cand[:, 2] <= dpatrol_max/2.) \
            * (lguide_star1_cand[:, 3] >= dpatrol_g1min) \
            * (lguide_star1_cand[:, 3] <= dpatrol_g1max)
    else:
        ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
            * (lguide_star1_cand[:, 2] <= dpatrol_max/2.) \
            * ((lguide_star1_cand[:, 3] >= dpatrol_g1min)
               + (lguide_star1_cand[:, 3] <= dpatrol_g1max))
    if debug:
        print("Number of guide_star1 stars in patrol region:", len(ii[ii]),
              lguide_star1_cand[ii, 2], dpatrol_min/2., dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]
    else:
        if debug:
            print("Insufficients stars in guide_star1 patrol region.")
        return

    if dpatrol_g2min < dpatrol_g2max:
        ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
            * (lguide_star2_cand[:, 2] <= dpatrol_max/2.) \
            * (lguide_star2_cand[:, 3] >= dpatrol_g2min) \
            * (lguide_star2_cand[:, 3] <= dpatrol_g2max)
    else:
        ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
            * (lguide_star2_cand[:, 2] <= dpatrol_max/2.) \
            * ((lguide_star2_cand[:, 3] >= dpatrol_g2min)
               + (lguide_star2_cand[:, 3] <= dpatrol_g2max))
    if debug:
        print("Number of guide_star2 stars in patrol region:", len(ii[ii]),
              lguide_star2_cand[ii, 2], dpatrol_min/2., dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]
    else:
        if debug:
            print("Insufficients stars in guide_star2 patrol region.")
        return

    if dpatrol_w1min < dpatrol_w1max:
        ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
            * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.) \
            * (lwfs_star1_cand[:, 3] >= dpatrol_w1min) \
            * (lwfs_star1_cand[:, 3] <= dpatrol_w1max)
    else:
        ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
            * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.) \
            * ((lwfs_star1_cand[:, 3] >= dpatrol_w1min)
               + (lwfs_star1_cand[:, 3] <= dpatrol_w1max))
    if debug:
        print("Number of wfs1 stars in patrol region:", len(ii[ii]),
              lwfs_star1_cand[ii, 2], dpatrol_min/2., dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]
    else:
        if debug:
            print("Insufficients stars in wfs_star1 patrol region.")
        return

    if dpatrol_w2min < dpatrol_w2max:
        ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
            * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.) \
            * (lwfs_star2_cand[:, 3] >= dpatrol_w2min) \
            * (lwfs_star2_cand[:, 3] <= dpatrol_w2max)
    else:
        ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
            * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.) \
            * ((lwfs_star2_cand[:, 3] >= dpatrol_w2min)
               + (lwfs_star2_cand[:, 3] <= dpatrol_w2max))
    if debug:
        print("Number of wfs2 stars in patrol region:", len(ii[ii]),
              lwfs_star2_cand[ii, 2], dpatrol_min/2., dpatrol_max/2.)
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]
    else:
        if debug:
            print("Insufficients stars in wfs_star2 patrol region.")
        return

    # if len(lguide_star1_cand) < 1 or len(lguide_star2_cand) < 1 \
    #     or len(lwfs_star1_cand) < 1 or len(wfs_star2_cand) < 1:
    #         if debug: print "Insufficient number of stars in patrol region."
    #         return None

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    #######################
    # Force gp1, gp2, wfs1, and/or wfs2 if SDSS ID is given
    #######################
    if config.getint("General", "force_gp1") > 0:
        lguide_star1_cand = force(config, "force_gp1", lguide_star1_cand)
    if config.getint("General", "force_gp2") > 0:
        lguide_star2_cand = force(config, "force_gp2", lguide_star2_cand)
    if config.getint("General", "force_wfs1") > 0:
        lwfs_star1_cand = force(config, "force_wfs1", lwfs_star1_cand)
    if config.getint("General", "force_wfs2") > 0:
        lwfs_star2_cand = force(config, "force_wfs2", lwfs_star2_cand)

    gl1 = len(lguide_star1_cand)
    gl2 = len(lguide_star2_cand)
    gp_pair = numpy.zeros((gl1*gl2, 6))
    igp = 0
    for i in range(gl1):
        for j in range(gl2):
            gp_pair[igp, 0] = i
            gp_pair[igp, 1] = j
            gp_pair[igp, 2] = lguide_star1_cand[i, 3]
            gp_pair[igp, 3] = lguide_star2_cand[j, 3]
            gp_pair[igp, 4] = (lguide_star1_cand[i, gppickcol]
                               + lguide_star2_cand[j, gppickcol])/2.
            ab = numpy.abs(lguide_star1_cand[i, 3]
                           - lguide_star2_cand[j, 3])
            if ab >= 180.:
                ab = 360. - ab
            gp_pair[igp, 5] = ab
            igp = igp + 1
    sortind = numpy.argsort(gp_pair[:, 4])
    # print(gp_pair[sortind])
    gp_pair = gp_pair[sortind]
    gp_pair = gp_pair[gp_pair[:, 5] > probemind0]

    wl1 = len(lwfs_star1_cand)
    wl2 = len(lwfs_star2_cand)
    wp_pair = numpy.zeros((wl1*wl2, 6))
    iwp = 0
    for i in range(wl1):
        for j in range(wl2):
            wp_pair[iwp, 0] = i
            wp_pair[iwp, 1] = j
            wp_pair[iwp, 2] = lwfs_star1_cand[i, 3]
            wp_pair[iwp, 3] = lwfs_star2_cand[j, 3]
            wp_pair[iwp, 4] = (lwfs_star1_cand[i, wfspickcol]
                               + lwfs_star2_cand[j, wfspickcol])/2.
            ab = numpy.abs(lwfs_star1_cand[i, 3] - lwfs_star2_cand[j, 3])
            if ab >= 180.:
                ab = 360. - ab
            wp_pair[iwp, 5] = ab
            iwp = iwp + 1
    sortind = numpy.argsort(wp_pair[:, 4])
    wp_pair = wp_pair[sortind]
    wp_pair = wp_pair[wp_pair[:, 5] > probemind0]

    for i in range(len(gp_pair)):
        for j in range(len(wp_pair)):
            dup = (gp_pair[i, 2]-wp_pair[j, 2])*(gp_pair[i, 2]-wp_pair[j, 3]) \
                * (gp_pair[i, 3]-wp_pair[j, 2])*(gp_pair[i, 3]-wp_pair[j, 3])
            if dup == 0:
                # if debug: print "Probe Angles for GP1, GP2, WF1, WF2: %.1f,
                # %.1f, %.1f, %.1f" %
                # (gp_pair[i,2],gp_pair[i,3],wp_pair[j,2],wp_pair[j,3])
                continue
            temp = numpy.zeros((4,))
            temp[0:2] = gp_pair[i, 2:4]
            temp[2:4] = wp_pair[j, 2:4]
            temp_order = numpy.zeros((4,))
            temp_order[2] = 1
            temp_order[3] = 1
            args = numpy.argsort(temp)
            # temp_ang = temp[args]
            # print(temp_ang)
            check = numpy.abs(numpy.diff(temp_order[args]))
            if numpy.sum(check) < 3:
                if debug:
                    print("Probe Angles for GP1, GP2, WF1, WF2: %.1f, %.1f,"
                          " %.1f, %.1f" % (gp_pair[i, 2], gp_pair[i, 3],
                                           wp_pair[j, 2], wp_pair[j, 3]))
                continue
            if numpy.abs(gp_pair[i, 2]-wp_pair[j, 2]) < probemind0:
                if debug:
                    print("Probe Angles between GP1 and WF1 < 22.1 degs: %.1f"
                          % (gp_pair[i, 2]-wp_pair[j, 2]))
                continue
            if numpy.abs(gp_pair[i, 3]-wp_pair[j, 2]) < probemind0:
                if debug:
                    print("Probe Angles between GP2 and WF1 < 22.1 degs: %.1f"
                          % (gp_pair[i, 3]-wp_pair[j, 2]))
                continue
            if numpy.abs(gp_pair[i, 2]-wp_pair[j, 3]) < probemind0:
                if debug:
                    print("Probe Angles between GP1 and WF2 < 22.1 degs: %.1f"
                          % (gp_pair[i, 2]-wp_pair[j, 3]))
                continue
            if numpy.abs(gp_pair[i, 3]-wp_pair[j, 3]) < probemind0:
                if debug:
                    print("Probe Angles between GP2 and WF2 < 22.1 degs: %.1f"
                          % (gp_pair[i, 3]-wp_pair[j, 3]))
                continue

            gs1 = lguide_star1_cand[gp_pair[i, 0]]
            gs2 = lguide_star2_cand[gp_pair[i, 1]]
            wfs1 = lwfs_star1_cand[wp_pair[j, 0]]
            wfs2 = lwfs_star2_cand[wp_pair[j, 1]]

            sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
            if debug:
                print("Guideprobe/WFS solution:")
                print("gp1  %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[0][0:9]))
                print("gp2  %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[1][0:9]))
                print("wfs1 %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[2][0:9]))
                print("wfs2 %d:%d RA, DEC = %f,%f"
                      " mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[3][0:9]))

            return sol

    # pick gp star 1 (pick uband brightest)
    # gs1 = lguide_star1_cand[ numpy.argmin( lguide_star1_cand[:,gppickcol]) ]
    # #print "gs1 theta = %.1f" % gs1[3], ": %d:%d" % tuple(gs1[0:2])
    # # calc angular distances to all guide star 2, wfs1 and wfs2  candidates
    # aad_gp2 =  aDist(gs1, lguide_star2_cand)
    # aad_wfs1 = aDist(gs1, lwfs_star1_cand)
    # aad_wfs2 = aDist(gs1, lwfs_star2_cand)
    #
    # if debug: print "Possible  thetas of gs2 stars:", lguide_star2_cand[:,3]
    # if debug: print "Possible Dthetas of gs2 stars:", aad_gp2[:]
    # # loop over list of guide star 2 candidates by ascendending angular distance
    # ll = numpy.argsort( -abs(aad_gp2) )
    #
    # for i in range( len(ll) ):
    #     l = ll[i]
    #     gs2 = lguide_star2_cand[l]
    #     if debug: print "gs2 theta = %.1f, dtheta = %.1f" % (gs2[3], aad_gp2[l]), ": %d:%d" % tuple(gs2[0:2])
    #     # check if for this combination of there exist wfs stars on either side
    #     # of the focal plane between those stars ... i.e. such that they can be
    #     # arranged on a circle in the order gp1 - wfs1 - gp2 - wfs2
    #     jj = (aad_wfs1 % 360.) < (aad_gp2[l] % 360. - 5.)
    #     jj *= (aad_wfs1 % 360.) > 5.
    #     if debug: print "Possible  thetas of wfs1 stars:", lwfs_star1_cand[jj,3]
    #     if debug: print "Possible Dthetas of wfs1 stars:", aad_wfs1[jj]
    #     kk  = ( (aad_wfs2 % 360.) > (aad_gp2[l] % 360. + 5.) )
    #     kk *= ( (aad_wfs2 % 360.) < 355. )
    #     if debug: print "Possible  thetas of wfs2 stars:", lwfs_star2_cand[kk,3]
    #     if debug: print "Possible Dthetas of wfs2 stars:", aad_wfs2[kk]
    #     if not any(jj) or not any(kk):
    #         # if not wfs1 and wfs2 stars exist between gp1 and gp2
    #         #  choose next possible gp2
    #         continue
    #     # if wfs1 and wfs2 stars do exist between gp1 and gp2
    #     # pick uband brightest one and call this problem solved!
    #     wfs1 = lwfs_star1_cand[jj][ numpy.argmin(lwfs_star1_cand[jj,wfspickcol]) ]
    #     wfs2 = lwfs_star2_cand[kk][ numpy.argmin(lwfs_star2_cand[kk,wfspickcol]) ]
    #     if debug: print "wfs1 theta = %.1f " % (wfs1[3])
    #     if debug: print "wfs2 theta = %.1f " % (wfs2[3])
    #
    #     sol = rvToGlobal(ra,dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]) )
    #     if debug:
    #         print "Guideprobe/WFS solution:"
    #         print "gp1  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[0][0:9])
    #         print "gp2  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[1][0:9])
    #         print "wfs1 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[2][0:9])
    #         print "wfs2 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f" % tuple(sol[3][0:9])
    #
    #     return sol

    return


def oldfindGuideWFSSol(ra, dec, pa, oo, config, debug=True):
    """
    Implements the logic of picking guide probe and wavefront sensor stars
    for one shot position.

    In short:
    1) reject all stars that do not fall into the magnitude limits.
    2) reject all stars that do not fall into the patrol region.
    3) pick brightest guide probe 1 candidate.
       (Band can be defined in config file)
    4) pick the guide probe 2 candidate which has the larges angular distance
    5) check if there are suitable wavefront sensor stars on either side of
       the two guide probes if not pick next brightest gp 2 candidate
       and repeat
    6) pick the two brightest suitable wavefront sensor stars
       (Band can be defined in config file).

    Note: We currently enforce a minimum angular distance (in focal plane)
    of 5 deg.
    """
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = findStarCand(oo, config, "guide1")
    # print "Number of suitable guidestars for probe 1:", len(guide_star1_cand)
    # find those stars which are suitable candidates for guiding in probe 2
    guide_star2_cand = findStarCand(oo, config, "guide2")
    # print "Number of suitable guidestars for probe 2:", len(guide_star2_cand)
    # find those stars which are suitable candidates for WFS in probe 1
    wfs_star1_cand = findStarCand(oo, config, "wfs1")
    # print "Number of suitable wfsstars for probe 1:", len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = findStarCand(oo, config, "wfs2")
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec, pa, guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec, pa, guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec, pa, wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec, pa, wfs_star2_cand)

    # pre-select all stars that fall into probe patrol region
    ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]

    ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]

    ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]

    ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]

    if len(lguide_star1_cand) < 1 or len(lguide_star2_cand) < 1 \
       or len(lwfs_star1_cand) < 1 or len(wfs_star2_cand) < 1:
        if debug:
            print("Insufficient number of stars in patrol region.")
        return None

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    # pick gp star 1 (pick uband brightest)
    gs1 = lguide_star1_cand[numpy.argmin(lguide_star1_cand[:, gppickcol])]
    # print "gs1 theta = %.1f" % gs1[3], ": %d:%d" % tuple(gs1[0:2])
    # calc angular distances to all guide star 2, wfs1 and wfs2  candidates
    aad_gp2 = aDist(gs1, lguide_star2_cand)
    aad_wfs1 = aDist(gs1, lwfs_star1_cand)
    aad_wfs2 = aDist(gs1, lwfs_star2_cand)

    if debug:
        print("Possible  thetas of gs2 stars:", lguide_star2_cand[:, 3])
        print("Possible Dthetas of gs2 stars:", aad_gp2[:])
    # loop over list of guide star 2 candidates by ascendending
    # angular distance
    ll = numpy.argsort(-abs(aad_gp2))

    for i in range(len(ll)):
        l = ll[i]
        gs2 = lguide_star2_cand[l]
        if debug:
            print("gs2 theta = %.1f, dtheta = %.1f" % (gs2[3],
                                                       aad_gp2[l]),
                  ": %d:%d" % tuple(gs2[0:2]))
        # check if for this combination of there exist wfs stars on either side
        # of the focal plane between those stars ... i.e. such that they can be
        # arranged on a circle in the order gp1 - wfs1 - gp2 - wfs2
        jj = (aad_wfs1 % 360.) < (aad_gp2[l] % 360. - 5.)
        jj *= (aad_wfs1 % 360.) > 5.
        if debug:
            print("Possible  thetas of wfs1 stars:", lwfs_star1_cand[jj, 3])
            print("Possible Dthetas of wfs1 stars:", aad_wfs1[jj])
        kk = ((aad_wfs2 % 360.) > (aad_gp2[l] % 360. + 5.))
        kk *= ((aad_wfs2 % 360.) < 355.)
        if debug:
            print("Possible  thetas of wfs2 stars:", lwfs_star2_cand[kk, 3])
            print("Possible Dthetas of wfs2 stars:", aad_wfs2[kk])
        if not any(jj) or not any(kk):
            # if not wfs1 and wfs2 stars exist between gp1 and gp2
            #  choose next possible gp2
            continue
        # if wfs1 and wfs2 stars do exist between gp1 and gp2
        # pick uband brightest one and call this problem solved!
        wfs1 = lwfs_star1_cand[jj][numpy.argmin(lwfs_star1_cand[jj,
                                                                wfspickcol])]
        wfs2 = lwfs_star2_cand[kk][numpy.argmin(lwfs_star2_cand[kk,
                                                                wfspickcol])]
        if debug:
            print("wfs1 theta = %.1f " % (wfs1[3]))
            print("wfs2 theta = %.1f " % (wfs2[3]))

        sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
        if debug:
            print("Guideprobe/WFS solution:")
            print("gp1  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[0][0:9]))
            print("gp2  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[1][0:9]))
            print("wfs1 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[2][0:9]))
            print("wfs2 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[3][0:9]))

        return sol

    return None


def findGuideWFSSol20(ra, dec, pa, oo, config, debug=True):
    """
    Implements the logic of picking guide probe and wavefront sensor
    stars for one shot position.

    In short:
    1) reject all stars that do not fall into the magnitude limits.
    2) reject all stars that do not fall into the patrol region.
    3) pick brightest guide probe 1 candidate. (Band can be defined in
       config file)
    4) pick the guide probe 2 candidate which has the larges angular distance
    5) check if there are suitable wavefront sensor stars on either side of
       the two guide probes if not pick next brightest gp 2 candidate
       and repeat
    6) pick the two brightest suitable wavefront sensor stars
       (Band can be defined in config file).

    Note: We currently enforce a minimum angular distance (in focal plane)
    of 5 deg.
    """
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = findStarCand20(oo, config, "guide1")
    # print "Number of suitable guidestars for probe 1:", len(guide_star1_cand)
    # find those stars which are suitable candidates for guiding in probe 2
    guide_star2_cand = findStarCand20(oo, config, "guide2")
    # print "Number of suitable guidestars for probe 2:", len(guide_star2_cand)
    # find those stars which are suitable candidates for WFS in probe 1
    wfs_star1_cand = findStarCand20(oo, config, "wfs1")
    # print "Number of suitable wfsstars for probe 1:", len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = findStarCand20(oo, config, "wfs2")
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec, pa, guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec, pa, guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec, pa, wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec, pa, wfs_star2_cand)

    # pre-select all stars that fall into probe patrol region
    ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]

    ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]

    ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]

    ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]

    if len(lguide_star1_cand) < 1 or len(lguide_star2_cand) < 1 \
       or len(lwfs_star1_cand) < 1 or len(wfs_star2_cand) < 1:
        if debug:
            print("Insufficient number of stars in patrol region.")
        return

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    # pick gp star 1 (pick uband brightest)
    gs1 = lguide_star1_cand[numpy.argmin(lguide_star1_cand[:, gppickcol])]
    # print "gs1 theta = %.1f" % gs1[3], ": %d:%d" % tuple(gs1[0:2])
    # calc angular distances to all guide star 2, wfs1 and wfs2  candidates
    aad_gp2 = aDist(gs1, lguide_star2_cand)
    aad_wfs1 = aDist(gs1, lwfs_star1_cand)
    aad_wfs2 = aDist(gs1, lwfs_star2_cand)

    if debug:
        print("Possible  thetas of gs2 stars:", lguide_star2_cand[:, 3])
        print("Possible Dthetas of gs2 stars:", aad_gp2[:])
    # loop over list of guide star 2 candidates by ascendending
    #  angular distance
    ll = numpy.argsort(-abs(aad_gp2))

    for i in range(len(ll)):
        l = ll[i]
        gs2 = lguide_star2_cand[l]
        if debug:
            print("gs2 theta = %.1f, dtheta = %.1f" % (gs2[3], aad_gp2[l]),
                  ": %d:%d" % tuple(gs2[0:2]))
        # check if for this combination of there exist wfs stars on either side
        # of the focal plane between those stars ... i.e. such that they can be
        # arranged on a circle in the order gp1 - wfs1 - gp2 - wfs2
        jj = (aad_wfs1 % 360.) < (aad_gp2[l] % 360. - 5.)
        jj *= (aad_wfs1 % 360.) > 5.
        if debug:
            print("Possible  thetas of wfs1 stars:", lwfs_star1_cand[jj, 3])
            print("Possible Dthetas of wfs1 stars:", aad_wfs1[jj])
        kk = ((aad_wfs2 % 360.) > (aad_gp2[l] % 360. + 5.))
        kk *= ((aad_wfs2 % 360.) < 355.)
        if debug:
            print("Possible  thetas of wfs2 stars:", lwfs_star2_cand[kk, 3])
            print("Possible Dthetas of wfs2 stars:", aad_wfs2[kk])
        if not any(jj) or not any(kk):
            # if not wfs1 and wfs2 stars exist between gp1 and gp2
            #  choose next possible gp2
            continue
        # if wfs1 and wfs2 stars do exist between gp1 and gp2
        # pick uband brightest one and call this problem solved!
        wfs1 = lwfs_star1_cand[jj][numpy.argmin(lwfs_star1_cand[jj,
                                                                wfspickcol])]
        wfs2 = lwfs_star2_cand[kk][numpy.argmin(lwfs_star2_cand[kk,
                                                                wfspickcol])]
        if debug:
            print("wfs1 theta = %.1f " % (wfs1[3]))
            print("wfs2 theta = %.1f " % (wfs2[3]))

        sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
        if debug:
            print("Guideprobe/WFS solution:")
            print("gp1  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[0][0:9]))
            print("gp2  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[1][0:9]))
            print("wfs1 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[2][0:9]))
            print("wfs2 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[3][0:9]))

        return sol

    return None


def findGuideWFSSol21(ra, dec, pa, oo, config, debug=True):
    """
    Implements the logic of picking guide probe and wavefront
    sensor stars for one shot position.

    In short:
    1) reject all stars that do not fall into the magnitude limits.
    2) reject all stars that do not fall into the patrol region.
    3) pick brightest guide probe 1 candidate.
       (Band can be defined in config file)
    4) pick the guide probe 2 candidate which has the larges
       angular distance
    5) check if there are suitable wavefront sensor stars on either
       side of the two guide probes if not pick next brightest gp 2
       candidate and repeat
    6) pick the two brightest suitable wavefront sensor stars
       (Band can be defined in config file).

    Note: We currently enforce a minimum angular distance (in focal plane)
          of 5 deg.
    """
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = findStarCand21(oo, config, "guide1")
    # print "Number of suitable guidestars for probe 1:", len(guide_star1_cand)
    # find those stars which are suitable candidates for guiding in probe 2
    guide_star2_cand = findStarCand21(oo, config, "guide2")
    # print "Number of suitable guidestars for probe 2:", len(guide_star2_cand)
    # find those stars which are suitable candidates for WFS in probe 1
    wfs_star1_cand = findStarCand21(oo, config, "wfs1")
    # print "Number of suitable wfsstars for probe 1:", len(wfs_star1_cand)
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = findStarCand21(oo, config, "wfs2")
    # print "Number of suitable wfsstars for probe 2:", len(wfs_star1_cand)

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec, pa, guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec, pa, guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec, pa, wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec, pa, wfs_star2_cand)

    # pre-select all stars that fall into probe patrol region
    ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]

    ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]

    ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]

    ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]

    if len(lguide_star1_cand) < 1 or len(lguide_star2_cand) < 1 \
       or len(lwfs_star1_cand) < 1 or len(wfs_star2_cand) < 1:
        if debug:
            print("Insufficient number of stars in patrol region.")
        return None

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    # pick gp star 1 (pick uband brightest)
    gs1 = lguide_star1_cand[numpy.argmin(lguide_star1_cand[:, gppickcol])]
    # print "gs1 theta = %.1f" % gs1[3], ": %d:%d" % tuple(gs1[0:2])
    # calc angular distances to all guide star 2, wfs1 and wfs2  candidates
    aad_gp2 = aDist(gs1, lguide_star2_cand)
    aad_wfs1 = aDist(gs1, lwfs_star1_cand)
    aad_wfs2 = aDist(gs1, lwfs_star2_cand)

    if debug:
        print("Possible  thetas of gs2 stars:", lguide_star2_cand[:, 3])
        print("Possible Dthetas of gs2 stars:", aad_gp2[:])
    # loop over list of guide star 2 candidates by ascendending
    # angular distance
    ll = numpy.argsort(-abs(aad_gp2))

    for i in range(len(ll)):
        l = ll[i]
        gs2 = lguide_star2_cand[l]
        if debug:
            print("gs2 theta = %.1f, dtheta = %.1f" % (gs2[3], aad_gp2[l]),
                  ": %d:%d" % tuple(gs2[0:2]))
        # check if for this combination of there exist wfs stars on either side
        # of the focal plane between those stars ... i.e. such that they can be
        # arranged on a circle in the order gp1 - wfs1 - gp2 - wfs2
        jj = (aad_wfs1 % 360.) < (aad_gp2[l] % 360. - 5.)
        jj *= (aad_wfs1 % 360.) > 5.
        if debug:
            print("Possible  thetas of wfs1 stars:", lwfs_star1_cand[jj, 3])
            print("Possible Dthetas of wfs1 stars:", aad_wfs1[jj])
        kk = ((aad_wfs2 % 360.) > (aad_gp2[l] % 360. + 5.))
        kk *= ((aad_wfs2 % 360.) < 355.)
        if debug:
            print("Possible  thetas of wfs2 stars:", lwfs_star2_cand[kk, 3])
            print("Possible Dthetas of wfs2 stars:", aad_wfs2[kk])
        if not any(jj) or not any(kk):
            # if not wfs1 and wfs2 stars exist between gp1 and gp2
            #  choose next possible gp2
            continue
        # if wfs1 and wfs2 stars do exist between gp1 and gp2
        # pick uband brightest one and call this problem solved!
        wfs1 = lwfs_star1_cand[jj][numpy.argmin(lwfs_star1_cand[jj,
                                                                wfspickcol])]
        wfs2 = lwfs_star2_cand[kk][numpy.argmin(lwfs_star2_cand[kk,
                                                                wfspickcol])]
        if debug:
            print("wfs1 theta = %.1f " % (wfs1[3]))
            print("wfs2 theta = %.1f " % (wfs2[3]))

        sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
        if debug:
            print("Guideprobe/WFS solution:")
            print("gp1  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[0][0:9]))
            print("gp2  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[1][0:9]))
            print("wfs1 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[2][0:9]))
            print("wfs2 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[3][0:9]))

        return sol

    return None


def findGuideWFSSol22(ra, dec, pa, oo, config, debug=True):
    """
    Implements the logic of picking guide probe and wavefront sensor stars
    for one shot position.

    In short:
    1) reject all stars that do not fall into the magnitude limits.
    2) reject all stars that do not fall into the patrol region.
    3) pick brightest guide probe 1 candidate. (Band can be defined in
       config file)
    4) pick the guide probe 2 candidate which has the larges angular
       distance
    5) check if there are suitable wavefront sensor stars on either
       side of the two guide probes if not pick next brightest gp 2
       candidate and repeat
    6) pick the two brightest suitable wavefront sensor stars
       (Band can be defined in config file).

    Note: We currently enforce a minimum angular distance (in focal plane)
          of 5 deg.
    """
    dpatrol_min = config.getfloat("General", "dpatrol_min")
    dpatrol_max = config.getfloat("General", "dpatrol_max")

    # find those stars which are suitable candidates for guiding in probe 1
    guide_star1_cand = findStarCand22(oo, config, "guide1")
    print("Number of suitable guidestars for probe 1:", len(guide_star1_cand))
    # find those stars which are suitable candidates for guiding in probe 2
    guide_star2_cand = findStarCand22(oo, config, "guide2")
    print("Number of suitable guidestars for probe 2:", len(guide_star2_cand))
    # find those stars which are suitable candidates for WFS in probe 1
    wfs_star1_cand = findStarCand22(oo, config, "wfs1")
    print("Number of suitable wfsstars for probe 1:", len(wfs_star1_cand))
    # find those stars which are suitable candidates for WFS in probe 2
    wfs_star2_cand = findStarCand22(oo, config, "wfs2")
    print("Number of suitable wfsstars for probe 2:", len(wfs_star1_cand))

    # translate star coordinates to focal plane coordinate system
    lguide_star1_cand = toLocalRV(ra, dec, pa, guide_star1_cand)
    lguide_star2_cand = toLocalRV(ra, dec, pa, guide_star2_cand)
    lwfs_star1_cand = toLocalRV(ra, dec, pa, wfs_star1_cand)
    lwfs_star2_cand = toLocalRV(ra, dec, pa, wfs_star2_cand)

    # pre-select all stars that fall into probe patrol region
    ii = (lguide_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star1_cand = lguide_star1_cand[ii]

    ii = (lguide_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lguide_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of guide_star2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lguide_star2_cand = lguide_star2_cand[ii]

    ii = (lwfs_star1_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star1_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs1 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star1_cand = lwfs_star1_cand[ii]

    ii = (lwfs_star2_cand[:, 2] >= dpatrol_min/2.) \
        * (lwfs_star2_cand[:, 2] <= dpatrol_max/2.)
    if debug:
        print("Number of wfs2 stars in patrol region:", len(ii[ii]))
    if len(ii[ii]) > 0:
        lwfs_star2_cand = lwfs_star2_cand[ii]

    if len(lguide_star1_cand) < 1 or len(lguide_star2_cand) < 1 \
       or len(lwfs_star1_cand) < 1 or len(wfs_star2_cand) < 1:
        if debug:
            print("Insufficient number of stars in patrol region.")
        return None

    gppickcol = config.getint("General", "gppickcol") + 3
    wfspickcol = config.getint("General", "wfspickcol") + 3

    # pick gp star 1 (pick uband brightest)
    gs1 = lguide_star1_cand[numpy.argmin(lguide_star1_cand[:, gppickcol])]
    # print "gs1 theta = %.1f" % gs1[3], ": %d:%d" % tuple(gs1[0:2])
    # calc angular distances to all guide star 2, wfs1 and wfs2  candidates
    aad_gp2 = aDist(gs1, lguide_star2_cand)
    aad_wfs1 = aDist(gs1, lwfs_star1_cand)
    aad_wfs2 = aDist(gs1, lwfs_star2_cand)

    if debug:
        print("Possible  thetas of gs2 stars:", lguide_star2_cand[:, 3])
        print("Possible Dthetas of gs2 stars:", aad_gp2[:])
    # loop over list of guide star 2 candidates by
    # ascendending angular distance
    ll = numpy.argsort(-abs(aad_gp2))

    for i in range(len(ll)):
        l = ll[i]
        gs2 = lguide_star2_cand[l]
        if debug:
            print("gs2 theta = %.1f, dtheta = %.1f" % (gs2[3], aad_gp2[l]),
                  ": %d:%d" % tuple(gs2[0:2]))
        # check if for this combination of there exist wfs stars on either side
        # of the focal plane between those stars ... i.e. such that they can be
        # arranged on a circle in the order gp1 - wfs1 - gp2 - wfs2
        jj = (aad_wfs1 % 360.) < (aad_gp2[l] % 360. - 5.)
        jj *= (aad_wfs1 % 360.) > 5.
        if debug:
            print("Possible  thetas of wfs1 stars:", lwfs_star1_cand[jj, 3])
            print("Possible Dthetas of wfs1 stars:", aad_wfs1[jj])
        kk = ((aad_wfs2 % 360.) > (aad_gp2[l] % 360. + 5.))
        kk *= ((aad_wfs2 % 360.) < 355.)
        if debug:
            print("Possible  thetas of wfs2 stars:", lwfs_star2_cand[kk, 3])
            print("Possible Dthetas of wfs2 stars:", aad_wfs2[kk])
        if not any(jj) or not any(kk):
            # if not wfs1 and wfs2 stars exist between gp1 and gp2
            #  choose next possible gp2
            continue
        # if wfs1 and wfs2 stars do exist between gp1 and gp2
        # pick uband brightest one and call this problem solved!
        wfs1 = lwfs_star1_cand[jj][numpy.argmin(lwfs_star1_cand[jj,
                                                                wfspickcol])]
        wfs2 = lwfs_star2_cand[kk][numpy.argmin(lwfs_star2_cand[kk,
                                                                wfspickcol])]
        if debug:
            print("wfs1 theta = %.1f " % (wfs1[3]))
            print("wfs2 theta = %.1f " % (wfs2[3]))

        sol = rvToGlobal(ra, dec, pa, numpy.vstack([gs1, gs2, wfs1, wfs2]))
        if debug:
            print("Guideprobe/WFS solution:")
            print("gp1  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[0][0:9]))
            print("gp2  %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[1][0:9]))
            print("wfs1 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[2][0:9]))
            print("wfs2 %d:%d RA, DEC = %f,%f mags=%.2f,%.2f,%.2f,%.2f,%.2f"
                  % tuple(sol[3][0:9]))
        return sol
    return None


def hasBrightStar(ra, dec, pa, config, oo):
    boo = findStarCand(oo, config, "fplane")
    if len(boo) == 0:
        # no bright star at all
        return False
    lboo = toLocalRV(ra, dec, pa, boo)
    dfplane = config.getfloat("General", "dfplane")
    return any(lboo[:, 2] <= dfplane/2.)
