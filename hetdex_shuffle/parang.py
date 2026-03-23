'''parang - routine to compute the parallactic angle for an input apparent
declination.

The trajectory calculations are derived by Frank Ray in HET Technical Reports
42, 43, and 90 (copies are available at
http://puck.as.utexas.edu/HET/references).

12SEP12  Modified by M. E. Cornell from the original by Max Fabricius.  Actual
         calculations are taken from POINT's newtraj.c.
'''

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy

# Begin HET constants from POINT's HET.h, dated 09APR12

PI = 4.0*numpy.arctan(1.0)            # a constant
PHI = 30.681436               # the latitude of the HET (deg)
R2D = 180.0/PI                # degrees per radian
D2R = PI/180.0                # Radians per degree
RMAX = 8.5                    # post-WFU tracker's angular freedom (deg)
ALT = 55.055223               # the altitude of the HET (deg)
P = numpy.cos(ALT/R2D)*numpy.cos(PHI/R2D)  # cos alt * cos latitude
Q = numpy.sin(ALT/R2D)*numpy.sin(PHI/R2D)  # sin alt * sin latitude
TDE_N = numpy.arcsin(P + Q)
TDE_S = numpy.arcsin(-P + Q)
# End HET constants from POINT's HET.h


def parang(app_dec, az_input=None, verbose=False):
    if verbose:
        print("Apparent declination: %7.3f deg\n" % (app_dec))
    app_dec /= R2D

    # algorithm for best (east) azimuth taken directly from tcsmon/rdscope.c
    # note that we use the apparent declination, rather than the observed
    if ((app_dec >= TDE_S) and (app_dec <= TDE_N)):
        if az_input is None:
            az = numpy.arccos((numpy.sin(app_dec)-Q)/P)  # TR43 eq. 5
        else:
            az = az_input*D2R
        tde = numpy.arcsin(P * numpy.cos(az) + Q)    # TR43 eq. 2
        hc = numpy.arcsin((numpy.cos(ALT/R2D)*numpy.sin(az)) /
                          (-numpy.cos(tde)))  # TR43 eq. 3
        # TR43 eq. 7, parallactic angle
        parang = numpy.arccos(numpy.cos(hc)*numpy.cos(az) +
                              numpy.sin(hc)*numpy.sin(az) *
                              numpy.sin(PHI/R2D))

        if az_input is None:
            az *= R2D
            az2 = 360 - az
        else:
            if az*R2D > 180.:
                az = 360. - az*R2D
                az2 = 360. - az
            else:
                az *= R2D
                az2 = 360 - az

        parang *= R2D
        parang2 = 360 - parang
        tracks = 2
        if (verbose):
            print("There will be 2 tracks:")
            print("  East Track: Az = %5.1f has a parallactic angle of %7.3f"
                  " deg" % (az, parang))
            print("  West Track: Az = %5.1f has a parallactic angle of %7.3f"
                  " deg" % (az2, parang2))
        if (az > 58 and az < 78):
            if (verbose):
                print("WARNING: THE CCAS TOWER will obscure some of the mirror"
                      " in the EAST TRACK.")
    elif (app_dec >= (TDE_S - RMAX/R2D) and (app_dec <= TDE_S)):
        az = 180
        az2 = 180
        parang = 180
        parang2 = 180
        tracks = 1
        if (verbose):
            print("There will be 1 track:")
            print("  Az = %5.1f has a parallactic angle of %7.3f deg" %
                  (az, parang))
    elif (app_dec <= (TDE_N + RMAX/R2D) and (app_dec >= TDE_N)):
        az = 0.0
        az2 = 0.0
        parang = 0
        parang2 = 0
        tracks = 1
        if (verbose):
            print("There will be 1 track:")
            print("  Az = %5.1f has a parallactic angle of %7.3f deg" %
                  (az, parang))
    else:
        az = 0
        az2 = 0
        parang = 0
        parang2 = 0
        tracks = 0
        if (verbose):
            print("The Declination for the object is out of range for the HET")

    return tracks, parang, parang2, az, az2
