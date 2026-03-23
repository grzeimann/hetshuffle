""" Deals with local FITS catalogues of
stars for Shuffle code

Author: DJ Farrow

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


class LocalCatalogue:
    """ Class to store local catalogue information
    """
    type = ""
    tbldata = []

    # Stores the required column headings for different types of catalogue
    _required_columns_dict = {"SDSS": ['objid', 'ra', 'dec', 'u', 'g', 'r',
                                       'i', 'z']}

    @classmethod
    def read_catalogue(self, filename, debug=True, type="SDSS"):
        """
        reads in a FITS catalogue
        stores tabledata in instance

        Parameters
        ----------
        filename (string)
            the name of the FITS catalogue

        debug (bool)
            sets debug mode or not

        type (string)
            the type of catalogue

        Returns
        -------
        True if successful, otherwise False
        """

        if debug:
            print("Reading {}.....".format(filename))

        try:
            hdus = fits.open(filename, mode="readonly")
            self.tbldata = hdus[1].data
        except Exception as e:
            print("ERROR when reading the input file: {}".format(e))
            return False

        column_names = self.tbldata.columns.names

        # Check input table has required columns
        missing = []
        for col in self._required_columns_dict[type]:
            if col not in column_names:
                missing.append(col)

        if len(missing) > 0:
            print("ERROR: input table of type {type} missing these columns:\n"
                  " {cols}".format(cols=missing, type=type))
            self.tbldata = []
            return False

        return True

    @classmethod
    def query_preloaded_catalogue(self, ra, dec, boxsize, debug=False):
        """Retrieves stars from preloaded catalogue.  ( call read_catalogue
        first )

        Parameters
        ----------
        ra  = center RA of field
        dec = center DEC of field
        boxsize = determines size of box around ra/dec for which sources should
        be retrieved

        Returns
        --------
        Array of stars with format
        IDa IDb RA DEC u g r i z
        """
        if len(self.tbldata) == 0:
            print("ERROR: No viable catalogue loaded into this instance of"
                  " LocalCatalogue, run read_catalogue method first")
            return []

        # Find rows in range
        delta_ra = boxsize/np.cos(dec/180.*np.pi)/2.
        delta_dec = boxsize/2.
        mask_of_rows_in_range = (self.tbldata['ra'] < (ra + delta_ra))
        mask_of_rows_in_range &= (self.tbldata['ra'] > (ra - delta_ra))
        mask_of_rows_in_range &= (self.tbldata['dec'] < (dec + delta_dec))
        mask_of_rows_in_range &= (self.tbldata['dec'] > (dec - delta_dec))
        selected_rows = self.tbldata[mask_of_rows_in_range][:]

        # For debugging, plot selected stars and centre of input shot
        if debug:
            plt.plot(selected_rows['ra'], selected_rows['dec'], 'k.',
                     markersize=0.9, rasterized=True)
            plt.plot(ra, dec, "r*", markersize=18.0)
            plt.text(ra, dec, "Shot centre", fontsize=18, color='r')
            plt.xlabel(" RA (Deg.) ")
            plt.ylabel(" Dec (Deg.) ")
            plt.show()

        # Ugly thing with IDs in order to fit in with
        # the rest of the Shuffle code - convert int64 to
        # two floats with 8 digits each
        oo = []
        for row in selected_rows:
            IDa = float(str(row[0])[:8])
            IDb = float(str(row[0])[8:])
            oo.append([IDa, IDb] + list(row[1:]))

        return np.array(oo)
