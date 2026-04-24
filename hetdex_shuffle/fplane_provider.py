from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from .astrometry import TangentPlane as TP

class FPlaneCentroids(object):
    """
    Authoritative provider for focal-plane IFU centroids and EN polygons.

    Builds a single TangentPlane at the shuffle center and reuses it to
    compute RA/Dec centroids for IFUs consistently across the pipeline.
    """
    def __init__(self, ra0, dec0, pa, fplane_offset_deg, ifu_centers_deg, ifu_ids):
        self.ra0 = float(ra0)
        self.dec0 = float(dec0)
        self.pa = float(pa)
        self.fplane_offset_deg = float(fplane_offset_deg)
        # Normalize ids to strings
        self.ifu_ids = [self._norm_id(x) for x in ifu_ids]
        self.ifu_centers_deg = np.asarray(ifu_centers_deg, dtype=float)
        # Single TP reused everywhere (EN plane definition)
        self._tp = TP(self.ra0, self.dec0, self.pa + self.fplane_offset_deg, 1.0, 1.0)
        # Cache dict of RA/Dec centroids
        self._sky_dict = None

    @staticmethod
    def _norm_id(x):
        try:
            s = str(x, 'utf-8') if isinstance(x, (bytes, bytearray)) else str(x)
        except Exception:
            s = str(x)
        return s.strip()

    def _build_cache(self):
        sky = {}
        for ifu_id, cen in zip(self.ifu_ids, self.ifu_centers_deg):
            x_ifu_deg = float(cen[0])
            y_ifu_deg = float(cen[1])
            ra_c, dec_c = self._tp.xy2raDec(3600.0 * x_ifu_deg, 3600.0 * y_ifu_deg)
            sky[ifu_id] = (float(ra_c), float(dec_c))
        self._sky_dict = sky

    def sky_centroids_dict(self):
        if self._sky_dict is None:
            self._build_cache()
        return dict(self._sky_dict)

    def sky_centroid(self, ifu_id):
        if self._sky_dict is None:
            self._build_cache()
        return self._sky_dict.get(self._norm_id(ifu_id))

    def ids(self):
        return list(self.ifu_ids)
