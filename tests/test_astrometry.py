from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from hetdex_shuffle.astrometry import FPlane, FPlaneCentroids, TangentPlane
from hetdex_shuffle.visualize import _build_acam_wcs, _sky_to_acam


def test_tangent_plane_round_trip_with_rotation():
    tangent = TangentPlane(180.0, 30.0, 23.0, 1.0, 1.0)
    ra = np.array([180.001, 179.998])
    dec = np.array([30.002, 29.999])

    x_arcsec, y_arcsec = tangent.raDec2xy(ra, dec)
    roundtrip_ra, roundtrip_dec = tangent.xy2raDec(x_arcsec, y_arcsec)

    assert np.allclose(roundtrip_ra, ra, atol=1e-9)
    assert np.allclose(roundtrip_dec, dec, atol=1e-9)


def test_acam_wcs_maps_shuffle_center_to_configured_origin():
    wcs = _build_acam_wcs(180.0, 30.0, 17.0, 2.128, 179.9, 388.96,
                          0.2709, 775, 773)
    x, y = _sky_to_acam(180.0, 30.0, 180.0, 30.0, 17.0, 2.128,
                        0.2709, 179.9, 388.96, wcs=wcs)

    assert x == pytest.approx(179.9, abs=1e-6)
    assert y == pytest.approx(388.96, abs=1e-6)


def test_fplane_centroids_normalize_ids_and_cache_results():
    centers = np.array([[0.001, 0.002], [-0.001, 0.0]])
    provider = FPlaneCentroids(180.0, 30.0, 17.0, 1.43,
                               centers, [b"058", " 059 "])

    assert provider.ids() == ["058", "059"]
    centroids = provider.sky_centroids_dict()
    assert set(centroids) == {"058", "059"}
    assert provider.sky_centroid(b"058") == pytest.approx(centroids["058"])
    assert provider.sky_centroid("missing") is None


def test_packaged_fplane_loads_known_ifu():
    fplane_path = Path(__file__).parents[1] / "hetdex_shuffle" / "configs" / "fplane.txt"
    fplane = FPlane(str(fplane_path))

    assert len(fplane.ifus) > 90
    assert fplane.by_ifuslot("058").ifuslot == "058"
