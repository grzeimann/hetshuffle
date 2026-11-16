from __future__ import annotations

import logging

import numpy as np
import pytest
from astropy.table import Table

from hetdex_shuffle import findStars
from hetdex_shuffle.astrometry import FPlane


def test_proper_motion_updates_positions_and_records_offsets():
    catalog = np.array([
        [1.0, 2.0, 180.0, 30.0, 0.0, 15.0, 16.0, 17.0, 18.0],
    ])
    gaia = np.array([
        [3.0, 4.0, 180.0005, 30.0003, 0.0, 15.2, 0.0, 0.0, 0.0,
         100.0, -50.0],
    ])

    result = findStars.update_coords(
        catalog.copy(), gaia, 2025.0, logging.getLogger("test"),
    )

    expected_ra_offset = 25.0 * 100.0 / 1e3 / 3600.0 / np.cos(np.deg2rad(30.0))
    expected_dec_offset = 25.0 * -50.0 / 1e3 / 3600.0
    assert result.shape == (1, 11)
    assert result[0, 2] == pytest.approx(180.0 + expected_ra_offset)
    assert result[0, 3] == pytest.approx(30.0 + expected_dec_offset)
    assert result[0, -2:] == pytest.approx(
        [expected_ra_offset, expected_dec_offset],
    )


def test_local_panstarrs_query_maps_rows_and_halves_search_radius(monkeypatch):
    table = Table({
        "objID": ["123456789012345678"],
        "raMean": [180.01],
        "decMean": [30.02],
        "gMeanPSFMag": [15.1],
        "rMeanPSFMag": [15.2],
        "iMeanPSFMag": [15.3],
        "zMeanPSFMag": [15.4],
    })
    call = {}

    def fake_catalog(ra, dec, radius, config):
        call.update(ra=ra, dec=dec, radius=radius, config=config)
        return table

    monkeypatch.setattr(findStars, "het_panstarrs_data", fake_catalog)
    config = object()

    result = findStars.queryPANSTARRS_LOCAL(180.0, 30.0, 0.2, config)

    assert call == {"ra": 180.0, "dec": 30.0, "radius": 0.1,
                    "config": config}
    assert result.tolist() == [[123456789.0, 12345678.0, 180.01, 30.02,
                                0.0, 15.1, 15.2, 15.3, 15.4]]


def test_fplane_excludes_slots_and_can_skip_empty_hardware(tmp_path):
    fplane_path = tmp_path / "fplane.txt"
    fplane_path.write_text(
        "# IFUSLOT X_FP Y_FP SPECID SPECSLOT IFUID IFUROT PLATESC\n"
        "000 0.0 0.0 00 900 000 0.0 1.0\n"
        "001 10.0 20.0 123 901 456 0.0 1.0\n"
    )

    populated = FPlane(str(fplane_path), exclude_ifuslot=["001"])
    assert len(populated.ifus) == 1
    assert populated.by_ifuslot("000").ifuid == "N01"
    with pytest.raises(KeyError):
        populated.by_ifuslot("001")

    skip_empty = FPlane(str(fplane_path), skip_empty=True)
    assert [ifu.ifuslot for ifu in skip_empty.ifus] == ["001"]
