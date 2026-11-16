from __future__ import annotations

import logging

import numpy as np
import pytest

from hetdex_shuffle import findStars


def test_ra_wrap_and_polar_coordinate_round_trip():
    wrapped = findStars.aDist2(np.array([359.0, 1.0]), 0.0)
    assert np.allclose(wrapped, [-1.0, 1.0])

    radius = np.array([1.0, 2.0])
    angle = np.array([45.0, 220.0])
    x, y = findStars.rv2xy(radius, angle)
    roundtrip_radius, roundtrip_angle = findStars.xy2rv(x, y)
    assert np.allclose(roundtrip_radius, radius)
    assert np.allclose(roundtrip_angle, angle % 360.0)


def test_local_global_coordinate_round_trip():
    sources = np.zeros((2, 9), dtype=float)
    sources[:, 2] = [179.99, 180.01]
    sources[:, 3] = [29.995, 30.005]

    local = findStars.toLocal(180.0, 30.0, 23.0, sources)
    global_sources = findStars.toGlobal(180.0, 30.0, 23.0, local)

    assert np.allclose(global_sources[:, 2:4], sources[:, 2:4], atol=1e-10)


def test_remove_duplicates_keeps_first_close_source_and_distant_source():
    sources = np.zeros((3, 9), dtype=float)
    sources[:, 2:4] = [[180.0, 30.0], [180.0001, 30.0], [180.01, 30.0]]

    result = findStars.removeDuplicates(sources, d=0.001)

    assert len(result) == 2
    assert np.allclose(result[:, 2], [180.0, 180.01])


def test_candidate_selection_applies_magnitude_limits(shuffle_config):
    sources = np.zeros((2, 9), dtype=float)
    sources[:, 2:4] = [[180.0, 30.0], [180.01, 30.0]]
    sources[:, 5:9] = [[15.0, 16.0, 17.0, 18.0],
                        [25.0, 25.0, 25.0, 25.0]]

    result = findStars.fastfindStarCand(sources, shuffle_config, "acam")

    assert len(result) == 1
    assert result[0, 2] == pytest.approx(180.0)


def test_acam_boundary_contains_center_but_not_distant_source(shuffle_config):
    sources = np.zeros((2, 9), dtype=float)
    sources[:, 2:4] = [[180.0, 30.0], [180.2, 30.0]]

    inside = findStars.stars_inside_acam(180.0, 30.0, sources, 0.0,
                                         shuffle_config)

    assert inside.tolist() == [True, False]


def test_empty_proper_motion_catalog_preserves_positions():
    catalog = np.array([[1.0, 2.0, 180.0, 30.0, 0.0, 15.0, 16.0, 17.0, 18.0]])
    result = findStars.update_coords(catalog.copy(), np.empty((0, 9)),
                                     2026.0, logging.getLogger("test"))

    assert np.array_equal(result, catalog)
