from __future__ import annotations

from io import BytesIO

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.collections import PatchCollection

from hetdex_shuffle import mast_cutouts, visualize
from scripts.experiment_mast_ps1 import (
    _acam_geometry,
    _config_for_experiment,
    _legacy_layer_for_dec,
    _required_fov_deg,
)


class FakeResponse:
    def __init__(self, status_code=200, text="", headers=None,
                 json_payload=None, content=b""):
        self.status_code = status_code
        self.text = text
        self.headers = headers or {}
        self._json_payload = json_payload
        self.content = content

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json_payload


def _fits_bytes():
    data = np.arange(30, dtype=np.float32).reshape(5, 6)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [180.0, 30.0]
    wcs.wcs.crpix = [2.0, 1.5]
    wcs.wcs.cd = [[-0.0001, 0.00002], [0.00003, 0.0001]]
    buffer = BytesIO()
    fits.PrimaryHDU(data=data, header=wcs.to_header()).writeto(buffer)
    return buffer.getvalue(), data, wcs


class SuccessfulSession:
    def __init__(self, fits_content):
        self.fits_content = fits_content
        self.phases = iter(["EXECUTING", "COMPLETED"])
        self.post_kwargs = None

    def post(self, url, **kwargs):
        self.post_kwargs = kwargs
        return FakeResponse(
            status_code=303,
            headers={"Location": "/api/v0/cutouts/async/test-job"},
        )

    def get(self, url, params=None, **kwargs):
        if url.endswith("/phase"):
            return FakeResponse(text=next(self.phases))
        if url.endswith("/results/summary"):
            return FakeResponse(json_payload={"total_jobs": 1, "completed_jobs": 1})
        if url.endswith("/results/cutouts"):
            return FakeResponse(json_payload={
                "results": [{
                    "mission": "ps1",
                    "position": [180.0, 30.0],
                    "size_px": 2048,
                    "filter": "g",
                    "mission_extras": {"survey": "3pi", "file_type": "data"},
                    "science": "/files/test.fits",
                }],
                "links": {},
            })
        if url.endswith("/files/test.fits"):
            return FakeResponse(content=self.fits_content)
        raise AssertionError(f"unexpected GET: {url}")

    def close(self):
        pass


def test_mast_success_preserves_wcs_and_return_conventions(monkeypatch):
    fits_content, original, original_wcs = _fits_bytes()
    session = SuccessfulSession(fits_content)
    monkeypatch.setattr(mast_cutouts.requests, "Session", lambda: session)
    monkeypatch.setattr(mast_cutouts, "MAST_POLL_INTERVAL_SEC", 0)

    result = mast_cutouts.retrieve_image_mast(
        180.0, 30.0, 0.01, yflip=True, scale=1.5, return_wcs=True
    )

    assert len(result) == 5
    image, cd, url, source, wcs = result
    assert image.shape == original.shape
    assert image.dtype.kind == original.dtype.kind
    assert image.dtype.itemsize == original.dtype.itemsize
    np.testing.assert_array_equal(image, original[::-1, :])
    assert source == "MAST/PS1/3pi/g"
    assert url.endswith("/files/test.fits")
    assert wcs.wcs.crval.tolist() == pytest.approx(original_wcs.wcs.crval)
    assert wcs.wcs.crpix[1] == pytest.approx(original.shape[0] + 1 - 1.5)
    assert np.asarray(cd)[:, 1].tolist() == pytest.approx(
        -np.asarray(original_wcs.pixel_scale_matrix)[:, 1]
    )
    fields = dict(session.post_kwargs["files"])
    assert fields["output_format"][1] == "fits"
    assert fields["ps1.filter"][1] == "g"
    assert fields["ps1.survey"][1] == "3pi"
    assert fields["ps1.stack_file"][1] == "data"

    session.phases = iter(["COMPLETED"])
    four_tuple = mast_cutouts.retrieve_image_mast(
        180.0, 30.0, 0.01, yflip=False, scale=1.5, return_wcs=False
    )
    assert len(four_tuple) == 4


def test_mast_timeout_is_bounded(monkeypatch):
    class NeverCompletes(SuccessfulSession):
        def get(self, url, params=None, **kwargs):
            if url.endswith("/phase"):
                return FakeResponse(text="EXECUTING")
            if url.endswith("/results/summary"):
                return FakeResponse(json_payload={})
            raise AssertionError("results should not be requested")

    fits_content, _, _ = _fits_bytes()
    session = NeverCompletes(fits_content)
    monkeypatch.setattr(mast_cutouts.requests, "Session", lambda: session)
    monkeypatch.setattr(mast_cutouts, "MAST_POLL_INTERVAL_SEC", 0)
    monkeypatch.setattr(mast_cutouts, "MAST_POLL_TIMEOUT_SEC", 0)

    with pytest.raises(TimeoutError):
        mast_cutouts.retrieve_image_mast(180.0, 30.0, 0.01)


def test_mast_failed_job_reports_service_error(monkeypatch):
    fits_content, _, _ = _fits_bytes()
    session = SuccessfulSession(fits_content)

    def failed_get(url, params=None, **kwargs):
        if url.endswith("/phase"):
            return FakeResponse(text="ERROR")
        if url.endswith("/results/summary"):
            return FakeResponse(json_payload={"failed_jobs": 1})
        if url.endswith("/error"):
            return FakeResponse(text="worker could not locate PS1 stack")
        raise AssertionError(f"unexpected GET: {url}")

    session.get = failed_get
    monkeypatch.setattr(mast_cutouts.requests, "Session", lambda: session)

    with pytest.raises(mast_cutouts.MastCutoutError, match="could not locate"):
        mast_cutouts.retrieve_image_mast(180.0, 30.0, 0.01)


def test_mast_rejects_invalid_wcs(monkeypatch):
    buffer = BytesIO()
    fits.PrimaryHDU(data=np.ones((8, 8), dtype=np.float32)).writeto(buffer)
    session = SuccessfulSession(buffer.getvalue())
    monkeypatch.setattr(mast_cutouts.requests, "Session", lambda: session)
    monkeypatch.setattr(mast_cutouts, "MAST_POLL_INTERVAL_SEC", 0)

    with pytest.raises(mast_cutouts.MastCutoutError, match="celestial WCS"):
        mast_cutouts.retrieve_image_mast(180.0, 30.0, 0.01)


def test_mast_rejects_multiple_products(monkeypatch):
    fits_content, _, _ = _fits_bytes()
    session = SuccessfulSession(fits_content)
    original_get = session.get

    def get_with_two_results(url, params=None, **kwargs):
        response = original_get(url, params=params, **kwargs)
        if url.endswith("/results/cutouts"):
            response._json_payload["results"].append(dict(response._json_payload["results"][0]))
        return response

    session.get = get_with_two_results
    monkeypatch.setattr(mast_cutouts.requests, "Session", lambda: session)
    monkeypatch.setattr(mast_cutouts, "MAST_POLL_INTERVAL_SEC", 0)

    with pytest.raises(mast_cutouts.MastCutoutError, match="exactly one"):
        mast_cutouts.retrieve_image_mast(180.0, 30.0, 0.01)


def test_mast_s3_science_uri_is_resolved_to_public_https():
    uri = "s3://stpubdata/panstarrs/test.fits"
    assert mast_cutouts._science_download_url(uri, "https://example/job") == (
        "https://stpubdata.s3.amazonaws.com/panstarrs/test.fits"
    )


def test_footprint_coverage_detects_insufficient_rotated_source_grid():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [180.0, 30.0]
    wcs.wcs.crpix = [50.5, 50.5]
    wcs.wcs.cd = [[-0.001, 0.0002], [0.0001, 0.001]]
    source_corners = wcs.wcs_pix2world(
        np.array([[0, 0], [99, 0], [0, 99], [99, 99]], dtype=float), 0
    )
    covered = mast_cutouts.assess_footprint_coverage(
        wcs, (100, 100), source_corners
    )
    assert covered["bbox_inside"]
    assert covered["minimum_edge_margin_pixels"] == pytest.approx(0.0)

    expanded = np.vstack([source_corners, [180.08, 30.08]])
    insufficient = mast_cutouts.assess_footprint_coverage(
        wcs, (100, 100), expanded
    )
    assert not insufficient["bbox_inside"]
    assert insufficient["minimum_edge_margin_pixels"] < 0


def test_mast_wcs_can_enter_acam_reprojection(shuffle_config, tmp_path, monkeypatch):
    _, _, wcs = _fits_bytes()
    data = np.ones((16, 16), dtype=float)
    image = fits.HDUList([fits.PrimaryHDU(data=data, header=wcs.to_header())])
    captured = {}

    def fake_reproject(im_hdul, output_wcs, shape_out, **kwargs):
        captured["ctype"] = [
            im_hdul[0].header["CTYPE1"],
            im_hdul[0].header["CTYPE2"],
        ]
        captured["shape_out"] = shape_out
        return np.ones(shape_out), np.ones(shape_out)

    monkeypatch.setattr(visualize, "reproject_interp", fake_reproject)
    monkeypatch.setattr(
        visualize,
        "plotFocalPlaneLRS_acam",
        lambda *args, **kwargs: PatchCollection([]),
    )
    output = tmp_path / "mast_acam.jpg"
    visualize.visualize_acam_clean(
        180.0, 30.0, 0.0, image, [], [], [], shuffle_config, str(output),
        img_source="MAST/PS1/3pi/g", target_name="mast-test",
    )
    assert captured["ctype"] == ["RA---TAN", "DEC--TAN"]
    assert captured["shape_out"] == (833, 835)
    assert output.is_file()


def test_experimental_field_geometry_preserves_pa_and_covers_acam(tmp_path):
    config = _config_for_experiment(tmp_path)
    ra = 320.53572192360383
    dec = 49.431000955430086
    pa = 72.61566216399552
    output_wcs, _, sky_corners = _acam_geometry(config, ra, dec, pa)
    nx = config.getint("General", "acam_x_length")
    ny = config.getint("General", "acam_y_length")
    request_center = output_wcs.wcs_pix2world(
        [[(nx - 1) / 2.0 + 30.0, (ny - 1) / 2.0 + 30.0]], 0
    )[0]
    fov = _required_fov_deg(request_center[0], request_center[1], sky_corners, 64)

    assert request_center[0] == pytest.approx(320.51251570593604, abs=1e-9)
    assert request_center[1] == pytest.approx(49.42711270998604, abs=1e-9)
    assert fov == pytest.approx(0.08040840757788704, abs=1e-12)
    assert output_wcs.wcs.crval.tolist() == pytest.approx([ra, dec])


def test_experimental_legacy_layer_hypothesis_follows_declination():
    assert _legacy_layer_for_dec(49.431) == "ls-dr9-north"
    assert _legacy_layer_for_dec(30.0) == "ls-dr11-south"
