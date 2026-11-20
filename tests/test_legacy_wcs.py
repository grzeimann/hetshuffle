from io import BytesIO

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from hetdex_shuffle import visualize


class _Response:
    url = "https://example.invalid/viewer/fits-cutout?test=1"

    def __init__(self, content):
        self.content = content

    def raise_for_status(self):
        return None


def _legacy_fits_bytes():
    data = np.arange(20, dtype=np.float32).reshape(4, 5)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [214.5, 25.1]
    wcs.wcs.crpix = [3.0, 2.0]
    wcs.wcs.cd = [[-0.0002, 0.00003], [0.00004, 0.0002]]
    buffer = BytesIO()
    header = wcs.to_header()
    # The production Legacy getter reads the service's explicit CD cards.
    # Make the synthetic response match that format rather than Astropy's
    # equivalent PC/CDELT serialization.
    for key in ("PC1_1", "PC1_2", "PC2_1", "PC2_2", "CDELT1", "CDELT2"):
        header.pop(key, None)
    header["CD1_1"] = -0.0002
    header["CD1_2"] = 0.00003
    header["CD2_1"] = 0.00004
    header["CD2_2"] = 0.0002
    fits.PrimaryHDU(data=data, header=header).writeto(buffer)
    return buffer.getvalue(), data, wcs


def test_legacy_yflip_keeps_array_and_celestial_wcs_consistent(monkeypatch):
    content, original, original_wcs = _legacy_fits_bytes()
    monkeypatch.setattr(
        visualize.requests,
        "get",
        lambda *args, **kwargs: _Response(content),
    )

    result = visualize.retrieve_image_legacy(
        214.5,
        25.1,
        4 * 1.5 / 3600.0,
        yflip=True,
        scale=1.5,
        return_wcs=True,
    )
    image, cd, url, source, flipped_wcs = result

    np.testing.assert_array_equal(image, original[::-1, :])
    assert source == "Legacy DR11"
    assert url == _Response.url

    # A raw pixel and its row-reversed counterpart must describe the same
    # sky position after the getter's array-Y transformation.
    raw_pixel = np.array([[1.0, 1.0]])
    flipped_pixel = np.array([[1.0, original.shape[0] - 1.0 - 1.0]])
    np.testing.assert_allclose(
        original_wcs.wcs_pix2world(raw_pixel, 0),
        flipped_wcs.wcs_pix2world(flipped_pixel, 0),
        rtol=0,
        atol=1e-10,
    )

    # The second pixel-axis column changes sign, including the off-diagonal
    # term, and the reference pixel is reflected about the array height.
    expected_cd = np.asarray(original_wcs.pixel_scale_matrix, dtype=float)
    expected_cd[:, 1] *= -1.0
    np.testing.assert_allclose(
        np.asarray(cd, dtype=float), expected_cd, rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(
        flipped_wcs.wcs.crpix,
        [original_wcs.wcs.crpix[0],
         original.shape[0] + 1.0 - original_wcs.wcs.crpix[1]],
        rtol=0,
        atol=1e-12,
    )
