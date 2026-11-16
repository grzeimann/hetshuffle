from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from matplotlib.collections import PatchCollection

from hetdex_shuffle import visualize


def test_acam_renderer_writes_synthetic_outputs(tmp_path, shuffle_config, monkeypatch):
    calls = {}

    def fake_reproject(im_hdul, output_wcs, shape_out, **kwargs):
        calls["shape_out"] = shape_out
        return np.ones(shape_out, dtype=float), np.ones(shape_out, dtype=float)

    monkeypatch.setattr(visualize, "reproject_interp", fake_reproject)
    monkeypatch.setattr(
        visualize,
        "plotFocalPlaneLRS_acam",
        lambda *args, **kwargs: PatchCollection([]),
    )

    output = tmp_path / "synthetic_acam.png"
    image = fits.HDUList([fits.PrimaryHDU(np.ones((16, 16), dtype=float))])
    candidates = [(180.0, 30.0, "center", 15.0, 16.0, 17.0)]

    visualize.visualize_acam_clean(
        180.0,
        30.0,
        17.0,
        image,
        candidates,
        [],
        [],
        shuffle_config,
        str(output),
        target_ra=180.0,
        target_dec=30.0,
        target_name="synthetic",
    )

    assert calls["shape_out"] == (833, 835)
    assert output.is_file() and output.stat().st_size > 0
    assert output.with_suffix(".csv").is_file()
    assert Path(shuffle_config.get("directories", "ds9_regions"),
                "synthetic_acam.reg").is_file()
    assert "#RA, Dec, acam_x" in output.with_suffix(".csv").read_text()
