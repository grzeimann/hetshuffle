#!/usr/bin/env python
"""Development-only Legacy Survey/DSS ACAM reference renders."""

from __future__ import print_function

import argparse
import json
from pathlib import Path

import numpy

from hetdex_shuffle import visualize
from experiment_mast_ps1 import (
    _acam_geometry,
    _config_for_experiment,
    _render_acam,
    _required_fov_deg,
    _hdul_from_processed,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ra", type=float, default=320.53572192360383)
    parser.add_argument("--dec", type=float, default=49.431000955430086)
    parser.add_argument("--pa", type=float, default=72.61566216399552)
    parser.add_argument("--scale", type=float, default=1.5)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("/private/tmp/hetshuffle-acam-references"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config = _config_for_experiment(args.output_dir)
    output_wcs, output_corners, sky_corners = _acam_geometry(
        config, args.ra, args.dec, args.pa
    )
    nx = config.getint("General", "acam_x_length")
    ny = config.getint("General", "acam_y_length")
    request_center = output_wcs.wcs_pix2world(
        [[(nx - 1) / 2.0 + 30.0, (ny - 1) / 2.0 + 30.0]], 0
    )[0]
    fov_deg = _required_fov_deg(
        request_center[0], request_center[1], sky_corners, 64
    )
    output_shape = (
        config.getint("General", "acam_y_length") + 60,
        config.getint("General", "acam_x_length") + 60,
    )
    diagnostics = {
        "retrieval_center_deg": request_center.tolist(),
        "acam_wcs_crval_deg": [args.ra, args.dec],
        "pa_deg": args.pa,
        "fov_deg": fov_deg,
        "sky_corners_deg": sky_corners.tolist(),
        "true_region_pixel_corners": output_corners.tolist(),
        "references": {},
    }
    comparison = []

    for name, getter, source in [
        (
            "legacy",
            lambda: visualize.retrieve_image_legacy(
                request_center[0], request_center[1], fov_deg, yflip=True,
                scale=args.scale, return_wcs=True,
            ),
            "Legacy DR11",
        ),
        (
            "dss",
            lambda: visualize.retrieve_image_DSS(
                request_center[0], request_center[1], fov_deg, yflip=True,
                return_wcs=True,
            ),
            "DSS",
        ),
    ]:
        try:
            result = getter()
            if result[4] is None:
                raise RuntimeError("provider returned no celestial WCS")
            hdul = _hdul_from_processed(result[0], result[4])
            fits_path = args.output_dir / (name + "_reference.fits")
            image_path = args.output_dir / (name + "_acam.jpg")
            hdul.writeto(fits_path, overwrite=True)
            coverage = _render_acam(
                hdul, config, image_path, args.ra, args.dec, args.pa, source,
                output_wcs, output_shape, sky_corners,
            )
            diagnostics["references"][name] = {
                "source": source,
                "shape": list(result[0].shape),
                "fits": str(fits_path),
                "acam": str(image_path),
                "coverage": coverage,
            }
            comparison.append((source, image_path))
        except Exception as exc:
            diagnostics["references"][name] = {"error": str(exc)}

    if comparison:
        from experiment_mast_ps1 import _side_by_side
        side_by_side = args.output_dir / "reference_side_by_side.png"
        _side_by_side(comparison, side_by_side)
        diagnostics["side_by_side"] = str(side_by_side)

    (args.output_dir / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2, default=str)
    )
    print(json.dumps(diagnostics, indent=2, default=str))


if __name__ == "__main__":
    main()
