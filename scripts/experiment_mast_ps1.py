#!/usr/bin/env python
"""Development-only MAST PS1 FITS/ACAM integration experiment.

Artifacts default to /private/tmp so this script does not populate the
repository with downloaded FITS files or rendered images.
"""

from __future__ import print_function

import argparse
import configparser
import json
import logging
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy
from astropy.io import fits

from hetdex_shuffle import visualize
from hetdex_shuffle.mast_cutouts import (
    _decode_mast_fits,
    assess_footprint_coverage,
    fetch_mast_cutout,
)


LEGACY_NORTH_SOUTH_BOUNDARY_DEG = 32.375


def _legacy_layer_for_dec(dec):
    """Return the documented Legacy layer for an experimental request.

    This selector is intentionally isolated from ``visualize.retrieve_image``:
    the production fallback remains unchanged while the layer hypothesis is
    tested against the raw service response.
    """
    if float(dec) >= LEGACY_NORTH_SOUTH_BOUNDARY_DEG:
        return "ls-dr9-north"
    return "ls-dr11-south"


def _roundtrip_metrics(wcs, shape, target_ra, target_dec):
    ny, nx = shape[:2]
    pixels = numpy.array([
        [0.0, 0.0],
        [nx / 2.0, ny / 2.0],
        [nx - 1.0, 0.0],
        [0.0, ny - 1.0],
        [nx - 1.0, ny - 1.0],
    ])
    world = wcs.wcs_pix2world(pixels, 0)
    recovered = wcs.wcs_world2pix(world, 0)
    center_world = world[1]
    dra = (center_world[0] - target_ra + 180.0) % 360.0 - 180.0
    return {
        "sample_pixels": pixels.tolist(),
        "sample_world_deg": world.tolist(),
        "max_roundtrip_pixel_error": float(
            numpy.max(numpy.abs(recovered - pixels))
        ),
        "center_offset_arcsec": [
            float(dra * 3600.0 * numpy.cos(numpy.deg2rad(target_dec))),
            float((center_world[1] - target_dec) * 3600.0),
        ],
    }


def _config_for_experiment(output_dir):
    config = configparser.ConfigParser(defaults={"xpa_method": "local"})
    config_path = Path(__file__).parents[1] / "hetdex_shuffle" / "configs" / "shuffle.cfg"
    if not config.read(config_path):
        raise RuntimeError("could not read packaged shuffle configuration")
    for option, path in {
        "images": output_dir / "images",
        "acam_images": output_dir / "acam_images",
        "ds9_regions": output_dir / "regions",
    }.items():
        config.set("directories", option, str(path))
    for path in [output_dir / "images", output_dir / "acam_images", output_dir / "regions"]:
        path.mkdir(parents=True, exist_ok=True)
    return config


def _hdul_from_processed(data, wcs):
    return fits.HDUList([fits.PrimaryHDU(data=data, header=wcs.to_header())])


def _acam_geometry(config, ra, dec, pa=0.0):
    """Return the established ACAM output WCS and true-image sky corners."""
    nx = config.getint("General", "acam_x_length")
    ny = config.getint("General", "acam_y_length")
    w_out = visualize._build_acam_wcs(
        ra,
        dec,
        pa,
        config.getfloat("offsets", "acam"),
        config.getfloat("General", "acam_x_origin"),
        config.getfloat("General", "acam_y_origin"),
        config.getfloat("General", "acam_pix_scale"),
        nx,
        ny,
    )
    pad = 30
    w_out_pad = w_out.deepcopy()
    w_out_pad.wcs.crpix = [w_out.wcs.crpix[0] + pad,
                           w_out.wcs.crpix[1] + pad]
    corners = numpy.array([
        [pad, pad],
        [pad + nx - 1, pad],
        [pad, pad + ny - 1],
        [pad + nx - 1, pad + ny - 1],
    ], dtype=float)
    return w_out_pad, corners, w_out_pad.wcs_pix2world(corners, 0)


def _required_fov_deg(ra, dec, sky_corners, margin_pixels, native_scale=0.25):
    """Derive a square native cutout FOV from the ACAM sky footprint."""
    sky_corners = numpy.asarray(sky_corners, dtype=float)
    east_arcsec = ((sky_corners[:, 0] - ra + 180.0) % 360.0 - 180.0)
    east_arcsec *= numpy.cos(numpy.deg2rad(dec)) * 3600.0
    north_arcsec = (sky_corners[:, 1] - dec) * 3600.0
    max_offset_arcsec = max(
        float(numpy.max(numpy.abs(east_arcsec))),
        float(numpy.max(numpy.abs(north_arcsec))),
    )
    # MAST requests a square cutout centered on the requested position. The
    # ACAM shuffle center is offset from the detector midpoint, so the square
    # must cover the largest absolute offset on either axis.
    return 2.0 * (max_offset_arcsec + margin_pixels * native_scale) / 3600.0


def _reproject_coverage(hdul, output_wcs, shape_out):
    """Measure finite/footprint coverage using the same reproject call."""
    image, footprint = visualize.reproject_interp(
        hdul,
        output_wcs,
        shape_out=shape_out,
        return_footprint=True,
        order="bilinear",
    )
    image = numpy.asarray(image, dtype=float)
    footprint = numpy.asarray(footprint, dtype=float)
    valid = (footprint > 0) & numpy.isfinite(image)
    return {
        "footprint_positive_fraction": float((footprint > 0).mean()),
        "finite_reprojected_fraction": float(numpy.isfinite(image).mean()),
        "valid_reprojected_fraction": float(valid.mean()),
        "nan_reprojected_fraction": float(numpy.isnan(image).mean()),
    }


def _render_acam(hdul, config, outfile, ra, dec, pa, source, output_wcs,
                 shape_out, sky_corners):
    coverage = _reproject_coverage(hdul, output_wcs, shape_out)
    started = time.monotonic()
    visualize.visualize_acam_clean(
        ra,
        dec,
        pa,
        hdul,
        [],
        [],
        [],
        config,
        str(outfile),
        img_source=source,
        target_name="mast_ps1_experiment",
    )
    coverage["acam_seconds"] = time.monotonic() - started
    coverage["source"] = source
    coverage["source_footprint"] = assess_footprint_coverage(
        _wcs_from_hdul(hdul), hdul[0].data.shape, sky_corners
    )
    return coverage


def _wcs_from_hdul(hdul):
    from astropy.wcs import WCS
    return WCS(hdul[0].header, naxis=2)


def _retrieve_reference(ra, dec, fov_deg, scale, include_dss):
    """Use Legacy first, then the established ESO DSS fallback."""
    failures = []
    try:
        result = visualize.retrieve_image_legacy(
            ra, dec, fov_deg, yflip=True, scale=scale, return_wcs=True
        )
        if result[4] is not None:
            return result, "Legacy DR11", failures
        failures.append("Legacy returned no WCS")
    except Exception as exc:
        failures.append("Legacy: {}".format(exc))

    if include_dss:
        try:
            result = visualize.retrieve_image_DSS(
                ra, dec, fov_deg, yflip=True, return_wcs=True
            )
            if result[4] is not None:
                return result, "DSS", failures
            failures.append("DSS returned no WCS")
        except Exception as exc:
            failures.append("DSS: {}".format(exc))
    raise RuntimeError("no reference image available: {}".format("; ".join(failures)))


def _side_by_side(paths, output_path):
    fig, axes = plt.subplots(1, len(paths), figsize=(8 * len(paths), 8))
    axes = numpy.atleast_1d(axes)
    for axis, (label, path) in zip(axes, paths):
        axis.imshow(plt.imread(path))
        axis.set_title(label)
        axis.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    # These are the resolved coordinates from the established J212209
    # do_shuffle command on this experimental branch.  Keep this script tied
    # to that real field; arbitrary coordinates belong in separate tests.
    parser.add_argument("--ra", type=float, default=320.53572192360383)
    parser.add_argument("--dec", type=float, default=49.431000955430086)
    parser.add_argument("--pa", type=float, default=72.61566216399552)
    parser.add_argument("--fov-deg", type=float, default=0.0)
    parser.add_argument("--scale", type=float, default=1.5)
    parser.add_argument(
        "--coverage-margin-pixels",
        type=int,
        default=64,
        help="native PS1 margin on each side of the ACAM sky footprint",
    )
    parser.add_argument(
        "--no-acam-coverage-sizing",
        action="store_true",
        help="reproduce the original nominal FOV request instead of sizing from ACAM",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/private/tmp/hetshuffle-mast-ps1"),
    )
    parser.add_argument(
        "--also-dss",
        dest="also_dss",
        action="store_true",
        default=True,
        help="render ESO DSS in addition to the Legacy attempt (default)",
    )
    parser.add_argument(
        "--no-dss",
        dest="also_dss",
        action="store_false",
        help="skip the ESO DSS request",
    )
    parser.add_argument(
        "--compare-cds",
        action="store_true",
        help="also retrieve the current CDS HiPS result and render it",
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
    derived_fov_deg = _required_fov_deg(
        request_center[0],
        request_center[1],
        sky_corners,
        args.coverage_margin_pixels,
    )
    request_fov_deg = args.fov_deg
    if not args.no_acam_coverage_sizing:
        request_fov_deg = max(request_fov_deg, derived_fov_deg)
    output_shape = (
        config.getint("General", "acam_y_length") + 60,
        config.getint("General", "acam_x_length") + 60,
    )

    raw, result, timing = fetch_mast_cutout(
        request_center[0], request_center[1], request_fov_deg, scale=args.scale
    )
    raw_path = args.output_dir / "mast_ps1_raw.fits"
    raw_path.write_bytes(raw)

    data, cd, image_wcs = _decode_mast_fits(raw, yflip=True)
    processed_hdul = _hdul_from_processed(data, image_wcs)
    processed_path = args.output_dir / "mast_ps1_yflipped.fits"
    processed_hdul.writeto(processed_path, overwrite=True)
    acam_seconds = _render_acam(
        processed_hdul,
        config,
        args.output_dir / "mast_ps1_acam.jpg",
        args.ra,
        args.dec,
        args.pa,
        "MAST/PS1/3pi/g",
        output_wcs,
        output_shape,
        sky_corners,
    )

    diagnostics = {
        "request": {
            "ra": request_center[0],
            "dec": request_center[1],
            "acam_wcs_crval": [args.ra, args.dec],
            "pa_deg": args.pa,
            "nominal_fov_deg": args.fov_deg,
            "derived_acam_fov_deg": derived_fov_deg,
            "request_fov_deg": request_fov_deg,
            "scale": args.scale,
            "coverage_margin_pixels": args.coverage_margin_pixels,
            "native_size_pixels": timing["native_size_pixels"],
            "filter": "g",
            "survey": "3pi",
            "stack_file": "data",
            "output_format": "fits",
        },
        "result": result,
        "timing": dict(timing, acam_seconds=acam_seconds["acam_seconds"]),
        "acam_footprint": {
            "true_region_pixel_corners": output_corners.tolist(),
            "true_region_sky_corners_deg": sky_corners.tolist(),
        },
        "image": {
            "shape": list(data.shape),
            "dtype": str(data.dtype),
            "pixel_scale_matrix_deg_per_pixel": numpy.asarray(
                image_wcs.pixel_scale_matrix
            ).tolist(),
            "ctype": list(image_wcs.wcs.ctype),
            "crval": list(image_wcs.wcs.crval),
            "crpix": list(image_wcs.wcs.crpix),
            "wcs_roundtrip": _roundtrip_metrics(
                image_wcs, data.shape, request_center[0], request_center[1]
            ),
        },
        "artifacts": {
            "raw_fits": str(raw_path),
            "yflipped_fits": str(processed_path),
            "mast_acam": str(args.output_dir / "mast_ps1_acam.jpg"),
        },
    }

    diagnostics["coverage"] = acam_seconds

    comparison_paths = [
        ("MAST PS1 g", args.output_dir / "mast_ps1_acam.jpg"),
    ]
    reference_failures = []
    try:
        reference, reference_source, reference_failures = _retrieve_reference(
            request_center[0], request_center[1], request_fov_deg, args.scale,
            include_dss=False
        )
        reference_data, _, reference_wcs = reference[0], reference[1], reference[4]
        reference_hdul = _hdul_from_processed(reference_data, reference_wcs)
        reference_slug = "legacy" if reference_source.startswith("Legacy") else "dss"
        reference_fits = args.output_dir / (reference_slug + "_reference.fits")
        reference_hdul.writeto(reference_fits, overwrite=True)
        reference_image = args.output_dir / (reference_slug + "_acam.jpg")
        reference_coverage = _render_acam(
            reference_hdul,
            config,
            reference_image,
            args.ra,
            args.dec,
            args.pa,
            reference_source,
            output_wcs,
            output_shape,
            sky_corners,
        )
        diagnostics["reference"] = {
            "source": reference_source,
            "shape": list(reference_data.shape),
            "fits": str(reference_fits),
            "acam": str(reference_image),
            "coverage": reference_coverage,
            "failures_before_success": reference_failures,
        }
        comparison_paths.append((reference_source, reference_image))
    except Exception as exc:
        diagnostics["reference"] = {
            "error": str(exc),
            "failures": reference_failures,
        }

    if args.also_dss:
        try:
            dss = visualize.retrieve_image_DSS(
                request_center[0], request_center[1], request_fov_deg,
                yflip=True, return_wcs=True
            )
            if dss[4] is None:
                raise RuntimeError("DSS returned no WCS")
            dss_hdul = _hdul_from_processed(dss[0], dss[4])
            dss_fits = args.output_dir / "dss_reference.fits"
            dss_hdul.writeto(dss_fits, overwrite=True)
            dss_image = args.output_dir / "dss_acam.jpg"
            dss_coverage = _render_acam(
                dss_hdul, config, dss_image, args.ra, args.dec, args.pa, "DSS",
                output_wcs, output_shape, sky_corners,
            )
            diagnostics["dss"] = {
                "fits": str(dss_fits),
                "acam": str(dss_image),
                "coverage": dss_coverage,
            }
            comparison_paths.append(("DSS", dss_image))
        except Exception as exc:
            diagnostics["dss"] = {"error": str(exc)}

    side_by_side = args.output_dir / "mast_reference_side_by_side.png"
    _side_by_side(comparison_paths, side_by_side)
    diagnostics["side_by_side"] = str(side_by_side)

    if args.compare_cds:
        cds_started = time.monotonic()
        cds = visualize.retrieve_image_hips(
            request_center[0],
            request_center[1],
            request_fov_deg,
            yflip=True,
            scale=args.scale,
            return_wcs=True,
        )
        cds_hdul = _hdul_from_processed(cds[0], cds[4]) if cds[4] is not None else None
        if cds_hdul is not None:
            cds_path = args.output_dir / "cds_acam.jpg"
            cds_acam_seconds = _render_acam(
                cds_hdul, config, cds_path, args.ra, args.dec, args.pa, cds[3],
                output_wcs, output_shape, sky_corners,
            )
            diagnostics["cds_comparison"] = {
                "source": cds[3],
                "shape": list(cds[0].shape),
                "elapsed_seconds": time.monotonic() - cds_started,
                "acam_seconds": cds_acam_seconds,
                "artifact": str(cds_path),
            }
        else:
            diagnostics["cds_comparison"] = {
                "source": cds[3],
                "error": "CDS result did not provide a WCS",
            }

    diagnostics_path = args.output_dir / "diagnostics.json"
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2, default=str))
    print(json.dumps(diagnostics, indent=2, default=str))


if __name__ == "__main__":
    main()
