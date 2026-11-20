#!/usr/bin/env python
"""Small sequential reliability/performance benchmark for image providers.

This is a development-only script.  It calls the existing retrieval functions
without changing production dispatch or timeout defaults, and writes all
results outside the repository by default.
"""

from __future__ import print_function

import argparse
import csv
import json
import math
import signal
import time
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS

from hetdex_shuffle import mast_cutouts, visualize
from hetdex_shuffle.mast_cutouts import assess_footprint_coverage
try:
    from experiment_mast_ps1 import (
        _acam_geometry,
        _config_for_experiment,
        _required_fov_deg,
    )
except ImportError:  # Importing this development script as a package module.
    from scripts.experiment_mast_ps1 import (
        _acam_geometry,
        _config_for_experiment,
        _required_fov_deg,
    )


PRIMARY_HIPS = "CDS/P/PanSTARRS/DR1/g"
PROVIDERS = ("hips2fits", "legacy", "mast", "dss")
SEQUENCES = {
    "A": ("hips2fits", "legacy", "mast", "dss"),
    "B": ("hips2fits", "legacy", "dss", "mast"),
    "C": ("hips2fits", "mast", "legacy", "dss"),
}


class BenchmarkError(RuntimeError):
    def __init__(self, category, detail, trace=None):
        super().__init__(detail)
        self.category = category
        self.trace = trace or {}


class HipsFallbackAbort(BaseException):
    """Stop the existing HiPS getter before it tries another survey."""


class AttemptTimeout(TimeoutError):
    pass


class _Deadline:
    def __init__(self, seconds):
        self.seconds = float(seconds)
        self.previous = None

    def _handler(self, signum, frame):
        raise AttemptTimeout(
            "attempt exceeded {:.1f}-second benchmark budget".format(
                self.seconds
            )
        )

    def __enter__(self):
        if not hasattr(signal, "SIGALRM"):
            raise RuntimeError("the benchmark requires Unix SIGALRM support")
        self.previous = signal.signal(signal.SIGALRM, self._handler)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, self.previous)
        return False


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _field_definitions():
    """Return established fields plus three declination-spanning fields."""
    return [
        {
            "id": "NGC5548",
            "input_ra_deg": 214.49808333333334,
            "input_dec_deg": 25.136833333333332,
            "ra_shuffle_deg": 214.5363528473741,
            "dec_shuffle_deg": 25.10984344789337,
            "pa_deg": 108.0591979935536,
            "ifu": "056",
            "selection": "Established LRS2-B NGC 5548 test field.",
        },
        {
            "id": "J212209",
            "input_ra_deg": 320.5385833333333,
            "input_dec_deg": 49.44005555555556,
            "ra_shuffle_deg": 320.53572192360383,
            "dec_shuffle_deg": 49.431000955430086,
            "pa_deg": 72.61566216399552,
            "ifu": "603",
            "selection": "Established IFU 603 development field.",
        },
        {
            "id": "DECM05_RA060",
            "input_ra_deg": 60.0,
            "input_dec_deg": -5.0,
            "ra_shuffle_deg": 60.0,
            "dec_shuffle_deg": -5.0,
            "pa_deg": 0.0,
            "ifu": "synthetic",
            "selection": "Uncrowded declination-span sample near -5 degrees.",
        },
        {
            "id": "M64",
            "input_ra_deg": 194.1821,
            "input_dec_deg": 21.6837,
            "ra_shuffle_deg": 194.1821,
            "dec_shuffle_deg": 21.6837,
            "pa_deg": 0.0,
            "ifu": "synthetic",
            "selection": "Recognizable, moderate-declination field with broad survey coverage.",
        },
        {
            "id": "NGC2403",
            "input_ra_deg": 114.214,
            "input_dec_deg": 65.603,
            "ra_shuffle_deg": 114.214,
            "dec_shuffle_deg": 65.603,
            "pa_deg": 0.0,
            "ifu": "synthetic",
            "selection": "Recognizable high-declination field near +65 degrees.",
        },
    ]


def _prepare_fields(config):
    prepared = []
    nx = config.getint("General", "acam_x_length")
    ny = config.getint("General", "acam_y_length")
    for field in _field_definitions():
        output_wcs, output_corners, sky_corners = _acam_geometry(
            config,
            field["ra_shuffle_deg"],
            field["dec_shuffle_deg"],
            field["pa_deg"],
        )
        request_center = output_wcs.wcs_pix2world(
            [[(nx - 1) / 2.0 + 30.0, (ny - 1) / 2.0 + 30.0]], 0
        )[0]
        fov_deg = _required_fov_deg(
            request_center[0], request_center[1], sky_corners, 64
        )
        item = dict(field)
        item.update({
            "request_ra_deg": float(request_center[0]),
            "request_dec_deg": float(request_center[1]),
            "fov_deg": float(fov_deg),
            "acam_sky_corners_deg": numpy.asarray(sky_corners).tolist(),
            "acam_output_corners": numpy.asarray(output_corners).tolist(),
        })
        prepared.append(item)
    return prepared


def _classify_exception(exc, trace):
    if isinstance(exc, AttemptTimeout):
        return "service_timeout"
    if isinstance(exc, BenchmarkError):
        return exc.category
    if trace.get("http_status", 0) >= 400:
        return "http_error"
    name = exc.__class__.__name__.lower()
    text = str(exc).lower()
    if "timeout" in name or "timed out" in text:
        return "service_timeout"
    if "http" in name or "status" in text:
        return "http_error"
    if any(token in name or token in text for token in (
        "connection", "dns", "nodename", "name resolution", "urlopen error",
        "max retries",
    )):
        return "connection_or_service_error"
    if "wcs" in text or "fits" in text or "image" in text:
        return "invalid_fits_or_wcs"
    if "coverage" in text or "outside" in text:
        return "incomplete_sky_coverage"
    return "connection_or_service_error"


def _validate_result(result, field, trace):
    if len(result) < 5:
        raise BenchmarkError("invalid_fits_or_wcs", "getter returned no WCS")
    data, _, url, source, image_wcs = result
    if not source:
        raise BenchmarkError(
            "service_unavailable",
            "getter returned the blank-image fallback",
            trace,
        )
    if data is None or numpy.asarray(data).ndim != 2:
        raise BenchmarkError("invalid_fits_or_wcs", "image is not 2-D", trace)
    data = numpy.asarray(data)
    finite = numpy.isfinite(data)
    if not finite.any() or not numpy.isfinite(data[finite]).all():
        raise BenchmarkError("invalid_fits_or_wcs", "image has no usable finite data", trace)
    if image_wcs is None or not image_wcs.has_celestial:
        raise BenchmarkError("missing_or_invalid_wcs", "missing celestial WCS", trace)
    coverage = assess_footprint_coverage(
        image_wcs, data.shape, field["acam_sky_corners_deg"]
    )
    if not coverage["bbox_inside"]:
        raise BenchmarkError(
            "incomplete_sky_coverage",
            "ACAM sky footprint is outside returned image: {}".format(coverage),
            trace,
        )
    trace["url"] = url
    trace["source"] = source
    trace["shape"] = list(data.shape)
    trace["dtype"] = str(data.dtype)
    trace["coverage"] = coverage
    trace["wcs"] = {
        "ctype": list(image_wcs.wcs.ctype),
        "crval_deg": image_wcs.wcs.crval.tolist(),
        "crpix": image_wcs.wcs.crpix.tolist(),
        "pixel_scale_matrix_deg_per_pixel": numpy.asarray(
            image_wcs.pixel_scale_matrix
        ).tolist(),
    }
    return data


def _call_hips_primary(field):
    trace = {"survey": PRIMARY_HIPS}
    real_get = visualize.requests.get

    def primary_only_get(url, params=None, **kwargs):
        if params and params.get("hips") != PRIMARY_HIPS:
            raise HipsFallbackAbort(
                "primary HiPS request failed; fallback survey suppressed"
            )
        try:
            response = real_get(url, params=params, **kwargs)
        except Exception as exc:
            trace["exception"] = "{}: {}".format(exc.__class__.__name__, exc)
            raise
        trace["http_status"] = getattr(response, "status_code", None)
        trace["download_bytes"] = len(getattr(response, "content", b""))
        return response

    visualize.requests.get = primary_only_get
    try:
        return visualize.retrieve_image_hips(
            field["request_ra_deg"], field["request_dec_deg"],
            field["fov_deg"], yflip=True, scale=1.5, return_wcs=True
        ), trace
    except HipsFallbackAbort as exc:
        raise BenchmarkError(
            "http_error" if trace.get("http_status", 0) else "connection_or_service_error",
            trace.get("exception", str(exc)),
            trace,
        )
    finally:
        visualize.requests.get = real_get


def _call_legacy(field):
    requested_pixels = int(round(field["fov_deg"] * 3600.0 / 0.25))
    size_pix = max(1, min(512, requested_pixels))
    pixscale = field["fov_deg"] * 3600.0 / size_pix
    params = {
        "ra": field["request_ra_deg"],
        "dec": field["request_dec_deg"],
        "layer": "ls-dr11",
        "pixscale": pixscale,
        "size": size_pix,
        "bands": "g",
    }
    trace = {
        "layer": "ls-dr11",
        "bands": "g",
        "format": "fits",
        "request_params": params,
        "url": "https://www.legacysurvey.org/viewer/fits-cutout?{}".format(
            urlencode(params)
        ),
    }
    original_get = visualize.requests.get

    def traced_get(url, params=None, **kwargs):
        response = original_get(url, params=params, **kwargs)
        trace["http_status"] = getattr(response, "status_code", None)
        trace["download_bytes"] = len(getattr(response, "content", b""))
        trace["response_url"] = getattr(response, "url", "")
        return response

    visualize.requests.get = traced_get
    try:
        result = visualize.retrieve_image_legacy(
            field["request_ra_deg"], field["request_dec_deg"],
            field["fov_deg"], yflip=True, scale=0.25, return_wcs=True
        )
    except BaseException as exc:
        if not isinstance(exc, KeyboardInterrupt):
            exc.benchmark_trace = trace
        raise
    finally:
        visualize.requests.get = original_get
    trace["http_status"] = trace.get("http_status", 200)
    trace["download_bytes"] = trace.get("download_bytes")
    return result, trace


def _call_dss(field):
    params = {
        "ra": field["request_ra_deg"],
        "dec": field["request_dec_deg"],
        "x": field["fov_deg"] * 60.0,
        "y": field["fov_deg"] * 60.0,
        "mime-type": "download-fits",
    }
    trace = {
        "format": "fits",
        "request_params": params,
        "url": "http://archive.eso.org/dss/dss/image?{}".format(
            urlencode(params)
        ),
    }
    original_open = visualize.fits.open
    original_urlopen = urlopen

    def traced_open(name, *args, **kwargs):
        if isinstance(name, str) and name.startswith("http"):
            trace["url"] = name
            with original_urlopen(name, timeout=20) as response:
                payload = response.read()
                trace["http_status"] = response.getcode()
            trace["download_bytes"] = len(payload)
            return original_open(BytesIO(payload), *args, **kwargs)
        return original_open(name, *args, **kwargs)

    visualize.fits.open = traced_open
    try:
        result = visualize.retrieve_image_DSS(
            field["request_ra_deg"], field["request_dec_deg"],
            field["fov_deg"], yflip=True, return_wcs=True
        )
    except BaseException as exc:
        if not isinstance(exc, KeyboardInterrupt):
            exc.benchmark_trace = trace
        raise
    finally:
        visualize.fits.open = original_open
    return result, trace


def _call_mast(field):
    trace = {
        "survey": "3pi",
        "filter": "g",
        "stack_file": "data",
        "request_position": [field["request_ra_deg"], field["request_dec_deg"]],
        "fov_deg": field["fov_deg"],
        "native_size_pixels": mast_cutouts._native_cutout_size(field["fov_deg"]),
    }
    original_fetch = mast_cutouts.fetch_mast_cutout
    original_phase = mast_cutouts._job_phase
    original_job_url = mast_cutouts._job_url_from_response

    def traced_job_url(response):
        job_url = original_job_url(response)
        trace["job_url"] = job_url
        trace["submission_http_status"] = response.status_code
        return job_url

    def traced_phase(session, job_url):
        phase = original_phase(session, job_url)
        trace["last_phase"] = phase
        trace["job_url"] = job_url
        return phase

    def traced_fetch(*args, **kwargs):
        content, record, timing = original_fetch(*args, **kwargs)
        trace.update({
            "job_url": timing.get("job_url"),
            "science_url": timing.get("science_url"),
            "job_result": record,
            "submit_seconds": timing.get("submit_seconds"),
            "poll_seconds": timing.get("poll_seconds"),
            "download_seconds": timing.get("download_seconds"),
            "download_bytes": timing.get("bytes"),
            "native_size_pixels": timing.get("native_size_pixels"),
            "raw_content": content,
        })
        return content, record, timing

    mast_cutouts.fetch_mast_cutout = traced_fetch
    mast_cutouts._job_phase = traced_phase
    mast_cutouts._job_url_from_response = traced_job_url
    try:
        result = visualize.retrieve_image_mast(
            field["request_ra_deg"], field["request_dec_deg"],
            field["fov_deg"], yflip=True, scale=0.25, return_wcs=True
        )
    except BaseException as exc:
        if not isinstance(exc, KeyboardInterrupt):
            exc.benchmark_trace = trace
        raise
    finally:
        mast_cutouts.fetch_mast_cutout = original_fetch
        mast_cutouts._job_phase = original_phase
        mast_cutouts._job_url_from_response = original_job_url
    return result, trace


def _call_provider(provider, field):
    return {
        "hips2fits": _call_hips_primary,
        "legacy": _call_legacy,
        "mast": _call_mast,
        "dss": _call_dss,
    }[provider](field)


def _write_product(output_dir, field, provider, repetition, result):
    data, _, _, _, image_wcs = result
    if image_wcs is None:
        return None
    path = output_dir / "products" / "{}_{}_{}.fits".format(
        field["id"], provider, repetition
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([
        fits.PrimaryHDU(data=numpy.asarray(data), header=image_wcs.to_header())
    ]).writeto(path, overwrite=True)
    return str(path)


def _run_attempt(field, provider, repetition, output_dir, timeout_seconds):
    started_wall = _utc_now()
    started = time.monotonic()
    trace = {}
    result = None
    error = None
    validation_started = None
    validation_elapsed = None
    product_path = None
    try:
        with _Deadline(timeout_seconds):
            result, trace = _call_provider(provider, field)
            validation_started = time.monotonic()
            _validate_result(result, field, trace)
            validation_elapsed = time.monotonic() - validation_started
            product_path = _write_product(
                output_dir, field, provider, repetition, result
            )
    except BaseException as exc:
        if isinstance(exc, KeyboardInterrupt):
            raise
        error = exc
        if hasattr(exc, "benchmark_trace"):
            trace.update(exc.benchmark_trace)
        if isinstance(exc, BenchmarkError):
            trace.update(exc.trace)
    elapsed = time.monotonic() - started
    raw_content = trace.pop("raw_content", None)
    if raw_content is not None:
        raw_path = output_dir / "products" / "{}_{}_{}_raw.fits".format(
            field["id"], provider, repetition
        )
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_bytes(raw_content)
        trace["raw_product"] = str(raw_path)

    success = error is None
    category = "success" if success else _classify_exception(error, trace)
    record = {
        "started_utc": started_wall,
        "provider": provider,
        "field_id": field["id"],
        "repetition": repetition,
        "request_ra_deg": field["request_ra_deg"],
        "request_dec_deg": field["request_dec_deg"],
        "fov_deg": field["fov_deg"],
        "pa_deg": field["pa_deg"],
        "success": success,
        "failure_category": category,
        "error": "" if error is None else str(error),
        "elapsed_seconds": elapsed,
        "validation_seconds": validation_elapsed,
        "product_path": product_path,
        "http_status": trace.get("http_status"),
        "download_bytes": trace.get("download_bytes"),
        "shape": json.dumps(trace.get("shape")) if trace.get("shape") else "",
        "source": trace.get("source", ""),
        "url": trace.get("url", trace.get("science_url", "")),
        "job_url": trace.get("job_url", ""),
        "last_phase": trace.get("last_phase", ""),
        "submit_seconds": trace.get("submit_seconds"),
        "poll_seconds": trace.get("poll_seconds"),
        "download_seconds": trace.get("download_seconds"),
        "coverage": json.dumps(trace.get("coverage")) if trace.get("coverage") else "",
        "wcs": json.dumps(trace.get("wcs")) if trace.get("wcs") else "",
    }
    return record


def _csv_fields():
    return [
        "started_utc", "provider", "field_id", "repetition",
        "request_ra_deg", "request_dec_deg", "fov_deg", "pa_deg",
        "success", "failure_category", "error", "elapsed_seconds",
        "validation_seconds", "product_path", "http_status",
        "download_bytes", "shape", "source", "url", "job_url",
        "last_phase", "submit_seconds", "poll_seconds", "download_seconds",
        "coverage", "wcs",
    ]


def _append_csv(path, record):
    fields = _csv_fields()
    write_header = not path.exists()
    with path.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if write_header:
            writer.writeheader()
        writer.writerow({key: record.get(key, "") for key in fields})


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _summarize(records):
    summary = {}
    for provider in PROVIDERS:
        rows = [r for r in records if r["provider"] == provider]
        successes = [r for r in rows if r["success"]]
        success_times = [float(r["elapsed_seconds"]) for r in successes]
        failure_rows = [r for r in rows if not r["success"]]
        failure_times = [float(r["elapsed_seconds"]) for r in failure_rows]
        failures = Counter(r["failure_category"] for r in rows if not r["success"])
        summary[provider] = {
            "attempts": len(rows),
            "successful_eligible": len(successes),
            "success_rate": len(successes) / float(len(rows)) if rows else 0.0,
            "minimum_seconds": min(success_times) if success_times else None,
            "median_seconds": float(numpy.median(success_times)) if success_times else None,
            "mean_seconds": float(numpy.mean(success_times)) if success_times else None,
            "maximum_seconds": max(success_times) if success_times else None,
            "failed_minimum_seconds": min(failure_times) if failure_times else None,
            "failed_median_seconds": float(numpy.median(failure_times)) if failure_times else None,
            "failed_mean_seconds": float(numpy.mean(failure_times)) if failure_times else None,
            "failed_maximum_seconds": max(failure_times) if failure_times else None,
            "failure_count": len(rows) - len(successes),
            "failure_categories": dict(failures),
            "covered_fields": sorted({r["field_id"] for r in successes}),
        }
    return summary


def _fallback_simulation(records):
    by_key = {(r["field_id"], r["repetition"], r["provider"]): r for r in records}
    result = {}
    fields = sorted({r["field_id"] for r in records})
    repetitions = sorted({r["repetition"] for r in records})
    for name, sequence in SEQUENCES.items():
        simulations = []
        for field_id in fields:
            for repetition in repetitions:
                total = 0.0
                final_provider = None
                success = False
                for provider in sequence:
                    row = by_key.get((field_id, repetition, provider))
                    if row is None:
                        continue
                    total += float(row["elapsed_seconds"])
                    if row["success"]:
                        success = True
                        final_provider = provider
                        break
                simulations.append({
                    "field_id": field_id,
                    "repetition": repetition,
                    "success": success,
                    "elapsed_seconds": total,
                    "final_provider": final_provider,
                })
        successful = [x for x in simulations if x["success"]]
        finals = Counter(x["final_provider"] for x in successful)
        times = [x["elapsed_seconds"] for x in successful]
        result[name] = {
            "sequence": list(sequence),
            "simulations": len(simulations),
            "successful": len(successful),
            "success_rate": len(successful) / float(len(simulations)) if simulations else 0.0,
            "exhausted": len(simulations) - len(successful),
            "median_time_to_first_usable_seconds": float(numpy.median(times)) if times else None,
            "mean_time_to_first_usable_seconds": float(numpy.mean(times)) if times else None,
            "maximum_time_to_first_usable_seconds": max(times) if times else None,
            "final_provider_counts": dict(finals),
            "simulations_detail": simulations,
        }
    return result


def _plot_results(records, output_dir):
    providers = list(PROVIDERS)
    labels = ["CDS HiPS", "Legacy", "MAST", "DSS"]
    fig, ax = plt.subplots(figsize=(10, 6))
    for index, provider in enumerate(providers):
        rows = [r for r in records if r["provider"] == provider]
        for row in rows:
            color = "tab:blue" if row["success"] else "tab:red"
            ax.scatter(index, row["elapsed_seconds"], color=color, alpha=0.75)
    ax.set_xticks(range(len(providers)), labels)
    ax.set_ylabel("Attempt elapsed time (s)")
    ax.set_title("Image retrieval timing; blue=eligible, red=failed")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "retrieval_times.png", dpi=150)
    plt.close(fig)

    counts = []
    for provider in providers:
        rows = [r for r in records if r["provider"] == provider]
        counts.append([
            sum(bool(r["success"]) for r in rows),
            sum(not bool(r["success"]) for r in rows),
        ])
    counts = numpy.asarray(counts)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels, counts[:, 0], label="eligible", color="tab:blue")
    ax.bar(labels, counts[:, 1], bottom=counts[:, 0], label="failed", color="tab:red")
    ax.set_ylabel("Attempts")
    ax.set_title("Retrieval outcomes")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "retrieval_outcomes.png", dpi=150)
    plt.close(fig)


def _markdown_report(output_dir, fields, records, summary, simulations):
    lines = [
        "# hetshuffle image retrieval benchmark",
        "",
        "Development-only sequential benchmark; production dispatch and timeout defaults were not changed.",
        "",
        "## Design",
        "",
        "- Fields: 5; repetitions per provider/field: 3; planned attempts: 60.",
        "- Attempt order: interleaved by repetition and field; no concurrency.",
        "- Per-attempt end-to-end budget: 30 seconds by default.",
        "- Retrieval timing excludes plotting and product serialization; validation is recorded separately.",
        "",
        "## Fields",
        "",
        "| Field | Retrieval center (RA, Dec) | FOV (deg) | PA (deg) | Selection |",
        "|---|---:|---:|---:|---|",
    ]
    for field in fields:
        lines.append(
            "| {} | {:.8f}, {:.8f} | {:.8f} | {:.6f} | {} |".format(
                field["id"], field["request_ra_deg"], field["request_dec_deg"],
                field["fov_deg"], field["pa_deg"], field["selection"],
            )
        )
    lines.extend([
        "",
        "## Provider results",
        "",
        "| Provider | Successful / eligible | Success rate | Median time (s) | Mean time (s) | Maximum time (s) | Failures |",
        "|---|---:|---:|---:|---:|---:|---|",
    ])
    names = {"hips2fits": "CDS hips2fits", "legacy": "Legacy Survey", "mast": "MAST PS1", "dss": "ESO DSS"}
    for provider in PROVIDERS:
        item = summary[provider]
        lines.append(
            "| {} | {}/{} | {:.1%} | {} | {} | {} | {} |".format(
                names[provider], item["successful_eligible"], item["attempts"],
                item["success_rate"],
                "-" if item["median_seconds"] is None else "{:.3f}".format(item["median_seconds"]),
                "-" if item["mean_seconds"] is None else "{:.3f}".format(item["mean_seconds"]),
                "-" if item["maximum_seconds"] is None else "{:.3f}".format(item["maximum_seconds"]),
                "; ".join("{}={}".format(k, v) for k, v in item["failure_categories"].items()) or "none",
            )
        )
    lines.extend([
        "",
        "Failed-attempt timing (separate from successful retrieval statistics):",
        "",
        "| Provider | Minimum failed time (s) | Median failed time (s) | Mean failed time (s) | Maximum failed time (s) |",
        "|---|---:|---:|---:|---:|",
    ])
    for provider in PROVIDERS:
        item = summary[provider]
        lines.append(
            "| {} | {} | {} | {} | {} |".format(
                names[provider],
                "-" if item["failed_minimum_seconds"] is None else "{:.3f}".format(item["failed_minimum_seconds"]),
                "-" if item["failed_median_seconds"] is None else "{:.3f}".format(item["failed_median_seconds"]),
                "-" if item["failed_mean_seconds"] is None else "{:.3f}".format(item["failed_mean_seconds"]),
                "-" if item["failed_maximum_seconds"] is None else "{:.3f}".format(item["failed_maximum_seconds"]),
            )
        )
    lines.extend(["", "## Exploratory fallback simulations", ""])
    lines.append("These simulations pair the same field/repetition observations across providers and include failed-provider elapsed time. They are not reliability forecasts.")
    lines.extend([
        "",
        "| Sequence | Success rate | Median first-usable time (s) | Mean (s) | Exhausted | Final provider counts |",
        "|---|---:|---:|---:|---:|---|",
    ])
    for name, item in simulations.items():
        lines.append(
            "| {} ({}) | {:.1%} | {} | {} | {} | {} |".format(
                name, " → ".join(item["sequence"]), item["success_rate"],
                "-" if item["median_time_to_first_usable_seconds"] is None else "{:.3f}".format(item["median_time_to_first_usable_seconds"]),
                "-" if item["mean_time_to_first_usable_seconds"] is None else "{:.3f}".format(item["mean_time_to_first_usable_seconds"]),
                item["exhausted"], item["final_provider_counts"],
            )
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        "This sample is suitable for comparing the observed run and failure behavior only. It is too small to establish long-term uptime or tail-latency guarantees. Coverage failures are classified separately from operational failures where the returned WCS shows that the ACAM footprint is outside the image.",
        "",
        "Artifacts: `retrieval_attempts.csv`, `retrieval_times.png`, `retrieval_outcomes.png`, and `benchmark_results.json`.",
        "",
    ])
    (output_dir / "benchmark_report.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("/private/tmp/hetshuffle-image-retrieval-benchmark"),
    )
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    args = parser.parse_args()
    if args.repetitions < 1 or args.timeout_seconds <= 0:
        parser.error("repetitions and timeout-seconds must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config = _config_for_experiment(args.output_dir / "config")
    fields = _prepare_fields(config)
    csv_path = args.output_dir / "retrieval_attempts.csv"
    records = []
    schedule = []
    for repetition in range(1, args.repetitions + 1):
        for field_index, field in enumerate(fields):
            start = (repetition + field_index) % len(PROVIDERS)
            for offset in range(len(PROVIDERS)):
                schedule.append((field, PROVIDERS[(start + offset) % len(PROVIDERS)], repetition))

    (args.output_dir / "field_definitions.json").write_text(json.dumps(fields, indent=2))
    (args.output_dir / "schedule.json").write_text(json.dumps([
        {"field_id": f["id"], "provider": p, "repetition": r}
        for f, p, r in schedule
    ], indent=2))

    for index, (field, provider, repetition) in enumerate(schedule, 1):
        print("[{}/{}] {} {} repetition {}".format(
            index, len(schedule), provider, field["id"], repetition
        ), flush=True)
        record = _run_attempt(
            field, provider, repetition, args.output_dir, args.timeout_seconds
        )
        records.append(record)
        _append_csv(csv_path, record)
        print("  {} {:.3f}s {}".format(
            "SUCCESS" if record["success"] else "FAIL",
            record["elapsed_seconds"], record["failure_category"],
        ), flush=True)

    summary = _summarize(records)
    simulations = _fallback_simulation(records)
    _plot_results(records, args.output_dir)
    results = {
        "planned_attempts": len(schedule),
        "completed_attempts": len(records),
        "timeout_seconds": args.timeout_seconds,
        "fields": fields,
        "summary": summary,
        "fallback_simulations": simulations,
    }
    (args.output_dir / "benchmark_results.json").write_text(json.dumps(results, indent=2))
    _markdown_report(args.output_dir, fields, records, summary, simulations)
    print(json.dumps({
        "completed_attempts": len(records),
        "planned_attempts": len(schedule),
        "summary": summary,
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
