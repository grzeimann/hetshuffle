"""Experimental client for the asynchronous MAST cutouts service.

This module is intentionally not connected to the production image-selection
path.  It provides a native PS1 FITS/WCS retrieval path for experiments and
tests.
"""

from __future__ import absolute_import, division, print_function

import math
import logging
import time
import warnings
from io import BytesIO
from urllib.parse import urljoin

import numpy
import requests
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import WCS
from astropy.wcs.wcs import FITSFixedWarning


MAST_CUTOUTS_BASE_URL = "https://cutouts.mast.stsci.edu"
MAST_ASYNC_URL = MAST_CUTOUTS_BASE_URL + "/api/v0/cutouts/async"
MAST_NATIVE_PIXEL_SCALE_ARCSEC = 0.25
MAST_MAX_SIZE_PIXELS = 6000
MAST_REQUEST_TIMEOUT_SEC = 60
MAST_POLL_INTERVAL_SEC = 0.5
MAST_POLL_TIMEOUT_SEC = 30 * 60
MAST_RESULTS_LIMIT = 100


class MastCutoutError(RuntimeError):
    """Raised when a MAST cutout cannot be retrieved or validated."""


def assess_footprint_coverage(image_wcs, image_shape, sky_points):
    """Assess whether sky points fall inside an image WCS footprint.

    ``sky_points`` is an array-like collection of ``(ra, dec)`` pairs in
    degrees. This is deliberately a geometric check only; it does not claim
    that every in-bounds pixel contains useful survey signal.
    """
    sky_points = numpy.asarray(sky_points, dtype=float)
    if sky_points.ndim != 2 or sky_points.shape[1] != 2:
        raise ValueError("sky_points must have shape (N, 2)")
    pixels = numpy.asarray(image_wcs.wcs_world2pix(sky_points, 0), dtype=float)
    height, width = image_shape[:2]
    x_min, x_max = float(pixels[:, 0].min()), float(pixels[:, 0].max())
    y_min, y_max = float(pixels[:, 1].min()), float(pixels[:, 1].max())
    margins = [x_min, (width - 1.0) - x_max,
               y_min, (height - 1.0) - y_max]
    return {
        "pixel_points": pixels.tolist(),
        "pixel_bbox": [[x_min, x_max], [y_min, y_max]],
        "image_shape": [int(height), int(width)],
        "all_points_inside": bool(numpy.all(pixels[:, 0] >= 0) and
                                   numpy.all(pixels[:, 0] <= width - 1) and
                                   numpy.all(pixels[:, 1] >= 0) and
                                   numpy.all(pixels[:, 1] <= height - 1)),
        "bbox_inside": bool(min(margins) >= 0),
        "minimum_edge_margin_pixels": float(min(margins)),
    }


def _native_cutout_size(fov_deg):
    """Return a native PS1 square size that covers ``fov_deg``."""
    if not math.isfinite(float(fov_deg)) or float(fov_deg) <= 0:
        raise ValueError("fov_deg must be a positive finite value")
    size = int(math.ceil(float(fov_deg) * 3600.0 /
                        MAST_NATIVE_PIXEL_SCALE_ARCSEC))
    if size < 1:
        size = 1
    if size > MAST_MAX_SIZE_PIXELS:
        raise MastCutoutError(
            "requested native PS1 cutout is {} pixels, exceeding MAST's {} "
            "pixel limit".format(size, MAST_MAX_SIZE_PIXELS)
        )
    return size


def _multipart_fields(ra, dec, size_pixels):
    """Build current documented multipart fields for one PS1 product."""
    # requests sends a multipart/form-data body when fields are supplied via
    # ``files``.  The (None, value) form prevents requests from treating the
    # values as file uploads while preserving repeated field names.
    values = [
        ("position", "{:.12g},{:.12g}".format(float(ra), float(dec))),
        ("size", str(int(size_pixels))),
        ("output_format", "fits"),
        ("ps1.filter", "g"),
        ("ps1.survey", "3pi"),
        ("ps1.stack_file", "data"),
    ]
    return [(name, (None, value)) for name, value in values]


def _job_url_from_response(response):
    if response.status_code not in (200, 303):
        raise MastCutoutError(
            "MAST submission failed (HTTP {}): {}".format(
                response.status_code, response.text[:500]
            )
        )
    location = response.headers.get("Location", "")
    if not location:
        raise MastCutoutError("MAST submission returned no Location header")
    return urljoin(MAST_ASYNC_URL + "/", location)


def _job_phase(session, job_url):
    response = session.get(
        job_url.rstrip("/") + "/phase",
        timeout=MAST_REQUEST_TIMEOUT_SEC,
    )
    response.raise_for_status()
    return response.text.strip().strip('"').upper()


def _job_summary(session, job_url):
    response = session.get(
        job_url.rstrip("/") + "/results/summary",
        timeout=MAST_REQUEST_TIMEOUT_SEC,
    )
    response.raise_for_status()
    return response.json()


def _job_error(session, job_url):
    try:
        response = session.get(
            job_url.rstrip("/") + "/error",
            timeout=MAST_REQUEST_TIMEOUT_SEC,
        )
        detail = response.text.strip()
        if detail:
            return detail[:1000]
    except Exception as exc:  # pragma: no cover - diagnostic fallback
        return "unable to read MAST error endpoint: {}".format(exc)
    return "MAST returned no error details"


def _science_download_url(science_uri, job_url):
    """Resolve a result's HTTP or public S3 science location."""
    if science_uri.startswith("s3://"):
        bucket_and_key = science_uri[len("s3://"):]
        bucket, separator, key = bucket_and_key.partition("/")
        if not separator or not bucket or not key:
            raise MastCutoutError(
                "invalid S3 science URI returned by MAST: {}".format(science_uri)
            )
        return "https://{}.s3.amazonaws.com/{}".format(bucket, key)
    return urljoin(job_url.rstrip("/") + "/", science_uri)


def _wait_for_job(session, job_url):
    terminal = {"COMPLETED", "ERROR", "ABORTED"}
    started = time.monotonic()
    last_summary = {}

    while True:
        phase = _job_phase(session, job_url)
        try:
            last_summary = _job_summary(session, job_url)
        except Exception:
            # The phase endpoint is authoritative for lifecycle state. A
            # summary failure should not hide the eventual phase or error.
            last_summary = {}

        total_jobs = last_summary.get("total_jobs", 0)
        finished_jobs = sum(
            last_summary.get(key, 0)
            for key in ("completed_jobs", "failed_jobs", "skipped_jobs")
        )
        if total_jobs and finished_jobs >= total_jobs:
            if last_summary.get("failed_jobs", 0) or last_summary.get("skipped_jobs", 0):
                raise MastCutoutError(
                    "MAST job completed with failed or skipped cutouts: {}; {}".format(
                        last_summary, _job_error(session, job_url)
                    )
                )
            return "COMPLETED", last_summary, time.monotonic() - started

        if phase in terminal:
            if phase != "COMPLETED":
                raise MastCutoutError(
                    "MAST job {} ended in {}: {}".format(
                        job_url.rstrip("/").split("/")[-1],
                        phase,
                        _job_error(session, job_url),
                    )
                )
            return phase, last_summary, time.monotonic() - started

        if time.monotonic() - started >= MAST_POLL_TIMEOUT_SEC:
            raise TimeoutError(
                "MAST job {} did not complete within {} seconds".format(
                    job_url.rstrip("/").split("/")[-1],
                    MAST_POLL_TIMEOUT_SEC,
                )
            )
        time.sleep(MAST_POLL_INTERVAL_SEC)


def _result_records(session, job_url):
    """Fetch all pages of result metadata for one asynchronous job."""
    records = []
    page_url = job_url.rstrip("/") + "/results/cutouts"
    params = {
        "output_format": "json",
        "page": 0,
        "limit": MAST_RESULTS_LIMIT,
    }

    while page_url:
        response = session.get(
            page_url,
            params=params,
            timeout=MAST_REQUEST_TIMEOUT_SEC,
        )
        response.raise_for_status()
        payload = response.json()
        records.extend(payload.get("results", []))
        links = payload.get("links") or {}
        next_url = links.get("next") if isinstance(links, dict) else None
        page_url = urljoin(page_url, next_url) if next_url else None
        params = None
    return records


def fetch_mast_cutout(ra, dec, fov_deg, scale=1.5, session=None):
    """Submit and download one native PS1 g-band FITS cutout.

    Returns ``(fits_bytes, result_record, timing)``.  This lower-level helper
    is used by the development integration script so it can preserve the raw
    MAST FITS product exactly as downloaded.
    """
    if not math.isfinite(float(scale)) or float(scale) <= 0:
        raise ValueError("scale must be a positive finite value")
    size_pixels = _native_cutout_size(fov_deg)
    own_session = session is None
    session = session or requests.Session()
    started = time.monotonic()
    try:
        submit_started = time.monotonic()
        response = session.post(
            MAST_ASYNC_URL,
            files=_multipart_fields(ra, dec, size_pixels),
            allow_redirects=False,
            timeout=MAST_REQUEST_TIMEOUT_SEC,
        )
        job_url = _job_url_from_response(response)
        submit_seconds = time.monotonic() - submit_started

        _, summary, poll_seconds = _wait_for_job(session, job_url)
        result_records = _result_records(session, job_url)
        if len(result_records) != 1:
            raise MastCutoutError(
                "expected exactly one MAST result for one g/data request, "
                "received {} (summary={!r})".format(
                    len(result_records), summary
                )
            )
        result = result_records[0]
        extras = result.get("mission_extras") or {}
        if result.get("mission") != "ps1" or result.get("filter") != "g":
            raise MastCutoutError(
                "MAST returned an unexpected mission/filter result: {!r}".format(
                    result
                )
            )
        if extras.get("survey") not in (None, "3pi"):
            raise MastCutoutError(
                "MAST returned a non-3pi result: {!r}".format(result)
            )
        if extras.get("file_type") not in (None, "data", "stack.data"):
            raise MastCutoutError(
                "MAST returned a non-science stack product: {!r}".format(result)
            )
        science = result.get("science")
        if not science:
            raise MastCutoutError(
                "MAST result has no science FITS URL: {!r}".format(result)
            )
        science_url = _science_download_url(science, job_url)

        download_started = time.monotonic()
        fits_response = session.get(
            science_url,
            timeout=MAST_REQUEST_TIMEOUT_SEC,
        )
        fits_response.raise_for_status()
        download_seconds = time.monotonic() - download_started
        timing = {
            "submit_seconds": submit_seconds,
            "poll_seconds": poll_seconds,
            "download_seconds": download_seconds,
            "total_seconds": time.monotonic() - started,
            "bytes": len(fits_response.content),
            "native_size_pixels": size_pixels,
            "science_uri": science,
            "science_url": science_url,
            "job_url": job_url,
        }
        return fits_response.content, result, timing
    finally:
        if own_session:
            session.close()


def _decode_mast_fits(content, yflip):
    """Decode and validate a MAST FITS image, preserving its celestial WCS."""
    with fits.open(BytesIO(content), memmap=False) as hdul:
        selected = None
        for hdu in hdul:
            if hdu.data is not None:
                selected = hdu
                break
        if selected is None:
            raise MastCutoutError("MAST FITS product contains no image data")
        header = selected.header.copy()
        imarray = numpy.asarray(selected.data).copy()

    if imarray.ndim != 2:
        raise MastCutoutError(
            "MAST PS1 science product is not a 2-D image: shape {}".format(
                imarray.shape
            )
        )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            warnings.simplefilter("ignore", VerifyWarning)
            image_wcs = WCS(header, naxis=2)
        if not image_wcs.has_celestial:
            raise ValueError("WCS has no celestial axes")
    except Exception as exc:
        raise MastCutoutError("MAST FITS has invalid celestial WCS: {}".format(exc))

    if yflip:
        imarray = imarray[::-1, :]
        image_wcs.wcs.crpix[1] = imarray.shape[0] + 1.0 - image_wcs.wcs.crpix[1]
        if image_wcs.wcs.has_cd():
            image_wcs.wcs.cd[:, 1] *= -1.0
        else:
            image_wcs.wcs.pc[:, 1] *= -1.0

    cd = numpy.asarray(image_wcs.pixel_scale_matrix, dtype=float)
    if cd.shape != (2, 2) or not numpy.all(numpy.isfinite(cd)):
        raise MastCutoutError("MAST FITS WCS has no finite 2-D pixel scale matrix")
    return imarray, numpy.matrix(cd), image_wcs


def retrieve_image_mast(ra, dec, fov_deg, yflip=True, scale=1.5,
                        return_wcs=False):
    """Experimentally retrieve a native PS1 g-band FITS image from MAST.

    The MAST service's ``size`` parameter is native PS1 pixels, so this
    backend requests enough 0.25 arcsec pixels to cover ``fov_deg``. It does
    not resample the image to the hips2fits pixel dimensions. ``scale`` is
    retained for interface compatibility and validation only.
    """
    started = time.monotonic()
    content, result, timing = fetch_mast_cutout(
        ra, dec, fov_deg, scale=scale
    )
    decode_started = time.monotonic()
    imarray, cd, image_wcs = _decode_mast_fits(content, yflip=yflip)
    decode_seconds = time.monotonic() - decode_started
    total_seconds = time.monotonic() - started

    source = "MAST/PS1/3pi/g"
    science_url = timing["science_url"]
    logging_message = (
        "MAST PS1 cutout source=%s shape=%s dtype=%s bytes=%d "
        "submit=%.3fs poll=%.3fs download=%.3fs decode=%.3fs total=%.3fs"
        % (
            source,
            imarray.shape,
            imarray.dtype,
            timing["bytes"],
            timing["submit_seconds"],
            timing["poll_seconds"],
            timing["download_seconds"],
            decode_seconds,
            total_seconds,
        )
    )
    logging.getLogger("shuffle").info(logging_message)

    result_tuple = (imarray, cd, science_url, source)
    if return_wcs:
        return result_tuple + (image_wcs,)
    return result_tuple
