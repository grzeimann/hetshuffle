import configparser

import numpy
import pytest

from hetdex_shuffle import visualize


def _config(image_service=None):
    config = configparser.ConfigParser()
    config.add_section("General")
    if image_service is not None:
        config.set("General", "image_service", image_service)
    return config


def _image(source):
    return numpy.ones((2, 2)), numpy.eye(2), "https://example.invalid/image", source


@pytest.mark.parametrize("service, function_name", [
    ("CDS", "retrieve_image_hips"),
    ("Legacy", "retrieve_image_legacy"),
    ("MAST", "retrieve_image_mast"),
    ("DSS", "retrieve_image_DSS"),
])
def test_explicit_image_service_selects_only_requested_provider(
        monkeypatch, service, function_name):
    calls = []

    def provider(*args, **kwargs):
        calls.append(function_name)
        return _image(service)

    monkeypatch.setattr(visualize, function_name, provider)
    result = visualize.retrieve_image(
        180.0, 30.0, 0.01, _config(service), True, scale=1.5
    )

    assert result[3] == service
    assert calls == [function_name]


def test_auto_tries_mast_before_dss_after_cds_and_legacy_fail(monkeypatch):
    calls = []

    def empty_cds(*args, **kwargs):
        calls.append("CDS")
        return visualize._retrieved_image(
            *visualize.return_blank_image(2, 1.5),
            return_wcs=kwargs.get("return_wcs", False),
        )

    def failed_legacy(*args, **kwargs):
        calls.append("Legacy")
        raise RuntimeError("legacy unavailable")

    def successful_mast(*args, **kwargs):
        calls.append("MAST")
        return _image("MAST/PS1/3pi/g")

    def unexpected_dss(*args, **kwargs):
        calls.append("DSS")
        raise AssertionError("DSS should not be reached")

    monkeypatch.setattr(visualize, "retrieve_image_hips", empty_cds)
    monkeypatch.setattr(visualize, "retrieve_image_legacy", failed_legacy)
    monkeypatch.setattr(visualize, "retrieve_image_mast", successful_mast)
    monkeypatch.setattr(visualize, "retrieve_image_DSS", unexpected_dss)

    result = visualize.retrieve_image(
        180.0, 30.0, 0.01, _config(), True, scale=1.5
    )

    assert result[3] == "MAST/PS1/3pi/g"
    assert calls == ["CDS", "Legacy", "MAST"]


def test_explicit_provider_failure_returns_blank_without_switching(monkeypatch):
    calls = []

    def failed_mast(*args, **kwargs):
        calls.append("MAST")
        raise RuntimeError("MAST unavailable")

    def unexpected_dss(*args, **kwargs):
        calls.append("DSS")
        raise AssertionError("explicit MAST must not switch to DSS")

    monkeypatch.setattr(visualize, "retrieve_image_mast", failed_mast)
    monkeypatch.setattr(visualize, "retrieve_image_DSS", unexpected_dss)

    result = visualize.retrieve_image(
        180.0, 30.0, 0.01, _config("mast"), True, scale=1.5
    )

    assert result[3] == ""
    assert calls == ["MAST"]


def test_unknown_image_service_is_rejected():
    with pytest.raises(ValueError, match="Unknown image service"):
        visualize.retrieve_image(
            180.0, 30.0, 0.01, _config("unknown"), True, scale=1.5
        )
