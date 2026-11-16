from __future__ import annotations

import importlib.metadata
from pathlib import Path
from types import SimpleNamespace

import pytest


def test_package_metadata_and_console_entry_points():
    import hetdex_shuffle
    from hetdex_shuffle import copy_config, do_shuffle_target, obj_dist_ifu_vis

    assert hetdex_shuffle.__version__ == importlib.metadata.version("hetdex-shuffle")
    assert callable(copy_config.main)
    assert callable(do_shuffle_target.main)
    assert callable(obj_dist_ifu_vis.main)


def test_default_config_and_argument_parsing():
    from hetdex_shuffle.do_shuffle_target import load_config, parse_args

    config = load_config("/path/that/does/not/exist/shuffle.cfg")
    assert config.getboolean("General", "visualize_ACAM")
    assert config.get("General", "catalog") == "PANSTARRS"

    config_path = Path(__file__).parents[1] / "hetdex_shuffle" / "configs" / "shuffle.cfg"
    args = parse_args([
        "-c",
        str(config_path),
        "12:34:56",
        "--",
        "-30:15:00",
        "0.01",
        "0",
        "58",
    ])
    assert args.ra == pytest.approx(12 + 34 / 60 + 56 / 3600)
    assert args.dec == pytest.approx(-30.25)
    assert args.radius == pytest.approx(0.01)
    assert args.ifuslot == "058"


def test_packaged_configuration_can_be_copied(tmp_path):
    from hetdex_shuffle.copy_config import COPY_FILES, copy

    copy(SimpleNamespace(to_dir=str(tmp_path), force=False))
    for filename in COPY_FILES:
        copied = tmp_path / filename
        assert copied.is_file()
        assert copied.stat().st_size > 0


def test_missing_config_warns_and_uses_packaged_default(tmp_path):
    from hetdex_shuffle.do_shuffle_target import load_config

    with pytest.warns(UserWarning):
        config = load_config(str(tmp_path / "missing.cfg"))
    assert config.has_section("General")
