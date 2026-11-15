from __future__ import annotations

import configparser
import os
from pathlib import Path

import pytest


os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture
def shuffle_config(tmp_path):
    """Return the packaged configuration with test output directories."""
    package_config = Path(__file__).parents[1] / "hetdex_shuffle" / "configs" / "shuffle.cfg"
    config = configparser.ConfigParser(defaults={"xpa_method": "local"})
    assert config.read(package_config)

    output_dirs = {
        "images": tmp_path / "images",
        "acam_images": tmp_path / "acam_images",
        "ds9_regions": tmp_path / "regions",
    }
    for option, path in output_dirs.items():
        config.set("directories", option, str(path))
    return config
