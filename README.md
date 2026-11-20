# HET Shuffle (hetdex-shuffle)

[![CI](https://github.com/grzeimann/hetshuffle/actions/workflows/python-tests.yml/badge.svg)](https://github.com/grzeimann/hetshuffle/actions/workflows/python-tests.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

HET Shuffle selects guide and wavefront-sensor stars and produces verified
target placements for HET observations. It is used on the mountain by
astronomers and telescope operators setting up HET tracks.

For contributor procedures, see the [Git workflow](git-workflow.md). For the
testing, CI, and scientific validation policy, see the
[testing and CI plan](testing-and-ci-plan.md).

## Install

hetshuffle supports Python 3.12 or newer. The standard mountain workflow
includes `reproject` and PostgreSQL catalog access through
`psycopg2-binary`.
`pyds9` is optional and is not required for ordinary operation or testing.

### Conda environment

```bash
conda env create --file environment.yml
conda activate hetshuffle
python -m pip install .
```

### Development installation

Install the development extra in the active environment to add pytest,
coverage, Ruff, and pre-commit:

```bash
python -m pip install --editable '.[dev]'
python -m pre_commit install --hook-type pre-commit
python -m pre_commit install --hook-type pre-push
```

Optional integrations can be installed separately:

```bash
python -m pip install '.[ds9]'     # DS9/XPA support; may require X11
python -m pip install '.[oracle]'  # Oracle catalog support
```

The package runtime dependencies are also listed in
[requirements.txt](requirements.txt). The Conda environment is the preferred
way to reproduce the supported mountain dependency boundary.

## Quickstart

Copy the packaged configuration files into a clean working directory:

```bash
mkdir -p hetshuffle-run
cd hetshuffle-run
shuffle_config
```

Verify the command-line installation:

```bash
do_shuffle --help
```

Example network-backed shuffle using the standard PANSTARRS catalog:

```bash
do_shuffle 21:22:09.26 '+49:26:24.20' 0 0 603 0 0 J212209 \
  --use_brightness true \
  --visualize_probestars false \
  --catalog PANSTARRS \
  -v
```

The example requires access to the configured catalog and image services.
Use `--config` and `--fplane_file` when running with site-specific files.
Right Ascension is supplied in hours; declination accepts a signed
degree-based sexagesimal value. Quote positive declinations containing `+`.

## Command-line usage

```text
do_shuffle RA DEC RADIUS TRACK IFUSLOT [X_OFFSET] [Y_OFFSET] [TARGET] [options]
```

- `RA`: target right ascension in hours, such as `08:00:16.08`.
- `DEC`: target declination in degrees, such as `+40:29:56.00`.
- `RADIUS`: search radius in degrees for guide and WFS stars.
- `TRACK`: `0` for East or `1` for West.
- `IFUSLOT`: three-digit target IFU slot; `000` is the boresight.
- `X_OFFSET`, `Y_OFFSET`: optional ACAM-plane offsets in arcseconds.
- `TARGET`: optional observation name.

Common options include:

- `--config PATH`: use a specific `shuffle.cfg`.
- `--fplane_file PATH`: use a specific focal-plane file.
- `--catalog NAME`: select a catalog such as `PANSTARRS`, `GAIA`, or
  `PANSTARRS_LOCAL`.
- `--image_service NAME`: override the configured image service
  (`auto`, `CDS`, `Legacy`, `MAST`, or `DSS`).
- `--visualize PATH`: write a sky-view visualization.
- `--visualize_ACAM PATH`: write an ACAM-frame visualization.
- `--visualize_probestars VALUE`: control probe-star visualization.
- `--epoch`, `--pmRA`, `--pmDEC`: provide target proper-motion information.
- `-v` or `--verbose`: increase logging verbosity; it can be repeated.

Image retrieval is controlled by `image_service` in the `[General]` section
of `shuffle.cfg`. Supported values are `auto`, `CDS`, `Legacy`, `MAST`,
and `DSS`. The default `auto` order is CDS HiPS, combined Legacy Survey
DR11, MAST Pan-STARRS, then ESO DSS. If all services fail, a blank fallback
image is produced.

## Configuration and data

`shuffle_config` copies the packaged `shuffle.cfg` and `fplane.txt` into the
current directory. The configuration controls catalog selection, patrol
geometry, magnitude limits, visualization, output names, and output
directories.

The standard configuration uses:

- `PANSTARRS` for network-based star selection.
- `fplane.txt` for the focal-plane geometry.
- ACAM and sky visualizations enabled by default.
- `PANSTARRS_LOCAL` support through the `[database]` settings when the
  mountain-local PostgreSQL catalog is used.

For local catalog operation, configure `db_user`, `db_pwd`, `db_host`, and
`db_database` in `shuffle.cfg`. Do not commit credentials or private
configuration files.

## Output products

Depending on the configuration and selected visualization options, a
successful shuffle can produce:

- `shout.result`: shuffled target coordinates and selected probe IDs.
- `shout.info`: shot metadata, displacement, calibration-star count, and
  failure code.
- `shout.probestars`, `shout.ifustars`, and `shout.acamstars`: selected star
  catalogs.
- `HETcommands/`: HET trajectory and probe-command files.
- `images/`: sky-view JPEG products.
- `ACAMimages/`: ACAM JPEG and CSV products.
- `regions/`: DS9 region files.
- `gp_wfs/`: guide and wavefront-sensor products when enabled.

The exact filenames and directories are controlled by the `[General]` and
`[directories]` sections of `shuffle.cfg`.

## Troubleshooting

- If catalog services time out, retry or select another supported catalog.
- If Gaia proper-motion lookup returns no data, hetshuffle retains the catalog
  coordinates and continues. A VizieR mirror can be selected with, for
  example, `HETSHUFFLE_VIZIER_SERVER=vizier.cfa.harvard.edu`.
- If no guide or WFS stars are found, increase `RADIUS` or adjust the
  magnitude and separation limits in `shuffle.cfg`.
- Keep `IFUSLOT` three digits, such as `058`; the command normalizes it
  internally.
- For target proper motion, provide `--epoch`, `--pmRA`, and `--pmDEC`.
- Image-service warnings do not necessarily mean the shuffle failed. Check
  the exit status, `shout.info`, and the generated image products. A fallback
  image service may still produce valid visualizations.
- DS9/XPA interaction is optional and is not required for ordinary tests.

## Development and testing

Run the local checks from an environment with the development extra installed:

```bash
python -m pre_commit run --all-files
make ci
```

The repository currently has 38 deterministic, network-free tests covering
configuration and CLI behavior, astrometry, ACAM geometry and rendering,
catalog selection helpers, proper-motion handling, focal-plane parsing, and
output formats. GitHub Actions runs the same Python 3.12 test workflow with
coverage reporting.

Live catalog, database, image-service, DS9/XPA, and full real-observation
checks remain manual or scheduled integration work. They should not be
required for ordinary pull requests.

## Contributing

New contributors should contact a repository maintainer and request GitHub
`Write` access before beginning work. Approved contributors clone the
canonical repository, work on feature branches, run the local hooks and
`make ci`, and submit pull requests against protected `main`.

See the [Git workflow](git-workflow.md),
[pull-request template](.github/pull_request_template.md), and
[branch-protection guidance](.github/branch-protection.md).

## License

HET Shuffle is distributed under the BSD 3-Clause License. See
[LICENSE](LICENSE).
