# hetshuffle Testing and CI Implementation Plan

Status: proposed

This plan adapts the useful code-quality parts of the Panacea repository for
hetshuffle. It intentionally excludes Panacea's documentation, PDF, link-check,
and documentation-deployment workflows.

The desired outcome is that every pull request is automatically checked with a
general test suite before merge, while the test suite also protects the
scientific behavior that matters to hetshuffle: coordinate transformations,
focal-plane geometry, star selection, configuration handling, and generated
outputs.

## What to carry over from Panacea

Panacea provides a good collaboration pattern:

- `.github/workflows/python-tests.yml` runs on pushes to the main branch and on
  pull requests targeting it.
- The workflow installs the package and development dependencies, runs Ruff,
  runs pytest under coverage, and uploads coverage and failure artifacts.
- `pyproject.toml` owns the development dependencies through a `dev` extra.
- `pytest.ini` makes the test environment predictable.
- The `Makefile` mirrors the commands used by CI, so contributors can run the
  same checks locally.
- `.pre-commit-config.yaml` catches whitespace, merge-conflict markers, and
  common Python issues before code is pushed.
- The pull-request template asks for a summary, tests run, and useful logs or
  screenshots.

Panacea's test suite is intentionally lightweight and deterministic. Its
20 tests cover numerical utilities, parsers, packaged resources, CLI help,
synthetic data, and a pipeline smoke path. Running it locally produced
`20 passed`.

## Current hetshuffle gaps to address first

The repository currently has no test directory, CI workflow, development
dependency extra, or pull-request template. Before relying on CI, reconcile
the package metadata:

- `pyproject.toml` lists `reproject`, but the legacy `setup.py` and
  `requirements.txt` do not list it consistently.
- `pyproject.toml` declares Python `>=3.8`, `environment.yml` uses Python 3.11,
  and `setup.py` still advertises a Python 2.7 classifier. Choose and document
  the supported Python versions, then make the metadata agree.
- Keep one authoritative runtime dependency list where practical. If
  `setup.py` must remain for compatibility, keep it synchronized with
  `pyproject.toml` until it can be retired.
- Do not make tests depend on live Astroquery, Pan-STARRS, SkyView, database,
  DS9, or external HTTP services. Those are integration/manual checks, not
  reliable per-commit tests.

## Primary deployment target: the mountain environment

The most important operational user is the HET mountain environment. The
compatibility target is therefore the currently deployed Python 3.12
environment, rather than the developer laptop or the older Python versions
listed in legacy metadata.

The environment currently provides the core scientific stack:

- Python 3.12.
- NumPy 1.x, SciPy, Astropy, Astroquery, Matplotlib, and Six.

The environment currently does not provide `reproject`, `pyds9`, or
`psycopg2-binary`. This creates three decisions that should be resolved before
making CI a merge gate:

1. `reproject` is imported at module load time by `visualize.py`, so it is
   effectively required for importing the current visualization path. Either
   install it as part of the mountain deployment and CI, or move the import to
   the specific reprojection path and provide a clear error when that feature
   is requested without it.
2. `psycopg2` is imported at module load time by `findStars.py`, although the
   database-backed catalog path is only one part of the workflow. Either
   install it on the mountain or make database access an explicit optional
   dependency with a lazy import and an actionable error. The ordinary local
   catalog and geometry tests should not require a database client.
3. `pyds9` is already imported inside the DS9-specific function, so it is a
   good candidate to remain optional. Tests should verify that normal
   non-DS9 operation works without it and that the DS9 path reports its
   requirement clearly.

The plan should support the mountain environment without depending on exact
patch versions. CI should use Python 3.12 and compatible dependency ranges,
while a small documented smoke check should be run in the actual mountain
environment after deployment. That smoke check should verify package import,
configuration access, CLI help/version, and a network-free synthetic path.
It should not require live catalog access, DS9, or database credentials.

## Proposed implementation phases

### Phase 0: Define the supported environment

1. Make Python 3.12 the first CI target because it matches the mountain
   deployment. Do not pin CI to the exact patch version.
2. Decide whether Python 3.8 remains a supported secondary target. If it does,
   add it to a CI matrix only after the mountain-compatible suite is stable;
   otherwise raise the minimum version and remove stale Python 2 metadata.
3. Decide which of `reproject`, `psycopg2-binary`, and `pyds9` are required for
   the ordinary mountain workflow and which are optional. Make imports and
   dependency extras reflect that decision.
4. Reconcile runtime dependencies in `pyproject.toml`, `setup.py`,
   `requirements.txt`, and `environment.yml`.
5. Add a small developer dependency group, initially containing:

   ```text
   pytest
   coverage
   pytest-cov
   ruff
   ```

   Black, mypy, and pre-commit can be added after the first CI loop is
   working. They should not block the initial scientific test effort.

### Phase 1: Add a deterministic general test suite

Create `tests/` with fixtures that use small arrays, temporary directories,
and synthetic configuration files. The suite should be fast enough to run on
every commit and must not require observatory data or network access.

The initial test groups should be:

1. **Import and packaging smoke tests**
   - Import the package and key modules.
   - Verify the three console entry points can be constructed.
   - Verify packaged `shuffle.cfg`, `fplane.txt`, and `shuffle_version` are
     readable from an installed package.

2. **Configuration and CLI tests**
   - Parse representative valid command-line arguments.
   - Verify `--help` and `--version` exit successfully.
   - Load the default configuration when no local configuration is supplied.
   - Copy configuration files into a temporary directory and verify their
     contents exist.
   - Exercise invalid or missing configuration paths where the code promises a
     useful error.

3. **Astrometry and geometry tests**
   - Test `TangentPlane.raDec2xy()` and `xy2raDec()` as a round trip for small
     offsets, including a non-zero position angle.
   - Test `FPlaneCentroids` ID normalization, centroid caching, and lookup.
   - Test `_build_acam_wcs()` and `_sky_to_acam()` at the shuffle center and
     for known small offsets.
   - Test ACAM footprint boundaries and the orientation/sign convention.
   - Prefer physical invariants and a few documented reference values over
     large golden files.

4. **Star-selection and coordinate-helper tests**
   - Test angular-distance helpers, local/global coordinate conversions, and
     duplicate/close-source filtering with synthetic catalogs.
   - Test candidate filtering at magnitude, separation, and ACAM-boundary
     limits.
   - Test fallback behavior using mocked service responses rather than making
     catalog queries.

5. **Visualization smoke tests**
   - Use Matplotlib's non-interactive `Agg` backend.
   - Supply a tiny synthetic FITS image and configuration.
   - Mock `reproject_interp` and all network image retrieval.
   - Verify an ACAM chart, CSV, and DS9 region output can be created.
   - For source annotations, verify the source label is propagated and mapped
     to the expected short name.
   - Keep pixel-level appearance checks minimal; require a manual screenshot
     review for major changes to chart layout, orientation, or labels.

6. **Output and data-format tests**
   - Test representative `obj_dist_ifu_vis` output headers and rows.
   - Test star/catalog writer output with small in-memory examples.
   - Verify empty-result and failure-result formats remain readable.

The scientifically minimal first suite does not need to exercise every legacy
algorithmic variant in `findStars.py` or run a complete real observing-night
shuffle. It should establish the coordinate, selection, configuration, and
output invariants that would make a regression scientifically visible.

### Phase 2: Make local commands match CI

Add a `Makefile` with commands equivalent to the Panacea pattern:

```make
lint:
	ruff check .

test:
	python -m pytest -q

coverage:
	python -m coverage run -m pytest -q
	python -m coverage report -m
	python -m coverage xml -o coverage.xml

ci: lint coverage
```

The exact tab characters should be used in the actual Makefile. Contributors
should be able to run `make ci` locally and get substantially the same result
as GitHub Actions.

### Phase 3: Add GitHub Actions

Add `.github/workflows/python-tests.yml` with this behavior:

- Trigger on `push` to `main` and on `pull_request` targeting `main`.
- Use `actions/checkout@v4` and `actions/setup-python@v5`.
- Start with Python 3.12 to match the mountain deployment; expand the matrix
  only after additional support is intentional.
- Enable pip caching.
- Install the package with its development extra.
- Run `ruff check .`.
- Run the complete pytest suite under coverage.
- Print the coverage summary.
- Upload `coverage.xml` and pytest/coverage diagnostics when useful,
  especially on failure.
- Set `permissions: contents: read` unless a later workflow needs more.
- Use `pull_request`, not `pull_request_target`, for code execution from
  forks. Do not expose repository secrets to tests.
- Add concurrency cancellation so superseded commits do not consume runners.

The first import and CLI smoke tests should be representative of the mountain
environment. If optional dependencies have been made lazy, include one test
job that installs only the ordinary mountain dependency set and verifies that
the package and non-DS9 CLI paths still import successfully.

Do not initially make a coverage percentage a hard gate. First establish a
stable baseline; then introduce a modest threshold that prevents the number
from declining. Coverage is a guide, not a substitute for scientifically
meaningful tests.

### Phase 4: Establish merge and commit policy

Configure the `main` branch so that:

- Direct pushes are discouraged or protected.
- A passing Python test workflow is required before merge.
- Pull requests require at least one review when practical.
- Conversations and requested changes are resolved before merge.
- The PR branch is up to date with `main` before merging when a conflict or
  behavior change makes that important.

The effective workflow becomes:

```text
local edit
  -> local make ci
  -> focused commit
  -> push branch
  -> GitHub PR
  -> automatic lint + general tests + coverage
  -> scientific/manual review when applicable
  -> merge to main
```

This means a trusted collaborator can contribute normally, but trust does not
replace reproducible verification. Fork PRs should still run in a clean,
non-secret CI environment.

### Phase 5: Add pre-commit hygiene

After the test suite and CI are stable, add a Panacea-like
`.pre-commit-config.yaml` with:

- end-of-file fixing,
- trailing-whitespace checks,
- merge-conflict checks,
- Ruff,
- optionally Black once formatting policy is agreed.

Run the full test suite on pre-push only if it remains fast enough for normal
development. CI remains authoritative; hooks are a convenience and should not
make emergency recovery or scientific exploration unnecessarily difficult.

## Pull-request expectations

Add a lightweight `.github/pull_request_template.md` containing:

- Summary and motivation.
- Files or behavior changed.
- Tests run locally, including the exact command.
- Whether the change affects scientific calculations, coordinate conventions,
  catalog selection, or output formats.
- Screenshots, logs, or representative output for visualization changes.
- Related issue or discussion, if any.

For a new or untested feature, use a draft PR. A PR should become ready for
review only after the general suite passes and the contributor has described
any manual or data-dependent verification that remains.

## Recommended test priorities

### P0: required before the first CI merge gate

- Package import and entry-point smoke tests.
- Configuration loading and copying.
- CLI help/version/argument parsing.
- Tangent-plane round trips.
- ACAM WCS center/orientation tests.
- Focal-plane centroid normalization and lookup.
- Deterministic star-distance and candidate-selection tests.
- Network-free visualization output smoke test.

### P1: add as the corresponding code is touched

- Proper-motion coordinate updates.
- Parallactic-angle boundary cases.
- F-plane file parsing and representative IFU geometry.
- Catalog fallback and malformed-response handling.
- CSV, DS9-region, and failed-shuffle output formats.
- Full CLI execution with a tiny synthetic input fixture.

### Manual or scheduled integration checks

Keep these outside the per-commit suite unless they can be made deterministic:

- Live Pan-STARRS, Gaia, SkyView, or other catalog/image retrieval.
- Database-backed catalog access.
- DS9/XPA interaction.
- A complete real-observation shuffle using large or private data.
- Visual inspection of a representative finder chart.

These can be run manually before a release or by a separately triggered,
credentialed workflow. They should never be required for ordinary fork PRs.

## Suggested implementation order

Use small, reviewable commits rather than one large infrastructure change:

1. Reconcile dependency metadata and define the supported Python version.
2. Resolve the mountain dependency boundary for reprojection, database
   catalogs, and DS9.
3. Add development dependencies and basic pytest configuration.
4. Add P0 deterministic tests and synthetic fixtures.
5. Add the Makefile commands and run them locally.
6. Add the GitHub Actions test workflow.
7. Add the PR template and branch-protection settings.
8. Add pre-commit hooks and P1 tests incrementally.

Each step should leave `main` installable and the general test suite passing.

## Implementation checklist

- [x] 1. Reconcile dependency metadata and define Python 3.12+ as the supported
  deployment floor. Verify the package metadata, environment definition, and
  requirements are consistent with the mountain target.
- [ ] 2. Resolve the mountain dependency boundary for `reproject`, database
  catalogs, and DS9.
- [ ] 3. Add development dependencies and basic pytest configuration.
- [ ] 4. Add P0 deterministic tests and synthetic fixtures.
- [ ] 5. Add Makefile commands that mirror CI.
- [ ] 6. Add the GitHub Actions test workflow.
- [ ] 7. Add the pull-request template and configure branch protection.
- [ ] 8. Add pre-commit hooks and P1 tests incrementally.

Current checkpoint: item 1 is complete and committed locally. The branch is
not being pushed until the remaining implementation items are finished.
