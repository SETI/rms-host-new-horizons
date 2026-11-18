# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`rms-host-new-horizons` is the New Horizons host plug-in for `oops` (rms-oops). It is a port
of `rms-oops/src/oops/hosts/newhorizons/` to the newer `oops.Host` / `oops.host` API, part
of splitting the mission hosts out of rms-oops (see the siblings `rms-host-cassini`,
`rms-host-galileo` and `rms-host-voyager`). LORRI is ported.

## Standards

Tooling and conventions follow **rms-polymath**. `.claude/rules/*.md`, which were copied from
it, are the authoritative coding, testing, documentation, dependency and environment
standards. The repo layout, CI and scripts come from `rms-devenv/repo_template`, adapted as
rms-polymath adapted them, by way of `rms-host-cassini` and `rms-host-galileo`. Model
changes to packaging, CI or scripts on rms-polymath's `pyproject.toml`, `.flake8`,
`.github/workflows/` and `scripts/`, and keep them in step with the sibling hosts. The
rules' references to the two published `.pyi` stubs and `stubtest` apply here to the stubs
described below.

- Hand-aligned columns (assignments, dict colons, `from x     import`) are deliberate. Never
  run `ruff format` or any other auto-formatter.
- Keep the file banners: a line of 90 `#` characters at the top and bottom of each module,
  with a `# host_new_horizons/<file>.py` header (`.pyi` for a stub).

## Current state

- `scripts/run-all-checks.sh` is the single source of truth for the checks: ruff, flake8
  `E12`/`E13`, mypy on `tests/`, pytest, pyroma, stubtest, pip-audit, the Sphinx build
  (`-W -n`), codespell and PyMarkdown. `ruff format` is disabled.
  `.github/workflows/run-tests.yml` runs the same set. `scripts/setup-venv.sh` builds the
  `venv/` that the scripts expect; set `VENV` to use another.
- The `oops.Host` API first appears in rms-oops 0.4.0, the minimum version in
  `pyproject.toml`.
- The type information is published in two stubs: `host_new_horizons/__init__.pyi`
  declares `NewHorizonsLORRI` in full, and `host_new_horizons/oops/__init__.pyi`
  re-exports it. `stubtest` checks them against the runtime in the check script and in CI,
  so adding, renaming or re-signing a public member means updating
  `host_new_horizons/__init__.pyi` in the same change. The `[tool.mypy] exclude`, the
  override list in `pyproject.toml` and `.stubtest-allowlist` keep it off the stub-less
  modules. The two `type: ignore[override]` comments are real mismatches with `oops.Host`:
  `from_file` names its first parameter `fileinfo`, and `from_index` returns
  `list[Snapshot]`.
- `docs/conf.py` resolves the short class names in docstrings (`FCPath`, `Observation`,
  ...) through intersphinx with `_SHORT_NAMES`. Add a name there when a nitpicky build
  reports a class from another RMS package.
- Coverage must stay at or above 65% (`fail_under` in `pyproject.toml`, and the project and
  patch targets in `codecov.yml`); it is at 100%. The ubuntu / Python 3.13 CI job uploads
  to Codecov and fails if the upload fails. The upload authenticates through GitHub OIDC
  (`use_oidc`), so it needs no `CODECOV_TOKEN` secret.
- Tests never touch SPICE kernels or spicedb. `tests/oops/conftest.py` replaces them with
  fakes for every test, so new tests get this automatically. Synthetic FITS files and
  labels are written with `write_fits()` and `write_label()` there. `test_files/` holds
  cut-down archive products, described in `tests/oops/test_lorri_files.py`.

## Architecture gotchas

- The oops support lives in the `host_new_horizons.oops` subpackage: `_NewHorizons` in
  `_new_horizons.py` and `NewHorizonsLORRI` in `lorri.py`. `NewHorizonsLORRI` calls
  `_register()` at import time, and `oops/__init__.py` imports it. `Host.from_file()`
  dispatches through the registered hosts' `_detect_in_*` methods, so the host is
  invisible until `host_new_horizons.oops` (or its module) is imported.
- `host_new_horizons/__init__.py` imports the hosts lazily, through a PEP 562
  `__getattr__`, so that importing `host_new_horizons` does not import oops. Keep any
  import of `host_new_horizons.oops` out of `host_new_horizons/__init__.py`.
- `_NewHorizons` is never instantiated; all its state lives in class attributes. Its
  `initialize()` furnishes every kernel through spicedb at once, for the whole period
  `START_TIME` to `STOP_TIME`, or a metakernel instead; there is no monthly loading.
- Every value of a `from_file()` Snapshot comes from the FITS file. Given a PDS3 label,
  the host only follows its `^IMAGE` pointer to the FITS file (`_fits_path()`, which tries
  the case of the label's own name first), so the two inputs give identical Snapshots and
  the path subfields always name the FITS file. The host opens FITS files itself
  (`_open_fits()`), because `oops.Host._read_fileinfo()` leaves a FITS file of another
  host open, and closes every HDUList once the arrays and headers are copied out of it,
  including one passed in.
- The clock counts come from the FITS header: `SPCSCLK0` for the start and `STOPMET`,
  converted by `_sclk_from_met()`, for the end. Both match the PDS3 label's spelling on
  every archive product checked.
- `from_index()` reads the team's volume index with the supplement beside it
  (`oops.Host.SUPPLEMENT`), which must exist. The tables write their shared columns
  differently (the index gives clock counts with a partition and rounds the stop count
  differently), so `_match_index_rows()` matches rows on `VOLUME_ID` and
  `FILE_SPECIFICATION_NAME` and copies the supplement's row into the index's before the
  merge, so that the supplement, which follows the labels, wins every shared column.
  `oops.Host.from_index()` reads the rows itself with the default validation, which
  compares every shared column, so it fails on these indexes.
- The camera frame `NH_LORRI` is a `Cmatrix` that flips the y and z axes of the SPICE
  frame, registered as `NH_LORRI_FLIPPED`, because the instrument kernel gives the
  boresight along -Z. `_define_camera_frame()` registers both lazily, and only when the
  pointing comes from the C kernels.
- The distortion model (`_LORRI_COEFF` and its numerical inverse `_LORRI_INV_COEFF`) is
  fixed in `lorri.py`; the instrument kernel supplies only the image shape.
- Running anything needs SPICE kernels through `spicedb` (`SPICE_PATH`, or
  `$OOPS_RESOURCES/SPICE`, plus `SPICE_SQLITE_DB_NAME`). `spicedb` currently ships inside
  rms-oops.

## Repo etiquette

- Commit subjects are a plain, capitalized, imperative sentence with no type prefix and no
  trailing period. PRs are squash-merged onto `main`.
- Branch names: `<initials>_<YYMMDD>_<topic>`.
- Versions come from `setuptools_scm`; never hand-write a version.
