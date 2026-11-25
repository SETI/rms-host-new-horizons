# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`rms-host-new-horizons` is the New Horizons host plug-in for `oops` (rms-oops). It reads
LORRI images and the Ring-Moon Systems Node metadata indexes through the `oops.Host` /
`oops.host` API. Its siblings are `rms-host-cassini`, `rms-host-galileo` and
`rms-host-voyager`, which read alike.

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
- Until rms-oops has a release with the host contract that this package uses (the
  streaming `_index_rows`, `_set_path_attributes`), `pyproject.toml` installs rms-oops
  from the branch of its open pull request, so that CI can run. That line becomes a
  versioned requirement once the release is out, and before this package is published,
  because PyPI refuses a package with a direct URL dependency.
- The type information is published in two stubs: `host_new_horizons/__init__.pyi`
  declares `NewHorizonsLORRI` in full, and `host_new_horizons/oops/__init__.pyi`
  re-exports it. `stubtest` checks them against the runtime in the check script and in CI,
  so adding, renaming or re-signing a public member means updating
  `host_new_horizons/__init__.pyi` in the same change. The `[tool.mypy] exclude`, the
  override list in `pyproject.toml` and `.stubtest-allowlist` keep it off the stub-less
  modules. The two `type: ignore[override]` comments are real mismatches with
  `oops.Host`: `from_file` names its first parameter `fileinfo`, and neither `from_file`
  nor `from_index` takes the `**kwargs` that `oops.Host` declares.
- `docs/conf.py` resolves the short class names in docstrings (`FCPath`, `Observation`,
  ...) through intersphinx with `_SHORT_NAMES`. Add a name there when a nitpicky build
  reports a class from another RMS package.
- Coverage must stay at or above 65% (`fail_under` in `pyproject.toml`, and the project and
  patch targets in `codecov.yml`); it is at 100%. The ubuntu / Python 3.13 CI job uploads
  to Codecov and fails if the upload fails. The upload authenticates through GitHub OIDC
  (`use_oidc`), so it needs no `CODECOV_TOKEN` secret.
- Tests never touch SPICE kernels or spicedb. `tests/oops/conftest.py` replaces them with
  fakes for every test, so new tests get this automatically. Each test also gets its own
  empty file cache under `tmp_path` (the `file_cache` fixture). Synthetic FITS files and
  labels are written with `write_fits()`, `write_label()` and `write_product()` there.
  `test_files/` holds cut-down archive products and a cut-down metadata index with its
  supplement.

## Architecture gotchas

- The oops support lives in the `host_new_horizons.oops` subpackage: `_NewHorizons` in
  `_new_horizons.py` and `NewHorizonsLORRI` in `lorri.py`. `NewHorizonsLORRI` calls
  `_register()` at import time, and `oops/__init__.py` imports it. `Host.from_file()`
  dispatches through the registered hosts' `_detect_in_*` methods, so the host is
  invisible until `host_new_horizons.oops` (or its module) is imported. Nothing loads
  kernels at import.
- `host_new_horizons/__init__.py` imports the hosts lazily, through a PEP 562
  `__getattr__`, so that importing `host_new_horizons` does not import oops. Keep any
  import of `host_new_horizons.oops` out of `host_new_horizons/__init__.py`.
- `_NewHorizons` is never instantiated; all its state lives in class attributes. Its
  `initialize()` furnishes every kernel through spicedb at once, for the whole period
  `START_TIME` to `STOP_TIME`, or a metakernel instead. `NewHorizonsLORRI.initialize()`
  (keywords `asof` and `meta` only) calls it, builds the FOVs, and registers the camera
  frame `NH_LORRI`, a `Cmatrix` that flips the y and z axes of the SPICE frame
  registered as `NH_LORRI_FLIPPED`. It is idempotent, and `from_file()` and
  `from_index()` call it.
- Every read of a label or a FITS file goes through one function, `_read()`, which
  reads the files (`_read_and_interpret()`), checks them (a FITS file that ends early,
  a primary array that is not one band, an array whose shape is not the label's IMAGE
  LINES and LINE_SAMPLES), extracts every metadata item and header value the
  observation needs, and turns any other failure (KeyError, ValueError, TypeError,
  IndexError, OSError other than FileNotFoundError, astropy's truncation warning) into
  `HostFormatError` naming the file passed. Code that calls SPICE stays outside it, so
  that a cspyce error is never relabeled. `_read()` returns a `_LorriFiles` namedtuple.
- What `from_file()` reads: a FITS input reads the FITS file alone; a label input reads
  the label and the FITS file its `^IMAGE` names, written in the case of the label's own
  name (`_fits_path()`), except that an astrometry read of a label reads the label
  alone. Labels are parsed with the strict method (`_DEFAULT_PDS3_OPTIONS`), which
  every current LORRI label passes as written, so the host keeps
  `oops.Host._correct_label_text()` as it is. There is no parser fallback.
- The metadata come from every source read. `_items_from_pds3()` (labels and index
  rows), `_items_from_label()` and `_items_from_header()` give each source's items
  under PDS3 names, leaving out missing values (`_is_missing()`: None,
  `oops.host.is_pds3_null()`, the FITS fill -1e32, and lists or sets of those), and
  `_merge_items()` merges them by `_WINNING_SOURCE`, whose comments cite the evidence
  and the label-differences issue. Never add a rule that looks at the values of one
  image. `dict` is the merged items (START_TIME as a `datetime`, the type the label
  gives) plus every other key of each source without missing values; `label_dict` is
  `Pds3Label.as_dict()` as it is, `header_dict` is `dict()` of the primary astropy
  Header, and `headers` is the list of astropy Headers. Nothing is reshaped for JSON.
- The I/F calibration (`_iof_calibrations()`) is a solar-normalized `FlatCalib` from
  RSOLAR, built outside `_read()` because it needs the target's ephemeris. It is
  omitted without RSOLAR, for a zero exposure, and for a target without a registered
  path (KeyError) or without SP kernel data at mid-exposure (cspyce's OSError). There
  is no spacecraft-distance fallback and no small-body kernel loading.
- No pointing check anywhere: a missing C kernel shows up when the frame is used.
- `from_index()` reads the metadata index and its supplement through
  `oops.Host._index_rows()` with `validate=()`; `MERGE_SUPPLEMENT_BY_FILE` is True, on
  the default key (VOLUME_ID, FILE_SPECIFICATION_NAME). The index value wins where
  both tables give a column. Rows are not repaired or checked: `dict` is the row as the
  tables give it, `_mask` entries included, and `filepath` is the relative
  `VOLUME_ID/FILE_SPECIFICATION_NAME`. BINNING_MODE is looked up without regard to
  case. The supplement is required (`supplement=None` raises `HostIndexError` unless
  `row_dicts` is given). The method returns a generator expression after the call-time
  work, so argument and label errors are raised at the call.
- `select` goes through `oops.Host._apply_select()` on the one-item tuple
  (`_selected()`); None returns the observation itself.
- Running anything needs SPICE kernels through `spicedb` (`SPICE_PATH`, or
  `$OOPS_RESOURCES/SPICE`, plus `SPICE_SQLITE_DB_NAME`). `spicedb` currently ships inside
  rms-oops.

## Repo etiquette

- Commit subjects are a plain, capitalized, imperative sentence without a type prefix or a
  trailing period. PRs are squash-merged onto `main`.
- Branch names: `<initials>_<YYMMDD>_<topic>`.
- Versions come from `setuptools_scm`; never hand-write a version.
