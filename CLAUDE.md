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
- The minimum rms-oops in `pyproject.toml` is 0.5.0, the first release with
  `oops.OopsSpiceError` (the `oops.Host` API first appears in 0.4.0). Until 0.5.0 is on
  PyPI, install rms-oops from its branch into the venv by hand.
- The type information is published in two stubs: `host_new_horizons/__init__.pyi`
  declares `NewHorizonsLORRI` in full, and `host_new_horizons/oops/__init__.pyi`
  re-exports it. `stubtest` checks them against the runtime in the check script and in CI,
  so adding, renaming or re-signing a public member means updating
  `host_new_horizons/__init__.pyi` in the same change. The `[tool.mypy] exclude`, the
  override list in `pyproject.toml` and `.stubtest-allowlist` keep it off the stub-less
  modules. The `type: ignore[override]` comment is a real mismatch with `oops.Host`:
  `from_file` names its first parameter `fileinfo`.
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
- What `from_file()` reads follows the operator's rules of 2026-11-18, shared by the four
  hosts: a FITS input reads the FITS file alone; a label input reads the label and the
  FITS file its `^IMAGE` names (`_fits_path()`, which tries the case of the label's own
  name first), except that an astrometry read of a label reads the label alone and needs
  no FITS file (and refuses `geom` or `pointing` from the header with `OopsValueError`). The
  path subfields always name the FITS file. The host opens FITS files itself
  (`_open_fits()`), because `oops.Host._read_fileinfo()` leaves a FITS file of another
  host open, and closes every HDUList once the arrays and headers are copied out of it,
  including one passed in. Labels are parsed by `_parse_label()`, which falls back from
  "fast" to "strict" to "loose".
- The metadata come from every source read. `_items_from_pds3()` and
  `_items_from_header()` give each source's items under PDS3 names (START_TIME in seconds
  TDB, the FITS start being SPCSCET - EXPTIME/2; STOPMET converted by `_sclk_from_met()`),
  leaving out missing and null values (`_is_missing()`, `_NULL_VALUES`), and
  `_merge_items()` merges them. Where both give an item, `_WINNING_SOURCE` names the one
  source, "label" or "header", whose value it takes for every image; each entry's comment
  cites the evidence and the issue that records the difference. Never add a rule that
  looks at the values of one image. The target body comes from the merged TARGET_NAME by
  the fixed translation in `_target_body()`. Every computed value (times, FOV, target)
  comes from the merged items. Nothing missing is None: a subfield or `dict` key with no
  value is absent, as are `error` and `quality` on a raw product. `dict` is the merged
  items plus every other key of each source; `label_dict`, `header_dict` and `headers`
  are the raw sources. All four, and `dict` of `from_index()`, must stay
  JSON-serializable (`_label_as_json()`, `_header_as_json()`, `_row_as_json()`). Labels
  are always parsed with `sources=False`; `label_dict` leaves out END, None values, and
  the parser's helper keys (lower-case suffixes such as `_fmt`, `_unit`, `_offset`, and
  `objects`); a repeated FITS keyword is keyed `NAME[n]` (`_header_cards()`), as the
  VICAR hosts key a repeated parameter. If no source gives the start time, it comes from
  SPACECRAFT_CLOCK_START_COUNT through the clock kernel (`_add_start_time_from_clock()`).
- By default, `calibrations` holds one `FlatCalib('I/F', ...)` when the FITS header gives
  RSOLAR (a calibrated product): pi r^2 / (texp RSOLAR 176), r the target's solar range
  in AU from the SP kernels at the Snapshot's midtime (`_iof_calibrations()`). The
  default is empty for a raw product, a label-only read, a zero exposure, the Sun (zero
  AU from itself), a target with no path or no ephemeris, and every index row. A
  `calibrations` override replaces it in `from_file()` and `from_index()` alike.
- With the pointing from the C kernels, `from_file()` evaluates the camera frame at
  mid-exposure (`_check_pointing()`) and raises `oops.OopsSpiceError` with the SPICE error
  code if the kernels give none; it catches only the classes cspyce raises. It does not
  check when `frame` is given or `pointing` is not "spice", and `from_index()` never
  checks, so that one image without pointing does not stop the whole index.
- `from_index()` follows the shared contract: `supplement` defaults to
  `oops.Host.SUPPLEMENT` (the one beside the index, which must then exist), and None
  raises `HostIndexError` unless `row_dicts` is given; `label=` is passed to
  `_read_index_rows()`; `select` other than None raises `OopsValueError`; any other
  keyword raises `OopsTypeError`. Each row's path subfields name the FITS file in full,
  from `_volumes_directory()` of the index's location, with `abspath` only when the file
  is local or cached. `_row_as_json()` does not apply the table reader's masks as they
  are, because the supplements declare zero as the null constant of vectors where a zero
  can be real (see the comment at `_TARGET_GEOMETRY`). Given
  `row_dicts`, it reads neither table. `oops.Host._read_index_rows()` pairs the rows by
  position when the two tables list the same files in the same order and otherwise by
  `VOLUME_ID` and `FILE_SPECIFICATION_NAME` (`file_key=_INDEX_KEYS`; the NHJULO_1001
  index lists `lor_0035089864_0x632_eng` last instead of at row 678, where the sorted
  supplement has it, and NHxxLO_1999 departs from its supplement's order at row 1924).
  It checks only those two columns (`validate=_INDEX_KEYS`), and the supplement, which
  follows the labels, wins every shared column (`prefer_supplement=True`; the index gives
  clock counts with a partition, rounds the stop count differently, gives
  PRODUCT_CREATION_TIME as a date alone, and shortens some TARGET_NAMEs).
  `oops.Host.from_index()` identifies the host from the index label and passes the call
  to `from_index()`, which reads the rows.
- `from_index()` returns an iterator. It checks its arguments and reads and filters the
  tables itself, then returns a generator that builds each row's Snapshot when the
  iteration reaches it, so only an error in one row waits for that row. It must not be a
  generator function, which would defer every error to the first `next()`.
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
