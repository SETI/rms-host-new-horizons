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
- Until rms-oops 0.5.0, the first release with `oops.OopsSpiceError` and the index
  reader the host uses, is on PyPI, `pyproject.toml` installs rms-oops from the branch
  of its open pull request, so that CI can run. That line becomes `rms-oops>=0.5.0`
  once the release is out, and before this package is published, because PyPI refuses
  a package with a direct URL dependency.
- The type information is published in two stubs: `host_new_horizons/__init__.pyi`
  declares `NewHorizonsLORRI` in full, and `host_new_horizons/oops/__init__.pyi`
  re-exports it. `stubtest` checks them against the runtime in the check script and in CI,
  so adding, renaming or re-signing a public member means updating
  `host_new_horizons/__init__.pyi` in the same change. The `[tool.mypy] exclude`, the
  override list in `pyproject.toml` and `.stubtest-allowlist` keep it off the stub-less
  modules. The two `type: ignore[override]` comments are real mismatches with
  `oops.Host`: `from_file` names its first parameter `fileinfo`, neither `from_file`
  nor `from_index` takes the `**kwargs` that `oops.Host` declares, and both type
  `select` as None alone, since any other value raises.
- `docs/conf.py` resolves the short class names in docstrings (`FCPath`, `Observation`,
  ...) through intersphinx with `_SHORT_NAMES`. Add a name there when a nitpicky build
  reports a class from another RMS package.
- Coverage must stay at or above 65% (`fail_under` in `pyproject.toml`, and the project and
  patch targets in `codecov.yml`); it is at 100%. The ubuntu / Python 3.13 CI job uploads
  to Codecov and fails if the upload fails. The upload authenticates through GitHub OIDC
  (`use_oidc`), so it needs no `CODECOV_TOKEN` secret.
- Tests never touch SPICE kernels or spicedb. `tests/oops/conftest.py` replaces them with
  fakes for every test, so new tests get this automatically. Each test also gets its own
  empty file cache under `tmp_path` (the `file_cache` fixture), so whether a remote file
  counts as cached never depends on the shared cache. Synthetic FITS files and
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
  path subfields always name the FITS file. A path is opened or parsed by
  `oops.Host._read_fileinfo()` (`formats='PF'`, `pds3_fallback=True`): a ".fit" or
  ".fits" path is a FITS file, closed by oops if it is not a LORRI file; a ".lbl" path,
  or another file that begins as a PDS3 label, is parsed by the requested method and then
  each later one of "fast", "strict" and "loose"; any other file is a WrongHostError. The
  host closes every HDUList once the arrays and headers are copied out of it, including
  one passed in. LORRI labels parse as they are written, so the host keeps
  `oops.Host._correct_label_text()` as it is.
- The metadata come from every source read. `_items_from_pds3()` and
  `_items_from_header()` give each source's items under PDS3 names (START_TIME in seconds
  TDB, the FITS start being SPCSCET - EXPTIME/2; STOPMET converted by `_sclk_from_met()`),
  leaving out missing and null values (`_is_missing()`, which calls
  `oops.host.is_pds3_null()` and adds the FITS fill value -1e32), and
  `_merge_items()` merges them. Where both give an item, `_WINNING_SOURCE` names the one
  source, "label" or "header", whose value it takes for every image; each entry's comment
  cites the evidence and the issue that records the difference. Never add a rule that
  looks at the values of one image. The target body comes from the merged TARGET_NAME by
  the fixed translation in `_target_body()`. Every computed value (times, FOV, target)
  comes from the merged items. A missing value is never None: a subfield or `dict` key
  whose value is missing is absent, as are `error` and `quality` on a raw product. `dict` is the merged
  items plus every other key of each source; `label_dict`, `header_dict` and `headers`
  are the raw sources. All four, and `dict` of `from_index()`, must stay
  JSON-serializable (`_label_as_json()`, `_header_as_json()`, `_row_as_json()`). Labels
  are always parsed with `sources=False`; `label_dict` leaves out END, None values, and
  the parser's helper keys (lower-case suffixes such as `_fmt`, `_unit`, `_offset`, and
  `objects`); a repeated FITS keyword is keyed `NAME[n]` (`_header_cards()`), as the
  VICAR hosts key a repeated parameter. If the sources lack the start time, it comes from
  SPACECRAFT_CLOCK_START_COUNT through the clock kernel (`_add_start_time_from_clock()`),
  and the `time_from_sclk` subfield, which every Snapshot has, is True.
- By default, `calibrations` holds one `FlatCalib('I/F', ...)` when the FITS header gives
  RSOLAR (a calibrated product): pi r^2 / (texp RSOLAR 176), r the target's solar range
  in AU from the SP kernels at the Snapshot's midtime (`_iof_calibrations()`). The
  default is empty for a raw product, a label-only read, a zero exposure, the Sun (zero
  AU from itself), a target without a path or an ephemeris, and every index row. A
  `calibrations` override replaces it in `from_file()` and `from_index()` alike.
- With the pointing from the C kernels, `from_file()` checks the camera frame before
  the overrides are applied (`_check_pointing()`, which passes the cadence and the
  overrides to `oops.Host._check_spice_frame()` with `defer=False`), and raises
  `oops.OopsSpiceError` with the SPICE error code if the kernels give none. oops builds
  the frame and the time as `_apply_overrides()` will and evaluates the frame relative
  to J2000, so the check reaches the SPICE frame `NH_LORRI_FLIPPED` beneath the camera's
  Cmatrix, at the midtime after the `cadence` override and a "T" `timeshift`. An "F"
  `timeshift` wraps the Cmatrix in a FrameShift, which shifts only the Cmatrix's own
  fixed rotation, so the kernels are still read at that midtime. A `tracker` override
  reads them once instead, at the start, midtime, or end that it selects (oops's
  `_tracker_epoch()`). Whether "F" should shift the whole chain is an open question in rms-oops (cited in the `#` comment in
  `_check_pointing()`); the host follows what oops does. It does not check when `frame`
  is given or `pointing` is not "spice". `from_index()` makes the same check for each
  row with `defer=True`: a row without pointing gets an `oops.frame.SpiceErrorFrame`,
  which raises `oops.OopsSpiceError` when it is evaluated, so that one image without
  pointing does not stop the whole index. In the tests, the fake SPICE frame is a class
  of its own (`FakeSpiceFrame` in `tests/oops/conftest.py`), so that a test can take its
  data away while the camera's Cmatrix still works, and `spice_check_times()` in
  `tests/oops/test_lorri_merge.py` records the times at which the check evaluates it.
- `from_file()` and `from_index()` both refuse a `select` other than None
  (`oops.Host._check_select()`, `OopsValueError`) at call time; neither takes
  `**kwargs`, so Python's own `TypeError` reports an unknown keyword. `abspath` is given
  only by a read of the FITS file's arrays: `from_file()` passes the HDUList's file
  (for a remote file, the path that oops's `FCPath.retrieve()` returned) as `abspath=` to
  `oops.Host._apply_overrides()`. `_hdulist_path()` gives that file for `abspath` and
  for the `filepath` of an HDUList passed in: `hdulist.filename()` if it names a file,
  and None for an HDUList built in memory or read from a stream, which astropy names
  after the stream's type ("<class '_io.BytesIO'>"). An astrometry read and `from_index()` leave it out, so
  their Snapshots lack `abspath`, even for a local or cached file. The host has no
  path code of its own (each row's file is `oops.Host._index_file_path()`
  with the FITS suffix), and reading never creates a directory. Parsed content that
  was not read from a file has no path subfields; a full read of such a Pds3Label raises
  `HostFormatError`.
- `from_index()` reads the index through `oops.Host._index_rows()`, the only code that
  reads a table. The call reads only the two labels (given the `.tab` path instead of
  the label, pdsparser reads the whole table looking for an attached label before it
  turns to the `.lbl` beside it); the iteration reads both tables 100
  rows at a time (oops's `CHUNK_ROWS`, through pdstable's `iter_dicts`) and builds each
  row's Snapshot when it reaches the row. Supplement rows read ahead of their index
  rows are held until their partners come; in the 16 metadata index pairs at most one
  is held at a time (2026-11-21 survey). Memory stays the same however long the index
  is only with the frame and path cache fix of SETI/rms-oops#252, which the CI does not
  install: with it, peak RSS was 163 MB over the 1,119 rows of NHJULO_1001, 166 MB over
  the 21,753 rows of NHxxLO_1999, and 166 MB over the 9,641 rows of NHKELO_2001, every
  one of which gets a SpiceErrorFrame (2026-11-19, rms-oops 191800d). Without it, each
  row keeps its camera frame in the oops caches, a few KB per row (about 13 KB for a
  row with a SpiceErrorFrame). It is
  not a generator function, which would defer the argument and label errors to the
  first `next()`.
  `supplement` defaults to `oops.Host.SUPPLEMENT` (the one beside the index, which must
  then exist), and None or an empty string, which oops would take as "no supplement",
  raises `HostIndexError` unless `row_dicts` is given.
  `_index_rows()` returns given `row_dicts` before it reads either label, so a direct
  call with `row_dicts` skips every label check (existence, parsing, tables, host,
  the "_index" name, `pds3_options`); only `select`, `fov_type`, `filter` (against the
  first row) and each row are checked. `oops.Host.from_index()` still parses the index
  label to pick the host.
  The host reads the Ring-Moon Systems Node metadata indexes and their supplements.
  Each row's path subfields name the FITS file in full, under
  `oops.Host._volumes_directory()` of the index's location.
- `MERGE_SUPPLEMENT_BY_FILE` is True: rows are paired by `VOLUME_ID` and
  `FILE_SPECIFICATION_NAME` (`file_key=_INDEX_KEYS`), because the NHJULO_1001 index
  lists `lor_0035089864_0x632_eng` out of place, last, instead of at row 679 (counting
  from 1), where the supplement, which is in file order, has it; in the cumulative
  NHxxLO_1999 the two tables differ from row 1,925 on for the same reason. The other 14 of the 16 metadata
  indexes are in step, and pairing by file gives them the same rows. The key must carry
  `VOLUME_ID`, which oops's default key (FILE_SPECIFICATION_NAME alone) does not: the
  cumulative NHxxLO_1999 and NHxxLO_2999 tables list four products twice, once for
  NHKELO and once for NHK2LO, under one FILE_SPECIFICATION_NAME. The explicit key also
  makes oops refuse a table without VOLUME_ID at call time. A row of `row_dicts` (or a
  null value) without either key column raises `HostFormatError` when the iteration
  reaches it (`_snapshot_from_row()`).
  Only the file key is compared (`validate=()`), and the supplement wins every shared column
  (`prefer_supplement=True`; the index gives clock counts with a partition, rounds the
  stop count differently, gives PRODUCT_CREATION_TIME as a date alone, and shortens some
  TARGET_NAMEs). The supplement is not the label either: it rounds EXPOSURE_DURATION to
  the millisecond and drops some SEQUENCE_ID prefixes.
- The index files can be corrected, unlike the product labels; their problems are listed
  in one issue in SETI/rms-data-projects, which the `#` comment beside each workaround
  cites by number and item. Every workaround for one must still give the right
  answer on a corrected table: pairing by file works in any order; `_row_as_json()`
  treats zeros specially only where the reader masks them as the declared null constant
  (see the comment at `_TARGET_GEOMETRY`), and `_VECTOR_COMPONENT` accepts QUATERNION as
  well as the supplement's QUATERNIAN; `_SUPPLEMENT_TARGET_NAMES` translates only the
  five bare moon names ("IO" to "J1 IO", as every label writes it), and leaves a name
  already in the label's form alone. Each has a test with a corrected value.
  `oops.Host.from_index()` identifies the host from the index label and passes the call
  to `from_index()`.
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

- Commit subjects are a plain, capitalized, imperative sentence without a type prefix or a
  trailing period. PRs are squash-merged onto `main`.
- Branch names: `<initials>_<YYMMDD>_<topic>`.
- Versions come from `setuptools_scm`; never hand-write a version.
