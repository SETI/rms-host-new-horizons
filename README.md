<!-- pyml disable MD025 -->
# rms-host-new-horizons

[![GitHub release; latest by date](https://img.shields.io/github/v/release/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/releases)
[![Test Status](https://img.shields.io/github/actions/workflow/status/SETI/rms-host-new-horizons/run-tests.yml?branch=main)](https://github.com/SETI/rms-host-new-horizons/actions)
[![Documentation Status](https://readthedocs.org/projects/rms-host-new-horizons/badge/?version=latest)](https://rms-host-new-horizons.readthedocs.io/en/latest/?badge=latest)
[![Code coverage](https://img.shields.io/codecov/c/github/SETI/rms-host-new-horizons/main?logo=codecov)](https://codecov.io/gh/SETI/rms-host-new-horizons)
<br />
[![PyPI - Version](https://img.shields.io/pypi/v/rms-host-new-horizons)](https://pypi.org/project/rms-host-new-horizons)
[![PyPI - Format](https://img.shields.io/pypi/format/rms-host-new-horizons)](https://pypi.org/project/rms-host-new-horizons)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/rms-host-new-horizons)](https://pypi.org/project/rms-host-new-horizons)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/rms-host-new-horizons)](https://pypi.org/project/rms-host-new-horizons)
<br />
[![GitHub commits since latest release](https://img.shields.io/github/commits-since/SETI/rms-host-new-horizons/latest)](https://github.com/SETI/rms-host-new-horizons/commits/main/)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/commits/main/)
[![GitHub last commit](https://img.shields.io/github/last-commit/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/commits/main/)
<br />
[![Number of GitHub open issues](https://img.shields.io/github/issues-raw/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/issues)
[![Number of GitHub closed issues](https://img.shields.io/github/issues-closed-raw/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/issues)
[![Number of GitHub open pull requests](https://img.shields.io/github/issues-pr-raw/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/pulls)
[![Number of GitHub closed pull requests](https://img.shields.io/github/issues-pr-closed-raw/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/pulls)
<br />
![GitHub License](https://img.shields.io/github/license/SETI/rms-host-new-horizons)
[![Number of GitHub stars](https://img.shields.io/github/stars/SETI/rms-host-new-horizons)](https://github.com/SETI/rms-host-new-horizons/stargazers)
![GitHub forks](https://img.shields.io/github/forks/SETI/rms-host-new-horizons)
<!-- start-after-point -->

# Introduction

`rms-host-new-horizons` is the New Horizons mission host for
[OOPS](https://github.com/SETI/rms-oops), the Object-Oriented Python SPICE toolkit for
observation geometry. It turns images from the Long Range Reconnaissance Imager (LORRI)
into OOPS `Snapshot` objects, with the frame, field of view, timing, target, and light
source that each image requires. It is a product of the
[PDS Ring-Moon Systems Node](https://pds-rings.seti.org).

# Installation

The `host_new_horizons` module is available via the `rms-host-new-horizons` package on
PyPI and can be installed with:

```sh
pip install rms-host-new-horizons
```

Geometry calculations need the New Horizons SPICE kernels. The host selects them through
the kernel database of `spicedb`, named by `SPICE_SQLITE_DB_NAME`, and reads them from
the directory named by the `SPICE_PATH` environment variable (or `$OOPS_RESOURCES/SPICE`).

# Getting Started

Importing `host_new_horizons.oops` registers the New Horizons host with `oops.Host`,
which can then identify and read any LORRI image, from either its FITS file or its PDS3
label:

```python
import oops
import host_new_horizons.oops

obs = oops.Host.from_file('NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit')
```

The host class can also be used directly. It reads a volume index together with the
supplemental index beside it (here `NHJULO_2001_supplemental_index.lbl`), which must
exist:

```python
from host_new_horizons import NewHorizonsLORRI

snapshot = NewHorizonsLORRI.from_file(
    'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.lbl')
snapshots = NewHorizonsLORRI.from_index('NHJULO_2001_index.lbl')
```

# Features

- `NewHorizonsLORRI.from_file()` makes a `Snapshot` from a raw or calibrated LORRI
  image, given its FITS file or its PDS3 label. Given the FITS file, it reads only that
  file. Given the label, it reads the label and the FITS file beside it, and takes the
  metadata from both. Its path names the FITS file either way. The image array is the
  one the file stores; a calibrated image also carries its error estimates and quality
  flags. `astrometry=True` reads no image data, and given a PDS3 label it reads the label
  alone, without the FITS file.
- Where the label and the FITS header disagree, each value always comes from the one of
  the two that is accurate for it. The start time, the exposure duration, and the image
  size come from the FITS header; the label gives the start time rounded to the
  millisecond. The target name and the spacecraft clock counts come from the label; the
  header's target name is the data pipeline's guess, such as "FJUPITER" for Jupiter. A
  value that neither gives is left out rather than set to None, and so are the archive's
  null values, such as "N/A" and "UNKNOWN".
- The target body comes from the target name by fixed translations of the label's forms
  of names: "J17 CALLIRRHOE" gives CALLIRRHOE, "136472 MAKEMAKE" gives MAKEMAKE, and
  "ASTEROID 486958 (2014 MU69)" gives ARROKOTH. Any other name is used as it is, whether
  or not SPICE knows it.
- The spacecraft clock counts at the start and end of the exposure are spelled as the
  PDS3 label spells them. They come from the label, or from the FITS header (its start
  count, and its end time converted to a count) when the label is not read or gives
  none. The metadata, the PDS3 label, and the FITS headers are published as dictionaries
  that can be written as JSON.
- An image whose pointing the C kernels do not cover raises `oops.OopsSpiceError`, with
  the SPICE error code in its message. Giving a frame, or taking the pointing from the
  FITS header, skips that check.
- The pointing comes from the C kernels, or from the FITS header with
  `pointing='fits'`; the spacecraft's position comes from the SP kernels, or from the
  FITS header with `geom='fits'`. The camera distortion is modeled by a polynomial,
  which `fov_type='flat'` replaces with an undistorted field of view.
- `NewHorizonsLORRI.from_index()` makes one `Snapshot` for each row of a volume index.
  It merges in the columns of the supplemental index, whose rows must list the same
  products, matched on VOLUME_ID and FILE_SPECIFICATION_NAME in any order. Where both
  tables give a column, the supplement's value is kept.

Details are available in the
[module documentation](https://rms-host-new-horizons.readthedocs.io/en/latest/module.html).

# Contributing

Information on contributing to this package can be found in the
[Contributing Guide](https://github.com/SETI/rms-host-new-horizons/blob/main/CONTRIBUTING.md).

# Links

- [Documentation](https://rms-host-new-horizons.readthedocs.io)
- [Repository](https://github.com/SETI/rms-host-new-horizons)
- [Issue tracker](https://github.com/SETI/rms-host-new-horizons/issues)
- [PyPi](https://pypi.org/project/rms-host-new-horizons)

# Licensing

This code is licensed under the [Apache License v2.0](https://github.com/SETI/rms-host-new-horizons/blob/main/LICENSE).
