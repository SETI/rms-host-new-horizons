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
/seti/newnav/capped-run.sh pip install rms-host-new-horizons
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
  image. Every value comes from the FITS file, so the FITS file and its PDS3 label give
  the same `Snapshot`, and its path names the FITS file either way. The image array is
  the one the file stores; a calibrated image also carries its error estimates and
  quality flags. `astrometry=True` reads no image data.
- The exposure is centered on the middle of the exposure that the FITS header gives. The
  spacecraft clock counts at the start and end of the exposure are given as they appear
  in the PDS3 label.
- The pointing comes from the C kernels, or from the FITS header with
  `pointing='fits'`; the spacecraft's position comes from the SP kernels, or from the
  FITS header with `geom='fits'`. The camera distortion is modeled by a polynomial,
  which `fov_type='flat'` replaces with an undistorted field of view.
- `NewHorizonsLORRI.from_index()` makes one `Snapshot` for each row of a volume index.
  It merges in the columns of the supplemental index, matching its rows to the index's
  by VOLUME_ID and FILE_SPECIFICATION_NAME, and keeps the supplement's value of a column
  that both give.

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
