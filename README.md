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
observation geometry. It reads images from the Long Range Reconnaissance Imager (LORRI)
of New Horizons, as the Planetary Data System (PDS) archives them, and turns each one
into an OOPS `Snapshot`: the image array, together with the time, field of view,
pointing, target, and light source that OOPS needs to compute the image's geometry. It
is a product of the [PDS Ring-Moon Systems Node](https://pds-rings.seti.org).

# Features

- **Either file of a product**: an image can be opened from its FITS file or from its
  PDS3 label, and the pixels are identical either way.
- **Metadata from both labels**: when a PDS3 label is read together with its image
  file, the values in the label and in the FITS header are merged under the PDS3 names,
  with one fixed choice for each value on which they disagree.
- **Whole volumes at once**: a PDS3 index and its supplemental index give one `Snapshot`
  per image, without opening the images.
- **Pointing from SPICE or from the FITS header**: the pointing comes from the New
  Horizons C kernels, from the FITS header, or from a frame that the caller supplies.
- **Calibration to I/F**: a calibrated image carries the factor that converts its data
  to I/F.

# Installation

`rms-host-new-horizons` needs Python 3.11 or later. Install it from PyPI with:

```sh
pip install rms-host-new-horizons
```

The host also needs the New Horizons SPICE kernels. The
[setup chapter](https://rms-host-new-horizons.readthedocs.io/en/latest/user_guide/user_guide_setup.html)
of the user guide describes the kernels and the environment variables that locate them.

# Quick Start

Importing `host_new_horizons.oops` registers the host with `oops.Host`, which can then
open any LORRI image:

```python
import oops
import host_new_horizons.oops

obs = oops.Host.from_file('NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit')
```

The host class can also be used directly, for one image or for every image in a volume:

```python
from host_new_horizons import NewHorizonsLORRI

snapshot = NewHorizonsLORRI.from_file(
    'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.lbl')
snapshots = list(NewHorizonsLORRI.from_index('NHJULO_2001_index.lbl'))
```

Replace each path with the location of the file in your copy of the archive, or with
its URL. `from_index` also reads the supplemental index beside the index, which must
exist. It returns an iterator, which builds each image's Snapshot when the loop reaches
it; `list()` builds them all at once.

# Documentation

The [documentation](https://rms-host-new-horizons.readthedocs.io) holds the user guide
and the API reference. To build it locally, run `make html` in the `docs` directory after
installing the package with its `docs` extra.

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
