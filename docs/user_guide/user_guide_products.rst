==================
Supported Products
==================

Instrument
==========

The Long Range Reconnaissance Imager (LORRI) of New Horizons is a telescope with a CCD
of 1024 by 1024 pixels. The host reads the LORRI images of the 14 PDS3 volumes of the
mission, from the launch through the Jupiter flyby, the Pluto encounter, and the
Kuiper Belt extended missions. Each mission phase has two volumes: one of raw images,
such as NHJULO_1001, and one of calibrated images, such as NHJULO_2001. The host reads
the current version of each volume.

Products
========

Each LORRI product is a FITS file, such as ``lor_0030713591_0x633_sci.fit``, with a
detached PDS3 label beside it, such as ``lor_0030713591_0x633_sci.lbl``. The FITS file
holds a header of its own, so a product carries two labels. :doc:`user_guide_metadata`
describes how the host combines them.

Every image becomes an OOPS :class:`~oops.observation.Snapshot`, with the image in its
``data`` subfield and the image's lines and samples as its *v* and *u* axes. The array
is the primary image exactly as the FITS file stores it. A raw image holds a few columns
beyond the field of view, and ``data`` keeps them. A calibrated image also holds its
error estimates and its quality flags, in extensions named ``LORRI Error image`` and
``LORRI Quality flag image``, which the host gives as the ``error`` and ``quality``
subfields. A raw image has neither, and an astrometry read leaves all three arrays
unread.

Binning modes and fields of view
================================

The number of lines of the image decides the field of view. An image of 1024 lines is
unbinned, with a field of view of 1024 by 1024 pixels. An image of 256 lines is binned
4 by 4: each pixel combines a 4 by 4 block of the CCD, so its field of view is 256 by
256 pixels. An index row gives the mode in its ``BINNING_MODE`` column, as ``1X1`` or
``4X4`` in either case.

The ``fov_type`` option chooses how the field of view models the camera's distortion:

``"fast"``
    The default. A polynomial converts angles to pixels, and a second, precomputed
    polynomial converts pixels back to angles.

``"slow"``
    The same polynomial, inverted by Newton's method where the inverse is needed.

``"flat"``
    A field of view without distortion, 4.96 microradians per unbinned pixel.

Calibration
===========

A calibrated image carries a solar-normalized I/F calibration: one
:class:`~oops.calibration.FlatCalib` named ``"I/F"`` for an extended source, which
converts its data to I/F. Its factor is π r² / (t RSOLAR F), in which RSOLAR is the
FITS header's solar spectral response, the conversion from DN per second to the radiance
of a solar spectrum; t is the exposure duration; F is 176 erg/cm²/s/Å, the solar flux
at 1 AU at the pivot wavelength of LORRI; and r is the target's distance from the Sun in
AU at the middle of the exposure, from the SPICE kernels.

The FITS header of a calibrated image also gives the responses to other spectra:
RPLUTO, RCHARON, RJUPITER, RPHOLUS, and RMU69 for extended sources, and the P keywords
for point sources. The host does not use them, and they stay in the ``header_dict``
subfield for a caller who wants a calibration of another kind.

The observation has no calibration when the FITS header gives no RSOLAR, as for a raw
image, or for an astrometry read of a PDS3 label, which leaves the FITS header unread.
It also has none when the exposure is zero, and when the target has no ephemeris: no
registered path in OOPS, or no position in the SPICE kernels at that time. The host
loads no kernels for small bodies, so the minor planets and the Kuiper Belt objects that
the mission imaged have no calibration unless the caller furnishes their kernels and
defines their paths in OOPS. Index rows have none, because the tables do not give
RSOLAR. The ``calibrations`` override supplies a calibration in every case, and
``calibrations=[]`` turns it off.
