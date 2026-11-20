==================
Supported Products
==================

Instrument
==========

The Long Range Reconnaissance Imager (LORRI) of New Horizons is a telescope with a CCD
of 1024 by 1024 pixels. The host reads the LORRI images of the 14 PDS3 volumes of the
mission, from the launch through the Jupiter flyby, the Pluto encounter, and the
Kuiper Belt extended missions. Each mission phase has two volumes: one of raw images,
such as NHJULO_1001, and one of calibrated images, such as NHJULO_2001.

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
error estimates and its quality flags, which the host gives as the ``error`` and
``quality`` subfields. A raw image has neither, and an astrometry read leaves all three
arrays unread.

Binning modes and fields of view
================================

The binning mode, given by the ``binning_mode`` subfield, decides the field of view. An
image of 1024 lines is unbinned, ``1X1``, with a field of view of 1024 by 1024 pixels.
An image of 256 lines is binned ``4X4``: each pixel combines a 4 by 4 block of the CCD,
so its field of view is 256 by 256 pixels. The host takes the mode from the number of
lines, and from the ``BINNING_MODE`` column for an index row.

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

A calibrated image carries one :class:`~oops.calibration.FlatCalib` named ``"I/F"``,
which converts its data to I/F. Its factor is π r² / (t RSOLAR F), in which RSOLAR is
the FITS header's conversion from DN per second to the radiance of a solar spectrum, t
is the exposure duration, F is 176 erg/cm²/s/Å, the solar flux at 1 AU at the pivot
wavelength of LORRI, and r is the target's distance from the Sun in AU at the middle of
the exposure, from the SPICE kernels.

The observation's list of calibrations is empty for a raw image, whose FITS header lacks
RSOLAR, and for an astrometry read of a PDS3 label, which leaves the FITS header unread.
It is also empty
when the exposure is zero, when the target is the Sun or ``NONE``, when the kernels
cannot place the target at that time, and for every index row. The caller can supply a
calibration with the ``calibrations`` override in every case.
