=============================
Time, Pointing, and Targets
=============================

Time
====

``START_TIME`` is the start of the exposure. The exposure lasts ``EXPOSURE_DURATION``,
and the Snapshot's default cadence covers exactly that interval. From the FITS header,
the start is ``SPCSCET``, the middle of the exposure in seconds TDB, less half of
``EXPTIME``. The PDS3 label gives the same start in UTC, rounded to the millisecond, and
the header's value is used when both are read (see :doc:`user_guide_metadata`). An
exposure of zero gives a Snapshot whose exposure is zero.

Pointing
========

By default the camera's pointing comes from the C kernels, through the SPICE frame
``NH_LORRI``. The instrument kernel puts the camera's boresight along the frame's -*z*
axis, while OOPS expects it along +*z*. The host therefore registers the SPICE frame as
``NH_LORRI_FLIPPED`` and defines the camera frame ``NH_LORRI`` from it, with its *y*
and *z* axes reversed. The ``spice_to_frame`` subfield holds that rotation.

The host does not check that the C kernels hold the pointing of an image. An image
without pointing still gives a Snapshot, and its frame raises the SPICE error, such as
``SPICE(NOFRAMECONNECT)``, when it is used.

The pointing can also come from the FITS header, with the ``pointing`` option (see
:doc:`user_guide_files`), and the spacecraft's position likewise with the ``geom``
option. A caller that has better pointing can pass it with the ``frame`` override, for
example as a :class:`~oops.frame.Cmatrix`, which is used as it is given.

Three subfields describe the pointing for a caller that works with the C kernels
directly: ``ck_object_id`` is -98000, the spacecraft whose attitude the C kernels give;
``ck_sclk_id`` is -98, the spacecraft clock of those kernels; and
``spice_frame_is_frozen`` is False, because the camera frame turns with the spacecraft
and is never fixed to one attitude.

Targets and light sources
=========================

The target body comes from ``TARGET_NAME``. The PDS3 label gives the name in the form
of the PDS catalog, which is not always the name that SPICE knows the body by, so the
host translates it by fixed rules:

.. list-table::
   :header-rows: 1

   * - ``TARGET_NAME``
     - Target
   * - A moon of Jupiter, such as ``J17 CALLIRRHOE``
     - The moon's name, ``CALLIRRHOE``
   * - A numbered minor planet with a name, such as ``136472 MAKEMAKE`` or
       ``50000 QUAOAR (2002 LM60)``
     - The name, ``MAKEMAKE`` or ``QUAOAR``
   * - ``ASTEROID 486958 (2014 MU69)``, and ``MU69`` in the FITS header
     - ``ARROKOTH``
   * - A name that begins with ``PLUTO``, such as ``PLUTO-CHARON``
     - ``PLUTO``
   * - ``CALIBRATION``, or no name at all
     - ``NONE``

Any other name is the target as it is, whether or not SPICE knows it. The light source
is the Sun, and the ``lightsource`` override replaces it.
