=============================
Time, Pointing, and Targets
=============================

Time
====

``START_TIME`` is the start of the exposure. The exposure lasts ``EXPOSURE_DURATION``,
and the Snapshot's default cadence covers exactly that interval. From the FITS header,
the start is ``SPCSCET``, the middle of the exposure in seconds TDB, less half of
``EXPTIME``. The PDS3 label gives the same start in UTC, rounded to the millisecond, and
the header's value is used when both are read (see :doc:`user_guide_metadata`).

If no label read gives a start time, the host converts
``SPACECRAFT_CLOCK_START_COUNT`` to a time with the SPICE clock kernel, and the
``time_from_sclk`` subfield is True. An image cannot be read if it has neither a start
time nor a clock count (see :doc:`user_guide_errors`).

Pointing
========

By default the camera's pointing comes from the C kernels, through the SPICE frame
``NH_LORRI``. The instrument kernel puts the camera's boresight along the frame's -*z*
axis, while OOPS expects it along +*z*. The host therefore registers the SPICE frame as
``NH_LORRI_FLIPPED`` and defines the camera frame ``NH_LORRI`` from it, with its *y*
and *z* axes reversed. The ``spice_to_frame`` subfield holds that rotation.

When :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` builds the frame
from the C kernels, it checks at once that they hold the pointing at the middle of the
exposure, after any ``cadence`` or ``timeshift`` override. If they do not, it raises
:exc:`oops.OopsSpiceError` (see :doc:`user_guide_errors`).
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` looks the pointing up
in the same way, and a row without pointing still gives a Snapshot, whose frame raises
:exc:`oops.OopsSpiceError` when it is evaluated.

The pointing can also come from the FITS header, with the ``pointing`` option (see
:doc:`user_guide_files`), and the spacecraft's position likewise with the ``geom``
option. A caller that has better pointing can pass it with the ``frame`` override, for
example as a :class:`~oops.frame.Cmatrix`. With either, the host does not consult the
C kernels, so a missing C kernel never raises an error.

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
is always the Sun.

The ``label_target`` subfield keeps ``TARGET_NAME`` as the label spells it.
