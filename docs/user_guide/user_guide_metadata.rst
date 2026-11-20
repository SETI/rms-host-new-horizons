=======================
Metadata and Subfields
=======================

The two labels
==============

A LORRI product describes itself twice. The *PDS3 label* is the detached ``.lbl`` file
beside the image. The *FITS header* is the label inside the FITS file itself: the
primary header, which describes the image, and the headers of the extensions that
follow it. The two were written at different times and do not always agree, and each
holds values that the other lacks.

The host takes its metadata from whichever of the two it read (see
:doc:`user_guide_files`): from the FITS header alone when given the FITS file, from the
PDS3 label alone in an astrometry read of the label, and from both when given the label
otherwise.

Merging the two labels
======================

The FITS header uses its own keywords. Where a header value has a PDS3 equivalent, the
host converts it to the PDS3 name and to the form that the PDS3 label uses:

- ``SPCSCET``, the middle of the exposure, less half of ``EXPTIME``, becomes
  ``START_TIME``.
- ``EXPTIME`` becomes ``EXPOSURE_DURATION``.
- ``SPCCBTNM`` becomes ``TARGET_NAME``.
- ``SPCSCLK0`` becomes ``SPACECRAFT_CLOCK_START_COUNT``.
- ``STOPMET``, a count of seconds, becomes ``SPACECRAFT_CLOCK_STOP_COUNT``, written as
  described below.
- ``NAXIS2`` gives the number of lines of the image, which the PDS3 label gives as the
  ``LINES`` of its ``IMAGE`` object.

A converted value appears in the merged metadata only under its PDS3 name, and the
original keyword stays in ``header_dict``. ``NAXIS2`` is the exception: it keeps its
name, because the PDS3 label gives the lines only inside its ``IMAGE`` object, which
keeps its own name. A header keyword that has no PDS3 equivalent keeps its own name.

Where only one label gives an item, its value is used. Where both give it and agree,
that value is used. Where both give it and disagree, the item always takes its value
from the same label, for every image:

.. list-table::
   :header-rows: 1

   * - Items
     - Source
     - Why
   * - ``START_TIME``
     - FITS header
     - The label gives the start time rounded to the millisecond. The header's start
       agrees to 0.0001 ms with the time of the label's own clock count. On an image
       taken while the spacecraft was spinning, the rounding can move the pointing by
       as much as 48.7 pixels.
   * - ``EXPOSURE_DURATION``
     - FITS header
     - The two differ only in the last binary digit, and the header's value keeps the
       start plus half of the exposure equal to ``SPCSCET``.
   * - The lines of the image
     - FITS header
     - The header describes the array that ``data`` holds, which the field of view must
       match.
   * - ``TARGET_NAME``
     - PDS3 label
     - The header's ``SPCCBTNM`` is the data pipeline's guess at the target, such as
       ``FJUPITER`` for Jupiter, ``ANSA`` for Jupiter's rings, or a pointing such as
       ``RADEC=270.66,-14.63=P4`` for Kerberos.
   * - ``SPACECRAFT_CLOCK_START_COUNT`` and ``SPACECRAFT_CLOCK_STOP_COUNT``
     - PDS3 label
     - The two agreed on every product checked.

If the chosen label gives no value for an item, the other label's value is used.

Missing values
==============

A value that no label read gives is left out: there is no ``dict`` key and no subfield
for it, and never a None or a placeholder in its place. These count as no value:

- the PDS3 null constants ``UNK``, ``N/A``, and ``NULL``, and the values ``---``,
  ``UNKNOWN``, and an empty string, which the FITS header writes for a body, a request,
  or a file it does not give, all compared in any case;
- the FITS header's fill value -1e32, which it writes for a value its pipeline cannot
  compute, such as the target geometry of an image that has no target;
- a list whose elements all count as no value.

A list in the PDS3 label with some missing elements among others is a vector, such as
``SC_TARGET_POSITION_VECTOR``, in which each position has a meaning. In ``dict`` it
keeps its length, with None for each missing element.

The dictionaries
================

Four subfields hold the metadata as dictionaries. All four can be written with
:func:`json.dump`: times are ISO 8601 text, sets are lists, and an item has the same
type whichever label it came from.

``dict``
    The merged metadata: every item of either label read, under the PDS3 name where the
    host converts it, and under its own name otherwise. ``START_TIME`` is written to the
    millisecond, as the label writes it. Where both labels have a key that the host
    does not convert, the PDS3 label's value is kept.

``label_dict``
    The PDS3 label as its parser reads it, present only when the label was read. Times
    are the label's own text, and sets are lists in the label's order. The ``END``
    keyword and the extra keys that the parser adds beside each keyword, such as
    ``START_TIME_fmt`` and ``^IMAGE_offset``, are left out.

``header_dict``
    Every card of the primary FITS header that has a value, present only when the FITS
    file was read. The blank, ``COMMENT``, and ``HISTORY`` cards are left out. A keyword
    that occurs more than once is keyed by its name and its occurrence counted from
    zero: ``REFFILE[0]``, ``REFFILE[1]``, and so on.

``headers``
    A list with the header of every part of the FITS file, the primary header first,
    present only when the FITS file was read. Each is keyed as ``header_dict`` is, and
    the text of its blank, ``COMMENT``, and ``HISTORY`` cards is kept as a list under
    each of those keywords.

Subfields
=========

Besides ``dict``, ``label_dict``, ``header_dict``, ``headers``, and ``data``, every
Snapshot has these subfields. One whose item no label gives is absent.

.. list-table::
   :header-rows: 1

   * - Subfield
     - Contents
   * - ``host``, ``instrument``
     - ``"New Horizons"`` and ``"LORRI"``.
   * - ``binning_mode``
     - ``"1X1"`` or ``"4X4"`` (see :doc:`user_guide_products`).
   * - ``label_target``
     - ``TARGET_NAME`` as the label spells it, or, without the label, as the FITS
       header spells it.
   * - ``spacecraft_clock_start_count``, ``spacecraft_clock_stop_count``
     - ``SPACECRAFT_CLOCK_START_COUNT`` and ``SPACECRAFT_CLOCK_STOP_COUNT``, described
       below.
   * - ``time_from_sclk``
     - True if no label read gives a start time, so that the time comes from the clock
       count (see :doc:`user_guide_geometry`); False otherwise.
   * - ``error``, ``quality``
     - The error estimates and quality flags of a calibrated image, except in an
       astrometry read (see :doc:`user_guide_products`).
   * - ``spice_kernels``
     - The SPICE kernels used for the observation.
   * - ``spice_frame_name``, ``spice_frame_id``, ``spice_to_frame``
     - The SPICE frame of the camera, ``NH_LORRI`` (-98300), and the rotation from it
       to the camera frame, which reverses its *y* and *z* axes (see
       :doc:`user_guide_geometry`).
   * - ``filepath``, ``basename``, ``abspath``, ``image_url``
     - The FITS file, absent for a parsed label or FITS file that was not read from a
       file (see :doc:`user_guide_files`).

The spacecraft clock counts
===========================

``spacecraft_clock_start_count`` and ``spacecraft_clock_stop_count`` are the counts at
the start and the end of the exposure, written as the PDS3 label writes them,
``SSSSSSSSSS:TTTTT``, such as ``0030713585:49000``: ten digits of whole seconds and five
of the ticks of 1/50000 s. From the FITS header, the start count is ``SPCSCLK0``, and
the stop count is ``STOPMET``, a number of seconds, rounded to the nearest tick and
written in the same form.
