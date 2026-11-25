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
original keyword stays in ``header_dict``. A header keyword that has no PDS3 equivalent keeps its own name.

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
     - The two agree on every product.

If the chosen label lacks a value for an item, the other label's value is used.

Missing values
==============

A value missing from every label read is left out of ``dict`` and the subfields, and
nothing takes its place, neither None nor a placeholder. These count as missing:

- the placeholders ``UNK``, ``N/A``, ``NULL``, ``UNKNOWN``, ``NOT APPLICABLE``, ``---``,
  and an empty string, compared in any case and without surrounding blanks;
- the FITS header's fill value -1e32, which it writes for a value its pipeline cannot
  compute, such as the target geometry of an image that has no target;
- a FITS card without a value;
- a list or set whose elements all count as missing.

A list in the PDS3 label with some missing elements among others is kept as the label
gives it. The missing values stay in ``label_dict``, ``header_dict``, and ``headers``,
which hold the labels as they were read.

The dictionaries
================

Four subfields hold the metadata:

``dict``
    The merged metadata: every item of either label read, under the PDS3 name where the
    host converts it, and under its own name otherwise, without the values that count
    as missing. ``START_TIME`` is a :class:`datetime.datetime` in UTC, as the PDS3
    label gives it, whichever label it came from. Where both labels have a key that the
    host does not convert, the PDS3 label's value is kept. The blank, ``COMMENT``, and
    ``HISTORY`` cards of the FITS header are left out.

``label_dict``
    The PDS3 label as its parser gives it, as a dictionary, present only when the label
    was read. It includes the extra keys that the parser adds beside a keyword, such as
    ``START_TIME_fmt``, the time as the label writes it.

``header_dict``
    The primary FITS header as a dictionary, with each value as astropy reads it,
    present only when the FITS file was read.

``headers``
    The astropy header of every part of the FITS file, the primary header first,
    present only when the FITS file was read.

Subfields
=========

Besides ``dict``, ``label_dict``, ``header_dict``, ``headers``, and ``data``, a Snapshot
read from a file has these subfields. A subfield whose item is missing from every label
read is absent.

.. list-table::
   :header-rows: 1

   * - Subfield
     - Contents
   * - ``spacecraft_clock_start_count``, ``spacecraft_clock_stop_count``
     - ``SPACECRAFT_CLOCK_START_COUNT`` and ``SPACECRAFT_CLOCK_STOP_COUNT``, described
       below.
   * - ``error``, ``quality``
     - The error estimates and quality flags of a calibrated image, except in an
       astrometry read (see :doc:`user_guide_products`).
   * - ``spice_frame_name``, ``spice_frame_id``, ``spice_to_frame``
     - The SPICE frame of the camera, ``NH_LORRI`` (-98300), and the rotation from it
       to the camera frame, which reverses its *y* and *z* axes (see
       :doc:`user_guide_geometry`).
   * - ``ck_object_id``, ``ck_sclk_id``, ``spice_frame_is_frozen``
     - -98000, -98, and False (see :doc:`user_guide_geometry`).
   * - ``filepath``, ``basename``, ``image_url``, ``image_abspath``, ``label_url``,
       ``label_abspath``
     - The FITS file and the PDS3 label (see :doc:`user_guide_files`).

The spacecraft clock counts
===========================

``spacecraft_clock_start_count`` and ``spacecraft_clock_stop_count`` are the counts at
the start and the end of the exposure, written as the PDS3 label writes them,
``SSSSSSSSSS:TTTTT``, such as ``0030713585:49000``: ten digits of whole seconds and five
of the ticks of 1/50000 s. From the FITS header, the start count is ``SPCSCLK0``, and
the stop count is ``STOPMET``, a number of seconds, rounded to the nearest tick and
written in the same form. The two give the same counts on every product.
