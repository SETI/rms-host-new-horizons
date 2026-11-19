====================
Reading a PDS3 Index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` makes one
:class:`~oops.observation.Snapshot` for each row of a PDS3 index. It reads only the
index tables and opens none of the images, which makes it much faster than opening a
volume's images one at a time:

.. code-block:: python

   from host_new_horizons import NewHorizonsLORRI

   snapshots = NewHorizonsLORRI.from_index(
       'metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl')

The index and its supplement
============================

The metadata tree of the Ring-Moon Systems Node holdings has, for each LORRI volume, a
PDS3 index, such as ``NHJULO_2001_index.lbl``, with one row per image, and a
supplemental index, ``NHJULO_2001_supplemental_index.lbl``, whose rows give further
values for the same images. The host needs both. Each row of the supplement is
merged into the matching row of the index; where both tables give a column, the
supplement's value is kept, because the supplement writes each value as the product's
label does. The index writes some of them differently: it gives each clock count with
its partition, as in ``1/0030598438:48850``, writes the product creation time as a date
alone, and shortens some target names.

The ``supplement`` argument chooses the supplemental index:

- By default, the host reads the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``, in upper or lower case. That file must
  exist, and an index whose file name holds no ``_index`` raises
  :exc:`~oops.host.HostIndexError`.
- A path or URL reads that file instead.
- None raises :exc:`~oops.host.HostIndexError`, unless ``row_dicts`` is given.

The two tables must list the same images, each once, but not necessarily in the same
order. The host matches the rows by their ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME``
and raises :exc:`~oops.host.HostIndexError` if an image is listed twice or appears in
only one of the tables. The Snapshots come in the order of the index.

Other arguments
===============

``label``
    The index's PDS3 label, if it has already been parsed.

``row_dicts``
    The rows themselves, one dictionary per row, if the tables have already been read.
    The host then reads neither table, and ``supplement`` and ``label`` are ignored.

``filter``
    A dictionary of column names and the values to keep. A row is kept when each column
    equals its value, is one of a set of values, or falls in a range given as a tuple of
    two values.

``select``
    Must be None. Each row is one image, so there is nothing to select, and any other
    value raises :exc:`~oops.OopsValueError`.

``fov_type``, ``return_all_planets``, ``pds3_options``, and the standard overrides
    The same as for :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`,
    applied to every row.

Any other keyword raises :exc:`~oops.OopsTypeError`, so that a misspelled argument is
not silently ignored.

What one row gives
==================

Each row gives a :class:`~oops.observation.Snapshot` with no ``data``. Its ``dict``
subfield holds the merged row, with the same treatment of missing values as a label (see
:doc:`user_guide_metadata`). The supplement declares a null constant for many columns,
and the host leaves out a value equal to it, with one exception: zero is the null
constant of the vectors and distances, and there a zero can be real. The target's
vectors and distances are left out only for an image that has no target geometry, whose
vector from the spacecraft to the target is zero. They are kept otherwise, so that an
image of the Sun keeps the Sun's distance of zero from itself. Any other vector is left
out only when all of its components are zero.

The subfields are those of an image read from its label (see
:doc:`user_guide_metadata`), except ``data``, ``error``, ``quality``, ``label_dict``,
``header_dict``, and ``headers``. The binning mode comes from the ``BINNING_MODE``
column. The pointing always comes from the C kernels unless the ``frame`` override is
given, and the position from the SPICE ephemerides. A row has no I/F calibration,
because it gives no RSOLAR.

The host does not check the C kernels for the pointing of an index row. If the kernels
hold no pointing for an image, its row still gives a Snapshot, and the error appears
only when the Snapshot's frame is evaluated.

File paths
==========

A row's ``FILE_SPECIFICATION_NAME`` is the path of the PDS3 label, but each Snapshot
refers to the FITS file, with ``.fit`` in place of the label's suffix, in the same case.
``filepath``, ``basename``, and ``image_url`` give it as a full path or URL, formed from
the location of the index:

- For an index in a volume's ``INDEX`` directory, the volume directories are those
  beside the volume that holds it.
- For an index in the ``metadata`` tree of the Ring-Moon Systems Node holdings, such as
  ``.../metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl``, they are in the
  matching ``volumes`` tree, ``.../volumes/NHxxLO_xxxx``.
- For an index anywhere else, they are in the index's own directory.

The row's ``VOLUME_ID`` and file specification complete the path. ``abspath`` is present
only when the FITS file is local or already in the local cache.

Reading an index through oops.Host
==================================

:meth:`oops.host.Host.from_index` cannot read a LORRI index. It requires the index and
the supplement to agree in every column they share, but the two tables write the clock
counts, the product creation time, and some target names differently, so the rows fail
that test. Call :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index`
instead.
