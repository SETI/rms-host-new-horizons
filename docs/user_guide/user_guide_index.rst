====================
Reading a PDS3 Index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` makes one
:class:`~oops.observation.Snapshot` for each row of a PDS3 index. It reads only the
index tables and leaves the images unopened, which makes it much faster than opening a
volume's images one at a time.

It returns an iterator, which gives the Snapshots row by row, in the order of the index.
The call itself checks the arguments and reads the PDS3 labels of the index and its
supplement, so a wrong argument or a label that cannot be read raises an error at once.
The rows are read as the iteration proceeds, 100 at a time, and each row's Snapshot is
built when the iteration reaches that row, so neither table is held in memory whole. An
error in a row is raised when the iteration reaches it, after the rows before it have
been given, and the iteration ends there. A caller that needs ``len()`` or indexing
passes the iterator to ``list()``:

.. code-block:: python

   from host_new_horizons import NewHorizonsLORRI

   index = 'metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl'
   for snapshot in NewHorizonsLORRI.from_index(index):
       print(snapshot.basename, snapshot.time)

   snapshots = list(NewHorizonsLORRI.from_index(index))

The index and its supplement
============================

The metadata tree of the Ring-Moon Systems Node holdings has, for each LORRI volume, a
PDS3 index, such as ``NHJULO_2001_index.lbl``, with one row per image, and a
supplemental index, ``NHJULO_2001_supplemental_index.lbl``, whose rows give further
values for the same images. :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index`
reads such a metadata index and its supplement, and it needs both: the binning mode and
the exposure duration are only in the supplement.

The ``supplement`` argument chooses the supplemental index:

- By default, the host reads the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``, in upper or lower case. That file must
  exist, and an index whose file name lacks ``_index`` raises
  :exc:`~oops.host.HostIndexError`.
- A path or URL reads that file instead, from whatever directory it names.
- None raises :exc:`~oops.host.HostIndexError`, unless ``row_dicts`` is given.

Each index row is merged with the supplement row that has the same ``VOLUME_ID`` and
``FILE_SPECIFICATION_NAME``, wherever that row is in the supplement, because an index
can list an image out of order. The cumulative indexes list a few images twice under one
file name, once for each of two volumes, and ``VOLUME_ID`` tells them apart. No other
column of the two rows is compared, because the two tables write several shared values
in different forms. Where both tables give a column, the merged row has the index's
value. For example, the index gives each clock count with its partition, as in
``1/0030598438:48850``, where the supplement gives ``0030598438:48850``, so the merged
row has the index's count.

The host raises :exc:`~oops.host.HostIndexError`, naming both tables, when the two do
not list the same images, at the row where it finds the difference.

Other arguments
===============

``label``
    The index's PDS3 label, if it has already been parsed.

``row_dicts``
    The rows themselves, one dictionary per row, each already merged with its
    supplement row, as a list or any other iterable. The host then reads neither table,
    and ``supplement`` and ``label`` are ignored.

``filter``
    A dictionary of column names and the values to keep. A row is kept when each column
    equals its value, is one of a set of values, or falls in a range given as a tuple of
    two values.

``pds3_options``
    Options for the PDS3 label parser. The index's label is parsed with all of them, and
    the supplement's label with the method alone. The tables' labels are parsed as they
    are written, with the fast method unless the options give another.

``fov_type``, ``select``, and the standard overrides
    The same as for :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`,
    applied to every row.

Any other keyword raises :exc:`TypeError`, so that a misspelled argument is not silently
ignored.

What one row gives
==================

Each row gives a :class:`~oops.observation.Snapshot` without ``data``, built from the
row's values as the tables give them:

- The exposure starts at ``START_TIME`` and lasts ``EXPOSURE_DURATION``. A row whose
  exposure is zero gives a Snapshot whose exposure is zero.
- The field of view is that of ``BINNING_MODE``, ``1X1`` or ``4X4``, read in either
  case; the supplement writes it in lower case.
- The target comes from ``TARGET_NAME`` as for an image (see :doc:`user_guide_geometry`).
- The pointing comes from the C kernels and the position from the SPICE ephemerides, and
  the ``frame`` and ``path`` overrides replace them. As for an image, a row whose
  pointing the C kernels do not hold still gives a Snapshot, and its frame raises the
  SPICE error when it is used, so one image without pointing does not stop the rest of
  the index.
- There is no calibration, because the tables do not give RSOLAR, unless the
  ``calibrations`` override gives one.

Its subfields are:

``dict``
    The merged row exactly as the tables give it, null values included, with the
    ``_mask`` entry that the table reader gives beside each column, which is True where
    the value is a declared null or could not be read.

``filepath`` and ``basename``
    The row's ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME`` joined by a slash, as the
    index gives them, such as
    ``NHJULO_2001/data/20070110_003070/lor_0030702839_0x630_sci.lbl``, and the last part
    of that path. The path names the product's label, relative to the volumes directory
    of the archive; a caller that wants a file joins it with its own copy or URL of that
    directory. Reading an index opens no image or label file, so a row has no
    ``image_url``, ``image_abspath``, ``label_url``, or ``label_abspath``.

``spacecraft_clock_start_count`` and ``spacecraft_clock_stop_count``
    The clock counts as the row gives them, with the index's partition.

``spice_frame_name``, ``spice_frame_id``, ``spice_to_frame``, ``ck_object_id``, ``ck_sclk_id``, and ``spice_frame_is_frozen``
    The same as for an image (see :doc:`user_guide_metadata`).

Reading an index through oops.Host
==================================

:meth:`oops.host.Host.from_index` recognizes a LORRI index from its label and passes it,
with every other argument, to
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index`, which reads the index
and its supplement itself. The two calls give the same Snapshots.
