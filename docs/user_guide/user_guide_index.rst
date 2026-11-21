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
built when the iteration reaches that row, so neither table is held in memory whole. The
memory used stays the same however long the index is when OOPS releases each row's
camera frame along with its Snapshot; an OOPS that keeps those frames adds a few
kilobytes for each row. An error in a row is raised when the iteration reaches it, after
the rows before it have been given. The one exception is a row longer than the label
says, which is found as the rows are read: its error comes at the first of the 100 rows
read with it. The iteration ends at the row that raises, and the rows after it are not built.
The ``filter`` argument, described below, can leave such a row out. To build
every row at once, pass the iterator to ``list()``:

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
values for the same images. ``from_index()`` reads such a metadata index and its
supplement, and the host needs both. Each row of the supplement is merged into the
matching row of the index; where both tables give a column, the
supplement's value is kept, because it is closer to the product's label. The index
gives each clock count with its partition, as in ``1/0030598438:48850``, rounds the stop
count differently, writes the product creation time as a date alone, and shortens some
target names, such as ``50000 QUAOAR (2002 LM)``.

The supplement does not match the product labels everywhere either. It gives the
exposure duration rounded to the millisecond, such as 0.501 s where the label gives
0.5006 s. It leaves out the number and colon at the start of some sequence IDs. And it
gives the target names of five moons of Jupiter, Io, Europa, Himalia, Elara, and
Callirrhoe, without the ``J1``, ``J2``, ``J6``, ``J7``, or ``J17`` that their labels put
in front, as in ``IO`` for ``J1 IO``. Every label of each of these moons writes its name
the same way, so until the supplement is corrected, the host writes these five names as
the labels do. Ganymede and Callisto have no number in front in any label, and the
supplement agrees.

The ``supplement`` argument chooses the supplemental index:

- By default, the host reads the file beside the index whose name has
  ``_supplemental_index`` in place of ``_index``, in upper or lower case. That file must
  exist, and an index whose file name lacks ``_index`` raises
  :exc:`~oops.host.HostIndexError`.
- A path or URL reads that file instead.
- None raises :exc:`~oops.host.HostIndexError`, unless ``row_dicts`` is given.

The two tables must list the same images, but not necessarily in the same order. Until
the index is corrected, the index of ``NHJULO_1001``, and so the cumulative
``NHxxLO_1999``, lists one image out of place, last, while the supplement lists every
image in file order. So each index row is merged with the supplement row
that has the same ``VOLUME_ID`` and ``FILE_SPECIFICATION_NAME``. The supplement is read
only as far ahead as the row that matches the current index row, and the rows read
ahead are kept until their own index rows come. The Snapshots come in the order of the
index.

The host raises :exc:`~oops.host.HostIndexError`, naming both tables, when the two do
not list the same images:

- When an index row's image is not in the supplement. The host finds this only once it
  has read the whole supplement, so the error comes after the Snapshots of the rows
  before it. An image that the index lists twice is found this way too.
- When the supplement lists an image a second time while its first row is still being
  kept.
- After the last index row, when the supplement lists images that the index does not, or
  lists an image more times than the index.

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
    two values. A row without one of the columns is left out. A column that neither
    table has raises :exc:`~oops.OopsValueError` when the method is called. With
    ``row_dicts``, the columns are checked against the first row instead: when the
    method is called if the rows are a list or another sequence, and when the iteration
    reads the first row otherwise. A range that cannot be compared with a row's value
    raises :exc:`~oops.OopsTypeError` when the iteration reaches that row.

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

Each row gives a :class:`~oops.observation.Snapshot` without ``data``. Its ``dict``
subfield holds the merged row, with the same treatment of missing values as a label (see
:doc:`user_guide_metadata`). The supplement declares a null constant for many columns,
and the host leaves out a value equal to it, with one exception: zero is the null
constant of the vectors and distances, and there a zero can be real. The target's
vectors and distances are left out only for an image that has no target geometry, whose
vector from the spacecraft to the target is zero. They are kept otherwise, so that an
image of the Sun keeps the Sun's distance of zero from itself. Any other vector is left
out only when all of its components are zero. These rules apply only until the
supplement is corrected: once it declares a null constant that is never a real value,
no zero is equal to it, and every zero is kept.

The subfields are those of an image read from its label (see
:doc:`user_guide_metadata`), except ``data``, ``error``, ``quality``, ``label_dict``,
``header_dict``, and ``headers``. The binning mode comes from the ``BINNING_MODE``
column, and the exposure from the supplement's ``EXPOSURE_DURATION``, which is rounded to
the millisecond. The target comes from ``TARGET_NAME`` as for a label (see
:doc:`user_guide_geometry`). For the five moons of Jupiter that the supplement lists
without their number, ``TARGET_NAME`` in ``dict`` and ``label_target`` give the label's
name, such as ``J1 IO``, as for an image read from its label. By default the pointing
comes from the C kernels and the position from the SPICE ephemerides, and the ``frame``
and ``path`` overrides replace them. A row's list of calibrations is empty, because the
tables lack RSOLAR, but the ``calibrations`` override supplies one.

The host looks up the pointing of each row at the middle of its exposure, but a row
whose pointing the C kernels do not hold still gives a Snapshot. Its frame raises
:exc:`oops.OopsSpiceError` when it is evaluated, so one image without pointing does not
stop the rest of the index.

File paths
==========

A row's ``FILE_SPECIFICATION_NAME`` is the path of the PDS3 label, but each Snapshot
refers to the FITS file, with ``.fit`` in place of the label's suffix, in the same case.
``filepath``, ``basename``, and ``image_url`` give it as a full path or URL, formed from
the location of the index:

- For an index in the ``metadata`` tree of the Ring-Moon Systems Node holdings, such as
  ``.../metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl``, the volume
  directories are in the matching ``volumes`` tree, ``.../volumes/NHxxLO_xxxx``.
- For an index anywhere else, they are in the index's own directory.

The row's ``VOLUME_ID`` and file specification complete the path. The Snapshots lack
``abspath``, because reading an index does not read the FITS files, even local ones.
Building the Snapshots leaves the file system as it was: it does not create directories.

Reading an index through oops.Host
==================================

:meth:`oops.host.Host.from_index` recognizes a LORRI index from its label and passes it,
with every other argument, to
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index`, which reads the index
and its supplement itself. The two calls give the same Snapshots and raise the same
errors.
