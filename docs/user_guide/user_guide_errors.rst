======
Errors
======

The host reports each kind of failure with one exception class. The classes come from
OOPS, so code that reads images from several missions can catch them the same way.

Reading an image file
=====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` raises:

:exc:`~oops.host.HostFormatError`
    The label or the FITS file exists but cannot be read, as when access to it is denied,
    or cannot be interpreted: the FITS file is not valid FITS or ends early; its primary
    array is not one image band; the array does not have the lines and samples that the
    label's ``IMAGE`` object gives; or the files lack an item that the Snapshot needs,
    such as the exposure duration, the start time, or lines of an image size that LORRI
    takes. The message names the file that was passed and the failure.

:exc:`~oops.host.PdsHostError`
    The parser cannot read the PDS3 label.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI label or FITS file.

:exc:`FileNotFoundError`
    The file does not exist, or the FITS file that a PDS3 label names does not exist, in
    a read that needs it.

:exc:`~oops.OopsValueError`
    ``geom``, ``pointing``, or ``fov_type`` is not one of its values; ``geom`` or
    ``pointing`` asks for the FITS header in an astrometry read of a PDS3 label; or
    ``pds3_options`` holds an option that the parser does not accept.

:exc:`~oops.OopsIndexError`
    ``select`` is an index other than 0 or -1.

:exc:`TypeError`
    A keyword argument is not one that the method accepts.

Missing pointing raises no error here: the Snapshot's frame raises the SPICE error when
it is used (see :doc:`user_guide_geometry`).

Reading a PDS3 index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` raises these when it
is called:

:exc:`~oops.host.HostIndexError`
    ``supplement`` is None and ``row_dicts`` is not given; the index's file name holds no
    ``_index`` from which to find the supplement; or a label does not describe a table.

:exc:`~oops.host.PdsHostError`
    The PDS3 label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The PDS3 label of the index or of the supplement does not exist.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, or ``pds3_options`` holds an option that the
    parser does not accept.

:exc:`TypeError`
    A keyword argument is not one that the method accepts.

The iteration raises these when it reaches the row concerned, after the Snapshots of
the rows before it, and ends there:

:exc:`~oops.host.HostIndexError`
    The two tables do not list the same images, or a table file holds a different number
    of rows than its label says.

:exc:`FileNotFoundError`
    A table file does not exist, which the iteration finds when it starts.

:exc:`KeyError`
    A row lacks a column that the Snapshot needs, or a column that ``filter`` names.

The rows are given as the tables give them, so a value that the Snapshot cannot use
raises the error of the code that reads it. Missing pointing never stops the iteration:
the Snapshot of a row without pointing raises the SPICE error only when its frame is
used.
