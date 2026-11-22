======
Errors
======

The host reports each kind of failure with one exception class. The classes come from
OOPS, so code that reads images from several missions can catch them the same way.

Reading an image file
=====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` raises:

:exc:`oops.OopsSpiceError`
    The C kernels lack the pointing of the camera at the time when the observation will
    read them, the middle of the exposure as a ``cadence`` override or a ``timeshift`` of
    the start and stop times sets it, or the start, middle, or end that a ``tracker``
    override selects (see :doc:`user_guide_geometry`). The message gives
    the time and the FITS file, and ends with the SPICE error code, such as
    ``SPICE(NOFRAMECONNECT)``. It is never raised for missing pointing when the
    caller supplies the pointing with the ``frame`` override, or takes it from the FITS
    header with the ``pointing`` option. Whatever the pointing, it is also raised when
    the labels read lack a start time and the clock kernel cannot convert the clock
    count.

:exc:`~oops.host.HostFormatError`
    The labels read lack the exposure duration, the lines of the image, or both a start
    time and a clock count; the image has neither 1024 nor 256 lines; the FITS file
    is not valid FITS; or a parsed label that was not read from a file is given for a
    normal read, so that no file path is known for its FITS file. The message gives the
    item, and the file if one is known. For an image passed in memory that was not read
    from a file, it says so instead.

:exc:`~oops.host.PdsHostError`
    None of the fast, strict, and loose methods of the parser can read the PDS3 label,
    the label has no ``^IMAGE`` pointer to a FITS file, or, in a read that opens the
    FITS file, the name it points to does not end in ``.fit`` or ``.fits``.

:exc:`~oops.host.WrongHostError`
    The file, or the parsed label or FITS file passed in, is not of a New Horizons
    LORRI image.

:exc:`FileNotFoundError`
    The FITS file that a PDS3 label points to does not exist, in a read that needs it.

:exc:`~oops.OopsValueError`
    ``geom``, ``pointing``, or ``fov_type`` is not one of its values; ``select`` is not
    None; ``geom`` or ``pointing`` asks for the FITS header in an astrometry read of a
    PDS3 label; or ``pds3_options`` holds an option that the parser does not accept. The
    options are checked only when the method parses a PDS3 label file, so they are not
    checked for a FITS file or for a label passed in already parsed.

:exc:`TypeError`
    A keyword argument is not one that the method accepts.

:exc:`OSError`
    A file cannot be read.

Reading a PDS3 index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` raises the errors
listed below. When the rows are given in ``row_dicts``, the call leaves the tables and
their labels unread, so the errors that come from reading a label or a table are never
raised: the labels need not exist or be readable, the index's file name need not
contain ``_index``, the index is not checked to be a New Horizons LORRI index, and
``pds3_options`` is not checked. The checks of ``select``, ``fov_type``, and ``filter``
are still made, as their entries say, and so are those of each row.

:exc:`~oops.host.HostIndexError`
    Raised by the call itself: ``supplement`` is None or an empty string, and
    ``row_dicts`` is not given; ``supplement`` is left at its default, and the index's
    file name lacks the ``_index`` from which the host finds the supplement; or either
    label does not describe an index table with the columns ``VOLUME_ID`` and
    ``FILE_SPECIFICATION_NAME``. Raised during the iteration: after the last row of a
    table file that holds fewer rows than its label says; for a row longer than its
    label's row length, at the first of the 100 rows read with it, so those 100 rows
    are not built; for an index row whose image the supplement does not list, once the
    whole supplement has been read; or for a supplement row that lists an image a second
    time while its first row is still being kept. Raised after the last row: a table
    file holds more rows than its label says, or the supplement lists images that the
    index does not, or lists an image more times than the index.

:exc:`~oops.host.PdsHostError`
    The PDS3 label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The supplement's PDS3 label does not exist, which the call itself finds; or a table
    file does not exist, which the iteration finds when it starts.

:exc:`~oops.host.HostFormatError`
    A row lacks ``VOLUME_ID`` or ``FILE_SPECIFICATION_NAME``, which only a row given in
    ``row_dicts`` or a null value can; a row's ``BINNING_MODE`` is missing or is not
    ``1X1`` or ``4X4``; or a row lacks the exposure duration, or lacks both a start time
    and a clock count. This error is raised when the iteration reaches that row.

:exc:`oops.OopsSpiceError`
    A row lacks a start time, and the clock kernel cannot convert its clock count. This
    error is raised when the iteration reaches that row.

:exc:`~oops.OopsTypeError`
    A range in ``filter`` cannot be compared with a row's value. This error is raised
    when the iteration reaches that row.

:exc:`TypeError`
    A keyword argument is not one that the method accepts, which the call itself finds.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, ``select`` is not None, ``filter`` asks for a
    column that neither table has, or ``pds3_options`` holds an option that the parser
    does not accept. The options are checked only when the call parses the index's
    label, so they are not checked when ``label`` is given. With ``row_dicts`` given as
    an iterable other than a list or another sequence, a ``filter`` column that the
    first row lacks raises this error when the iteration reads that row.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`OSError`
    The PDS3 label of the index or of a named supplement does not exist or cannot be
    read.

The iteration ends at the first row that raises, after the Snapshots of the rows before
it. Missing pointing never stops the iteration: the Snapshot of a row without pointing
raises :exc:`oops.OopsSpiceError` only when its frame is evaluated.
