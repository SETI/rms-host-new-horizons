======
Errors
======

The host reports each kind of failure with one exception class. The classes come from
OOPS, so code that reads images from several missions can catch them the same way.

Reading an image file
=====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` raises:

:exc:`oops.OopsSpiceError`
    The C kernels hold no pointing for the camera at the middle of the exposure. The
    message gives the time and the FITS file, and ends with the SPICE error code, such
    as ``SPICE(NOFRAMECONNECT)``. It is never raised for missing pointing when the
    caller supplies the pointing with the ``frame`` override, or takes it from the FITS
    header with the ``pointing`` option. Whatever the pointing, it is also raised when
    no label read gives a start time and the clock kernel cannot convert the clock
    count.

:exc:`~oops.host.HostFormatError`
    No label read gives the exposure duration, the lines of the image, or either a
    start time or a clock count; the image has neither 1024 nor 256 lines; the FITS file
    is not valid FITS; or a parsed label that was not read from a file is given for a
    normal read, so that no file path is known for its FITS file. The message gives the
    item, and the file if one is known. For an image passed in memory that was not read
    from a file, it says so instead.

:exc:`~oops.host.PdsHostError`
    None of the fast, strict, and loose methods of the parser can read the PDS3 label,
    or the label has no ``^IMAGE`` pointer to a FITS file.

:exc:`~oops.host.WrongHostError`
    The file, or the parsed label or FITS file passed in, is not of a New Horizons
    LORRI image.

:exc:`FileNotFoundError`
    The FITS file that a PDS3 label points to does not exist, in a read that needs it.

:exc:`~oops.OopsValueError`
    ``geom``, ``pointing``, or ``fov_type`` is not one of its values; ``select`` is not
    None; ``geom`` or ``pointing`` asks for the FITS header in an astrometry read of a
    PDS3 label; or ``pds3_options`` holds an option that the parser does not accept.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts.

:exc:`OSError`
    A file cannot be read.

Reading a PDS3 index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` raises:

:exc:`~oops.host.HostIndexError`
    Raised by the call itself: the index is a volume's own index, in a directory named
    ``INDEX``; ``supplement`` is None and ``row_dicts`` is not given; the index's file
    name holds no ``_index`` from which to find the supplement; or either label does not
    describe an index table with the columns ``VOLUME_ID`` and
    ``FILE_SPECIFICATION_NAME``. Raised when the iteration reaches the row concerned:
    a row that the table file lacks or that does not match its label; an index row whose
    image the supplement does not list, found once the whole supplement has been read;
    or a supplement row that lists an image a second time while its first row is still
    being kept. Raised after the last row: the supplement lists images that the index
    does not, or lists an image more times than the index.

:exc:`~oops.host.PdsHostError`
    The PDS3 label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The supplement's PDS3 label does not exist, which the call itself finds; or a table
    file does not exist, which the iteration finds when it starts.

:exc:`~oops.host.HostFormatError`
    A row's ``BINNING_MODE`` is missing or is not ``1X1`` or ``4X4``, or a row lacks the
    exposure duration, or lacks both a start time and a clock count. This error is
    raised when the iteration reaches that row.

:exc:`oops.OopsSpiceError`
    A row lacks a start time, and the clock kernel cannot convert its clock count. This
    error is raised when the iteration reaches that row.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts. A range in ``filter`` that
    cannot be compared with a row's value raises this error when the iteration reaches
    that row.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, ``select`` is not None, ``filter`` asks for a
    column that neither table has, or ``pds3_options`` holds an option that the parser
    does not accept.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`OSError`
    The PDS3 label of the index or of a named supplement does not exist or cannot be
    read.

The iteration ends at the first row that raises, after the Snapshots of the rows before
it. Missing pointing never stops the iteration: the Snapshot of a row without pointing
raises :exc:`oops.OopsSpiceError` only when its frame is evaluated.
