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
    item and the file.

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
    The index is a volume's own index, in a directory named ``INDEX``; ``supplement`` is
    None and ``row_dicts`` is not given; an image appears in only one of the index and
    the supplement; the two tables list the images in different orders and one of them
    lists an image twice; the index's file name holds no ``_index`` from which to find
    the supplement; or either table cannot be read as an index table.

:exc:`~oops.host.PdsHostError`
    The PDS3 label of the index or of the supplement cannot be parsed.

:exc:`FileNotFoundError`
    The supplemental index does not exist.

:exc:`~oops.host.HostFormatError`
    A row's ``BINNING_MODE`` is missing or is not ``1X1`` or ``4X4``, or a row gives no
    exposure duration, or neither a start time nor a clock count. This error is raised
    when the iteration reaches that row.

:exc:`oops.OopsSpiceError`
    A row gives no start time, and the clock kernel cannot convert its clock count. This
    error is raised when the iteration reaches that row.

:exc:`~oops.OopsTypeError`
    A keyword argument is not one that the method accepts, or a range in ``filter``
    cannot be compared with the values of its column.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, ``select`` is not None, ``filter`` asks for a
    column that no row has, or ``pds3_options`` holds an option that the parser does not
    accept.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`OSError`
    The index or a named supplement does not exist or cannot be read.

Every error except those of a row is raised by the call itself, before the iteration
starts. The iteration ends at the first row that raises. Missing pointing never stops
the iteration: the Snapshot of a row without pointing raises
:exc:`oops.OopsSpiceError` only when its frame is evaluated.
