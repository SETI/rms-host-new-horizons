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
    is not valid FITS; or a parsed label or FITS file passed in was not read from a
    file. The message gives the item and the file.

:exc:`~oops.host.PdsHostError`
    None of the fast, strict, and loose methods of the parser can read the PDS3 label,
    or the label has no ``^IMAGE`` pointer to a FITS file.

:exc:`~oops.host.WrongHostError`
    The file, or the parsed label or FITS file passed in, is not of a New Horizons
    LORRI image.

:exc:`FileNotFoundError`
    The FITS file that a PDS3 label points to does not exist, in a read that needs it.

:exc:`~oops.OopsValueError`
    ``geom``, ``pointing``, or ``fov_type`` is not one of its values; ``geom`` or
    ``pointing`` asks for the FITS header in an astrometry read of a PDS3 label; or
    ``pds3_options`` holds an option that the parser does not accept.

:exc:`OSError`
    A file cannot be read.

Reading a PDS3 index
====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` raises:

:exc:`~oops.host.HostIndexError`
    ``supplement`` is None and ``row_dicts`` is not given; the supplement and the index
    do not list the same images, each once; or the index's file name holds no
    ``_index`` from which to find the supplement.

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
    A keyword argument is not one that the method accepts.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, or ``select`` is not None.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`OSError`
    The index or a named supplement does not exist or cannot be read.

The method does not check the pointing, and never raises :exc:`oops.OopsSpiceError` for
missing pointing.
