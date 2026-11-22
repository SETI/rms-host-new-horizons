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
    When the call is made, any of these:

    - ``supplement`` is left at its default, and the index's file name lacks the
      ``_index`` from which the host finds the supplement. This is not checked when
      ``row_dicts`` is given.
    - ``supplement`` is None or an empty string, and ``row_dicts`` is not given.
    - A table cannot be read as an index, because its label lacks the table it points
      to.
    - A table lacks the column ``VOLUME_ID`` or ``FILE_SPECIFICATION_NAME``.

    During the iteration:

    - At an index row whose image the supplement does not list, once the whole
      supplement has been read. An image that the index lists more times than the
      supplement is found this way too.
    - At a supplement row that lists an image a second time while its first row is
      still being kept.
    - After the last row of a table that holds fewer rows than its label says. Every
      row that the table holds is given first.
    - At a row that is longer than its label says. The rows are read 100 at a time, and
      none of the 100 read with that row is given.

    After the last row: a table holds more rows than its label says, or the supplement
    lists images that the index does not, or lists an image more times than the index.

    The tables are checked only for these errors. A value that cannot be read as the
    number its column describes gives a warning and is left out of the row, as a missing
    value is, so a row whose exposure cannot be read raises
    :exc:`~oops.host.HostFormatError`. A row that is shorter than its label says is read
    as it is.

:exc:`FileNotFoundError`
    The PDS3 label of the supplement does not exist, when the call is made; or a table
    file does not exist, when the iteration starts.

:exc:`~oops.host.HostFormatError`
    A row lacks ``VOLUME_ID`` or ``FILE_SPECIFICATION_NAME``, which only a row given in
    ``row_dicts`` or a null value can; a row's ``BINNING_MODE`` is missing or is not
    ``1X1`` or ``4X4``; or a row lacks the exposure duration, or lacks both a start time
    and a clock count. The message gives the item and the FITS file, or the item alone
    for a row without ``VOLUME_ID`` or ``FILE_SPECIFICATION_NAME``, whose file is not
    known. This error is raised when the iteration reaches that row.

:exc:`oops.OopsSpiceError`
    A row lacks a start time, and the clock kernel cannot convert its clock count. This
    error is raised when the iteration reaches that row.

:exc:`TypeError`
    A keyword argument is not one that the method accepts, when the call is made. This
    is Python's own error for an unexpected keyword argument.

:exc:`~oops.OopsTypeError`
    A range in ``filter`` cannot be compared with a row's value. This error is raised
    when the iteration reaches that row.

:exc:`~oops.OopsValueError`
    ``fov_type`` is not one of its values, ``select`` is not None, ``filter`` gives a
    column that neither table has, or ``pds3_options`` names a method the parser does
    not have or an option it does not take. The call checks the options when it is
    made, unless ``row_dicts`` is given. When ``row_dicts`` is given, ``filter`` is
    checked against the first row: when the call is made if the rows are a list or
    another sequence, and otherwise when the iteration reads that row.

:exc:`~oops.host.PdsHostError`
    The PDS3 label of the index or of the supplement cannot be parsed.

:exc:`~oops.host.WrongHostError`
    The file is not a New Horizons LORRI index.

:exc:`OSError`
    The PDS3 label of the index or of a named supplement does not exist or cannot be
    read.

The iteration ends at the first row that raises, after the Snapshots of the rows before
it. The ``filter`` argument (see :doc:`user_guide_index`) can leave out a row whose
Snapshot would raise, but it cannot avoid an error in reading or pairing the rows of the
tables, because the rows are filtered only after they have been read and paired.
Missing pointing never stops the iteration: the Snapshot of a row without pointing
raises :exc:`oops.OopsSpiceError` only when its frame is evaluated.
