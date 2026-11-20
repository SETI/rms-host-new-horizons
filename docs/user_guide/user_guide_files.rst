=====================
Reading an Image File
=====================

:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` makes one
:class:`~oops.observation.Snapshot` from one LORRI image. It accepts either file of the
product, the FITS file or its PDS3 label:

.. code-block:: python

   from host_new_horizons import NewHorizonsLORRI

   snapshot = NewHorizonsLORRI.from_file(
       'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit')
   snapshot = NewHorizonsLORRI.from_file(
       'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.lbl')

What is read
============

A path that ends in ``.fit`` or ``.fits``, in any case, is read as a FITS file. A path
that ends in ``.lbl``, in any case, is read as a PDS3 label, and so is any other file
whose first bytes hold ``PDS_VERSION_ID``. Any other file raises
:exc:`~oops.host.WrongHostError`. What the host reads depends on the file it is given
and on the ``astrometry`` option. An astrometry read is for work that needs only the
image's time and geometry, so it leaves the image arrays unread.

.. list-table::
   :header-rows: 1

   * - Input
     - Normal read
     - With ``astrometry=True``
   * - FITS file
     - The FITS file alone: its headers and its arrays.
     - The FITS file's headers alone.
   * - PDS3 label
     - The label, and the whole FITS file beside it.
     - The label alone. The FITS file need not exist.

The metadata come from whatever was read, as :doc:`user_guide_metadata` describes. The
``data`` subfield, and the ``error`` and ``quality`` subfields of a calibrated image,
hold the image arrays except in an astrometry read.

An astrometry read of a PDS3 label reads the label alone, without the FITS header. It
accepts only the default ``geom`` and ``pointing`` options, which take the position and
the pointing from the SPICE kernels, and its list of calibrations is empty.

Finding the FITS file
=====================

Given a PDS3 label, the host looks in the label's directory for the FITS file that the
label's ``^IMAGE`` pointer gives. That name must end in ``.fit`` or ``.fits``, in any
case, as the name of every archive product does; a label that points to a file with
another suffix raises :exc:`~oops.host.PdsHostError`. The archive spells the FITS file's
name in the same case as the label's, so the host tries that case first, then the name
as the label spells it, then the name in lower case, and last in upper case. An
astrometry read does not look for the file, so its path subfields give the name in the
case of the label's own name, whatever its suffix.

The ``filepath``, ``basename``, ``abspath``, and ``image_url`` subfields always give the
FITS file, also when the input is the label. ``abspath`` is the file's local path, and
it is present only when that file exists: a local file, or a remote one that a read has
put in the local cache. After an astrometry read of a remote label, the FITS file has
not been retrieved, so ``abspath`` is absent unless an earlier read put the file in the
cache.

Parsed input
============

A label or FITS file that has already been opened can be passed in its place: a
:class:`~pdsparser.Pds3Label` is read as a PDS3 label, and an
:class:`~astropy.io.fits.HDUList` as the FITS file. The host closes an
:class:`~astropy.io.fits.HDUList` once it has copied out the headers and arrays.

A label or FITS file parsed from text or built in memory, rather than read from a file,
has no location, so the Snapshot has no ``filepath``, ``basename``, ``abspath``, or
``image_url``. An :class:`~astropy.io.fits.HDUList` of that kind holds everything the
host needs. A :class:`~pdsparser.Pds3Label` of that kind can be read with
``astrometry=True``, but a normal read raises :exc:`~oops.host.HostFormatError`,
because the host cannot find the FITS file beside a label that has no file path.

Options
=======

These options are particular to the New Horizons host:

``geom``
    ``"spice"``, the default, takes the spacecraft's position from the SPICE ephemerides.
    ``"fits"`` takes it from the FITS header, as a path in a straight line through the
    position and velocity that the header gives at the middle of the exposure.

``pointing``
    ``"spice"``, the default, takes the pointing from the C kernels (see
    :doc:`user_guide_geometry`). ``"fits"`` takes it from the FITS header's world
    coordinates and north clock angle. ``"fitsapp"`` treats the header's pointing as
    apparent, corrected for aberration, rather than geometric. ``"fits90"`` uses the
    header's pointing and adds 90° to the north clock angle of an image taken before
    2013, which corrects an error in the first archive delivery.

``fov_type``
    ``"fast"``, ``"slow"``, or ``"flat"``, as :doc:`user_guide_products` describes.

``return_all_planets``
    True lists in the ``spice_kernels`` subfield the kernels of every planet and of
    Pluto, and not only those of the target body.

``pds3_options``
    Options for the PDS3 label parser. The host starts with the parser's fast method and,
    if that cannot read a label, tries the strict method and then the loose one. A method
    given in ``pds3_options`` starts that sequence at the method given, and the compound
    method, if given, is tried alone.

``astrometry``
    True for an astrometry read, as described above.

``select``
    Must be None. Each file holds one image, so there is nothing to select, and any
    other value raises :exc:`~oops.OopsValueError`.

The standard OOPS overrides, ``target``, ``lightsource``, ``path``, ``frame``,
``cadence``, ``fov``, ``calibrations``, ``timeshift``, ``navigation``, ``tracker``, and
``parallel``, replace what the host would otherwise build. They are listed with
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`. Any other keyword
raises :exc:`~oops.OopsTypeError`, so that a misspelled argument is not silently
ignored.

Reading a label through oops.Host
=================================

:meth:`oops.host.Host.from_file` recognizes a LORRI image from its PDS3 label or its
FITS header and passes the path, with every other argument, to
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`. The host then parses
the label itself, trying the fast, strict, and loose methods in turn, so the two calls
read the same labels and give the same Snapshot, I/F calibration included.
