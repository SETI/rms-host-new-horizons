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
the pointing from the SPICE kernels, and it has no calibration.

Finding the FITS file
=====================

Given a PDS3 label, the host reads the FITS file that the label's ``^IMAGE`` pointer
names, in the label's directory. The labels write that name in upper case, such as
``LOR_0030713591_0X633_SCI.FIT``, while the archive stores each file in the same case as
its label, which is lower case. The host therefore writes the name in the case of the
label's own name: in lower case if the label's name is all lower case, in upper case if
it is all upper case, and as the pointer writes it otherwise.

The file subfields name the FITS file, also when the input is the label:

``filepath`` and ``basename``
    The path of the FITS file and its name.

``image_url``
    The absolute path or URL of the FITS file.

``image_abspath``
    The path of the FITS file on the local disk: the file itself if it is local, or the
    copy that the read downloaded into the local cache if it is remote. It is present
    whenever the read opened the FITS file, which an astrometry read of a label does
    not.

``label_url`` and ``label_abspath``
    The absolute path or URL of the PDS3 label and its path on the local disk, present
    only when the input is the label.

A FITS file that has already been opened as an :class:`~astropy.io.fits.HDUList` can be
passed in place of its path, which is how :meth:`oops.host.Host.from_file` passes it on.
The file subfields then name the file that it was opened from, and the host closes it
once it has copied out the headers and arrays.

Options
=======

These options are particular to the New Horizons host:

``geom``
    ``"spice"``, the default, takes the spacecraft's position from the SPICE ephemerides.
    ``"fits"`` takes it from the FITS header, as a path in a straight line about the Sun
    through the position and velocity that the header gives at the middle of the
    exposure.

``pointing``
    ``"spice"``, the default, takes the pointing from the C kernels (see
    :doc:`user_guide_geometry`). ``"fits"`` takes it from the FITS header's world
    coordinates and north clock angle. ``"fitsapp"`` treats the header's pointing as
    apparent, corrected for aberration, rather than geometric. ``"fits90"`` uses the
    header's pointing and adds 90° to the north clock angle of an image taken before
    2013, which corrects an error in the first archive delivery.

``fov_type``
    ``"fast"``, ``"slow"``, or ``"flat"``, as :doc:`user_guide_products` describes.

``pds3_options``
    Options for the PDS3 label parser. The host parses every label with the parser's
    strict method unless these options give another. Every LORRI label in the archive
    parses with the strict method as it is written.

``astrometry``
    True for an astrometry read, as described above.

``select``
    Chooses among the observations of the file, which holds one. None, the default, and
    0 return that observation, and a slice returns a tuple of it.

The standard OOPS overrides, ``target``, ``lightsource``, ``path``, ``frame``,
``cadence``, ``fov``, ``calibrations``, ``timeshift``, ``navigation``, ``tracker``, and
``parallel``, replace what the host would otherwise build. They are listed with
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`. A ``frame`` override,
such as a :class:`~oops.frame.Cmatrix`, is used as it is given. Any other keyword raises
:exc:`TypeError`, so that a misspelled argument is not silently ignored.

Reading a file through oops.Host
================================

:meth:`oops.host.Host.from_file` recognizes a LORRI image from its PDS3 label or its
FITS header and passes it, with every other argument, to
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file`. A label is passed as
its path, and the host parses it itself with the strict method, so the two calls read
the same labels.

The two calls differ in one default. :meth:`oops.host.Host.from_file` passes
``calibrations=[]`` unless the caller gives calibrations, and an empty list turns the
host's I/F calibration off. To have the I/F calibration, call
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` directly.
