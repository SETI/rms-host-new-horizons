##########################################################################################
# host_new_horizons/oops/lorri.py
##########################################################################################

import pathlib

import numpy as np

import julian

from astropy.io import fits as pyfits
from filecache  import FCPath
from pdsparser  import Pds3Label

import oops
from ._new_horizons import _NewHorizons

__all__ = ['NewHorizonsLORRI']

# The instrument kernel gives the boresight along -Z, so the camera frame is the SPICE
# frame with its y and z axes flipped
_SPICE_FRAME_NAME = 'NH_LORRI'
_SPICE_FRAME_ID = -98300
_FLIPPED_FRAME_ID = 'NH_LORRI_FLIPPED'
_SPICE_TO_FRAME = oops.Matrix3([[1, 0, 0],
                                [0,-1, 0],
                                [0, 0,-1]])

# The name of the full-resolution detector in the instrument kernel
_IK_NAME = 'NH_LORRI_1X1'

# The dimensions of an unbinned LORRI image in pixels, as the instrument kernel must give
# them
_LORRI_PIXELS = 1024

# The binning mode of an image, keyed by its number of lines
_BINNING_MODES = {1024: '1X1', 256: '4X4'}

# The pixel scale of the flat FOV in radians, the IFOV of the instrument kernel rounded
_FLAT_SCALE = 4.96e-6

# The suffixes of a FITS file; any other file is read as a PDS3 label
_FITS_SUFFIXES = ('.fit', '.fits')

_GEOMS = ('spice', 'fits')
_POINTINGS = ('spice', 'fits', 'fitsapp', 'fits90')
_FOV_TYPES = ('fast', 'slow', 'flat')

# Before 2013, the FITS headers of the first PDS delivery give a north clock angle 90
# degrees short; pointing="fits90" adds it back for an image taken in or before this year
_FITS90_LAST_YEAR = 2012

# The names of the extensions of a calibrated product that hold its error estimates and
# quality flags
_ERROR_EXTNAME = 'LORRI ERROR IMAGE'
_QUALITY_EXTNAME = 'LORRI QUALITY FLAG IMAGE'

# The target names of a FITS header that has no target. The calibration images of the
# launch phase, whose labels give the TARGET_NAME "CALIBRATION", give "---".
_NO_TARGET_NAMES = ('---', 'N/A')

# The ticks of the New Horizons spacecraft clock per second
_SCLK_TICKS = 50000

# The columns that identify the same image in an index and its supplement
_INDEX_KEYS = ['VOLUME_ID', 'FILE_SPECIFICATION_NAME']

# FITS header keywords that only carry text, left out of the `dict` subfield
_COMMENTARY_KEYWORDS = ('', 'COMMENT', 'HISTORY')

##########################################################################################
# The distortion model of
#   Owen Jr., W.M., 2011. New Horizons LORRI Geometric Calibration of August 2006. JPL
#   IOM 343L-11-002.
# These polynomials convert from X,Y (radians) to U,V (pixels).
##########################################################################################

_LORRI_F  =  2619.008   # mm
_LORRI_E2 =  2.696e-5   # / mm2
_LORRI_E5 =  1.988e-5   # / mm
_LORRI_E6 = -2.864e-5   # / mm
_LORRI_KX =  76.9231    # samples/mm
_LORRI_KY = -76.9231    # lines/mm

_LORRI_COEFF = np.zeros((4, 4, 2))
_LORRI_COEFF[1,0,0] = _LORRI_KX           * _LORRI_F
_LORRI_COEFF[3,0,0] = _LORRI_KX*_LORRI_E2 * _LORRI_F**3
_LORRI_COEFF[1,2,0] = _LORRI_KX*_LORRI_E2 * _LORRI_F**3
_LORRI_COEFF[1,1,0] = _LORRI_KX*_LORRI_E5 * _LORRI_F**2
_LORRI_COEFF[2,0,0] = _LORRI_KX*_LORRI_E6 * _LORRI_F**2

_LORRI_COEFF[0,1,1] = _LORRI_KY           * _LORRI_F
_LORRI_COEFF[2,1,1] = _LORRI_KY*_LORRI_E2 * _LORRI_F**3
_LORRI_COEFF[0,3,1] = _LORRI_KY*_LORRI_E2 * _LORRI_F**3
_LORRI_COEFF[0,2,1] = _LORRI_KY*_LORRI_E5 * _LORRI_F**2
_LORRI_COEFF[1,1,1] = _LORRI_KY*_LORRI_E6 * _LORRI_F**2

# The inverse distortion model, from U,V (pixels) to X,Y (radians). These coefficients
# were computed by numerically solving the above polynomials. Maximum errors from applying
# the original distortion model and then inverting:
#   U DIFF MIN MAX -0.00261017184761 0.0018410196501
#   V DIFF MIN MAX -0.00263250108583 0.00209864359397
_LORRI_INV_COEFF = np.zeros((4, 4, 2))
_LORRI_INV_COEFF[:,:,0] = [
    [  5.62144901e-10,  1.80741669e-15, -3.62872755e-15, -7.65201036e-21],
    [  4.96369475e-06,  1.27370808e-12, -2.24498555e-14,  0.00000000e+00],
    [  1.83505816e-12,  1.87593608e-18,  0.00000000e+00,  0.00000000e+00],
    [ -2.24714913e-14,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]]
_LORRI_INV_COEFF[:,:,1] = [
    [ -3.90196341e-10, -4.96369475e-06, -1.27377633e-12,  2.24721956e-14],
    [ -1.96689309e-15, -1.83495964e-12, -1.87499462e-18,  0.00000000e+00],
    [  2.51861242e-15,  2.24491504e-14,  0.00000000e+00,  0.00000000e+00],
    [  7.85425602e-21,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00]]


class NewHorizonsLORRI(oops.Host):
    """The New Horizons LORRI host: Snapshot constructors for Long Range Reconnaissance
    Imager images."""

    NAME = 'New Horizons LORRI'

    _INSTRUMENT_KERNEL = None
    _FOVS = {}
    _initialized = False

    @staticmethod
    def from_file(fileinfo, *, geom='spice', pointing='spice', fov_type='fast',
                  pds3_options=None, astrometry=False, target=None, lightsource=None,
                  path=None, frame=None, cadence=None, fov=None, calibrations=None,
                  timeshift=None, navigation=None, tracker=None, parallel=None, **kwargs):
        """A Snapshot object based on a given New Horizons LORRI image file.

        Every value comes from the FITS file. Given a PDS3 label, the host reads the FITS
        file that the label's ``^IMAGE`` pointer names, so the label and its FITS file
        give the same Snapshot.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | HDUList): The path to the
                FITS file or its PDS3 label, or its parsed content as a Pds3Label or an
                HDUList. The FITS file is closed once it has been read, including an
                HDUList that was passed in.
            geom (str, optional): "spice" to take the spacecraft's position from the SP
                kernels; "fits" to take it from the FITS header, as a linear path relative
                to the Sun.
            pointing (str, optional): "spice" to take the pointing from the C kernels;
                "fits" to take it from the FITS header; "fitsapp" to treat the pointing in
                the FITS header as apparent rather than geometric; "fits90" to use the
                FITS header's pointing, adding 90 degrees to the north clock angle of an
                image taken before 2013, which corrects an error in the first PDS
                delivery.
            fov_type (str, optional): "fast" to model the camera distortion with a
                polynomial and its numerical inverse; "slow" to invert the polynomial by
                Newton's method; "flat" for a field of view without distortion.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser.
            astrometry (bool, optional):
                True to specify that the returned Observations will only be used for
                timing and astrometry. In this case, no data arrays are read and returned
                in the Observation.
            target (str | oops.Body, optional): Override for the default target Body of
                the observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of the observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument.
            cadence (oops.Cadence, optional): Override for the Cadence of the observation.
            fov (oops.FOV, optional): Override for the default FOV of the observing
                instrument.
            calibrations (oops.Calibration | list[oops.Calibration], optional): Override
                for the calibration or list of calibrations.
            timeshift (tuple[float, str], optional): Assign a Fittable time shift to one
                or more attributes of the Observation. The first item is the initial time
                shift in seconds. The second item is a string indicating which attributes
                are to be time-shifted: "T" for the start/stop time, "P" for the path,
                and/or "F" for the frame. If multiple letters are given, the time shifts
                are coupled to whichever attribute is listed first.
            navigation (tuple[float, ...], optional): Two or three rotation angles in
                radians to apply to the frame in order to align the calculated geometry
                with the observed geometry. In this case, a fittable Navigation frame
                wraps the default frame. Specify three angles to include a rotation about
                the optic axis or two for a pointing offset without rotation. Use (0,0) or
                (0,0,0) as the input if you have no starting guess.
            tracker (str, optional): Assign a TrackerFrame to the Observation, to ensure
                that the target body remains at a fixed position within the FOV. Use one
                of "start", "midtime", and "end", indicating the time within the
                Observation at which the target body's calculated position within the FOV
                is accurate.
            parallel (Observation, optional): A parallel Observation (same origin and
                time, different frame and FOV) relative to which this Observation has a
                fixed offset. Note that no more than one of `parallel` and `tracker` can
                be specified.
            **kwargs: Additional keyword arguments; they are accepted and ignored. This
                includes `select` and `calibration`. The data array is always the image
                as the file stores it.

        Returns:
            Snapshot: The observation, with subfields `dict`, `headers`, `host`,
            `instrument`, `binning_mode`, `label_target`,
            `spacecraft_clock_start_count`, `spacecraft_clock_stop_count`,
            `spice_kernels`, `filepath`, `basename`, `spice_to_frame`, `spice_frame_name`,
            `spice_frame_id`, `abspath` and `image_url` inserted, and `data`, `error` and
            `quality` unless `astrometry` is True. `dict` is the primary FITS header
            without its COMMENT and HISTORY cards, and `headers` is the list of the
            headers of every HDU. `data` is the primary image as stored, which in a raw
            product includes the columns beyond the field of view; `error` and `quality`
            are the error estimates and quality flags of a calibrated product, and None
            for a raw one. `label_target` is the FITS header's target name as written.
            The two clock counts are the start and end of the exposure, in the form
            "SSSSSSSSSS:TTTTT" that the PDS3 label uses; either is None if the header
            does not give it. `filepath` and the other path subfields name the FITS file
            even when a label was given.

        Raises:
            FileNotFoundError: If the FITS file that a PDS3 label names does not exist.
            oops.OopsValueError: If `geom`, `pointing` or `fov_type` is not one of the
                values listed above.
            oops.host.HostFormatError: If the image is not 1024 or 256 lines, if the
                FITS file is not valid FITS, or if a PDS3 label names no FITS file.
            oops.host.PdsHostError: If the PDS3 label cannot be parsed.
            oops.host.WrongHostError: If `fileinfo` is not a New Horizons LORRI PDS3 label
                or FITS file.
            OSError: If the file cannot be read.
        """

        NewHorizonsLORRI._check_option('geom', geom, _GEOMS)
        NewHorizonsLORRI._check_option('pointing', pointing, _POINTINGS)
        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)

        NewHorizonsLORRI.initialize()   # define everything the first time through; use
                                        # defaults unless initialize() was called

        # Default Pds3Label options are method="fast", sources=True
        pds3_options = oops.Host._pds3_options_with_defaults(pds3_options)

        # A FITS file is opened here rather than by oops.Host._read_fileinfo(), which
        # leaves the file open when it belongs to another host
        if isinstance(fileinfo, (str, pathlib.Path, FCPath)):
            filepath = FCPath(fileinfo)
            if filepath.suffix.lower() in _FITS_SUFFIXES:
                parsed = NewHorizonsLORRI._open_fits(filepath)
            else:
                parsed = NewHorizonsLORRI._read_fileinfo(filepath, formats='P',
                                                         pds3_options=pds3_options)
        else:
            filepath = None
            parsed = fileinfo

        if isinstance(parsed, Pds3Label):
            filepath = NewHorizonsLORRI._fits_path(parsed)
            hdulist = NewHorizonsLORRI._open_fits(filepath)
        else:
            hdulist = parsed
            if filepath is None:
                filepath = FCPath(hdulist.filename())

        # Everything needed is copied out of the file, which is then closed, including
        # an HDUList that oops.Host.from_file() or the caller opened
        try:
            headers = [hdu.header.copy() for hdu in hdulist]
            arrays = None
            if not astrometry:
                arrays = NewHorizonsLORRI._read_arrays(hdulist)
        finally:
            hdulist.close()

        overrides = {
            'target'      : target,
            'lightsource' : lightsource,
            'path'        : path,
            'frame'       : frame,
            'cadence'     : cadence,
            'fov'         : fov,
            'calibrations': calibrations,
            'timeshift'   : timeshift,
            'navigation'  : navigation,
            'tracker'     : tracker,
            'parallel'    : parallel,
        }

        return NewHorizonsLORRI._make_snapshot(headers, filepath, arrays, geom=geom,
                                               pointing=pointing, fov_type=fov_type,
                                               overrides=overrides)

    @staticmethod
    def from_index(filepath, *, supplement=None, row_dicts=None, pds3_options=None,
                   filter=None, fov_type='fast', select=None, target=None,
                   lightsource=None, path=None, frame=None, cadence=None, fov=None,
                   calibrations=None, timeshift=None, navigation=None, tracker=None,
                   parallel=None, **kwargs):
        """A list of Snapshot objects, one for each row of a New Horizons LORRI index.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The full path to a New Horizons LORRI
                index file or its PDS label, such as ``NHJULO_2001_index.lbl``.
            supplement (str | pathlib.Path | FCPath, optional): The path to the PDS label
                of the supplemental index, whose columns are merged into each row. Its
                rows must name the same files in the same order as the index. Where both
                tables give a column, the supplement's value is kept, because the
                supplement is spelled as the product labels are. If not specified, the
                directory of `filepath` is searched for a file whose name matches except
                that "_index" is replaced by "_supplemental_index" in upper or lower case.
            row_dicts (list[dict[str, Any]], optional): A list of dictionaries, one per
                row of the index. Include if `filepath` and `supplement` have already been
                read. If provided, `supplement` is ignored.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser.
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, will be
                included in the returned list.
            fov_type (str, optional): "fast", "slow" or "flat", as for :meth:`from_file`.
            select (int | slice | str | tuple[int | slice | str, ...]): Accepted for
                compatibility with :meth:`oops.Host.from_index` and ignored, because each
                row describes a single observation.
            target (str | oops.Body, optional): Override for the default target Body of
                every observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default LightSource
                of every observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument.
            cadence (oops.Cadence, optional): Override for the Cadence of the observation.
            fov (oops.FOV, optional): Override for the FOV of every observation,
                replacing the FOV that each row's binning mode would otherwise define.
            calibrations (oops.Calibration | list[oops.Calibration], optional): Override
                for the calibration or list of calibrations.
            timeshift (tuple[float, str], optional): Assign a Fittable time shift to one
                or more attributes of the Observation. The first item is the initial time
                shift in seconds. The second item is a string indicating which attributes
                are to be time-shifted: "T" for the start/stop time, "P" for the path,
                and/or "F" for the frame. If multiple letters are given, the time shifts
                are coupled to whichever attribute is listed first.
            navigation (tuple[float, ...], optional): Two or three rotation angles in
                radians to apply to the frame in order to align the calculated geometry
                with the observed geometry. In this case, a fittable Navigation frame
                wraps the default frame. Specify three angles to include a rotation about
                the optic axis or two for a pointing offset without rotation. Use (0,0) or
                (0,0,0) as the input if you have no starting guess.
            tracker (str, optional): Assign a TrackerFrame to the Observation, to ensure
                that the target body remains at a fixed position within the FOV. Use one
                of "start", "midtime", and "end", indicating the time within the
                Observation at which the target body's calculated position within the FOV
                is accurate.
            parallel (Observation, optional): A parallel Observation (same origin and
                time, different frame and FOV) relative to which this Observation has a
                fixed offset. Note that no more than one of `parallel` and `tracker` can
                be specified.
            **kwargs: Additional keyword arguments; they are accepted and ignored.

        Returns:
            list[Snapshot]: One observation per row of the index, each with subfields
            `dict`, `host`, `instrument`, `binning_mode`, `label_target`,
            `spacecraft_clock_start_count`, `spacecraft_clock_stop_count`,
            `spice_kernels`, `filepath`, `basename`, `spice_to_frame`, `spice_frame_name`
            and `spice_frame_id` inserted. `dict` is the row, merged with its row of the
            supplement. `filepath` is the VOLUME_ID and the FILE_SPECIFICATION_NAME, with
            the label's suffix replaced by the FITS file's, joined by a slash.
            `label_target` is the row's TARGET_NAME, and the target is that name repaired
            as for :meth:`from_file`.

        Raises:
            FileNotFoundError: If the supplement does not exist.
            oops.OopsValueError: If `fov_type` is not "fast", "slow" or "flat".
            oops.host.HostFormatError: If a row's BINNING_MODE is not "1X1" or "4X4".
            oops.host.HostIndexError: If the supplement does not have the same number of
                rows as the index, if a row's VOLUME_ID or FILE_SPECIFICATION_NAME
                differs from the index's, or if `supplement` is not specified and the
                name of `filepath` contains no "_index".
            oops.host.WrongHostError: If `filepath` is not a New Horizons LORRI index.
            OSError: If `filepath`, or a `supplement` that is specified, does not exist or
                cannot be read.
        """

        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)
        NewHorizonsLORRI.initialize()

        # Read and filter the index table
        if row_dicts is None:
            supplement = supplement or oops.Host.SUPPLEMENT
            row_dicts = NewHorizonsLORRI._read_index_rows(
                filepath, supplement=supplement,
                validate=NewHorizonsLORRI._match_index_rows, pds3_options=pds3_options)
        row_dicts = oops.Host._filter_index_rows(row_dicts, filter)

        overrides = {
            'target'      : target,
            'lightsource' : lightsource,
            'path'        : path,
            'frame'       : frame,
            'cadence'     : cadence,
            'fov'         : fov,
            'calibrations': calibrations,
            'timeshift'   : timeshift,
            'navigation'  : navigation,
            'tracker'     : tracker,
            'parallel'    : parallel,
        }

        # Create the list of Snapshot objects
        snapshots = []
        for row_dict in row_dicts:
            obs = NewHorizonsLORRI._snapshot_from_row(row_dict, fov_type=fov_type,
                                                      overrides=overrides)
            snapshots.append(obs)

        return snapshots

    ######################################################################################
    # Support methods
    ######################################################################################

    @staticmethod
    def _check_option(name, value, allowed):
        """Raise OopsValueError unless `value` is one of `allowed`."""

        if value not in allowed:
            raise oops.OopsValueError(f'invalid {name} {value!r}; use one of '
                                      + ', '.join(repr(item) for item in allowed))

    @staticmethod
    def _fits_path(label):
        """The FITS file that a PDS3 label's ``^IMAGE`` pointer names.

        The label may spell the name in a different case from the file system. The
        archive names the FITS file in the same case as the label file, so the name is
        tried first in the case of the label's own name, then as written, then in lower
        and in upper case. On a file system that ignores case, the first spelling tried
        is the one returned.

        Parameters:
            label (Pds3Label): The parsed label.

        Returns:
            FCPath: The path to the FITS file, in the label's directory.

        Raises:
            FileNotFoundError: If no spelling of the name exists.
            oops.host.HostFormatError: If the label has no ``^IMAGE`` pointer to a
                separate file.
        """

        name = label.get('^IMAGE')
        if not isinstance(name, str):
            raise oops.host.HostFormatError(f'PDS3 label names no FITS file: '
                                            f'{label.filepath}')

        label_path = FCPath(label.filepath)
        if label_path.name == label_path.name.lower():
            first = name.lower()
        elif label_path.name == label_path.name.upper():
            first = name.upper()
        else:
            first = name

        parent = label_path.parent
        for option in dict.fromkeys([first, name, name.lower(), name.upper()]):
            path = parent / option
            if path.exists():
                return path

        raise FileNotFoundError(f'FITS file {name} of PDS3 label not found: '
                                f'{label.filepath}')

    @staticmethod
    def _open_fits(filepath):
        """Open a LORRI FITS file.

        Parameters:
            filepath (FCPath): The path to the file.

        Returns:
            HDUList: The open file.

        Raises:
            oops.host.HostFormatError: If the file is not valid FITS.
            oops.host.WrongHostError: If the file is not a New Horizons LORRI FITS file;
                the file is closed first.
            OSError: If the file does not exist or cannot be read.
        """

        local_path = filepath.retrieve()    # a missing file is not a format error
        try:
            hdulist = pyfits.open(local_path)
        except OSError as err:
            raise oops.host.HostFormatError(f'not a valid FITS file: {filepath}') from err

        if not NewHorizonsLORRI._detect_in_fits(hdulist):
            hdulist.close()
            raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} FITS file: '
                                           f'{filepath}')
        return hdulist

    @staticmethod
    def _read_arrays(hdulist):
        """The image, error and quality arrays of an open LORRI FITS file.

        Parameters:
            hdulist (HDUList): The FITS file.

        Returns:
            dict[str, numpy.ndarray | None]: The arrays keyed by subfield name: `data`,
            the primary image, and `error` and `quality`, the extensions of a calibrated
            product that hold them, or None where the file has no such extension. Each
            array is a copy, so it outlives the file.
        """

        arrays = {'data': np.array(hdulist[0].data), 'error': None, 'quality': None}
        for hdu in hdulist[1:]:
            extname = str(hdu.header.get('EXTNAME', '')).strip().upper()
            if extname == _ERROR_EXTNAME:
                arrays['error'] = np.array(hdu.data)
            elif extname == _QUALITY_EXTNAME:
                arrays['quality'] = np.array(hdu.data)

        return arrays

    @staticmethod
    def _make_snapshot(headers, filepath, arrays, *, geom, pointing, fov_type,
                       overrides):
        """A Snapshot built from the headers and arrays of a LORRI FITS file.

        Parameters:
            headers (list[astropy.io.fits.Header]): The header of every HDU.
            filepath (FCPath): The path to the FITS file.
            arrays (dict[str, numpy.ndarray | None] | None): The arrays as
                :meth:`_read_arrays` returns them; None to omit them.
            geom (str): "spice" or "fits", as for :meth:`from_file`.
            pointing (str): "spice", "fits", "fitsapp" or "fits90", as for
                :meth:`from_file`.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.host.HostFormatError: If the image is not 1024 or 256 lines.
        """

        header = headers[0]

        # Key information from the header
        texp = header['EXPTIME']
        tdb_midtime = header['SPCSCET']
        tstart = tdb_midtime - texp / 2.

        lines = header['NAXIS2']
        if lines not in _BINNING_MODES:
            raise oops.host.HostFormatError(f'unrecognized LORRI image shape: '
                                            f'{lines} lines in {filepath}')
        binning_mode = _BINNING_MODES[lines]

        label_target = header['SPCCBTNM']
        target_name = NewHorizonsLORRI._repair_target(label_target)

        # Define the path and frame
        if geom == 'spice':
            default_path = _NewHorizons.PATH_ID
        else:
            default_path = NewHorizonsLORRI._fits_path_from_header(header, filepath,
                                                                   tdb_midtime)

        using_cks = overrides['frame'] is None and pointing == 'spice'
        if pointing == 'spice':
            if overrides['frame'] is None:
                NewHorizonsLORRI._define_camera_frame()
            default_frame = _SPICE_FRAME_NAME
        else:
            default_frame = NewHorizonsLORRI._fits_frame_from_header(
                header, filepath, pointing=pointing, path=default_path,
                tdb_midtime=tdb_midtime)

        # Define the Snapshot parameters
        params = {
            'cadence'     : oops.cadence.SnapCadence(tstart, texp),
            'fov'         : NewHorizonsLORRI._FOVS[binning_mode, fov_type],
            'path'        : default_path,
            'frame'       : default_frame,
            'target'      : target_name,
            'lightsource' : 'SUN',
            'calibrations': [],         # TBD
        }

        # Apply standard overrides
        oops.Host._apply_overrides(params, filepath, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if arrays is not None:
            obs.insert_subfields(arrays)

        # Extra attributes: use standardized names wherever practical
        dict_ = {key: value for (key, value) in header.items()
                 if key not in _COMMENTARY_KEYWORDS}
        obs.insert_subfield('dict', dict_)
        obs.insert_subfield('headers', headers)
        stop_count = NewHorizonsLORRI._sclk_from_met(header.get('STOPMET'))
        NewHorizonsLORRI._insert_common_subfields(
            obs, binning_mode=binning_mode, label_target=label_target, target=target_name,
            start_count=header.get('SPCSCLK0'), stop_count=stop_count, ck=using_cks)
        return obs

    @staticmethod
    def _snapshot_from_row(row_dict, *, fov_type, overrides):
        """A Snapshot built from a row of an index merged with its supplement.

        Parameters:
            row_dict (dict[str, Any]): The row.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_index`.

        Raises:
            oops.host.HostFormatError: If the row's BINNING_MODE is not "1X1" or "4X4".
        """

        tstart = julian.tdb_from_tai(julian.tai_from_iso(row_dict['START_TIME']))
        texp = row_dict['EXPOSURE_DURATION']

        binning_mode = row_dict['BINNING_MODE'].upper()
        if binning_mode not in _BINNING_MODES.values():
            raise oops.host.HostFormatError(f'unrecognized LORRI BINNING_MODE: '
                                            f'{row_dict["BINNING_MODE"]!r}')

        # The index names the label; the observation names the FITS file beside it
        label_name = row_dict['FILE_SPECIFICATION_NAME']
        stem, _, suffix = label_name.rpartition('.')
        fits_suffix = 'FIT' if suffix.isupper() else 'fit'
        filepath = row_dict['VOLUME_ID'] + '/' + stem + '.' + fits_suffix

        using_cks = overrides['frame'] is None
        if using_cks:
            NewHorizonsLORRI._define_camera_frame()

        # The table reader returns NumPy strings; the subfields are Python strings, as
        # from a FITS file
        label_target = str(row_dict['TARGET_NAME'])
        target_name = NewHorizonsLORRI._repair_target(label_target)
        params = {
            'cadence'     : oops.cadence.SnapCadence(tstart, texp),
            'fov'         : NewHorizonsLORRI._FOVS[binning_mode, fov_type],
            'path'        : _NewHorizons.PATH_ID,
            'frame'       : _SPICE_FRAME_NAME,
            'target'      : target_name,
            'lightsource' : 'SUN',
            'calibrations': [],         # TBD
        }
        oops.Host._apply_overrides(params, filepath, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        obs.insert_subfield('dict', row_dict)
        NewHorizonsLORRI._insert_common_subfields(
            obs, binning_mode=binning_mode, label_target=label_target, target=target_name,
            start_count=str(row_dict['SPACECRAFT_CLOCK_START_COUNT']),
            stop_count=str(row_dict['SPACECRAFT_CLOCK_STOP_COUNT']), ck=using_cks)
        return obs

    @staticmethod
    def _insert_common_subfields(obs, *, binning_mode, label_target, target, start_count,
                                 stop_count, ck):
        """Insert the subfields that file and index observations share.

        Parameters:
            obs (Snapshot): The observation.
            binning_mode (str): "1X1" or "4X4".
            label_target (str): The target name as the file or index writes it.
            target (str): The name of the default target, whose kernels are included in
                `spice_kernels`.
            start_count (str | None): The spacecraft clock count at the start of the
                exposure.
            stop_count (str | None): The spacecraft clock count at the end of the
                exposure.
            ck (bool): True if the pointing comes from the C kernels, so that they are
                included in `spice_kernels`.
        """

        obs.insert_subfield('host', 'New Horizons')
        obs.insert_subfield('instrument', 'LORRI')
        obs.insert_subfield('binning_mode', binning_mode)
        obs.insert_subfield('label_target', label_target)
        obs.insert_subfield('spacecraft_clock_start_count', start_count)
        obs.insert_subfield('spacecraft_clock_stop_count', stop_count)

        # With a custom frame, pointing never came from a CK, so any CK that happens to be
        # furnished is unrelated to this observation and must not be reported as used.
        kernels = _NewHorizons.used_kernels(obs.time, 'lorri', target=target, ck=ck)
        obs.insert_subfield('spice_kernels', kernels)
        obs.insert_subfield('spice_to_frame', _SPICE_TO_FRAME)
        obs.insert_subfield('spice_frame_name', _SPICE_FRAME_NAME)
        obs.insert_subfield('spice_frame_id', _SPICE_FRAME_ID)

    @staticmethod
    def _repair_target(name):
        """The target name of a FITS header, repaired.

        A header gives "---" or "N/A" when the image has no target, as the calibration
        images of the launch phase do, and some images of Pluto give a name that begins
        with "PLUTO".

        Parameters:
            name (str): The header's SPCCBTNM.

        Returns:
            str: The name of the target body: "NONE" for an image with no target, "PLUTO"
            for a name beginning with "PLUTO", and otherwise `name` itself.
        """

        if name.strip() in _NO_TARGET_NAMES:
            return 'NONE'
        if name.startswith('PLUTO'):
            return 'PLUTO'
        return name

    @staticmethod
    def _sclk_from_met(met):
        """A mission elapsed time converted to a clock count in the PDS3 label's form.

        Parameters:
            met (float | None): Seconds of the spacecraft clock, as the FITS header's
                STARTMET and STOPMET give them.

        Returns:
            str | None: The count as "SSSSSSSSSS:TTTTT", ten digits of whole seconds and
            five of the 1/50000-second ticks past them; None if `met` is None.
        """

        if met is None:
            return None

        (seconds, ticks) = divmod(round(met * _SCLK_TICKS), _SCLK_TICKS)
        return f'{seconds:010d}:{ticks:05d}'

    @staticmethod
    def _fits_path_from_header(header, filepath, tdb_midtime):
        """The spacecraft's path as the FITS header gives it.

        Parameters:
            header (astropy.io.fits.Header): The primary header.
            filepath (FCPath): The path to the FITS file, whose name makes the path's ID
                unique.
            tdb_midtime (float): The middle of the exposure in seconds TDB.

        Returns:
            oops.Path: A linear path relative to the Sun, through the spacecraft's
            position and velocity at mid-exposure.
        """

        # The header gives the vectors from the spacecraft to the Sun
        pos = oops.Vector3([-header['SPCSSCX'], -header['SPCSSCY'], -header['SPCSSCZ']])
        vel = oops.Vector3([-header['SPCSSCVX'], -header['SPCSSCVY'],
                            -header['SPCSSCVZ']])

        # The path ID has to be unique to this observation
        sc_path = oops.path.LinearPath((pos, vel), tdb_midtime, 'SUN', frame='J2000',
                                       path_id='.NH_PATH_' + filepath.name)
        return oops.Path.as_waypoint(sc_path)

    @staticmethod
    def _fits_frame_from_header(header, filepath, *, pointing, path, tdb_midtime):
        """The camera frame as the FITS header gives it.

        Parameters:
            header (astropy.io.fits.Header): The primary header.
            filepath (FCPath): The path to the FITS file, whose name makes the frame's ID
                unique.
            pointing (str): "fits", "fitsapp" or "fits90", as for :meth:`from_file`.
            path (str | oops.Path): The spacecraft's path, needed for "fitsapp".
            tdb_midtime (float): The middle of the exposure in seconds TDB.

        Returns:
            oops.Frame: A fixed frame relative to J2000.
        """

        # The boresight; FITS pixel (1,1) is the center of the lower left pixel
        u_center = header['NAXIS1'] / 2. + 0.5
        v_center = header['NAXIS2'] / 2. + 0.5
        (ra_deg, dec_deg) = NewHorizonsLORRI._radec_from_uv(u_center, v_center, header)
        north_clock_deg = header['SPCEMEN']

        if pointing == 'fitsapp':
            # Insert the apparent vector as the actual one to reverse the aberration
            event = oops.Path.as_path(path).event_at_time(tdb_midtime)
            event.neg_arr_j2000 = oops.Vector3.from_ra_dec_length(
                ra_deg * oops.RPD, dec_deg * oops.RPD, recursive=False)
            (ra, dec) = event.ra_and_dec(apparent=True)
            ra_deg  = ra.vals  * oops.DPR
            dec_deg = dec.vals * oops.DPR

        elif pointing == 'fits90':
            year = int(header['SPCUTCID'][:4])
            if year <= _FITS90_LAST_YEAR:
                north_clock_deg += 90.

        frame = oops.frame.Cmatrix.from_ra_dec(ra_deg, dec_deg, north_clock_deg,
                                               oops.Frame.J2000,
                                               frame_id='.NH_FRAME_' + filepath.name)
        return oops.Frame.as_wayframe(frame)

    @staticmethod
    def _radec_from_uv(u, v, header):
        """The right ascension and declination at FITS pixel coordinates `(u,v)`.

        The gnomonic projection of the header's WCS parameters, adapted from STScI
        source. Here `(u,v) = (1,1)` is the center of the lower left pixel, unlike oops,
        in which `(0,0)` is the corner of the image and pixel centers have half-integer
        values.

        Parameters:
            u (float): The FITS horizontal pixel coordinate.
            v (float): The FITS vertical pixel coordinate.
            header (astropy.io.fits.Header): The FITS header, which provides the WCS
                parameters.

        Returns:
            tuple[float, float]: The right ascension and declination in degrees.
        """

        dx = u - header['CRPIX1']
        dy = v - header['CRPIX2']

        xi  = (header['CD1_1'] * dx + header['CD1_2'] * dy) / oops.DPR
        eta = (header['CD2_1'] * dx + header['CD2_2'] * dy) / oops.DPR

        ra0  = header['CRVAL1'] / oops.DPR
        dec0 = header['CRVAL2'] / oops.DPR

        denom = np.cos(dec0) - eta * np.sin(dec0)

        ra  = oops.DPR * (np.arctan(xi / denom) + ra0)
        dec = oops.DPR * np.arctan((eta * np.cos(dec0) + np.sin(dec0))
                                   / np.sqrt(denom**2 + xi**2))

        return (ra, dec)

    @staticmethod
    def _match_index_rows(index_row, supplement_row):
        """True if a row of the index and a row of its supplement describe one image.

        The rows match on :data:`_INDEX_KEYS`. The other columns the two tables share are
        written differently: the index gives each clock count with its partition, as in
        "1/0030598438:48850", and rounds the stop count differently from the label,
        whereas the supplement copies the label. :meth:`oops.Host._read_index_rows` keeps
        the index's value of a shared column, so this function copies the supplement's
        row into the index's before that merge, so that every shared column is spelled as
        the product labels spell it.

        Parameters:
            index_row (dict[str, Any]): The row of the index, updated in place.
            supplement_row (dict[str, Any]): The row of the supplement.

        Returns:
            bool: True if the rows match.

        Raises:
            oops.host.HostIndexError: If the rows differ in a column of
                :data:`_INDEX_KEYS`.
        """

        for key in _INDEX_KEYS:
            if index_row[key] != supplement_row[key]:
                raise oops.host.HostIndexError(f'{key} value mismatch: '
                                               f'{index_row[key]!r}, '
                                               f'{supplement_row[key]!r}')

        index_row.update(supplement_row)
        return True

    ######################################################################################
    # Detectors
    ######################################################################################

    @staticmethod
    def _detect_in_pds3(label):
        """True if the given parsed PDS3 label describes this host's data."""

        return (label.get('INSTRUMENT_HOST_NAME', '') == 'NEW HORIZONS'
                and label.get('INSTRUMENT_ID', '') == 'LORRI')

    @staticmethod
    def _detect_in_fits(hdulist):
        """True if the given FITS HDUList describes this host's data."""

        header = hdulist[0].header
        return (str(header.get('MISSION', '')).strip().upper() == 'NEW HORIZONS'
                and str(header.get('INSTRU', '')).strip().upper() == 'LOR')

    @staticmethod
    def _detect_in_index(label):
        """True if the given index label refers to this host/instrument.

        A New Horizons LORRI index label names the host and instrument at its top level,
        so unlike the other detectors, this one never needs to look at the rows.
        """

        return NewHorizonsLORRI._detect_in_pds3(label)

    @staticmethod
    def _detect_in_row(row_dict):
        """True if the given row of an index file refers to this host/instrument."""

        return (row_dict.get('INSTRUMENT_HOST', '') == 'NEW HORIZONS'
                and row_dict.get('INSTRUMENT_ID', '') == 'LORRI')

    ######################################################################################
    # Initialization
    ######################################################################################

    @staticmethod
    def initialize(*, asof=None, meta=None, planets=None, irregulars=True):
        """Initialize key information about the LORRI instrument.

        Furnishes the mission's kernels and builds the fields of view. Must be called
        first; :meth:`from_file` and :meth:`from_index` call it with the defaults. After
        the first call, later calls to this function are ignored.

        Parameters:
            asof (str, optional): Only use SPICE kernels that existed before this date;
                None to ignore.
            meta (str, optional): The name of a metakernel to furnish instead of the
                kernels selected from the database.
            planets (list, optional): A list of planets to pass to
                :meth:`~oops.Body.define_solar_system`. None or 0 means all.
            irregulars (bool, optional): True to include the irregular satellites; False
                otherwise.

        Raises:
            oops.host.HostError: If the instrument kernel does not give a 1024x1024 image.
        """

        # Quick exit after first call
        if NewHorizonsLORRI._initialized:
            return

        _NewHorizons.initialize(asof=asof, meta=meta, planets=planets,
                                irregulars=irregulars)

        # Load the instrument kernel
        (NewHorizonsLORRI._INSTRUMENT_KERNEL, _) = _NewHorizons.spice_instrument_kernel(
            'LORRI', asof=asof)
        info = NewHorizonsLORRI._INSTRUMENT_KERNEL['INS'][_IK_NAME]

        shape = (info['PIXEL_SAMPLES'], info['PIXEL_LINES'])
        if shape != (_LORRI_PIXELS, _LORRI_PIXELS):
            raise oops.host.HostError(f'unsupported LORRI image shape {shape} in the '
                                      'LORRI kernel')

        # Display directions: [u,v] = [right,down]
        full_fov_slow = oops.fov.PolynomialFOV(shape, coefft_uv_from_xy=_LORRI_COEFF,
                                               coefft_xy_from_uv=None)
        full_fov_fast = oops.fov.PolynomialFOV(shape, coefft_uv_from_xy=_LORRI_COEFF,
                                               coefft_xy_from_uv=_LORRI_INV_COEFF)
        full_fov_flat = oops.fov.FlatFOV((_FLAT_SCALE, -_FLAT_SCALE), shape)

        # The FOVs, keyed by the binning mode and the FOV type
        fovs = {}
        for (fov_type, full_fov) in [('slow', full_fov_slow), ('fast', full_fov_fast),
                                     ('flat', full_fov_flat)]:
            fovs['1X1', fov_type] = full_fov
            fovs['4X4', fov_type] = oops.fov.SubsampledFOV(full_fov, 4)

        NewHorizonsLORRI._FOVS = fovs
        NewHorizonsLORRI._initialized = True

    @staticmethod
    def _define_camera_frame():
        """Register the camera frame, NH_LORRI, and the SPICE frame beneath it.

        :meth:`initialize` must have been called first. Built lazily (and only once) so
        that observations using a custom frame never construct or depend on the SPICE
        camera frame. The SPICE frame NH_LORRI is registered as NH_LORRI_FLIPPED, and the
        camera frame NH_LORRI flips its y and z axes, because the instrument kernel gives
        the boresight along -Z.
        """

        # The frame is tested against the Frame registry rather than against a flag on
        # this class. Registering a second frame under an ID that already exists does not
        # fail; it is renamed with a suffix and never found again, so the test has to ask
        # the registry what it actually holds.
        if oops.Frame.frame_id_exists(_SPICE_FRAME_NAME):
            return

        flipped = oops.frame.SpiceFrame(_SPICE_FRAME_NAME, frame_id=_FLIPPED_FRAME_ID)
        _ = oops.frame.Cmatrix(_SPICE_TO_FRAME, flipped, frame_id=_SPICE_FRAME_NAME)

    @staticmethod
    def _reset():
        """Reset the internal New Horizons LORRI parameters.

        Can be useful for debugging.
        """

        NewHorizonsLORRI._INSTRUMENT_KERNEL = None
        NewHorizonsLORRI._FOVS = {}
        NewHorizonsLORRI._initialized = False
        _NewHorizons.reset()


NewHorizonsLORRI._register()

##########################################################################################
