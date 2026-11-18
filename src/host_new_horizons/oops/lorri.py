##########################################################################################
# host_new_horizons/oops/lorri.py
##########################################################################################

import datetime as dt
import pathlib

import numpy as np

import cspyce
import julian

from astropy.io import fits as pyfits
from filecache  import FCPath
from pdsparser  import Pds3Label, PdsError

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

# The target name of the calibration images of the launch phase in the label and index,
# which means the image has no target. Their FITS header gives "---", a null value.
_NO_TARGET_NAMES = ('CALIBRATION',)

# Values that mean a metadata item is missing: the PDS3 null constants, and "---", which
# the FITS header writes for a body it does not name. They still appear in `label_dict`
# and `header_dict`, but nowhere else.
_NULL_VALUES = ('UNK', 'N/A', 'NULL', '---')

# The ticks of the New Horizons spacecraft clock per second
_SCLK_TICKS = 50000

# The columns that identify the same image in an index and its supplement
_INDEX_KEYS = ('VOLUME_ID', 'FILE_SPECIFICATION_NAME')

# FITS header keywords that only carry text, left out of `header_dict` and `dict`
_COMMENTARY_KEYWORDS = ('', 'COMMENT', 'HISTORY')

# The PDS3 label parser's methods, in the order they are tried. "fast" rejects some valid
# labels, such as a quoted value that starts at the beginning of the next line; "strict"
# reads every valid label; "loose" also accepts some invalid ones, such as a keyword with
# no value.
_PDS3_METHODS = ('fast', 'strict', 'loose')

# The parser's other method, which is used only when it is asked for
_OTHER_PDS3_METHODS = ('compound',)

# The PDS3 names of the metadata items that both the detached label and the FITS header
# give, and the FITS keywords each is converted from. The detached label gives the lines
# of the image in its IMAGE object, which keeps its own name in `dict`; LINES is the
# item's name only inside the host.
_FITS_ITEM_KEYWORDS = {
    'START_TIME'                  : ('SPCSCET', 'EXPTIME'),
    'EXPOSURE_DURATION'           : ('EXPTIME',),
    'TARGET_NAME'                 : ('SPCCBTNM',),
    'SPACECRAFT_CLOCK_START_COUNT': ('SPCSCLK0',),
    'SPACECRAFT_CLOCK_STOP_COUNT' : ('STOPMET',),
    'LINES'                       : ('NAXIS2',),
}

# The precision of the detached label's START_TIME (milliseconds) and EXPOSURE_DURATION.
# Within half of it, the label's value is the FITS header's value rounded.
_TIME_DIGITS = 3
_LABEL_ROUNDING = 0.5e-3 + 1.e-6        # seconds, with room for floating-point error

##########################################################################################
# The resolution rules. Where the detached label and the FITS header both give an item
# and the values differ, the rule for that item picks one. Evidence: a survey of the
# label and header of 7,897 products, a random sample of all 14 volumes (2026-11-18),
# recorded with every other difference between the two in #13.
##########################################################################################


def _rounded_label_loses(label_value, fits_value):
    """The FITS header's value if the label's is the same number rounded; otherwise the
    label's."""

    if abs(label_value - fits_value) <= _LABEL_ROUNDING:
        return fits_value
    return label_value


def _label_wins(label_value, fits_value):
    """The detached label's value, the default rule."""

    return label_value


def _fits_wins(label_value, fits_value):
    """The FITS header's value."""

    return fits_value


def _is_spice_body(name):
    """True if the kernel pool knows `name` as the name of a body."""

    try:
        cspyce.bodn2c(name)
    except (KeyError, ValueError):
        return False
    return True


def _spice_body_wins(label_value, fits_value):
    """The FITS header's value if it names a body that SPICE knows and the label's does
    not; otherwise the label's."""

    if _is_spice_body(fits_value) and not _is_spice_body(label_value):
        return fits_value
    return label_value


_RESOLUTION_RULES = {
    # START_TIME: the label writes the FITS header's start (SPCSCET - EXPTIME/2) rounded
    # to the millisecond on all 7,897 products surveyed; the header's value is the
    # precise one (#13).
    'START_TIME': _rounded_label_loses,

    # EXPOSURE_DURATION: equal on 7,170 products; on the other 727 they differ only in
    # the last bit of the header's binary value, 9.9676 against 9.967600000000001 (#13).
    'EXPOSURE_DURATION': _rounded_label_loses,

    # TARGET_NAME: both give a name and they differ on 1,005 products, and neither source
    # is consistently the better one. The label writes a catalog name ("J1 IO", "J17
    # CALLIRRHOE") where the header writes the SPICE name ("IO") on 100; the header writes
    # a pointing ("RADEC=270.66,-14.63=P4") or a misspelling ("FJUPITER") where the label
    # names the body ("KERBEROS", "JUPITER") on 93; on 807 SPICE knows neither name. The
    # target body must be one that SPICE knows, so a name it knows wins; otherwise the
    # label's, the default (#13).
    'TARGET_NAME': _spice_body_wins,

    # The clock counts agree on all 7,897 products surveyed, so the default applies.
    'SPACECRAFT_CLOCK_START_COUNT': _label_wins,
    'SPACECRAFT_CLOCK_STOP_COUNT' : _label_wins,

    # LINES: they agree on all 7,897 products surveyed. The header describes the array
    # that `data` holds, which the field of view must match.
    'LINES': _fits_wins,
}

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

        What is read depends on the input:

        - A FITS file: the FITS file alone, with its pixels unless `astrometry` is True.
        - A PDS3 label: the label, and the FITS file that its ``^IMAGE`` pointer names
          (in the label's directory, in the case of the label's own name first), with
          its pixels. With `astrometry`, the label alone; the FITS file need not exist.

        The metadata come from whichever of the two was read: the detached label, and
        the label attached to the image, which is the FITS header. Each item of the FITS
        header is converted to the name and form that the PDS3 label gives it:

        - START_TIME, the start of the exposure: SPCSCET, the middle of the exposure,
          less half of EXPTIME.
        - EXPOSURE_DURATION: EXPTIME.
        - TARGET_NAME: SPCCBTNM.
        - SPACECRAFT_CLOCK_START_COUNT: SPCSCLK0.
        - SPACECRAFT_CLOCK_STOP_COUNT: STOPMET, a count of seconds, converted to the
          form "SSSSSSSSSS:TTTTT" (ten digits of seconds and five of the 1/50000-second
          ticks).
        - The lines of the image: NAXIS2, the LINES of the label's IMAGE object.

        A value of "UNK", "N/A", "NULL" or "---" counts as missing. Where only one source
        gives an item, its value is used. Where both give it and the values differ, a
        rule for each item chooses:

        - START_TIME and EXPOSURE_DURATION: the FITS header's if the label's is the
          same value rounded to the millisecond; otherwise the label's.
        - TARGET_NAME: the one that the SPICE kernel pool knows as a body, if only one
          is, as with the header's "CALLIRRHOE" rather than the label's "J17
          CALLIRRHOE", or the label's "KERBEROS" rather than the header's
          "RADEC=270.66,-14.63=P4"; otherwise the label's.
        - The lines of the image: the FITS header's, which describe the array in `data`.
        - The clock counts: the label's.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | HDUList): The path to the
                FITS file or its PDS3 label, or its parsed content as a Pds3Label or an
                HDUList, which is read as the label or the FITS file would be. The FITS
                file is closed once it has been read, including an HDUList that was
                passed in.
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
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. If the
                parser's method, by default "fast", cannot read the label, the methods
                after it in the order "fast", "strict", "loose" are tried in turn.
            astrometry (bool, optional):
                True to specify that the returned Observations will only be used for
                timing and astrometry. In this case, no data arrays are read and returned
                in the Observation, and a PDS3 label is read without its FITS file.
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
            Snapshot: The observation, with subfields `dict`, `host`, `instrument`,
            `binning_mode`, `spice_kernels`, `filepath`, `basename`, `spice_to_frame`,
            `spice_frame_name` and `spice_frame_id` inserted, and these where they apply:

            - `label_target`, `spacecraft_clock_start_count` and
              `spacecraft_clock_stop_count`: the merged TARGET_NAME and clock counts,
              absent where no source read gives them. The counts are the start and end
              of the exposure, spelled as the PDS3 label spells them, as in
              "0030713585:49000".
            - `label_dict`: the PDS3 label as a dictionary, if it was read.
            - `header_dict`: the primary FITS header without its blank, COMMENT and
              HISTORY cards, and `headers`, the list of the headers of every HDU, if the
              FITS file was read.
            - `data`: the primary image as stored, which in a raw product includes the
              columns beyond the field of view, unless `astrometry` is True.
            - `error` and `quality`: the error estimates and quality flags of a
              calibrated product, unless `astrometry` is True; absent for a raw product.
            - `abspath` and `image_url`, unless `astrometry` is True and the label is
              remote, so that the FITS file was never fetched (`image_url` is then
              still given).

            `dict` holds the merged items under their PDS3 names, START_TIME as an ISO
            string to the millisecond, and every other key of each source read under its
            own name; the FITS keywords converted to the items are left out, and so is
            any value that counts as missing. The target is the merged TARGET_NAME,
            repaired: "NONE" if it is missing or "CALIBRATION", which mean the image has
            no target, and "PLUTO" for a name beginning with "PLUTO". `filepath` and the
            other path subfields name the FITS file even when a label was given.

        Raises:
            FileNotFoundError: If the FITS file that a PDS3 label names does not exist,
                except in an astrometry read.
            oops.OopsValueError: If `geom`, `pointing` or `fov_type` is not one of the
                values listed above.
            oops.host.HostError: If `geom` or `pointing` is not "spice" in an astrometry
                read of a PDS3 label, which reads no FITS header, or if the pointing
                comes from the C kernels and they give none at mid-exposure; the message
                then ends with the SPICE error, which begins with its short code, such as
                "SPICE(NOFRAMECONNECT)".
            oops.host.HostFormatError: If the image is not 1024 or 256 lines, if the
                FITS file is not valid FITS, if a PDS3 label names no FITS file, if an
                HDUList passed in was not read from a file, or if no source read gives
                the start time, the exposure duration or the lines of the image.
            oops.host.PdsHostError: If no method of the PDS3 label parser can read the
                label.
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
                parsed = NewHorizonsLORRI._parse_label(filepath, pds3_options)
        else:
            filepath = None
            parsed = fileinfo
            if isinstance(parsed, Pds3Label) and not NewHorizonsLORRI._detect_in_pds3(
                    parsed):
                raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} PDS3 '
                                               f'label: {parsed.filepath}')
            if not isinstance(parsed, (Pds3Label, pyfits.HDUList)):
                raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} PDS3 '
                                               f'label or FITS file: {parsed!r}')

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

        label = None
        hdulist = None
        if isinstance(parsed, Pds3Label):
            label = parsed
            if astrometry:
                # An astrometry read of a label reads the label alone
                if geom != 'spice' or pointing != 'spice':
                    raise oops.host.HostError('geom and pointing from the FITS header '
                                              'need the FITS file; an astrometry read of '
                                              'a PDS3 label reads only the label')
                filepath = NewHorizonsLORRI._fits_path(label, check=False)
            else:
                filepath = NewHorizonsLORRI._fits_path(label)
                hdulist = NewHorizonsLORRI._open_fits(filepath)
        else:
            hdulist = parsed
            if filepath is None:
                # The path subfields name the file, so an HDUList must come from one
                if hdulist.filename() is None:
                    hdulist.close()
                    raise oops.host.HostFormatError('an HDUList without a file name '
                                                    'cannot name its image')
                filepath = FCPath(hdulist.filename())
                if not NewHorizonsLORRI._detect_in_fits(hdulist):
                    hdulist.close()
                    raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} '
                                                   f'FITS file: {filepath}')

        # Everything needed is copied out of the file, which is then closed, including
        # an HDUList that oops.Host.from_file() or the caller opened
        headers = None
        arrays = None
        if hdulist is not None:
            try:
                headers = [hdu.header.copy() for hdu in hdulist]
                if not astrometry:
                    arrays = NewHorizonsLORRI._read_arrays(hdulist)
            finally:
                hdulist.close()

        return NewHorizonsLORRI._make_snapshot(filepath, label=label, headers=headers,
                                               arrays=arrays, geom=geom,
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
                rows must name the same files as the index's, matched on VOLUME_ID and
                FILE_SPECIFICATION_NAME, in any order. Where both tables give a column,
                the supplement's value is kept, because the supplement is spelled as the
                product labels are. If not specified, the
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
            list[Snapshot]: One observation per row of the index, in the index's order,
            each with subfields `dict`, `host`, `instrument`, `binning_mode`,
            `spice_kernels`, `filepath`, `basename`, `spice_to_frame`, `spice_frame_name`
            and `spice_frame_id` inserted, and `label_target`,
            `spacecraft_clock_start_count` and `spacecraft_clock_stop_count` where the
            row gives TARGET_NAME and the clock counts other than "UNK", "N/A", "NULL"
            or "---". `dict` is the row, merged with its row of the supplement.
            `filepath` is the VOLUME_ID and the FILE_SPECIFICATION_NAME, with the label's
            suffix replaced by the FITS file's, joined by a slash. The target is the
            row's TARGET_NAME repaired as for :meth:`from_file`.

        Raises:
            FileNotFoundError: If the supplement does not exist.
            oops.OopsValueError: If `fov_type` is not "fast", "slow" or "flat".
            oops.host.HostFormatError: If a row's BINNING_MODE is not "1X1" or "4X4", or
                if a row gives no START_TIME or EXPOSURE_DURATION.
            oops.host.HostIndexError: If the supplement and the index do not list the
                same files, each once, or if `supplement` is not specified and the name
                of `filepath` contains no "_index".
            oops.host.WrongHostError: If `filepath` is not a New Horizons LORRI index.
            OSError: If `filepath`, or a `supplement` that is specified, does not exist or
                cannot be read.
        """

        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)
        NewHorizonsLORRI.initialize()

        # Read and filter the index table
        if row_dicts is None:
            row_dicts = NewHorizonsLORRI._read_index_and_supplement(
                filepath, supplement, pds3_options=pds3_options)
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
    def _parse_label(filepath, pds3_options):
        """A LORRI PDS3 label, parsed by the first method of the parser that reads it.

        Parameters:
            filepath (FCPath): The path to the label.
            pds3_options (dict[str, Any]): The options for the parser, including its
                method. If it is one of :data:`_PDS3_METHODS`, the methods after it are
                tried in turn if it fails; a method of :data:`_OTHER_PDS3_METHODS` is
                tried alone.

        Returns:
            Pds3Label: The parsed label.

        Raises:
            oops.OopsValueError: If the method in `pds3_options` is not one of
                :data:`_PDS3_METHODS` or :data:`_OTHER_PDS3_METHODS`, or if another
                option is not one the parser accepts.
            oops.host.PdsHostError: If no method can read the label; the message gives
                each method's error.
            oops.host.WrongHostError: If the label does not describe a New Horizons
                LORRI image.
            OSError: If the file does not exist or cannot be read.
        """

        first = pds3_options['method']
        if first in _PDS3_METHODS:
            methods = _PDS3_METHODS[_PDS3_METHODS.index(first):]
        elif first in _OTHER_PDS3_METHODS:
            methods = (first,)
        else:
            raise oops.OopsValueError(f'invalid PDS3 parser method {first!r}; use one '
                                      'of ' + ', '.join(repr(m) for m in _PDS3_METHODS
                                                        + _OTHER_PDS3_METHODS))

        errors = []
        for method in methods:
            try:
                label = Pds3Label(filepath, **(pds3_options | {'method': method}))
            except TypeError as err:
                raise oops.OopsValueError(f'invalid PDS3 parser option: {err}') from err
            # Besides PdsError, the fast method raises a ValueError for a value it cannot
            # convert, such as a time with an hour of 25 (JulianValidateFailure)
            except (PdsError, ValueError) as err:
                errors.append(f'method {method}: {err}')
                continue

            if not NewHorizonsLORRI._detect_in_pds3(label):
                raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} PDS3 '
                                               f'label: {filepath}')
            return label

        raise oops.host.PdsHostError(f'PDS3 label cannot be parsed: {filepath}; '
                                     + '; '.join(errors))

    @staticmethod
    def _fits_path(label, *, check=True):
        """The FITS file that a PDS3 label's ``^IMAGE`` pointer names.

        The label may spell the name in a different case from the file system. The
        archive names the FITS file in the same case as the label file, so the name is
        tried first in the case of the label's own name, then as written, then in lower
        and in upper case. On a file system that ignores case, the first spelling tried
        is the one returned.

        Parameters:
            label (Pds3Label): The parsed label.
            check (bool, optional): True to look for the file; False to return the first
                spelling without looking, so that no FITS file is needed.

        Returns:
            FCPath: The path to the FITS file, in the label's directory.

        Raises:
            FileNotFoundError: If `check` is True and no spelling of the name exists.
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
        if not check:
            return parent / first

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
            dict[str, numpy.ndarray]: The arrays keyed by subfield name: `data`, the
            primary image, and `error` and `quality`, the extensions of a calibrated
            product that hold them, each only if the file has it. Each array is a copy,
            so it outlives the file.
        """

        arrays = {'data': np.array(hdulist[0].data)}
        for hdu in hdulist[1:]:
            extname = str(hdu.header.get('EXTNAME', '')).strip().upper()
            if extname == _ERROR_EXTNAME:
                arrays['error'] = np.array(hdu.data)
            elif extname == _QUALITY_EXTNAME:
                arrays['quality'] = np.array(hdu.data)

        return arrays

    @staticmethod
    def _make_snapshot(filepath, *, label, headers, arrays, geom, pointing, fov_type,
                       overrides):
        """A Snapshot built from a PDS3 label, the headers and arrays of a FITS file, or
        both.

        Parameters:
            filepath (FCPath): The path to the FITS file.
            label (Pds3Label | None): The detached label; None if it was not read.
            headers (list[astropy.io.fits.Header] | None): The header of every HDU; None
                if the FITS file was not read.
            arrays (dict[str, numpy.ndarray] | None): The arrays as :meth:`_read_arrays`
                returns them; None to omit them.
            geom (str): "spice" or "fits", as for :meth:`from_file`.
            pointing (str): "spice", "fits", "fitsapp" or "fits90", as for
                :meth:`from_file`.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.host.HostFormatError: If the image is not 1024 or 256 lines, or if no
                source gives the start time, the exposure duration or the lines.
        """

        header = None if headers is None else headers[0]
        source = filepath if label is None else label.filepath

        label_dict = None if label is None else label.as_dict()
        label_items = {}
        if label_dict is not None:
            label_items = NewHorizonsLORRI._items_from_pds3(label_dict)
            if 'LINES' in label_dict.get('IMAGE', {}):
                label_items['LINES'] = label_dict['IMAGE']['LINES']

        header_dict = None
        fits_items = {}
        if header is not None:
            header_dict = {key: value for (key, value) in header.items()
                           if key not in _COMMENTARY_KEYWORDS}
            fits_items = NewHorizonsLORRI._items_from_header(header_dict)

        items = NewHorizonsLORRI._merge_items(label_items, fits_items)
        for name in ('EXPOSURE_DURATION', 'START_TIME', 'LINES'):
            NewHorizonsLORRI._require(items, name, source)

        lines = items['LINES']
        if lines not in _BINNING_MODES:
            raise oops.host.HostFormatError(f'unrecognized LORRI image shape: '
                                            f'{lines} lines in {source}')
        binning_mode = _BINNING_MODES[lines]

        tstart = items['START_TIME']
        texp = items['EXPOSURE_DURATION']
        tdb_midtime = tstart + texp / 2.
        target_name = NewHorizonsLORRI._repair_target(items.get('TARGET_NAME'))

        # Define the path and frame. A default that an override replaces is not built,
        # so the header need not give its values, and the apparent pointing uses the
        # observer that the Snapshot will have. geom and pointing other than "spice"
        # are refused before this point when the header was not read. The header's
        # position, velocity and pointing are those of its own mid-exposure, SPCSCET.
        if overrides['path'] is not None:
            default_path = overrides['path']
        elif geom == 'spice':
            default_path = _NewHorizons.PATH_ID
        else:
            default_path = NewHorizonsLORRI._fits_path_from_header(header, filepath,
                                                                   header['SPCSCET'])

        using_cks = overrides['frame'] is None and pointing == 'spice'
        if overrides['frame'] is not None:
            default_frame = overrides['frame']
        elif pointing == 'spice':
            NewHorizonsLORRI._define_camera_frame()
            check_time = (tdb_midtime if overrides['cadence'] is None
                          else overrides['cadence'].midtime)
            NewHorizonsLORRI._check_pointing(check_time, filepath)
            default_frame = _SPICE_FRAME_NAME
        else:
            default_frame = NewHorizonsLORRI._fits_frame_from_header(
                header, filepath, pointing=pointing, path=default_path,
                tdb_midtime=header['SPCSCET'])

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

        # Apply standard overrides. A remote FITS file that was not read was not
        # fetched, so no local copy of it exists to name.
        oops.Host._apply_overrides(params, filepath, **overrides)
        if header is None and not filepath.is_local():
            del params['abspath']

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if arrays is not None:
            obs.insert_subfields(arrays)

        # The merged items under their PDS3 names, then every other key of each source
        # without missing values. A null in one source never hides the other's value.
        dict_ = {}
        if label_dict is not None:
            dict_.update({key: value for (key, value) in label_dict.items()
                          if not NewHorizonsLORRI._is_missing(value)})
        if header_dict is not None:
            # NAXIS2 stays, because LINES is not a key of `dict`
            converted = {keyword for (name, keywords) in _FITS_ITEM_KEYWORDS.items()
                         if name != 'LINES' for keyword in keywords}
            dict_.update({key: value for (key, value) in header_dict.items()
                          if key not in converted and key not in dict_
                          and not NewHorizonsLORRI._is_missing(value)})
        dict_.update(NewHorizonsLORRI._dict_items(items))

        obs.insert_subfield('dict', dict_)
        if label_dict is not None:
            obs.insert_subfield('label_dict', label_dict)
        if header_dict is not None:
            obs.insert_subfield('header_dict', header_dict)
            obs.insert_subfield('headers', headers)
        NewHorizonsLORRI._insert_common_subfields(obs, items, binning_mode=binning_mode,
                                                  target=target_name, ck=using_cks)
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
            oops.host.HostFormatError: If the row's BINNING_MODE is not "1X1" or "4X4",
                or if the row gives no START_TIME or EXPOSURE_DURATION.
        """

        binning_mode = row_dict['BINNING_MODE'].upper()
        if binning_mode not in _BINNING_MODES.values():
            raise oops.host.HostFormatError(f'unrecognized LORRI BINNING_MODE: '
                                            f'{row_dict["BINNING_MODE"]!r}')

        # The index names the label; the observation names the FITS file beside it
        label_name = row_dict['FILE_SPECIFICATION_NAME']
        stem, _, suffix = label_name.rpartition('.')
        fits_suffix = 'FIT' if suffix.isupper() else 'fit'
        filepath = row_dict['VOLUME_ID'] + '/' + stem + '.' + fits_suffix

        items = NewHorizonsLORRI._items_from_pds3(row_dict)
        for name in ('EXPOSURE_DURATION', 'START_TIME'):
            NewHorizonsLORRI._require(items, name, filepath)

        tstart = items['START_TIME']
        texp = items['EXPOSURE_DURATION']

        using_cks = overrides['frame'] is None
        if using_cks:
            NewHorizonsLORRI._define_camera_frame()

        target_name = NewHorizonsLORRI._repair_target(items.get('TARGET_NAME'))
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
        NewHorizonsLORRI._insert_common_subfields(obs, items, binning_mode=binning_mode,
                                                  target=target_name, ck=using_cks)
        return obs

    @staticmethod
    def _check_pointing(tdb, filepath):
        """Raise HostError if the C kernels give no pointing of the camera at a time.

        Parameters:
            tdb (float): The time in seconds TDB.
            filepath (FCPath): The path to the FITS file, for the message.

        Raises:
            oops.host.HostError: If the camera frame cannot be evaluated at `tdb`. The
                message ends with the SPICE error, which begins with its short code,
                such as "SPICE(NOFRAMECONNECT)", so that missing SPICE data can be told
                from other failures.
        """

        frame = oops.Frame.as_frame(_SPICE_FRAME_NAME).wrt(oops.Frame.J2000)
        try:
            frame.transform_at_time(tdb)
        except Exception as err:    # cspyce maps each SPICE error to a builtin exception
            if 'SPICE(' not in str(err):
                raise
            raise oops.host.HostError(f'no {NewHorizonsLORRI.NAME} pointing for '
                                      f'{filepath}: {err}') from err

    @staticmethod
    def _require(items, name, source):
        """Raise HostFormatError if a metadata item the Snapshot needs is missing.

        Parameters:
            items (dict[str, Any]): The metadata items.
            name (str): The name of the item.
            source (str | FCPath): The file or index row the items came from, for the
                message.

        Raises:
            oops.host.HostFormatError: If `items` does not contain `name`.
        """

        if name not in items:
            description = {'LINES': 'lines of the image (IMAGE LINES or NAXIS2)'}
            raise oops.host.HostFormatError(f'no {description.get(name, name)} for '
                                            f'LORRI image {source}')

    @staticmethod
    def _is_missing(value):
        """True if a value means the item is missing.

        Parameters:
            value (Any): The value.

        Returns:
            bool: True if `value` is None, as for an undefined FITS card, or one of
            :data:`_NULL_VALUES`, ignoring case and surrounding blanks.
        """

        if value is None:
            return True
        return isinstance(value, str) and value.strip().upper() in _NULL_VALUES

    @staticmethod
    def _items_from_pds3(dict_):
        """The metadata items that a PDS3 label or an index row gives.

        Parameters:
            dict_ (dict[str, Any]): The label as a dictionary, or the index row.

        Returns:
            dict[str, Any]: The items, keyed by PDS3 name, leaving out any that is
            missing, counts as missing, or cannot be read: START_TIME, in seconds TDB;
            EXPOSURE_DURATION, in seconds; and TARGET_NAME,
            SPACECRAFT_CLOCK_START_COUNT and SPACECRAFT_CLOCK_STOP_COUNT, as Python
            strings.
        """

        values = {key: dict_[key] for key in _FITS_ITEM_KEYWORDS
                  if key in dict_ and not NewHorizonsLORRI._is_missing(dict_[key])}

        items = {}
        if 'START_TIME' in values:
            # A label gives START_TIME as a datetime, an index row as an ISO string. A
            # time or an exposure that does not convert counts as missing, so that the
            # other source's value is used if there is one.
            start_time = values['START_TIME']
            if isinstance(start_time, dt.datetime):
                start_time = start_time.isoformat()
            try:
                items['START_TIME'] = julian.tdb_from_tai(julian.tai_from_iso(
                    str(start_time)))
            except ValueError:
                pass
        if 'EXPOSURE_DURATION' in values:
            try:
                items['EXPOSURE_DURATION'] = float(values['EXPOSURE_DURATION'])
            except (TypeError, ValueError):
                pass

        # The table reader returns NumPy strings; the subfields are Python strings
        for key in ('TARGET_NAME', 'SPACECRAFT_CLOCK_START_COUNT',
                    'SPACECRAFT_CLOCK_STOP_COUNT'):
            if key in values:
                items[key] = str(values[key])

        return items

    @staticmethod
    def _items_from_header(header_dict):
        """The metadata items that a FITS header gives, converted to the PDS3 forms.

        Parameters:
            header_dict (dict[str, Any]): The primary header as a dictionary.

        Returns:
            dict[str, Any]: The items as for :meth:`_items_from_pds3`, converted from the
            keywords of :data:`_FITS_ITEM_KEYWORDS`, and LINES, the lines of the image.
            An item is left out if a keyword it needs is missing, undefined, or counts
            as missing.
        """

        def given(keyword):
            return (keyword in header_dict
                    and not NewHorizonsLORRI._is_missing(header_dict[keyword]))

        items = {}
        if given('EXPTIME'):
            items['EXPOSURE_DURATION'] = header_dict['EXPTIME']
            if given('SPCSCET'):
                items['START_TIME'] = header_dict['SPCSCET'] - header_dict['EXPTIME'] / 2.
        if given('SPCCBTNM'):
            items['TARGET_NAME'] = header_dict['SPCCBTNM']
        if given('SPCSCLK0'):
            items['SPACECRAFT_CLOCK_START_COUNT'] = header_dict['SPCSCLK0']
        if given('STOPMET'):
            items['SPACECRAFT_CLOCK_STOP_COUNT'] = NewHorizonsLORRI._sclk_from_met(
                header_dict['STOPMET'])
        if given('NAXIS2'):
            items['LINES'] = header_dict['NAXIS2']

        return items

    @staticmethod
    def _merge_items(label_items, fits_items):
        """The metadata items of the detached label and the FITS header, merged.

        Parameters:
            label_items (dict[str, Any]): The items of the label.
            fits_items (dict[str, Any]): The items of the FITS header.

        Returns:
            dict[str, Any]: Every item of either source. Where both give an item, the
            values are compared, and where they differ, the item's rule in
            :data:`_RESOLUTION_RULES` chooses.
        """

        items = fits_items | label_items
        for name in label_items.keys() & fits_items.keys():
            if label_items[name] != fits_items[name]:
                items[name] = _RESOLUTION_RULES[name](label_items[name],
                                                      fits_items[name])
        return items

    @staticmethod
    def _dict_items(items):
        """The merged items as `dict` gives them.

        Parameters:
            items (dict[str, Any]): The merged items, which include START_TIME.

        Returns:
            dict[str, Any]: The items under their PDS3 names, without LINES, which the
            label gives inside its IMAGE object; START_TIME is an ISO string to the
            millisecond, as the label writes it.
        """

        dict_ = {key: value for (key, value) in items.items() if key != 'LINES'}
        tai = julian.tai_from_tdb(dict_['START_TIME'])
        dict_['START_TIME'] = julian.iso_from_tai(tai, digits=_TIME_DIGITS)
        return dict_

    @staticmethod
    def _insert_common_subfields(obs, items, *, binning_mode, target, ck):
        """Insert the subfields that file and index observations share.

        Parameters:
            obs (Snapshot): The observation.
            items (dict[str, Any]): The metadata items; TARGET_NAME and the clock counts
                become subfields where they are given.
            binning_mode (str): "1X1" or "4X4".
            target (str): The name of the default target, whose kernels are included in
                `spice_kernels`.
            ck (bool): True if the pointing comes from the C kernels, so that they are
                included in `spice_kernels`.
        """

        obs.insert_subfield('host', 'New Horizons')
        obs.insert_subfield('instrument', 'LORRI')
        obs.insert_subfield('binning_mode', binning_mode)
        for (subfield, name) in [('label_target', 'TARGET_NAME'),
                                 ('spacecraft_clock_start_count',
                                  'SPACECRAFT_CLOCK_START_COUNT'),
                                 ('spacecraft_clock_stop_count',
                                  'SPACECRAFT_CLOCK_STOP_COUNT')]:
            if name in items:
                obs.insert_subfield(subfield, items[name])

        # With a custom frame, pointing never came from a CK, so any CK that happens to be
        # furnished is unrelated to this observation and must not be reported as used.
        kernels = _NewHorizons.used_kernels(obs.time, 'lorri', target=target, ck=ck)
        obs.insert_subfield('spice_kernels', kernels)
        obs.insert_subfield('spice_to_frame', _SPICE_TO_FRAME)
        obs.insert_subfield('spice_frame_name', _SPICE_FRAME_NAME)
        obs.insert_subfield('spice_frame_id', _SPICE_FRAME_ID)

    @staticmethod
    def _repair_target(name):
        """The name of the target body for a TARGET_NAME.

        A label or index gives "CALIBRATION" for the calibration images of the launch
        phase, whose FITS header gives no target. Some images of Pluto give a name that
        begins with "PLUTO".

        Parameters:
            name (str | None): The merged TARGET_NAME; None if no source gives one.

        Returns:
            str: The name of the target body: "NONE" for an image with no target, "PLUTO"
            for a name beginning with "PLUTO", and otherwise `name` itself.
        """

        if name is None or name.strip() in _NO_TARGET_NAMES:
            return 'NONE'
        if name.startswith('PLUTO'):
            return 'PLUTO'
        return name

    @staticmethod
    def _sclk_from_met(met):
        """A mission elapsed time converted to a clock count in the PDS3 label's form.

        Parameters:
            met (float): Seconds of the spacecraft clock, as the FITS header's STARTMET
                and STOPMET give them.

        Returns:
            str: The count as "SSSSSSSSSS:TTTTT", ten digits of whole seconds and five of
            the 1/50000-second ticks past them.
        """

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

        # The two-argument arctangent keeps the quadrant of a pixel beyond a pole
        ra  = oops.DPR * (np.arctan2(xi, denom) + ra0)
        dec = oops.DPR * np.arctan((eta * np.cos(dec0) + np.sin(dec0))
                                   / np.sqrt(denom**2 + xi**2))

        return (ra, dec)

    @staticmethod
    def _read_index_and_supplement(filepath, supplement, *, pds3_options):
        """The rows of an index, each merged with the row of its supplement for the same
        image.

        The rows are matched on :data:`_INDEX_KEYS` rather than by position, because the
        supplement of NHJULO_1001 lists two rows in a different order from its index.
        The other columns the two tables share are written differently: the index gives
        each clock count with its partition, as in "1/0030598438:48850", and rounds the
        stop count differently from the label, whereas the supplement copies the label.
        So the supplement's value of a shared column is kept.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The path to the label of the index.
            supplement (str | pathlib.Path | FCPath | object | None): The path to the
                label of the supplement; None or :data:`oops.Host.SUPPLEMENT` for the
                file beside the index whose name has "_supplemental_index" in place of
                "_index".
            pds3_options (dict[str, Any] | None): The options for the PDS3 label parser.

        Returns:
            list[dict[str, Any]]: One dictionary per row of the index, in its order.

        Raises:
            FileNotFoundError: If the supplement does not exist.
            oops.host.HostIndexError: If the two tables do not list the same images, each
                once, or if `supplement` is not specified and the name of `filepath`
                contains no "_index".
            oops.host.WrongHostError: If either table is not a New Horizons LORRI index.
        """

        filepath = FCPath(filepath)
        if supplement is None or supplement is oops.Host.SUPPLEMENT:
            name = (filepath.name.replace('_index', '_supplemental_index')
                                 .replace('_INDEX', '_SUPPLEMENTAL_INDEX'))
            if name == filepath.name:
                raise oops.host.HostIndexError(f'supplemental index required for '
                                               f'{filepath}')
            supplement = filepath.with_name(name)

        row_dicts = NewHorizonsLORRI._read_index_rows(filepath, supplement=None,
                                                      pds3_options=pds3_options)
        supp_dicts = NewHorizonsLORRI._read_index_rows(supplement, supplement=None,
                                                       pds3_options=pds3_options)

        def key(row_dict):
            return tuple(str(row_dict[name]) for name in _INDEX_KEYS)

        supp_by_key = {}
        for supp_dict in supp_dicts:
            if key(supp_dict) in supp_by_key:
                raise oops.host.HostIndexError(f'supplemental index lists '
                                               f'{key(supp_dict)} twice: {supplement}')
            supp_by_key[key(supp_dict)] = supp_dict

        index_keys = {key(row_dict) for row_dict in row_dicts}
        if len(index_keys) != len(row_dicts):
            raise oops.host.HostIndexError(f'index lists an image twice: {filepath}')
        if index_keys != supp_by_key.keys():
            missing = sorted(index_keys - supp_by_key.keys())
            extra = sorted(supp_by_key.keys() - index_keys)
            raise oops.host.HostIndexError(
                f'index files list different images: {filepath}, {supplement}; '
                f'only in the index: {missing[:3]}; only in the supplement: {extra[:3]}')

        return [row_dict | supp_by_key[key(row_dict)] for row_dict in row_dicts]

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
