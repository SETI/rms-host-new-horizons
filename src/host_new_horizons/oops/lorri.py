##########################################################################################
# host_new_horizons/oops/lorri.py
##########################################################################################
"""The New Horizons LORRI host, which makes a Snapshot from each LORRI image.

The :doc:`user guide </user_guide/user_guide>` describes the products it reads, how it
merges the two labels of each product, and the errors it raises.
"""

import collections
import datetime as dt
import re
import warnings

import numpy as np

import julian

from astropy.io               import fits as pyfits
from astropy.utils.exceptions import AstropyUserWarning
from filecache                import FCPath
from pdsparser                import Pds3Label

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

# The C kernel IDs of the pointing: the spacecraft bus, whose attitude the C kernels give,
# and its clock. The camera frame turns with the spacecraft, so it is never frozen.
_CK_OBJECT_ID = -98000
_CK_SCLK_ID = -98
_SPICE_FRAME_IS_FROZEN = False

# The dimensions of an unbinned LORRI image in pixels
_LORRI_PIXELS = 1024

# The binning mode of an image, keyed by its number of lines
_BINNING_MODES = {1024: '1X1', 256: '4X4'}

# The pixel scale of the flat FOV in radians, the IFOV of the instrument kernel rounded
_FLAT_SCALE = 4.96e-6

_GEOMS = ('spice', 'fits')
_POINTINGS = ('spice', 'fits', 'fitsapp', 'fits90')
_FOV_TYPES = ('fast', 'slow', 'flat')

# The parser method of the PDS3 labels. Every current LORRI product label parses under
# strict as it is written, so the host corrects no label text.
_DEFAULT_PDS3_OPTIONS = {'method': 'strict'}

# Before 2013, the FITS headers of the first PDS delivery give a north clock angle 90
# degrees short; pointing="fits90" adds it back for an image taken in or before this year
_FITS90_LAST_YEAR = 2012

# The names of the extensions of a calibrated product that hold its error estimates and
# quality flags, compared in upper case. Raw products have other extensions instead.
_ERROR_EXTNAME = 'LORRI ERROR IMAGE'
_QUALITY_EXTNAME = 'LORRI QUALITY FLAG IMAGE'

# The text with which astropy warns of a FITS file that ends early
_TRUNCATED_WARNING = 'File may have been truncated'

# The target name of the calibration images of the launch phase in the label and index
# (NHLALO_2001), which means the image has no target. Their FITS header gives "---", a
# null value.
_NO_TARGET_NAMES = ('CALIBRATION',)

# The target body of a target name is the name that SPICE knows the body by, translated
# from the label's form of the name by fixed rules (_target_body()). A moon of Jupiter is
# "J<n> <NAME>", as in "J1 IO" and "J17 CALLIRRHOE"; a numbered minor planet with a name
# is "<number> <NAME>", optionally followed by its provisional designation, as in "136472
# MAKEMAKE" and "50000 QUAOAR (2002 LM60)"; a numbered minor planet without a name is
# "ASTEROID <number> (<designation>)". The name of a named minor planet needs at least
# three letters, so that a provisional designation such as "2014 MU" is not read as one.
_JUPITER_MOON_NAME = re.compile(r'J\d+ (.+)')
_NAMED_MINOR_PLANET = re.compile(r'\d+ ([A-Z][A-Z\' -]+[A-Z])(?: \(.+\))?')
_NUMBERED_MINOR_PLANET = re.compile(r'ASTEROID (\d+) \(.+\)')

# The names that SPICE knows numbered minor planets by, where the label gives only their
# number and designation, keyed by number
_MINOR_PLANET_NAMES = {
    '486958': 'ARROKOTH',       # "ASTEROID 486958 (2014 MU69)", NAIF ID 2486958 (NHKELO)
}

# The bodies whose name in the FITS header of an image is one SPICE does not know, and
# the name it knows. The header names Arrokoth "MU69" on every image of it (1,078 in a
# random sample of 7,195 products), so an image read without its label still has a
# target body.
_HEADER_TARGET_NAMES = {
    'MU69': 'ARROKOTH',
}

# The number that the FITS header writes for a value its pipeline cannot compute, such as
# the target geometry of an image without a target (SPCTSORN, SPCTRANG, ...); the label
# gives "N/A" for the same values. It counts as missing, and stays in `header_dict` and
# `headers`.
_FITS_FILL = -1.e32

# The ticks of the New Horizons spacecraft clock per second
_SCLK_TICKS = 50000

# The solar flux at 1 AU at the pivot wavelength of LORRI, 6076.2 Angstroms, in
# erg/cm^2/s/Angstrom, and the astronomical unit in km. The FITS header of a calibrated
# product gives RSOLAR, the conversion from DN/s/pixel to the radiance of a solar
# spectrum in erg/cm^2/s/Angstrom/sr, from which I/F = pi r^2 DN / (texp RSOLAR F) for a
# target r AU from the Sun.
_SOLAR_FLUX = 176.
_AU = 149597870.7

# FITS header keywords that only carry text, left out of `dict`
_COMMENTARY_KEYWORDS = ('', 'COMMENT', 'HISTORY')

# The FITS keywords that are converted to the PDS3 items START_TIME, EXPOSURE_DURATION,
# TARGET_NAME and the clock counts, and so appear in `dict` only under those names. NAXIS2
# gives the lines of the image, which the label gives inside its IMAGE object, so it
# keeps its own name in `dict`.
_CONVERTED_KEYWORDS = ('SPCSCET', 'EXPTIME', 'SPCCBTNM', 'SPCSCLK0', 'STOPMET')

##########################################################################################
# The resolution rules. Where the detached label and the FITS header both give an item
# and the values differ, the item's winning source gives it, the same source for every
# image: "label" for the detached label, "header" for the FITS header. The winner is the
# source that is accurate for the item. A value missing from the winning source is not a
# difference; the other source's value is used. The evidence comes from comparisons of
# the label and the header of 7,897 and 7,195 products, random samples of all 14
# volumes; every difference between the two is listed in SETI/rms-host-new-horizons#13.
##########################################################################################

_WINNING_SOURCE = {
    # START_TIME: the header's, SPCSCET - EXPTIME/2. On 2,342 products, it equals the
    # start that the label's own SPACECRAFT_CLOCK_START_COUNT gives through the spacecraft
    # clock kernel to 0.0001 ms; the label's START_TIME is that time rounded to the
    # millisecond, from 0.5 ms early to 0.5 ms late on every product. On frames taken
    # while the spacecraft spun, in cruise as well as at launch, the rounding moves the
    # pointing at mid-exposure by as much as 48.7 pixels (lor_0055952951_0x630_eng).
    'START_TIME': 'header',

    # EXPOSURE_DURATION: the header's. Both are accurate, matching the difference of the
    # clock counts to 0.0008 ms, and they differ only in the last binary digit (238 of
    # 7,195 products, all in NHKELO). The header's keeps START_TIME plus half of the
    # exposure equal to SPCSCET.
    'EXPOSURE_DURATION': 'header',

    # TARGET_NAME: the label's. The header's SPCCBTNM is what the label's NOTE calls the
    # data pipeline's best guess at the target, outside the peer review: "FJUPITER" for
    # JUPITER (64 of 7,195 products), "ANSA" for J RINGS (1,035), pointings such as
    # "RADEC=270.66,-14.63=P4" for KERBEROS (46), and "---" for CALIBRATION. Where the
    # label's name is not the one SPICE knows ("J1 IO", "ASTEROID 486958 (2014 MU69)"),
    # it is still accurate; the target body comes from it by _target_body().
    'TARGET_NAME': 'label',

    # The clock counts: the label's. They agree on every product compared, and all 43,441
    # current labels give exactly the stop count that STOPMET gives.
    'SPACECRAFT_CLOCK_START_COUNT': 'label',
    'SPACECRAFT_CLOCK_STOP_COUNT' : 'label',

    # LINES: the header's. They agree on every product compared, and the header describes
    # the array that `data` holds, which the field of view must match.
    'LINES': 'header',
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

# What the one read function returns: the paths of the image file and of the label (None
# if no label was read), the label as a dictionary and the FITS headers and arrays (each
# None if not read), the merged metadata items, the start and duration of the exposure,
# the FOV's binning mode, `dict`, and the geometry and pointing from the FITS header
# (None unless requested)
_LorriFiles = collections.namedtuple('_LorriFiles',
                                     ['image_path', 'label_path', 'label_dict', 'headers',
                                      'header_dict', 'arrays', 'items', 'tstart', 'texp',
                                      'binning_mode', 'dict_', 'sun_state', 'radec'])


class NewHorizonsLORRI(oops.Host):
    """The New Horizons LORRI host: Snapshot constructors for Long Range Reconnaissance
    Imager images."""

    NAME = 'New Horizons LORRI'

    # The rows of an index and its supplement are paired by file, on the default key
    # (VOLUME_ID, FILE_SPECIFICATION_NAME). Each of the 16 LORRI metadata indexes lists
    # the same files as its supplement, and 14 list them in the same order. In
    # NHJULO_1001 the index lists lor_0035089864_0x632_eng out of place, last, where the
    # supplement lists it at row 679 (counting from 1); in the cumulative NHxxLO_1999 the
    # two differ from row 1,925 on for the same reason. The cumulative
    # NHxxLO_1999 and NHxxLO_2999 list four products twice under one
    # FILE_SPECIFICATION_NAME, once for NHKELO and once for NHK2LO, which the VOLUME_ID
    # of the key tells apart. SETI/rms-data-projects#262 item 1
    #: Each row of an index is paired with the row of its supplement that gives the same
    #: VOLUME_ID and FILE_SPECIFICATION_NAME, wherever that row is.
    MERGE_SUPPLEMENT_BY_FILE = True

    _FOVS = {}
    _initialized = False

    @staticmethod
    def from_file(fileinfo, *, geom='spice', pointing='spice', fov_type='fast',
                  pds3_options=None, astrometry=False, select=None, target=None,
                  lightsource=None, path=None, frame=None, cadence=None, fov=None,
                  calibrations=None, timeshift=None, navigation=None, tracker=None,
                  parallel=None):
        """A Snapshot object based on a given New Horizons LORRI image file.

        What is read depends on the input:

        - A FITS file: the FITS file alone, with its arrays unless `astrometry` is True.
        - A PDS3 label: the label, and the FITS file that its ``^IMAGE`` pointer names,
          with its arrays. The FITS file is in the label's directory, and its name is
          written in the case of the label's own name, upper or lower. With `astrometry`,
          the label alone; the FITS file need not exist.

        The metadata come from whichever of the two was read: the detached label, and
        the label attached to the image, which is the FITS header. These FITS keywords
        are converted to the name and form that the PDS3 label gives the item:

        - START_TIME, the start of the exposure: SPCSCET, the middle of the exposure,
          less half of EXPTIME.
        - EXPOSURE_DURATION: EXPTIME.
        - TARGET_NAME: SPCCBTNM.
        - SPACECRAFT_CLOCK_START_COUNT: SPCSCLK0.
        - SPACECRAFT_CLOCK_STOP_COUNT: STOPMET, a count of seconds, converted to the
          form "SSSSSSSSSS:TTTTT" (ten digits of seconds and five of the 1/50000-second
          ticks, rounded to the nearest tick).
        - The lines of the image: NAXIS2, the LINES of the label's IMAGE object.

        A value counts as missing if it is "UNK", "N/A", "NULL", "UNKNOWN", "NOT
        APPLICABLE", "---" or blank, ignoring case and surrounding blanks, the FITS
        header's fill value -1e32, a FITS card without a value, or a list or set of such
        values. Where only one source gives an item, its value is used. Where both give
        it, each item takes the value of one source, the same for every image:

        - START_TIME, EXPOSURE_DURATION and the lines of the image: the FITS header's.
          The label gives the start time rounded to the millisecond.
        - TARGET_NAME and the clock counts: the label's. The FITS header's target name
          is the data pipeline's guess, such as "FJUPITER" for JUPITER, or a pointing
          such as "RADEC=270.66,-14.63=P4" for KERBEROS.

        The FOV follows the lines of the image: 1024 lines for the full-resolution mode
        and 256 for the 4x4 binned mode.

        The target body is the name that SPICE knows the body by, translated from
        TARGET_NAME: "J17 CALLIRRHOE" gives CALLIRRHOE, "136472 MAKEMAKE" gives
        MAKEMAKE, "ASTEROID 486958 (2014 MU69)" and the FITS header's "MU69" give
        ARROKOTH, a name beginning with "PLUTO" gives PLUTO, and a missing name or
        "CALIBRATION", which mean the image has no target, give "NONE". Any other name
        is the target body as it is.

        The light source is the Sun.

        The calibration is a solar-normalized I/F: one
        :class:`~oops.calibration.FlatCalib` named "I/F" for an extended source, whose
        factor is pi r^2 / (texp RSOLAR F), in which texp is the exposure duration,
        RSOLAR is the FITS header's solar spectral response, the conversion from DN/s
        per pixel to the radiance of a solar spectrum, F = 176 erg/cm^2/s/Angstrom is
        the solar flux at 1 AU at the pivot wavelength of LORRI, and r is the distance
        of the target (after the `target` override) from the Sun at mid-exposure (the
        middle of the `cadence` override, if it is given) in AU, from the SP kernels.
        The responses to the other spectra that a calibrated header gives (RPLUTO,
        RCHARON, RJUPITER, RPHOLUS, RMU69) and the P keywords of point sources are not
        used; they stay in `header_dict`. There is no calibration if the header gives no
        RSOLAR, as for a raw product or a read of the label alone, if the exposure is
        zero, or if the target has no ephemeris: it has no registered path, or the SP
        kernels do not give its position at that time.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | HDUList): The path to the FITS file
                or its PDS3 label, or the FITS file opened as an HDUList, which is read
                as the file it was opened from. The FITS file is closed once it has
                been read, including an HDUList that was passed in.
            geom (str, optional): "spice" to take the spacecraft's position from the SP
                kernels; "fits" to take it from the FITS header, as a linear path
                relative to the Sun through the position and velocity that SPCSSCX,
                SPCSSCY, SPCSSCZ, SPCSSCVX, SPCSSCVY and SPCSSCVZ give at SPCSCET.
            pointing (str, optional): "spice" to take the pointing from the C kernels;
                "fits" to take it from the WCS parameters and the north clock angle
                SPCEMEN of the FITS header; "fitsapp" to treat the pointing in the FITS
                header as apparent rather than geometric; "fits90" to use the FITS
                header's pointing, adding 90 degrees to the north clock angle of an
                image taken before 2013, which corrects an error in the first PDS
                delivery.
            fov_type (str, optional): "fast" to model the camera distortion with a
                polynomial and its numerical inverse; "slow" to invert the polynomial by
                Newton's method; "flat" for a field of view without distortion.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. The
                method is "strict" unless they give another.
            astrometry (bool, optional): True to specify that the returned Observation
                will only be used for timing and astrometry. In this case, no data
                arrays are read, and a PDS3 label is read without its FITS file.
            select (int | slice | tuple[int | slice, ...], optional): An index, slice,
                or tuple of them that chooses among the observations of the file, which
                holds one; None, the default, returns that observation itself. `0` also
                returns it, and a slice returns a tuple.
            target (str | oops.Body, optional): Override for the default target Body of
                the observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default
                LightSource of the observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument, such as a :class:`~oops.frame.Cmatrix`, used as it is given.
            cadence (oops.Cadence, optional): Override for the Cadence of the
                observation.
            fov (oops.FOV, optional): Override for the default FOV of the observing
                instrument.
            calibrations (oops.Calibration | list[oops.Calibration], optional): Override
                for the calibration or list of calibrations; an empty list for none.
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

        Returns:
            Snapshot | tuple[Snapshot, ...]: The observation, or the tuple that a slice in
            `select` chooses, with the subfields `target`, `lightsource`, `calibrations`,
            `dict`, `spice_to_frame` (diag(1,-1,-1), from the SPICE frame to the camera
            frame), `spice_frame_name` ("NH_LORRI"), `spice_frame_id` (-98300),
            `ck_object_id` (-98000), `ck_sclk_id` (-98) and `spice_frame_is_frozen`
            (False), and these where they apply:

            - `filepath`, `basename` and `image_url`: the FITS file, as a path, its name,
              and its absolute path or URL, even when a label was given.
            - `image_abspath`: the FITS file on the local disk, if it was read, which an
              astrometry read of a label leaves out.
            - `label_url` and `label_abspath`: the label, as an absolute path or URL and
              on the local disk, if a label was given.
            - `spacecraft_clock_start_count` and `spacecraft_clock_stop_count`: the
              merged clock counts of the start and end of the exposure, as in
              "0030713585:49000". Each is absent where no source read gives it.
            - `label_dict`: the PDS3 label as the parser's dictionary, with the keys the
              parser adds, if the label was read.
            - `header_dict`: the primary FITS header as a dictionary, and `headers`, the
              astropy Header of every HDU of the file, if the FITS file was read.
            - `data`: the primary image as stored, which in a raw product includes the
              columns beyond the field of view, unless `astrometry` is True.
            - `error` and `quality`: the error estimates and quality flags of a
              calibrated product, from the extensions named "LORRI Error image" and
              "LORRI Quality flag image", unless `astrometry` is True; absent for a raw
              product, which has no such extensions.

            `dict` holds every key of each source read under its own name, the label's
            where both give a key, without the values that count as missing, and the
            merged items under their PDS3 names, with START_TIME as a
            :class:`datetime.datetime` in UTC, as the label gives it. The FITS keywords
            converted to those items, and the blank, COMMENT and HISTORY cards, are
            left out.

        Raises:
            FileNotFoundError: If the file does not exist, or if the FITS file that a
                PDS3 label names does not exist, except in an astrometry read.
            oops.OopsValueError: If `geom`, `pointing` or `fov_type` is not one of the
                values listed above; if `geom` or `pointing` is not "spice" in an
                astrometry read of a PDS3 label, which leaves the FITS header unread; or
                if an option in `pds3_options` is not one the parser accepts, or the
                parser cannot convert a value in the label.
            oops.OopsIndexError: If `select` is an index other than 0 or -1.
            oops.host.HostFormatError: If the FITS file is not valid FITS; if the label
                or the FITS file exists but cannot be read, as when access to it is
                denied; or if either cannot be interpreted: the file ends early, the
                primary array is not one image band, the array does not have the lines
                and samples that the label's IMAGE object gives, or the files lack an
                item the observation needs, such as the exposure duration, the start
                time, or lines of an image size that LORRI takes. The message names the
                file that was passed.
            oops.host.PdsHostError: If the PDS3 label cannot be parsed.
            oops.host.WrongHostError: If `fileinfo` is not a New Horizons LORRI PDS3
                label or FITS file.
            TypeError: If a keyword argument is not one of those listed.
        """

        NewHorizonsLORRI._check_option('geom', geom, _GEOMS)
        NewHorizonsLORRI._check_option('pointing', pointing, _POINTINGS)
        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)

        NewHorizonsLORRI.initialize()

        files = NewHorizonsLORRI._read(fileinfo, astrometry=astrometry, geom=geom,
                                       pointing=pointing, pds3_options=pds3_options)

        # The IDs of a path and frame from the header are made unique by the file's name
        id_name = files.image_path.name
        if geom == 'spice':
            default_path = _NewHorizons.PATH_ID
        else:
            default_path = NewHorizonsLORRI._path_from_sun_state(files.sun_state,
                                                                 id_name)

        if pointing == 'spice':
            default_frame = _SPICE_FRAME_NAME
        else:
            default_frame = NewHorizonsLORRI._frame_from_radec(files.radec, id_name,
                                                               pointing=pointing,
                                                               path=default_path)

        target_name = NewHorizonsLORRI._target_body(files.items.get('TARGET_NAME'))
        params = {
            'cadence'     : oops.cadence.SnapCadence(files.tstart, files.texp),
            'fov'         : NewHorizonsLORRI._FOVS[files.binning_mode, fov_type],
            'path'        : default_path,
            'frame'       : default_frame,
            'target'      : target_name,
            'lightsource' : 'SUN',
            'calibrations': [],
        }
        oops.Host._apply_overrides(params, target=target, lightsource=lightsource,
                                   path=path, frame=frame, cadence=cadence, fov=fov,
                                   calibrations=calibrations, timeshift=timeshift,
                                   navigation=navigation, tracker=tracker,
                                   parallel=parallel)

        # The I/F calibration, for the target and at the midtime that the overrides leave
        if calibrations is None and files.header_dict is not None:
            params['calibrations'] = NewHorizonsLORRI._iof_calibrations(
                files.header_dict, files.texp, target=params['target'],
                tdb=params['cadence'].midtime)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if files.arrays is not None:
            obs.insert_subfields(files.arrays)

        obs.insert_subfield('dict', files.dict_)
        if files.label_dict is not None:
            obs.insert_subfield('label_dict', files.label_dict)
        if files.header_dict is not None:
            obs.insert_subfield('header_dict', files.header_dict)
            obs.insert_subfield('headers', files.headers)
        NewHorizonsLORRI._insert_common_subfields(obs, files.items)
        oops.Host._set_path_attributes(obs, files.image_path, files.label_path)

        return NewHorizonsLORRI._selected(obs, select)

    @staticmethod
    def from_index(filepath, *, supplement=oops.Host.SUPPLEMENT, label=None,
                   row_dicts=None, pds3_options=None, filter=None, fov_type='fast',
                   select=None, target=None, lightsource=None, path=None, frame=None,
                   cadence=None, fov=None, calibrations=None, timeshift=None,
                   navigation=None, tracker=None, parallel=None):
        """An iterator over Snapshots, one for each row of a New Horizons LORRI index.

        It reads a metadata index of the Ring-Moon Systems Node, such as
        ``NHJULO_2001_index.lbl``, and its supplement,
        ``NHJULO_2001_supplemental_index.lbl``, which gives the binning mode and the
        exposure duration that the index lacks. When this method is called, it reads
        only the labels of the two tables. The iteration then reads the tables 100 rows
        at a time and builds each row's Snapshot when it reaches that row.

        Each index row is paired with the supplement row that gives the same VOLUME_ID
        and FILE_SPECIFICATION_NAME (:attr:`MERGE_SUPPLEMENT_BY_FILE` is True), wherever
        that row is, because an index can list a file out of order. A supplement row
        read ahead of its partner is held until the partner is reached, so the rows
        held grow with how far the two tables are out of step. No other column of the
        two rows is compared. Where both tables give a column, the merged row has the
        index's value.

        Every row of the index gives an observation, in the index's order, with the
        values of the row as the tables give them:

        - The start of the exposure is START_TIME, a UTC time, and the exposure is
          EXPOSURE_DURATION. An exposure of zero gives an observation with an exposure
          of zero.
        - The FOV is that of BINNING_MODE, "1X1" or "4X4" in either case.
        - The target body comes from TARGET_NAME as for :meth:`from_file`.
        - The pointing comes from the C kernels and the position from the SP kernels.
          A row whose pointing the C kernels cannot give still gives an observation;
          its frame raises the SPICE error when it is used.
        - There is no calibration, because the tables do not give RSOLAR.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The path to the label of a New
                Horizons LORRI metadata index.
            supplement (str | pathlib.Path | FCPath | object | None, optional): The path
                to the PDS label of the supplemental index, in any directory. The
                default, `oops.Host.SUPPLEMENT`, reads the file in the directory of
                `filepath` whose name matches except that "_index" is replaced by
                "_supplemental_index" in upper or lower case. The index needs its
                supplement, so None is refused unless `row_dicts` is given.
            label (Pds3Label, optional): The PDS label of `filepath` if it has already
                been read.
            row_dicts (collections.abc.Iterable[dict[str, Any]], optional): The rows of
                the index, each already merged with its supplement row. If they are
                given, the tables are not read, and `supplement` and `label` are
                ignored.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. The
                index's label is parsed with all of them, and the supplement's label
                with the method alone, "strict" unless they give another; for example,
                ``pds3_options={'method': 'fast'}`` parses both with "fast".
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, bounds
                included, are yielded.
            fov_type (str, optional): "fast", "slow" or "flat", as for :meth:`from_file`.
            select (int | slice | tuple[int | slice, ...], optional): Applied to the
                observation of each row as for :meth:`from_file`.
            target (str | oops.Body, optional): Override for the default target Body of
                every observation. Use "NONE" for inertial pointing, indicating that no
                Solar System body was tracked.
            lightsource (oops.LightSource, optional): Override for the default
                LightSource of every observation.
            path (str | oops.Path, optional): Override for the Path of the observer.
            frame (str | oops.Frame, optional): Override for the Frame of the observing
                instrument, such as a :class:`~oops.frame.Cmatrix`, used as it is given.
            cadence (oops.Cadence, optional): Override for the Cadence of every
                observation.
            fov (oops.FOV, optional): Override for the FOV of every observation.
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

        Returns:
            collections.abc.Iterator[Snapshot | tuple[Snapshot, ...]]: One observation
            per filtered row of the index, in the index's order, or the tuple that a
            slice in `select` chooses from it. Each has the subfields `target`,
            `lightsource`, `calibrations`, `dict`, `filepath`, `basename`,
            `spice_to_frame`, `spice_frame_name`, `spice_frame_id`, `ck_object_id`,
            `ck_sclk_id` and `spice_frame_is_frozen`, as for :meth:`from_file`, and
            `spacecraft_clock_start_count` and `spacecraft_clock_stop_count` as the
            row gives them, with the index's partition, as in "1/0030713585:49000".
            `dict` is the merged row exactly as the tables give it, declared null values
            included, with the "_mask" entry of each column. `filepath` is the row's
            VOLUME_ID and FILE_SPECIFICATION_NAME joined by "/", as the index gives
            them, which names the product's label; `basename` is its last component. A
            caller that wants the file joins `filepath` with its own root directory. A
            caller that needs ``len()`` or indexing wraps the call in ``list()``.

        Raises:
            FileNotFoundError: When this method is called, if the label of `filepath`
                or of the supplement does not exist; when the iteration starts, if a
                table file does not exist.
            oops.OopsValueError: When this method is called, if `fov_type` is not
                "fast", "slow" or "flat", or if a value in `pds3_options` is invalid.
            oops.host.HostIndexError: When this method is called, if `supplement` is
                None while `row_dicts` is not given, if `supplement` is
                `oops.Host.SUPPLEMENT` and the name of `filepath` lacks "_index", or if
                a label does not describe a table. During the iteration, at the row
                concerned, if the tables do not give the same files, or if a table file
                holds a different number of rows than its label says.
            oops.host.PdsHostError: When this method is called, if the PDS label of
                `filepath` or `supplement` cannot be parsed.
            oops.host.WrongHostError: When this method is called, if `filepath` is not a
                New Horizons LORRI index.
            KeyError: During the iteration, at a row that lacks a column the
                observation needs, or a column that `filter` names.
            TypeError: When this method is called, if a keyword argument is not one of
                those listed.
        """

        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)

        if supplement is None and row_dicts is None:
            raise oops.host.HostIndexError(f'{NewHorizonsLORRI.NAME} index needs its '
                                           f'supplemental index: {filepath}')

        # Only the file key is compared (validate=()): the other columns that the two
        # tables share are written differently, so that a check of them fails on every
        # row of 10 of the 16 metadata indexes. The index gives each clock count with its
        # partition, as in "1/0030598438:48850", rounds the stop count differently, gives
        # PRODUCT_CREATION_TIME as a date without a time, and shortens some target names,
        # as in "50000 QUAOAR (2002 LM)". SETI/rms-data-projects#262 items 2 and 3
        rows = NewHorizonsLORRI._index_rows(filepath, supplement=supplement, label=label,
                                            row_dicts=row_dicts,
                                            pds3_options=pds3_options, filter=filter,
                                            validate=())

        NewHorizonsLORRI.initialize()

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

        def snapshot(row):
            """The observation of one row, as `select` chooses it."""
            obs = NewHorizonsLORRI._snapshot_from_row(row, fov_type=fov_type,
                                                      overrides=overrides)
            return NewHorizonsLORRI._selected(obs, select)

        # Build each row's Snapshot when the iteration reaches it
        return (snapshot(row) for row in rows)

    ######################################################################################
    # Reading the files
    ######################################################################################

    @staticmethod
    def _read(fileinfo, *, astrometry, geom, pointing, pds3_options):
        """Read and interpret the label, the FITS file, or both: the one function
        through which this host reads its files.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | HDUList): The file, as for
                :meth:`from_file`.
            astrometry (bool): True to read no arrays, and a PDS3 label without its FITS
                file.
            geom (str): "spice" or "fits", as for :meth:`from_file`.
            pointing (str): "spice", "fits", "fitsapp" or "fits90", as for
                :meth:`from_file`.
            pds3_options (dict[str, Any] | None): The caller's PDS3 parser options.

        Returns:
            _LorriFiles: Everything the observation is built from.

        Raises:
            oops.host.HostFormatError: For any failure to read or interpret the files
                other than those listed here, naming the file that was passed.
            FileNotFoundError: If a file does not exist.
            oops.OopsValueError: As described for :meth:`from_file`.
            oops.host.HostError: As raised by oops for a file that cannot be parsed or
                is not this host's.
        """

        name = (fileinfo.filename() if isinstance(fileinfo, pyfits.HDUList)
                else fileinfo)
        try:
            return NewHorizonsLORRI._read_and_interpret(fileinfo, astrometry=astrometry,
                                                        geom=geom, pointing=pointing,
                                                        pds3_options=pds3_options)
        except (FileNotFoundError, oops.host.HostError, oops.OopsValueError):
            raise
        except (KeyError, IndexError, TypeError, ValueError, OSError,
                AstropyUserWarning) as err:
            raise oops.host.HostFormatError(f'cannot read the {NewHorizonsLORRI.NAME} '
                                            f'file {name}: {type(err).__name__}: '
                                            f'{err}') from err

    @staticmethod
    def _read_and_interpret(fileinfo, *, astrometry, geom, pointing, pds3_options):
        """The work of :meth:`_read`, which converts the errors it raises.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | HDUList): The file, as for
                :meth:`from_file`.
            astrometry (bool): True to read no arrays, and a PDS3 label without its FITS
                file.
            geom (str): "spice" or "fits", as for :meth:`from_file`.
            pointing (str): "spice", "fits", "fitsapp" or "fits90", as for
                :meth:`from_file`.
            pds3_options (dict[str, Any] | None): The caller's PDS3 parser options.

        Returns:
            _LorriFiles: Everything the observation is built from.
        """

        options = _DEFAULT_PDS3_OPTIONS | (pds3_options or {})
        parsed = NewHorizonsLORRI._read_fileinfo(fileinfo, formats='PF',
                                                 pds3_options=options)

        label = None
        label_path = None
        hdulist = None
        if isinstance(parsed, Pds3Label):
            label = parsed.as_dict()
            label_path = FCPath(parsed.filepath)
            image_path = NewHorizonsLORRI._fits_path(label, label_path)
            if astrometry:
                if geom != 'spice' or pointing != 'spice':
                    raise oops.OopsValueError('geom and pointing from the FITS header '
                                              'need the FITS file; an astrometry read of '
                                              'a PDS3 label reads only the label')
            else:
                hdulist = NewHorizonsLORRI._read_fileinfo(image_path, formats='F')
        else:
            # oops opens a FITS file from its local copy, so the caller's path, which may
            # be a URL, names it; an HDUList passed in names only the file it was opened
            # from
            hdulist = parsed
            image_path = FCPath(fileinfo if hdulist is not fileinfo
                                else hdulist.filename())

        headers = None
        header_dict = None
        arrays = None
        if hdulist is not None:
            (headers, arrays) = NewHorizonsLORRI._read_fits(hdulist,
                                                            astrometry=astrometry)
            header_dict = dict(headers[0])

        # The label describes the array as its IMAGE object's LINES and LINE_SAMPLES
        if label is not None and arrays is not None:
            shape = (label['IMAGE']['LINES'], label['IMAGE']['LINE_SAMPLES'])
            if arrays['data'].shape != shape:
                raise ValueError(f'the image array has shape {arrays["data"].shape}, '
                                 f'but its PDS3 label gives {shape}')

        label_items = {} if label is None else NewHorizonsLORRI._items_from_label(label)
        fits_items = ({} if header_dict is None
                      else NewHorizonsLORRI._items_from_header(header_dict))
        items = NewHorizonsLORRI._merge_items(label_items, fits_items)

        dict_ = {}
        if header_dict is not None:
            dict_.update({key: value for (key, value)
                          in NewHorizonsLORRI._without_missing(header_dict).items()
                          if key not in _COMMENTARY_KEYWORDS
                          and key not in _CONVERTED_KEYWORDS})
        if label is not None:
            dict_.update(NewHorizonsLORRI._without_missing(label))
        dict_.update(NewHorizonsLORRI._dict_items(items))

        sun_state = None
        radec = None
        if geom == 'fits':
            sun_state = NewHorizonsLORRI._sun_state_from_header(header_dict)
        if pointing != 'spice':
            radec = NewHorizonsLORRI._radec_from_header(header_dict, pointing)

        return _LorriFiles(image_path=image_path, label_path=label_path,
                           label_dict=label, headers=headers, header_dict=header_dict,
                           arrays=arrays, items=items, tstart=items['START_TIME'],
                           texp=items['EXPOSURE_DURATION'],
                           binning_mode=_BINNING_MODES[items['LINES']], dict_=dict_,
                           sun_state=sun_state, radec=radec)

    @staticmethod
    def _fits_path(label, label_path):
        """The FITS file that a PDS3 label's ``^IMAGE`` pointer names.

        The labels write the name of the FITS file in upper case, as in
        "LOR_0030713591_0X633_SCI.FIT", while the archive stores the file in the same
        case as the label file, which is lower case. The name is therefore written in
        the case of the label's own name: lower case if that name is all lower case,
        upper case if it is all upper case, and as the label writes it otherwise.

        Parameters:
            label (dict[str, Any]): The parsed label as a dictionary.
            label_path (FCPath): The path to the label.

        Returns:
            FCPath: The path to the FITS file, in the label's directory.
        """

        name = label['^IMAGE']
        if label_path.name == label_path.name.lower():
            name = name.lower()
        elif label_path.name == label_path.name.upper():
            name = name.upper()
        return label_path.parent / name

    @staticmethod
    def _read_fits(hdulist, *, astrometry):
        """The headers and arrays of an open LORRI FITS file, which is then closed.

        Parameters:
            hdulist (HDUList): The FITS file.
            astrometry (bool): True to read the headers alone.

        Returns:
            tuple[list[astropy.io.fits.Header], dict[str, numpy.ndarray] | None]: A copy
            of the header of every HDU, and the arrays keyed by subfield name, None if
            `astrometry` is True: `data`, the primary image, and `error` and `quality`,
            the extensions of a calibrated product that hold them, each only if the file
            has it. Each array is a copy, so it outlives the file.

        Raises:
            AstropyUserWarning: If the file ends before the end that its headers give.
            ValueError: If the primary array is not two-dimensional.
        """

        try:
            # astropy only warns of a file that ends early, giving the arrays it can,
            # and silently leaves out an extension of which only part is present
            with warnings.catch_warnings():
                warnings.filterwarnings('error', message=_TRUNCATED_WARNING,
                                        category=AstropyUserWarning)
                headers = [hdu.header.copy() for hdu in hdulist]
                if astrometry:
                    return (headers, None)

                arrays = {'data': np.array(hdulist[0].data)}
                for hdu in hdulist[1:]:
                    extname = str(hdu.header.get('EXTNAME', '')).strip().upper()
                    if extname == _ERROR_EXTNAME:
                        arrays['error'] = np.array(hdu.data)
                    elif extname == _QUALITY_EXTNAME:
                        arrays['quality'] = np.array(hdu.data)
        finally:
            hdulist.close()

        if arrays['data'].ndim != 2:
            raise ValueError(f'the primary array is not one image band: shape '
                             f'{arrays["data"].shape}')
        return (headers, arrays)

    ######################################################################################
    # Interpreting the metadata
    ######################################################################################

    @staticmethod
    def _is_missing(value):
        """True if a value means the item is missing.

        Parameters:
            value (Any): The value.

        Returns:
            bool: True if `value` is None, as for a FITS card without a value; one of
            :data:`oops.host.PDS3_NULL_TEXTS`, ignoring case and surrounding blanks; the
            FITS fill value -1e32; or a list, tuple or set whose elements are all
            missing.
        """

        if value is None or oops.host.is_pds3_null(value):
            return True
        if isinstance(value, float):
            return value == _FITS_FILL
        if isinstance(value, (list, tuple, set, frozenset)):
            return all(NewHorizonsLORRI._is_missing(item) for item in value)
        return False

    @staticmethod
    def _without_missing(dict_):
        """A dictionary without the values that count as missing.

        Parameters:
            dict_ (dict[str, Any]): A label or a FITS header as a dictionary, which may
                hold dictionaries in turn.

        Returns:
            dict[str, Any]: A copy of `dict_` without any key whose value is missing (see
            :meth:`_is_missing`), at any depth.
        """

        return {key: (NewHorizonsLORRI._without_missing(value)
                      if isinstance(value, dict) else value)
                for (key, value) in dict_.items()
                if not NewHorizonsLORRI._is_missing(value)}

    @staticmethod
    def _items_from_pds3(dict_):
        """The metadata items that a PDS3 label or an index row gives.

        Parameters:
            dict_ (dict[str, Any]): The label as a dictionary, or the index row.

        Returns:
            dict[str, Any]: The items, keyed by PDS3 name, leaving out any that is
            missing or counts as missing: START_TIME, in seconds TDB; EXPOSURE_DURATION,
            in seconds; and TARGET_NAME, SPACECRAFT_CLOCK_START_COUNT and
            SPACECRAFT_CLOCK_STOP_COUNT, as Python strings.
        """

        values = {key: dict_[key] for key in _WINNING_SOURCE
                  if key in dict_ and not NewHorizonsLORRI._is_missing(dict_[key])}

        # A label gives START_TIME as a datetime, an index row as an ISO string
        items = {}
        if 'START_TIME' in values:
            start_time = values['START_TIME']
            if isinstance(start_time, dt.datetime):
                start_time = start_time.isoformat()
            items['START_TIME'] = julian.tdb_from_tai(julian.tai_from_iso(
                str(start_time)))
        if 'EXPOSURE_DURATION' in values:
            items['EXPOSURE_DURATION'] = float(values['EXPOSURE_DURATION'])

        # The table reader returns NumPy strings; the subfields are Python strings
        for key in ('TARGET_NAME', 'SPACECRAFT_CLOCK_START_COUNT',
                    'SPACECRAFT_CLOCK_STOP_COUNT'):
            if key in values:
                items[key] = str(values[key])

        return items

    @staticmethod
    def _items_from_label(label):
        """The metadata items that a PDS3 label gives, with LINES from its IMAGE object.

        Parameters:
            label (dict[str, Any]): The parsed label as a dictionary.

        Returns:
            dict[str, Any]: The items as for :meth:`_items_from_pds3`, and LINES where
            the label's IMAGE object gives it.
        """

        items = NewHorizonsLORRI._items_from_pds3(label)
        lines = label.get('IMAGE', {}).get('LINES')
        if not NewHorizonsLORRI._is_missing(lines):
            items['LINES'] = lines
        return items

    @staticmethod
    def _items_from_header(header_dict):
        """The metadata items that a FITS header gives, converted to the PDS3 forms.

        Parameters:
            header_dict (dict[str, Any]): The primary header as a dictionary.

        Returns:
            dict[str, Any]: The items as for :meth:`_items_from_pds3`, and LINES, the
            lines of the image. An item is left out if a keyword it needs is missing or
            counts as missing.
        """

        def given(keyword):
            """True if the header gives the keyword a value that is not missing."""
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
            value is that of the item's winning source in :data:`_WINNING_SOURCE`.
        """

        items = fits_items | label_items
        for name in label_items.keys() & fits_items.keys():
            winner = label_items if _WINNING_SOURCE[name] == 'label' else fits_items
            items[name] = winner[name]
        return items

    @staticmethod
    def _dict_items(items):
        """The merged items as `dict` gives them.

        Parameters:
            items (dict[str, Any]): The merged items, which include START_TIME.

        Returns:
            dict[str, Any]: The items under their PDS3 names, without LINES, which the
            label gives inside its IMAGE object. START_TIME is a naive
            :class:`datetime.datetime` in UTC, to the microsecond, the type that the
            label gives it.
        """

        dict_ = {key: value for (key, value) in items.items() if key != 'LINES'}
        iso = julian.iso_from_tai(julian.tai_from_tdb(items['START_TIME']), digits=6)
        dict_['START_TIME'] = dt.datetime.fromisoformat(iso)
        return dict_

    @staticmethod
    def _sclk_from_met(met):
        """A mission elapsed time converted to a clock count in the PDS3 label's form.

        Parameters:
            met (float): Seconds of the spacecraft clock, as the FITS header's STOPMET
                gives them.

        Returns:
            str: The count as "SSSSSSSSSS:TTTTT", ten digits of whole seconds and five of
            the 1/50000-second ticks past them, rounded to the nearest tick.
        """

        (seconds, ticks) = divmod(round(met * _SCLK_TICKS), _SCLK_TICKS)
        return f'{seconds:010d}:{ticks:05d}'

    @staticmethod
    def _target_body(name):
        """The name of the target body for a TARGET_NAME, by a fixed translation.

        Parameters:
            name (str | None): The merged TARGET_NAME; None if the sources lack one.

        Returns:
            str: The name of the target body: "NONE" for an image without a target;
            "PLUTO" for a name beginning with "PLUTO"; the name of a moon of Jupiter
            without its "J<n> " prefix; the name of a named minor planet without its
            number and designation; the name in :data:`_MINOR_PLANET_NAMES` of a
            numbered minor planet without one, and in :data:`_HEADER_TARGET_NAMES` of a
            FITS header's name; and otherwise `name` itself.
        """

        if name is None or name.strip() in _NO_TARGET_NAMES:
            return 'NONE'
        if name.startswith('PLUTO'):
            return 'PLUTO'
        if match := _JUPITER_MOON_NAME.fullmatch(name):
            return match[1]
        if match := _NAMED_MINOR_PLANET.fullmatch(name):
            return match[1]
        if match := _NUMBERED_MINOR_PLANET.fullmatch(name):
            return _MINOR_PLANET_NAMES.get(match[1], name)
        return _HEADER_TARGET_NAMES.get(name, name)

    @staticmethod
    def _sun_state_from_header(header_dict):
        """The spacecraft's position and velocity relative to the Sun from the FITS
        header.

        Parameters:
            header_dict (dict[str, Any]): The primary header as a dictionary.

        Returns:
            tuple[oops.Vector3, oops.Vector3, float]: The position in km and the velocity
            in km/s in J2000, and their epoch, SPCSCET, in seconds TDB.
        """

        # The header gives the vectors from the spacecraft to the Sun
        pos = oops.Vector3([-header_dict['SPCSSCX'], -header_dict['SPCSSCY'],
                            -header_dict['SPCSSCZ']])
        vel = oops.Vector3([-header_dict['SPCSSCVX'], -header_dict['SPCSSCVY'],
                            -header_dict['SPCSSCVZ']])
        return (pos, vel, header_dict['SPCSCET'])

    @staticmethod
    def _radec_from_header(header_dict, pointing):
        """The pointing of the camera as the FITS header gives it.

        Parameters:
            header_dict (dict[str, Any]): The primary header as a dictionary.
            pointing (str): "fits", "fitsapp" or "fits90", as for :meth:`from_file`.

        Returns:
            tuple[float, float, float, float]: The right ascension and declination of
            the center of the image and the north clock angle, in degrees, and their
            epoch, SPCSCET, in seconds TDB. For "fits90", 90 degrees are added to the
            clock angle of an image taken in or before 2012.
        """

        # The boresight; FITS pixel (1,1) is the center of the lower left pixel
        u_center = header_dict['NAXIS1'] / 2. + 0.5
        v_center = header_dict['NAXIS2'] / 2. + 0.5
        (ra_deg, dec_deg) = NewHorizonsLORRI._radec_from_uv(u_center, v_center,
                                                            header_dict)
        north_clock_deg = header_dict['SPCEMEN']
        if pointing == 'fits90' and int(header_dict['SPCUTCID'][:4]) <= _FITS90_LAST_YEAR:
            north_clock_deg += 90.
        return (ra_deg, dec_deg, north_clock_deg, header_dict['SPCSCET'])

    @staticmethod
    def _radec_from_uv(u, v, header_dict):
        """The right ascension and declination at FITS pixel coordinates `(u,v)`.

        The gnomonic projection of the header's WCS parameters, adapted from STScI
        source. Here `(u,v) = (1,1)` is the center of the lower left pixel, unlike oops,
        in which `(0,0)` is the corner of the image and pixel centers have half-integer
        values.

        Parameters:
            u (float): The FITS horizontal pixel coordinate.
            v (float): The FITS vertical pixel coordinate.
            header_dict (dict[str, Any]): The FITS header, which provides the WCS
                parameters.

        Returns:
            tuple[float, float]: The right ascension and declination in degrees.
        """

        dx = u - header_dict['CRPIX1']
        dy = v - header_dict['CRPIX2']

        xi  = (header_dict['CD1_1'] * dx + header_dict['CD1_2'] * dy) / oops.DPR
        eta = (header_dict['CD2_1'] * dx + header_dict['CD2_2'] * dy) / oops.DPR

        ra0  = header_dict['CRVAL1'] / oops.DPR
        dec0 = header_dict['CRVAL2'] / oops.DPR

        denom = np.cos(dec0) - eta * np.sin(dec0)

        # The two-argument arctangent keeps the quadrant of a point beyond a celestial
        # pole, where denom is negative and arctan(xi / denom) is 180 degrees off
        ra  = oops.DPR * (np.arctan2(xi, denom) + ra0)
        dec = oops.DPR * np.arctan((eta * np.cos(dec0) + np.sin(dec0))
                                   / np.sqrt(denom**2 + xi**2))

        return (ra, dec)

    ######################################################################################
    # Building the observation
    ######################################################################################

    @staticmethod
    def _check_option(name, value, allowed):
        """Raise OopsValueError unless `value` is one of `allowed`."""

        if value not in allowed:
            raise oops.OopsValueError(f'invalid {name} {value!r}; use one of '
                                      + ', '.join(repr(item) for item in allowed))

    @staticmethod
    def _path_from_sun_state(sun_state, id_name):
        """The spacecraft's path as the FITS header gives it.

        Parameters:
            sun_state (tuple[oops.Vector3, oops.Vector3, float]): The position, velocity
                and epoch, as :meth:`_sun_state_from_header` gives them.
            id_name (str): The name that makes the path's ID unique.

        Returns:
            oops.Path: A linear path relative to the Sun, through the spacecraft's
            position and velocity at the epoch.
        """

        (pos, vel, epoch) = sun_state
        sc_path = oops.path.LinearPath((pos, vel), epoch, 'SUN', frame='J2000',
                                       path_id='.NH_PATH_' + id_name)
        return oops.Path.as_waypoint(sc_path)

    @staticmethod
    def _frame_from_radec(radec, id_name, *, pointing, path):
        """The camera frame as the FITS header gives it.

        Parameters:
            radec (tuple[float, float, float, float]): The right ascension, declination
                and north clock angle in degrees, and their epoch, as
                :meth:`_radec_from_header` gives them.
            id_name (str): The name that makes the frame's ID unique.
            pointing (str): "fits", "fitsapp" or "fits90", as for :meth:`from_file`.
            path (str | oops.Path): The spacecraft's path, needed for "fitsapp".

        Returns:
            oops.Frame: A fixed frame relative to J2000.
        """

        (ra_deg, dec_deg, north_clock_deg, epoch) = radec
        if pointing == 'fitsapp':
            # Insert the apparent vector as the actual one to reverse the aberration
            event = oops.Path.as_path(path).event_at_time(epoch)
            event.neg_arr_j2000 = oops.Vector3.from_ra_dec_length(
                ra_deg * oops.RPD, dec_deg * oops.RPD, recursive=False)
            (ra, dec) = event.ra_and_dec(apparent=True)
            ra_deg  = ra.vals  * oops.DPR
            dec_deg = dec.vals * oops.DPR

        frame = oops.frame.Cmatrix.from_ra_dec(ra_deg, dec_deg, north_clock_deg,
                                               oops.Frame.J2000,
                                               frame_id='.NH_FRAME_' + id_name)
        return oops.Frame.as_wayframe(frame)

    @staticmethod
    def _iof_calibrations(header_dict, texp, *, target, tdb):
        """The solar-normalized I/F calibration of a calibrated product, if it can be
        built.

        Parameters:
            header_dict (dict[str, Any]): The primary FITS header as a dictionary.
            texp (float): The exposure duration in seconds.
            target (str | oops.Body): The target of the observation.
            tdb (float): The middle of the exposure in seconds TDB.

        Returns:
            list[oops.Calibration]: A FlatCalib named "I/F" whose factor is
            pi r^2 / (texp RSOLAR F), in which r is the distance of the target from the
            Sun at `tdb` in AU from the SP kernels and F the solar flux at 1 AU; empty
            if the header gives no RSOLAR, as for a raw product, if the exposure is
            zero, or if the target has no ephemeris: no registered path, or no position
            in the SP kernels at `tdb`.
        """

        rsolar = header_dict.get('RSOLAR')
        if NewHorizonsLORRI._is_missing(rsolar) or texp == 0.:
            return []

        try:
            path = oops.Path.as_waypoint(target.path if isinstance(target, oops.Body)
                                         else target)
        except KeyError:            # a target without a registered path
            return []

        try:
            event = path.wrt('SUN').event_at_time(tdb)
        except OSError:             # cspyce's error for SPKINSUFFDATA
            return []

        solar_range = event.pos.norm().vals / _AU
        factor = np.pi * solar_range**2 / (texp * rsolar * _SOLAR_FLUX)
        return [oops.calibration.FlatCalib('I/F', factor)]

    @staticmethod
    def _snapshot_from_row(row, *, fov_type, overrides):
        """A Snapshot built from a row of an index merged with its supplement.

        Parameters:
            row (dict[str, Any]): The row, as the table reader gives it.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_index`.
        """

        items = NewHorizonsLORRI._items_from_pds3(row)
        params = {
            'cadence'     : oops.cadence.SnapCadence(items['START_TIME'],
                                                     items['EXPOSURE_DURATION']),
            'fov'         : NewHorizonsLORRI._FOVS[str(row['BINNING_MODE']).upper(),
                                                   fov_type],
            'path'        : _NewHorizons.PATH_ID,
            'frame'       : _SPICE_FRAME_NAME,
            'target'      : NewHorizonsLORRI._target_body(items.get('TARGET_NAME')),
            'lightsource' : 'SUN',
            'calibrations': [],
        }
        oops.Host._apply_overrides(params, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        obs.insert_subfield('dict', row)
        NewHorizonsLORRI._insert_common_subfields(obs, items)
        oops.Host._set_path_attributes(obs, f'{row["VOLUME_ID"]}/'
                                            f'{row["FILE_SPECIFICATION_NAME"]}', None)
        return obs

    @staticmethod
    def _insert_common_subfields(obs, items):
        """Insert the subfields that file and index observations share.

        Parameters:
            obs (Snapshot): The observation.
            items (dict[str, Any]): The metadata items; the clock counts become subfields
                where they are given.
        """

        for (subfield, name) in [('spacecraft_clock_start_count',
                                  'SPACECRAFT_CLOCK_START_COUNT'),
                                 ('spacecraft_clock_stop_count',
                                  'SPACECRAFT_CLOCK_STOP_COUNT')]:
            if name in items:
                obs.insert_subfield(subfield, items[name])

        obs.insert_subfield('spice_to_frame', _SPICE_TO_FRAME)
        obs.insert_subfield('spice_frame_name', _SPICE_FRAME_NAME)
        obs.insert_subfield('spice_frame_id', _SPICE_FRAME_ID)
        obs.insert_subfield('ck_object_id', _CK_OBJECT_ID)
        obs.insert_subfield('ck_sclk_id', _CK_SCLK_ID)
        obs.insert_subfield('spice_frame_is_frozen', _SPICE_FRAME_IS_FROZEN)

    @staticmethod
    def _selected(obs, select):
        """The observation that `select` chooses from a file or row that holds one.

        Parameters:
            obs (Snapshot): The observation.
            select (int | slice | tuple[int | slice, ...] | None): The selection.

        Returns:
            Snapshot | tuple[Snapshot, ...]: `obs` itself if `select` is None; otherwise
            what :meth:`oops.host.Host._apply_select` chooses from the one-item tuple
            holding it.
        """

        if select is None:
            return obs
        return oops.Host._apply_select((obs,), select)

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
        so the rows are never needed to recognize it.
        """

        return NewHorizonsLORRI._detect_in_pds3(label)

    ######################################################################################
    # Initialization
    ######################################################################################

    @staticmethod
    def initialize(*, asof=None, meta=None):
        """Initialize key information about the LORRI instrument.

        Furnishes the mission's kernels, builds the fields of view, and registers the
        camera frame NH_LORRI. :meth:`from_file` and :meth:`from_index` call it with
        the defaults if it has not been called. After the first call, later calls to
        this function are ignored.

        Parameters:
            asof (str, optional): Only use SPICE kernels that existed before this date;
                None to ignore.
            meta (str, optional): The name of a metakernel to furnish instead of the
                kernels selected from the database.
        """

        if NewHorizonsLORRI._initialized:
            return

        _NewHorizons.initialize(asof=asof, meta=meta)

        # Display directions: [u,v] = [right,down]
        shape = (_LORRI_PIXELS, _LORRI_PIXELS)
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

        # The SPICE frame NH_LORRI is registered as NH_LORRI_FLIPPED, and the camera
        # frame NH_LORRI flips its y and z axes
        flipped = oops.frame.SpiceFrame(_SPICE_FRAME_NAME, frame_id=_FLIPPED_FRAME_ID)
        _ = oops.frame.Cmatrix(_SPICE_TO_FRAME, flipped, frame_id=_SPICE_FRAME_NAME)

        NewHorizonsLORRI._initialized = True

    @staticmethod
    def _reset():
        """Reset the internal New Horizons LORRI parameters, so that the next
        :meth:`initialize` runs again.
        """

        NewHorizonsLORRI._FOVS = {}
        NewHorizonsLORRI._initialized = False
        _NewHorizons._reset()


NewHorizonsLORRI._register()

##########################################################################################
