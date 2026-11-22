##########################################################################################
# host_new_horizons/oops/lorri.py
##########################################################################################
"""The New Horizons LORRI host, which makes a Snapshot from each LORRI image.

The :doc:`user guide </user_guide/user_guide>` describes the products it reads, how it
merges the two labels of each product, and the errors it raises.
"""

import collections
import datetime as dt
import pathlib
import re

import numpy as np

import cspyce
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

# What OopsSpiceError reports when the C kernels lack the pointing of the camera
_NO_POINTING = 'no New Horizons LORRI pointing'

# The name of the full-resolution detector in the instrument kernel
_IK_NAME = 'NH_LORRI_1X1'

# The dimensions of an unbinned LORRI image in pixels, as the instrument kernel must give
# them
_LORRI_PIXELS = 1024

# The binning mode of an image, keyed by its number of lines
_BINNING_MODES = {1024: '1X1', 256: '4X4'}

# The pixel scale of the flat FOV in radians, the IFOV of the instrument kernel rounded
_FLAT_SCALE = 4.96e-6

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
    '486958': 'ARROKOTH',       # "ASTEROID 486958 (2014 MU69)", NAIF ID 2486958
}

# The bodies whose name in the FITS header of an image is one SPICE does not know, and
# the name it knows. The header names Arrokoth "MU69" on all 1,078 of its images in a
# scan of 7,195 products, so an image read without its label still has a target body.
_HEADER_TARGET_NAMES = {
    'MU69': 'ARROKOTH',
}

# The texts that mean a metadata item is missing are oops.host.PDS3_NULL_TEXTS, compared
# ignoring case and surrounding blanks (oops.host.is_pds3_null()): the PDS3 null
# constants; "---" and "UNKNOWN", which the FITS header writes for a body or a request it
# does not name (SPCCBTNM, SPCTCB, SPC1CB, REQCOMM, ...) and the label for an observation
# it does not describe (NEWHORIZONS:OBSERVATION_DESC); the empty string, which the header
# writes for a file it does not name (REFSUBIM); and "NOT APPLICABLE", which is absent
# from every LORRI source scanned (2026-11-19: as a value in all 43,441 product labels,
# anywhere in the 16 metadata index and supplement tables and their labels, and anywhere
# in the FITS headers of 1,400 random products), so counting it as missing leaves the
# metadata of every current product as it was. A list whose elements are all missing is
# missing too. These values still appear in `label_dict`, `header_dict` and `headers`,
# but nowhere else.

# The number that the FITS header writes for a value its pipeline cannot compute, such as
# the target geometry of an image without a target (SPCTSORN, SPCTRANG, ...) or the
# body-fixed coordinates of a body without a rotation model (SPCTSCLA, ...); the label
# gives "N/A" for the same values. It too counts as missing, and stays in `header_dict`
# and `headers`.
_FITS_FILL = -1.e32

# The ticks of the New Horizons spacecraft clock per second
_SCLK_TICKS = 50000

# The columns that identify the same image in an index and its supplement, by which the
# two tables are paired (see MERGE_SUPPLEMENT_BY_FILE). The key needs VOLUME_ID, which
# oops leaves out of its default key: the cumulative tables of NHxxLO_1999 and
# NHxxLO_2999 list four products twice, once for NHKELO and once for NHK2LO, under the
# same FILE_SPECIFICATION_NAME. The host also needs VOLUME_ID for the path of each row's
# file. Giving the key makes oops refuse a table that lacks either column when
# from_index() is called; a row of `row_dicts` that lacks either raises HostFormatError
# when the iteration reaches it.
_INDEX_KEYS = ('VOLUME_ID', 'FILE_SPECIFICATION_NAME')

# The names that the index and the supplement give five moons of Jupiter, and the names
# that their product labels give them. Each of these moons has one name in all 2,237
# labels of NHJULO_1001 and NHJULO_2001: "J1 IO" in all 358 labels of Io, and likewise for
# the other four. GANYMEDE and CALLISTO have no "J<n> " in front in any label, and the
# tables agree with the labels for them. An index row's TARGET_NAME is translated to the
# label's form, so that from_file() and from_index() give the same TARGET_NAME and
# `label_target`. A corrected table gives the label's form, which is left as it is.
# SETI/rms-data-projects#262 item 6
_SUPPLEMENT_TARGET_NAMES = {
    'IO'        : 'J1 IO',
    'EUROPA'    : 'J2 EUROPA',
    'HIMALIA'   : 'J6 HIMALIA',
    'ELARA'     : 'J7 ELARA',
    'CALLIRRHOE': 'J17 CALLIRRHOE',
}

# The null constants that the label of a supplement declares are not all null values.
# The table reader masks a value equal to its column's constant. A masked -999 (angles)
# or -999999 (the approximate pixel of the target) is never a real value, and neither is
# a masked zero distance; but zero is also the constant of the components of vectors, and
# there a zero can be real. A scan of the 43,441 rows of the 14 volume supplements
# (2026-11-19) found:
# - The vectors from the spacecraft to the target and from the target to the Sun, the
#   distance of each, and the target's velocity are all zero in the 5,207 rows of images
#   without target geometry, whose labels give "N/A" for them. Where the spacecraft-target
#   vector is not zero, a zero is real: the Sun's distance from itself is zero in the 288
#   rows whose target is the Sun, and the velocity of a star (672 rows), or one component
#   of it (90 rows), rounds to zero in the table's four decimal places. So these columns
#   are missing exactly when the spacecraft-target vector is zero.
# - Any other vector is missing when all of its components are zero: the quaternion in 48
#   rows, and the spacecraft's velocity relative to the Earth, which is missing from
#   every label and which every row gives as zeros, although only its X and Z declare
#   zero as null.
# These rules apply only to a value that the reader masks, which is one equal to the null
# constant that the supplement's label declares. Once the supplement is corrected to
# declare a fill value that is never real, a zero is no longer masked and is kept.
# SETI/rms-data-projects#262 item 4
_TARGET_GEOMETRY = re.compile(r'(SC_TARGET|TARGET_SUN)_(POSITION|VELOCITY)_VECTOR_[XYZ]'
                              r'|TARGET_CENTER_DISTANCE|SOLAR_DISTANCE')
_TARGET_POSITION = tuple(f'SC_TARGET_POSITION_VECTOR_{axis}' for axis in 'XYZ')

# A column that is one component of a vector, and the name of the vector. The supplements
# spell the quaternion "QUATERNIAN"; the label's spelling, which a corrected supplement
# would use, is accepted too. SETI/rms-data-projects#262 item 5
_VECTOR_COMPONENT = re.compile(r'(.+_VECTOR)_[XYZ]|(QUATERNI[AO]N)_[1-4]')

# The solar flux at 1 AU at the pivot wavelength of LORRI, 6076.2 Angstroms, in
# erg/cm^2/s/Angstrom, and the astronomical unit in km. The FITS header of a calibrated
# product gives RSOLAR, the conversion from DN/s/pixel to the radiance of a solar
# spectrum in erg/cm^2/s/Angstrom/sr, from which I/F = pi r^2 DN / (texp RSOLAR F) for a
# target r AU from the Sun.
_SOLAR_FLUX = 176.
_AU = 149597870.7

# FITS header keywords that only carry text, left out of `header_dict` and `dict`
_COMMENTARY_KEYWORDS = ('', 'COMMENT', 'HISTORY')

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

# The number of decimal places of seconds in the times of a PDS3 label
_TIME_DIGITS = 3

##########################################################################################
# The resolution rules. Where the detached label and the FITS header both give an item
# and the values differ, the item's winning source gives it, the same source for every
# image: "label" for the detached label, "header" for the FITS header. The winner is the
# source that is accurate for the item. A value missing from the winning source is not a
# difference; the other source's value is used. Evidence: two scans comparing the label
# and the header of 7,897 and 7,195 products, random samples of all 14 volumes
# (2026-11-18), recorded with every other difference between the two in #13.
##########################################################################################

_WINNING_SOURCE = {
    # START_TIME: the header's, SPCSCET - EXPTIME/2. On 2,342 products, it equals the
    # start that the label's own SPACECRAFT_CLOCK_START_COUNT gives through the spacecraft
    # clock kernel to 0.0001 ms; the label's START_TIME is that time rounded to the
    # millisecond, from 0.5 ms early to 0.5 ms late on every product. On frames taken
    # while the spacecraft spun, in cruise as well as at launch, the rounding moves the
    # pointing at mid-exposure by as much as 48.7 pixels (lor_0055952951_0x630_eng,
    # #13).
    'START_TIME': 'header',

    # EXPOSURE_DURATION: the header's. Both are accurate, matching the difference of the
    # clock counts to 0.0008 ms, and they differ only in the last binary digit (238 of
    # 7,195 products, all in NHKELO). The header's keeps START_TIME plus half of the
    # exposure equal to SPCSCET (#13).
    'EXPOSURE_DURATION': 'header',

    # TARGET_NAME: the label's. The header's SPCCBTNM is what the label's NOTE calls the
    # data pipeline's best guess at the target, outside the peer review: "FJUPITER" for
    # JUPITER (64 of 7,195 products), "ANSA" for J RINGS (1,035), pointings such as
    # "RADEC=270.66,-14.63=P4" for KERBEROS (46), and "---" for CALIBRATION. Where the
    # label's name is not the one SPICE knows ("J1 IO", "ASTEROID 486958 (2014 MU69)"),
    # it is still accurate; the target body comes from it by _target_body() (#13).
    'TARGET_NAME': 'label',

    # The clock counts: the label's, the default. They agree on every product scanned.
    'SPACECRAFT_CLOCK_START_COUNT': 'label',
    'SPACECRAFT_CLOCK_STOP_COUNT' : 'label',

    # LINES: the header's. They agree on every product scanned, and the header describes
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


class NewHorizonsLORRI(oops.Host):
    """The New Horizons LORRI host: Snapshot constructors for Long Range Reconnaissance
    Imager images."""

    NAME = 'New Horizons LORRI'

    # The rows of an index and its supplement are paired by file. A survey of all 16 LORRI
    # metadata indexes (2026-11-19, the VOLUME_ID and FILE_SPECIFICATION_NAME of every
    # row) found both tables listing the same files, each once, in each index. In 14 the
    # order is the same. In NHJULO_1001 the index lists lor_0035089864_0x632_eng out of
    # place, last, where the supplement, which is in file order, lists it at row 679
    # (counting from 1), so the two tables differ from row 679 on; in the cumulative
    # NHxxLO_1999 they differ from row 1,925 on for the same reason. Pairing by file gives
    # the same rows on tables in the same order, so it stays right once the index is
    # corrected. SETI/rms-data-projects#262 item 1
    MERGE_SUPPLEMENT_BY_FILE = True

    _INSTRUMENT_KERNEL = None
    _FOVS = {}
    _initialized = False

    @staticmethod
    def from_file(fileinfo, *, geom='spice', pointing='spice', fov_type='fast',
                  return_all_planets=False, pds3_options=None, astrometry=False,
                  select=None, target=None, lightsource=None, path=None, frame=None,
                  cadence=None, fov=None, calibrations=None, timeshift=None,
                  navigation=None, tracker=None, parallel=None):
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
          ticks, rounded to the nearest tick).
        - The lines of the image: NAXIS2, the LINES of the label's IMAGE object.

        A value counts as missing if it is "UNK", "N/A", "NULL", "NOT APPLICABLE",
        "---", "UNKNOWN" or blank, ignoring case, the FITS header's fill value -1e32, or
        a list of such values. Where only one source gives an item, its value is used.
        Where both give it and the values differ, each item takes the value of one
        source, the same for every image:

        - START_TIME, EXPOSURE_DURATION and the lines of the image: the FITS header's.
          The label gives the start time rounded to the millisecond.
        - TARGET_NAME and the clock counts: the label's. The FITS header's target name
          is the data pipeline's guess, such as "FJUPITER" for JUPITER, or a pointing
          such as "RADEC=270.66,-14.63=P4" for KERBEROS.

        If neither source gives the start time, it is SPACECRAFT_CLOCK_START_COUNT
        converted by the spacecraft clock kernel.

        The target body is the name that SPICE knows the body by, translated from
        TARGET_NAME: "J17 CALLIRRHOE" gives CALLIRRHOE, "136472 MAKEMAKE" gives
        MAKEMAKE, "ASTEROID 486958 (2014 MU69)" and the FITS header's "MU69" give
        ARROKOTH, a name beginning with "PLUTO" gives PLUTO, and a missing name or
        "CALIBRATION", which mean the image has no target, give "NONE". Any other name
        is the target body as it is.

        Parameters:
            fileinfo (str | pathlib.Path | FCPath | Pds3Label | HDUList): The path to the
                FITS file or its PDS3 label, or its parsed content as a Pds3Label or an
                HDUList, which is read as the label or the FITS file would be. The FITS
                file is closed once it has been read, including an HDUList that was
                passed in. Parsed content that was not read from a file has no path
                subfields; a Pds3Label of that kind can be read only with `astrometry`,
                since the FITS file beside it cannot be found.
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
            return_all_planets (bool, optional): Include in `spice_kernels` the kernels
                of every planet and of Pluto, each with its system barycenter, besides
                those of the target body.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. If the
                parser's method, by default "fast", cannot read the label, the methods
                after it in the order "fast", "strict", "loose" are tried in turn. The
                option `sources` is always False.
            astrometry (bool, optional):
                True to specify that the returned Observations will only be used for
                timing and astrometry. In this case, no data arrays are read and returned
                in the Observation, and a PDS3 label is read without its FITS file.
            select (None): Must be None. Each file holds a single observation, so there
                is nothing to select from; any other value raises an error.
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

        Returns:
            Snapshot: The observation, with subfields `dict`, `host`, `instrument`,
            `binning_mode`, `time_from_sclk`, `spice_kernels`, `spice_to_frame`,
            `spice_frame_name` and `spice_frame_id` inserted, and these where they apply:

            - `filepath`, `basename` and `image_url`: the FITS file, as a path or URL,
              even when a label was given; absent for parsed content that was not read
              from a file.
            - `abspath`: the local path of the FITS file whose pixels this call read:
              the file itself if it is local, or the copy that the read retrieved into
              the local cache if it is remote. An astrometry read leaves the pixels
              unread, so `abspath` is then absent, as it is for an HDUList built in
              memory.
            - `time_from_sclk`: True if the sources lack the start time, so that it comes
              from SPACECRAFT_CLOCK_START_COUNT through the spacecraft clock kernel;
              False otherwise.
            - `label_target`: the target name as the labels give it, TARGET_NAME as
              merged: the detached label's where it gives one, otherwise the FITS
              header's SPCCBTNM. Absent where neither gives one.
            - `spacecraft_clock_start_count` and `spacecraft_clock_stop_count`: the
              merged clock counts of the start and end of the exposure, spelled as the
              PDS3 label spells them, as in "0030713585:49000". Each is absent where no
              source read gives it.
            - `label_dict`: the PDS3 label as a dictionary, if it was read, with every
              time written as the label writes it and every set as a list in the
              label's order, and without the END key, any value of None, or the keys
              that the parser adds (such as "START_TIME_fmt", "EXPOSURE_DURATION_unit",
              "^IMAGE_offset" and "objects").
            - `header_dict`: the primary FITS header as a dictionary without its blank,
              COMMENT and HISTORY cards and its cards without a value, and `headers`, a
              list of the header of every HDU as a dictionary, with the text of the
              blank, COMMENT and HISTORY cards as a list under each keyword, if the
              FITS file was read. Any other keyword that a header repeats is given at
              each occurrence as "NAME[0]", "NAME[1]", and so on.
            - `calibrations`: for a calibrated product, whose FITS header gives RSOLAR,
              a list of one `oops.calibration.FlatCalib` named "I/F" that converts the
              image's DN to I/F, pi r^2 DN / (texp RSOLAR F). Here RSOLAR converts DN/s
              per pixel to the radiance of a solar spectrum, F = 176 erg/cm^2/s/A is the
              solar flux at 1 AU at the pivot wavelength of LORRI, and r is the distance
              of the target from the Sun at mid-exposure (the middle of `cadence`, if it
              is given) in AU, from the SP kernels.
              The list is empty for a raw product, for an astrometry read of a label,
              for an exposure of zero, and if the target (after the `target` override)
              is "NONE" or the Sun, has no path, or has no position in the kernels at
              that time.
            - `data`: the primary image as stored, which in a raw product includes the
              columns beyond the field of view, unless `astrometry` is True.
            - `error` and `quality`: the error estimates and quality flags of a
              calibrated product, unless `astrometry` is True; absent for a raw product.

            `dict` holds the merged items under their PDS3 names, START_TIME as an ISO
            string to the millisecond, and every other key of each source read under its
            own name, the label's where both sources have a key; the FITS keywords
            converted to the items are left out, and so is any value that counts as
            missing. A key of the label whose value is a list with missing elements
            among others, which is a vector, keeps the list with None for each missing
            element. Every value in `dict`, `label_dict`, `header_dict` and `headers`
            can be written as JSON.

        Raises:
            FileNotFoundError: If the FITS file that a PDS3 label names does not exist,
                except in an astrometry read.
            oops.OopsSpiceError: If the pointing comes from the C kernels and they give
                none at the time when the observation will read them: mid-exposure, or
                the middle of `cadence` if it is given, moved by a `timeshift` of the
                start/stop time ("T"). A `timeshift` of the frame ("F") shifts only the
                camera's fixed rotation above the SPICE frame, so the kernels are read
                at that same time. With a `tracker` override, the kernels are read once,
                at the start, midtime, or end of that exposure, as `tracker` selects.
                Its message gives the time and the FITS file and ends with the SPICE
                error code, such as "SPICE(NOFRAMECONNECT)". Never raised when `frame`
                is given or `pointing` is not "spice". Also raised if no source gives the
                start time and the spacecraft clock kernel cannot convert the clock
                count.
            oops.OopsValueError: If `geom`, `pointing` or `fov_type` is not one of the
                values listed above, if `select` is not None, if `geom` or `pointing` is
                not "spice" in an astrometry read of a PDS3 label, which leaves the FITS
                header unread, or if an option in `pds3_options` is not one the parser
                accepts while the method parses a PDS3 label file; the options are not
                checked for a FITS file or for a Pds3Label passed in.
            oops.host.HostFormatError: If the image is not 1024 or 256 lines, if the
                FITS file is not valid FITS, if a Pds3Label passed in was not read from a
                file and `astrometry` is False, so that no file path is known for its
                FITS file, or if the sources read lack the exposure duration, the lines
                of the image, or both the start time and the clock count.
            oops.host.PdsHostError: If no method of the PDS3 label parser can read the
                label, if the label lacks an ``^IMAGE`` pointer to a FITS file, or, except
                in an astrometry read, if the name it points to does not end in ".fit" or
                ".fits" in either case.
            oops.host.WrongHostError: If `fileinfo` is not a New Horizons LORRI PDS3 label
                or FITS file.
            OSError: If the file cannot be read.
            TypeError: If a keyword argument is not one of those listed.
        """

        oops.Host._check_select(select, unit='New Horizons LORRI file')
        NewHorizonsLORRI._check_option('geom', geom, _GEOMS)
        NewHorizonsLORRI._check_option('pointing', pointing, _POINTINGS)
        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)

        NewHorizonsLORRI.initialize()   # define everything the first time through; use
                                        # defaults unless initialize() was called

        # Default Pds3Label options are method="fast"; the parser's source text is never
        # wanted, in label_dict or anywhere else
        pds3_options = oops.Host._pds3_options_with_defaults(pds3_options)
        pds3_options['sources'] = False

        # A label that the requested parser method cannot read is parsed again by each
        # later method of "fast", "strict" and "loose"
        if isinstance(fileinfo, (str, pathlib.Path, FCPath)):
            filepath = FCPath(fileinfo)
            parsed = NewHorizonsLORRI._read_fileinfo(filepath, formats='PF',
                                                     pds3_options=pds3_options,
                                                     pds3_fallback=True)
        else:
            filepath = None
            parsed = fileinfo
            if isinstance(parsed, Pds3Label) and not NewHorizonsLORRI._detect_in_pds3(
                    parsed):
                source = ('Pds3Label without a file path' if parsed.filepath is None
                          else parsed.filepath)
                raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} PDS3 '
                                               f'label: {source}')
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

        # The path subfields name the FITS file. Parsed content that was not read from a
        # file has none, and a full read of such a label cannot find its FITS file.
        label = None
        hdulist = None
        if isinstance(parsed, Pds3Label):
            label = parsed
            if astrometry:
                # An astrometry read of a label reads the label alone
                if geom != 'spice' or pointing != 'spice':
                    raise oops.OopsValueError('geom and pointing from the FITS header '
                                              'need the FITS file; an astrometry read of '
                                              'a PDS3 label reads only the label')
                if label.filepath is not None:
                    filepath = NewHorizonsLORRI._fits_path(label, check=False)
            elif label.filepath is None:
                raise oops.host.HostFormatError('no file path is known for this '
                                                'Pds3Label, so its FITS file cannot be '
                                                'found')
            else:
                filepath = NewHorizonsLORRI._fits_path(label)
                hdulist = NewHorizonsLORRI._read_fileinfo(filepath, formats='F')
        else:
            hdulist = parsed
            if filepath is None:
                source = NewHorizonsLORRI._hdulist_path(hdulist)
                if source is not None:
                    filepath = FCPath(source)
                if not NewHorizonsLORRI._detect_in_fits(hdulist):
                    hdulist.close()
                    name = ('HDUList without a file path' if filepath is None
                            else filepath)
                    raise oops.host.WrongHostError(f'not a {NewHorizonsLORRI.NAME} '
                                                   f'FITS file: {name}')

        # Everything needed is copied out of the file, which is then closed, including
        # an HDUList that oops.Host.from_file() or the caller opened. A read of the
        # arrays has the FITS file locally, at the path that the HDUList was opened from
        # (for a remote file, the path that FCPath.retrieve() returned), and only that
        # path is given as `abspath`.
        headers = None
        arrays = None
        abspath = None
        if hdulist is not None:
            try:
                headers = [hdu.header.copy() for hdu in hdulist]
                if not astrometry:
                    arrays = NewHorizonsLORRI._read_arrays(hdulist)
                    abspath = NewHorizonsLORRI._hdulist_path(hdulist)
            finally:
                hdulist.close()

        return NewHorizonsLORRI._make_snapshot(filepath, abspath=abspath, label=label,
                                               headers=headers, arrays=arrays, geom=geom,
                                               pointing=pointing, fov_type=fov_type,
                                               return_all_planets=return_all_planets,
                                               overrides=overrides)

    @staticmethod
    def from_index(filepath, *, supplement=oops.Host.SUPPLEMENT, label=None,
                   row_dicts=None, pds3_options=None, filter=None, fov_type='fast',
                   return_all_planets=False, select=None, target=None, lightsource=None,
                   path=None, frame=None, cadence=None, fov=None, calibrations=None,
                   timeshift=None, navigation=None, tracker=None, parallel=None):
        """An iterator over Snapshots, one for each row of a New Horizons LORRI index.

        It reads a Ring-Moon Systems Node metadata index and its supplement. When this
        method is called, it checks the arguments and reads only the labels of the index
        and the supplement. The iteration then reads the two tables 100 rows at a time,
        through ``oops.Host._index_rows()``, and builds each row's Snapshot when it
        reaches that row, so the tables are never held in memory whole. The memory used
        stays the same however long the index is only with an rms-oops whose frame and
        path caches release each row's camera frame along with its Snapshot; with one
        whose caches keep them, each row adds a few kilobytes.

        With `row_dicts`, the given rows take the place of the tables, and the call leaves
        both labels unread. The errors below that come from reading a label or a table
        are then never raised: the labels need not exist, the index need not be New
        Horizons LORRI, the name of `filepath` need not contain "_index", and
        `pds3_options` is not checked. What is still checked is `select`, `fov_type`,
        `filter` against the first row, and each row when the iteration reaches it.

        Each index row is paired with the supplement row that gives the same VOLUME_ID
        and FILE_SPECIFICATION_NAME (:attr:`MERGE_SUPPLEMENT_BY_FILE` is True), because
        the index of NHJULO_1001, and so the cumulative index NHxxLO_1999, lists one file
        out of place, while the supplement lists the files in file order. The supplement
        is read ahead only as far as the partner of the current row, and each supplement
        row read ahead of its partner is kept until the partner is reached. The 100-row
        chunks bound the rows parsed at a time, but the rows kept grow with how far the
        supplement is out of step with the index; in the archive's tables, one row at
        most is kept at a time.

        Parameters:
            filepath (str | pathlib.Path | FCPath): The full path to a New Horizons LORRI
                index file or its PDS label, such as ``NHJULO_2001_index.lbl``. Given the
                table, the PDS3 parser reads the whole table, looking for a label in it,
                before it reads the label beside it, so the label is the faster input.
            supplement (str | pathlib.Path | FCPath | object, optional): The path to the
                PDS label of the supplemental index, whose columns are merged into each
                row; a LORRI index needs its supplement. Its rows must give the same
                files as the index's, in any order. Where both tables give a column, the
                supplement's value is kept, because it is closer to the product label:
                the index gives each clock count with its partition, rounds the stop
                count differently, gives PRODUCT_CREATION_TIME as a date alone, and
                shortens some target names. The default, `oops.Host.SUPPLEMENT`, reads
                the file in the directory of `filepath` whose name matches except that
                "_index" is replaced by "_supplemental_index" in upper or lower case.
            label (Pds3Label, optional): The PDS label of `filepath` if it has already
                been read.
            row_dicts (collections.abc.Iterable[dict[str, Any]], optional): The rows of
                the index, each already merged with its supplement row, as a list or any
                other iterable, which is read as the iteration proceeds. If they are
                given, they are used in place of the tables, `supplement`, `label` and
                `pds3_options` are ignored, and the tables and their labels are left
                unread.
            pds3_options (dict[str, Any], optional): Any keyword arguments and their
                values to be passed to the :class:`~pdsparser.Pds3Label` parser. The
                options are checked when the call is made, unless `row_dicts` is given,
                when they are not used. The index's label is parsed with all of them, and
                the supplement's label with the method alone. Each label is parsed with
                one method, "fast" unless they give another, and no later method is
                tried.
            filter (dict[str, Any], optional): A dictionary of index column names and
                their values or a set or tuple of values. If provided, only rows of the
                index in which each named column equals the value, falls in the set of
                values, or falls in a range defined by a tuple of two values, will be
                included. A row without a named column is left out.
            fov_type (str, optional): "fast", "slow" or "flat", as for :meth:`from_file`.
            return_all_planets (bool, optional): Include in `spice_kernels` the kernels
                of every planet and of Pluto, as for :meth:`from_file`.
            select (None): Must be None. Each row describes a single observation, so
                there is nothing to select from; any other value raises an error.
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

        Returns:
            collections.abc.Iterator[Snapshot]: One observation per filtered row of the
            index, in the index's order, each with subfields `dict`, `host`,
            `instrument`, `binning_mode`, `time_from_sclk`, `spice_kernels`, `filepath`,
            `basename`, `image_url`, `spice_to_frame`, `spice_frame_name` and
            `spice_frame_id` inserted, and `label_target`,
            `spacecraft_clock_start_count` and `spacecraft_clock_stop_count` where the
            row gives TARGET_NAME and the clock counts with values that do not count as
            missing, as for :meth:`from_file`. `abspath` is absent, because the FITS
            file is not read. The tables lack RSOLAR, which the I/F calibration needs,
            so the calibrations are empty unless the `calibrations` override gives
            some. The iteration stops at the first row that raises an error. `filter`
            can leave out a row whose Snapshot would raise, but it cannot avoid an error
            in reading or pairing the rows of the tables, because the rows are filtered
            only after they have been read and paired.

            `dict` is the row, merged with its row of the supplement,
            without the masks of the table reader or NumPy types, so that it can be
            written as JSON, and without any value that counts as missing: a PDS3 null
            such as "N/A", and a value that the table reader masks as its column's null
            constant, except a zero that is real. A zero is the null constant of the
            vectors and distances of the target, which are left out where the vector
            from the spacecraft to the target is zero (an image without target
            geometry) and kept elsewhere, as for an image of the Sun; of the other
            vectors, which are left out where all of their components are zero; and of
            the other distances, which are left out. The start of the exposure is
            START_TIME, or, if the row lacks it, SPACECRAFT_CLOCK_START_COUNT
            converted by the spacecraft clock kernel, and `time_from_sclk` says which.
            The exposure is the supplement's EXPOSURE_DURATION, which is rounded to the
            millisecond. The supplement gives five moons of Jupiter a name without the
            "J<n> " that every one of their product labels puts in front, such as "IO"
            for "J1 IO". These five names are written as the labels write them, in
            TARGET_NAME and in `label_target`, so both are the same as from
            :meth:`from_file`. The target body comes from TARGET_NAME as for
            :meth:`from_file`. The C kernels are consulted at the time at which
            :meth:`from_file` would read them, the row's midtime moved by a `cadence`
            override and by a `timeshift` of the start/stop time ("T"), or the start,
            midtime, or end that a `tracker` override selects, but a row without
            pointing still gives an observation, and its frame raises
            `oops.OopsSpiceError` when it is evaluated.
            `filepath`, `basename`, and `image_url` name the FITS file, as a full path
            or URL: the row's FILE_SPECIFICATION_NAME, with the label's
            suffix replaced by the FITS file's in the same case, in the volume directory
            named by its VOLUME_ID. The volume directories are found from the location
            of `filepath`: under "volumes" in place of "metadata" for an index in an RMS
            Node metadata tree
            (".../metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl"), and
            otherwise in the index's own directory.

        Raises:
            FileNotFoundError: When this method is called, if the supplement's label does
                not exist; when the iteration starts, if a table file does not exist.
            oops.OopsSpiceError: When the iteration reaches a row that lacks
                START_TIME, if the spacecraft clock kernel cannot convert its
                SPACECRAFT_CLOCK_START_COUNT.
            oops.OopsTypeError: When the iteration reaches a row, if a range in `filter`
                cannot be compared with that row's value.
            oops.OopsValueError: When this method is called, if `fov_type` is not
                "fast", "slow" or "flat", if `select` is not None, if `filter` uses a
                column that neither table has, or if `pds3_options` names a method the
                parser does not have or an option it does not take. The call checks the
                options when it is made, unless `row_dicts` is given. With `row_dicts`,
                `filter` is checked against the first row: when this method is called if
                `row_dicts` is a sequence, such as a list, and otherwise when the
                iteration reads the first row.
            oops.host.HostFormatError: When the iteration reaches a row that lacks
                VOLUME_ID or FILE_SPECIFICATION_NAME, which only a row of `row_dicts` or a
                null value can, whose BINNING_MODE is missing or not "1X1" or "4X4" in
                either case, or that lacks EXPOSURE_DURATION, or lacks both START_TIME
                and SPACECRAFT_CLOCK_START_COUNT.
            oops.host.HostIndexError: When this method is called, if `supplement` is None
                or an empty string while `row_dicts` is not given; if `supplement` is
                `oops.Host.SUPPLEMENT` and the name of `filepath` lacks "_index" while
                `row_dicts` is not given; or if either label does not describe an index
                table with the columns VOLUME_ID and FILE_SPECIFICATION_NAME. During the
                iteration: after the last row that a table file holds, if it holds fewer
                rows than its label says; for a row longer than the label's row length,
                when the iteration reaches the first of the 100 rows read with it, so
                those 100 rows are not given; when the supplement gives a file a second
                time while the row with its first copy is still unpaired; when the whole
                supplement has been read without finding the file of an index row; and,
                after the last index row, if a table file holds more rows than its label
                says, or if the supplement gives files that the index does not, or gives
                a file more times than the index.
            oops.host.PdsHostError: When this method is called, if the PDS label of
                `filepath` or `supplement` cannot be parsed.
            oops.host.WrongHostError: When this method is called, if `filepath` is not a
                New Horizons LORRI index.
            OSError: When this method is called, if the label of `filepath`, or of a
                `supplement` that is specified, does not exist or cannot be read.
            TypeError: When this method is called, if a keyword argument is not one of
                those listed.
        """

        oops.Host._check_select(select, unit='row of a New Horizons LORRI index')
        NewHorizonsLORRI._check_option('fov_type', fov_type, _FOV_TYPES)

        # The supplement is required, so None or an empty path, either of which the index
        # reader takes to mean "read the index alone", is refused
        if row_dicts is None and not supplement:
            raise oops.host.HostIndexError(f'{NewHorizonsLORRI.NAME} index needs its '
                                           f'supplemental index: {filepath}')

        # Only the labels are read here; the rows are read as the iteration reaches them.
        # Given rows are used as they are: _index_rows() returns them before it reads
        # either label, so the label checks are not made.
        # The other columns that the two tables share, besides the file key, are written
        # differently: the index gives each clock count with its partition, as in
        # "1/0030598438:48850", rounds the stop count differently from the label, gives
        # PRODUCT_CREATION_TIME as a date without a time, and shortens some target names,
        # as in "50000 QUAOAR (2002 LM)", where the supplement gives the label's value.
        # Until the index is corrected, only the file key is compared (validate=()) and
        # the supplement's value of a shared column is kept. A corrected index gives the
        # same values as the supplement, so the rows are the same either way. The
        # supplement still differs from the labels elsewhere: it rounds
        # EXPOSURE_DURATION to the millisecond (F7.3, 0.501 where the label gives
        # 0.5006), leaves out the "12150:" prefix of some SEQUENCE_IDs, and names five
        # moons of Jupiter in another way (see _SUPPLEMENT_TARGET_NAMES).
        # SETI/rms-data-projects#262 items 2, 3 and 6
        rows = NewHorizonsLORRI._index_rows(filepath, supplement=supplement, label=label,
                                            row_dicts=row_dicts,
                                            pds3_options=pds3_options, filter=filter,
                                            validate=(), file_key=_INDEX_KEYS,
                                            prefer_supplement=True)

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

        # Build each row's Snapshot when the iteration reaches it
        volumes = oops.Host._volumes_directory(filepath)
        return (NewHorizonsLORRI._snapshot_from_row(row_dict, volumes=volumes,
                                                    fov_type=fov_type,
                                                    return_all_planets=return_all_planets,
                                                    overrides=overrides)
                for row_dict in rows)

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
    def _fits_path(label, *, check=True):
        """The FITS file that a PDS3 label's ``^IMAGE`` pointer names.

        The label may spell the name in a different case from the file system. The
        archive names the FITS file in the same case as the label file, so the name is
        tried first in the case of the label's own name, then as written, then in lower
        and in upper case. On a file system that ignores case, the first spelling tried
        is the one returned.

        Parameters:
            label (Pds3Label): The parsed label.
            check (bool, optional): True to look for the file, which must then be named
                with the suffix ".fit" or ".fits" in either case; False to return the
                first spelling without looking, so that no FITS file is needed.

        Returns:
            FCPath: The path to the FITS file, in the label's directory.

        Raises:
            FileNotFoundError: If `check` is True and no spelling of the name exists.
            oops.host.PdsHostError: If the label has no ``^IMAGE`` pointer to a
                separate file, or if `check` is True and the name it points to does not
                end in ".fit" or ".fits".
        """

        name = label.get('^IMAGE')
        if not isinstance(name, str):
            raise oops.host.PdsHostError(f'PDS3 label lacks an ^IMAGE pointer to a FITS '
                                         f'file: {label.filepath}')

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

        # oops opens a FITS file only by its suffix, which every archive product has
        if not name.lower().endswith(('.fit', '.fits')):
            raise oops.host.PdsHostError(f'PDS3 label points to {name}, whose name does '
                                         f'not end in .fit or .fits: {label.filepath}')

        for option in dict.fromkeys([first, name, name.lower(), name.upper()]):
            path = parent / option
            if path.exists():
                return path

        raise FileNotFoundError(f'FITS file {name} of PDS3 label not found: '
                                f'{label.filepath}')

    @staticmethod
    def _hdulist_path(hdulist):
        """The path of the local file that an HDUList was read from.

        astropy names an HDUList read from a stream, such as an ``io.BytesIO``, after
        the stream's type, as in "<class '_io.BytesIO'>", so its name counts as a path
        only if it names a file.

        Parameters:
            hdulist (HDUList): The FITS file.

        Returns:
            str | None: The path of the file, or None if the HDUList was built in memory
            or read from a stream.
        """

        name = hdulist.filename()
        if name is None or not pathlib.Path(name).is_file():
            return None
        return name

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
    def _make_snapshot(filepath, *, abspath, label, headers, arrays, geom, pointing,
                       fov_type, return_all_planets, overrides):
        """A Snapshot built from a PDS3 label, the headers and arrays of a FITS file, or
        both.

        Parameters:
            filepath (FCPath | None): The path to the FITS file; None for parsed content
                that was not read from a file.
            abspath (str | None): The local path of the FITS file whose arrays were read;
                None if they were not read from a local file.
            label (Pds3Label | None): The detached label; None if it was not read.
            headers (list[astropy.io.fits.Header] | None): The header of every HDU; None
                if the FITS file was not read.
            arrays (dict[str, numpy.ndarray] | None): The arrays as :meth:`_read_arrays`
                returns them; None to omit them.
            geom (str): "spice" or "fits", as for :meth:`from_file`.
            pointing (str): "spice", "fits", "fitsapp" or "fits90", as for
                :meth:`from_file`.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            return_all_planets (bool): True to include the kernels of every planet in
                `spice_kernels`, as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_file`.

        Raises:
            oops.OopsSpiceError: If the pointing comes from the C kernels and they give
                none at the pointing time, or if the sources lack the start time and the
                spacecraft clock kernel cannot convert the clock count.
            oops.host.HostFormatError: If the image is not 1024 or 256 lines, or if the
                sources lack the exposure duration, the lines, or both the start time and
                the clock count.
        """

        header = None if headers is None else headers[0]

        label_dict = None
        label_items = {}
        if label is not None:
            raw_label_dict = label.as_dict()
            label_dict = NewHorizonsLORRI._label_as_json(raw_label_dict)
            label_items = NewHorizonsLORRI._items_from_pds3(raw_label_dict)
            lines = raw_label_dict.get('IMAGE', {}).get('LINES')
            if not NewHorizonsLORRI._is_missing(lines):
                label_items['LINES'] = lines

        header_dict = None
        fits_items = {}
        if header is not None:
            header_dict = {key: value for (key, value)
                           in NewHorizonsLORRI._header_cards(header)
                           if key not in _COMMENTARY_KEYWORDS and value is not None}
            fits_items = NewHorizonsLORRI._items_from_header(header_dict)

        items = NewHorizonsLORRI._merge_items(label_items, fits_items)
        NewHorizonsLORRI._require(items, 'EXPOSURE_DURATION', filepath)
        time_from_sclk = NewHorizonsLORRI._add_start_time_from_clock(items, filepath)
        for name in ('START_TIME', 'LINES'):
            NewHorizonsLORRI._require(items, name, filepath)

        lines = items['LINES']
        if lines not in _BINNING_MODES:
            raise oops.host.HostFormatError(f'unrecognized LORRI image shape: '
                                            f'{lines} lines in '
                                            f'{NewHorizonsLORRI._source(filepath)}')
        binning_mode = _BINNING_MODES[lines]

        tstart = items['START_TIME']
        texp = items['EXPOSURE_DURATION']
        cadence = oops.cadence.SnapCadence(tstart, texp)
        target_name = NewHorizonsLORRI._target_body(items.get('TARGET_NAME'))

        # Define the path and frame. A default that an override replaces is not built,
        # so the header need not give its values, and the apparent pointing uses the
        # observer that the Snapshot will have. geom and pointing other than "spice"
        # are refused before this point when the header was not read. The header's
        # position, velocity and pointing are those of its own mid-exposure, SPCSCET.
        # The IDs of a path and frame from the header are made unique by the file's
        # name, or by the clock count for an HDUList that was not read from a file.
        if filepath is not None:
            id_name = filepath.name
        elif header is not None:
            id_name = str(header.get('SPCSCLK0'))
        if overrides['path'] is not None:
            default_path = overrides['path']
        elif geom == 'spice':
            default_path = _NewHorizons.PATH_ID
        else:
            default_path = NewHorizonsLORRI._fits_path_from_header(header, id_name,
                                                                   header['SPCSCET'])

        using_cks = overrides['frame'] is None and pointing == 'spice'
        if overrides['frame'] is not None:
            default_frame = overrides['frame']
        elif pointing == 'spice':
            NewHorizonsLORRI._define_camera_frame()
            NewHorizonsLORRI._check_pointing(cadence, overrides, filepath)
            default_frame = _SPICE_FRAME_NAME
        else:
            # Only "fitsapp" needs the epoch, to evaluate the observer's path
            default_frame = NewHorizonsLORRI._fits_frame_from_header(
                header, id_name, pointing=pointing, path=default_path,
                tdb_midtime=header.get('SPCSCET'))

        # Define the Snapshot parameters
        params = {
            'cadence'     : cadence,
            'fov'         : NewHorizonsLORRI._FOVS[binning_mode, fov_type],
            'path'        : default_path,
            'frame'       : default_frame,
            'target'      : target_name,
            'lightsource' : 'SUN',
            'calibrations': [],         # TBD
        }

        oops.Host._apply_overrides(params, filepath, abspath=abspath, **overrides)

        # The I/F calibration of a calibrated product, for the target and at the midtime
        # that the overrides leave
        if overrides['calibrations'] is None and header_dict is not None:
            params['calibrations'] = NewHorizonsLORRI._iof_calibrations(
                header_dict, texp, target=params['target'],
                tdb=params['cadence'].midtime)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        if arrays is not None:
            obs.insert_subfields(arrays)

        # The merged items under their PDS3 names, then every other key of each source
        # without missing values, the label's where both have a key. A null in one
        # source never hides the other's value.
        dict_ = {}
        if label_dict is not None:
            dict_.update(NewHorizonsLORRI._without_missing(label_dict))
        if header_dict is not None:
            # NAXIS2 stays, because LINES is not a key of `dict`
            converted = {keyword for (name, keywords) in _FITS_ITEM_KEYWORDS.items()
                         if name != 'LINES' for keyword in keywords}
            dict_.update({key: value for (key, value)
                          in NewHorizonsLORRI._without_missing(header_dict).items()
                          if key not in converted and key not in dict_})
        dict_.update(NewHorizonsLORRI._dict_items(items))

        obs.insert_subfield('dict', dict_)
        if label_dict is not None:
            obs.insert_subfield('label_dict', label_dict)
        if header_dict is not None:
            obs.insert_subfield('header_dict', header_dict)
            obs.insert_subfield('headers', [NewHorizonsLORRI._header_as_json(header)
                                            for header in headers])
        NewHorizonsLORRI._insert_common_subfields(obs, items, binning_mode=binning_mode,
                                                  target=params['target'], ck=using_cks,
                                                  time_from_sclk=time_from_sclk,
                                                  return_all_planets=return_all_planets)
        return obs

    @staticmethod
    def _snapshot_from_row(row_dict, *, volumes, fov_type, return_all_planets, overrides):
        """A Snapshot built from a row of an index merged with its supplement.

        Parameters:
            row_dict (dict[str, Any]): The row, as the table reader gives it.
            volumes (FCPath): The directory that holds the volume directories, as
                :meth:`oops.Host._volumes_directory` gives it.
            fov_type (str): "fast", "slow" or "flat", as for :meth:`from_file`.
            return_all_planets (bool): True to include the kernels of every planet in
                `spice_kernels`, as for :meth:`from_file`.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.

        Returns:
            Snapshot: The observation, with the subfields listed for :meth:`from_index`.

        Raises:
            oops.OopsSpiceError: If the row lacks START_TIME and the spacecraft clock
                kernel cannot convert its SPACECRAFT_CLOCK_START_COUNT. Missing pointing
                is raised only when the observation's frame is evaluated.
            oops.host.HostFormatError: If the row lacks VOLUME_ID or
                FILE_SPECIFICATION_NAME, if its BINNING_MODE is missing or not "1X1" or
                "4X4", or if it lacks EXPOSURE_DURATION, or lacks both START_TIME and
                SPACECRAFT_CLOCK_START_COUNT.
        """

        row_dict = NewHorizonsLORRI._row_as_json(row_dict)

        # Five moons of Jupiter take the name that their product labels give them
        target_name = row_dict.get('TARGET_NAME')
        if isinstance(target_name, str) and target_name in _SUPPLEMENT_TARGET_NAMES:
            row_dict['TARGET_NAME'] = _SUPPLEMENT_TARGET_NAMES[target_name]

        # The index names the label; the observation names the FITS file beside it. A
        # table without either column is refused before any row is read, but a row of
        # `row_dicts`, or a null value, can still lack one.
        for name in _INDEX_KEYS:
            if name not in row_dict:
                raise oops.host.HostFormatError(f'no {name} in a LORRI index row')
        label_path = oops.Host._index_file_path(volumes, row_dict)
        stem, _, suffix = label_path.name.rpartition('.')
        fits_suffix = 'FIT' if suffix.isupper() else 'fit'
        filepath = label_path.parent / (stem + '.' + fits_suffix)

        # A null binning mode is missing from the row, and so unrecognized
        binning_mode = str(row_dict.get('BINNING_MODE', '')).upper()
        if binning_mode not in _BINNING_MODES.values():
            given = repr(row_dict['BINNING_MODE']) if binning_mode else 'missing'
            raise oops.host.HostFormatError(f'unrecognized LORRI BINNING_MODE: {given} '
                                            f'in {filepath}')

        items = NewHorizonsLORRI._items_from_pds3(row_dict)
        NewHorizonsLORRI._require(items, 'EXPOSURE_DURATION', filepath)
        time_from_sclk = NewHorizonsLORRI._add_start_time_from_clock(items, filepath)
        NewHorizonsLORRI._require(items, 'START_TIME', filepath)

        cadence = oops.cadence.SnapCadence(items['START_TIME'],
                                           items['EXPOSURE_DURATION'])

        # A row whose pointing the C kernels cannot supply is kept, with a frame that
        # raises OopsSpiceError when it is evaluated, so that one image without
        # pointing does not stop the whole index
        using_cks = overrides['frame'] is None
        default_frame = _SPICE_FRAME_NAME
        if using_cks:
            NewHorizonsLORRI._define_camera_frame()
            default_frame = NewHorizonsLORRI._check_pointing(cadence, overrides, filepath,
                                                             defer=True)

        target_name = NewHorizonsLORRI._target_body(items.get('TARGET_NAME'))
        params = {
            'cadence'     : cadence,
            'fov'         : NewHorizonsLORRI._FOVS[binning_mode, fov_type],
            'path'        : _NewHorizons.PATH_ID,
            'frame'       : default_frame,
            'target'      : target_name,
            'lightsource' : 'SUN',
            'calibrations': [],         # TBD
        }
        oops.Host._apply_overrides(params, filepath, **overrides)

        obs = oops.observation.Snapshot(('v', 'u'), **params)
        obs.insert_subfield('dict', row_dict)
        NewHorizonsLORRI._insert_common_subfields(obs, items, binning_mode=binning_mode,
                                                  target=params['target'], ck=using_cks,
                                                  time_from_sclk=time_from_sclk,
                                                  return_all_planets=return_all_planets)
        return obs

    @staticmethod
    def _check_pointing(cadence, overrides, filepath, *, defer=False):
        """Check that the C kernels give the pointing of the camera where the observation
        will read it.

        :meth:`_define_camera_frame` must have been called first, and the overrides must
        not yet be applied. :meth:`oops.Host._check_spice_frame` builds the frame and the
        time as the overrides will, and evaluates the camera frame relative to J2000,
        which evaluates every frame in its chain, and so the SPICE frame beneath it. The
        camera frame itself is a fixed rotation that never consults the kernels, so a
        time shift of the frame ("F") shifts only that rotation, and the kernels are read
        at the observation's midtime, moved by a `cadence` override and by a time shift
        of the cadence ("T"). With a `tracker` override, the kernels are read once, when
        the tracking frame is built, at the start, midtime, or end of that cadence that
        the override selects.

        Parameters:
            cadence (oops.cadence.Cadence): The cadence of the observation before the
                overrides.
            overrides (dict[str, Any]): The values of the standard overrides, as described
                for :meth:`from_file`.
            filepath (FCPath | None): The path to the FITS file, for the message.
            defer (bool, optional): True to return a frame that raises the error when it
                is evaluated, rather than raising it now.

        Returns:
            str | oops.frame.SpiceErrorFrame: The camera frame's ID, "NH_LORRI", if the
            kernels give the pointing; otherwise, with `defer`, a SpiceErrorFrame that
            stands in for the camera frame before the overrides.

        Raises:
            oops.OopsSpiceError: If `defer` is False and the kernels lack the pointing.
                The message names the time and the file and ends with the SPICE error
                code, such as "SPICE(NOFRAMECONNECT)".
        """

        # What "F" should mean for a camera frame that is a fixed rotation over a SPICE
        # frame is an open question, SETI/rms-oops#251; the host follows oops
        return oops.Host._check_spice_frame(_SPICE_FRAME_NAME, cadence, overrides,
                                            description=_NO_POINTING, filepath=filepath,
                                            defer=defer)

    @staticmethod
    def _require(items, name, filepath):
        """Raise HostFormatError if a metadata item the Snapshot needs is missing.

        Parameters:
            items (dict[str, Any]): The metadata items.
            name (str): The name of the item.
            filepath (FCPath | None): The image file, for the message.

        Raises:
            oops.host.HostFormatError: If `items` does not contain `name`.
        """

        if name not in items:
            description = {'LINES': 'lines of the image (IMAGE LINES or NAXIS2)'}
            raise oops.host.HostFormatError(f'no {description.get(name, name)} for '
                                            f'LORRI image '
                                            f'{NewHorizonsLORRI._source(filepath)}')

    @staticmethod
    def _source(filepath):
        """The image file, or a description of parsed content without one, for a message.

        Parameters:
            filepath (FCPath | None): The image file; None for parsed content that was
                not read from a file.

        Returns:
            str | FCPath: `filepath`, or "without a file path" if it is None.
        """

        return 'without a file path' if filepath is None else filepath

    @staticmethod
    def _is_missing(value):
        """True if a value means the item is missing.

        Parameters:
            value (Any): The value.

        Returns:
            bool: True if `value` is None, as for an undefined FITS card; one of
            :data:`oops.host.PDS3_NULL_TEXTS`, ignoring case and surrounding blanks; the
            FITS fill value :data:`_FITS_FILL`; or a list or tuple whose elements are all
            missing.
        """

        if oops.host.is_pds3_null(value):
            return True
        if isinstance(value, float):
            return value == _FITS_FILL
        if isinstance(value, (list, tuple)):
            return all(NewHorizonsLORRI._is_missing(item) for item in value)
        return False

    @staticmethod
    def _without_missing(dict_):
        """A dictionary without the values that count as missing.

        Parameters:
            dict_ (dict[str, Any]): A label, a FITS header or an index row as a
                dictionary, which may hold dictionaries in turn.

        Returns:
            dict[str, Any]: A copy of `dict_` without any key whose value is missing (see
            :meth:`_is_missing`), at any depth. The lists of a LORRI label that have
            missing elements among others are vectors (SC_TARGET_POSITION_VECTOR,
            QUATERNION, ...), in which each position has a meaning, so each missing
            element becomes None in place.
        """

        result = {}
        for (key, value) in dict_.items():
            if NewHorizonsLORRI._is_missing(value):
                continue
            if isinstance(value, dict):
                value = NewHorizonsLORRI._without_missing(value)
            elif isinstance(value, list):
                value = [None if NewHorizonsLORRI._is_missing(item) else item
                         for item in value]
            result[key] = value
        return result

    @staticmethod
    def _label_as_json(label_dict):
        """A PDS3 label as a dictionary that can be written as JSON.

        Parameters:
            label_dict (dict[str, Any]): The label as :meth:`Pds3Label.as_dict` gives it.

        Returns:
            dict[str, Any]: A copy of `label_dict` without the END key, any key whose
            value is None, or the keys that the parser adds (those with a lower-case
            suffix, such as "_fmt", "_unit" and "_offset", and "objects"), in which
            every date or time is the text the label writes for it (its "_fmt" key) and
            every set is a list in the label's order (its "_list" key), at any depth.
        """

        result = {}
        for (key, value) in label_dict.items():
            if key == 'END' or value is None or key != key.upper():
                continue
            result[key] = NewHorizonsLORRI._json_value(value, key, label_dict)
        return result

    @staticmethod
    def _json_value(value, key, siblings):
        """One value of a PDS3 label in a form that can be written as JSON.

        Parameters:
            value (Any): The value.
            key (str): Its key.
            siblings (dict[str, Any]): The dictionary that holds it, whose "_fmt" and
                "_list" keys give the text of a time and the order of a set.

        Returns:
            Any: The value, with dates and times as text, sets and tuples as lists, and
            dictionaries converted by :meth:`_label_as_json`.
        """

        if isinstance(value, dict):
            return NewHorizonsLORRI._label_as_json(value)
        if isinstance(value, (dt.date, dt.time)):           # a datetime is also a date
            return siblings.get(key + '_fmt', value.isoformat())
        if isinstance(value, (set, frozenset)):
            return list(siblings.get(key + '_list', sorted(value, key=str)))
        if isinstance(value, (list, tuple)):
            return [NewHorizonsLORRI._json_value(item, '', {}) for item in value]
        return value

    @staticmethod
    def _header_as_json(header):
        """A FITS header as a dictionary that can be written as JSON.

        Parameters:
            header (astropy.io.fits.Header): The header.

        Returns:
            dict[str, Any]: Every card with a value, keyed as :meth:`_header_cards` keys
            it, and the text of the blank, COMMENT and HISTORY cards as a list under
            each of those keywords.
        """

        result = {}
        for (key, value) in NewHorizonsLORRI._header_cards(header):
            if key in _COMMENTARY_KEYWORDS:
                result.setdefault(key, []).append(str(value))
            elif value is not None:
                result[key] = value
        return result

    @staticmethod
    def _header_cards(header):
        """The cards of a FITS header, keyed so that a repeated keyword keeps every value.

        Parameters:
            header (astropy.io.fits.Header): The header.

        Returns:
            list[tuple[str, Any]]: The keyword and value of each card in header order. A
            keyword other than blank, COMMENT and HISTORY that occurs more than once is
            given at each occurrence as the keyword followed by its occurrence counted
            from zero in brackets, such as "REFFILE[1]", the form the VICAR hosts use for
            a repeated parameter.
        """

        counts = collections.Counter(header.keys())
        seen = collections.Counter()
        cards = []
        for (key, value) in header.items():
            if key not in _COMMENTARY_KEYWORDS and counts[key] > 1:
                occurrence = seen[key]
                seen[key] += 1
                key = f'{key}[{occurrence}]'
            cards.append((key, value))
        return cards

    @staticmethod
    def _row_as_json(row_dict):
        """A row of an index as a dictionary that can be written as JSON.

        Parameters:
            row_dict (dict[str, Any]): The row, as the table reader gives it, with a key
                ending in "_mask" beside a column that is True where the column's value
                equals the null constant that the column declares.

        Returns:
            dict[str, Any]: The row without its masks, with each NumPy value converted to
            the Python type that holds it, and without any value that counts as missing
            or is a null value that the reader masks. A masked value is a null value
            unless it is a zero that is real: the target's vectors and distances are
            missing only where the vector from the spacecraft to the target is zero
            (:data:`_TARGET_GEOMETRY`), and the components of any other vector only
            where all of them are zero, including a component whose column does not
            declare a null constant.
        """

        values = {key: (value.item() if isinstance(value, np.generic) else value)
                  for (key, value) in row_dict.items() if not key.endswith('_mask')}
        no_target = all(values.get(key, 0.) == 0. for key in _TARGET_POSITION)

        vectors = {}
        for key in values:
            if match := _VECTOR_COMPONENT.fullmatch(key):
                vectors.setdefault(match[1] or match[2], []).append(key)

        missing = set()
        for key in values:
            if not row_dict.get(key + '_mask', False):
                continue
            if values[key] != 0:                    # -999 or -999999
                missing.add(key)
            elif _TARGET_GEOMETRY.fullmatch(key):
                if no_target:
                    missing.add(key)
            elif match := _VECTOR_COMPONENT.fullmatch(key):
                components = vectors[match[1] or match[2]]
                if all(values[component] == 0 for component in components):
                    missing.update(components)
            else:
                missing.add(key)

        return {key: value for (key, value) in values.items()
                if key not in missing and not NewHorizonsLORRI._is_missing(value)}

    @staticmethod
    def _add_start_time_from_clock(items, filepath):
        """Give the merged items a START_TIME from the clock count if the sources lack it.

        Parameters:
            items (dict[str, Any]): The merged metadata items, changed in place.
            filepath (FCPath | None): The path to the FITS file, for the message.

        Returns:
            bool: True if START_TIME was converted from the clock count; False if a
            source gives it, or if the sources lack the clock count too.

        Raises:
            oops.OopsSpiceError: If the spacecraft clock kernel cannot convert the
                SPACECRAFT_CLOCK_START_COUNT.
        """

        if 'START_TIME' in items or 'SPACECRAFT_CLOCK_START_COUNT' not in items:
            return False

        count = items['SPACECRAFT_CLOCK_START_COUNT']
        try:
            items['START_TIME'] = cspyce.scs2e(_NewHorizons.SPACECRAFT_ID, count)
        except oops.OopsSpiceError.CSPYCE_ERRORS as err:
            if not oops.OopsSpiceError.is_spice_error(err):
                raise
            raise oops.OopsSpiceError.from_spice_error(
                err, description=f'no time for the LORRI clock count {count}',
                filepath=filepath) from err
        return True

    @staticmethod
    def _iof_calibrations(header_dict, texp, *, target, tdb):
        """The I/F calibration of a calibrated product, if it can be built.

        Parameters:
            header_dict (dict[str, Any]): The primary FITS header as a dictionary.
            texp (float): The exposure duration in seconds.
            target (str | oops.Body): The target of the observation.
            tdb (float): The middle of the exposure in seconds TDB.

        Returns:
            list[oops.Calibration]: A FlatCalib named "I/F" whose factor is
            pi r^2 / (texp RSOLAR F), in which r is the distance of the target from the
            Sun at `tdb` in AU from the SP kernels and F the solar flux at 1 AU; empty
            if the header lacks RSOLAR, as for a raw product, if the exposure is
            zero, if the target has no path or the kernels do not give its position at
            `tdb`, or if the target is the Sun, for which I/F has no meaning.
        """

        rsolar = header_dict.get('RSOLAR')
        if not isinstance(rsolar, (int, float)) or rsolar <= 0. or texp <= 0.:
            return []

        try:
            path = oops.Path.as_waypoint(target.path if isinstance(target, oops.Body)
                                         else target)
        except KeyError:            # a target that lacks a path
            return []

        try:
            event = path.wrt('SUN').event_at_time(tdb)
        except oops.OopsSpiceError.CSPYCE_ERRORS as err:
            if not oops.OopsSpiceError.is_spice_error(err):
                raise
            return []

        solar_range = event.pos.norm().vals / _AU
        if solar_range == 0.:
            return []
        factor = np.pi * solar_range**2 / (texp * rsolar * _SOLAR_FLUX)
        return [oops.calibration.FlatCalib('I/F', factor)]

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
            label gives inside its IMAGE object; START_TIME is an ISO string to the
            millisecond, as the label writes it.
        """

        dict_ = {key: value for (key, value) in items.items() if key != 'LINES'}
        tai = julian.tai_from_tdb(dict_['START_TIME'])
        dict_['START_TIME'] = julian.iso_from_tai(tai, digits=_TIME_DIGITS)
        return dict_

    @staticmethod
    def _insert_common_subfields(obs, items, *, binning_mode, target, ck, time_from_sclk,
                                 return_all_planets):
        """Insert the subfields that file and index observations share.

        Parameters:
            obs (Snapshot): The observation.
            items (dict[str, Any]): The metadata items; TARGET_NAME and the clock counts
                become subfields where they are given.
            binning_mode (str): "1X1" or "4X4".
            target (str | oops.Body): The target of the observation, after any `target`
                override, whose kernels are included in `spice_kernels`.
            ck (bool): True if the pointing comes from the C kernels, so that they are
                included in `spice_kernels`.
            time_from_sclk (bool): True if the start time was converted from the clock
                count, which becomes the `time_from_sclk` subfield.
            return_all_planets (bool): True to include the kernels of every planet and
                of Pluto in `spice_kernels`.
        """

        obs.insert_subfield('host', 'New Horizons')
        obs.insert_subfield('instrument', 'LORRI')
        obs.insert_subfield('binning_mode', binning_mode)
        obs.insert_subfield('time_from_sclk', time_from_sclk)
        for (subfield, name) in [('label_target', 'TARGET_NAME'),
                                 ('spacecraft_clock_start_count',
                                  'SPACECRAFT_CLOCK_START_COUNT'),
                                 ('spacecraft_clock_stop_count',
                                  'SPACECRAFT_CLOCK_STOP_COUNT')]:
            if name in items:
                obs.insert_subfield(subfield, items[name])

        # With a custom frame, pointing never came from a CK, so any CK that happens to be
        # furnished is unrelated to this observation and must not be reported as used.
        # A `target` override can be a Body, whose kernels are found by its name.
        kernels = _NewHorizons.used_kernels(obs.time, 'lorri',
                                            target=getattr(target, 'name', target), ck=ck,
                                            return_all_planets=return_all_planets)
        obs.insert_subfield('spice_kernels', kernels)
        obs.insert_subfield('spice_to_frame', _SPICE_TO_FRAME)
        obs.insert_subfield('spice_frame_name', _SPICE_FRAME_NAME)
        obs.insert_subfield('spice_frame_id', _SPICE_FRAME_ID)

    @staticmethod
    def _target_body(name):
        """The name of the target body for a TARGET_NAME, by a fixed translation.

        The label gives the PDS catalog's form of a name, which is not always the name
        SPICE knows the body by: "J17 CALLIRRHOE" for CALLIRRHOE, "136472 MAKEMAKE" for
        MAKEMAKE, "ASTEROID 486958 (2014 MU69)" for ARROKOTH. A label or index gives
        "CALIBRATION" for the calibration images of the launch phase, whose FITS header
        gives the null value "---" as the target. Some images of Pluto give a name that
        begins with "PLUTO".

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
    def _fits_path_from_header(header, id_name, tdb_midtime):
        """The spacecraft's path as the FITS header gives it.

        Parameters:
            header (astropy.io.fits.Header): The primary header.
            id_name (str): The name that makes the path's ID unique, as
                :meth:`_make_snapshot` gives it.
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
                                       path_id='.NH_PATH_' + id_name)
        return oops.Path.as_waypoint(sc_path)

    @staticmethod
    def _fits_frame_from_header(header, id_name, *, pointing, path, tdb_midtime):
        """The camera frame as the FITS header gives it.

        Parameters:
            header (astropy.io.fits.Header): The primary header.
            id_name (str): The name that makes the frame's ID unique, as
                :meth:`_make_snapshot` gives it.
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
                                               frame_id='.NH_FRAME_' + id_name)
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
            asof (str, optional): A date in ISO format. Each kernel is chosen from
                those released no later than the start (00:00) of that date, the leap
                seconds kernel included, and a query of which every kernel was released
                after the start of the date silently uses the kernel it would choose
                without a date. A set of SP or C kernels keeps only those released by
                the start of the date, so the later part of the mission can lack
                coverage. The C kernels and a metakernel named in the database never fall
                back, and ValueError is raised instead. The LORRI instrument kernel is
                chosen from the database by this date, also with `meta`. None to choose
                the newest kernels.
            meta (str, optional): The name of a metakernel in the database, or the path
                to a metakernel file, to furnish instead of the kernels selected from the
                database. A metakernel named in the database is the newest one released
                by the start of `asof`; a metakernel file is used whatever the date. The
                kernels that it lists are furnished whatever their release dates, after
                the Solar System's, so a leap seconds kernel among them takes the place
                of the one that `asof` chose.
            planets (int | collections.abc.Sequence[int], optional): The planet number,
                1 to 9, or a sequence of planet numbers, such as (5, 9) for Jupiter and
                Pluto, whose bodies :meth:`~oops.Body.define_solar_system` defines with
                their moons. None or 0 means all nine.
            irregulars (bool, optional): True to include the irregular satellites; False
                otherwise.

        Raises:
            oops.host.HostError: If the instrument kernel does not give a 1024x1024 image.
            RuntimeError: If Saturn is among the planets and `asof` is a date from
                2003-11-05 through 2004-04-26, because the frames kernel of Saturn's
                small moons chosen for those dates has a malformed line that SPICE
                rejects.
            ValueError: If `meta` is None and none of the spacecraft's C kernels was
                released by the start of `asof`, as for any date up to and including
                2015-03-12 outside that window; or if `meta` names a metakernel in the
                database and none of that name was released by the start of `asof`.
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
