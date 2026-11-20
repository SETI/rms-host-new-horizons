##########################################################################################
# tests/oops/test_lorri.py
##########################################################################################
"""Tests for the New Horizons LORRI host, host_new_horizons.oops.lorri."""

from collections.abc import Iterator
import json
import pathlib
import tempfile
from typing import Any

import cspyce
import julian
import numpy as np
import oops
import pytest
from astropy.io import fits
from conftest import (FITS_HEADER, TEST_FILES, FakeSpice, FakeSpiceDB, cache_file,
                      write_fits, write_label)
from filecache import FCPath
from pdsparser import Pds3Label

import host_new_horizons
import host_new_horizons.oops
from host_new_horizons.oops.lorri import NewHorizonsLORRI

INDEX = TEST_FILES / 'NHJULO_2001_index_edited.lbl'
SUPPLEMENT = TEST_FILES / 'NHJULO_2001_supplemental_index_edited.lbl'

# The middle and length of the exposure in the synthetic FITS header
MIDTIME = FITS_HEADER['SPCSCET']
TEXP = FITS_HEADER['EXPTIME']


def load(tmp_path: pathlib.Path, *, fits_keywords: dict[str, Any] | None = None,
         **options: Any) -> Any:
    """The Snapshot of a synthetic FITS file written with the given header keywords.

    The file is written in a new directory under `tmp_path`, so one test can load more
    than one.
    """

    directory = pathlib.Path(tempfile.mkdtemp(dir=tmp_path))
    path = write_fits(directory / 'lor_test.fit', **(fits_keywords or {}))
    return NewHorizonsLORRI.from_file(path, **options)


def vals(value: Any) -> Any:
    """The values of a polymath object as nested lists."""

    return np.asarray(value.vals).tolist()


def boresight(obs: Any) -> Any:
    """The J2000 direction of the optic axis at mid-exposure."""

    xform = obs.frame.wrt(oops.Frame.J2000).transform_at_time(obs.midtime)
    return xform.unrotate(oops.Vector3.ZAXIS)


def row(**columns: Any) -> dict[str, Any]:
    """An index row merged with its supplement, with the given columns replaced."""

    row_dict = {
        'VOLUME_ID'                   : 'NHJULO_2001',
        'FILE_SPECIFICATION_NAME'     : 'data/20070110_003071/'
                                        'lor_0030713591_0x633_sci.lbl',
        'START_TIME'                  : '2007-01-10T05:41:08.363',
        'EXPOSURE_DURATION'           : 5.,
        'BINNING_MODE'                : '4x4',
        'TARGET_NAME'                 : 'J17 CALLIRRHOE',
        'SPACECRAFT_CLOCK_START_COUNT': '0030713585:49000',
        'SPACECRAFT_CLOCK_STOP_COUNT' : '0030713590:49000',
    }
    row_dict.update(columns)
    return row_dict

##########################################################################################
# Registration and exports
##########################################################################################


def test_lorri_is_registered_with_oops() -> None:
    assert NewHorizonsLORRI in oops.Host._HOSTS


def test_module_exports() -> None:
    assert host_new_horizons.oops.__all__ == ['NewHorizonsLORRI']


@pytest.mark.parametrize('module', [host_new_horizons, host_new_horizons.oops])
def test_package_exposes_host(module: Any) -> None:
    assert module.NewHorizonsLORRI is NewHorizonsLORRI

##########################################################################################
# Detectors
##########################################################################################


@pytest.mark.parametrize(('label', 'expected'), [
    ({'INSTRUMENT_HOST_NAME': 'NEW HORIZONS', 'INSTRUMENT_ID': 'LORRI'}, True),
    ({'INSTRUMENT_HOST_NAME': 'NEW HORIZONS', 'INSTRUMENT_ID': 'MVIC'}, False),
    ({'INSTRUMENT_HOST_NAME': 'CASSINI ORBITER', 'INSTRUMENT_ID': 'LORRI'}, False),
    ({}, False),
])
def test_detect_in_pds3(label: dict[str, str], expected: bool) -> None:
    assert NewHorizonsLORRI._detect_in_pds3(label) is expected


@pytest.mark.parametrize(('label', 'expected'), [
    ({'INSTRUMENT_HOST_NAME': 'NEW HORIZONS', 'INSTRUMENT_ID': 'LORRI'}, True),
    ({'INSTRUMENT_HOST_NAME': 'GALILEO ORBITER', 'INSTRUMENT_ID': 'SSI'}, False),
])
def test_detect_in_index(label: dict[str, str], expected: bool) -> None:
    assert NewHorizonsLORRI._detect_in_index(label) is expected


@pytest.mark.parametrize(('keywords', 'expected'), [
    ({'MISSION': 'New Horizons', 'INSTRU': 'lor'}, True),
    ({'MISSION': 'NEW HORIZONS', 'INSTRU': 'LOR     '}, True),
    ({'MISSION': 'New Horizons', 'INSTRU': 'mvi'}, False),
    ({'MISSION': 'Cassini', 'INSTRU': 'lor'}, False),
    ({}, False),
])
def test_detect_in_fits(keywords: dict[str, str], expected: bool) -> None:
    primary = fits.PrimaryHDU()
    for (key, value) in keywords.items():
        primary.header[key] = value
    assert NewHorizonsLORRI._detect_in_fits(fits.HDUList([primary])) is expected


@pytest.mark.parametrize(('row_dict', 'expected'), [
    ({'INSTRUMENT_HOST': 'NEW HORIZONS', 'INSTRUMENT_ID': 'LORRI'}, True),
    ({'INSTRUMENT_HOST': 'NEW HORIZONS', 'INSTRUMENT_ID': 'ALICE'}, False),
    ({}, False),
])
def test_detect_in_row(row_dict: dict[str, str], expected: bool) -> None:
    assert NewHorizonsLORRI._detect_in_row(row_dict) is expected

##########################################################################################
# initialize()
##########################################################################################


def test_initialize_loads_the_mission(fake_spice: FakeSpice,
                                      fake_spicedb: FakeSpiceDB) -> None:
    NewHorizonsLORRI.initialize(asof='2020-01-01', planets=[9], irregulars=False)
    assert fake_spice.solar_system[0][1] == {'asof': '2020-01-01T00:00:00',
                                             'planets': [9], 'irregulars': False}
    assert fake_spicedb.named('select_inst')[0][1]['asof'] == '2020-01-01T00:00:00'


def test_initialize_with_a_metakernel(fake_spicedb: FakeSpiceDB) -> None:
    NewHorizonsLORRI.initialize(meta='nh.tm')
    assert fake_spicedb.named('furnish_by_metafile')[0][0] == ('nh.tm',)


def test_initialize_rejects_an_unexpected_image_shape(fake_spicedb: FakeSpiceDB) -> None:
    fake_spicedb.kernel_dict = {'INS': {'NH_LORRI_1X1': {'PIXEL_SAMPLES': 1024,
                                                         'PIXEL_LINES': 512}}}
    with pytest.raises(oops.host.HostError,
                       match=r'unsupported LORRI image shape \(1024, 512\)'):
        NewHorizonsLORRI.initialize()


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    NewHorizonsLORRI.initialize()
    NewHorizonsLORRI.initialize()
    assert len(fake_spice.solar_system) == 1


def test_initialize_can_repeat_after_reset(fake_spice: FakeSpice) -> None:
    NewHorizonsLORRI.initialize()
    NewHorizonsLORRI._reset()
    NewHorizonsLORRI.initialize()
    assert len(fake_spice.solar_system) == 2


def test_reset_clears_state() -> None:
    NewHorizonsLORRI.initialize()
    NewHorizonsLORRI._reset()
    assert NewHorizonsLORRI._INSTRUMENT_KERNEL is None
    assert NewHorizonsLORRI._FOVS == {}
    assert not NewHorizonsLORRI._initialized


def test_fov_keys() -> None:
    NewHorizonsLORRI.initialize()
    assert set(NewHorizonsLORRI._FOVS) == {(binning, fov_type)
                                           for binning in ('1X1', '4X4')
                                           for fov_type in ('fast', 'slow', 'flat')}


@pytest.mark.parametrize(('fov_type', 'inverse'), [('fast', True), ('slow', False)])
def test_polynomial_fovs(fov_type: str, inverse: bool) -> None:
    NewHorizonsLORRI.initialize()
    fov = NewHorizonsLORRI._FOVS['1X1', fov_type]
    assert isinstance(fov, oops.fov.PolynomialFOV)
    assert vals(fov.uv_shape) == [1024, 1024]
    assert (fov.coefft_xy_from_uv is not None) is inverse


def test_polynomial_fov_distortion_coefficients() -> None:
    NewHorizonsLORRI.initialize()
    coefft = NewHorizonsLORRI._FOVS['1X1', 'fast'].coefft_uv_from_xy
    assert coefft[1, 0, 0] == pytest.approx(76.9231 * 2619.008)
    assert coefft[0, 1, 1] == pytest.approx(-76.9231 * 2619.008)
    assert coefft[3, 0, 0] == pytest.approx(76.9231 * 2.696e-5 * 2619.008**3)
    assert coefft[1, 1, 1] == pytest.approx(-76.9231 * -2.864e-5 * 2619.008**2)


def test_fast_fov_inverts_the_distortion() -> None:
    NewHorizonsLORRI.initialize()
    fov = NewHorizonsLORRI._FOVS['1X1', 'fast']
    uv = oops.Pair([(100., 200.), (900., 50.), (512., 512.)])
    roundtrip = fov.uv_from_los(fov.los_from_uv(uv))
    assert roundtrip.vals == pytest.approx(uv.vals, abs=0.003)


def test_flat_fov() -> None:
    NewHorizonsLORRI.initialize()
    fov = NewHorizonsLORRI._FOVS['1X1', 'flat']
    assert isinstance(fov, oops.fov.FlatFOV)
    assert vals(fov.uv_scale) == [4.96e-6, -4.96e-6]


@pytest.mark.parametrize('fov_type', ['fast', 'slow', 'flat'])
def test_binned_fovs(fov_type: str) -> None:
    NewHorizonsLORRI.initialize()
    fov = NewHorizonsLORRI._FOVS['4X4', fov_type]
    assert isinstance(fov, oops.fov.SubsampledFOV)
    assert fov.fov is NewHorizonsLORRI._FOVS['1X1', fov_type]
    assert vals(fov.uv_shape) == [256, 256]

##########################################################################################
# The camera frame
##########################################################################################


def test_camera_frame_flips_the_spice_frame(fake_spice: FakeSpice) -> None:
    NewHorizonsLORRI._define_camera_frame()
    assert fake_spice.spice_frames == [('NH_LORRI', 'NH_LORRI_FLIPPED')]
    xform = oops.Frame.as_frame('NH_LORRI').wrt(oops.Frame.J2000).transform_at_time(0.)
    assert vals(xform.unrotate(oops.Vector3.ZAXIS)) == [0., 0., -1.]


def test_camera_frame_is_defined_once(fake_spice: FakeSpice) -> None:
    NewHorizonsLORRI._define_camera_frame()
    NewHorizonsLORRI._define_camera_frame()
    assert len(fake_spice.spice_frames) == 1

##########################################################################################
# from_file() on a FITS file
##########################################################################################


def test_from_file_time(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    assert obs.time == pytest.approx((MIDTIME - TEXP / 2., MIDTIME + TEXP / 2.))
    assert obs.texp == TEXP


@pytest.mark.parametrize(('lines', 'binning'), [(1024, '1X1'), (256, '4X4')])
def test_from_file_binning(tmp_path: pathlib.Path, lines: int, binning: str) -> None:
    path = write_fits(tmp_path / 'lor_test.fit', lines=lines, samples=lines)
    obs = NewHorizonsLORRI.from_file(path)
    assert obs.binning_mode == binning
    assert obs.fov is NewHorizonsLORRI._FOVS[binning, 'fast']


def test_from_file_rejects_an_unknown_shape(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit', lines=512, samples=512)
    with pytest.raises(oops.host.HostFormatError,
                       match='unrecognized LORRI image shape: 512 lines'):
        NewHorizonsLORRI.from_file(path)


@pytest.mark.parametrize('fov_type', ['slow', 'flat'])
def test_from_file_fov_type(tmp_path: pathlib.Path, fov_type: str) -> None:
    obs = load(tmp_path, fov_type=fov_type)
    assert obs.fov is NewHorizonsLORRI._FOVS['1X1', fov_type]


@pytest.mark.parametrize(('option', 'value'), [
    ('geom', 'header'),
    ('pointing', 'wcs'),
    ('fov_type', 'fastest'),
])
def test_from_file_rejects_an_invalid_option(tmp_path: pathlib.Path, option: str,
                                             value: str) -> None:
    options: dict[str, Any] = {option: value}
    with pytest.raises(oops.OopsValueError, match=f'invalid {option} {value!r}; use one'):
        load(tmp_path, **options)


# The arguments are checked before the file is opened, so the file need not exist
@pytest.mark.parametrize('select', [0, slice(None), 'data'])
def test_from_file_refuses_select(tmp_path: pathlib.Path, select: Any) -> None:
    with pytest.raises(oops.OopsValueError,
                       match=r'each New Horizons LORRI file is one observation; select '
                             r'must be None, not '):
        NewHorizonsLORRI.from_file(tmp_path / 'missing.fit', select=select)


def test_from_file_rejects_an_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.OopsTypeError,
                       match=r"unexpected keyword argument for "
                             r"NewHorizonsLORRI\.from_file: 'calibration'"):
        NewHorizonsLORRI.from_file(tmp_path / 'missing.fit', calibration=[])


def test_from_file_data(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    expected = np.arange(1024 * 1024, dtype='>f4').reshape(1024, 1024)
    assert obs.data.dtype == np.dtype('>f4')
    assert np.array_equal(obs.data, expected)


def test_from_file_error_and_quality(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    assert np.array_equal(obs.error, obs.data + 0.5)
    assert np.array_equal(obs.quality, np.ones((1024, 1024)))


def test_from_file_raw_product(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit', lines=256, samples=257,
                      calibrated=False)
    obs = NewHorizonsLORRI.from_file(path)
    assert obs.data.shape == (256, 257)
    assert 'error' not in obs.subfields
    assert 'quality' not in obs.subfields


def test_from_file_astrometry(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, astrometry=True)
    assert 'data' not in obs.subfields
    assert 'error' not in obs.subfields
    assert 'quality' not in obs.subfields
    assert obs.headers[0]['SPCCBTNM'] == 'JUPITER'


def test_from_file_header_dict_is_the_primary_header(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    with fits.open(path, mode='update') as hdulist:
        hdulist[0].header['COMMENT'] = 'a comment'
        hdulist[0].header['HISTORY'] = 'a history'
    obs = NewHorizonsLORRI.from_file(path)
    assert obs.header_dict['SPCCBTNM'] == 'JUPITER'
    assert 'COMMENT' not in obs.header_dict
    assert 'HISTORY' not in obs.header_dict


@pytest.mark.parametrize(('key', 'value'), [
    ('START_TIME', '2007-01-10T05:41:08.363'),
    ('EXPOSURE_DURATION', 5.0),
    ('TARGET_NAME', 'JUPITER'),
    ('SPACECRAFT_CLOCK_START_COUNT', '0030713585:49000'),
    ('SPACECRAFT_CLOCK_STOP_COUNT', '0030713590:49000'),
    ('NAXIS2', 1024),
    ('MISSION', 'New Horizons'),
])
def test_from_file_dict_gives_pds3_names(tmp_path: pathlib.Path, key: str,
                                         value: Any) -> None:
    assert load(tmp_path).dict[key] == value


@pytest.mark.parametrize('keyword', ['SPCSCET', 'EXPTIME', 'SPCCBTNM', 'SPCSCLK0',
                                     'STOPMET', 'COMMENT'])
def test_from_file_dict_leaves_out_converted_keywords(tmp_path: pathlib.Path,
                                                      keyword: str) -> None:
    assert keyword not in load(tmp_path).dict


def test_from_file_dict_leaves_out_null_values(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, fits_keywords={'CB3LOAD': '---'})
    assert 'CB3LOAD' not in obs.dict
    assert obs.header_dict['CB3LOAD'] == '---'


# The FITS header writes -1e32 for a value it cannot compute, such as the Sun's distance
# from a target it has no ephemeris for
def test_from_file_dict_leaves_out_the_fits_fill_value(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, fits_keywords={'SPCTSORN': -1.e32})
    assert 'SPCTSORN' not in obs.dict
    assert obs.header_dict['SPCTSORN'] == -1.e32


def test_from_file_headers(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    assert [header.get('EXTNAME') for header in obs.headers] == [
        None, 'LORRI Error image', 'LORRI Quality flag image']


def test_from_file_subfields(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    assert obs.host == 'New Horizons'
    assert obs.instrument == 'LORRI'
    assert obs.label_target == 'JUPITER'
    assert obs.target == 'JUPITER'
    assert obs.lightsource == 'SUN'
    assert obs.calibrations == []

##########################################################################################
# Calibration
##########################################################################################

# The conversion of a calibrated image to radiance for a solar spectrum, as the FITS
# header of a calibrated 4x4 image gives it, in (DN/s/pixel)/(erg/cm^2/s/A/sr)
RSOLAR = 5114880.0

# Jupiter's distance from the Sun in the jupiter_path and sun_path fixtures, in AU
JUPITER_RANGE = 1.e9 / 149597870.7


def test_from_file_iof_calibration(tmp_path: pathlib.Path, sun_path: oops.Path,
                                   jupiter_path: oops.Path) -> None:
    """I/F is pi r^2 DN / (texp RSOLAR F), where F is the solar flux at 1 AU, 176
    erg/cm^2/s/A at the pivot wavelength of LORRI, and r the target's solar range in
    AU."""

    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR})
    (calibration,) = obs.calibrations
    assert calibration.name == 'I/F'
    assert calibration.factor.vals == pytest.approx(
        np.pi * JUPITER_RANGE**2 / (TEXP * RSOLAR * 176.), rel=1.e-12)


def test_from_file_calibration_follows_the_target(tmp_path: pathlib.Path,
                                                  sun_path: oops.Path,
                                                  jupiter_path: oops.Path) -> None:
    oops.path.FixedPath((2.e9, 0., 0.), oops.Path.SSB, oops.Frame.J2000, path_id='SATURN')
    jupiter = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR})
    saturn = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR}, target='SATURN')
    assert saturn.calibrations[0].factor.vals == pytest.approx(
        4. * jupiter.calibrations[0].factor.vals)


# The target recedes from the Sun at 1e6 km/s, so it is twice as far 1000 s later
def test_from_file_calibration_at_the_cadence_override(tmp_path: pathlib.Path,
                                                       sun_path: oops.Path) -> None:
    oops.path.LinearPath((oops.Vector3((1.e9, 0., 0.)), oops.Vector3((1.e6, 0., 0.))),
                         MIDTIME, 'SSB', frame='J2000', path_id='JUPITER')
    default = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR})
    later = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR},
                 cadence=oops.cadence.SnapCadence(MIDTIME + 1000. - TEXP / 2., TEXP))
    assert later.calibrations[0].factor.vals == pytest.approx(
        4. * default.calibrations[0].factor.vals)


def test_raw_product_has_no_calibration(tmp_path: pathlib.Path, sun_path: oops.Path,
                                        jupiter_path: oops.Path) -> None:
    obs = load(tmp_path, fits_keywords={'RSOLAR': None})
    assert obs.calibrations == []


def test_label_astrometry_has_no_calibration(tmp_path: pathlib.Path, sun_path: oops.Path,
                                             jupiter_path: oops.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit', RSOLAR=RSOLAR)
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    assert NewHorizonsLORRI.from_file(label, astrometry=True).calibrations == []


def test_zero_exposure_has_no_calibration(tmp_path: pathlib.Path, sun_path: oops.Path,
                                          jupiter_path: oops.Path) -> None:
    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR, 'EXPTIME': 0.})
    assert obs.calibrations == []


# I/F is not defined for the Sun itself, which is zero AU from the Sun
def test_image_of_the_sun_has_no_calibration(tmp_path: pathlib.Path,
                                             sun_path: oops.Path) -> None:
    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR, 'SPCCBTNM': 'SUN'})
    assert obs.calibrations == []


# The target's solar range cannot be found for a missing target, or for one without a path
@pytest.mark.parametrize('target', ['---', 'CALLIRRHOE'])
def test_target_without_a_path_has_no_calibration(tmp_path: pathlib.Path,
                                                  sun_path: oops.Path,
                                                  target: str) -> None:
    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR, 'SPCCBTNM': target})
    assert obs.calibrations == []


def test_target_without_spice_coverage_has_no_calibration(
        tmp_path: pathlib.Path, sun_path: oops.Path, jupiter_path: oops.Path,
        monkeypatch: pytest.MonkeyPatch) -> None:
    def event_at_time(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError('SPICE(SPKINSUFFDATA) -- Insufficient ephemeris data')

    monkeypatch.setattr(oops.path.FixedPath, 'event_at_time', event_at_time)
    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR})
    assert obs.calibrations == []


def test_calibration_error_without_a_spice_code_passes_through(
        tmp_path: pathlib.Path, sun_path: oops.Path, jupiter_path: oops.Path,
        monkeypatch: pytest.MonkeyPatch) -> None:
    def event_at_time(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError('something else')

    monkeypatch.setattr(oops.path.FixedPath, 'event_at_time', event_at_time)
    with pytest.raises(RuntimeError, match=r'^something else$'):
        load(tmp_path, fits_keywords={'RSOLAR': RSOLAR})


def test_calibrations_override(tmp_path: pathlib.Path, sun_path: oops.Path,
                               jupiter_path: oops.Path) -> None:
    calibration = oops.calibration.FlatCalib('DN', 1.)
    obs = load(tmp_path, fits_keywords={'RSOLAR': RSOLAR}, calibrations=[calibration])
    assert obs.calibrations == [calibration]


def test_from_index_has_no_calibration() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(RSOLAR=RSOLAR)])
    assert obs.calibrations == []


@pytest.mark.parametrize(('label_target', 'target'), [
    ('PLUTO', 'PLUTO'),
    ('PLUTO-CHARON', 'PLUTO'),
    ('CHARON', 'CHARON'),
])
def test_from_file_target(tmp_path: pathlib.Path, label_target: str, target: str) -> None:
    obs = load(tmp_path, fits_keywords={'SPCCBTNM': label_target})
    assert obs.target == target
    assert obs.label_target == label_target


# The target body of each form of name that the labels give
@pytest.mark.parametrize(('name', 'body'), [
    ('J1 IO', 'IO'),
    ('J17 CALLIRRHOE', 'CALLIRRHOE'),
    ('J RINGS', 'J RINGS'),
    ('ASTEROID 486958 (2014 MU69)', 'ARROKOTH'),
    ('MU69', 'ARROKOTH'),
    ('ASTEROID 516977 (2012 HZ84)', 'ASTEROID 516977 (2012 HZ84)'),
    ('136472 MAKEMAKE', 'MAKEMAKE'),
    ('50000 QUAOAR (2002 LM60)', 'QUAOAR'),
    ('15810 ARAWN (1994 JR1)', 'ARAWN'),
    ('2060 CHIRON', 'CHIRON'),
    ('2011 HJ103', '2011 HJ103'),
    ('2020 KR11', '2020 KR11'),
    ('2014 MU', '2014 MU'),
    ('NGC 3532', 'NGC 3532'),
    ('WOLF 359', 'WOLF 359'),
    ('RADEC=270.66,-14.63=P4', 'RADEC=270.66,-14.63=P4'),
    ('PLUTO', 'PLUTO'),
    ('PLUTO-CHARON', 'PLUTO'),
    ('CALIBRATION', 'NONE'),
    (None, 'NONE'),
])
def test_target_body(name: str | None, body: str) -> None:
    """Each form of name translates by its fixed rule; any other name stays as it is."""

    assert NewHorizonsLORRI._target_body(name) == body


def test_from_file_target_comes_from_the_label_name(tmp_path: pathlib.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit', SPCCBTNM='MU69')
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT',
                        TARGET_NAME='"ASTEROID 486958 (2014 MU69)"')
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.target == 'ARROKOTH'
    assert obs.label_target == 'ASTEROID 486958 (2014 MU69)'


def test_from_file_spice_kernels_of_the_target_body(tmp_path: pathlib.Path,
                                                    fake_spicedb: FakeSpiceDB) -> None:
    write_fits(tmp_path / 'lor_test.fit', SPCCBTNM='CALLIRRHOE')
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT',
                        TARGET_NAME='"J17 CALLIRRHOE"')
    _ = NewHorizonsLORRI.from_file(label)
    (_, kwargs) = fake_spicedb.named('used_basenames')[0]
    assert kwargs['bodies'] == [5, 517]


@pytest.mark.parametrize('null', ['---', 'N/A', 'UNK', 'NULL', 'UNKNOWN', ''])
def test_from_file_null_target(tmp_path: pathlib.Path, null: str) -> None:
    obs = load(tmp_path, fits_keywords={'SPCCBTNM': null})
    assert obs.target == 'NONE'
    assert 'label_target' not in obs.subfields
    assert 'TARGET_NAME' not in obs.dict


def test_from_file_clock_counts(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, fits_keywords={'STOPMET': 30713586.48})
    assert obs.spacecraft_clock_start_count == '0030713585:49000'
    assert obs.spacecraft_clock_stop_count == '0030713586:24000'


@pytest.mark.parametrize(('keyword', 'subfield', 'key'), [
    ('SPCSCLK0', 'spacecraft_clock_start_count', 'SPACECRAFT_CLOCK_START_COUNT'),
    ('STOPMET', 'spacecraft_clock_stop_count', 'SPACECRAFT_CLOCK_STOP_COUNT'),
])
def test_from_file_missing_clock_count(tmp_path: pathlib.Path, keyword: str,
                                       subfield: str, key: str) -> None:
    obs = load(tmp_path, fits_keywords={keyword: None})
    assert subfield not in obs.subfields
    assert key not in obs.dict


@pytest.mark.parametrize(('met', 'count'), [
    (30713585.98, '0030713585:49000'),
    (30598438.977, '0030598438:48850'),
    (408709280.99999999, '0408709281:00000'),
    (0., '0000000000:00000'),
    # 299193745.002 * 50000 is 14959687250099.998 in floating point: rounded, not
    # truncated
    (299193745.002, '0299193745:00100'),
])
def test_sclk_from_met(met: float, count: str) -> None:
    assert NewHorizonsLORRI._sclk_from_met(met) == count


def test_from_file_path_subfields(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    obs = NewHorizonsLORRI.from_file(path)
    assert obs.filepath == FCPath(path)
    assert obs.basename == 'lor_test.fit'
    assert obs.abspath == path.resolve()
    assert obs.image_url == path.absolute().as_posix()


def test_from_file_spice_subfields(tmp_path: pathlib.Path,
                                   fake_spicedb: FakeSpiceDB) -> None:
    obs = load(tmp_path)
    assert obs.spice_frame_name == 'NH_LORRI'
    assert obs.spice_frame_id == -98300
    assert vals(obs.spice_to_frame) == [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    assert obs.spice_kernels == fake_spicedb.used


def test_from_file_spice_kernels_query(tmp_path: pathlib.Path,
                                       fake_spicedb: FakeSpiceDB) -> None:
    _ = load(tmp_path)
    (_, kwargs) = fake_spicedb.named('used_basenames')[0]
    assert kwargs['types'] is None
    assert kwargs['bodies'] == [5, 599]


# The planets, Pluto included, each with its system barycenter
ALL_PLANETS = [1, 199, 2, 299, 3, 399, 4, 499, 5, 599, 6, 699, 7, 799, 8, 899, 9, 999]


def test_from_file_return_all_planets(tmp_path: pathlib.Path,
                                      fake_spicedb: FakeSpiceDB) -> None:
    """The planets come first, and the target, a moon, is kept after them."""

    _ = load(tmp_path, fits_keywords={'SPCCBTNM': 'CALLIRRHOE'}, return_all_planets=True)
    (_, kwargs) = fake_spicedb.named('used_basenames')[0]
    assert kwargs['bodies'] == ALL_PLANETS + [517]


def test_from_index_return_all_planets(fake_spicedb: FakeSpiceDB) -> None:
    """A target that is one of the planets is not listed twice."""

    _ = list(NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(TARGET_NAME='JUPITER')],
                                         return_all_planets=True))
    (_, kwargs) = fake_spicedb.named('used_basenames')[0]
    assert kwargs['bodies'] == ALL_PLANETS


def test_from_file_default_geometry(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path)
    assert obs.path.path_id == 'NEW HORIZONS'
    assert obs.frame.frame_id == 'NH_LORRI'


def test_from_file_custom_frame_reports_no_cks(tmp_path: pathlib.Path,
                                               fake_spice: FakeSpice,
                                               fake_spicedb: FakeSpiceDB) -> None:
    obs = load(tmp_path, frame='J2000')
    assert obs.frame.frame_id == 'J2000'
    assert fake_spice.spice_frames == []
    assert 'CK' not in fake_spicedb.named('used_basenames')[0][1]['types']


def test_from_file_overrides(tmp_path: pathlib.Path, jupiter_path: oops.Path) -> None:
    lightsource = oops.lightsource.star_lookup('SIRIUS')
    obs = load(tmp_path, target='CHARON', lightsource=lightsource, path='JUPITER')
    assert obs.target == 'CHARON'
    assert obs.lightsource is lightsource
    assert obs.path.path_id == 'JUPITER'


def test_from_file_cadence_override(tmp_path: pathlib.Path) -> None:
    cadence = oops.cadence.SnapCadence(100., 2.)
    obs = load(tmp_path, cadence=cadence)
    assert obs.time == (100., 102.)


def test_from_file_tracker(tmp_path: pathlib.Path, jupiter_path: oops.Path) -> None:
    obs = load(tmp_path, tracker='midtime')
    assert isinstance(obs.frame, oops.frame.TrackerFrame)


def test_from_file_navigation(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, navigation=(0., 0.))
    assert isinstance(obs.frame, oops.frame.Navigation)


def test_from_file_timeshift(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, timeshift=(1., 'T'))
    assert obs.time == pytest.approx((MIDTIME - TEXP / 2. + 1., MIDTIME + TEXP / 2. + 1.))

##########################################################################################
# from_file() with the geometry and pointing of the FITS header
##########################################################################################


def test_fits_geometry_position(tmp_path: pathlib.Path, sun_path: oops.Path) -> None:
    obs = load(tmp_path, geom='fits')
    event = obs.path.wrt('SUN').event_at_time(obs.midtime)
    assert event.pos.vals == pytest.approx([-1.e8, -2.e8, -3.e8])
    assert event.vel.vals == pytest.approx([-1., -2., -3.])


def test_fits_geometry_moves_linearly(tmp_path: pathlib.Path,
                                      sun_path: oops.Path) -> None:
    obs = load(tmp_path, geom='fits')
    event = obs.path.wrt('SUN').event_at_time(obs.midtime + 10.)
    assert event.pos.vals == pytest.approx([-1.e8 - 10., -2.e8 - 20., -3.e8 - 30.])


def test_fits_geometry_path_id(tmp_path: pathlib.Path, sun_path: oops.Path) -> None:
    obs = load(tmp_path, geom='fits')
    assert obs.path.path_id == '.NH_PATH_lor_test.fit'


def test_radec_from_uv_at_the_reference_pixel() -> None:
    (ra, dec) = NewHorizonsLORRI._radec_from_uv(512.5, 512.5, FITS_HEADER)
    assert ra == pytest.approx(237.5)
    assert dec == pytest.approx(-33.8)


def test_radec_from_uv_off_the_reference_pixel() -> None:
    header = FITS_HEADER | {'CRVAL2': 0., 'CD1_2': 0., 'CD2_1': 0.}
    (ra, dec) = NewHorizonsLORRI._radec_from_uv(612.5, 512.5, header)
    assert ra == pytest.approx(237.5 + np.degrees(np.arctan(np.radians(-0.028))))
    assert dec == pytest.approx(0.)


def test_radec_from_uv_beyond_the_pole() -> None:
    header = FITS_HEADER | {'CRVAL2': 89.99, 'CRPIX2': 412.5, 'CD1_1': 2.8e-4,
                            'CD1_2': 0., 'CD2_1': 0., 'CD2_2': 2.8e-4}
    (ra, dec) = NewHorizonsLORRI._radec_from_uv(512.5, 512.5, header)
    assert ra == pytest.approx(237.5 + 180.)
    assert dec == pytest.approx(90. - 0.018)


def test_fits_geometry_with_a_path_override(tmp_path: pathlib.Path,
                                            jupiter_path: oops.Path) -> None:
    keywords = dict.fromkeys(('SPCSSCX', 'SPCSSCY', 'SPCSSCZ'))
    obs = load(tmp_path, geom='fits', path='JUPITER', fits_keywords=keywords)
    assert obs.path.path_id == 'JUPITER'


def test_fits_pointing_with_a_frame_override(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, pointing='fits', frame='J2000', fits_keywords={'SPCEMEN': None})
    assert obs.frame.frame_id == 'J2000'


def test_fitsapp_pointing_uses_the_path_override(tmp_path: pathlib.Path,
                                                 sun_path: oops.Path) -> None:
    # The overriding path moves at 30 km/s, so its aberration turns the boresight; the
    # default path from the header is at rest
    keywords = {'SPCSSCVX': 0., 'SPCSSCVY': 0., 'SPCSSCVZ': 0.}
    moving = oops.path.LinearPath((oops.Vector3((1.e8, 0., 0.)),
                                   oops.Vector3((30., 0., 0.))), MIDTIME, 'SSB',
                                  frame='J2000', path_id='MOVING')
    fits_obs = load(tmp_path, geom='fits', pointing='fits', fits_keywords=keywords)
    app_obs = load(tmp_path, geom='fits', pointing='fitsapp', path=moving,
                   fits_keywords=keywords)
    separation = boresight(app_obs).sep(boresight(fits_obs)).vals
    assert 1.e-5 < separation < 1.e-4


def test_fits_pointing_boresight(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, pointing='fits')
    (ra, dec) = NewHorizonsLORRI._radec_from_uv(512.5, 512.5, FITS_HEADER)
    expected = oops.Vector3.from_ra_dec_length(np.radians(ra), np.radians(dec))
    assert boresight(obs).vals == pytest.approx(expected.vals)


def test_fits_pointing_frame_id(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, pointing='fits')
    assert obs.frame.frame_id == '.NH_FRAME_lor_test.fit'


def test_fits_pointing_reports_no_cks(tmp_path: pathlib.Path, fake_spice: FakeSpice,
                                      fake_spicedb: FakeSpiceDB) -> None:
    _ = load(tmp_path, pointing='fits')
    assert fake_spice.spice_frames == []
    assert 'CK' not in fake_spicedb.named('used_basenames')[0][1]['types']


def fits_frame_matrix(clock: float) -> Any:
    """The matrix of the frame defined by the synthetic header's pointing and `clock`."""

    (ra, dec) = NewHorizonsLORRI._radec_from_uv(512.5, 512.5, FITS_HEADER)
    frame = oops.frame.Cmatrix.from_ra_dec(ra, dec, clock)
    return frame.transform_at_time(0.).matrix.vals


@pytest.mark.parametrize(('utc', 'clock'), [
    ('2007-010T05:41:10.863', 275. + 90.),
    ('2012-365T23:59:59.000', 275. + 90.),
    ('2013-001T00:00:00.000', 275.),
])
def test_fits90_pointing(tmp_path: pathlib.Path, utc: str, clock: float) -> None:
    obs = load(tmp_path, pointing='fits90', fits_keywords={'SPCUTCID': utc})
    matrix = obs.frame.wrt(oops.Frame.J2000).transform_at_time(obs.midtime).matrix.vals
    assert matrix == pytest.approx(fits_frame_matrix(clock))


def test_fits_pointing_keeps_the_clock_angle(tmp_path: pathlib.Path) -> None:
    obs = load(tmp_path, pointing='fits')
    matrix = obs.frame.wrt(oops.Frame.J2000).transform_at_time(obs.midtime).matrix.vals
    assert matrix == pytest.approx(fits_frame_matrix(275.))


def test_fitsapp_pointing_without_motion(tmp_path: pathlib.Path) -> None:
    # The spacecraft path of the fake SPICE is at rest, so there is no aberration
    fits_obs = load(tmp_path, pointing='fits')
    app_obs = load(tmp_path, pointing='fitsapp')
    assert boresight(app_obs).vals == pytest.approx(boresight(fits_obs).vals, abs=1.e-12)


def test_fitsapp_pointing_with_motion(tmp_path: pathlib.Path,
                                      sun_path: oops.Path) -> None:
    # Moving at 30 km/s along x shifts the apparent direction by up to 1e-4 radians
    keywords = {'SPCSSCVX': -30., 'SPCSSCVY': 0., 'SPCSSCVZ': 0.}
    fits_obs = load(tmp_path, geom='fits', pointing='fits', fits_keywords=keywords)
    app_obs = load(tmp_path, geom='fits', pointing='fitsapp', fits_keywords=keywords)
    separation = boresight(app_obs).sep(boresight(fits_obs)).vals
    assert 1.e-5 < separation < 1.e-4

##########################################################################################
# from_file() on a PDS3 label
##########################################################################################


def ignores_case(directory: pathlib.Path) -> bool:
    """True if the file system of `directory` ignores the case of file names."""

    probe = directory / 'case_probe'
    probe.write_text('')
    return (directory / 'CASE_PROBE').exists()


# On a file system that ignores case, the first spelling tried is found whatever the case
# of the file, so the last two cases find the name as the label writes it
@pytest.mark.parametrize(('label_name', 'pointer', 'fits_name', 'found_ignoring_case'), [
    ('lor_test.lbl', 'LOR_TEST.FIT', 'lor_test.fit', 'lor_test.fit'),
    ('LOR_TEST.LBL', 'lor_test.fit', 'LOR_TEST.FIT', 'LOR_TEST.FIT'),
    ('Lor_Test.lbl', 'Lor_Test.fit', 'Lor_Test.fit', 'Lor_Test.fit'),
    ('Lor_Test.lbl', 'LOR_TEST.FIT', 'lor_test.fit', 'LOR_TEST.FIT'),
    ('Lor_Test.lbl', 'lor_test.fit', 'LOR_TEST.FIT', 'lor_test.fit'),
])
def test_from_file_label_finds_the_fits_file(tmp_path: pathlib.Path, label_name: str,
                                             pointer: str, fits_name: str,
                                             found_ignoring_case: str) -> None:
    expected = found_ignoring_case if ignores_case(tmp_path) else fits_name
    write_fits(tmp_path / fits_name)
    label = write_label(tmp_path / label_name, pointer)
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.basename == expected
    assert obs.filepath == FCPath(tmp_path / expected)


def test_from_file_label_gives_the_fits_snapshot(tmp_path: pathlib.Path) -> None:
    fits_path = write_fits(tmp_path / 'lor_test.fit')
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    from_label = NewHorizonsLORRI.from_file(label)
    from_fits = NewHorizonsLORRI.from_file(fits_path)
    assert from_label.time == from_fits.time
    assert np.array_equal(from_label.data, from_fits.data)
    assert from_label.image_url == from_fits.image_url
    assert from_label.abspath == from_fits.abspath


def test_from_file_parsed_label(tmp_path: pathlib.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit')
    label = Pds3Label(write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT'))
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.basename == 'lor_test.fit'


def test_from_file_label_with_a_missing_fits_file(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    with pytest.raises(FileNotFoundError,
                       match=r'FITS file LOR_TEST\.FIT of PDS3 label not found'):
        NewHorizonsLORRI.from_file(label)


def test_from_file_label_without_an_image(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', None)
    with pytest.raises(oops.host.PdsHostError,
                       match=r'PDS3 label lacks an \^IMAGE pointer to a FITS file: '
                             r'.*lor_test\.lbl'):
        NewHorizonsLORRI.from_file(label)


def test_from_file_label_of_another_instrument(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'mvic_test.lbl', 'MVIC_TEST.FIT',
                        INSTRUMENT_ID='"MVIC"')
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI PDS3'):
        NewHorizonsLORRI.from_file(label)


def test_from_file_label_astrometry_needs_no_fits_file(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    assert obs.filepath == FCPath(tmp_path / 'lor_test.fit')
    assert 'data' not in obs.subfields
    assert 'headers' not in obs.subfields


def test_from_file_label_astrometry_reads_no_fits_file(tmp_path: pathlib.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit', SPCCBTNM='EUROPA', EXPTIME=0.1)
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    assert obs.target == 'JUPITER'
    assert obs.texp == 5.


def test_from_file_label_astrometry_values(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=256)
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    tstart = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T05:41:08.363'))
    assert obs.time == pytest.approx((tstart, tstart + 5.), abs=1.e-6)
    assert obs.binning_mode == '4X4'
    assert obs.fov is NewHorizonsLORRI._FOVS['4X4', 'fast']
    assert obs.label_target == 'JUPITER'
    assert obs.spacecraft_clock_start_count == '0030713585:49000'
    assert obs.spacecraft_clock_stop_count == '0030713590:49000'
    assert obs.dict['INSTRUMENT_ID'] == 'LORRI'
    assert obs.frame.frame_id == 'NH_LORRI'


def test_from_file_label_astrometry_target(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT',
                        TARGET_NAME='"CALIBRATION"')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    assert obs.target == 'NONE'


def test_from_file_label_astrometry_rejects_an_unknown_shape(
        tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=512)
    with pytest.raises(oops.host.HostFormatError,
                       match='unrecognized LORRI image shape: 512 lines'):
        NewHorizonsLORRI.from_file(label, astrometry=True)


@pytest.mark.parametrize(('option', 'value'), [('geom', 'fits'), ('pointing', 'fits')])
def test_from_file_label_astrometry_rejects_fits_geometry(
        tmp_path: pathlib.Path, option: str, value: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    options: dict[str, Any] = {option: value}
    with pytest.raises(oops.OopsValueError, match='need the FITS file'):
        NewHorizonsLORRI.from_file(label, astrometry=True, **options)

##########################################################################################
# from_file() on other inputs
##########################################################################################


def test_from_file_hdulist(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    hdulist = fits.open(path)
    obs = NewHorizonsLORRI.from_file(hdulist)
    assert obs.filepath == FCPath(path)
    assert np.array_equal(obs.data, np.arange(1024 * 1024).reshape(1024, 1024))


def test_from_file_rejects_an_hdulist_of_another_instrument(
        tmp_path: pathlib.Path) -> None:
    hdulist = fits.open(write_fits(tmp_path / 'mvic_test.fit', INSTRU='mvi'))
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI FITS'):
        NewHorizonsLORRI.from_file(hdulist)
    assert hdulist._file.closed


def test_from_file_closes_an_hdulist(tmp_path: pathlib.Path) -> None:
    hdulist = fits.open(write_fits(tmp_path / 'lor_test.fit'))
    _ = NewHorizonsLORRI.from_file(hdulist)
    assert hdulist._file.closed


def test_from_file_closes_the_fits_file(tmp_path: pathlib.Path,
                                        monkeypatch: pytest.MonkeyPatch) -> None:
    closed = []
    original = fits.HDUList.close

    def close(self: fits.HDUList, *args: Any, **kwargs: Any) -> None:
        closed.append(self.filename())
        original(self, *args, **kwargs)

    path = write_fits(tmp_path / 'lor_test.fit')
    monkeypatch.setattr(fits.HDUList, 'close', close)
    _ = NewHorizonsLORRI.from_file(path)
    assert closed == [str(path)]


def test_from_file_fits_of_another_instrument(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'mvic_test.fit', INSTRU='mvi')
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI FITS'):
        NewHorizonsLORRI.from_file(path)


def test_from_file_invalid_fits(tmp_path: pathlib.Path) -> None:
    path = tmp_path / 'lor_test.fit'
    path.write_bytes(b'not a FITS file')
    with pytest.raises(oops.host.HostFormatError, match='not a valid FITS file'):
        NewHorizonsLORRI.from_file(path)


def test_from_file_initializes_lorri(tmp_path: pathlib.Path) -> None:
    _ = load(tmp_path)
    assert NewHorizonsLORRI._initialized


@pytest.mark.parametrize('suffix', ['fit', 'lbl'])
def test_host_from_file_dispatches_to_lorri(tmp_path: pathlib.Path, suffix: str) -> None:
    write_fits(tmp_path / 'lor_test.fit')
    write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = oops.Host.from_file(tmp_path / f'lor_test.{suffix}')
    assert obs.instrument == 'LORRI'
    assert obs.basename == 'lor_test.fit'

##########################################################################################
# from_index()
##########################################################################################


def test_from_index_rows(fake_spicedb: FakeSpiceDB) -> None:
    obs = list(NewHorizonsLORRI.from_index(INDEX))
    assert len(obs) == 30


def test_from_index_builds_each_row_when_reached(monkeypatch: pytest.MonkeyPatch) -> None:
    """The result is an iterator that builds each row's Snapshot when it reaches the row
    and yields the rows in index order."""

    built: list[Any] = []
    snapshot_from_row = NewHorizonsLORRI._snapshot_from_row

    def spy(row_dict: dict[str, Any], **kwargs: Any) -> Any:
        built.append(row_dict['FILE_SPECIFICATION_NAME'])
        return snapshot_from_row(row_dict, **kwargs)

    monkeypatch.setattr(NewHorizonsLORRI, '_snapshot_from_row', spy)
    snapshots = NewHorizonsLORRI.from_index(INDEX)
    assert isinstance(snapshots, Iterator)
    assert built == []
    first = next(snapshots)
    assert len(built) == 1
    rest = list(snapshots)
    assert len(rest) == 29
    assert first.basename == 'lor_0030702839_0x630_sci.fit'


def test_from_index_uses_the_supplement_values() -> None:
    obs = next(NewHorizonsLORRI.from_index(INDEX))
    assert obs.dict['SPACECRAFT_CLOCK_START_COUNT'] == '0030702838:48850'
    assert obs.dict['SPACECRAFT_CLOCK_STOP_COUNT'] == '0030702838:49000'


def test_from_index_merges_the_supplement() -> None:
    obs = next(NewHorizonsLORRI.from_index(INDEX))
    assert obs.dict['REDUCTION_LEVEL'] == 'LEVEL2'
    assert obs.dict['BINNING_MODE'] == '1x1'


def test_from_index_with_supplement() -> None:
    obs = list(NewHorizonsLORRI.from_index(INDEX, supplement=SUPPLEMENT))
    assert len(obs) == 30


def test_oops_host_from_index_passes_the_index_to_this_host() -> None:
    dispatched = oops.Host.from_index(INDEX, supplement=SUPPLEMENT)
    direct = NewHorizonsLORRI.from_index(INDEX, supplement=SUPPLEMENT)
    assert [obs.dict for obs in dispatched] == [obs.dict for obs in direct]


def copy_index(tmp_path: pathlib.Path, *, index: bool = True,
               supplement: bool = True) -> pathlib.Path:
    """Copy the edited index and supplement into `tmp_path`; return the index label."""

    names = []
    if index:
        names += ['NHJULO_2001_index_edited.lbl', 'NHJULO_2001_index_edited.tab']
    if supplement:
        names += ['NHJULO_2001_supplemental_index_edited.lbl',
                  'NHJULO_2001_supplemental_index_edited.tab']
    for name in names:
        (tmp_path / name).write_bytes((TEST_FILES / name).read_bytes())
    return tmp_path / 'NHJULO_2001_index_edited.lbl'


def test_from_index_requires_the_supplement_beside_it(tmp_path: pathlib.Path) -> None:
    index = copy_index(tmp_path, supplement=False)
    with pytest.raises(FileNotFoundError, match='supplemental_index_edited'):
        NewHorizonsLORRI.from_index(index)


def test_from_index_requires_a_named_supplement(tmp_path: pathlib.Path) -> None:
    index = copy_index(tmp_path)
    renamed = index.rename(tmp_path / 'NHJULO_2001_rows.lbl')
    with pytest.raises(oops.host.HostIndexError, match='supplemental index required'):
        NewHorizonsLORRI.from_index(renamed)


def test_from_index_requires_the_volume_column(tmp_path: pathlib.Path) -> None:
    """A supplement without VOLUME_ID is refused when from_index() is called, rather than
    paired by FILE_SPECIFICATION_NAME alone."""

    index = copy_index(tmp_path)
    supplement = tmp_path / 'NHJULO_2001_supplemental_index_edited.lbl'
    text = supplement.read_text()
    supplement.write_text(text.replace('= VOLUME_ID\n', '= VOLUME_XX\n'))
    with pytest.raises(oops.host.HostIndexError,
                       match=r'^the supplement lacks the column VOLUME_ID: '):
        NewHorizonsLORRI.from_index(index)


def test_from_index_rejects_another_host(tmp_path: pathlib.Path) -> None:
    index = copy_index(tmp_path)
    index.write_text(index.read_text().replace('INSTRUMENT_ID         = "LORRI"',
                                               'INSTRUMENT_ID         = "MVIC" '))
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI index'):
        NewHorizonsLORRI.from_index(index)


def rewrite_rows(tmp_path: pathlib.Path, edit: Any, *,
                 table: str = 'supplemental_index') -> pathlib.Path:
    """Copy the index and supplement into `tmp_path`, apply `edit` to the list of the
    rows of one table as bytes, and return the index label."""

    index = copy_index(tmp_path)
    path = tmp_path / f'NHJULO_2001_{table}_edited.tab'
    rows = path.read_bytes().splitlines(keepends=True)
    path.write_bytes(b''.join(edit(rows)))
    return index


def built_before(snapshots: Iterator[Any], error: type[Exception], match: str) -> int:
    """The number of Snapshots that an iteration yields before it raises `error`, whose
    message must match `match`."""

    built: list[Any] = []
    with pytest.raises(error, match=match):
        built.extend(snapshots)     # keeps each Snapshot given before the error
    return len(built)


def test_from_index_pairs_rows_by_file() -> None:
    assert NewHorizonsLORRI.MERGE_SUPPLEMENT_BY_FILE is True


def test_from_index_matches_supplement_rows_in_another_order(
        tmp_path: pathlib.Path) -> None:
    index = rewrite_rows(tmp_path, lambda rows: rows[::-1])
    in_order = NewHorizonsLORRI.from_index(INDEX)
    reordered = NewHorizonsLORRI.from_index(index)
    assert [obs.dict for obs in reordered] == [obs.dict for obs in in_order]


def test_from_index_reads_no_table_when_called(tmp_path: pathlib.Path) -> None:
    """The call reads only the two labels; the tables are first read when the iteration
    starts."""

    index = copy_index(tmp_path)
    (tmp_path / 'NHJULO_2001_index_edited.tab').unlink()
    snapshots = NewHorizonsLORRI.from_index(index)
    with pytest.raises(FileNotFoundError):
        next(snapshots)


def test_from_index_builds_rows_before_a_missing_row(tmp_path: pathlib.Path) -> None:
    """Every row that the index table holds gives its Snapshot before the iteration
    reaches the first row that the label lists but the table lacks (row 20 of 30)."""

    index = rewrite_rows(tmp_path, lambda rows: rows[:20], table='index')
    assert built_before(NewHorizonsLORRI.from_index(index), oops.host.HostIndexError,
                        'invalid index table') == 20


def test_from_index_rejects_a_mismatched_supplement(tmp_path: pathlib.Path) -> None:
    """An index row whose file the supplement does not list raises once the whole
    supplement has been read, after the rows before it."""

    index = rewrite_rows(tmp_path, lambda rows: [row.replace(b'lor_0030713591_0x633',
                                                             b'lor_0030713591_0x999')
                                                 for row in rows])
    match = (r"the supplement does not list the file of index row 21, "
             r"'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.lbl'")
    assert built_before(NewHorizonsLORRI.from_index(index), oops.host.HostIndexError,
                        match) == 21


def test_from_index_rejects_a_supplement_row_given_twice(tmp_path: pathlib.Path) -> None:
    """A supplement row read ahead that gives the same file as another row read ahead
    raises at once."""

    index = rewrite_rows(tmp_path, lambda rows: [rows[1], rows[1], rows[0]] + rows[3:])
    assert built_before(NewHorizonsLORRI.from_index(index), oops.host.HostIndexError,
                        r'a second time, at row 1') == 0


def test_from_index_rejects_an_index_row_given_twice(tmp_path: pathlib.Path) -> None:
    """An index row that repeats an earlier row's file is left without a partner, because
    the earlier row used it, and raises after the rows before it."""

    index = rewrite_rows(tmp_path, lambda rows: rows[:-1] + rows[:1], table='index')
    assert built_before(NewHorizonsLORRI.from_index(index), oops.host.HostIndexError,
                        'the supplement does not list the file of index row 29') == 29


def test_from_index_filter() -> None:
    obs = NewHorizonsLORRI.from_index(INDEX, filter={'BINNING_MODE': '4x4'})
    assert {o.binning_mode for o in obs} == {'4X4'}


def test_from_index_rejects_an_unknown_filter_column() -> None:
    """The call itself raises, before any iteration."""

    with pytest.raises(oops.OopsValueError, match='NOPE'):
        NewHorizonsLORRI.from_index(INDEX, filter={'NOPE': 1})


def test_from_index_rejects_an_invalid_fov_type() -> None:
    with pytest.raises(oops.OopsValueError, match="invalid fov_type 'wide'"):
        NewHorizonsLORRI.from_index(INDEX, fov_type='wide')


@pytest.mark.parametrize(('label_target', 'target'), [
    ('CALIBRATION', 'NONE'),
    ('PLUTO-CHARON', 'PLUTO'),
])
def test_from_index_row_target(label_target: str, target: str) -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(TARGET_NAME=label_target)])
    assert obs.target == target


# The supplement gives five moons of Jupiter without the "J<n> " of their labels; a
# corrected supplement, and the names of the other moons, are left as they are
@pytest.mark.parametrize(('row_target', 'label_target'), [
    ('IO', 'J1 IO'),
    ('CALLIRRHOE', 'J17 CALLIRRHOE'),
    ('J1 IO', 'J1 IO'),
    ('GANYMEDE', 'GANYMEDE'),
])
def test_from_index_row_label_target(row_target: str, label_target: str) -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(TARGET_NAME=row_target)])
    assert obs.label_target == label_target


def test_from_index_row_dict_gives_the_label_target_name() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(TARGET_NAME='EUROPA')])
    assert obs.dict['TARGET_NAME'] == 'J2 EUROPA'


def test_from_index_row_snapshot() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()])
    tstart = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T05:41:08.363'))
    assert obs.time == pytest.approx((tstart, tstart + 5.))
    assert obs.fov is NewHorizonsLORRI._FOVS['4X4', 'fast']
    assert obs.target == 'CALLIRRHOE'
    assert obs.frame.frame_id == 'NH_LORRI'
    assert obs.path.path_id == 'NEW HORIZONS'


def test_from_index_row_subfields() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()])
    assert obs.binning_mode == '4X4'
    assert obs.label_target == 'J17 CALLIRRHOE'
    assert obs.spacecraft_clock_start_count == '0030713585:49000'
    assert obs.spacecraft_clock_stop_count == '0030713590:49000'
    assert obs.spice_frame_id == -98300


@pytest.mark.parametrize(('spec', 'filepath'), [
    ('data/20070110_003071/lor_0030713591_0x633_sci.lbl',
     'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit'),
    ('DATA/20070110_003071/LOR_0030713591_0X633_SCI.LBL',
     'NHJULO_2001/DATA/20070110_003071/LOR_0030713591_0X633_SCI.FIT'),
])
def test_from_index_filepath_names_the_fits_file(spec: str, filepath: str) -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX,
                                         row_dicts=[row(FILE_SPECIFICATION_NAME=spec)])
    assert str(obs.filepath) == (TEST_FILES / filepath).as_posix()
    assert obs.basename == filepath.rpartition('/')[2]


def test_from_index_names_the_fits_file_in_full() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()])
    image = TEST_FILES / 'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit'
    assert obs.image_url == image.as_posix()


# The FITS file of row() and an index in a metadata tree, relative to the holdings
ROW_IMAGE = ('volumes/NHxxLO_xxxx/NHJULO_2001/data/20070110_003071/'
             'lor_0030713591_0x633_sci.fit')
ROW_INDEX = 'metadata/NHxxLO_xxxx/NHJULO_2001/NHJULO_2001_index.lbl'


def test_from_index_local_image_has_abspath(tmp_path: pathlib.Path) -> None:
    image = tmp_path / ROW_IMAGE
    image.parent.mkdir(parents=True)
    image.write_bytes(b'')
    (obs,) = NewHorizonsLORRI.from_index(tmp_path / ROW_INDEX, row_dicts=[row()])
    assert obs.abspath == image.resolve()


def test_from_index_missing_local_image_has_no_abspath(tmp_path: pathlib.Path) -> None:
    (obs,) = NewHorizonsLORRI.from_index(tmp_path / ROW_INDEX, row_dicts=[row()])
    assert 'abspath' not in obs.subfields


def test_from_index_creates_no_directories(tmp_path: pathlib.Path) -> None:
    holdings = tmp_path / 'holdings'
    (holdings / ROW_INDEX).parent.mkdir(parents=True)
    list(NewHorizonsLORRI.from_index(holdings / ROW_INDEX, row_dicts=[row()]))
    assert not (holdings / 'volumes').exists()


# An index read leaves the FITS file unfetched, so only an earlier read can have cached it
def test_from_index_remote_image_has_no_abspath() -> None:
    (obs,) = NewHorizonsLORRI.from_index('https://h.org/' + ROW_INDEX, row_dicts=[row()])
    assert 'abspath' not in obs.subfields


def test_from_index_remote_image_in_the_cache_has_abspath() -> None:
    cached = cache_file('https://h.org/' + ROW_IMAGE)
    (obs,) = NewHorizonsLORRI.from_index('https://h.org/' + ROW_INDEX, row_dicts=[row()])
    assert obs.abspath == cached.resolve()



# A volume's own index lives in its INDEX directory, in either case
@pytest.mark.parametrize('index', ['/h/volumes/NHxxLO_xxxx/NHK2LO_1001/index/index.lbl',
                                   '/h/NHK2LO_1001/INDEX/INDEX.LBL'])
def test_from_index_refuses_a_volume_index(index: str) -> None:
    with pytest.raises(oops.host.HostIndexError,
                       match=r'New Horizons LORRI reads the Ring-Moon Systems Node '
                             r'metadata index, such as "NHJULO_2001_index\.lbl", not a '
                             r"volume's own index: .*(index|INDEX)\.(lbl|LBL)$"):
        NewHorizonsLORRI.from_index(index, row_dicts=[row()])


# The supplement is required, so None is refused rather than taken as "search for it"
def test_from_index_refuses_no_supplement() -> None:
    with pytest.raises(oops.host.HostIndexError,
                       match=r'New Horizons LORRI index needs its supplemental index: '
                             r'.*NHJULO_2001_index_edited\.lbl'):
        NewHorizonsLORRI.from_index(INDEX, supplement=None)


def test_from_index_rows_need_no_supplement() -> None:
    assert len(list(NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()],
                                                supplement=None))) == 1


def test_from_index_uses_a_given_label(monkeypatch: pytest.MonkeyPatch) -> None:
    label = Pds3Label(INDEX)
    labels: list[Any] = []
    index_rows = NewHorizonsLORRI._index_rows

    def spy(filepath: Any, **kwargs: Any) -> Any:
        labels.append(kwargs.get('label'))
        return index_rows(filepath, **kwargs)

    monkeypatch.setattr(NewHorizonsLORRI, '_index_rows', spy)
    NewHorizonsLORRI.from_index(INDEX, label=label)
    assert labels == [label]


def test_from_index_rejects_an_unknown_keyword() -> None:
    with pytest.raises(oops.OopsTypeError,
                       match=r"unexpected keyword argument for "
                             r"NewHorizonsLORRI\.from_index: 'supplementary'"):
        NewHorizonsLORRI.from_index(INDEX, supplementary=SUPPLEMENT)


@pytest.mark.parametrize('select', [0, slice(None), 'data'])
def test_from_index_refuses_select(select: Any) -> None:
    """The call itself raises, before any iteration."""

    with pytest.raises(oops.OopsValueError,
                       match='each row of a New Horizons LORRI index is one observation'):
        NewHorizonsLORRI.from_index(INDEX, select=select)


def test_from_index_rejects_an_unknown_binning_mode() -> None:
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(BINNING_MODE='2x2')])
    with pytest.raises(oops.host.HostFormatError,
                       match="unrecognized LORRI BINNING_MODE: '2x2'"):
        next(snapshots)


@pytest.mark.parametrize('null', ['N/A', None])
def test_from_index_rejects_a_missing_binning_mode(null: str | None) -> None:
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(BINNING_MODE=null)])
    with pytest.raises(oops.host.HostFormatError,
                       match=r'unrecognized LORRI BINNING_MODE: missing in '
                             r'.*NHJULO_2001/'):
        next(snapshots)


def test_from_file_fits_pointing_needs_no_header_epoch(tmp_path: pathlib.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit', SPCSCET=None)
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label, pointing='fits')
    assert obs.frame.frame_id == '.NH_FRAME_lor_test.fit'


def test_from_index_custom_frame_reports_no_cks(fake_spice: FakeSpice,
                                                fake_spicedb: FakeSpiceDB) -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()], frame='J2000')
    assert obs.frame.frame_id == 'J2000'
    assert fake_spice.spice_frames == []
    assert 'CK' not in fake_spicedb.named('used_basenames')[0][1]['types']


@pytest.mark.parametrize('null', ['N/A', 'UNK', None])
def test_from_index_row_null_values(null: str | None) -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(
        TARGET_NAME=null, SPACECRAFT_CLOCK_START_COUNT=null)])
    assert obs.target == 'NONE'
    assert 'label_target' not in obs.subfields
    assert 'spacecraft_clock_start_count' not in obs.subfields
    assert 'TARGET_NAME' not in obs.dict


def masked(**columns: float) -> dict[str, Any]:
    """Columns of a supplement row, each masked by the table reader, as its value equals
    the null constant that the column declares."""

    result: dict[str, Any] = {}
    for (key, value) in columns.items():
        result[key] = np.float64(value)
        result[key + '_mask'] = np.True_
    return result


def vector(name: str, *values: float, mask: bool = True) -> dict[str, Any]:
    """The components of a vector of a supplement row, masked where they are zero."""

    result: dict[str, Any] = {}
    suffixes = 'XYZ' if len(values) == 3 else '1234'
    for (suffix, value) in zip(suffixes, values, strict=True):
        key = f'{name}_{suffix}'
        result[key] = np.float64(value)
        result[key + '_mask'] = np.bool_(mask and value == 0.)
    return result


def row_dict_of(**columns: Any) -> dict[str, Any]:
    """The `dict` of the Snapshot of a row with the given columns."""

    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(**columns)])
    row_dict: dict[str, Any] = obs.dict
    return row_dict


def test_from_index_row_mask_keys_are_left_out() -> None:
    row_dict = row_dict_of(SOLAR_ELONGATION=np.float64(9.4),
                           SOLAR_ELONGATION_mask=np.False_)
    assert [key for key in row_dict if key.endswith('_mask')] == []


def test_from_index_row_value_of_a_nonzero_null_constant_is_left_out() -> None:
    assert 'PHASE_ANGLE' not in row_dict_of(**masked(PHASE_ANGLE=-999.))


def test_from_index_row_zero_null_constant_of_a_scalar_is_left_out() -> None:
    row_dict = row_dict_of(**masked(SC_GEOCENTRIC_DISTANCE=0.))
    assert 'SC_GEOCENTRIC_DISTANCE' not in row_dict


# A zero component of a vector that is not zero is a real value, rounded to zero
def test_from_index_row_zero_component_of_a_vector_is_kept() -> None:
    row_dict = row_dict_of(**vector('QUATERNIAN', 0.5, 0., 0.5, 0.7))
    assert row_dict['QUATERNIAN_2'] == 0.


# The supplement writes zeros for the spacecraft's velocity relative to the Earth, which
# is missing from every label; the column of the Y component does not declare a null
# constant
def test_from_index_row_zero_vector_is_left_out() -> None:
    columns = vector('SC_EARTH_VELOCITY_VECTOR', 0., 0., 0.)
    columns['SC_EARTH_VELOCITY_VECTOR_Y_mask'] = np.False_
    row_dict = row_dict_of(**columns)
    assert [key for key in row_dict if key.startswith('SC_EARTH_VELOCITY')] == []


# A corrected supplement declares a fill value that is never real, so a zero is not
# masked, and a zero vector is a real value
def test_from_index_row_zero_vector_that_is_not_masked_is_kept() -> None:
    row_dict = row_dict_of(**vector('SC_EARTH_VELOCITY_VECTOR', 0., 0., 0., mask=False))
    assert row_dict['SC_EARTH_VELOCITY_VECTOR_Y'] == 0.


# A corrected supplement spells the quaternion as the label does
def test_from_index_row_zero_component_with_the_label_spelling_is_kept() -> None:
    row_dict = row_dict_of(**vector('QUATERNION', 0.5, 0., 0.5, 0.7))
    assert row_dict['QUATERNION_2'] == 0.


# An image without target geometry has a zero vector from the spacecraft to the target
def test_from_index_row_target_geometry_without_a_target_is_left_out() -> None:
    row_dict = row_dict_of(**vector('SC_TARGET_POSITION_VECTOR', 0., 0., 0.),
                           **vector('TARGET_SUN_POSITION_VECTOR', 0., 0., 0.),
                           **masked(SOLAR_DISTANCE=0.))
    assert [key for key in row_dict
            if key.startswith(('SC_TARGET_', 'TARGET_SUN_', 'SOLAR_'))] == []


# The Sun is zero kilometers from itself
def test_from_index_row_zero_target_geometry_of_a_target_is_kept() -> None:
    row_dict = row_dict_of(**vector('SC_TARGET_POSITION_VECTOR', 2.e8, 1.e9, 4.e8),
                           **vector('TARGET_SUN_POSITION_VECTOR', 0., 0., 0.),
                           **masked(SOLAR_DISTANCE=0.))
    assert [row_dict[key] for key in ('TARGET_SUN_POSITION_VECTOR_X', 'SOLAR_DISTANCE')
            ] == [0., 0.]


def test_from_index_dict_is_json() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(
        APPROX_TARGET_LINE=np.int64(527), SEQUENCE_ID=np.str_('JELR_JOBSATM01'),
        SOLAR_ELONGATION_mask=np.False_)])
    assert json.loads(json.dumps(obs.dict)) == obs.dict


def test_from_index_dict_of_the_test_index_is_json() -> None:
    obs = next(NewHorizonsLORRI.from_index(INDEX))
    assert json.loads(json.dumps(obs.dict)) == obs.dict


@pytest.mark.parametrize(('missing', 'column'), [
    (('START_TIME', 'SPACECRAFT_CLOCK_START_COUNT'), 'START_TIME'),
    (('EXPOSURE_DURATION',), 'EXPOSURE_DURATION'),
])
def test_from_index_row_without_a_needed_column(missing: tuple[str, ...],
                                                column: str) -> None:
    """The row raises when the iteration reaches it, after the row before it."""

    row_dict = row()
    for key in missing:
        del row_dict[key]
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row(), row_dict])
    next(snapshots)
    with pytest.raises(oops.host.HostFormatError,
                       match=f'no {column} for LORRI image .*NHJULO_2001/data/'):
        next(snapshots)


def fake_scs2e(monkeypatch: pytest.MonkeyPatch) -> list[tuple[int, str]]:
    """Replace cspyce.scs2e with a fake that gives 1000 seconds plus the count's seconds;
    return the list of its calls."""

    calls: list[tuple[int, str]] = []

    def scs2e(sc: int, count: str) -> float:
        calls.append((sc, count))
        return 1000. + float(count.partition(':')[0])

    monkeypatch.setattr(cspyce, 'scs2e', scs2e)
    return calls


# Without START_TIME, the start comes from the clock count through the clock kernel
def test_from_index_row_time_from_the_clock_count(
        monkeypatch: pytest.MonkeyPatch) -> None:
    calls = fake_scs2e(monkeypatch)
    row_dict = row()
    del row_dict['START_TIME']
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row_dict])
    assert calls == [(-98, '0030713585:49000')]
    assert obs.time == (1000. + 30713585., 1000. + 30713590.)
    assert obs.time_from_sclk is True


def test_from_file_time_from_the_clock_count(tmp_path: pathlib.Path,
                                             monkeypatch: pytest.MonkeyPatch) -> None:
    fake_scs2e(monkeypatch)
    obs = load(tmp_path, fits_keywords={'SPCSCET': None})
    assert obs.time[0] == 1000. + 30713585.
    assert obs.time_from_sclk is True


def test_from_file_time_from_the_header(tmp_path: pathlib.Path) -> None:
    assert load(tmp_path).time_from_sclk is False


def test_from_index_row_time_from_its_start_time() -> None:
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row()])
    assert obs.time_from_sclk is False


def test_time_from_an_unreadable_clock_count(monkeypatch: pytest.MonkeyPatch) -> None:
    def scs2e(sc: int, count: str) -> float:
        raise ValueError('SPICE(INVALIDSCLKSTRING) -- Invalid SCLK string')

    monkeypatch.setattr(cspyce, 'scs2e', scs2e)
    row_dict = row()
    del row_dict['START_TIME']
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row_dict])
    with pytest.raises(oops.OopsSpiceError,
                       match=r'no time for the LORRI clock count 0030713585:49000 for '
                             r'.*lor_0030713591_0x633_sci\.fit: '
                             r'SPICE\(INVALIDSCLKSTRING\)'):
        next(snapshots)


def test_clock_error_without_a_spice_code_passes_through(
        monkeypatch: pytest.MonkeyPatch) -> None:
    def scs2e(sc: int, count: str) -> float:
        raise ValueError('something else')

    monkeypatch.setattr(cspyce, 'scs2e', scs2e)
    row_dict = row()
    del row_dict['START_TIME']
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row_dict])
    with pytest.raises(ValueError, match=r'^something else$'):
        next(snapshots)


# Without a clock count, which would give the start time instead
@pytest.mark.parametrize(('column', 'value'), [
    ('START_TIME', '2007-01-10T25:41:08'),
    ('EXPOSURE_DURATION', 'five seconds'),
])
def test_from_index_row_with_an_invalid_value(column: str, value: str) -> None:
    row_dict = row(**{column: value})
    del row_dict['SPACECRAFT_CLOCK_START_COUNT']
    snapshots = NewHorizonsLORRI.from_index(INDEX, row_dicts=[row_dict])
    with pytest.raises(oops.host.HostFormatError, match=f'no {column} for LORRI image'):
        next(snapshots)

##########################################################################################
