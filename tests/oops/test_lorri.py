##########################################################################################
# tests/oops/test_lorri.py
##########################################################################################
"""Tests of NewHorizonsLORRI.from_file(): what it reads, the errors it raises, and the
observation it builds."""

import os
import pathlib
import shutil
import subprocess
import sys
from typing import Any

import cspyce
import numpy as np
import oops
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from conftest import (FITS_HEADER, TEST_FILES, FakeSpice, FakeSpiceFrame, write_fits,
                      write_label, write_product)
from filecache import FCPath

from host_new_horizons.oops.lorri import NewHorizonsLORRI

SCI = 'lor_0030713591_0x633_sci_edited'
ENG = 'lor_0030713591_0x633_eng_edited'


def single(obs: Any) -> oops.observation.Snapshot:
    """The observation, checked to be a single Snapshot."""

    assert isinstance(obs, oops.observation.Snapshot)
    return obs

##########################################################################################
# Options and initialization
##########################################################################################


@pytest.mark.parametrize('option', ['geom', 'pointing', 'fov_type'])
def test_invalid_option(tmp_path: pathlib.Path, option: str) -> None:
    """An option value that is not listed raises OopsValueError naming it."""

    with pytest.raises(oops.OopsValueError, match=f"invalid {option} 'bad'"):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), **{option: 'bad'})


def test_calibration_keyword_is_refused(tmp_path: pathlib.Path) -> None:
    """The calibration keyword is an unknown keyword."""

    with pytest.raises(TypeError, match="unexpected keyword argument 'calibration'"):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                   calibration=False)


def test_initialize_builds_the_fovs() -> None:
    """The FOVs cover both binning modes and the three FOV types."""

    NewHorizonsLORRI.initialize()
    assert sorted(NewHorizonsLORRI._FOVS) == [('1X1', 'fast'), ('1X1', 'flat'),
                                              ('1X1', 'slow'), ('4X4', 'fast'),
                                              ('4X4', 'flat'), ('4X4', 'slow')]


@pytest.mark.parametrize(('fov_type', 'fov_class'), [
    ('fast', oops.fov.PolynomialFOV),
    ('slow', oops.fov.PolynomialFOV),
    ('flat', oops.fov.FlatFOV),
])
def test_full_resolution_fov(tmp_path: pathlib.Path, fov_type: str,
                             fov_class: type) -> None:
    """A 1024-line image has the full-resolution FOV of the requested type."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            fov_type=fov_type))
    assert isinstance(obs.fov, fov_class)


def test_slow_fov_has_no_inverse(tmp_path: pathlib.Path) -> None:
    """The "slow" FOV inverts the polynomial by Newton's method."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            fov_type='slow'))
    assert obs.fov.coefft_xy_from_uv is None


def test_binned_fov(tmp_path: pathlib.Path) -> None:
    """A 256-line image has the 4x4 binned FOV."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', lines=256,
                                                       samples=256)))
    assert tuple(obs.fov.uv_shape.vals) == (256, 256)


def test_initialize_registers_the_camera_frame(fake_spice: FakeSpice) -> None:
    """The SPICE frame NH_LORRI is registered as NH_LORRI_FLIPPED."""

    NewHorizonsLORRI.initialize()
    assert fake_spice.spice_frames == [('NH_LORRI', 'NH_LORRI_FLIPPED')]


def test_camera_frame_flips_y_and_z(tmp_path: pathlib.Path) -> None:
    """The camera frame NH_LORRI is the SPICE frame with its y and z axes flipped."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    matrix = oops.Frame.as_frame(obs.frame).transform_at_time(0.).matrix.vals
    assert np.array_equal(matrix, np.diag([1., -1., -1.]))


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    """A second call does nothing."""

    NewHorizonsLORRI.initialize()
    NewHorizonsLORRI.initialize(asof='2020-01-01')
    assert len(fake_spice.spice_frames) == 1


def test_import_loads_no_kernels() -> None:
    """Importing the host module opens no SPICE database."""

    src = pathlib.Path(__file__).parents[2] / 'src'
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join([str(src), env.get('PYTHONPATH', '')])
    code = ('import spicedb\n'
            'def refuse(*args, **kwargs):\n'
            '    raise RuntimeError("kernels loaded at import")\n'
            'spicedb.open_db = refuse\n'
            'import host_new_horizons.oops\n'
            'print("imported")\n')
    result = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True,
                            text=True, check=True)
    assert result.stdout.strip() == 'imported'

##########################################################################################
# What is read
##########################################################################################


def test_fits_file_alone(tmp_path: pathlib.Path) -> None:
    """A FITS file is read alone, with its arrays, headers and no label."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert sorted(set(obs.subfields) & {'data', 'error', 'quality', 'headers',
                                        'header_dict', 'label_dict', 'label_url',
                                        'label_abspath'}) == [
        'data', 'error', 'header_dict', 'headers', 'quality']


def test_fits_data(tmp_path: pathlib.Path) -> None:
    """`data` is the primary image."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.data[1, 2] == 1026.


def test_error_and_quality_by_extname(tmp_path: pathlib.Path) -> None:
    """`error` and `quality` come from the extensions with those names."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.error[0, 0] == 0.5
    assert obs.quality[0, 0] == 1


def test_raw_product_has_no_error_or_quality(tmp_path: pathlib.Path) -> None:
    """A raw product's extensions are neither `error` nor `quality`."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit',
                                                       calibrated=False)))
    assert {'error', 'quality'} & set(obs.subfields) == set()


def test_headers_of_every_hdu(tmp_path: pathlib.Path) -> None:
    """`headers` holds the astropy Header of every HDU."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert [header.get('EXTNAME') for header in obs.headers] == [
        None, 'LORRI Error image', 'LORRI Quality flag image']


def test_fits_astrometry_reads_no_arrays(tmp_path: pathlib.Path) -> None:
    """An astrometry read of a FITS file reads its headers and no arrays."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            astrometry=True))
    assert sorted(set(obs.subfields) & {'data', 'error', 'quality', 'headers'}) == [
        'headers']


def test_label_reads_the_fits_file(tmp_path: pathlib.Path) -> None:
    """A label is read with the FITS file that it names, in the label's case."""

    obs = single(NewHorizonsLORRI.from_file(write_product(tmp_path)))
    assert obs.filepath == FCPath(tmp_path / 'lor_0030713591_0x633_sci.fit')


def test_label_gives_label_dict_and_arrays(tmp_path: pathlib.Path) -> None:
    """A label read gives the label and the FITS file's content."""

    obs = single(NewHorizonsLORRI.from_file(write_product(tmp_path)))
    assert {'label_dict', 'header_dict', 'headers', 'data'} <= set(obs.subfields)


def test_upper_case_label_names_an_upper_case_fits_file(tmp_path: pathlib.Path) -> None:
    """A label named in upper case finds its FITS file in upper case."""

    write_fits(tmp_path / 'LOR_1.FIT')
    obs = single(NewHorizonsLORRI.from_file(write_label(tmp_path / 'LOR_1.LBL',
                                                        'lor_1.fit')))
    assert obs.basename == 'LOR_1.FIT'


def test_mixed_case_label_keeps_the_pointer_case(tmp_path: pathlib.Path) -> None:
    """A label named in mixed case finds its FITS file as the pointer writes it."""

    write_fits(tmp_path / 'Lor_1.Fit')
    obs = single(NewHorizonsLORRI.from_file(write_label(tmp_path / 'Lor_1.lbl',
                                                        'Lor_1.Fit')))
    assert obs.basename == 'Lor_1.Fit'


def test_label_astrometry_needs_no_fits_file(tmp_path: pathlib.Path) -> None:
    """An astrometry read of a label reads the label alone."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT')
    obs = single(NewHorizonsLORRI.from_file(label, astrometry=True))
    assert sorted(set(obs.subfields) & {'data', 'headers', 'header_dict', 'label_dict',
                                        'image_abspath'}) == ['label_dict']


def test_label_astrometry_names_the_fits_file(tmp_path: pathlib.Path) -> None:
    """An astrometry read of a label still names the FITS file."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT')
    obs = single(NewHorizonsLORRI.from_file(label, astrometry=True))
    assert obs.image_url == (tmp_path / 'lor_1.fit').as_posix()


@pytest.mark.parametrize('option', ['geom', 'pointing'])
def test_label_astrometry_refuses_header_geometry(tmp_path: pathlib.Path,
                                                  option: str) -> None:
    """Geometry from the FITS header needs the FITS file."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT')
    with pytest.raises(oops.OopsValueError, match='need the FITS file'):
        NewHorizonsLORRI.from_file(label, astrometry=True, **{option: 'fits'})


def test_missing_fits_file_of_a_label(tmp_path: pathlib.Path) -> None:
    """A label whose FITS file does not exist raises FileNotFoundError."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT')
    with pytest.raises(FileNotFoundError, match=r'lor_1\.fit'):
        NewHorizonsLORRI.from_file(label)


def test_hdulist_is_read_and_closed(tmp_path: pathlib.Path) -> None:
    """An HDUList is read as its file and closed afterward."""

    path = write_fits(tmp_path / 'a.fit')
    hdulist = fits.open(path)
    obs = single(NewHorizonsLORRI.from_file(hdulist))
    assert obs.filepath == FCPath(path)
    assert hdulist._file.closed


def test_dispatcher_reads_a_label(tmp_path: pathlib.Path) -> None:
    """oops.Host.from_file hands a LORRI label to this host."""

    obs = single(oops.Host.from_file(write_product(tmp_path)))
    assert obs.label_url == (tmp_path / 'lor_0030713591_0x633_sci.lbl').as_posix()


def test_dispatcher_reads_a_fits_file(tmp_path: pathlib.Path) -> None:
    """oops.Host.from_file hands a LORRI FITS file to this host."""

    path = write_fits(tmp_path / 'a.fit')
    obs = single(oops.Host.from_file(path))
    assert obs.filepath == FCPath(path)


def test_default_parse_method_is_strict(tmp_path: pathlib.Path,
                                        monkeypatch: pytest.MonkeyPatch) -> None:
    """A label is parsed with the strict method unless the caller gives another."""

    options = []
    original = NewHorizonsLORRI._read_fileinfo

    def spy(fileinfo: Any, **kwargs: Any) -> Any:
        options.append(kwargs.get('pds3_options'))
        return original(fileinfo, **kwargs)

    monkeypatch.setattr(NewHorizonsLORRI, '_read_fileinfo', spy)
    NewHorizonsLORRI.from_file(write_product(tmp_path))
    assert options[0] == {'method': 'strict'}


def test_caller_parse_method(tmp_path: pathlib.Path,
                             monkeypatch: pytest.MonkeyPatch) -> None:
    """The caller's parser options replace the default method."""

    options = []
    original = NewHorizonsLORRI._read_fileinfo

    def spy(fileinfo: Any, **kwargs: Any) -> Any:
        options.append(kwargs.get('pds3_options'))
        return original(fileinfo, **kwargs)

    monkeypatch.setattr(NewHorizonsLORRI, '_read_fileinfo', spy)
    NewHorizonsLORRI.from_file(write_product(tmp_path), pds3_options={'method': 'fast'})
    assert options[0] == {'method': 'fast'}


@pytest.mark.parametrize(('name', 'shape'), [(SCI, (256, 256)), (ENG, (256, 257))])
def test_archive_products(tmp_path: pathlib.Path, name: str,
                          shape: tuple[int, int]) -> None:
    """The cut-down archive products, whose labels name their FITS files in upper case,
    are read with their arrays."""

    obs = single(NewHorizonsLORRI.from_file(TEST_FILES / f'{name}.lbl'))
    assert obs.data.shape == shape


def test_archive_calibrated_product_extensions() -> None:
    """The archive's calibrated product gives its error and quality arrays."""

    obs = single(NewHorizonsLORRI.from_file(TEST_FILES / f'{SCI}.fit'))
    assert obs.error.shape == obs.quality.shape

##########################################################################################
# Files that cannot be read
##########################################################################################


def test_truncated_fits_file(tmp_path: pathlib.Path) -> None:
    """A FITS file that ends early raises HostFormatError naming the file."""

    path = write_fits(tmp_path / 'a.fit')
    data = path.read_bytes()
    path.write_bytes(data[:len(data) // 2])
    with pytest.raises(oops.host.HostFormatError, match=r'cannot read .*a\.fit'):
        NewHorizonsLORRI.from_file(path)


def test_truncated_fits_file_in_astrometry(tmp_path: pathlib.Path) -> None:
    """A FITS file that ends early raises HostFormatError even when no arrays are read."""

    path = write_fits(tmp_path / 'a.fit')
    data = path.read_bytes()
    path.write_bytes(data[:len(data) - 4000])
    with pytest.raises(oops.host.HostFormatError, match='truncated'):
        NewHorizonsLORRI.from_file(path, astrometry=True)


def test_multiband_image(tmp_path: pathlib.Path) -> None:
    """A primary array of more than one band raises HostFormatError."""

    path = write_fits(tmp_path / 'a.fit', data=np.zeros((2, 1024, 1024), dtype='>f4'),
                      calibrated=False)
    with pytest.raises(oops.host.HostFormatError, match='not one image band'):
        NewHorizonsLORRI.from_file(path)


def test_array_that_does_not_match_its_label(tmp_path: pathlib.Path) -> None:
    """An array whose shape differs from its label's raises HostFormatError."""

    write_fits(tmp_path / 'lor_1.fit')
    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT', samples=1028)
    with pytest.raises(oops.host.HostFormatError,
                       match=r'PDS3 label gives \(1024, 1028\)'):
        NewHorizonsLORRI.from_file(label)


def test_missing_exposure(tmp_path: pathlib.Path) -> None:
    """A file without an exposure duration raises HostFormatError."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT', EXPOSURE_DURATION=None)
    with pytest.raises(oops.host.HostFormatError, match="KeyError: 'EXPOSURE_DURATION'"):
        NewHorizonsLORRI.from_file(label, astrometry=True)


def test_missing_start_time(tmp_path: pathlib.Path) -> None:
    """A FITS file without SPCSCET, and so without a start time, raises
    HostFormatError."""

    with pytest.raises(oops.host.HostFormatError, match="KeyError: 'START_TIME'"):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', SPCSCET=None))


def test_unsupported_image_size(tmp_path: pathlib.Path) -> None:
    """An image of a size that LORRI does not take raises HostFormatError."""

    with pytest.raises(oops.host.HostFormatError, match='KeyError: 512'):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', lines=512,
                                              samples=512))


def test_label_without_image_object(tmp_path: pathlib.Path) -> None:
    """A label without the IMAGE object that describes the array raises
    HostFormatError."""

    write_fits(tmp_path / 'lor_1.fit')
    label = tmp_path / 'lor_1.lbl'
    text = write_label(label, 'LOR_1.FIT').read_text()
    label.write_text(text.split('OBJECT = IMAGE')[0] + 'END\r\n')
    with pytest.raises(oops.host.HostFormatError, match="KeyError: 'IMAGE'"):
        NewHorizonsLORRI.from_file(label)


def test_invalid_fits_file(tmp_path: pathlib.Path) -> None:
    """A ".fit" file that is not FITS raises HostFormatError."""

    path = tmp_path / 'a.fit'
    path.write_text('not FITS')
    with pytest.raises(oops.host.HostFormatError, match='not a valid FITS file'):
        NewHorizonsLORRI.from_file(path)


def test_unparsable_label(tmp_path: pathlib.Path) -> None:
    """A label that the parser refuses raises PdsHostError."""

    label = tmp_path / 'lor_1.lbl'
    label.write_text('PDS_VERSION_ID = PDS3\r\nTHIS IS NOT PDS3\r\nEND\r\n')
    with pytest.raises(oops.host.PdsHostError, match='cannot parse the PDS3 label'):
        NewHorizonsLORRI.from_file(label)


# oops leaves open the FITS file that it refuses for this host, which only the garbage
# collector closes
@pytest.mark.filterwarnings('ignore::ResourceWarning',
                            'ignore::pytest.PytestUnraisableExceptionWarning')
def test_fits_file_of_another_instrument(tmp_path: pathlib.Path) -> None:
    """A FITS file of another instrument raises WrongHostError."""

    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI'):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', INSTRU='mvc'))


def test_missing_file(tmp_path: pathlib.Path) -> None:
    """A file that does not exist raises FileNotFoundError."""

    with pytest.raises(FileNotFoundError):
        NewHorizonsLORRI.from_file(tmp_path / 'a.fit')

##########################################################################################
# The observation
##########################################################################################


def test_times_from_the_header(tmp_path: pathlib.Path) -> None:
    """The exposure starts EXPTIME/2 before SPCSCET and lasts EXPTIME."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.time == (FITS_HEADER['SPCSCET'] - 2.5, FITS_HEADER['SPCSCET'] + 2.5)


def test_zero_exposure_loads(tmp_path: pathlib.Path) -> None:
    """An exposure of zero gives an observation with an exposure of zero."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', EXPTIME=0.)))
    assert obs.texp == 0.


def test_lightsource_is_the_sun(tmp_path: pathlib.Path) -> None:
    """The light source is the Sun."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.lightsource == 'SUN'


def test_path_from_the_sp_kernels(tmp_path: pathlib.Path) -> None:
    """With geom="spice", the observer is the registered spacecraft path."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.path == oops.Path.as_waypoint('NEW HORIZONS')


def test_path_from_the_header(tmp_path: pathlib.Path, sun_path: oops.Path) -> None:
    """With geom="fits", the observer moves linearly about the Sun through the header's
    position at SPCSCET."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), geom='fits'))
    event = oops.Path.as_path(obs.path).event_at_time(FITS_HEADER['SPCSCET'])
    assert np.array_equal(event.pos.vals, [-1.e8, -2.e8, -3.e8])


def test_path_from_the_header_moves(tmp_path: pathlib.Path, sun_path: oops.Path) -> None:
    """With geom="fits", the observer moves with the header's velocity."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), geom='fits'))
    event = oops.Path.as_path(obs.path).event_at_time(FITS_HEADER['SPCSCET'] + 10.)
    assert np.allclose(event.pos.vals, [-1.e8 - 10., -2.e8 - 20., -3.e8 - 30.])


def boresight(obs: oops.observation.Snapshot) -> tuple[float, float]:
    """The right ascension and declination of the camera's optic axis in degrees."""

    matrix = oops.Frame.as_frame(obs.frame).transform_at_time(0.).matrix
    (ra, dec, _) = (matrix.unrotate(oops.Vector3.ZAXIS)).to_ra_dec_length()
    return (ra.vals * oops.DPR % 360., dec.vals * oops.DPR)


def frame_matrix(obs: oops.observation.Snapshot) -> np.ndarray:
    """The rotation matrix from J2000 to the observation's frame."""

    return np.asarray(oops.Frame.as_frame(obs.frame).transform_at_time(0.).matrix.vals)


def wcs_frame_matrix(path: pathlib.Path, clock: float) -> np.ndarray:
    """The rotation matrix to the frame whose optic axis is at the WCS position of the
    image center, with celestial north at the given clock angle in degrees."""

    (ra, dec) = WCS(fits.getheader(path)).all_pix2world(512.5, 512.5, 1)
    cmatrix = oops.frame.Cmatrix.from_ra_dec(float(ra), float(dec), clock)
    return np.asarray(cmatrix.transform_at_time(0.).matrix.vals)


@pytest.mark.parametrize(('pointing', 'utc', 'clock'), [
    ('fits', '2007-010T05:41:10.863', 275.),
    ('fits90', '2007-010T05:41:10.863', 365.),
    ('fits90', '2013-001T00:00:00.000', 275.),
], ids=['fits', 'fits90 before 2013', 'fits90 after 2012'])
def test_north_clock_angle(tmp_path: pathlib.Path, pointing: str, utc: str,
                           clock: float) -> None:
    """The frame from the header has north at SPCEMEN, plus 90 degrees for "fits90" on
    an image taken before 2013."""

    path = write_fits(tmp_path / 'a.fit', SPCUTCID=utc)
    obs = single(NewHorizonsLORRI.from_file(path, pointing=pointing))
    assert np.allclose(frame_matrix(obs), wcs_frame_matrix(path, clock), atol=1.e-12)


def test_fitsapp_removes_the_aberration(tmp_path: pathlib.Path,
                                        sun_path: oops.Path) -> None:
    """With "fitsapp", the optic axis moves from the header's by the aberration of the
    spacecraft's motion, a few arcseconds."""

    path = write_fits(tmp_path / 'a.fit', SPCSSCVX=-30., SPCSSCVY=0., SPCSSCVZ=0.)
    fits_axis = boresight(single(NewHorizonsLORRI.from_file(path, geom='fits',
                                                            pointing='fits')))
    app_axis = boresight(single(NewHorizonsLORRI.from_file(path, geom='fits',
                                                           pointing='fitsapp')))
    shift = np.hypot((app_axis[0] - fits_axis[0]) * np.cos(fits_axis[1] * oops.RPD),
                     app_axis[1] - fits_axis[1])
    assert 1.e-5 < shift < 1.e-2


@pytest.mark.parametrize(('crval', 'pixel'), [
    ((237.5, -33.8), (700., 300.)),
    ((40., 89.95), (512.5, 1024.)),
], ids=['away from the pole', 'beyond the pole'])
def test_radec_from_uv_matches_the_wcs(crval: tuple[float, float],
                                       pixel: tuple[float, float]) -> None:
    """The gnomonic projection matches astropy's, also for a point beyond a celestial
    pole, where denom is negative."""

    header = dict(FITS_HEADER, CRVAL1=crval[0], CRVAL2=crval[1])
    (ra, dec) = NewHorizonsLORRI._radec_from_uv(pixel[0], pixel[1], header)
    wcs = WCS(fits.Header(dict(header, CTYPE1='RA---TAN', CTYPE2='DEC--TAN')))
    (wcs_ra, wcs_dec) = wcs.all_pix2world(pixel[0], pixel[1], 1)
    assert ra % 360. == pytest.approx(float(wcs_ra), abs=1.e-9)
    assert dec == pytest.approx(float(wcs_dec), abs=1.e-9)


def test_no_pointing_check(tmp_path: pathlib.Path,
                           monkeypatch: pytest.MonkeyPatch) -> None:
    """An image without C-kernel pointing loads; its frame raises the SPICE error when
    it is used."""

    def no_pointing(self: FakeSpiceFrame, time: Any, *, quick: Any = None) -> None:
        raise KeyError('SPICE(NOFRAMECONNECT)')

    monkeypatch.setattr(FakeSpiceFrame, 'transform_at_time', no_pointing)
    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    with pytest.raises(KeyError, match='NOFRAMECONNECT'):
        oops.Frame.as_frame(obs.frame).wrt(oops.Frame.J2000).transform_at_time(
            obs.midtime)


def test_frame_override_is_used_as_given(tmp_path: pathlib.Path) -> None:
    """A Cmatrix given as the frame is the observation's frame."""

    cmatrix = oops.frame.Cmatrix(oops.Matrix3.IDENTITY)
    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            frame=cmatrix))
    assert obs.frame is cmatrix


def test_navigation_override_wraps_the_frame(tmp_path: pathlib.Path) -> None:
    """A navigation override wraps the default frame in a Navigation frame."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            navigation=(0., 0.)))
    assert isinstance(obs.frame, oops.frame.Navigation)


def test_target_override(tmp_path: pathlib.Path) -> None:
    """A target override replaces the target body."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            target='PLUTO'))
    assert obs.target == 'PLUTO'

##########################################################################################
# The I/F calibration
##########################################################################################

AU = 149597870.7


def test_iof_calibration(tmp_path: pathlib.Path, jupiter_path: oops.Path) -> None:
    """A calibrated product of a target with an ephemeris gets one solar I/F
    calibration, pi r^2 / (texp RSOLAR 176)."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    (calib,) = obs.calibrations
    assert calib.name == 'I/F'
    assert calib.factor.vals == pytest.approx(np.pi * 25. / (5. * 5.e6 * 176.))


def test_iof_calibration_of_a_label_read(tmp_path: pathlib.Path,
                                         jupiter_path: oops.Path) -> None:
    """A label read gets the calibration from the FITS header that it reads."""

    obs = single(NewHorizonsLORRI.from_file(write_product(tmp_path)))
    assert [calib.name for calib in obs.calibrations] == ['I/F']


def test_iof_calibration_of_a_target_body(tmp_path: pathlib.Path,
                                          jupiter_path: oops.Path) -> None:
    """A target override given as a Body is placed by its path."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            target=oops.Body.lookup('JUPITER')))
    assert [calib.name for calib in obs.calibrations] == ['I/F']


def test_iof_calibration_at_the_cadence_override(tmp_path: pathlib.Path,
                                                 sun_path: oops.Path) -> None:
    """The target's distance is taken at the middle of a cadence override."""

    oops.path.LinearPath((oops.Vector3((AU, 0., 0.)), oops.Vector3((AU, 0., 0.))), 0.,
                         'SUN', path_id='JUPITER')
    cadence = oops.cadence.SnapCadence(1., 0.)
    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            cadence=cadence))
    assert obs.calibrations[0].factor.vals == pytest.approx(
        np.pi * 4. / (5. * 5.e6 * 176.))


def test_no_iof_calibration_for_a_raw_product(tmp_path: pathlib.Path,
                                              jupiter_path: oops.Path) -> None:
    """A header without RSOLAR gives no calibration."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit',
                                                       RSOLAR=None)))
    assert obs.calibrations == []


def test_no_iof_calibration_for_a_zero_exposure(tmp_path: pathlib.Path,
                                                jupiter_path: oops.Path) -> None:
    """An exposure of zero gives no calibration."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit', EXPTIME=0.)))
    assert obs.calibrations == []


def test_no_iof_calibration_without_a_path(tmp_path: pathlib.Path) -> None:
    """A target without a registered path gives no calibration."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.calibrations == []


def test_no_iof_calibration_without_an_ephemeris(tmp_path: pathlib.Path,
                                                 jupiter_path: oops.Path,
                                                 monkeypatch: pytest.MonkeyPatch) -> None:
    """A target whose position the SP kernels do not give gives no calibration."""

    def no_ephemeris(self: oops.Path, time: Any, **kwargs: Any) -> None:
        raise OSError('SPICE(SPKINSUFFDATA)')

    monkeypatch.setattr(oops.path.FixedPath, 'event_at_time', no_ephemeris)
    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.calibrations == []


def test_no_iof_calibration_for_a_label_alone(tmp_path: pathlib.Path,
                                              jupiter_path: oops.Path) -> None:
    """An astrometry read of a label reads no RSOLAR, so it gives no calibration."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT')
    obs = single(NewHorizonsLORRI.from_file(label, astrometry=True))
    assert obs.calibrations == []


def test_empty_calibrations_turn_calibration_off(tmp_path: pathlib.Path,
                                                 jupiter_path: oops.Path) -> None:
    """calibrations=[] gives no calibration."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'),
                                            calibrations=[]))
    assert obs.calibrations == []

##########################################################################################
# Subfields
##########################################################################################


def test_spice_frame_subfields(tmp_path: pathlib.Path) -> None:
    """The SPICE frame of the camera and its rotation to the camera frame."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.spice_frame_name == 'NH_LORRI'
    assert obs.spice_frame_id == -98300


def test_spice_to_frame(tmp_path: pathlib.Path) -> None:
    """spice_to_frame flips the y and z axes."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert np.array_equal(obs.spice_to_frame.vals, np.diag([1, -1, -1]))


def test_ck_subfields(tmp_path: pathlib.Path) -> None:
    """The C kernel IDs of the pointing, whose frame is never frozen."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert obs.ck_object_id == -98000
    assert obs.ck_sclk_id == -98
    assert obs.spice_frame_is_frozen is False


def test_ck_sclk_id_is_the_clock_of_the_ck_object(tmp_path: pathlib.Path) -> None:
    """SPICE gives ck_sclk_id as the clock of ck_object_id."""

    obs = single(NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit')))
    assert cspyce.ckmeta(obs.ck_object_id, 'SCLK') == obs.ck_sclk_id


def test_path_subfields_of_a_fits_read(tmp_path: pathlib.Path) -> None:
    """A FITS read names the FITS file and no label."""

    path = write_fits(tmp_path / 'a.fit')
    obs = single(NewHorizonsLORRI.from_file(path))
    assert obs.basename == 'a.fit'
    assert obs.image_url == path.as_posix()
    assert obs.image_abspath == path.resolve()
    assert {'label_url', 'label_abspath'} & set(obs.subfields) == set()


def test_path_subfields_of_a_label_read(tmp_path: pathlib.Path) -> None:
    """A label read names the FITS file, which it read, and the label."""

    label = write_product(tmp_path)
    obs = single(NewHorizonsLORRI.from_file(label))
    assert obs.image_abspath == (tmp_path / 'lor_0030713591_0x633_sci.fit').resolve()
    assert obs.label_url == label.as_posix()
    assert obs.label_abspath == label.resolve()


def test_archive_label_path_subfields(tmp_path: pathlib.Path) -> None:
    """The FITS file of an archive label is found beside it in the label's case."""

    for suffix in ('lbl', 'fit'):
        shutil.copy(TEST_FILES / f'{SCI}.{suffix}', tmp_path / f'{SCI}.{suffix}')
    obs = single(NewHorizonsLORRI.from_file(tmp_path / f'{SCI}.lbl'))
    assert obs.basename == f'{SCI}.fit'

##########################################################################################
# select
##########################################################################################


def test_select_the_observation(tmp_path: pathlib.Path) -> None:
    """select=0 chooses the file's one observation."""

    obs = NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), select=0)
    assert isinstance(obs, oops.observation.Snapshot)


def test_select_a_slice(tmp_path: pathlib.Path) -> None:
    """A slice chooses a tuple of the observations."""

    obs = NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), select=slice(None))
    assert isinstance(obs, tuple)


def test_select_out_of_range(tmp_path: pathlib.Path) -> None:
    """An index beyond the one observation raises OopsIndexError."""

    with pytest.raises(oops.OopsIndexError, match='select index out of range: 1'):
        NewHorizonsLORRI.from_file(write_fits(tmp_path / 'a.fit'), select=1)

##########################################################################################
