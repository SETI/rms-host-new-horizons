##########################################################################################
# tests/oops/test_lorri_files.py
##########################################################################################
"""Tests that a LORRI FITS file, its PDS3 label, and its index row give the same Snapshot
wherever the label and the FITS header agree, and that a label's Snapshot follows the
resolution rules where they do not.

The files in ``test_files`` are cut down from NHJULO_1001 and NHJULO_2001,
data/20070110_003071/lor_0030713591_0x633: the raw (``eng``) product is unchanged, the
calibrated (``sci``) product keeps the first 16 rows of each array and zeros the rest,
and the index and supplement of NHJULO_2001 keep rows 21-50. Each file is renamed with
an ``_edited`` suffix, and the labels name the renamed files.
"""

import pathlib
from typing import Any

import numpy as np
import oops
import pytest
from astropy.io import fits
from conftest import TEST_FILES

from host_new_horizons.oops.lorri import NewHorizonsLORRI

STEM = 'lor_0030713591_0x633'
INDEX = TEST_FILES / 'NHJULO_2001_index_edited.lbl'

PRODUCTS = ['eng', 'sci']
INPUTS = ['lbl', 'host_fit', 'host_lbl']

# The index and the label give the start of the exposure to the millisecond
INDEX_TIME_TOLERANCE = 1.e-3    # seconds

# The corners and center of the 256x256 field of view
CORNERS = oops.Pair([(0, 0), (256, 0), (0, 256), (256, 256), (128, 128)])


def load(product: str, source: str) -> Any:
    """The Snapshot of the product ('eng' or 'sci') from one source.

    The source is 'fit' or 'lbl' through NewHorizonsLORRI.from_file(), or 'host_fit' or
    'host_lbl' through oops.Host.from_file().
    """

    suffix = source.rpartition('_')[2]
    path = TEST_FILES / f'{STEM}_{product}_edited.{suffix}'
    if source.startswith('host'):
        return oops.Host.from_file(path)
    return NewHorizonsLORRI.from_file(path)


def label_value(product: str, keyword: str) -> str:
    """The quoted value of a keyword in the product's PDS3 label."""

    text = (TEST_FILES / f'{STEM}_{product}_edited.lbl').read_text()
    line = next(line for line in text.splitlines() if line.startswith(keyword + ' '))
    return line.partition('"')[2].partition('"')[0]


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
@pytest.mark.parametrize('attribute', ['time', 'texp', 'target', 'lightsource',
                                       'binning_mode', 'spacecraft_clock_start_count',
                                       'spacecraft_clock_stop_count', 'filepath',
                                       'basename', 'abspath', 'image_url',
                                       'header_dict', 'spice_kernels', 'spice_frame_name',
                                       'spice_frame_id', 'calibrations', 'shape'])
def test_attribute_matches_the_fits_file(product: str, source: str,
                                         attribute: str) -> None:
    assert getattr(load(product, source), attribute) == getattr(load(product, 'fit'),
                                                                attribute)


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
def test_data_matches_the_fits_file(product: str, source: str) -> None:
    assert np.array_equal(load(product, source).data, load(product, 'fit').data)


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
def test_headers_match_the_fits_file(product: str, source: str) -> None:
    headers = [str(header) for header in load(product, source).headers]
    assert headers == [str(header) for header in load(product, 'fit').headers]


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
def test_fov_matches_the_fits_file(product: str, source: str) -> None:
    los = load(product, source).fov.los_from_uv(CORNERS).vals
    assert np.array_equal(los, load(product, 'fit').fov.los_from_uv(CORNERS).vals)


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
def test_path_and_frame_match_the_fits_file(product: str, source: str) -> None:
    obs = load(product, source)
    assert obs.path.path_id == 'NEW HORIZONS'
    assert obs.frame.frame_id == 'NH_LORRI'


@pytest.mark.parametrize('product', PRODUCTS)
def test_path_subfields_name_the_fits_file(product: str) -> None:
    obs = load(product, 'lbl')
    assert obs.basename == f'{STEM}_{product}_edited.fit'


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize(('subfield', 'keyword'), [
    ('spacecraft_clock_start_count', 'SPACECRAFT_CLOCK_START_COUNT'),
    ('spacecraft_clock_stop_count', 'SPACECRAFT_CLOCK_STOP_COUNT'),
])
def test_clock_counts_match_the_label(product: str, subfield: str, keyword: str) -> None:
    assert load(product, 'fit').subfields[subfield] == label_value(product, keyword)


@pytest.mark.parametrize('product', PRODUCTS)
def test_data_is_the_primary_array(product: str) -> None:
    with fits.open(TEST_FILES / f'{STEM}_{product}_edited.fit') as hdulist:
        expected = hdulist[0].data.copy()
    assert np.array_equal(load(product, 'fit').data, expected)


def test_raw_product_shape() -> None:
    obs = load('eng', 'fit')
    assert obs.data.shape == (256, 257)
    assert obs.data.dtype == np.dtype('>i2')


def test_raw_product_has_no_error_or_quality() -> None:
    obs = load('eng', 'fit')
    assert 'error' not in obs.subfields
    assert 'quality' not in obs.subfields


@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize('source', INPUTS)
def test_dict_matches_the_fits_file(product: str, source: str) -> None:
    from_fits = load(product, 'fit').dict
    dict_ = load(product, source).dict
    del from_fits['TARGET_NAME']
    assert {key: dict_[key] for key in from_fits} == from_fits


# The two labels name the target differently ("J17 CALLIRRHOE" and "CALLIRRHOE"); the
# detached label's name wins, and both give the same target body
@pytest.mark.parametrize('product', PRODUCTS)
@pytest.mark.parametrize(('source', 'name'), [
    ('fit', 'CALLIRRHOE'), ('host_fit', 'CALLIRRHOE'),
    ('lbl', 'J17 CALLIRRHOE'), ('host_lbl', 'J17 CALLIRRHOE')])
def test_target_name_of_each_input(product: str, source: str, name: str) -> None:
    obs = load(product, source)
    assert obs.label_target == name
    assert obs.dict['TARGET_NAME'] == name


@pytest.mark.parametrize('product', PRODUCTS)
def test_dict_of_the_label_adds_the_label(product: str) -> None:
    obs = load(product, 'lbl')
    assert obs.dict['PRODUCT_ID'] == f'LOR_0030713591_0X633_{product.upper()}'
    assert obs.label_dict['TARGET_NAME'] == 'J17 CALLIRRHOE'


def test_calibrated_product_error_and_quality() -> None:
    obs = load('sci', 'fit')
    with fits.open(TEST_FILES / f'{STEM}_sci_edited.fit') as hdulist:
        assert np.array_equal(obs.error, hdulist[1].data)
        assert np.array_equal(obs.quality, hdulist[2].data)


# The label's catalog name "J17 CALLIRRHOE" translates to the body name the header gives
@pytest.mark.parametrize('source', ['fit', 'lbl'])
def test_target(source: str) -> None:
    obs = load('sci', source)
    assert obs.target == 'CALLIRRHOE'
    assert label_value('sci', 'TARGET_NAME') == 'J17 CALLIRRHOE'

##########################################################################################
# The index row of the same image
##########################################################################################


def index_row() -> Any:
    """The Snapshot of the calibrated product from its row of the index."""

    return next(obs for obs in NewHorizonsLORRI.from_index(INDEX)
                if obs.basename == f'{STEM}_sci.fit')


def test_index_row_names_the_fits_file() -> None:
    assert index_row().filepath == f'NHJULO_2001/data/20070110_003071/{STEM}_sci.fit'


@pytest.mark.parametrize('attribute', ['binning_mode', 'spacecraft_clock_start_count',
                                       'spacecraft_clock_stop_count', 'texp'])
def test_index_row_matches_the_fits_file(attribute: str) -> None:
    assert getattr(index_row(), attribute) == getattr(load('sci', 'fit'), attribute)


def test_index_row_time_matches_the_fits_file() -> None:
    assert index_row().time == pytest.approx(load('sci', 'fit').time,
                                             abs=INDEX_TIME_TOLERANCE)


def test_index_row_fov_matches_the_fits_file() -> None:
    assert index_row().fov is load('sci', 'fit').fov


def test_index_row_target_matches_the_fits_file() -> None:
    assert index_row().target == load('sci', 'fit').target


def test_index_row_target_is_a_python_string() -> None:
    assert type(index_row().target) is str

##########################################################################################
# An astrometry read of the label alone
##########################################################################################


@pytest.fixture
def label_alone(tmp_path: pathlib.Path) -> pathlib.Path:
    """A copy of the calibrated product's label, in a directory without its FITS file."""

    path = tmp_path / f'{STEM}_sci_edited.lbl'
    path.write_bytes((TEST_FILES / path.name).read_bytes())
    return path


def test_label_astrometry_needs_no_fits_file(label_alone: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(label_alone, astrometry=True)
    assert obs.basename == f'{STEM}_sci_edited.fit'


@pytest.mark.parametrize('attribute', ['binning_mode', 'spacecraft_clock_start_count',
                                       'spacecraft_clock_stop_count', 'fov'])
def test_label_astrometry_matches_the_fits_file(label_alone: pathlib.Path,
                                                attribute: str) -> None:
    obs = NewHorizonsLORRI.from_file(label_alone, astrometry=True)
    assert getattr(obs, attribute) == getattr(load('sci', 'fit'), attribute)


def test_label_astrometry_time_matches_the_fits_file(label_alone: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(label_alone, astrometry=True)
    assert obs.time == pytest.approx(load('sci', 'fit').time, abs=INDEX_TIME_TOLERANCE)


def test_label_astrometry_target_is_the_label_target(label_alone: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(label_alone, astrometry=True)
    assert obs.target == 'CALLIRRHOE'

##########################################################################################
