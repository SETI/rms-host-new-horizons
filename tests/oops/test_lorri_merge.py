##########################################################################################
# tests/oops/test_lorri_merge.py
##########################################################################################
"""Tests of what NewHorizonsLORRI.from_file() reads from each input, and of how it merges
the detached PDS3 label with the FITS header."""

import pathlib
from typing import Any

import julian
import oops
import pytest
from conftest import FITS_HEADER, write_fits, write_label
from filecache import FCPath
from pdsparser import Pds3Label

from host_new_horizons.oops.lorri import NewHorizonsLORRI

# The start of the exposure in the synthetic FITS header, and in the label, which gives
# it rounded to the millisecond
FITS_START = FITS_HEADER['SPCSCET'] - FITS_HEADER['EXPTIME'] / 2.
LABEL_START = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T05:41:08.363'))


def label_and_fits(tmp_path: pathlib.Path, *, fits_keywords: dict[str, Any] | None = None,
                   lines: int | None = 1024, **label_keywords: Any) -> pathlib.Path:
    """Write a synthetic FITS file and its label in `tmp_path`; return the label's path.

    The label keywords override those of write_label(); `fits_keywords` override the
    FITS header.
    """

    write_fits(tmp_path / 'lor_test.fit', **(fits_keywords or {}))
    return write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=lines,
                       **label_keywords)


def from_label(tmp_path: pathlib.Path, **kwargs: Any) -> Any:
    """The Snapshot of a full read of a label written by label_and_fits()."""

    return NewHorizonsLORRI.from_file(label_and_fits(tmp_path, **kwargs))

##########################################################################################
# What each input reads
##########################################################################################


# The subfields present after each of the four reads
@pytest.mark.parametrize(('suffix', 'astrometry', 'present'), [
    ('fit', False, {'data', 'error', 'quality', 'header_dict', 'headers'}),
    ('fit', True, {'header_dict', 'headers'}),
    ('lbl', False, {'data', 'error', 'quality', 'header_dict', 'headers', 'label_dict'}),
    ('lbl', True, {'label_dict'}),
])
def test_subfields_of_each_read(tmp_path: pathlib.Path, suffix: str, astrometry: bool,
                                present: set[str]) -> None:
    label_and_fits(tmp_path)
    obs = NewHorizonsLORRI.from_file(tmp_path / f'lor_test.{suffix}',
                                     astrometry=astrometry)
    optional = {'data', 'error', 'quality', 'header_dict', 'headers', 'label_dict'}
    assert optional & set(obs.subfields) == present


def test_label_dict_is_the_label(tmp_path: pathlib.Path) -> None:
    label = label_and_fits(tmp_path)
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.label_dict == Pds3Label(label, sources=True).as_dict()


def test_dict_of_a_label_holds_both_sources(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path)
    assert obs.dict['INSTRUMENT_ID'] == 'LORRI'
    assert obs.dict['MISSION'] == 'New Horizons'


def test_dict_of_a_label_matches_the_fits_file_where_they_agree(
        tmp_path: pathlib.Path) -> None:
    label = label_and_fits(tmp_path)
    from_fits = NewHorizonsLORRI.from_file(tmp_path / 'lor_test.fit').dict
    dict_ = NewHorizonsLORRI.from_file(label).dict
    assert {key: dict_[key] for key in from_fits} == from_fits


def test_dict_start_time_is_an_iso_string(tmp_path: pathlib.Path) -> None:
    assert from_label(tmp_path).dict['START_TIME'] == '2007-01-10T05:41:08.363'


def test_label_astrometry_on_a_remote_label_names_no_local_file(
        tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT'))
    label.filepath = FCPath('https://example.org/volume/lor_test.lbl')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    assert 'abspath' not in obs.subfields
    assert obs.image_url == 'https://example.org/volume/lor_test.fit'


def test_parsed_label_of_another_instrument(tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_label(tmp_path / 'mvic.lbl', 'MVIC.FIT',
                                  INSTRUMENT_ID='"MVIC"'))
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI PDS3'):
        NewHorizonsLORRI.from_file(label)


def test_parsed_object_of_another_kind() -> None:
    with pytest.raises(oops.host.WrongHostError, match='label or FITS file: 42'):
        NewHorizonsLORRI.from_file(42)

##########################################################################################
# Merging the label and the FITS header
##########################################################################################


@pytest.mark.parametrize(('fits_start', 'expected'), [
    (LABEL_START + 0.0004, LABEL_START + 0.0004),     # the label rounds it: the header's
    (LABEL_START - 0.0004, LABEL_START - 0.0004),
    (LABEL_START + 0.002, LABEL_START),                # they differ: the label's
])
def test_start_time_rule(tmp_path: pathlib.Path, fits_start: float,
                         expected: float) -> None:
    spcscet = fits_start + FITS_HEADER['EXPTIME'] / 2.
    obs = from_label(tmp_path, fits_keywords={'SPCSCET': spcscet})
    assert obs.time[0] == pytest.approx(expected, abs=1.e-7)


@pytest.mark.parametrize(('exptime', 'expected'), [
    (4.9996, 4.9996),       # the label rounds it: the header's
    (4.9, 5.0),             # they differ: the label's
])
def test_exposure_duration_rule(tmp_path: pathlib.Path, exptime: float,
                                expected: float) -> None:
    obs = from_label(tmp_path, fits_keywords={'EXPTIME': exptime})
    assert obs.texp == expected


@pytest.mark.parametrize(('label_target', 'fits_target', 'expected'), [
    ('J17 CALLIRRHOE', 'CALLIRRHOE', 'CALLIRRHOE'),    # only the header's is a body
    ('JUPITER', 'FJUPITER', 'JUPITER'),                # only the label's is a body
    ('J RINGS', 'ANSA', 'J RINGS'),                    # neither: the label's
    ('PLUTO', 'JUPITER', 'PLUTO'),                     # both: the label's
])
def test_target_name_rule(tmp_path: pathlib.Path, label_target: str, fits_target: str,
                          expected: str) -> None:
    obs = from_label(tmp_path, TARGET_NAME=f'"{label_target}"',
                     fits_keywords={'SPCCBTNM': fits_target})
    assert obs.label_target == expected
    assert obs.dict['TARGET_NAME'] == expected


@pytest.mark.parametrize(('keyword', 'value', 'subfield'), [
    ('SPCSCLK0', '0030713585:00000', 'spacecraft_clock_start_count'),
    ('STOPMET', 30713599.0, 'spacecraft_clock_stop_count'),
])
def test_clock_count_rule(tmp_path: pathlib.Path, keyword: str, value: Any,
                          subfield: str) -> None:
    obs = from_label(tmp_path, fits_keywords={keyword: value})
    label_counts = {'spacecraft_clock_start_count': '0030713585:49000',
                    'spacecraft_clock_stop_count': '0030713590:49000'}
    assert obs.subfields[subfield] == label_counts[subfield]


def test_lines_rule(tmp_path: pathlib.Path) -> None:
    write_fits(tmp_path / 'lor_test.fit', lines=256, samples=256)
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=1024)
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.binning_mode == '4X4'


@pytest.mark.parametrize(('keyword', 'subfield', 'value'), [
    ('SPACECRAFT_CLOCK_START_COUNT', 'spacecraft_clock_start_count', '0030713585:49000'),
    ('SPACECRAFT_CLOCK_STOP_COUNT', 'spacecraft_clock_stop_count', '0030713590:49000'),
    ('TARGET_NAME', 'label_target', 'JUPITER'),
])
@pytest.mark.parametrize('label_value', [None, '"UNK"', '"N/A"'])
def test_item_missing_from_the_label_comes_from_the_header(
        tmp_path: pathlib.Path, keyword: str, subfield: str, value: str,
        label_value: str | None) -> None:
    obs = from_label(tmp_path, **{keyword: label_value})
    assert obs.subfields[subfield] == value


@pytest.mark.parametrize(('keyword', 'subfield', 'value'), [
    ('SPCSCLK0', 'spacecraft_clock_start_count', '0030713585:49000'),
    ('STOPMET', 'spacecraft_clock_stop_count', '0030713590:49000'),
    ('SPCCBTNM', 'label_target', 'JUPITER'),
])
def test_item_missing_from_the_header_comes_from_the_label(
        tmp_path: pathlib.Path, keyword: str, subfield: str, value: str) -> None:
    obs = from_label(tmp_path, fits_keywords={keyword: None})
    assert obs.subfields[subfield] == value


def test_times_missing_from_the_label_come_from_the_header(
        tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, START_TIME=None, EXPOSURE_DURATION=None, lines=None)
    assert obs.time[0] == FITS_START


def test_item_missing_from_both_sources_is_absent(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, SPACECRAFT_CLOCK_STOP_COUNT='"N/A"',
                     fits_keywords={'STOPMET': None})
    assert 'spacecraft_clock_stop_count' not in obs.subfields
    assert 'SPACECRAFT_CLOCK_STOP_COUNT' not in obs.dict


def test_null_values_stay_in_the_raw_dicts(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, TARGET_NAME='"N/A"', fits_keywords={'SPCCBTNM': '---'})
    assert obs.label_dict['TARGET_NAME'] == 'N/A'
    assert obs.header_dict['SPCCBTNM'] == '---'
    assert obs.target == 'NONE'

##########################################################################################
# Values the Snapshot needs
##########################################################################################


@pytest.mark.parametrize(('keyword', 'message'), [
    ('EXPTIME', 'no EXPOSURE_DURATION for LORRI image .*lor_test.fit'),
    ('SPCSCET', 'no START_TIME for LORRI image .*lor_test.fit'),
])
def test_fits_file_without_a_needed_item(tmp_path: pathlib.Path, keyword: str,
                                         message: str) -> None:
    keywords: dict[str, Any] = {keyword: None}
    path = write_fits(tmp_path / 'lor_test.fit', **keywords)
    with pytest.raises(oops.host.HostFormatError, match=message):
        NewHorizonsLORRI.from_file(path)


@pytest.mark.parametrize(('keywords', 'lines', 'message'), [
    ({'START_TIME': None}, 1024, 'no START_TIME for LORRI image .*lor_test.lbl'),
    ({'EXPOSURE_DURATION': None}, 1024, 'no EXPOSURE_DURATION for LORRI image'),
    ({}, None, 'no lines of the image'),
])
def test_label_astrometry_without_a_needed_item(tmp_path: pathlib.Path,
                                                keywords: dict[str, Any],
                                                lines: int | None, message: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=lines,
                        **keywords)
    with pytest.raises(oops.host.HostFormatError, match=message):
        NewHorizonsLORRI.from_file(label, astrometry=True)

##########################################################################################
# Parsing the label
##########################################################################################


# "fast" rejects a quoted value that starts at the beginning of the next line; "strict"
# rejects a set of unquoted words with a blank inside, which "loose" accepts
@pytest.mark.parametrize('extra', [
    ('NOTE =', '"a note that starts', 'on the next line"'),
    ('NOTE =', '"a note that starts', 'on the next line"', 'WORDS = {a, b c}'),
])
def test_label_parser_falls_back(tmp_path: pathlib.Path, extra: tuple[str, ...]) -> None:
    write_fits(tmp_path / 'lor_test.fit')
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', extra=extra)
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.basename == 'lor_test.fit'


def test_label_parser_starts_at_the_given_method(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT',
                        extra=('NOTE = "a note" junk',))
    with pytest.raises(oops.host.PdsHostError, match=r'method strict: .*method loose: '):
        NewHorizonsLORRI.from_file(label, astrometry=True,
                                   pds3_options={'method': 'strict'})


# An hour of 25 makes the fast method raise a ValueError rather than a PdsError
@pytest.mark.parametrize('start_time', ['2007-01-10T25:41:08.363', '"unterminated'])
def test_unreadable_label(tmp_path: pathlib.Path, start_time: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', START_TIME=start_time)
    with pytest.raises(oops.host.PdsHostError,
                       match=r'PDS3 label cannot be parsed: .*method fast: '):
        NewHorizonsLORRI.from_file(label, astrometry=True)


@pytest.mark.parametrize(('options', 'message'), [
    ({'method': 'quick'}, "invalid PDS3 parser method 'quick'"),
    ({'bogus': 1}, 'invalid PDS3 parser option'),
])
def test_invalid_label_parser_options(tmp_path: pathlib.Path, options: dict[str, Any],
                                      message: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    with pytest.raises(oops.OopsValueError, match=message):
        NewHorizonsLORRI.from_file(label, pds3_options=options)

##########################################################################################
# Missing pointing
##########################################################################################


def no_pointing(monkeypatch: pytest.MonkeyPatch, message: str) -> None:
    """Make the camera frame raise RuntimeError(message) when it is evaluated."""

    def transform_at_time(self: Any, *args: Any, **kwargs: Any) -> Any:
        raise RuntimeError(message)

    monkeypatch.setattr(oops.frame.Cmatrix, 'transform_at_time', transform_at_time)


def test_missing_pointing_raises_host_error(tmp_path: pathlib.Path,
                                            monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT) -- sxform -- insufficient data')
    with pytest.raises(oops.host.HostError,
                       match=r'no New Horizons LORRI pointing for .*lor_test.fit: '
                             r'SPICE\(NOFRAMECONNECT\)'):
        NewHorizonsLORRI.from_file(path)


def test_other_pointing_errors_pass_through(tmp_path: pathlib.Path,
                                            monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'something else')
    with pytest.raises(RuntimeError, match='something else'):
        NewHorizonsLORRI.from_file(path)


@pytest.mark.parametrize('options', [{'frame': 'J2000'}, {'pointing': 'fits'}])
def test_pointing_is_not_checked_without_the_c_kernels(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch,
        options: dict[str, Any]) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT)')
    obs = NewHorizonsLORRI.from_file(path, **options)
    assert obs.basename == 'lor_test.fit'

##########################################################################################
