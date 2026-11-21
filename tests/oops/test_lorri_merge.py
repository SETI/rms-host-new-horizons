##########################################################################################
# tests/oops/test_lorri_merge.py
##########################################################################################
"""Tests of what NewHorizonsLORRI.from_file() reads from each input, and of how it merges
the detached PDS3 label with the FITS header."""

import json
import pathlib
from typing import Any

import julian
import oops
import pytest
from astropy.io import fits
from conftest import FITS_HEADER, FakeSpiceFrame, write_fits, write_label
from filecache import FCPath
from pdsparser import Pds3Label

from host_new_horizons.oops.lorri import (_FITS_ITEM_KEYWORDS, _WINNING_SOURCE,
                                          NewHorizonsLORRI)

# The start of the exposure in the synthetic FITS header, and in the label, which gives
# it rounded to the millisecond
FITS_START = FITS_HEADER['SPCSCET'] - FITS_HEADER['EXPTIME'] / 2.
LABEL_START = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T05:41:08.363'))

# The middle of the exposure in the synthetic FITS header
MIDTIME = FITS_HEADER['SPCSCET']


def label_and_fits(tmp_path: pathlib.Path, *, fits_keywords: dict[str, Any] | None = None,
                   lines: int | None = 1024, **label_keywords: Any) -> pathlib.Path:
    """Write a synthetic FITS file and its label in `tmp_path`; return the label's path.

    The label keywords override those of write_label(); `fits_keywords` override the
    FITS header.
    """

    write_fits(tmp_path / 'lor_test.fit', **(fits_keywords or {}))
    return write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=lines,
                       **label_keywords)


def from_label(tmp_path: pathlib.Path, *, geom: str = 'spice', **kwargs: Any) -> Any:
    """The Snapshot of a full read of a label written by label_and_fits()."""

    return NewHorizonsLORRI.from_file(label_and_fits(tmp_path, **kwargs), geom=geom)

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
    expected = Pds3Label(label, sources=False).as_dict()
    assert obs.label_dict.keys() == {key for key in expected
                                     if key == key.upper()} - {'END'}


def helper_keys(dict_: dict[str, Any]) -> list[str]:
    """The keys that the PDS3 parser adds, at any depth: those with a lower-case suffix,
    such as "_fmt", "_unit" and "_offset", and "objects"."""

    keys = []
    for (key, value) in dict_.items():
        if key != key.upper():
            keys.append(key)
        if isinstance(value, dict):
            keys += helper_keys(value)
    return keys


# The label gives a unit, a time, a pointer, an object, and a NOTE that the parser unwraps
@pytest.mark.parametrize('subfield', ['dict', 'label_dict'])
def test_parser_helper_keys_are_left_out(tmp_path: pathlib.Path, subfield: str) -> None:
    obs = from_label(tmp_path, extra=('NOTE = "a note"',))
    assert helper_keys(obs.subfields[subfield]) == []


# A caller's sources=True is overruled, so that no "_source" key appears anywhere
def test_label_is_parsed_without_source_text(tmp_path: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(label_and_fits(tmp_path),
                                     pds3_options={'sources': True})
    keys = {*obs.dict, *obs.label_dict}
    assert [key for key in keys if key.endswith('_source')] == []


# A keyword that a header repeats is kept at each occurrence, as "NAME[n]"
@pytest.mark.parametrize('subfield', ['header_dict', 'headers'])
def test_header_keeps_a_repeated_keyword(tmp_path: pathlib.Path, subfield: str) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    with fits.open(path, mode='update') as hdulist:
        hdulist[0].header.append(('REFFILE', 'a.fit'))
        hdulist[0].header.append(('REFFILE', 'b.fit'))
    obs = NewHorizonsLORRI.from_file(path)
    header = obs.header_dict if subfield == 'header_dict' else obs.headers[0]
    assert [header.get(key) for key in ('REFFILE', 'REFFILE[0]', 'REFFILE[1]')] == [
        None, 'a.fit', 'b.fit']


def test_label_dict_writes_times_as_the_label_does(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, PRODUCT_CREATION_TIME='2016-10-31T00:00:00')
    assert obs.label_dict['PRODUCT_CREATION_TIME'] == '2016-10-31T00:00:00'


def test_label_dict_writes_a_set_as_a_list(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, extra=('SPICE_FILE_NAME = {"b.bc", "a.bsp"}',))
    assert obs.label_dict['SPICE_FILE_NAME'] == ['b.bc', 'a.bsp']


@pytest.mark.parametrize(('suffix', 'astrometry'), [
    ('fit', False), ('fit', True), ('lbl', False), ('lbl', True)])
@pytest.mark.parametrize('subfield', ['dict', 'label_dict', 'header_dict', 'headers'])
def test_metadata_subfields_are_json(tmp_path: pathlib.Path, suffix: str,
                                     astrometry: bool, subfield: str) -> None:
    """Each metadata subfield of each read survives a JSON round trip; a subfield that a
    read does not give is checked as empty."""

    write_fits(tmp_path / 'lor_test.fit', COMMENT='a comment')
    write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT',
                PRODUCT_CREATION_TIME='2016-10-31T00:00:00',
                extra=('SPICE_FILE_NAME = {"b.bc", "a.bsp"}',))
    obs = NewHorizonsLORRI.from_file(tmp_path / f'lor_test.{suffix}',
                                     astrometry=astrometry)
    value = obs.subfields.get(subfield, {})
    assert json.loads(json.dumps(value)) == value


@pytest.mark.parametrize('subfield', ['dict', 'label_dict'])
def test_end_of_the_label_is_left_out(tmp_path: pathlib.Path, subfield: str) -> None:
    assert 'END' not in from_label(tmp_path).subfields[subfield]


def test_headers_keep_the_commentary_cards(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    with fits.open(path, mode='update') as hdulist:
        hdulist[0].header['COMMENT'] = 'one'
        hdulist[0].header['COMMENT'] = 'two'
    obs = NewHorizonsLORRI.from_file(path)
    assert obs.headers[0]['COMMENT'] == ['one', 'two']


def test_undefined_header_card_is_left_out(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    with fits.open(path, mode='update') as hdulist:
        hdulist[0].header['UNDEF'] = None
    obs = NewHorizonsLORRI.from_file(path)
    assert 'UNDEF' not in obs.header_dict


def test_dict_of_a_label_holds_both_sources(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path)
    assert obs.dict['INSTRUMENT_ID'] == 'LORRI'
    assert obs.dict['MISSION'] == 'New Horizons'


def test_label_key_wins_over_the_same_header_key_in_dict(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, MISSION='"LABEL MISSION"')
    assert obs.dict['MISSION'] == 'LABEL MISSION'


def test_dict_of_a_label_matches_the_fits_file_where_they_agree(
        tmp_path: pathlib.Path) -> None:
    label = label_and_fits(tmp_path)
    from_fits = NewHorizonsLORRI.from_file(tmp_path / 'lor_test.fit').dict
    dict_ = NewHorizonsLORRI.from_file(label).dict
    assert {key: dict_[key] for key in from_fits} == from_fits


def test_dict_start_time_is_an_iso_string(tmp_path: pathlib.Path) -> None:
    assert from_label(tmp_path).dict['START_TIME'] == '2007-01-10T05:41:08.363'


def remote_label(tmp_path: pathlib.Path) -> Pds3Label:
    """A parsed synthetic label that claims to come from a remote volume."""

    label = Pds3Label(write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT'))
    label.filepath = FCPath('https://example.org/volume/lor_test.lbl')
    return label


def test_label_astrometry_on_a_remote_label_names_no_local_file(
        tmp_path: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(remote_label(tmp_path), astrometry=True)
    assert 'abspath' not in obs.subfields
    assert obs.image_url == 'https://example.org/volume/lor_test.fit'


def test_parsed_label_of_another_instrument(tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_label(tmp_path / 'mvic.lbl', 'MVIC.FIT',
                                  INSTRUMENT_ID='"MVIC"'))
    with pytest.raises(oops.host.WrongHostError, match='not a New Horizons LORRI PDS3'):
        NewHorizonsLORRI.from_file(label)


PATH_SUBFIELDS = {'filepath', 'basename', 'abspath', 'image_url'}


def label_without_a_file(tmp_path: pathlib.Path) -> Pds3Label:
    """A synthetic label parsed from its text, so that it has no file path."""

    return Pds3Label(write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT').read_text())


def test_astrometry_of_a_label_without_a_file(tmp_path: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(label_without_a_file(tmp_path), astrometry=True)
    assert PATH_SUBFIELDS & set(obs.subfields) == set()


def test_full_read_of_a_label_without_a_file(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match='no file path is known for this Pds3Label, so its FITS '
                             'file cannot be found'):
        NewHorizonsLORRI.from_file(label_without_a_file(tmp_path))


def hdulist_without_a_file(tmp_path: pathlib.Path, **keywords: Any) -> fits.HDUList:
    """A synthetic FITS file, written with the given header keywords and copied into
    memory, so that it has no file name."""

    with fits.open(write_fits(tmp_path / 'lor_test.fit', **keywords)) as hdulist:
        return fits.HDUList([hdu.copy() for hdu in hdulist])


def test_message_for_an_hdulist_without_a_file(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.HostFormatError,
                       match=r'^no EXPOSURE_DURATION for LORRI image without a file '
                             r'path$'):
        NewHorizonsLORRI.from_file(hdulist_without_a_file(tmp_path, EXPTIME=None))


def test_hdulist_of_another_instrument_without_a_file(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.host.WrongHostError,
                       match=r'^not a New Horizons LORRI FITS file: HDUList without a '
                             r'file path$'):
        NewHorizonsLORRI.from_file(hdulist_without_a_file(tmp_path, INSTRU='mvi'))


def test_label_of_another_instrument_without_a_file(tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_label(tmp_path / 'mvic.lbl', 'MVIC.FIT',
                                  INSTRUMENT_ID='"MVIC"').read_text())
    with pytest.raises(oops.host.WrongHostError,
                       match=r'^not a New Horizons LORRI PDS3 label: Pds3Label without '
                             r'a file path$'):
        NewHorizonsLORRI.from_file(label)


def test_hdulist_without_a_file(tmp_path: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(hdulist_without_a_file(tmp_path))
    assert PATH_SUBFIELDS & set(obs.subfields) == set()


# The clock count, not a file name, makes the IDs of the header's frame and path unique
def test_hdulist_without_a_file_names_its_fits_frame_by_its_clock_count(
        tmp_path: pathlib.Path) -> None:
    obs = NewHorizonsLORRI.from_file(hdulist_without_a_file(tmp_path), pointing='fits')
    assert obs.frame.frame_id == '.NH_FRAME_' + FITS_HEADER['SPCSCLK0']


def test_fits_file_without_image_lines(tmp_path: pathlib.Path) -> None:
    path = tmp_path / 'lor_test.fit'
    primary = fits.PrimaryHDU()
    for (key, value) in FITS_HEADER.items():
        primary.header[key] = value
    primary.writeto(path)
    with pytest.raises(oops.host.HostFormatError, match='no lines of the image'):
        NewHorizonsLORRI.from_file(path, astrometry=True)


def test_parsed_object_of_another_kind() -> None:
    with pytest.raises(oops.host.WrongHostError, match='label or FITS file: 42'):
        NewHorizonsLORRI.from_file(42)

##########################################################################################
# Merging the label and the FITS header
##########################################################################################


def test_every_item_has_one_winning_source() -> None:
    assert set(_WINNING_SOURCE) == set(_FITS_ITEM_KEYWORDS)
    assert set(_WINNING_SOURCE.values()) <= {'label', 'header'}


# The header's start wins whether or not the label's is the same time rounded
@pytest.mark.parametrize('fits_start', [LABEL_START + 0.0004, LABEL_START - 0.0004,
                                        LABEL_START + 0.002, LABEL_START - 1.])
def test_start_time_rule(tmp_path: pathlib.Path, fits_start: float) -> None:
    spcscet = fits_start + FITS_HEADER['EXPTIME'] / 2.
    obs = from_label(tmp_path, fits_keywords={'SPCSCET': spcscet})
    assert obs.time[0] == fits_start


def test_dict_start_time_follows_the_header(tmp_path: pathlib.Path) -> None:
    spcscet = LABEL_START + 0.002 + FITS_HEADER['EXPTIME'] / 2.
    obs = from_label(tmp_path, fits_keywords={'SPCSCET': spcscet})
    assert obs.dict['START_TIME'] == '2007-01-10T05:41:08.365'


@pytest.mark.parametrize('exptime', [4.9996, 4.9, 6.0])
def test_exposure_duration_rule(tmp_path: pathlib.Path, exptime: float) -> None:
    obs = from_label(tmp_path, fits_keywords={'EXPTIME': exptime})
    assert obs.texp == exptime


@pytest.mark.parametrize(('label_target', 'fits_target', 'target'), [
    ('J17 CALLIRRHOE', 'CALLIRRHOE', 'CALLIRRHOE'),
    ('JUPITER', 'FJUPITER', 'JUPITER'),
    ('J RINGS', 'ANSA', 'J RINGS'),
    ('KERBEROS', 'RADEC=270.66,-14.63=P4', 'KERBEROS'),
    ('PLUTO', 'JUPITER', 'PLUTO'),
])
def test_target_name_rule(tmp_path: pathlib.Path, label_target: str, fits_target: str,
                          target: str) -> None:
    obs = from_label(tmp_path, TARGET_NAME=f'"{label_target}"',
                     fits_keywords={'SPCCBTNM': fits_target})
    assert obs.label_target == label_target
    assert obs.dict['TARGET_NAME'] == label_target
    assert obs.target == target


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
    keywords: dict[str, Any] = {keyword: label_value}
    obs = from_label(tmp_path, **keywords)
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


def test_null_in_the_label_does_not_hide_the_header_value(
        tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, MISSION='"N/A"')
    assert obs.dict['MISSION'] == 'New Horizons'


def test_undefined_header_card_counts_as_missing(tmp_path: pathlib.Path) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    with fits.open(path, mode='update') as hdulist:
        hdulist[0].header['EXPTIME'] = None
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label)
    assert obs.texp == 5.0
    assert obs.time[0] == pytest.approx(LABEL_START, abs=1.e-6)


@pytest.mark.parametrize('dt', [0., 1.])
def test_fits_geometry_uses_the_header_epoch(tmp_path: pathlib.Path,
                                             sun_path: oops.Path, dt: float) -> None:
    """The header's position and velocity are those at SPCSCET, not at the observation's
    midtime; a velocity of 1000 km/s makes an epoch even 1 microsecond off show as 1 m.

    Without EXPTIME the header lacks a start time, so the label's, 2 s later, sets a
    midtime that differs from SPCSCET.
    """

    keywords = {'SPCSSCVX': 1000., 'SPCSSCVY': 0., 'SPCSSCVZ': 0., 'EXPTIME': None}
    obs = from_label(tmp_path, START_TIME='2007-01-10T05:41:10.363', geom='fits',
                     fits_keywords=keywords)
    assert obs.midtime == pytest.approx(FITS_HEADER['SPCSCET'] + 2., abs=1.e-3)
    event = obs.path.wrt('SUN').event_at_time(FITS_HEADER['SPCSCET'] + dt)
    assert event.pos.vals == pytest.approx([-1.e8 - 1000. * dt, -2.e8, -3.e8],
                                           rel=0., abs=1.e-4)


def test_null_values_stay_in_the_raw_dicts(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, TARGET_NAME='"N/A"', fits_keywords={'SPCCBTNM': '---'})
    assert obs.label_dict['TARGET_NAME'] == 'N/A'
    assert obs.header_dict['SPCCBTNM'] == '---'
    assert obs.target == 'NONE'


@pytest.mark.parametrize('null', ['UNKNOWN', '', '  ', 'unk', '---'])
def test_header_null_values_are_left_out_of_dict(tmp_path: pathlib.Path,
                                                 null: str) -> None:
    obs = from_label(tmp_path, fits_keywords={'SPCTCB': null})
    assert 'SPCTCB' not in obs.dict


@pytest.mark.parametrize('null', ['"UNKNOWN"', '""', '"N/A"', '"NOT APPLICABLE"'])
def test_label_null_values_are_left_out_of_dict(tmp_path: pathlib.Path,
                                                null: str) -> None:
    obs = from_label(tmp_path, extra=(f'NEWHORIZONS:OBSERVATION_DESC = {null}',))
    assert 'NEWHORIZONS:OBSERVATION_DESC' not in obs.dict


def test_list_of_nulls_is_left_out_of_dict(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path,
                     extra=('SC_TARGET_POSITION_VECTOR = ("N/A", "N/A", "N/A")',))
    assert 'SC_TARGET_POSITION_VECTOR' not in obs.dict


def test_null_inside_a_vector_is_none_in_dict(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, extra=('SC_TARGET_POSITION_VECTOR = (1.0, "N/A", 3.0)',))
    assert obs.dict['SC_TARGET_POSITION_VECTOR'] == [1.0, None, 3.0]


def test_null_inside_an_object_is_left_out_of_dict(tmp_path: pathlib.Path) -> None:
    obs = from_label(tmp_path, extra=('OBJECT = THING', '  NAME = "N/A"', '  SIZE = 3',
                                      'END_OBJECT = THING'))
    assert 'NAME' not in obs.dict['THING']
    assert obs.dict['THING']['SIZE'] == 3

##########################################################################################
# Values the Snapshot needs
##########################################################################################


# Without SPCSCET, the start would come from the clock count, SPCSCLK0
@pytest.mark.parametrize(('missing', 'message'), [
    (('EXPTIME',), 'no EXPOSURE_DURATION for LORRI image .*lor_test.fit'),
    (('SPCSCET', 'SPCSCLK0'), 'no START_TIME for LORRI image .*lor_test.fit'),
])
def test_fits_file_without_a_needed_item(tmp_path: pathlib.Path, missing: tuple[str, ...],
                                         message: str) -> None:
    keywords: dict[str, Any] = dict.fromkeys(missing)
    path = write_fits(tmp_path / 'lor_test.fit', **keywords)
    with pytest.raises(oops.host.HostFormatError, match=message):
        NewHorizonsLORRI.from_file(path)


@pytest.mark.parametrize(('keywords', 'lines', 'message'), [
    ({'START_TIME': None, 'SPACECRAFT_CLOCK_START_COUNT': None}, 1024,
     'no START_TIME for LORRI image .*lor_test.fit'),
    ({'EXPOSURE_DURATION': None}, 1024, 'no EXPOSURE_DURATION for LORRI image'),
    ({}, None, 'no lines of the image'),
    ({}, '"N/A"', 'no lines of the image'),
])
def test_label_astrometry_without_a_needed_item(tmp_path: pathlib.Path,
                                                keywords: dict[str, Any],
                                                lines: int | str | None,
                                                message: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', lines=lines,
                        **keywords)
    with pytest.raises(oops.host.HostFormatError, match=message):
        NewHorizonsLORRI.from_file(label, astrometry=True)

##########################################################################################
# Parsing the label
##########################################################################################


# "fast" rejects a quoted value that starts at the beginning of the next line; "strict"
# rejects a set of unquoted words with a blank inside, which "loose" accepts
NOTE = ('NOTE =', '"a note that starts', 'on the next line"')


@pytest.mark.parametrize(('extra', 'expected'), [
    (NOTE, {'NOTE': 'a note that starts\non the next line'}),
    (NOTE + ('WORDS = {a, b c}',), {'NOTE': 'a note that starts\non the next line',
                                    'WORDS': ['a', 'b', 'c']}),
])
def test_label_parser_falls_back(tmp_path: pathlib.Path, extra: tuple[str, ...],
                                 expected: dict[str, Any]) -> None:
    write_fits(tmp_path / 'lor_test.fit')
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', extra=extra)
    obs = NewHorizonsLORRI.from_file(label)
    assert {key: obs.label_dict[key] for key in expected} == expected


def test_label_with_another_suffix_is_read(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.txt', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True)
    assert obs.label_target == 'JUPITER'


def test_file_that_is_neither_label_nor_fits(tmp_path: pathlib.Path) -> None:
    path = tmp_path / 'lor_test.txt'
    path.write_text('not a label\n')
    with pytest.raises(oops.host.WrongHostError,
                       match=r'^not a New Horizons LORRI file: .*lor_test\.txt$'):
        NewHorizonsLORRI.from_file(path)


def test_label_parser_compound_method(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    obs = NewHorizonsLORRI.from_file(label, astrometry=True,
                                     pds3_options={'method': 'compound'})
    assert obs.label_target == 'JUPITER'


def test_label_parser_compound_method_is_tried_alone(tmp_path: pathlib.Path) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT', START_TIME='"x')
    with pytest.raises(oops.host.PdsHostError,
                       match=r'^cannot parse the PDS3 label .*lor_test\.lbl: [^;]*$'):
        NewHorizonsLORRI.from_file(label, astrometry=True,
                                   pds3_options={'method': 'compound'})


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
                       match=r'^cannot parse the PDS3 label .*lor_test\.lbl: '
                             r'method fast: .*; method strict: .*; method loose: '):
        NewHorizonsLORRI.from_file(label, astrometry=True)


@pytest.mark.parametrize(('options', 'message'), [
    ({'method': 'quick'}, "^invalid method: 'quick'$"),
    ({'bogus': 1}, "unexpected keyword argument 'bogus'$"),
])
def test_invalid_label_parser_options(tmp_path: pathlib.Path, options: dict[str, Any],
                                      message: str) -> None:
    label = write_label(tmp_path / 'lor_test.lbl', 'LOR_TEST.FIT')
    with pytest.raises(oops.OopsValueError, match=message):
        NewHorizonsLORRI.from_file(label, pds3_options=options)

##########################################################################################
# Missing pointing
##########################################################################################


def no_pointing(monkeypatch: pytest.MonkeyPatch, message: str,
                error: type[Exception] = RuntimeError) -> None:
    """Make the SPICE frame beneath the camera frame raise `error(message)` when it is
    evaluated. Only that frame lacks data: the camera frame, a Cmatrix on it, still gives
    its own fixed rotation."""

    def transform_at_time(self: Any, *args: Any, **kwargs: Any) -> Any:
        raise error(message)

    monkeypatch.setattr(FakeSpiceFrame, 'transform_at_time', transform_at_time)


# The synthetic image's midtime, SPCSCET, is 2007-01-10T05:41:10.863 UTC
@pytest.mark.parametrize('suffix', ['fit', 'lbl'])
def test_missing_pointing_raises_oops_spice_error(tmp_path: pathlib.Path,
                                                  monkeypatch: pytest.MonkeyPatch,
                                                  suffix: str) -> None:
    label_and_fits(tmp_path)
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT) -- sxform -- insufficient data')
    with pytest.raises(oops.OopsSpiceError,
                       match=r'^no New Horizons LORRI pointing at '
                             r'2007-01-10T05:41:10\.863 for .*lor_test\.fit: '
                             r'SPICE\(NOFRAMECONNECT\)$'):
        NewHorizonsLORRI.from_file(tmp_path / f'lor_test.{suffix}')


def test_missing_pointing_is_found_beneath_the_camera_frame(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The camera frame's own fixed rotation still works when only the SPICE frame lacks
    data, and from_file() raises all the same."""

    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(CKINSUFFDATA) -- no data')
    NewHorizonsLORRI.initialize()
    NewHorizonsLORRI._define_camera_frame()
    camera = oops.Frame.as_frame('NH_LORRI')
    assert isinstance(camera.transform_at_time(MIDTIME), oops.Transform)
    with pytest.raises(oops.OopsSpiceError, match=r'SPICE\(CKINSUFFDATA\)$'):
        NewHorizonsLORRI.from_file(path)


def test_missing_pointing_error_attributes(tmp_path: pathlib.Path,
                                           monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(CKINSUFFDATA) -- no data')
    with pytest.raises(oops.OopsSpiceError) as excinfo:
        NewHorizonsLORRI.from_file(path)
    assert excinfo.value.spice_code == 'SPICE(CKINSUFFDATA)'
    assert excinfo.value.filepath == FCPath(path)
    assert not isinstance(excinfo.value, oops.host.HostError)


# The classes of exception that cspyce raises
@pytest.mark.parametrize('error', [LookupError, KeyError, OSError, RuntimeError,
                                   ValueError])
def test_each_spice_exception_class_is_missing_pointing(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch,
        error: type[Exception]) -> None:
    """Every class of exception that cspyce raises becomes OopsSpiceError."""

    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT)', error)
    with pytest.raises(oops.OopsSpiceError, match=r'SPICE\(NOFRAMECONNECT\)$'):
        NewHorizonsLORRI.from_file(path)


def test_other_pointing_errors_pass_through(tmp_path: pathlib.Path,
                                            monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'something else', TypeError)
    with pytest.raises(TypeError, match='something else'):
        NewHorizonsLORRI.from_file(path)


def test_pointing_is_not_checked_with_a_cmatrix(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    frame = oops.frame.Cmatrix(oops.Matrix3([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
                               oops.Frame.J2000, frame_id='GIVEN_CMATRIX')
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT)')
    obs = NewHorizonsLORRI.from_file(path, frame=frame)
    assert obs.frame.frame_id == frame.frame_id


# An index row merged with its supplement, and the middle of its exposure
INDEX_ROW = {'VOLUME_ID': 'NHJULO_2001', 'FILE_SPECIFICATION_NAME': 'data/lor_test.lbl',
             'START_TIME': '2007-01-10T05:41:08.363', 'EXPOSURE_DURATION': 5.,
             'BINNING_MODE': '1x1'}
ROW_MIDTIME = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T05:41:08.363')) + 2.5


def row_without_pointing(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    """The Snapshot of an index row whose pointing the C kernels cannot supply."""

    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT) -- sxform -- insufficient data')
    (obs,) = NewHorizonsLORRI.from_index(tmp_path / 'NHJULO_2001_index.lbl',
                                         row_dicts=[dict(INDEX_ROW)])
    return obs


def test_from_index_keeps_a_row_without_pointing(tmp_path: pathlib.Path,
                                                 monkeypatch: pytest.MonkeyPatch) -> None:
    assert row_without_pointing(tmp_path, monkeypatch).basename == 'lor_test.fit'


# The row's start, 2007-01-10T05:41:08.363 UTC, is when its frame is evaluated here
def test_frame_of_a_row_without_pointing_raises_oops_spice_error(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    obs = row_without_pointing(tmp_path, monkeypatch)
    with pytest.raises(oops.OopsSpiceError,
                       match=r'^no New Horizons LORRI pointing at '
                             r'2007-01-10T05:41:08\.363 for .*NHJULO_2001/data/'
                             r'lor_test\.fit: SPICE\(NOFRAMECONNECT\)$'):
        obs.frame.transform_at_time(obs.time[0])


def spice_check_times(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Record the times at which the pointing check evaluates the SPICE frame beneath the
    camera frame; return the list that the times are appended to."""

    times: list[float] = []
    checking: list[bool] = []
    original_check = oops.Host._check_spice_frame
    original_transform = FakeSpiceFrame.transform_at_time

    def check_spice_frame(*args: Any, **kwargs: Any) -> Any:
        checking.append(True)
        try:
            return original_check(*args, **kwargs)
        finally:
            checking.pop()

    def transform_at_time(self: Any, time: Any, *, quick: Any = None) -> Any:
        if checking:
            times.append(float(oops.Scalar.as_scalar(time).vals))
        return original_transform(self, time, quick=quick)

    monkeypatch.setattr(oops.Host, '_check_spice_frame', staticmethod(check_spice_frame))
    monkeypatch.setattr(FakeSpiceFrame, 'transform_at_time', transform_at_time)
    return times


# What a time shift of the frame should mean for a camera frame that is a fixed rotation
# over a SPICE frame is an open question, SETI/rms-oops#251; these tests follow what oops
# does now
@pytest.mark.parametrize(('options', 'expected'), [
    ({'timeshift': (3., 'T')}, MIDTIME + 3.),
    ({'timeshift': (3., 'PT')}, MIDTIME + 3.),
    ({'timeshift': (3., 'F')}, MIDTIME),
    ({'timeshift': (3., 'TF')}, MIDTIME + 3.),
    ({'timeshift': (3., 't'), 'cadence': oops.cadence.SnapCadence(1000., 10.)}, 1008.),
])
def test_pointing_is_checked_at_the_shifted_midtime(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, options: dict[str, Any],
        expected: float) -> None:
    """A time shift of the cadence ("T") moves the time at which the C kernels are read.
    One of the frame ("F") shifts only the camera's fixed rotation, so the SPICE frame
    beneath it is read at the observation's own time."""

    times = spice_check_times(monkeypatch)
    NewHorizonsLORRI.from_file(write_fits(tmp_path / 'lor_test.fit'), **options)
    assert times == [pytest.approx(expected, abs=1.e-6)]


def pointing_until(monkeypatch: pytest.MonkeyPatch, last: float) -> None:
    """Make the SPICE frame beneath the camera frame lack data after the time `last`."""

    original = FakeSpiceFrame.transform_at_time

    def transform_at_time(self: Any, time: Any, *, quick: Any = None) -> Any:
        if float(oops.Scalar.as_scalar(time).vals) > last:
            raise RuntimeError('SPICE(CKINSUFFDATA) -- no data')
        return original(self, time, quick=quick)

    monkeypatch.setattr(FakeSpiceFrame, 'transform_at_time', transform_at_time)


def test_frame_time_shift_loads_where_the_pointing_ends(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With the pointing ending one second after the midtime, a three-second shift of the
    frame alone loads, because the observation reads the SPICE frame at its midtime."""

    path = write_fits(tmp_path / 'lor_test.fit')
    pointing_until(monkeypatch, MIDTIME + 1.)
    obs = NewHorizonsLORRI.from_file(path, timeshift=(3., 'F'))
    transform = oops.Frame.as_frame(obs.frame).wrt(oops.Frame.J2000).transform_at_time(
        obs.midtime)
    assert isinstance(transform, oops.Transform)


def test_cadence_time_shift_is_checked_where_the_pointing_ends(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With the pointing ending one second after the midtime, a three-second shift of the
    cadence moves the midtime past the end and so is refused at once."""

    path = write_fits(tmp_path / 'lor_test.fit')
    pointing_until(monkeypatch, MIDTIME + 1.)
    with pytest.raises(oops.OopsSpiceError, match=r'SPICE\(CKINSUFFDATA\)$'):
        NewHorizonsLORRI.from_file(path, timeshift=(3., 'T'))


def test_pointing_is_checked_at_the_cadence_override(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    times = spice_check_times(monkeypatch)
    cadence = oops.cadence.SnapCadence(1000., 10.)
    NewHorizonsLORRI.from_file(write_fits(tmp_path / 'lor_test.fit'), cadence=cadence)
    assert times == [1005.]


@pytest.mark.parametrize(('options', 'expected'), [
    ({'cadence': oops.cadence.SnapCadence(1000., 10.)}, 1005.),
    ({'timeshift': (3., 'T')}, ROW_MIDTIME + 3.),
    ({'timeshift': (3., 'f')}, ROW_MIDTIME),
])
def test_from_index_pointing_is_checked_at_the_overridden_midtime(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, options: dict[str, Any],
        expected: float) -> None:
    times = spice_check_times(monkeypatch)
    list(NewHorizonsLORRI.from_index(tmp_path / 'NHJULO_2001_index.lbl',
                                     row_dicts=[dict(INDEX_ROW)], **options))
    assert times == [pytest.approx(expected, abs=1.e-6)]


@pytest.mark.parametrize('options', [{'frame': 'J2000'}, {'pointing': 'fits'}])
def test_pointing_is_not_checked_without_the_c_kernels(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch,
        options: dict[str, Any]) -> None:
    path = write_fits(tmp_path / 'lor_test.fit')
    no_pointing(monkeypatch, 'SPICE(NOFRAMECONNECT)')
    obs = NewHorizonsLORRI.from_file(path, **options)
    assert obs.basename == 'lor_test.fit'

##########################################################################################
