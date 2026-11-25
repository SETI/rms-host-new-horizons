##########################################################################################
# tests/oops/test_lorri_merge.py
##########################################################################################
"""Tests of the metadata of NewHorizonsLORRI.from_file(): the items of the detached label
and the FITS header, which source wins each, and `dict`, `label_dict` and
`header_dict`."""

import datetime as dt
import pathlib
from typing import Any

import oops
import pytest
from astropy.io.fits import Header, PrimaryHDU
from conftest import FITS_HEADER, write_fits, write_label

from host_new_horizons.oops.lorri import NewHorizonsLORRI


def read(directory: pathlib.Path, *, fits: dict[str, Any] | None = None,
         label: dict[str, Any] | None = None,
         extra: tuple[str, ...] = ()) -> oops.observation.Snapshot:
    """The observation of a synthetic label and its FITS file, with the given keywords."""

    write_fits(directory / 'lor_1.fit', **(fits or {}))
    obs = NewHorizonsLORRI.from_file(write_label(directory / 'lor_1.lbl', 'LOR_1.FIT',
                                                 extra=extra, **(label or {})))
    assert isinstance(obs, oops.observation.Snapshot)
    return obs


def read_fits(directory: pathlib.Path, **keywords: Any) -> oops.observation.Snapshot:
    """The observation of a synthetic FITS file read alone."""

    obs = NewHorizonsLORRI.from_file(write_fits(directory / 'lor_1.fit', **keywords))
    assert isinstance(obs, oops.observation.Snapshot)
    return obs


def read_label(directory: pathlib.Path, **keywords: Any) -> oops.observation.Snapshot:
    """The observation of a synthetic label read alone, for astrometry."""

    obs = NewHorizonsLORRI.from_file(write_label(directory / 'lor_1.lbl', 'LOR_1.FIT',
                                                 **keywords), astrometry=True)
    assert isinstance(obs, oops.observation.Snapshot)
    return obs

##########################################################################################
# The winning source of each item
##########################################################################################


def test_start_time_from_the_header(tmp_path: pathlib.Path) -> None:
    """The header's start time wins over the label's, which is rounded."""

    obs = read(tmp_path)
    assert obs.time[0] == 221679736.0467595 - 2.5


def test_start_time_from_a_label_alone(tmp_path: pathlib.Path) -> None:
    """A label read alone gives its own start time, in UTC."""

    obs = read_label(tmp_path)
    assert obs.dict['START_TIME'] == dt.datetime(2007, 1, 10, 5, 41, 8, 363000)


def test_exposure_from_the_header(tmp_path: pathlib.Path) -> None:
    """The header's exposure wins over the label's."""

    obs = read(tmp_path, label={'EXPOSURE_DURATION': '4.999 <s>'})
    assert obs.texp == 5.


def test_exposure_from_a_label_alone(tmp_path: pathlib.Path) -> None:
    """A label read alone gives its own exposure."""

    obs = read_label(tmp_path, EXPOSURE_DURATION='4.999 <s>')
    assert obs.texp == 4.999


def test_target_from_the_label(tmp_path: pathlib.Path) -> None:
    """The label's target name wins over the header's."""

    obs = read(tmp_path, fits={'SPCCBTNM': 'FJUPITER'})
    assert obs.dict['TARGET_NAME'] == 'JUPITER'


def test_target_from_the_header_when_the_label_has_none(tmp_path: pathlib.Path) -> None:
    """A null target name in the label falls to the header's."""

    obs = read(tmp_path, fits={'SPCCBTNM': 'PLUTO'}, label={'TARGET_NAME': '"N/A"'})
    assert obs.target == 'PLUTO'


def test_start_time_from_the_header_when_the_label_has_none(tmp_path: pathlib.Path
                                                            ) -> None:
    """A label without START_TIME takes the header's."""

    obs = read(tmp_path, label={'START_TIME': None})
    assert obs.time[0] == 221679736.0467595 - 2.5


def test_times_from_the_label_when_the_header_has_no_exposure(tmp_path: pathlib.Path
                                                              ) -> None:
    """A header without EXPTIME gives neither the exposure nor the start, and the
    label's are used."""

    obs = read(tmp_path, fits={'EXPTIME': None}, label={'EXPOSURE_DURATION': '4.999 <s>'})
    assert obs.texp == 4.999
    assert obs.dict['START_TIME'] == dt.datetime(2007, 1, 10, 5, 41, 8, 363000)


def test_null_label_lines(tmp_path: pathlib.Path) -> None:
    """A label read alone whose image lines are null gives no FOV, and raises
    HostFormatError."""

    label = write_label(tmp_path / 'lor_1.lbl', 'LOR_1.FIT', lines='"N/A"')
    with pytest.raises(oops.host.HostFormatError, match="KeyError: 'LINES'"):
        NewHorizonsLORRI.from_file(label, astrometry=True)


def test_header_without_lines(tmp_path: pathlib.Path) -> None:
    """A FITS header without NAXIS2 gives no FOV, and raises HostFormatError."""

    path = tmp_path / 'lor_1.fit'
    PrimaryHDU(header=Header(FITS_HEADER)).writeto(path)
    with pytest.raises(oops.host.HostFormatError, match="KeyError: 'LINES'"):
        NewHorizonsLORRI.from_file(path, astrometry=True)


def test_clock_counts_from_the_label(tmp_path: pathlib.Path) -> None:
    """The label's clock counts win over the header's."""

    obs = read(tmp_path, fits={'SPCSCLK0': '0030713585:48999'})
    assert obs.spacecraft_clock_start_count == '0030713585:49000'


def test_start_count_of_a_fits_file_alone(tmp_path: pathlib.Path) -> None:
    """A FITS file read alone gives the start count from SPCSCLK0."""

    obs = read_fits(tmp_path, SPCSCLK0='0030713585:12345')
    assert obs.spacecraft_clock_start_count == '0030713585:12345'


def test_stop_count_of_a_fits_file_alone(tmp_path: pathlib.Path) -> None:
    """A FITS file read alone gives the stop count from STOPMET, rounded to the nearest
    1/50000-second tick."""

    obs = read_fits(tmp_path, STOPMET=30713590.123456)
    assert obs.spacecraft_clock_stop_count == '0030713590:06173'


def test_no_clock_counts(tmp_path: pathlib.Path) -> None:
    """A clock count that no source gives is absent."""

    obs = read_fits(tmp_path, SPCSCLK0=None, STOPMET=None)
    assert {'spacecraft_clock_start_count', 'spacecraft_clock_stop_count'} & set(
        obs.subfields) == set()

##########################################################################################
# Missing values
##########################################################################################


def test_fits_fill_value_is_missing(tmp_path: pathlib.Path) -> None:
    """The FITS fill value -1e32 is left out of `dict`."""

    obs = read_fits(tmp_path, SPCTSORN=-1.e32)
    assert 'SPCTSORN' not in obs.dict


def test_fits_fill_value_stays_in_header_dict(tmp_path: pathlib.Path) -> None:
    """The FITS fill value stays in `header_dict`."""

    obs = read_fits(tmp_path, SPCTSORN=-1.e32)
    assert obs.header_dict['SPCTSORN'] == -1.e32


def test_null_header_text_is_missing(tmp_path: pathlib.Path) -> None:
    """A null text in the header, such as "---", is left out of `dict`."""

    obs = read_fits(tmp_path, SPCTCB='---')
    assert 'SPCTCB' not in obs.dict


def test_null_header_target_gives_no_target(tmp_path: pathlib.Path) -> None:
    """A FITS file read alone whose target is "---" has no target."""

    obs = read_fits(tmp_path, SPCCBTNM='---')
    assert obs.target == 'NONE'


def test_null_label_value_is_missing(tmp_path: pathlib.Path) -> None:
    """A null label value is left out of `dict`."""

    obs = read(tmp_path, extra=('OBSERVATION_DESC = "N/A"',))
    assert 'OBSERVATION_DESC' not in obs.dict


def test_null_label_value_stays_in_label_dict(tmp_path: pathlib.Path) -> None:
    """A null label value stays in `label_dict`."""

    obs = read(tmp_path, extra=('OBSERVATION_DESC = "N/A"',))
    assert obs.label_dict['OBSERVATION_DESC'] == 'N/A'


def test_list_of_null_values_is_missing(tmp_path: pathlib.Path) -> None:
    """A label list whose elements are all null is left out of `dict`."""

    obs = read(tmp_path, extra=('SC_TARGET_POSITION_VECTOR = ("N/A", "N/A", "N/A")',))
    assert 'SC_TARGET_POSITION_VECTOR' not in obs.dict


def test_list_with_some_null_values_is_kept(tmp_path: pathlib.Path) -> None:
    """A label list with only some null elements is kept as the parser gives it."""

    obs = read(tmp_path, extra=('QUATERNION = (1.0, "N/A", 0.0, 0.0)',))
    assert obs.dict['QUATERNION'] == [1.0, 'N/A', 0.0, 0.0]


def test_null_value_in_a_label_object_is_missing(tmp_path: pathlib.Path) -> None:
    """A null value inside a label object is left out of that object in `dict`."""

    obs = read(tmp_path, extra=('OBJECT = HISTOGRAM', '  NOTE = "UNK"',
                                '  ITEMS = 32', 'END_OBJECT = HISTOGRAM'))
    assert 'NOTE' not in obs.dict['HISTOGRAM']

##########################################################################################
# dict, label_dict and header_dict
##########################################################################################


@pytest.mark.parametrize('reader', [read, read_fits, read_label])
def test_dict_start_time_is_a_datetime(tmp_path: pathlib.Path, reader: Any) -> None:
    """START_TIME in `dict` is a datetime whichever file is read."""

    obs = reader(tmp_path)
    assert isinstance(obs.dict['START_TIME'], dt.datetime)


def test_dict_start_time_of_the_header(tmp_path: pathlib.Path) -> None:
    """The header's START_TIME in `dict` is the UTC time of SPCSCET - EXPTIME/2, to the
    microsecond."""

    obs = read_fits(tmp_path)
    assert obs.dict['START_TIME'] == dt.datetime(2007, 1, 10, 5, 41, 8, 362574)


def test_converted_keywords_are_left_out_of_dict(tmp_path: pathlib.Path) -> None:
    """A FITS keyword converted to a PDS3 item appears only under the item's name."""

    obs = read_fits(tmp_path)
    assert {'SPCSCET', 'EXPTIME', 'SPCCBTNM', 'SPCSCLK0', 'STOPMET', 'NAXIS2'} & set(
        obs.dict) == set()


def test_converted_items_in_dict(tmp_path: pathlib.Path) -> None:
    """The converted items appear under their PDS3 names."""

    obs = read_fits(tmp_path)
    assert {'START_TIME', 'EXPOSURE_DURATION', 'TARGET_NAME',
            'SPACECRAFT_CLOCK_START_COUNT', 'SPACECRAFT_CLOCK_STOP_COUNT'} <= set(
        obs.dict)


def test_other_header_keywords_in_dict(tmp_path: pathlib.Path) -> None:
    """A header keyword that is not converted keeps its own name in `dict`."""

    obs = read_fits(tmp_path)
    assert obs.dict['SPCEMEN'] == 275.0


def test_naxis2_stays_in_the_header_dict(tmp_path: pathlib.Path) -> None:
    """The converted NAXIS2 is in `header_dict` under its own name."""

    obs = read_fits(tmp_path)
    assert obs.header_dict['NAXIS2'] == 1024


def test_integer_exposure_time_is_a_float(tmp_path: pathlib.Path) -> None:
    """An integer EXPTIME gives a float EXPOSURE_DURATION in `dict`."""

    obs = read_fits(tmp_path, EXPTIME=5)
    assert isinstance(obs.dict['EXPOSURE_DURATION'], float)


def test_commentary_cards_are_left_out_of_dict(tmp_path: pathlib.Path) -> None:
    """The COMMENT and HISTORY cards are left out of `dict`."""

    obs = read_fits(tmp_path, COMMENT='a comment', HISTORY='a history')
    assert {'COMMENT', 'HISTORY'} & set(obs.dict) == set()


def test_label_key_wins_over_header_keyword(tmp_path: pathlib.Path) -> None:
    """Where the label and the header both have a key, `dict` has the label's value."""

    obs = read(tmp_path, extra=('MISSION = "NEW HORIZONS MISSION"',))
    assert obs.dict['MISSION'] == 'NEW HORIZONS MISSION'


def test_label_dict_is_the_parsers(tmp_path: pathlib.Path) -> None:
    """`label_dict` keeps the keys that the parser adds."""

    obs = read(tmp_path)
    assert obs.label_dict['START_TIME_fmt'] == '2007-01-10T05:41:08.363'


def test_header_dict_is_the_primary_header(tmp_path: pathlib.Path) -> None:
    """`header_dict` keeps the keywords that `dict` converts."""

    obs = read_fits(tmp_path)
    assert obs.header_dict['SPCSCET'] == 221679736.0467595

##########################################################################################
# The target body
##########################################################################################


@pytest.mark.parametrize(('name', 'body'), [
    ('JUPITER', 'JUPITER'),
    ('J17 CALLIRRHOE', 'CALLIRRHOE'),
    ('136472 MAKEMAKE', 'MAKEMAKE'),
    ('50000 QUAOAR (2002 LM60)', 'QUAOAR'),
    ('ASTEROID 486958 (2014 MU69)', 'ARROKOTH'),
    ('ASTEROID 15810 (1994 JR1)', 'ASTEROID 15810 (1994 JR1)'),
    ('MU69', 'ARROKOTH'),
    ('PLUTO-CHARON', 'PLUTO'),
    ('CALIBRATION', 'NONE'),
    (None, 'NONE'),
])
def test_target_body(name: str | None, body: str) -> None:
    """Each form of target name gives the name SPICE knows the body by."""

    assert NewHorizonsLORRI._target_body(name) == body


def test_calibration_frames_have_no_target(tmp_path: pathlib.Path) -> None:
    """The label target CALIBRATION, with the header's "---", gives no target."""

    obs = read(tmp_path, fits={'SPCCBTNM': '---'}, label={'TARGET_NAME': '"CALIBRATION"'})
    assert obs.target == 'NONE'

##########################################################################################
