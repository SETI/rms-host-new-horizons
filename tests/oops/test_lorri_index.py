##########################################################################################
# tests/oops/test_lorri_index.py
##########################################################################################
"""Tests of NewHorizonsLORRI.from_index(), on the cut-down metadata index NHJULO_2001 and
its supplement, and on rows given directly."""

import pathlib
import shutil
from typing import Any

import julian
import oops
import pytest
from conftest import TEST_FILES

from host_new_horizons.oops.lorri import NewHorizonsLORRI

INDEX = TEST_FILES / 'NHJULO_2001_index_edited.lbl'
SUPPLEMENT_NAME = 'NHJULO_2001_supplemental_index_edited'

# A row of an index merged with its supplement, as a caller might give it
ROW = {
    'VOLUME_ID'                   : 'NHJULO_2001',
    'FILE_SPECIFICATION_NAME'     : 'data/20070110_003070/lor_0030702839_0x630_sci.lbl',
    'START_TIME'                  : '2007-01-10T02:42:01.360',
    'EXPOSURE_DURATION'           : 0.5,
    'BINNING_MODE'                : '4x4',
    'TARGET_NAME'                 : 'JUPITER',
    'SPACECRAFT_CLOCK_START_COUNT': '1/0030702838:48850',
    'SPACECRAFT_CLOCK_STOP_COUNT' : '1/0030702839:00000',
}


def index_files() -> list[str]:
    """The FILE_SPECIFICATION_NAME of every row of the index table, in order."""

    lines = INDEX.with_suffix('.tab').read_text().splitlines()
    return [line.split(',')[1].strip('" ') for line in lines]


def observations(**kwargs: Any) -> list[oops.observation.Snapshot]:
    """The observations of the cut-down index."""

    return list(NewHorizonsLORRI.from_index(INDEX, **kwargs))


def row_observation(**changes: Any) -> oops.observation.Snapshot:
    """The observation of :data:`ROW` with the given changes, given as row_dicts."""

    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[ROW | changes])
    assert isinstance(obs, oops.observation.Snapshot)
    return obs


def test_returns_an_iterator() -> None:
    """from_index returns an iterator, which a caller wraps in list() for len()."""

    rows = NewHorizonsLORRI.from_index(INDEX)
    assert iter(rows) is rows


def test_every_row_in_index_order() -> None:
    """Every row gives an observation, in the index's order."""

    names = [obs.basename for obs in observations()]
    assert names == [name.rpartition('/')[2] for name in index_files()]


def test_filepath_is_the_relative_label_path() -> None:
    """`filepath` joins VOLUME_ID and FILE_SPECIFICATION_NAME as the index gives them."""

    obs = observations()[0]
    assert obs.filepath == 'NHJULO_2001/data/20070110_003070/lor_0030702839_0x630_sci.lbl'


def test_row_has_no_file_urls() -> None:
    """A row names no file beyond its relative path."""

    obs = observations()[0]
    assert {'image_url', 'image_abspath', 'label_url', 'label_abspath'} & set(
        obs.subfields) == set()


def test_dict_is_the_row_with_its_masks() -> None:
    """`dict` is the merged row as the tables give it, with the "_mask" entries."""

    obs = observations()[0]
    assert 'BINNING_MODE_mask' in obs.dict


def test_index_value_wins_over_the_supplements() -> None:
    """Where both tables give a column, the row has the index's value."""

    obs = observations()[0]
    assert obs.dict['SPACECRAFT_CLOCK_START_COUNT'] == '1/0030702838:48850'


def test_a_declared_null_stays_in_dict_as_the_table_gives_it() -> None:
    """A column that holds a declared null text appears in `dict` with that text."""

    obs = row_observation(FILTER_NAME='N/A')
    assert obs.dict['FILTER_NAME'] == 'N/A'


def test_clock_counts_as_the_row_gives_them() -> None:
    """The clock-count subfields are the row's, with the index's partition."""

    obs = observations()[0]
    assert obs.spacecraft_clock_stop_count == '1/0030702839:00000'


def test_times_from_the_row() -> None:
    """The exposure starts at START_TIME and lasts EXPOSURE_DURATION."""

    obs = observations()[0]
    start = julian.tdb_from_tai(julian.tai_from_iso('2007-01-10T02:42:01.360'))
    assert obs.time == pytest.approx((start, start + 0.003))


def test_fov_from_a_lower_case_binning_mode() -> None:
    """BINNING_MODE is looked up without regard to case."""

    obs = row_observation(BINNING_MODE='4x4')
    assert tuple(obs.fov.uv_shape.vals) == (256, 256)


def test_fov_type() -> None:
    """The FOV type applies to every row."""

    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[ROW], fov_type='flat')
    assert isinstance(obs.fov.fov, oops.fov.FlatFOV)


def test_zero_exposure_row_loads() -> None:
    """A row with an exposure of zero gives an observation with an exposure of zero."""

    obs = row_observation(EXPOSURE_DURATION=0.)
    assert obs.texp == 0.


def test_no_calibration() -> None:
    """A row gets no calibration, because the tables do not give RSOLAR."""

    obs = row_observation()
    assert obs.calibrations == []


def test_target_from_the_row() -> None:
    """The target body comes from TARGET_NAME by the same rules as from_file."""

    obs = row_observation(TARGET_NAME='J1 IO')
    assert obs.target == 'IO'


def test_row_subfields() -> None:
    """A row has the SPICE and C kernel subfields of a file."""

    obs = row_observation()
    assert obs.ck_object_id == -98000


def test_row_pointing_from_the_c_kernels() -> None:
    """A row's frame is the camera frame NH_LORRI."""

    obs = row_observation()
    assert obs.frame == oops.Frame.as_wayframe('NH_LORRI')


def test_frame_override_on_every_row() -> None:
    """A frame override is used as given on every row."""

    cmatrix = oops.frame.Cmatrix(oops.Matrix3.IDENTITY)
    (obs,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[ROW], frame=cmatrix)
    assert obs.frame is cmatrix


def test_select_on_every_row() -> None:
    """select applies to each row's observation."""

    (chosen,) = NewHorizonsLORRI.from_index(INDEX, row_dicts=[ROW], select=slice(None))
    assert isinstance(chosen, tuple)


def test_filter() -> None:
    """filter keeps the rows whose column has the value."""

    rows = observations(filter={'BINNING_MODE': '4x4'})
    assert {obs.dict['BINNING_MODE'] for obs in rows} == {'4x4'}


def test_no_supplement_is_refused() -> None:
    """supplement=None raises HostIndexError, because the index needs its supplement."""

    with pytest.raises(oops.host.HostIndexError, match='needs its supplemental index'):
        NewHorizonsLORRI.from_index(INDEX, supplement=None)


def test_no_supplement_with_rows() -> None:
    """supplement=None is accepted with row_dicts, which are already merged."""

    rows = list(NewHorizonsLORRI.from_index(INDEX, supplement=None, row_dicts=[ROW]))
    assert len(rows) == 1


def test_supplement_in_another_directory(tmp_path: pathlib.Path) -> None:
    """An explicit supplement is read from the directory it names."""

    for suffix in ('.lbl', '.tab'):
        shutil.copy(TEST_FILES / (SUPPLEMENT_NAME + suffix), tmp_path)
    rows = observations(supplement=tmp_path / (SUPPLEMENT_NAME + '.lbl'))
    assert len(rows) == len(index_files())


def test_missing_explicit_supplement(tmp_path: pathlib.Path) -> None:
    """An explicit supplement that does not exist raises FileNotFoundError."""

    with pytest.raises(FileNotFoundError):
        NewHorizonsLORRI.from_index(INDEX, supplement=tmp_path / 'missing.lbl')


def test_supplement_out_of_order(tmp_path: pathlib.Path) -> None:
    """Each index row is paired with the supplement row of the same file, wherever it
    is."""

    for name in ('NHJULO_2001_index_edited', SUPPLEMENT_NAME):
        for suffix in ('.lbl', '.tab'):
            shutil.copy(TEST_FILES / (name + suffix), tmp_path)
    table = tmp_path / (SUPPLEMENT_NAME + '.tab')
    rows = table.read_bytes().split(b'\r\n')[:-1]
    table.write_bytes(b'\r\n'.join([*reversed(rows), b'']))

    shuffled = list(NewHorizonsLORRI.from_index(tmp_path / INDEX.name))
    assert [obs.texp for obs in shuffled] == [obs.texp for obs in observations()]


def test_invalid_fov_type() -> None:
    """An unknown FOV type raises OopsValueError."""

    with pytest.raises(oops.OopsValueError, match="invalid fov_type 'bad'"):
        NewHorizonsLORRI.from_index(INDEX, fov_type='bad')


def test_unknown_keyword() -> None:
    """A keyword that is not listed raises TypeError."""

    with pytest.raises(TypeError, match="unexpected keyword argument 'asof'"):
        NewHorizonsLORRI.from_index(INDEX, asof='2020-01-01')


def test_dispatcher_reads_the_index() -> None:
    """oops.Host.from_index hands a LORRI index to this host."""

    rows = list(oops.Host.from_index(INDEX))
    assert len(rows) == len(index_files())

##########################################################################################
