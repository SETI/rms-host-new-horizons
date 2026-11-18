##########################################################################################
# tests/oops/test_new_horizons.py
##########################################################################################
"""Tests for the New Horizons mission support, host_new_horizons.oops._new_horizons."""

import julian
import pytest
import spicedb
from conftest import FakeSpice, FakeSpiceDB

from host_new_horizons.oops._new_horizons import _NewHorizons

MISSION = ('2006-01-19', '2022-01-01')

##########################################################################################
# initialize() and reset()
##########################################################################################


def test_initialize_defines_the_solar_system(fake_spice: FakeSpice) -> None:
    _NewHorizons.initialize(planets=[9], asof='2020-01-01', irregulars=False)
    assert fake_spice.solar_system == [(MISSION, {'asof': '2020-01-01T00:00:00',
                                                  'planets': [9],
                                                  'irregulars': False})]


def test_initialize_defines_the_spacecraft_path(fake_spice: FakeSpice) -> None:
    _NewHorizons.initialize()
    assert fake_spice.spice_paths == [('NEW HORIZONS', 'SSB')]


def test_initialize_furnishes_the_database_kernels(fake_spicedb: FakeSpiceDB) -> None:
    _NewHorizons.initialize(asof='2020-01-01')
    asof = '2020-01-01T00:00:00'
    calls = [(name, args, kwargs) for (name, args, kwargs) in fake_spicedb.calls
             if name.startswith('furnish')]
    assert calls == [
        ('furnish_lsk', (), {'asof': asof}),
        ('furnish_pck', (_NewHorizons.PCK_BODIES,), {'asof': asof}),
        ('furnish_inst', (-98,), {'asof': asof}),
        ('furnish_spk', (-98,), {'time': MISSION, 'name': 'NH-SPK-PREDICTED%',
                                 'asof': asof}),
        ('furnish_spk', (-98,), {'time': MISSION, 'name': 'NH-SPK-RECONSTRUCTED%',
                                 'asof': asof}),
        ('furnish_ck', (-98,), {'time': MISSION, 'asof': asof}),
    ]


def test_initialize_records_the_kernel_names() -> None:
    _NewHorizons.initialize()
    assert _NewHorizons.names == ['furnish_lsk.kernel', 'furnish_pck.kernel',
                                  'furnish_inst.kernel', 'furnish_spk.kernel',
                                  'furnish_spk.kernel', 'furnish_ck.kernel']


def test_pck_bodies_are_the_jupiter_and_pluto_systems() -> None:
    assert _NewHorizons.PCK_BODIES == [5, 9, 501, 502, 503, 504, 514, 515, 516, 901, 902,
                                       903, 904, 905]


def test_initialize_with_a_metakernel(fake_spicedb: FakeSpiceDB) -> None:
    _NewHorizons.initialize(meta='nh.tm')
    calls = [name for (name, _, _) in fake_spicedb.calls if name.startswith('furnish')]
    assert calls == ['furnish_by_metafile']
    assert fake_spicedb.named('furnish_by_metafile') == [(('nh.tm',),
                                                          {'time': MISSION,
                                                           'asof': None})]


def test_initialize_closes_the_database(fake_spicedb: FakeSpiceDB) -> None:
    _NewHorizons.initialize()
    names = [name for (name, _, _) in fake_spicedb.calls]
    assert names[0] == 'open_db'
    assert names[-1] == 'close_db'


def test_initialize_closes_the_database_after_a_failure(fake_spicedb: FakeSpiceDB,
                                                        monkeypatch: pytest.MonkeyPatch
                                                        ) -> None:
    def fail(*args: object, **kwargs: object) -> list[str]:
        raise OSError('no database')

    monkeypatch.setattr(spicedb, 'furnish_lsk', fail)
    with pytest.raises(OSError, match='no database'):
        _NewHorizons.initialize()
    assert fake_spicedb.calls[-1][0] == 'close_db'


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    _NewHorizons.initialize()
    _NewHorizons.initialize()
    assert len(fake_spice.solar_system) == 1
    assert _NewHorizons.initialized


def test_reset_clears_state() -> None:
    _NewHorizons.initialize()
    _NewHorizons.reset()
    assert not _NewHorizons.initialized
    assert _NewHorizons.names == []

##########################################################################################
# spice_instrument_kernel()
##########################################################################################


def test_spice_instrument_kernel(fake_spicedb: FakeSpiceDB) -> None:
    (kernel, name) = _NewHorizons.spice_instrument_kernel('LORRI', asof='2020-01-01')
    assert kernel == fake_spicedb.kernel_dict
    assert name == 'nh_lorri_v201.ti'
    assert fake_spicedb.named('select_inst') == [
        ((-98,), {'inst': 'lorri', 'types': 'IK', 'asof': '2020-01-01T00:00:00'})]
    assert fake_spicedb.named('furnish_kernels') == [((fake_spicedb.inst_kernels,),
                                                      {'fast': True})]

##########################################################################################
# used_kernels()
##########################################################################################

TIME = (julian.tdb_from_iso('2007-01-10T05:41:08'),
        julian.tdb_from_iso('2007-01-10T05:41:13'))


@pytest.mark.parametrize(('target', 'bodies'), [
    ('CALLIRRHOE', [5, 517]),
    ('JUPITER', [5, 599]),
    ('ARROKOTH', [2486958]),
    ('HD 1234', []),
    (None, []),
])
def test_used_kernels_bodies(fake_spicedb: FakeSpiceDB, target: str | None,
                             bodies: list[int]) -> None:
    _ = _NewHorizons.used_kernels(TIME, 'lorri', target=target)
    assert fake_spicedb.named('used_basenames')[0][1]['bodies'] == bodies


def test_used_kernels_query(fake_spicedb: FakeSpiceDB) -> None:
    kernels = _NewHorizons.used_kernels(TIME, 'lorri')
    assert kernels == fake_spicedb.used
    (_, kwargs) = fake_spicedb.named('used_basenames')[0]
    assert kwargs['types'] is None
    assert kwargs['inst'] == 'lorri'
    assert kwargs['sc'] == -98


def test_used_kernels_converts_the_times_to_tai(fake_spicedb: FakeSpiceDB) -> None:
    _ = _NewHorizons.used_kernels(TIME, 'lorri')
    time = fake_spicedb.named('used_basenames')[0][1]['time']
    assert time == pytest.approx(tuple(julian.tai_from_tdb(t) for t in TIME))


def test_used_kernels_without_cks(fake_spicedb: FakeSpiceDB) -> None:
    _ = _NewHorizons.used_kernels(TIME, 'lorri', ck=False)
    types = fake_spicedb.named('used_basenames')[0][1]['types']
    assert types == [t for t in spicedb.KERNEL_TYPE_SORT_ORDER if t != 'CK']

##########################################################################################
