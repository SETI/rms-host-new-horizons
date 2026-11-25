##########################################################################################
# tests/oops/test_new_horizons.py
##########################################################################################
"""Tests for the New Horizons mission support, host_new_horizons.oops._new_horizons."""

import pytest
import spicedb
from conftest import FakeSpice, FakeSpiceDB

from host_new_horizons.oops._new_horizons import _NewHorizons

MISSION = ('2006-01-19', '2022-01-01')


def test_initialize_defines_the_solar_system(fake_spice: FakeSpice) -> None:
    """The bodies are defined for the whole mission, as of the given date."""

    _NewHorizons.initialize(asof='2020-01-01')
    assert fake_spice.solar_system == [(MISSION, {'asof': '2020-01-01'})]


def test_initialize_defines_the_spacecraft_path(fake_spice: FakeSpice) -> None:
    """The spacecraft path is registered relative to the barycenter."""

    _NewHorizons.initialize()
    assert fake_spice.spice_paths == [('NEW HORIZONS', 'SSB')]


def test_initialize_furnishes_the_database_kernels(fake_spicedb: FakeSpiceDB) -> None:
    """The database kernels are furnished for the whole mission, as of the given date."""

    _NewHorizons.initialize(asof='2020-01-01')
    asof = '2020-01-01'
    assert fake_spicedb.calls == [
        ('open_db', (), {}),
        ('furnish_lsk', (), {'asof': asof}),
        ('furnish_pck', (_NewHorizons.PCK_BODIES,), {'asof': asof}),
        ('furnish_inst', (-98,), {'asof': asof}),
        ('furnish_spk', (-98,), {'time': MISSION, 'name': 'NH-SPK-PREDICTED%',
                                 'asof': asof}),
        ('furnish_spk', (-98,), {'time': MISSION, 'name': 'NH-SPK-RECONSTRUCTED%',
                                 'asof': asof}),
        ('furnish_ck', (-98,), {'time': MISSION, 'asof': asof}),
        ('close_db', (), {}),
    ]


def test_pck_bodies_are_the_jupiter_and_pluto_systems() -> None:
    """The PC kernels cover the Jupiter and Pluto systems."""

    assert _NewHorizons.PCK_BODIES == [5, 9, 501, 502, 503, 504, 514, 515, 516, 901, 902,
                                       903, 904, 905]


def test_initialize_with_a_metakernel(fake_spicedb: FakeSpiceDB) -> None:
    """A metakernel replaces the database selection."""

    _NewHorizons.initialize(meta='nh.tm')
    assert fake_spicedb.calls == [('open_db', (), {}),
                                  ('furnish_by_metafile', ('nh.tm',),
                                   {'time': MISSION, 'asof': None}),
                                  ('close_db', (), {})]


def test_initialize_closes_the_database_after_a_failure(fake_spicedb: FakeSpiceDB,
                                                        monkeypatch: pytest.MonkeyPatch
                                                        ) -> None:
    """The database is closed even when furnishing fails."""

    def fail(*args: object, **kwargs: object) -> list[str]:
        """Raise the error of a database that cannot be read."""

        raise OSError('no database')

    monkeypatch.setattr(spicedb, 'furnish_lsk', fail)
    with pytest.raises(OSError, match='no database'):
        _NewHorizons.initialize()
    assert fake_spicedb.calls[-1][0] == 'close_db'


def test_initialize_runs_once(fake_spice: FakeSpice) -> None:
    """A second call does nothing."""

    _NewHorizons.initialize()
    _NewHorizons.initialize(asof='2020-01-01')
    assert len(fake_spice.solar_system) == 1


def test_reset_lets_initialize_run_again(fake_spice: FakeSpice) -> None:
    """After a reset, the next call furnishes again."""

    _NewHorizons.initialize()
    _NewHorizons._reset()
    _NewHorizons.initialize()
    assert len(fake_spice.solar_system) == 2

##########################################################################################
