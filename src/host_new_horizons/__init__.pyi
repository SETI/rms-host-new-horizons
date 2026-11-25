##########################################################################################
# host_new_horizons/__init__.pyi
##########################################################################################
"""Type stub for :mod:`host_new_horizons`.

The source types its methods in their docstrings rather than in their signatures, so the
type information for public symbols is published here instead. Only package stubs exist,
so a name is annotated when it is imported from the package that exports it and not when
it is imported from the module that defines it. The stub describes the shape of the API
exactly: every public name, its parameters, which of them are keyword-only, and which have
defaults. Types are given where they are unambiguous and are `Any` elsewhere.

The detectors and support methods that a Host subclass overrides or calls, such as
`_register` and `_detect_in_fits`, are private and are not published here.
"""

import pathlib
from collections.abc import Iterable, Iterator
from typing import Any, overload

from filecache import FCPath
from pdsparser import Pds3Label
from oops import Body
from oops.cadence import Cadence
from oops.calibration import Calibration
from oops.fov import FOV
from oops.frame import Frame
from oops.host import Host
from oops.lightsource import LightSource
from oops.observation import Observation, Snapshot
from oops.path import Path

__all__ = ['NewHorizonsLORRI']

_Select = int | slice | tuple[int | slice, ...]

class NewHorizonsLORRI(Host):
    # The runtime names the first parameter `fileinfo`, where Host names it `filepath`,
    # and does not take **kwargs, where Host does, so that Python refuses a keyword not
    # listed with TypeError. A missing file raises FileNotFoundError, and a file that
    # cannot be read or interpreted raises oops.host.HostFormatError, or its subclass
    # PdsHostError for a label that cannot be parsed. select chooses among the one
    # observation of the file: None or 0 returns it, and a slice returns a tuple.
    @overload  # type: ignore[override]
    @staticmethod
    def from_file(
        fileinfo: str | pathlib.Path | FCPath, *,
        geom: str = 'spice', pointing: str = 'spice', fov_type: str = 'fast',
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        select: None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None
    ) -> Snapshot: ...
    @overload
    @staticmethod
    def from_file(
        fileinfo: str | pathlib.Path | FCPath, *,
        geom: str = 'spice', pointing: str = 'spice', fov_type: str = 'fast',
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        select: _Select,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None
    ) -> Snapshot | tuple[Snapshot, ...]: ...
    # Every row is one Snapshot of a Ring-Moon Systems Node metadata index, read with its
    # supplement; supplement=None without row_dicts raises oops.host.HostIndexError. The
    # call reads only the labels; the iterator reads the index and its supplement 100
    # rows at a time, pairs each index row with the supplement row of the same file, and
    # builds each row's Snapshot when it reaches the row. Host takes **kwargs here too,
    # and the runtime does not, so that Python refuses a keyword not listed.
    @overload  # type: ignore[override]
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: Iterable[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, fov_type: str = 'fast',
        select: None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None
    ) -> Iterator[Snapshot]: ...
    @overload
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: Iterable[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, fov_type: str = 'fast',
        select: _Select,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None
    ) -> Iterator[Snapshot | tuple[Snapshot, ...]]: ...
    @staticmethod
    def initialize(*, asof: str | None = None, meta: str | None = None) -> None: ...

##########################################################################################
