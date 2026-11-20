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
from collections.abc import Iterable, Iterator, Sequence
from typing import Any

from astropy.io.fits import HDUList
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

class NewHorizonsLORRI(Host):
    # The runtime names the first parameter `fileinfo`, where Host names it `filepath`.
    # A missing or unreadable file raises OSError (FileNotFoundError for a FITS file that
    # a label names but that does not exist). A file whose content cannot be parsed or
    # used raises oops.host.HostFormatError, or its subclass PdsHostError for a label
    # (also for a label that lacks an ^IMAGE pointer). geom or pointing other than "spice"
    # in an astrometry read of a PDS3 label raises oops.OopsValueError. An image whose
    # C-kernel pointing is missing raises oops.OopsSpiceError, whose message ends with the
    # SPICE error code, such as "SPICE(NOFRAMECONNECT)"; never when `frame` is given or
    # `pointing` is not "spice". A select other than None raises oops.OopsValueError, and
    # any keyword not listed raises oops.OopsTypeError. from_index checks the pointing of
    # each row only to give a row without it a frame that raises oops.OopsSpiceError when
    # it is evaluated.
    @staticmethod
    def from_file(  # type: ignore[override]
        fileinfo: str | pathlib.Path | FCPath | Pds3Label | HDUList, *,
        geom: str = 'spice', pointing: str = 'spice', fov_type: str = 'fast',
        return_all_planets: bool = False,
        pds3_options: dict[str, Any] | None = None, astrometry: bool = False,
        select: int | slice | str | tuple[int | slice | str, ...] | None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None, **kwargs: Any) -> Snapshot: ...
    # Every row is one Snapshot. supplement=None without row_dicts raises
    # oops.host.HostIndexError, a select other than None raises oops.OopsValueError, and
    # any keyword not listed raises oops.OopsTypeError. The call reads only the labels;
    # the iterator reads the index and its supplement a few rows at a time, pairs each
    # index row with the supplement row of the same file, and builds each row's Snapshot
    # when it reaches the row. An index in a volume's own INDEX directory raises
    # oops.host.HostIndexError at call time.
    @staticmethod
    def from_index(
        filepath: str | pathlib.Path | FCPath, *,
        supplement: str | pathlib.Path | FCPath | object | None = ...,
        label: Pds3Label | None = None,
        row_dicts: Iterable[dict[str, Any]] | None = None,
        pds3_options: dict[str, Any] | None = None,
        filter: dict[str, Any] | None = None, fov_type: str = 'fast',
        return_all_planets: bool = False,
        select: int | slice | str | tuple[int | slice | str, ...] | None = None,
        target: str | Body | None = None, lightsource: LightSource | None = None,
        path: str | Path | None = None, frame: str | Frame | None = None,
        cadence: Cadence | None = None, fov: FOV | None = None,
        calibrations: Calibration | list[Calibration] | None = None,
        timeshift: tuple[float, str] | None = None,
        navigation: tuple[float, ...] | None = None, tracker: str | None = None,
        parallel: Observation | None = None, **kwargs: Any) -> Iterator[Snapshot]: ...
    @staticmethod
    def initialize(*, asof: str | None = None, meta: str | None = None,
        planets: int | Sequence[int] | None = None,
        irregulars: bool = True) -> None: ...

##########################################################################################
