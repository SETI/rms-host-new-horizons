##########################################################################################
# tests/oops/conftest.py
##########################################################################################
"""Shared fixtures that let the New Horizons host run without SPICE kernels or spicedb.

Every spicedb function the package calls is replaced by :class:`FakeSpiceDB`, and the
oops and cspyce calls that would read SPICE kernels are replaced by kernel-free
equivalents. Everything else, including oops itself, runs for real.
"""

import pathlib
from dataclasses import dataclass, field
from typing import Any

import cspyce
import numpy as np
import oops
import pytest
import spicedb
from astropy.io import fits

from host_new_horizons.oops.lorri import NewHorizonsLORRI

TEST_FILES = pathlib.Path(__file__).parents[2] / 'test_files'

# The NAIF IDs that the fake cspyce.bodn2c knows
BODY_IDS = {
    'JUPITER'   : 599,
    'CALLIRRHOE': 517,
    'PLUTO'     : 999,
    'ARROKOTH'  : 2486958,
}

# The primary header of a synthetic LORRI FITS file, before any overrides. Its values
# agree with those of write_label(), as the archive's do: the label's START_TIME is the
# start of the exposure, SPCSCET - EXPTIME/2, rounded to the millisecond.
FITS_HEADER: dict[str, Any] = {
    'MISSION' : 'New Horizons',
    'INSTRU'  : 'lor',
    'EXPTIME' : 5.0,
    'SPCSCET' : 221679736.0467595,
    'SPCUTCID': '2007-010T05:41:10.863',
    'SPCCBTNM': 'JUPITER',
    'SPCSCLK0': '0030713585:49000',
    'STARTMET': 30713585.98,
    'STOPMET' : 30713590.98,
    'SPCEMEN' : 275.0,
    'CRPIX1'  : 512.5,
    'CRPIX2'  : 512.5,
    'CRVAL1'  : 237.5,
    'CRVAL2'  : -33.8,
    'CD1_1'   : -2.8e-4,
    'CD1_2'   : 1.0e-6,
    'CD2_1'   : 1.0e-6,
    'CD2_2'   : 2.8e-4,
    'SPCSSCX' : 1.0e8,
    'SPCSSCY' : 2.0e8,
    'SPCSSCZ' : 3.0e8,
    'SPCSSCVX': 1.0,
    'SPCSSCVY': 2.0,
    'SPCSSCVZ': 3.0,
}


def write_fits(path: pathlib.Path, *, lines: int = 1024, samples: int = 1024,
               calibrated: bool = True, **keywords: Any) -> pathlib.Path:
    """Write a synthetic LORRI FITS file and return its path.

    The primary image counts up from zero. A calibrated file adds the error and quality
    extensions; otherwise, a histogram extension stands in for a raw product's. Keywords
    override :data:`FITS_HEADER`; a value of None deletes the keyword.
    """

    data = np.arange(lines * samples, dtype='>f4').reshape(lines, samples)
    primary = fits.PrimaryHDU(data)
    for (key, value) in (FITS_HEADER | keywords).items():
        if value is not None:
            primary.header[key] = value

    hdus: list[Any] = [primary]
    # The extension names are set in the header, because the name argument of an HDU
    # would write them in upper case, unlike the archive's
    if calibrated:
        extensions = [(data + 0.5, 'LORRI Error image'),
                      (np.ones((lines, samples), dtype='int16'),
                       'LORRI Quality flag image')]
    else:
        extensions = [(np.zeros(32, dtype='int32'), 'HISTOGRAM')]
    for (array, extname) in extensions:
        hdu = fits.ImageHDU(array)
        hdu.header['EXTNAME'] = extname
        hdus.append(hdu)

    fits.HDUList(hdus).writeto(path)
    return path


# The keywords of a synthetic detached label, before any overrides, as they are written
LABEL_KEYWORDS: dict[str, str] = {
    'INSTRUMENT_HOST_NAME'        : '"NEW HORIZONS"',
    'INSTRUMENT_ID'               : '"LORRI"',
    'TARGET_NAME'                 : '"JUPITER"',
    'START_TIME'                  : '2007-01-10T05:41:08.363',
    'STOP_TIME'                   : '2007-01-10T05:41:13.363',
    'SPACECRAFT_CLOCK_START_COUNT': '"0030713585:49000"',
    'SPACECRAFT_CLOCK_STOP_COUNT' : '"0030713590:49000"',
    'EXPOSURE_DURATION'           : '5.000 <s>',
}


def write_label(path: pathlib.Path, image: str | None, *, lines: int | str | None = 1024,
                extra: tuple[str, ...] = (), **keywords: str | None) -> pathlib.Path:
    """Write a detached PDS3 label pointing at `image` and return its path.

    The label gives :data:`LABEL_KEYWORDS`, which `keywords` override; a value is written
    as given, so a string needs its quotes, and None leaves the keyword out. `extra`
    lines follow them, and then an IMAGE object of the given number of lines (no LINES
    if `lines` is None). With `image` None, the label has no ``^IMAGE`` pointer.
    """

    text = ['PDS_VERSION_ID = PDS3',
            'RECORD_TYPE = FIXED_LENGTH',
            'RECORD_BYTES = 2880']
    if image is not None:
        text.append(f'^IMAGE = ("{image}", 1)')
    text += [f'{key} = {value}' for (key, value) in (LABEL_KEYWORDS | keywords).items()
             if value is not None]
    text += list(extra)
    text += ['OBJECT = IMAGE']
    if lines is not None:
        text += [f'  LINES = {lines}']
    text += ['  LINE_SAMPLES = 1024',
             'END_OBJECT = IMAGE',
             'END', '']
    path.write_text('\r\n'.join(text), encoding='utf-8', newline='')
    return path


def ik_dict() -> dict[str, Any]:
    """The part of the LORRI instrument kernel that the host reads as a dictionary."""

    return {'INS': {'NH_LORRI_1X1': {'PIXEL_SAMPLES': 1024, 'PIXEL_LINES': 1024}}}


@dataclass
class FakeKernel:
    """A stand-in for a spicedb KernelInfo object."""

    filespec: str


@dataclass
class FakeSpiceDB:
    """Records every spicedb call and returns configurable results.

    Attributes:
        calls (list[tuple[str, tuple, dict]]): Each call as (name, args, kwargs).
        inst_kernels (list[FakeKernel]): Returned by `select_inst`.
        kernel_dict (dict): Returned by `as_dict`.
        used (list[str]): Returned by `used_basenames`.
    """

    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = field(default_factory=list)
    inst_kernels: list[FakeKernel] = field(
        default_factory=lambda: [FakeKernel('nh_lorri_v201.ti')])
    kernel_dict: dict[str, Any] = field(default_factory=ik_dict)
    used: list[str] = field(default_factory=lambda: ['nh_lorri_v201.ti', 'nh_v220.tf'])

    def _record(self, func: str, /, *args: Any, **kwargs: Any) -> None:
        self.calls.append((func, args, kwargs))

    def named(self, name: str) -> list[tuple[tuple[Any, ...], dict[str, Any]]]:
        """The (args, kwargs) of every recorded call to the named function."""

        return [(args, kwargs) for (n, args, kwargs) in self.calls if n == name]

    def open_db(self, *args: Any, **kwargs: Any) -> None:
        self._record('open_db', *args, **kwargs)

    def close_db(self, *args: Any, **kwargs: Any) -> None:
        self._record('close_db', *args, **kwargs)

    def _furnish(self, func: str, *args: Any, **kwargs: Any) -> list[str]:
        self._record(func, *args, **kwargs)
        return [func + '.kernel']

    def furnish_lsk(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_lsk', *args, **kwargs)

    def furnish_pck(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_pck', *args, **kwargs)

    def furnish_inst(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_inst', *args, **kwargs)

    def furnish_spk(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_spk', *args, **kwargs)

    def furnish_ck(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_ck', *args, **kwargs)

    def furnish_by_metafile(self, *args: Any, **kwargs: Any) -> list[str]:
        return self._furnish('furnish_by_metafile', *args, **kwargs)

    def select_inst(self, *args: Any, **kwargs: Any) -> list[FakeKernel]:
        self._record('select_inst', *args, **kwargs)
        return self.inst_kernels

    def furnish_kernels(self, *args: Any, **kwargs: Any) -> list[str]:
        self._record('furnish_kernels', *args, **kwargs)
        return []

    def as_dict(self, *args: Any, **kwargs: Any) -> dict[str, Any]:
        self._record('as_dict', *args, **kwargs)
        return self.kernel_dict

    def as_names(self, kernels: list[FakeKernel]) -> list[str]:
        self._record('as_names', kernels)
        return [kernel.filespec for kernel in kernels]

    def used_basenames(self, *args: Any, **kwargs: Any) -> list[str]:
        self._record('used_basenames', *args, **kwargs)
        return self.used


_SPICEDB_FUNCTIONS = ('open_db', 'close_db', 'furnish_lsk', 'furnish_pck', 'furnish_inst',
                      'furnish_spk', 'furnish_ck', 'furnish_by_metafile', 'select_inst',
                      'furnish_kernels', 'as_dict', 'as_names', 'used_basenames')


@dataclass
class FakeSpice:
    """Records the oops calls that would otherwise read SPICE kernels.

    Attributes:
        solar_system (list[tuple[tuple, dict]]): The (args, kwargs) of each call to
            `Body.define_solar_system`.
        spice_paths (list[tuple[str, str]]): The (path ID, origin) of each SpicePath.
        spice_frames (list[tuple[str, str]]): The (SPICE frame, frame ID) of each
            SpiceFrame.
    """

    solar_system: list[tuple[tuple[Any, ...], dict[str, Any]]] = field(
        default_factory=list)
    spice_paths: list[tuple[str, str]] = field(default_factory=list)
    spice_frames: list[tuple[str, str]] = field(default_factory=list)


@pytest.fixture(autouse=True)
def fake_spicedb(monkeypatch: pytest.MonkeyPatch) -> FakeSpiceDB:
    """Replace every spicedb function the package calls with a recording fake."""

    fake = FakeSpiceDB()
    for name in _SPICEDB_FUNCTIONS:
        monkeypatch.setattr(spicedb, name, getattr(fake, name))
    return fake


@pytest.fixture(autouse=True)
def fake_spice(monkeypatch: pytest.MonkeyPatch) -> FakeSpice:
    """Keep SPICE kernels out of every test and isolate the oops registries.

    Kernel-reading calls are replaced: `cspyce.bodn2c` looks names up in
    :data:`BODY_IDS`, a SpicePath becomes a FixedPath at the origin, and a SpiceFrame
    becomes an identity Cmatrix on J2000. The Path and Frame registries are copied so that
    paths and frames registered by one test never reach another.
    """

    fake = FakeSpice()

    def define_solar_system(*args: Any, **kwargs: Any) -> None:
        fake.solar_system.append((args, kwargs))

    def spice_path(spice_id: str, origin: str) -> oops.Path:
        fake.spice_paths.append((spice_id, origin))
        return oops.path.FixedPath((0., 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                                   path_id=spice_id)

    def spice_frame(spice_frame: str, *, frame_id: str) -> oops.Frame:
        fake.spice_frames.append((spice_frame, frame_id))
        return oops.frame.Cmatrix(oops.Matrix3.IDENTITY, oops.Frame.J2000,
                                  frame_id=frame_id)

    def bodn2c(name: str) -> int:
        if name not in BODY_IDS:
            raise KeyError(f'SPICE(BODYNAMENOTFOUND) -- body name "{name}" not found')
        return BODY_IDS[name]

    monkeypatch.setattr(cspyce, 'bodn2c', bodn2c)
    monkeypatch.setattr(oops.Body, 'define_solar_system',
                        staticmethod(define_solar_system))
    monkeypatch.setattr(oops.path, 'SpicePath', spice_path)
    monkeypatch.setattr(oops.frame, 'SpiceFrame', spice_frame)

    monkeypatch.setattr(oops.Path, '_PATH_REGISTRY', dict(oops.Path._PATH_REGISTRY))
    monkeypatch.setattr(oops.Path, '_PATH_CACHE', dict(oops.Path._PATH_CACHE))
    monkeypatch.setattr(oops.Frame, '_FRAME_REGISTRY', dict(oops.Frame._FRAME_REGISTRY))
    monkeypatch.setattr(oops.Frame, '_FRAME_CACHE', dict(oops.Frame._FRAME_CACHE))
    return fake


@pytest.fixture(autouse=True)
def reset_host() -> Any:
    """Start and end every test with the New Horizons and LORRI class states reset."""

    NewHorizonsLORRI._reset()
    yield
    NewHorizonsLORRI._reset()


@pytest.fixture
def sun_path() -> oops.Path:
    """A fixed path registered as SUN, at the Solar System barycenter."""

    return oops.path.FixedPath((0., 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                               path_id='SUN')


@pytest.fixture
def jupiter_path(monkeypatch: pytest.MonkeyPatch) -> oops.Path:
    """A fixed path registered as JUPITER, with a Body of that name that a TrackerFrame
    can track; the Body registry is restored afterward."""

    monkeypatch.setattr(oops.Body, 'BODY_REGISTRY', dict(oops.Body.BODY_REGISTRY))
    path = oops.path.FixedPath((1.e9, 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                               path_id='JUPITER')
    oops.Body('JUPITER', path, oops.Frame.J2000)
    return path

##########################################################################################
