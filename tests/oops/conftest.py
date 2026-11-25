##########################################################################################
# tests/oops/conftest.py
##########################################################################################
"""Shared fixtures that let the New Horizons host run without SPICE kernels or spicedb.

Every spicedb function the package calls is replaced by :class:`FakeSpiceDB`, and the
oops calls that would read SPICE kernels are replaced by kernel-free equivalents.
Everything else, including oops itself, runs for real.
"""

import pathlib
from dataclasses import dataclass, field
from typing import Any

import filecache
import numpy as np
import oops
import pytest
import spicedb
from astropy.io import fits
from filecache import FileCache

from host_new_horizons.oops.lorri import NewHorizonsLORRI

TEST_FILES = pathlib.Path(__file__).parents[2] / 'test_files'

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
    'RSOLAR'  : 5.0e6,
}


def write_fits(path: pathlib.Path, *, lines: int = 1024, samples: int = 1024,
               calibrated: bool = True, data: np.ndarray | None = None,
               **keywords: Any) -> pathlib.Path:
    """Write a synthetic LORRI FITS file and return its path.

    The primary image counts up from zero, unless `data` gives it. A calibrated file adds
    the error and quality extensions; otherwise, a histogram extension stands in for a
    raw product's. Keywords override :data:`FITS_HEADER`; a value of None deletes the
    keyword.
    """

    if data is None:
        data = np.arange(lines * samples, dtype='>f4').reshape(lines, samples)
    primary = fits.PrimaryHDU(data)
    for (key, value) in (FITS_HEADER | keywords).items():
        if value is not None:
            primary.header[key] = value

    hdus: list[Any] = [primary]
    # The extension names are set in the header, because the name argument of an HDU
    # would write them in upper case, unlike the archive's
    extensions: list[tuple[np.ndarray, str]]
    if calibrated:
        extensions = [(np.full((lines, samples), 0.5, dtype='>f4'), 'LORRI Error image'),
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


def write_label(path: pathlib.Path, image: str, *, lines: int | str = 1024,
                samples: int = 1024, extra: tuple[str, ...] = (),
                **keywords: str | None) -> pathlib.Path:
    """Write a detached PDS3 label pointing at `image` and return its path.

    The label gives :data:`LABEL_KEYWORDS`, which `keywords` override; a value is written
    as given, so a string needs its quotes, and None leaves the keyword out. `extra`
    lines follow them, and then an IMAGE object of the given lines and samples.
    """

    text = ['PDS_VERSION_ID = PDS3',
            'RECORD_TYPE = FIXED_LENGTH',
            'RECORD_BYTES = 2880',
            f'^IMAGE = ("{image}", 1)']
    text += [f'{key} = {value}' for (key, value) in (LABEL_KEYWORDS | keywords).items()
             if value is not None]
    text += list(extra)
    text += ['OBJECT = IMAGE',
             f'  LINES = {lines}',
             f'  LINE_SAMPLES = {samples}',
             'END_OBJECT = IMAGE',
             'END', '']
    path.write_text('\r\n'.join(text), encoding='utf-8', newline='')
    return path


def write_product(directory: pathlib.Path, **keywords: Any) -> pathlib.Path:
    """Write a label and its FITS file in the archive's lower case and return the label.

    The FITS file is written by :func:`write_fits` with `keywords`; the label names it
    in upper case, as the archive's labels do.
    """

    write_fits(directory / 'lor_0030713591_0x633_sci.fit', **keywords)
    return write_label(directory / 'lor_0030713591_0x633_sci.lbl',
                       'LOR_0030713591_0X633_SCI.FIT')


@dataclass
class FakeSpiceDB:
    """Records every spicedb call.

    Attributes:
        calls (list[tuple[str, tuple, dict]]): Each call as (name, args, kwargs).
    """

    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = field(default_factory=list)

    def recorder(self, name: str) -> Any:
        """A function that records each call under `name` and returns an empty list."""

        def record(*args: Any, **kwargs: Any) -> list[str]:
            """Record the call and return an empty list."""

            self.calls.append((name, args, kwargs))
            return []

        return record


_SPICEDB_FUNCTIONS = ('open_db', 'close_db', 'furnish_lsk', 'furnish_pck', 'furnish_inst',
                      'furnish_spk', 'furnish_ck', 'furnish_by_metafile')


class FakeSpiceFrame(oops.Frame):
    """A stand-in for a SpiceFrame: the identity rotation from J2000."""

    def __init__(self, frame_id: str) -> None:
        """Register the frame under `frame_id`.

        Parameters:
            frame_id (str): The ID of the frame.
        """

        self._reference = oops.Frame.J2000
        self._origin = None
        self._shape = ()
        self._register(frame_id)

    def transform_at_time(self, time: Any, *, quick: Any = None) -> oops.Transform:
        """The identity transform, at any time."""

        return oops.Transform(oops.Matrix3.IDENTITY, oops.Vector3.ZERO, self._wayframe,
                              self._reference)


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
        monkeypatch.setattr(spicedb, name, fake.recorder(name))
    return fake


@pytest.fixture(autouse=True)
def fake_spice(monkeypatch: pytest.MonkeyPatch) -> FakeSpice:
    """Keep SPICE kernels out of every test and isolate the oops registries.

    Kernel-reading calls are replaced: a SpicePath becomes a FixedPath at the origin, and
    a SpiceFrame becomes a :class:`FakeSpiceFrame`. The Path and Frame registries,
    caches, and waypoint and wayframe pools are copied so that paths and frames made by
    one test never reach another.
    """

    fake = FakeSpice()

    def define_solar_system(*args: Any, **kwargs: Any) -> None:
        """Record the call without defining any body."""

        fake.solar_system.append((args, kwargs))

    def spice_path(spice_id: str, origin: str) -> oops.Path:
        """Record the call and return a fixed path at the origin."""

        fake.spice_paths.append((spice_id, origin))
        return oops.path.FixedPath((0., 0., 0.), oops.Path.SSB, oops.Frame.J2000,
                                   path_id=spice_id)

    def spice_frame(spice_frame: str, *, frame_id: str) -> oops.Frame:
        """Record the call and return a :class:`FakeSpiceFrame`."""

        fake.spice_frames.append((spice_frame, frame_id))
        return FakeSpiceFrame(frame_id)

    monkeypatch.setattr(oops.Body, 'define_solar_system',
                        staticmethod(define_solar_system))
    monkeypatch.setattr(oops.path, 'SpicePath', spice_path)
    monkeypatch.setattr(oops.frame, 'SpiceFrame', spice_frame)

    # Each copy is made from items(), which works whether a pool is a plain dictionary
    # or a WeakValueDictionary, from which garbage collection can remove a member during
    # a copy made by dict()
    monkeypatch.setattr(oops.Path, '_PATH_REGISTRY',
                        dict(oops.Path._PATH_REGISTRY.items()))
    monkeypatch.setattr(oops.Path, '_PATH_CACHE', dict(oops.Path._PATH_CACHE.items()))
    monkeypatch.setattr(oops.Frame, '_FRAME_REGISTRY',
                        dict(oops.Frame._FRAME_REGISTRY.items()))
    monkeypatch.setattr(oops.Frame, '_FRAME_CACHE', dict(oops.Frame._FRAME_CACHE.items()))

    # A path or frame defined as an earlier one was shares its waypoint or wayframe, and
    # so its ID, so each class's pool is copied too; otherwise a fixed SUN path made by
    # one test would be the NEW HORIZONS path of the next, both being fixed at the origin
    for path_class in oops.Path._PATH_SUBCLASSES:
        if '_WAYPOINTS' in vars(path_class):
            monkeypatch.setattr(path_class, '_WAYPOINTS',
                                dict(path_class._WAYPOINTS.items()))
    for frame_class in oops.Frame._FRAME_SUBCLASSES:
        for pool in ('_WAYFRAMES', '_FOR_NAME', '_FRAME_LOOKUP'):
            if pool in vars(frame_class):
                monkeypatch.setattr(frame_class, pool,
                                    dict(getattr(frame_class, pool).items()))
    return fake


@pytest.fixture(autouse=True)
def file_cache(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> FileCache:
    """Give every test its own empty file cache, so that what an earlier test or another
    program put in the shared cache never decides whether a remote file is cached."""

    cache = FileCache(cache_root=tmp_path / 'filecache')
    monkeypatch.setattr(filecache.file_cache_path, '_DEFAULT_FILECACHE', cache)
    return cache


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
def jupiter_path(monkeypatch: pytest.MonkeyPatch, sun_path: oops.Path) -> oops.Path:
    """A fixed path registered as JUPITER, 5 AU from the Sun, with a Body of that name
    that a TrackerFrame can track; the Body registry is restored afterward."""

    monkeypatch.setattr(oops.Body, 'BODY_REGISTRY', dict(oops.Body.BODY_REGISTRY))
    path = oops.path.FixedPath((5. * 149597870.7, 0., 0.), oops.Path.SSB,
                               oops.Frame.J2000, path_id='JUPITER')
    oops.Body('JUPITER', path, oops.Frame.J2000)
    return path

##########################################################################################
