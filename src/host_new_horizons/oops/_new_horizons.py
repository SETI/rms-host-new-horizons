##########################################################################################
# host_new_horizons/oops/_new_horizons.py
##########################################################################################

import julian
import spicedb
import cspyce
import oops

from oops.body import Body


class _NewHorizons:
    """An instance-free class to hold New Horizons-specific parameters.

    Attributes:
        START_TIME (str): The start of the period whose kernels are furnished, as an ISO
            date; the launch, 2006 January 19.
        STOP_TIME (str): The end of the period whose kernels are furnished, as an ISO
            date.
        SPACECRAFT_ID (int): The NAIF ID of the New Horizons spacecraft.
        PATH_ID (str): The ID of the spacecraft's registered Path.
        PCK_BODIES (list[int]): The NAIF IDs of the bodies whose PC kernels are
            furnished: the Jupiter and Pluto systems.
        initialized (bool): True after :meth:`initialize` has been called.
        names (list[str]): The names of the kernels that :meth:`initialize` furnished.
    """

    START_TIME = '2006-01-19'
    STOP_TIME  = '2022-01-01'

    SPACECRAFT_ID = -98
    PATH_ID = 'NEW HORIZONS'

    PCK_BODIES = ([5, 9] + list(range(501, 505)) + list(range(514, 517))
                  + list(range(901, 906)))

    initialized = False
    names = []

    @staticmethod
    def initialize(*, asof=None, meta=None, planets=None, irregulars=True):
        """Furnish the SPICE kernels of the New Horizons mission.

        Defines the Solar System bodies and the Path of the spacecraft,
        :attr:`PATH_ID`. Unless `meta` names a metakernel, the kernels come from the
        spicedb database: the leap seconds, the PC kernels of :attr:`PCK_BODIES`, the
        instrument, frame and clock kernels of the spacecraft, its predicted and
        reconstructed SP kernels, and its C kernels, each for the whole period from
        :attr:`START_TIME` to :attr:`STOP_TIME`. After the first call, later calls to
        this function are ignored.

        Parameters:
            asof (str, optional): A date in ISO format, which every database query
                receives as the date-time that :meth:`_format_asof` returns. Each kernel
                is chosen from those released no later than the start (00:00) of that
                date. A query of which every kernel was released after the start of the
                date silently uses the kernel it would choose without a date; a set of
                SP or C kernels keeps only those released by the start of the date, so
                the later part of the mission can lack coverage. The C kernels and a
                metakernel named in the database never fall back, and ValueError is
                raised instead. The kernels of the Solar System, the leap seconds kernel
                among them, are chosen by this date also with `meta`. None to choose the
                newest kernels.
            meta (str, optional): The name of a metakernel in the database, or the path
                to a metakernel file, to furnish instead of the kernels selected from the
                database. A metakernel named in the database is the newest one released
                by the start of `asof`; a metakernel file is used whatever the date. The
                kernels that it lists are furnished whatever their release dates, after
                the Solar System's, so a leap seconds kernel among them takes the place
                of the one that `asof` chose.
            planets (int | collections.abc.Sequence[int], optional): The planet number,
                1 to 9, or a sequence of planet numbers, such as (5, 9) for Jupiter and
                Pluto, whose bodies :meth:`~oops.Body.define_solar_system` defines with
                their moons. None or 0 means all nine.
            irregulars (bool, optional): True to include the irregular satellites; False
                otherwise.

        Raises:
            RuntimeError: If Saturn is among the planets and `asof` is a date from
                2003-11-05 through 2004-04-26, because the frames kernel of Saturn's
                small moons chosen for those dates has a malformed line that SPICE
                rejects.
            ValueError: If `meta` is None and none of the spacecraft's C kernels was
                released by the start of `asof`, as for any date up to and including
                2015-03-12 outside that window; or if `meta` names a metakernel in the
                database and none of that name was released by the start of `asof`.
        """

        if _NewHorizons.initialized:
            return

        asof = _NewHorizons._format_asof(asof)
        time = (_NewHorizons.START_TIME, _NewHorizons.STOP_TIME)
        sc = _NewHorizons.SPACECRAFT_ID

        Body.define_solar_system(*time, asof=asof, planets=planets,
                                 irregulars=irregulars)

        spicedb.open_db()
        try:
            if meta is None:
                names  = spicedb.furnish_lsk(asof=asof)
                names += spicedb.furnish_pck(_NewHorizons.PCK_BODIES, asof=asof)
                names += spicedb.furnish_inst(sc, asof=asof)
                names += spicedb.furnish_spk(sc, time=time, name='NH-SPK-PREDICTED%',
                                             asof=asof)
                names += spicedb.furnish_spk(sc, time=time,
                                             name='NH-SPK-RECONSTRUCTED%', asof=asof)
                names += spicedb.furnish_ck(sc, time=time, asof=asof)
            else:
                names = spicedb.furnish_by_metafile(meta, time=time, asof=asof)
        finally:
            spicedb.close_db()

        _ = oops.path.SpicePath(_NewHorizons.PATH_ID, 'SSB')

        _NewHorizons.names = names
        _NewHorizons.initialized = True

    @staticmethod
    def reset():
        """Reset the internal parameters.

        Can be useful for debugging. The kernels stay furnished.
        """

        _NewHorizons.names = []
        _NewHorizons.initialized = False

    @staticmethod
    def spice_instrument_kernel(inst, asof=None):
        """A dictionary containing the Instrument Kernel information.

        Also furnishes it for use by the SPICE tools.

        Parameters:
            inst (str): The instrument name, e.g., "LORRI"; case is ignored.
            asof (str, optional): A date in ISO date or date-time format. The instrument
                kernel is chosen from those released no later than it, at 00:00 for a
                date without a time, or, if every one was released after it, the kernel
                chosen without a date is used, silently. None to choose the newest.

        Returns:
            tuple[dict, str]: The dictionary generated by :func:`textkernel.from_file`
            and the name of the kernel.
        """

        spicedb.open_db()
        try:
            kernel_info = spicedb.select_inst(_NewHorizons.SPACECRAFT_ID,
                                              inst=inst.lower(), types='IK',
                                              asof=_NewHorizons._format_asof(asof))
            spicedb.furnish_kernels(kernel_info, fast=True)
        finally:
            spicedb.close_db()

        return (spicedb.as_dict(kernel_info), spicedb.as_names(kernel_info)[0])

    @staticmethod
    def used_kernels(time, inst, *, target=None, ck=True, return_all_planets=False):
        """The kernels associated with a New Horizons observation.

        Parameters:
            time (tuple[float, float]): A (start, stop) tuple of times in seconds TDB.
            inst (str): A fragment of the instrument kernel's file name that identifies
                the instrument, such as "lorri"; case is ignored.
            target (str, optional): The name of the target body. If SPICE knows it, the
                kernels of the body and of its system barycenter are included along with
                those of the spacecraft; otherwise, only the spacecraft's.
            ck (bool, optional): True (default) to include CK (pointing) kernels; False to
                exclude them, e.g. for an observation whose pointing came from a custom
                frame rather than SPICE, where any furnished CK is unrelated to how the
                observation was actually pointed.
            return_all_planets (bool, optional): True to include the kernels of every
                planet and of Pluto, each with its system barycenter, besides those of
                the target.

        Returns:
            list[str]: The basenames of the furnished kernels that apply to the
            observation.
        """

        # A planet or a moon, with an ID from 100 to 999, orbits its system barycenter,
        # whose ID is its own divided by 100; any other body stands alone
        bodies = []
        if return_all_planets:
            for barycenter in range(1, 10):
                bodies += [barycenter, barycenter * 100 + 99]
        if target is not None:
            try:
                body = cspyce.bodn2c(target)
            except (KeyError, ValueError):
                pass
            else:
                if 100 <= body < 1000:
                    bodies.append(body // 100)
                bodies.append(body)

        types = None
        if not ck:
            types = [t for t in spicedb.KERNEL_TYPE_SORT_ORDER if t != 'CK']

        # spicedb compares times in seconds TAI
        time_tai = tuple(julian.tai_from_tdb(t) for t in time)
        return spicedb.used_basenames(types=types, time=time_tai, inst=inst,
                                      sc=_NewHorizons.SPACECRAFT_ID,
                                      bodies=list(dict.fromkeys(bodies)))

    @staticmethod
    def _format_asof(asof):
        """An "as of" date reformatted as the ISO date-time spicedb expects, or None."""

        if asof is None:
            return None
        (day, sec) = julian.day_sec_from_iso(asof)
        return julian.ymdhms_format_from_day_sec(day, sec)

##########################################################################################
