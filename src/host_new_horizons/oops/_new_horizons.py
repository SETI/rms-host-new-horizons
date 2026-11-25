##########################################################################################
# host_new_horizons/oops/_new_horizons.py
##########################################################################################
"""The New Horizons mission support that the LORRI host builds on.

It furnishes the mission's SPICE kernels once, for the whole mission.
"""

import spicedb

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
    """

    START_TIME = '2006-01-19'
    STOP_TIME  = '2022-01-01'

    SPACECRAFT_ID = -98
    PATH_ID = 'NEW HORIZONS'

    PCK_BODIES = ([5, 9] + list(range(501, 505)) + list(range(514, 517))
                  + list(range(901, 906)))

    initialized = False

    @staticmethod
    def initialize(*, asof=None, meta=None):
        """Furnish the SPICE kernels of the New Horizons mission.

        Defines the Solar System bodies and the Path of the spacecraft,
        :attr:`PATH_ID`. Unless `meta` names a metakernel, the kernels come from the
        spicedb database: the leap seconds, the PC kernels of :attr:`PCK_BODIES`, the
        instrument, frame and clock kernels of the spacecraft, its predicted and
        reconstructed SP kernels, and its C kernels, each for the whole period from
        :attr:`START_TIME` to :attr:`STOP_TIME`. After the first call, later calls to
        this function are ignored.

        Parameters:
            asof (str, optional): Only use SPICE kernels that existed before this date;
                None to ignore.
            meta (str, optional): The name of a metakernel to furnish instead of the
                kernels selected from the database.
        """

        if _NewHorizons.initialized:
            return

        time = (_NewHorizons.START_TIME, _NewHorizons.STOP_TIME)
        sc = _NewHorizons.SPACECRAFT_ID

        Body.define_solar_system(*time, asof=asof)

        spicedb.open_db()
        try:
            if meta is None:
                spicedb.furnish_lsk(asof=asof)
                spicedb.furnish_pck(_NewHorizons.PCK_BODIES, asof=asof)
                spicedb.furnish_inst(sc, asof=asof)
                spicedb.furnish_spk(sc, time=time, name='NH-SPK-PREDICTED%', asof=asof)
                spicedb.furnish_spk(sc, time=time, name='NH-SPK-RECONSTRUCTED%',
                                    asof=asof)
                spicedb.furnish_ck(sc, time=time, asof=asof)
            else:
                spicedb.furnish_by_metafile(meta, time=time, asof=asof)
        finally:
            spicedb.close_db()

        _ = oops.path.SpicePath(_NewHorizons.PATH_ID, 'SSB')

        _NewHorizons.initialized = True

    @staticmethod
    def _reset():
        """Reset the internal parameters, so that the next :meth:`initialize` runs again.

        The kernels stay furnished.
        """

        _NewHorizons.initialized = False

##########################################################################################
