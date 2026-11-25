##########################################################################################
# host_new_horizons/oops/__init__.py
##########################################################################################
"""The New Horizons hosts for oops.

The host is :class:`~host_new_horizons.oops.lorri.NewHorizonsLORRI`, documented with its
module below. Importing this package registers it with :class:`~oops.host.Host`, so that
its ``from_file()`` and ``from_index()`` recognize New Horizons LORRI images and
indexes.
"""

# Importing each instrument module registers its host with oops.Host
from .lorri import NewHorizonsLORRI

__all__ = ['NewHorizonsLORRI']

##########################################################################################
