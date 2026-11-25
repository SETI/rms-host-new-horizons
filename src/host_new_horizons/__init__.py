##########################################################################################
# host_new_horizons/__init__.py
##########################################################################################
"""The New Horizons mission hosts for the OOPS observation geometry toolkit.

The package exports :class:`~host_new_horizons.oops.lorri.NewHorizonsLORRI`. Importing
``host_new_horizons`` does not import oops; the first access to the class imports
:mod:`host_new_horizons.oops`, which registers the host with :class:`~oops.host.Host`.
"""

import importlib

__all__ = ['NewHorizonsLORRI']

# The oops hosts are imported on first access, so that importing host_new_horizons does
# not import oops. Importing host_new_horizons.oops registers every host with oops.Host.
_LAZY_NAMES = {
    'NewHorizonsLORRI': 'host_new_horizons.oops',
}


def __getattr__(name):
    """Import an oops host class on first access (PEP 562)."""

    if name in _LAZY_NAMES:
        return getattr(importlib.import_module(_LAZY_NAMES[name]), name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

##########################################################################################
