##########################################################################################
# host_new_horizons/__init__.py
##########################################################################################

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
