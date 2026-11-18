##########################################################################################
# host_new_horizons/oops/__init__.py
##########################################################################################

# Importing each instrument module registers its host with oops.Host
from .lorri import NewHorizonsLORRI

__all__ = ['NewHorizonsLORRI']

##########################################################################################
