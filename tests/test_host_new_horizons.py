##########################################################################################
# tests/test_host_new_horizons.py
##########################################################################################
"""Tests of the top-level :mod:`host_new_horizons` package."""

import os
import pathlib
import subprocess
import sys

import pytest

import host_new_horizons


def test_unknown_attribute_raises() -> None:
    """An unexported name raises AttributeError naming the module and the attribute."""

    with pytest.raises(AttributeError,
                       match="module 'host_new_horizons' has no attribute 'NoSuchHost'"):
        _ = host_new_horizons.NoSuchHost  # type: ignore[attr-defined]  # absent by design


def test_import_does_not_import_oops() -> None:
    """Importing host_new_horizons leaves oops unimported."""

    src = pathlib.Path(__file__).parents[1] / 'src'
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join([str(src), env.get('PYTHONPATH', '')])
    code = 'import sys, host_new_horizons; print("oops" in sys.modules)'
    result = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True,
                            text=True, check=True)
    assert result.stdout.strip() == 'False'

##########################################################################################
