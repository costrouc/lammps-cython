"""This file exists to understand the behavior of commands

These tests do not really show off new features of api
"""
import pytest
import numpy as np

from lammps.core import LammpsNormalError

def test_command_atom_style(lmp):
    assert lmp.system.style == 'atomic'
    lmp.command('atom_style charge')
    assert lmp.system.style == 'charge'


def test_unknown_command(lmp):
    unknown_command = 'asdfasdf'
    with pytest.raises(LammpsNormalError) as error:
        lmp.command(unknown_command)
    assert 'ERROR: Unknown command: %s' % unknown_command in error.value.message
