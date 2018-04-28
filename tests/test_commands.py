"""This file exists to understand the behavior of commands

These tests do not really show off new features of api
"""

import numpy as np


def test_command_atom_style(lmp):
    assert lmp.system.style == 'atomic'
    lmp.command('atom_style charge')
    assert lmp.system.style == 'charge'
