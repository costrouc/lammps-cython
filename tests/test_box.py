import os
from math import pi
import functools

import pytest
import numpy as np

import lammps

# Redefine Lammps command-line args so no annoying logs or stdout
Lammps = functools.partial(lammps.Lammps, args=[
    '-log', 'none',
    '-screen', 'none'
])


def test_box():
    pass


def test_lammps_initial_box():
    lammps = Lammps()
    assert lammps.box.lengths == (1., 1., 1.)
    assert np.all(np.isclose(lammps.box.angles, (pi/2., pi/2., pi/2.)))
