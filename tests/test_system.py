import pytest
import numpy as np


def test_system_lammps_initial(lmp):
    num_atoms = 0
    assert lmp.system.total == num_atoms
    # assert lmp.system.types == None
    # assert lmp.system.positions == None
    # assert lmp.system.velocities == None
    # assert lmp.system.forces == None


def test_system_melt(lmp_melt):
    num_atoms = 16
    assert lmp_melt.system.total == num_atoms
    assert np.all(lmp_melt.system.types == [1] * num_atoms)
    assert len(lmp_melt.system.atom_types) == 1
    assert np.isclose(lmp_melt.system.atom_types[0].mass, 1.)
    assert lmp_melt.system.atom_types[0].index == 1
    assert lmp_melt.system.total == num_atoms
    assert len(lmp_melt.system.tags) == num_atoms
    assert np.all(lmp_melt.system.tags.reshape(-1,) == np.arange(num_atoms)+1)
    assert lmp_melt.system.positions.shape == (num_atoms, 3)
    assert lmp_melt.system.velocities.shape == (num_atoms, 3)
    assert lmp_melt.system.forces.shape == (num_atoms, 3)
