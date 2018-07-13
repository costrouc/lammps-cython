"""Testing package integration

pymatgen
ase
"""
import math

import numpy as np

import lammps


def test_pymatgen_add_structure():
    import pymatgen as pmg

    a = 4.2
    alpha = 90
    alpha_r = math.radians(alpha)
    lattice = pmg.Lattice.from_parameters(a, a, a, alpha, alpha, alpha)
    symbols = ['Mg', 'O']
    positions = [(0, 0, 0), (0.5, 0.5, 0.5)]

    structure = pmg.Structure.from_spacegroup(1, lattice, symbols, positions)
    velocities = np.random.random((len(structure), 3))
    structure.add_site_property('velocities', velocities)
    elements = list(set(structure.species))

    lmp = lammps.Lammps()
    lmp.system.add_pymatgen_structure(structure, elements=elements)
    assert lmp.system.total == len(structure)
    assert np.all(np.isclose(lmp.box.lengths, (a, a, a)))
    assert np.all(np.isclose(lmp.box.angles, (alpha_r, alpha_r, alpha_r)))
    assert np.all(np.isclose(lmp.box.origin, (0, 0, 0)))
    print(velocities)
    print(np.array(structure.site_properties['velocities']))
    print(lmp.system.velocities)
    assert np.all(np.isclose(structure.site_properties['velocities'], velocities))
    assert np.all(np.isclose(lmp.system.velocities, velocities))



def test_ase_add_structure():
    pass
