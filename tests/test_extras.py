"""Testing package integration

 - pymatgen
 - ase
"""
import math

import numpy as np


def test_pymatgen_add_and_get_structure(lmp):
    import pymatgen as pmg

    a = 4.2
    alpha = 90
    alpha_r = math.radians(alpha)
    lattice = pmg.Lattice.from_parameters(a, a, a, alpha, alpha, alpha)
    symbols = ['Mg', 'O']
    positions = [(0, 0, 0), (0.5, 0.5, 0.5)]

    structure = pmg.Structure.from_spacegroup(225, lattice, symbols, positions)
    velocities = np.random.random((len(structure), 3))
    structure.add_site_property('velocities', velocities)
    elements = list(set(structure.species))

    lmp.system.add_pymatgen_structure(structure, elements=elements)
    assert lmp.system.total == len(structure)
    assert np.all(np.isclose(lmp.box.lengths, (a, a, a)))
    assert np.all(np.isclose(lmp.box.angles, (alpha_r, alpha_r, alpha_r)))
    assert np.all(np.isclose(lmp.box.origin, (0, 0, 0)))
    assert np.all(np.isclose(lmp.system.positions, structure.cart_coords, atol=1e-6))
    assert np.all(np.isclose(lmp.system.velocities, velocities))


def test_pymatgen_snapshot(lmp_melt):
    symbols = ['Ne']
    num_atoms = 16
    new_structure = lmp_melt.system.snapshot(symbols, format='pymatgen')
    assert len(new_structure) == num_atoms


def test_gsd_snapshot(lmp_melt):
    symbols = ['Ne']
    num_atoms = 16
    snapshot = lmp_melt.system.snapshot(symbols, format='snapshot')
    assert snapshot.particles.N == num_atoms
    snapshot.validate()


def test_ase_add_structure(lmp):
    from ase import spacegroup

    a = 4.2
    alpha = 90
    alpha_r = math.radians(alpha)
    symbols = ['Mg', 'O']
    positions = [(0, 0, 0), (0.5, 0.5, 0.5)]

    structure = spacegroup.crystal(
        symbols, basis=positions, spacegroup=225,
        cellpar=[a, a, a, alpha, alpha, alpha])
    velocities = np.random.random((len(structure), 3))
    structure.set_velocities(velocities)
    elements = list({(atom.symbol, atom.mass) for atom in structure})

    lmp.system.add_ase_structure(structure, elements=elements)
    assert lmp.system.total == len(structure)
    assert np.all(np.isclose(lmp.box.lengths, (a, a, a)))
    assert np.all(np.isclose(lmp.box.angles, (alpha_r, alpha_r, alpha_r)))
    assert np.all(np.isclose(lmp.box.origin, (0, 0, 0)))
    print(lmp.system.velocities)
    print(structure.get_velocities())
    print(lmp.system.velocities - structure.get_velocities())
    assert np.all(np.isclose(lmp.system.positions, structure.positions, atol=1e-6))
    assert np.all(np.isclose(lmp.system.velocities, velocities, atol=1e-6))
