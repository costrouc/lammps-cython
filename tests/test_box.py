from math import pi

import pytest
import numpy as np

import lammps


def test_lattice_const_to_lammps_box_cubic():
    lengths = (5, 5, 5)
    angles = (pi/2, pi/2, pi/2)
    origin = (0, 0, 0)
    a, b, c = lengths
    xlo, ylo, zlo = origin
    bounds, tilts, rotation_matrix = lammps.core.lattice_const_to_lammps_box(lengths, angles)
    assert np.all(np.isclose(bounds, [[xlo, xlo+a], [ylo, ylo+b], [zlo, zlo+c]]))
    assert np.all(np.isclose(tilts, (0, 0, 0)))
    assert np.all(np.isclose(rotation_matrix, np.eye(3)))


def test_lattice_const_to_lammps_box_cubic_offset_origin():
    lengths = (5, 5, 5)
    angles = (pi/2, pi/2, pi/2)
    origin = (4, 3, 2)
    a, b, c = lengths
    xlo, ylo, zlo = origin
    bounds, tilts, rotation_matrix = lammps.core.lattice_const_to_lammps_box(lengths, angles, origin=origin)
    assert np.all(np.isclose(bounds, [[xlo, xlo+a], [ylo, ylo+b], [zlo, zlo+c]]))
    assert np.all(np.isclose(tilts, (0, 0, 0)))
    assert np.all(np.isclose(rotation_matrix, np.eye(3)))


def test_lattice_const_to_lammps_box_rhomb():
    # 3C-SiC
    lengths = (3.0968, 3.0968, 3.0968)
    angles = (pi/3, pi/3, pi/3)
    bounds, tilts, rotation_matrix = lammps.core.lattice_const_to_lammps_box(lengths, angles)
    assert np.all(np.isclose(bounds, ((0, 3.0968), (0, 2.6819074704396493), (0, 2.528526611816982)), atol=1e-3))
    assert np.all(np.isclose(tilts, (1.5484000000000004, 1.5484000000000004, 0.8939691568132165)))


def test_lammps_box_to_lattice_const_cubic():
    bounds = [[0, 5], [0, 5], [0, 5]]
    tilts = (0, 0, 0)
    origin = (0, 0, 0)
    lengths, angles, origin = lammps.core.lammps_box_to_lattice_const(bounds, tilts)
    assert np.all(np.isclose(lengths, (5, 5, 5)))
    assert np.all(np.isclose(angles, (pi/2, pi/2, pi/2)))


def test_lammps_box_orthogonal_reversible():
    lengths = (4, 4, 4)
    angles = (pi/2, pi/2, pi/2)
    origin = (1, 2, 3)
    bounds, tilts, rotation_matrix = lammps.core.lattice_const_to_lammps_box(lengths, angles, origin=origin)
    lengths_r, angles_r, origin_r = lammps.core.lammps_box_to_lattice_const(bounds, tilts)
    assert np.all(np.isclose(lengths, lengths_r))
    assert np.all(np.isclose(angles, angles_r))
    assert np.all(np.isclose(origin, origin_r))


def test_lammps_box_tetrahedral_reversible():
    # LiTaO3
    lengths = (5.5338, 5.5338, 5.5338)
    angles = (56.14486291 * pi/180, 56.14486291 * pi/180, 56.14486291 * pi/180)
    origin = (1, 2, 3)
    bounds, tilts, rotation_matrix = lammps.core.lattice_const_to_lammps_box(lengths, angles, origin=origin)
    lengths_r, angles_r, origin_r = lammps.core.lammps_box_to_lattice_const(bounds, tilts)
    assert np.all(np.isclose(lengths, lengths_r))
    assert np.all(np.isclose(angles, angles_r))
    assert np.all(np.isclose(origin, origin_r))



def test_lammps_initial_box(lmp):
    assert lmp.box.dimension == 3
    assert np.all(np.isclose(lmp.box.lengths, (1., 1., 1.)))
    assert np.all(np.isclose(lmp.box.angles, (pi/2., pi/2., pi/2.)))
    assert np.all(np.isclose(lmp.box.bounds, [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]))
    assert np.all(np.isclose(lmp.box.tilts, [0, 0, 0]))
    assert np.all(np.isclose(lmp.box.lengths_angles, [[1, 1, 1], [pi/2, pi/2, pi/2]]))
    # lammps has some seriously weird initial behavior
    # has unit cell 1x1x1 with volume 0 ???
    # actually has non-deterministic behavior 0 or inf
    # assert np.isclose(lmp.box.volume, 0.)


def test_lammps_set_box_from_lattice_const(lmp):
    atom_types = 5
    lengths = (10, 10, 10)
    angles = (pi/2., pi/2., pi/2.)
    lmp.box.from_lattice_const(atom_types, lengths, angles)
    assert np.all(np.isclose(lmp.box.lengths, lengths))
    assert np.all(np.isclose(lmp.box.angles, angles))
    assert lmp.system.total == 0
    assert len(lmp.system.atom_types) == atom_types
    assert np.isclose(lmp.box.volume, 10**3)


def test_lammps_update_lattice_const(lmp):
    lengths = (10, 10, 10)
    angles = (pi/2., pi/2., pi/2.)
    lmp.box.update_lattice_const(lengths, angles)
    assert np.all(np.isclose(lmp.box.lengths, lengths))
    assert np.all(np.isclose(lmp.box.angles, angles))
    assert np.isclose(lmp.box.volume, 10**3)
