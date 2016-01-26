import pytest

from lammps.core import Lammps

import functools

# TODO
#   file
#   run(int)
#   time_step
#   time

# Redefine Lammps command-line args so no annoying 
# logs or stdout
Lammps_partial = functools.partial(Lammps, args=[
    '-log', 'none',
    '-screen', 'none'
])

# Lammps command line arguments
def test_lammps_init_default_units():
    lmp = Lammps_partial()
    assert lmp.units == 'lj'


def test_lammps_init_set_units():
    lmp = Lammps_partial(units='metal')
    assert lmp.units == 'metal'


def test_lammps_init_invalid_units():
    with pytest.raises(ValueError):
        Lammps_partial(units='invalid_units')


# MPI_Comm comm (don't want to do mpi tests)
def test_lammps_init_invalid_comm():
    with pytest.raises(TypeError):
        Lammps_partial(comm='invalid_comm')


# Style
def test_lammps_init_default_style():
    lmp = Lammps_partial()
    assert lmp.system.style == 'atomic'


def test_lammps_init_set_style():
    lmp = Lammps_partial(style='full')
    assert lmp.system.style == 'full'


def test_lammps_init_invalid_style():
    with pytest.raises(ValueError):
        Lammps_partial(style='invalid_style')


# version
def test_lammps_version():
    lmp = Lammps_partial()
    assert isinstance(lmp.__version__, str)


# command
def test_lammps_command():
    lmp = Lammps_partial()
    lmp.command('timestep 2.0')
    assert lmp.dt == 2.0


# file
def test_lammps_file(tmpdir):
    tmpfile = tmpdir.join("test_file.in")
    tmpfile.write("timestep 1.0\n")

    lmp = Lammps_partial()
    lmp.file(str(tmpfile))
    assert lmp.dt == 1.0


def test_lammps_file_twice(tmpdir):
    tmpfile1 = tmpdir.join("test_file1.in")
    tmpfile1.write("timestep 1.0\n")

    tmpfile2 = tmpdir.join("test_file2.in")
    tmpfile2.write("units full")

    lmp = Lammps_partial()

    lmp.file(str(tmpfile1))
    assert lmp.dt == 1.0
    assert lmp.units == 'lj'

    lmp.file(str(tmpfile2))
    assert lmp.dt == 1.0
    assert lmp.units == 'full'

# def test_run():
#     lmp = Lammps_partial()
#     lmp.file('in.test')
#     lmp.run(100)
#     .assertEqual(lmp.time_step, 100)
#     .assertEqual(lmp.time, 100 * 0.006)


# dt
def test_lammps_default_dt():
    lmp = Lammps_partial()
    assert lmp.dt == 0.005


def test_lammps_set_dt():
    lmp = Lammps_partial()
    lmp.dt = 13.0
    assert lmp.dt == 13.0


# reset
def test_lammps_reset():
    lmp = Lammps_partial()
    lmp.dt = 13.0
    lmp.reset()
    assert lmp.dt == 0.005


if __name__ == '__main__':
    pytest.main()
