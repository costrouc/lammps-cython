import os
import functools

import pytest
import lammps


# Redefine Lammps command-line args so no annoying logs or stdout
Lammps = functools.partial(lammps.Lammps, args=[
    '-log', 'none',
    '-screen', 'none'
])


# Lammps command line arguments
def test_lammps_init_default_units():
    lmp = Lammps()
    assert lmp.units == 'lj'


def test_lammps_init_set_get_units():
    lmp = Lammps(units='metal')
    assert lmp.units == 'metal'


def test_lammps_init_invalid_units():
    with pytest.raises(ValueError):
        Lammps(units='invalid_units')


# MPI_Comm comm (don't want to do mpi tests)
def test_lammps_init_invalid_comm():
    with pytest.raises(TypeError):
        Lammps(comm='invalid_comm')


# Style
def test_lammps_init_default_style():
    lmp = Lammps()
    assert lmp.system.style == 'atomic'


@pytest.mark.skipif(
    not pytest.config.getoption('--pkg-molecule'),
    reason='requires molecule package')
def test_lammps_init_set_get_style():
    lmp = Lammps(style='full')
    assert lmp.system.style == 'full'


def test_lammps_init_invalid_style():
    with pytest.raises(ValueError):
        Lammps(style='invalid_style')


# version
def test_lammps_version():
    lmp = Lammps()
    assert isinstance(lmp.__version__, str)


# command
def test_lammps_command():
    lmp = Lammps()
    lmp.command('timestep 2.0')
    assert lmp.dt == 2.0


# file
def test_lammps_file(tmpdir):
    tmpfile = tmpdir.join("test_file.in")
    tmpfile.write("timestep 1.0\n")

    lmp = Lammps()
    lmp.file(str(tmpfile))
    assert lmp.dt == 1.0


def test_lammps_file_twice(tmpdir):
    tmpfile1 = tmpdir.join("test_file1.in")
    tmpfile1.write("timestep 1.0\n")

    tmpfile2 = tmpdir.join("test_file2.in")
    tmpfile2.write("timestep 2.0\n")

    lmp = Lammps()

    lmp.file(str(tmpfile1))
    assert lmp.dt == 1.0

    lmp.file(str(tmpfile2))
    assert lmp.dt == 2.0


# Run
def test_lammps_run():
    # This tests has a dependency of the
    # LAMMPS example melt
    # dt tested
    # time step tested
    # time tested
    # This is hardly a unit test... (a better way?)

    lmp = Lammps()
    lmp.file(os.path.join(lammps.__path__[0], 'data', 'melt.in'))

    assert lmp.dt == 0.005
    assert lmp.time_step == 100
    assert lmp.time == lmp.time_step * lmp.dt


# time_step
def test_lammps_default_time_step():
    lmp = Lammps()
    assert lmp.time_step == 0


def test_lammps_set_get_time_step():
    lmp = Lammps()
    lmp.time_step = 100
    assert lmp.time_step == 100


# dt
def test_lammps_default_dt():
    lmp = Lammps()
    assert lmp.dt == 0.005


def test_lammps_set_get_dt():
    lmp = Lammps()
    lmp.dt = 13.0
    assert lmp.dt == 13.0


# time
def test_lammps_default_time():
    lmp = Lammps()
    assert lmp.time == 0.0


# reset
def test_lammps_reset():
    lmp = Lammps()
    lmp.dt = 13.0
    lmp.reset()
    assert lmp.dt == 0.005
