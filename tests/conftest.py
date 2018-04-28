import pytest
import functools


def pytest_addoption(parser):
    parser.addoption('--pkg-molecule', action='store_true', default=False, help='test molecule package')


@pytest.fixture
def lmp():
    import lammps
    # Redefine Lammps command-line args so no annoying logs or stdout
    Lammps = functools.partial(lammps.Lammps, args=[
        '-log', 'none',
        '-screen', 'none'
    ])
    return Lammps()


@pytest.fixture
def lmp_melt(lmp):
    melt_in = """
    lattice         bcc 1
    region          box block 0 2 0 2 0 2
    create_box      1 box
    create_atoms    1 box
    mass            1 1.0

    velocity        all create 10.0 123

    pair_style	lj/cut 1.5
    pair_coeff	1 1 1.0 1.0 2.5

    fix		1 all nve

    run 5
    """
    for line in melt_in.split('\n'):
        lmp.command(line)
    return lmp
