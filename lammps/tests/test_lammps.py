import unittest

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
    '-screen', 'none',
])


class LammpsBaseTest(unittest.TestCase):

    #Input Arguments
    ## Units
    def test_default_units(self):
        lmp = Lammps_partial()
        self.assertEqual(lmp.units, 'lj')

    def test_set_units(self):
        lmp = Lammps_partial(units='metal')
        self.assertEqual(lmp.units, 'metal')

    def test_bad_units(self):
        self.assertRaises(ValueError, Lammps_partial, units='bad_units')

    ## MPI_Comm comm (don't want to do mpi tests)
    def test_bad_comm(self):
        self.assertRaises(TypeError, Lammps_partial, comm='bad_comm')

    ## Style
    def test_default_style(self):
        lmp = Lammps_partial()
        self.assertEqual(lmp.system.style, 'atomic')

    def test_set_style(self):
        lmp = Lammps_partial(style='full')
        self.assertEqual(lmp.system.style, 'full')

    def test_bad_style(self):
        self.assertRaises(ValueError, Lammps_partial, style='bad_style')

    # version
    def test_version(self):
        lmp = Lammps_partial()
        self.assertIsInstance(lmp.__version__, str)

    # command
    def test_command(self):
        lmp = Lammps_partial()
        lmp.command('timestep 2.0')
        self.assertEqual(lmp.dt, 2.0)

    # file
    def test_file(self):
        lmp = Lammps_partial()
        lmp.file('in.test')
        self.assertEqual(lmp.dt, 0.006)
        self.assertEqual(lmp.units, 'metal')
        self.assertEqual(lmp.system.style, 'charge')

    # def test_run(self):
    #     lmp = Lammps_partial()
    #     lmp.file('in.test')
    #     lmp.run(100)
    #     self.assertEqual(lmp.time_step, 100)
    #     self.assertEqual(lmp.time, 100 * 0.006)

    # dt
    def test_default_dt(self):
        lmp = Lammps_partial()
        self.assertEqual(lmp.dt, 0.005)

    def test_set_dt(self):
        lmp = Lammps_partial()
        lmp.dt = 13.0
        self.assertEqual(lmp.dt, 13.0)

    # reset
    def test_reset(self):
        lmp = Lammps_partial()
        lmp.dt = 13.0
        lmp.reset()
        self.assertEqual(lmp.dt, 0.005)


    
    
if __name__ == '__main__':
    unittest.main()
