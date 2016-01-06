from Cython.Build import cythonize

import mpi4py
import numpy

from distutils.core import setup
from distutils.extension import Extension


sourcefiles = ['lammps/lammps.pyx']

extension =  Extension('lammps',
                       sources=sourcefiles,
                       include_dirs=['/home/costrouc/software/lammps-7Dec15/src',
                                     '/usr/lib/mpich/include',
                                     mpi4py.get_include(),
                                     numpy.get_include()],
                       libraries=['lammps'],
                       library_dirs=['/home/costrouc/.local/lib'],
                       language='c++')

setup(
    name='lammps',
    version='0.0.1',
    author='Christopher Ostrouchov',
    description='A Cython Interface to LAMMPS',
    ext_modules=cythonize(extension),
)
