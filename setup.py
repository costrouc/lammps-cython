import mpi4py
import numpy

from Cython.Build import cythonize

from distutils.core import setup
from distutils.extension import Extension

from configparser import ConfigParser

lammps_config = ConfigParser()
lammps_config.read('lammps.cfg')

include_dirs = [
    mpi4py.get_include(),
    numpy.get_include(),
    lammps_config.get('lammps', 'lammps_include_dir'),
    lammps_config.get('mpi', 'mpi_include_dir')
]

# TODO: Should maybe include mpi_cxx, mpi, python3.4m
libraries = [lammps_config.get('lammps', 'lammps_library')]
library_dirs = [lammps_config.get('lammps', 'lammps_library_dir')]

extension = Extension(
    'lammps.core',
    sources=['lammps/core.pyx'],
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs,
    language='c++'
)

setup(
    name='lammps',
    packages=['lammps'],
    version='0.0.1',
    description='Pythonic Wrapper to LAMMPS',
    author='Christopher Ostrouchov',
    author_email='chris.ostrouchov+lammps@gmail.com',
    url='https://github.com/costrouc/lammps-python',
    download_url='https://github.com/costrouc/lammps-python/tarball/master',
    keywords=['lammps', 'molecular dynamics', 'cython', 'wrapper'],
    ext_modules=cythonize(extension)
)
