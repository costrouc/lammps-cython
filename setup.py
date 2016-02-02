from configparser import ConfigParser

import mpi4py
import numpy

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

lammps_config = ConfigParser()
lammps_config.read('lammps.cfg')

include_dirs = [
    mpi4py.get_include(),
    numpy.get_include(),
    lammps_config.get('lammps', 'lammps_include_dir'),
    lammps_config.get('mpi', 'mpi_include_dir')
]

# TODO: Should maybe include mpi_cxx, mpi, python3.4m
libraries = [lammps_config.get('lammps', 'lammps_library'), 
             lammps_config.get('mpi', 'mpi_library')]
library_dirs = [lammps_config.get('lammps', 'lammps_library_dir')]

extensions = [
    Extension(
        'lammps.core',
        sources=['lammps/core.pyx'],
        include_dirs=include_dirs,
        libraries=libraries,
        library_dirs=library_dirs,
        language='c++'
    )
]

setup(
    name='lammps',
    version='0.1.0',
    packages=find_packages(),
    package_data={
        'lammps': ['data/*.in']
    },
    description='Pythonic Wrapper to LAMMPS',
    long_description='Pythonic Wrapper to LAMMPS (LONG)',
    author='Christopher Ostrouchov',
    author_email='chris.ostrouchov+lammps@gmail.com',
    url='https://github.com/costrouc/lammps-python',
    download_url='https://github.com/costrouc/lammps-python/tarball/master',
    keywords=['lammps', 'molecular dynamics', 'cython', 'wrapper'],
    ext_modules=cythonize(extensions),
    scripts=['scripts/pylammps'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
