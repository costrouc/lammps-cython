from Cython.Build import cythonize

import mpi4py
import numpy

from distutils.core import setup
from distutils.extension import Extension
from configparser import ConfigParser

lammps_config_file = 'lammps.cfg'
lammps_config = ConfigParser()
lammps_config.read(lammps_config_file)

include_dirs = [
    mpi4py.get_include(),
    numpy.get_include()
]
include_dirs.append(lammps_config.get('lammps', 'lammps_include_dir'))
include_dirs.append(lammps_config.get('mpi', 'mpi_include_dir'))

libraries = [lammps_config.get('lammps', 'lammps_library')]
library_dirs = [lammps_config.get('lammps', 'lammps_library_dir')]

source_files = ['lammps/lammps.pyx']

extensions =  Extension('lammps',
                        sources=source_files,
                        include_dirs=include_dirs,
                        libraries=libraries,
                        library_dirs=['/home/costrouc/.local/lib'],
                        language='c++')

setup(
    name='lammps',
    version='0.0.1',
    author='Christopher Ostrouchov',
    description='Pythonic Wrapper to LAMMPS',
    ext_modules=cythonize(extensions, compiler_directives={'embedsignature': True}),
)
