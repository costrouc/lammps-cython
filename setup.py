from configparser import ConfigParser

import mpi4py
import numpy

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
# To use a consistent encoding
from codecs import open
from os import path


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

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='lammps',
    version='0.1.0',
    description='Pythonic Wrapper to LAMMPS',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Christopher Ostrouchov',
    author_email='chris.ostrouchov+lammps@gmail.com',
    license="MIT",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.3'
        'Programming Language :: Python :: 3.4'
        'Programming Language :: Python :: 3.5'
        'Programming Language :: Python :: 3.6'
    ],
    url='https://gitlab.com/costrouc/lammps-cython',
    download_url='https://gitlab.com/costrouc/lammps-cython/-/archive/master/lammps-cython-master.zip',
    keywords=['lammps', 'molecular dynamics', 'cython', 'wrapper', 'mpi'],
    ext_modules=cythonize(extensions),
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={
        'lammps': ['data/*.in']
    },
    entry_points={
        'console_scripts': [
            'pylammps=lammps.__main__:main'
        ]
    },
    scripts=['scripts/pylammps'],
    setup_requires=['pytest-runner', 'setuptools>=38.6.0'],
    tests_require=['pytest', 'pytest-cov'],
)
