from setuptools import setup, find_packages
from setuptools.extension import Extension
# To use a consistent encoding
from codecs import open
from os import path

# https://stackoverflow.com/questions/2379898/make-distutils-look-for-numpy-header-files-in-the-correct-place
try:
    from Cython.Distutils import build_ext as _build_ext
except:
    # If we couldn't import Cython, use the normal setuptools
    # and look for a pre-compiled .c file instead of a .pyx file
    from setuptools.command.build_ext import build_ext as _build_ext
    ext_modules = [Extension("lammps.core", sources=['lammps/core.cpp'], language='c++')]
else:
    # If we successfully imported Cython, look for a .pyx file
    ext_modules = [Extension("lammps.core", sources=['lammps/core.pyx'], language='c++')]

# https://stackoverflow.com/questions/2379898/make-distutils-look-for-numpy-header-files-in-the-correct-place
class build_ext(_build_ext):
    """needed to import numpy and mpi4py includes"""
    def run(self):
        import numpy
        import mpi4py

        # not sure if necissary
        from configparser import ConfigParser
        import os

        lammps_config = ConfigParser()
        # Give precedence to config file in user home config directory
        if os.path.isfile(os.path.expanduser('~/.config/lammps-site.cfg')):
            lammps_config.read(os.path.expanduser('~/.config/lammps-site.cfg'))
        else:
            lammps_config.read('lammps.cfg')

        def config_to_list(key1, key2):
            return [s.strip() for s in lammps_config.get(key1, key2).split(',')]

        # Add mpi4py, numpy, and custom headers to include_dirs
        self.include_dirs.extend([
            numpy.get_include(),
            mpi4py.get_include()
        ])
        self.include_dirs.extend(config_to_list('lammps', 'lammps_include_dir'))
        self.include_dirs.extend(config_to_list('mpi', 'mpi_include_dir'))

        self.library_dirs.extend(config_to_list('lammps', 'lammps_library_dir'))
        self.library_dirs.extend(config_to_list('mpi', 'mpi_library_dir'))

        self.libraries.extend(config_to_list('lammps', 'lammps_library'))
        self.libraries.extend(config_to_list('mpi', 'mpi_library'))

        # Call original build_ext command
        _build_ext.run(self)


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='lammps-cython',
    version='0.5.8',
    description='Pythonic Wrapper to LAMMPS using cython',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Christopher Ostrouchov',
    author_email='chris.ostrouchov+lammps@gmail.com',
    license="MIT",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    url='https://gitlab.com/costrouc/lammps-cython',
    download_url='https://gitlab.com/costrouc/lammps-cython/-/archive/master/lammps-cython-master.zip',
    keywords=['lammps', 'molecular dynamics', 'cython', 'wrapper', 'mpi'],
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    packages=find_packages(exclude=['docs', 'tests']),
    entry_points={
        'console_scripts': [
            'pylammps=lammps.__main__:main'
        ]
    },
    setup_requires=['pytest-runner', 'setuptools>=38.6.0'],
    tests_require=['pytest', 'pytest-cov'],
    install_requires=['mpi4py', 'numpy'],
    extras_require={
        'all': ['pymatgen', 'ase'],
        'pymatgen': 'pymatgen',
        'ase': 'ase'
    }
)
