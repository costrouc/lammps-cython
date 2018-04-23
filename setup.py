from setuptools import setup, find_packages
from setuptools.extension import Extension
# To use a consistent encoding
from codecs import open
from os import path

# https://stackoverflow.com/questions/2379898/make-distutils-look-for-numpy-header-files-in-the-correct-place
try:
    from Cython.setuptools import build_ext
except:
    # If we couldn't import Cython, use the normal setuptools
    # and look for a pre-compiled .c file instead of a .pyx file
    from setuptools.command.build_ext import build_ext
    ext_modules = [Extension("lammps.core", sources=['lammps/core.cpp'], language='c++')]
else:
    # If we successfully imported Cython, look for a .pyx file
    ext_modules = [Extension("lammps.core", sources=['lammps/core.pyx'], language='c++')]

# https://stackoverflow.com/questions/2379898/make-distutils-look-for-numpy-header-files-in-the-correct-place
class CustomBuildExtCommand(build_ext):
    """build_ext command for use when numpy and mpi4py headers are needed."""
    def run(self):
        import numpy
        import mpi4py

        # Add mpi4py, numpy, and custom headers to include_dirs
        self.include_dirs.append(numpy.get_include())
        self.include_dirs.append(mpi4py.get_include())

        # not sure if necissary
        from configparser import ConfigParser
        lammps_config = ConfigParser()
        lammps_config.read('lammps.cfg')
        self.include_dirs.append(lammps_config.get('lammps', 'lammps_include_dir'))
        self.include_dirs.append(lammps_config.get('mpi', 'mpi_include_dir'))

        self.libraries = [
            lammps_config.get('lammps', 'lammps_library'),
            lammps_config.get('mpi', 'mpi_library'),
        ]

        self.library_dirs = [
            lammps_config.get('lammps', 'lammps_library_dir')
        ]

        # Call original build_ext command
        build_ext.run(self)


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
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=Extension("lammps.core", sources=['lammps/core.pyx'], language='c++'),
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={
        'lammps': ['data/*.in']
    },
    entry_points={
        'console_scripts': [
            'pylammps=lammps.__main__:main'
        ]
    },
    setup_requires=['pytest-runner', 'setuptools>=38.6.0'],
    tests_require=['pytest', 'pytest-cov'],
    install_requires=['mpi4py', 'numpy'],
)
