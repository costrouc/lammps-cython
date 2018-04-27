# lammps-cython

A high-performance cython wrapper around LAMMPS. Lammps is a great
molecular dynamics package that has an unmatched set of potentials and
fixes. This package offers unique features such as minimizing I/O by
allowing direct access to thermostats and atom properties and allowing
interactive lammps within python interpreters such a `ipython`.  The
goal of this project is to put an opinionated wrapper around LAMMPS
(the good parts) and allow the user to easily extend it's
functionality in python. The api should feel very similar to
[HOOMD](https://codeblue.umich.edu/hoomd-blue/) and is being actively
developed.

<table>
<tr>
  <td>Latest Release</td>
  <td><img src="https://img.shields.io/pypi/v/lammps-cython.svg" alt="latest release"/></td>
</tr>
<tr>
  <td>Package Status</td>
  <td><img src="https://img.shields.io/pypi/status/lammps-cython.svg" alt="status" /></td>
</tr>
<tr>
  <td>License</td>
  <td><img src="https://img.shields.io/pypi/l/lammps-cython.svg" alt="license" /></td>
</tr>
<tr>
  <td>Build Status</td>
  <td> <a href="https://gitlab.com/costrouc/lammps-cython/pipelines"> <img
src="https://gitlab.com/costrouc/lammps-cython/badges/master/pipeline.svg"
alt="gitlab pipeline status" /> </a> </td>
</tr>
</table>


# Documentation

Full documentation can be found at
[lammps-cython](https://costrouc.gitlab.io/lammps-cython/).

# Features

 - Full MPI support
 - Pythonic API inspired by
 [HOOMD](https://codeblue.umich.edu/hoomd-blue/)
 - Supports Python 2 and 3
 - Heavily documented and tested
 - Elimination of unnecessary file I/O for thermostats and atoms properties

A neat feature of the wrapper is that lammps can be run regularly
using the following script (use "-i" instead of stdin). This is the
command `pylammps` when the package is installed.

```python
from lammps import Lammps
import sys
Lammps(args=sys.args)
```

# Installation

`lammps-cython` has several options for installation. The easiest way
is using the provided docker containter image
[costrouc/lammps-cython](https://hub.docker.com/r/costrouc/lammps-cython/). There
are plans to support conda and pip wheels. However currently other
methods require manual installation of lammps. Detailed installation
are provieded in the
[documentation](https://costrouc.gitlab.io/lammps-cython/installation.html). If
you have any issues with installation be submit an issue at the
[gitlab repository](https://gitlab.com/costrouc/lammps-cython/).

The general path to installation is install [LAMMPS as a shared
library](http://lammps.sandia.gov/doc/Section_start.html#start-4) then
edit `~/.config/lammps-site.cfg` to include the paths of necissary
libraries. See example below.

``` ini
[lammps]
lammps_include_dir = /usr/local/include/lammps/
lammps_library_dir = /usr/local/lib/
# true library filename is liblammps.so notice lib and .so are removed
lammps_library = lammps

# use mpic++ -showme to list libraries and includes
[mpi]
mpi_include_dir = /usr/lib/x86_64-linux-gnu/openmpi/include
mpi_library_dir = /usr/lib/x86_64-linux-gnu/openmpi/lib
# no necissarily needed (default are mpi, mpi_cxx)
mpi_library     = mpi
```

Then `pip install lammps-cython` should just work.

## Docker Image

The docker image
[costrouc/lammps-cython](https://hub.docker.com/r/costrouc/lammps-cython/)
uses `python3.5` and has the library preinstalled with the executables
`pylammps` and `lammps` available.


# Tutorials

Work is being done to show how to use the features of `lammps-cython`
for now just visit the [tutorial page](https://costrouc.gitlab.io/lammps-cython/tutorial.html).

These will turn to links when the tutorial exists.

  - basic usage
  - modify atom positions
  - get forces and velocity for each atom and compute potential energy

# Contributing

All contributions, bug reports, bug fixes, documentation improvements,
enhancements and ideas are welcome!

Contributors:

  - [Chris Ostrouchov](https://gitlab.com/costrouc) (maintainer)

# License

MIT
