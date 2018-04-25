Installation
============

Installing lammps-cython can be a complex process. But from my
experience not too hard on most standard linux distributions. If
installation is hard for you I would recommend using the docker images
or conda (not done yet) which do not require manual installation of
`LAMMPS <http://lammps.sandia.gov/>`_.

Docker
------

This method should work for Linux, OSX, and windows. The docker image
includes the lammps shared library, python3.5, and lammps
binary. Everything you need to start experimenting with the package.

.. code-block:: bash

   docker pull costrouc/lammps-cython

Conda
-----

**Does Not Work Yet**

`Conda <https://github.com/conda/conda>`_ is an OS agnostic package
manager. It is developed by `contiuum analytics
<https://anaconda.io>`_. `Conda installation
<https://docs.anaconda.com/anaconda/install/>`_ is by far the easier
way to get started with python and installing complex scientific
codes. Once you have conda installed it will be as simple as

.. code-block:: bash

   conda install lammps-cython

Pip
---

Installation via pip requires building the LAMMPS shared library. This
is because pypi does not handle binaries in wheels well. Hopefully
this changes with `manylinux <https://github.com/pypa/manylinux>`_.

Building LAMMPS Shared Library (Linux and OSX)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I will attempt to keep this current on how to install in
Ubuntu. Should work similarly for all Linux distributions. I am not
knowledgeable on how to install on OSX or Windows (maybe this will
change). Before installing lammps make sure that you have all the
necessary dependencies. If you are working on a cluster there is a
good chance that you will not need to install the dependencies.

Linux Dependencies
""""""""""""""""""

For linux derivatives you should have packages very similarly named to
the ones below. I have installed on both debian and ubuntu and these
packages should work.

.. code-block:: bash

   apt install build-essential libopenblas-dev libfftw3-dev libopenmpi-dev

OSX Dependecies
"""""""""""""""

`Homebrew <http://brewformulas.org/>`_ is the choice of developer for
the package manager of OSX. See `homebrew installation instructions
<https://brew.sh/>`_ if you do not already have the package
manager. While I have minimal experience using homebrew these
instructions should work.

.. code-block:: bash

   brew install openblas mpich fftw

To install the equivalent of ``build-essential`` in linux you will
need to have xcode installed. This can be done in the terminal via
``xcode-select --install``. Not sure if necessary and probably already
installed on your machine.

Installation
^^^^^^^^^^^^

Make sure that your mpi implementation is MPI3. Next download the
latest LAMMPS `release <http://lammps.sandia.gov/download.html>`_ and
untar the download. I would actually recommend using the `LAMMPS
github repo <https://github.com/lammps/lammps/releases>`_ the archives
are smaller, download faster, and are easier to find. Replace
``<version>`` with your version you are downloading. ``wget`` is used
for downloading the archive from the terminal but there are of course
other ways.

.. code-block:: bash

   wget https://github.com/lammps/lammps/archive/<version>.tar.gz
   cd lammps-<version>/src

Inside of the lammps directory you will see a directory called
**src**. :command:`cd` to the src directory. Next make the LAMMPS
shared library. Building the shared library should take around 5
minutes without adding on packages. If you need to add packages to
LAMMPS this can easily be done. ``make package-status`` will list the
available packages and if they are activated. You can activate any
package by typing ``make yes-<package>``. To deactivate a package type
``make no-<package>``. As an example if you need to use the `manybody
<http://lammps.sandia.gov/doc/Section_packages.html#manybody-package>`_
package this can be installed via ``make yes-manybody``. Note that
some of the packages require additional complex dependencies.

.. code-block:: bash

   make mode=shlib mpi -j4
   cp liblammps_mpi.so /usr/local/lib/liblammps.so
   mkdir /usr/local/include/lammps/; cp *.h /usr/local/include/lammps/


Now all include files has been coppied to
``/usr/local/include/lammps/`` and the shared library has been coppied
to ``/usr/local/lib``. This will require root permissions. Another
great location that does not require root is ``$HOME/.local/lib/`` and
``$HOME/.local/include/``.

.. note::

   Ensure that the LAMMPS library is installed without any undefined
   symbols. To check the shared library file run :command:`ldd -r
   liblammps_mpi.so`. You should not see any undefined symbols. One
   reported issue is that the gfortran library is not included by
   default with LAMMPS when installing additional add-on packages. To
   fix this issue simply build the LAMMPS library with the gfortran
   shared library included :commmand:`-lgfortran`.

Installating lammps-cython
^^^^^^^^^^^^^^^^^^^^^^^^^^

lammps-cython installation should be easy if you have exactly followed
the steps above. If the LAMMPS shared library is located in
``/usr/local/lib``, the include files are located in
``/usr/local/include/lammps`` your installation should just work with
pip. Soon I will move to PyPi but it is not packaged well enough right
now.

.. code-block:: bash

   pip install numpy mpi4py cython
   pip install lammps-cython


If it does not you can do the simple manual installation.

.. code-block:: bash

   pip install numpy mpi4py
   pip download lammps-cython
   <not sure here>

I hope to soon support easy installation via the pip install
lammps. Download the lammps-python `source
<http://github.com/costrouc/lammps-python/tarbal/stable>`_. Untar
folder and cd into folder. Next install all the dependencies.

Common Installation Errors
--------------------------

There are some common errors that should be checked before submitting
an issue on the github repository.

.. code-block:: python

   >>> import lammps
   from .core import Lammps

   ImportError: liblammps.so: cannot open shared object file: No such file or directory


This error results because python cannot find the LAMMPS
library. Meaning that the lammps library is the not in the standard
library search path. On a typical linux system these paths are
:command:`/usr/lib` and :command:`/usr/local/lib`. If you would like
to have the LAMMPS library in another directory not in the standard
path you must modify the environment variable
:command:`LD_LIBRARY_PATH`.

For any other errors PLEASE add an `issue to the gitlab page
<https://gitlab.com/costrouc/lammps-cython>`_. I check gitlab often
and really want to make this a long-term supported addition to the
LAMMPS community!
