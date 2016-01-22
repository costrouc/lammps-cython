Installation
============

Installing lammpy-python is an easy process. But from my experience
**easy** is always relative. lammps-python is Python 2 and 3 compatible!

Dependencies
------------

 - Some MPI implementation (preferably implementing the MPI3 api)  
 - `mpi4py <https://bitbucket.org/mpi4py/mpi4py/>`_
 - `numpy <http://www.numpy.org/>`_
 - `cython <http://cython.org/>`_

Creating LAMMPS Library
-----------------------

I will attempt to keep this current (1/10/2016) on how to install in
Ubuntu. Should work similarly for all Linux distributions. I am not
knowledgeable on how to install on OSX or Windows (maybe this will
change). Before installing lammps make sure that you have all the
necessary dependencies. If you are working on a cluster there is a
good chance that you will not need to install the dependencies.

.. code-block:: bash

   $ sudo apt install mpich libfftw3-dev libpng-dev libjpeg-dev

Make sure that your mpi implementation is MPI3. Next download the
latest LAMMPS `release <http://lammps.sandia.gov/download.html>`_ and
untar the download.

.. code-block:: bash

   $ tar -xf <lammps-download>.tar.gz

Inside of the lammps directory you will see a directory called
**src**. :command:`cd` to the src directory. Next make the lammps
shared library. Building the shared library should take around 5
minutes.

.. code-block:: bash

   $ cd <lammps_download>/src
   $ make mode=shlib ubuntu

Remember two important paths.
 - path to lammps src directory <lammps_include_dir>
 - name of the lammps shared library <lammps_library>

.. note::

   Ensure that the LAMMPS library is installed without any undefined
   symbols. To check the shared library file run :command:`ldd -r
   <lammps_library>`. You should not see any undefined symbols. One
   reported issue is that the gfortran library is not included by
   default with LAMMPS when installing additional packages. To fix
   this issue simply build the LAMMPS library with the gfortran shared
   library included :commmand:`-lgfortran`.


Installing lammps-python
------------------------

I hope to soon support easy installation via the pip install
lammps. Download the lammps-python `source
<http://github.com/costrouc/lammps-python/tarbal/stable>`_. Untar
folder and cd into folder. Next install all the dependencies.

.. code-block:: bash

   $ pip install -r requirements.txt

If :command:`pip install` does not work simply make sure that mpi4py,
numpy, and cython are installed and available to python.

Edit the lammps.cfg to have the correct directories and
filenames. Often times the lammps.cfg does not require much
editing. Note the previously mentioned <lammps_include_dir> and
<lammps_library>.

.. code-block:: bash

   $ `python setup.py install`  

You now have lammps-python installed! You can easily check

.. code-block:: bash

   $ python
   Python 3.4.3+ (default, Oct 14 2015, 16:03:50)
   [GCC 5.2.1 20151010] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import lammps

Next see how to use lammps-python in the :doc:`tutorial`.

There are some common errors that should be checked before looking
at the mailing list.

.. code-block::
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

.. code-block::
   >>> import lammps
   
   ImportError: core.cpython.so undefined symbol *****

This error is my fault for improperly writting the setup.py install
file. First check that the <lammps_library> has no undefined symbols
(see warning above). Next run :command:`ldd -r core.cpython.so`. You
can easily find this library in the lammps directory when you build
lammps-python with the command :command:`python setup.py build_ext
-i`. There should be no undefined symbols. For a quick fix simply
modify the setup.py file such that it includes a shared library where
the symbol is defined (variable libraries). If you do not have the
expertise please submit an issue on the github page.

For any other errors PLEASE add an issue to the github page. I check
github often and really want to make this a long-term supported
addition to the LAMMPS community!
