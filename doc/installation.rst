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
