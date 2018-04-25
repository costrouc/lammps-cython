.. lammps-cython documentation master file, created by
   sphinx-quickstart on Wed Apr 25 12:10:05 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to lammps-cython's documentation!
=========================================

Welcome to the lammps-cython wrapper documentation. To get started I
recommend you checkout the :doc:`installation` instructions and then
walk through the :doc:`tutorial`. The installation can be as easy as
``docker pull`` of ``conda install``. Full reference documentation is
available in the :doc:`modindex` section. The pythonic api heavily
borrows from `HOOMD <https://codeblue.umich.edu/hoomd-blue/>`_ but
hopefully the community settles on a universal api for molecular
dynamics.

 - Full MPI support
 - Pythonic API inspired by `HOOMD <https://codeblue.umich.edu/hoomd-blue/>`_
 - Supports Python 2 and 3
 - Heavily documented and tested
 - Elimination of unnecessary file I/O

.. code-block:: python

   from lammps import Lammps
   import sys
   Lammps(args=sys.args)


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
