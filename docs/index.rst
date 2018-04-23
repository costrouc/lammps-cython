Welcome to LAMMPS Python Wrapper
================================

Welcome to the LAMMPS's Python Wrapper documentation. Being a
relatively new (1 week old) project the documentation and api are in
flux. To get started I recommend you checkout the :doc:`installation`
instructions and then walk through the :doc:`tutorial`. Full reference
documentation is available in the :doc:`api` section. The pythonic api
heavily borrows from `HOOMD <https://codeblue.umich.edu/hoomd-blue/>`_
but hopefully the community settles on a universal api for molecular
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

Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   tutorial
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

