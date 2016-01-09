# lammps-python

A cython wrapper around lammps. Lammps is a great molecular dynamics
package that currently does not have a convenient way to run. The goal
of this project is to put an opinionated wrapper around LAMMPS (the
good parts) and allow the user to easily extend it's functionality in
python.

# Features

 - Full MPI support  
 - Pythonic MD api inspired by
[HOOMD](https://codeblue.umich.edu/hoomd-blue/)  
 - 
 - Run lammps regularly (use "-i" instead of stdin) 
```python 
from lammps import Lammps
import sys
Lammps(args=sys.args)
```

# Install
First install all the dependencies. You must install an MPI
implementation separately.
`pip install -r requirements.txt`

python setup.py install

# Dependencies

- Some MPI implementation (preferably implementing the MPI3 api)  
- [mpi4py](https://bitbucket.org/mpi4py/mpi4py/)  
- [numpy](http://www.numpy.org/)  
- [cython](http://cython.org/)  
