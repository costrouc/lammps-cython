# lammps-python

A cython wrapper around lammps. Lammps is a great molecular dynamics
package that currently does not have a convenient way to run. The goal
of this project is to put an opinionated wrapper around LAMMPS (the
good parts) and allow the user to easily extend it's functionality in
python.

## Features

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


## Dependencies

- Some MPI implementation (preferably implementing the MPI3 api)  
- [mpi4py](https://bitbucket.org/mpi4py/mpi4py/)  
- [numpy](http://www.numpy.org/)  
- [cython](http://cython.org/)  


# Documentation for now
## Creation of Box and Atoms
This is a step that I was very confused with at first. It is not the
#most elegant method but again it the lammps code does not make it easy.

I will use the method

`region <region_id> prism xlo xhi ylo yhi zlo zhi xy xz yz`
`create_box <num_atom_types> <region_id>`

Repeat for as many atoms as there are
`create_atoms <num_atom_types> single`
