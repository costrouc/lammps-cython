# lammps-python

A cython wrapper around lammps. Lammps is a great molecular dynamics
package that currently does not have a convenient way to run. The goal
of this project is to put an opinionated wrapper around LAMMPS (the
good parts) and allow the user to easily extend it's functionality in
python. The api should feel very similar to 
[HOOMD](https://codeblue.umich.edu/hoomd-blue/)  

# Features

 - Full MPI support  
 - Pythonic MD api inspired by
[HOOMD](https://codeblue.umich.edu/hoomd-blue/)  
 - Run lammps regularly (use "-i" instead of stdin) 
```python 
from lammps import Lammps
import sys
Lammps(args=sys.args)
```

# Install
First install all the dependencies. You must install an MPI
implementation separately.  
> `pip install -r requirements.txt`  

Edit the lammps.cfg to have the correct directories and
filenames. Often times the lammps.cfg does not require much editing.  
> `python setup.py install`

## Installing LAMMPS Ubuntu
I will attempt to keep this current (1/10/2016) on how to install in
Ubuntu. Should work similarly for Linux distributions. I am not
knowledgeable on how to install on OSX or Windows (maybe this will
change). Run these commands most likely as a super user. If you are
working on a cluster it is a good chance that you will not need to
install the dependencies.

```bash
sudo apt install mpich libfftw3-dev libpng-dev libjpeg-dev
```

Assuming the dependencies are installed lammps should be __very__ easy
to install.

### Building Library
Since we are using python to execute lammps we do not need to make the
LAMMPS executable.

```bash
cd <lammps_download_folder>/src
make mode=shlib ubuntu
```

Honestly this make configuration should work for most systems.

# Dependencies

- Some MPI implementation (preferably implementing the MPI3 api)  
- [mpi4py](https://bitbucket.org/mpi4py/mpi4py/)  
- [numpy](http://www.numpy.org/)  
- [cython](http://cython.org/)  
