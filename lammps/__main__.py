from lammps import Lammps
import sys


def main():
    lmp = Lammps(args=sys.argv[1:])


if __name__ == "__main__":
    main()
