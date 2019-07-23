{ nixpkgs ? <nixpkgs> }:

let pkgs = import nixpkgs { };
    build = import ./build.nix {
      inherit pkgs;
      pythonPackages = pkgs.python3Packages;
    };
in {
  lammps-cython = build.package;
  lammps-cython-docs = build.docs;
  lammps-cython-docker = build.docker;
}
