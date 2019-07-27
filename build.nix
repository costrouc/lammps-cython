{ pkgs ? import <nixpkgs> {}, pythonPackages ? "python3Packages" }:

rec {
  package = pythonPackages.buildPythonPackage rec {
    pname = "lammps-cython";
    version = "master";
    disabled = pythonPackages.isPy27;

    src = builtins.filterSource
      (path: _: !builtins.elem  (builtins.baseNameOf path) [".git" "result" "docs"])
      ./.;

    buildInputs = with pythonPackages; [
      pytestrunner
      pkgs.lammps-mpi
      cython
    ];

    checkInputs = with pythonPackages; [
      pytest
      pytestcov
      pkgs.openssh
      pymatgen
      ase
      gsd
     ];

    propagatedBuildInputs = with pythonPackages; [
      numpy
      mpi4py
    ];

    postConfigure = ''
      echo "Creating lammps.cfg file..."
      cat << EOF > lammps.cfg
      [lammps]
      lammps_include_dir = ${pkgs.lammps-mpi}/include
      lammps_library_dir = ${pkgs.lammps-mpi}/lib
      lammps_library = lammps_mpi

      [mpi]
      mpi_include_dir = ${pythonPackages.mpi4py.mpi}/include
      mpi_library_dir = ${pythonPackages.mpi4py.mpi}/lib
      mpi_library     = mpi
      EOF
    '';

    checkPhase = ''
      pytest --pkg-molecule
    '';
  };

  sdist = pkgs.stdenv.mkDerivation {
    name = "lammps-cython-sdist";

    buildInputs = [ package ];

    src = builtins.filterSource
        (path: _: !builtins.elem  (builtins.baseNameOf path) [".git" "result"])
        ./.;

    buildPhase = ''
      python setup.py sdist
    '';

    installPhase = ''
      mkdir -p $out
      cp dist/* $out
    '';
  };

  docs = pkgs.stdenv.mkDerivation {
    name = "docs";
    version = "master";

    src = builtins.filterSource
        (path: _: !builtins.elem  (builtins.baseNameOf path) [".git" "result"])
        ./.;

    buildInputs = with pythonPackages; [
      package
      pkgs.openssh
      sphinx
      sphinx_rtd_theme
    ];

    buildPhase = ''
      cd docs;
      sphinx-apidoc -o source/ ../lammps
      sphinx-build -b html -d build/doctrees . build/html
    '';

    installPhase = ''
     mkdir -p $out
     cp -r build/html/* $out
     touch $out/.nojekyll
    '';
  };

  docker = let
      pythonEnv = pythonPackages.python.withPackages
                   (ps: with ps; [ jupyterlab package ase pymatgen gsd ]);
    in pkgs.dockerTools.buildLayeredImage {
      name = "lammps-cython";
      tag = "latest";
      config.Cmd = [ "${pythonEnv.interpreter}" ];
      maxLayers = 120;
    };
}
