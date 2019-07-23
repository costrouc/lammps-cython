{ pkgs ? import <nixpkgs> {}, pythonPackages ? "python3Packages" }:

{
  package = pythonPackages.buildPythonPackage rec {
    pname = "lammps-cython";
    version = "master";
    disabled = pythonPackages.isPy27;

    src = ./.;

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

  docs = pkgs.stdenv.mkDerivation {
    name = "docs";
    version = "master";

    src = ./.;

    buildInputs = with pythonPackages; [
      lammps-cython
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
    '';
  };

  docker = let
      pythonEnv = pythonPackages.python.withPackages
                   (ps: with ps; [ jupyterlab lammps-cython ase pymatgen gsd ]);
    in pkgs.dockerTools.buildLayeredImage {
      name = "lammps-cython";
      tag = "latest";
      config.Cmd = [ "${pythonEnv.interpreter}" ];
      maxLayers = 120;
    };
}
