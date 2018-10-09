{ pkgs ? import <nixpkgs> {}, pythonPackages ? "python36Packages" }:

let
  elem = builtins.elem;
  basename = path: with pkgs.lib; last (splitString "/" path);
  startsWith = prefix: full: let
    actualPrefix = builtins.substring 0 (builtins.stringLength prefix) full;
  in actualPrefix == prefix;

  src-filter = path: type: with pkgs.lib;
    let
      ext = last (splitString "." path);
    in
      !elem (basename path) [".git" "__pycache__" ".eggs"] &&
      !elem ext ["egg-info" "pyc"] &&
      !startsWith "result" path;

   basePythonPackages = if builtins.isAttrs pythonPackages
     then pythonPackages
     else builtins.getAttr pythonPackages pkgs;

   mpi = basePythonPackages.mpi4py.mpi;
in
basePythonPackages.buildPythonPackage rec {
  pname = "lammps-cython";
  version = "0.5.8";
  disabled = (!basePythonPackages.isPy3k);

  src = builtins.filterSource src-filter ./.;

  buildInputs = with basePythonPackages; [ pytestrunner pkgs.lammps-mpi cython ];
  checkInputs = with basePythonPackages; [ pytest pytestcov pkgs.openssh pymatgen ase gsd ];
  propagatedBuildInputs = with basePythonPackages; [ numpy mpi4py ];

  preBuild = ''
    echo "Creating lammps.cfg file..."
    cat << EOF > lammps.cfg
    [lammps]
    lammps_include_dir = ${pkgs.lammps-mpi}/include
    lammps_library_dir = ${pkgs.lammps-mpi}/lib
    lammps_library = lammps_mpi

    [mpi]
    mpi_include_dir = ${mpi}/include
    mpi_library_dir = ${mpi}/lib
    mpi_library     = mpi
    EOF
  '';

  meta = with pkgs; {
    description = "Pythonic Wrapper to LAMMPS using cython";
    homepage = https://gitlab.com/costrouc/lammps-cython;
    license = lib.licenses.mit;
    maintainers = with lib.maintainers; [ costrouc ];
  };
}
