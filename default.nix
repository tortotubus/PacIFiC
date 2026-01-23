{ pkgs ? import <nixpkgs> {} }:

let
  pname   = "pacific";
  version = "0.0.1";

  src = pkgs.fetchFromGitHub {
    owner = "tortotubus";
    repo  = "PacIFiC";
    fetchSubmodules = true;
  };

  package = pkgs.stdenv.mkDerivation {
    inherit pname version;
    src = ./.;

    nativeBuildInputs = with pkgs; [
      gcc
      cmake
      gnumake
      pkg-config
      (hdf5-mpi.override { mpi = openmpi; })
    ];

    buildInputs = with pkgs; [
      zlib
      xercesc
    ];

    cmakeFlags = [
      "-DCMAKE_BUILD_TYPE=Release"
      "-DUSE_SUBMODULES=ON"
      "-DOCTREE_BASILISK_PROVIDER=VENDORED"
      "-DCMAKE_INSTALL_PREFIX=$out"
    ];
  };

  devShell = pkgs.mkShell {
    inputsFrom = [ package ];
    packages = with pkgs; [
      gdb
      strace
      doxygen
      ninja
      apptainer
    ];
  };
in
{
  inherit package devShell;

  # Nice defaults for common commands:
  default = package;
}
