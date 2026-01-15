{ pkgs ? import <nixpkgs> {} }:

let
  lib = pkgs.lib;

  pname = "pacific";
  version = "0.1.0";

  package = pkgs.stdenv.mkDerivation {
    inherit pname version;
    src = ./.;

    nativeBuildInputs = with pkgs; [ 
      gcc
      cmake 
      gnumake
      mpi
      zlib
      xercesc
      hdf5-mpi
    ];

    buildInputs = with pkgs; [ 
    ];

    configurePhase = ''
      cmake -S . -B build -G Ninja -DCMAKE_INSTALL_PREFIX=$out
    '';

    buildPhase = ''cmake --build build'';
    
    installPhase = ''cmake --install build'';
  };

  shell = pkgs.mkShell {
    inputsFrom = [ package ];
    packages = with pkgs; [ 
      gdb 
      strace
      doxygen
    ];
  };

in
# nix-shell sets IN_NIX_SHELL, nix-build doesn't.
if lib.inNixShell then shell else package
