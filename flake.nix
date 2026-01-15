{
  description = "PacIFiC (pacific) build + dev shell";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    beam-basilisk-src.url = "github:tortotubus/beam-basilisk";
  };

  outputs = { self, nixpkgs, flake-utils, beam-basilisk-src }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        pname = "pacific";
        version = "0.1.0";

        pacific = pkgs.stdenv.mkDerivation {
          inherit pname version;
          src = self;

          nativeBuildInputs = with pkgs; [
            cmake
            gnumake
          ];

          buildInputs = with pkgs; [
            gcc
            mpi
            xercesc
            hdf5-mpi
            zlib
          ];

          cmakeGenerator = "Unix Makefiles";
          cmakeBuildType = "Release";

          cmakeFlags = [
            "-DCMAKE_INSTALL_PREFIX=${placeholder "out"}"
            "-DBASILISK_SRC=${beam-basilisk-src}"
          ];

          # If your CMake really needs MPI wrapper compilers:
          # preConfigure = ''
          #   export CC=mpicc
          #   export CXX=mpicxx
          # '';
        };
      in {
        packages.default = pacific;
        packages.pacific = pacific;

        devShells.default = pkgs.mkShell {
          inputsFrom = [ pacific ];
          packages = with pkgs; [
            gdb
            strace
            doxygen
          ];
        };
      }
    );
}
