{
  description = "PacIFiC (source-tree)";

  # Optional: helps CI/users if they haven't enabled flakes globally
  nixConfig = {
    experimental-features = [
      "nix-command"
      "flakes"
    ];
  };

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
      cudaToolchainFor = pkgs:
        let
          cudaPackages =
            if pkgs.stdenv.hostPlatform.isLinux && pkgs ? cudaPackages then
              pkgs.cudaPackages
            else
              null;
        in
        {
          cudaToolkit =
            if cudaPackages != null && cudaPackages ? cudatoolkit then
              cudaPackages.cudatoolkit
            else
              null;
          cudaNvcc =
            if cudaPackages != null && cudaPackages ? cuda_nvcc then
              cudaPackages.cuda_nvcc
            else
              null;
          # nixpkgs selects this stdenv for the CUDA package set. Using it for PacIFiC as well
          # keeps C, C++, and nvcc's host compiler in one supported compiler family.
          stdenv =
            if cudaPackages != null && cudaPackages ? backendStdenv then
              cudaPackages.backendStdenv
            else
              pkgs.stdenv;
        };
    in
    {
      overlays.default =
        final: prev:
        let
          toolchain = cudaToolchainFor final;
        in
        {
          pacific = final.callPackage ./pacific.nix {
            stdenv = toolchain.stdenv;
            cudaToolkit = toolchain.cudaToolkit;
            cudaNvcc = toolchain.cudaNvcc;
          };
        };

      packages = forAllSystems (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
            config.allowUnfree = true;
          };
        in
        {
          pacific = pkgs.pacific;
          default = pkgs.pacific;
        }
      );

      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
            config.allowUnfree = true;
          };
          toolchain = cudaToolchainFor pkgs;
          hostCC = toolchain.stdenv.cc;
        in
        {
          default = (pkgs.mkShell.override { stdenv = toolchain.stdenv; }) {
            inputsFrom = [ pkgs.pacific ];
            packages =
              (with pkgs; [
                gdb
                petsc
                doxygen
                python313Packages.mkdocs
                apptainer
              ])
              ++ pkgs.lib.optionals (toolchain.cudaToolkit != null) [ toolchain.cudaToolkit ]
              ++ pkgs.lib.optionals (toolchain.cudaNvcc != null) [ toolchain.cudaNvcc ];
            shellHook = ''
              export CC="${hostCC}/bin/gcc"
              export CXX="${hostCC}/bin/g++"
              export PETSC_DIR="${pkgs.petsc}"
              export HYPRE_DIR="${pkgs.hypre}"
              # Avoid relocatable Eigen configs reached through the system profile,
              # which may expose share/eigen3 without its corresponding headers.
              export Eigen3_ROOT="${pkgs.eigen}"
            ''
            + pkgs.lib.optionalString (toolchain.cudaToolkit != null) ''
              export CUDA_PATH="${toolchain.cudaToolkit}"
              export CUDA_HOME="${toolchain.cudaToolkit}"
            ''
            + pkgs.lib.optionalString (toolchain.cudaNvcc != null) ''
              export CUDACXX="${toolchain.cudaNvcc}/bin/nvcc"
              export CUDAHOSTCXX="${hostCC}/bin/g++"
            '';
          };
        }
      );

      # Nice: `nix flake check` builds it
      checks = forAllSystems (system: {
        pacific = self.packages.${system}.pacific;
      });
    };
}
