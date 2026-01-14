# Welcome to PacIFiC!

PacIFiC is a high-performance MPI parallel c/c++ software to compute particle-laden flows at the particle scale. PacIFiC stands for "PArtiCles In FluId Computations".

PacIFiC is open-sourced under the MIT license, and is developed by the research group of Prof. Anthony Wachs at the University of British Columbia, Vancouver, Canada with the support of IFP Energies nouvelles, France.

```mermaid
graph TD;
  pacific-->Octree;
  pacific-->Cartesian;
  pacific-->Grains3D;
  Octree-.->basilisk;
  Octree-->FDM1[DLMFD];
  Octree-->LBM;
  Octree-->MovingCutCell;
  Octree-->VTKHyperTree;
  Octree-->eulerian_caps;
  Octree-->lagrangian_caps
  Cartesian-->FLUID;
  FLUID-->FDM2[DLMFD];
  FLUID-->DirectionSplitting;
  Cartesian-->MacWorld;
  MacWorld-->MAC;
  MacWorld-.->HYPRE;
  MacWorld-.->PETSc;
```

Documentation will come soon.

## CMake Installation

### Requirements

For RPM-based distributions
```
  sudo dnf install 
```

### Build
```
mkdir -p build
cd build
cmake .. && make
```

