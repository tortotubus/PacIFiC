# Grains3D 

## Introduction

Grains3D is a high-performance DEM library. If this project is being built as a part of PacIFiC project, please refer to the root-level README.md. However, if you wish to build the Grains3D project alone, you may follow the instructions provided here.

## Requirements

The following depenencies are required:

 * OpenMPI
 * Xerces-C++ 
 * zlib

On rpm-based distributions, these can be obtained using 

```bash
sudo dnf install -y openmpi-devel xerces-c-devel zlib-ng-devel
```

## Building 

The Grains library requires the CMake meta-build tool and a build tool (e.g. ninja or gnumake). We will use ninja here. 

First, create a build directory where you would like to build the system. Here, we will create a build folder in the Grains3D root folder: To do so, configure and generate the build system by running

```bash
cmake -S . -B build -G "Ninja" # Alternatively, "Unix Makefiles"
```

To actually build the grains library, finally run

```bash
cmake --build build 
```

## Install

If you wish to install to your system, run

```bash
sudo cmake --build build --target install
```