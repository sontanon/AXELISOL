# Elliptic-Solver
## Purpose
This elliptic solver can solve two types of equations:
1. A flat laplacian in axisymmetric space, i.e. the linear elliptic equation

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(\mathring{\nabla}^2\,&plus;s(\rho,&space;z)\right)\,u(\rho,z)&space;=&space;f(\rho,&space;z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(\mathring{\nabla}^2\,&plus;s(\rho,&space;z)\right)\,u(\rho,z)&space;=&space;f(\rho,&space;z)" title="\left(\mathring{\nabla}^2\,+s(\rho, z)\right)\,u(\rho,z) = f(\rho, z)" /></a>

where

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathring{\nabla}^2&space;=&space;\frac{\partial^2}{\partial\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial&space;\rho}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathring{\nabla}^2&space;=&space;\frac{\partial^2}{\partial\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial&space;\rho}" title="\mathring{\nabla}^2 = \frac{\partial^2}{\partial\rho^2} + \frac{\partial^2}{\partial z^2} + \frac{1}{\rho}\,\frac{\partial}{\partial \rho}" /></a>

is the flat laplacian in cyllindrical coordinates.

2. A general linear elliptic equation with variable coefficients
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;\left(a(\rho,z)\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b(\rho,&space;z)\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c(\rho,&space;z)\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d(\rho,&space;z)\,\frac{\partial}{\partial&space;\rho}&space;&plus;&space;e(\rho,&space;z)\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)=f(\rho,&space;z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;\left(a(\rho,z)\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b(\rho,&space;z)\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c(\rho,&space;z)\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d(\rho,&space;z)\,\frac{\partial}{\partial&space;\rho}&space;&plus;&space;e(\rho,&space;z)\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)=f(\rho,&space;z)" title="\left(a(\rho,z)\,\frac{\partial^2}{\partial \rho^2} + b(\rho, z)\,\frac{\partial^2}{\partial \rho\,\partial z} + c(\rho, z)\,\frac{\partial^2}{\partial z^2} + d(\rho, z)\,\frac{\partial}{\partial \rho} + e(\rho, z)\,\frac{\partial}{\partial z} + s(\rho, z) \right)\,u(\rho, z)=f(\rho, z)" /></a>

Both equations are solved in an axisymmetric space which also has equatorial symmetry.

## Sparse Direct Solver
This solver discretizes the partial differential equation via finite differences into a linear system of equations. Such as system results to be sparse and can thus be solved with a specialized sparse solver. As such, **PARDISO** is the chosen direct solver. Throughout this project, two different versions of PARDISO may be utilized

1. The oficial and up-to-date version available in the [PARDISO project page](https://pardiso-project.org/). This version will be referenced as `VANILLA PARDISO`. For more information, please consult the [User Guide](https://pardiso-project.org/manual/manual.pdf).
2. INTEL MKL's version optimized for INTEL processors and INTEL compilers. This version, which will be referenced as `INTEL PARDISO`, can be aquired by installing [MKL](https://software.intel.com/en-us/mkl) and it has its own [User Guide](https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface).

For more information on how these two version differ in their calls, see the *Usage* section.

## Prerequisites
To compile `ELLSOLVE` you must have two main components:
1. **INTEL compilers `icc`, `ifort` and an installation of MKL** for SPARSE BLAS. 
  All compiler and enviroment variables must be properly setup. After installing the compilers and MKL, this can usually be done by including the following lines in `~/.bashrc`.
```  
export INTEL=/opt/intel
export PATH=$INTEL/bin:$PATH
export LD_LIBRARY_PATH=$INTEL/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$INTEL/lib:$LIBRARY_PATH
source $INTEL/bin/ifortvars.sh intel64
source $INTEL/bin/iccvars.sh intel64
source $INTEL/bin/compilervars.sh intel64
source $INTEL/mkl/bin/mklvars.sh intel64
```
2. **PARDSIO**. INTEL MKL ships with its own version of PARDISO, however you can choose to use the ofical version in [here](https://pardiso-project.org/). If you choose this version, be sure to aquire the version compatible with INTEL's compilers (libpardiso500-INTEL1301-X86-64.so) and to add it to `LD_LIBRARY_PATH` and `LIBRARY_PATH`. For example, if you install the shared library `.so` into `~/PARDISO`, everything is taken care by including these lines inside your `~/.bashrc`.
```
export PARDISO=$HOME/PARDISO
export PARDISO_LIC_PATH=$PARDISO
export LD_LIBRARY_PATH=$PARDISO:$LD_LIBRARY_PATH
export LIBRARY_PATH=$PARDISO:$LIBRARY_PATH
```

## Usage
`ELLSOLVE` has three fundamental subroutines.

* `flat_laplacian`
* `pardiso_start`
* `pardiso_stop`

## Compilation
Assuming you have all the prerequisites, compilation is as simple as choosing one of the four versions and simply calling `make` which will create a `bin` directory for the object files and produce the executable.
