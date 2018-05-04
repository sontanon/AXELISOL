# Elliptic-Solver
## Purpose
Solves a flat laplacian in axisymmetric space. Specifically, the linear elliptic equation

<a href="https://www.codecogs.com/eqnedit.php?latex=\left(\mathring{\nabla}^2\,&plus;s(\rho,&space;z)\right)\,u(\rho,z)&space;=&space;f(\rho,&space;z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left(\mathring{\nabla}^2\,&plus;s(\rho,&space;z)\right)\,u(\rho,z)&space;=&space;f(\rho,&space;z)" title="\left(\mathring{\nabla}^2\,+s(\rho, z)\right)\,u(\rho,z) = f(\rho, z)" /></a>

where

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathring{\nabla}^2&space;=&space;\frac{\partial^2}{\partial\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial&space;\rho}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathring{\nabla}^2&space;=&space;\frac{\partial^2}{\partial\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial&space;\rho}" title="\mathring{\nabla}^2 = \frac{\partial^2}{\partial\rho^2} + \frac{\partial^2}{\partial z^2} + \frac{1}{\rho}\,\frac{\partial}{\partial \rho}" /></a>

is the flat laplacian in cyllindrical coordinates.

## Sparse Direct Solver
This solver utilizes [PARDISO](https://pardiso-project.org/) as the main solver. Although, this solver does considerable effort to simplify how the user calls it, please consult PARDISO's User Guide for more information.

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
2. **PARDSIO**. INTEL MKL ships with its own version of PARDISO, however you can choose to use the separtedly distributted version in [here](https://pardiso-project.org/). If you choose this version, be sure to aquire the INTEL version and to add it to your `LIBRARY_PATH`. For example, if you install the shared libraries `.so` into a directory inside your home directory, everything is taken care by including these lines inside your `~/.bashrc`.
```
export PARDISO=$HOME/PARDISO
export PARDISO_LIC_PATH=$PARDISO
export LD_LIBRARY_PATH=$PARDISO:$LD_LIBRARY_PATH
export LIBRARY_PATH=$PARDISO:$LIBRARY_PATH
```
`ELLSOLVE` has four working version. Two are written in C and two in FORTRAN. Each of these two versions corresponds to using INTEL MKL's PARDISO or the official version. 

## Usage
`ELLSOLVE` has three fundamental subroutines.

* `flat_laplacian`
* `pardiso_start`
* `pardiso_stop`
