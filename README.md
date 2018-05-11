# `AXELISOL`
## Purpose.
`AXELISOL` (Axisymmetric Elliptic Solver) is an elliptic solver that can solve two types of elliptic equations in an axisymmetric space.
* A flat Laplacian in axisymmetric space.
* A general elliptic equation with variable coefficients.

### Flat Laplacian.
A flat Laplacian takes the following form in axisymmetric space:

<a href="http://www.codecogs.com/eqnedit.php?latex=\left(\frac{\partial^2}{\partial\,\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial\,z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)&space;=&space;f(\rho,&space;z)\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?\left(\frac{\partial^2}{\partial\,\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial\,z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)&space;=&space;f(\rho,&space;z)\,," title="\left(\frac{\partial^2}{\partial\,\rho^2} + \frac{\partial^2}{\partial\,z^2} + \frac{1}{\rho}\,\frac{\partial}{\partial\,\rho} + s(\rho, z) \right)\,u(\rho, z) = f(\rho, z)\,," /></a>

where `s` is a linear source, `f` is a right hand side and `u` is sought-after solution.

### General Elliptic Equation.
A general elliptic equation with variable coefficients has the form

<a href="http://www.codecogs.com/eqnedit.php?latex=\left(a\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;e\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s\right)\,u&space;=&space;f\,." target="_blank"><img src="http://latex.codecogs.com/gif.latex?\left(a\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;e\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s\right)\,u&space;=&space;f\,." title="\left(a\,\frac{\partial^2}{\partial \rho^2} + b\,\frac{\partial^2}{\partial \rho\,\partial z} + c\,\frac{\partial^2}{\partial z^2} + d\,\frac{\partial}{\partial\,\rho} + e\,\frac{\partial}{\partial z} + s\right)\,u = f\,." /></a>

Now `a`, `b`, `c`, `d`, `e`, `f` and `s` are all functions of ρ, z, and is in the flat case, `u` is the solution.

## Boundary Conditions.

`AXELISOL` uses a cartesian grid in ρ, z and thus requires four boundary conditions corresponding to the four edges of the grid.

* Symmetry in the ρ, z axes.

All functions in axisymmetric space must have a definite parity about the ρ axis and if the space has equatorial symmetry, also about the z axis. 

* Robin boundary condition at grid edge.

On the external boundary, `AXELISOL` uses a Robin-type boundary condition:

<a href="http://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;r}&space;&plus;&space;\frac{(u&space;-&space;u_\infty)}{r}&space;=&space;0\,." target="_blank"><img src="http://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;r}&space;&plus;&space;\frac{(u&space;-&space;u_\infty)}{r}&space;=&space;0\,." title="\frac{\partial u}{\partial r} + \frac{(u - u_\infty)}{r} = 0\,." /></a>

Notice that r here stands for the radial coordinate <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;r&space;=&space;\sqrt{\rho^2&space;&plus;&space;z^2}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;r&space;=&space;\sqrt{\rho^2&space;&plus;&space;z^2}" title="r = \sqrt{\rho^2 + z^2}" /></a>. This condition is equivalent to having a solution `u` behave at large r like 

<a href="http://www.codecogs.com/eqnedit.php?latex=u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}\,," title="u = u_\infty + \frac{M}{r}\,," /></a>

where <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;u_\infty" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;u_\infty" title="u_\infty" /></a> the solution’s value at infinity and <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;M" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;M" title="M" /></a> is a constant .

## Sparse Direct Solver.
This solver uses finite differences to discretize the partial differential equation into a linear system of equations. Such a system is sparse and can thus be solved with a specialized solver. **PARDISO** has been chosen to fulfill this purpose. 
There are at least two versions of PARDISO. The first is the official version available at PARDISO’s [project page](https://pardiso-project.org/). However, this program utilizes [Intel MKL’s version](https://software.intel.com/en-us/mkl-developer-reference-c-intel-mkl-pardiso-parallel-direct-sparse-solver-interface).

## Prerequisites.
As already mentioned, this program requires Intel MKL’s version of PARDISO. However, a full installation of MKL is also necessary. MKL can be downloaded for free in the following [link](https://software.intel.com/en-us/mkl).
To build this application, the user may choose to use Intel’s proprietary compilers `icc` and `ifort` or the GNU alternatives `gcc` and `gfortran`.
In summary, to compile `AXELISOL` you need.
* Intel MKL: For PARDISO, BLAS and Sparse BLAS.
* A C compiler such as `gcc` or `icc`.
* (Optional) A FORTRAN compiler `gfortran` or `ifort` to build a FORTRAN-based executable.

A common way of knowing if MKL is properly installed is if `MKLROOT` is valid enviroment variable.

## Building `AXELISOL`.
In the `src` directory there are two main files. One in C and another in FORTRAN. All other source files are in C, including the CSR matrix generators and the PARDISO wrappers and related subroutines.  To build `AXELISOL` you must first choose to build using a C main file and C style calls or, alternatively, use a FORTRAN main file and FORTRAN style subroutine calls.
To build using C, type:
```
make C
```
To build using FORTRAN:
```
make FORTRAN
```
By default, `make` will try to build using Intel’s compilers. However, if you wish to specify which type of compiler to utilize, pass the compiler option to `make`. For example, to build a C target with GNU use
```
make C compiler=gnu
```
And to build FORTRAN with Intel:
```
make FORTRAN compiler=intel
```
Type `make help` to see a summary of building options.

## Usage.
`ELLSOLVE` has four fundamental subroutines.

* `flat_laplacian`
* `pardiso_start`
* `pardiso_stop`

## Preconditioning



