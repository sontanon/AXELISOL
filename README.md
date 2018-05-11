# `AXELISOL`
## Purpose.
`AXELISOL` (Axisymmetric Elliptic Solver) is an elliptic solver that can solve two types of elliptic equations in an axisymmetric space *with equatorial symmetry*.

* A flat Laplacian in axisymmetric space.
* A general elliptic equation with variable coefficients.

The solution is solved at either **2nd order** or **4th order** finite differences.

### A note on equatorial symmetry.
Please note that for the moment, `AXELISOL` only supports equatorial symmetry. In other words, all functions must have a definite parity not only about the ρ axis, *but also about the z axis*.

### Flat Laplacian.
A flat Laplacian takes the following form in axisymmetric space:

<a href="http://www.codecogs.com/eqnedit.php?latex=\left(\frac{\partial^2}{\partial\,\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial\,z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)&space;=&space;f(\rho,&space;z)\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?\left(\frac{\partial^2}{\partial\,\rho^2}&space;&plus;&space;\frac{\partial^2}{\partial\,z^2}&space;&plus;&space;\frac{1}{\rho}\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;s(\rho,&space;z)&space;\right)\,u(\rho,&space;z)&space;=&space;f(\rho,&space;z)\,," title="\left(\frac{\partial^2}{\partial\,\rho^2} + \frac{\partial^2}{\partial\,z^2} + \frac{1}{\rho}\,\frac{\partial}{\partial\,\rho} + s(\rho, z) \right)\,u(\rho, z) = f(\rho, z)\,," /></a>

where `s` is a linear source, `f` is a right hand side and `u` is the sought-after solution.

### General Elliptic Equation.
A general elliptic equation with variable coefficients has the form

<a href="http://www.codecogs.com/eqnedit.php?latex=\left(a\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;e\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s\right)\,u&space;=&space;f\,." target="_blank"><img src="http://latex.codecogs.com/gif.latex?\left(a\,\frac{\partial^2}{\partial&space;\rho^2}&space;&plus;&space;b\,\frac{\partial^2}{\partial&space;\rho\,\partial&space;z}&space;&plus;&space;c\,\frac{\partial^2}{\partial&space;z^2}&space;&plus;&space;d\,\frac{\partial}{\partial\,\rho}&space;&plus;&space;e\,\frac{\partial}{\partial&space;z}&space;&plus;&space;s\right)\,u&space;=&space;f\,." title="\left(a\,\frac{\partial^2}{\partial \rho^2} + b\,\frac{\partial^2}{\partial \rho\,\partial z} + c\,\frac{\partial^2}{\partial z^2} + d\,\frac{\partial}{\partial\,\rho} + e\,\frac{\partial}{\partial z} + s\right)\,u = f\,." /></a>

Now `a`, `b`, `c`, `d`, `e`, `f` and `s` are all functions of ρ, z, and is in the flat case, `u` is the solution.

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

A common way of knowing if MKL is properly installed is if `MKLROOT` is a valid enviroment variable.

## Building `AXELISOL`.
In the `src` directory there are two main files. One in C and another in FORTRAN. All other source files are in C, including the CSR matrix generators and the PARDISO wrappers and related subroutines.  To build `AXELISOL` you must first choose to build using a C main file and C style calls or, alternatively, use a FORTRAN main file and FORTRAN style subroutine calls.
To build using C, type:
```console
$ make C
```
To build using FORTRAN:
```console
$ make FORTRAN
```
By default, `make` will try to build using Intel’s compilers. However, if you wish to specify which type of compiler to utilize, pass the compiler option to `make`. For example, to build a C target with GNU use
```console
$ make C compiler=gnu
```
And to build FORTRAN with Intel:
```console
$ make FORTRAN compiler=intel
```
Type `make help` to see a summary of building options. Please note that `make` uses the preprocessor to generate FORTRAN subroutine headers. In other words, if you want FORTRAN headers make sure that you define the preprocessor macro `FORTRAN`. This is done in `make` using the `CFLAG` `-D FORTRAN`.

## Boundary Conditions.

`AXELISOL` uses a cartesian grid in ρ, z and thus requires four boundary conditions corresponding to the four edges of the grid.

### Symmetry about the ρ, z axes.

All functions in axisymmetric space must have a definite parity about the ρ axis and if the space has equatorial symmetry, also about the z axis (for the moment, equatorial symmetry is required).

### Robin boundary condition at grid edge.

On the external boundary, `AXELISOL` uses a Robin-type boundary condition:

<a href="http://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;r}&space;&plus;&space;\frac{(u&space;-&space;u_\infty)}{r}&space;=&space;0\,." target="_blank"><img src="http://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;r}&space;&plus;&space;\frac{(u&space;-&space;u_\infty)}{r}&space;=&space;0\,." title="\frac{\partial u}{\partial r} + \frac{(u - u_\infty)}{r} = 0\,." /></a>

Notice that r here stands for the radial coordinate <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;r&space;=&space;\sqrt{\rho^2&space;&plus;&space;z^2}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;r&space;=&space;\sqrt{\rho^2&space;&plus;&space;z^2}" title="r = \sqrt{\rho^2 + z^2}" /></a>. This condition is equivalent to having the solution `u` behave at large r like 

<a href="http://www.codecogs.com/eqnedit.php?latex=u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-2})\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-2})\,," title="u = u_\infty + \frac{M}{r} + \mathcal{O}(r^{-2})\,," /></a>

where <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;u_\infty" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;u_\infty" title="u_\infty" /></a> is the solution’s value at infinity and <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;M" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;M" title="M" /></a> is a constant.

`AXELISOL` can implement three types of Robin boundary conditions:

| Robin BC Type | Equation equivalent |
|:-------------:|:--------------------|
| 1 | <a href="http://www.codecogs.com/eqnedit.php?latex=u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-2})\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-2})\,," title="u = u_\infty + \frac{M}{r} + \mathcal{O}(r^{-2})\,," /></a> |
| 2 | <a href="http://www.codecogs.com/eqnedit.php?latex=u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-3})\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-3})\,," title="u = u_\infty + \frac{M}{r} + \mathcal{O}(r^{-3})\,," /></a> |
| 3 | <a href="http://www.codecogs.com/eqnedit.php?latex=u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-4})\,," target="_blank"><img src="http://latex.codecogs.com/gif.latex?u&space;=&space;u_\infty&space;&plus;&space;\frac{M}{r}&space;&plus;&space;\mathcal{O}(r^{-4})\,," title="u = u_\infty + \frac{M}{r} + \mathcal{O}(r^{-4})\,," /></a> |

## Grid and variables.
As already mentioned, this program uses a cartesian ρ, z grid. The grid structure is drawn in the follwing figure:

![Grid Strucutre](https://github.com/sontanon/AXELISOL/blob/master/figures/grid.png)

* The red points are the ghost zones used for symmetry about the ρ, z axes. If you wish to use second order finite differences, there must be at least one ghost zones. On the other hand, fourth order finite differences requires at least two ghost zones. The number of ghost zones is reffered as `ghost` in subroutine calls. In this example there are two ghost zones.
* The black points are the actual interior points of the grid and are called `NrInterior` and `NzInterior` corresponding to ρ, z respectively. The figure above hast 7 interior points in both directions.
* Finally, the blue points are the external boundary which are just a single strip in both directions.

Note that the step sizes in each direction Δρ, Δz (which are called `dr` and `dz` in the code) are uniform but can be different between themselves.

The following table contains a summmary of the grid variables 

| Grid variable        | Description           |
|:---------------------|:----------------------|
| `NrInterior`         | Number of interior points in ρ. | 
| `NzInterior`         | Number of interior points in z. |
| `ghost`              | Number of ghost zones. |
| `NrTotal`            | Total grid extension in ρ: `ghost` + `NrInterior` + 1. |
| `NzTotal`            | Total grid extension in z: `ghost` + `NzInterior` + 1. |
| `ARRAY_DIM`          | Total number of points: `NrTotal * NzTotal` |
| `dr`                 | Spatial step size in ρ. |
| `dz`                 | Spatial step size in z. |

### Memory access.
All grid functions (such as the solution `u`, RHS `f`, and various coefficients) have the same geoemtric structure. However, all functions are stored in **linear memory** where data is ρ-major ordered, i.e. z is the fast index.

## Usage.
`AXELISOL` has four fundamental subroutines.

### `pardiso_start`
Before anything PARDISO must be initialized with some default parameters and memory allocation. As such, be sure to call from C 

```C
pardiso_start(NrInterior, NzInterior);
```
Or from FORTRAN
```FORTRAN
CALL PARDISO_START(NRINTERIOR, NZINTERIOR)
```

### `pardiso_stop`
By the same token, at the end of the program or when no more elliptic equations are to be solved, this subroutine will clear all memory structures associated to PARDISO and the solver.

```C
pardiso_stop();
```
```FORTRAN
CALL PARDISO_STOP()
```

### `flat_laplacian`
This is the first of the elliptic solver subroutines. The argument list is as follows (for C and FORTRAN):

```C
flat_laplacian(u, res, s, f, 
                u_inf, robin, r_sym, z_sym, 
                NrInterior, NzInterior, ghost, dr, dz, order, 
                lr_use, precond_use);
```
```FORTRAN
CALL FLAT_LAPLACIAN(U, RES, S, F,&
                    U_INF, ROBIN, R_SYM, Z_SYM,& 
                    NRINTERIOR, NZINTERIOR, GHOST, DR, DZ, ORDER,& 
                    LR_USE, PRECOND_USE)
```

The variables are described in the following table:

| Argument      | Input/Ouput | Type | Description           | README section |
|:--------------|:-----------:|:----:|:----------------------|----------------|
| `u`           | Output | Double precision array of size `ARRAY_DIM` | Solution. | Flat Laplacian |
| `res`         | Output | Double precision array of size `ARRAY_DIM` | Residual. | Flat Laplacian |
| `s`           | Input  | Double precision array of size `ARRAY_DIM` | Linear source. | Flat Laplacian |
| `f`           | Input  | Double precision array of size `ARRAY_DIM` | Right-hand side.| Flat Laplacian |
| `u_inf`       | Input  | Double  | Solution value at infinity for Robin BC. | Boundary Conditions |
| `robin`       | Input  | Integer | Robin BC type: 1, 2, or 3. | Boundary Conditions |
| `r_sym`       | Input  | Integer | Parity in ρ axis: -1(odd), +1(even) | Boundary Conditions |
| `z_sym`       | Input  | Integer | Parity in z axis: -1(odd), +1(even) | Boundary Conditions |
| `NrInterior`  | Input  | Integer | Number of interior points in ρ. | Grid and variables |
| `NzInterior`  | Input  | Integer | Number of interior points in z. | Grid and variables |
| `ghost`       | Input  | Integer | Number of ghost zones. | Grid and variables |
| `dr`          | Input  | Double  | Spatial step size in ρ. | Grid and variables |
| `dz`          | Input  | Double  | Spatial step size in z. | Grid and variables |
| `order`       | Input  | Integer | Finite difference order: 2 or 4 | Grid and variables |
| `lr_use`      | Input  | Integer | Use low rank update: 1(on), 0(off) | Low Rank Update and Preconditioning |
| `precond_use` | Input  | Integer | Use CGS preconditioner: 1(on), 0(off) | Low Rank Update and Preconditioning |

### `general_elliptic`
This subroutine solve the general linear elliptic equation. The argument list is as follows (for C and FORTRAN):

```C
general_elliptic(u, res, a, b, c, d, e, s, f, 
                u_inf, robin, r_sym, z_sym, 
                NrInterior, NzInterior, ghost, dr, dz, order, 
                lr_use, precond_use);
```
```FORTRAN
CALL GENERAL_ELLITPTIC(U, RES, A, B, C, D, E, S, F,&
                    U_INF, ROBIN, R_SYM, Z_SYM,& 
                    NRINTERIOR, NZINTERIOR, GHOST, DR, DZ, ORDER,& 
                    LR_USE, PRECOND_USE)
```
Most arguments are the sames as in the flat case, the only difference lies in the variable coefficients.

| Argument      | Input/Ouput | Type | Description           | README section |
|:--------------|:-----------:|:----:|:----------------------|----------------|
| `a`           | Input  | Double precision array of size `ARRAY_DIM` | <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial^2_{\rho&space;\rho}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;\partial^2_{\rho&space;\rho}" title="\partial^2_{\rho \rho}" /></a> coefficient. | General Elliptic Equation |
| `b`           | Input  | Double precision array of size `ARRAY_DIM` | <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial^2_{\rho&space;z}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;\partial^2_{\rho&space;z}" title="\partial^2_{\rho z}" /></a> coefficient. | General Elliptic Equation |
| `c`           | Input  | Double precision array of size `ARRAY_DIM` | <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial^2_{z&space;z}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;\partial^2_{z&space;z}" title="\partial^2_{z z}" /></a> coefficient. | General Elliptic Equation |
| `d`           | Input  | Double precision array of size `ARRAY_DIM` | <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial_{\rho}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;\partial_{\rho}" title="\partial_{\rho}" /></a> coefficient. | General Elliptic Equation |
| `e`           | Input  | Double precision array of size `ARRAY_DIM` | <a href="http://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial_{z}" target="_blank"><img src="http://latex.codecogs.com/gif.latex?\inline&space;\partial_{z}" title="\partial_{z}" /></a> coefficient. | General Elliptic Equation |
| `s`           | Input  | Double precision array of size `ARRAY_DIM` | Linear source. | General Elliptic Equation |
| `f`           | Input  | Double precision array of size `ARRAY_DIM` | Right-hand side.| General Elliptic Equation |

## Low Rank Update and Preconditioning



