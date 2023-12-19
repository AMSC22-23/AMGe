# AMGe
Multigrid implementation for the Advanced Modeling for Scientific Computing @polimi

# Description of the Problem
Given the Poisson equation expressed in a domanin  Ω ⊂ R^n with a boundary ∂Ω as

             −∇^2φ = f in Ω

              φ = g on ∂Ω (Dirichlet Condition)

where ∇^2 is the Laplacian operator, φ = φ(x) is our scalar variable that is a function of space, f = f(x ) is
a forcing function and g = g(x) is a prescribed function along the boundary.

Ω is discretized using a 2nd-order finite difference approximation on a cartesian mesh composed of a number N nodes in the x- and y-directions with uniform spacing h. The result is a system of (N − 2) x (N − 2) linear equations for the unknown values of φi,j in the interior of the domain. Given a natural ordering of the unknowns, it is possible to express the linear system:

               A Φ = f
where the matrix A is banded as a result of the structured mesh discretization.

A multigrid (MG) method is an iterative algorithm of the form x(k+1) = MG(x(k)), k ≥ 0, for solving the (typically) sparse linear systems of equations like the previous one.
MG methods are based on a hierarchy of levels (associated with a hierarchy of discretisations).
MG-cycle reduces all error components by a fixed amount (bounded well below
one) independent of the dimension n of the system.
The main idea of multigrid is to accelerate the convergence of a basic iterative
method by a global correction of the fine grid solution approximation
accomplished by solving a coarse problem.




# Overview
The main purpose of this project is studying the convergence and scalability of Multigrid poisson equation solver. 


* The smoothers used are Jacobi and Gauss-Seidel. 

* We use tests with different number of iterations, different size of the mesh problem and different levels of multigrid.

* Jacobi is parallelized with OpenMP.

* The method is only applicable to structured square meshes having a number of nodes in each dimension equal to 2^n + 1.

# Build with Make
The instructions for using it:

## Prerequisities
* Make
* OpenMP
* Google Benchmark

## Build the project
```
cd build
make
```
## Run the executable
```

```


