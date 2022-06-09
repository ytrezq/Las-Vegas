# ECDLP-LasVegas
## Oblique Elimination

### Introduction
This folder contains code for the paper titled **Oblique Elimination to solve the elliptic cirve discrete logarithm problem**.

Algorithm 2 from this paper is implemented in this code base.
The code can be compiled using the make file. 
An executable named lasVeagas is geneted after successful compilation.

### Compiling/Running Code.
	1. Compiling : make 
	2. Running : mpirun -n 2 lasVegas


The following dependencies are required for the code to work.
### Dependencies:
	1. C/C++ compiler
	2. OpemMPI
	3. Number Theory Library (NTL) and GNU-GMP
	4. Sage (genetate input)
