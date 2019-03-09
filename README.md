# HPC19
Spring 2019: Advanced Topics in Numerical Analysis: High Performance Computing

## Running hw2 on CIMS server
```
module purge
module load gcc-8.1
cd hw2/
make
./MMult1
./val_test01_solved
./val_test02_solved
./omp_solved2
...
./omp_solved6
./jacobi2D-omp
./gs2D-omp
```