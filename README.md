# MPI_sorting
Sort key value pair using MPI in cpp <br />
The objective of this experiment is to sort a structure by usingmultiple processes using MPI. We are
using OpenMPI with C++ <br/>
We have implemented library sort.h with functions for merge sort, quick sort, radix sort, and heap
sort. Each process with its own data calls the sort function and selects the numbers which need to be
kept in memory. In this experiment, we are not using any additional memory and taken minimum
communication with other processes which has a great impact on speed and scalability.
## setup
### Create executable file using:
1. mpiCC -c -o sort.o sort.cpp
2. ar -rsclibpsort.a sort.o
### Run the file using:
mpirun -n 64 a.out
