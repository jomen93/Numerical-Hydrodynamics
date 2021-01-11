mpicc hello_mpi.cpp -O3 -o hello
mpiexec --oversubscribe -np 5 ./hello

rm hello