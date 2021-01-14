mpicc mensajes_paralelo.cpp -O3 -o message
mpiexec --oversubscribe -np 6 ./message

rm message