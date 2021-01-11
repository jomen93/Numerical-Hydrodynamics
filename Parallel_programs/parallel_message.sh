mpicc mensajes_paralelo.cpp -O3 -o message
mpiexec --oversubscribe -np 4 ./message

rm message