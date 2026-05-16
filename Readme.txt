make -f Makefile_PC clean

make -f Makefile_PC
or 
make -f Makefile_PC all

mpirun -np <number_of_processes> ./grit.x

mpirun -np 1 ./grit.x
mpirun -np 2 ./grit.x