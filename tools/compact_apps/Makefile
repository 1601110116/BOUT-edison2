#CC=mpiCC
CC=CC
RM=rm -f

default: tdma

tdma: tdma_parallel.o mpi.o initialize.o assemble_matrix.o \
	print.o forward_sweep.o back_sweep.o destroy.o dump_solution.o 
	$(CC) -o tdma tdma_parallel.o mpi.o initialize.o assemble_matrix.o \
	 print.o forward_sweep.o back_sweep.o destroy.o dump_solution.o $(IPM)

tdma_parallel.o: tdma_parallel.cxx 
	$(CC) -c tdma_parallel.cxx 

mpi.o : mpi.cxx
	$(CC) -c mpi.cxx

initialize.o : initialize.cxx
	$(CC) -c initialize.cxx

assemble_matrix.o: assemble_matrix.cxx
	$(CC) -c assemble_matrix.cxx 

forward_sweep.o: forward_sweep.cxx
	$(CC) -c forward_sweep.cxx

back_sweep.o: back_sweep.cxx
	$(CC) -c back_sweep.cxx

print.o : print.cxx
	$(CC) -c print.cxx

dump_solution.o : dump_solution.cxx
	$(CC) -c dump_solution.cxx

destroy.o : destroy.cxx
	$(CC) -c destroy.cxx

clean: *.o tdma
	$(RM) *.o tdma
