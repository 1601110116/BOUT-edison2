#include <iostream>
#include <fstream.h>
#include "mpi_header.hxx"

int MPI_init()
{
	int nproc;
	int myid;

	ofstream fout;
	
	MPI::Init();

	nproc = MPI::COMM_WORLD.Get_size();
	myid  = MPI:: COMM_WORLD.Get_rank();

	if(myid ==0)
	{
		fout.open("nproc.txt");
		fout << "nproc" << "\t" << nproc << endl;
	}	
return 0;

}

int MPI_get_rank()
{
	return MPI::COMM_WORLD.Get_rank();
}


int MPI_finalize()
{
	MPI::Finalize();
	return 0;
}
