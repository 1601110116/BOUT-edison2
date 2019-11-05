#include "tdma.hxx"
#include <math.h>
#include <iostream>
#include <mpi.h>
int TDMA::assemble_matrix(TDMA *tdma)
{

int i;
int myrank;
int nproc;

	myrank = MPI::COMM_WORLD.Get_rank();
	nproc  = MPI::COMM_WORLD.Get_size();

	for (i=0;i<=tdma->N-1;i++)
		{
		  tdma->f[i] = -1*tdma->h*tdma->h;
		  tdma->a[i] =  1;
		  tdma->b[i] = -2;
		  tdma->c[i] =  1;
		}


	if(myrank==0)
	tdma->f[0] = tdma->f_left*tdma->h*tdma->h - tdma->u_left;

      	if(myrank==nproc-1)
	tdma->f[tdma->N-1] = tdma->f_right*tdma->h*tdma->h - tdma->u_right;

return 0;

}
