#include "tdma.hxx"
#include <iostream>
#include <mpi.h>
using namespace std;
int TDMA::back_sweep(TDMA *tdma)
{
  double *a, *b, *c, *f, *u;
  int N;
  int i;
  MPI::Status status;
  MPI::Request send_req, recv_req;
  int myrank;
  int nproc;
  double val_right;


  a=tdma->a;
  b=tdma->b;
  c=tdma->c;
  f=tdma->f;
  u=tdma->u;
  N=tdma->N;

  
  myrank=MPI::COMM_WORLD.Get_rank();
  nproc=MPI::COMM_WORLD.Get_size();

  if(myrank==nproc-1)
  u[N-1] = f[N-1]/b[N-1];

  
   if(myrank!=nproc-1)
      {
	MPI::COMM_WORLD.Recv(&val_right,1,MPI::DOUBLE,myrank+1,MPI::ANY_TAG,status);
 		u[N-1]=(f[N-1]-c[N-1]*val_right)/b[N-1];
      }
 	

	for (i=N-2;i>=0;i--)
	{
		u[i]=(f[i]-c[i]*u[i+1])/b[i];
	}

  if(myrank!=0)
  {
  send_req=MPI::COMM_WORLD.Isend(&u[0],1,MPI::DOUBLE,myrank-1,0);
  send_req.Wait(status);
  }
  

return 0;
}

