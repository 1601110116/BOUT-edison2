#include "tdma.hxx"
#include "mpi_header.hxx"
#include <iostream>
using namespace std;

int TDMA::forward_sweep(TDMA *tdma)
{
  double *a,*b, *c, *f;
  double a_left, b_left, c_left, f_left;
  int i;
  double alpha;
  double alpha_left;
  int myrank;
  int nproc;
  double inittime;
  double endtime;
  double mpitime;
  

  MPI::Status status;
  MPI::Request send_req, recv_req;

  inittime=MPI::Wtime();  

  a=tdma->a;  
  b=tdma->b;
  c=tdma->c;
  f=tdma->f;


  myrank=MPI::COMM_WORLD.Get_rank();
  nproc=MPI::COMM_WORLD.Get_size(); 

  if(myrank!=0){
    MPI::COMM_WORLD.Recv(&a_left,1, MPI::DOUBLE,myrank-1,MPI::ANY_TAG,status);
    MPI::COMM_WORLD.Recv(&b_left,1, MPI::DOUBLE,myrank-1,MPI::ANY_TAG,status);
    MPI::COMM_WORLD.Recv(&c_left,1, MPI::DOUBLE,myrank-1,MPI::ANY_TAG,status);
    MPI::COMM_WORLD.Recv(&f_left,1, MPI::DOUBLE,myrank-1,MPI::ANY_TAG,status);
    std::cout <<"Receiving at \t"<<myrank<<endl;
   }

  if(myrank!=0)
   {
      alpha_left = a[0]/b_left;
      b[0]=b[0]-c_left*alpha_left;
      f[0]=f[0]-f_left*alpha_left;
   }

  for (i=1;i<=tdma->N-1;i++)
    {
      alpha=a[i]/b[i-1];
      b[i]=b[i] - c[i-1]*alpha;
      f[i]=f[i] -f[i-1]*alpha;
    }


if(myrank!=nproc-1)
    {
	send_req=MPI::COMM_WORLD.Isend(&a[tdma->N-1],1,MPI::DOUBLE,myrank+1,0);
        send_req.Wait(status);
	send_req=MPI::COMM_WORLD.Isend(&b[tdma->N-1],1,MPI::DOUBLE,myrank+1,0);
        send_req.Wait(status);
	send_req=MPI::COMM_WORLD.Isend(&c[tdma->N-1],1,MPI::DOUBLE,myrank+1,0);
        send_req.Wait(status);
	send_req=MPI::COMM_WORLD.Isend(&f[tdma->N-1],1,MPI::DOUBLE,myrank+1,0);
        send_req.Wait(status);

// might have to get rid of this for performance 
std::cout << "Sending from\t"<<myrank <<endl;
    }
  
  endtime=MPI::Wtime();

  mpitime=endtime-inittime;

    //cout << endl << "******" << endl;
//  std::cout << myrank <<"\t" << mpitime << std::endl;



  return 0;
  


}


	
