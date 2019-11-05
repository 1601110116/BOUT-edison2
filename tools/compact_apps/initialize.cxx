#include <iostream>
#include "tdma.hxx"
#include "mpi.h"
#include "fstream.h"
using namespace std;
int TDMA::initialize(TDMA *tdma)
{
  int i;
  int myrank;
  int N;
  ifstream fin;
  ofstream out;
  int nproc;

myrank=MPI::COMM_WORLD.Get_rank();
nproc=MPI::COMM_WORLD.Get_size();


if(myrank==0)
{
//cout << "Enter the value of N"<< endl;
//cin >> tdma->N;
  fin.open("tdma.in");
  fin >> N;
  fin.close();

}


MPI::COMM_WORLD.Bcast(&N,1, MPI::INT,0);
tdma->N=N;

ofstream fout;

 //tdma->N=N;
 tdma->a=new double[tdma->N];
 tdma->b=new double[tdma->N];
 tdma->c=new double[tdma->N];
 tdma->f=new double[tdma->N];
 tdma->u=new double[tdma->N]; 

 tdma->f_left = -1;
 tdma->f_right = -1;
 tdma->u_left = 0;
 tdma->u_right = 0;

 tdma->L = 1;

//weak scaling
 tdma->h = tdma->L/(double(tdma->N)-1);

 for (i=0;i<=tdma->N-1;i++)
   {
	  tdma->a[i]=0;
	  tdma->b[i]=0;
	  tdma->c[i]=0;
	  tdma->f[i]=0;
          tdma->u[i] = 0;
   }




 return 0;
}
