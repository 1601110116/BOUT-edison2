#include <iostream>
#include "tdma.hxx"
#include "fstream.h"
#include "mpi.h"
#include "string.h"
#include <sstream>

using namespace std;

int TDMA::print(TDMA *tdma)
{

  int i;
  ofstream fout;
  stringstream out;
  string procdir; 
  int myrank;
  int nproc;
  char *file = new char[procdir.size()+1];

  myrank=MPI::COMM_WORLD.Get_rank(); 
  nproc = MPI::COMM_WORLD.Get_size();
  
  out << myrank;
  procdir = "dump/"+out.str();	
  
  std::copy(procdir.begin(),procdir.end(),file);
  file[procdir.size()]='\0';

  
  fout.open(file);
  
  for(i=0;i<=tdma->N-1;i++)
    {
      fout << tdma->a[i]<< "\t" << tdma->b[i] <<"\t" << tdma->c[i] << "\t" 
	<< tdma->f[i] << endl;
    }

  fout.close();
  delete[] file;	


  return 0;
}
