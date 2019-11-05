#include <mpi.h>
#include "tdma.hxx"
#include <fstream.h>
#include <iostream>
#include <sstream>
using namespace std;
int TDMA::dump_solution(TDMA *tdma)
{
ofstream fout;
int myrank;
int nproc;
stringstream out;
string procdir;
char *file = new char[procdir.size()+1];
int i;


myrank=MPI::COMM_WORLD.Get_rank();
nproc=MPI::COMM_WORLD.Get_size();

out << myrank;
procdir = "dump/solution."+out.str();

  std::copy(procdir.begin(),procdir.end(),file);
  file[procdir.size()]='\0';

  fout.open(file);
 
  for(i=0;i<=tdma->N-1;i++)
 	fout << u[i] <<endl;

  fout.close();
  delete[] file;

  return 0;

}
 
