#include <iostream>
#include "tdma.hxx"
#include "mpi_header.hxx"
#include "mpi.h"
#include <fstream.h>

using namespace std;
int main()
{
//Solves the tri-diagonal system Au=f
//where 
//A_i= (a_i b_i c_i)
//In expanded form, we attempt to solve
//a_i u_{i-1} +b_i u_i +c_i u_{i+1} = f_i

//get coefficients a_i, b_i, c_i 
//after elementary row transformation
//apply transformation to rows 2 .. N to convert to upper triangular system
//Do back substitution from rows N .. 1 to get solution.


double starttime;
double walltime;
int myrank;

ofstream outfile;

MPI_init();

starttime=MPI::Wtime();

TDMA tdma;
tdma.initialize(&tdma);
tdma.assemble_matrix(&tdma);
//prints general diagnostics
tdma.print(&tdma);
tdma.forward_sweep(&tdma);
tdma.back_sweep(&tdma);
//tdma.print(&tdma);
tdma.dump_solution(&tdma);
//tdma.destroy(&tdma);

walltime=MPI::Wtime()-starttime;

myrank = MPI::COMM_WORLD.Get_rank();
if(myrank==0)
{
outfile.open("dump/walltime.txt");
outfile <<"Walltime \t"<< walltime <<endl;
outfile.close();
}

MPI_finalize();

return 0;
}
