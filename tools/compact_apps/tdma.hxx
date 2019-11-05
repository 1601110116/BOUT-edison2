#include "mpi.h"
class TDMA{

public:

  double *a;
  double *b;
  double *c;
  double *f;
  double *u; //solution vector

  int  N;

  double f_left;
  double f_right;
  double u_left ;
  double u_right ;
  double h;
  double  L;


  int initialize(TDMA *);
  int assemble_matrix(TDMA *);
  int print(TDMA *);
  int forward_sweep(TDMA *);
  int back_sweep(TDMA *);
  int dump_solution(TDMA *);
  int destroy(TDMA *);

};
