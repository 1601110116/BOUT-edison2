#include "tdma.hxx"

int TDMA::destroy(TDMA *tdma)
{
  delete [] tdma->a;
  delete [] tdma->b;
  delete [] tdma->c;
  delete [] tdma->f;
  delete [] tdma->u;
  
  return 0;
}
