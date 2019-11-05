#include <iostream>
using namespace std;

int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy) {
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

int main() {
	int n=4;
	double *vx, *vy, x,y;

	vx = new double [n];
	vy = new double [n];
	vx[0]=0.0;
	vx[1]=1.0;
	vx[2]=1.0;
	vx[3]=0.0;
	vy[0]=0.0;
	vy[1]=0.0;
	vy[2]=1.0;
	vy[3]=1.0;

	cout << pnpoly(n,vx,vy,0.5,0.5) << endl;
	cout << pnpoly(n,vx,vy,0.5,1.5) << endl;
	return 0;
	
}
