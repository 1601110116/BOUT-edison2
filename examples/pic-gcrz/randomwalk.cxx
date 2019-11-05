#include <cmath>
#include <istream>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

const double dl_tb = 0.1;

int main() {

	double gaussian();
	int npoints;
	//cout << "Enter the number of points." << endl;
	npoints = 10000;
	ofstream outfile("gauss.dat");
	for(int i=0;i<npoints;i++) {
		outfile << gaussian() << endl;
		//cout << gaussian() << endl;
	}
	outfile.close();

	return 0;
}

double gaussian() {
	int static iset=0;
	double static res2;
	double rsq,v1,v2,fac;

	if(iset ==0) {
		rsq=1.1;
		while (rsq >= 1.0 or rsq==0.0) {
			v1 = rand()/(double)(RAND_MAX);
			v1=2.0*v1-1.0;
			v2 = rand()/(double)(RAND_MAX);
			v2=2.0*v2-1.0;
			rsq=v1*v1+v2*v2;
		}
		fac = sqrt(-2.0*log(rsq)/rsq);
		iset = 1;
		res2 = v2*fac;
		return v1*fac;
	}
	else {
		iset=0;
		return res2;
	}
}


void randomwalk(double dr, double dz, double dr_tb, double dz_tb) {
	double  r,dl;

	r=sqrt(dr*dr+dz*dz);
	dl = gaussian()+dl_tb;

	dr_tb=dz/r*dl;
	dz_tb=-dr/r*dl;
}
