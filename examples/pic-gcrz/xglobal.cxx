/*******************************************************************************
*guiding-center module for divertor configuration(with X-points) in RZ coordinate
 *******************************************************************************/
#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <iostream>
#include <fstream>
#include <invert_parderiv.hxx>
#include <sourcex.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>
#include <bout/globalfield.hxx>
//#include "interp2d.cxx"
//#include <bout/globalfield.hxx>
#include "cylinder.cxx"
#include <bout/globalfield.hxx>

//Field3D fake;
//Field2D ndes;
//GlobalField2D *gndes;
//constant
const BoutReal Mp = 1.6726e-27; // proton mass
const BoutReal E_charge = 1.6022e-19; // elementary charge
const BoutReal PI=3.141592653589793238462643383279502884197;
const BoutReal TWOPI=2.0*PI;
const int neq=4;
const int kr=4,kz=4,kp=4,kpsi=4,ktheta=4;
//parameters
BoutReal Ti,v_th,bmag,rmag_rho,rmag,epsilon_l,v_para,v_perp,currenttime=0.,vmag,tmag,R_limiter, Z_limiter,t;
BoutReal yperiod;

//interpolation schemeGlobalField2D
int interpn;
//Electric field inclued or not
int flagEle,flagphi;
//particle's mass, charge, magnetic moment. unchanged
BoutReal AA,ZZ,mu,energy,y[4],Pzeta;

//particle's rho_para, position vector. Which will be evolved.
BoutReal rho_para,rr,rz,rp,orbitR,orbitZ,dt;
int nresult;

//field quantities
int nx,ny,ngy,ngz,zperiod;
int jyseps1_1,jyseps2_2,ixseps1,MYPE;
BoutReal zlength;
Field2D Rxy,Zxy;
BoutReal psia,psib;

Field3D fake;

	mesh1D meshpsi;
	scal1D fpol,pres,workk1,workk2,qpsi;
	meshRZ gmesh,coremesh;
	scalRZ psip,bphi,scal1,scal2,scal3,scal4,rcore,zcore;
	vecRZ vec_b,vec_gradb,vec_curlb,vec_bxgradb;
	int npsi_core,ntheta_core;
	BoutReal *gridpsi_core,*gridtheta_core,dtheta,dpsi;
	
const int npic=10,ntime=10,ndensity=100;
BoutReal nres[npic][5],density[ndensity][ntime];
//int locate[npic]; //recored particle's region. 0,core; 1,sol; 2, pf; 9,out of domain.
int indexrun=0;	
BoutReal ppsi[npic];
particle impurity[npic];
Field2D ndes;
GlobalField2D *gndes;
int nlcfs,nwall,ndomain;
BoutReal *r_lcfs,*z_lcfs,*r_wall,*z_wall,*r_domain, *z_domain;

BoutReal interp3d ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, BoutReal *zarray, int nz, Field3D &fxyz, BoutReal x, BoutReal y, BoutReal z);
void gcf(BoutReal t, BoutReal *y, int rkneq, BoutReal *dydt);
void rk4(BoutReal *y, int rkneq, BoutReal t, BoutReal dt, BoutReal *yout);
BoutReal gaussian();
int pnpoly(int nvert, BoutReal *vertx, BoutReal *verty, BoutReal testx, BoutReal testy);
void checkpsi(GlobalField2D &gpsixy);
void domain(GlobalField2D &gRxy, GlobalField2D &gZxy, int nv, BoutReal *vr, BoutReal *vz);

void region_y(BoutReal x, BoutReal &y, int indexregion);
void region_locate(BoutReal x, BoutReal y, int &indexregion);
void map_core(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek);

BoutReal dot (BoutReal v1r, BoutReal v1p, BoutReal v1z, BoutReal v2r,BoutReal v2p,BoutReal v2z) {
	return (v1r*v2r+v1p*v2p+v1z*v2z);
}

BoutReal dot (BoutReal v1r, BoutReal v1p, BoutReal v1z) {
	return (v1r*v1r+v1p*v1p+v1z*v1z);
}

cylinder cross (BoutReal v1r, BoutReal v1p, BoutReal v1z, BoutReal v2r,BoutReal v2p,BoutReal v2z) {
	cylinder res;
	res.r = v1p*v2z-v1z*v2p;
	res.p = v1z*v2r-v1r*v2z;
	res.z = v1r*v2p-v1p*v2r;
	return res;
}

BoutReal *gridr,*gridz,*gridpsi;
int nxefit,nyefit,maxstep;
GlobalField2D *gRxy, *gZxy, *gpsixy;
BoutReal simagx,sibdry;

int physics_init ( bool restarting ) {
	Field2D Bxy,Bpxy,Btxy,hthe,I,psixy;
	int xindex,yindex,MYPE;

	// Load metrics
	mesh->get(nx, "nx");
	mesh->get(ny, "ny");
	ngy = mesh->ngy;
	ngz = mesh->ngz;
	mesh->get(Rxy, "Rxy");
	mesh->get(Zxy, "Zxy");
	mesh->get(Bxy, "Bxy");
	mesh->get(Bpxy, "Bpxy");
	mesh->get(Btxy, "Btxy");
	mesh->get(hthe, "hthe");
	mesh->get(I,    "sinty");
	mesh->get(psixy, "psixy");
	mesh->get(jyseps1_1, "jyseps1_1");
	mesh->get(jyseps2_2, "jyseps2_2");
	mesh->get(ixseps1, "ixseps1");

	// Load normalisation values
	mesh->get(bmag, "bmag");
	mesh->get(rmag, "rmag");

	//
	Options *globalOptions = Options::getRoot();
	Options *options = globalOptions->getSection("gc");
	OPTION(options, AA, 1.0);
	OPTION(options, ZZ, 1.0);
	OPTION(options, Ti, 1.0e3);
	OPTION(options, v_para, 0.5);
	OPTION(options, rr, 0.8);
	OPTION(options, rz, 3.1415926);
	OPTION(options, rp, 0.0);
	OPTION(options, dt, 0.01);
	OPTION(options, maxstep, 1000);
	OPTION(options, interpn, 4);
	OPTION(options, flagEle, 0);
	OPTION(options, flagphi, 0);
	globalOptions->get("zperiod", zperiod, 1);



	//normalise
	Rxy = Rxy / rmag;
	Zxy = Zxy / rmag;
	Bxy = Bxy / bmag;
	mesh->Bxy = Bxy;
	Bpxy = Bpxy / bmag;
	Btxy = Btxy / bmag;
	hthe = hthe / rmag;
	I = I*bmag*rmag*rmag;
	psixy = psixy/bmag/rmag/rmag;
	mesh->dx = mesh->dx /(rmag*rmag*bmag);
	mesh->get(psia,"psi_axis");
	psia = psia /bmag/rmag/rmag;
	mesh->get(psib,"psi_bndry");
	psib = psib /bmag/rmag/rmag;
	//psixy,psia,psib is normalised to (0,1) futher.
	//psixy = (psixy-psia)/(psib-psia);
	
	//CALCULATE METRICS
	mesh->g11 = (Rxy*Bpxy)^2;
	mesh->g22 = 1.0 / (hthe^2);
	mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
	mesh->g12 = 0.0;
	mesh->g13 = -I*mesh->g11;
	mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
	mesh->J = hthe / Bpxy;
	mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
	mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
	mesh->g_33 = Rxy*Rxy;
	mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
	mesh->g_13 = I*Rxy*Rxy;
	mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
	mesh->geometry();


	ndes = 0.0;
	gndes = new GlobalField2D(mesh);
	gndes->gather(ndes);
	MPI_Comm_rank(BoutComm::get(), &MYPE);
	if (MYPE==0) {
		xindex=125;
		yindex=39;
		for(int i=0;i<10;i++) {
			(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			output.write("%5i%20.10e\n",i,(*gndes)(xindex,yindex));
		}
	}
	ndes = gndes->scatter();
	dump.add(ndes,"ndes",0);

	SOLVE_FOR(fake);
	return 0;

}

int physics_run(BoutReal t) {
	ddt(fake) = 0.;
	return 0;
}
