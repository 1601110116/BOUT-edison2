/*******************************************************************************
*guiding-center module for ciucular geometry. run one sigle processor now.
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

#include <iostream>
#include <fstream>

#include <invert_parderiv.hxx>
#include <sourcex.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>
#include "interp2d.cxx"
#include "Csplines.cxx"

//constant
const BoutReal Mp = 1.6726e-27; // proton mass
const BoutReal E_charge = 1.6022e-19; // elementary charge
const BoutReal PI=3.141592653589793238462643383279502884197;
const BoutReal TWOPI=2.0*PI;
const int neq=4,kx=4,ky=4;
//parameters
BoutReal Ti,v_th,bmag,rmag_rho,rmag,epsilon_l,v_para,v_perp,currenttime=0.;

//interpolation scheme
int interpn;
//Electric field inclued or not
int flagEle,flagphi;
//particle's mass, charge, magnetic moment. unchanged
BoutReal AA,ZZ,mu;

//particle's rho_para, position vector. Which will be evolved.
BoutReal rho_para,rvecx,rvecy,rvecz,orbitR,orbitZ,dt;
int nresult;
BoutReal tmpb,energy,Pzeta,tmpbt,tmpphi,tmpphi0,tmpphi1;
BoutReal result[4],result2[4];

//field quantities
int nx,ny,ngy,ngz,zperiod;
BoutReal zlength;
Field2D scaf1,scaf2,scaf3,scaf4,Rxy,Zxy;
Field2D phi0scaf1,phi0scaf2;
Field3D phi1,phi1scaf1,phi1scaf2;
Vector2D vecf1,vecf2,vecf3;
Vector2D tmpvecf2d,phi0vecf1;
Vector3D tmpvecf3d,phi1vecf1;
BoutReal *xgrid, *knotx, *ygrid, *knoty, *zgrid;
BoutReal *shiftangle, *bcoefshiftangle;
BoutReal psia,psib;
BoutReal ***phi1t;
BoutReal ****phi1t3d;
BoutReal **bcoef2ds1, **bcoef2ds2, **bcoef2ds3, **bcoef2ds4, **bcoef2dR, **bcoef2dZ, **bcoef2dBt;
BoutReal **bcoef2dp0s1, **bcoef2dp0s2,  **bcoef2dphi0;
BoutReal **bcoef2dv1x, **bcoef2dv1y, **bcoef2dv1z;
BoutReal **bcoef2dv2x, **bcoef2dv2y, **bcoef2dv2z;
BoutReal **bcoef2dv3x, **bcoef2dv3y, **bcoef2dv3z;
BoutReal **bcoef2dp0v1x, **bcoef2dp0v1y, **bcoef2dp0v1z;

Field3D fake;

int locate (BoutReal *xvec, int n, BoutReal x);
BoutReal interp1dlin ( BoutReal *xarray, int nx, BoutReal *fx, BoutReal x);
BoutReal interp2dlin ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y);
BoutReal interpqua ( BoutReal *xarray, int nx, BoutReal *yarray, BoutReal x);
BoutReal interp2dqua ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y);
BoutReal interpcub ( BoutReal *xarray, int nx, BoutReal *yarray, BoutReal x);
BoutReal interp2dcub ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y);
void cdbbsnak(BoutReal *x, int n, BoutReal *knot, int nk);
void cdbbscoef2d(BoutReal *x, BoutReal *y, Field2D &fxy, int nx, int ny, BoutReal *knotx, BoutReal *knoty, int kx, int ky, BoutReal **coef2d);
void cdbbscoef2d(BoutReal *x, BoutReal *y, BoutReal **fxy, int nx, int ny, BoutReal *knotx, BoutReal *knoty, int kx, int ky, BoutReal **coef2d);
BoutReal cdbbsval2d(BoutReal *knotx, int nkx, BoutReal *knoty, int nky, BoutReal **coef2d, int nx, int ny, BoutReal x, BoutReal y, int derivx, int derivy);
void gcf(BoutReal t, BoutReal *y, int rkneq, BoutReal *dydt);
void rk4(BoutReal *y, int rkneq, BoutReal t, BoutReal dt, BoutReal *yout);
void periodicboundary(Field2D &fxy);
void periodicboundary(Vector2D &fxy);
void periodicboundary(Field3D &fxyz);
void periodicboundary(Vector3D &fxyz);
void smoothboundary(Field2D &fxy);
void smoothboundary(Vector2D &fxy);
void smoothboundary(Field3D &fxyz);
void smoothboundary(Vector3D &fxyz);
BoutReal interp2d ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y) {
  if (interpn == 2)
      return interp2dlin ( xarray, nx, yarray, ny, fxy, x, y);
   else if (interpn == 3)
      return interp2dqua ( xarray, nx, yarray, ny, fxy, x, y);
   else if (interpn == 4)
       return interp2dcub ( xarray, nx, yarray, ny, fxy, x, y);
}
BoutReal interp3d ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, BoutReal *zarray, int nz, Field3D &fxyz, BoutReal x, BoutReal y, BoutReal z);
Field2D phi0;

int physics_init ( bool restarting )
{
	Field2D Bxy,Bpxy,Btxy,hthe,I,psixy,shiftangle2d;
	Vector2D B0vec,unit_B0vec;

	ifstream inFile;

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

	shiftangle = new BoutReal [nx];
	mesh->get(shiftangle2d, "shiftangle2d");
	for(int i=0;i<nx;i++)
		shiftangle[i]=shiftangle2d[i][0];

	phi0 = 0.0;
/*	inFile.open ("data/sheathphi0.dat");
	if(inFile.is_open()) {
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
				{inFile >> phi0[i][2+j];}
	}
	else output.write("can't open sheathphi0.dat\n");
	inFile.close();
	for(int i=0;i<nx;i++)
		for(int j=0;j<2;j++)
			{phi0[i][ny+2+j]=phi0[i][2+j];phi0[i][j]=phi0[i][ny+j];}
*/
	/*int nphi1t = 72;
	phi1t = new BoutReal **[nx];
	for (int i=0;i<nx;i++) {
		phi1t[i] = new BoutReal *[ngy];
		for (int j=0;j<ngy;j++)
			{phi1t[i][j] = new BoutReal [nphi1t];}
	}
	inFile.open ("data/sheathphi1.dat");
	if(inFile.is_open()) {
		for(int k=0;k<nphi1t;k++)
			for(int j=0;j<ny;j++)
				for(int i=0;i<nx;i++)
					{inFile >> phi1t[i][2+j][k];}
	}
	else output.write("can't open sheathphi1.dat\n");
	inFile.close(); 
	for(int k=0;k<nphi1t;k++)
		for(int i=0;i<nx;i++)
			for(int j=0;j<2;j++)
				{phi1t[i][ny+2+j][k]=phi1t[i][2+j][k];phi1t[i][j][k]=phi1t[i][ny+j][k];} */

/*	int nphi1t3d = 172;
	phi1t3d = new BoutReal ***[nx];
	for (int i=0;i<nx;i++) {
			phi1t3d[i] = new BoutReal **[ngy];
			for (int j=0;j<ngy;j++) {
				phi1t3d[i][j] = new BoutReal *[ngz];
				for (int k=0;k<ngz;k++)
					phi1t3d[i][j][k] = new BoutReal [nphi1t3d];
			}
	}
	inFile.open ("data/sheathphi1-3D.dat");
	if(inFile.is_open()) {
		for(int l=0;l<nphi1t3d;l++)
			for(int k=0;k<ngz-1;k++)
				for(int j=0;j<ny;j++)
					for(int i=0;i<nx;i++)
						{inFile >> phi1t3d[i][2+j][k][l];}
	}
	else output.write("can't open sheathphi1-3D.dat\n");
	inFile.close();
*/
	phi1 =0.0;
/*	for(int k=0;k<ngz-1;k++)
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
				{phi1[i][2+j][k]=phi1t3d[i][2+j][k][0];}
*/
	// Load normalisation values
	mesh->get(bmag, "bmag");
	mesh->get(rmag, "rmag");

	/*************** READ OPTIONS *************************/
	Options *globalOptions = Options::getRoot();
	Options *options = globalOptions->getSection("gc");
	OPTION(options, AA, 1.0);
	OPTION(options, ZZ, 1.0);
	OPTION(options, Ti, 1.0e3);
	OPTION(options, v_para, 0.5);
	OPTION(options, rvecx, 0.8);
	OPTION(options, rvecy, 3.1415926);
	OPTION(options, rvecz, 0.0);
	OPTION(options, dt, 0.01);
	OPTION(options, nresult, 1000);
	OPTION(options, interpn, 4);
	OPTION(options, flagEle, 0);
	OPTION(options, flagphi, 0);
	globalOptions->get("zperiod", zperiod, 1);

	/*for (int i=0;i<nx;i++) {
		for (int j=0;j<ngy;j++)
			{delete phi1t[i][j];}
		delete phi1t[i];
	}
	delete phi1t;*/

	v_th = sqrt(2.0*E_charge*Ti/Mp);
	rmag_rho = Mp*v_th/E_charge/bmag;
	epsilon_l=rmag_rho/rmag;
	v_perp = sqrt(1.0-v_para*v_para);

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
	//psixy is normalised to (0,1) futher.
	psixy = (psixy-psia)/(psib-psia);
	//denormalisation and re-normalisation of phi
	BoutReal va=9.468397e06;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ngy;j++)
			{phi0[i][j] = va*rmag*bmag/Ti*(phi0[i][j]*Bxy[i][j]);}
	for(int i=0;i<nx;i++)
		for(int j=0;j<ngy;j++)
			for(int k=0;k<ngz;k++)
				{phi1[i][j][k] = va*rmag*bmag/Ti*(phi1[i][j][k]*Bxy[i][j]);}
	//phi0 is generated analytically
	BoutReal *vx,vxmin,vxmax;
	vx = new BoutReal [nx];
	for (int i=0;i<nx;i++)
		{vx[i]=psixy(i,0);}
	vxmin=vx[0];
	vxmax=vx[nx-1];
	for (int i=0;i<nx;i++)
		{vx[i]=(vx[i]-vxmin)/(vxmax-vxmin)*PI;}
	for(int i=0;i<nx;i++)
		for(int j=0;j<ngy;j++)
			{phi0[i][j]=cos(vx[i]);}
	for(int i=0;i<nx;i++)
		for(int j=0;j<2;j++)
			{phi0[i][ny+2+j]=phi0[i][2+j];phi0[i][j]=phi0[i][ny+j];}
	//

	// Set B field
	B0vec.covariant = false;
	B0vec.x = 0.;
	B0vec.y = Bpxy / hthe;
	B0vec.z = 0.;
	unit_B0vec.covariant = false;
	unit_B0vec.x = 0.;
	unit_B0vec.y = B0vec.y / Bxy;
	unit_B0vec.z = 0.;
	

	/**************** CALCULATE METRICS ******************/
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

	//set arrays of grids in x,y,z dimensions
	xgrid = new BoutReal [nx];
	ygrid = new BoutReal [ngy];
	zgrid = new BoutReal [ngz];
	for(int i=0;i<nx;i++)
		{xgrid[i]=psixy[i][3];}
	BoutReal dy = 2.0*PI/ny;
	ygrid[0] = -2.0*dy;
	for(int i=1;i<ngy;i++)
		{ygrid[i] = ygrid[i-1]+dy;}
	zlength = TWOPI/(BoutReal(zperiod));
	BoutReal dz = zlength/(ngz-1);
	zgrid[0]=0.0;
	for(int i=1;i<ngz;i++)
		{zgrid[i]=zgrid[i-1]+dz;}

	//fill boundary guard cell of phi1 with values using periodic condition. notice the zshift.
	 periodicboundary(phi1); // periodicboundary(Field3D &fxyz), this function fill guard cell.

	//scale "Bxy, unit_B0vec*Curl(unit_B0vec)" and vector "Grad(Bxy), B0vec, Curl(B0vec), Grad(phi)" 
	// are needed to be interplated when advancing particle.

	vecf1.covariant = false;
	vecf2.covariant = false;
	vecf3.covariant = false;
	phi0vecf1.covariant = false;
	tmpvecf2d.covariant = false;
	phi1vecf1.covariant = false;
	tmpvecf3d.covariant = false;
/*	vecf1.setBoundary("periodic");
	vecf2.setBoundary("periodic");
	vecf3.setBoundary("periodic");
	tmpvecf2d.setBoundary("periodic");
	tmpvecf3d.setBoundary("periodic");
	scaf2.setBoundary("periodic");
	scaf3.setBoundary("periodic");
	scaf4.setBoundary("periodic");
	phi0scaf1.setBoundary("periodic");
	phi0scaf2.setBoundary("periodic");
	phi0vecf1.setBoundary("periodic");
	phi1scaf1.setBoundary("periodic");
	phi1scaf2.setBoundary("periodic");
	phi1vecf1.setBoundary("periodic");*/

	scaf1 = Bxy;

	tmpvecf2d = Curl(unit_B0vec);
//	tmpvecf2d.applyBoundary();
	periodicboundary( tmpvecf2d );
//	smoothboundary( tmpvecf2d );

	scaf2 = unit_B0vec*tmpvecf2d;
//	scaf2.applyBoundary();
	periodicboundary( scaf2 );
//	smoothboundary( scaf2 );

	vecf1 = B0vec;

	vecf2 = Curl(B0vec);
//	vecf2.applyBoundary();
	periodicboundary( vecf2 );
//	smoothboundary( vecf2 );

	vecf3 = Grad(Bxy);
//	vecf3.applyBoundary();
	periodicboundary( vecf3 );
//	smoothboundary( vecf3 );

	scaf3 = vecf3*vecf1;
//	scaf3.applyBoundary();
	periodicboundary( scaf3 );
//	smoothboundary( scaf3 );

	scaf4 = vecf3*vecf2;
//	scaf4.applyBoundary();
	periodicboundary( scaf4 );
//	smoothboundary( scaf4 );

	vecf3 = vecf1^vecf3;
//	vecf3.applyBoundary();
	periodicboundary( vecf3 );
//	smoothboundary( vecf3 );
	
	tmpvecf2d = Grad(phi0);
//	tmpvecf2d.applyBoundary();
	periodicboundary( tmpvecf2d );
//	smoothboundary( tmpvecf2d );

	phi0scaf1 = tmpvecf2d*vecf1;
//	phi0scaf1.applyBoundary();
	periodicboundary( phi0scaf1 );
//	smoothboundary( phi0scaf1 );

	phi0scaf2 = tmpvecf2d*vecf2;
//	phi0scaf2.applyBoundary();
	periodicboundary( phi0scaf2 );
//	smoothboundary( phi0scaf2 );
	//phi0scaf2 = 0.0;

	phi0vecf1 = vecf1^tmpvecf2d;
//	phi0vecf1.applyBoundary();
	periodicboundary( phi0vecf1 );
//	smoothboundary( phi0vecf1 );

	tmpvecf3d = Grad(phi1);
//	tmpvecf3d.applyBoundary();
	periodicboundary( tmpvecf3d );
//	smoothboundary( tmpvecf3d );

	phi1scaf1 = tmpvecf3d*vecf1;
//	phi1scaf1.applyBoundary();
	periodicboundary( phi1scaf1 );
//	smoothboundary( phi1scaf1 );

	phi1scaf2 = tmpvecf3d*vecf2;
//	phi1scaf2.applyBoundary();
	periodicboundary( phi1scaf2 );
//	smoothboundary( phi1scaf2 );

	phi1vecf1 = vecf1^tmpvecf3d;
//	phi1vecf1.applyBoundary();
	periodicboundary( phi1vecf1 );
//	smoothboundary( phi1vecf1 );

	dump.add(phi0scaf1,"phi0scaf1",0);
	dump.add(phi0scaf2,"phi0scaf2",0);
	dump.add(phi0vecf1.x,"phi0vecf1x",0);
	dump.add(phi0vecf1.y,"phi0vecf1y",0);
	dump.add(phi0vecf1.z,"phi0vecf1z",0);
	dump.add(phi1scaf1,"phi1scaf1",0);
	dump.add(phi1scaf2,"phi1scaf2",0);
	dump.add(phi1vecf1.x,"phi1vecf1x",0);
	dump.add(phi1vecf1.y,"phi1vecf1y",0);
	dump.add(phi1vecf1.z,"phi1vecf1z",0);

//Bspline
	knotx = new BoutReal [nx+kx];
	knoty = new BoutReal [ngy+ky];

	bcoefshiftangle = new BoutReal [nx];

	bcoef2ds1 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2ds1[i] = new BoutReal [ngy];

	bcoef2ds2 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2ds2[i] = new BoutReal [ngy];

	bcoef2ds3 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2ds3[i] = new BoutReal [ngy];

	bcoef2ds4 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2ds4[i] = new BoutReal [ngy];
	
   bcoef2dphi0 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dphi0[i] = new BoutReal [ngy];

	bcoef2dp0s1 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dp0s1[i] = new BoutReal [ngy];
   
	bcoef2dp0s2 = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dp0s2[i] = new BoutReal [ngy];

	bcoef2dR = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dR[i] = new BoutReal [ngy];

	bcoef2dZ = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dZ[i] = new BoutReal [ngy];

	bcoef2dBt = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dBt[i] = new BoutReal [ngy];

	bcoef2dv1x = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv1x[i] = new BoutReal [ngy];

	bcoef2dv1y = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv1y[i] = new BoutReal [ngy];

	bcoef2dv1z = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv1z[i] = new BoutReal [ngy];

	bcoef2dv2x = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv2x[i] = new BoutReal [ngy];

	bcoef2dv2y = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv2y[i] = new BoutReal [ngy];

	bcoef2dv2z = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv2z[i] = new BoutReal [ngy];

	bcoef2dv3x = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv3x[i] = new BoutReal [ngy];

	bcoef2dv3y = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv3y[i] = new BoutReal [ngy];

	bcoef2dv3z = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dv3z[i] = new BoutReal [ngy];

	bcoef2dp0v1x = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dp0v1x[i] = new BoutReal [ngy];

	bcoef2dp0v1y = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dp0v1y[i] = new BoutReal [ngy];

	bcoef2dp0v1z = new BoutReal *[nx];
	for (int i=0;i<nx;i++)
		bcoef2dp0v1z[i] = new BoutReal [ngy];

	cdbbsnak(xgrid,nx,knotx,nx+kx);
	cdbbsnak(ygrid,ngy,knoty,ngy+ky);

	cdbbscoef(xgrid, shiftangle, nx, knotx, kx, bcoefshiftangle);
	cdbbscoef2d(xgrid, ygrid, scaf1, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2ds1);
	cdbbscoef2d(xgrid, ygrid, scaf2, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2ds2);
	cdbbscoef2d(xgrid, ygrid, scaf3, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2ds3);
	cdbbscoef2d(xgrid, ygrid, scaf4, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2ds4);
	cdbbscoef2d(xgrid, ygrid, phi0, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dphi0);
	cdbbscoef2d(xgrid, ygrid, phi0scaf1, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dp0s1);
	cdbbscoef2d(xgrid, ygrid, phi0scaf2, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dp0s2);
	cdbbscoef2d(xgrid, ygrid, Rxy, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dR);
	cdbbscoef2d(xgrid, ygrid, Zxy, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dZ);
	cdbbscoef2d(xgrid, ygrid, Btxy, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dBt);
	cdbbscoef2d(xgrid, ygrid, vecf1.x, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv1x);
	cdbbscoef2d(xgrid, ygrid, vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv1y);
	cdbbscoef2d(xgrid, ygrid, vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv1z);
	cdbbscoef2d(xgrid, ygrid, vecf2.x, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv2x);
	cdbbscoef2d(xgrid, ygrid, vecf2.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv2y);
	cdbbscoef2d(xgrid, ygrid, vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv2z);
	cdbbscoef2d(xgrid, ygrid, vecf3.x, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv3x);
	cdbbscoef2d(xgrid, ygrid, vecf3.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv3y);
	cdbbscoef2d(xgrid, ygrid, vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dv3z);
	cdbbscoef2d(xgrid, ygrid, phi0vecf1.x, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dp0v1x);
	cdbbscoef2d(xgrid, ygrid, phi0vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dp0v1y);
	cdbbscoef2d(xgrid, ygrid, phi0vecf1.y, nx, ngy,knotx, knoty,kx, ky, (BoutReal**)bcoef2dp0v1z);

	output.write("normalisation: \n");
	output.write("bmag: \n");
	output.write("%lf\n",bmag);
	output.write("rmag: \n");
	output.write("%lf\n",rmag);
	output.write("rmag_rho: \n");
	output.write("%lf\n",rmag_rho);
	output.write("epsilon_l: \n");
	output.write("%lf\n",epsilon_l);
	output.write("v_th: \n");
	output.write("%lf\n",v_th);
	output.write("T_gyro: \n");
	output.write("%15.5e\n",Mp/E_charge/bmag);
	output.write("T_ploidal: \n");
	output.write("%15.5e\n",rmag/v_th);
    
	SOLVE_FOR(fake);

	dump.add(phi0,"phi0",0);

	int tmpn;
	BoutReal tmpy,tmpshiftangle,tmpz;
	//BoutReal t,result[4],result2[4];
	v_perp = sqrt(1.0-v_para*v_para);
	tmpn=floor(rvecy / TWOPI);
	tmpy = rvecy -  tmpn * TWOPI;
	tmpshiftangle = cdbbsval(knotx, nx+kx, bcoefshiftangle, nx, rvecx, 0);
	tmpz = rvecz - tmpn * tmpshiftangle;
	tmpz = tmpz - floor(tmpz/zlength)*zlength;
	BoutReal B_init;
	//B_init = interp2d(xgrid,nx,ygrid,ngy,scaf1,rvecx,tmpy);
	B_init = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds1, nx, ngy,rvecx,tmpy,0,0);
	mu=AA*v_perp*v_perp/B_init;  
	rho_para=AA*v_para/ZZ/B_init;
	//t=0.0;
	result[0]=rvecx;
	result[1]=rvecy;
	result[2]=rvecz;
	result[3]=rho_para;
	output.write("result: \n");
	//orbitR = interp2d(xgrid,nx,ygrid,ngy,Rxy,result[0],tmpy);
	orbitR = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dR, nx, ngy,result[0],tmpy,0,0);
	//orbitZ = interp2d(xgrid,nx,ygrid,ngy,Zxy,result[0],tmpy);
	orbitZ = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dZ, nx, ngy,result[0],tmpy,0,0);
	tmpb = B_init;
	tmpphi0 =  cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dphi0, nx, ngy,result[0],tmpy,0,0);
	tmpphi1 = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1,result[0],tmpy,tmpz);
	
	tmpphi = tmpphi0 + tmpphi1;
	energy = mu*tmpb + (ZZ*result[3]*tmpb)*(ZZ*result[3]*tmpb)/AA + tmpphi*ZZ;
	//tmpbt = interp2d(xgrid,nx,ygrid,ngy,Btxy,result[0],tmpy);
	tmpbt = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dBt, nx, ngy,result[0],tmpy,0,0);
	//Pzeta = result[0]-epsilon_l*orbitR*tmpbt*result[3];
	Pzeta = result[0]*(psib-psia)+psia-epsilon_l*orbitR*tmpbt*result[3];

	dump.add(orbitR,"R",1);
	dump.add(orbitZ,"Z",1);
	dump.add(mu,"mu",1);
	dump.add(energy,"energy",1);
	dump.add(Pzeta,"pzeta",1);
	return 0;
}


int physics_run(BoutReal t) {
	int tmpn;
	BoutReal tmpy,tmpshiftangle,tmpz;


	if( t>0.0 and fmod(t,1.0) ==0.0) {

		rk4(result, neq, t, dt, result2);
		for(int j=0;j<4;j++)
			{result[j]=result2[j];}
		tmpn=floor(result[1] / TWOPI);
		tmpy = result[1] -  tmpn * TWOPI;
		tmpshiftangle = cdbbsval(knotx, nx+kx, bcoefshiftangle, nx, result[0], 0);
		tmpz = result[2] - tmpn * tmpshiftangle;
		tmpz = tmpz - floor(tmpz/zlength)*zlength;
		// orbitR = interp2d(xgrid,nx,ygrid,ngy,Rxy,result[0],tmpy);
		orbitR = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dR, nx, ngy,result[0],tmpy,0,0);
		//orbitZ = interp2d(xgrid,nx,ygrid,ngy,Zxy,result[0],tmpy);
		orbitZ = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dZ, nx, ngy,result[0],tmpy,0,0);
		// tmpb = interp2d(xgrid,nx,ygrid,ngy,scaf1,result[0],tmpy);
		tmpb = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds1, nx, ngy,result[0],tmpy,0,0);
		tmpphi0 =  cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dphi0, nx, ngy,result[0],tmpy,0,0);
		tmpphi1 = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1,result[0],tmpy,tmpz);
		tmpphi = tmpphi0 + tmpphi1;
		energy = mu*tmpb + (ZZ*result[3]*tmpb)*(ZZ*result[3]*tmpb)/AA + tmpphi*ZZ;
		//tmpbt = interp2d(xgrid,nx,ygrid,ngy,Btxy,result[0],tmpy);
		tmpbt = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dBt, nx, ngy,result[0],tmpy,0,0);
		//Pzeta = result[0] - epsilon_l*orbitR*tmpbt*result[3];
		Pzeta = result[0]*(psib-psia)+psia-epsilon_l*orbitR*tmpbt*result[3];
		
		
	}

	ddt(fake) = 0.;
	return 0;
}

void gcf(BoutReal t, BoutReal *y, int rkneq, BoutReal *dydt) {
	BoutReal tmpsca1,tmpsca2,tmpsca3,tmpsca4,tmpphisca1,tmpphisca2,D,tmp1,tmpx,tmpy,tmpz,tmprho_para;
	BoutReal tmpvec1[3],tmpvec2[3],tmpvec3[3],tmpphivec1[3],tmpshiftangle;
	BoutReal tmpphi0sca1,tmpphi1sca1,tmpphi0sca2,tmpphi1sca2,tmpphi0vec1[3],tmpphi1vec1[3];
	int tmpn;
	tmpx = y[0];
	tmpn=floor(y[1] / TWOPI);
	tmpy = y[1] -  tmpn * TWOPI;
	tmpshiftangle = cdbbsval(knotx, nx+kx, bcoefshiftangle, nx, y[0], 0);
	tmpz = y[2] - tmpn * tmpshiftangle;
	tmpz = tmpz - floor(tmpz/zlength)*zlength;
	tmprho_para = y[3];

	//tmpsca1 = interp2d(xgrid,nx,ygrid,ngy,scaf1,tmpx,tmpy);
	tmpsca1 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds1, nx, ngy,tmpx,tmpy,0,0);
	//tmpsca2 = interp2d(xgrid,nx,ygrid,ngy,scaf2,tmpx,tmpy);
	tmpsca2 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds2, nx, ngy,tmpx,tmpy,0,0);
	//tmpsca3 = interp2d(xgrid,nx,ygrid,ngy,scaf3,tmpx,tmpy);
	tmpsca3 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds3, nx, ngy,tmpx,tmpy,0,0);
	//tmpsca4 = interp2d(xgrid,nx,ygrid,ngy,scaf4,tmpx,tmpy);
	tmpsca4 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2ds4, nx, ngy,tmpx,tmpy,0,0);
	//tmpphi0sca1 =  interp2d(xgrid,nx,ygrid,ngy,phi0scaf1,tmpx,tmpy);
	tmpphi0sca1 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dp0s1, nx, ngy,tmpx,tmpy,0,0);
	tmpphi1sca1 = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1scaf1,tmpx,tmpy,tmpz);
	tmpphisca1 = tmpphi0sca1 + tmpphi1sca1;
	//tmpphi0sca2 =  interp2d(xgrid,nx,ygrid,ngy,phiscaf2,tmpx,tmpy);
	tmpphi0sca2 = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dp0s2, nx, ngy,tmpx,tmpy,0,0);
	tmpphi1sca2 = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1scaf2,tmpx,tmpy,tmpz);
	tmpphisca2 = tmpphi0sca2 + tmpphi1sca2;
	//tmpvec1[0] = interp2d(xgrid,nx,ygrid,ngy,(vecf1.x),tmpx,tmpy);
	tmpvec1[0] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv1x, nx, ngy,tmpx,tmpy,0,0);
	//tmpvec1[1] = interp2d(xgrid,nx,ygrid,ngy,(vecf1.y),tmpx,tmpy);
	tmpvec1[1] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv1y, nx, ngy,tmpx,tmpy,0,0);
	tmpvec1[2] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv1z, nx, ngy,tmpx,tmpy,0,0);
	//tmpvec2[0] = interp2d(xgrid,nx,ygrid,ngy,(vecf2.x),tmpx,tmpy);
	tmpvec2[0] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv2x, nx, ngy,tmpx,tmpy,0,0);
	//tmpvec2[1] = interp2d(xgrid,nx,ygrid,ngy,(vecf2.y),tmpx,tmpy);
	tmpvec2[1] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv2y, nx, ngy,tmpx,tmpy,0,0);
	tmpvec2[2] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv2z, nx, ngy,tmpx,tmpy,0,0);
	//tmpvec3[0] = interp2d(xgrid,nx,ygrid,ngy,(vecf3.x),tmpx,tmpy);
	tmpvec3[0] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv3x, nx, ngy,tmpx,tmpy,0,0);
	//tmpvec3[1] = interp2d(xgrid,nx,ygrid,ngy,(vecf3.y),tmpx,tmpy);
	tmpvec3[1] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv3y, nx, ngy,tmpx,tmpy,0,0);
	tmpvec3[2] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dv3z, nx, ngy,tmpx,tmpy,0,0);
	//tmpphi0vec1[0] = interp2d(xgrid,nx,ygrid,ngy,(phi0vecf1.x),tmpx,tmpy);
	tmpphi0vec1[0] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dp0v1x, nx, ngy,tmpx,tmpy,0,0);
	tmpphi1vec1[0] = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1vecf1.x,tmpx,tmpy,tmpz);
	tmpphivec1[0] = tmpphi0vec1[0] + tmpphi1vec1[0];
	//tmpphivec1[1] = interp2d(xgrid,nx,ygrid,ngy,(phi0vecf1.y),tmpx,tmpy);
	tmpphi0vec1[1] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dp0v1y, nx, ngy,tmpx,tmpy,0,0);
	tmpphi1vec1[1] = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1vecf1.y,tmpx,tmpy,tmpz);
	tmpphivec1[1] = tmpphi0vec1[1] + tmpphi1vec1[1];
	tmpphi0vec1[2] = cdbbsval2d(knotx,nx+kx,knoty,ngy+ky,(BoutReal**)bcoef2dp0v1z, nx, ngy,tmpx,tmpy,0,0);
	tmpphi1vec1[2] = interp3d(xgrid,nx,ygrid,ngy,zgrid,ngz,phi1vecf1.z,tmpx,tmpy,tmpz);
	tmpphivec1[2] = tmpphi0vec1[2] + tmpphi1vec1[2];

	D = 1.0+epsilon_l*tmprho_para*(tmpsca2);
	tmp1 = 2.0*ZZ*tmprho_para*tmprho_para/AA*tmpsca1+mu/ZZ;
	dydt[0] = 1.0/D*(ZZ*tmprho_para/AA*tmpvec1[0]+epsilon_l*tmpvec2[0]+0.5*epsilon_l/tmpsca1/tmpsca1*(tmp1*tmpvec3[0]+tmpphivec1[0]));
	dydt[0] = dydt[0]/(psib-psia);
	dydt[1] = 1.0/D*(ZZ*tmprho_para/AA*tmpvec1[1]+epsilon_l*tmpvec2[1]+0.5*epsilon_l/tmpsca1/tmpsca1*(tmp1*tmpvec3[1]+tmpphivec1[1]));
	dydt[2] = 1.0/D*(ZZ*tmprho_para/AA*tmpvec1[2]+epsilon_l*tmpvec2[2]+0.5*epsilon_l/tmpsca1/tmpsca1*(tmp1*tmpvec3[2]+tmpphivec1[2]));
	dydt[3] = -0.5/(tmpsca1*tmpsca1*D)*(tmp1*(tmpsca3+epsilon_l*tmprho_para*tmpsca4)+tmpphisca1 + epsilon_l*tmprho_para*tmpphisca2);
}

void rk4(BoutReal *y, int rkneq, BoutReal t, BoutReal dt, BoutReal *yout) {
	BoutReal dydx[rkneq],yt[rkneq],dyt[rkneq],dym[rkneq],dth,dt6,th;

	dth=dt*0.5;
	dt6=dt/6.0;
	th=t+dth;

	gcf(t,y,rkneq,dydx);
	for(int i=0;i < rkneq;i++)
	{yt[i] = y[i]+dth*dydx[i];}


	gcf(th,yt,rkneq,dyt);
	for(int i=0;i < rkneq;i++)
	{yt[i] = y[i]+dth*dyt[i];}

	gcf(th,yt,rkneq,dym);
	for(int i=0;i < rkneq;i++)
	{yt[i]=y[i]+dt*dym[i];dym[i]=dyt[i]+dym[i];}

	gcf(t+dt,yt,rkneq,dyt);
	for(int i=0;i < rkneq;i++)
	{yout[i]=y[i]+dt6*(dydx[i]+dyt[i]+2.0*dym[i]);}
}

void periodicboundary(Field2D &fxy) {
	for ( int j=0; j < 2; j++)
		for ( int i=0; i < nx ; i++)
			{fxy[i][j]=fxy[i][j+ny];fxy[i][2+j+ny]=fxy[i][2+j];}
}

void periodicboundary(Vector2D &fxy) {
	periodicboundary( (fxy.x) );
	periodicboundary( (fxy.y) );
	periodicboundary( (fxy.z) );
}

void periodicboundary(Field3D &fxyz) {
	//fxyz(x,y,zlength)=fxyz(x,y,0), for any x,y.
	for(int j=0;j<ny;j++)
		for(int i=0;i<nx;i++)
			fxyz[i][j+2][ngz-1] = fxyz[i][j+2][0];
	//
	BoutReal tmpz, *tmpf;
	tmpf = new BoutReal [ngz];
	for(int k=0;k<ngz;k++) {
		for(int i=0;i<nx;i++) {
			//fxyz[i][ny+2][k]=fxyz[i][2][k];
			tmpz=zgrid[k]-shiftangle[i];
			tmpz=tmpz - floor(tmpz/zlength)*zlength;
			for(int j=0;j<ngz;j++)
				tmpf[j]=fxyz[i][2][j];
			fxyz[i][ny+2][k]=interp1dlin (zgrid, ngz, tmpf, tmpz);
			//fxyz[i][ny+3][k]=fxyz[i][3][k];
			for(int j=0;j<ngz;j++)
				tmpf[j]=fxyz[i][3][j];
			fxyz[i][ny+3][k]=interp1dlin (zgrid, ngz, tmpf, tmpz);
			//fxyz[i][0][k]=fxyz[i][ny][k];
			tmpz=zgrid[k]+shiftangle[i];
			tmpz=tmpz - floor(tmpz/zlength)*zlength;
			for(int j=0;j<ngz;j++)
				tmpf[j]=fxyz[i][ny][j];
			fxyz[i][0][k]=interp1dlin (zgrid, ngz, tmpf, tmpz);
			//fxyz[i][1][k]=fxyz[i][ny+1][k];
			for(int j=0;j<ngz;j++)
				tmpf[j]=fxyz[i][ny+1][j];
			fxyz[i][1][k]=interp1dlin (zgrid, ngz, tmpf, tmpz);
		}
	}
}

void periodicboundary(Vector3D &fxyz) {
	periodicboundary( (fxyz.x) );
	periodicboundary( (fxyz.y) );
	periodicboundary( (fxyz.z) );
}

void smoothboundary(Field2D &fxy) {
	int nbound=5;
	BoutReal nearzero = 1.0e-12;
	for(int i=0;i<nx;i++)
	{for(int j=0;j<nbound;j++)
		{if (fabs(fxy[i][j]/(fxy[i][nbound]+nearzero)) > 10.0)
			{fxy[i][j] = fxy[i][nbound];}
		if (fabs(fxy[i][ngy-1-j]/(fxy[i][ngy-1-nbound]+nearzero)) > 10.0)
			{fxy[i][ngy-1-j] = fxy[i][ngy-1-nbound];}
		}
	}
}

void smoothboundary(Vector2D &fxy) {
	smoothboundary( (fxy.x) );
	smoothboundary( (fxy.y) );
	smoothboundary( (fxy.z) );
}
 
void smoothboundary(Field3D &fxyz) {
   ;
}

void smoothboundary(Vector3D &fxyz) {
	smoothboundary( (fxyz.x) );
	smoothboundary( (fxyz.y) );
	smoothboundary( (fxyz.z) );
}
