#include "cylinder.h"
#include "Csplines.cxx"
int pnpoly(int nvert, BoutReal *vertx, BoutReal *verty, BoutReal testx, BoutReal testy) {
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

BoutReal cylinder::dot(cylinder vec) {
	return (r*vec.r+p*vec.p+z*vec.z);
}

cylinder cylinder::cross(cylinder vec) {
	cylinder res;
	res.r = p*vec.z-z*vec.p;
	res.p = z*vec.r-r*vec.z;
	res.z = r*vec.p-p*vec.r;
	return res;
}

mesh1D::mesh1D() {
}

mesh1D::mesh1D(int in, BoutReal *x, int ik) {
	n = in;
	k = ik;
	grid = new BoutReal[n];
	knot = new BoutReal[n+k];
	for(int i=0;i<n;i++) {
		grid[i] = x[i];
	}
	cdbbsnak(grid,n,knot,n+k);
}

void mesh1D::init(int in, BoutReal *x, int ik) {
	n = in;
	k = ik;
	grid = new BoutReal[n];
	knot = new BoutReal[n+k];
	for(int i=0;i<n;i++) {
		grid[i] = x[i];
	}
	cdbbsnak(grid,n,knot,n+k);
}

scal1D::scal1D() {
}

/*
void particle::map() {
		if( pnpoly(nlcfs, r_lcfs, z_lcfs, y[0], y[1]) == 1) {
			region = 0;
			psi = psip.val(y[0],y[1],0,0,gmesh);
		}
		else if ( pnpoly(ndomain, r_domain, z_domain, y[0],y[1]) == 1) {
			psi = psip.val(y[0],y[1],0,0,gmesh);
			if (psi >= psib) locate[i] = 1;
			else region = 2;
		}
		else {
			region = 9;
		}
}
*/
scal1D::scal1D(mesh1D m) {
	f = new BoutReal[m.n];
}

void scal1D::init(mesh1D m) {
	f = new BoutReal[m.n];
}

void scal1D::coef(mesh1D m) {
	bscoef = new BoutReal[m.n];
	cdbbscoef(m.grid,f,m.n,m.knot,m.k,bscoef);
}

BoutReal scal1D::val(BoutReal x, int xderiv,mesh1D m) {
	return cdbbsval(m.knot, m.n + m.k, bscoef, m.n, x, xderiv);
}

meshRZ::meshRZ() {
}


meshRZ::meshRZ(int ir, BoutReal *vr, int iz, BoutReal *vz, int ikr, int ikz) {
	nr=ir;
	nz=iz;
	kr=ikr;
	kz=ikz;
	gridr = new BoutReal[nr];
	gridz = new BoutReal[nz];
	knotr = new BoutReal[nr+kr];
	knotz = new BoutReal[nz+kz];
	for(int i=0;i<nr;i++) {
		gridr[i] = vr[i];
	}
	for(int i=0;i<nz;i++) {
		gridz[i] = vz[i];
	}
	cdbbsnak(gridr,nr,knotr,nr+kr);
	cdbbsnak(gridz,nz,knotz,nz+kz);
}

void meshRZ::init(int ir, int iz, int ikr, int ikz) {
	nr=ir;
	nz=iz;
	kr=ikr;
	kz=ikz;
	gridr = new BoutReal[nr];
	gridz = new BoutReal[nz];
	knotr = new BoutReal[nr+kr];
	knotz = new BoutReal[nz+kz];
}


void meshRZ::init(int ir, BoutReal *vr, int iz, BoutReal *vz, int ikr, int ikz) {
	nr=ir;
	nz=iz;
	kr=ikr;
	kz=ikz;
	gridr = new BoutReal[nr];
	gridz = new BoutReal[nz];
	knotr = new BoutReal[nr+kr];
	knotz = new BoutReal[nz+kz];
	for(int i=0;i<nr;i++) {
		gridr[i] = vr[i];
	}
	for(int i=0;i<nz;i++) {
		gridz[i] = vz[i];
	}
	cdbbsnak(gridr,nr,knotr,nr+kr);
	cdbbsnak(gridz,nz,knotz,nz+kz);
}

scalRZ::scalRZ() {
}
/*
scalRZ::~scalRZ() {
	clearf();
	clearcoef();
}
*/
scalRZ::scalRZ(meshRZ mrz) {
	int nr=mrz.nr, nz=mrz.nz;

	f = new BoutReal*[nr];
	f[0] = new BoutReal[nr*nz];
	for (long i=1;i<nr;i++) {
		f[i] = f[i-1] + nz;
	}
}

void scalRZ::init(meshRZ mrz) {
	int nr=mrz.nr, nz=mrz.nz;

	f = new BoutReal*[nr];
	f[0] = new BoutReal[nr*nz];
	for (long i=1;i<nr;i++) {
		f[i] = f[i-1] + nz;
	}

	bscoef = new BoutReal*[mrz.nr];
	bscoef[0] = new BoutReal[mrz.nr*mrz.nz];
	for (long i=1;i<mrz.nr;i++) {
		bscoef[i] = bscoef[i-1] + mrz.nz;
	}
}

void scalRZ::coef(meshRZ mrz) {
	cdbbscoef2d(mrz.gridr,mrz.gridz,f,mrz.nr,mrz.nz,mrz.knotr,mrz.knotz,mrz.kr,mrz.kz,bscoef);
}

BoutReal scalRZ::val(BoutReal r, BoutReal z, int derivx, int derivy, meshRZ mrz) {
	return cdbbsval2d(mrz.knotr,mrz.nr+mrz.kr,mrz.knotz,mrz.nz+mrz.kz,bscoef,mrz.nr,mrz.nz,r,z,derivx,derivy);
}

void scalRZ::clearf() {
	delete[] f[0];
	delete f;
}

void scalRZ::clearcoef() {
	delete[] bscoef[0];
	delete bscoef;
}

vecRZ::vecRZ() {
}

vecRZ::vecRZ(meshRZ mrz) : vr(mrz),vp(mrz),vz(mrz) {
}
/*
vecRZ::~vecRZ() {
	clearf();
	clearcoef();
}
*/
void vecRZ::init(meshRZ mrz) {
	vr.init(mrz);
	vp.init(mrz);
	vz.init(mrz);
}

void vecRZ::coef(meshRZ mrz) {
	(vr.coef)(mrz);
	(vp.coef)(mrz);
	(vz.coef)(mrz);
}

cylinder vecRZ::val(BoutReal r, BoutReal z, int derivx, int derivy, meshRZ mrz) {
	cylinder res;
	res.r = (vr.val)(r, z, derivx, derivy, mrz);
	res.p = (vp.val)(r, z, derivx, derivy, mrz);
	res.z = (vz.val)(r, z, derivx, derivy, mrz);
	return res;
}

void vecRZ::clearf() {
	vr.clearf();
	vp.clearf();
	vz.clearf();
}

void vecRZ::clearcoef() {
	vr.clearcoef();
	vp.clearcoef();
	vz.clearcoef();
}
