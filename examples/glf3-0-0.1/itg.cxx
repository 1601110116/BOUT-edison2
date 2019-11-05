/****************************************************************
 '3+0' field gyro-fluid model for ITG simulation
 ****************************************************************/
#include <bout.hxx>
#include <boutmain.hxx>
#include <gyro_average.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <cmath>

// Solved Variables
Field3D Vor;  // Generalized Vorticity
Field3D Vpar; // Parallel Flow Velocity
Field3D Prs;  // Ion Pressure
Field3D Ti;   // Ion Temperature

// Derived quantities
Field3D Phi;  // Electrostatic potential
Field3D Ni;   // Ion Density

// Flux Surface Averages
Field2D VorHat,VparHat,PrsHat,TiHat;
Field2D NiHat,NeHat,TeHat;

// Equilibrium quantities
Vector2D B0vec; // Equilibrium B field vector
Vector2D b0xcv,b0xcv2; // curvature term
Field2D qsafe,rho;
Field2D qprofBout;

// [itg] options
bool nonlinear;    // Include nonlinear terms
bool pressure_eq;  // select pressure or temperature equation
bool filter_z;
int filter_z_mode;
int phi_flags;  
int profile;       // select cyclone (=0) or ITER (=1) case 
int prsrhs_term;   // select RHS term of pressure or Ti equation 
BoutReal viscosity;
BoutReal hypervis;   

BRACKET_METHOD bm = BRACKET_STD;
//BRACKET_METHOD bm = BRACKET_ARAKAWA;

// Fundamental constants
const BoutReal qe = 1.602e-19;      // Electron charge
const BoutReal Mp = 1.67262158e-27; // Proton mass
BoutReal AA = 1.0;                  // Ion atomic mass (1.0:Proton)
BoutReal ZZ = 1.0;                  // Ion charge

// Normalisation factors
BoutReal Lbar;   // Scale length = minor radius
BoutReal Tbar;   // Timescale = Lbar / Cs
BoutReal Ti0;    // Ion temperature at r=0 
BoutReal Bbar;   // Magnetic field at r=0

// Basic parameters
BoutReal Rmajor,aminor; // Major and minor radii
BoutReal epsilon; // Inversed aspect ratio
BoutReal Te0;     // Electron temperature at r=0
BoutReal Cs;      // Ion thermal speed at r=0 : sqrt(Ti0/Mi)
BoutReal rho_i;   // Ion Larmor Radius
BoutReal delta;   // rho*=rho_i/a in this model
BoutReal mu1,mu2; // Viscosity, hyper-viscosity

FieldGroup comms; // Communications

#define CYCLONE 0
#define ITER 1

const Field3D Curv(const Field3D &f);
int nx,ny;


////////////////////////////////////////////////////////////////////////
int physics_init(bool restarting)
{
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("itg");

  // Read [itg] options
  OPTION(options, phi_flags,  0); 
  OPTION(options, nonlinear, false);
  OPTION(options, pressure_eq, true);
  OPTION(options, filter_z, false);
  OPTION(options, filter_z_mode, 1);
  OPTION(options, profile, 0);
  OPTION(options, prsrhs_term, 1); 
  OPTION(options, viscosity, -1.0);  
  OPTION(options, hypervis, -1.0);  

  // Read grid file
  Field2D Rxy, Bpxy, Btxy, Bxy, hthe, I_shear;

  GRID_LOAD(Rxy);    // Major radius [m]
  GRID_LOAD(hthe);   // Poloidal arc length [m/radian]

  GRID_LOAD(Bpxy);   // Poloidal B field [T]
  GRID_LOAD(Btxy);   // Toroidal B field [T]
  GRID_LOAD(Bxy);    // Total B field [T]

  mesh->get(I_shear, "sinty");  // dnu/dpsi*theta

  // q profile 
  mesh->get(qsafe, "ShiftAngle");  // q*2*pi
  qsafe=qsafe/(2*M_PI);
  dump.add(qsafe,  "qsafe",  0);

  // local pitch
  qprofBout=Btxy*hthe/(Bpxy*Rxy);
  dump.add(qprofBout, "qprofBout", 0);

  // number of nodes
  mesh->get(nx, "nx"); 
  mesh->get(ny, "ny"); 

  // b0 x kappa
  b0xcv.covariant=false;
  mesh->get(b0xcv.x,  "bxcvx");
  mesh->get(b0xcv.y,  "bxcvy");
  mesh->get(b0xcv.z,  "bxcvz");

  // Assign 2D profiles directly
  BoutReal R_Lne,R_Lti,R_Lte,Te_T0,rmajor,rtmp,xx,deltar,r1,r2;
  BoutReal factorni,factorti,factorte;

    if(profile==CYCLONE){ // cyclone base case (SI units, Ti[eV])
       R_Lne=2.22; R_Lti=6.92;  R_Lte=R_Lti;  
       Rmajor=1.3; aminor=0.48; Bbar=1.9;
       Ti0=4.1e3; Te0=Ti0; deltar=0.3;

    }else if(profile==ITER){ // ITER-like case (SI units, Ti[eV])
       R_Lne=1.0/0.54; R_Lti=1.0/0.11; R_Lte=R_Lne; 
       Rmajor=2.6; aminor=0.94; Bbar=4.6;
       Ti0=26.e3; Te0=8.e3; deltar=0.3;
    }

    epsilon=aminor/Rmajor;
    rmajor=1.0/epsilon;  // normalized major radius, R/a

    mesh->get(PrsHat, "pressure");
    rho  =0.0*PrsHat;
    NiHat=0.0*PrsHat;
    TiHat=0.0*PrsHat;
    TeHat=0.0*PrsHat;

    r1=0.1; r2=0.9;
    for(int i=0;i<mesh->ngx;i++)
      for(int j=0;j<mesh->ngy;j++) {
        rtmp=(r2-r1)*((float)(mesh->XGLOBAL(i))/(float)nx)+r1;
        rho[i][j]=rtmp;

        xx=(rtmp-0.5)/deltar;
        NiHat[i][j]=exp(-deltar*R_Lne/rmajor*tanh(xx));
        TiHat[i][j]=exp(-deltar*R_Lti/rmajor*tanh(xx));
        TeHat[i][j]=exp(-deltar*R_Lte/rmajor*tanh(xx));
      }

    xx=-0.5/deltar;
    factorni=exp(-deltar*R_Lne/rmajor*tanh(xx)); // value at r=0
    factorti=exp(-deltar*R_Lti/rmajor*tanh(xx)); 
    factorte=exp(-deltar*R_Lte/rmajor*tanh(xx)); 

    NiHat=NiHat/factorni; // density is normalized by value at r=0
    TiHat=TiHat/factorti; // Ti is normalized by value at r=0
    TeHat=TeHat/factorte; 

    Te_T0=Te0/Ti0;
    TeHat=Te_T0*TeHat;    // Te is normalized by Ti0 

    PrsHat  = NiHat*TiHat;  
    NeHat   = NiHat;        

    dump.add(PrsHat,  "PrsHat", 0);
    dump.add(TiHat,   "TiHat",  0);
    dump.add(TeHat,   "TeHat",  0);
    dump.add(NiHat,   "NiHat",  0);
    dump.add(rho,     "rho", 0);

  // Basic parameters (SI units)
  Lbar=aminor;                  // Scale length 
  Cs = sqrt(qe*Ti0 / (AA*Mp));  // Ion thermal speed at r=0
  Tbar = Lbar / Cs;             // Timescale
  rho_i = Cs / (ZZ*qe*Bbar/(AA*Mp)); // Gyroradius 
  delta = rho_i / Lbar;  // Here, delta means rho*

  output << "\n\tParameters\n";
  output.write("\tLbar = %e [m], rho_i = %e [m]\n",Lbar, rho_i);
  output.write("\tTbar = %e [s], Cs = %e [m/s]\n", Tbar,Cs);
  output.write("\trho_star = %e, epsilon = %e \n", delta,epsilon);
  output.write("\tBbar = %e [T], Ti0 = %e [eV]\n", Bbar,Ti0);

  dump.add(rho_i,   "rho_i", 0);
  dump.add(delta,   "delta", 0);
  dump.add(Lbar,    "Lbar", 0);
  dump.add(Tbar,    "Tbar", 0);
  dump.add(Bbar,    "Bbar", 0);
  dump.add(Cs,      "Cs", 0); 

  // Normalize
  b0xcv.x /= Bbar;
  b0xcv.y *= Lbar*Lbar;
  b0xcv.z *= Lbar*Lbar;
  dump.add(b0xcv.x, "b0xcvx", 0);
  dump.add(b0xcv.y, "b0xcvy", 0);
  dump.add(b0xcv.z, "b0xcvz", 0);

  hthe /= Lbar; 
  Bpxy /= Bbar;
  Btxy /= Bbar;
  Bxy  /= Bbar;

  Rxy  /= Lbar; 
  mesh->dx /= Lbar*Lbar*Bbar;
  I_shear *= Lbar*Lbar*Bbar;

  // Metric components
  mesh->g11 = (Rxy^2)*(Bpxy^2);
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I_shear^2)*mesh->g11+(Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I_shear*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);

  mesh->J = hthe / Bpxy;
  mesh->Bxy = Bxy;

  mesh->g_11 = (I_shear^2)*(Rxy^2)+1.0/mesh->g11;
  mesh->g_22 = (Bxy^2)*(hthe^2)/(Bpxy^2);
  mesh->g_33 = Rxy*Bxy;
  mesh->g_23 = (Btxy*hthe*Rxy)/Bpxy;
  mesh->g_12 = I_shear*mesh->g_23;
  mesh->g_13 = I_shear*(Rxy^2);

  mesh->geometry();

  // Set B field vector  
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;

  // For comparison bet. kappa & grad lnB
  b0xcv2.covariant = false;
  b0xcv2=B0vec^Grad(log(mesh->Bxy))/mesh->Bxy;
  dump.add(b0xcv2.x, "b0xcv2x", 0);
  dump.add(b0xcv2.y, "b0xcv2y", 0);
  dump.add(b0xcv2.z, "b0xcv2z", 0);

  // Viscosity, hyper-viscosity
  mu1=viscosity*pow(delta,2.0);
  mu2=hypervis*pow(delta,4.0);

  // Add Equations to be solved
    SOLVE_FOR(Vor);   comms.add(Vor);
    SOLVE_FOR(Vpar);  comms.add(Vpar);
  if(pressure_eq){
    SOLVE_FOR(Prs);   comms.add(Prs);
  }else{
    SOLVE_FOR(Ti);    comms.add(Ti);
  }

  dump.add(Phi,"Phi",1);
  comms.add(Phi);

  Phi.setBoundary("Phi");

  cout << "init done. \n";
  return 0;
}


////////////////////////////////////////////////////////////////////////
int physics_run(BoutReal t)
{
  // Laplacian Inversion
  Field3D bInvLap;
  Field2D aInvLap;

  bInvLap = -Vor/NiHat/delta/delta;
  aInvLap = -1.0/TeHat/delta/delta;
  Phi = invert_laplace(bInvLap, phi_flags, &aInvLap);
  Phi.applyBoundary();

  mesh->communicate(comms);

  Ni = NeHat / TeHat * Phi ; // Adiabatic electron 

  if(pressure_eq)
    Ti = (Prs  - (Ni * TiHat)) / NiHat;
  else
    Prs  = Ni*TiHat + NiHat*Ti;

  // Vorticity Eq.
    ddt(Vor)  = NiHat*(-Div_par(Vpar) + Curv(Phi) + Curv(Prs)/NiHat);
    ddt(Vor) += delta*bracket(PrsHat,(Ni-Vor)/NiHat,bm);
    ddt(Vor) += delta*bracket(NiHat,Phi,bm);

  // Parallel Velocity Eq.
    ddt(Vpar) = -Grad_par(Phi)-Grad_par(Prs)/NiHat;

  // Pressure or Temperature Eq.
  if(pressure_eq){
    if(prsrhs_term==1)
      ddt(Prs)  = 5.0/3.0*PrsHat*(-Div_par(Vpar));
    else if(prsrhs_term==2)
      ddt(Prs)  = 5.0/3.0*PrsHat*(-Div_par(Vpar) + Curv(Ti));
    else if(prsrhs_term==3)
      ddt(Prs)  = 5.0/3.0*PrsHat*(-Div_par(Vpar) + Curv(Phi+Ti) + Curv(Prs)/NiHat);
  }else{
    if(prsrhs_term==1)
      ddt(Ti)   = 2.0/3.0*TiHat*(-Div_par(Vpar));
    else if(prsrhs_term==2)
      ddt(Ti)   = 2.0/3.0*TiHat*(-Div_par(Vpar) + Curv(2.5*Ti));
    else if(prsrhs_term==3)
      ddt(Ti)   = 2.0/3.0*TiHat*(-Div_par(Vpar) + Curv(Phi+2.5*Ti) + Curv(Prs)/NiHat);
  }

  // Viscous damping 
  if(viscosity > 0.0){
      ddt(Vor) += mu1 * Delp2(Vor);
      ddt(Vpar)+= mu1 * Delp2(Vpar);
    if(pressure_eq)
      ddt(Prs) += mu1 * Delp2(Prs);
    else
      ddt(Ti)  += mu1 * Delp2(Ti);
  }

  // Hyper-viscous damping 
  if(hypervis > 0.0){
      ddt(Vor) -= mu2 * Delp2(Delp2(Vor));
      ddt(Vpar)-= mu2 * Delp2(Delp2(Vpar));
    if(pressure_eq)
      ddt(Prs) -= mu2 * Delp2(Delp2(Prs));
    else
      ddt(Ti)  -= mu2 * Delp2(Delp2(Ti));
  }

  // ExB advection
  if(nonlinear){
      ddt(Vor) -= delta*delta*bracket(Phi,Vor,bm);
      ddt(Vpar)-= delta*delta*bracket(Phi,Vpar,bm);
    if(pressure_eq)
      ddt(Prs) -= delta*delta*bracket(Phi,Prs,bm);
    else
      ddt(Ti)  -= delta*delta*bracket(Phi,Ti,bm);
  }else{
    if(pressure_eq)
      ddt(Prs) += delta*bracket(PrsHat,Phi,bm); // source for linear ITG 
    else
      ddt(Ti)  += delta*bracket(TiHat,Phi,bm);  // source for linear ITG
  }

  // Filtering
  if(filter_z){
      ddt(Vor)  = filter(ddt(Vor),  filter_z_mode);
      ddt(Vpar) = filter(ddt(Vpar), filter_z_mode);
    if(pressure_eq)
      ddt(Prs)  = filter(ddt(Prs),  filter_z_mode);
    else
      ddt(Ti)   = filter(ddt(Ti),   filter_z_mode);
  }

  return 0;
}

// Curvature operator 
const Field3D Curv(const Field3D &f)
{
  return -delta*(b0xcv*Grad(f)/mesh->Bxy+bracket(log(mesh->Bxy), f, bm));
}

