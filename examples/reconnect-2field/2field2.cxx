/*****************************************************************************
 * 2 field (Apar, vorticity) model for benchmarking 
 * simple slab reconnection model
 *****************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <bout/constants.hxx>

#include <invert_laplace.hxx>
#include <initialprofiles.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Jpar0, Te0, Ni0;

// 3D evolving fields
Field3D Upar, Apar;

// Derived 3D variables
Field3D Phi, Jpar;

// External coil field
Field3D Apar_ext, Jpar_ext, Phi0_ext, Upar0_ext;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe;

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.60217646e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.67262158e-27; // Ion mass
const BoutReal Me = 9.1093816e-31;  // Electron mass
const BoutReal Me_Mi = Me / Mi; // Electron mass / Ion mass

// normalisation parameters
int NOUT;
BoutReal TIMESTEP;
BoutReal length_norm, time_norm, speed_norm, timestep_norm;
BoutReal Te_norm, Ne_norm, B0_norm;
BoutReal Cs, rho_s, wci, beta_hat;

BoutReal etab0x, etab2x, muix2x, muix4b;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

int phi_flags; // Inversion flags

bool nonlinear;
bool parallel_lc;
bool include_jpar0;
int  jpar_bndry;

// Function Declarations
void smooth_bndry(Field3D f, int bndry);
int physics_init(bool restarting);
int physics_run(BoutReal t);

// Function Definitions

int physics_init(bool restarting) {

  // Load 2D profiles
  GRID_LOAD3(Jpar0, Te0, Ni0);
  Ni0 *= 1e20; // To m^-3
  
  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->Bxy,  "Bxy");
  
  // Read some parameters
  Options *globalOptions = Options::getRoot();

  // Global Options
  OPTION(globalOptions, NOUT, 1);        // Total # of time steps
  OPTION(globalOptions, TIMESTEP, 1.0);  // Normalised time/step  
  SAVE_ONCE2(NOUT,TIMESTEP);

  Options *options = globalOptions->getSection("2field");  
  // normalisation values
  OPTION(options, nonlinear, false);
  OPTION(options, parallel_lc, true);
  OPTION(options, include_jpar0, true);
  OPTION(options, jpar_bndry, 0);

  OPTION(options, etab0x,  1.e-3); // Normalised parallel resistivity             
  OPTION(options, etab2x,  0.e0 ); // Normalised par hyper-resistivity      Delp2
  OPTION(options, muix2x,  1.e-3); // Normalised perp ion viscosity         Delp2
  OPTION(options, muix4b,  0.e0 ); // Normalised perp ion hyper-viscosity   Dpar4
  SAVE_ONCE4(etab0x,etab2x,muix2x,muix4b);

  OPTION(options, phi_flags,   0);
  
  int bracket_method;
  OPTION(options, bracket_method, 0);
  switch(bracket_method) {
  case 0: {
    bm = BRACKET_STD; 
    output << "\tBrackets: default differencing\n";
    break;
  }
  case 1: {
    bm = BRACKET_SIMPLE; 
    output << "\tBrackets: simplified operator\n";
    break;
  }
  case 2: {
    bm = BRACKET_ARAKAWA; 
    output << "\tBrackets: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm = BRACKET_CTU; 
    output << "\tBrackets: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

  ///////////////////////////////////////////////////
  // Normalisation
  
  Te_norm = max(Te0, true);
  if(Te_norm < 1)
    Te_norm = 1000;
  Ne_norm = max(Ni0, true);
  if(Ne_norm < 1)
    Ne_norm = 1.e19;
  B0_norm  = max(mesh->Bxy, true);
  
  // Sound speed in m/s
  Cs = sqrt(Charge*Te_norm / Mi);
  speed_norm = Cs;

  // drift scale
  rho_s = Cs * Mi / (Charge * B0_norm);
  length_norm = rho_s;

  // Ion cyclotron frequency
  wci = Charge * B0_norm / Mi;
  time_norm = 1/wci;
  timestep_norm = TIMESTEP*time_norm;

  beta_hat = MU0 * Charge*Te_norm * Ne_norm / (B0_norm*B0_norm);
  
  output << "\tNormalisations:" << endl;
  output << "\tCs       = " << Cs << endl;
  output << "\trho_s    = " << rho_s << endl;
  output << "\twci      = " << wci << endl; 
  output << "\tbeta_hat = " << beta_hat << endl; 
  
  SAVE_ONCE3(Te_norm, Ne_norm, B0_norm);
  SAVE_ONCE4(Cs, rho_s, wci, beta_hat);
  SAVE_ONCE4(speed_norm, length_norm, time_norm, timestep_norm);

  // Normalise geometry 
  Rxy  /= rho_s;
  hthe /= rho_s;
  mesh->dx /= rho_s*rho_s*B0_norm;

  // Normalise magnetic field
  Bpxy /= B0_norm;
  Btxy /= B0_norm;
  mesh->Bxy /= B0_norm;
  
  // Plasma quantities
  Jpar0 /= Ne_norm*Charge*Cs;

  // CALCULATE METRICS

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = 0.;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11;
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = 0.;
  mesh->g_13 = 0.;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;

  mesh->geometry();

  // Tell BOUT++ which variables to evolve
  SOLVE_FOR2(Upar, Apar);
  
  // Set boundary conditions
  Jpar.setBoundary("Jpar");
  Phi.setBoundary("Phi");
  
  // Add any other variables to be dumped to file
  SAVE_REPEAT2(Phi, Jpar);
  SAVE_ONCE(Jpar0);
  
  // Generate external field
  
  initial_profile("Apar_ext", Apar_ext);
  Jpar_ext = -Delp2(Apar_ext);
  SAVE_ONCE2(Apar_ext,Jpar_ext);

  initial_profile("Phi0_ext", Phi0_ext);
  Upar0_ext = Delp2(Phi0_ext)/mesh->Bxy;
  SAVE_ONCE2(Phi0_ext,Upar0_ext);
  
  return 0;
}

const Field3D Grad_parP_LtoC(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_LtoC(f);
//    if(nonlinear) {
      result -= beta_hat * bracket(Apar_ext+Apar, f, BRACKET_ARAKAWA);
//    }else {
//      result -= beta_hat * bracket(Apar_ext, f, BRACKET_ARAKAWA);
//    }
  }else {
//    if(nonlinear) {
      result = Grad_parP((Apar+Apar_ext)*beta_hat, f);
//    }else {
//      result = Grad_parP(Apar_ext*beta_hat, f);
//    }
  }
  return result;
}

const Field3D Grad_parP_CtoL(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_CtoL(f);
//    if(nonlinear) {
      result -= beta_hat * bracket(Apar+Apar_ext, f, BRACKET_ARAKAWA);
//    }else {
//      result -= beta_hat * bracket(Apar_ext, f, BRACKET_ARAKAWA);
//    }
  }else {
//    if(nonlinear) {
        result = Grad_parP((Apar+Apar_ext)*beta_hat, f);
//    }else {
//      result = Grad_parP(Apar_ext*beta_hat, f);
//    }
  }
  return result;
}

int physics_run(BoutReal t) {
  // Solve EM fields

  // Upar = (1/B) * Delp2(Phi)
  Phi = invert_laplace(mesh->Bxy*Upar, phi_flags);
  Phi.applyBoundary(); // For target plates only
  Upar.applyBoundary();

  mesh->communicate(Upar, Phi, Apar);
  
  Jpar = -Delp2(Apar+Apar_ext);
  Jpar.applyBoundary();
  mesh->communicate(Jpar);
  if(jpar_bndry > 0) 
    smooth_bndry(Jpar,jpar_bndry);
/*  if(jpar_bndry > 0) {
    // Boundary in jpar
    if(mesh->firstX()) {
      for(int i=jpar_bndry;i>=0;i--)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    Jpar[i][j][k] = Jpar[i+1][j][k];
	  }
    }
    if(mesh->lastX()) {
      for(int i=mesh->ngx-jpar_bndry-1;i<mesh->ngx;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    Jpar[i][j][k] = Jpar[i-1][j][k];
	  }
    }
  }
*/

  // VORTICITY
  ddt(Upar) = 0.0;

    // Parallel current
  if (!nonlinear)
    ddt(Upar) += SQ(mesh->Bxy)*Grad_par_LtoC(Jpar/mesh->Bxy);
  else
    ddt(Upar) += SQ(mesh->Bxy)*Grad_parP_LtoC(Jpar/mesh->Bxy);

  if(include_jpar0) {
   ddt(Upar) -= SQ(mesh->Bxy)*beta_hat * bracket(Apar+Apar_ext, Jpar0/mesh->Bxy, BRACKET_ARAKAWA);
  }
 
    // ExB advection  
    ddt(Upar) -= bracket(Phi0_ext, Upar, bm); 
    //ddt(Upar) -= bracket(Phi, Upar0_ext, bm);   
  if(nonlinear) {
    ddt(Upar) -= bracket(Phi, Upar, bm);  
  }

    // Dissipation
  if(muix2x > 0.)
    ddt(Upar) += muix2x*Delp2(Upar);

/*
  if(numvisc > 0) {
    //mu4 = 2.0*mesh->DY*hthe/rhos;
    //mu4*=mu4; //2nd power
    //mu4*=mu4; //4th power
    Field3D Dpar2U;
    Dpar2U = Grad2par2(U);
    Dpar2U.applyBoundary('neumann');
    ddt(Upar) -= muix4b*Grad2par2(Dpar2U);
  }
*/

  // APAR
  ddt(Apar) = 0.0;

    // E_parallel
    ddt(Apar) -= Grad_parP_CtoL(Phi0_ext) / beta_hat;
  if (!nonlinear)
    ddt(Apar) -= Grad_par_CtoL(Phi) / beta_hat;
  else
    ddt(Apar) -= Grad_parP_CtoL(Phi) / beta_hat;

    // Dissipation
  if(etab0x > 0.)
    ddt(Apar) -= etab0x*Jpar / beta_hat;

  if(etab2x > 0.) {
    Field3D Delp2J;
    mesh->communicate(Jpar);
    Delp2J = Delp2(Jpar);
    ddt(Apar) += etab2x*Delp2J / beta_hat;
  }

  return 0;
}

void smooth_bndry(Field3D f, int bndry)
{
  // Smooth boundary values of f
  if(mesh->firstX()) {
    for(int i=bndry;i>=0;i--)
	for(int j=0;j<mesh->ngy;j++)
	for(int k=0;k<mesh->ngz-1;k++)
	    f[i][j][k] = f[i+1][j][k];
  }
  if(mesh->lastX()) {
    for(int i=mesh->ngx-bndry-1;i<mesh->ngx;i++)
	for(int j=0;j<mesh->ngy;j++)
    for(int k=0;k<mesh->ngz-1;k++)
      f[i][j][k] = f[i-1][j][k];
  }
}
