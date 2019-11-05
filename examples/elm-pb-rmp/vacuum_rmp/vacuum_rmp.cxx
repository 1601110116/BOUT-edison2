





/*******************************************************************************
 * High-Beta Flute-Reduced MHD
 * see elm_reduced.pdf
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * Can also include the Vpar compressional term
 *******************************************************************************/

#include "bout.hxx"
#include "initialprofiles.hxx"
#include "invert_laplace.hxx"
#include "interpolation.hxx"
#include "derivs.hxx"
#include <math.h>
#include "sourcex.hxx"
#include <boutmain.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>


//*******************************************gui2
#include "mpi.h"
#include <fstream>
#include <stdio.h>
#include <bout_types.hxx>
#include <globals.hxx>
#include <options.hxx>
#include <fft.hxx>
#include <utils.hxx>
#include <dcomplex.hxx>
#include <cmath>
#include <output.hxx>
#include <lapack_routines.hxx> // Tridiagonal & band inversion routines
#include <boutexception.hxx>
#include <bout/sys/timer.hxx>

// This was defined in nvector.h
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
int MYPE_IN_CORE=0;
int laplace_maxmode; ///< The maximum Z mode to solve for
bool invert_async_send; ///< If true, use asyncronous send in parallel algorithms
bool invert_use_pdd; ///< If true, use PDD algorithm
bool invert_low_mem;    ///< If true, reduce the amount of memory used
bool laplace_all_terms; // applies to Delp2 operator and laplacian inversion
bool laplace_nonuniform; // Non-uniform mesh correction
//*******************************************end gui2




// 2D inital profiles
Field2D J0, P0; // Current and pressure
Vector2D b0xcv; // Curvature term
Field2D beta, gradparB;   // Used for Vpar terms
Field2D phi0;   // When diamagnetic terms used
Field2D U0, Psixy, x;   //0th vorticity of equilibrium flow, 
                       //radial flux coordinate, normalized radial flux coordinate

// B field vectors
Vector2D B0vec; // B0 field vector

// V0 field vectors
Vector2D V0net; //net flow 

// 3D evolving variables
Field3D U, Psi, P, Vpar;

// Derived 3D variables
Field3D Jpar, phi; // Parallel current, electric potential

Field3D Jpar2; //  Delp2 of Parallel current

Field3D tmpP2; // Grad2_par2new of pressure
Field3D tmpU2; // Grad2_par2new of Parallel vorticity
Field3D tmpA2; // Grad2_par2new of Parallel vector potential

// Constraint
Field3D C_phi;

// Parameters
BoutReal density; // Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
BoutReal dnorm; // For diamagnetic terms: 1 / (2. * wci * Tbar)
BoutReal dia_fact; // Multiply diamagnetic term by this
BoutReal delta_i; // Normalized ion skin depth

BoutReal diffusion_p4;   //xqx: parallel hyper-viscous diffusion for pressure
BoutReal diffusion_u4;   //xqx: parallel hyper-viscous diffusion for vorticity
BoutReal diffusion_a4;   //xqx: parallel hyper-viscous diffusion for vector potential

BoutReal diffusion_par;  // Parallel pressure diffusion
BoutReal heating_P;  // heating power in pressure
BoutReal hp_width;  // heating profile radial width in pressure
BoutReal hp_length;  // heating radial domain in pressure
BoutReal sink_P;     // sink in pressure
BoutReal sp_width;   // sink profile radial width in pressure
BoutReal sp_length;  // sink radial domain in pressure

BoutReal sink_Ul;     // left edge sink in vorticity
BoutReal su_widthl;   // left edge sink profile radial width in vorticity
BoutReal su_lengthl;  // left edge sink radial domain in vorticity

BoutReal sink_Ur;     // right edge sink in vorticity
BoutReal su_widthr;   // right edge sink profile radial width in vorticity
BoutReal su_lengthr;  // right edge sink radial domain in vorticity

BoutReal viscos_par;  // Parallel viscosity
BoutReal viscos_perp; // Perpendicular viscosity
BoutReal hyperviscos; // Hyper-viscosity (radial)
Field3D hyper_mu_x; // Hyper-viscosity coefficient

// options
bool include_curvature, include_jpar0, compress0;
bool evolve_pressure;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar; 
BoutReal g; // Only if compressible
bool phi_curv;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE      
/*
 * Bracket method 
 *
 * BRACKET_STD      - Same as b0xGrad_dot_Grad, methods in BOUT.inp
 * BRACKET_SIMPLE   - Subset of terms, used in BOUT-06
 * BRACKET_ARAKAWA  - Arakawa central differencing (2nd order)
 * BRACKET_CTU      - 1st order upwind method
 *
 */

// Bracket method for advection terms 
BRACKET_METHOD bm_exb;
BRACKET_METHOD bm_mag;
int bm_exb_flag;
int bm_mag_flag;
/* BRACKET_METHOD bm_ExB = BRACKET_STD;
   BRACKET_METHOD bm_mflutter = BRACKET_STD; */

bool diamag;
bool diamag_grad_t; // Grad_par(Te) term in Psi equation
bool diamag_phi0;   // Include the diamagnetic equilibrium phi0

bool eHall;
BoutReal AA; // ion mass in units of the proton mass; AA=Mi/Mp

//net flow, Er=-R*Bp*Dphi0,Dphi0=-D_min-0.5*D_0*(1.0-tanh(D_s*(x-x0)))
Field2D V0; //net flow amplitude
Field2D Dphi0; //differential potential to flux
BoutReal D_0;   // potential amplitude
BoutReal D_s;   // shear parameter
BoutReal x0;   //velocity peak location
BoutReal sign; //direction of flow
BoutReal Psiaxis, Psibndry;  
bool withflow; 
bool K_H_term; //Kelvin-Holmhotz term
Field2D perp;     //for test    
BoutReal D_min; //constant in flow

//for C_mod
bool experiment_Er;  //read in total Er from experiment

bool nogradparj;
bool filter_z;
int filter_z_mode;
int low_pass_z;
int zonal_flow;
int zonal_field;
int zonal_bkgd;
bool relax_j_vac;
BoutReal relax_j_tconst; // Time-constant for j relax
Field3D Psitarget;   // The (moving) target to relax to

bool smooth_j_x;  // Smooth Jpar in the x direction

int jpar_bndry_width; // Zero jpar in a boundary region

bool parallel_lr_diff; // Use left and right shifted stencils for parallel differences

bool parallel_lagrange; // Use (semi-) Lagrangian method for parallel derivatives
bool parallel_project;  // Use Apar to project field-lines

Field3D Xip_x, Xip_z;     // Displacement of y+1 (in cell index space)

Field3D Xim_x, Xim_z;     // Displacement of y-1 (in cell index space)

bool phi_constraint; // Solver for phi using a solver constraint 

bool include_rmp; // Include RMP coil perturbation
bool simple_rmp;  // Just use a simple form for the perturbation
int rmp_n, rmp_m; // toroidal and poloidal mode numbers
BoutReal rmp_polwid;  // Poloidal width (-ve -> full, fraction of 2pi)
BoutReal rmp_polpeak; // Peak poloidal location (fraction of 2pi)
BoutReal rmp_factor;  // Multiply amplitude by this factor
BoutReal rmp_ramp;    // Ramp-up time for RMP [s]. negative -> instant
BoutReal rmp_freq;    // Amplitude oscillation frequency [Hz] (negative -> no oscillation)
BoutReal rmp_rotate;  // Rotation rate [Hz]
Field3D rmp_Psi0; // Parallel vector potential from Resonant Magnetic Perturbation (RMP) coils
Field3D rmp_Psi;  // Value used in calculations
Field3D rmp_Psix;
Field3D rmp_Psib;
Field3D rmp_dApdt; // Time variation
Vector3D rmp_B; // RMP field

BoutReal vac_lund, core_lund;       // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist,  core_resist;  // The resistivities (just 1 / S)
Field3D eta;                    // Resistivity profile (1 / S)
bool spitzer_resist;  // Use Spitzer formula for resistivity
BoutReal Zeff;            // Z effective for resistivity formula

BoutReal hyperresist;    // Hyper-resistivity coefficient (in core only)
BoutReal ehyperviscos;   // electron Hyper-viscosity coefficient
Field3D hyper_eta_x; // Radial resistivity profile
Field3D hyper_eta_z; // Toroidal resistivity profile

int damp_width;     // Width of inner damped region
BoutReal damp_t_const;  // Timescale of damping

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, B0, hthe;
Field2D I; // Shear factor

const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass

// Communication objects
FieldGroup comms;

int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner
int jacobian(BoutReal t); // Jacobian-vector multiply

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner with phi constraint

void advect_tracer(const Field3D &p,  // phi (input)
		   const Field3D &delta_x, const Field3D &delta_z, // Current location (input) 
		   Field3D &F_dx, Field3D &F_dz); // Time-derivative of location

const Field3D Grad2_par2new(const Field3D &f); //for 4th order diffusion


//*********************gui2
const Field3D invert_laplace_gui2(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c = NULL, const Field2D *d = NULL);

int invert_laplace_gui3(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c = NULL, const Field2D *d = NULL);

int invert_laplace_gui4(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c = NULL, const Field2D *d = NULL);

int invert_laplace_ser(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a,
                       const Field2D *ccoef=NULL, const Field2D *d=NULL);

void laplace_tridag_coefs_gui2(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c,
			       const Field2D *ccoef, const Field2D *d);
//*********************end_gui2

const Field3D Grad2_par2new(const Field3D &f)
{
  /*
   * This function implements d2/dy2 where y is the poloidal coordinate theta
   */


#ifdef CHECK
  int msg_pos = msg_stack.push("Grad2_par2new( Field3D )");
#endif



  Field3D result = D2DY2(f);
  

#ifdef TRACK
  result.name = "Grad2_par2new("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

int physics_init(bool restarting)
{
  bool noshear;
  
  output.write("Solving high-beta flute reduced equations\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  mesh->get(P0, "pressure"); // Pascals

  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Load metrics
  if(mesh->get(Rxy,  "Rxy")) { // m
    output.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1
  mesh->get(Psixy, "psixy");//get Psi             
  mesh->get(Psiaxis,"psi_axis");//axis flux      
  mesh->get(Psibndry,"psi_bndry");//edge flux   

  //////////////////////////////////////////////////////////////
  // Read parameters from the options file
  // 
  // Options.get ( NAME,    VARIABLE,    DEFAULT VALUE)
  //
  // or if NAME = "VARIABLE" then just
  //
  // OPTION(VARIABLE, DEFAULT VALUE) 
  //
  // Prints out what values are assigned
  /////////////////////////////////////////////////////////////

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("highbeta");

  OPTION(options, density,           1.0e19); // Number density [m^-3]

  OPTION(options, evolve_jpar,       false);  // If true, evolve J raher than Psi
  OPTION(options, phi_constraint,    false);  // Use solver constraint for phi

  // Effects to include/exclude
  OPTION(options, include_curvature, true);
  OPTION(options, include_jpar0,     true);
  OPTION(options, nogradparj,       false);  // exclude (B0^2)*Grad_par(Jpar)?
  OPTION(options, evolve_pressure,   true);
  
  OPTION(options, compress0,          false);
  OPTION(options, nonlinear,         false);

  // option for ExB Poisson Bracket 
  OPTION(options, bm_exb_flag,         0);
  switch(bm_exb_flag) {
  case 0: {
    bm_exb = BRACKET_STD;
    output << "\tBrackets for ExB: default differencing\n";
    break;
  }
  case 1: {
    bm_exb = BRACKET_SIMPLE;
    output << "\tBrackets for ExB: simplified operator\n";
    break;
  }
  case 2: {
    bm_exb = BRACKET_ARAKAWA;
    output << "\tBrackets for ExB: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm_exb = BRACKET_CTU;
    output << "\tBrackets for ExB: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

  // option for magnetic flutter Poisson Bracket
  OPTION(options, bm_mag_flag,         0);
  switch(bm_mag_flag) {
  case 0: {
    bm_mag = BRACKET_STD;
    output << "\tBrackets: default differencing\n";
    break;
  }
  case 1: {
    bm_mag = BRACKET_SIMPLE;
    output << "\tBrackets: simplified operator\n";
    break;
  }
  case 2: {
    bm_mag = BRACKET_ARAKAWA;
    output << "\tBrackets: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm_mag = BRACKET_CTU;
    output << "\tBrackets: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

  OPTION(options, eHall,            false);  // electron Hall or electron parallel pressue gradient effects? 
  OPTION(options, AA,               1.0);    // ion mass in units of proton mass

  OPTION(options, diamag,            false);  // Diamagnetic effects? 
  OPTION(options, diamag_grad_t,     diamag); // Grad_par(Te) term in Psi equation
  OPTION(options, diamag_phi0,       diamag); // Include equilibrium phi0
  OPTION(options, dia_fact,          1.0);    // Scale diamagnetic effects by this factor
  
  OPTION(options, withflow,          false);    //withflow or not 
  OPTION(options, K_H_term,          true);    //keep K-H term 
  OPTION(options, D_0,                0.0);    // velocity magnitude  
  OPTION(options, D_s,                0.0);    // flowshear    
  OPTION(options, x0,                 0.0);    //flow location   
  OPTION(options, sign,               1.0);    //flow direction, -1 means negative electric field   
  OPTION(options, D_min,           3000.0);    //a constant   

  OPTION(options, experiment_Er,     false);

  OPTION(options, noshear,           false);

  OPTION(options, relax_j_vac,       false); // Relax vacuum current to zero
  OPTION(options, relax_j_tconst,    0.1);
  
  // Toroidal filtering
  OPTION(options, filter_z,          false);  // Filter a single n
  OPTION(options, filter_z_mode,     1);
  OPTION(options, low_pass_z,       -1);      // Low-pass filter
  OPTION(options, zonal_flow,       -1);      // zonal flow filter
  OPTION(options, zonal_field,      -1);      // zonal field filter
  OPTION(options, zonal_bkgd,       -1);      // zonal background P filter

  // Radial smoothing
  OPTION(options, smooth_j_x,       false);  // Smooth Jpar in x

  // Jpar boundary region
  OPTION(options, jpar_bndry_width, -1);

  // Parallel differencing
  OPTION(options, parallel_lr_diff, false);
  OPTION(options, parallel_lagrange, false); // Use a (semi-) Lagrangian method for Grad_parP
  OPTION(options, parallel_project, false);

  // RMP-related options
  OPTION(options, include_rmp,       false);  // Read RMP data from grid
  
  OPTION(options, simple_rmp,        false);  // Include a simple RMP model
  OPTION(options, rmp_factor,         1.0);
  OPTION(options, rmp_ramp,          -1.0);
  OPTION(options, rmp_freq,          -1.0);
  OPTION(options, rmp_rotate,         0.0);

  // Vacuum region control
  OPTION(options, vacuum_pressure,   0.02);   // Fraction of peak pressure
  OPTION(options, vacuum_trans,      0.005);  // Transition width in pressure
  
  // Resistivity and hyper-resistivity options
  OPTION(options, vac_lund,          0.0);    // Lundquist number in vacuum region
  OPTION(options, core_lund,         0.0);    // Lundquist number in core region
  OPTION(options, hyperresist,       -1.0);
  OPTION(options, ehyperviscos,      -1.0);
  OPTION(options, spitzer_resist,    false);  // Use Spitzer resistivity
  OPTION(options, Zeff,              2.0);    // Z effective

  // Inner boundary damping
  OPTION(options, damp_width,        0);
  OPTION(options, damp_t_const,      0.1);

  // Viscosity and hyper-viscosity
  OPTION(options, viscos_par,        -1.0);  // Parallel viscosity
  OPTION(options, viscos_perp,       -1.0);  // Perpendicular viscosity
  OPTION(options, hyperviscos,       -1.0);  // Radial hyperviscosity
  
  // parallel pressure diffusion
  OPTION(options, diffusion_par,        -1.0);  // Parallel pressure diffusion
  OPTION(options, diffusion_p4,        -1.0);   //xqx: parallel hyper-viscous diffusion for pressure
  OPTION(options, diffusion_u4,        -1.0);   //xqx: parallel hyper-viscous diffusion for vorticity
  OPTION(options, diffusion_a4,        -1.0);   //xqx: parallel hyper-viscous diffusion for vector potential

  // heating factor in pressure
  OPTION(options, heating_P,        -1.0);  //  heating power in pressure
  OPTION(options, hp_width,         0.1);  //  the percentage of radial grid points for heating profile radial width in pressure
  OPTION(options, hp_length,        0.04);  //  the percentage of radial grid points for heating profile radial domain in pressure

  // sink factor in pressure
  OPTION(options, sink_P,           -1.0);  //  sink in pressure
  OPTION(options, sp_width,         0.05);  //  the percentage of radial grid points for sink profile radial width in pressure
  OPTION(options, sp_length,        0.04);  //  the percentage of radial grid points for sink profile radial domain in pressure


  // left edge sink factor in vorticity
  OPTION(options, sink_Ul,           -1.0);  //  left edge sink in vorticity
  OPTION(options, su_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in vorticity
  OPTION(options, su_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in vorticity

  // right edge sink factor in vorticity
  OPTION(options, sink_Ur,           -1.0);  //  right edge sink in vorticity
  OPTION(options, su_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in vorticity
  OPTION(options, su_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in vorticity

  // Compressional terms
  OPTION(options, phi_curv,          true);
  options->get("gamma",             g,                 5.0/3.0);
  
  // Field inversion flags
  OPTION(options, phi_flags,         0);
  OPTION(options, apar_flags,        0);

  x=(Psixy-Psiaxis)/(Psibndry-Psiaxis);     
 
  if(experiment_Er){          //get er from experiment
    mesh->get(Dphi0,"Epsi");  
    diamag_phi0=false;
    K_H_term=false;
  } else {
    Dphi0=-D_min-0.5*D_0*(1.0-tanh(D_s*(x-x0)));
  }

  if(sign<0)                //change flow direction
    Dphi0*=-1;    

  V0=-Rxy*Bpxy*Dphi0/B0;      
 
































  if(!include_curvature)
    b0xcv = 0.0;
  
  if(!include_jpar0)
    J0 = 0.0;

  if(noshear) {
    if(include_curvature)
      b0xcv.z += I*b0xcv.x;
    mesh->ShiftXderivs = false;
    I = 0.0;
  }
  
  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      if(include_curvature)
	b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
      output.write("\t****** I=0\n");
    }
  }

  //////////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES
  
  if(mesh->get(Bbar, "bmag")) // Typical magnetic field
    Bbar = 1.0;
  if(mesh->get(Lbar, "rmag")) // Typical length scale
    Lbar = 1.0;
  
  Va = sqrt(Bbar*Bbar / (MU0*density*Mi));

  Tbar = Lbar / Va;

  dnorm = dia_fact * Mi / (2.*1.602e-19*Bbar*Tbar);

  delta_i = AA*60.67*5.31e5/sqrt(density/1e6)/(Lbar*100.0);

  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Va = %e m/s   Tbar = %e s\n", Va, Tbar);
  output.write("                dnorm = %e\n", dnorm);
  output.write("    Resistivity\n");

  if(eHall)
    output.write("                delta_i = %e   AA = %e \n", delta_i, AA);

  if(vac_lund > 0.0) {
    output.write("        Vacuum  Tau_R = %e s   eta = %e Ohm m\n", vac_lund * Tbar, 
		 MU0 * Lbar * Lbar / (vac_lund * Tbar));
    vac_resist = 1. / vac_lund;
  }else {
    output.write("        Vacuum  - Zero resistivity -\n");
    vac_resist = 0.0;
  }
  if(core_lund > 0.0) {
    output.write("        Core    Tau_R = %e s   eta = %e Ohm m\n", core_lund * Tbar,
		 MU0 * Lbar * Lbar / (core_lund * Tbar));
    core_resist = 1. / core_lund;
  }else {
    output.write("        Core    - Zero resistivity -\n");
    core_resist = 0.0;
  }

  if(hyperresist > 0.0) {
    output.write("    Hyper-resistivity coefficient: %e\n", hyperresist);
    dump.add(hyper_eta_x, "hyper_eta_x", 1);
    dump.add(hyper_eta_z, "hyper_eta_z", 1);
  }

  if(ehyperviscos > 0.0) {
    output.write("    electron Hyper-viscosity coefficient: %e\n", ehyperviscos);
  }

  if(hyperviscos > 0.0) {
    output.write("    Hyper-viscosity coefficient: %e\n", hyperviscos);
    dump.add(hyper_mu_x, "hyper_mu_x", 1);
  }

  if(diffusion_par > 0.0) {
    output.write("    diffusion_par: %e\n", diffusion_par);
    dump.add(diffusion_par, "diffusion_par", 1);
  }

  //xqx: parallel hyper-viscous diffusion for pressure
  if(diffusion_p4 > 0.0) {
    output.write("    diffusion_p4: %e\n", diffusion_p4);
    dump.add(diffusion_p4, "diffusion_p4", 1);
  }

  //xqx: parallel hyper-viscous diffusion for vorticity
  if(diffusion_u4 > 0.0) {
    output.write("    diffusion_u4: %e\n", diffusion_u4);
    dump.add(diffusion_u4, "diffusion_u4", 1);
  }

  //xqx: parallel hyper-viscous diffusion for vector potential
  if(diffusion_a4 > 0.0) {
    output.write("    diffusion_a4: %e\n", diffusion_a4);
    dump.add(diffusion_a4, "diffusion_a4", 1);
  }

  if(heating_P > 0.0) {
    output.write("    heating_P(watts): %e\n", heating_P);
    dump.add(heating_P, "heating_P", 1);

    output.write("    hp_width(%): %e\n",hp_width);
    dump.add(hp_width, "hp_width", 1);

    output.write("    hp_length(%): %e\n",hp_length);
    dump.add(hp_length, "hp_length", 1);
  }

  if(sink_P > 0.0) {
    output.write("    sink_P(rate): %e\n", sink_P);
    dump.add(sink_P, "sink_P", 1);

    output.write("    sp_width(%): %e\n",sp_width);
    dump.add(sp_width, "sp_width", 1);

    output.write("    sp_length(%): %e\n",sp_length);
    dump.add(sp_length, "sp_length", 1);
  }

  if(K_H_term)                      
    output.write("    keep K-H term\n");
  else
    output.write("   drop K-H term\n"); 

  Field2D Te;
  Te = P0 / (2.0*density * 1.602e-19); // Temperature in eV

  J0 = - MU0*Lbar * J0 / B0;           
  P0 = 2.0*MU0 * P0 / (Bbar*Bbar);
  V0=V0/Va;        
  Dphi0*=Tbar;     
 
  b0xcv.x /= Bbar;
  b0xcv.y *= Lbar*Lbar;
  b0xcv.z *= Lbar*Lbar;

  Rxy  /= Lbar;
  Bpxy /= Bbar;
  Btxy /= Bbar;
  B0   /= Bbar;
  hthe /= Lbar;
  mesh->dx   /= Lbar*Lbar*Bbar;
  I *= Lbar*Lbar*Bbar;
  
  BoutReal pnorm = max(P0, true); // Maximum over all processors
  
  vacuum_pressure *= pnorm; // Get pressure from fraction
  vacuum_trans *= pnorm;

  // Transitions from 0 in core to 1 in vacuum
  vac_mask = (1.0 - tanh( (P0 - vacuum_pressure) / vacuum_trans )) / 2.0;

  if(spitzer_resist) {
    // Use Spitzer resistivity 
    output.write("\tTemperature: %e -> %e [eV]\n", min(Te), max(Te));
    eta = 0.51*1.03e-4*Zeff*20.*(Te^(-1.5)); // eta in Ohm-m. NOTE: ln(Lambda) = 20
    output.write("\tSpitzer resistivity: %e -> %e [Ohm m]\n", min(eta), max(eta));
    eta /= MU0 * Va * Lbar;
    output.write("\t -> Lundquist %e -> %e\n", 1.0/max(eta), 1.0/min(eta));
  }else {
    // transition from 0 for large P0 to resistivity for small P0
    eta = core_resist + (vac_resist - core_resist) * vac_mask;
  }

  dump.add(eta, "eta", 0);






	
  /**************** CALCULATE METRICS ******************/

    mesh->g11 = (Rxy*Bpxy)^2;
    mesh->g22 = 1.0 / (hthe^2);
    mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -I*mesh->g11;
    mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);

  mesh->J = hthe / Bpxy;
  mesh->Bxy = B0;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (B0*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry(); // Calculate quantities from metric tensor

  // Set B field vector
  
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;
  
  V0net.covariant = false;                              //presentation for net flow
  V0net.x = 0.;
  V0net.y = Rxy*Btxy*Bpxy/(hthe*B0*B0)*Dphi0;
  V0net.z = -Dphi0;

  U0=B0vec*Curl(V0net)/B0;     //get 0th vorticity for Kelvin-Holmholtz term





  if(simple_rmp)
    include_rmp = true;

  if(include_rmp) {
    // Including external field coils. 
    if(simple_rmp) {
      // Use a fairly simple form for the perturbation

      Field2D pol_angle;
      if(mesh->get(pol_angle, "dy")) {
	output.write("     ***WARNING: need poloidal angle for simple RMP\n");
	include_rmp = false;
      }else {
	OPTION(options, rmp_n,  3);
	OPTION(options, rmp_m,  9);
	OPTION(options, rmp_polwid, -1.0);
	OPTION(options, rmp_polpeak, 0.5);
	// Divide n by the size of the domain
        int zperiod;
        globalOptions->get("zperiod", zperiod, 1);
	if((rmp_n % zperiod) != 0)
	  output.write("     ***WARNING: rmp_n (%d) not a multiple of zperiod (%d)\n", rmp_n, zperiod);

	output.write("\tMagnetic perturbation: n = %d, m = %d, magnitude %e Tm\n", 
		     rmp_n, rmp_m, rmp_factor);
	
	rmp_Psi0 = 0.0;

	BoutReal ***d = rmp_Psi0.getData();

	// Set the outer boundary
	for(int jx=mesh->ngx-4;jx<mesh->ngx;jx++)
	  for(int jy=0;jy<mesh->ngy;jy++)
	    for(int jz=0;jz<mesh->ngz;jz++) {
	      
	      BoutReal angle = rmp_m * pol_angle[jx][jy] * ((BoutReal) mesh->YGLOBAL (jy)) + rmp_n * ((BoutReal) jz) * mesh->dz;

	      d[jx][jy][jz] = (((BoutReal)(jx - 4)) / ((BoutReal)(mesh->ngx - 5))) * rmp_factor * cos(angle);
		//	      d[jx][jy][jz] = (((BoutReal)(jx - 4)) / ((BoutReal)(mesh->ngx - 5))) * rmp_factor;

              if(rmp_polwid > 0.0) {
                // Multiply by a Gaussian in poloidal angle
output.write("\tCheck point, if rmp_polwid: %e\n", rmp_polwid);
                BoutReal gx = ((pol_angle[jx][jy] / (2.*PI)) - rmp_polpeak) / rmp_polwid;
                d[jx][jy][jz] *= exp(-gx*gx);
              }
	    }
	output.write("\tCheck point 101: %e ,%e\n", max(rmp_Psi0), min(rmp_Psi0) );
	
        rmp_Psib = rmp_Psi0;
	dump.add(rmp_Psib, "rmp_Psib", 0);
	rmp_Psi0 = rmp_Psi0.shiftZ(false); // Shift into field-aligned coords

	output.write("\tCheck point 1: %e ,%e\n", max(rmp_Psi0), min(rmp_Psi0) );

	// Now have a simple model for Psi due to coils at the outer boundary
	// Need to calculate Psi inside the domain, enforcing j = 0

        Jpar = 0.0;
output.write("\tCheck point 2\n");
        rmp_Psix = rmp_Psi0;
	dump.add(rmp_Psix, "rmp_Psix", 0);
	output.write("\tCheck point 3 %e\n", max(Jpar));

	  // rmp_Psi0 = invert_laplace_gui2(Jpar, rmp_Psi0,INVERT_AC_IN_LAP | INVERT_4TH_ORDER | INVERT_OUT_SET | INVERT_AC_OUT_LAP, NULL);
	  rmp_Psi0 = invert_laplace_gui2(Jpar, rmp_Psi0,INVERT_AC_IN_LAP | INVERT_4TH_ORDER | INVERT_OUT_SET, NULL);

	// Currently only implemented for 4th-order serial inversion (NXPE = 1)
      }
}else {
      // Load perturbation from grid file.
      include_rmp = mesh->get(rmp_Psi0, "rmp_A"); // Only include if found
      
      // Multiply by factor
      rmp_Psi0 *= rmp_factor;
    }
  }


  if(include_rmp) {
    // Normalise RMP quantities

    rmp_Psi0 /= Bbar * Lbar;
    
    rmp_ramp /= Tbar;
    rmp_freq *= Tbar;
    rmp_rotate *= Tbar;
    
    rmp_Psi = rmp_Psi0;
    rmp_dApdt = 0.0;

    bool apar_changing = false;

    output.write("Including magnetic perturbation\n");
    if(rmp_ramp > 0.0) {
      output.write("\tRamping up over period t = %e (%e ms)\n", rmp_ramp, rmp_ramp*Tbar*1000.);
      apar_changing = true;
    }
    if(rmp_freq > 0.0) {
      output.write("\tOscillating with frequency f = %e (%e kHz)\n", rmp_freq, rmp_freq/Tbar/1000.);
      apar_changing = true;
    }
    if(rmp_rotate != 0.0) {
      output.write("\tRotating with a frequency f = %e (%e kHz)\n", rmp_rotate, rmp_rotate/Tbar/1000.);
      apar_changing = true;
    }

    if(apar_changing) {
      dump.add(rmp_Psi, "rmp_Psi", 1);
      dump.add(rmp_dApdt, "rmp_dApdt", 1);
    }else {
      dump.add(rmp_Psi, "rmp_Psi", 0);
    }
  }else
    rmp_Psi = 0.0;



    rmp_B = Grad(rmp_Psi0) ^ B0vec;
    // Want contravariant components: B^x is change in psi, B^z is change in toroidal angle
    rmp_B.toContravariant();
  dump.add(rmp_B.x, "rmpbx", 0);
  dump.add(rmp_B.y, "rmpby", 0);
  dump.add(rmp_B.z, "rmpbz", 0);






  /**************** SET VARIABLE LOCATIONS *************/

  P.setLocation(CELL_CENTRE);
  U.setLocation(CELL_CENTRE);
  phi.setLocation(CELL_CENTRE);
  Psi.setLocation(CELL_YLOW);
  Jpar.setLocation(CELL_YLOW);
  Vpar.setLocation(CELL_YLOW);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(U);
  SOLVE_FOR(P);

  if(evolve_jpar) {
    output.write("Solving for jpar: Inverting to get Psi\n");
    SOLVE_FOR(Jpar);
    dump.add(Psi, "Psi", 1);
  }else {
    output.write("Solving for Psi, Differentiating to get jpar\n");
    SOLVE_FOR(Psi);
    dump.add(Jpar, "jpar", 1);
  }
  
  if(parallel_lagrange) {
    // Evolving the distortion of the flux surfaces (Ideal-MHD only!)
    
    bout_solve(Xip_x, "Xip_x");
    bout_solve(Xip_z, "Xip_z");
    
    bout_solve(Xim_x, "Xim_x");
    bout_solve(Xim_z, "Xim_z");
  }
  
  if(parallel_project) {
    // Add Xi to the dump file
    dump.add(Xip_x, "Xip_x", 1);
    dump.add(Xip_z, "Xip_z", 1);
    
    dump.add(Xim_x, "Xim_x", 1);
    dump.add(Xim_z, "Xim_z", 1);
  }

  if(compress0) {
    output.write("Including compression (Vpar) effects\n");
    
    SOLVE_FOR(Vpar);

    beta = B0*B0 / ( 0.5 + (B0*B0 / (g*P0)));
    gradparB = Grad_par(B0) / B0;

    output.write("Beta in range %e -> %e\n",
      min(beta), max(beta));
  }

  if(phi_constraint) {
    // Implicit Phi solve using IDA
    
    if(!bout_constrain(phi, C_phi, "phi")) {
      output.write("ERROR: Cannot constrain. Run again with phi_constraint=false\n");
      bout_error("Aborting.\n");
    }
    
    // Set preconditioner
    solver->setPrecon(precon_phi);

  }else {
    // Phi solved in RHS (explicitly)
    dump.add(phi, "phi", 1);

    // Set preconditioner
    solver->setPrecon(precon);

    // Set Jacobian
    solver->setJacobian(jacobian);
  }

  // Diamagnetic phi0
  if(diamag_phi0) {
    // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
    phi0 = -0.5*dnorm*P0/B0;
    SAVE_ONCE(phi0);
  }
  
 
  // Add some equilibrium quantities and normalisations
  // everything needed to recover physical units
  SAVE_ONCE2(J0, P0);
  SAVE_ONCE4(density, Lbar, Bbar, Tbar);
  SAVE_ONCE2(Va, B0);
  SAVE_ONCE2(Dphi0, U0); 
  SAVE_ONCE(V0);       
  
  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero
  
  if(!restarting) {
    // Only if not restarting: Check initial perturbation

    // Set U to zero where P0 < vacuum_pressure
    U = where(P0 - vacuum_pressure, U, 0.0);

    // Phi should be consistent with U
    phi = invert_laplace(U, phi_flags, NULL);
    
    if(diamag) {
      phi -= 0.5*dnorm * P / B0;
    }
  }

  /************** SETUP COMMUNICATIONS **************/
  
  comms.add(U);
  comms.add(Psi);
  comms.add(P);
  comms.add(phi);

  phi.setBoundary("phi"); // Set boundary conditions
  tmpU2.setBoundary("U");
  tmpP2.setBoundary("P");
  tmpA2.setBoundary("J");
  
  if(evolve_jpar) {
    comms.add(Jpar);
  }else {
    // otherwise Need to communicate Jpar separately
    Jpar.setBoundary("J");
  }
  Jpar2.setBoundary("J");

  return 0;
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  
  if(parallel_lagrange || parallel_project) {
    // Moving stencil locations
    
    Field3D fp, fm; // Interpolated on + and - y locations
    
    fp = interpolate(f, Xip_x, Xip_z);
    fm = interpolate(f, Xim_x, Xim_z);
    
    result.allocate();
    for(int i=0;i<mesh->ngx;i++)
      for(int j=1;j<mesh->ngy-1;j++)
	for(int k=0;k<mesh->ngz-1;k++) {
	  result[i][j][k] = (fp[i][j+1][k] - fm[i][j-1][k])/(2.*mesh->dy[i][j]*sqrt(mesh->g_22[i][j]));
	}
  }else {
    if(parallel_lr_diff) {
      // Use left/right biased stencils. NOTE: First order only!
      if(loc == CELL_YLOW) {
	result = Grad_par_CtoL(f);
      }else
	result = Grad_par_LtoC(f);
    }else
      result = Grad_par(f, loc);
    
    if(nonlinear) {
      result -= bracket(Psi + rmp_Psi, f, bm_mag)*B0;
    }
  }

  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t)
{
  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  if(nonlinear) {
    vac_mask = (1.0 - tanh( ((P0 + P) - vacuum_pressure) / vacuum_trans )) / 2.0;
    
    // Update resistivity
    if(spitzer_resist) {
      // Use Spitzer formula
      Field3D Te;
      Te = (P0+P)*Bbar*Bbar/(4.*MU0) / (density * 1.602e-19); // eV
      eta = 0.51*1.03e-4*Zeff*20.*(Te^(-1.5)); // eta in Ohm-m. ln(Lambda) = 20
      eta /= MU0 * Va * Lbar; // Normalised eta
    }else {
      // Use specified core and vacuum Lundquist numbers
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
    }
  }

  ////////////////////////////////////////////
  // Resonant Magnetic Perturbation code

  if(include_rmp) {

    if( (rmp_ramp > 0.0) || (rmp_freq > 0.0) || (rmp_rotate != 0.0) ) {
      // Need to update the RMP terms
    
      if((rmp_ramp > 0.0) && (t < rmp_ramp)) {
	// Still in ramp phase
	
	rmp_Psi = (t / rmp_ramp) * rmp_Psi0 ; // Linear ramp
	
	rmp_dApdt = rmp_Psi0 / rmp_ramp;
      }else {
	rmp_Psi = rmp_Psi0;
	rmp_dApdt = 0.0;
      }
      
      if(rmp_freq > 0.0) {
	// Oscillating the amplitude
	
	rmp_dApdt = rmp_dApdt * sin(2.*PI * rmp_freq * t)
	  + rmp_Psi * (2.*PI * rmp_freq) * cos(2.*PI * rmp_freq * t);
	
	rmp_Psi *= sin(2.*PI * rmp_freq * t);
	
      }
      
      if(rmp_rotate != 0.0) {
	// Rotate toroidally at given frequency
	
	rmp_Psi = rmp_Psi.shiftZ(2*PI * rmp_rotate * t);
	
	rmp_dApdt = rmp_dApdt.shiftZ(2*PI * rmp_rotate * t);
	
	// Add toroidal rotation term. CHECK SIGN
	
	rmp_dApdt += DDZ(rmp_Psi) * 2*PI * rmp_rotate;
      }
      
      // Set to zero in the core
	rmp_Psi *= 1.0 ;//vac_mask;
    }else {
      // Set to zero in the core region
      rmp_Psi = rmp_Psi0; // * vac_mask;  // Only in vacuum -> skin current -> diffuses inwards
    }
  }
  
  ////////////////////////////////////////////
  // Inversion

  if(evolve_jpar) {
    // Invert laplacian for Psi
    Psi = invert_laplace(Jpar, apar_flags, NULL);
  }

  if(phi_constraint) {
    // Phi being solved as a constraint
    
    Field3D Ctmp = phi;
    Ctmp.setBoundary("phi"); // Look up boundary conditions for phi
    Ctmp.applyBoundary();
    Ctmp -= phi; // Now contains error in the boundary
    
    C_phi = Delp2(phi) - U; // Error in the bulk
    C_phi.setBoundaryTo(Ctmp);

  }else {
    // Invert laplacian for phi
    phi = invert_laplace(U, phi_flags, NULL);
    
    if(diamag) {
      phi -= 0.5*dnorm * P / B0;
    }
    
    // Apply a boundary condition on phi for target plates
    phi.applyBoundary();
  }

  // Perform communications
  mesh->communicate(comms);

  if(!evolve_jpar) {
    // Get J from Psi
    Jpar = Delp2(Psi + rmp_Psi);

    Jpar.applyBoundary();
    mesh->communicate(Jpar);

    if(jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary
      
      for(int i=0;i<jpar_bndry_width;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    if(mesh->firstX())
	      Jpar[i][j][k] = 0.0;
	    if(mesh->lastX())
	      Jpar[mesh->ngx-1-i][j][k] = 0.0;
	  }
    }

    // Smooth j in x
    if(smooth_j_x)
      Jpar = smooth_x(Jpar);

    //xqx begin
    // Get Delp2(J) from J
    Jpar2 = Delp2(Jpar);

    Jpar2.applyBoundary();
    mesh->communicate(Jpar2);

    if(jpar_bndry_width > 0) {
      // Zero jpar2 in boundary regions. Prevents vorticity drive
      // at the boundary
      
      for(int i=0;i<jpar_bndry_width;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    if(mesh->firstX())
	      Jpar2[i][j][k] = 0.0;
	    if(mesh->lastX())
	      Jpar2[mesh->ngx-1-i][j][k] = 0.0;
	  }
    }
    //xqx end
  }
  ////////////////////////////////////////////////////
  // Project perturbed field-lines

  if(parallel_project) {
    Vector3D Btilde; // Perturbed field
    Btilde = Grad(Psi) ^ B0vec;
    // Want contravariant components: B^x is change in psi, B^z is change in toroidal angle
    Btilde.toContravariant();
    
    // Calculate displacements in index space
    Xip_x = 0.;
    Xip_z = 0.;
    Xim_x = 0.;
    Xim_z = 0.;
    Field2D sgx = sqrt(mesh->g_11); // 1/(R*Bp)
    Field2D sgz = sqrt(mesh->g_33); // 1/R
    for(int i=0;i<mesh->ngx;i++)
      for(int j=1;j<mesh->ngy-1;j++)
	for(int k=0;k<mesh->ngz-1;k++) {
	  BoutReal bxratio = (Btilde.x[i][j][k]*sgx[i][j])/B0[i][j]; // |Bx| / B
	  bxratio *= (hthe[i][j] * B0[i][j] / Bpxy[i][j]); // Multiply by parallel length
	  Xip_x[i][j+1][k] = bxratio / (sgx[i][j+1] * mesh->dx[i][j+1]); // Convert to psi, then index
	  Xim_x[i][j-1][k] = bxratio / (sgx[i][j-1] * mesh->dx[i][j-1]); 
	  
	  BoutReal bzratio = (Btilde.z[i][j][k]*sgz[i][j])/B0[i][j]; // |Bz| / B
	  bzratio *= (hthe[i][j] * B0[i][j] / Bpxy[i][j]); // Multiply by parallel length
	  Xip_z[i][j+1][k] = bzratio / (sgz[i][j+1] * mesh->dz); // Convert to psi, then index
	  Xim_z[i][j-1][k] = bzratio / (sgz[i][j-1] * mesh->dz); 
	}
  }

  ////////////////////////////////////////////////////
  // Parallel electric field

  if(evolve_jpar) {
    // Jpar
    
    ddt(Jpar) = -Grad_parP(B0*U, CELL_YLOW) / B0 + eta*Delp2(Jpar);

    /*
      ddt(Psi) = -Grad_parP(B0*phi, CELL_YLOW) / B0 + eta*Jpar;
      // Apply boundary condition on Psi
      apply_boundary(ddt(Psi), "Psi");
      mesh->communicate(ddt(psi));
      ddt(Jpar) = Delp2(ddt(Psi));
    */

    if(relax_j_vac) {
      // Make ddt(Jpar) relax to zero.
      
      ddt(Jpar) -= vac_mask * Jpar / relax_j_tconst;
    }
  }else {
    // Vector potential
      
    ddt(Psi) = -Grad_parP(B0*phi, CELL_CENTRE) / B0 + eta*Jpar;
    //xqx      ddt(Psi) = -Grad_parP(B0*phi, CELL_YLOW) / B0 + eta*Jpar;

    if(eHall) {
      ddt(Psi) +=  0.25*delta_i*(Grad_parP(B0*P, CELL_CENTRE) / B0 
                                 +b0xGrad_dot_Grad(P0, Psi));   // electron parallel pressure
    }

    if(diamag_phi0)
      ddt(Psi) -= b0xGrad_dot_Grad(phi0, Psi);   // Equilibrium flow

    if(withflow)                                //net flow 
      ddt(Psi)-= V_dot_Grad(V0net, Psi);

    //F_Psi = -Grad_par_CtoL(B0*phi) / B0 + eta*Jpar;
    
    if(diamag_grad_t) {
      // grad_par(T_e) correction

      ddt(Psi) += 1.71 * dnorm * 0.5 * Grad_parP(P, CELL_YLOW) / B0;
    }

    // Hyper-resistivity
    if(hyperresist > 0.0) {
      ddt(Psi) -= eta*hyperresist * Delp2(Jpar);
    }

    // electron Hyper-viscosity coefficient
    if(ehyperviscos > 0.0) {
      ddt(Psi) -= eta*ehyperviscos * Delp2(Jpar2);
    }

    //xqx: parallel hyper-viscous diffusion for vector potential
    if(diffusion_a4 > 0.0){
      tmpA2 = Grad2_par2new(Psi);
      mesh->communicate(tmpA2);
      tmpA2.applyBoundary();
      ddt(Psi) -= diffusion_a4 * Grad2_par2new(tmpA2);}

    // Vacuum solution
    if(relax_j_vac) {
      // Calculate the J and Psi profile we're aiming for
      Field3D Jtarget = Jpar * (1.0 - vac_mask); // Zero in vacuum
      
      // Invert laplacian for Psi
      Psitarget = invert_laplace(Jtarget, apar_flags, NULL);
      
      // Add a relaxation term in the vacuum
      ddt(Psi) = ddt(Psi)*(1. - vac_mask) - (Psi - Psitarget)*vac_mask / relax_j_tconst;
    }
  }

  if(parallel_lagrange) {
    // Move tracers with the plasma
    //advect_tracer(phi, Xi_x, Xi_z, F_Xi_x, F_Xi_z);

    // Follow intersection of fieldlines with neighbouring
    // poloidal planes

    // Calculate the ExB velocity on the grid-points
    Field3D vx0, vz0, vxp, vzp, vxm, vzm;
    vx0 = DDZ(phi);
    vz0 = -DDX(phi);
    
    // Interpolate to get velocity at intersections
    vxp = interpolate(vx0, Xip_x, Xip_z);
    vzp = interpolate(vz0, Xip_x, Xip_z);
    
    vxm = interpolate(vx0, Xim_x, Xim_z);
    vzm = interpolate(vz0, Xim_x, Xim_z);
    
    // Convert into relative velocity of
    // grid points and intersections.
    // NOTE: Should be Christoffel terms here
  }

  ////////////////////////////////////////////////////
  // Vorticity equation

  ddt(U) = (B0^2) * b0xGrad_dot_Grad(Psi + rmp_Psi, J0, CELL_CENTRE); // Grad j term

  ddt(U) += b0xcv*Grad(P);  // curvature term

  if(!nogradparj) {
    // Parallel current term
    ddt(U) -= (B0^2)*Grad_parP(Jpar, CELL_CENTRE); // b dot grad j
    //ddt(U) -= (B0^2)*Grad_par_LtoC(Jpar);
  }

  if(withflow&&K_H_term)                         //K_H_term
    ddt(U) -=b0xGrad_dot_Grad(phi,U0); 

  if(diamag_phi0)
    ddt(U) -= b0xGrad_dot_Grad(phi0, U);   // Equilibrium flow

  if(withflow)                            // net flow
    ddt(U) -= V_dot_Grad(V0net, U);

  if(nonlinear) {
    ddt(U) -= bracket(phi, U, bm_exb)*B0;    // Advection
  }

  // Viscosity terms 
  if(viscos_par > 0.0)
    ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
  
  //xqx: parallel hyper-viscous diffusion for vorticity
  if(diffusion_u4 > 0.0){
    tmpU2 = Grad2_par2new(U);
    mesh->communicate(tmpU2);
    tmpU2.applyBoundary();
    //    tmpU2.applyBoundary("neumann");
    ddt(U) -= diffusion_u4 * Grad2_par2new(tmpU2);}
  
  if(viscos_perp > 0.0)
    ddt(U) += viscos_perp * Delp2(U);     // Perpendicular viscosity

  // Hyper-viscosity
  if(hyperviscos > 0.0) {
    // Calculate coefficient.
    
    hyper_mu_x = hyperviscos * mesh->g_11*SQ(mesh->dx) * abs(mesh->g11*D2DX2(U)) / (abs(U) + 1e-3);
    hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries
    
    ddt(U) += hyper_mu_x * mesh->g11*D2DX2(U);
    
    if(first_run) { // Print out maximum values of viscosity used on this processor
      output.write("   Hyper-viscosity values:\n");
      output.write("      Max mu_x = %e, Max_DC mu_x = %e\n", max(hyper_mu_x), max(hyper_mu_x.DC()));
    }
  }

  // left edge sink terms 
  if(sink_Ul > 0.0){
    ddt(U) -=  sink_Ul*sink_tanhxl(P0,U,su_widthl,su_lengthl); // core sink
  }

  // right edge sink terms 
  if(sink_Ur > 0.0){
    ddt(U) -=  sink_Ur*sink_tanhxr(P0,U,su_widthr,su_lengthr); //  sol sink
  }

  ////////////////////////////////////////////////////
  // Pressure equation

  ddt(P) = 0.0;
  if(evolve_pressure) {
    ddt(P) -= b0xGrad_dot_Grad(phi, P0);
    
    if(diamag_phi0)
      ddt(P) -= b0xGrad_dot_Grad(phi0, P);   // Equilibrium flow

    if(withflow)                              //net flow
      ddt(P) -= V_dot_Grad(V0net, P);  

    if(nonlinear)
      ddt(P) -= bracket(phi, P, bm_exb)*B0;    // Advection      
  }

  // Parallel diffusion terms 
  if(diffusion_par > 0.0)
    ddt(P) += diffusion_par * Grad2_par2(P); // Parallel diffusion

  //xqx: parallel hyper-viscous diffusion for pressure
  if(diffusion_p4 > 0.0){
    tmpP2 = Grad2_par2new(P);
    mesh->communicate(tmpP2);
    tmpP2.applyBoundary();
    ddt(P) = diffusion_p4 * Grad2_par2new(tmpP2);}

  // heating source terms 
  if(heating_P > 0.0){
    BoutReal pnorm = P0[0][0];
    ddt(P) += heating_P*source_expx2(P0,2.*hp_width,0.5*hp_length)*(Tbar/pnorm); // heat source
    ddt(P) += (100.*source_tanhx(P0,hp_width,hp_length)+0.01) * mesh->g11 * D2DX2(P) * (Tbar/Lbar/Lbar) ;     // radial diffusion
  }

  // sink terms 
  if(sink_P > 0.0){
    ddt(P) -= sink_P*sink_tanhxr(P0,P,sp_width,sp_length)*Tbar; // sink
  }

  ////////////////////////////////////////////////////
  // Compressional effects
  
  if(compress0) {
    
    //ddt(P) += beta*( - Grad_parP(Vpar, CELL_CENTRE) + Vpar*gradparB );
    ddt(P) -= beta*Div_par_CtoL(Vpar);
    
    if(phi_curv) {
      ddt(P) -= 2.*beta*b0xcv*Grad(phi);
    }

    // Vpar equation

    //ddt(Vpar) = -0.5*Grad_parP(P + P0, CELL_YLOW);
    ddt(Vpar) = -0.5*Grad_par_LtoC(P + P0);

    if(nonlinear)
      ddt(Vpar) -= bracket(phi, Vpar, bm_exb)*B0; // Advection
  }
  
  if(filter_z) {
    // Filter out all except filter_z_mode
    
    if(evolve_jpar) {
      ddt(Jpar) = filter(ddt(Jpar), filter_z_mode);
    }else
      ddt(Psi) = filter(ddt(Psi), filter_z_mode);

    ddt(U) = filter(ddt(U), filter_z_mode);
    ddt(P) = filter(ddt(P), filter_z_mode);
  }

  if(low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
    if(evolve_jpar) {
      ddt(Jpar) = lowPass(ddt(Jpar), low_pass_z, zonal_field);
    }else
      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

    ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
    ddt(P) = lowPass(ddt(P), low_pass_z, zonal_bkgd);
  }

  if(damp_width > 0) {
    for(int i=0;i<damp_width;i++) {
      for(int j=0;j<mesh->ngy;j++)
	for(int k=0;k<mesh->ngz;k++) {
	  if(mesh->firstX())
	    ddt(U)[i][j][k] -= U[i][j][k] / damp_t_const;
	  if(mesh->lastX())
	    ddt(U)[mesh->ngx-1-i][j][k] -= U[mesh->ngx-1-i][j][k] / damp_t_const;
	}
    }
  }

  first_run = false;

  return 0;
}


/*******************************************************************************
 * Preconditioner described in elm_reduced.cpp
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 * 
 * o Return values should be in vars (overwriting system state)
 *
 * NOTE: EXPERIMENTAL
 * enable by setting solver / use_precon = true in BOUT.inp
 *******************************************************************************/

const Field3D Lp(const Field3D &f)
{
  return b0xcv*Grad(f);
}

const Field3D Lpsi(const Field3D &f)
{
  Field3D jpre = Delp2(f);
  jpre.setBoundary("J");
  jpre.applyBoundary();

  mesh->communicate(jpre);

  Field3D result = b0xGrad_dot_Grad(f, J0);
  
  return (B0^2)*(result - Grad_parP(jpre) );
}

const Field3D Pschur(const Field3D &f)
{
  Field3D phitmp = invert_laplace(f, phitmp, phi_flags, NULL);
  // Need to communicate phi
  mesh->communicate(phitmp);

  Field3D dP = b0xGrad_dot_Grad(phitmp, P0);
  Field3D dPsi = Grad_parP(B0*phitmp) / B0;
  dP.setBoundary("P"); dP.applyBoundary();
  dPsi.setBoundary("Psi"); dPsi.applyBoundary();
  
  Field3D dJ = Delp2(dPsi);
  dJ.setBoundary("J"); dJ.applyBoundary();

  mesh->communicate(dP, dPsi, dJ);
  
  Field3D result = b0xcv*Grad(dP) + (B0^2)*(Grad_par(dJ) - b0xGrad_dot_Grad(dPsi, J0));
  result.setBoundary("U"); result.applyBoundary();
  
  return result;
}

int precon(BoutReal t, BoutReal gamma, BoutReal delta)
{
  Field3D U1;
  Field3D P2, Psi2, U2;
  Field3D P3, Psi3, U3;

  //output.write("precon t = %e, gamma = %e\n", t, gamma);

  mesh->communicate(ddt(P), ddt(Psi), ddt(U));

  // First matrix. Only modifies vorticity
  
  U1 = ddt(U) + gamma*(Lp(ddt(P)) + Lpsi(ddt(Psi)));
  U1.setBoundary("U"); U1.applyBoundary();
  
  // Second matrix. If linear, only modify vorticity
  // NB: This is the key step, inverting Pschur
  
  P2 = ddt(P);
  Psi2 = ddt(Psi);

  U2=U1;
  //  U2 = invert_parderiv(1.0, -gamma*gamma*B0*B0, U1);
  U2 = U1 + gamma*gamma*Pschur(U1); // Binomial expansion of P^-1
  
  // Third matrix
  
  Field3D phitmp = invert_laplace(U2, phi_flags, NULL);
  mesh->communicate(phitmp);
  
  P3 = P2 + gamma*b0xGrad_dot_Grad(phitmp, P0);
  Psi3 = Psi2 + gamma*Grad_par(B0*phitmp) / B0;
  U3 = U2;
  
  // Put result into system state
  
  P = P3;
  Psi = Psi3;
  U = U3;
  
  P.applyBoundary();
  Psi.applyBoundary();
  U.applyBoundary();
  
  return 0;
}

/*******************************************************************************
 * Jacobian-vector multiply
 *
 * Input
 *   System state is in (P, Psi, U)
 *   Vector v is in (F_P, F_Psi, F_U)
 * Output
 *   Jacobian-vector multiplied Jv should be in (P, Psi, U)
 * 
 * NOTE: EXPERIMENTAL
 * enable by setting solver / use_jacobian = true in BOUT.inp
 *******************************************************************************/

int jacobian(BoutReal t)
{
  // NOTE: LINEAR ONLY!
  
  // Communicate
  mesh->communicate(ddt(P), ddt(Psi), ddt(U));

  phi = invert_laplace(ddt(U), phi_flags, NULL);

  Jpar = Delp2(ddt(Psi));

  mesh->communicate(phi, Jpar);

  Field3D JP = -b0xGrad_dot_Grad(phi, P0);
  JP.setBoundary("P"); JP.applyBoundary();
  
  Field3D JPsi = -Grad_par(B0*phi, CELL_YLOW) / B0;
  JPsi.setBoundary("Psi"); JPsi.applyBoundary();

  Field3D JU = b0xcv*Grad(ddt(P))
    - (B0^2)*Grad_par(Jpar, CELL_CENTRE)
    + (B0^2) * b0xGrad_dot_Grad(ddt(Psi), J0, CELL_CENTRE);
  JU.setBoundary("U"); JU.applyBoundary();

  // Put result into vars

  P = JP;
  Psi = JPsi;
  U = JU;

  return 0;
}

/*******************************************************************************
 * Preconditioner for when phi solved as a constraint
 * Currently only possible with the IDA solver
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 * 
 * o Return values should be in vars (overwriting system state)
 *******************************************************************************/

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta)
{
  P = ddt(P);
  Psi = ddt(Psi);
  /*
  invert_laplace(F_U, phi, phi_flags, NULL);
  phi = C_phi + phi;
  */
  phi = invert_laplace(C_phi - ddt(U), phi_flags, NULL);
  U = ddt(U);
  return 0;
}

/*******************************************************************************
 * Tracer advection code
 * 
 * May be used to track field-lines and calculate Grad_par along
 * highly perturbed field-lines.
 *
 * Aug 2009
 *******************************************************************************/

void advect_tracer(const Field3D &p,  // phi (input)
		   const Field3D &delta_x, const Field3D &delta_z, // Current location (input) 
		   Field3D &F_dx, Field3D &F_dz)
{
  // Calculate the ExB velocity
  Field3D vx, vz;
  
  vx = DDZ(p);
  vz = -DDX(p);

  // Interpolate these velocities onto the current locations
  vx = interpolate(vx, delta_x, delta_z);
  vz = interpolate(vz, delta_x, delta_z);
  
  // Time derivative of the cell index:
  F_dx = vx / mesh->dx;
  F_dz = vz / mesh->dz;
}






//**********************************************************************************gui2






/// Laplacian inversion initialisation. Called once at the start to get settings
int invert_init()
{
  BoutReal filter; ///< Fraction of Z modes to filter out. Between 0 and 1

  output.write("Initialising Laplacian inversion routines\n");

  Options *options = Options::getRoot();

  // Communication options
  Options *commOpts = options->getSection("comms");
  commOpts->get("async", invert_async_send, true);

  // Inversion options
  Options *lapOpts = options->getSection("laplace");
  OPTION(lapOpts, filter, 0.2);
  lapOpts->get("low_mem", invert_low_mem, false);
  lapOpts->get("use_pdd", invert_use_pdd, false);
  lapOpts->get("all_terms", laplace_all_terms, false);
  OPTION(lapOpts, laplace_nonuniform, true);

  if(mesh->firstX() && mesh->lastX()) {
    // This processor is both the first and the last in X
    // -> No parallelisation needed
    output.write("\tUsing serial algorithm\n");

  }else {
    // Need to use a parallel algorithm
    if(invert_use_pdd) {
      output.write("\tUsing PDD algorithm\n");
    }else
      output.write("\tUsing parallel Thomas algorithm\n");
  }

  int ncz = mesh->ngz-1;

  // convert filtering into an integer number of modes
  laplace_maxmode = ROUND((1.0 - filter) * ((double) (ncz / 2)));

  // Can be overriden by max_mode option
  lapOpts->get("max_mode", laplace_maxmode, laplace_maxmode);

  if(laplace_maxmode < 0) laplace_maxmode = 0;
  if(laplace_maxmode > ncz/2) laplace_maxmode = ncz/2;

  // Broadcast this since rounding errors could cause mismatch across processors
  // THIS LINE CAUSES SEGFAULT ON LLNL GRENDEL
  //MPI_Bcast(&laplace_maxmode, 1, MPI_INT, 0, BoutComm::get());

  return 0;
}

/// Returns the coefficients for a tridiagonal matrix for laplace. Used by Delp2 too

/**********************************************************************************
 *                                 SERIAL CODE
 **********************************************************************************/

/// Perpendicular laplacian inversion (serial)
/*!
 * Inverts an X-Z slice (FieldPerp) using band-diagonal solvers
 * This code is only for serial i.e. mesh->NXPE == 1
 */
int invert_laplace_ser(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a,
                       const Field2D *ccoef, const Field2D *d) {
  int ncx = mesh->ngx-1;
  int ncz = mesh->ngz-1;

  static dcomplex **bk = NULL, *bk1d;
  static dcomplex **xk, *xk1d;
  int xbndry; // Width of the x boundary

  if(!mesh->firstX() || !mesh->lastX()) {
    output.write("Error: invert_laplace only works for mesh->NXPE = 1\n");
    return 1;
  }

  x.allocate();

  int jy = b.getIndex();
  x.setIndex(jy);

  if(bk == NULL) {
    // Allocate memory

    bk = cmatrix(mesh->ngx, ncz/2 + 1);
    bk1d = new dcomplex[mesh->ngx];

    xk = cmatrix(mesh->ngx, ncz/2 + 1);
    xk1d = new dcomplex[mesh->ngx];
  }

  xbndry = 2;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;

  #pragma omp parallel for
  for(int ix=0;ix<mesh->ngx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)

    ZFFT(b[ix], mesh->zShift[ix][jy], bk[ix]);
  }

  if(!mesh->periodicX) {
    if(flags & INVERT_IN_SET) {
      // Setting the inner boundary from x
      #pragma omp parallel for
      for(int ix=0;ix<xbndry;ix++)
	ZFFT(x[ix], mesh->zShift[ix][jy], xk[ix]);
    }

    if(flags & INVERT_OUT_SET) {
      // Setting the outer boundary from x
      #pragma omp parallel for
      for(int ix=0;ix<xbndry;ix++)
	ZFFT(x[ncx-ix], mesh->zShift[ncx-ix][jy], xk[ncx-ix]);
    }
  }

  if((flags & INVERT_4TH_ORDER) && (!mesh->periodicX)) { // Not implemented for parallel calculations or periodic X
    // Use band solver - 4th order

    static dcomplex **A = (dcomplex**) NULL;
    int xstart, xend;

    if(A == (dcomplex**) NULL)
      A = cmatrix(mesh->ngx, 5);

    // Get range for 4th order: Need at least 2 each side
    if(xbndry > 1) {
      xstart = xbndry;
      xend = ncx-xbndry;
    }else {
      xstart = 2;
      xend = mesh->ngx-2;
    }

    for(int iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x

      BoutReal coef1=0.0, coef2=0.0, coef3=0.0, coef4=0.0,
        coef5=0.0, coef6=0.0, kwave, flt;
      ///////// PERFORM INVERSION /////////

      // shift freqs according to FFT convention
      kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]

      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;

      // set bk1d
      for(int ix=0;ix<mesh->ngx;ix++)
	bk1d[ix] = bk[ix][iz]*flt;

      // Fill in interior points

      for(int ix=xstart;ix<=xend;ix++) {

	// Set coefficients
	coef1 = mesh->g11[ix][jy];  // X 2nd derivative
	coef2 = mesh->g33[ix][jy];  // Z 2nd derivative
	coef3 = mesh->g13[ix][jy];  // X-Z mixed derivatives
	coef4 = 0.0;          // X 1st derivative
	coef5 = 0.0;          // Z 1st derivative
	coef6 = 0.0;          // Constant

        if(d != (Field2D*) NULL) {
          // Multiply Delp2 component by a factor
          coef1 *= (*d)[ix][jy];
          coef2 *= (*d)[ix][jy];
          coef3 *= (*d)[ix][jy];
        }

	if(a != (Field2D*) NULL)
	  coef6 = (*a)[ix][jy];

	if(laplace_all_terms) {
	  coef4 = mesh->G1[ix][jy];
	  coef5 = mesh->G3[ix][jy];
	}
       	laplace_nonuniform = true;
	if(laplace_nonuniform) {
	  //    output.write("\t************************************check point 11: laplace_nonuniform\n");
	  // non-uniform mesh correction
	  if((ix != 0) && (ix != ncx))
	    //	    coef4 += mesh->g11[ix][jy]*( (1.0/mesh->dx[ix+1][jy]) - (1.0/mesh->dx[ix-1][jy]) )/(2.0*mesh->dx[ix][jy]);
      coef4 -= 0.5*((mesh->dx[ix+1][jy] - mesh->dx[ix-1][jy])/SQ(mesh->dx[ix][jy]))*coef1; // BOUT-06 term
	}

	if(ccoef != NULL) {
	  // A first order derivative term (1/c)\nabla_perp c\cdot\nabla_\perp x

	  if((ix > 1) && (ix < (mesh->ngx-2)))
	    coef4 += mesh->g11[ix][jy] * ((*ccoef)[ix-2][jy] - 8.*(*ccoef)[ix-1][jy] + 8.*(*ccoef)[ix+1][jy] - (*ccoef)[ix+2][jy]) / (12.*mesh->dx[ix][jy]*((*ccoef)[ix][jy]));
	}

	// Put into matrix
	coef1 /= 12.* SQ(mesh->dx[ix][jy]);
	coef2 *= SQ(kwave);
	coef3 *= kwave / (12. * mesh->dx[ix][jy]);
	coef4 /= 12. * mesh->dx[ix][jy];
	coef5 *= kwave;

	A[ix][0] = dcomplex(    -coef1 +   coef4 ,     coef3 );
	A[ix][1] = dcomplex( 16.*coef1 - 8*coef4 , -8.*coef3 );
	A[ix][2] = dcomplex(-30.*coef1 - coef2 + coef6, coef5);
	A[ix][3] = dcomplex( 16.*coef1 + 8*coef4 ,  8.*coef3 );
	A[ix][4] = dcomplex(    -coef1 -   coef4 ,    -coef3 );
      }

      if(xbndry < 2) {
	// Use 2nd order near edges

	int ix = 1;

	coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

        if(d != (Field2D*) NULL) {
          // Multiply Delp2 component by a factor
          coef1 *= (*d)[ix][jy];
          coef2 *= (*d)[ix][jy];
          coef3 *= (*d)[ix][jy];
        }

	A[ix][0] = 0.0; // Should never be used
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;

	ix = ncx-1;

	coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	A[ix][0] = 0.0;
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;  // Should never be used
      }

      // Boundary conditions

      for(int ix=0;ix<xbndry;ix++) {
	// Set zero-value. Change to zero-gradient if needed

        if(!(flags & INVERT_IN_RHS))
          bk1d[ix] = 0.0;
        if(!(flags & INVERT_OUT_RHS))
          bk1d[ncx-ix] = 0.0;

	A[ix][0] = A[ix][1] = A[ix][3] = A[ix][4] = 0.0;
	A[ix][2] = 1.0;

	A[ncx-ix][0] = A[ncx-ix][1] = A[ncx-ix][3] = A[ncx-ix][4] = 0.0;
	A[ncx-ix][2] = 1.0;
      }

      if(flags & INVERT_IN_SET) {
	// Set values of inner boundary from X
	for(int ix=0;ix<xbndry;ix++)
	  bk1d[ix] = xk[ix][iz];
      }

      if(flags & INVERT_OUT_SET) {
	// Set values of outer boundary from X
	for(int ix=0;ix<xbndry;ix++)
	  bk1d[ncx-ix] = xk[ncx-ix][iz];
	//	output.write("\t**************************************Used the 4th order invert\n");
      }

      if(iz == 0) {
	// DC

	// Inner boundary
	if(flags & INVERT_DC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}

	// Outer boundary
	if(flags & INVERT_DC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}

      }else {
	// AC

	// Inner boundarySQ(kwave)*coef2
	if(flags & INVERT_AC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}else if(flags & INVERT_AC_IN_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  //	output.write("\t*************************************Enforce zero laplacian for 2nd and 4th-orde\n");
	  int ix = 1;

	  coef1=mesh->g11[ix][jy]/(12.* SQ(mesh->dx[ix][jy]));

	  coef2=mesh->g33[ix][jy];

	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];

	  // Combine 4th order at 1 with 2nd order at 0
	  A[1][0] = 0.0; // Not used
	  A[1][1] = dcomplex( (14. - SQ(mesh->dx[0][jy]*kwave)*mesh->g33[0][jy]/mesh->g11[0][jy])*coef1  ,  -coef3 );
	  A[1][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[1][3] = dcomplex( 16.*coef1  , coef3 );
	  A[1][4] = dcomplex(    -coef1  ,     0.0 );

	  coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	  // Use 2nd order at 1
	  A[0][0] = 0.0;  // Should never be used
	  A[0][1] = 0.0;
	  A[0][2] = dcomplex(coef1, -coef3);
	  A[0][3] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	  A[0][4] = dcomplex(coef1,  coef3);

	}

	// Outer boundary
	if(flags & INVERT_AC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}else if(flags & INVERT_AC_OUT_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  // NOTE: Currently ignoring XZ term and coef4 assumed zero on boundary
	  // FIX THIS IF IT WORKS

	  int ix = ncx-1;

	  coef1=mesh->g11[ix][jy]/(12.* SQ(mesh->dx[ix][jy]));

	  coef2=mesh->g33[ix][jy];

	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];

	  // Combine 4th order at ncx-1 with 2nd order at ncx
	  A[ix][0] = dcomplex(    -coef1  ,     0.0 );
	  A[ix][1] = dcomplex( 16.*coef1  , -coef3 );
	  A[ix][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[ix][3] = dcomplex( (14. - SQ(mesh->dx[ncx][jy]*kwave)*mesh->g33[ncx][jy]/mesh->g11[ncx][jy])*coef1  ,  coef3 );
	  A[ix][4] = 0.0; // Not used

	  coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	  // Use 2nd order at ncx - 1
	  A[ncx][0] = dcomplex(coef1, -coef3);
	  A[ncx][1] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	  A[ncx][2] = dcomplex(coef1,  coef3);
	  A[ncx][3] = 0.0;  // Should never be used
	  A[ncx][4] = 0.0;
	}
      }

      // Perform inversion
      cband_solve(A, mesh->ngx, 2, 2, bk1d);

      if((flags & INVERT_KX_ZERO) && (iz == 0)) {
        // Set the Kx = 0, n = 0 component to zero. For now just subtract
        // Should do in the inversion e.g. Sherman-Morrison formula

        dcomplex offset(0.0);
        for(int ix=0;ix<=ncx;ix++)
          offset += bk1d[ix];
        offset /= (BoutReal) (ncx+1);
        for(int ix=0;ix<=ncx;ix++)
          bk1d[ix] -= offset;
      }

      // Fill xk
      for (int ix=0; ix<=ncx; ix++)
	xk[ix][iz]=bk1d[ix];

    }
  }else {
    // Use tridiagonal system in x - 2nd order

    static dcomplex *avec = (dcomplex*) NULL, *bvec, *cvec;

    if(avec == (dcomplex*) NULL) {
      avec = new dcomplex[mesh->ngx];
      bvec = new dcomplex[mesh->ngx];
      cvec = new dcomplex[mesh->ngx];
    }
    //output.write("\t**************************************Used the 2th order invert\n");
    for(int iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x

      // set bk1d
      BoutReal flt;
      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;

      for(int ix=0;ix<=ncx;ix++)
	bk1d[ix] = bk[ix][iz] * flt;

      ///////// PERFORM INVERSION /////////

      laplace_nonuniform = true;
      for(int ix=xbndry;ix<=ncx-xbndry;ix++) {
	laplace_tridag_coefs_gui2(ix, jy, iz, avec[ix], bvec[ix], cvec[ix], ccoef, d);

	if(a != (Field2D*) NULL)
	  bvec[ix] += (*a)[ix][jy];
      }

      if(!mesh->periodicX) {
	// Need boundary conditions
	/// By default, set RHS to zero, unless INVERT_*_RHS set
	if(!(flags & INVERT_IN_RHS)) {
	  for(int ix=0;ix<xbndry;ix++)
	    bk1d[ix] = 0.;
	}
	if(!(flags & INVERT_OUT_RHS)) {
	  for(int ix=mesh->ngx-xbndry;ix<mesh->ngx;ix++)
	    bk1d[ix] = 0.;
	}

	// Set boundary conditions

	if(iz == 0) {
	  // DC

	  // Inner boundary
	  if(flags & INVERT_DC_IN_GRAD) {
	    // Zero gradient at inner boundary

	    if((flags & INVERT_IN_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient

	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }
	      // Symmetric on last point
	      avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.0,0.0);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_IN_SET) {
	    for(int ix=0;ix<xbndry;ix++) {
	      avec[ix] = 0.0;
	      bvec[ix] = 1.0;
	      cvec[ix] = 0.0;
	      bk1d[ix] = xk[ix][iz];
	    }
	  }else {
	    // Zero value at inner boundary
	    if(flags & INVERT_IN_SYM) {
	      // Use anti-symmetric boundary to set zero-value

	      // Zero-gradient for first point(s)
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }

	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      }else {
		// Antisymmetric across boundary between cells
		avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      }

	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(0.,0.);
	      }
	    }
	  }

	  // Outer boundary
	  if(flags & INVERT_DC_OUT_GRAD) {
	    // Zero gradient at outer boundary

	    if((flags & INVERT_OUT_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient

	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      // Symmetric on last point
	      int ix = xbndry-1;
	      avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;

	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_OUT_SET) {
	    // Setting the values in the outer boundary
	    for(int ix=0;ix<xbndry;ix++) {
	      avec[ncx-ix] = 0.0;
	      bvec[ncx-ix] = 1.0;
	      cvec[ncx-ix] = 0.0;
	      bk1d[ncx-ix] = xk[ncx-ix][iz];
	    }
	  }else {
	    // Zero value at outer boundary
	    if(flags & INVERT_OUT_SYM) {
	      // Use anti-symmetric boundary to set zero-value

	      // Zero-gradient for first point(s)
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      int ix = xbndry-1;
	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      }else {
		// Antisymmetric across boundary between cells
		avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	}else {
	  // AC

	  // Inner boundary
	  if(flags & INVERT_AC_IN_GRAD) {
	    // Zero gradient at inner boundary

	    if((flags & INVERT_IN_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient

	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }
	      // Symmetric on last point
	      avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_IN_SET) {
	    // Setting the values in the boundary
	    for(int ix=0;ix<xbndry;ix++) {
	      avec[ix] = 0.0;
	      bvec[ix] = 1.0;
	      cvec[ix] = 0.0;
	      bk1d[ix] = xk[ix][iz];
	    }
	  }else if(flags & INVERT_AC_IN_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (int ix=0;ix<xbndry;ix++) {
	      avec[ix] = 0.0;
	      bvec[ix] = -1.0;
	      cvec[ix] = exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
	    }
	  }else {
	    // Zero value at inner boundary

	    if(flags & INVERT_IN_SYM) {
	      // Use anti-symmetric boundary to set zero-value

	      // Zero-gradient for first point(s)
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }

	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      }else {
		// Antisymmetric across boundary between cells
		avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      }

	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(0.,0.);
	      }
	    }
	  }

	  // Outer boundary
	  if(flags & INVERT_AC_OUT_GRAD) {
	    // Zero gradient at outer boundary

	    if((flags & INVERT_OUT_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient

	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      // Symmetric on last point
	      int ix = xbndry-1;
	      avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;

	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_AC_OUT_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (int ix=0;ix<xbndry;ix++) {
	      avec[ncx-ix] = exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*mesh->dx[ncx-ix][jy]);;
	      bvec[ncx-ix] = -1.0;
	      cvec[ncx-ix] = 0.0;
	    }
	  }else if(flags & INVERT_OUT_SET) {
	    // Setting the values in the outer boundary
	    for(int ix=0;ix<xbndry;ix++) {
	      avec[ncx-ix] = 0.0;
	      bvec[ncx-ix] = 1.0;
	      cvec[ncx-ix] = 0.0;
	      bk1d[ncx-ix] = xk[ncx-ix][iz];
	    }
	  }else {
	    // Zero value at outer boundary

	    if(flags & INVERT_OUT_SYM) {
	      // Use anti-symmetric boundary to set zero-value

	      // Zero-gradient for first point(s)
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      int ix = xbndry-1;
	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      }else {
		// Antisymmetric across boundary between cells
		avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	}

	// Call tridiagonal solver
	tridag(avec, bvec, cvec, bk1d, xk1d, mesh->ngx);

	if((flags & INVERT_IN_SYM) && (xbndry > 1)) {
	  // (Anti-)symmetry on inner boundary. Nothing to do if only one boundary cell
	  int xloc = 2*xbndry;
	  if(!mesh->BoundaryOnCell)
	    xloc--;

	  if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	    // Inner gradient zero - symmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ix] = xk1d[xloc-ix];
	  }else {
	    // Inner value zero - antisymmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ix] = -xk1d[xloc-ix];
	  }
	}
	if((flags & INVERT_OUT_SYM) && (xbndry > 1)) {
	  // (Anti-)symmetry on outer boundary. Nothing to do if only one boundary cell

	  int xloc =  mesh->ngx - 2*xbndry;
	  if(mesh->BoundaryOnCell)
	    xloc--;

	  if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	    // Outer gradient zero - symmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ncx-ix] = xk1d[xloc + ix];
	  }else {
	    // Outer value zero - antisymmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ncx-ix] = -xk1d[xloc + ix];
	  }
	}
      } else {
	// Periodic in X, so no boundaries
	cyclic_tridag(avec+2, bvec+2, cvec+2, bk1d+2, xk1d+2, mesh->ngx-4);

	// Copy boundary regions
	for(int ix=0;ix<2;ix++) {
	  xk1d[ix] = xk1d[mesh->ngx-4+ix];
	  xk1d[mesh->ngx-2+ix] = xk1d[2+ix];
	}
      }

      if((flags & INVERT_KX_ZERO) && (iz == 0)) {
        dcomplex offset(0.0);
        for(int ix=0;ix<=ncx;ix++)
          offset += bk1d[ix];
        offset /= (BoutReal) (ncx+1);
        for(int ix=0;ix<=ncx;ix++)
          bk1d[ix] -= offset;
      }

      // Fill xk

      for (int ix=0; ix<=ncx; ix++){
	xk[ix][iz]=xk1d[ix];
      }
    }
  }

  // Done inversion, transform back

  for(int ix=0; ix<=ncx; ix++){

    if(flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;

    ZFFT_rev(xk[ix], mesh->zShift[ix][jy], x[ix]);

    x[ix][mesh->ngz-1] = x[ix][0]; // enforce periodicity
  }

  return 0;
}




int invert_laplace_gui4(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d)
{
  if(mesh->NXPE == 1) {
    // Just use the serial code
    return invert_laplace_ser(b, x, flags, a, c, d);
  }


  return 0;
}


int invert_laplace_gui3(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  int jy, jy2;
  FieldPerp xperp;
  int ret;

  Timer timer("invert");

  x.allocate();

  int ys = mesh->ystart, ye = mesh->yend;

  if(MYPE_IN_CORE == 0) {
    // NOTE: REFINE THIS TO ONLY SOLVE IN BOUNDARY Y CELLS
    ys = 0;
    ye = mesh->ngy-1;
  }

  if((mesh->NXPE == 1) || invert_low_mem) {

    for(jy=ys; jy <= ye; jy++) {
      if((flags & INVERT_IN_SET) || (flags & INVERT_OUT_SET))
	xperp = x.slice(jy); // Using boundary values

      if((ret = invert_laplace_gui4(b.slice(jy), xperp, flags, a, c, d)))
	return(ret);
      x = xperp;
    }
  }


  x.setLocation(b.getLocation());

  return 0;
}

const Field3D invert_laplace_gui2(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  output.write("\tCheck point 5\n");
  invert_laplace_gui3(b, x, flags, a, c, d);
  return x;
}







void laplace_tridag_coefs_gui2(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c,
                          const Field2D *ccoef, const Field2D *d)
{
  BoutReal coef1, coef2, coef3, coef4, coef5, kwave;

  kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]

  coef1=mesh->g11[jx][jy];     ///< X 2nd derivative coefficient
  coef2=mesh->g33[jx][jy];     ///< Z 2nd derivative coefficient
  coef3=2.*mesh->g13[jx][jy];  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  if(laplace_all_terms) {
    coef4 = mesh->G1[jx][jy]; // X 1st derivative
    coef5 = mesh->G3[jx][jy]; // Z 1st derivative
  }

  if(d != (Field2D*) NULL) {
    // Multiply Delp2 component by a factor
    coef1 *= (*d)[jx][jy];
    coef2 *= (*d)[jx][jy];
    coef3 *= (*d)[jx][jy];
    coef4 *= (*d)[jx][jy];
    coef5 *= (*d)[jx][jy];
  }

  laplace_nonuniform = true;
  // output.write("\t&&&&&&&&&&&&&&&&&&&Check point 15, nonuniform =true\n");
  if(laplace_nonuniform) {
    // non-uniform mesh correction
    //  output.write("\t&&&&&&&&&&&&&&&&&&&Check point 15, nonuniform =true\n");
    if((jx != 0) && (jx != (mesh->ngx-1))) {
      //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
      coef4 -= 0.5*((mesh->dx[jx+1][jy] - mesh->dx[jx-1][jy])/SQ(mesh->dx[jx][jy]))*coef1; // BOUT-06 term
    }
  }

  if(ccoef != NULL) {
    // A first order derivative term

    if((jx > 0) && (jx < (mesh->ngx-1)))
      coef4 += mesh->g11[jx][jy] * ((*ccoef)[jx+1][jy] - (*ccoef)[jx-1][jy]) / (2.*mesh->dx[jx][jy]*((*ccoef)[jx][jy]));
  }

  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[jx][jy] * mesh->IntShiftTorsion[jx][jy] * mesh->IntShiftTorsion[jx][jy];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }

  coef1 /= SQ(mesh->dx[jx][jy]);
  coef3 /= 2.*mesh->dx[jx][jy];
  coef4 /= 2.*mesh->dx[jx][jy];

  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}


//**********************************************************************************end gui2
