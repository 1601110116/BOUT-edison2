/*******************************************************************************
 * High-Beta gyrofluid 4-field
 * see glf-4field.pdf
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * Based on previous elm-pb code
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <interpolation.hxx>
#include <derivs.hxx>
#include <math.h>
#include <msg_stack.hxx>

// 2D inital profiles
Field2D J0, P0, n0, T0; // Current and pressure
Vector2D b0xcv; // Curvature term
Field2D beta, gradparB;   // Used for Vpar terms
Field2D phi0;   // When diamagnetic terms used
Field2D U0, Psixy, x;   //0th vorticity of equilibrium flow,
                       //radial flux coordinate, normalized radial flux coordinate
BoutReal ny;

// B field vectors
Vector2D B0vec; // B0 field vector

// V0 field vectors
Vector2D V0net; //net flow

// 3D evolving variables
Field3D U, Psi, P, Vpar, Lambdai;

// Derived 3D variables
Field3D Un, Jpar, phi, Vipar; // Parallel current, electric potential

Field3D Jpar2,Grad2p; //  Delp2 of Parallel current

// Parameters
BoutReal density; // Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
BoutReal dnorm; // For diamagnetic terms: 1 / (2. * wci * Tbar)
BoutReal dia_fact; // Multiply diamagnetic term by this
BoutReal delta_i; // Normalized ion skin depth

BoutReal diffusion_p4;   //M: 4th Parallel pressure diffusion

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

BoutReal diffusion_u4;  //4th order perp diffusion in vorticity
Field3D Delp2u;

// options
bool include_curvature, include_jpar0, compress0;
bool evolve_pressure;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar;
int NX,NY,NZ,MZ;
bool lasttimestep;

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

//gyrofluid and Parallel ion motion
Field3D  gyrophi, gyron;   // gyroaveraged potential, gyroaveraged ion density
Field2D P_G0, n_iG0, U0_sour, gyroU0;
Field3D  gyropar, gyrosour, diasour,ni1, Pe1, ns; //invert parameter,perturbed ion denisty, perturbed electron pressure
Field3D gyrophi0, gyrophi_d, gyrophi0_d; // gyroaveraged equilibrium potential, phi-gyrophi, phi0-gyrophi0
Field2D rhoi, invertpar;
BoutReal gynorm; //Normalisation: e*Va*Lbar*Bbar/(2T0)
BoutReal omegai; //temperature,Ti_perp=Ti_parallel=Te_perp=Te_parallel, ion gyroradius and frequency
bool gyrofluid;  //use gyrofluid or not
bool electron_pressure; //include electro pressure in Ohm's law
bool FLR_terms;         //keep gyroaveraged phi and phi0
Field2D gyroa,gyrod, gyroe, gyrof; //invert parameters
int gyron_flags, gyrophi_flags, gyrophi0_flags,gyroPsi_flags, gyroVipar_flags, U0_flags, Ugdc_flags; //invert flags
Vector2D sour1, Grad_n0, Grad_phi0;
Vector3D sour2;
Field3D sour3;
Field3D gyroPsi, gyroPsi_d;   //gyro-averaged Psi and gyroPsi_d = gyroPsi - Psi
bool GradB_terms;       //Keep Grad(B) terms
bool Parallel_ion;      //with parallel ion motion or not
Field3D gyroVipar, gyroVipar_d; //gyro-averaged ion parallel velocity and gyroVipar = gyroVipar - Vipar
bool current_P;     //parallel current in pressure term

//gyroviscosity
bool gyroviscos;
Field2D delp2P0, delp2phi0;
Field3D delp2phi,delp2P, exblinear0, exblinear, exbnonlinear;

//density and temperature control
bool profile_control;
bool fit_pressure;
BoutReal P_min,PDped,Pdel,P_x0,Pbottom,Pa;
bool constant_density;  //using constant density or in vorticity
BoutReal T0_const, n0_const;                //T0=Ti=Te, constant temperature or constant density
bool T0_profile;
BoutReal T_0, T_s, T_min, T_x0;

//electron inertial
bool electron_inertial;  // keep electron inertial term in Ohm's law
Field3D myLambda;        // Lambda = Psi-Me/Mi*Jpar/(gynorm*gynorm*n0)
BoutReal massratio;      // Me/Mi
Field2D Psia,Psid;      //inversion parameters for Psi
int Psi_flags;

//calculate dc part of phi separatively
bool separate_phi;
Field3D phi_dc, p_dc, ni1_dc, Ug_sour, Ug_dc;  //dc part of phi and P

//for pellet simulation
bool include_pellet, finite_perturbation;
Field3D P_pellet, P_pelletx,P_pellety,P_pelletz;
BoutReal plt_x0,plt_xw,plt_y0,plt_yw,plt_z0,plt_zw,plt_P0;
bool DIIID_varyped5;
Field2D newp0;

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

BoutReal filter_weight;
bool smooth_j_x;  // Smooth Jpar in the x direction

int jpar_bndry_width; // Zero jpar in a boundary region
int U_bndry_width; // Zero jpar in a boundary region

bool parallel_lr_diff; // Use left and right shifted stencils for parallel differences
bool parallel_lagrange; // Use (semi-) Lagrangian method for parallel derivatives
bool parallel_project;  // Use Apar to project field-lines

Field3D Xip_x, Xip_z;     // Displacement of y+1 (in cell index space)
Field3D Xim_x, Xim_z;     // Displacement of y-1 (in cell index space)

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

const BoutReal PI = 3.1415927;
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass
const BoutReal Me = 9.1094e-31;     // Electron mass

// Communication objects
FieldGroup comms;

const Field3D Grad2_par2new(const Field3D &f); //for 4th order diffusion

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

  output.write("Solving high-beta gyrofluid 4-field equations\n");
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
  mesh->get(NX,    "nx");
  mesh->get(NY,    "ny");

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

  OPTION(options, gyrofluid,          false);
  OPTION(options, electron_pressure,  false);
  OPTION(options, FLR_terms,          false);
  OPTION(options, GradB_terms,        false);
  OPTION(options, Parallel_ion,       false);
  OPTION(options, current_P,          false); //parallel current term in presure equation

  OPTION(options, gyroviscos,         false);
  OPTION(options, electron_inertial,  false);  //electron inertial


  OPTION(options, profile_control,    false);  //change density and temperature profiles
  OPTION(options, constant_density,   false);
  OPTION(options, n0_const,           10e19);   //density m^-3
  OPTION(options, T0_const,            1000);   //temperature in eV

  OPTION(options, fit_pressure,       false);  //using fit pressure or not, only for cbm girds
  OPTION(options, P_min,              950.0);
  OPTION(options, Pa,                 0.002);
  OPTION(options, Pdel,                22.0);
  OPTION(options, P_x0,                0.85);
  OPTION(options, PDped,              0.062);
  OPTION(options, Pbottom,            948.0);

  OPTION(options, T0_profile,         false);  //non-constant temperature
  OPTION(options, T_0,                    0);
  OPTION(options, T_s,                    0);
  OPTION(options, T_min,                  0);
  OPTION(options, T_x0,                   0);

  OPTION(options, include_pellet,     false);  //include pellet modification
  OPTION(options, finite_perturbation,false);
  OPTION(options, plt_P0,               1.0);
  OPTION(options, plt_x0,               0.5);
  OPTION(options, plt_xw,               1.0);
  OPTION(options, plt_y0,               0.5);
  OPTION(options, plt_yw,               1.0);
  OPTION(options, plt_z0,               0.5);
  OPTION(options, plt_zw,               1.0);
  OPTION(options, DIIID_varyped5,     false);

  OPTION(options, evolve_jpar,       false);  // If true, evolve J raher than Psi

  // Effects to include/exclude
  OPTION(options, include_curvature, true);
  OPTION(options, include_jpar0,     true);
  OPTION(options, nogradparj,       false);  // exclude (B0^2)*Grad_par(Jpar)?
  OPTION(options, evolve_pressure,   true);

  OPTION(options, compress0,         false);
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

  //net EXB flow
  OPTION(options, withflow,         false);    //withflow or not
  OPTION(options, K_H_term,         false);    //keep K-H term
  OPTION(options, D_0,                0.0);    // velocity magnitude
  OPTION(options, D_s,                0.0);    // flowshear
  OPTION(options, x0,                 0.0);    //flow location
  OPTION(options, sign,               1.0);    //flow direction, -1 means negative electric field
  OPTION(options, D_min,           3000.0);    //a constant

  OPTION(options, noshear,           false);

  OPTION(options, relax_j_vac,       false); // Relax vacuum current to zero
  OPTION(options, relax_j_tconst,    0.1);

  OPTION(options, separate_phi,     false);

  // Toroidal filtering
  OPTION(options, filter_z,          false);  // Filter a single n
  OPTION(options, filter_z_mode,     1);
  OPTION(options, low_pass_z,       -1);      // Low-pass filter
  OPTION(options, zonal_flow,       -1);      // zonal flow filter
  OPTION(options, zonal_field,      -1);      // zonal field filter
  OPTION(options, zonal_bkgd,       -1);      // zonal background P filter

  // Radial smoothing
  OPTION(options, filter_weight,       -1);
  OPTION(options, smooth_j_x,       false);  // Smooth Jpar in x

  // Jpar boundary region
  OPTION(options, jpar_bndry_width, -1);
  OPTION(options, U_bndry_width,    -1);

  // Parallel differencing
  OPTION(options, parallel_lr_diff, false);
  OPTION(options, parallel_lagrange, false); // Use a (semi-) Lagrangian method for Grad_parP
  OPTION(options, parallel_project, false);

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
  OPTION(options, diffusion_par,      -1.0);  // Parallel pressure diffusion
  OPTION(options, diffusion_p4,       -1.0);  // M: 4th Parallel pressure diffusion

  OPTION(options, diffusion_u4,      -1.0);

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


  // Field inversion flags
  OPTION(options, phi_flags,         0);
  OPTION(options, apar_flags,        0);
  OPTION(options, gyron_flags,       0);
  OPTION(options, gyrophi_flags,     0);
  OPTION(options, gyrophi0_flags,    0);
  OPTION(options, Psi_flags,         0);
  OPTION(options, gyroPsi_flags,     0);
  OPTION(options, gyroVipar_flags,   0);
  OPTION(options, U0_flags,          0);
  OPTION(options, Ugdc_flags,        0);

  OPTION(options, lasttimestep,  false);
  OPTION(globalOptions, MZ,          1);

  x=(Psixy-Psiaxis)/(Psibndry-Psiaxis);

  Dphi0=-D_min-0.5*D_0*(1.0-tanh(D_s*(x-x0)));

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

  //density and temperature profiles control
  if(profile_control){

    if(fit_pressure){                  //only for cbm grids
      //  P0 = 0.5*P_0*(1-tanh(P_s*(x-P_x0)))+P_min-21000*log(x+1.7)-31000*log(5.5-x)+67600;
      P0=P_min*(exp(-(P_x0-x)/PDped)+(1.0+Pa*(P_x0-x)/PDped)*exp((P_x0-x)/PDped)*Pdel)/(exp(-(P_x0-x)/PDped)+exp((P_x0-x)/PDped))-Pbottom;
      output.write("warning: using a fit function for pressure!!!!!!!!!!!!!!!!!!!\n");
      }

    if(DIIID_varyped5){  //only for DIIID varyped5 equilibrium,ny=128
      newp0=0.0;
      SAVE_ONCE(newp0);
      for(int jx=0;jx<mesh->ngx;jx++)
	for(int jy=0;jy<mesh->ngy;jy++){
          if(mesh->YGLOBAL(jy)>7&&mesh->YGLOBAL(jy)<120){
             newp0[jx][jy] = 37441800.0*x[jx][jy]*x[jx][jy]*x[jx][jy]-107703720.0*x[jx][jy]*x[jx][jy]+103023464.0*x[jx][jy]-32761272.0+180.0*exp(10.8*tanh(19.9*(x[jx][jy]-1.0127)));
             if(mesh->XGLOBAL(jx)>208)
	       newp0[jx][jy] = 134.08;
          }
          else
             newp0[jx][jy] = 134.07;
     	  mesh->communicate(newp0);
        }
      P0 = newp0;
    }


    if(constant_density){
      n0 = n0_const;
      T0 = 0.5*P0/n0;       //K
      output.write("constant n0 = %e m^3\n  ", n0_const);
      }else{
	if(T0_profile){
	  T0 = 1.0;
	  T0 = 0.5*T_0*(1-tanh(T_s*(x-T_x0)))+T_min;//-21000*log(x+1.7)-31000*log(5.5-x)+67600;
          output.write("warning: using a fit function for temperature!!!!!!!!!!!!!!!!!!!\n");
	}else{
         T0 = T0_const;        //eV
         output.write("constant T0 = %e eV\n ", T0_const);}

         T0*=1.602e-19;          //K
         n0 = 0.5*P0/T0;
      }
  }else{                     //temperature and density from grid
    mesh->get(n0, "Ni0");    // *10^20 m^-3
    mesh->get(T0, "Ti0");    //eV
    n0*=1.0e20;              // m^-3
    T0*=1.602e-19;           //K
    output.write("read density and temperature from grid\n");
  }


  if(include_pellet){
    P_pellet = 0.0;
    P_pelletx = 0.0;
    P_pellety = 0.0;
    P_pelletz = 0.0;

    P_pelletx = exp(-0.5*(x-plt_x0)*(x-plt_x0)/(plt_xw*plt_xw));
    Field2D pol_angle;
    mesh->get(pol_angle, "dy");
    int zperiod;
    globalOptions->get("zperiod", zperiod, 1);
    for(int jx=0;jx<mesh->ngx;jx++)
      for(int jz=0;jz<mesh->ngz;jz++)
        for(int jy=0;jy<mesh->ngy;jy++){
          BoutReal p_angle=pol_angle[jx][jy] * ((BoutReal) mesh->YGLOBAL (jy));
	  BoutReal y = sin(0.5*p_angle-plt_y0*PI);
	   P_pellety[jx][jy][jz] = exp(-0.5*y*y/(plt_yw*plt_yw));
    }

    for(int jx=0;jx<mesh->ngx;jx++)
      for(int jy=0;jy<mesh->ngy;jy++)
        for(int jz=0;jz<mesh->ngz;jz++){
          BoutReal t_angle = ((BoutReal) jz) * mesh->dz*zperiod;
	  BoutReal z = sin(0.5*t_angle-plt_z0*PI);
	  P_pelletz[jx][jy][jz] = exp(-0.5*z*z/(plt_zw*plt_zw));
     }

    P_pellet = plt_P0*P_pelletx*P_pellety*P_pelletz;
    P_pellet = P_pellet.shiftZ(false);
    P_pellet = 2.0*MU0 * P_pellet / (Bbar*Bbar);
    SAVE_ONCE4(P_pelletx,P_pellety,P_pelletz,P_pellet);
  }

  T0=T0/(1.602e-19);
  SAVE_ONCE(T0);    //in eV

  T0*=1.602e-19;
  rhoi = Mi*sqrt(T0/Mi)/(1.602e-19*B0);   //local rhoi
  omegai = 1.602e-19*Bbar/Mi;
  T0=2.*T0/(Mi*Va*Va);                       //normalized T0
  n0=n0/density;                              //normalized n0
  gynorm = omegai*Tbar;
  invertpar = rhoi*rhoi/(Lbar*Lbar);
  gyroa = 1.0;                            //invert parameters
  gyrod = -0.5*invertpar;
  gyroe = 2./invertpar;
  gyrof = -0.25*invertpar;
  SAVE_ONCE(gyrod);
  SAVE_ONCE(gyroe);

  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Va = %e m/s   Tbar = %e s\n", Va, Tbar);
  output.write("                dnorm = %e\n", dnorm);
  output.write("omegai = %e Hz\n",omegai/(2*3.1415));
  output.write("gynorm = %e   ", gynorm);

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

  //M: 4th order diffusion of p
  if(diffusion_p4 > 0.0) {
    output.write("    diffusion_p4: %e\n", diffusion_p4);
    dump.add(diffusion_p4, "diffusion_p4", 1);
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

  if(withflow&&K_H_term)
    output.write("    keep K-H term\n");
  else
    output.write("   drop K-H term\n");

  if(gyrofluid)
    output.write("    using gyrofluid\n");

  if(electron_pressure)
    output.write("    keep electron pressure\n");

  J0 = - MU0*Lbar * J0 / B0;
  P0 = 2.0*MU0 * P0 / (Bbar*Bbar);
  V0=V0/Va;
  Dphi0*=Tbar;

  massratio = Me/Mi;                     //mass ratio
  Psia = 1.0;                            //inversion parameters
  Psid = -massratio/(gynorm*gynorm*n0);  //
  output.write("massratio = %e\n", massratio);

  b0xcv.x /= Bbar;
  b0xcv.y *= Lbar*Lbar;
  b0xcv.z *= Lbar*Lbar;

  Rxy  /= Lbar;
  Bpxy /= Bbar;
  Btxy /= Bbar;
  B0   /= Bbar;
  hthe /= Lbar;
  mesh->dx   /= Lbar*Lbar*Bbar;
  I    *= Lbar*Lbar*Bbar;

  BoutReal pnorm = max(P0, true); // Maximum over all processors

  vacuum_pressure *= pnorm; // Get pressure from fraction
  vacuum_trans *= pnorm;

  // Transitions from 0 in core to 1 in vacuum
  vac_mask = (1.0 - tanh( (P0 - vacuum_pressure) / vacuum_trans )) / 2.0;

  eta = core_resist + (vac_resist - core_resist) * vac_mask;
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

  /**************** SET VARIABLE LOCATIONS *************/

  P.setLocation(CELL_CENTRE);
  U.setLocation(CELL_CENTRE);
  phi.setLocation(CELL_CENTRE);
  Psi.setLocation(CELL_YLOW);
  Lambdai.setLocation(CELL_CENTRE);
  Jpar.setLocation(CELL_YLOW);
  Vpar.setLocation(CELL_YLOW);
  gyropar.setLocation(CELL_YLOW);
  gyrophi.setLocation(CELL_YLOW);
  gyron.setLocation(CELL_YLOW);
  gyrophi_d.setLocation(CELL_YLOW);
  gyrophi0_d.setLocation(CELL_YLOW);
  gyrophi0.setLocation(CELL_YLOW);
  sour2.setLocation(CELL_YLOW);
  sour3.setLocation(CELL_YLOW);
  myLambda.setLocation(CELL_CENTRE);
  gyrosour.setLocation(CELL_YLOW);
  Vipar.setLocation(CELL_YLOW);
  gyroVipar.setLocation(CELL_YLOW);
  gyroVipar_d.setLocation(CELL_YLOW);
  gyroPsi.setLocation(CELL_YLOW);
  gyroPsi_d.setLocation(CELL_YLOW);

  delp2phi0.setLocation(CELL_YLOW);
  delp2phi.setLocation(CELL_YLOW);
  delp2P0.setLocation(CELL_YLOW);
  delp2P.setLocation(CELL_YLOW);
  exblinear0.setLocation(CELL_YLOW);
  exblinear.setLocation(CELL_YLOW);
  exbnonlinear.setLocation(CELL_YLOW);
  Delp2u.setLocation(CELL_YLOW);

  U0_sour.setLocation(CELL_YLOW);
  gyroU0.setLocation(CELL_YLOW);
  P_G0.setLocation(CELL_YLOW);
  n_iG0.setLocation(CELL_YLOW);
  phi_dc.setLocation(CELL_YLOW);
  p_dc.setLocation(CELL_YLOW);
  ni1_dc.setLocation(CELL_YLOW);
  Ug_sour.setLocation(CELL_YLOW);
  Ug_dc.setLocation(CELL_YLOW);

  P_pellet.setLocation(CELL_YLOW);
  P_pelletx.setLocation(CELL_YLOW);
  P_pellety.setLocation(CELL_YLOW);
  P_pelletz.setLocation(CELL_YLOW);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(U);
  SOLVE_FOR(P);

  if(electron_inertial){
    output.write("Keep electron inertial in Ohm's law\n");
    SOLVE_FOR(myLambda);
    dump.add(Psi, "Psi", 1);
    dump.add(Jpar, "jpar", 1);
  }
  else{
    output.write("Solving for Psi, Differentiating to get jpar\n");
    SOLVE_FOR(Psi);
    dump.add(Jpar, "jpar", 1);
  }

  if(Parallel_ion){
     SOLVE_FOR(Lambdai);
     output.write("4-field with ion parallel motion\n");
     dump.add(Vipar, "Vipar", 1);
     if(gyrofluid){
       dump.add(gyroVipar, "gyroVipar", 1);
       dump.add(gyroPsi,   "gyroPsi",   1);
      }
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

    dump.add(phi, "phi", 1);

  // Diamagnetic phi0,  not devided by B0!!!!!!!!
    if(constant_density)      //constant density
      phi0 = -0.5*dnorm*P0;
    else                     //non-constant density, valid for isothermal model
      phi0 = -0.25*T0/(gynorm)*log(P0);
    SAVE_ONCE(phi0);

  delp2phi0.setBoundary("delp2phi0");
  delp2P0.setBoundary("delp2P0");

  delp2P0 = 0.0;
  delp2phi0 = 0.0;
  U0_sour = 0.0;
  P_G0 = 0.0;
  n_iG0 = 0.0;
  gyroU0 = 0.0;

  delp2P0 = Delp2(P0);
  delp2P0.applyBoundary();
  mesh->communicate(delp2P0);

  delp2phi0 = Delp2(phi0);
  delp2phi0.applyBoundary();
  mesh->communicate(delp2phi0);

  if(gyrofluid){
    dump.add(gyrophi, "gyrophi", 1);
    dump.add(gyrophi0, "gyrophi0", 1);
    dump.add(gyron,   "gyron",   1);
    dump.add(ni1,     "ni1",     1);
    SAVE_ONCE4(P_G0, gyroU0, U0_sour, n_iG0);

    U0_sour = -0.5/T0*delp2P0;
    //  gyroU0 = invert_laplace(U0_sour,U0_flags,&gyroe,NULL,NULL);
    gyroU0 = U0_sour*gynorm*B0/gyroe;
    mesh->communicate(gyroU0);

    P_G0 = P0-T0*gyroU0/(gynorm*B0);
    n_iG0 = 0.5*(P_G0/T0-gyroU0/(gynorm*B0));
  }

  if(electron_pressure)
  dump.add(Pe1,     "Pe1",     1);


  phi_dc = 0.0;
  Ug_dc = 0.0;
  Ug_sour = 0.0;
  p_dc = 0.0;
  ni1_dc = 0.0;
  if(separate_phi){
   dump.add(p_dc, "p_dc", 1);
   dump.add(phi_dc, "phi_dc", 1);
   dump.add(Ug_dc, "Ug_dc", 1);
   dump.add(Ug_sour, "Ug_sour", 1);
   p_dc.setBoundary("p_dc");
   Ug_sour.setBoundary("Ug_sour");
   phi_flags = 784;       //remove dc phi in invert_laplace
   output.write("caculate dc part of phi separately, and phi_flags = %d\n", phi_flags);
  }


  // Add some equilibrium quantities and normalisations
  // everything needed to recover physical units
  SAVE_ONCE2(J0, P0);
  SAVE_ONCE4(density, Lbar, Bbar, Tbar);
  SAVE_ONCE2(Va, B0);
  SAVE_ONCE2(Dphi0, V0);
  SAVE_ONCE2(rhoi, n0);

  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero
      gyropar = 0.0;
    ns = 0.0;
    gyrosour = 0.0;

  if(!restarting&&!lasttimestep) {
    // Only if not restarting: Check initial perturbation

    if(include_pellet&&finite_perturbation)
      P = P_pellet;

    //Psi should be consistent with myLambda
   if(electron_inertial)
     Psi = invert_laplace(myLambda,Psi_flags,&Psia,NULL,&Psid);

    // Set U to zero where P0 < vacuum_pressure
   // U = where(P0 - vacuum_pressure, U, 0.0);

    // Phi should be consistent with U
  if(electron_pressure)
    Pe1 = 0.5*(P+T0*U/(gynorm*B0));
  if(gyrofluid){
    if(nonlinear&&separate_phi){
      p_dc = filter(P, 0);        //get dc part of pressure
      mesh->communicate(p_dc);
      p_dc.applyBoundary();
      p_dc = smooth_x(p_dc);

      Ug_sour = -0.25/T0*invertpar*Delp2(p_dc);   //get dc part of gyrofluid vorticity
      mesh->communicate(Ug_sour);
      Ug_sour.applyBoundary();
      Ug_sour = smooth_x(Ug_sour);
      Ug_dc = invert_laplace(Ug_sour, Ugdc_flags, &gyroa, NULL, &gyrof);
      mesh->communicate(Ug_dc);
      Ug_dc *= (B0*gynorm);
    }

    ni1 = 0.5*(P/T0-(U+Ug_dc)/(gynorm*B0));

    //first step invert: get gyroaveraged density
    gyron = invert_laplace(ni1,gyron_flags,&gyroa,NULL,&gyrod);
    mesh->communicate(gyron);

   //sencond step invert: get potential
    if(constant_density){
      ns = T0*(gyron-ni1-U/(gynorm*B0))/n0;
      gyrosour = -ns/invertpar;
      gyropar = invert_laplace(gyrosour,phi_flags,NULL);
      mesh->communicate(gyropar);
      phi=0.5*gyropar/gynorm;
      mesh->communicate(phi);
    }
    else{
      ns = T0*(gyron-ni1-(U+Ug_dc)/(gynorm*B0))/n0;
      ns = smooth_x(ns);
      sour1 = Grad(log(n0), CELL_YLOW);
      sour2 = Grad(ns, CELL_YLOW);
      gyrosour = -ns-invertpar*sour1*sour2;
      mesh->communicate(gyrosour);
      gyropar = invert_laplace(gyrosour, gyrophi_flags, NULL, &n0, NULL);
      mesh->communicate(gyropar);
      phi = 0.5*(gyropar/invertpar+ns)/gynorm;    //phi include dc part in nonlinear simulations

      if(nonlinear&&separate_phi){
        phi_dc = filter(phi, 0);    //linear dc part of phi, need to be removed
        mesh->communicate(phi_dc);
        phi -= phi_dc;              //remove linear dc phi

        ni1_dc = filter(ni1, 0);
        mesh->communicate(ni1_dc);
        phi_dc = -0.25*T0/gynorm*log(1.+2.*T0*ni1_dc/P_G0);  //get nonlinear dc phi
      }
      phi += phi_dc;
    }

    //third step invert: get gyroaveraged potential
    if(FLR_terms){
      gyrophi = invert_laplace(phi,gyrophi_flags,&gyroa,NULL,&gyrod);
      mesh->communicate(gyrophi);
      gyrophi_d = gyrophi-phi;

      if(diamag_phi0){
        gyrophi0 = invert_laplace(phi0,gyrophi0_flags,&gyroa,NULL,&gyrod);
        mesh->communicate(gyrophi0);
        gyrophi0_d = gyrophi0-phi0;
      }
    }
 }
  else{
    if(constant_density){
     Un = B0*U/n0;
     phi = invert_laplace(Un, phi_flags, NULL);
     mesh->communicate(phi);
     if(diamag){
       phi -= 0.5*dnorm * P / n0;
      }
    }
    else{
      if(nonlinear&&separate_phi){
        p_dc = filter(P, 0);        //get dc part of pressure
        mesh->communicate(p_dc);
        p_dc.applyBoundary();
        phi_dc = -0.25*T0/gynorm*log(1.+p_dc/P0);
       }

      if(diamag){
        sour3 = 0.0;
	sour3 = Delp2(P);
	mesh->communicate(sour3);
        diasour = B0*U/n0-0.5*dnorm/n0*sour3;
        phi = invert_laplace(diasour, phi_flags, NULL, &n0, NULL);
        mesh->communicate(phi);
     }
     else{
       diasour = B0*U/n0;
       phi = invert_laplace(diasour, phi_flags, NULL, &n0, NULL);
       mesh->communicate(phi);
     }

      phi += phi_dc;  //get dc part back

    }
  }

  if(Parallel_ion){
    Vipar = Lambdai;
    if(gyrofluid){
      gyroPsi = invert_laplace(Psi,gyroPsi_flags,&gyroa,NULL,&gyrod);
      mesh->communicate(gyroPsi);
      gyroPsi_d = gyroPsi-Psi;
      Vipar -= gynorm*B0*gyroPsi_d;
    }
  }
}


if(lasttimestep&&!restarting)
{
  NZ=MZ-1;
  output.write("Initial Profiles are loaded from .txt files of all evolving quatities at last time step \n");
  BoutReal lstimeP[NX][NY][NZ],lstimeU[NX][NY][NZ],lstimePsi[NX][NY][NZ];
  ifstream pFile1,pFile2,pFile3;
  pFile1.open ("data/lstime_P.txt", ios::in );
  pFile2.open ("data/lstime_U.txt", ios::in );
  pFile3.open ("data/lstime_Psi.txt", ios::in );

  for(int jx=0;jx<NX;jx++)
    for(int jy=0;jy<NY;jy++)
      for(int jz=0;jz<NZ;jz++)
         {
	   pFile1 >> lstimeP[jx][jy][jz];
	   pFile2 >> lstimeU[jx][jy][jz];
	   pFile3 >> lstimePsi[jx][jy][jz];
	 }
  pFile1.close();
  pFile2.close();
  pFile3.close();

  for(int jx=0;jx<mesh->ngx;jx++)
  for(int jy=0;jy<mesh->ngy;jy++)
  for(int jz=0;jz<mesh->ngz;jz++)
  {
     P[jx][jy][jz] = lstimeP[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
     U[jx][jy][jz] = lstimeU[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
     if(electron_inertial)
       myLambda[jx][jy][jz] = lstimePsi[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
     else
       Psi[jx][jy][jz] = lstimePsi[mesh->XGLOBAL (jx)][mesh->YGLOBAL (jy)][jz];
  }
}



  /************** SETUP COMMUNICATIONS **************/

  comms.add(U);
  comms.add(P);

   // Set boundary conditions
  phi.setBoundary("phi");
  gyrophi.setBoundary("gyrophi");
  gyrophi0.setBoundary("gyrophi0");
  gyropar.setBoundary("gyropar");
  sour1.setBoundary("sour1");
  sour2.setBoundary("sour2");
  sour3.setBoundary("sour3");
  diasour.setBoundary("diasour");
  gyron.setBoundary("gyron");
  Grad2p.setBoundary("Grad2p");
  Jpar.setBoundary("J");
  Jpar2.setBoundary("J");
  gyroVipar.setBoundary("gyroVipar");
  gyroPsi.setBoundary("gyroPsi");


  delp2phi.setBoundary("delp2phi");
  delp2P.setBoundary("delp2P");
  exblinear0.setBoundary("exblinear0");
  exblinear.setBoundary("exblinear");
  exbnonlinear.setBoundary("exbnonlinear");

  Delp2u.setBoundary("Delp2u");

  if(electron_inertial){
    comms.add(myLambda);
    Psi.setBoundary("Psi");
  }
  else
    comms.add(Psi);

  if(Parallel_ion)
    comms.add(Lambdai);

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
      result -= bracket(Psi*B0 , f, bm_mag);
    }
  }

  return result;
}

const Field3D Grad_parP_gyro(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
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
      result -= bracket(gyroPsi*B0, f, bm_mag);
    }
  }

  return result;
}

const Field3D Grad_parP_back(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
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
	  result[i][j][k] = (fp[i][j][k] - fm[i][j-2][k])/(2.*mesh->dy[i][j]*sqrt(mesh->g_22[i][j]));
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
      result -= bracket(Psi*B0, f, bm_mag);
    }
  }

  return result;
}

const Field3D Grad_parP_forward(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
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
	  result[i][j][k] = (fp[i][j+2][k] - fm[i][j][k])/(2.*mesh->dy[i][j]*sqrt(mesh->g_22[i][j]));
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
      result -= bracket(Psi*B0, f, bm_mag);
    }
  }

  return result;
}

const Field3D Grad_parP0(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
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
  }

  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t)
{

   mesh->communicate(comms);

  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  if(nonlinear) {
    vac_mask = (1.0 - tanh( ((P0 + P) - vacuum_pressure) / vacuum_trans )) / 2.0;
        eta = core_resist + (vac_resist - core_resist) * vac_mask;
   }

   // Zero U in boundary regions. Prevents vorticity drive at the boundary
  if(U_bndry_width > 0) {
      for(int i=0;i<U_bndry_width;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	     if(mesh->firstX())
	      U[i][j][k] = 0.0;
	    if(mesh->lastX())
	      U[mesh->ngx-1-i][j][k] = 0.0;
	  }
   }

  //get Psi
  if(electron_inertial)
    {
      Psi = invert_laplace(myLambda,Psi_flags,&Psia,NULL,&Psid);
      mesh->communicate(Psi);
    }

  //Three step inversion using Pade approximation
  if(electron_pressure)
    Pe1 = 0.5*(P+T0*U/(gynorm*B0));

  if(gyrofluid){
    if(nonlinear&&separate_phi){
      p_dc = filter(P, 0);        //get dc part of pressure
      mesh->communicate(p_dc);
      p_dc.applyBoundary();
      p_dc = smooth_x(p_dc);

      Ug_sour = -0.25/T0*invertpar*Delp2(p_dc);   //get dc part of gyrofluid vorticity
      mesh->communicate(Ug_sour);
      Ug_sour.applyBoundary();
      Ug_sour = smooth_x(Ug_sour);
      Ug_dc = invert_laplace(Ug_sour, Ugdc_flags, &gyroa, NULL, &gyrof);
      mesh->communicate(Ug_dc);
      Ug_dc *= (B0*gynorm);
     }

    ni1 = 0.5*(P/T0-(U+Ug_dc)/(gynorm*B0));

    //first step invert: get gyroaveraged density
    gyron = invert_laplace(ni1,gyron_flags,&gyroa,NULL,&gyrod);
    mesh->communicate(gyron);

    //sencond step invert: get potential
    if(constant_density){
      ns = T0*(gyron-ni1-U/(gynorm*B0))/n0;
      gyrosour = -ns/invertpar;
      mesh->communicate(gyrosour);
      gyropar = invert_laplace(gyrosour,phi_flags,NULL);
      mesh->communicate(gyropar);
      phi=0.5*gyropar/gynorm;
      mesh->communicate(phi);
    }
    else{
      ns = T0*(gyron-ni1-(U+Ug_dc)/(gynorm*B0))/n0;
      ns = smooth_x(ns);
      sour1 = Grad(log(n0),CELL_YLOW);
      sour2 = Grad(ns,CELL_YLOW);
      sour1.applyBoundary();
      sour2.applyBoundary();
      mesh->communicate(sour1);
      mesh->communicate(sour2);
      gyrosour = -ns-invertpar*sour1*sour2;
      gyrosour = smooth_x(gyrosour);
      gyropar = invert_laplace(gyrosour, phi_flags, NULL, &n0, NULL);
      mesh->communicate(gyropar);
      gyropar = smooth_x(gyropar);
      phi = 0.5*(gyropar/invertpar+ns)/gynorm;
      phi = smooth_x(phi);

      if(nonlinear&&separate_phi){
        phi_dc = filter(phi, 0);    //linear dc part of phi, need to be removed
        mesh->communicate(phi_dc);
        phi -= phi_dc;              //remove linear dc phi

        ni1_dc = filter(ni1, 0);
        mesh->communicate(ni1_dc);
        phi_dc = -0.25*T0/(gynorm*B0)*log(1.+2.*T0*ni1_dc/P_G0);  //get nonlinear dc phi
      }
      phi += phi_dc;
    }

    //third step invert: get gyroaveraged potential
    if(FLR_terms){
      gyrophi = invert_laplace(phi,gyrophi_flags,&gyroa,NULL,&gyrod);
      mesh->communicate(gyrophi);
      gyrophi_d = gyrophi-phi;

      if(diamag_phi0){
        gyrophi0 = invert_laplace(phi0,gyrophi0_flags,&gyroa,NULL,&gyrod);
        mesh->communicate(gyrophi0);
        gyrophi0_d = gyrophi0-phi0;
      }
    }
}
else{
  if(constant_density){
    Un = B0*U/n0;
    phi = invert_laplace(Un, phi_flags, NULL);
    if(diamag) {
      phi -= 0.5*dnorm * P / n0;
     }
   }
   else{
     if(nonlinear&&separate_phi){
       p_dc = filter(P, 0);        //get dc part of pressure
       mesh->communicate(p_dc);
       p_dc.applyBoundary();
       phi_dc = -0.25*T0/gynorm*log(1.+p_dc/P0);   //get dc part of phi
     }

     if(diamag){
       sour3 = Delp2(P);
       sour3.applyBoundary();
       mesh->communicate(sour3);
       if(include_pellet&&!finite_perturbation)
	  diasour =  B0*U/(n0+0.5*P_pellet/T0)-0.5*dnorm/(n0+0.5*P_pellet/T0)*sour3;
       else
          diasour = B0*U/n0-0.5*dnorm/n0*sour3;
       phi = invert_laplace(diasour, phi_flags, NULL, &n0, NULL);
     }
     else{
       diasour = B0*U/n0;
       phi = invert_laplace(diasour, phi_flags, NULL, &n0, NULL);
     }
     phi += phi_dc;          //get dc part back
   }

   mesh->communicate(phi);
}

  // Get J from Psi
  Jpar = Delp2(Psi*B0)/B0;
  Jpar.applyBoundary();
  mesh->communicate(Jpar);

  // Zero j in boundary regions. Prevents vorticity drive at the boundary
  if(jpar_bndry_width > 0) {
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


   // Zero jpar2 in boundary regions. Prevents vorticity drive at the boundary
   if(jpar_bndry_width > 0) {
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
  ////////////////////////////////////////////////////
  // Project perturbed field-lines

  if(parallel_project) {
    Vector3D Btilde; // Perturbed field
    Btilde = Grad(Psi*B0) ^ B0vec/B0;
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

  if(electron_inertial) {

    ddt(myLambda) = 0;
    ddt(myLambda) = massratio/(gynorm*gynorm)*bracket(phi,J0/n0, bm_exb) -Grad_parP0(phi, CELL_CENTRE)/B0 + eta*Jpar;

    if(diamag_phi0)
      ddt(myLambda) += bracket(Psi*B0, phi0, bm_mag)/B0;

    if(nonlinear){
      ddt(myLambda) -= bracket(phi, myLambda, bm_exb);

    }
    if(electron_pressure){
      ddt(myLambda)+=1.0/(4.0*gynorm*n0*B0*B0)*b0xGrad_dot_Grad(P0,Psi*B0);
      ddt(myLambda)+=1.0/(2.0*gynorm*n0*B0)*Grad_parP_back(Pe1, CELL_YLOW);
    }

    // Hyper-resistivity
    if(hyperresist > 0.0)
      ddt(myLambda) -= eta*hyperresist * Delp2(Jpar);

    // electron Hyper-viscosity coefficient
    if(ehyperviscos > 0.0)
      ddt(myLambda) -= eta*ehyperviscos * Delp2(Jpar2);

  }else {   // Evolving vector potential

    ddt(Psi) = -Grad_parP(phi, CELL_CENTRE) / B0  + eta*Jpar;

    if(eHall) {
      ddt(Psi) +=  0.25*delta_i*(Grad_parP(B0*P, CELL_CENTRE) / B0
                                 +b0xGrad_dot_Grad(P0, Psi*B0)/B0);   // electron parallel pressure
    }

    if(diamag_phi0)
      ddt(Psi) += bracket(Psi*B0, phi0, bm_mag)/B0;   // Equilibrium flow

    if(withflow)                                //net flow
      ddt(Psi)-= V_dot_Grad(V0net, Psi*B0)/B0;

    if(electron_pressure){
         ddt(Psi)+=1.0/(4.0*gynorm*n0*B0*B0)*b0xGrad_dot_Grad(P0,Psi*B0);
         ddt(Psi)+=1.0/(2.0*gynorm*n0*B0)*Grad_parP_back(Pe1, CELL_YLOW);
    }


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

  //Inversion for gyroPsi
  if(Parallel_ion&&gyrofluid){
    gyroPsi = invert_laplace(Psi,gyroPsi_flags,&gyroa,NULL,&gyrod);
    gyroPsi.applyBoundary();
    mesh->communicate(gyroPsi);

    gyroPsi_d = gyroPsi-Psi;
    mesh->communicate(gyroPsi_d);
}


  ////////////////////////////////////////////////////
  //ion parallel motion equaiton

  if(Parallel_ion){
    ddt(Lambdai) = 0;

    ddt(Lambdai) -= 0.5/n0*Grad_parP(P, CELL_CENTRE);
    ddt(Lambdai) -= 0.5/n0*b0xGrad_dot_Grad(P0, Psi*B0)/B0;


    if(nonlinear)
      ddt(Lambdai) -= bracket(phi, Lambdai, bm_exb);

    if(gyrofluid&&FLR_terms){
      ddt(Lambdai) -= gynorm*Grad_parP_gyro(gyrophi_d, CELL_YLOW);
      ddt(Lambdai) -= 0.25/n0*b0xGrad_dot_Grad(P0, gyroPsi_d*B0)/B0;
    }

   if(gyrofluid&&FLR_terms&&nonlinear){
      ddt(Lambdai) -= bracket(gyrophi_d, Lambdai, bm_exb);
      ddt(Lambdai) += gynorm*bracket(gyrophi, gyroPsi_d*B0, bm_exb);
      ddt(Lambdai) -= gynorm*bracket(phi, gyroPsi_d*B0, bm_exb);
      ddt(Lambdai) += 0.5*T0/n0*bracket(gyroPsi_d*B0, ni1, bm_mag);
    }

   if(GradB_terms)
      ddt(Lambdai) -= 1./(gynorm*B0*n0)*b0xGrad_dot_Grad(B0, P0*Lambdai*B0);

   if(GradB_terms&&gyrofluid&&FLR_terms)
      ddt(Lambdai) += T0/n0*b0xGrad_dot_Grad(B0, n0*gyroPsi_d*B0);

   Vipar = Lambdai;              //get ion parallel veloctiy

   if(gyrofluid){
     Vipar -= gynorm*B0*gyroPsi_d;

     gyroVipar = invert_laplace(Vipar,gyroVipar_flags,&gyroa,NULL,&gyrod);
     gyroVipar.applyBoundary();
     mesh->communicate(gyroVipar);

     gyroVipar_d = gyroVipar - Vipar;
     mesh->communicate(gyroVipar_d);
   }

}

  ////////////////////////////////////////////////////
  // Vorticity equation

  ddt(U) = B0* b0xGrad_dot_Grad(Psi*B0, J0, CELL_CENTRE); // Grad j term

  ddt(U) += b0xcv*Grad(P);  // curvature term

  if(!nogradparj)
    ddt(U) -= (B0^2)*Grad_parP_forward(Jpar, CELL_CENTRE); // b dot grad j

  if(withflow&&K_H_term)                         //K_H_term
    ddt(U) -= bracket(phi,U0, bm_exb);

  if(diamag_phi0)
    ddt(U) -= bracket(phi0, U, bm_exb);   // Equilibrium flow

  if(withflow)                            // net flow
    ddt(U) -= V_dot_Grad(V0net, U);

  if(nonlinear)
    ddt(U) -= bracket(phi, U, bm_exb);    // Advection

  if(gyrofluid&&FLR_terms)
    ddt(U) += gynorm*B0*bracket(gyrophi_d,n_iG0, bm_exb);

  if(gyrofluid&&diamag_phi0)
    ddt(U) -= bracket(phi, gyroU0, bm_exb);

  if(gyrofluid&&diamag_phi0&&FLR_terms)
    ddt(U) += gynorm*B0*bracket(gyrophi0_d, ni1, bm_exb);

  if(gyrofluid&&nonlinear&&FLR_terms)
    ddt(U) += gynorm*B0*bracket(gyrophi_d,ni1, bm_exb);

  if(gyrofluid&&nonlinear&&separate_phi){
     ddt(U) -= bracket(phi0, Ug_dc, bm_exb);
     ddt(U) -= bracket(phi, Ug_dc, bm_exb);
   }

  if(GradB_terms)
    ddt(U) += 2.*J0*b0xGrad_dot_Grad(B0, Psi*B0);

  if(GradB_terms&&FLR_terms)
    ddt(U) += gynorm*n0*b0xGrad_dot_Grad(B0, gyrophi_d);

  if(Parallel_ion&&gyrofluid&&FLR_terms)
    ddt(U) -= gynorm*B0*B0*Grad_parP(n0*gyroVipar_d/B0, CELL_CENTRE);

  if(Parallel_ion&&gyrofluid&&FLR_terms&&nonlinear)
    ddt(U) -= gynorm*B0*B0*bracket(gyroPsi_d*B0, n0*Vipar/B0, bm_mag);

  //gyroviscosity
  if(diamag&&gyroviscos){
    delp2phi = Delp2(phi);
    delp2phi.applyBoundary();
    mesh->communicate(delp2phi);
    ddt(U) += 1./(8.*gynorm*B0)*bracket(delp2phi, P0, bm_exb);

    ddt(U) += 1./(8.*gynorm*B0)*bracket(phi, delp2P0, bm_exb);

    exblinear = bracket(phi, P0, bm_exb);
    exblinear.applyBoundary();
    mesh->communicate(exblinear);
    ddt(U) -= 1./(8.*gynorm*B0)*Delp2(exblinear);

    if(diamag_phi0){
      ddt(U) += 1./(8.*gynorm*B0)*bracket(delp2phi0, P, bm_exb);

      delp2P = Delp2(P);
      delp2P.applyBoundary();
      mesh->communicate(delp2P);
      ddt(U) += 1./(8.*gynorm*B0)*bracket(phi0, delp2P, bm_exb);

      exblinear0 = bracket(phi0, P, bm_exb);
      exblinear0.applyBoundary();
      mesh->communicate(exblinear0);
      ddt(U) -= 1./(8.*gynorm*B0)*Delp2(exblinear0);
    }

    if(nonlinear){
      ddt(U) += 1./(8.*gynorm*B0)*bracket(delp2phi, P, bm_exb);

      delp2P = Delp2(P);
      delp2P.applyBoundary();
      mesh->communicate(delp2P);
      ddt(U) += 1./(8.*gynorm*B0)*bracket(phi, delp2P, bm_exb);

      exbnonlinear = bracket(phi, P, bm_exb);
      exbnonlinear.applyBoundary();
      mesh->communicate(exbnonlinear);
      ddt(U) -= 1./(8.*gynorm*B0)*Delp2(exbnonlinear);
    }
   }

  if(diffusion_u4 > 0.0){         //4th order perp diffsion, 1/sqrt(n0)
    Delp2u = Grad2_par2new(U);
    Delp2u.applyBoundary();
    mesh->communicate(Delp2u);
    ddt(U) -= diffusion_u4*Grad2_par2new(Delp2u);
   }

  // Viscosity terms
  if(viscos_par > 0.0)
    ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity

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
    if(gyrofluid)
        ddt(P) -= bracket(phi, P_G0);
    else
        ddt(P) -= bracket(phi, P0);

    if(include_pellet&&!finite_perturbation)
      ddt(P) -= bracket(phi, P_pellet);

    if(diamag_phi0)
      ddt(P) -= bracket(phi0, P, bm_exb);   // Equilibrium flow

    if(withflow)                              //net flow
      ddt(P) -= V_dot_Grad(V0net, P);

    if(nonlinear)
      ddt(P) -= bracket(phi, P, bm_exb);    // Advection

    if(gyrofluid&&FLR_terms)
      ddt(P) -= T0*bracket(gyrophi_d,n_iG0, bm_exb);

    if(gyrofluid&&diamag_phi0&&FLR_terms)
      ddt(P) -= T0*bracket(gyrophi0_d, ni1, bm_exb);

    if(gyrofluid&&nonlinear&&FLR_terms)
      ddt(P) -= T0*bracket(gyrophi_d, ni1, bm_exb);

    if(GradB_terms){
      ddt(P) += 2.*T0*J0/gynorm*b0xGrad_dot_Grad(B0, Psi*B0)/B0;
      ddt(P) -= 2.*P0/(B0*B0)*b0xGrad_dot_Grad(B0, phi);
    }

    if(GradB_terms&&FLR_terms)
      ddt(P) -= P0/B0*b0xGrad_dot_Grad(B0, gyrophi_d);

    if(Parallel_ion)
      ddt(P) -= B0*Grad_parP(P0*Vipar/B0);

    if(Parallel_ion&&gyrofluid)
      ddt(P) -= B0*Grad_parP(0.5*P0*gyroVipar_d/B0, CELL_CENTRE);

    if(Parallel_ion&&current_P)
      ddt(P) -= 1./gynorm*b0xGrad_dot_Grad(0.5*P0*J0/n0, Psi*B0)+B0*Grad_parP(0.5*T0*Jpar/gynorm, CELL_YLOW);

    if(Parallel_ion&&nonlinear&&gyrofluid&&FLR_terms)
      ddt(P) += bracket(gyroPsi_d*B0, 0.5*P0*Vipar/B0, bm_mag);

  }

  // Parallel diffusion terms
  if(diffusion_par > 0.0)
    ddt(P) += diffusion_par * Grad2_par2(P); // Parallel diffusion

  //M: 4th order Parallel diffusion terms
   if(diffusion_p4 > 0.0){
    Grad2p=Grad2_par2new(P);
    Grad2p.applyBoundary();
    mesh->communicate(Grad2p);
    ddt(P) -= diffusion_p4 * Grad2_par2new(Grad2p);
    }

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

 if(filter_weight > 0.0){
   ddt(U) = nl_filter_x(ddt(U), filter_weight);
   ddt(U) = nl_filter_y(ddt(U), filter_weight);
  }

  if(filter_z) {
    // Filter out all except filter_z_mode

    if(electron_inertial) {
      ddt(myLambda) = filter(ddt(myLambda), filter_z_mode);
    }else
      ddt(Psi) = filter(ddt(Psi), filter_z_mode);

    ddt(U) = filter(ddt(U), filter_z_mode);
    ddt(P) = filter(ddt(P), filter_z_mode);

    if(Parallel_ion)
      ddt(Lambdai) = filter(ddt(Lambdai), filter_z_mode);

  }

  if(low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
    if(electron_inertial) {
      ddt(myLambda) = lowPass(ddt(myLambda), low_pass_z, zonal_field);
    }else
      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

    ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
    ddt(P) = lowPass(ddt(P), low_pass_z, zonal_bkgd);

    if(Parallel_ion)
      ddt(Lambdai) = lowPass(ddt(Lambdai), low_pass_z, zonal_flow);

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
