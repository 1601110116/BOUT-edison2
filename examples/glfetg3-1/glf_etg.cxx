/*******************************************************************************
 * gyrofluid 3+1 model for ETG simulation 
 * see glf-etg.pdf
 * Based on P.Synder thesis
 * Developed by P.W.Xi and X.Q.Xu
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

// 2D inital profiles
Field2D P0, Pe0, n0, T0, test; // pressure, electron density, temperature
Field2D rhoe, omegae; //electron gyroradius and frequency
Vector2D b0xcv; // Curvature term
Vector2D B0vec; // B0 field vector

// 3D evolving variables
Field3D ne, Lambda, Ppar, Pperp;

// Derived 3D variables
Field3D phi, Vpar, Psi, Tperp, Jpar; // potential, electron parallel velocity, vector potential, perpendicular temperature, parallel current
Field3D gyrophi, gyrophi_a, gyrophi_b, gyroPsi, gyroPsi_a, gyrone; //gyroavered quantities

//Intermediate variables
Field3D ns, Ts, Ts1, phi_sour, phi_i;    // used in quasi-neutralility
Field3D gyroVpar;       // used in Ampere's law
Field3D gyrophi_as, gyrophi_ai, gyrophi_bs, gyrophi_bs1, gyrophi_bs2;   //used for gyrophi_a and gyrophi_b
Field3D Psi_sour, gyroPsi_as, gyroPsi_ai;    //for Psi
Field2D gyroa, gyrob, gyroc, gyrod, gyroe;   //inversion parameters


// Parameters
BoutReal density, temperature; // electron Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Vbar; // Normalisation constants
BoutReal Zion, massratio;    //ion mass number, Mi/Me
BoutReal etae, tau;   //  Ln/Lt, Ti/Te
BoutReal omegae_bar, rhoe_bar, V_alfven; 
BoutReal L_tem;  //temperature length scale

// options
bool include_curvature, include_jpar0;
bool nonlinear;
bool electromagnetic;

int general_flags; //inversion flags for any quantities

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

bool Te_fit;
BoutReal Te_0, Te_s, Te_x0, Te_C;
BoutReal Psiaxis, Psibndry;  
Field2D x, Psixy;

//fileter
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

BoutReal diffusion_ne4;   //4th order diffusion for smoothing electron density
Field3D Grad2_ne;     

bool parallel_lr_diff;
bool parallel_project;  // Use Apar to project field-lines

Field3D Xip_x, Xip_z;     // Displacement of y+1 (in cell index space)

Field3D Xim_x, Xim_z;     // Displacement of y-1 (in cell index space)

BoutReal vac_lund, core_lund;       // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist,  core_resist;  // The resistivities (just 1 / S)
Field3D eta;                    // Resistivity profile (1 / S)

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, B0, hthe;
Field2D I; // Shear factor

const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass
const BoutReal Me = 9.1094e-31;     // Electron mass

// Communication objects
FieldGroup comms;

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


const Field3D Grad2_par2new(const Field3D &f)    //4th order diffusion for smoothing parallel perturbation
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
 
  output.write("Solving gyrofluid 3+1 equations for ETG\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  //  mesh->get(n0, "Ne0");    //10^20*m^-3
  //  mesh->get(T0, "Te0"); // eV
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
  Options *options = globalOptions->getSection("glfetg");

  OPTION(options, density,            1.0e19); // Number density [m^-3]

  // Effects to include/exclude
  OPTION(options, include_curvature,   true);
  OPTION(options, include_jpar0,      false);

  
  OPTION(options, nonlinear,         false);
  OPTION(options, electromagnetic,    true);    //include electromagnetic perturbation
  OPTION(options, Zion,                1.0);    // ion mass in units of proton mass
  OPTION(options, etae,                1.0);    //Ln/Lt
  OPTION(options, tau,                 1.0);    //Ti/Te
  OPTION(options, L_tem,               1.0);    //temperature length scale

  OPTION(options, relax_j_vac,       false);       // Relax vacuum current to zero
  OPTION(options, relax_j_tconst,      0.1);

  OPTION(options, Te_fit,            false);
  OPTION(options, Te_0,                0.0); 
  OPTION(options, Te_s,                0.0);
  OPTION(options, Te_x0,               0.0);
  OPTION(options, Te_C,                0.0);
  
  // Toroidal filtering
  OPTION(options, filter_z,          false);       // Filter a single n
  OPTION(options, filter_z_mode,         1);
  OPTION(options, low_pass_z,           -1);      // Low-pass filter
  OPTION(options, zonal_flow,           -1);      // zonal flow filter
  OPTION(options, zonal_field,          -1);      // zonal field filter
  OPTION(options, zonal_bkgd,           -1);      // zonal background P filter

  // Radial smoothing
  OPTION(options, smooth_j_x,        false);  // Smooth Jpar in x
  // Jpar boundary region
  OPTION(options, jpar_bndry_width,     -1);
 
  OPTION(options, diffusion_ne4,        -1);

  // Parallel differencing
  OPTION(options, parallel_lr_diff,  false);
  OPTION(options, parallel_project,  false);

  // Vacuum region control
  OPTION(options, vacuum_pressure,   0.02);   // Fraction of peak pressure
  OPTION(options, vacuum_trans,      0.005);  // Transition width in pressure
  
  // Resistivity and hyper-resistivity options
  OPTION(options, vac_lund,          0.0);    // Lundquist number in vacuum region
  OPTION(options, core_lund,         0.0);    // Lundquist number in core region
  
  // Field inversion flags
  OPTION(options, general_flags,       0);


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

   if(!include_curvature)
    b0xcv = 0.0;
  
 
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
  //  if(mesh->get(Lbar, "L_T")) // Typical length scale
    Lbar = L_tem;
  // if(mesh->get(temperature, "Te_x")) //Typical temperature
    temperature = 0.5*Te_0; 
  
  temperature = temperature*1.602e-19;      //convert from  eV to J
  Vbar = sqrt(temperature/Me);
  Tbar = Lbar / Vbar;
  
  x=(Psixy-Psiaxis)/(Psibndry-Psiaxis);   //normalized radial coordinate     
   
  if(Te_fit)
    T0 = 0.5*Te_0*(1-tanh(Te_s*(x-Te_x0)))+Te_C;

  //definitions
  T0=T0*1.602e-19;                           //convert from  eV to J
  rhoe = Me*sqrt(T0/Me)/(1.602e-19*B0);     //2D, equilibrium
  omegae = -1.602e-19*B0/Me;                //2D
  omegae_bar = 1.602e-19*Bbar/Me;           //constant
  rhoe_bar = Me*Vbar/(1.602e-19*Bbar);      //constant
  V_alfven = Bbar*Bbar/(MU0*Me*density);    //constant

  Pe0 =  P0/(1+tau/Zion);   //electron equilibrium pressure 
  
  if(Te_fit)                //for fit electron temperature
    n0 = Pe0/T0;     
 
  gyroa = 1.0;
  gyrob = -0.5*rhoe*rhoe/(Lbar*Lbar);
  gyroc = 2.0;
  gyrod = 0.5*gyrob;
  gyroe = -(1+tau)*rhoe*rhoe/(Lbar*Lbar);
            
  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Vbar = %e m/s   Tbar = %e s\n", Vbar, Tbar);
  output.write("omegae_bar = %e rad/s, rhoe_bar = %e m, temperatrue = %e eV\n", omegae_bar, rhoe_bar, temperature/(1.602e-19)); 

  //normalize equlibrium quantities
  n0 = n0/density;    
  T0=T0/(Me*Vbar*Vbar);     
  P0 = P0/(density*Me*Vbar*Vbar);
  Pe0 = Pe0/(density*Me*Vbar*Vbar);
  omegae = -omegae*Me/(Bbar*1.602e-19);  

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

  // Transitions from 0 in core to 1 in vacuum
  vac_mask = (1.0 - tanh( (P0 - vacuum_pressure) / vacuum_trans )) / 2.0;

    // transition from 0 for large P0 to resistivity for small P0
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

  /**************** SET VARIABLE LOCATIONS *************/
  ne.setLocation(CELL_CENTRE);
  Lambda.setLocation(CELL_CENTRE);
  Ppar.setLocation(CELL_CENTRE);
  Pperp.setLocation(CELL_CENTRE);
  Vpar.setLocation(CELL_CENTRE);
 
  phi.setLocation(CELL_CENTRE);
  Psi.setLocation(CELL_YLOW);
  gyrophi.setLocation(CELL_YLOW);
  gyrophi_a.setLocation(CELL_YLOW);
  gyrophi_b.setLocation(CELL_YLOW);
  gyroPsi.setLocation(CELL_YLOW);
  gyroPsi_a.setLocation(CELL_YLOW);
  gyrone.setLocation(CELL_YLOW);
  ns.setLocation(CELL_CENTRE);
  Ts1.setLocation(CELL_YLOW);
  Ts.setLocation(CELL_CENTRE);
  Grad2_ne.setLocation(CELL_CENTRE);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(ne);
  SOLVE_FOR(Lambda);
  SOLVE_FOR(Ppar);
  SOLVE_FOR(Pperp); 
  SOLVE_FOR(Vpar);  //not real solve Vpar,only for restart 
 
  if(parallel_project) {
    // Add Xi to the dump file
    dump.add(Xip_x, "Xip_x", 1);
    dump.add(Xip_z, "Xip_z", 1);
    
    dump.add(Xim_x, "Xim_x", 1);
    dump.add(Xim_z, "Xim_z", 1);
  }

    dump.add(phi, "phi", 1);
    dump.add(Psi, "Psi", 1);
    dump.add(gyrophi, "gyrophi", 1);
    dump.add(gyroPsi, "gyroPsi", 1);
    dump.add(gyrone, "gyrone", 1);
    dump.add(gyroVpar, "gyroVpar", 1);

  // everything needed to recover physical units
      SAVE_ONCE4(n0, P0, T0, rhoe);
      SAVE_ONCE4(Lbar, Bbar, Tbar, Vbar);
      SAVE_ONCE2(omegae, B0);   
  
  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero
  
  if(!restarting) {
  
    Tperp = (Pperp-T0*ne)/n0;
    //electric field   
    //get gyroaveraged density 
    ns = invert_laplace(ne,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(ns); 
    Ts1 = invert_laplace(-4.*Tperp,general_flags,&gyroc,NULL,&gyrod);
    mesh->communicate(Ts1); 
    Ts = invert_laplace(Ts1+2.0*Tperp,general_flags,&gyroc,NULL,&gyrod);
    mesh->communicate(Ts); 
    gyrone = ns-n0*Ts/T0;
   
    //get potential      
    phi_sour = -tau*tau/(1.0+tau)*rhoe_bar*T0/(Lbar*n0)*gyrone;
    phi_i = invert_laplace(phi_sour,general_flags,&gyroa,NULL,&gyroe); 
    mesh->communicate(phi_i);
    phi = (phi_i-tau/(tau+1.0)*rhoe_bar*T0/(Lbar*n0)*gyrone)/B0;
      
    //get gyroaveraged potential
    gyrophi = invert_laplace(phi,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi);
   
    gyrophi_as = invert_laplace(2.*phi,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi_as);
    gyrophi_ai= 2.*phi-gyrophi_as;
    gyrophi_a = invert_laplace(gyrophi_ai,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi_a);

    gyrophi_bs2 = invert_laplace(0.5*phi,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi_bs2);
    gyrophi_bs1 = invert_laplace(gyrophi_bs2,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi_bs1);
    gyrophi_bs = invert_laplace(gyrophi_bs1,general_flags,&gyroa,NULL,&gyrob);
    mesh->communicate(gyrophi_bs);
    gyrophi_b = gyrophi_a + gyrophi_bs;
  
    Vpar = Lambda;
   }

  /************** SETUP COMMUNICATIONS **************/
  
  comms.add(ne);
  comms.add(Lambda);
  comms.add(Ppar);
  comms.add(Pperp);
  comms.add(Vpar);

  ns.setBoundary("ns");
  Ts1.setBoundary("Ts1");
  Ts.setBoundary("Ts");
  phi_i.setBoundary("phi_i");

  gyroVpar.setBoundary("gyroVpar");
  Psi.setBoundary("Psi");

  gyrophi.setBoundary("gyrophi");
  gyrophi_as.setBoundary("gyrophi_as");
  gyrophi_a.setBoundary("gyrophi_a");
  gyrophi_bs2.setBoundary("gyrophi_bs2");
  gyrophi_bs1.setBoundary("gyrophi_bs1");
  gyrophi_bs.setBoundary("gyrophi_bs");

  gyroPsi.setBoundary("gyroPsi");
  gyroPsi_as.setBoundary("gyroPsi_as");
  gyroPsi_a.setBoundary("gyroPsi_a");
  
  Jpar.setBoundary("Jpar");
 
  Grad2_ne.setBoundary("Grad2_ne");

  return 0;
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  if(parallel_lr_diff) {
   // Use left/right biased stencils. NOTE: First order only!
    if(loc == CELL_YLOW) {
	result = Grad_par_CtoL(f);
      }else
	result = Grad_par_LtoC(f);
    }else
      result = Grad_par(f, loc);
    
    if(nonlinear&&electromagnetic) {
      result -= bracket(Psi, f, bm_mag)*B0;
    }
  return result;
}

const Field3D Grad_parP_gyro(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  if(parallel_lr_diff) {
   // Use left/right biased stencils. NOTE: First order only!
    if(loc == CELL_YLOW) {
	result = Grad_par_CtoL(f);
      }else
	result = Grad_par_LtoC(f);
    }else
      result = Grad_par(f, loc);
    
    if(nonlinear&&electromagnetic) {
      result -= bracket(gyroPsi, f, bm_mag)*B0;
    }
  return result;
}

const Field3D Grad_parP0(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  if(parallel_lr_diff) {
   // Use left/right biased stencils. NOTE: First order only!
    if(loc == CELL_YLOW) {
	result = Grad_par_CtoL(f);
      }else
	result = Grad_par_LtoC(f);
    }else
      result = Grad_par(f, loc);
  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t)
{
  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  if(nonlinear) {
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
   }

  Tperp = (Pperp-T0*ne)/n0;

    //electric field   
    //get gyroaveraged density 
    ns = invert_laplace(ne,general_flags,&gyroa,NULL,&gyrob);
    ns.applyBoundary();
    mesh->communicate(ns); 
    Ts1 = invert_laplace(-4.*Tperp,general_flags,&gyroc,NULL,&gyrod);
    Ts1.applyBoundary();
    mesh->communicate(Ts1); 
    Ts = invert_laplace(Ts1+2.0*Tperp,general_flags,&gyroc,NULL,&gyrod);
    Ts.applyBoundary();
    mesh->communicate(Ts); 
    gyrone = ns-n0*Ts/T0;
   
    //get potential
    
    phi_sour = -tau*tau/(1.0+tau)*rhoe_bar*T0/(Lbar*n0)*gyrone;
    phi_i = invert_laplace(phi_sour,general_flags,&gyroa,NULL,&gyroe); 
    mesh->communicate(phi_i);
    phi_i.applyBoundary(); 
    phi = (phi_i-tau/(tau+1.0)*rhoe_bar*T0/(Lbar*n0)*gyrone)/B0;   

    //get gyroaveraged potential
    gyrophi = invert_laplace(phi,general_flags,&gyroa,NULL,&gyrob);
    gyrophi.applyBoundary();
    mesh->communicate(gyrophi);
   
    gyrophi_as = invert_laplace(2.*phi,general_flags,&gyroa,NULL,&gyrob);
    gyrophi_as.applyBoundary();
    mesh->communicate(gyrophi_as);
    gyrophi_ai= 2.*phi-gyrophi_as;
    gyrophi_a = invert_laplace(gyrophi_ai,general_flags,&gyroa,NULL,&gyrob);
    gyrophi_a.applyBoundary();
    mesh->communicate(gyrophi_a);

    gyrophi_bs2 = invert_laplace(0.5*phi,general_flags,&gyroa,NULL,&gyrob);
    gyrophi_bs2.applyBoundary();
    mesh->communicate(gyrophi_bs2);
    gyrophi_bs1 = invert_laplace(gyrophi_bs2,general_flags,&gyroa,NULL,&gyrob);
    gyrophi_bs1.applyBoundary();
    mesh->communicate(gyrophi_bs1);
    gyrophi_bs = invert_laplace(gyrophi_bs1,general_flags,&gyroa,NULL,&gyrob);
    gyrophi_bs.applyBoundary();
    mesh->communicate(gyrophi_bs);
    gyrophi_b = gyrophi_a + gyrophi_bs;

     //get magnetic field 
    if(electromagnetic){  
      gyroVpar = invert_laplace(Vpar,general_flags,&gyroa,NULL,&gyrob);
      gyroVpar.applyBoundary();
      mesh->communicate(gyroVpar); 

      Psi_sour = -Vbar*omegae_bar*Lbar*n0*gyroVpar/(B0*V_alfven*V_alfven);    
      Psi = invert_laplace(Psi_sour,general_flags,NULL,NULL); 
      Psi.applyBoundary();
      mesh->communicate(Psi);

      //get gyroaveraged vector potential
      gyroPsi = invert_laplace(Psi,general_flags,&gyroa,NULL,&gyrob);
      gyroPsi.applyBoundary();
      mesh->communicate(gyroPsi);
  
      gyroPsi_as = invert_laplace(2.*Psi,general_flags,&gyroa,NULL,&gyrob);
      gyroPsi_as.applyBoundary();
      mesh->communicate(gyroPsi_as);
      gyroPsi_ai= 2.*Psi-gyroPsi_as;
      gyroPsi_a = invert_laplace(gyroPsi_ai,general_flags,&gyroa,NULL,&gyrob);
      gyroPsi_a.applyBoundary();
      mesh->communicate(gyroPsi_a);
    }else 
      Psi = 0.0;

  mesh->communicate(comms);

   
  if(electromagnetic){                                         // electromagnetic get J from Psi
    Jpar = Delp2(Psi);
    Jpar.applyBoundary();
    mesh->communicate(Jpar);
  }else                                                        // electrostatic get J from Vpar
    Jpar = Lbar*Vbar*omegae_bar/(V_alfven*V_alfven)*n0*Vpar/B0;
 
 
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

  ///////////////////////////////////////////////////
  //electron guilding center density
  ddt(ne) = 0.;
  
  ddt(ne) -= B0*Grad_parP_gyro(n0*Vpar/B0, CELL_YLOW);
  ddt(ne) += b0xGrad_dot_Grad(n0, gyrophi);
  ddt(ne) += 0.5*etae*b0xGrad_dot_Grad(n0, gyrophi_a);
  ddt(ne) -= 2.*n0/omegae*b0xGrad_dot_Grad(B0, gyrophi);
  ddt(ne) -= 0.5*n0/omegae*b0xGrad_dot_Grad(B0, gyrophi_a);
  ddt(ne) += 1./(Tbar*omegae_bar*omegae*B0)*b0xGrad_dot_Grad(B0, Ppar+Pperp);

  if(nonlinear){
    ddt(ne) -= bracket(gyrophi, ne, bm_exb)*B0;
    ddt(ne) -= 0.5*n0/T0*bracket(gyrophi, Tperp, bm_exb)*B0;  
  }

  if(diffusion_ne4 > 0.0){
    Grad2_ne = Grad2_par2new(ne);
    Grad2_ne.applyBoundary();
    mesh->communicate(Grad2_ne);
    ddt(ne) -= diffusion_ne4 * Grad2_par2new(Grad2_ne);
  }
    

  ////////////////////////////////////////////////////
  //electron parallel motion 

  ddt(Lambda) = 0.;
  
  ddt(Lambda) += omegae_bar*Tbar*Grad_parP_gyro(gyrophi*B0, CELL_YLOW);
  ddt(Lambda) -= B0/n0*Grad_parP_gyro(Ppar/B0, CELL_CENTRE);
  ddt(Lambda) -= Pperp/(n0*B0)*Grad_parP0(B0, CELL_CENTRE);
  ddt(Lambda) += 0.5*omegae_bar*Tbar*gyrophi_a*Grad_parP0(B0, CELL_CENTRE);
  ddt(Lambda) += 4./(omegae_bar*Tbar*n0*omegae*B0)*b0xGrad_dot_Grad(B0, Pe0*Vpar);    

  if(electromagnetic){
   ddt(Lambda) -= T0/n0*(1.+etae)*b0xGrad_dot_Grad(n0, gyroPsi);
   ddt(Lambda) -= 0.5*T0*etae/n0*b0xGrad_dot_Grad(n0, gyroPsi_a);
  }
  
  if(core_lund>0)
    ddt(Lambda) -= omegae_bar*Tbar*B0*eta*Jpar;
 
  if(nonlinear){
    ddt(Lambda) -= bracket(gyrophi, Lambda, bm_exb)*B0;
    if(electromagnetic){
      ddt(Lambda) -= omegae_bar*Tbar*bracket(gyrophi, B0*gyroPsi, bm_exb)*B0;
      ddt(Lambda) += 0.5*bracket(gyroPsi_a, Tperp, bm_mag)*B0; 
    } 
   }
    
  Vpar = Lambda;
    //if electromagnetic, need to get Vpar
  if(electromagnetic)
    Vpar += Lbar*omegae_bar/Vbar*B0*gyroPsi;     

  ////////////////////////////////////////////////////
  //papallel pressure 
  ddt(Ppar) = 0.;

  ddt(Ppar) -= 3.*B0*Grad_parP_gyro(Pe0*Vpar/B0, CELL_CENTRE);
  ddt(Ppar) -= 2.*Pe0*Vpar/B0*Grad_parP0(B0, CELL_CENTRE);
  ddt(Ppar) += (1.+etae)*T0*b0xGrad_dot_Grad(n0, gyrophi);
  ddt(Ppar) += 0.5*etae*T0*b0xGrad_dot_Grad(n0, gyrophi_a);
  ddt(Ppar) -= 4.*n0*T0/B0*b0xGrad_dot_Grad(B0, gyrophi);
  ddt(Ppar) -= 0.5*n0*T0/B0*b0xGrad_dot_Grad(B0, gyrophi_a);
 
  if(nonlinear){
    ddt(Ppar) -= bracket(gyrophi, Ppar, bm_exb)*B0;
    ddt(Ppar) -= 0.5*n0*bracket(gyrophi_a, Tperp, bm_exb)*B0;  
  } 

  ////////////////////////////////////////////////////
  // perpendicular pressure
  ddt(Pperp) = 0.;

  ddt(Pperp) -= B0*B0*Grad_parP_gyro(Pe0*Vpar/(B0*B0), CELL_YLOW);
  ddt(Pperp) += (1.+etae)*T0*b0xGrad_dot_Grad(n0, gyrophi);
  ddt(Pperp) += (0.5+0.5*etae)*T0*b0xGrad_dot_Grad(n0, gyrophi_a);
  ddt(Pperp) += etae*T0*b0xGrad_dot_Grad(n0, gyrophi_b); 
  ddt(Pperp) -= 6.*n0*T0/B0*b0xGrad_dot_Grad(B0, gyrophi);
  ddt(Pperp) -= 3.*n0*T0/B0*b0xGrad_dot_Grad(B0, gyrophi_a);
  ddt(Pperp) -= 2.*n0*T0/B0*b0xGrad_dot_Grad(B0, gyrophi_b);

  if(nonlinear){
    ddt(Pperp) -= bracket(gyrophi, Pperp, bm_exb)*B0;
    ddt(Pperp) -= 0.5*bracket(gyrophi_a, Pperp, bm_exb)*B0;
    ddt(Pperp) += n0*bracket(gyrophi_b, Tperp, bm_exb)*B0;
   
    if(electromagnetic)
    ddt(Pperp) += 0.5*Pe0*bracket(gyroPsi_a, Vpar, bm_mag)*B0;
  }
  

  ddt(Vpar) = 0;

  ///////////////////////////////////////////////////
  //toroidal filter
  if(filter_z) {
    // Filter out all except filter_z_mode
    ddt(ne) = filter(ddt(ne), filter_z_mode);
    ddt(Lambda) = filter(ddt(Lambda), filter_z_mode);
    ddt(Ppar) = filter(ddt(Ppar), filter_z_mode);
    ddt(Pperp) = filter(ddt(Pperp), filter_z_mode);   
  }

  if(low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
    ddt(ne) = lowPass(ddt(ne), low_pass_z, zonal_field);
    ddt(Lambda) = lowPass(ddt(Lambda), low_pass_z, zonal_flow);
    ddt(Ppar) = lowPass(ddt(Ppar), low_pass_z, zonal_bkgd);
    ddt(Pperp) = lowPass(ddt(Pperp), low_pass_z, zonal_bkgd);
  }

  first_run = false;

  return 0;
}




