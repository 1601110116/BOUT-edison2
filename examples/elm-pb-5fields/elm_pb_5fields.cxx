/*******************************************************************************
 * High-Beta Flute-Reduced MHD with 5-field of (N_i, T_e, T_i, U, Psi)
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * T. Xia
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

BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x; //the total height, average width and center of profile of N0, the  start of flat region of N0
BoutReal Tconst; //the ampitude of congstant temperature

BoutReal laplace_alpha; //test the effect of first order term of invert Laplace function

// 2D inital profiles
Field2D  J0, P0; // Current and pressure
Vector2D b0xcv; // Curvature term
Field2D phi0;   // When diamagnetic terms used

Field2D N0,Ti0,Te0;  // number density and temperature

// B field vectors
Vector2D B0vec; // B0 field vector

// V0 field vectors
Vector2D V0vec; // V0 field vector in convection
Vector2D V0eff; // effective V0 field vector in Ohm's law

// 3D evolving variables
Field3D U, Psi, P;
Field3D N, Te, Ti;

// Derived 3D variables
Field3D Jpar, phi; // Parallel current, electric potential

Field3D Jpar2; //  Delp2 of Parallel current

// Constraint
Field3D C_phi;

// Parameters
BoutReal density; // Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
BoutReal Nbar, Tibar, Tebar;
BoutReal dnorm; // For diamagnetic terms: 1 / (2. * wci * Tbar)
BoutReal unorm1, unorm2;  // for normalization coefficients in vorticity equation
BoutReal jnorm;

BoutReal dia_fact; // Multiply diamagnetic term by this
BoutReal delta_i; // Normalized ion skin depth

BoutReal diffusion_par;  // Parallel pressure diffusion
BoutReal diffusion_p4;   //M: 4th Parallel pressure diffusion

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
bool evolve_pressure, continuity;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar; 
BoutReal g; // Only if compressible
bool phi_curv;

bool diamag;
bool diamag_grad_t; // Grad_par(Te) term in Psi equation
bool diamag_phi0;   // Include the diamagnetic equilibrium phi0

bool eHall;
BoutReal AA; // ion mass in units of the proton mass; AA=Mi/Mp

BoutReal Vt0; // equilibrium toroidal flow normalized to Alfven velocity
BoutReal Vp0; // equilibrium poloidal flow normalized to Alfven velocity

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



//********************


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
Field3D rmp_dApdt; // Time variation

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
const BoutReal KB = 1.38065e-23;     // Boltamann constant

// Communication objects
FieldGroup comms;

void advect_tracer(const Field3D &p,  // phi (input)
		   const Field3D &delta_x, const Field3D &delta_z, // Current location (input) 
		   Field3D &F_dx, Field3D &F_dz); // Time-derivative of location

const Field3D Grad2_par2new(const Field3D &f); //for 4th order diffusion

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x);

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

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x)
{
  Field2D result;
  result.allocate();

  BoutReal Grid_NX, Grid_NXlimit;    //the grid number on x, and the 
  mesh->get(Grid_NX, "NX");
  Grid_NXlimit = n0_bottom_x * Grid_NX;

  for(int jx=0;jx<mesh->ngx;jx++)
    {
      BoutReal mgx = mesh->GlobalX(jx);
      BoutReal xgrid_num = Grid_NXlimit/Grid_NX;
      if (mgx > xgrid_num)
	mgx = xgrid_num;
      BoutReal rlx = mgx - n0_center;
      BoutReal temp = exp(rlx/n0_width);
      BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
      for(int jy=0;jy<mesh->ngy;jy++)
	result[jx][jy] = 0.5*(1.0 - dampr) * n0_height + n0_ave;  
    }
  
  mesh->communicate(result);

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

  if(mesh->get(N0,  "Ni0")) { // T                                            
    output.write("Error: Cannot read Ni0 from grid\n");
    return 1;
  }  


  if(mesh->get(Ti0,  "Ti0")) { // T                                            
    output.write("Error: Cannot read Ti0 from grid\n");
    return 1;
  }

  if(mesh->get(Te0,  "Te0")) { // T                                            
    output.write("Error: Cannot read Te0 from grid\n");
    return 1;
  }


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

  //  options.setSection("highbeta");    // The section header for all these settings

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("highbeta");

  OPTION(options, n0_height,         0.4);   //the total height of profile of N0, in percentage of Ni_x
  OPTION(options, n0_ave,            0.01);   //the center or average of N0, in percentage of Ni_x
  OPTION(options, n0_width,          0.1);   //the width of the gradient of N0,in percentage of x
  OPTION(options, n0_center,         0.633);   //the grid number of the center of N0, in percentage of x
  OPTION(options, n0_bottom_x,       0.81);  //the start of flat region of N0 on SOL side, in percentage of x  
  OPTION(options, Tconst,            -1.0);   //the amplitude of constant temperature, in percentage

  OPTION(options, laplace_alpha,     1.0);   //test parameter for invert Lapalace

  OPTION(options, density,           1.0e19); // Number density [m^-3]

  OPTION(options, continuity,         false);  // use continuity equation

  OPTION(options, evolve_jpar,       false);  // If true, evolve J raher than Psi
  OPTION(options, phi_constraint,    false);  // Use solver constraint for phi

  // Effects to include/exclude
  OPTION(options, include_curvature, true);
  OPTION(options, include_jpar0,     true);
  OPTION(options, evolve_pressure,   true);
  
  OPTION(options, compress0,          false);
  OPTION(options, nonlinear,         false);

  OPTION(options, AA,               1.0);    // ion mass in units of proton mass

  OPTION(options, diamag,            false);  // Diamagnetic effects? 
  OPTION(options, diamag_grad_t,     diamag); // Grad_par(Te) term in Psi equation
  OPTION(options, diamag_phi0,       diamag); // Include equilibrium phi0
  OPTION(options, dia_fact,          1.0);    // Scale diamagnetic effects by this factor

  OPTION(options, Vt0,          0.0);    // equilibrium toroidal flow normalized to Alfven velocity
  OPTION(options, Vp0,          0.0);    // equilibrium poloidal flow normalized to Alfven velocity

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
  OPTION(options, diffusion_p4,        -1.0);  // M: 4th Parallel pressure diffusion


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

  if(mesh->get(Tibar, "Ti_x")) // Typical ion temperature scale
    Tibar = 1.0;

  if(mesh->get(Tebar, "Te_x")) // Typical electron temperature scale
    Tebar = 1.0;

    Nbar = 1.e20/density;
  
  Va = sqrt(Bbar*Bbar / (MU0*Mi*density));

  Tbar = Lbar / Va;

  dnorm = dia_fact * KB*Tbar*Tibar*11604 / (1.602e-19*Lbar*Lbar*Bbar);

  unorm1 = 1.0;
  unorm2 = KB*(Tibar+Tebar)/2.*11604.*Tbar*Tbar/(Mi*Lbar*Lbar);

  delta_i = AA*60.67*5.31e5/sqrt(density/1e6)/(Lbar*100.0);

  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Va = %e m/s   Tbar = %e s\n", Va, Tbar);
  output.write("                Nbar = %e * %e m^-3\n",Nbar,density);
  output.write("                Tibar = %e eV   Tebar = %e eV\n", Tibar, Tebar);
  output.write("                dnorm = %e\n", dnorm);
  output.write("                unorm1 = %e    unorm2 = %e\n", unorm1,unorm2);
  output.write("    Resistivity\n");


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


  J0 = MU0*Lbar * J0 / B0;
  P0 = P0/(KB * (Tibar+Tebar)/2. * density * 11604);

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
  if(Tconst<0.0 && n0_ave>0.0)
    {
      N0 = N0tanh(n0_height*Nbar, n0_ave*Nbar, n0_width, n0_center, n0_bottom_x);

      Ti0 = P0/N0/2.0;
      Te0 = Ti0;
    }
  else if (Tconst>0.0)
    {
      Ti0 = Tconst;
      Te0 = Ti0;
      N0 = P0/(Ti0+Te0);
    }
  else
    {
      if(mesh->get(N0,  "Ni0")) { // N_i0 
	output.write("Error: Cannot read Ni0 from grid\n");
	return 1;
      }  
      
      if(mesh->get(Ti0,  "Ti0")) { // T_i0                                        
	output.write("Error: Cannot read Ti0 from grid\n");
	return 1;
      }

      if(mesh->get(Te0,  "Te0")) { // T_e0
	output.write("Error: Cannot read Te0 from grid\n");
	return 1;
      }
      N0 /= density*Nbar;
      Ti0 /= Tibar;
      Te0 /= Tebar;
    }


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

  // Set V0vec field vector
  
  V0vec.covariant = false;
  V0vec.x = 0.;
  V0vec.y = Vp0 / hthe;
  V0vec.z = Vt0 / Rxy;

  // Set V0eff field vector

  V0eff.covariant = false;
  V0eff.x = 0.;
  V0eff.y = -(Btxy/(B0*B0))*(Vp0*Btxy-Vt0*Bpxy) / hthe;
  V0eff.z =  (Bpxy/(B0*B0))*(Vp0*Btxy-Vt0*Bpxy) / Rxy;

  /**************** SET VARIABLE LOCATIONS *************/

  P.setLocation(CELL_CENTRE);
  U.setLocation(CELL_CENTRE);
  phi.setLocation(CELL_CENTRE);
  Psi.setLocation(CELL_YLOW);
  Jpar.setLocation(CELL_YLOW);

  N.setLocation(CELL_CENTRE);
  Ti.setLocation(CELL_CENTRE);
  Te.setLocation(CELL_CENTRE);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(U);
  SOLVE_FOR(N);
  SOLVE_FOR(Ti);
  SOLVE_FOR(Te);


  output.write("Solving for Psi, Differentiating to get jpar\n");
  SOLVE_FOR(Psi);
  dump.add(Jpar, "jpar", 1);

  dump.add(P, "P", 1);
  
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

  if(phi_constraint) {
    // Implicit Phi solve using IDA
    
    if(!bout_constrain(phi, C_phi, "phi")) {
      output.write("ERROR: Cannot constrain. Run again with phi_constraint=false\n");
      bout_error("Aborting.\n");
    }
    
  }else {
    // Phi solved in RHS (explicitly)
    dump.add(phi, "phi", 1);
  }

  // Diamagnetic phi0
  if(diamag && diamag_phi0) {
    // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
    phi0 = -0.5*dnorm*P0/B0/N0;
    SAVE_ONCE(phi0);
  }

  // Add some equilibrium quantities and normalisations
  // everything needed to recover physical units
  SAVE_ONCE2(J0, P0);
  SAVE_ONCE4(Nbar, Lbar, Bbar, Tbar);
  SAVE_ONCE2(Va, B0);
  SAVE_ONCE3(Ti0, Te0, N0);

  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero
  
  if(!restarting) {
    // Only if not restarting: Check initial perturbation

    // Set U to zero where P0 < vacuum_pressure
    U = where(P0 - vacuum_pressure, U, 0.0);

    P = N*(Ti0 + Te0) + N0*(Ti + Te) +  N*(Ti + Te);

    //    Field2D lap_temp = 0.0;
    Field2D logn0 = laplace_alpha * N0;
    Field3D ubyn = U*B0/N0;
    // Phi should be consistent with U
    if(laplace_alpha <= 0.0)
      phi = invert_laplace(ubyn, phi_flags, NULL)/B0;
    else
      phi = invert_laplace(ubyn, phi_flags, NULL, &logn0, NULL)/B0;
    //phi = invert_laplace(U, phi_flags, NULL);    

    if(diamag) {
      phi -= Ti0/(Ti0+Te0)*dnorm * P / B0;
    }
  }

  /************** SETUP COMMUNICATIONS **************/
  
  comms.add(U);
  comms.add(Psi);
  comms.add(phi);
  comms.add(N);
  comms.add(Ti);
  comms.add(Te);

  phi.setBoundary("phi"); // Set boundary conditions

  P.setBoundary("P");
  Jpar.setBoundary("J");
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
      //      result -= b0xGrad_dot_Grad(Psi + rmp_Psi, f);
      result -= b0xGrad_dot_Grad(Psi, f);
    }
  }

  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t)
{
  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  //  if(nonlinear) {
  //    vac_mask = (1.0 - tanh( ((P0 + P) - vacuum_pressure) / vacuum_trans )) / 2.0;
  
  ////////////////////////////////////////////
  // Inversion

  P = N*(Ti0 + Te0) + N0*(Ti + Te) +  N*(Ti + Te);
  P.applyBoundary();
  mesh->communicate(P);

  Field2D logn0 = laplace_alpha * N0;
  Field3D ubyn = U*B0/N0;

  // Invert laplacian for phi
  if(laplace_alpha <= 0.0)
    phi = invert_laplace(ubyn, phi_flags, NULL)/B0;
  else
    phi = invert_laplace(ubyn, phi_flags, NULL, &logn0, NULL)/B0;
  
  if(diamag) {
    phi -= Ti0/(Ti0+Te0)*dnorm * P / B0;
  }
    
    // Apply a boundary condition on phi for target plates
  phi.applyBoundary();
    //  }

  // Perform communications
  mesh->communicate(comms);


  Jpar = -Delp2(Psi);
  
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
  Jpar2 = -Delp2(Jpar);

  Jpar2.applyBoundary();
  mesh->communicate(Jpar2);

  if(jpar_bndry_width > 0) 
    {
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




  ////////////////////////////////////////////////////
  // Parallel electric field

      
  ddt(Psi) = -Grad_parP(B0*phi, CELL_CENTRE) / B0 - eta*Jpar;

  if(diamag && diamag_phi0)
    ddt(Psi) -= b0xGrad_dot_Grad(B0*phi0, Psi)/B0;   // Equilibrium flow


    // Hyper-resistivity
  if(hyperresist > 0.0) {
    ddt(Psi) += eta*hyperresist * Delp2(Jpar);
  }

  ////////////////////////////////////////////////////
  // Vorticity equation

  ddt(U) = -(B0^2) * b0xGrad_dot_Grad(Psi, J0, CELL_CENTRE); // Grad j term

 
  ddt(U) += 2.0* unorm2 * b0xcv*Grad(P);  // curvature term


  ddt(U) += (B0^2)*Grad_parP(Jpar, CELL_CENTRE); // b dot grad j

  
  if(diamag && diamag_phi0)
    ddt(U) -= b0xGrad_dot_Grad(phi0, U);   // Equilibrium flow


  if(nonlinear) {
    ddt(U) -= b0xGrad_dot_Grad(phi, U);    // Advection
  }

  // Viscosity terms 
  if(viscos_par > 0.0)
    ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
  
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



  ///////////////////////////////////////////////
  // number density equation

  ddt(N) = 0.0;

  if (continuity)
    {
      if (diamag && diamag_phi0)
	ddt(N) -= 2.0*N0/B0 *  b0xcv*Grad(phi0*B0);
      ddt(N) -= 2.0*N0/B0 *  b0xcv*Grad(phi*B0);
    }
  ddt(N) -=  b0xGrad_dot_Grad(phi, N0);

  if(diamag && diamag_phi0)
    ddt(N) -= b0xGrad_dot_Grad(phi0, N);   // Equilibrium flow

  if(nonlinear) {
    ddt(N) -= b0xGrad_dot_Grad(phi, N);    // Advection                         
  }


  ///////////////////////////////////////////////                                // ion temperature equation                                                    

  ddt(Ti) = 0.0;

  if (continuity)
    {
      if (diamag && diamag_phi0)
	ddt(Ti) += 2.0*Ti0/B0 *  b0xcv*Grad(phi0*B0);
      ddt(Ti) += 2.0*Ti0/B0 *  b0xcv*Grad(phi*B0);
    }

  ddt(Ti) -= b0xGrad_dot_Grad(phi, Ti0);

  if(diamag && diamag_phi0)
    ddt(Ti) -= b0xGrad_dot_Grad(phi0, Ti);   // Equilibrium flow

  if(nonlinear) {
    ddt(Ti) -= b0xGrad_dot_Grad(phi, Ti);    // Advection
  }


  ///////////////////////////////////////////////                                // electron temperature equation                                                    
  ddt(Te) = 0.0;

  if (continuity)
    {
      if (diamag && diamag_phi0)
	ddt(Te) += 2.0*Ti0/B0 *  b0xcv*Grad(phi0*B0);
      ddt(Te) += 2.0*Ti0/B0 *  b0xcv*Grad(phi*B0);
    }

  ddt(Te) -= b0xGrad_dot_Grad(phi, Te0);

  if(diamag && diamag_phi0)
    ddt(Te) -= b0xGrad_dot_Grad(phi0, Te);   // Equilibrium flow

  if(nonlinear) {
    ddt(Te) -= b0xGrad_dot_Grad(phi, Te);    // Advection
  }

  if(filter_z) {
    // Filter out all except filter_z_mode
    
 
    ddt(Psi) = filter(ddt(Psi), filter_z_mode);
    
    ddt(U) = filter(ddt(U), filter_z_mode);
    //    ddt(P) = filter(ddt(P), filter_z_mode);

    ddt(N) = filter(ddt(N), filter_z_mode);
    ddt(Ti) = filter(ddt(Ti), filter_z_mode);
    ddt(Te) = filter(ddt(Te), filter_z_mode);

  }

  if(low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

      ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
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


