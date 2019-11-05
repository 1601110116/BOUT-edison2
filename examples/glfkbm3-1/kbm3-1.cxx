/*******************************************************************************
 * gyrofluid 3+1 model for KBM simulation 
 * see KBM3-1.pdf
 * Based on P.Synder thesis
 * Developed by C.H.Ma, P.W.Xi and X.Q.Xu
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
// #include <integrops.hxx>
// #include <inttorops.hxx>

// 2D inital profiles
Field2D P0, n0, T0, Ti0, Te0, nG0, J0, U0, beta; // pressure, electron density, temperature
Field2D Pi0, Pe0;
Field2D rhoi, Vthermal; //ion gyroradius, thermal velocity
Vector2D b0xcv; // Curvature term
Vector2D B0vec; // B0 field vector
BoutReal Rmajor, P00;
int timestep = 0;
string path;
bool fakerun;

// 3D evolving variables
Field3D ni, Lambda_i, Ppar_i, Pperp_i, U, Psi, Ppar_e, Pperp_e,P;
Field3D Ppar, Pperp, Pi, Pe;

// Derived 3D variables
Field3D Tperp_i, Tpar_i, ne, Tperp_e, Tpar_e, Te;  
Field3D Ni_tmp, Ne_tmp;
Field3D gyroni, phi, phi_f;
Field3D gyrophi, gyrophi2, gyrophi3, gyrophi_a, gyrophi_b, gyroPsi, gyroPsi_a;
Field3D Jpar, Vpar_i, gyroVpar_i, Vpar_e;

//Intermediate variables
Field3D ns, Ts, Ts1, nmid, gyrosour, gyropar, sour3;    // used in quasi-neutralility
Vector2D sour1;
Vector3D sour2;
Field3D gyrophi_as, gyrophi_ai, gyrophi_bi, gyrophi_bs1, gyrophi_bs2, gyroPsi_as, gyroPsi_ai;   //used for gyrophi_a, gyrophi_b, gyroPsi_a

// Parameters
BoutReal density; //normalize density [m^-3] 
BoutReal Bbar, Lbar, Tbar, Vbar; // Normalisation constants
BoutReal Va, T0_max;   //  Ln/Lt, Alfven velocity
BoutReal omega_bar;  //e*Bbar/Mi 
BoutReal Cnor; //Tbar*omega_bar

// options
bool include_curvature;
int curv_model;
bool nonlinear;
int nonlinear_terms;
bool shear_Alfven_wave, shear_Alfven_advection, eHall, eHall_adiabatic, quasineutral, electrostatic;
BoutReal beta_coeff, kpar, Low_limit, Low_limit_n0;
int electrostatic_poisson_model;
bool Zero_Te, GLF_terms, solve_lambda, gyroaverage, FLR_effect, evolve_ne, evolve_Te;
bool compression, continuity, energy_flux, isotropic;

//Landau damping
bool Landau_damping_i, Landau_damping_e;
Field3D Qperp_i, Qpar_i, Qperp_e, Qpar_e;    //perpendicular and parallel heat flux
BoutReal Lpar, k0_landau;

//toroidal closure
bool toroidal_closure1,toroidal_closure2,toroidal_closure3;  //include toroidal closure in pressure equation
Field2D vgbxhat, vgbyhat, vgbzhat, modb0xglB0, modvgb; // grad B velocity components
Vector2D b0xglB0, vgb,vgbhat, vcurv; // grad-B drift
Vector2D bxgradp, bxkappa;
BoutReal k0_toroidal, zperiod;
Field3D toroidal1,toroidal2,toroidal3,toroidal4,toroidal5;
Field3D mod_wd_Vpar,mod_wd_Tpar,mod_wd_Tperp;
BoutReal nu1r, nu1i, nu2r, nu2i, nu3r, nu3i, nu4r, nu4i, nu5r, nu5i;

//inversion 
int ns_flags, phi_flags, Psi_flags, Vpar_flags; //inversion flags for any quantities
Field2D gyroa, gyrob, gyroc, gyrod, gyroe, gyrop;   //inversion parameters
Field2D phia, phid;   //inversion parameters for electrostatic phi

/***************** profile control *********************/

int Equilibrium_case;

BoutReal P_0,s_p,P_min,PDped,Pdel,P_x0,Pbottom,Pa, T_0, s_t,etai;
BoutReal N_0, N_d, N_w, N_m, N_a, N_b, N_s;
BoutReal T_d, T_w, T_m, T_a, T_b, T_s;

BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x, n0_topgrad, n0_bottomgrad; //the total height, average width and center of profile of N0
BoutReal T0_const, const_Te0;

BoutReal n0_cyclone; 

bool fit_pressure, phi_constant_density, constant_density, n0_fake_prof;
bool fit_all, fit_density, fit_temperature;

/******************************************************/

//filter
bool filter_z;
int filter_z_mode;
int low_pass_z;
int zonal_flow;
int zonal_field;
int zonal_bkgd;

BoutReal nl_y;

//resistivity
BoutReal vac_lund, core_lund, hyperresist;       // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist,  core_resist;  // The resistivities (just 1 / S)
Field3D eta;                    // Resistivity profile (1 / S)

bool smooth_j_x, smooth_U0_x;  // Smooth Jpar in the x direction
int jpar_bndry_width; // Zero jpar in a boundary region
BoutReal diffusion4_Ppar_e;   //4th order diffusion for smoothing electron density
Field3D Grad2_Ppar_e;     
BoutReal diffusion_Ppar4;   //4th order diffusion for smoothing electron density
Field3D Grad2_Ppar;   

BoutReal viscos_par, viscos_Ppar_e;  //parallel viscosity in vorticity
BoutReal viscos_perp; // Perpendicular viscosity
BoutReal hyperviscos; // Hyper-viscosity (radial)
Field3D hyper_mu_x; // Hyper-viscosity coefficient
BoutReal n0_max;    

BoutReal sink_Ul;     // left edge sink in vorticity
BoutReal su_widthl;   // left edge sink profile radial width in vorticity
BoutReal su_lengthl;  // left edge sink radial domain in vorticity

BoutReal sink_Ur;     // right edge sink in vorticity
BoutReal su_widthr;   // right edge sink profile radial width in vorticity
BoutReal su_lengthr;  // right edge sink radial domain in vorticity

BoutReal sink_nil;     // left edge sink in density
BoutReal sni_widthl;   // left edge sink profile radial width in density
BoutReal sni_lengthl;  // left edge sink radial domain in vorticity

BoutReal sink_nir;     // right edge sink in density
BoutReal sni_widthr;   // right edge sink profile radial width in density
BoutReal sni_lengthr;  // right edge sink radial domain in density

BoutReal sink_Ppar_il;     // left edge sink in pressure
BoutReal sPpar_i_widthl;   // left edge sink profile radial width in pressure
BoutReal sPpar_i_lengthl;  // left edge sink radial domain in pressure

BoutReal sink_Ppar_ir;     // right edge sink in pressure
BoutReal sPpar_i_widthr;   // right edge sink profile radial width in pressure
BoutReal sPpar_i_lengthr;  // right edge sink radial domain in pressure

BoutReal sink_Pperp_il;     // left edge sink in pressure
BoutReal sPperp_i_widthl;   // left edge sink profile radial width in pressure
BoutReal sPperp_i_lengthl;  // left edge sink radial domain in pressure

BoutReal sink_Pperp_ir;     // right edge sink in pressure
BoutReal sPperp_i_widthr;   // right edge sink profile radial width in pressure
BoutReal sPperp_i_lengthr;  // right edge sink radial domain in pressure

BoutReal sink_Ppar_el;     // left edge sink in pressure
BoutReal sPpar_e_widthl;   // left edge sink profile radial width in pressure
BoutReal sPpar_e_lengthl;  // left edge sink radial domain in pressure

BoutReal sink_Ppar_er;     // right edge sink in pressure
BoutReal sPpar_e_widthr;   // right edge sink profile radial width in pressure
BoutReal sPpar_e_lengthr;  // right edge sink radial domain in pressure

BoutReal sink_Pperp_el;     // left edge sink in pressure
BoutReal sPperp_e_widthl;   // left edge sink profile radial width in pressure
BoutReal sPperp_e_lengthl;  // left edge sink radial domain in pressure

BoutReal sink_Pperp_er;     // right edge sink in pressure
BoutReal sPperp_e_widthr;   // right edge sink profile radial width in pressure
BoutReal sPperp_e_lengthr;  // right edge sink radial domain in pressure

BoutReal sink_Lambda_ir;     // right edge sink in pressure
BoutReal sLambda_i_widthr;   // right edge sink profile radial width in pressure
BoutReal sLambda_i_lengthr;  // right edge sink radial domain in pressure

bool parallel_lr_diff;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

/************* for code testing *********************/
bool code_test;
Field3D mytest,mytest3,mytesta,mytest3a;
Vector2D mytest2;

bool ddtU_Terms_test;
Field3D ddtU_phiU0, ddtU_phi_fnG0, ddtU_GradJpar, ddtU_PsiJ0, ddtU_GradP, ddtU_Vpar, ddtU_B0phiT0;
bool ddtPsi_Terms_test;
Field3D ddtPsi_Gradphi, ddtPsi_Jpar, ddtPsi_Psin0, ddtPsi_PsiT0, ddtPsi_eHall, ddtPsi_hyperresist;
bool GradparJ_test;
Field3D GradparJ_C, GradparJ_U;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, B0, hthe;
Field2D I; // Shear factor
BoutReal Psiaxis, Psibndry;  
Field2D x, Psixy;

const BoutReal PI = 3.14159265;
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


const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x, BoutReal n0_topgrad = 0, BoutReal n0_bottomgrad = 0);
const Field3D Curv(const Field3D &f);
const Field3D Curv(const Field2D &g, const Field3D &f);

const Field3D field_larger(const Field3D &f, const BoutReal limit)
{
  Field3D result;
  result.allocate();

//  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
      {
        if(f[jx][jy][jz] >= limit)
	  result[jx][jy][jz] = f[jx][jy][jz];
	else
	  result[jx][jy][jz] = limit;
      }
  mesh->communicate(result);
  return(result);
}

const Field2D field_larger(const Field2D &f, const BoutReal limit)
{
  Field2D result;
  result.allocate();

//  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
     for(int jy=0;jy<mesh->ngy;jy++)
     {
        if(f[jx][jy] >= limit)
           result[jx][jy] = f[jx][jy];
        else
           result[jx][jy] = limit;
     }
  mesh->communicate(result);
  return(result);
}

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

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x, BoutReal n0_topgrad, BoutReal n0_bottomgrad)
{
   Field2D result;
   result.allocate();

   BoutReal Grid_NX, Grid_NXlimit;    //the grid number on x, and the 
   BoutReal Jysep;
   mesh->get(Grid_NX, "nx");
   mesh->get(Jysep, "jyseps1_1");
   Grid_NXlimit = n0_bottom_x * Grid_NX; 
   output.write("Jysep1_1 = %i   Grid number = %e\n", int(Jysep), Grid_NX);

   if (Jysep > 0.) //for single null geometry
   {
      BoutReal Jxsep, Jysep2;
      mesh->get(Jxsep, "ixseps1");
      mesh->get(Jysep2, "jyseps2_2");
      //output.write("Jysep2_2 = %i   Ixsep1 = %i\n", int(Jysep2), int(Jxsep));

      for(int jx=0;jx<mesh->ngx;jx++)
      {
         BoutReal mgx = mesh->GlobalX(jx);
         BoutReal xgrid_num = (Jxsep+1.)/Grid_NX;
         //output.write("mgx = %e xgrid_num = %e\n", mgx);
         for (int jy=0;jy<mesh->ngy;jy++)
         {
            int globaly = mesh->YGLOBAL(jy);
            //output.write("local y = %i;   global y: %i\n", jy, globaly);
            if ( mgx > xgrid_num || (globaly<=int(Jysep)-4) || (globaly>int(Jysep2)) )
               mgx = xgrid_num;
            BoutReal rlx = mgx - n0_center;
            BoutReal temp = exp(rlx/n0_width);
            BoutReal dampr = (((1. + n0_bottomgrad * rlx) * temp - (1. - n0_topgrad * rlx) / temp) / (temp + 1.0 / temp));
            result[jx][jy] = 0.5*(1.0 - dampr) * n0_height + n0_ave;  
         }
      }
   }
   else //circular geometry
   {
      for(int jx=0;jx<mesh->ngx;jx++)
      {
         BoutReal mgx = mesh->GlobalX(jx);
         BoutReal xgrid_num = Grid_NXlimit/Grid_NX;
         if (mgx > xgrid_num)
            mgx = xgrid_num;
         BoutReal rlx = mgx - n0_center;
         BoutReal temp = exp(rlx/n0_width);
         BoutReal dampr = (((1. + n0_bottomgrad * rlx) * temp - (1. - n0_topgrad * rlx) / temp) / (temp + 1.0 / temp));
         for(int jy=0;jy<mesh->ngy;jy++)
            result[jx][jy] = 0.5*(1.0 - dampr) * n0_height + n0_ave;  
      }
   }

   mesh->communicate(result);

   return result;
}

int physics_init(bool restarting)
{

   output.write("Solving gyrofluid 3+1 equations for KBM\n");
   output.write("\tFile    : %s\n", __FILE__);
   output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

   //////////////////////////////////////////////////////////////
   // Load data from the grid

   // Load 2D profiles
   mesh->get(J0, "Jpar0"); 
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

   // Btxy = smooth_x(smooth_x(smooth_x(smooth_x(smooth_x(Btxy)))));
   //  hthe = smooth_x(smooth_x(smooth_x(smooth_x(smooth_x(hthe)))));
   // B0 = smooth_x(smooth_x(smooth_x(smooth_x(smooth_x(B0)))));
   // Bpxy = smooth_x(smooth_x(smooth_x(smooth_x(smooth_x(Bpxy)))));

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
   Options *options = globalOptions->getSection("glfkbm");
   globalOptions->get("zperiod", zperiod, 1);

   OPTION(options, density,          1.0e19); // Number density [m^-3]

   // Effects to include/exclude
   OPTION(options, include_curvature,   true);
   OPTION(options, curv_model,             1);
   OPTION(options, nonlinear,          false);
   OPTION(options, nonlinear_terms,      100);
   OPTION(options, shear_Alfven_wave,   true);
   OPTION(options, shear_Alfven_advection, false);
   OPTION(options, eHall,             false);
   OPTION(options, eHall_adiabatic,   false);
   OPTION(options, beta_coeff,           1.);
   OPTION(options, kpar,                -1.);
   OPTION(options, Zero_Te,            false);
   OPTION(options, const_Te0,          -2000);
   OPTION(options, quasineutral,        true);
   OPTION(options, phi_constant_density, false);
   OPTION(options, electrostatic,       false);
   OPTION(options, electrostatic_poisson_model, 0);
   OPTION(options, GLF_terms,           true);
   OPTION(options, solve_lambda,        false);
   OPTION(options, gyroaverage,           true);   // Using gyroaveraged phi
   OPTION(options, FLR_effect,             false);   // Include FLR effect (phi_f terms)
   OPTION(options, continuity,             false);   // use continuity equation 
   OPTION(options, compression,            false);
   OPTION(options, energy_flux,            false);
   OPTION(options, isotropic, false);
   OPTION(options, Low_limit, 1.e-10);
   OPTION(options, Low_limit_n0, 1.e-2);

   OPTION(options, evolve_ne,              false);
   OPTION(options, evolve_Te,              false);

   OPTION(options, fakerun,                false);
   OPTION(options, path,                      "./");

   //Landau damping
   OPTION(options, Landau_damping_i,  false);    //Landau damping
   OPTION(options, Landau_damping_e,  false);    //Landau damping
   OPTION(options, Lpar,                1.0);

   //toroidal colsure
   OPTION(options, toroidal_closure1, false);    //toroidal closure in pressure equation
   OPTION(options, toroidal_closure2, false);
   OPTION(options, toroidal_closure3, false);
   OPTION(options, nu1r,                0.0);
   OPTION(options, nu1i,                0.0);
   OPTION(options, nu2r,                0.0);
   OPTION(options, nu2i,                0.0);
   OPTION(options, nu3r,                0.0);
   OPTION(options, nu3i,                0.0);
   OPTION(options, nu4r,                0.0);
   OPTION(options, nu4i,                0.0);
   OPTION(options, nu5r,                0.0);
   OPTION(options, nu5i,                0.0);

   /*********************** profile control ***************/
   OPTION(options, Equilibrium_case,      -1);

   ///////// case 1:

   OPTION(options, fit_all,          false);
   OPTION(options, fit_pressure,     false);
   OPTION(options, fit_density,      false);
   OPTION(options, fit_temperature,  false);
   OPTION(options, P_0,                 8675);
   OPTION(options, s_p,                0.225);
   OPTION(options, P_min,               1.14);
   OPTION(options, P_x0,                0.92);
   OPTION(options, PDped,              0.075);
   OPTION(options, Pa,                0.0028);
   OPTION(options, Pdel,                22.0); 
   OPTION(options, Pbottom,            -21.5);
   OPTION(options, etai,                 1.0);    //Ln/Lt
   OPTION(options, N_0,        1.13905322872);
   OPTION(options, N_d,       0.398383003076);
   OPTION(options, N_w,       0.654934553111);
   OPTION(options, N_m,       0.124230472429);
   OPTION(options, N_a,       0.3);
   OPTION(options, N_b,       0.2);
   OPTION(options, N_s,       0.34);
   OPTION(options, T_0,        1.13905322872);
   OPTION(options, T_d,       0.398383003076);
   OPTION(options, T_w,       0.654934553111);
   OPTION(options, T_m,       0.124230472429);
   OPTION(options, T_a,       0.3);
   OPTION(options, T_b,       0.2);
   OPTION(options, T_s,       0.34);

   ///////// case 2:

   OPTION(options, n0_cyclone,        1.0e19);

   ///////// case 4:

   OPTION(options, n0_fake_prof,     true);
   OPTION(options, n0_height,         0.4);   //the total height of profile of N0, in percentage of Ni_x
   OPTION(options, n0_ave,           0.01);   //the center or average of N0, in percentage of Ni_x
   OPTION(options, n0_width,          0.1);   //the width of the gradient of N0,in percentage of x
   OPTION(options, n0_center,       0.633);   //the grid number of the center of N0, in percentage of x
   OPTION(options, n0_bottom_x,      0.81);  //the start of flat region of N0 on SOL side, in percentage of x 
   OPTION(options, n0_topgrad,       0.0);   //the gradient on the top
   OPTION(options, n0_bottomgrad,     0.0);   //the gradient on the bottom
   OPTION(options, T0_const,         1000);   // Constant T0 in eV

   ///////// case 5:
   
   OPTION(options, P00,           23072.3);        

   /********************************************************/  

   // Toroidal filtering
   OPTION(options, filter_z,          false);       // Filter a single n
   OPTION(options, filter_z_mode,         1);
   OPTION(options, low_pass_z,           -1);      // Low-pass filter
   OPTION(options, zonal_flow,           -1);      // zonal flow filter
   OPTION(options, zonal_field,          -1);      // zonal field filter
   OPTION(options, zonal_bkgd,           -1);      // zonal background P filter
   OPTION(options, nl_y,                 -1);      // nonlinear filter in y 

   // Radial smoothing
   OPTION(options, smooth_j_x,        false);  // Smooth Jpar in x 
   OPTION(options, smooth_U0_x,         false);
   OPTION(options, jpar_bndry_width,     -1);  // Jpar boundary region

   // Resistivity and hyper-resistivity options
   OPTION(options, vac_lund,          0.0);    // Lundquist number in vacuum region
   OPTION(options, core_lund,         0.0);    // Lundquist number in core region   
   OPTION(options, hyperresist,       0.0);   //hyper-resistivity

   OPTION(options, diffusion4_Ppar_e,    -1);
   OPTION(options, diffusion_Ppar4,      -1);

   OPTION(options, viscos_par,        -1.0); 
   OPTION(options, viscos_perp,       -1.0);  // Perpendicular viscosity
   OPTION(options, hyperviscos,       -1.0);  // Radial hyperviscosity
   OPTION(options, viscos_Ppar_e,     -1.0); 

   // left edge sink factor in vorticity
   OPTION(options, sink_Ul,           -1.0);  //  left edge sink in vorticity
   OPTION(options, su_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in vorticity
   OPTION(options, su_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in vorticity

   // right edge sink factor in vorticity
   OPTION(options, sink_Ur,           -1.0);  //  right edge sink in vorticity
   OPTION(options, su_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in vorticity
   OPTION(options, su_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in vorticity

   // left edge sink factor in density
   OPTION(options, sink_nil,           -1.0);  //  left edge sink in density
   OPTION(options, sni_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in density
   OPTION(options, sni_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in density

   // right edge sink factor in density
   OPTION(options, sink_nir,           -1.0);  //  right edge sink in density
   OPTION(options, sni_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in density
   OPTION(options, sni_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in density

   // left edge sink factor in pressure
   OPTION(options, sink_Ppar_il,           -1.0);  //  left edge sink in pressure
   OPTION(options, sPpar_i_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPpar_i_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in pressure

   // right edge sink factor in pressure
   OPTION(options, sink_Ppar_ir,           -1.0);  //  right edge sink in pressure
   OPTION(options, sPpar_i_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPpar_i_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in pressure

   // left edge sink factor in pressure
   OPTION(options, sink_Pperp_il,           -1.0);  //  left edge sink in pressure
   OPTION(options, sPperp_i_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPperp_i_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in pressure

   // right edge sink factor in pressure
   OPTION(options, sink_Pperp_ir,           -1.0);  //  right edge sink in pressure
   OPTION(options, sPperp_i_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPperp_i_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in pressure

   // left edge sink factor in pressure
   OPTION(options, sink_Ppar_el,           -1.0);  //  left edge sink in pressure
   OPTION(options, sPpar_e_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPpar_e_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in pressure

   // right edge sink factor in pressure
   OPTION(options, sink_Ppar_er,           -1.0);  //  right edge sink in pressure
   OPTION(options, sPpar_e_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPpar_e_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in pressure

   // left edge sink factor in pressure
   OPTION(options, sink_Pperp_el,           -1.0);  //  left edge sink in pressure
   OPTION(options, sPperp_e_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPperp_e_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in pressure

   // right edge sink factor in pressure
   OPTION(options, sink_Pperp_er,           -1.0);  //  right edge sink in pressure
   OPTION(options, sPperp_e_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in pressure
   OPTION(options, sPperp_e_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in pressure

   OPTION(options, sink_Lambda_ir,           -1.0);  //  right edge sink in pressure
   OPTION(options, sLambda_i_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in pressure
   OPTION(options, sLambda_i_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in pressure


   // Vacuum region control
   OPTION(options, vacuum_pressure,   0.02);   // Fraction of peak pressure
   OPTION(options, vacuum_trans,      0.005);  // Transition width in pressure

   OPTION(options, code_test,        false);
   OPTION(options, ddtU_Terms_test,  false);
   OPTION(options, ddtPsi_Terms_test,false);
   OPTION(options, GradparJ_test,    false);

   // Field inversion flags
   OPTION(options, Psi_flags,           0);
   OPTION(options, phi_flags,           0);
   OPTION(options, ns_flags,            0);
   OPTION(options, Vpar_flags,          0);

   OPTION(options, parallel_lr_diff,  false);
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
   // NORMALISE QUANTITIE  
   x=(Psixy-Psiaxis)/(Psibndry-Psiaxis);   //normalized radial coordinate     
   for(int jx=0; jx<mesh->ngx; jx++)
      for(int jy=0; jy<mesh->ngy; jy++)
         if(x[jx][jy] > 1.)
            x[jx][jy] = 1.;
   mesh->communicate(x);
   SAVE_ONCE(x);

   /******************************* profile control ************************************/
   constant_density = false;

   switch(Equilibrium_case){
      case 1: 
         P0=P_0*exp(s_p*(P_min*(exp(-(P_x0-x)/PDped)+(1.0+Pa*(P_x0-x)/PDped)*exp((P_x0-x)/PDped)*Pdel)/(exp(-(P_x0-x)/PDped)+exp((P_x0-x)/PDped))+Pbottom));
         if(etai>0)
            s_t = etai*s_p/(1+etai); 
         if(etai>100){ //constant density
            s_t = s_p;
            constant_density = true;}

         if(T_0<0)   //if negtive, T_0 = n_0, the unit is 10^19 m-3
            T_0 = -0.5*P_0/(T_0*1.e19)/(1.602e-19);   //get real T_0 in eV

         T0=T_0*exp(s_t*(P_min*(exp(-(P_x0-x)/PDped)+(1.0+Pa*(P_x0-x)/PDped)*exp((P_x0-x)/PDped)*Pdel)/(exp(-(P_x0-x)/PDped)+exp((P_x0-x)/PDped))+Pbottom));

         T0 *= 1.602e-19;          //K
         Ti0 = T0;
         if(const_Te0 > 0)
            Te0 = const_Te0 * 1.602e-19;
         else
            Te0 = T0;
         n0 = 0.5*P0/T0;                 

         if(mesh->get(Bbar, "bmag")) // Typical magnetic field
            Bbar = 1.0;

         Lbar = 0.021;
         T0_max = max(T0, true); // Maximum over all processors
         Vbar = sqrt(T0_max/Mi); 
         output.write("using cbm18 grids\n");
         break;


      case 2: 
         if(fit_density)
            n0 = N_0 * exp(-N_d * tanh((x - N_m) / N_w));
         else
            mesh->get(n0, "NiHat");  //normalized
         mesh->get(T0, "TeHat");  //normalized
         n0 *= n0_cyclone;
         T0 *= 2000.*1.602e-19;
         Ti0 = T0;
         if(const_Te0 > 0)
            Te0 = const_Te0 * 1.602e-19;
         else
            Te0 = T0;
         P0 = 2 * n0 * T0;

         // density = n0_cyclone;
         mesh->get(Bbar, "Bbar"); 
         Lbar = 0.585;
         output.write("   Cs = %e m/s \n", sqrt(2000.*1.602e-19/Mi));   
         // Vbar = sqrt(2000. * 1.602e-19 / Mi) * 10;
         Vbar = sqrt(Bbar*Bbar/(MU0*Mi*density));
         mesh->get(Rmajor, "Rmajor");
         Lpar = 2.0*PI*(0.854+2.184*0.5*0.5)*Rmajor; //2*pi*q*R

         output.write("using cyclone grids\n");
         break;

      case 4:

         if(fit_pressure)
            P0=P_0*exp(s_p*(P_min*(exp(-(P_x0-x)/PDped)+(1.0+Pa*(P_x0-x)/PDped)*exp((P_x0-x)/PDped)*Pdel)/(exp(-(P_x0-x)/PDped)+exp((P_x0-x)/PDped))+Pbottom));
         else
            mesh->get(P0, "pressure"); // Pascals

         if(fit_all)
         {
            n0 = N_a * ((1 + N_s * x) * exp(-(x - N_m) / N_w) - exp((x - N_m) / N_w)) / (exp(-(x - N_m) / N_w) + exp((x - N_m) / N_w)) + N_b;
            // n0 = field_larger(n0, Low_limit_n0);
            n0 *= 1.e20;
            T0 = T_a * ((1 + T_s * x) * exp(-(x - T_m) / T_w) - exp((x - T_m) / T_w)) / (exp(-(x - T_m) / T_w) + exp((x - T_m) / T_w)) + T_b;
            T0 *= 1.6e-19;
            P0 = 2. * n0 * T0;
         }
         else if(fit_density)
         {
            n0 = N_a * ((1 + N_s * x) * exp(-(x - N_m) / N_w) - exp((x - N_m) / N_w)) / (exp(-(x - N_m) / N_w) + exp((x - N_m) / N_w)) + N_b;
            n0 = field_larger(n0, Low_limit_n0);
            n0 *= 1.e20;
            T0 = P0 / (2. * n0);
         }
         else if(n0_fake_prof)
         {
            n0 = N0tanh(n0_height, n0_ave, n0_width, n0_center, n0_bottom_x, n0_topgrad, n0_bottomgrad);
            n0 *= 1.e20;
            T0 = P0 / (2. * n0);
         }
         else
         {
            T0 = T0_const * 1.602e-19;
            n0 = P0 / (2. * T0);
         }
         Ti0 = T0;
         Te0 = T0;

         if(mesh->get(Bbar, "bmag")) // Typical magnetic field
            Bbar = 1.0;

         if(mesh->get(Lbar, "rmag")) // Typical length scale
            Lbar = 0.021;

         T0_max = max(T0, true); // Maximum over all processors
         // Vbar = sqrt(T0_max/Mi); 
         Vbar = sqrt(Bbar*Bbar/(MU0*Mi*density));
         output.write("using cbm18 grids\n");
         break;

      case 5:
         mesh->get(P0, "pressure"); // Pascals
         n0 = (P0 / P00) ^ (0.3);
         n0 *= n0_height * 1.e20;
         T0 = P0 / (2. * n0);
         Ti0 = T0;
         Te0 = T0;

         if(mesh->get(Bbar, "bmag")) // Typical magnetic field
            Bbar = 1.0;

         if(mesh->get(Lbar, "rmag")) // Typical length scale
            Lbar = 0.021;

         T0_max = max(T0, true); // Maximum over all processors
         // Vbar = sqrt(T0_max/Mi); 
         Vbar = sqrt(Bbar*Bbar/(MU0*Mi*density));
         output.write("using cbm18 w/ BS grids\n");
         break;

      case 6:

         output.write("using cmod grid\n");
         if (mesh->get(P0, "pressure_s")) // Pascals
         {
            mesh->get(P0, "pressure");
            output.write("Using pressure as P0.\n");
         }
         else
            output.write("Using pressure_s as P0.\n");
         if(mesh->get(n0,  "Niexp")) { // N_i0                                          
            output.write("Error: Cannot read Ni0 from grid\n");
            return 1;
         } 

         if(mesh->get(Ti0,  "Tiexp")) { // T_i0                                         
            output.write("Error: Cannot read Ti0 from grid\n");
            return 1;
         }

         if(mesh->get(Te0,  "Teexp")) { // T_e0  
            output.write("Error: Cannot read Te0 from grid\n");
            return 1;
         }
         n0 *= 1e20;
         Ti0 *= 1.6e-19;
         Te0 *= 1.6e-19;
         T0 = Ti0;

         if(mesh->get(Bbar, "bmag")) // Typical magnetic field
            Bbar = 1.0;

         if(mesh->get(Lbar, "rmag")) // Typical length scale
            Lbar = 0.021;

         T0_max = max(T0, true); // Maximum over all processors
         // Vbar = sqrt(T0_max/Mi); 
         Vbar = sqrt(Bbar*Bbar/(MU0*Mi*density));
         break;

      default : 
         output.write("need to give equilibrium_case\n");
         return 1; 
   }
   /************************************************************************************/

   Tbar = Lbar/Vbar;
   Va = sqrt(Bbar*Bbar/(MU0*Mi*density));        //Using Alfven velocity

   rhoi = sqrt(T0*Mi)/(1.602e-19*B0);     //2D, equilibrium
   omega_bar = 1.602e-19*Bbar/Mi;     
   Cnor = omega_bar*Tbar;                    //normalization coefficient
   Vthermal = sqrt(T0/Mi);                 //thermal velcocity 

   gyroa = 1.0;
   gyrob = -0.5*rhoi*rhoi/(Lbar*Lbar);
   gyroc = 2.0;
   gyrod = 2.0*gyrob;
   gyroe = 4.0*gyrob;
   gyrop = -Lbar * Lbar / (2 * rhoi * rhoi);

   output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
   output.write("                Vbar = %e m/s   Tbar = %e s\n", Vbar, Tbar);
   output.write("Cnor = %e,  Va = %e\n", Cnor, Va); 

   //normalize equlibrium quantities
   n0 = n0/density;    
   T0=T0/(Mi*Vbar*Vbar);   
   Ti0 = Ti0 / (Mi * Vbar * Vbar);
   Te0 = Te0 / (Mi * Vbar * Vbar);
   beta = 2*MU0*P0/Bbar/Bbar;  
   P0 = P0/(density*Mi*Vbar*Vbar);
   J0 = MU0 * Lbar * J0 / B0;
   Vthermal = Vthermal/Vbar;

   Lpar /= Lbar;                               //landau damping parameter 

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

   Pi0 = n0 * Ti0;
   Pe0 = n0 * Te0;
   SAVE_ONCE2(Pi0, Pe0);


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

   if(hyperviscos > 0.0) {
      output.write("    Hyper-viscosity coefficient: %e\n", hyperviscos);
      dump.add(hyper_mu_x, "hyper_mu_x", 1);
   }


   n0_max = max(n0, true);


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


   SAVE_ONCE4(mesh->g_22, hthe, Bpxy, Btxy);
   // Set B field vector

   B0vec.covariant = false;
   B0vec.x = 0.;
   B0vec.y = Bpxy / hthe;
   B0vec.z = 0.;

   //with out phi0

   if(false && !electrostatic)
   {
      nG0.setBoundary("nG0");
      nG0 = gyrob*Delp2(n0);
      nG0.applyBoundary();
      mesh->communicate(nG0);
      nG0 += n0;
      U0 = Cnor*B0*(n0-nG0);  
      nG0 = n0;
   }
   else
   {
      nG0 = n0;
      U0 = 0;
   }
   SAVE_ONCE2(U0, beta);

   /**************** SET VARIABLE LOCATIONS *************/
   ni.setLocation(CELL_YLOW);
   Lambda_i.setLocation(CELL_YLOW);
   Ppar_i.setLocation(CELL_CENTRE);
   Pperp_i.setLocation(CELL_CENTRE);
   U.setLocation(CELL_CENTRE);
   Psi.setLocation(CELL_YLOW);
   Te.setLocation(CELL_CENTRE);
   Ppar_e.setLocation(CELL_CENTRE);
   Pperp_e.setLocation(CELL_CENTRE); 
   P.setLocation(CELL_CENTRE);
   Ppar.setLocation(CELL_CENTRE);
   Pperp.setLocation(CELL_CENTRE);
   Pi.setLocation(CELL_CENTRE);
   Pe.setLocation(CELL_CENTRE);

   Tperp_i.setLocation(CELL_CENTRE);
   Tpar_i.setLocation(CELL_CENTRE);
   ne.setLocation(CELL_YLOW);
   Tperp_e.setLocation(CELL_CENTRE);
   Tpar_e.setLocation(CELL_CENTRE);
   phi.setLocation(CELL_CENTRE);
   gyroni.setLocation(CELL_YLOW);
   phi_f.setLocation(CELL_YLOW);
   gyrophi.setLocation(CELL_CENTRE);
   gyrophi2.setLocation(CELL_CENTRE);
   gyrophi_a.setLocation(CELL_CENTRE);
   gyrophi_b.setLocation(CELL_CENTRE);
   gyroPsi.setLocation(CELL_YLOW);
   gyroPsi_a.setLocation(CELL_YLOW); 
   Jpar.setLocation(CELL_YLOW);
   gyroVpar_i.setLocation(CELL_YLOW);
   Vpar_i.setLocation(CELL_YLOW);
   Vpar_e.setLocation(CELL_YLOW);

   ns.setLocation(CELL_CENTRE);
   Ts.setLocation(CELL_CENTRE);
   Ts1.setLocation(CELL_CENTRE);
   nmid.setLocation(CELL_CENTRE);
   sour2.setLocation(CELL_YLOW);
   gyrosour.setLocation(CELL_YLOW);
   gyropar.setLocation(CELL_CENTRE);
   gyrophi_as.setLocation(CELL_CENTRE);
   gyrophi_ai.setLocation(CELL_CENTRE);
   gyrophi_bi.setLocation(CELL_CENTRE);
   gyrophi_bs1.setLocation(CELL_CENTRE);
   gyrophi_bs2.setLocation(CELL_CENTRE);
   gyroPsi_as.setLocation(CELL_CENTRE);
   gyroPsi_ai.setLocation(CELL_CENTRE);

   Qperp_i.setLocation(CELL_CENTRE);
   Qpar_i.setLocation(CELL_CENTRE);
   Qperp_e.setLocation(CELL_CENTRE);
   Qpar_e.setLocation(CELL_CENTRE);

   mytest.setLocation(CELL_CENTRE);
   mytest3.setLocation(CELL_CENTRE);
   mytesta.setLocation(CELL_CENTRE);
   mytest3a.setLocation(CELL_CENTRE);

   Grad2_Ppar_e.setLocation(CELL_CENTRE);
   /**************** SET EVOLVING VARIABLES *************/

   // Tell BOUT which variables to evolve

   if(compression)
      SOLVE_FOR(Lambda_i);
   if(!electrostatic)
   {
      if(evolve_ne)
         SOLVE_FOR(ne);
      else
         SOLVE_FOR(U);
      SOLVE_FOR(Psi);
      if(!Zero_Te)
      {
         if(evolve_Te)
            SOLVE_FOR(Te);
         else
         {
            SOLVE_FOR(Ppar_e);
            SOLVE_FOR(Pperp_e);
         }
      }
   }
   SOLVE_FOR(ni);
   SOLVE_FOR(Ppar_i);
   SOLVE_FOR(Pperp_i);  

   if(evolve_ne)
      dump.add(U, "U", 1);
   else
      dump.add(ne, "ne", 1);

   dump.add(phi, "phi", 1);
   dump.add(Jpar, "Jpar", 1);
   dump.add(P, "P", 1);
   dump.add(Vpar_i, "Vpar_i", 1);
   dump.add(Vpar_e, "Vpar_e", 1);
   dump.add(Tpar_i, "Tpar_i", 1);
   dump.add(Tperp_i, "Tperp_i", 1);
   dump.add(Tpar_e, "Tpar_e", 1);
   dump.add(Tperp_e, "Tperp_e", 1);

   if(evolve_Te)
   {
      dump.add(Ppar_e, "Ppar_e", 1);
      dump.add(Pperp_e, "Pperp_e", 1);
   }

   if(Landau_damping_i){
      Qperp_i=0.0;
      Qpar_i=0.0;
      dump.add(Qperp_i,  "Qperp_i", 1);
      dump.add(Qpar_i,  "Qpar_i", 1);

      Qpar_i.setBoundary("Qpar_i");
      Qperp_i.setBoundary("Qperp_i");
      Qpar_i.setLocation(CELL_CENTRE);
      Qperp_i.setLocation(CELL_CENTRE);
   }

   if(Landau_damping_e){
      Qperp_e=0.0;
      Qpar_e=0.0;
      dump.add(Qperp_e,  "Qperp_e", 1);
      dump.add(Qpar_e,  "Qpar_e", 1);

      Qpar_e.setBoundary("Qpar_e");
      Qperp_e.setBoundary("Qperp_e");
      Qpar_e.setLocation(CELL_CENTRE);
      Qperp_e.setLocation(CELL_CENTRE);
   }


   if(code_test){

      mytest.setBoundary("ni");
      mytest3.setBoundary("Psi");

      mytesta.setBoundary("mytesta");
      mytest3a.setBoundary("mytest3a");

      mytest2.setBoundary("mytest2");  
      mytest2 = Grad(log(n0),CELL_YLOW);
      mytest2.applyBoundary();
      mesh->communicate(mytest2);

      dump.add(mytest2, "mytest2", 1);
      dump.add(mytest, "mytest", 1);
      dump.add(mytest3, "mytest3", 1);
      dump.add(mytesta, "mytesta", 1);
      dump.add(mytest3a, "mytest3a", 1);



      dump.add(ns, "ns", 1);
      dump.add(gyroni, "gyroni", 1);
      dump.add(nmid, "nmid", 1);
      dump.add(sour1, "sour1", 1);
      dump.add(sour2, "sour2", 1);
      dump.add(sour3, "sour3", 1);
      dump.add(gyrosour, "gyrosour", 1);
      dump.add(gyropar,  "gyropar", 1); 
      dump.add(gyrophi,  "gyrophi", 1); 
      dump.add(gyrophi2,  "gyrophi2", 1); 
      dump.add(gyroPsi, "gyroPsi", 1);
      dump.add(phi_f, "phi_f", 1);
   }

   if(ddtU_Terms_test)
   {
      ddtU_phiU0.setBoundary("U");
      ddtU_phi_fnG0.setBoundary("U");
      ddtU_GradJpar.setBoundary("U");
      ddtU_PsiJ0.setBoundary("U");
      ddtU_GradP.setBoundary("U");
      ddtU_Vpar.setBoundary("U");
      ddtU_B0phiT0.setBoundary("U");

      dump.add(ddtU_phiU0, "ddtU_phiU0", 1);
      dump.add(ddtU_phi_fnG0, "ddtU_phi_fnG0", 1);
      dump.add(ddtU_GradJpar, "ddtU_GradJpar", 1);
      dump.add(ddtU_PsiJ0, "ddtU_PsiJ0", 1);
      dump.add(ddtU_GradP, "ddtU_GradP", 1);
      dump.add(ddtU_Vpar, "ddtU_Vpar", 1);
      dump.add(ddtU_B0phiT0, "ddtU_B0phiT0", 1);
   }

   if(ddtPsi_Terms_test)
   {
      ddtPsi_Gradphi.setBoundary("Psi");
      ddtPsi_Jpar.setBoundary("Psi");
      ddtPsi_Psin0.setBoundary("Psi");
      ddtPsi_PsiT0.setBoundary("Psi");
      ddtPsi_eHall.setBoundary("Psi");
      ddtPsi_hyperresist.setBoundary("Psi");

      dump.add(ddtPsi_Gradphi, "ddtPsi_Gradphi", 1);
      dump.add(ddtPsi_Jpar, "ddtPsi_Jpar", 1);
      dump.add(ddtPsi_Psin0, "ddtPsi_Psin0", 1);
      dump.add(ddtPsi_PsiT0, "ddtPsi_PsiT0", 1);
      dump.add(ddtPsi_eHall, "ddtPsi_eHall", 1);
      dump.add(ddtPsi_hyperresist, "ddtPsi_hyperresist", 1);
   }

   if(GradparJ_test)
   {
      GradparJ_C.setBoundary("U");
      GradparJ_U.setBoundary("U");

      dump.add(GradparJ_C, "GradparJ_C", 1);
      dump.add(GradparJ_U, "GradparJ_U", 1);
   }

   // everything needed to recover physical units
   SAVE_ONCE6(n0, P0, T0, Ti0, Te0, rhoi);
   SAVE_ONCE4(B0, Vthermal, nG0, J0);   

   /********************* toroidal closure  *********************************/ 
   /********
   if(true || toroidal_closure3)
   {
      b0xgradlogB0(B0, modb0xglB0, vgbxhat, vgbyhat, vgbzhat, b0xglB0);
      vgbxhat.setBoundary("modb0xglB0");
      vgbxhat.applyBoundary();
      vgbyhat.setBoundary("modb0xglB0");
      vgbyhat.applyBoundary();
      vgbzhat.setBoundary("modb0xglB0");
      vgbzhat.applyBoundary();

      mod_wd_Vpar.setBoundary("toroidal1");
      mod_wd_Vpar.applyBoundary();
      mod_wd_Tpar.setBoundary("toroidal1");
      mod_wd_Tpar.applyBoundary();
      mod_wd_Tperp.setBoundary("toroidal1");
      mod_wd_Tperp.applyBoundary();

      dump.add(mod_wd_Vpar, "mod_wd_Vpar", 1);
      dump.add(mod_wd_Tpar, "mod_wd_Tpar", 1);
      dump.add(mod_wd_Tperp, "mod_wd_Tperp", 1);

      vgbhat.covariant = false;
      vgbhat.x = vgbxhat;
      vgbhat.y = vgbyhat;
      vgbhat.z = vgbzhat;
      modb0xglB0.setBoundary("modb0xglB0");
      modb0xglB0.applyBoundary();

      vgb = T0/(B0*Cnor)*modb0xglB0*vgbhat;
      k0_toroidal = 0.05;// *zperiod;

      bxgradp = B0vec ^ Grad(P0) / (B0 * B0 * B0);
      bxkappa = B0vec ^ Grad(P0 + B0 * B0 / 2.) / (B0 * B0 * B0);
      mesh->communicate(bxgradp);
      mesh->communicate(bxkappa);

      switch(curv_model)
      {
         case 1:   
         default:
            vcurv = T0/(B0*Cnor)*b0xcv;
            break;
         case 2:
            vcurv = T0 / (B0 * Cnor) * b0xglB0;
            break;
         case 3:
            vcurv = .5 * T0 / (B0 * Cnor) * (b0xcv + b0xglB0);
            break;
         case 4:
            vcurv = T0 / (B0 * Cnor) * bxkappa;
            break;
         case 5:
            vcurv = .5 * T0 / (B0 * Cnor) * (bxkappa + b0xglB0);
            break;
      }

      (vcurv.x).applyBoundary("neumann");
      (vcurv.y).applyBoundary("neumann");
      (vcurv.z).applyBoundary("neumann");

      toroidal1=0;
      toroidal2=0;
      toroidal3=0;
      toroidal4=0;
      toroidal5=0;

      SAVE_ONCE2(vgb,vcurv);
      SAVE_ONCE4(modb0xglB0, vgbxhat, vgbzhat, k0_toroidal);
      SAVE_ONCE4(b0xcv, b0xglB0, bxgradp, bxkappa);
      dump.add(toroidal1,"toroidal1",1);
      dump.add(toroidal2,"toroidal2",1);
      dump.add(toroidal3,"toroidal3",1);
      dump.add(toroidal4,"toroidal4",1); 
      dump.add(toroidal5,"toroidal5",1);   
   }
**********/

   /****************************************************************************/ 
   /////////////// CHECK VACUUM ///////////////////////
   // In vacuum region, initial vorticity should equal zero
   sour1.setBoundary("sour1");
   sour2.setBoundary("sour2");  
   gyropar.setBoundary("phi");
   gyrosour.setBoundary("phi");
   phi.setBoundary("phi");
   gyrophi.setBoundary("phi");
   gyrophi2.setBoundary("phi");
   gyrophi3.setBoundary("phi");
   nmid.setBoundary("phi");


   ns=0.0;
   gyropar=0.0;
   gyrosour=0.0;

   sour1 = 0;
   sour2 = 0;

   if(electrostatic)
   {
      Ppar_e = 0.;
      Ppar_i = 0.;
      U = 0.;
      Psi = 0.;
   }

   if(!restarting) { 
      if(!electrostatic && !evolve_ne)
         ni = -U/(Cnor*B0);      //ne is 0 initially!
      if(nonlinear)
      {
         Tperp_i = (Pperp_i - Ti0 * ni) / n0;
         Tpar_i = (Ppar_i - Ti0 * ni)/ n0;
      }
      else
      {
         Tperp_i = (Pperp_i-Ti0*ni)/n0;
         // mesh->communicate(Tperp_i);
         Tpar_i = (Ppar_i-Ti0*ni)/n0;
         // mesh->communicate(Tpar_i);
      }

      if(!electrostatic && !evolve_ne)
      {
         ne = ni;
         if(quasineutral)
            ne += U/(Cnor*B0);
      }

      if(!Zero_Te && !electrostatic){
         if(evolve_Te)
         {
            Tperp_e = Te;
            mesh->communicate(Tperp_e);
            Tpar_e = Te;
            mesh->communicate(Tpar_e);
         }
         else
         {
            Tperp_e = (Pperp_e-Te0*ne)/n0;
            mesh->communicate(Tperp_e);
            Tpar_e = (Ppar_e-Te0*ne)/n0;
            mesh->communicate(Tpar_e);
         }
      }else{
         Tperp_e = 0;
         Tpar_e = 0;
      }

      /********************* electric field *************************/   
      if(gyroaverage)
      {
         //get gyroaveraged density 
         ns = invert_laplace(ni,ns_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(ns); 
         Ts1 = invert_laplace(-Tperp_i,ns_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(Ts1);  
         Ts = invert_laplace(Ts1+Tperp_i,ns_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(Ts); 
         gyroni = ns -n0*Ts/T0;

         //get potential

         /***/
         if(electrostatic)
         {
            if(electrostatic_poisson_model == 0)
            {
               /***** Old way ******/
               nmid = gyroni / Cnor / n0 * T0;
               mesh->communicate(nmid);
               nmid.applyBoundary();
               gyrosour = nmid;
               gyrosour -= rhoi * rhoi / (Lbar * Lbar) * Delp2(nmid);
               mesh->communicate(gyrosour);
               gyrosour.applyBoundary();
               phid = - 2. * rhoi * rhoi / (Lbar * Lbar);
               phia = 1.;
               phi = invert_laplace(gyrosour, phi_flags, &phia, NULL, &phid);
               mesh->communicate(phi);
               phi.applyBoundary();
               /******/
            }
            else
            {
               /***** Kim's ******/
               nmid = 0.5 * gyroni / Cnor / n0 * T0;
               mesh->communicate(nmid);
               nmid.applyBoundary();
               gyrosour = nmid;
               phid = - 2. * rhoi * rhoi / (Lbar * Lbar);
               phia = 1.;
               gyropar = invert_laplace(gyrosour, phi_flags, &phia, NULL, &phid);
               mesh->communicate(gyropar);
               gyropar.applyBoundary();
               phi = gyropar + gyrosour;
               mesh->communicate(phi);
               phi.applyBoundary();
               /*****/
            }
         }
         else
         {
            if(evolve_ne)
               nmid = (gyroni - ne) / Cnor / n0 * T0;
            else
               nmid = (gyroni - ni - U / (Cnor * B0)) / Cnor / n0 * T0;
            mesh->communicate(nmid);
            gyrosour = -Lbar * Lbar / rhoi / rhoi * nmid;
            gyrosour += Delp2(nmid);
            mesh->communicate(gyrosour);
            gyrosour.applyBoundary();
            if(phi_constant_density)
               phi = invert_laplace(gyrosour, phi_flags, NULL);
            else
               phi = invert_laplace(gyrosour, phi_flags, NULL, &n0, NULL);
            mesh->communicate(phi);
            phi.applyBoundary();
         }

         /********   Old way *****
           nmid = (gyroni-ni-U/(Cnor*B0))/Cnor/n0*T0;
           sour2 = Grad(nmid,CELL_YLOW);   
           sour2.applyBoundary();   
           mesh->communicate(sour2);
           gyrosour = -Lbar*Lbar/rhoi/rhoi*nmid-sour1*sour2;

           if(constant_density)
           gyropar = invert_laplace(gyrosour, phi_flags, NULL);  
           else
           gyropar = invert_laplace(gyrosour, phi_flags, NULL, &n0, NULL);  
           mesh->communicate(gyropar);
           phi = gyropar+nmid;
          *********/

         //get gyroaveraged potential         

         gyrophi = invert_laplace(phi,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi); 
         gyrophi.applyBoundary();
         gyrophi2 = invert_laplace(gyrophi, phi_flags, &gyroa, NULL, &gyrob);
         mesh->communicate(gyrophi2);
         gyrophi2.applyBoundary();
         gyrophi3 = invert_laplace(gyrophi2, phi_flags, &gyroa, NULL, &gyrob);
         mesh->communicate(gyrophi3);
         gyrophi3.applyBoundary();

         phi_f = gyrophi - phi;
         gyrophi_a = 2 * (gyrophi2 - gyrophi);

         /****
         gyrophi_as = invert_laplace(2.*phi,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi_as);
         gyrophi_ai= gyrophi_as-2.*phi;
         gyrophi_a = invert_laplace(gyrophi_ai,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi_a);

         gyrophi_bs2 = invert_laplace(2.*phi,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi_bs2);
         gyrophi_bi=gyrophi_bs2-2.*phi;
         gyrophi_bs1 = invert_laplace(gyrophi_bi,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi_bs1);
         gyrophi_b = invert_laplace(gyrophi_bs1,phi_flags,&gyroa,NULL,&gyrob);
         mesh->communicate(gyrophi_b);     
         ****/
      }
      else
      {
         U = where(P0 - vacuum_pressure, U, 0.0);
         phi = invert_laplace(U * B0 / n0, phi_flags, NULL, &n0, NULL);
         mesh->communicate(phi);
         phi.applyBoundary();

         gyrophi = phi;
         gyrophi2 = phi;
         gyrophi3 = phi;


         phi_f = 0;

         gyrophi_as = 2. * phi;
         gyrophi_ai = 0;
         gyrophi_a = 0;

         gyrophi_bs2 = 2. * phi;
         gyrophi_bi = 0;
         gyrophi_bs1 = 0;
         gyrophi_b = 0;
      }

   }

   /************** SETUP COMMUNICATIONS **************/

   comms.add(ni);
   if(compression)
   {
      comms.add(Lambda_i);
      Vpar_i.setBoundary("Lambda_i");
      Vpar_e.setBoundary("Lambda_i");
   }
   comms.add(Ppar_i);
   comms.add(Pperp_i);
   if(!electrostatic)
   {
      if(evolve_ne)
         comms.add(ne);
      else
         comms.add(U);
      comms.add(Psi);

      if(!Zero_Te){
         if(evolve_Te)
            comms.add(Te);
         else
         {
            comms.add(Ppar_e);
            comms.add(Pperp_e);
         }
      }
   }

   Jpar.setBoundary("Jpar");

   /***/
   toroidal1.setBoundary("toroidal1");
   toroidal2.setBoundary("toroidal2");
   toroidal3.setBoundary("toroidal3");
   toroidal4.setBoundary("toroidal4");
   toroidal5.setBoundary("toroidal5");
   /***/

   sour3.setBoundary("sour3");  

   Grad2_Ppar_e.setBoundary("Grad2_Ppar_e");
   Grad2_Ppar.setBoundary("Grad2_Ppar");

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

   if(nonlinear) {
      //result -= bracket(Psi*B0, f, bm_mag)/B0;
      result -= bracket(Psi, f, bm_mag) * B0; // corrected as 6-field
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

   if(nonlinear) {
      // if(nonlinear_terms > 4)
      //    result -= bracket(gyroPsi, f, bm_mag)*B0;
      // else
         result -= bracket(Psi, f, bm_mag)*B0;
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

   // output.write("t = %f\n", t);
   if(fakerun)
   {
      int rank;
      char filename[200];
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      sprintf(filename, "%sBOUT.dmp.%d.nc", path.c_str(), rank);
      timestep = (int)floor(t);

      if(false)
      {
         output.write("t = %f, timestep = %d from process %d\n", t, timestep, rank);
         output.write("Reading from file %s\n", filename);
      }

      DataFormat *file = data_format(filename);
      file->openr(filename);
      file->setRecord(timestep);
      ni.allocate();
      file->read_rec(**(ni.getData()), "ni", mesh->ngx, mesh->ngy, mesh->ngz);
      U.allocate();
      file->read_rec(**(U.getData()), "U", mesh->ngx, mesh->ngy, mesh->ngz);
      Psi.allocate();
      file->read_rec(**(Psi.getData()), "Psi", mesh->ngx, mesh->ngy, mesh->ngz);
      Lambda_i.allocate();
      file->read_rec(**(Lambda_i.getData()), "Lambda_i", mesh->ngx, mesh->ngy, mesh->ngz);
      Ppar_i.allocate();
      file->read_rec(**(Ppar_i.getData()), "Ppar_i", mesh->ngx, mesh->ngy, mesh->ngz);
      Pperp_i.allocate();
      file->read_rec(**(Pperp_i.getData()), "Pperp_i", mesh->ngx, mesh->ngy, mesh->ngz);
      Ppar_e.allocate();
      file->read_rec(**(Ppar_e.getData()), "Ppar_e", mesh->ngx, mesh->ngy, mesh->ngz);
      Pperp_e.allocate();
      file->read_rec(**(Pperp_e.getData()), "Pperp_e", mesh->ngx, mesh->ngy, mesh->ngz);
      file->close();
   }

   /* 
      mytest = Grad_parP(phi, CELL_CENTRE);
      mytest.applyBoundary();
      mesh->communicate(mytest);

      mytest3 = Delp2(mytest);
      mytest3.applyBoundary();
      mesh->communicate(mytest3);

      mytesta = Delp2(phi);
      mytesta.applyBoundary();
      mesh->communicate(mytesta);

      mytest3a = Grad_parP(mytesta, CELL_CENTRE);
      mytest3a.applyBoundary();
      mesh->communicate(mytest3a);
      */

   mesh->communicate(comms);

   ////////////////////////////////////////////
   // Transitions from 0 in core to 1 in vacuum
   if(nonlinear) 
      eta = core_resist + (vac_resist - core_resist) * vac_mask;

   if(!electrostatic && !evolve_ne)
   {
      ne = ni;
      if(quasineutral)
         ne += U/(Cnor*B0);
   }

   if(nonlinear)
   {
      Ni_tmp = field_larger(n0 + ni, Low_limit);
      Ne_tmp = field_larger(n0 + ne, Low_limit);
   }

   if(nonlinear)
   {
      Tperp_i = (Pperp_i-Ti0*ni)/n0;
      mesh->communicate(Tperp_i);
      Tpar_i = (Ppar_i-Ti0*ni)/n0;
      mesh->communicate(Tpar_i);
   }
   else
   {
      Tperp_i = (Pperp_i-Ti0*ni)/n0;
      mesh->communicate(Tperp_i);
      Tpar_i = (Ppar_i-Ti0*ni)/n0;
      mesh->communicate(Tpar_i);
   }

   if(!Zero_Te && !electrostatic){
      if(evolve_Te)
      {
         Tpar_e = Te;
         Tperp_e = Te;
         Ppar_e = ne * Te0 + n0 * Te;
         Pperp_e = ne * Te0 + n0 * Te;
         if(nonlinear)
         {
            Ppar_e += ne * Te;
            Pperp_e += ne * Te;
         }
      }
      else
      {
         if(nonlinear)
         {
            Tperp_e = (Pperp_e-Te0*ne)/Ne_tmp;
            mesh->communicate(Tperp_e);
            Tpar_e = (Ppar_e-Te0*ne)/Ne_tmp;     
            mesh->communicate(Tpar_e);
         }
         else
         {
            Tperp_e = (Pperp_e-Te0*ne)/n0;
            mesh->communicate(Tperp_e);
            Tpar_e = (Ppar_e-Te0*ne)/n0;     
            mesh->communicate(Tpar_e);
         }
      }
   }else{
      Tperp_e = 0;
      Tpar_e = 0;
      Ppar_e = Te0*ni;
      Pperp_e = Te0*ni;
   }

   Pi = Ppar_i + Pperp_i;
   // mesh->communicate(Pi);
   Pe = Ppar_e + Pperp_e;
   // mesh->communicate(Pe);
   Ppar = Ppar_i + Ppar_e;
   // mesh->communicate(Ppar);
   Pperp = Pperp_i + Pperp_e;
   // mesh->communicate(Pperp);
   P = Ppar_i+Ppar_e+Pperp_i+Pperp_e;  
   // mesh->communicate(P);

   /********************* electric field *************************/   
   if(gyroaverage)
   {
      //get gyroaveraged density 
      ns = invert_laplace(ni,ns_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(ns); 
      Ts1 = invert_laplace(-Tperp_i,ns_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(Ts1);  
      Ts = invert_laplace(Ts1+Tperp_i,ns_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(Ts); 
      gyroni = ns -n0*Ts/T0; 

      //get potential

      /**********/
      if(electrostatic)
      {
         if(electrostatic_poisson_model == 0)
         {
            /*********/
            nmid = gyroni / Cnor / n0 * T0;
            mesh->communicate(nmid);
            nmid.applyBoundary();
            gyrosour = nmid;
            gyrosour -= rhoi * rhoi / (Lbar * Lbar) * Delp2(nmid);
            mesh->communicate(gyrosour);
            gyrosour.applyBoundary();
            phid = - 2. * rhoi * rhoi / (Lbar * Lbar);
            phia = 1.;
            phi = invert_laplace(gyrosour, phi_flags, &phia, NULL, &phid);
            mesh->communicate(phi);
            phi.applyBoundary();
            /*********/
         }
         else
         {
            /***** Kim's ******/
            nmid = 0.5 * gyroni / Cnor / n0 * T0;
            mesh->communicate(nmid);
            nmid.applyBoundary();
            gyrosour = nmid;
            phid = - 2. * rhoi * rhoi / (Lbar * Lbar);
            phia = 1.;
            gyropar = invert_laplace(gyrosour, phi_flags, &phia, NULL, &phid);
            mesh->communicate(gyropar);
            gyropar.applyBoundary();
            phi = gyropar + gyrosour;
            mesh->communicate(phi);
            phi.applyBoundary();
            /*********/
         }
      }
      else
      {
         if(evolve_ne)
            nmid = (gyroni - ne) / Cnor / n0 * T0;
         else
            nmid = (gyroni - ni - U / (Cnor * B0)) / Cnor / n0 * T0;
         mesh->communicate(nmid);
         gyrosour = -Lbar * Lbar / rhoi / rhoi * nmid;
         gyrosour += Delp2(nmid);
         mesh->communicate(gyrosour);
         gyrosour.applyBoundary();
         if(phi_constant_density)
            phi = invert_laplace(gyrosour, phi_flags, NULL);
         else
            phi = invert_laplace(gyrosour, phi_flags, NULL, &n0, NULL);
         mesh->communicate(phi);
         phi.applyBoundary();
      }
      /*****/

      /***********       Old way   ************
        nmid = (gyroni-ni-U/(Cnor*B0))/Cnor/n0*T0;
        sour2 = Grad(nmid,CELL_YLOW);
        sour2.applyBoundary();
        mesh->communicate(sour2);
        gyrosour = -Lbar*Lbar/rhoi/rhoi*nmid-sour1*sour2;

        if(constant_density)
        gyropar = invert_laplace(gyrosour, phi_flags, NULL);  
        else
        gyropar = invert_laplace(gyrosour, phi_flags, NULL, &n0, NULL);  

        mesh->communicate(gyropar);
        phi=gyropar+nmid;                       //not divided by B0 
       ***********/

      //get gyroaveraged potential

      gyrophi = invert_laplace(phi,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi); 
      gyrophi.applyBoundary();
      gyrophi2 = invert_laplace(gyrophi, phi_flags, &gyroa, NULL, &gyrob);
      mesh->communicate(gyrophi2);
      gyrophi2.applyBoundary();
      gyrophi3 = invert_laplace(gyrophi2, phi_flags, &gyroa, NULL, &gyrob);
      mesh->communicate(gyrophi3);
      gyrophi3.applyBoundary();

      phi_f = gyrophi - phi;
      gyrophi_a = 2 * (gyrophi2 - gyrophi);

      /***
      gyrophi_as = invert_laplace(2.*phi,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi_as);
      gyrophi_ai= gyrophi_as-2.*phi;
      gyrophi_a = invert_laplace(gyrophi_ai,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi_a);

      gyrophi_bs2 = invert_laplace(2.*phi,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi_bs2);
      gyrophi_bi=gyrophi_bs2-2.*phi;
      gyrophi_bs1 = invert_laplace(gyrophi_bi,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi_bs1);
      gyrophi_b = invert_laplace(gyrophi_bs1,phi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyrophi_b);  
      ***/
   }
   else
   {
      U = where(P0 - vacuum_pressure, U, 0.0);
      phi = invert_laplace(U * B0 / n0, phi_flags, NULL, &n0, NULL);
      mesh->communicate(phi);
      phi.applyBoundary();

      gyrophi = phi;
      gyrophi2 = phi;
      gyrophi3 = phi;

      phi_f = 0;

      gyrophi_as = 2. * phi;
      gyrophi_ai = 0;
      gyrophi_a = 0;

      gyrophi_bs2 = 2. * phi;
      gyrophi_bi = 0;
      gyrophi_bs1 = 0;
      gyrophi_b = 0;
   }

   /*************************** Psi and current *********************************/

   if(gyroaverage)
   {
      //get gyroaveraged vector potential
      gyroPsi = invert_laplace(Psi,Psi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyroPsi);

      gyroPsi_as = invert_laplace(2.*Psi,Psi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyroPsi_as);
      gyroPsi_ai= gyroPsi_as-2.*Psi;
      gyroPsi_a = invert_laplace(gyroPsi_ai,Psi_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyroPsi_a);   
   }
   else
   {
      //get gyroaveraged vector potential
      gyroPsi = Psi;

      gyroPsi_as = 2 * Psi;
      gyroPsi_ai = 0;
      gyroPsi_a = 0;
   }

   Jpar = Delp2(Psi);
   Jpar.applyBoundary();
   mesh->communicate(Jpar);   

   // Smooth j in x
   if(smooth_j_x)
      Jpar = smooth_x(Jpar);

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

   /*************** ion and electron parallel velocity ***********************/

   if(compression)
   {
      Vpar_i = Lambda_i;
      if(solve_lambda)
         Vpar_i -= Cnor*B0*gyroPsi;
      gyroVpar_i = invert_laplace(Vpar_i,Vpar_flags,&gyroa,NULL,&gyrob);
      mesh->communicate(gyroVpar_i);

      Vpar_e = gyroVpar_i;
      if(shear_Alfven_wave)
      {
         if(nonlinear)
            Vpar_e += Va*Va/Vbar/Vbar*B0*Jpar/(Cnor*Ni_tmp);   
         else
            Vpar_e += Va*Va/Vbar/Vbar*B0*Jpar/(Cnor*n0);   
      }
      Vpar_e.applyBoundary();
      mesh->communicate(Vpar_e);
   }

   /*********************** Landau_damping ***********************************/
   // k0_landau = 0.314/Lpar;
   if(Equilibrium_case == 2)
      k0_landau = 0.05*2.*PI / (2.0*PI*(0.854+2.184*0.5*0.5)*Rmajor);
   else
      k0_landau = 2. / (20. * (2 * PI * 5.));

   /***
   if(Landau_damping_i){
      if(nonlinear)
      {
         Qpar_i = -sqrt(8./PI)*Vthermal*iSign_kpar(Ppar_i, k0_landau);
         Qperp_i = -sqrt(2./PI)*Vthermal*iSign_kpar(Pperp_i, k0_landau);
      }
      else
      {
         Qpar_i = -n0*sqrt(8./PI)*Vthermal*iSign_kpar(Tpar_i, k0_landau);
         Qperp_i = -n0*sqrt(2./PI)*Vthermal*iSign_kpar(Tperp_i, k0_landau);
         Qperp_i -= n0*sqrt(2./PI)*Vthermal*iSign_kpar(Cnor * (gyrophi2 - gyrophi), k0_landau);
      }
      mesh->communicate(Qpar_i);
      mesh->communicate(Qperp_i);
      Qpar_i.applyBoundary();
      Qperp_i.applyBoundary();
   }

   if(Landau_damping_e){
      // Qpar_e = -sqrt(8./PI)*sqrt(1836)*Vthermal*iSign_kpar(Ppar_e, k0_landau);
      // Qperp_e = -sqrt(2./PI)*sqrt(1836)*Vthermal*iSign_kpar(Pperp_e, k0_landau);
      Qpar_e = -n0*sqrt(8./PI)*sqrt(1836)*Vthermal*iSign_kpar(Tpar_e, k0_landau);
      Qperp_e = -n0*sqrt(2./PI)*sqrt(1836)*Vthermal*iSign_kpar(Tperp_e, k0_landau);
      mesh->communicate(Qpar_e);
      mesh->communicate(Qperp_e);
      Qpar_e.applyBoundary();
      Qperp_e.applyBoundary();
   }  
   ****/

   ///////////// Mod omegad closures

   /**********
   if(toroidal_closure3)
   {
      if(compression)
      {
      mod_wd_Vpar = mod_wd(Vpar_i, &vcurv, k0_toroidal);
         mesh->communicate(mod_wd_Vpar);
         mod_wd_Vpar.applyBoundary();
      }
      mod_wd_Tpar = mod_wd(Tpar_i, &vcurv, k0_toroidal);
      mod_wd_Tperp = mod_wd(Tperp_i, &vcurv, k0_toroidal);
      mesh->communicate(mod_wd_Tpar);
      mesh->communicate(mod_wd_Tperp);
      mod_wd_Tpar.applyBoundary();
      mod_wd_Tperp.applyBoundary();
   }
   *****/


   /****************************** evolving equations ********************************************/

   //ion density equation
   ddt(ni) = 0.; 

   ddt(ni) -= bracket(gyrophi, nG0, bm_exb);

   if(electrostatic)
   {
      switch(curv_model)
      {
         case 1:
         default:
            ddt(ni) -= 1./(Cnor*B0)* b0xcv * Grad(Pi);
            break;
         case 2:
            ddt(ni) -= 1./(Cnor*B0) * bracket(B0, Pi, bm_exb);
            break;
         case 3:
            ddt(ni) -= .5/(Cnor*B0)* b0xcv * Grad(Pi);
            ddt(ni) -= .5/(Cnor*B0) * bracket(B0, Pi, bm_exb);
            break;
         case 4:
            ddt(ni) -= 1./(Cnor*B0)* bxkappa * Grad(Pi);
            break;
         case 5:
            ddt(ni) -= .5/(Cnor*B0)* bxkappa * Grad(Pi);
            ddt(ni) -= .5/(Cnor*B0) * bracket(B0, Pi, bm_exb);
            break;
      }
   }

   if(continuity)
   {
      if(isotropic)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(ni) -= 2. * n0 / B0 * b0xcv * Grad(phi);
               break;
            case 2:
               ddt(ni) -= 2. * n0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 3:
               ddt(ni) -= 1. * n0 / B0 * b0xcv * Grad(phi);
               ddt(ni) -= 1. * n0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 4:
               ddt(ni) -= 2. * n0 / B0 * bxkappa * Grad(phi);
               break;
            case 5:
               ddt(ni) -= 1. * n0 / B0 * bxkappa * Grad(phi);
               ddt(ni) -= 1. * n0 / B0 * bracket(B0, phi, bm_exb);
               break;
         }

         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(ni) -= 2. * ni / B0 * b0xcv * Grad(phi);
                  break;
               case 2:
                  ddt(ni) -= 2. * ni / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 3:
                  ddt(ni) -= 1. * ni / B0 * b0xcv * Grad(phi);
                  ddt(ni) -= 1. * ni / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 4:
                  ddt(ni) -= 2. * ni / B0 * bxkappa * Grad(phi);
                  break;
               case 5:
                  ddt(ni) -= 1. * ni / B0 * bxkappa * Grad(phi);
                  ddt(ni) -= 1. * ni / B0 * bracket(B0, phi, bm_exb);
                  break;
            }
         }
      }
      else
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(ni) -= n0/B0*b0xcv * Grad(gyrophi);
               ddt(ni) -= n0/B0*b0xcv * Grad(gyrophi2);
               break;
            case 2:
               ddt(ni) -= n0/B0* bracket(B0, gyrophi, bm_exb);
               ddt(ni) -= n0/B0* bracket(B0, gyrophi2, bm_exb);
               break;
            case 3:
               ddt(ni) -= .5 * n0/B0*b0xcv * Grad(gyrophi);
               ddt(ni) -= .5 * n0/B0*b0xcv * Grad(gyrophi2);
               ddt(ni) -= .5 * n0/B0* bracket(B0, gyrophi, bm_exb);
               ddt(ni) -= .5 * n0/B0* bracket(B0, gyrophi2, bm_exb);
               break;
            case 4:
               ddt(ni) -= n0/B0*bxkappa * Grad(gyrophi);
               ddt(ni) -= n0/B0*bxkappa * Grad(gyrophi2);
               break;
            case 5:
               ddt(ni) -= .5 * n0/B0*bxkappa * Grad(gyrophi);
               ddt(ni) -= .5 * n0/B0*bxkappa * Grad(gyrophi2);
               ddt(ni) -= .5 * n0/B0* bracket(B0, gyrophi, bm_exb);
               ddt(ni) -= .5 * n0/B0* bracket(B0, gyrophi2, bm_exb);
               break;
         }

         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(ni) -= ni/B0*b0xcv * Grad(gyrophi);
                  ddt(ni) -= ni/B0*b0xcv * Grad(gyrophi2);
                  break;
               case 2:
                  ddt(ni) -= ni/B0* bracket(B0, gyrophi, bm_exb);
                  ddt(ni) -= ni/B0* bracket(B0, gyrophi2, bm_exb);
                  break;
               case 3:
                  ddt(ni) -= .5 * ni/B0*b0xcv * Grad(gyrophi);
                  ddt(ni) -= .5 * ni/B0*b0xcv * Grad(gyrophi2);
                  ddt(ni) -= .5 * ni/B0* bracket(B0, gyrophi, bm_exb);
                  ddt(ni) -= .5 * ni/B0* bracket(B0, gyrophi2, bm_exb);
                  break;
               case 4:
                  ddt(ni) -= ni/B0*bxkappa * Grad(gyrophi);
                  ddt(ni) -= ni/B0*bxkappa * Grad(gyrophi2);
                  break;
               case 5:
                  ddt(ni) -= .5 * ni/B0*bxkappa * Grad(gyrophi);
                  ddt(ni) -= .5 * ni/B0*bxkappa * Grad(gyrophi2);
                  ddt(ni) -= .5 * ni/B0* bracket(B0, gyrophi, bm_exb);
                  ddt(ni) -= .5 * ni/B0* bracket(B0, gyrophi2, bm_exb);
                  break;
            }
         }
      }

      if(compression)
      {
         ddt(ni) -= n0 * B0 * Grad_parP(Vpar_i/B0, CELL_CENTRE); // Kim line 1
         
         if(nonlinear && nonlinear_terms > 2)
            ddt(ni) -= ni * B0 * Grad_par(Vpar_i/B0, CELL_CENTRE); // Kim line 1
      }
   }


   if(FLR_effect)
      ddt(ni) -= nG0 / Ti0 * bracket(gyrophi2 - gyrophi, Ti0, bm_exb); // Kim line 3

   if(nonlinear && nonlinear_terms > 0){
      ddt(ni) -= bracket(gyrophi, ni, bm_exb);
      if(nonlinear_terms > 4)
         ddt(ni) -= bracket(gyrophi2 - gyrophi, n0 * Tperp_i / Ti0, bm_exb);
      // ddt(ni) -= 0.5*n0/Ti0*bracket(gyrophi_a, Tperp_i, bm_exb);  
   } 

   if(viscos_par > 0.0)
      ddt(ni) += viscos_par*Grad2_par2(ni); 

   if(viscos_perp > 0.0)
      ddt(ni) += viscos_perp * Delp2(ni);

   if(sink_nil > 0.0){
      ddt(ni) -=  sink_nil*sink_tanhxl(n0,ni,sni_widthl,sni_lengthl); // core sink
   }

   if(sink_nir > 0.0){
      ddt(ni) -=  sink_nir*sink_tanhxr(n0,ni,sni_widthr,sni_lengthr); //  sol sink
   }


   ////////////////////////////////////////////////////
   //ion parallel motion 

   if(compression)
   {
      ddt(Lambda_i) = 0.; 

      if(!electrostatic)
      {
         ddt(Lambda_i) -= Grad_parP(Ppar, CELL_YLOW) / n0;
         ddt(Lambda_i) -= bracket(-Psi, P0, bm_mag) * B0 / n0;

         switch(curv_model)
         {
            case 1:
            default:
               ddt(Lambda_i) -= 4./(Cnor*B0)*V_dot_Grad(Ti0 * b0xcv, Vpar_i);   // Kim line 3 
               break;
            case 2:
               ddt(Lambda_i) -= 4.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
            case 3:
               ddt(Lambda_i) -= 2./(Cnor*B0)*V_dot_Grad(Ti0 * b0xcv, Vpar_i);   // Kim line 3 
               ddt(Lambda_i) -= 2.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
            case 4:
               ddt(Lambda_i) -= 4./(Cnor*B0)*V_dot_Grad(Ti0 * bxkappa, Vpar_i);   // Kim line 3 
               break;
            case 5:
               ddt(Lambda_i) -= 2./(Cnor*B0)*V_dot_Grad(Ti0 * bxkappa, Vpar_i);   // Kim line 3 
               ddt(Lambda_i) -= 2.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
         }
         ddt(Lambda_i) -= Tperp_i/B0*Grad_par(B0, CELL_CENTRE);  // Kim line 2
         ddt(Lambda_i) -= Cnor*(gyrophi2 - gyrophi)/B0*Grad_parP0(B0, CELL_CENTRE); // Kim line 2
      }
      else
      {
         ddt(Lambda_i) -= Grad_parP(Ppar_i, CELL_YLOW) / n0;
         ddt(Lambda_i) -= Cnor*Grad_parP0(gyrophi, CELL_YLOW);             // Kim line 1
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Lambda_i) -= 4./(Cnor*B0)*V_dot_Grad(Ti0 * b0xcv, Vpar_i);   // Kim line 3 
               break;
            case 2:
               ddt(Lambda_i) -= 4.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
            case 3:
               ddt(Lambda_i) -= 2./(Cnor*B0)*V_dot_Grad(Ti0 * b0xcv, Vpar_i);   // Kim line 3 
               ddt(Lambda_i) -= 2.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
            case 4:
               ddt(Lambda_i) -= 4./(Cnor*B0)*V_dot_Grad(Ti0 * bxkappa, Vpar_i);   // Kim line 3 
               break;
            case 5:
               ddt(Lambda_i) -= 2./(Cnor*B0)*V_dot_Grad(Ti0 * bxkappa, Vpar_i);   // Kim line 3 
               ddt(Lambda_i) -= 2.*Ti0/(Cnor*B0)* bracket(B0, Vpar_i, bm_exb);   // Kim line 3 
               break;
         }
         ddt(Lambda_i) -= Tperp_i/B0*Grad_par(B0, CELL_CENTRE);  // Kim line 2
         ddt(Lambda_i) -= Cnor*(gyrophi2 - gyrophi)/B0*Grad_parP0(B0, CELL_CENTRE); // Kim line 2
      }

      // ddt(Lambda_i) -= B0/n0*Grad_parP_gyro(Ppar_i/B0, CELL_CENTRE);  // Kim line 1

      if(false && !solve_lambda && !electrostatic)
      {
         ddt(Lambda_i) -= B0 / n0 * Grad_parP_gyro(Ppar_e / B0, CELL_CENTRE);

         ddt(Lambda_i) -= 2.*Tperp_i/B0*Grad_parP0(B0, CELL_CENTRE);  // Kim line 2
         ddt(Lambda_i) -= 0.5*Cnor*(gyrophi2 - gyrophi)/B0*Grad_parP0(B0, CELL_CENTRE); // Kim line 2
      }

      if(false && !electrostatic)
         ddt(Lambda_i) += Ti0 / n0 * bracket(Psi * B0, nG0, bm_mag);

      if(solve_lambda)
      {
         ddt(Lambda_i) += T0 / n0 * bracket(gyroPsi * B0, nG0, bm_mag);
         ddt(Lambda_i) += 1. / B0 * b0xGrad_dot_Grad(gyroPsi*B0, T0);
         ddt(Lambda_i) += 0.5 * B0 * bracket(gyroPsi_a, T0, bm_mag);
      }

      if(nonlinear && nonlinear_terms > 0){
         ddt(Lambda_i) -= bracket(gyrophi, Lambda_i, bm_exb);
         // ddt(Lambda_i) += 0.5*bracket(gyroPsi_a, Tperp_i, bm_mag)*B0;  
      }  

      if(toroidal_closure2)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Lambda_i) -= 2./(Cnor*B0)* V_dot_Grad(Ti0 * nu5i * b0xcv, Vpar_i);
               break;
            case 2:
               ddt(Lambda_i) -= 2.*Ti0/(Cnor*B0)* bracket(nu5i * B0, Vpar_i, bm_exb);
               break;
            case 3:
               ddt(Lambda_i) -= 1./(Cnor*B0)* V_dot_Grad(Ti0 * nu5i * b0xcv, Vpar_i);
               ddt(Lambda_i) -= 1.*Ti0/(Cnor*B0)* bracket(nu5i*B0, Vpar_i, bm_exb);
               break;
            case 4:
               ddt(Lambda_i) -= 2./(Cnor*B0)* V_dot_Grad(Ti0 * nu5i*bxkappa, Vpar_i);
               break;
            case 5:
               ddt(Lambda_i) -= 1./(Cnor*B0)* V_dot_Grad(Ti0 * nu5i*bxkappa, Vpar_i);
               ddt(Lambda_i) -= 1.*Ti0/(Cnor*B0)* bracket(nu5i*B0, Vpar_i, bm_exb);
               break;
         }
      }

      if(toroidal_closure3){
         /****/
         ddt(Lambda_i) -= 2. * nu5r * mod_wd_Vpar;
         /****
           toroidal5 = 1.*nu5r*(mod_wd(Vpar_i,&vgb,k0_toroidal)+mod_wd(Vpar_i,&vcurv,k0_toroidal));
           toroidal5.applyBoundary();
           mesh->communicate(toroidal5);
           ddt(Lambda_i) -= toroidal5;
          ***/
      }

      if(sink_Lambda_ir > 0.0){
         ddt(Lambda_i) -=  sink_Lambda_ir*sink_tanhxr(P0,Pperp_e,sLambda_i_widthr,sLambda_i_lengthr); //  sol sink
      }
   }

   ////////////////////////////////////////////////////
   //papallel pressure  
   ddt(Ppar_i) = 0.;

   // ddt(Ppar_i) -= bracket(phi, P0 / 2., bm_exb);

   ddt(Ppar_i) -= Ti0 * bracket(gyrophi, nG0, bm_exb);
   ddt(Ppar_i) -= nG0 * bracket(gyrophi, Ti0, bm_exb);

   if(FLR_effect)
      ddt(Ppar_i) -= nG0 * bracket(gyrophi2 - gyrophi, Ti0, bm_exb);  // Kim line 1

   if(continuity)
   {
      if(isotropic)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Ppar_i) -= 10. / 3. * Pi0 / B0 * b0xcv * Grad(phi);
               break;
            case 2:
               ddt(Ppar_i) -= 10. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 3:
               ddt(Ppar_i) -= 5. / 3. * Pi0 / B0 * b0xcv * Grad(phi);
               ddt(Ppar_i) -= 5. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 4:
               ddt(Ppar_i) -= 10. / 3. * Pi0 / B0 * bxkappa * Grad(phi);
               break;
            case 5:
               ddt(Ppar_i) -= 5. / 3. * Pi0 / B0 * bxkappa * Grad(phi);
               ddt(Ppar_i) -= 5. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
         }

         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Ppar_i) -= 10. / 3. * Ppar_i / B0 * b0xcv * Grad(phi);
                  break;
               case 2:
                  ddt(Ppar_i) -= 10. / 3. * Ppar_i / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 3:
                  ddt(Ppar_i) -= 5. / 3. * Ppar_i / B0 * b0xcv * Grad(phi);
                  ddt(Ppar_i) -= 5. / 3. * Ppar_i / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 4:
                  ddt(Ppar_i) -= 10. / 3. * Ppar_i / B0 * bxkappa * Grad(phi);
                  break;
               case 5:
                  ddt(Ppar_i) -= 5. / 3. * Ppar_i / B0 * bxkappa * Grad(phi);
                  ddt(Ppar_i) -= 5. / 3. * Ppar_i / B0 * bracket(B0, phi, bm_exb);
                  break;
            }
         }
      }
      else
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Ppar_i) -= 3. * Pi0 / B0 * b0xcv * Grad(gyrophi);    // Kim line 3
               ddt(Ppar_i) -= n0 * Ti0 / B0 * b0xcv * Grad(gyrophi2); // Kim line 3
               break;
            case 2:
               ddt(Ppar_i) -= 3. * Pi0 / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
               ddt(Ppar_i) -= n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
               break;
            case 3:
               ddt(Ppar_i) -= 1.5 * Pi0 / B0 * b0xcv * Grad(gyrophi);    // Kim line 3
               ddt(Ppar_i) -= .5 * n0 * Ti0 / B0 * b0xcv * Grad(gyrophi2); // Kim line 3
               ddt(Ppar_i) -= 1.5 * Pi0 / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
               ddt(Ppar_i) -= .5 * n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
               break;
            case 4:
               ddt(Ppar_i) -= 3. * Pi0 / B0 * bxkappa * Grad(gyrophi);    // Kim line 3
               ddt(Ppar_i) -= n0 * Ti0 / B0 * bxkappa * Grad(gyrophi2); // Kim line 3
               break;
            case 5:
               ddt(Ppar_i) -= 1.5 * Pi0 / B0 * bxkappa * Grad(gyrophi);    // Kim line 3
               ddt(Ppar_i) -= .5 * n0 * Ti0 / B0 * bxkappa * Grad(gyrophi2); // Kim line 3
               ddt(Ppar_i) -= 1.5 * Pi0 / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
               ddt(Ppar_i) -= .5 * n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
               break;
         }

         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Ppar_i) -= 3. * Ppar_i / B0 * b0xcv * Grad(gyrophi);    // Kim line 3
                  ddt(Ppar_i) -= Ppar_i / B0 * b0xcv * Grad(gyrophi2); // Kim line 3
                  break;
               case 2:
                  ddt(Ppar_i) -= 3. * Ppar_i / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
                  ddt(Ppar_i) -= Ppar_i / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
                  break;
               case 3:
                  ddt(Ppar_i) -= 1.5 * Ppar_i / B0 * b0xcv * Grad(gyrophi);    // Kim line 3
                  ddt(Ppar_i) -= .5 * Ppar_i / B0 * b0xcv * Grad(gyrophi2); // Kim line 3
                  ddt(Ppar_i) -= 1.5 * Ppar_i / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
                  ddt(Ppar_i) -= .5 * Ppar_i / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
                  break;
               case 4:
                  ddt(Ppar_i) -= 3. * Ppar_i / B0 * bxkappa * Grad(gyrophi);    // Kim line 3
                  ddt(Ppar_i) -= Ppar_i / B0 * bxkappa * Grad(gyrophi2); // Kim line 3
                  break;
               case 5:
                  ddt(Ppar_i) -= 1.5 * Ppar_i / B0 * bxkappa * Grad(gyrophi);    // Kim line 3
                  ddt(Ppar_i) -= .5 * Ppar_i / B0 * bxkappa * Grad(gyrophi2); // Kim line 3
                  ddt(Ppar_i) -= 1.5 * Ppar_i / B0 * bracket(B0, gyrophi, bm_exb);    // Kim line 3
                  ddt(Ppar_i) -= .5 * Ppar_i / B0 * bracket(B0, gyrophi2, bm_exb); // Kim line 3
                  break;
            }
         }
      }

      if(compression)
      {
         if(isotropic)
         {
            ddt(Ppar_i) -= 2. / 3. * Pi0 * B0 * Grad_parP(Vpar_i / B0, CELL_CENTRE);
            if(nonlinear && nonlinear_terms > 2)
               ddt(Ppar_i) -= 2. / 3. * Ppar_i * B0 * Grad_par(Vpar_i / B0, CELL_CENTRE);
         }
         else
         {
            ddt(Ppar_i) -= 3. * Pi0 * B0 * Grad_parP_gyro(Vpar_i / B0, CELL_CENTRE); // Kim line 2
            if(nonlinear && nonlinear_terms > 2)
               ddt(Ppar_i) -= 3. * Ppar_i * B0 * Grad_par(Vpar_i / B0, CELL_CENTRE); // Kim line 2
         }
         // ddt(Ppar_i) -= 2.*T0*n0*Vpar_i/B0*Grad_parP0(B0, CELL_CENTRE);
      }
   }

   if(energy_flux)
   {
      if(isotropic)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pi0 * b0xcv, Tpar_i + Tperp_i);
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * b0xcv * Grad(Ti0);
               break;
            case 2:
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
            case 3:
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pi0 * b0xcv, Tpar_i + Tperp_i);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * b0xcv * Grad(Ti0);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
            case 4:
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pi0 * bxkappa, Tpar_i + Tperp_i);
               ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bxkappa * Grad(Ti0);
               break;
            case 5:
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pi0 * bxkappa, Tpar_i + Tperp_i);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bxkappa * Grad(Ti0);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 3)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Ppar_i * b0xcv, Tpar_i + Tperp_i);
                  break;
               case 2:
                  ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * Ppar_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
               case 3:
                  ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Ppar_i * b0xcv, Tpar_i + Tperp_i);
                  ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * Ppar_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
               case 4:
                  ddt(Ppar_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Ppar_i * bxkappa, Tpar_i + Tperp_i);
                  break;
               case 5:
                  ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Ppar_i * bxkappa, Tpar_i + Tperp_i);
                  ddt(Ppar_i) -= 5. / 6. / (Cnor * B0) * Ppar_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
            }
         }
      }
      else
      {
         switch(curv_model)
         {
            case 1:
            default:
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 4.* b0xcv * Grad(Ti0 * ni);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv * 7., Ppar_i);
               ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv , Pperp_i);
               break;
            case 2:
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 4.* bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= 1. /(Cnor*B0) * Ti0 * bracket(B0, 7. * Ppar_i, bm_exb);
               ddt(Ppar_i) -= 1. /(Cnor*B0) * Ti0 * bracket(B0, Pperp_i, bm_exb);
               break;
            case 3:
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 2.* b0xcv * Grad(Ti0 * ni);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv * 7., Ppar_i);
               ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv , Pperp_i);
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 2.* bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0, 7. * Ppar_i, bm_exb);
               ddt(Ppar_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0, Pperp_i, bm_exb);
               break;
            case 4:
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 4.* bxkappa * Grad(Ti0 * ni);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa * 7., Ppar_i);
               ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa , Pperp_i);
               break;
            case 5:
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 2.* bxkappa * Grad(Ti0 * ni);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa * 7., Ppar_i);
               ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa , Pperp_i);
               // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
               ddt(Ppar_i) += Ti0/(Cnor*B0) * 2.* bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
               ddt(Ppar_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0, 7. * Ppar_i, bm_exb);
               ddt(Ppar_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0, Pperp_i, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 3)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 4.* b0xcv * Grad(Ti0 * ni);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tpar_i * b0xcv * 7., Ppar_i);
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tpar_i * b0xcv , Pperp_i);
                  break;
               case 2:
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 4.* bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * Tpar_i * bracket(B0, 7. * Ppar_i, bm_exb);
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * Tpar_i * bracket(B0, Pperp_i, bm_exb);
                  break;
               case 3:
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 2.* b0xcv * Grad(Ti0 * ni);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tpar_i * b0xcv * 7., Ppar_i);
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tpar_i * b0xcv , Pperp_i);
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 2.* bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * Tpar_i * bracket(B0, 7. * Ppar_i, bm_exb);
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * Tpar_i * bracket(B0, Pperp_i, bm_exb);
                  break;
               case 4:
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 4.* bxkappa * Grad(Ti0 * ni);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tpar_i * bxkappa * 7., Ppar_i);
                  ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tpar_i * bxkappa , Pperp_i);
                  break;
               case 5:
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 2.* bxkappa * Grad(Ti0 * ni);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tpar_i * bxkappa * 7., Ppar_i);
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tpar_i * bxkappa , Pperp_i);
                  // ddt(Ppar_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  ddt(Ppar_i) += Tpar_i/(Cnor*B0) * 2.* bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Ppar_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * Tpar_i * bracket(B0, 7. * Ppar_i, bm_exb);
                  ddt(Ppar_i) -= .5 /(Cnor*B0) * Tpar_i * bracket(B0, Pperp_i, bm_exb);
                  break;
            }
         }
      }
   }

   if(nonlinear && nonlinear_terms > 0){
      ddt(Ppar_i) -= bracket(gyrophi, Ppar_i, bm_exb);
      if(nonlinear_terms > 4)
         ddt(Ppar_i) -= nG0 * bracket(gyrophi2 - gyrophi, Tperp_i, bm_exb);  // Kim line 1
      // ddt(Ppar_i) -= 0.5*n0*bracket(gyrophi_a, Tperp_i, bm_exb);  
   } 

   if(Landau_damping_i){
      if(!nonlinear)
      {
         ddt(Ppar_i) -= 2.*Qperp_i/B0*Grad_par(B0, CELL_CENTRE);
         // ddt(Ppar_i) -= B0*Grad_parP_gyro(Qpar_i/B0);
         ddt(Ppar_i) -= B0 * Grad_par(Qpar_i / B0);
      }
      else
         ddt(Ppar_i) -= Grad_par(Qpar_i);
   }

   if(toroidal_closure2){
      switch(curv_model)
      {
         case 1:
         default:
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu1i * Pi0 * b0xcv, Tpar_i);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu2i * Pi0 * b0xcv, Tperp_i);
            break;
         case 2:
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 2. /(Cnor*B0) * Pi0 * bracket(nu1i * B0, Tpar_i, bm_exb);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 2. /(Cnor*B0) * Pi0 * bracket(nu2i * B0, Tperp_i, bm_exb);
            break;
         case 3:
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu1i * Pi0 * b0xcv, Tpar_i);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu2i * Pi0 * b0xcv, Tperp_i);
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu1i * B0, Tpar_i, bm_exb);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu2i * B0, Tperp_i, bm_exb);
            break;
         case 4:
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu1i * Pi0 * bxkappa, Tpar_i);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu2i * Pi0 * bxkappa, Tperp_i);
            break;
         case 5:
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu1i * Pi0 * bxkappa, Tpar_i);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu2i * Pi0 * bxkappa, Tperp_i);
            // ddt(Ppar_i) -= nu1i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Ppar_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu1i * B0, Tpar_i, bm_exb);
            // ddt(Ppar_i) -= nu2i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Ppar_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu2i * B0, Tperp_i, bm_exb);
            break;
      }
   }

   if(toroidal_closure3){
      /****/
      ddt(Ppar_i) -= 2. * n0 * (nu1r * mod_wd_Tpar + nu2r * mod_wd_Tperp);
      /***
        toroidal1=1*nu1r*n0*(mod_wd(Tpar_i,&vgb,k0_toroidal)+mod_wd(Tpar_i,&vcurv,k0_toroidal));
        toroidal2=1*nu2r*n0*(mod_wd(Tperp_i,&vgb,k0_toroidal)+mod_wd(Tperp_i,&vcurv,k0_toroidal));
        toroidal1.applyBoundary();
        toroidal2.applyBoundary();
        mesh->communicate(toroidal1);
        mesh->communicate(toroidal2);
        ddt(Ppar_i) -= toroidal1+toroidal2;  
       ***/
   }

   if(sink_Ppar_il > 0.0){
      ddt(Ppar_i) -=  sink_Ppar_il*sink_tanhxl(P0,Ppar_i,sPpar_i_widthl,sPpar_i_lengthl); // core sink
   }

   if(sink_Ppar_ir > 0.0){
      ddt(Ppar_i) -=  sink_Ppar_ir*sink_tanhxr(P0,Ppar_i,sPpar_i_widthr,sPpar_i_lengthr); //  sol sink
   }

   /*
      if(diffusion_Ppar_i4 > 0.0){
      Grad2_Ppar_i = Grad2_par2new(Ppar_i);
      Grad2_Ppar_i.applyBoundary();
      mesh->communicate(Grad2_Ppar_i);
      ddt(Ppar_i) -= diffusion_Ppar_i4 * Grad2_par2new(Grad2_Ppar_i);
      } */

   ////////////////////////////////////////////////////
   // perpendicular ion pressure
   ddt(Pperp_i) = 0.;

   // ddt(Pperp_i) -= bracket(phi, P0 / 2., bm_exb);

   ddt(Pperp_i) -= Ti0 * bracket(gyrophi2, nG0, bm_exb);
   ddt(Pperp_i) -= nG0 * bracket(gyrophi2, Ti0, bm_exb);

   if(FLR_effect)
   {
      // ddt(Pperp_i) -= 0.5 * T0 * bracket(gyrophi_a, nG0, bm_exb);
      // ddt(Pperp_i) -= 0.5 * nG0  * bracket(gyrophi_a, T0, bm_exb);
      ddt(Pperp_i) -= 2 * nG0 * bracket(gyrophi3 - gyrophi2, Ti0, bm_exb);        // Kim line 1
   }

   if(continuity)
   {
      if(isotropic)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Pperp_i) -= 10. / 3. * Pi0 / B0 * b0xcv * Grad(phi);
               break;
            case 2:
               ddt(Pperp_i) -= 10. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 3:
               ddt(Pperp_i) -= 5. / 3. * Pi0 / B0 * b0xcv * Grad(phi);
               ddt(Pperp_i) -= 5. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
            case 4:
               ddt(Pperp_i) -= 10. / 3. * Pi0 / B0 * bxkappa * Grad(phi);
               break;
            case 5:
               ddt(Pperp_i) -= 5. / 3. * Pi0 / B0 * bxkappa * Grad(phi);
               ddt(Pperp_i) -= 5. / 3. * Pi0 / B0 * bracket(B0, phi, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Pperp_i) -= 10. / 3. * Pperp_i / B0 * b0xcv * Grad(phi);
                  break;
               case 2:
                  ddt(Pperp_i) -= 10. / 3. * Pperp_i / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 3:
                  ddt(Pperp_i) -= 5. / 3. * Pperp_i / B0 * b0xcv * Grad(phi);
                  ddt(Pperp_i) -= 5. / 3. * Pperp_i / B0 * bracket(B0, phi, bm_exb);
                  break;
               case 4:
                  ddt(Pperp_i) -= 10. / 3. * Pperp_i / B0 * bxkappa * Grad(phi);
                  break;
               case 5:
                  ddt(Pperp_i) -= 5. / 3. * Pperp_i / B0 * bxkappa * Grad(phi);
                  ddt(Pperp_i) -= 5. / 3. * Pperp_i / B0 * bracket(B0, phi, bm_exb);
                  break;
            }
         }
      }
      else
      {
         // ddt(Pperp_i) -= 3. * n0 * T0 / B0 * b0xcv * Grad(gyrophi);
         // ddt(Pperp_i) -= 1.5 * n0 * T0 / B0 * b0xcv * Grad(gyrophi_a);   // Kim line 4
         // ddt(Pperp_i) -= 3. * n0 * T0 / B0 * b0xcv * Grad(gyrophi);
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Pperp_i) -= 2 * n0 * Ti0 / B0 * b0xcv * Grad(gyrophi3);   // Kim line 4
               ddt(Pperp_i) -= n0 * Ti0 / B0 * b0xcv * Grad(gyrophi2);
               break;
            case 2:
               ddt(Pperp_i) -= 2 * n0 * Ti0 / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
               ddt(Pperp_i) -= n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb);
               break;
            case 3:
               ddt(Pperp_i) -= 1. * n0 * Ti0 / B0 * b0xcv * Grad(gyrophi3);   // Kim line 4
               ddt(Pperp_i) -= 0.5 * n0 * Ti0 / B0 * b0xcv * Grad(gyrophi2);
               ddt(Pperp_i) -= 1. * n0 * Ti0 / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
               ddt(Pperp_i) -= 0.5 * n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb);
               break;
            case 4:
               ddt(Pperp_i) -= 2 * n0 * Ti0 / B0 * bxkappa * Grad(gyrophi3);   // Kim line 4
               ddt(Pperp_i) -= n0 * Ti0 / B0 * bxkappa * Grad(gyrophi2);
               break;
            case 5:
               ddt(Pperp_i) -= 1. * n0 * Ti0 / B0 * bxkappa * Grad(gyrophi3);   // Kim line 4
               ddt(Pperp_i) -= 0.5 * n0 * Ti0 / B0 * bxkappa * Grad(gyrophi2);
               ddt(Pperp_i) -= 1. * n0 * Ti0 / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
               ddt(Pperp_i) -= 0.5 * n0 * Ti0 / B0 * bracket(B0, gyrophi2, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 1)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Pperp_i) -= 2 * Pperp_i / B0 * b0xcv * Grad(gyrophi3);   // Kim line 4
                  ddt(Pperp_i) -= Pperp_i / B0 * b0xcv * Grad(gyrophi2);
                  break;
               case 2:
                  ddt(Pperp_i) -= 2 * Pperp_i / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
                  ddt(Pperp_i) -= Pperp_i / B0 * bracket(B0, gyrophi2, bm_exb);
                  break;
               case 3:
                  ddt(Pperp_i) -= 1. * Pperp_i / B0 * b0xcv * Grad(gyrophi3);   // Kim line 4
                  ddt(Pperp_i) -= 0.5 * Pperp_i / B0 * b0xcv * Grad(gyrophi2);
                  ddt(Pperp_i) -= 1. * Pperp_i / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
                  ddt(Pperp_i) -= 0.5 * Pperp_i / B0 * bracket(B0, gyrophi2, bm_exb);
                  break;
               case 4:
                  ddt(Pperp_i) -= 2 * Pperp_i / B0 * bxkappa * Grad(gyrophi3);   // Kim line 4
                  ddt(Pperp_i) -= Pperp_i / B0 * bxkappa * Grad(gyrophi2);
                  break;
               case 5:
                  ddt(Pperp_i) -= 1. * Pperp_i / B0 * bxkappa * Grad(gyrophi3);   // Kim line 4
                  ddt(Pperp_i) -= 0.5 * Pperp_i / B0 * bxkappa * Grad(gyrophi2);
                  ddt(Pperp_i) -= 1. * Pperp_i / B0 * bracket(B0, gyrophi3, bm_exb);   // Kim line 4
                  ddt(Pperp_i) -= 0.5 * Pperp_i / B0 * bracket(B0, gyrophi2, bm_exb);
                  break;
            }
         }
      }

      if(compression)
      {
         if(isotropic)
         {
            ddt(Pperp_i) -= 2. / 3. * Pi0 * B0 * Grad_parP(Vpar_i / B0, CELL_CENTRE);
            if(nonlinear && nonlinear_terms > 2)
               ddt(Pperp_i) -= 2. / 3. * Pperp_i * B0 * Grad_par(Vpar_i / B0, CELL_CENTRE);
         }
         else
         {
            ddt(Pperp_i) -= Pi0 * B0 * Grad_parP_gyro(Vpar_i / B0, CELL_CENTRE); // Kim line 2
            ddt(Pperp_i) += Pi0 * Vpar_i /B0 * Grad_par(B0, CELL_CENTRE);
            if(nonlinear && nonlinear_terms > 2)
            {
               ddt(Pperp_i) -= Pperp_i * B0 * Grad_par(Vpar_i / B0, CELL_CENTRE); // Kim line 2
               ddt(Pperp_i) += Pperp_i * Vpar_i /B0 * Grad_par(B0, CELL_CENTRE);
            }
         }
      }
   }

   if(energy_flux)
   {
      if(isotropic)
      {
         switch(curv_model)
         {
            case 1:
            default:
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pi0 * b0xcv, Tpar_i + Tperp_i);
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * b0xcv * Grad(Ti0);
               break;
            case 2:
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
            case 3:
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pi0 * b0xcv, Tpar_i + Tperp_i);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * b0xcv * Grad(Ti0);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
            case 4:
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pi0 * bxkappa, Tpar_i + Tperp_i);
               ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bxkappa * Grad(Ti0);
               break;
            case 5:
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pi0 * bxkappa, Tpar_i + Tperp_i);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bxkappa * Grad(Ti0);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * Pi0 * bracket(B0, Tpar_i + Tperp_i, bm_exb);
               ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * n0 * (Tpar_i + Tperp_i) * bracket(B0, Ti0, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 3)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pperp_i * b0xcv, Tpar_i + Tperp_i);
                  break;
               case 2:
                  ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * Pperp_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
               case 3:
                  ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pperp_i * b0xcv, Tpar_i + Tperp_i);
                  ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * Pperp_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
               case 4:
                  ddt(Pperp_i) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(Pperp_i * bxkappa, Tpar_i + Tperp_i);
                  break;
               case 5:
                  ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * V_dot_Grad(Pperp_i * bxkappa, Tpar_i + Tperp_i);
                  ddt(Pperp_i) -= 5. / 6. / (Cnor * B0) * Pperp_i * bracket(B0, Tpar_i + Tperp_i, bm_exb);
                  break;
            }
         }
      }
      else
      {
         switch(curv_model)
         {
            case 1:
            default:
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += Ti0 / (Cnor*B0) * 3. * b0xcv * Grad(Ti0 * ni);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv , Ppar_i);
               ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(5. * Ti0 * b0xcv , Pperp_i);
               break;
            case 2:
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += Ti0 / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= 1. /(Cnor*B0) * Ti0 * bracket(B0 , Ppar_i, bm_exb);
               ddt(Pperp_i) -= 1. /(Cnor*B0) * 5. * Ti0 * bracket(B0 , Pperp_i, bm_exb);
               break;
            case 3:
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += .5 * Ti0 / (Cnor*B0) * 3. * b0xcv * Grad(Ti0 * ni);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * b0xcv , Ppar_i);
               ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(5. * Ti0 * b0xcv , Pperp_i);
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += .5 * Ti0 / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0 , Ppar_i, bm_exb);
               ddt(Pperp_i) -= .5 /(Cnor*B0) * 5. * Ti0 * bracket(B0 , Pperp_i, bm_exb);
               break;
            case 4:
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += Ti0 / (Cnor*B0) * 3. * bxkappa * Grad(Ti0 * ni);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa , Ppar_i);
               ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(5. * Ti0 * bxkappa , Pperp_i);
               break;
            case 5:
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += .5 * Ti0 / (Cnor*B0) * 3. * bxkappa * Grad(Ti0 * ni);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(Ti0 * bxkappa , Ppar_i);
               ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(5. * Ti0 * bxkappa , Pperp_i);
               // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
               ddt(Pperp_i) += .5 * Ti0 / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
               // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
               ddt(Pperp_i) -= .5 /(Cnor*B0) * Ti0 * bracket(B0 , Ppar_i, bm_exb);
               ddt(Pperp_i) -= .5 /(Cnor*B0) * 5. * Ti0 * bracket(B0 , Pperp_i, bm_exb);
               break;
         }
         if(nonlinear && nonlinear_terms > 3)
         {
            switch(curv_model)
            {
               case 1:
               default:
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += Tperp_i / (Cnor*B0) * 3. * b0xcv * Grad(Ti0 * ni);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tperp_i * b0xcv , Ppar_i);
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(5. * Tperp_i * b0xcv , Pperp_i);
                  break;
               case 2:
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += Tperp_i / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * Tperp_i * bracket(B0 , Ppar_i, bm_exb);
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * 5. * Tperp_i * bracket(B0 , Pperp_i, bm_exb);
                  break;
               case 3:
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += .5 * Tperp_i / (Cnor*B0) * 3. * b0xcv * Grad(Ti0 * ni);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tperp_i * b0xcv , Ppar_i);
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(5. * Tperp_i * b0xcv , Pperp_i);
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += .5 * Tperp_i / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * Tperp_i * bracket(B0 , Ppar_i, bm_exb);
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * 5. * Tperp_i * bracket(B0 , Pperp_i, bm_exb);
                  break;
               case 4:
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += Tperp_i / (Cnor*B0) * 3. * bxkappa * Grad(Ti0 * ni);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(Tperp_i * bxkappa , Ppar_i);
                  ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(5. * Tperp_i * bxkappa , Pperp_i);
                  break;
               case 5:
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += .5 * Tperp_i / (Cnor*B0) * 3. * bxkappa * Grad(Ti0 * ni);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(Tperp_i * bxkappa , Ppar_i);
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * V_dot_Grad(5. * Tperp_i * bxkappa , Pperp_i);
                  // ddt(Pperp_i) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  ddt(Pperp_i) += .5 * Tperp_i / (Cnor*B0) * 3. * bracket(B0, Ti0 * ni, bm_exb);
                  // ddt(Pperp_i) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * Tperp_i * bracket(B0 , Ppar_i, bm_exb);
                  ddt(Pperp_i) -= .5 /(Cnor*B0) * 5. * Tperp_i * bracket(B0 , Pperp_i, bm_exb);
                  break;
            }
         }
      }
   }

   if(nonlinear && nonlinear_terms > 0){
      ddt(Pperp_i) -= bracket(gyrophi, Pperp_i, bm_exb);
      if(nonlinear_terms > 4)
         ddt(Pperp_i) -= 2 * nG0 * bracket(gyrophi3 - gyrophi2, Tperp_i, bm_exb);        // Kim line 1
      // ddt(Pperp_i) -= 0.5*bracket(gyrophi_a, Pperp_i, bm_exb);
      // ddt(Pperp_i) += 0.5*Ti0*n0*bracket(gyroPsi_a*B0, Vpar_i, bm_mag);    
   }

   if(Landau_damping_i)
   {
      if(!nonlinear)
      {
         // ddt(Pperp_i) -= B0*B0*Grad_parP_gyro(Qperp_i/(B0*B0));
         ddt(Pperp_i) -= B0 * B0 * Grad_par(Qperp_i / (B0 * B0));
      }
      else
         ddt(Pperp_i) -= Grad_par(Qperp_i);
   }

   if(toroidal_closure2){
      switch(curv_model)
      {
         case 1:
         default:
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu3i * Pi0 * b0xcv, Tpar_i);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu4i * Pi0 * b0xcv, Tperp_i);
            break;
         case 2:
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 2. /(Cnor*B0) * Pi0 * bracket(nu3i * B0, Tpar_i, bm_exb);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 2. /(Cnor*B0) * Pi0 * bracket(nu4i * B0, Tperp_i, bm_exb);
            break;
         case 3:
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu3i * Pi0 * b0xcv, Tpar_i);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu4i * Pi0 * b0xcv, Tperp_i);
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu3i * B0, Tpar_i, bm_exb);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu4i * B0, Tperp_i, bm_exb);
            break;
         case 4:
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu3i * Pi0 * bxkappa, Tpar_i);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 2. /(Cnor*B0) * V_dot_Grad(nu4i * Pi0 * bxkappa, Tperp_i);
            break;
         case 5:
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu3i * Pi0 * bxkappa, Tpar_i);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 1. /(Cnor*B0) * V_dot_Grad(nu4i * Pi0 * bxkappa, Tperp_i);
            // ddt(Pperp_i) -= nu3i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,Tpar_i);
            ddt(Pperp_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu3i * B0, Tpar_i, bm_exb);
            // ddt(Pperp_i) -= nu4i*n0*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Tperp_i);   
            ddt(Pperp_i) -= 1. /(Cnor*B0) * Pi0 * bracket(nu4i * B0, Tperp_i, bm_exb);
            break;
      }
   }

   if(toroidal_closure3){
      /***/
      ddt(Pperp_i) -= 2.0 * n0 * (nu3r * mod_wd_Tpar + nu4r * mod_wd_Tperp);
      /***
        toroidal3=1*nu3r*n0*(mod_wd(Tpar_i,&vgb,k0_toroidal)+mod_wd(Tpar_i,&vcurv,k0_toroidal));
        toroidal4=1*nu4r*n0*(mod_wd(Tperp_i,&vgb,k0_toroidal)+mod_wd(Tperp_i,&vcurv,k0_toroidal));
        toroidal3.applyBoundary();
        toroidal4.applyBoundary();
        mesh->communicate(toroidal3);
        mesh->communicate(toroidal4);
        ddt(Pperp_i) -= toroidal3+toroidal4;  
       ***/
   } 

   if(sink_Pperp_il > 0.0){
      ddt(Pperp_i) -=  sink_Pperp_il*sink_tanhxl(P0,Pperp_i,sPperp_i_widthl,sPperp_i_lengthl); // core sink
   }

   if(sink_Pperp_ir > 0.0){
      ddt(Pperp_i) -=  sink_Pperp_ir*sink_tanhxr(P0,Pperp_i,sPperp_i_widthr,sPperp_i_lengthr); //  sol sink
   }

   ///////////////////////////////////////////////////
   // vorticity
   if(evolve_ne)
   {
      ddt(ne) = -bracket(phi, ne, bm_exb);
      ddt(ne) += 1./(Cnor*B0)*V_dot_Grad(b0xcv, Ppar_e);
      ddt(ne) += 1./(Cnor*B0)*V_dot_Grad(b0xcv, Pperp_e);
      ddt(ne) -= Va * Va / Vbar / Vbar * B0 / Cnor * Grad_parP(Jpar, CELL_CENTRE);
   }
   else
   {


      if(!electrostatic)
      {
         ddt(U) = 0.0;
         ddt(U) -= bracket(phi, U0, bm_exb);

         if(GradparJ_test)
         {
            Vector2D Jparpara;
            Jparpara.covariant = false;
            Jparpara.x = 0;
            Jparpara.y = -Va * Va / Vbar / Vbar;
            Jparpara.z = 0;
            GradparJ_U = B0 * B0 / sqrt(mesh->g_22) * V_dot_Grad(Jparpara, Jpar);
            GradparJ_C = -Va * Va / Vbar / Vbar * B0 * B0 * Grad_parP(Jpar, CELL_CENTRE);
         }

         if(shear_Alfven_wave)
         {
            if(kpar < 0)
            {
               if(shear_Alfven_advection)
                  ddt(U) -= Vpar_Grad_par(Va * Va / Vbar / Vbar * B0 * B0, Jpar, CELL_CENTRE);
               else
                  ddt(U) -= Va * Va / Vbar / Vbar * B0 * B0 * Grad_parP(Jpar, CELL_CENTRE);
            }
            else
               ddt(U) += Va * Va / Vbar / Vbar * B0 * B0 * Jpar / kpar;
            ddt(U) -= B0 * B0 * B0 * bracket(Psi, J0, bm_mag);
         }

         switch(curv_model)
         {
            case 1:
            default:
               ddt(U) += b0xcv * Grad(P);
               break;
            case 2:
               ddt(U) += bracket(B0, P, bm_exb);
               break;
            case 3:
               ddt(U) += .5 * b0xcv * Grad(P);
               ddt(U) += .5 * bracket(B0, P, bm_exb);
               break;
            case 4:
               ddt(U) += bxkappa * Grad(P);
               break;
            case 5:
               ddt(U) += .5 * bxkappa * Grad(P);
               ddt(U) += .5 * bracket(B0, P, bm_exb);
               break;
         }



         if(FLR_effect)
         {
            ddt(U) -= Cnor * B0 * bracket(-phi_f, nG0, bm_exb);
            switch(curv_model)
            {
               case 1:
               default:
                  ddt(U) += 2. * Cnor * n0 * b0xcv * Grad(phi_f);
                  ddt(U) += 0.5 * Cnor * n0 * b0xcv * Grad(gyrophi_a);
                  break;
               case 2:
                  ddt(U) += 2. * Cnor * n0 * bracket(B0, phi_f, bm_exb);
                  ddt(U) += 0.5 * Cnor * n0 * bracket(B0, gyrophi_a, bm_exb);
                  break;
               case 3:
                  ddt(U) += 1. * Cnor * n0 * b0xcv * Grad(phi_f);
                  ddt(U) += 0.25 * Cnor * n0 * b0xcv * Grad(gyrophi_a);
                  ddt(U) += 1. * Cnor * n0 * bracket(B0, phi_f, bm_exb);
                  ddt(U) += 0.25 * Cnor * n0 * bracket(B0, gyrophi_a, bm_exb);
                  break;
               case 4:
                  ddt(U) += 2. * Cnor * n0 * bxkappa * Grad(phi_f);
                  ddt(U) += 0.5 * Cnor * n0 * bxkappa * Grad(gyrophi_a);
                  break;
               case 5:
                  ddt(U) += 1. * Cnor * n0 * bxkappa * Grad(phi_f);
                  ddt(U) += 0.25 * Cnor * n0 * bxkappa * Grad(gyrophi_a);
                  ddt(U) += 1. * Cnor * n0 * bracket(B0, phi_f, bm_exb);
                  ddt(U) += 0.25 * Cnor * n0 * bracket(B0, gyrophi_a, bm_exb);
                  break;
            }
            ddt(U) -= 0.5 * Cnor * n0 / (Ti0 * B0) * bracket(-gyrophi_a, Ti0, bm_exb);

            if(false && compression)
               ddt(U) -= Cnor * B0 * B0 * Grad_parP(n0 * (gyroVpar_i - Vpar_i) / B0, CELL_CENTRE);
         }

         if(ddtU_Terms_test)
         {
            ddtU_phiU0 = 0;
            ddtU_phi_fnG0 = 0;
            ddtU_GradJpar = 0;
            ddtU_PsiJ0 = 0;
            ddtU_GradP = 0;
            ddtU_Vpar = 0;
            ddtU_B0phiT0 = 0;

            ddtU_phiU0 -= bracket(phi, U0, bm_exb);
            ddtU_phi_fnG0 += Cnor * B0 * bracket(phi_f, nG0, bm_exb);
            ddtU_GradJpar -= Va * Va / Vbar / Vbar * B0 * B0 * Grad_parP(Jpar, CELL_CENTRE);
            ddtU_PsiJ0 -= Cnor * B0 * B0 * B0 * bracket(Psi, J0 / B0, bm_mag);
            ddtU_GradP += b0xcv * Grad(Ppar_e + Ppar_i + Pperp_e + Pperp_i);
            ddtU_Vpar -= Cnor * B0 * B0 * Grad_parP(n0 * (gyroVpar_i - Vpar_i) / B0, CELL_CENTRE);
            ddtU_B0phiT0 += 0.5 * Cnor * n0 * b0xcv * Grad(gyrophi_a);
            ddtU_B0phiT0 += 0.5 * Cnor * n0 / (T0 * B0) * bracket(gyrophi_a, T0, bm_exb);

            ddtU_phiU0.applyBoundary();
            ddtU_phi_fnG0.applyBoundary();
            ddtU_GradJpar.applyBoundary();
            ddtU_PsiJ0.applyBoundary();
            ddtU_GradP.applyBoundary();
            ddtU_Vpar.applyBoundary();
            ddtU_B0phiT0.applyBoundary();
         }

         if(nonlinear && nonlinear_terms > 0) {
            ddt(U) -= bracket(phi, U, bm_exb);
            if(nonlinear_terms > 4)
            {
               ddt(U) -= Cnor*B0*bracket(-phi_f, ni, bm_exb);
               if(compression)
                  ddt(U) -= Cnor*B0*B0*B0*bracket(gyroPsi-Psi, n0*Vpar_i/B0, bm_mag);  
               ddt(U) -= 0.5*Cnor*B0*bracket(-gyrophi_a, n0 / T0 * Tperp_i, bm_exb);
            }
         }  

         if(viscos_par > 0.0)
            ddt(U) += viscos_par*Grad2_par2(U); 

         if(viscos_perp > 0.0)
            ddt(U) += viscos_perp * Delp2(U);

         // Hyper-viscosity
         if(hyperviscos > 0.0) {
            // Calculate coefficient.

            hyper_mu_x = hyperviscos * mesh->g_11*SQ(mesh->dx) * abs(mesh->g11*D2DX2(U)) / (abs(U) + 1e-3);
            hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries

            ddt(U) -= hyperviscos * Delp2(Delp2(U));

            if(first_run) { // Print out maximum values of viscosity used on this processor
               output.write("   Hyper-viscosity values:\n");
               output.write("      Max mu_x = %e, Max_DC mu_x = %e\n", max(hyper_mu_x), max(hyper_mu_x.DC()));
            }
         }

         if(sink_Ul > 0.0){
            ddt(U) -=  sink_Ul*sink_tanhxl(P0,U,su_widthl,su_lengthl); // core sink
         }

         if(sink_Ur > 0.0){
            ddt(U) -=  sink_Ur*sink_tanhxr(P0,U,su_widthr,su_lengthr); //  sol sink
         }
      }
   }

   //////////////////////////////////////////////////
   //Ohm's Law

   if(!electrostatic)
   {
      ddt(Psi) = 0.0; 
      ddt(Psi) += eta * Jpar;
      if(kpar < 0.)
      {
         if(shear_Alfven_advection)
            ddt(Psi) -= Vpar_Grad_par(1. / (B0 * beta_coeff), phi, CELL_CENTRE); 
         else
            ddt(Psi) -= 1. / (B0 * beta_coeff) * Grad_parP(phi, CELL_CENTRE); 
      }
      else
         ddt(Psi) -= 1. / (B0 * beta_coeff * kpar) * phi; 
      // ddt(Psi) -= T0 / (Cnor * B0 * n0) * bracket(B0 * Psi, n0, bm_mag);
      // ddt(Psi) -= 1. / (Cnor * B0) * bracket(B0 * Psi, T0, bm_mag);

      if(false && compression)
      {
         ddt(Psi) -= 1. / Cnor * (Tpar_e - Tperp_e) / (B0 * B0) * Grad_parP(B0, CELL_CENTRE);
      }

      if(eHall)
      {
         if(eHall_adiabatic)
         {
            if(kpar < 0.)
            {
               if(shear_Alfven_advection)
                  ddt(Psi) -= Vpar_Grad_par(-1. / (Cnor * n0 * beta_coeff * B0) * Te0, ne, CELL_YLOW);  
               else
                  ddt(Psi) += 1. / (Cnor * n0 * beta_coeff) * Te0 * Grad_parP(ne, CELL_YLOW) / B0;  
            }
            else
               ddt(Psi) += 1. / (Cnor * n0 * beta_coeff * kpar) * Te0 * ne / B0;
            ddt(Psi) -= 1. / (Cnor * n0) * bracket(Psi, Pe0, bm_mag);
         }
         else
         {
            if(kpar < 0.)
            {
               if(shear_Alfven_advection)
                  ddt(Psi) -= Vpar_Grad_par(-1. / (Cnor * n0 * beta_coeff * B0), Ppar_e, CELL_YLOW);  
               else
                  ddt(Psi) += 1. / (Cnor * n0 * beta_coeff) * Grad_parP(Ppar_e, CELL_YLOW) / B0;  
            }
            else
               ddt(Psi) += 1. / (Cnor * n0 * beta_coeff * kpar) * Ppar_e / B0;  
            ddt(Psi) -= 1. / (Cnor * n0) * bracket(Psi, Pe0, bm_mag);  
         }
      }

      if(hyperresist > 0.0)
         ddt(Psi) -= eta*hyperresist*Delp2(Jpar);
      if(ddtPsi_Terms_test)
      {
         ddtPsi_Gradphi = 0;
         ddtPsi_Jpar = 0;
         ddtPsi_Psin0 = 0;
         ddtPsi_PsiT0 = 0;
         ddtPsi_eHall = 0;
         ddtPsi_hyperresist = 0;

         ddtPsi_Gradphi = -1./B0*Grad_parP(phi, CELL_CENTRE);
         ddtPsi_Jpar = eta*Jpar; 
         ddtPsi_Psin0 -= T0/(Cnor*B0*n0)*bracket(B0*Psi, n0, bm_mag);
         ddtPsi_PsiT0 -= 1./(Cnor*B0)*bracket(B0*Psi, T0, bm_mag);

         if(eHall)
            ddtPsi_eHall += 1./(Cnor*n0)*Grad_parP(T0*ne/B0, CELL_YLOW);  

         if(hyperresist > 0.0)
            ddtPsi_hyperresist -= eta*hyperresist*Delp2(Jpar);

         ddtPsi_Gradphi.applyBoundary();
         ddtPsi_Jpar.applyBoundary();
         ddtPsi_Psin0.applyBoundary();
         ddtPsi_PsiT0.applyBoundary();
         ddtPsi_eHall.applyBoundary();
         ddtPsi_hyperresist.applyBoundary();
      }
   }

   // if(sink_Ur > 0.0){
   //   ddt(Psi) -=  sink_Ur*sink_tanhxr(P0,Psi,su_widthr,su_lengthr); //  sol sink
   // }

   //////////////////////////////////////////////////
   //parallel electron pressure 
   if(!Zero_Te){


      if(evolve_Te)
      {
         ddt(Te) = 0.0;
         
         ddt(Te) -= bracket(phi, Te0, bm_exb);

         if(continuity)
         {
            ddt(Te) -= 4. / 3. * Te0 / B0 * b0xcv * Grad(phi);
            if(nonlinear)
            {
               ddt(Te) -= 4. / 3. * Te / B0 * b0xcv * Grad(phi);
            }
         }

         if(compression)
         {
            ddt(Te) -= 2. / 3. * Te0 * B0 * Grad_parP(Vpar_e/B0, CELL_CENTRE);
            if(nonlinear)
               ddt(Te) -= 2. / 3. * Te * B0 * Grad_par(Vpar_e/B0, CELL_CENTRE);
         }

         if(nonlinear)
            ddt(Te) -= bracket(phi, Te, bm_exb);

         if(Landau_damping_e)
            ddt(Te) -= Grad_par(Qpar_e) * 2. / (3. * n0);

         if(energy_flux)
         {
            ddt(Te) += 4. / 3. / (Cnor * B0) * Te0 / n0 * b0xcv * Grad(Ppar_e);
            ddt(Te) -= 10. / 3. / (Cnor * B0) * V_dot_Grad(-Te0 * b0xcv, Te);
            ddt(Te) += 10. / 3. / (Cnor * B0) * Te * b0xcv * Grad(Te0);
            if(nonlinear)
            {
               ddt(Te) += 4. / 3. / (Cnor * B0) * Te / n0 * b0xcv * Grad(Ppar_e);
               ddt(Te) -= 10. / 3. / (Cnor * B0) * V_dot_Grad(-Te * b0xcv, Te);
            }
         }
      }
      else
      {
         if(!electrostatic)
         {
            ddt(Ppar_e) = 0.0;

            ddt(Ppar_e) -= Te0 * bracket(phi, n0, bm_exb);
            ddt(Ppar_e) -= n0 * bracket(phi, Te0, bm_exb);

            if(continuity)
            {
               if(isotropic) 
               {
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Ppar_e) -= 10. / 3. * Pe0 / B0 * b0xcv * Grad(phi);
                        break;
                     case 2:
                        ddt(Ppar_e) -= 10. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 3:
                        ddt(Ppar_e) -= 5. / 3. * Pe0 / B0 * b0xcv * Grad(phi);
                        ddt(Ppar_e) -= 5. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 4:
                        ddt(Ppar_e) -= 10. / 3. * Pe0 / B0 * bxkappa * Grad(phi);
                        break;
                     case 5:
                        ddt(Ppar_e) -= 5. / 3. * Pe0 / B0 * bxkappa * Grad(phi);
                        ddt(Ppar_e) -= 5. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 1)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Ppar_e) -= 10. / 3. * Ppar_e / B0 * b0xcv * Grad(phi);
                           break;
                        case 2:
                           ddt(Ppar_e) -= 10. / 3. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 3:
                           ddt(Ppar_e) -= 5. / 3. * Ppar_e / B0 * b0xcv * Grad(phi);
                           ddt(Ppar_e) -= 5. / 3. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 4:
                           ddt(Ppar_e) -= 10. / 3. * Ppar_e / B0 * bxkappa * Grad(phi);
                           break;
                        case 5:
                           ddt(Ppar_e) -= 5. / 3. * Ppar_e / B0 * bxkappa * Grad(phi);
                           ddt(Ppar_e) -= 5. / 3. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                     }
                  }
               }
               else
               {
                  // ddt(Ppar_e) -= 4. * n0 * Te0 / B0 * b0xcv * Grad(phi);  
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Ppar_e) -= 4. * Pe0 / B0 * b0xcv * Grad(phi);
                        break;
                     case 2:
                        ddt(Ppar_e) -= 4. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 3:
                        ddt(Ppar_e) -= 2. * Pe0 / B0 * b0xcv * Grad(phi);
                        ddt(Ppar_e) -= 2. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 4:
                        ddt(Ppar_e) -= 4. * Pe0 / B0 * bxkappa * Grad(phi);
                        break;
                     case 5:
                        ddt(Ppar_e) -= 2. * Pe0 / B0 * bxkappa * Grad(phi);
                        ddt(Ppar_e) -= 2. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 1)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Ppar_e) -= 4. * Ppar_e / B0 * b0xcv * Grad(phi);
                           break;
                        case 2:
                           ddt(Ppar_e) -= 4. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 3:
                           ddt(Ppar_e) -= 2. * Ppar_e / B0 * b0xcv * Grad(phi);
                           ddt(Ppar_e) -= 2. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 4:
                           ddt(Ppar_e) -= 4. * Ppar_e / B0 * bxkappa * Grad(phi);
                           break;
                        case 5:
                           ddt(Ppar_e) -= 2. * Ppar_e / B0 * bxkappa * Grad(phi);
                           ddt(Ppar_e) -= 2. * Ppar_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                     }
                  }
               }
            }

            if(compression)
            {
               if(true || isotropic)
               {
                  ddt(Ppar_e) -= 2. / 3. * Pe0 * B0 * Grad_parP(Vpar_e / B0, CELL_CENTRE);
                  if(nonlinear && nonlinear_terms > 2)
                     ddt(Ppar_e) -= 2. / 3. * Ppar_e * B0 * Grad_par(Vpar_e / B0, CELL_CENTRE);
               }
               else
               {
                  ddt(Ppar_e) -= 5. / 3. * B0 * Te0 * n0 * Grad_parP(Vpar_e / B0, CELL_YLOW);
                  if(nonlinear && nonlinear_terms > 2)
                     ddt(Ppar_e) -= 5. / 3. * B0 * Ppar_e * Grad_par(Vpar_e / B0, CELL_YLOW);
               }
            }

            if(nonlinear && nonlinear_terms > 0){
               ddt(Ppar_e) -= bracket(phi, Ppar_e, bm_exb);
            }

            if(Landau_damping_e)
            {
               // ddt(Ppar_e) -= B0*Grad_parP(Qpar_e/B0);  
               ddt(Ppar_e) -= Grad_par(Qpar_e);  
            }

            if(energy_flux)
            {
               if(isotropic)
               {
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * b0xcv, Tpar_e + Tperp_e);
                        ddt(Ppar_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * b0xcv * Grad(Te0);
                        break;
                     case 2:
                        ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Ppar_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                     case 3:
                        ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * b0xcv, Tpar_e + Tperp_e);
                        ddt(Ppar_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * b0xcv * Grad(Te0);
                        ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Ppar_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                     case 4:
                        ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * bxkappa, Tpar_e + Tperp_e);
                        ddt(Ppar_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bxkappa * Grad(Te0);
                        break;
                     case 5:
                        ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * bxkappa, Tpar_e + Tperp_e);
                        ddt(Ppar_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bxkappa * Grad(Te0);
                        ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Ppar_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 3)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Ppar_e * b0xcv, Tpar_e + Tperp_e);
                           break;
                        case 2:
                           ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * Ppar_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                        case 3:
                           ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Ppar_e * b0xcv, Tpar_e + Tperp_e);
                           ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * Ppar_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                        case 4:
                           ddt(Ppar_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Ppar_e * bxkappa, Tpar_e + Tperp_e);
                           break;
                        case 5:
                           ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Ppar_e * bxkappa, Tpar_e + Tperp_e);
                           ddt(Ppar_e) -= .5 * 5. / 3. / (Cnor * B0) * Ppar_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                     }
                  }
               }
               else
               {
                  // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ni);
                  // ddt(Ppar_e) -= Te0/(Cnor*B0)*b0xcv*Grad(4.*Te0*ne);
                  // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_i+Pperp_i);    
                  // ddt(Ppar_e) += Te0/(Cnor*B0)*b0xcv*Grad(7.*Ppar_e+Pperp_e);
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 4.* b0xcv * Grad(Te0 * ne);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv * 7., Ppar_e);
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv , Pperp_e);
                        break;
                     case 2:
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 4.* bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * Te0 * bracket(-B0, 7. * Ppar_e, bm_exb);
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * Te0 * bracket(-B0, Pperp_e, bm_exb);
                        break;
                     case 3:
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 2.* b0xcv * Grad(Te0 * ne);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv * 7., Ppar_e);
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv , Pperp_e);
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 2.* bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0, 7. * Ppar_e, bm_exb);
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0, Pperp_e, bm_exb);
                        break;
                     case 4:
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 4.* bxkappa * Grad(Te0 * ne);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa * 7., Ppar_e);
                        ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa , Pperp_e);
                        break;
                     case 5:
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 2.* bxkappa * Grad(Te0 * ne);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa * 7., Ppar_e);
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa , Pperp_e);
                        // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                        ddt(Ppar_e) -= Te0/(Cnor*B0) * 2.* bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0, 7. * Ppar_e, bm_exb);
                        ddt(Ppar_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0, Pperp_e, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 3)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 4.* b0xcv * Grad(Te0 * ne);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tpar_e * b0xcv * 7., Ppar_e);
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tpar_e * b0xcv , Pperp_e);
                           break;
                        case 2:
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 4.* bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * Tpar_e * bracket(-B0, 7. * Ppar_e, bm_exb);
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * Tpar_e * bracket(-B0, Pperp_e, bm_exb);
                           break;
                        case 3:
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 2.* b0xcv * Grad(Te0 * ne);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tpar_e * b0xcv * 7., Ppar_e);
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tpar_e * b0xcv , Pperp_e);
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 2.* bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * Tpar_e * bracket(-B0, 7. * Ppar_e, bm_exb);
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * Tpar_e * bracket(-B0, Pperp_e, bm_exb);
                           break;
                        case 4:
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 4.* bxkappa * Grad(Te0 * ne);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tpar_e * bxkappa * 7., Ppar_e);
                           ddt(Ppar_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tpar_e * bxkappa , Pperp_e);
                           break;
                        case 5:
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 2.* bxkappa * Grad(Te0 * ne);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tpar_e * bxkappa * 7., Ppar_e);
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tpar_e * bxkappa , Pperp_e);
                           // ddt(Ppar_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,4.*T0*ne);
                           ddt(Ppar_e) -= Tpar_e/(Cnor*B0) * 2.* bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Ppar_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, 7.*Ppar_e+Pperp_e);    
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * Tpar_e * bracket(-B0, 7. * Ppar_e, bm_exb);
                           ddt(Ppar_e) -= .5 /(Cnor*B0) * Tpar_e * bracket(-B0, Pperp_e, bm_exb);
                           break;
                     }
                  }
               }
            }

            if(viscos_Ppar_e > 0.0)
               ddt(Ppar_e) += viscos_Ppar_e*Grad2_par2(Ppar_e); 

            if(diffusion4_Ppar_e > 0.0){
               Grad2_Ppar_e = 0.0;
               Grad2_Ppar_e = Grad2_par2new(Ppar_e);
               Grad2_Ppar_e.applyBoundary();
               mesh->communicate(Grad2_Ppar_e);
               ddt(Ppar_e) -= diffusion4_Ppar_e * Grad2_par2new(Grad2_Ppar_e);
            }  

            if(sink_Ppar_el > 0.0){
               ddt(Ppar_e) -=  sink_Ppar_el*sink_tanhxl(P0,Ppar_e,sPpar_e_widthl,sPpar_e_lengthl); // core sink
            }

            if(sink_Ppar_er > 0.0){
               ddt(Ppar_e) -=  sink_Ppar_er*sink_tanhxr(P0,Ppar_e,sPpar_e_widthr,sPpar_e_lengthr); //  sol sink
            }

         }

         ///////////////////////////////////////////////////
         // perpendicular electron pressure

         if(!electrostatic)
         {
            ddt(Pperp_e) = 0.0;

            ddt(Pperp_e) -= Te0 * bracket(phi, n0, bm_exb);
            ddt(Pperp_e) -= n0 * bracket(phi, Te0, bm_exb);

            if(continuity)
            {
               if(isotropic)
               {
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Pperp_e) -= 10. / 3. * Pe0 / B0 * b0xcv * Grad(phi);
                        break;
                     case 2:
                        ddt(Pperp_e) -= 10. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 3:
                        ddt(Pperp_e) -= 5. / 3. * Pe0 / B0 * b0xcv * Grad(phi);
                        ddt(Pperp_e) -= 5. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 4:
                        ddt(Pperp_e) -= 10. / 3. * Pe0 / B0 * bxkappa * Grad(phi);
                        break;
                     case 5:
                        ddt(Pperp_e) -= 5. / 3. * Pe0 / B0 * bxkappa * Grad(phi);
                        ddt(Pperp_e) -= 5. / 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 1)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Pperp_e) -= 10. / 3. * Pperp_e / B0 * b0xcv * Grad(phi);
                           break;
                        case 2:
                           ddt(Pperp_e) -= 10. / 3. * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 3:
                           ddt(Pperp_e) -= 5. / 3. * Pperp_e / B0 * b0xcv * Grad(phi);
                           ddt(Pperp_e) -= 5. / 3. * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 4:
                           ddt(Pperp_e) -= 10. / 3. * Pperp_e / B0 * bxkappa * Grad(phi);
                           break;
                        case 5:
                           ddt(Pperp_e) -= 5. / 3. * Pperp_e / B0 * bxkappa * Grad(phi);
                           ddt(Pperp_e) -= 5. / 3. * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                     }
                  }
               }
               else
               {
                  // ddt(Pperp_e) -= 3. * n0 * Te0 / B0 * b0xcv * Grad(phi);  
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Pperp_e) -= 3. * Pe0 / B0 * b0xcv * Grad(phi);
                        break;
                     case 2:
                        ddt(Pperp_e) -= 3. * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 3:
                        ddt(Pperp_e) -= 1.5 * Pe0 / B0 * b0xcv * Grad(phi);
                        ddt(Pperp_e) -= 1.5 * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                     case 4:
                        ddt(Pperp_e) -= 3. * Pe0 / B0 * bxkappa * Grad(phi);
                        break;
                     case 5:
                        ddt(Pperp_e) -= 1.5 * Pe0 / B0 * bxkappa * Grad(phi);
                        ddt(Pperp_e) -= 1.5 * Pe0 / B0 * bracket(B0, phi, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 1)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Pperp_e) -= 3. * Pperp_e / B0 * b0xcv * Grad(phi);
                           break;
                        case 2:
                           ddt(Pperp_e) -= 3. * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 3:
                           ddt(Pperp_e) -= 1.5 * Pperp_e / B0 * b0xcv * Grad(phi);
                           ddt(Pperp_e) -= 1.5 * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                        case 4:
                           ddt(Pperp_e) -= 3. * Pperp_e / B0 * bxkappa * Grad(phi);
                           break;
                        case 5:
                           ddt(Pperp_e) -= 1.5 * Pperp_e / B0 * bxkappa * Grad(phi);
                           ddt(Pperp_e) -= 1.5 * Pperp_e / B0 * bracket(B0, phi, bm_exb);
                           break;
                     }
                  }
               }
            }

            if(compression)
            {
               if(true || isotropic)
               {
                  ddt(Pperp_e) -= 2. / 3. * Pe0 * B0 * Grad_parP(Vpar_e / B0, CELL_CENTRE);
                  if(nonlinear && nonlinear_terms > 2)
                     ddt(Pperp_e) -= 2. / 3. * Pperp_e * B0 * Grad_par(Vpar_e / B0, CELL_CENTRE);
               }
               else
               {
                  ddt(Pperp_e) -= 5. / 3. * B0*n0*Te0*Grad_parP(Vpar_e/B0, CELL_YLOW);
                  if(nonlinear && nonlinear_terms > 2)
                     ddt(Pperp_e) -= 5. / 3. * B0*n0*Pperp_e*Grad_par(Vpar_e/B0, CELL_YLOW);
               }
            }

            if(nonlinear && nonlinear_terms > 0){
               ddt(Pperp_e) -= bracket(phi, Pperp_e, bm_exb);
            }

            if(Landau_damping_e)
            {
               // ddt(Pperp_e) -= B0*B0*Grad_parP(Qperp_e/(B0*B0)); 
               ddt(Pperp_e) -= Grad_par(Qperp_e); 
            }

            if(energy_flux)
            {
               if(isotropic)
               {
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * b0xcv, Tpar_e + Tperp_e);
                        ddt(Pperp_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * b0xcv * Grad(Te0);
                        break;
                     case 2:
                        ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Pperp_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                     case 3:
                        ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * b0xcv, Tpar_e + Tperp_e);
                        ddt(Pperp_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * b0xcv * Grad(Te0);
                        ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Pperp_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                     case 4:
                        ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * bxkappa, Tpar_e + Tperp_e);
                        ddt(Pperp_e) += 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bxkappa * Grad(Te0);
                        break;
                     case 5:
                        ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pe0 * bxkappa, Tpar_e + Tperp_e);
                        ddt(Pperp_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bxkappa * Grad(Te0);
                        ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * Pe0 * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                        ddt(Pperp_e) += .5 * 5. / 3. / (Cnor * B0) * n0 * (Tpar_e + Tperp_e) * bracket(B0, Te0, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 3)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pperp_e * b0xcv, Tpar_e + Tperp_e);
                           break;
                        case 2:
                           ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * Pperp_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                        case 3:
                           ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pperp_e * b0xcv, Tpar_e + Tperp_e);
                           ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * Pperp_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                        case 4:
                           ddt(Pperp_e) -= 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pperp_e * bxkappa, Tpar_e + Tperp_e);
                           break;
                        case 5:
                           ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * V_dot_Grad(-Pperp_e * bxkappa, Tpar_e + Tperp_e);
                           ddt(Pperp_e) -= .5 * 5. / 3. / (Cnor * B0) * Pperp_e * bracket(-B0, Tpar_e + Tperp_e, bm_exb);
                           break;
                     }
                  }
               }
               else
               {
                  // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ni);
                  // ddt(Pperp_e) -= Curv(Te0/(Cnor*B0), 3.*Te0*ne);
                  // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_i+5.*Pperp_i);    
                  // ddt(Pperp_e) += Curv(Te0/(Cnor*B0), Ppar_e+5.*Pperp_e);
                  switch(curv_model)
                  {
                     case 1:
                     default:
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= Te0 / (Cnor*B0) * 3. * b0xcv * Grad(Te0 * ne);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv , Ppar_e);
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-5. * Te0 * b0xcv , Pperp_e);
                        break;
                     case 2:
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= Te0 / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * Te0 * bracket(-B0 , Ppar_e, bm_exb);
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * 5. * Te0 * bracket(-B0 , Pperp_e, bm_exb);
                        break;
                     case 3:
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= .5 * Te0 / (Cnor*B0) * 3. * b0xcv * Grad(Te0 * ne);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * b0xcv , Ppar_e);
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-5. * Te0 * b0xcv , Pperp_e);
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= .5 * Te0 / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0 , Ppar_e, bm_exb);
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * 5. * Te0 * bracket(-B0 , Pperp_e, bm_exb);
                        break;
                     case 4:
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= Te0 / (Cnor*B0) * 3. * bxkappa * Grad(Te0 * ne);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa , Ppar_e);
                        ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-5. * Te0 * bxkappa , Pperp_e);
                        break;
                     case 5:
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= .5 * Te0 / (Cnor*B0) * 3. * bxkappa * Grad(Te0 * ne);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Te0 * bxkappa , Ppar_e);
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-5. * Te0 * bxkappa , Pperp_e);
                        // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                        ddt(Pperp_e) -= .5 * Te0 / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                        // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * Te0 * bracket(-B0 , Ppar_e, bm_exb);
                        ddt(Pperp_e) -= .5 /(Cnor*B0) * 5. * Te0 * bracket(-B0 , Pperp_e, bm_exb);
                        break;
                  }
                  if(nonlinear && nonlinear_terms > 3)
                  {
                     switch(curv_model)
                     {
                        case 1:
                        default:
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= Tperp_e / (Cnor*B0) * 3. * b0xcv * Grad(Te0 * ne);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tperp_e * b0xcv , Ppar_e);
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-5. * Tperp_e * b0xcv , Pperp_e);
                           break;
                        case 2:
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= Tperp_e / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * Tperp_e * bracket(-B0 , Ppar_e, bm_exb);
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * 5. * Tperp_e * bracket(-B0 , Pperp_e, bm_exb);
                           break;
                        case 3:
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= .5 * Tperp_e / (Cnor*B0) * 3. * b0xcv * Grad(Te0 * ne);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tperp_e * b0xcv , Ppar_e);
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-5. * Tperp_e * b0xcv , Pperp_e);
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= .5 * Tperp_e / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * Tperp_e * bracket(-B0 , Ppar_e, bm_exb);
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * 5. * Tperp_e * bracket(-B0 , Pperp_e, bm_exb);
                           break;
                        case 4:
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= Tperp_e / (Cnor*B0) * 3. * bxkappa * Grad(Te0 * ne);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-Tperp_e * bxkappa , Ppar_e);
                           ddt(Pperp_e) -= 1. /(Cnor*B0) * V_dot_Grad(-5. * Tperp_e * bxkappa , Pperp_e);
                           break;
                        case 5:
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= .5 * Tperp_e / (Cnor*B0) * 3. * bxkappa * Grad(Te0 * ne);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-Tperp_e * bxkappa , Ppar_e);
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * V_dot_Grad(-5. * Tperp_e * bxkappa , Pperp_e);
                           // ddt(Pperp_e) += 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0,3.*T0*ne);
                           ddt(Pperp_e) -= .5 * Tperp_e / (Cnor*B0) * 3. * bracket(B0, Te0 * ne, bm_exb);
                           // ddt(Pperp_e) -= 0.5*T0/(Cnor*B0*B0)*b0xGrad_dot_Grad(B0, Ppar_e+5.*Pperp_e);    
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * Tperp_e * bracket(-B0 , Ppar_e, bm_exb);
                           ddt(Pperp_e) -= .5 /(Cnor*B0) * 5. * Tperp_e * bracket(-B0 , Pperp_e, bm_exb);
                           break;
                     }
                  }
               }
            }

            if(sink_Pperp_el > 0.0){
               ddt(Pperp_e) -=  sink_Pperp_el*sink_tanhxl(P0,Pperp_e,sPperp_e_widthl,sPperp_e_lengthl); // core sink
            }

            if(sink_Pperp_er > 0.0){
               ddt(Pperp_e) -=  sink_Pperp_er*sink_tanhxr(P0,Pperp_e,sPperp_e_widthr,sPperp_e_lengthr); //  sol sink
            }

         }
      }

   }

   ///////////////////////////////////////////////////
   //toroidal filter
   if(filter_z) {
      // Filter out all except filter_z_mode
      ddt(ni) = filter(ddt(ni), filter_z_mode);
      if(compression)
         ddt(Lambda_i) = filter(ddt(Lambda_i), filter_z_mode);
      ddt(Ppar_i) = filter(ddt(Ppar_i), filter_z_mode);
      ddt(Pperp_i) = filter(ddt(Pperp_i), filter_z_mode);   

      if(!electrostatic)
      {
         if(evolve_ne)
            ddt(ne) = filter(ddt(ne), filter_z_mode);
         else
            ddt(U) = filter(ddt(U), filter_z_mode);
         ddt(Psi) = filter(ddt(Psi), filter_z_mode);

         if(!Zero_Te)
         {
            if(evolve_Te)
               ddt(Te) = filter(ddt(Te), filter_z_mode);
            else
            {
               ddt(Ppar_e) = filter(ddt(Ppar_e), filter_z_mode);
               ddt(Pperp_e) = filter(ddt(Pperp_e), filter_z_mode);
            }
         }  
      }
   }

   if(low_pass_z > 0) {
      // Low-pass filter, keeping n up to low_pass_z
      ddt(ni) = lowPass(ddt(ni), low_pass_z, zonal_bkgd);
      if(compression)
         ddt(Lambda_i) = lowPass(ddt(Lambda_i), low_pass_z, zonal_bkgd);
      ddt(Ppar_i) = lowPass(ddt(Ppar_i), low_pass_z, zonal_bkgd);
      ddt(Pperp_i) = lowPass(ddt(Pperp_i), low_pass_z, zonal_bkgd);

      if(!electrostatic)
      {
         if(evolve_ne)
            ddt(ne) = lowPass(ddt(ne), low_pass_z, zonal_flow);
         else
            ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
         ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

         if(!Zero_Te){
            if(evolve_Te)
               ddt(Te) = lowPass(ddt(Te), low_pass_z, zonal_bkgd);
            else
            {
               ddt(Ppar_e) = lowPass(ddt(Ppar_e), low_pass_z, zonal_bkgd);
               ddt(Pperp_e) = lowPass(ddt(Pperp_e), low_pass_z, zonal_bkgd);
            }
         }
      }
   }

   if(nl_y > 0){
      // ddt(Lambda_i) = nl_filter_y(ddt(Lambda_i), nl_y);
      // ddt(U) = nl_filter_y(ddt(U), nl_y);
      ddt(ni) = nl_filter(ddt(ni), nl_y);
      // ddt(Ppar_i) = nl_filter_y(ddt(Ppar_i), nl_y);
      // ddt(Pperp_i) = nl_filter_y(ddt(Pperp_i), nl_y);
      // ddt(Ppar_e) = nl_filter_y(ddt(Ppar_e), nl_y);
      // ddt(Pperp_e) = nl_filter_y(ddt(Pperp_e), nl_y);
   }

   
   if(fakerun)
   {
      ddt(Lambda_i) = 0.0;
      ddt(U) = 0.0;
      ddt(ni) = 0.0;
      ddt(Psi) = 0.0;
      ddt(Ppar_i) = 0.0;
      ddt(Pperp_i) = 0.0;
      ddt(Ppar_e) = 0.0;
      ddt(Pperp_e) = 0.0;
   }

   first_run = false;

   return 0;
}

/****
// Curvature operator 
const Field3D Curv(const Field3D &f)
{
switch(curv_model) 
return V_dot_Grad(b0xcv, f);
else
return bracket(B0, f, BRACKET_STD);
}

// Curvature operator 
const Field3D Curv(const Field2D &g, const Field3D &f)
{
switch(curv_model) 
return V_dot_Grad(g * b0xcv, f);
else
return g * bracket(B0, f, BRACKET_STD);
}
 *****/

