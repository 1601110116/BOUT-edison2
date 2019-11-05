/****************************************************************
 '3+1' field gyro-fluid model for ITG simulation
  Developed by S. S. Kim and X. Q. Xu
  Beer's model [Phys. Plasmas 3, 4046(1996)] was used 
  Nonlinear part has not been completed yet.
 ****************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <cmath>

#include <initialprofiles.hxx>
#include <invert_parderiv.hxx>
#include <field_factory.hxx>
#include <utils.hxx>
#include <interpolation.hxx>
#include <math.h>
#include <msg_stack.hxx>

// #include <integrops.hxx>
// #include <inttorops.hxx>

// Solved Variables
Field3D Ni;   // guiding center density
Field3D Vpar; // parallel flow velocity
Field3D Pperp,Ppar; // perpendicular, parallel pressure 

// Derived quantities
Field3D Phi;  // electrostatic potential
Field3D gyroPhi;  // gyro-averaged potential
Field3D Vor;  // generalized vorticity
Field3D gyroPhi2,gyroPhi3;  
Field3D Tperp,Tpar;   
Field3D gyroNi,gyroTperp;
Field3D delgyroPhi,delgyroPhi2,sumgyroPhi;
Field3D sumPrs;
Field3D qpar_nl,qperp_nl,Dglf_par,Dglf_perp,Tperp1;   
Field3D r11,r12,r21,r22;   
  
Field2D phia,phid;
Field2D gyroa,gyrod;
Field2D tau;

Field3D F1;
Field2D F0;
Field3D qvar_nl;
BoutReal Lpar,k0;
Field2D gradpar_lnB,lnB;

// Flux Surface Averages
Field2D NiHat,TeHat;
Field2D TiHat0,PrsHat0;
Field2D TiHat_tot,PrsHat_tot;

// Equilibrium quantities
Vector2D b0xcv,b0xglB0; // curvature term
Field2D Rxy, Bpxy, Btxy, Bxy, hthe, I_shear;
Field2D rho,pol_angle;
Field2D nu,qprof;
Vector2D B0vec; // equilibrium B field vector

// For the calculation of modomegad 
Vector2D vd; // grad-B drift
Field2D modb0xglB0, b0xglB0xhat, b0xglB0yhat, b0xglB0zhat; 
Field3D mod_wd_Vpar,mod_wd_Tpar,mod_wd_Tperp;
BoutReal nu1r,nu2r,nu3r,nu4r,nu5r;
BoutReal nu1i,nu2i,nu3i,nu4i,nu5i;
const BoutReal zperiod = 1.0;

// [itg] options
bool nonlinear;  
bool FLR_nonlinearity;  
bool add_qperp_FLR_nonlin;
bool glf;  
bool nonFourier;
int nytimes;
bool add_q_gradparlnB;
bool add_qperp_gradparlnB;
bool add_vpar_gradparlnB;
bool filter_z;
int filter_z_mode;
bool filter_y_zonal;
int mmax;
int phi_flags,gyrophi_flags;  
bool cylind;
bool toroidal_closure;
bool include_modomegad;
BoutReal curv_Ti;   
int curv_model;
int low_pass_z;
int zonal_flow;
int zonal_field;
int zonal_bkgd;

BRACKET_METHOD bm = BRACKET_STD;
//BRACKET_METHOD bm = BRACKET_ARAKAWA;

// Fundamental constants
const BoutReal qe = 1.602e-19;      // Electron charge
const BoutReal Mp = 1.67262158e-27; // Proton mass
BoutReal AA = 2.0;                  // Ion atomic mass (1.0:Proton)
BoutReal ZZ = 1.0;                  // Ion charge

// Normalisation factors
BoutReal Lbar;   // Scale length = minor radius
BoutReal Tbar;   // Timescale = Lbar / Cs
BoutReal Ti0;    // Ion temperature 
BoutReal Bbar;   // Magnetic field 

// Basic parameters
BoutReal Rmajor,aminor; // Major and minor radii
BoutReal epsilon; // Inversed aspect ratio
BoutReal Te0;     // Electron temperature 
BoutReal Cs;      // Ion thermal speed : sqrt(Ti0/Mi)
BoutReal rho_i;   // Ion Larmor Radius
BoutReal delta;   // rho*=rho_i/a in this model

FieldGroup comms; // Communications

const Field3D Curv(const Field3D &f);
const Field3D Landau_damping(const Field3D &tvar);
const Field3D Gyro(const Field3D &f);
const Field3D mod_wd_F(const Field3D &f);
const Field3D LowPass_y(const Field3D &f, int m0);
const Field2D F0BC(const Field2D &f);
const Field3D F1BC(const Field3D &f);
const Field2D ZeroZero(const Field3D &f);

int i;
Field3D f1;
Field2D Qi,Pi_Rey,gradTi,Chi_i,I_phi;
Vector2D b0vec,rhat;
Field3D Vr;
Field2D Lti,R_Lti,gradNi,Lne,R_Lne,Chi_i_conv,time_conv;

////////////////////////////////////////////////////////////////////////
int physics_init(bool restarting)
{
    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("itg");

    // Read [itg] options
    OPTION(options, nonlinear, false);
    OPTION(options, FLR_nonlinearity, false);
    OPTION(options, add_qperp_FLR_nonlin, false);
    OPTION(options, filter_z, true);
    OPTION(options, filter_z_mode, 1);
    OPTION(options, low_pass_z, 16);
    OPTION(options, zonal_flow, 0);      // zonal flow filter
    OPTION(options, zonal_field, 0);     // zonal field filter
    OPTION(options, zonal_bkgd, 0);      // zonal background P filter
    OPTION(options, filter_y_zonal, false);
    OPTION(options, mmax, 3);
    OPTION(options, glf, true); 
    OPTION(options, nonFourier, true); 
    OPTION(options, nytimes, 9); 
    OPTION(options, add_q_gradparlnB, true); 
    OPTION(options, add_qperp_gradparlnB, true); 
    OPTION(options, add_vpar_gradparlnB, false); 
    OPTION(options, toroidal_closure, true);  
    OPTION(options, include_modomegad, true);  
    OPTION(options, curv_model, 1);  
    OPTION(options, curv_Ti, 1.0);  
    OPTION(options, cylind, false);
    OPTION(options, phi_flags,  1); 
    OPTION(options, gyrophi_flags,  1); 

    // Read grid file
    GRID_LOAD(Rxy);    // Major radius [m]
    GRID_LOAD(hthe);   // Poloidal arc length [m/radian]
    GRID_LOAD(Bpxy);   // Poloidal B field [T]
    GRID_LOAD(Btxy);   // Toroidal B field [T]
    GRID_LOAD(Bxy);    // Total B field [T]

    GRID_LOAD(PrsHat0);
    GRID_LOAD(TiHat0);
    GRID_LOAD(TeHat);
    GRID_LOAD(NiHat);
    GRID_LOAD(rho);
    GRID_LOAD(pol_angle);

    mesh->get(I_shear, "sinty");  // dnu/dpsi*theta

    // b0 x kappa
    b0xcv.covariant=false;
    mesh->get(b0xcv.x,  "bxcvx");
    mesh->get(b0xcv.y,  "bxcvy");
    mesh->get(b0xcv.z,  "bxcvz");

    // Basic parameters (SI units)
    mesh->get(Ti0, "Ti0");
    mesh->get(Bbar, "Bbar");
    mesh->get(aminor, "aminor");
    mesh->get(Rmajor, "Rmajor");

    Cs = sqrt(qe*Ti0 / (AA*Mp));  // Ion thermal speed 
    Lbar = aminor;
    Tbar = Lbar / Cs;             // Timescale
    rho_i = Cs / (ZZ*qe*Bbar/(AA*Mp)); // Gyroradius 
    delta = rho_i / Lbar;  // Here, delta means rho*
    epsilon = aminor/Rmajor;

    Lpar=2.0*M_PI*(0.854+2.184*0.5*0.5)*Rmajor; // =2*pi*q*R
    k0=0.05*2.*M_PI/Lpar;

    output << "\n\tParameters\n";
    output.write("\tLbar = %e [m], rho_i = %e [m]\n",Lbar, rho_i);
    output.write("\tTbar = %e [s], Cs = %e [m/s]\n", Tbar,Cs);
    output.write("\trho_star = %e, epsilon = %e \n", delta,epsilon);
    output.write("\tBbar = %e [T], Ti0 = %e [eV]\n", Bbar,Ti0);
    output.write("\tLpar = %e [m], k0 = %e [1/m]\n", Lpar,k0);

    // Normalize
    b0xcv.x /= Bbar;
    b0xcv.y *= Lbar*Lbar;
    b0xcv.z *= Lbar*Lbar;
    hthe /= Lbar; 
    Bpxy /= Bbar;
    Btxy /= Bbar;
    Bxy  /= Bbar;

    if(cylind) Rxy = Rmajor;
    Rxy  /= Lbar; 

    mesh->dx /= Lbar*Lbar*Bbar;
    if(mesh->non_uniform)
      mesh->d1_dx *= (Lbar*Lbar*Bbar);
    I_shear *= Lbar*Lbar*Bbar;

    Lpar /= aminor;
    k0 *= aminor;

    // local pitch
    nu=Btxy*hthe/(Bpxy*Rxy);
    dump.add(nu, "nu", 0);
    // q profile 
    qprof=mesh->averageY(nu);
    dump.add(qprof,  "qprof",  0);

    // SHIFTED RADIAL COORDINATES
    if(mesh->ShiftXderivs) {
      if(mesh->IncIntShear) {
        // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
        mesh->IntShiftTorsion = I_shear;
      }else {
        // Dimits style, using local coordinate system
        b0xcv.z += I_shear*b0xcv.x;
        I_shear = 0.0;  // I disappears from metric
      }
    }

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
    mesh->g_33 = Rxy*Rxy;
    mesh->g_23 = (Btxy*hthe*Rxy)/Bpxy;
    mesh->g_12 = I_shear*mesh->g_23;
    mesh->g_13 = I_shear*(Rxy^2);

    mesh->geometry();

    // Set B field vector  
    B0vec.covariant = false;
    B0vec.x = 0.;
    B0vec.y = Bpxy / hthe;
    B0vec.z = 0.;

    // Add Equations to be solved
    SOLVE_FOR(Ni);    comms.add(Ni);
    SOLVE_FOR(Vpar);  comms.add(Vpar);
    SOLVE_FOR(Ppar);  comms.add(Ppar);
    SOLVE_FOR(Pperp); comms.add(Pperp);
    dump.add(Tperp,"Tperp",1);
    dump.add(Tpar,"Tpar",1);

    // dump output
    dump.add(Phi,"Phi",1);

    dump.add(PrsHat0,"PrsHat0",0);
    dump.add(TiHat0, "TiHat0", 0);
    dump.add(TeHat,  "TeHat",  0);
    dump.add(NiHat,  "NiHat",  0);

    dump.add(rho_i,"rho_i", 0);
    dump.add(delta,"delta", 0);
    dump.add(Lbar, "Lbar", 0);
    dump.add(Tbar, "Tbar", 0);
    dump.add(Bbar, "Bbar", 0);
    dump.add(Cs,   "Cs", 0); 

    if(nonlinear){
      dump.add(Qi,"Qi",1);
      dump.add(Pi_Rey,"Pi_Rey",1);
      dump.add(gradTi,"gradTi",1);
      dump.add(Chi_i,"Chi_i",1);
      dump.add(I_phi,"I_phi",1);
    }

    b0vec=B0vec/Bxy;
    mesh->communicate(rho);
    rhat=Grad(rho);

    Rmajor/=Lbar;
    aminor/=Lbar;

    if(nonlinear){
      gradTi=mesh->averageY(-rhat*Grad(TiHat0))+1.e-30; // [T0/a]
      Lti = abs(TiHat0/gradTi)+1.e-30; // [a]
      R_Lti=Rmajor/Lti;
      gradNi=mesh->averageY(-rhat*Grad(NiHat))+1.e-30; // [n0/a]
      Lne = abs(NiHat/gradNi)+1.e-30; // [a]
      R_Lne=Rmajor/Lne;
      Chi_i_conv = Lne/aminor/TiHat0; // Chi_i --> [rhoi^2*vthi/Lne], vthi=sqrt(Ti/m)
      time_conv = aminor/Lti*sqrt(TiHat0); // [a/Cs0] --> [Lti/vthi]

      dump.add(gradNi,"gradNi",0);
      dump.add(Lti,"Lti",0);
      dump.add(Lne,"Lne",0);
      dump.add(R_Lti,"R_Lti",0);
      dump.add(R_Lne,"R_Lne",0);
      dump.add(Chi_i_conv,"Chi_i_conv",0);
      dump.add(time_conv,"time_conv",0);
    }

    F1.setBoundary("F1");         
    F0.setBoundary("F0");         

    nu1r= 1.232;  nu1i= 0.437;
    nu2r=-0.912;  nu2i= 0.362;
    nu3r=-1.164;  nu3i= 0.294;
    nu4r= 0.478;  nu4i=-1.926;
    nu5r= 0.515;  nu5i=-0.958;

    lnB = log(Bxy);
    mesh->communicate(lnB);
    gradpar_lnB=Grad_par(lnB);

    // b0xgradlogB0(Bxy,modb0xglB0,b0xglB0xhat,b0xglB0yhat,b0xglB0zhat,b0xglB0);

    cout << "init done. \n";

    return 0;
}


////////////////////////////////////////////////////////////////////////
int physics_run(BoutReal t)
{
    PrsHat_tot = PrsHat0; // Pressure profile is fixed
    TiHat_tot  = TiHat0;  

    Tperp = F1BC((Pperp  - (Ni * TiHat_tot))/NiHat);
    Tpar = F1BC((Ppar  - (Ni * TiHat_tot))/NiHat);

    tau = TiHat_tot/TeHat;

    phia = tau/(1.0+tau);
    phid = -TiHat_tot*delta*delta/Bxy/Bxy;

    gyroa = 1.0;               
    gyrod = 0.5*phid;

  // calculate vorticity from guiding center density and Ti
    gyroNi = Gyro(Ni);
    gyroTperp = Gyro(Tperp);
    Vor = gyroNi+NiHat/TiHat_tot*Gyro(gyroTperp-Tperp);

  // calculate potential from vorticity 
    Vor = F1BC(TiHat_tot/NiHat/(1.0+tau)*Vor);
    Phi = invert_laplace(Vor, phi_flags, &phia, NULL, &phid);
    Phi = F1BC(Phi/(1.0+tau) + Vor);

    if(nonlinear && filter_y_zonal)
      Phi = LowPass_y(Phi, mmax);

  // calculate gyroaveraged-potential from potential 
    gyroPhi  = Gyro(Phi);
    if(nonlinear && filter_y_zonal)
      gyroPhi  = LowPass_y(gyroPhi, mmax);

    gyroPhi2 = Gyro(gyroPhi);
    if(nonlinear && filter_y_zonal)
      gyroPhi2 = LowPass_y(gyroPhi2, mmax);

    gyroPhi3 = Gyro(gyroPhi2);
    if(nonlinear && filter_y_zonal)
      gyroPhi3 = LowPass_y(gyroPhi3, mmax);

    delgyroPhi =gyroPhi2-gyroPhi;
    delgyroPhi2=gyroPhi3-gyroPhi2;
    sumgyroPhi =gyroPhi2+gyroPhi;

    r11 = F1BC(3.0*(Ppar  + curv_Ti*NiHat*Tpar));
    r12 = F1BC(Ppar  + curv_Ti*NiHat*Tperp); 
    r21 = F1BC(Pperp + curv_Ti*NiHat*Tpar);
    r22 = F1BC(2.0*(Pperp + curv_Ti*NiHat*Tperp));

    mesh->communicate(comms);

    sumPrs = Ppar + Pperp;

  // Guidinig center density Eq.
    ddt(Ni)  =-Div_par(NiHat*Vpar) ;
    ddt(Ni) += 0.5*NiHat*(Curv(sumgyroPhi) + Curv(sumPrs)/NiHat); 
    ddt(Ni) += NiHat/TiHat_tot*delta*bracket(TiHat_tot,delgyroPhi,bm);
    ddt(Ni) += delta*bracket(NiHat,gyroPhi,bm);

  // Parallel Velocity Eq.
    ddt(Vpar)  = -Grad_par(gyroPhi)-Div_par(Ppar)/NiHat;
    ddt(Vpar) += -(Tperp + delgyroPhi)*gradpar_lnB;
    ddt(Vpar) += 2.0*TiHat_tot*Curv(Vpar);

  // Parallel Pressure Eq.
    ddt(Ppar)  = NiHat*delta*bracket(TiHat_tot,delgyroPhi,bm);
    ddt(Ppar) -= Div_par(PrsHat_tot*Vpar) + 2.0*PrsHat_tot*Div_par(Vpar);
    ddt(Ppar) += PrsHat_tot*Curv(1.5*gyroPhi+0.5*gyroPhi2);
    if(toroidal_closure)
      ddt(Ppar) += 0.5*TiHat_tot*Curv(r11+r12);

  // Perpendicular Pressure Eq.
    ddt(Pperp)  = 2.0*NiHat*delta*bracket(TiHat_tot,delgyroPhi2,bm);
    ddt(Pperp) -= Div_par(PrsHat_tot*Vpar);
    ddt(Pperp) += PrsHat_tot*Vpar*gradpar_lnB; 
    ddt(Pperp) += PrsHat_tot*Curv(gyroPhi3 + 0.5*gyroPhi2);
    if(toroidal_closure)
      ddt(Pperp) += 0.5*TiHat_tot*Curv(r21+r22);

    if(toroidal_closure && include_modomegad){

      if(curv_model==1) // use 0.5*(b0xcv + b0xGrad(lnB))
        vd = 0.5*TiHat_tot/Bxy*(b0xcv+b0xglB0); 
      else              // use b0xGrad(lnB)
        vd = TiHat_tot/Bxy*b0xglB0; 

      (vd.x).applyBoundary("neumann");
      (vd.y).applyBoundary("neumann");
      (vd.z).applyBoundary("neumann");

      mod_wd_Vpar = delta*mod_wd_F(Vpar);
      mod_wd_Tpar = delta*mod_wd_F(Tpar);
      mod_wd_Tperp= delta*mod_wd_F(Tperp);

      ddt(Vpar) -= 2.0*nu5r*mod_wd_Vpar; 
      ddt(Vpar) += nu5i*TiHat_tot*Curv(Vpar); 
      ddt(Ppar) -= 2.0*NiHat*(nu1r*mod_wd_Tpar+nu2r*mod_wd_Tperp); 
      ddt(Ppar) += PrsHat_tot*(nu1i*Curv(Tpar)+nu2i*Curv(Tperp));
      ddt(Pperp)-= 2.0*NiHat*(nu3r*mod_wd_Tpar+nu4r*mod_wd_Tperp);
      ddt(Pperp)+= PrsHat_tot*(nu3i*Curv(Tpar)+nu4i*Curv(Tperp));
    }

    if(glf){
     Tperp1=Tperp+delgyroPhi;
  // Tperp1=Tperp;

     Dglf_par  = sqrt(8.0/M_PI*TiHat_tot)*NiHat*Landau_damping(Tpar);
     Dglf_perp = sqrt(2.0/M_PI*TiHat_tot)*NiHat*Landau_damping(Tperp1);

     qperp_nl = sqrt(2.0/M_PI*TiHat_tot)*NiHat*qvar_nl;
     if(add_qperp_gradparlnB){ 
       Dglf_par -= 2.0*qperp_nl*gradpar_lnB;
       Dglf_perp += qperp_nl*gradpar_lnB;
     }
     if(add_vpar_gradparlnB){ 
       Dglf_par -= 2.0*PrsHat_tot*Vpar*gradpar_lnB;
     }

     Dglf_par=F1BC(Dglf_par);
     Dglf_perp=F1BC(Dglf_perp);

     ddt(Ppar) += Dglf_par;
     ddt(Pperp) += Dglf_perp;
    }

  // ExB advection
    ddt(Ppar)  += delta*bracket(PrsHat_tot,gyroPhi,bm); // source for linear ITG 
    ddt(Pperp) += delta*bracket(PrsHat_tot,gyroPhi2,bm);// source for linear ITG

    if(nonlinear){
      ddt(Ni)    -= delta*delta*bracket(gyroPhi,Ni,bm);
      ddt(Vpar)  -= delta*delta*bracket(gyroPhi,Vpar,bm);
      ddt(Ppar)  -= delta*delta*bracket(gyroPhi,Ppar,bm);
      ddt(Pperp) -= delta*delta*bracket(gyroPhi,Pperp,bm);
      if(FLR_nonlinearity){
        ddt(Ni)    -= NiHat/TiHat_tot*delta*delta*bracket(delgyroPhi,Tperp,bm);
        if(add_qperp_FLR_nonlin){
          ddt(Vpar) -= delta*delta*bracket(delgyroPhi,F1BC(qperp_nl/TiHat_tot),bm);
        }
        ddt(Ppar)  -= NiHat*delta*delta*bracket(delgyroPhi,Tperp,bm);
        ddt(Pperp) -= delta*delta*bracket(delgyroPhi,Pperp,bm);
        ddt(Pperp) -= 2.0*NiHat*delta*delta*bracket(delgyroPhi2,Tperp,bm);
      }
    }

  // Filtering
    if(filter_z){
      ddt(Ni)    = filter(ddt(Ni),    filter_z_mode);
      ddt(Vpar)  = filter(ddt(Vpar),  filter_z_mode);
      ddt(Ppar)  = filter(ddt(Ppar),  filter_z_mode);
      ddt(Pperp) = filter(ddt(Pperp), filter_z_mode);
    }

    if(nonlinear && filter_y_zonal){
      ddt(Ni)    = LowPass_y(ddt(Ni),    mmax);
      ddt(Vpar)  = LowPass_y(ddt(Vpar),  mmax);
      ddt(Ppar)  = LowPass_y(ddt(Ppar),  mmax);
      ddt(Pperp) = LowPass_y(ddt(Pperp), mmax);
    }

    if(nonlinear && (low_pass_z > 0)) {
      // Low-pass filter, keeping n up to low_pass_z
      ddt(Ni)    = lowPass(ddt(Ni),    low_pass_z, zonal_flow);
      ddt(Vpar)  = lowPass(ddt(Vpar),  low_pass_z, zonal_field);
      ddt(Ppar)  = lowPass(ddt(Ppar),  low_pass_z, zonal_bkgd);
      ddt(Pperp) = lowPass(ddt(Pperp), low_pass_z, zonal_bkgd);
    } 

    if(nonlinear){
      Vr=delta*(rhat*(b0vec^Grad(Phi)))/mesh->Bxy; // ExB velocity
      Qi=1.5*ZeroZero(Vr*(Ppar+2.0*Pperp)/3.0);  // Turbulent heat flux
      Pi_Rey=ZeroZero(Vr*Vpar); // Parallel Reynolds stress
      gradTi=mesh->averageY(-rhat*Grad(TiHat_tot))+1.e-30; // Ti gradient 
      Chi_i =Qi/gradTi/NiHat; // Thermal diffusivity
      I_phi=ZeroZero(Phi^2);  // Turbulence intensity
    }

  return 0;
}

// Curvature operator 
const Field3D Curv(const Field3D &f) // = -2*i*omega_d(f)/Ti
{
  if(curv_model==1)  // use 0.5*(b0xcv + b0xGrad(lnB))
    return -delta*(b0xcv*Grad(f)+bracket(Bxy, f, bm))/Bxy; 
  else               // use b0xGrad(lnB)
    return -2.0*delta*bracket(Bxy,f,bm)/Bxy;
}

// Landau damping operator 
const Field3D Landau_damping(const Field3D &tvar)
{
/*******
    if(nonFourier){
      qvar_nl = F1BC(-iSign_kpar(tvar,k0));
    }else{
      qvar_nl = F1BC(-mesh->iSign_kpar_FFTy(tvar,nytimes));
    }

    if(add_q_gradparlnB){
      return -Div_par(qvar_nl);
    }else{
      return -Grad_par(qvar_nl);
    }
******/
   return Field3D(0.);
}

// Gyro-averaging operator 
const Field3D Gyro(const Field3D &f)
{
  return F1BC(invert_laplace(f, gyrophi_flags, &gyroa, NULL, &gyrod));
}

// |omega_d| operator
const Field3D mod_wd_F(const Field3D &f)
{
   /***
  if(nonlinear){
   F1=0.0;
   for(i=1;i<=low_pass_z;i++){
     f1 = filter(f,i);
     F1+= mod_wd(f1, &vd,i/20.);
   }
  }else
   F1= mod_wd(f, &vd,1/20.);
  F1 = F1BC(F1);
  return F1;
  ***/
  return Field3D(0.);
}

// Lowpass filter for n=0 components
const Field3D LowPass_y(const Field3D &f, int m0)
{
  F0 = F0BC(f.DC());
  F1 = f - F0; // extract n=0 components
  F0 = mesh->lowPass_poloidal(F0,m0); // filter (m>m0,n=0) modes
  F1 += F0;    // add (m<=m0,n=0) modes 
  F1 = F1BC(F1);
  return F1;
}

const Field2D F0BC(const Field2D &f)
{
  F0=f;
  F0.applyBoundary();
  mesh->communicate(F0);
  return F0;
}

const Field3D F1BC(const Field3D &f)
{
  F1=f;
  F1.applyBoundary();
  mesh->communicate(F1);
  return F1;
}

const Field2D ZeroZero(const Field3D &f)
{   
  return mesh->averageY(F0BC(f.DC()));
}

