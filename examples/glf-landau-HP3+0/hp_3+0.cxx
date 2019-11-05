/*******************************************************************************
 * Parallel Poisson inversion test
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <invert_parderiv.hxx>
#include <field_factory.hxx>
#include <utils.hxx>

//-2D initial profiles
Field2D phi0, rho0;


//-3D evolving fields
Field3D rho;
Field3D nvar, uvar, tvar;


//-Derived 3D variables
Field3D phi;
Field3D qvar_l;    //-local heat flux
Field3D qvar_nl; //-nonlocal heat flux


//-Non-linear coefficients
Field3D kapa_Te;


//-Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe, I;


//-Parameters
BoutReal acoef, bcoef;


//-Settings
bool evolve_nvar, evolve_uvar, evolve_tvar;
bool evolve_rho;
int phi_flags; // Inversion flags


//-Group fields together for communication
FieldGroup comms;



InvertPar *invpar;

//FieldFactory f;
//-coefs for parallel poisson operator
//Field2D A, B;




int physics_init(bool restarting)
{
  output << "Poisson inversion\n";

  invpar = InvertPar::Create();


  /************* LOAD DATA FROM GRID FILE ****************/

  //-Load 2D profiles (set to zero if not found)
  GRID_LOAD(phi0);
  GRID_LOAD(rho0);


  //-Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->dx,   "dpsi");
  mesh->get(I,    "sinty");
  mesh->get(mesh->zShift, "qinty");

  //-Load normalisation values



  /*************** READ OPTIONS *************************/

  // Read some parameters
  Options *globalOptions = Options::getRoot();

  Options *options = globalOptions->getSection("invpar");

  OPTION(options, acoef,   1.0);
  OPTION(options, bcoef,   0.0);
  output.write("Using acoef=%f, bcoef=%f\n", acoef, bcoef);

  OPTION(options,phi_flags,  0);

  (globalOptions->getSection("rho"))->get("evolve", evolve_rho,   true);

  (globalOptions->getSection("nvar"))->get("evolve", evolve_nvar,   true);
  (globalOptions->getSection("uvar"))->get("evolve", evolve_uvar,   true);
  (globalOptions->getSection("tvar"))->get("evolve", evolve_tvar,   true);



  //- SHIFTED RADIAL COORDINATES -//
  //- SHIFTED RADIAL COORDINATES -//
  //- CALCULATE PARAMETERS -//
  //- PRINT Z INFORMATION -//
  //- SHIFTED GRIDS LOCATION -//
  //- NORMALISE QUANTITIES -//

  //- CALCULATE METRICS -//
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



  /**************** SET EVOLVING VARIABLES *************/
  if(evolve_rho) {
    bout_solve(rho, "rho");
    comms.add(rho);
    output.write("rho\n");
  }else
    initial_profile("rho", rho);



  if(evolve_nvar) {
    bout_solve(nvar, "nvar");
    comms.add(nvar);
    output.write("nvar\n");
  }else
    initial_profile("nvar", nvar);

  if(evolve_uvar) {
    bout_solve(uvar, "uvar");
    comms.add(uvar);
    output.write("uvar\n");
  }else
    initial_profile("uvar", uvar);

  if(evolve_tvar) {
    bout_solve(tvar, "tvar");
    comms.add(tvar);
    output.write("tvar\n");
  }else
    initial_profile("tvar", tvar);


  /************** SETUP COMMUNICATIONS **************/
  // add extra variables to communication
  comms.add(phi);
  comms.add(qvar_l);
  comms.add(qvar_nl);


  // Add any other variables to be dumped to file
  dump.add(qvar_l,  "qvar_l",  1);
  dump.add(qvar_nl, "qvar_nl",  1);
  dump.add(phi,     "phi",  1);
  //dump.add(deriv,   "deriv",  1);

  dump.add(acoef,  "acoef",  0);
  dump.add(bcoef,  "bcoef",  0);
  
  return(0);
}



int get_qpar_local(Field3D &tin, Field3D &qout)
{
  //-local heat diffusion
  qout = - 1e-6*Grad_par(tin); //-arbitrary scaling factor to make it run faster
  return(0);
}


int get_qpar_nonlocal(Field3D &tin, Field3D &qout)
{

  qout = 0.0;

  static const BoutReal alp=5.0;
  static const BoutReal bet=1.04;
  static Field3D dummy3D;


  //-range of Lorentzian terms
  static const int nMin=0;
  static const int nMax=7; //7;


  for (int n=nMin; n<=nMax; n++)
    {

      BoutReal alpn=pow(alp,n);
      BoutReal alpn2=alpn*alpn;

      //-invert parallel poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(alpn2); //-free term
      invpar->setCoefB(-1.0);  //-d2/dz2 term
      dummy3D = invpar->solve(-bet*alpn*Grad_par(tin));      

      qout = qout + 1e-4*dummy3D; //-arbitrary scaling factor to make it run faster
           
    }


  return(0);
}


//*********************************************************************************************//


int physics_run(BoutReal t)
{
  //-Solve for EM fields

  //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
  invpar->setCoefA(acoef);
  invpar->setCoefB(bcoef);
  phi = invpar->solve(rho);


  get_qpar_local(tvar, qvar_l);
  get_qpar_nonlocal(tvar, qvar_nl);


  //-Communicate variables
  mesh->communicate(comms);



  //-----------------set RHS of time-evolution equations

  // VORTICITY
  ddt(rho) = 0.0;


  //-1D fluid equations
  ddt(nvar) = -Div_par(uvar);
  ddt(uvar) = -Grad_par(nvar) -Grad_par(tvar);


  ddt(tvar) = 0.0;
  if(evolve_tvar) {
    ddt(tvar) += -2.0*Div_par(uvar);
    //ddt(tvar) +=  -Div_par(qvar_l);
    ddt(tvar) +=  -Div_par(qvar_nl);
  }

  

  return(0);
}

//*********************************************************************************************//
