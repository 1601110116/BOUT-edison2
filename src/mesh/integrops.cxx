/**************************************************************************
 * Various differential operators defined on BOUT grid
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <globals.hxx>
#include <bout.hxx>
#include <boutexception.hxx>
#include <integrops.hxx>
#include <difops.hxx>
#include <vecops.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <invert_parderiv.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

#include <interpolation.hxx>

#include <math.h>
#include <stdlib.h>

/*******************************************************************************
 * iSign_kpar()
 * Integral operator calculates non-local parallel heat flux for given
 * parallel temperature profile
 *******************************************************************************/

const Field3D iSign_kpar(const Field3D &vin, BoutReal k0, const int N)
{

  Field3D qout = 0.0;

  const BoutReal alp=5.0;
  const BoutReal bet=1.04;

  // InvertPar *invpar;
  // invpar = InvertPar::Create();
  // much better if static:

  static InvertPar* invpar = InvertPar::Create();

  //-range of Lorentzian terms
  const int nMin=0;
  const int nMax=N;

  //-spectral shift parameter
  //replace the default by the right value here
  //-use empirical rule k0=min(kpar)/20., min(kpar)=2PI/Lpar, where Lpar is parallel domain size
  // Note: Lpar is in same units as your equations, if normalized then Lpar normalized same way

  //const BoutReal k0=1.0; //-this is the default
  //const BoutReal k0=(2*PI/Lpar)/20.;

  if ( nMax < 7)
    throw BoutException("ERROR: At least 7 Lorentzians required for collisionless Landau closure fitting!\n");

  for (int n=nMin; n<nMax; n++)
  //for (int n=1; n<=1; n++) /*just for testing*/
    {

      BoutReal alpn=k0*pow(alp,n);
      BoutReal alpn2=alpn*alpn;

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(alpn2); //-free term
      invpar->setCoefB(-1.0);  //-d2/dz2 term

      qout += invpar->solve(bet*alpn*Grad_par(vin));

    }

#if 0
  delete invpar; //-deallocate the memory associated with the pointer
  //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

  return(qout);
}

/*******************************************************************************/

/*******************************************************************************
 * iSign_kpar_wcoll()
 * Integral operator calculates non-local parallel heat flux with collisions for given
 * parallel temperature profile
 *******************************************************************************/

const Field3D iSign_kpar_wcoll(const Field3D &vin, const Field3D &k0, const int N)
{

   Field3D qout = 0.0;
   Field3D ddyvin;

   ddyvin = DDY(vin);
   mesh->communicate(ddyvin);

   // InvertPar *invpar;
   // invpar = InvertPar::Create();
   // much better if static:

   static InvertPar* invpar = InvertPar::Create();

   //-range of Lorentzian terms
   const int nMin=0;
   const int nMax= (N == 0 ? 3 : N);

   // set coefficient for Lorentzian fitting
   const BoutReal *alpha;
   const BoutReal *beta;

   // N = 3:
   // Old coef. from  M.V. Umansky, et. al.,  Journal of Nuclear Materials 463, 506 (2015)
   // WARNING: deprecated
   const BoutReal alpha0[] = {0.0018, 0.0769, 2.4498};
   const BoutReal beta0[] = {0.1192, 0.4913, 2.1495};
   // new coef
   // N = 3:
   const BoutReal alpha3[] = {0.01315, 0.924, 14.1365};
   const BoutReal beta3[] = {0.2044, 1.3587, 8.9643};
   // N = 7:
   const BoutReal alpha7[] = { \
            0.007438, 0.6161, 5.9804, 37.9822, 234.3654, \
            1466.4331, 14981.4634};
   const BoutReal beta7[] = {0.1678, 1.1106, 5.6457, 33.1536, \
            202.738, 1254.2144, 9275.3323};
   // N = 12:
   const BoutReal alpha12[] = {0.001424, 0.20736, 2.5653, 14.927, \
             79.3050, 419.2399, 2215.7233, 11709.7857, \
             61885.2763, 327392.6096, 1773350.1566, 16903628.3745};
   const BoutReal beta12[] = {0.09419, 0.6741, 2.9628, 14.43958, \
             75.1106, 395.8293, 2090.8877, 11049.1471, \
             58392.0969, 308695.7371, 1645460.1472, 10794779.4293};
   switch (N) {
     case 0: {
       // WARNING: deprecated
       alpha = alpha0;
       beta = beta0;
       break;
       }
     case 3: {
       alpha = alpha3;
       beta = beta3;
       break;
     }
     case 7: {
       alpha = alpha7;
       beta = beta7;
       break;
     }
     case 12: {
       alpha = alpha12;
       beta = beta12;
       break;
     }
     default: {
       throw BoutException("ERROR: Invalid choice of number of Lorentzian terms. Must be in [3, 7, 12]\n");
     }
   }

   for (int n=nMin; n<nMax; n++)
      //for (int n=1; n<=1; n++) /*just for testing*/
   {

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(beta[n] * beta[n]); //-free term
      invpar->setCoefB(-1.0 / (k0 * k0));  //-d2/dz2 term

      qout += invpar->solve(alpha[n] / k0 * ddyvin);

   }

#if 0
   delete invpar; //-deallocate the memory associated with the pointer
   //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

   return(qout);
}

const Field2D iSign_kpar_wcoll(const Field2D &vin, const Field2D &k0, const int N)
{

   Field2D qout = 0.0;
   Field2D ddyvin;

   ddyvin = DDY(vin);
   mesh->communicate(ddyvin);

   // InvertPar *invpar;
   // invpar = InvertPar::Create();
   // much better if static:

   static InvertPar* invpar = InvertPar::Create();

   //-range of Lorentzian terms
   const int nMin=0;
   const int nMax= (N == 0 ? 3 : N);

   // set coefficient for Lorentzian fitting
   const BoutReal *alpha;
   const BoutReal *beta;

   // N = 3:
   // Old coef. from  M.V. Umansky, et. al.,  Journal of Nuclear Materials 463, 506 (2015)
   const BoutReal alpha0[] = {0.0018, 0.0769, 2.4498};
   const BoutReal beta0[] = {0.1192, 0.4913, 2.1495};
   // new coef
   // N = 3:
   const BoutReal alpha3[] = {0.01315, 0.924, 14.1365};
   const BoutReal beta3[] = {0.2044, 1.3587, 8.9643};
   // N = 7:
   const BoutReal alpha7[] = { \
            0.007438, 0.6161, 5.9804, 37.9822, 234.3654, \
            1466.4331, 14981.4634};
   const BoutReal beta7[] = {0.1678, 1.1106, 5.6457, 33.1536, \
            202.738, 1254.2144, 9275.3323};
   // N = 12:
   const BoutReal alpha12[] = {0.001424, 0.20736, 2.5653, 14.927, \
             79.3050, 419.2399, 2215.7233, 11709.7857, \
             61885.2763, 327392.6096, 1773350.1566, 16903628.3745};
   const BoutReal beta12[] = {0.09419, 0.6741, 2.9628, 14.43958, \
             75.1106, 395.8293, 2090.8877, 11049.1471, \
             58392.0969, 308695.7371, 1645460.1472, 10794779.4293};
   switch (N) {
     case 0: {
       alpha = alpha0;
       beta = beta0;
       break;
       }
     case 3: {
       alpha = alpha3;
       beta = beta3;
       break;
     }
     case 7: {
       alpha = alpha7;
       beta = beta7;
       break;
     }
     case 12: {
       alpha = alpha12;
       beta = beta12;
       break;
     }
     default: {
       throw BoutException("ERROR: Invalid choice of number of Lorentzian terms. Must be in [3, 7, 12]\n");
     }
   }

   for (int n=nMin; n<nMax; n++)
      //for (int n=1; n<=1; n++) /*just for testing*/
   {

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(beta[n] * beta[n]); //-free term
      invpar->setCoefB(-1.0 / (k0 * k0));  //-d2/dz2 term

      qout += invpar->solve(alpha[n] / k0 * ddyvin);

   }

#if 0
   delete invpar; //-deallocate the memory associated with the pointer
   //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

   return(qout);
}

