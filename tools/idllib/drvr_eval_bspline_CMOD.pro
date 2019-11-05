; Evaluate pedestal profiles from DIII-D mdsplus database, where profiles
; are obtained via Osborne's python tools.
; In particular, evaluate a b-spline 
; Written: 03/22/12 - rjg

pro drvr_eval_bspline_CMOD

   ; Example data for 146394, 2250, 'ML02', 'vpolsplpsi'

   ; The fit was with a b-spline function. List spline knots and coefficients in
   ; arrays, exactly as obtained from the database.
   
   ;; Knot location vector
;   tt =  [ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.17300000,    $
;           0.40900000, 0.80000000, 0.96000000, 1.0000000, 1.0000000,  $
;           1.0000000, 1.0000000]
   
   ;; Coefficient array
;   coef = [ -0.0000000, -0.0025227561, -0.0084869598, -1.0037094, 1.3987460,  $
;            -18.679731, 14.745304, 12.000000 ]
   
  ;; Knot location vector for CMOD
   tt =  [0.696491,    .696491,    0.696491,    0.696491,    0.696491,    $
          0.783518,   0.818329,    0.853140,    0.887951,    0.922762,    $
          1.013270,    1.013270,   1.013270,    1.013270,    1.013270]

   ;; Coefficient array
   coef = [84.302376,   99.595856,   38.382866,   110.242325,  53.770931, $ 
           68.373230,   7.870729,    165.294556,  -315.809784, 206.881226] 


   ;; To evaluate the spline, we will need to make a call external. For this to
   ;; work, we expect to have some routines in a directory that we can find.
   ;; Need to set environmental variable to point to this.
   
   ; do NOT appear to need(???)
   setenv,'IDLUTILS_DIR=./idlutils'

   
   ; Set up a psi array where we evaluate the fit.  Will evaluate from psi = 0.8 - 1.0.
   ;  Need to be careful about going beyond psi=1.0.  Some fits are not done past
   ; psi=1.0.  You cannot tell that from the fit coefficients.
;   psi = 0.001 * indgen(201) + 0.8

   file_path = '/global/u1/e/emd/cmod_jrt2011/cmod.1080321020.01200/1080321020/1200/v02/grids/cmod.1080321020.01200.x516_y64.grd.nc'
   g=file_import(file_path)
   ;g=file_import("./bout_Cmod_EDA1110201023_t900_x516y128.grd.nc")

   psi=(g.psixy[*,32]-g.PSI_AXIS)/(g.PSI_BNDRY-g.PSI_AXIS)  

   yy =  slatec_bvalu (psi, tt, coef )
     
   ;; Plot the results
   Plot, psi, yy, psym=5
   
end




   
   
   
   
   
   

