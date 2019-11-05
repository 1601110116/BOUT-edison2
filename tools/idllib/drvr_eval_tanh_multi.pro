
; Evaluate pedestal profiles from DIII-D mdsplus database, where profiles
; are obtained via Osborne's python tools.
; Written: 03/22/12 - rjg

pro drvr_eval_tanh_multi

   ; Example data for 146394, 2250, 'ML02', 'netanhpsi'
   ; The fit was with a tanh_multi function. List coefficients in an array,
   ;  exactly as obtained from the database.
   
   cc = [ 0.98309097, 0.050708881, 0.45031508, 0.037430270, -0.0016290952,  $
         0.0019897103, -3.4721258e-05 ]
   
   ; Set up a psi array where we evaluate the fit.  Will evaluate from psi = 0.8 - 1.0.
   ;  Need to be careful about going beyond psi=1.0.  Some fits are not done past
   ; psi=1.0.  You cannot tell that from the fit coefficients.
   psi = 0.001 * indgen(201) + 0.8
   
   
   ; Evaluate the coefficient array on the psi array
   param = 'none'
   yy = tanh_multi (cc, psi, param)
   
   
   ;; Plot the results
   Plot, psi, yy, psym=5
   
end




   
   
   
   
   
   

