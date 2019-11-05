; This function is used to fit the derivative of the profiles with
; Gaussian function.

function prof_fit, var, gname=gname, nterms=nterms
  
  ON_ERROR, 2
  
  nx = n_elements(var)
  
  if Not keyword_set(gname) then begin
     print,"No input grid, use index for x coordinate."
     rxy = findgen(nx)
  endif else begin
     if ( size(gname,/type) ne 7 ) then begin
        print, "Wrong input BOUT++ grid file."
        return,0
     endif else begin
        print, '*****************',gname
        g=file_import(gname)
        help,g,/str
        rmax=max(g.rxy[nx-1,*], ypeak)
        rxy = g.rxy[*,ypeak]
     endelse
  endelse
  
  if not keyword_set(nterms) then begin
     print,"Keyword nterms is not set. Use default value: 3."
     nterms=3
  endif else begin
     if (nterms ne 3) and  (nterms ne 4) and  (nterms ne 5) and (nterms ne 6) then begin
        print,"Keyword nterms is not in the range [3,6]. Use default value: 3."
        nterms=3
     endif
  endelse

  dvar = deriv(rxy,var)

  fit1 = gaussfit(rxy,dvar,coef1,nterms=nterms)
  
  print,"The coefficients for Gaussian function are:",coef1
  win,0
  plot, rxy, dvar, thick=3, chars=2, title='1st order derivative'
  oplot, rxy, fit1, col=2, thick=3

  opt = get_yesno("Is this fitting OK?")
  if opt ne 1 then begin
     coef2 = coef1
     min_dvar = min(dvar, x_min)
     case nterms of
        3: coef2[0] = min_dvar
        4: coef2[0] = min_dvar - coef1[3]
        5: coef2[0] = min_dvar - coef1[3] - rxy[x_min]*coef1[4]
        6: coef2[0] = min_dvar - coef1[3] - rxy[x_min]*coef1[4] - rxy[x_min]^2.*coef1[5]
     endcase
     coef2[2] = sqrt(coef1[0]/coef2[0])*coef1[2]
     rz = (rxy-coef2[1])/coef2[2]
;     print,nterms
     case nterms of
        6: fit1 = coef2[0]*exp(-0.5*rz^2.)+coef2[3]+coef2[4]*rxy+coef2[5]*rxy^2.
        5: fit1 = coef2[0]*exp(-0.5*rz^2.)+coef2[3]+coef2[4]*rxy
        4: fit1 = coef2[0]*exp(-0.5*rz^2.)+coef2[3]
        3: fit1 = coef2[0]*exp(-0.5*rz^2.)
     endcase
     print,"The coefficients for Gaussian function are:",coef2
     win,1
     plot, rxy, dvar, thick=3, chars=2, title='1st order derivative'
     oplot, rxy, fit1, col=2, thick=3

     opt = get_yesno("Is this fitting OK?")
     if opt ne 1 then begin
        print,"I quit!"
        return,0
     endif
;        repeat begin
;           for i = 0, nterms-1 do begin
;              print,"The coefficient",i+1," is:"
  endif
          
  varfit = fltarr(nx)
  for i=1,nx-1 do begin
     varfit[i] = varfit[i-1]+fit1[i]*(rxy[i]-rxy[i-1])
  endfor
  result = fltarr(nx)
  help,varfit
  result = varfit - varfit[nx-1] + var[nx-1]
  
  win,2
  plot, rxy, var, thick=3, chars=2, title='varible'
  oplot, rxy, result, col=2, thick=3
  
  return,result
end
