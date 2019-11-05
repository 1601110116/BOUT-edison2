function fprime_smooth, ff, xstart=xstart, xsep=xsep

  ON_ERROR, 2

  pt = ff
  result = ff

  nd = size(ff,/n_dim)

  if not (nd eq 2) then begin
     print,"Error! Input should be 2D."
     return, 0
  endif

  nx = n_elements(pt[*,0])
  ny = n_elements(pt[0,*])

  if not keyword_set(xstart) then begin
     tmp = max(abs(pt[*,ny/4]), xt)
     xstart = xt
     if (xstart le 0) then begin
        xstart = nx-1
     endif
     print,"Set xstart to be the peak of the 1st derivative of input: ", xstart
  endif

  if not keyword_set(xsep) then begin
     xsep = (xstart + nx-1)/2
     print,"Set xsep to be: ", xsep
  endif

  for j=0, ny-1 do begin
;     dpt = deriv(pt[*,j])
     dpt = pt[*,j]
     xind = where (dpt[xstart:*] eq 0.)
;     print,"xind is ",xind
;     xrange1 = xind[0]+xstart-nx/52
     xrange1 = (xind[0]+2*xstart)/2
     if (xrange1 le xstart) then begin
        xrange1 = xstart
     endif
     xrange2 = xsep + nx/20
     if (xrange2 ge nx-1) then begin
        xrange2 = nx-1
     endif
     dpt2 = dpt
     dpt_t = dpt[xrange1-1]
     deltn = xrange2-xrange1+1
     ddpt = dpt_t/(deltn*(deltn+1)/2.)
;     ddpt = dpt_t/deltn

;     print,dpt_t,ddpt
     for i=xrange1, xrange2 do begin
;        dpt2[i] = dpt_t-ddpt*(i-xrange1)
        dpt2[i] = ddpt*(xrange2-i)*(xrange2-i+1)/2.
;        print,i
     endfor

     result[*,j] = dpt2
;     result[0,j] = ff[0,j]
;     for i=1,nx-1 do begin
;        result[i,j] = result[i-1,j]+dpt2[i]
;     endfor

;     plot,dpt[*],chars=1.5,thick=3
;     oplot,dpt2[*],thick=3,col=2   
  endfor

;  plot,deriv(ff[*,ny/2]),chars=1.5,thick=3
;  oplot,deriv(result[*,ny/2]),thick=3,col=2

  return, result
end
