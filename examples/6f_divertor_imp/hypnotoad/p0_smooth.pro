function p0_smooth, p0, xstart=xstart, xend=xend

  ON_ERROR, 2

  pt = p0

  nd = size(p0,/n_dim)

  if ( (nd gt 2) or (nd eq 0) ) then begin
     print,"Error! P0 should be 1D or 2D."
     return, 0
  endif
  
  if (nd eq 1) then begin
     nx = n_elements(p0)
  endif else begin
     nx = n_elements(p0[*,0])
     ny = n_elements(p0[0,*])
     maxdp = 0.
     yy = 0
     for jy=0, ny-1 do begin
        dptmp = max(-deriv(p0[*,jy]),xt)
        if (dptmp gt maxdp) then begin
           maxdp = dptmp
           yy = jy
;           help,yy
        endif
     endfor
  endelse
  
  if not keyword_set(xstart) then begin
;  xstart = get_integer("please enter the index of xstart:")
     if (nd eq 1) then begin
        tmp = max(abs(deriv(deriv(p0))), xt)
     endif else begin
        tmp = max(abs(deriv(deriv(p0[*,yy]))), xt)
     endelse
     xstart = xt-nx/50
     print,"Set xstart to be the peak of the 2nd derivative of P0: ", xstart 
  endif

  if not keyword_set(xend) then begin
;  xend = get_integer("please enter the index of xend:")
     xend = nx-1
     print, "Set xend to be the last index:", xend
  endif
  
  if (xend gt nx-1) then begin
     print,"Warning: xend is larger than the length of input array."
     xend = nx-1
     print, "Set xend to be the last index:", xend
  endif

  if (nd eq 1) then begin
     dp = deriv(p0)
     dpt=dp[xstart-1]*.9
     for i=xstart,xend do begin
        pt[i] = pt[i-1] + dpt
        dpt = dpt*.9
     endfor
  endif else begin
     for jy=0, ny-1 do begin
        dp = deriv(p0[*,jy])
        dpt = dp[xstart-1]*.9
        for i=xstart,xend do begin
           pt[i,jy] = pt[i-1,jy] + dpt
           dpt = dpt*.9
        endfor
     endfor
  endelse

  return, pt
end
