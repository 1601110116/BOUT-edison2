;Import the experimentally measured temperature profiles, inarr, into the
;grid file, g. The input profile should be interpolated to the grid
;before this code.(ver 0.1). This function is only for single null geometry

function smoothtwo, inarr, xmin=xmin, xmax=xmax, width=width
  ON_ERROR, 2
  
  nx=n_elements(inarr)
  p=inarr
  ptemp=p[xmin:xmax]
  pt5=smooth(ptemp,width)
  pt5=smooth(pt5,5)
  for x=xmin,xmax do begin
     p[x]=pt5[x-xmin]
  endfor
  return,p
end

pro Ni2Gridall_im, filep0, filene, fileni, filete, fileti, fileer, filename;, smooth=smooth 
  
  ON_ERROR, 2

  ON_ERROR, 2


  if ( size(filep0,/type) ne 7 ) then begin
     print, "Wrong pressure input file."
     return
  endif

  if ( size(filene,/type) ne 7 ) then begin
     print, "Wrong electron density input file."
     return
  endif

  if ( size(fileni,/type) ne 7 ) then begin
     print, "Wrong ion density input file."
     return
  endif

  if ( size(filete,/type) ne 7 ) then begin
     print, "Wrong electron temperature input file."
     return
  endif

  if ( size(fileti,/type) ne 7 ) then begin
     print, "Wrong ion temperature input file."
     return
  endif

  if ( size(fileer,/type) ne 7 ) then begin
     print, "Wrong radial electric field input file."
     return
  endif

  if ( size(filename,/type) ne 7 ) then begin
     print, "Wrong grid file."
     return
  endif 

  g=file_import(filename)
  rm_temp=max(g.rxy[g.nx-1,*],ypeak)
  psn=(g.psixy[*,ypeak]-g.psi_axis)/(g.psi_bndry-g.psi_axis)

  restore,filep0
  restore,filene
  restore,fileni
  restore,filete
  restore,fileti
  restore,fileer  
;  help

  opt = get_yesno("Is the 6th file Er?(y/n):")
  if opt ne 1 then begin
     print,"Using Omgeb instead of Er."
  endif

  nx = g.nx
  ny = g.ny

  kb = 1.38e-23
;  profp0 = fltarr(nx,ny)
  profte = fltarr(nx,ny)
  profne = fltarr(nx,ny)
  profni = fltarr(nx,ny)
  profti = fltarr(nx,ny)
  profer = fltarr(nx,ny)
  profp = fltarr(nx,ny)
  profphi = fltarr(nx,ny)
  xsep = floor(g.ixseps1)
  if (xsep lt 0.) then begin
     xsep = nx-1
  endif
  if (xsep ge nx) then begin
     xsep = nx-1
  endif
  xsep2 = floor(g.ixseps2)
  if (xsep2 lt 0.) then begin
     xsep2 = nx-1
  endif
  if (xsep2 ge nx) then begin
     xsep2 = nx-1
  endif

  phi0 = fltarr(nx)

;  profp = p0
  tmp = where(te gt 0.01)
  nte = n_elements(tmp) 
  tet = te
  for i=nte, nx-1 do begin
     tet[i] = te[nte-1]
  endfor

  if opt eq 1 then begin
     for i=0, nx-1 do begin
        phi0[i] = total(er[0:i]*g.dx[0:i,ypeak])
     endfor
  endif else begin
     print,'Undevelpoed phi0 for Omgeb!'
  endelse

  for i=0, xsep do begin
     for j=0, ny-1 do begin
        if ((j gt g.jyseps1_1) and (j le g.jyseps2_1)) or ((j gt g.jyseps1_2) and (j le g.jyseps2_2)) then begin
           profte[i,j] = tet[i]
           profti[i,j] = ti[i]
           profni[i,j] = ni[i]
           profne[i,j] = nee[i]
           profp[i,j] = p0[i]
           if opt ne 1 then begin
              profer[i,j] = om[i]*g.rxy[i,j]*g.bpxy[i,j]*1000.
           endif else begin
              profer[i,j] = er[i]*g.rxy[i,j]*g.bpxy[i,j]/(g.rxy[i,ypeak]*g.bpxy[i,ypeak])*1000.
           endelse
           profphi[i,j] = phi0[i]*1000.
        endif else begin
           if (j le g.jyseps1_1) or (j gt g.jyseps2_2) then begin
              profte[i,j] = tet[xsep]
              profti[i,j] = ti[xsep]
              profni[i,j] = ni[xsep]
              profne[i,j] = nee[xsep]
              profp[i,j] = p0[xsep]
              if opt ne 1 then begin
                 profer[i,j] = om[xsep]*g.rxy[xsep,j]*g.bpxy[xsep,j]*1000.
              endif else begin
                 profer[i,j] = er[xsep]*g.rxy[xsep,j]*g.bpxy[xsep,j]/(g.rxy[xsep,ypeak]*g.bpxy[xsep,ypeak])*1000.
              endelse
              profphi[i,j] = phi0[xsep]*1000.
           endif else begin
              profte[i,j] = tet[xsep2]
              profti[i,j] = ti[xsep2]
              profni[i,j] = ni[xsep2]
              profne[i,j] = nee[xsep2]
              profp[i,j] = p0[xsep2]
              if opt ne 1 then begin
                 profer[i,j] = om[xsep2]*g.rxy[xsep2,j]*g.bpxy[xsep2,j]*1000.
              endif else begin
                 profer[i,j] = er[xsep2]*g.rxy[xsep2,j]*g.bpxy[xsep2,j]/(g.rxy[xsep2,ypeak]*g.bpxy[xsep2,ypeak])*1000.
              endelse
              profphi[i,j] = phi0[xsep2]*1000. 
           endelse
        endelse
     endfor
  endfor

  for i=xsep, nx-1 do begin
     for j=0, ny-1 do begin
        profte[i,j] = tet[i]
        profti[i,j] = ti[i]
        profni[i,j] = ni[i]
        profne[i,j] = nee[i]
        profp[i,j] = p0[i]
        if opt ne 1 then begin
           profer[i,j] = om[i]*g.rxy[i,j]*g.bpxy[i,j]*1000.
        endif else begin
           profer[i,j] = er[i]*g.rxy[i,j]*g.bpxy[i,j]/(g.rxy[i,ypeak]*g.bpxy[i,ypeak])*1000.
        endelse
        profphi[i,j] = phi0[i]*1000.
     endfor
  endfor

;  profni = profp/(profte+profti)/16.02
;  window,0
;  plot,profti[309,*]
;  window,1
;  plot,profte[309,*]

;hel, profni
  for i=xsep, nx-1 do begin
     for j=0, g.jyseps1_1 do begin
        profte[i,j] = profte[i,g.jyseps1_1+1]
        profti[i,j] = profti[i,g.jyseps1_1+1]
        profni[i,j] = profni[i,g.jyseps1_1+1]
        profne[i,j] = profne[i,g.jyseps1_1+1]
        profp[i,j] = profp[i,g.jyseps1_1+1]
        profer[i,j] = profer[i,g.jyseps1_1+1]
        profphi[i,j] = profphi[i,g.jyseps1_1+1]
     endfor
     for j = g.jyseps1_2+1, g.ny_inner-1 do begin
        profte[i,j] = profte[i,g.jyseps1_2]
        profti[i,j] = profti[i,g.jyseps1_2]
        profni[i,j] = profni[i,g.jyseps1_2]
        profne[i,j] = profne[i,g.jyseps1_2]
        profp[i,j] = profp[i,g.jyseps1_2]
        profer[i,j] = profer[i,g.jyseps1_2]
        profphi[i,j] = profphi[i,g.jyseps1_2]
     endfor
     for j = g.ny_inner, g.jyseps2_1 do begin
        profte[i,j] = profte[i,g.jyseps2_1+1]
        profti[i,j] = profti[i,g.jyseps2_1+1]
        profni[i,j] = profni[i,g.jyseps2_1+1]
        profne[i,j] = profne[i,g.jyseps2_1+1]
        profp[i,j] = profp[i,g.jyseps2_1+1]
        profer[i,j] = profer[i,g.jyseps2_1+1]
        profphi[i,j] = profphi[i,g.jyseps2_1+1]
     endfor
     for j = g.jyseps2_2+1, ny-1 do begin
        profte[i,j] = profte[i,g.jyseps2_2]
        profti[i,j] = profti[i,g.jyseps2_2]
        profni[i,j] = profni[i,g.jyseps2_2]
        profne[i,j] = profne[i,g.jyseps2_2]
        profp[i,j] = profp[i,g.jyseps2_2]
        profer[i,j] = profer[i,g.jyseps2_2]
        profphi[i,j] = profphi[i,g.jyseps2_2]
     endfor
  endfor

  for i=0, xsep-1 do begin
     for j=0, g.jyseps1_1 do begin
        profte[i,j] = profte[xsep,j]
        profti[i,j] = profti[xsep,j]
        profni[i,j] = profni[xsep,j]
        profne[i,j] = profne[xsep,j]
        profp[i,j] = profp[xsep,j]
        profer[i,j] = profer[xsep,j]
        profphi[i,j] = profphi[xsep,j]
     endfor
     for j = g.jyseps1_2+1, g.jyseps2_1 do begin
        profte[i,j] = profte[xsep,j]
        profti[i,j] = profti[xsep,j]
        profni[i,j] = profni[xsep,j]
        profne[i,j] = profne[xsep,j]
        profp[i,j] = profp[xsep,j]
        profer[i,j] = profer[xsep,j]
        profphi[i,j] = profphi[xsep,j]
     endfor
     for j = g.jyseps2_2+1, ny-1 do begin
        profte[i,j] = profte[xsep,j]
        profti[i,j] = profti[xsep,j]
        profni[i,j] = profni[xsep,j]
        profne[i,j] = profne[xsep,j]
        profp[i,j] = profp[xsep,j]
        profer[i,j] = profer[xsep,j]
        profphi[i,j] = profphi[xsep,j]
     endfor
  endfor  

;  if keyword_set(smooth) then begin
;     xmin1=xsep-20
;     xmax1=xsep+20
;     for j=0,ny-1 do begin
;        profp[*,j]=smoothtwo(profp[*,j],xmin=xmin1,xmax=xmax1,width=20)
;        profni[*,j]=smoothtwo(profni[*,j],xmin=xmin1,xmax=xmax1,width=20)
;        profte[*,j]=smoothtwo(profte[*,j],xmin=xmin1,xmax=xmax1,width=20)
;        profer[*,j]=smoothtwo(profer[*,j],xmin=xmin1,xmax=xmax1,width=20)
;        profphi[*,j]=smoothtwo(profphi[*,j],xmin=xmin1,xmax=xmax1,width=20)
;        profti[*,j]=smoothtwo(profti[*,j],xmin=xmin1,xmax=xmax1,width=20)
;     endfor
;  endif

;  f=file_export(filename, gnew)
  handle = file_open(filename, /write)
  s = file_write(handle, 'Niexp', profni)
  s = file_write(handle, 'Neexp', profne)
  s = file_write(handle, 'Tiexp', profti)
  s = file_write(handle, 'Teexp', profte)
  s = file_write(handle, 'pressure_s', profp)
  s = file_write(handle, 'E_r', profer)
  s = file_write(handle, 'Phi_0', profphi)
  s = file_write(handle, 'Nixexp', profni[0,ypeak])
  file_close, handle

end

  
