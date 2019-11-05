;***************************
; write_polslice.pro
;***************************
;
; Created:  
;   03/26/2012 I. Joseph
;
; Modified:
;   09/19/2012 I. Joseph
;
;******************************

;***************************
; PRO write_polslice 
;***************************
;
; This function maps 3d data on a BOUT grid to 3d data at constant toroidal angle 
; evaluated at the Rxy and Zxy mesh points for each slice of the "poloidal" plane.
; This data is written to NETCDF output file, along with the toroidal harmonics 
; (1st sine and cosine, as of now), and some equilibrium information.
;
; Required Input:
;   g = BOUT grid data structure
;   f = 3D data (nx,ny,nz)
;   nzperiod = toroidal mode # of z-domain of 3D data
;   nzmode = toroidal harmonic of interest within z-domain; default = 1
;   these two inputs determine the toroidal mode number NMODE = nzperiod*nzmode 
;
; Optional Input:
;   nxout = resize to nxout radial points
;   nyout = resize to nyout poloidal points
;   file = name of NETCDF output file; default is 'bout.polslice.nc'
;   time = option to print time taken
;   memory = use less memory for output; only write harmonic of interest
;
; Optional output: 
;   data = returns calculated output
;           NX, NY, NZ, NMODE, RXY, ZXY, NE0, TE0, NI0, TI0 = equilibrium info
;           FZ = data mapped to constant toroidal plane
;           FC = cosine harmonics: nth index = nth harmonic (* nmode)
;           FS = sine harmonics: nth index = nth harmonic (* nmode)
;           FC1, FS1 = first harmonics: fundamental (nmode)
;   gout = returns resized BOUT grid structure
;   savefile = option to write IDL save file
;
; Needs:
;  .r  resize_bgrid ; in order to resize the grid
;
; Tips:
;  Increasing the poloidal grid size is very important for removing numerical artifacts (Moire patterns).
;   Artifacts are generated by the strong magnetic shearing of the data which occurs near X-points.
;
;
;******************************
; Created:  
;  03/26/2012 I. Joseph
;
; Modified:
;   09/19/2012 I. Joseph
;     - added full fourier harmonic output
;     - added timing calculation
;     - changed defaults & variable names
;     - added memory "safe" mode
;
;******************************



pro write_polslice, g, f, nzperiod=nzperiod, nzmode=nzmode, nxout=nxout, nyout=nyout, gout=gout, data=data, file=file, savefile=savefile, time=time, memory=memory
  if keyword_set(time) then time0=systime(/sec)

  if not keyword_set(nzperiod) then begin
    print, 'write_polslice: nzperiod must be specified'
    return
  endif

  d = size(f)
  nx = d[1]
  ny = d[2]
  nz = d[3] 
  if not keyword_set(nzmode) then begin
    nzmode = 1
  endif
  if nzmode ge nz/2 or nzmode lt 1 then begin
    print, 'write_polslice: requires 0 < nzmode < nz/2'
    return
  endif
  nmode = nzperiod*nzmode

  if not keyword_set(nxout) and not keyword_set(nyout) then begin
    Rxy = g.Rxy
    Zxy = g.Zxy
    Ne0 = g.Ni0
    Ni0 = g.Ni0
    Te0 = g.Te0
    Ti0 = g.Ti0

    f0 = f/max(abs(f))
    fz0 = zshift(f0, -g.zshift, period=nzperiod)
    fz0 = fz0/max(abs(fz0))

  endif else begin
    if not keyword_set(nxout) then nxout=g.nx
    if not keyword_set(nyout) then nyout=g.ny
    nx = nxout
    ny = nyout

    gout = resize_bgrid(g,nxout=nxout, nyout=nyout)
    Rxy = gout.Rxy
    Zxy = gout.Zxy
    Ne0 = gout.Ni0
    Te0 = gout.Te0
    Ni0 = gout.Ni0
    Ti0 = gout.Ti0

    f0 = resize_bgrid_var(g, f, nxout=nxout, nyout=nyout)  
    f0 = f0/max(abs(f0))
    fz0 = zshift(f0, -gout.zshift, period=nzperiod)
    fz0 = fz0/max(abs(fz0))


  endelse

; Calculate fourier harmonics
;   For cosine & sine terms: nth index = nth harmonic


if not keyword_set(memory) then begin
  fkz = fft(fz0,dim=3)
  f00 = reform(real_part(fkz[*,*,0]))
  f1c = reform(real_part(fkz[*,*,nz-nzmode]+fkz[*,*,nzmode]))
  f1s = reform(imaginary(fkz[*,*,nz-nzmode]-fkz[*,*,nzmode]))


  nz2=(nz+2)/2
  fc=fltarr(nx,ny,nz2) 
  fs=fc
  fc[*,*,0] = reform(real_part(fkz[*,*,0]))
  fs[*,*,0] = reform(imaginary(fkz[*,*,0]))
  fc[*,*,nz2-1] = reform(real_part(fkz[*,*,nz2-1]))
  fs[*,*,nz2-1] = reform(imaginary(fkz[*,*,nz2-1]))
  for iz = 1, nz2-1 do begin
    fc[*,*,iz] = reform(real_part(fkz[*,*,nz-iz]+fkz[*,*,iz]))
    fs[*,*,iz] = reform(imaginary(fkz[*,*,nz-iz]-fkz[*,*,iz]))
  endfor 

; Write data
  data = {NX:nx, NY:ny, NZ:nz, RXY:Rxy, ZXY:Zxy, $
     NE0:Ne0, TE0:Te0, NI0:Ni0, TI0:Ti0, $
	 NZPERIOD:nzperiod, NZMODE:nzmode, NMODE:nmode, $
     PZ:fz0, PC:fc, PS:fs, P00:f00, P1C:f1c, P1S:f1s}

endif else begin

  f00 = total(fz0,3)/float(nz)
  f1c = 0*f00
  f1s = f1c
  for iz=0,nz-1 do begin
    f1c=f1c+fz0[*,*,iz]*cos(2.*!pi*float(iz*nzmode)/float(nz))
    f1s=f1s+fz0[*,*,iz]*sin(2.*!pi*float(iz*nzmode)/float(nz))
  endfor
  f1c=f1c/float(nz/2.)
  f1s=f1s/float(nz/2.)

; Write data
  data = {NX:nx, NY:ny, NZ:nz, RXY:Rxy, ZXY:Zxy, $
     NE0:Ne0, TE0:Te0, NI0:Ni0, TI0:Ti0, $
	 NZPERIOD:nzperiod, NZMODE:nzmode, NMODE:nmode, $
     P00:f00, P1C:f1c, P1S:f1s}

endelse


  if not keyword_set(file) then begin
    file='bout.polslice_'+strtrim(nx,2)+'x'+strtrim(ny,2)+'.nc'
  endif
  print, 'write_polslice: writing to '+file
  s=file_export(file, data)

  if keyword_set(savefile) then begin
    save, nx,ny,nz,nzperiod, nzmode, nmode, $
          Rxy,Zxy,Ne0,Te0,Ni0,Ti0, $
          fz0, f00, f1c, f1s
  end

  if keyword_set(time) then begin
    time1=systime(/sec)
    print, 'write_polslice: finished in '+strtrim(time1-time0,2)+' sec'
  endif 

end
