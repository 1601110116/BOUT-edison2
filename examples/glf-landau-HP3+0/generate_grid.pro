; Create an input file for 2D diffusion equation example
;=======================================================;

@pdb2idl

PRO generate_grid, $
nx=nx, ny=ny, rxy=rxy, zxy=zxy, metric=metric, core=core,$
file=file, save=save, plot=plot, debug=debug

;
;
;

  IF NOT KEYWORD_SET(nx) THEN nx = 10
  IF NOT KEYWORD_SET(ny) THEN ny = 64
  IF NOT KEYWORD_SET(file) THEN file="slab.grd.nc"


  ;;-magnetic field
  Bp0=1e0 ;-poloidal field [T]
  Bt0=1e1 ;-toroidal field [T] at reference location Rt0
  Rt0=1.0 ;-reference location [m] for Bt0
  eCharge=1.6e-19 ;-electron charge [MKS]


  ;;-domain boundaries
  Rmin=Rt0-0.01 ;-left boundary  [m]
  Rmax=Rt0+0.01 ;-right boundary [m]
  Zmin=0.00 ;-bottom boundary [m]
  Zmax=0.01 ;-top boundary [m]

  dR=(Rmax-Rmin)/(nx-1)
  dZ=(Zmax-Zmin)/(ny-1)

  hthe0=(Zmax-Zmin)/(2*!PI)
  hthe = FLTARR(nx, ny) + hthe0



  ;;-actual coordinates are needed only for plotting
  Rxy=fltarr(nx,ny)
  Zxy=fltarr(nx,ny)
  Bxy=fltarr(nx,ny)
  Bpxy=fltarr(nx,ny)
  Btxy=fltarr(nx,ny)
  Psixy=fltarr(nx,ny)
  dPsi=fltarr(nx,ny)

  dx=fltarr(nx,ny)
  dy=fltarr(nx,ny)


  ;;-curvature components
  bxcvx=fltarr(nx,ny)
  bxcvy=fltarr(nx,ny)
  bxcvz=fltarr(nx,ny)


  ;;-shear quantities
  sinty=fltarr(nx,ny);;+1.
  qinty=fltarr(nx,ny);;+1.



  ;;-set geometry and magnetic field on the grid
  for ix=0,nx-1 do begin
      for jy=0,ny-1 do begin

          Rxy[ix,jy]=Rmin+dR*ix
          Zxy[ix,jy]=Zmin+dZ*jy

          Bpxy[ix,jy]=Bp0 ;;-constant value
          Btxy[ix,jy]=Bt0*Rt0/Rxy[ix,jy]
          Bxy[ix,jy]=SQRT((Bpxy[ix,jy])^2+(Btxy[ix,jy])^2)

          Psixy[ix,jy]=0.5*Bp0*Rxy[ix,jy]^2
          dPsi[ix,jy]=Bp0*Rxy[ix,jy]*dR
          dy[ix,jy] = 2.0*!PI/ny

          bxcvz[ix,jy] = (Btxy[ix,jy]/Bxy[ix,jy])^2*(Bxy[ix,jy]/Bpxy[ix,jy])/(Rxy[ix,jy])^2

      endfor
  endfor



  ;;-set background 2D profiles
  Te0 = FLTARR(nx, ny)
  Ti0 = FLTARR(nx, ny)
  Ni0 = FLTARR(nx, ny)
  Vi0 = FLTARR(nx, ny)
  Ve0 = FLTARR(nx, ny)
  phi0= FLTARR(nx, ny)
  rho0= FLTARR(nx, ny)
  Ajpar0= FLTARR(nx, ny)


  Ln=1e10 ;;-density decay length [m]

  for ix=0,nx-1 do begin
      for jy=0,ny-1 do begin

          Te0[ix,jy] = 1e0 ;;-[eV]
          Ti0[ix,jy] = 1e-10 ;;-[eV]

          Vi0[ix,jy] = 0e0 ;;-[m/s]
          Ve0[ix,jy] = 0e0 ;;-[m/s]
          
          Ni0[ix,jy]  = 1e20*exp(-(Rxy[ix,jy]-Rmin)/Ln) ;;-[m-3]

          phi0[ix,jy] = 0e0 ;;-[V]
          rho0[ix,jy] = 0e0 ;;-[?]

          Ajpar0[ix,jy] = 0e0 ;;-[?]

      endfor
  endfor


  ;;-set elements of metric tensor
  g11=fltarr(nx,ny)
  g22=fltarr(nx,ny)
  g33=fltarr(nx,ny)
  g12=fltarr(nx,ny)
  g13=fltarr(nx,ny)
  g23=fltarr(nx,ny)

  g11[*,*] = 1.0
  g22[*,*] = 1.0
  g33[*,*] = 1.0
  g12[*,*] = 0.0
  g13[*,*] = 0.0
  g23[*,*] = 0.0
  

  if keyword_set(CORE) then begin
      ;; entire domain inside 'core' - periodic
      ixseps1 = nx
      ixseps2 = nx
      jyseps1_1 = -1
      jyseps2_2 = ny-1
  endif else begin
      ;; Topology: Set all points outside separatrix
      ;; so NOT periodic in Y
      ixseps1 = 0
      ixseps2 = 0
      jyseps1_1 = -1
      jyseps2_2 = ny-1
  endelse
  jyseps1_2=ny/2
  jyseps2_1=ny/2


  if keyword_set(SAVE) then begin

      print, "Writing data to ", file

      f = file_open(file, /create)

      status = file_write(f, "nx", nx)
      status = file_write(f, "ny", ny)


      status = file_write(f, "Rxy", Rxy)
      status = file_write(f, "Zxy", Zxy)

      status = file_write(f, "Bxy", Bxy)
      status = file_write(f, "Bpxy", Bpxy)
      status = file_write(f, "Btxy", Btxy)

      status = file_write(f, "bxcvx", bxcvx)
      status = file_write(f, "bxcvy", bxcvy)
      status = file_write(f, "bxcvz", bxcvz)

      ;;-supply metric tensor if not calculating internally
      ;;status = file_write(f, "g11", g11)
      ;;status = file_write(f, "g22", g22)
      ;;status = file_write(f, "g33", g33)
      ;;status = file_write(f, "g12", g12)
      ;;status = file_write(f, "g13", g13)
      ;;status = file_write(f, "g23", g23)

      status = file_write(f, "sinty", sinty)
      status = file_write(f, "qinty", qinty)      

      status = file_write(f, "ixseps1", ixseps1)
      status = file_write(f, "ixseps2", ixseps2)
      status = file_write(f, "jyseps1_1", jyseps1_1)
      status = file_write(f, "jyseps2_2", jyseps2_2)
      status = file_write(f, "jyseps2_1", jyseps2_1)
      status = file_write(f, "jyseps1_2", jyseps1_2)

      status = file_write(f, "Te0", Te0)
      status = file_write(f, "Ti0", Ti0)
      status = file_write(f, "Ni0", Ni0)
      status = file_write(f, "Vi0", Vi0)
      status = file_write(f, "Ve0", Ve0)
      status = file_write(f, "Ajpar0", Ajpar0)
      status = file_write(f, "phi0", phi0)
      status = file_write(f, "rho0", rho0)


      ;;-some auxiliary quantities
      status = file_write(f, "Ni_x", max(Ni0))
      status = file_write(f, "Vi_x", 1.0)
      status = file_write(f, "Te_x", max(Te0))
      status = file_write(f, "Ti_x", max(Ti0))
      status = file_write(f, "bmag", max(Bxy))


      status = file_write(f, "hthe", hthe)
      status = file_write(f, "hthe0", hthe0)
      status = file_write(f, "dpsi", dpsi)
      status = file_write(f, "dy", dy)


      gjy0=0
      status = file_write(f,"gjy0",gjy0)

      INIXNORM=0
      status = file_write(f,"INIXNORM",INIXNORM)

      IXLB2=34
      status = file_write(f,"IXLB2",IXLB2)

      JNIXNORM=0
      status = file_write(f,"JNIXNORM",JNIXNORM)


      jpar0=fltarr(nx,ny)+1.0
      status = file_write(f, "jpar0", jpar0)

      file_close, f

  endif


;
;
;
if keyword_set(DEBUG) then STOP
END
