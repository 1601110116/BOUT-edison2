scriptID("$Id: tokamak.bas,v 1.6 2007/08/03 19:30:47 pearls Exp $")

if (~exists("graphics_bas_id")) then
  read graphics.bas
endif

# Allow debug stops...
mdef RUN =
  run
  if (debug == yes) then
    stdout << return << "Stopping at RUN #" << nprob << "..."
    pb
    call debugger
  endif
mend

function tokamak(;arg) # Display list of tokamak.bas routines

  default(arg) = no

  stdout << return \
         << "User-callable functions defined in tokamak.bas..." << return
  call basisexe("grep '^function' "//tokamak_bas_pn)
  stdout << " "

  if (arg <> no) then
    stdout << "Internal functions defined in tokamak.bas..." << return
    call basisexe("grep '^  function' "//tokamak_bas_pn)
    stdout << " "
  endif

endf # tokamak

function tokamak_ds(;fname1,fname2) # Dead-start a tokamak equilibrium

  default(fname1) = "./tokamak.inp"
  default(fname2) = "./pfcoil.inp"

  if (fname1 == "help") then
    remark " "
    remark "Function ""tokamak_ds(;file_name,coil_specs)"" (or its alias, ""ds"") creates a"
    remark "free-boundary tokamak equilibrim from parameters in the file named by argument"
    remark "file_name [default: ""tokamak.inp""]."
    remark " "
    remark "PF coil specifications are read from a file named by the 2nd argument"
    remark "[default: ""pfcoil.inp""]. If this file does not exist, then the coil"
    remark "configuration will be generated."
    remark " "
    remark "The following sample input files may be copied from ~bulmer/Corsica/Share and"
    remark "used as a starting point:"
    remark " "
    remark "    lim_tokamak.inp   # limited tokamak"
    remark "    dn_tokamak.inp    # diverted (DN) tokamak"
    remark "    kstar_tokamak.inp # creates nominal KSTAR equilibrium"
    remark "    kstar_pfcoil.inp  # for use with kstar_tokamak.inp"
    remark " "
    remark "Limited configuration example..."
    remark " "
    remark "    caltrans tokamak.bas"
    remark "    ds(""lim_tokamak.inp"")"
    remark " "
    remark "Diverted example..."
    remark " "
    remark "    caltrans tokamak.bas"
    remark "    ds(""dn_tokamak.inp"")"
    remark " "
    remark "KSTAR example with ""real"" coil set..."
    remark " "
    remark "    caltrans tokamak.bas"
    remark "    ds(""kstar_tokamak.inp"",""kstar_pfcoil.inp"")"
    remark " "
    remark "Upon successful completion of the dead-start procedure, a free-boundary"
    remark "equilibrium savefile (.sav) file be created."
    remark " "
    return
  endif

  real r_major, r_minor, z_axis, kappa_95, delta_95
  real delta_sep, beta_p, l_i, psi_ext
  real rtor, r_grid1, r_grid(2), z_grid(2)

  if (nprob > 0) then
    warning("*** Equilibrium already exists: dead start needs new session ***")
    return
  endif

  # Read plasma specifications from disk...
  integer io = basopen(fname1, "i")
  if (io <> -1) then
    io = basopen(fname1, "r")
    noisy = yes
    chameleon plasma_id
    io >> plasma_id >> return
    probid = plasma_id // " Starting Case"
    noisy = no
    io >> plcm         >> return
    io >> r_major      >> return
    io >> r_minor      >> return
    io >> z_axis       >> return
    io >> kappa_95     >> return
    io >> delta_95     >> return
    io >> delta_sep    >> return
    io >> beta_p       >> return
    io >> l_i          >> return
    io >> psi_ext      >> return
    io >> btor >> rtor >> return
    io >> r_grid1
    if (eof <> 1) then
      r_grid(1) = r_grid1
      io >> r_grid(2)
      io >> z_grid       >> return
      io >> jm >> km     >> return
      io >> rclmin       >> return
      io >> rclmax       >> return
      io >> zclmin       >> return
      io >> zclmax       >> return
      io >> crdn         >> return
    else
      r_grid(1) = rtor - 1.2*r_minor
      r_grid(2) = rtor + 1.2*r_minor
      z_grid(2) = 0.5*kappa_95*(r_grid(2) - r_grid(1))
      z_grid(1) = -z_grid(2)
      jm = 33
      km = 65
    endif
    call basclose(io)
  else
    if (nprob > 0) then # An equilibrium must already exist
      return
    else
      warning("*** Cannot open input file: "//trim(fname1)//" ***")
      return (1)
    endif
  endif

  # Construct name for save-file...
  global chameleon saveName = tolower(trim(plasma_id))
  character*1 nextc
  integer i
  do i = 2,strlen(trim(plasma_id))
    if (substr(trim(plasma_id),i,1) == " " |
        substr(trim(plasma_id),i,1) == "/" ) then
      nextc = "_"
    else
      nextc = tolower(substr(trim(plasma_id),i,1))
    endif
    saveName = substr(saveName,1,i-1) // nextc
  enddo
  saveName = saveName // ".sav"

  # Convert from S.I. to "Corsica" units...
  r_major   = r_major * 100
  r_minor   = r_minor * 100
  z_axis    = z_axis * 100
  btor      = btor * 1.0e+04
  rtor      = rtor * 100
  r_grid    = r_grid * 100
  z_grid    = z_grid * 100
  rclmin    = rclmin * 100
  rclmax    = rclmax * 100
  zclmin    = zclmin * 100
  zclmax    = zclmax * 100
  crdn      = crdn * 1000

  # Set-up the computational grid parameters...

  ngp=0
  liml = (jm - 1) / 8 + 1
  rn = r_grid(1)
  rx = r_grid(2)
  zx = z_grid(2)
  zn = z_grid(1)

  if (r_grid(1) <= 0) then
    r_grid(1) = errt
    r_grid(2) = r_grid(2) + errt
  endif
  btor=btor*rtor/ro
  ro = 0.5 * (r_grid(1) + r_grid(2))
  if (abs(rtor - ro)/ro > 0.001) then
    chameleon query = "Center-of-grid ""ro"" not same as ""rtor"", continue? "
    if(substr(sbasget(query),1,1) <> "y") then
      call kaboom(0)
    endif
  endif
  dr = float(r_grid(2) - r_grid(1)) / (jm - 1)
  rb = (dr / 1.0017469212) * ((jm + 1) / 2 - liml)
  rfix = 0
  el = float(z_grid(2) - z_grid(1)) / float(r_grid(2) - r_grid(1))
  rlim(0) = ro - rb
  limlo = 0
  # Set analytic coil parameters...
  if (rbcoil <= 0) then
    rbcoil = (3*ro + rb)/4
  endif
  if (rocoil <= 0) then
    rocoil = (5*ro - rb)/4
  endif
  if (rbcoil > 3*rb) then
    rocoil = ro
    rbcoil = 3*rb
  endif
  if (elcoil <= 0) then
    elcoil = max(1,0.9*kappa_95)
    dcoil  = delta_95
  endif
  if (ncoil == 0) then
    ncoil  = 64
  endif

  # Load the tokamak configuration...
  call read_pfcoil(fname2)

  # Set plot parameters...
  if (rclmax == 0) then
    rclmin = 0
    rclmax = int(100*(max(rc) + 2*max(drc)))
    rclmax = 10*int(0.1*rclmax + 1)
    zclmax = int(100*(max(zc) + 2*max(dzc)))
    zclmax = 10*int(0.1*zclmax + 1)
    zclmin = -zclmax
    crdn   = 2000
  endif

  # X-point search box...
  rxpr(1) = max(ro-rb,100*(min(rc)+2*min(drc)))
  rxpr(2) = min(ro+rb,100*(max(rc)+2*min(drc)))
  zxpr(2) = min(z_grid(2),100*(max(zc)-3*min(drc)))
  zxpr(1) = -zxpr(2)

  # Set markers for plasma boundary...
  nsym  = 2              # For plasma up/down asymmetry
  ibdry = 1              # Use analytic boundary generator initially
  nbd   = min(7,odd(nc)) # No. "hard" marker points, but
  ebdu  = kappa_95       #   use user's shape parameters
  dbdu  = delta_95
  ebdl  = ebdu
  dbdl  = dbdu

  # Get the PF coil specifications...
  nsymc = 2         # Assume coil up/down asymmetry
  ircwt = 3         # Use J**2 minimization
  cc = -1           # Initialize coil currents
  ic = iota(nc)     # Use all circuits independently
  vltf = 0          # No constraint on external flux linkage until later

  # Simple profiles, to start...
  betaj = 0.1
  ipp = 1; alfa(1) = -1
  ipf = 1; alfa(0) = -1

  # Get initial equilibrium...
  ixpt     = 1              # Either limited or diverted
  keqic    = 1              # Start from analytic psi
  psix     = plcm * 1.0e+08 # Initial guess---may need to adjust
  nl       = 200            # Iteration limit on GS
  ni       = 1              # Inner-loop iteration limit (not needed)
  map      = 64             # 1/2 No. points in poloidal direction
  alphac   = 0.5            # Relaxation factor
  mls      = 2 * map +1     # Number of points around surface.
  note = "Initial equilibrium with simple profiles"
  chameleon nfbd_file = nfbd
  if (nfbd_file > 0) then
    chameleon rfbd_file = rfbd
    chameleon zfbd_file = zfbd
  endif
  nfbd = 0
  call gchange("Teq_input",0)
  call gchange("TeqGS_input",0)
  teq_interface
  separatrix_interface
  teqGS_interface
  RUN
  keqic = 0 # From now on, use last psi as initial guess
  ibdry = 0 # From now on, user specifies r,zbd points

  note = "with fuzzy markers"
  nbd = 0
  if (nfbd_file == 0) then # Analytic description
    call dshape(64,r_major,z_axis,r_minor,kappa_95,delta_95)
    rlim(0) = rfbd(1)
    zlim(0) = zfbd(1)
  else # Assume specified shape doesn't need limiter
    nfbd = nfbd_file
    rfbd = rfbd_file
    zfbd = zfbd_file
    call set_limiter_pt()
  endif
  alfbd = (pi * 1.0e+09) / (plcm * ro**2) # Empirical scaling
  alfbd = 2**(int(log(alfbd(1)) / log(2.0) + 1))
  RUN

  # Add separatrix-separation constraint...
  if (limiterd == 0 | delta_sep <> 0 & abs(delta_sep) < 0.001) then
    nbd = 1
    rbd(1) = min(rfbd)
    zbd(1) = 0
    note = "with delta_sep constraint for ~DN"
    if (limiterd ==1) then
      if(rlim(0) > rcntr) then
        rlim(0)=rlim(0)+dr
      else
        rlim(0)=rlim(0)-dr
      endif
    endif
    irl = 0
    if (delta_sep == 0) then
      rl = -1.0e-06
    else
      rl = sign(1.0e-06,delta_sep)
    endif
    RUN
  endif

  # Apply desired separation for SN configurations...

  if (abs(delta_sep) >= 0.001) then
    if (limiterd ==1) then
      if(rlim(0) > rcntr) then
        rlim(0)=rlim(0)+dr
      else
        rlim(0)=rlim(0)-dr
      endif
    endif
    rlimm=rlim;zlimm=zlim
    separatrix
#    rfbd=rls(::2);zfbd=zls(::2)
#    note = "moving fuzzy boundary point to separatrix"
    RUN
    if (abs(delta_sep) > 0) then
      note = "increasing separatrix separation"
      while (dsep < 10 * abs(delta_sep))
        rl = 10 * rl
        run
      endwhile
      if (delta_sep > 0) then # SNT
        note = "with delta_sep constraint for SNT"
        irl = -1
      elseif (delta_sep < 0) then # SNB
        note = "with delta_sep constraint for SNB"
        irl = 1
      endif
      rl = 100 * abs(delta_sep)
      RUN
      rfbd=rls(::2);zfbd=zls(::2)
    endif
  endif

  # Switch to "Ohmic" profiles...
  
  rl=0                      # Do not constrain distance between separatrices
  zxpr(2)=-1                # Avoid upper X-points
  aap = 0                   # No flat portion in oOhmic profile
  msrf = 51                 # No. flux surfaces
  lsrf = -1                 # Expand near axis & edge
  ipj = 2                   # Now specifying <J.B>, not FF'
  ipp = 3                   # Use parabolic profile form for p
  alfa(1) = 2; betp(1) = 1  # Pressure profile exponents
  ipf = 3                   # Use parabolic profile form for <J.B>
  alfa(0) = 2; betp(0) = 2  # Parallel current profile exponents
  betaj = 0.3
  note = "with Ohmic profiles"

  RUN
  # Converge to nominal profile parameters with HYBRD (ceq package)...
  # NOTE: li(3) is the ITER definition, li(1) is the GA definition
  note = "HYBRD solution for desired profiles"
  package ceq
  nbd = 0
  vltf = psi_ext
  epsfcn = 10 * epsj; tol = 10 * epsfcn; factor = 0.1
  nctot = 2
  vo  = ["li(3)","betap(1)"]
  vo0 = [l_i, beta_p]
  vi  = ["betp(0)","betaj"]
  x0  = [betp(0),betaj]
  if (debug == yes) then
    ihy = 0; run
  else
    ihy = 99; run; ihy = 0; note = " "
    saveq(saveName)
    stdout << return << "Saved to disk in `"//saveName//"'"
  endif
endf # tokamak_ds

function ds(;fname1,fname2) # Synonym for tokamak_ds

  default(fname1) = "./tokamak.inp"
  default(fname2) = "./pfcoil.inp"

  call tokamak_ds(fname1,fname2)

endf # ds

function tokamak_rd # (Re-) read tokamak configuration files

  call read_pfcoil
  call read_tfcoil
  call read_shape
  call read_fwall
  call read_passive

endf # tokamak_rd

function set_limiter # Map limiter coordinates to Corsica's rlimw,zlimw arrays

  integer i, k, n
  real f, ds, dist = 0

  if (lmax <> 0) then
    do i = idx_lim(1),idx_lim(2)-1
      dist = dist + sqrt( (rfw(i+1) - rfw(i))**2 + (zfw(i+1) - zfw(i))**2 )
    enddo
    ds = dist/(lmax - 1)
    rlimw = 0
    zlimw = 0
    integer nlimw = 1
    rlimw(nlimw) = rfw(idx_lim(1))
    zlimw(nlimw) = zfw(idx_lim(1))
    do i = idx_lim(1),idx_lim(2)-1
      n = sqrt( (rfw(i+1) - rfw(i))**2 + (zfw(i+1) - zfw(i))**2) / ds + 2
      do k = 1,n
        f = float(k)/float(n)
        if (nlimw < lmax) then
          nlimw = nlimw + 1
          rlimw(nlimw) = rfw(i) + f * (rfw(i+1) - rfw(i))
          zlimw(nlimw) = zfw(i) + f * (zfw(i+1) - zfw(i))
        endif
      enddo
      if (i < (idx_lim(2) - 1)) then
        dist = dist - sqrt( (rfw(i+1) - rfw(i))**2 + (zfw(i+1) - zfw(i))**2 )
        ds = dist/(lmax - nlimw - 1)
      endif
    enddo
    rlimw = 100 * rlimw
    zlimw = 100 * zlimw
  endif
  stdout << "*** Start-up limiter, using rfw/zfw(" << idx_lim(1) << ":" \
         << idx_lim(2) << "), installed ***"

endf # set_limiter

  function vs(;n_harmonics,n_nodes,n_surf,f_thetac) # Calculate growth rate
  
    default(n_harmonics) =  10
    default(n_nodes)     =  10
    default(n_surf)      = 101
    default(f_thetac)    =   0.01
  
    # Preserve user's settings
    character note_save = note
    integer ic_save = ic
    integer ihy_save = ihy
  
    # Preserve coil connections
    integer i, j
    ic = 0
    ic(1) = 1
    do i = 2,nc
      do j = 1,i-1
        if (ic_save(i) = ic_save(j)) then # These are in series
          ic(i) = ic(j)
          break
        endif
      enddo
      if(ic(i) == 0) ic(i) = max(ic) + 1
    enddo
  
    # Get growth rate & stability factor
    fe_ed(2)
    ic = ic_save
    disp(n_harmonics,n_nodes)
    thetac = f_thetac
    msrf = n_surf
    ihy = 0
    note = "For vertical stability, msrf & thetac modified!"
    run
    vst
    real results(2) = [max(gamma), stabf]
  
    # Restore user's settings (but NOT msrf & thetac)
    note = note_save
    ihy = ihy_save
  
    return results
  
  endf # vs

function set_symmetric # Make plasma and PF coils up/down symmetric

  if (nsym == 1 & nsymc == 1) then
    return
  endif

  # Select only coils with Zc > 0...
  integer n, nc_new = 0
  do n = 1,nc
    if (zc(n) > 0) then
      nc_new        = nc_new + 1
      ic(nc_new)    = nc_new
      pfid(nc_new)  = pfid(n)
      rc(nc_new)    = rc(n)
      zc(nc_new)    = zc(n)
      drc(nc_new)   = drc(n)
      dzc(nc_new)   = dzc(n)
      cc(nc_new)    = cc(n)
      ntc(nc_new)   = ntc(n)
      cccap(nc_new) = cccap(n)
      bccap(nc_new) = bccap(n)
    endif
  enddo
  nc     = 2 * nc_new
  ncplot = nc

  # Select hard markers where Z >= 0...
  irl = 0; rl = 0
  integer nbd_new = 0
  do n = 1,nbd
    if (zbd(n) >= 0) then
      nbd_new      = nbd_new + 1
      rbd(nbd_new) = rbd(n)
      zbd(nbd_new) = zbd(n)
    endif
  enddo
  nbd = nbd_new
  if (nbd == 2) zbd = 0

  nsymc = 1 # Enforce coil up/down symmetry
  nsym  = 1 # Enforce plasma up/down symmetry
  if (debug == yes) then
    ihy = 0; run
  else
    ihy = 99; run
    chameleon symSaveName = substr(saveName,1,strlen(saveName)-4) // "-sym.sav"
    saveq(symSaveName)
    stdout << return << "Saved to disk in `"//symSaveName//"'"
  endif

endf # set_symmetric

function read_fwall # Read first-wall, divertor, limiter coordinates

  chameleon fname = "./fwall.inp" # Generic file name
  integer io = basopen(fname, "i")
  if (io <> -1) then
    io = basopen(fname, "r")
    noisy = yes
    chameleon fwall_id
    io >> fwall_id >> return
    noisy = no
    global integer nfw # No. points for entire FW
    io >> nfw >> return >> return
nfw
    real bufr(3,nfw)
    io >> bufr
    global real rfw(nfw) = transpose(bufr(1,))
    global real zfw(nfw) = transpose(bufr(2,))
    integer kfw(nfw) = transpose(int(bufr(3,)))
    call basclose(io)

    # Get nominal limiter point...
    if (min(kfw) == -1) then
      global real rlim_nom = 100 * rfw(mnx(kfw))
      global real zlim_nom = 100 * zfw(mnx(kfw))
      kfw = where(kfw < 0, 1, kfw)
    endif

    # Map first-wall to Corsica's rplate,zplate arrays...
    nsegs = 10
    nplates = nfw
    call gchange("Strk", 0)
    integer i
    do i = 1,nfw-1
#?      if (kfw(i+1) == kfw(i)) then
        rplate(i,1) = rfw(i); rplate(i,2) = rfw(i+1)
        zplate(i,1) = zfw(i); zplate(i,2) = zfw(i+1)
#?      endif
    enddo
    rplate(nfw,1) = rfw(nfw)
    rplate(nfw,2) = rplate(1,1)
    zplate(nfw,1) = zfw(nfw)
    zplate(nfw,2) = zplate(1,1)
    rplate = 100*rplate
    zplate = 100*zplate
  
    # Get indices of limiter (==1) and divertor (==2) regions...
    global integer idx_lim(2), idx_div(2)
    do i = 1,nfw
      if (kfw(i) == 1) then
        idx_lim(1) = i
        break
      endif
    enddo
    if (idx_lim(1) > 0) then
      do i = idx_lim(1),nfw
        if (kfw(i) <> 1) then
          idx_lim(2) = i - 1
          break
        endif
      enddo
    endif
    do i = 1,nfw
      if (kfw(i) == 2) then
        idx_div(1) = i
        break
      endif
    enddo
    if (idx_div(1) > 0) then
      do i = idx_div(1),nfw
        if (kfw(i) <> 2) then
          idx_div(2) = i - 1
          break
        endif
      enddo
    endif
    stdout << "Loaded FW/DIV/LIM model:        " << fwall_id
  endif

endf # read_fwall

function read_passive # Read passive structure definition

  chameleon fname = "./passive.inp"

  # Globals used by psm
  global integer n_ps, n_bp, n_vv

  real t
  integer i, io = basopen(fname, "i")
  if (io <> -1) then
    io = basopen(fname,"r")
    noisy = yes
    chameleon passive_id
    io >> passive_id >> return
    noisy = no
    integer n_ps_save
    io >> n_ps >> return >> return
    n_ps_save = n_ps
    real bufr(6,n_ps)
    io >> bufr
    real r1_ps  = transpose(bufr(1,))
    real z1_ps  = transpose(bufr(2,))
    real r2_ps  = transpose(bufr(3,))
    real z2_ps  = transpose(bufr(4,))
    real t_ps   = transpose(bufr(5,))
    real rho_ps = transpose(bufr(6,))
    call basclose(io)

    nwires = 0
    do i = 1,n_ps
      nwires = nwires + 1
      rwires(nwires) = 0.5 * (r1_ps(i) + r2_ps(i))
      zwires(nwires) = 0.5 * (z1_ps(i) + z2_ps(i))
      t = atan2(z2_ps(i) - z1_ps(i),r2_ps(i) - r1_ps(i))
      if (t < 0) t = t + 2 * pi
      if (t > 2*pi) t = t - 2 * pi
      if (t > pi/4 & t < 3*pi/4 | t > 5*pi/4 & t < 7*pi/4) then # Type 2
        drwires(nwires) = t_ps(i) / abs(sin(t))
        dzwires(nwires) = abs(z2_ps(i) - z1_ps(i))
        awires(nwires)  = 0
        awires2(nwires) = t
        if (awires2(nwires) > pi) awires2(nwires) = awires2(nwires) - 2 * pi
      else # Type 1
        drwires(nwires) = abs(r2_ps(i) - r1_ps(i))
        dzwires(nwires) = t_ps(i) / abs(cos(t))
        awires(nwires) = t
        if(awires(nwires) > pi) awires(nwires) = awires(nwires) - 2 * pi
        awires2(nwires) = 0
      endif
      rhwires(nwires) = rho_ps(i)
    enddo
    stdout << "Loaded passive structure model: " << passive_id
  endif

endf # read_passive

function read_pfcoil(;fname) # Read PF coil specifications

  default(fname) = "./pfcoil.inp"

  # Filament spacing parameter...
  real spacing = (0.01 * ro) / 25 # Filament spacing

  integer io = basopen(fname, "i")
  if (io == -1) then # Generate coils
    nc = ncoil
    ncoil = 0
    real t = 2*pi*float(iota(nc) - 0.5)/nc
    rc = rocoil + rbcoil*cos(t + dcoil*sin(t))
    zc = rbcoil*elcoil*sin(t)
    drc = 0.5*(zc(2) - zc(1))
    dzc = drc
    rc = 0.01*rc
    zc = 0.01*zc
    drc = 0.01*drc
    dzc = 0.01*dzc
    ic = iota(nc)
    integer i
    do i = 1,nc
      pfid(i) = "PF"//format(i,-2)
    enddo
  else
    io = basopen(fname, "r")
    noisy = yes
    chameleon pfc_id
    io >> pfc_id >> return
    noisy = no
    io >> nc >> return >> return
    ncplot = nc
    real buffer(7,nc)
    character*8 tempc, bufferc(nc)
    integer i
    do i = 1,nc
      noisy = yes
      io >> tempc
      noisy = no
      bufferc(i) = tempc
      io >> buffer(1:7,i)
    enddo
    call basclose(io)
    
    global real ntc(nc), cccap(nc), bccap(nc)
    do i = 1,nc
      pfid(i)  = bufferc(i)
      rc(i)    = buffer(1,i)
      zc(i)    = buffer(2,i)
      drc(i)   = buffer(3,i)
      dzc(i)   = buffer(4,i)
      ntc(i)   = buffer(5,i)
      cccap(i) = buffer(6,i)
      bccap(i) = buffer(7,i)
    enddo

  endif

  # Filament distribution...
  nrc = int(drc / spacing + 0.5)
  nzc = int(dzc / spacing + 0.5)
  nrc = where(nrc < 1, 1, nrc)
  nzc = where(nzc < 1, 1, nzc)

  # Default circuit arrangement (if currents aren't set)...
  if (sum(cc) == 0) then
    cc = 1
  endif

endf # read_pfcoil

function read_shape # Read plasma shape coordinates

  chameleon fname = "./shape.inp" # Generic file name
  integer io = basopen(fname, "i")
  if (io <> -1) then
    io = basopen(fname, "r")
    noisy = yes
    chameleon shape_id
    io >> shape_id >> return
    noisy = no
    io >> nfbd >> return >> return
    real bufr(2,nfbd)
    io >> bufr
    rfbd = 100 * transpose(bufr(1,))
    zfbd = 100 * transpose(bufr(2,))
    call basclose(io)
    stdout << "Loaded reference plasma shape:  " << shape_id
    return (0)
  else
    return (1)
  endif

endf # read_shape

function read_tfcoil # Read TF coil description

  chameleon fname = "./tfcoil.inp"
  integer io = basopen(fname, "i")
  if (io <> -1) then
    io = basopen(fname, "r")
    noisy = yes
    chameleon tfcoil_id
    io >> tfcoil_id >> return >> return
    noisy = no

    # Inner periphery...
    global integer ntfi2
    io >> ntfi2
    real bufr(2,ntfi2)
    io >> bufr
    global real rtfi2 = transpose(bufr(1,))
    global real ztfi2 = transpose(bufr(2,))

    # Outer periphery...
    # (Can't use eq.ntfo, eq.rtfo,... due to array limit)
    global integer ntfo
    io >> return >> ntfo
    real bufr(2,ntfo)
    io >> bufr
    global real rtfo = transpose(bufr(1,))
    global real ztfo = transpose(bufr(2,))
  
    # Current centerline (for OOP forces)...
    global integer ntfc      # No. of TF coils
    global integer ntfcpts
    io >> return >> ntfc >> ntfcpts
    if (ntfc*ntfcpts > 0) then
      real bufr(2,ntfcpts)
      io >> bufr
      global real rtfcpts = transpose(bufr(1,))
      global real ztfcpts = transpose(bufr(2,))
    endif
    call basclose(io)

    rtfo = 100 * rtfo
    ztfo = 100 * ztfo
    rtfi2 = 100 * rtfi2
    ztfi2 = 100 * ztfi2

    stdout << "Loaded TF coil specifications:  " << tfcoil_id

  endif

endf # read_tfcoil

  function dshape(n,r0,z0,aminor,elongation,triangularity) # Generate D shape

    # Generate analytic D-shape coordinates and put into
    # fuzzy marker arrays
    nfbd = n
    real t, dt = 2 * pi / nfbd
    integer i
    do i = 1,nfbd
      t = dt * (i - 1)
      rfbd(i) = r0 + aminor * cos(t + triangularity * sin(t))
      zfbd(i) = z0 + elongation * aminor * sin(t)
    enddo

  endf # dshape

  function odd(n) # Return n truncated down to nearest odd integer

    return (2 * (int(n + 1) / 2) - 1)

  endf # odd

  function set_limiter_pt # Set rlim(0),zlim(0)

    # Set limiter point...
    if (exists("rlim_nom")) then
      rlim(0) = rlim_nom
      zlim(0) = zlim_nom
    else # (Placeholder)
      rlim(0) = min(r(3),0.5 * (r(1) + rbin))
      zlim(0) = 0
    endif

  endf # set_limiter_pt
  
  function limw_from_plates

     real slimw(lmax),splate(nplates),dplate(nplates)
     integer i
     dplate(2:)=sqrt((rplate(2:)-rplate(:nplates-1))**2+ \
       (zplate(2:)-zplate(:nplates-1))**2)
     do i=2,nplates
       splate(i)=splate(i-1)+dplate(i)
     enddo
     slimw=(iota(lmax)-1)*splate(nplates)/(lmax-1)
     interp(splate,rplate,nplates,splbdry,splbdry,slimw,&rlimw,rlimw,lmax)
     interp(splate,zplate,nplates,splbdry,splbdry,slimw,&zlimw,rlimw,lmax)

  endf
     
