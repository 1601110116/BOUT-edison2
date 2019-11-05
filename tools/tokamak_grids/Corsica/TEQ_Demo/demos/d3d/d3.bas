scriptID("$Id: d3.bas,v 1.67 2012/03/05 19:19:11 bulmer Exp $")

corsica   # put corsica on top

# The followg parameters control the disposition of output...
integer d3_plot = [ 1, 1, 1, 1, 1, 1, 1, 1] # == 0 to turn off plots
integer d3_save = [ 1, 1] # == 0 to disable save-file writes
logical d3_debug = false # If true, call debugger before each direct-eq run
ezccgmc = 10000 # Allow big ncgm files
integer which_way=0 # 2: move boundary by thetac using direct solve or
                    # 1: move boundary by thetac directly and
                    #    then first inverse equilibrium using FF' or.
                    # 0: move boundary by thetac directly and
                    #    then first inverse equilibrium using jparsave.
                    #-1: move boundary by thetac directly and
                    #    then first inverse equilibrium using qsave.

                    #    Final inverse equilibrium uses inv_k and inv_p

logical use_eqdsk_psi = true
logical use_zero_der =false
integer imth=11

# -1: toroidal flux coordinate and qsave

# Inverse solver convergence criteria...
integer nht_default = 400
real epsrk_default = 1.0e-06

# Defaults for d3() arguments...
integer confirm_default = -1 
integer make_inverse_default = 1
integer msrf_default = 81
real thetac_default = 0.01
real fix_edge_default = 0.9
real hsrf_default=0.5
integer ixpt_default=1
integer inv_k_default = 0
integer inv_p_default = -1


# Internal parameter(s)...
real dsep_flag = 39 # EFIT's value for DSEP "undefined"
real cctol_default = 25 # kAt
real ccefit_low = 0.1 # max(cc(F-coils))/Ip
real alfbd_default(2) = [1e5,1] # Defaults for fuzzy marker weights
real alphac_default = 0.7 # Relaxation parameter
real setbox_default = 5
real dr_factor = 2
real dz_factor = 2
integer upper_divertor_shot = 100771 # Shot number where upper divertor installed

MDEF Comment()=
  << trim($1)
MEND
MDEF StatusReport()=
  if (debug <> 0 | d3_debug) then
    << "<"//trim($1)//">"
  endif
MEND

MDEF d3() =

  character*128 gname
  integer confirm, make_inverse, new_msrf, exceptions = 0
  integer new_inv_k, new_inv_p
  real fix_edge, new_thetac, new_hsrf

  IFELSE($1,) (gname        = "help",               gname        = $1)
  IFELSE($2,) (confirm      = confirm_default,      confirm      = $2)
  IFELSE($3,) (fix_edge     = fix_edge_default,     fix_edge     = $3)
  IFELSE($4,) (make_inverse = make_inverse_default, make_inverse = $4)
  IFELSE($5,) (new_msrf     = msrf_default,         new_msrf     = $5)
  IFELSE($6,) (new_thetac   = thetac_default,       new_thetac   = $6)
  IFELSE($7,) (new_inv_k    = inv_k_default,        new_inv_k    = $7)
  IFELSE($8,) (new_inv_p    = inv_p_default,        new_inv_p    = $8)
  IFELSE($9,) (new_hsrf     = hsrf_default,         new_hsrf     = $9)

  if (trim(gname) == "help") then
    << " "
    << "Generate a Corsica equilibrium (and inverse equilibrium) to match a DIII-D EFIT"
    << "equilibrium as defined in EQDSK g- and a-file pairs."
    << " "
    << "   d3(eqdsk_g-file;confirm,fix_edge,make_inverse,msrf,thetac,inv_k,inv_p)"
    << " "
    << "where..."
    << "   eqdsk_g-file == EQDSK file name (g-file);"
    << "   confirm       > 0 pause time (in seconds) after frame display"
    << "                == 0 don't show plot frames, or"
    << "                <  0 to pause until user responds, [default: " \
                           << confirm_default << "];"
    << "   fix_edge     == 0 to use the EQDSK profiles as-is, or"
    << "                >  0 to modify the edge where psibar > fix_edge, [" \
                           << fix_edge << "];"
    << "   make_inverse == 1 to create an inverse equilibrium, or"
    << "                == 0 to make only a direct-solve equilibrium, [" \
                           << make_inverse << "];"
    << "   msrf         == number of flux surfaces [" \
                           << msrf_default << "];"
    << "   thetac       == separatrix proximity thetac [" \
                           << thetac_default << "]."
    << "   inv_k        == profile choice for inverse solve, inv_k [" \
                           << inv_k_default << "];"
    << "   inv_p        == profile choice for inverse solve, inv_p [" \
                           << inv_p_default << "];"
    << " "
    << "The a-file name is constructed by changing ""g"" to ""a"" in the first character"
    << "position of argument eqdsk_g-file."
    << " "
    << "Logical global variable ""d3_debug"" can be used to ""call debugger"" before each"
    << "direct-solve run."
    << " "
    exceptions = exceptions + 1
  elseif (basopen(gname,"i") == -1) then
    warning("Cannot open EQDSK file: "//trim(gname))
    exceptions = exceptions + 1
  else
    integer magic_number = [17, 33, 65, 129, 257] # For grid size ferret
    global integer i, inumber(32)
    global integer nw_eqdsk = 0, nh_eqdsk = 0

    # This kluge ferrets out the number of mesh points from any other
    # "numbers" that may exist in the 1st line of the file
    integer io = basopen(gname,"r")
    noisy = no
    do i = 1,length(inumber)
      io >> inumber(i)
      if (i > 1) then
        if (min(abs(inumber(i-1) - magic_number)) == 0) then
          if (min(abs(inumber(i) - magic_number)) == 0) then
            nw_eqdsk = inumber(i-1)
            nh_eqdsk = inumber(i)
            break
          endif
        endif
      endif
    enddo
    call basclose(io)
    forget i, io, inumber, magic_number
    # Verify grid dimensions were found...
    if (nw_eqdsk * nh_eqdsk == 0) then
      warning("Cannot interpret NW & NH in 1st line of file: "//trim(gname))
      exceptions = exceptions + 1
    elseif (nprob == 0) # All ok, load generic equilibrium for this grid size...
      character*128 restore_file = trim("d3d"//format(nw_eqdsk,0)//"x" \
                                    //format(nh_eqdsk,0)//".sav")
      Comment("Restoring generic equilibrium: "//trim(restore_file))
      if(basopen(restore_file,"i")<>0) restore_file="d3d.sav"
      forget nw_eqdsk, nh_eqdsk
    endif
  endif
  logical restore_now = true
  if (exceptions > 0 | nprob > 0) then
    restore_now = false
  endif
  if (restore_now) then
    read {restore}.bas # And restore if restore_now == true
  endif
  qmax=1e99  # not to limit q_axis
  inv_k=new_inv_k # reset inv_k to its value after restore
  inv_p=new_inv_p # reset inv_p to its value after restore
  hsrf =new_hsrf  # reset inv_p to its value after restore
  if (exceptions == 0) then
    call d3f(gname,confirm,fix_edge,make_inverse,new_msrf,new_thetac)
  endif
  forget gname, confirm, fix_edge, make_inverse, new_msrf, new_thetac, exceptions

MEND # d3

function d3f(gname,confirm,fix_edge,make_inverse,new_msrf,new_thetac)

  # Change number of flux surfaces.
  if (new_msrf <> msrf) then
    note = "Changing msrf to " // format(new_msrf,0)
    msrf = new_msrf
    liml=-1
    run
  endif

  call pause(int(confirm))
  
  if (filter_wp==-2) filter_wp= -1 # Must use spline coefs and not polys
  if (filter_wp== 2) filter_wp=  1 # Must use spline coefs and not polys
  if (filter_p==-2)  filter_p=  -1 # Must use spline coefs and not polys
  if (filter_p== 2)  filter_p=   1 # Must use spline coefs and not polys

  integer regrid_flag = 0
  character*32 shotid = "shot.time"
  character*32 fileid = "shot_time"
  character*32 ccname = "ccshot.time"
  chameleon gfile = gname
  integer i = index(gfile,"/")
  integer n = strlen(trim(gfile))
  while (i <> 0)
    gfile = substr(trim(gfile),i+1,n-i)
    i = index(gfile,"/")
    n = strlen(trim(gfile))
  endwhile
  chameleon aname = "a"//substr(gfile,2,n-1)
  character*32 save_file = "teq.sav"
  character*32 inv_save_file = "teq_inv.sav"

  # Define the circuit indices for the F and E-coils...
  global integer ic_free = [
     1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,20,
    19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,
    19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,
    20,20,19,21,21,21,21,24,24,24,24,24,24,24,21,21,20,19,19,20,20,
    19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,
    19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,20,20,19,19,
    20,22,22,22,22,23,23,23,23,23,23,23,22,22]

  character*8 shot_date
  character*6 shot_number
  character*5 time_slice
  real shottime 
  integer magic_number = [17, 33, 65, 129, 257] # For grid size ferret
  integer i, inumber(32)
  integer g_nw = 0, g_nh = 0

  if (basopen(gname,"i") == -1) then
    warning("Cannot open EQDSK file: "//trim(gname))
  else

    integer io = basopen(gname,"r")
    # This kluge ferrets out the number of mesh points from any other
    # "numbers" that may exist in the 1st line of the file
    noisy = no
    do i = 1,length(inumber)
      io >> inumber(i)
      if (i > 1) then
        if (min(abs(inumber(i-1) - magic_number)) == 0) then
          if (min(abs(inumber(i) - magic_number)) == 0) then
            g_nw = inumber(i-1)
            g_nh = inumber(i)
            break
          endif
        endif
      endif
    enddo
    shottime = float(inumber(5)) * 1.0e-03 # Shot time, seconds

    # Verify grid dimensions were found...
    if (g_nw * g_nh == 0) then
      warning("Cannot interpret NW & NH in 1st line of file: "//trim(gname))
      call basclose(io)
      return
    endif

    # Get shot number and time point...
    shot_date = format(inumber(1),-2) // "/" \
             // format(inumber(2),-2) // "/" \
             // format(inumber(3),-2)
    shot_number = format(inumber(4),-6)
    time_slice = format(inumber(5),5)
    shotName = shot_number
    shotTime = inumber(5) * 1.0e-03

    # Load passive structure model...
    restore d3d_wires.pfb
    parsestr("$n = "//shotName)
    if ($n < upper_divertor_shot) then # Omit upper divertor supports
      nwires = 25
    endif

    # Set shot id string and file names...
    shotid = format(inumber(4),-6) // "." // format(inumber(5),-5)
    ccname = "cc" // shotid
    fileid = format(inumber(4),-6) // "_" // format(inumber(5),-4)
    save_file = trim(fileid) // ".sav"
    inv_save_file = trim(fileid) // "_inv.sav"

    # Get directory name...
    character*256 directory = dirname(trim(gname))

    # Read data from A EQDSK file (if it exists)...
    if (directory <> " ") then
      aname = trim(directory) // aname
    endif
    call read_aeqdsk(aname)

    # Set problem identification...
    call setd3id(shotName,shotTime)

    # Re-grid, if necessary...
    if (jm <> g_nw) regrid_flag = regrid_flag + 1
    if (km <> g_nh) regrid_flag = regrid_flag + 1
    if (regrid_flag > 0 ) then
      if(liml<>-1) then
        integer ldiff = int((log(jm - 1) - log(liml - 1) + 0.001) / log(2))
        liml = 2**(int((log(jm - 1) + 0.001) / log(2)) - ldiff) + 1
        limlo = 0
      endif
      jm = g_nw
      km = g_nh
      note="Re-gridding with jm = " // format(jm,0) \
           // ", km = " // format(km,0) // ", liml = " // format(liml,0)
      run
    endif

    mdef SkipData()=
      real junk($1)
      io >> junk
      forget junk
    mend

    # Close the file and re-open to get the correct position for either
    # EFIT or Corsica-generated EQDSKs...
    call basclose(io)
    io = basopen(gname,"r")
    io >> return # Now, skip the header line
    # Read the data (names per EFIT convention with "g_")...
    real g_xdim, g_zdim, g_rzero, g_rgrid1, g_zmid
    io >> g_xdim >> g_zdim >> g_rzero >> g_rgrid1 >> g_zmid
    SkipData(4)
    real g_bcentr, g_cpasma
    real g_ssimag, g_ssibry, g_zmaxis
    io >> g_bcentr
    io >> g_cpasma
    global real sign_conv = sign(1.0,g_cpasma)
    if (g_cpasma / placur < 0) then
      note = "Change sign of current, flux, etc..."
      cursign
    endif
    io >> g_ssimag
    SkipData(3)
    io >> g_zmaxis
    SkipData(1)
    io >> g_ssibry
    SkipData(2)
    real g_fpol(g_nw), g_pres(g_nw)
    io >> g_fpol
    io >> g_pres
    global real g_ffprim(g_nw), g_pprime(g_nw)
    io >> g_ffprim
    io >> g_pprime
    real g_psirz(g_nw,g_nh)
    io >> g_psirz
    global real g_qpsi(g_nw)
    io >> g_qpsi
    integer g_nbbbs, g_limitr
    io >> g_nbbbs >> g_limitr
    real g_rzbbbs(2*g_nbbbs)
    io >> g_rzbbbs # Boundary coordinates
    real g_xylim(2*g_limitr)
    io >> g_xylim
    call basclose(io)

    # Convert units...
    btor = g_bcentr * 1.0e+04 # Gauss
    plcm = g_cpasma * 1.0e-06 # MA
    real delta_psi = (g_ssimag - g_ssibry) * 1.0e+08
    g_fpol = g_fpol * 1.0e+06
    g_pres = g_pres * 10

    # Construct EFIT grid dimensions...
    real g_r(g_nw) = g_rgrid1 + g_xdim*float(iota(g_nw) - 1)/(g_nw - 1)
    real g_z(g_nh) = g_zmid + g_zdim*(float(iota(g_nh) - 1)/(g_nh - 1) - 0.5)

    # Map boundary to fuzzy marker arrays...
    integer n, nskip
    nskip = g_nbbbs/200 + 1
    nfbd = g_nbbbs/nskip
    alfbd = alfbd_default(1)
    n = 0
    do i = 1, g_nbbbs, nskip
      n = n + 1
      rfbd(n) = 100 * g_rzbbbs(2*i-1)
      zfbd(n) = 100 * g_rzbbbs(2*i)
    enddo

    # Map limiter coordinates to Corsica r,zplate...
    real g_xlim = fromone(g_xylim(1::2))
    real g_ylim = fromone(g_xylim(2::2))
    real dsq
    nplates = length(g_xlim) - 1
    call teq_sync
    nplates = 0
    do i = 1, length(g_xlim) - 1
      dsq = (g_xlim(i+1) - g_xlim(i))**2 + (g_ylim(i+1) - g_ylim(i))**2
      if (dsq > 0) then
        nplates = nplates + 1
        rplate(nplates,1) = g_xlim(i)
        zplate(nplates,1) = g_ylim(i)
        rplate(nplates,2) = g_xlim(i+1)
        zplate(nplates,2) = g_ylim(i+1)
      endif
    enddo
    rplate = 100*rplate
    zplate = 100*zplate
    call teq_sync
    call limw_plates

    if (use_eqdsk_psi) then
      StatusReport("Mapping EQDSK psi(R,Z) to Corsica R,Z grid")
      cps(,1) = -sign_conv*shape(g_psirz,jmkm)*1.0e+08
      real r_=r,z_=z
      r=100*g_r;z=100*g_z
      #   get spline coefficents on EQDSK grid
      integer imth=11;call bispline(r,z,g_nw,g_nh,&cps,imth,n4)
      #   get psi on Corsica grid
      call boxm(&rm,&zm,&psi,jmkm,0,1,n4)
      #   set flux
      psi00  = psi
      cps(,1)= psi
      r=r_;z=z_
      #   get spline coefficients on Corisca grid
      integer imth=11;call bispline(r,z,jm,km,&cps,imth,n4)
      #   evaluate topology
      wvert = 1
      if (exists("limloc")) then
        if (limloc == "DN") rl = 1
      endif
      call separatrix
      rl = 0
      if (limiterd <> 0) then
        real dsep_eqdsk = sign(dsep,zxpt)
      endif
    endif

    # Load Corsica's v, w arrays & get spline coefficients...
    real v_eqdsk(g_nw), w_eqdsk(g_nw)
    global real x_eqdsk(g_nw)
    x_eqdsk = float(iota(g_nw) - 1) / float(g_nw - 1)
    betaj = (0.01 * ro)**2 / (1 + (0.01 * ro)**2)
    v_eqdsk = g_pres / betaj / (sign_conv * delta_psi)
    w_eqdsk = (g_fpol**2 - g_fpol(g_nw)**2) / (8 * pi * (sign_conv * delta_psi))
    call interp(x_eqdsk,v_eqdsk,g_nw,splbdry,splbdry,psibar,&v,&vp,msrf)
    call interp(x_eqdsk,w_eqdsk,g_nw,splbdry,splbdry,psibar,&w,&wp,msrf)
    call spline(psibar,v,msrf,splbdry,splbdry,&vp,&vpp,&vppp,-999)
    call spline(psibar,w,msrf,splbdry,splbdry,&wp,&wpp,&wppp,-999)
    ipscl = 0 # Use EFIT plasma current

    # Don's edge modification...
    if (fix_edge > 0 & fix_edge < 1) then
      ipscl = -1 # To let plasma current float since profiles modified
      entropy_flag=0
      integer mfix = mnx(abs(psibar - fix_edge))
      if(~use_zero_der) then
        filter_p=0
        cos0 = 1
        cos0(mfix:) = (1-(psibar(mfix:)-psibar(mfix))**10/(1-psibar(mfix))**10)
        vp = vp*cos0
        call spline(psibar,vp,msrf,splbdry,splbdry,v,v,v,&v)
        v = v-v(msrf)
        cos0 = 1
        cos0(mfix:) = (1-(psibar(mfix:)-psibar(mfix))**10/(1-psibar(mfix))**10)
        wp = wp*cos0
        call spline(psibar,wp,msrf,splbdry,splbdry,w,w,w,&w)
        w = w-w(msrf)
        call spline(psibar,v,msrf,splbdry,0,&vp,&vpp,&vppp,-999)
        call spline(psibar,w,msrf,splbdry,0,&wp,&wpp,&wppp,-999)
      else
        mfix=msrf-mfix
        v=zero_der(v,mfix)
        call spline(psibar,v,msrf,splbdry,0,&vp,&vpp,&vppp,-999)
        w=zero_der(w,mfix)
        call spline(psibar,w,msrf,splbdry,0,&wp,&wpp,&wppp,-999)
      endif
    endif
    ipj = 0; ipp = 999; ipf = 999
    prsx = 1.0001
    filter_wp=0;filter_p=0
    package eq
    crdn = 2200 # Appropriate for DIII-D
    nbd = 0
    ircwt = 2 # To minimize sum(NIcoil**2)
    ic = ic_free
    if (epsb <> 1) then
      alphac = alphac_default; nl = 301; ni = 1; epsj = 1.0e-06; epsb = 1
    endif
    kaxis = 1; zaxis0 = 100 * g_zmaxis
    kaxis = 0 # For better convergence (in general)

    # Nominal (inboard) limiter location...
    rlim(0) = 101.6
    zlim(0) = zfbd(mnx(rfbd))

    # Extract file name from path name...
    i = index(gname,"/")
    chameleon fname = gname
    while (i > 0)
      fname = substr(gname,i+1,strlen(trim(gname))-i)
      gname = fname
      i = index(gname,"/")
    endwhile
    if (exists("limloc")) then
      note = "Constructed from a," // trim(fname)
    else
      note = "Constructed from " // trim(fname)
    endif
    wvert=1
    # Refine constraints...
    rxpr = [min(rplate)+dr_factor*dr,160]
    zxpr = [max(min(zfbd)-20,min(zplate)),min(max(zfbd)+20,max(zplate))]
    irl = 0; rl = 0; ixpt = 1; alphac = alphac_default
    if (exists("limloc")) then # Use info from EQDSK-A
      StatusReport("EFIT LIMLOC: "//trim(limloc))
      if (exists("rxpt_efit")) then
        StatusReport("Adding 2 EFIT R,Zxpt constraints")
        integer nxpt = mnx((rxpt_efit - rfbd)**2 + (zxpt_efit - zfbd)**2)
        nbd = 2
        rbd = [rfbd(nxpt-2),rfbd(nxpt+2)]
        zbd = [zfbd(nxpt-2),zfbd(nxpt+2)]
        StatusReport("Restricting vertical extent of x-point search box")
        if (zxpt_efit < 0) then
          zxpr(1) = zxpt_efit - 10
          zxpr(2) = max(zplate)
        else
          zxpr(1) = min(zplate)
          zxpr(2) = zxpt_efit + 10
        endif
        if (limloc == "DN") then
          zxpr = [min(zfbd)-10,max(zfbd)+10]
        endif
      else # Limited
        StatusReport("Fictitious limiter point near inboard midplane")
        ixpt = 0
        alphac = 0.5 * (1 + alphac_default)
        real dist(nfbd)
        dist = rfbd**2 + zfbd**2
        rlim(0) = rfbd(mnx(dist))
        zlim(0) = zfbd(mnx(dist))
        if (rlim(0) < 0.5*(raxis+max(rfbd))) then
          StatusReport("Adding fixed marker at outboard edge")
          nbd = 1
          rbd(1) = rfbd(mxx(rfbd))
          zbd(1) = zfbd(mxx(rfbd))
        else
          StatusReport("Adding fixed marker at inboard edge")
          nbd = 1
          rbd(1) = rfbd(mnx(rfbd))
          zbd(1) = zfbd(mnx(rfbd))
        endif
      endif
      StatusReport("Installing EFIT E-coil currents")
      if (max(abs(ccefit)) == 0) then # Dont' use any
        StatusReport("EFIT coil currents not available")
        ircwt = 2
        alfbd = alfbd_default(1)
      else
        if (max(abs(ccefit(1:nfcc)))/abs(plcm) < ccefit_low) then
          StatusReport("Converting EFIT coil currents from A to A-turns")
          real turns_fc(nfcc) = [58,58,58,58,58,55,55,58,55]
          turns_fc(nfcc/2+1:nfcc) = turns_fc(1:nfcc/2)
          ccefit(1:nfcc) = turns_fc(1:nfcc)*ccefit(1:nfcc)
        endif
        StatusReport("Trying to match EFIT F-coil currents")
        cc = ccefit; cc0 = ccefit
        ic(nfcc+1:) = 0 # Fix E-coil currents
        ic = iota(nfcc)
        ircwt = 1
        alfbd = alfbd_default(2)
      endif
    else # Add a hard marker at outboard edge if none exist...
      if (nbd == 0) then
        StatusReport("Adding fixed marker at outboard edge")
        nbd = 1
        rbd(1) = rfbd(mxx(rfbd))
        zbd(1) = zfbd(mxx(rfbd)) - 0.1 # to avoid limiter point
      endif
    endif
    if (exists("limloc")) then
      if (limloc == "DN") alfbd = alfbd_default(1)
    endif
    if (d3_debug) call debugger
    run

    # Refine separatrix separation/weights...
    if (exists("dsep_efit")) then
      if (limloc == "DN") then
        if(exists("dsep_eqdsk")) then
          dsep_efit = dsep_eqdsk
        else
          dsep_efit = sign(1.0e-03,zxpt)
        endif
      endif
      if (dsep <> 0 & abs(dsep_efit) < dsep_flag) then
        rxpr = [min([rxpt_efit-dr,rxpt2_efit-dr,rxpr(1),rxpr(2)]), \
	       max([rxpt_efit+dr,rxpt2_efit+dr,rxpr(1),rxpr(2)])]
        zxpr = [min([zxpt_efit-dz,zxpt2_efit-dz,zxpr(1),zxpr(2)]), \
	       max([zxpt_efit+dz,zxpt2_efit+dz,zxpr(1),zxpr(2)])]
        if (limloc == "DN") then
          StatusReport("Adding DN constraint")
          rl = abs(dsep_efit)
          irl = int(sign(1.0,-dsep_efit))
          alfbd = alfbd_default(2)
          if (d3_debug) call debugger
          run
        else # SN with dsep constraint
          StatusReport("Setting nbd = 0, adding EFIT DSEP constraint")
          rl = abs(dsep_efit)
          irl = int(sign(1.0,-dsep_efit))
        endif
        nbd = 0
        alfbd = alfbd_default(2)
        if (d3_debug) call debugger
        run
        if (residj > epsj) then
          StatusReport("Refining x-point search box with setbox()")
          call setbox
          if (d3_debug) call debugger
          run
        endif
      endif
    endif

    if (limiterd <> 0) then
      StatusReport("Moving limiter point to realistic location")
      real dist(nfbd)
      do i = 1,nfbd
        dist(i) = min((rfbd(i) - rplate(,1))**2 + (zfbd(i) - zplate(,1))**2)
      enddo
      rlim(0) = rfbd(mnx(dist))
      zlim(0) = zfbd(mnx(dist))
      if (d3_debug) call debugger
      run
    endif

    if (exists("limloc")) then # Check coil current variance
      if (1000 * max(abs(cc - ccefit)) > cctol_default) then
        warning(trim(fileid)//" has significant F-coil current" \
                            //" variance from EFIT")
      endif
    endif
    if (residj > epsj) then
      if (residj < 50 * epsj) then
        warning(trim(fileid)//" not fully converged")
      else
        warning(trim(fileid)//" has poor convergence")
      endif
    endif

    if (ipscl <> 0) then # Plasma currents will differ...
      remark "Plasma current comparison (ipscl==-1)..."
      remark "EFIT:    "//format(int(g_cpasma),8)//" A"
      remark "Corsica: "//format(int(10*plc),8)//" A  (" \
       //format(100*(10*plc-g_cpasma)/g_cpasma,0,2,1)//" %)"
    endif
    if (sum(d3_plot(1:8)) <> 0 & confirm <>0) then
      Comment("Check free-boundary equilibrium...")
      win on d3-comparison
    endif
    if (d3_plot(1) <> 0) then
      pb; pause("Plasma boundary")
    endif
    if (d3_plot(2) <> 0) then
      ppprime; pause("Pressure")
    endif
    if (d3_plot(3) <> 0) then
      pffprime; pause("Toroidal flux function")
    endif
    if (d3_plot(4) <> 0) then
      pq; pause("q-profile")
    endif
    if (d3_plot(5) <> 0) then
      pd3cc; pause("Coil currents")
    endif
    if (d3_plot(6) <> 0) then
      layout; pause("Configuration")
    endif
    # Save the direct-solve equilibrium...
    if (d3_save(1) <> 0) then
      saveq(save_file)
      Comment("Created equilibrium save-file: "//trim(save_file))
    endif

    # Use these to restore direct equilibrium solution...
    real v_save = [v, vp, vpp, vppp]
    real w_save = [w, wp, wpp, wppp]
    if (make_inverse <> 0) then # Make an inverse equilibrium save-file...
      thetac = new_thetac
      if(which_way==2) then # direct equilibrium with thetac
        Comment(char(10)//"Generating direct equilibrium with thetac=" \
        //format(thetac,0,2,0))
        run
      endif

      if (limiterd <> 0 & (pxpt-psil)/dpsi0>thetac) thetac=0
      if (confirm <> 0) then
        Comment(char(10)//"Construct inverse equilibrium with thetac=" \
                //format(thetac,0,2,0)//"...")
      endif

      nht = nht_default
      epsrk = epsrk_default
      morph_inv(which_way)

      if (d3_plot(3) <> 0) then
        pffprime; pause("FF'-profile")
      endif
      if (d3_plot(4) <> 0) then
        pq; pause("q-profile")
      endif
      if(d3_plot(7) <> 0) then
        contour; pause("Configuration")
      endif
      if(d3_plot(8) <> 0) then
        testGS; pause("GS error estimate")  
      endif

      if (d3_save(2) <> 0) then
        saveq(inv_save_file)
        Comment("Created inverse-equilibrium save-file: "//trim(inv_save_file))
      endif

      Comment("Recovering direct equilibrium...")
      thetac = 0
      inv_eq = 0
      v = v_save(,1); vp = v_save(,2); vpp = v_save(,3); vppp = v_save(,4)
      w = w_save(,1); wp = w_save(,2); wpp = w_save(,3); wppp = w_save(,4)
      recoup
    endif
  endif

  if (sum(d3_plot(1:7)) <> 0 & confirm <>0) then
    win close d3-comparison
  endif

endf # d3f

function d3plots(;p1,p2,p3,p4,p5,p6,p7)

  default(p1) = d3_plot(1)
  default(p2) = d3_plot(2)
  default(p3) = d3_plot(3)
  default(p4) = d3_plot(4)
  default(p5) = d3_plot(5)
  default(p6) = d3_plot(6)
  default(p7) = d3_plot(7)

  d3_plot = [p1,p2,p3,p4,p5,p6,p7]
  d3_plot = where(d3_plot <> 0, 1, d3_plot)

  remark "PLOT SELECTION     ON/OFF"
  remark " (1) Equilibrium   "//format(d3_plot(1),3)
  remark " (2) p' profile    "//format(d3_plot(2),3)
  remark " (3) FF' profile   "//format(d3_plot(3),3)
  remark " (4) q profile     "//format(d3_plot(4),3)
  remark " (5) Coil currents "//format(d3_plot(5),3)
  remark " (6) Configuration "//format(d3_plot(6),3)
  remark " (7) Inverse eq.   "//format(d3_plot(7),3)

endf # d3plots

function setbox(;xptmargin)

  default(xptmargin) = setbox_default

  if (xptmargin == 0) then
    rxpr = [min(rplate)+dr_factor*dr, 160]
    zxpr = [min(zplate), max(zplate)]
    return
  endif

  call reorder_xps

  # Dominant and secondary x-points...
  real px_pair(2) = pxps(1)
  real rx_pair(2) = rxps(1)
  real zx_pair(2) = zxps(1)

  # Nearest x-point at opposite Z...
  real dxtoa, dxtoa_min = r1mach(2)
  integer n
  do n = 2,nxps
    if (zxps(n) / zx_pair(1) < 0) then
      dxtoa = (rxps(n) - raxis)**2 + (zxps(n) - zaxis)**2
      if (dxtoa < dxtoa_min) then
        dxtoa_min = dxtoa
        px_pair(2) = pxps(n)
        rx_pair(2) = rxps(n)
        zx_pair(2) = zxps(n)
      endif
    endif
  enddo

  rxpr = [min(rx_pair) - xptmargin, max(rx_pair) + xptmargin]
  zxpr = [min(zx_pair) - xptmargin, max(zx_pair) + xptmargin]
  if (setbox_default <> 0) then
    if (rxpr(1) < min(rplate)) rxpr(1) = min(rplate) + dr_factor*dr
    if (zxpr(1) < min(zplate)) zxpr(1) = min(zplate)
    if (zxpr(2) > max(zplate)) zxpr(2) = max(zplate)
  endif

endf # setbox

function reorder_xps

  if (nxps < 3) return

  real pxps_t = pxps
  real rxps_t = rxps
  real zxps_t = zxps

  real dxps = sqrt((rxps - raxis)**2 + (zxps - zaxis)**2)

  integer i, im
  do i = 1,length(dxps)
    im = mnx(dxps)
    pxps(i) = pxps_t(im)
    rxps(i) = rxps_t(im)
    zxps(i) = zxps_t(im)
    dxps(im) = r1mach(2)
  enddo

endf # reorder_xps

# Use this function to map the (2|6) EC circuit currents to the E-coils...
function map_ecc(eccefit,necc,nfcc)

  real c1, c2, c3, c4, c5, c6
  integer i
  integer nc_ec = nc - nfcc # No. of E-coils
  integer ic_ec(nc_ec)

  # Circuit configuration for 6 EC circuits...
  ic_ec = [
    1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1,
    1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1,
    1, 2, 2, 1, 1, 2, 2, 1, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 3,
    3, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1,
    2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1,
    2, 2, 1, 1, 2, 2, 1, 1, 2, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
    4, 4]
  # Mapping for 2 EC circuits...
  #   3,5 -> 1
  #   4,6 -> 2

  # Copy the EFIT currents...
  c1 = eccefit(1)
  c2 = eccefit(2)
  if (necc == 2) then
    c3 = eccefit(1)
    c4 = eccefit(2)
    c5 = eccefit(1)
    c6 = eccefit(2)
  else
    c3 = eccefit(3)
    c4 = eccefit(4)
    c5 = eccefit(5)
    c6 = eccefit(6)
  endif
  do i=1,nc-nfcc
    if (ic_ec(i) == 1) ccefit(nfcc+i) = c1
    if (ic_ec(i) == 2) ccefit(nfcc+i) = c2
    if (ic_ec(i) == 3) ccefit(nfcc+i) = c3
    if (ic_ec(i) == 4) ccefit(nfcc+i) = c4
    if (ic_ec(i) == 5) ccefit(nfcc+i) = c5
    if (ic_ec(i) == 6) ccefit(nfcc+i) = c6
  enddo

endf # map_ecc

function pd3cc(;ordering) # Plot coil current comparison
 default(ordering) = 1 # == 0 for poloidal ordering of coil currents
                       # == 1 for DIII-D ordering of coil currents
 if (~exists("nfcc")) then
   Comment("Coil currents from EFIT are not available")
   return
 endif

 integer n,icd3(nfcc)
 real ccd3(nfcc),ccd3_efit(nfcc)

 teqplotdefs; nf
 ezctitfr = 0.06

 # F-coil currents...
 ezcquad(12)
 frame 0 nfcc+1
 ezcextra  =  0.05 # extra vertical space beyond data
 if (ordering  ==  0) then
   titles "TEQ and EFIT (X) currents" \
          "F-coil ordering"//colon//" 1,2,3,4,5,8,9,7,6...","NIcoil   [kAt]"
   icd3 = [ic(1:5),ic(8:9),ic(7),ic(6),ic(10:14),ic(17:nfcc),ic(16),ic(15)]
   ccd3 = 1000*[cc(1:5),cc(8:9),cc(7),cc(6),cc(10:14),cc(17:nfcc), \
                cc(16),cc(15)]
   ccd3_efit = 1000*[ccefit(1:5),ccefit(8:9),ccefit(7),ccefit(6), \
                     ccefit(10:14),ccefit(17:nfcc),ccefit(16),ccefit(15)]
 else
   titles "TEQ (curve) and EFIT (X) F-coil currents" \
          "DIII-D F-coil number","NIcoil   [kAt]"
   icd3 = ic(1:nfcc)
   ccd3 = 1000*cc(1:nfcc)
   ccd3_efit = 1000*ccefit(1:nfcc)
 endif
 plot ccd3 thick=2 color=orangered
 plot ccd3_efit mark cross marksize 2
 plot 0.0*ones(nfcc) style dotted
 do n = 1,nfcc
   if (icd3(n) == 0) plot ccd3(n),n mark=circle color=yellow
 enddo

 # E-coil currents...
 integer iced3(6)
 real cc_ecoil(6), cc_efit_ecoil(6)
 iced3(1) = ic( 19); cc_ecoil(1) = cc( 19); cc_efit_ecoil(1) = ccefit( 19)
 iced3(2) = ic( 20); cc_ecoil(2) = cc( 20); cc_efit_ecoil(2) = ccefit( 20)
 iced3(3) = ic( 67); cc_ecoil(3) = cc( 67); cc_efit_ecoil(3) = ccefit( 67)
 iced3(4) = ic(128); cc_ecoil(4) = cc(128); cc_efit_ecoil(4) = ccefit(128)
 iced3(5) = ic(132); cc_ecoil(5) = cc(132); cc_efit_ecoil(5) = ccefit(132)
 iced3(6) = ic( 71); cc_ecoil(6) = cc( 71); cc_efit_ecoil(6) = ccefit( 71)
 ezcquad(34)
 frame 0 7
 titles "TEQ (curve) and EFIT (X) E-coil circuit currents" \
        "DIII-D E-coil circuit number","NIcoil   [kAt]"
 plot 1000*cc_ecoil,iota(6) thick=2 color=orangered
 plot 1000*cc_efit_ecoil,iota(6) mark cross marksize 2
 plot 0.0*ones(6),iota(6) style dotted
 do n = 1,6
   if(iced3(n) == 0) plot 1000*cc_ecoil(n),n mark=circle color=yellow
 enddo

 nf; teqplotreset

endf # pd3cc

function ppprime # Plot p' in mks units

  StartFrame
    titles probid "normalized psi" "p'"
    plot -sign_conv*g_pprime,x_eqdsk mark circle marksize 0.5
    plot pdsrf*1.0e+07,psibar color red
    plot 0.0*ones(msrf),psibar color gray70
  EndFrame

endf # ppprime

function pffprime # Plot FF' in mks units

 teqplotdefs; nf
 titles probid "normalized psi" "FF'"
 plot -sign_conv*g_ffprim,x_eqdsk mark circle marksize 0.5
 plot frsrf*fpsrf*1.0e-04,psibar color red
 plot 0.0*ones(msrf),psibar color gray70
 nf; teqplotreset

endf # pffprime

function pq # Plot q-profile

 teqplotdefs; nf
 titles probid "normalized psi" "q"
 frame ,,0
 plot g_qpsi,x_eqdsk mark circle marksize 0.5
 plot qsrf,psibar color red
 nf; teqplotreset

endf # pq

function dirname(pathname)

  character*256 dir_name = " "
  character*256 str1 = pathname
  character*256 str2
  integer ix

  while(true)
    ix = index(str1,"/")
    if (ix > 0) then
      str2 = substr(str1,ix+1,strlen(str1))
      dir_name = trim(dir_name) // substr(str1,1,ix)
      str1 = str2
    else
      break
    endif
  endwhile
  return trim(dir_name)

endf # dirname

function read_aeqdsk(fname)

  ## Read an A-EQDSK file and extract the EFIT configuration code,
  ## returned in global character variable `limloc'.  If the
  ## configuration is diverted, set the global variables `rxpt_efit'
  ## and `zxpt_efit' to contain the coordinates of the dominate
  ## X-point, and `dsep_efit' to the separatrix separation.
  ##
  ## EFIT A-EQDSK Limiter-Location Codes (LIMLOC)...
  ##
  ## Limited:             "IN", "OUT", "TOP", "BOT"
  ## Single-null:         "SNT", "SNB"
  ## Double-null:         "DN"
  ## Marginally diverted: "MAR" (treat as limited)

  ## So we can existence-check after call..
  forget rxpt_efit
  forget zxpt_efit
  forget rxpt2_efit
  forget zxpt2_efit
  forget dsep_efit
  forget limloc

  ## If the A-EQDSK file does not exist, do nothing...
  if (basopen(fname,"i") == -1) then
    return
  endif

  ## Since the file exists, create these variables...
  global character*3 limloc
  global real rxpt_efit, zxpt_efit, dsep_efit
  global real rxpt2_efit, zxpt2_efit
  integer mco2v, mco2r
  integer nsilop0, magpri0, nfcoil0, nesum0

  ## Read the contents of the file...
  real rnum1, rnum2, rnum3, rnum4, rnum5, junk
  character*8 str1
  integer n
  mdef SkipLines()=
    do n = 1,$1
      io >> return
    enddo
  mend
  mdef SkipItems()=
    if ($1 > 0) then
      do n = 1,$1
        io >> junk
      enddo
    endif
  mend
  integer io = basopen(fname,"r")
  SkipLines(3)
  noisy = yes
  io >> str1 >> rnum1 >> rnum2 >> rnum3 >> limloc
  io >> mco2v >> mco2r >> return
  noisy = no
  SkipLines(6)
   SkipItems(mco2v)
  SkipItems(mco2v)
  SkipItems(mco2r)
  SkipItems(mco2r)
  SkipLines(5)
  io >> rnum1 >> rnum2 >> rnum3 >> rnum4
  SkipLines(5)
  io >> nsilop0 >> magpri0 >> nfcoil0 >> nesum0
  global real cc_fcoils(nfcoil0), cc_ecoils(nesum0)
  SkipItems(nsilop0)
  SkipItems(magpri0)
  if (nfcoil0 > 0) then
    io >> cc_fcoils
  endif
  if (nesum0 > 0) then
    io >> cc_ecoils
  endif
  SkipLines(3)
  SkipItems(2)
  io >> rnum5
  call basclose(io)

  # Put EFIT coil currents into ccefit...
  global integer nfcc, necc
  nfcc = nfcoil0
  necc = nesum0
  if (nfcoil0 > 0 & nc > 0) then
    global real ccefit(nc)
    ccefit(1:nfcc) = cc_fcoils
    call map_ecc(cc_ecoils,necc,nfcc)
    ccefit = ccefit*1.0e-06
    ccefit = where(ccefit == 0, 1.0e-10, ccefit)
  endif

  ## Set the dominate X-point coordinates and separation...
  ## Also set the secondary from EFIT for possible to
  ## set the search box

  if (limloc == "SNB") then
    rxpt_efit = rnum1
    zxpt_efit = rnum2
    rxpt2_efit = rnum3
    zxpt2_efit = rnum4
    dsep_efit = rnum5
  elseif (limloc == "SNT") then
    rxpt_efit = rnum3
    zxpt_efit = rnum4
    rxpt2_efit = rnum1
    zxpt2_efit = rnum2
    dsep_efit = rnum5
  elseif (limloc == "DN") then
    if (rnum5 < 0) then
      rxpt_efit = rnum1
      zxpt_efit = rnum2
      rxpt2_efit = rnum3
      zxpt2_efit = rnum4
    else
      rxpt_efit = rnum3
      zxpt_efit = rnum4
      rxpt2_efit = rnum1
      zxpt2_efit = rnum2
    endif
    dsep_efit = rnum5
  else # So we can existence-test after call
    forget rxpt_efit
    forget zxpt_efit
    forget rxpt2_efit
    forget zxpt2_efit
    forget dsep_efit
  endif

endf # read_aeqdsk

  function setd3id(name,time)

    probid = "DIII-D   shot #"//substr(name,1,6)//" @ " \
           // format(time,0,3,1) // " s"
    if (exists("limloc")) then
      probid = trim(probid)//" ["//trim(limloc)//"]"
    endif

  endf # setd3id

function makebd(;regexp) # Make time-dependent plasma boundary database

  default(regexp) = "*_inv.sav"

  if (substr(regexp,1,1) == "h") then
    remark " "
    remark "Function `makebd(;regexp)' restores multiple save-files which match the regular"
    remark "expression `regexp' [default: ""*_inv.sav""] and records the time-point, boundary"
    remark "coordinates and shape parameters.  The data is saved in file `<shot>_bd.pfb'"
    remark "which is used by the self-contained function `getbd(t)'."
    remark " "
    return
  endif

  chameleon files = filelist(regexp)
  if (files(1) == "None") then
    fatal("makebd: no files match """//regexp//""".")
  endif

  integer i, n = length(files)
  if (n < 2) then
    fatal("makebd: need 2 or more inverse save-files")
  endif

  # Load equilibria...
  do i = 1,n
    <<format(i,int(log10(n)+1))<<" "<<files(i)
    restore ^files(i)
    if (inv_eq <> 0) then
      teq_inv
    else
      corsica; run
    endif
    if (i == 1) then
      if (residj > epsj) then
        fatal("makebd: first equilibrium """//trim(files(i))//""" not converged")
      else
        global character*16 shot_bd = shotName
        global real time_bd = shotTime
        global real rls_bd(mls,1) = rls
        global real zls_bd(mls,1) = zls
        global real rcntr_bd = rcntr
        global real rbore_bd = rbore
        global real elongn_bd = elongn
        global real trian_bd = trian
        global real thetac_bd = thetac
        global integer inv_eq_bd = inv_eq
        global integer mls_bd = mls
      endif
    endif
    if (shotName <> shot_bd) then
      fatal("makebd: """//trim(files(i))//""" has an inconsistent shot number")
    endif
    if (inv_eq <> inv_eq_bd) then
      fatal("makebd: """//trim(files(i))//""" has an inconsistent inv_eq")
    endif
    if (thetac <> thetac_bd) then
      fatal("makebd: """//trim(files(i))//""" has an inconsistent thetac")
    endif
    if (mls <> mls_bd) then
      fatal("makebd: """//trim(files(i))//""" has an inconsistent mls")
    endif
    if (residj > epsj) then
      warning("makebd: """//files(i)//""" not converged---will skip")
    elseif (i > 1) then
      time_bd   := shotTime
      rcntr_bd  := rcntr
      rbore_bd  := rbore
      elongn_bd := elongn
      trian_bd  := trian
      rls_bd    := rls
      zls_bd    := zls
    endif
  enddo

  chameleon pfbname = trim(shot_bd)//"_bd.pfb"
  create ^pfbname
    write inv_eq_bd, thetac_bd, mls_bd
    write time_bd, rcntr_bd, rbore_bd, elongn_bd, trian_bd
    write rls_bd, zls_bd
    writef getbd
  close
  << "Created " << pfbname

endf # makebd


function getbd(;t,usage) # Get time-dependent plasma boundary from database

  default(t) = 0
  default(usage) = 0

  if (type(t) == "character") then
    remark " "
    remark "Function `getbd(;t,usage)' interpolates for the plasma boundary at time `t' [s]"
    remark "using data stored in file ""<shot>_bd.pfb"" (created by function `makebd'). If"
    remark "the requested time is zero [default], the data in the PFB file will be graphed."
    remark " "
    remark "Using the interpolated boundary results when t > 0..."
    remark "   if usage == 0, the boundary will be returned in variables `rlst' & `zlst';"
    remark "   if usage  > 0, the boundary will be installed in (uk,vk) or (rfbd,zfbd); or"
    remark "   if usage  < 0, the boundary will be graphed."
    remark "If usage > 0 and inv_eq == 0, the fuzzy marker weighting vector `alfbd' will be"
    remark "assigned with ""alfbd=usage""."
    remark " "
    remark "This function returns the following integer status codes:"
    remark " "
    remark "   0 : no error"
    remark "   1 : cannot find file ""<shot>_bd.pfb"""
    remark "   2 : inv_eq inconsistent with ""<shot>_bd.pfb"""
    remark "   3 : mls inconsistent with ""<shot>_bd.pfb"""
    remark "   4 : thetac inconsistent with ""<shot>_bd.pfb"""
    remark "   5 : requested time out-of-range"
    remark " "
    return
  endif

  # Load the plasma boundary database...
  chameleon pfbname = trim(shotName)//"_bd.pfb"
  if (basopen(pfbname,"i") <> 0) then
    warning("getbd: cannot open """//pfbname//""".")
    return 1
  else
    if (~exists("pfbname_last")) then
      global chameleon pfbname_last = " "
    endif
    if (pfbname_last <> pfbname) then
      restore ^pfbname
      pfbname_last = pfbname
    endif
  endif
  integer i, n = length(time_bd)
  integer mls_bd = shape(rls_bd)(1)
  chameleon trange = format(min(time_bd),0,3,1)
  trange = trange//" to "//format(max(time_bd),0,3,1)//" s"

  if (t <= 0) then # just plot the boundaries
    <<return<<"makebd database "//trim(pfbname) \
    //", for: inv_eq="//format(inv_eq_bd,0) \
    //", mls="//format(mls_bd,0) \
    //", thetac="//format(thetac_bd,0,3,2) << return
    chameleon titlex = format(n,0)//" boundaries, from "//trange
    StartFrame
      attr scale = equal
      titles ":F2:Boundary database "//pfbname ":F2:"//titlex
      plot zls_bd,rls_bd # color=rainbow
    EndFrame
    pause
    StartFrame
      ezctitfr=0.06
      call agseti("Y/NICE.",yes)
      ezcquad(12)
        titles ":F2:Boundary shape parameters from "//pfbname " " ":F2:R:B:0:E:  & a"
        plot rcntr_bd, time_bd thick=3 color=orangered
        plot rbore_bd, time_bd thick=3 color=springgreen
        plot rcntr_bd, time_bd thick=3 mark=plus
        plot rbore_bd, time_bd thick=3 mark=plus
      ezcquad(34)
        titles " " ":F2:time [s]" ":F8:k:F2: & :F8:d"
        plot elongn_bd, time_bd thick=3 color=orangered
        plot trian_bd, time_bd thick=3 color=springgreen
        plot elongn_bd, time_bd thick=3 mark=plus
        plot trian_bd, time_bd thick=3 mark=plus
    EndFrame
    call agseti("Y/NICE.",no)
    return
  endif

  # Consistency checks...
  if (inv_eq <> inv_eq_bd) then
    warning("getbd: inv_eq inconsistent with data in "//pfbname)
    return 2
  elseif (mls <> mls_bd) then
    warning("getbd: mls inconsistent with data in "//pfbname)
    return 3
  elseif (thetac <> thetac_bd) then
    warning("getbd: thetac inconsistent with data in "//pfbname)
    return 4
  elseif (t < min(time_bd) | t > max(time_bd)) then
    warning("getbd: requested time is out of range: "//trange)
    return 5
  endif

  # Interpolate for boundary at time t...
  global real rlst(mls_bd), zlst(mls_bd)
  i = mnx(abs(t - time_bd))
  if (min(abs(t - time_bd)) == 0) then # No interpolation required
    rlst = rls_bd(,i)
    zlst = zls_bd(,i)
  else # Linearly interpolate
    real xt = (t - time_bd(1:n-1))/(t - time_bd(2:n))
    i = mnx(xt)
    real ft = (t - time_bd(i))/(time_bd(i+1) - time_bd(i))
    rlst = rls_bd(,i) + ft*(rls_bd(,i+1) - rls_bd(,i))
    zlst = zls_bd(,i) + ft*(zls_bd(,i+1) - zls_bd(,i))
  endif

  # Plot it or install in inverse solver arrays...
  if (usage < 0) then
    chameleon titlex = "Boundary (rlst,zlst) at t="//format(float(t),0,3,1)//" s"
    StartFrame
      attr scale = equal
      titles ":F2:Boundary database "//pfbname ":F2:"//titlex
      plot zlst,rlst thick=3
    EndFrame
  elseif (usage > 0) then
    if (inv_eq == 0) then
      nfbd = mls
      rfbd = rlst
      zfbd = zlst
      alfbd = usage
    else
      uk = rlst
      vk = zlst
    endif
  endif

  return 0

endf # getbd
