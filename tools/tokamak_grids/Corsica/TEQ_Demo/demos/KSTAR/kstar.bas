scriptID("$Id: kstar.bas,v 1.1 2011/05/04 21:32:07 bulmer Exp $")

# See MAINTENANCE NOTE, below.

corsica  # Put corsica on top

# Initialize defaults for kstar() arguments...
integer confirm_default      =  -1
real    roll_off_default     =   0
integer make_inverse_default =   0
integer msrf_default         = 101
real    thetac_default       =   0.01

# Initialize global control parameters...
logical use_eqdsk_psi  = false
logical use_eqdsk_ssep = true
logical use_eqdsk_cc   = true
real fbd_wts(2)        = [1000, 1]

logical     showMarkers = false
integer     kstar_plot  = [ 1, 1, 1, 1, 1, 1, 1, 1]  # 0 to turn off
integer     kstar_save  = [ 1, 1]  # == 0 to disable save-file writes
chameleon   cCORSICA    = "orangered"
chameleon   dCORSICA    = " (Corsica red)"
ezccgmc                 = 10000  # Allow big ncgm files

# Internal parameter(s)...
real dsep_flag = 39  # EQDSK value for SSEP "undefined"
real dr_factor = 2
integer which_way = 0
character*1 a_prefix = "a"

# Internal macros...
  MDEF Comment()=
    << trim($1)
  MEND
  MDEF StatusReport()=
    if (debug > 0) then
      << "<"//trim($1)//">"
    endif
  MEND

MDEF kstar() = # Interface to function kstarf()

  character*128 gname
  integer confirm, make_inverse, new_msrf, exceptions = 0
  real roll_off, new_thetac

  IFELSE($1,) (gname        = "help",               gname        = $1)
  IFELSE($2,) (confirm      = confirm_default,      confirm      = $2)
  IFELSE($3,) (roll_off     = roll_off_default,     roll_off     = $3)
  IFELSE($4,) (make_inverse = make_inverse_default, make_inverse = $4)
  IFELSE($5,) (new_msrf     = msrf_default,         new_msrf     = $5)
  IFELSE($6,) (new_thetac   = thetac_default,       new_thetac   = $6)

  if (trim(gname) == "help") then
    if (roll_off_default == 0) then
      chameleon roll_off_str = "0"
    elseif (roll_off_default == 1) then
      chameleon roll_off_str = "1"
    else
      integer nnn = int(-log10(1 - roll_off_default - errt))
      chameleon roll_off_str = format(roll_off_default,0,nnn,1)
    endif
    <<"==============================================================================="
    <<" "
    <<"Generate a Corsica equilibrium to match a KSTAR EFIT equilibrium defined in a"
    <<"pair (both g-file and a-file) of EQDSK files. Calling sequence:"
    <<" "
    <<"   kstar(g-file; confirm, roll_off, make_inverse, msrf, theta)"
    <<" "
    <<"where "";"" signifies the remaining arguments are optional. The argument"
    <<"definitions are..."
    <<" "
    <<"   g-file       : EQDSK file name"
    <<"   confirm      : integer, controls hold-time after graphics frame display..."
    <<"                  >  0 pause time (in seconds)"
    <<"                  == 0 do not display graphics"
    <<"                  <  0 pause until user responds [default: " \
                         <<confirm_default <<"]"
    <<"   roll_off     : 0 (or 1) to use the EQDSK profiles as-is, or roll-off"
    <<"                  profiles smoothly to zero where psibar={roll_off:1} [" \
                          <<roll_off_str <<"]"
    <<"   make_inverse : non-zero to create an inverse equilibrium in addition to"
    <<"                  direct-solve equilibrium [" \
                          <<make_inverse_default <<"]"
    <<"   msrf         : specify number of flux surfaces [" \
                          <<msrf_default <<"]"
    <<"   thetac       : specify thetac for inverse solver [" \
                          <<format(thetac_default,0,3,2) <<"]."
    <<" "
    <<"Default argument values may be altered in a session by specifying a value for"
    <<"<argname>_default, e.g., ""confirm_default = 0""."
    <<" "
    <<"Other controls available to the user are logical (true | false) variables"
    <<"""use_eqdsk_cc"", ""use_eqdsk_psi"" and ""use_eqdsk_ssep"" and real array"
    <<"""fbd_wts(2)"" which holds the values of the weighting factor applied to the"
    <<"plasma boundary for (1) no EQDSK coil currents available and (2) EQDSK coil"
    <<"currents available."
    <<" "
    <<"==============================================================================="
    exceptions = exceptions + 1
  elseif (basopen(gname,"i") == -1) then
    warning("Cannot open EQDSK file: "//trim(gname))
    exceptions = exceptions + 1
  else # Read EQDSK file(s)
    call reqdskg(gname, "none")
    chameleon gname = trim(gname)
    chameleon aname = a_prefix//substr(gname,2,strlen(trim(gname))-1)
    if (basopen(aname, "i") == 0) then
      call reqdska(aname, "none")
    else
      forget aname
    endif
  endif
  logical restore_now = true
  if (exceptions > 0 | nprob > 0) then  # User started-up with save-file
    restore_now = false
  else  # Choose generic reference case
    # MAINTENANCE NOTE: If a new KSTAR configuration is issued, add coding in
    # this section to restore the appropriate Corsica save-file.
    if (g_shot < 1938) then
      fatal("No Corsica save-file for shot numbers less than 1938")
    elseif (g_shot >= 1938 & g_shot < 4015) then
      chameleon restore_file = "kstar20091016.sav"
    elseif (g_shot >= 4015) then
      chameleon restore_file = "kstar20100815.sav"
    endif
  endif
  if (restore_now) then
    <<return<<"Restoring reference case: "<<restore_file
  endif
  if (restore_now) read {restore}.bas
  if (exceptions == 0) then
    call kstarf(gname,confirm,roll_off,make_inverse,new_msrf,new_thetac)
  endif
  forget i, confirm, roll_off, make_inverse, new_msrf, new_thetac, exceptions

MEND # kstar

function kstarf(gname;confirm,roll_off,make_inverse,new_msrf,new_thetac)

  default(confirm)      = confirm_default
  default(roll_off)     = roll_off_default
  default(make_inverse) = make_inverse_default
  default(new_msrf)     = msrf_default
  default(new_thetac)   = thetac_default

  # Initialize graphics frame hold-time...
  call pause(confirm)

  # Change number of flux surfaces.
  if (new_msrf <> msrf) then
    note = "Changing msrf to "//format(new_msrf,0)
    msrf = new_msrf
    run
  endif
  
  # Must use spline coefficients and not polys...
  if (filter_wp == -2) filter_wp = -1
  if (filter_wp == 2) filter_wp =  1
  if (filter_p == -2) filter_p = -1
  if (filter_p == 2) filter_p =  1

  # Get shot number and time point...
  shotName = format(g_shot,-6)
  shotTime = g_time*1.0e-03  # Convert ms to s

  # Set shot id string and file names...
  chameleon shotid = format(g_shot, -6)//"."//format(g_time,-5)
  chameleon ccname = "cc"//shotid
  chameleon fileid = format(g_shot, -6)//"_"//format(g_time,-4)
  chameleon save_file = trim(fileid)//".sav"
  chameleon inv_save_file = trim(fileid)//"_inv.sav"

  # Set problem identification...
  call set_kstar_id(shotName,shotTime)

  # Set-up for running...
  package eq
  filter_wp=0; filter_p=0; nht=200; epsrk=1.e-6

  # Re-grid, if necessary...
  integer regrid_flag = 0
  if (jm <> g_nw) regrid_flag = regrid_flag + 1
  if (km <> g_nh) regrid_flag = regrid_flag + 1
  if (regrid_flag > 0) then
    if (liml > 0) then
      integer ldiff = int((log(jm - 1) - log(liml - 1) + 0.001)/log(2))
      liml = 2**(int((log(jm - 1) + 0.001)/log(2)) - ldiff) + 1
      limlo = 0
    endif
    jm = g_nw
    km = g_nh
    note = "Re-gridding with jm = "//format(jm,0)//", km = "//format(km,0)
    run
  endif

  # Convert units...
  g_fpol = g_fpol*(0.01*ro*g_bcentr)/(g_fpol(g_nw))
  btor = g_bcentr*1.0e+04  # Convert T to G
  plcm = g_cpasma*1.0e-06  # Convert A to MA
  real dpsi_eqdsk = (g_ssimag - g_ssibry)*1.0e+08  # Convert Wb to Maxwells

  # Construct EFIT grid dimensions...
  real g_r(g_nw) = g_rgrid1 + g_xdim*float(iota(g_nw) - 1)/(g_nw - 1)
  real g_z(g_nh) = g_zmid + g_zdim*(float(iota(g_nh) - 1)/(g_nh - 1) - 0.5)

  # Map EFIT boundary to Corsica fuzzy marker arrays...
  nfbd = g_nbbbs
  alfbd = max(fbd_wts)
  rfbd = 100*g_rbbbs  # Convert m to cm
  zfbd = 100*g_zbbbs  # Convert m to cm

  # Map limiter coordinates to Corsica r,zplate...
  nplates = g_limitr
  call teq_sync
  nplates = 0
  integer i
  do i = 2,g_limitr
    if ((g_xlim(i) - g_xlim(i-1))**2 + (g_ylim(i) - g_ylim(i-1))**2 > 0) then
      nplates = nplates + 1
      rplate(nplates,1) = g_xlim(i-1)
      zplate(nplates,1) = g_ylim(i-1)
      rplate(nplates,2) = g_xlim(i)
      zplate(nplates,2) = g_ylim(i)
    endif
  enddo
  call teq_sync
  rplate = 100*rplate  # Convert m to cm
  zplate = 100*zplate  # Convert m to cm

  # Use consistent sign convention...
  global real sign_conv = sign(1.0, g_cpasma)
  if (g_cpasma/plcm < 0) then
    note = "Change sign of current, flux, etc..."
    cursign
  endif

  # Load profiles from EQDSK file...
  StatusReport("Mapping EQDSK profiles to Corsica arrays")
  global real x_eqdsk = float(iota(g_nw) - 1)/float(g_nw - 1)
  real v_eqdsk(g_nw), w_eqdsk(g_nw)
  betaj = (0.01*ro)**2/(1 + (0.01*ro)**2)
  v_eqdsk = 10*g_pres/betaj/(sign_conv*dpsi_eqdsk)
  w_eqdsk = (g_fpol**2 - g_fpol(g_nw)**2)*1e12/(8*pi*(sign_conv*dpsi_eqdsk))
  call interp(x_eqdsk,v_eqdsk,g_nw,splbdry,splbdry,psibar,&v,&vp,msrf)
  call interp(x_eqdsk,w_eqdsk,g_nw,splbdry,splbdry,psibar,&w,&wp,msrf)
  call spline(psibar,v,msrf,splbdry,splbdry,&vp,&vpp,&vppp,-999)
  call spline(psibar,w,msrf,splbdry,splbdry,&wp,&wpp,&wppp,-999)
  ipscl = 0 # Use EQDSK plasma current
  ipj = 0; ipp = 999; ipf = 999

  global real dsep_eqdsk
  if (use_eqdsk_psi) then
    StatusReport("Mapping EQDSK psi(R,Z) to Corsica R,Z grid")
    cps(,1) = -sign_conv*shape(g_psirz,jmkm)*1.0e+08
    # Get spline coefficents on EQDSK grid
    integer imeth = 11
    call bispline(100*g_r,100*g_z,g_nw,g_nh,&cps,imeth,n4)
    # Get psi on Corsica grid
    call boxm(&rm,&zm,&psi,jmkm,0,1,n4)
    # Set flux
    psi00  = psi
    cps(,1)= psi
    # Get spline coefficients on Corisca grid
    call bispline(r,z,jm,km,&cps,imeth,n4)
    # Evaluate topology...
    wvert = 1
    if (exists("a_limloc")) then
      if (a_limloc == "DN") rl = 1
    endif
    call separatrix
    rl = 0
    if (limiterd == 0) then
      dsep_eqdsk = sign(dsep, zxpt)
    endif
  endif

  # Fix profiles near edge...
  if (roll_off > 0 & roll_off < 1) then
    StatusReport("Modifying profiles near edge")
    filter_p=0
    integer mroll = mnx(abs(psibar - roll_off))
    cos0 = 1
    cos0(mroll:) = (1-(psibar(mroll:)-psibar(mroll))**10/(1-psibar(mroll))**10)
    vp = vp*cos0
    call spline(psibar,vp,msrf,splbdry,splbdry,v,v,v,&v)
    v = v - v(msrf) + v_eqdsk(g_nw)
    cos0 = 1
    cos0(mroll:) = (1-(psibar(mroll:)-psibar(mroll))**10/(1-psibar(mroll))**10)
    wp = wp*cos0
    call spline(psibar,wp,msrf,splbdry,splbdry,w,w,w,&w)
    w = w - w(msrf)
    call spline(psibar,v,msrf,splbdry,0,&vp,&vpp,&vppp,-999)
    call spline(psibar,w,msrf,splbdry,0,&wp,&wpp,&wppp,-999)
    ipscl = -1  # Let plasma current float since profiles modified
    entropy_flag=0
  endif

  prsx = 1.0001

  # Initialize coil currents...
  StatusReport("Initializing coil currents")
  if (min(ic) == 0) then
    fatal("Array ic has one or more zeros")
  endif
  ircwt = 1
  cc0 = 0
  cc = ic*1e-3  # To preserve coils in series
  alfbd = max(fbd_wts)

  nbd = 2
  rbd = [rfbd(mxx(rfbd)),rfbd(mnx(rfbd))]
  zbd = [zfbd(mxx(rfbd)),zfbd(mnx(rfbd))]
  if (epsb <> 1) then
    ixpt = 1; alphac = 0.7; nl = 301; ni = 1; epsj = 1.0e-06; epsb = 1
  endif
  kaxis = 1
  zaxis0 = 100*g_zmaxis  # Convert m to cm
  raxis0 = 100*g_rmaxis  # Convert m to cm
  kaxis = 0  # For better convergence (in general)
  rlim(0) = min(rplate)  # In case of limited equilibrium
  zlim(0) = zfbd(mnx(rfbd))
  if (exists("aname")) then
    note = "Constructed from "//a_prefix//","//trim(gname)
  else
    note = "Constructed from "//trim(gname)
  endif
  wvert=1

  # Refine constraints...
  irl = 0; rl = 0
  if (exists("aname")) then # Use info from a-file
    StatusReport("EQDSK limloc: "//trim(a_limloc))
    if (a_limloc == "SNB") then
      zxpr(2) = max(zplate)
    elseif (a_limloc == "SNT") then
      zxpr(1) = min(zplate)
    endif
    if (exists("a_rseps")) then
      StatusReport("Adding 2 r,zbd constraints near midplane")
      nbd = 2
      rbd = [rfbd(mxx(rfbd)),rfbd(mnx(rfbd))]
      zbd = [zfbd(mxx(rfbd)),zfbd(mnx(rfbd))]
    else # Limited
      StatusReport("Adding fixed marker at OB edge")
      nbd = 1
      rbd(1) = rfbd(mxx(rfbd))
      zbd(1) = zfbd(mxx(rfbd))
    endif
    if (use_eqdsk_cc) then
      StatusReport("Installing EQDSK coil currents as target values")
      if (getcc) then
        alfbd = min(fbd_wts)
      endif
    endif
    if (dsep > 0 & abs(a_ssep) < dsep_flag) then
      if (a_limloc == "DN") then
        StatusReport("Adding DN constraint")
        irl = 0
        rl = sign(1.0e-06,-a_ssep)
      elseif (use_eqdsk_ssep) then
        StatusReport("Adding EQDSK SSEP constraint")
        rl = abs(a_ssep)
        irl = sign(1, -a_ssep)
        if (dsep_eqdsk <> 0) then
          if (max(a_ssep/dsep_eqdsk, dsep_eqdsk/a_ssep) > 2) then
            chameleon msg = "Large difference in a_ssep and dsep_eqdsk"
            msg = msg//", ignoring constraint"
            warning(msg)
            rl = abs(dsep_eqdsk)
            if (rl <> 0) irl =sign(1, -dsep_eqdsk)
          endif
        endif
      endif
    endif
  else # Add a hard marker at inboard edge if none exist...
    if (nbd == 0) then
      StatusReport("Adding fixed marker at outboard edge")
      nbd = 1
      rbd(1) = rfbd(mxx(rfbd))
      zbd(1) = zfbd(mxx(rfbd)) - 0.1 # to avoid limiter point
    endif
  endif

  run

  if (residj > epsj) then
    if (residj < 50*epsj) then
      warning(trim(fileid)//" not fully converged")
    else
      warning(trim(fileid)//" has poor convergence")
    endif
  endif

  if (ipscl <> 0) then # Plasma currents will differ...
    remark " "
    remark "Plasma current comparison (ipscl == -1)..."
    remark "EQDSK:   "//format(int(g_cpasma),8)//" A"
    remark "Corsica: "//format(int(10*plc),8)//" A  (" \
    //format(100*(10*plc-g_cpasma)/g_cpasma,0,2,1)//" %)"
  endif

  if (which_way == -1) then # direct equilibrium with thetac
    thetac=new_thetac
    Comment(char(10)//"Generating direct equilibrium with thetac=" \
    //format(thetac,0,2,0))
    run
  endif

  if (sum(kstar_plot(1:8)) <> 0 & confirm <> 0) then
    Comment(char(10)//"Check equilibrium...")
    win on "kstar-comparison"
  endif
  if (kstar_plot(1) <> 0) then
    pbg; pause("Plasma boundary")
  endif
  if (kstar_plot(2) <> 0) then
    ppprime; pause("Pressure")
  endif
  if (kstar_plot(3) <> 0) then
    pffprime; pause("Toroidal flux function")
  endif
  if (kstar_plot(4) <> 0) then
    pq; pause("q-profile")
  endif
  if (kstar_plot(5) <> 0) then
    if (cc0!cc0 > 0) then
      pcoils; pause("Coil currents")
    endif
  endif
  if (kstar_plot(6) <> 0) then
    layout; pause("Configuration")
  endif

  # Save the direct-solve equilibrium...
  if (kstar_save(1) <> 0) then
    saveq(save_file)
    Comment(char(10)//"Created equilibrium save-file: "//trim(save_file))
  endif

  # Use these to restore direct equilibrium solution...
  real v_save = [v, vp, vpp, vppp]
  real w_save = [w, wp, wpp, wppp]

  if (make_inverse <> 0) then # Make an inverse equilibrium save-file...
    thetac=new_thetac
    if (confirm <> 0) then
      Comment(char(10)//"Generating inverse equilibrium with thetac=" \
      //format(thetac,0,2,0))
    endif
    morph_inv(which_way)
    
    if (kstar_plot(3) <> 0) then
      pffprime; pause("Toroidal flux function")
    endif
    if (kstar_plot(4) <> 0) then
      pq; pause("q-profile")
    endif
    if(kstar_plot(7) <> 0) then
      contour; pause("Configuration")
    endif
    if(kstar_plot(8) <> 0) then
      if (which_way == 0) then
        testGS(3); pause("GS error estimate")
      else
        testGS(1); pause("GS error estimate")
      endif
    endif
    if(kstar_plot(8) <> 0) then
      Comment("Testing q-profile input to inverse solver")
      teqinv(0,0)
      testGS(0); pause("GS error estimate")
    endif
    if (kstar_save(2) <> 0) then
      saveq(inv_save_file)
      Comment("Created inverse-equilibrium save-file: "//trim(inv_save_file))
    endif
    
    Comment("Restoring direct equilibrium save-file...")
    restore ^save_file
    inv_eq = 0
    v = v_save(,1); vp = v_save(,2); vpp = v_save(,3); vppp = v_save(,4)
    w = w_save(,1); wp = w_save(,2); wpp = w_save(,3); wppp = w_save(,4)
    run

  endif
  
  if (sum(kstar_plot) <> 0 & confirm <>0) then
    win close "kstar-comparison"
  endif

endf  # kstarf

function kstarplots(;p1,p2,p3,p4,p5,p6,p7,p8)  # Plot selection for morphing

  default(p1) = kstar_plot(1)
  default(p2) = kstar_plot(2)
  default(p3) = kstar_plot(3)
  default(p4) = kstar_plot(4)
  default(p5) = kstar_plot(5)
  default(p6) = kstar_plot(6)
  default(p7) = kstar_plot(7)
  default(p8) = kstar_plot(8)

  if (type(p1) <> "integer") then
    remark " "
    remark "Script function ""kstarplots(i1,...,in)"" sets and displays the current plot"
    remark "activation list for the morphing routine, kstar. The arguments are a list of"
    remark "integers: 1 to turn on the ith plot and 0 to turn off the ith plot. Arguments"
    remark "default to their present value. Execute kstarplots with no arguments to get the"
    remark "current settings."
    remark " "
    return
  endif

  kstar_plot = [p1,p2,p3,p4,p5,p6,p7,p8]
  kstar_plot = where(kstar_plot <> 0, 1, kstar_plot)

  remark " "
  remark "Plot activation via kstar_plot(1:8)..."
  remark " i content       kstar_plot(i)"
  remark " "
  remark " 1 Shape           "//format(kstar_plot(1),3)
  remark " 2 p & p' profiles "//format(kstar_plot(2),3)
  remark " 3 F & FF' profiles"//format(kstar_plot(3),3)
  remark " 4 q profile       "//format(kstar_plot(4),3)
  remark " 5 Coil currents   "//format(kstar_plot(5),3)
  remark " 6 Configuration   "//format(kstar_plot(6),3)
  remark " 7 Inverse eq.     "//format(kstar_plot(7),3)
  remark " 8 Inv. details    "//format(kstar_plot(8),3)
  remark " "
  remark "1/0 = ON/OFF"
  remark " "

endf  # kstarplots

##################    I N T E R N A L    R O U T I N E S    ###################

  function getcc  # Get EQDSK coil currents

    if (exists("aname") & a_nesum > 0 & a_nfcoil > 0) then
      if (nc <> a_nfcoil) then
        fatal("EQDSK circuit currents incompatible with device configuration")
      endif
      cc(1:nc) = ntc(1:nc)*a_ccbrsp*1e-6
      cc0 = cc
      return true
    else
      warning("*** EQDSK coil currents not available ***")
      return false
    endif

  endf  # getcc

  function pcoils(;num)  # Plot coil currents

    default(num) = npfc

    integer i
    real xtext = 0.21, dxtext = 0.0543, ytext = 0.07
    real y_min = min(min(cc(1:num)),min(cc0(1:num)))
    real y_max = max(max(cc(1:num)),max(cc0(1:num)))
    real y_rng = y_max - y_min
    y_min = y_min - 0.05*y_rng
    y_max = y_max + 0.10*y_rng

    real cc_dev = (cc(1:num) - cc0(1:num))/max(abs(cc0(1:num)))
    chameleon cnote = format(100*rmsdv(cc_dev),0,2,1)
    cnote = cnote//"% Corsica-EFIT r.m.s. deviation"

    chameleon pnote = "(Passive structure currents range from "
    pnote = pnote//format(min(cc(17:))*1e6,0,0,1)//" to "
    pnote = pnote//format(max(cc(17:))*1e6,0,0,1)//" A)"
    
    oframe
      call ezcstxqu(1)
      titles probid " " "Coil current [MAt]"//dCORSICA
      frame 0 num+1 y_min y_max
      plot [0,0],[0,num+1] color=gray70
      plot cc0(1:num),iota(num) thick=3
      plot cc0(1:num),iota(num) mark=circle marksize=2
      plot cc(1:num),iota(num) thick=3 style=dotted color=cCORSICA
      if (num == npfc) then
        ftext "Coil Index (driven coils only)" 0.57 ytext 10 0 0
      else
        ftext "Coil Index" 0.57 ytext 10 0 0
      endif
      ftext trim(pnote) 0.57 0.04 2 0 0
      ftext trim(cnote) 0.57 0.90 2 0 0
    cframe
    call ezcstxqu(0)

  endf  # pcoils

  function ppprime  # Plot p and p' comparison

    real y_min, y_max, y_rng
    real y_factor = 0.05
    oframe
      call ezcstxqu(1)
      ezctitfr = 0.06
      call ezcquad(12)
      y_min = min(min(g_pres),min(0.1*prsrf))
      y_max = max(max(g_pres),max(0.1*prsrf))
      y_rng = y_max - y_min
      y_min = y_min - y_factor*y_rng
      y_max = y_max + y_factor*y_rng
      titles probid " " "p  "//dCORSICA
      frame 0 1.1 y_min y_max
      plot 0.0*ones(msrf),psibar color=gray70
      plot [y_min,y_max],[1,1] color=gray70
      plot g_pres,x_eqdsk mark=cross
      if (showMarkers) then
        plot 0.1*prsrf,psibar color=cCORSICA mark=circle
      else
        plot 0.1*prsrf,psibar color=cCORSICA
      endif

      call ezcquad(34)
      y_min = min(min(-sign_conv*g_pprime),min(pdsrf*1e7))
      y_max = max(max(-sign_conv*g_pprime),max(pdsrf*1e7))
      y_rng = y_max - y_min
      y_min = y_min - y_factor*y_rng
      y_max = y_max + y_factor*y_rng
      titles " " "psibar" "p'  "//dCORSICA
      frame 0 1.1 y_min y_max
      plot 0.0*ones(msrf),psibar color=gray70
      plot [y_min,y_max],[1,1] color=gray70
      plot -sign_conv*g_pprime,x_eqdsk mark=cross
      if (showMarkers) then
        plot pdsrf*1.0e+07,psibar color=cCORSICA mark=circle
      else
        plot pdsrf*1.0e+07,psibar color=cCORSICA
      endif
    cframe
    call ezcstxqu(0)

  endf  # ppprime

  function pffprime  # Plot F & FF' comparison

    real y_min, y_max, y_rng

    real f_corsica = sign(1.0,btor)*frsrf*1.0e-06
    if (btor < 1) then
      chameleon y_label = "-F "//dCORSICA
    else
      chameleon y_label = "F  "//dCORSICA
    endif

    oframe
      call ezcstxqu(1)
      ezctitfr = 0.06
      call ezcquad(12)
      titles probid " " "F  "//dCORSICA
      real y_min = min(min(g_fpol),min(f_corsica))
      real y_max = max(max(g_fpol),max(f_corsica))
      real y_rng = y_max - y_min
      y_min = y_min - 0.05*y_rng
      y_max = y_max + 0.05*y_rng
      frame 0 1.1 y_min y_max
      plot 0.0*ones(msrf),psibar color=gray70
      plot [y_min,y_max],[1,1] color=gray70
      plot g_fpol,x_eqdsk mark=cross
      if (showMarkers) then
        plot f_corsica,psibar color=cCORSICA mark=circle
      else
        plot f_corsica,psibar color=cCORSICA
      endif

      call ezcquad(34)
      titles " " "psibar" "FF'  "//dCORSICA
      real y_min = min(min(-sign_conv*g_ffprim),min(frsrf*fpsrf*1e-4))
      real y_max = max(max(-sign_conv*g_ffprim),max(frsrf*fpsrf*1e-4))
      real y_rng = y_max - y_min
      y_min = y_min - 0.05*y_rng
      y_max = y_max + 0.05*y_rng
      frame 0 1.1 y_min y_max
      plot 0.0*ones(msrf),psibar color=gray70
      plot [y_min,y_max],[1,1] color=gray70
      plot -sign_conv*g_ffprim,x_eqdsk mark=cross
      if (showMarkers) then
        plot frsrf*fpsrf*1.0e-04,psibar color=cCORSICA mark=circle
      else
        plot frsrf*fpsrf*1.0e-04,psibar color=cCORSICA
      endif
    cframe
    call ezcstxqu(0)

  endf  # pffprime

  function pq  # Plot q-profile comparison

    oframe
      call ezcstxqu(1)
      titles probid "psibar" "q  "//dCORSICA
      frame 0 1.1 0 1.05*max(max(g_qpsi),max(qsrf))
      plot [0, 1.05*max(max(g_qpsi),max(qsrf))],[1,1] color=gray70
      plot g_qpsi,x_eqdsk mark=cross
      if (showMarkers) then
        plot qsrf,psibar color=cCORSICA mark=circle
      else
        plot qsrf,psibar color=cCORSICA
      endif
    cframe
    call ezcstxqu(0)

  endf  # pq

  function set_kstar_id(name, time)

    probid = "KSTAR #"//trim(name)//" @ "// format(time,0,3,1)//" s"
    if (exists("a_limloc")) then
      probid = trim(probid)//" ["//trim(a_limloc)//"]"
    endif

  endf  # set_kstar_id
