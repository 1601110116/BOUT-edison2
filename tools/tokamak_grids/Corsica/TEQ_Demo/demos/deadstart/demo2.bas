echo = yes
###############################################################################
#
# Starting with the up/down symmetric solution obtained with the deadstart
# procedure executed in demo1.bas, adjust the equilibrium to match desired
# properties (in this case from Garofalo's paper). Start-up with:
#
#     caltrans fdf2.sav demo2.bas
#
# (0) Display initial solution...
      win
      layout(0, 0)


# (1) Next, look at shape constraints with pb or pbg...
      pause

      pbg

      # Note boundary is defined by limiter point L0 and 62 marker points
      # in the +Z domain.


# (2) Next, add two constraints for R0 = 2.49, a = 0.71 and retract limiter
#     point to get an x-point limited solution...
      pause
      
      nbd = 2
      rbd = [249+71, 249-71]
      zbd = 0
      rlim(0) = rlim(0) + 5
      run
      pbg


# (3) Next, install reverse-shear profiles with betaN = 3.5 and (just a guess)
#     qmin = 1.8
      pause

      alfa(0:1) = [1, 2]
      betp(0:1) = 1
      epf = 0.6
      npf = -0.8
      # Redefine ceq problem:
      vo = ["qmin", "ctroy"]
      vo0 = [1.8, 3.5]
      vi = ["epf", "betaj"]
      x0 = [epf, betaj]
      ihy = 20
      run
      profiles


# (4) Next, add x-point constraint...
      pause

      kxp = 1
      rxpt0 = 213.5
      zxpt0 = 164
      alfbd = 1    # Since r,zfbd inappropriate
      ihy = 0
      run
      pbg


# (5) Next, reconverge profiles...
      pause

      ihy = 20  
      run
      profiles


# (6) Next, reset fuzzy boundary points...
      pause

      pls
      nfbd
      mlm
      rfbd = rls(2:mlm-1)
      zfbd = zls(2:mlm-1)
      pbg


# (7) Make some plots and save final state...
      pause

      probid = "FPF at Corsica Summer School"
      layout(0)
      pause

      kextflux = -1  # show external flux
      layout(1,0)  # draw coil sections in proportion to their current
      pause

      profiles
      pause

      saveq("fdf_css.sav")


# (8) Finally, refine the R-Z grid...

      gridup
      run

      gridup
      run

      zoom
      layout

      saveq("fdf_129x257.sav")

      pause

      quit
