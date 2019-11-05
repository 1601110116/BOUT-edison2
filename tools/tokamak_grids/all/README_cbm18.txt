IDL>  pdb2bout, "../gato/cbm18_dens8.dskgato.pdb", output="cbm18_dens8.grid.nc"
Loading input file...
**Last poloidal point duplicates the first. Removing...
Calculating poloidal field using SVD...
Difference in Bpxy computed two different ways:     0.096065331
Maximum percentage difference:        92.059741
Using given Bp, Br and Bz values
Calculating toroidal field from input f
Difference in Bt computed two different ways:   9.5367432e-07
Maximum percentage difference:    2.3790789e-05
Using Bt from input file
***Maximum mu0p is       30762.3
Is this pressure (not mu0*pressure)?y
Grid contains vacuum region
PSI normalised range: 0 to 1.44000
Number of radial grid points  : 259
Number of poloidal grid points: 512
Edge Q: 2.96531
Range of parallel current density [A/m^2]:       -774721.75
       359325.12
Inner Psi boundary:0.4
Outer Psi boundary:1.2
Psi range      0.400000      1.20000
Number of radial grid points:          100
Of which inside the plasma:           80
Is this range ok?y
====== SETTING PLASMA PROFILES =======
Some plasma parameters given in input file
Use given parameters?y
Got density, setting temperature (Ti = Te)
Maximum temperature (eV):      960.120
Generating plasma profiles:
  1. Flat temperature profile
  2. Flat density profile
  3. Te proportional to density
Profile option:1
Setting flat temperature profile
Temperature (eV):1000
Maximum density (10^20 m^-3):     0.960120
Is this ok?y
Increase radial resolution?y
Number of radial points:516
======== GENERATING ORTHOGONAL COORDINATES FOR BOUT ========= 
Number of poloidal grid points:64
Enter x index of equal hthe [0, 515] :317
Lines done:64 of 64 
Interpolating Rxy
Interpolating Zxy
Is this ok?y
Interpolating values onto new grid
Jpar
Bxy
Bpxy
Btxy
Density 0.0000000 -> 0.88833028
Maximum pressure [Pa]:        28462.103
Add vacuum region?no
Equilibrium correction options:
  0  No correction
  1  RBt using force balance
  2  hthe and RBt using force balance and q (FAILS)
  3  hthe and RBt using force balance and jpar
Enter option:1
Correcting f = R*Bt using radial force balance
% Program caused arithmetic error: Floating underflow
       0force imbalance:      0.035070791->    0.00011982640
      63force imbalance:      0.035480647->    0.00021981727
Maximum change in Bt =     0.0049090046
Maximum percentage change =       0.23831240
Calculating poloidal arc length dthe
Maximum difference in hthe:       0.24193937
Maximum percentage difference:        17.584302
Use new hthe?y
Checking parallel current
****Equilibrium has -ve toroidal field
Maximum difference in jpar0:        138848.07
Maximum percentage difference:        69776.691
Use new Jpar?n
Performing integrals
q is negative. Reversing values from equilibrium
Use new qsafe?y
****Minimum pressure is very small:       0.0000000
****Setting minimum pressure to 1% of maximum
Cannot write 2 dimensional double bxcvz. Writing as float
DONE
