FUNCTION GKV_RHS, X, Y, gamma_d=gamma__d, nu_z=nu__z
;
; COMPUTE the RHS for integrating the system of 
; differential equations describing 4-wave 
; problem from L.Chen, Z. Lin, and R. White, 
; Phys. Plasmas <7>, 3129 (Aug. 2000)
;
; Written by W.M. Nevins
;	5/4/04
COMMON LSODE, gamma_d, nu_z, delta

ydot = FLTARR(4)

ydot[0] =           y[0] - 2.*y[1]*y[2]*COS(y[3])
ydot[1] = - gamma_d*y[1] +    y[0]*y[2]*COS(y[3])
ydot[2] =    - nu_z*y[2] + 2.*y[0]*y[1]*COS(y[3])
ydot[3] =      delta   - y[0]*y[2]/y[1]*SIN(y[3])
RETURN, ydot
END ; ****** GKV_RHS ****** ;

FUNCTION GKV_LSODE, 	tmax=t_max, P0=P_0, S0=S_0, Z0=Z_0, PSI0=PSI_0,	$
			gamma_d=gamma__d, nu_z=nu__z, delta=delta_in,	$
			nSteps = n_steps
;
; Use LSODE to integrate system of differential equations describing 4-wave 
; problem from L.Chen, Z. Lin, and R. White, Phys. Plasmas <7>, 3129 (Aug. 2000)
;
; Input keywords:
;
;	tmax	Equations are integrtated from t=0 to t=tmax.
;		Defaults to tmax=100. (Optional)
;
;	nSteps	Values of P, S, Z, and PSI will be reported in
;		output objects at "nSteps" equally spaced points
;		between 0 and tmax (actual time step is chosen by
;		LSODE).  Defaults to 1000. (Optional)
;
;	gamma_d	Parameter of the same name in Chen et al, corresponding
;		to the ratio of the side-band damping rate to the linear
;		growth rate of  the pump.  In my experinece, the system
;		'runs-away' if gamma_d is less than (or equal to) 1. 
;		However, Chen et al suggest that lower values of gamma_d
;		can lead to stable fixed point solutions (see Fig. 4 of 
;		Chen et al).  Defaults to 2.0 (Optional)
;
;	delta	Parameter of the same name in Chen et al, corresponding
;		to the ratio of the frequency mis-match (between the pump
;		and the sideband) to the linear growth rate of  the pump.  
;		Defaults to 2.0 (Optional)
;
;	nu_z	Parameter of the same name in Chen et al, corresponding
;		to the ratio of the damping rate of the zonal flow to the
;		linear growth rate of the pump.  Defaults to 1.e-2, but
;		larger values should be used if you desire to integrate
;		in time until the solution reaches the fixed point.
;		(Optional)
;
;	P0	Initial value of pump amplitude.  Defaults to 0.1.
;		(Optional)
;
;	S0	Initial value of the side-band amplitude.
;		Defaults to 0.01. (Optional)
;
;	Z0	Initial value of the zonal flow amplitude.
;		Defaults to 0.01. (Optional)
;
;	PSI0	Initial value of the relative phase betwen the 
;		pump and the side-band.  Defaults to 0s. (Optional)
;
; Written by W.M. Nevins
;	5/4/04
; 
; Create common block to pass parameters to GKV_RHS,
; set default values, and read command line
;
COMMON LSODE, gamma_d, nu_z, delta

tmax=100.
IF(N_ELEMENTS(t_max) EQ 1) THEN tmax = t_max

P0=1.e-1
IF(N_ELEMENTS(P_0)   EQ 1) THEN P0 = P_0

S0=1.e-2
IF(N_ELEMENTS(S_0)   EQ 1) THEN S0 = S_0

Z0=1.e-2
IF(N_ELEMENTS(Z_0)   EQ 1) THEN Z0 = Z_0

PSI0=0.
IF(N_ELEMENTS(PSI_0)   EQ 1) THEN PSI0 = PSI_0

gamma_d=2.
IF(N_ELEMENTS(gamma__d)   EQ 1) THEN gamma_d=gamma__d

nu_z = 1.e-2
IF(N_ELEMENTS(nu__z)     EQ 1) THEN nu_z = nu__z

delta = 2.
IF(N_ELEMENTS(delta_in)     EQ 1) THEN delta = delta_in

nSteps = 1000
IF(N_ELEMENTS(n_steps) EQ 1) THEN nSteps = n_Steps > 10
;
; Create time vector
;
dt = tmax/(nSteps-1)
t = dt*FINDGEN(nSteps)
;
; Create output vectors
;
P   = FLTARR(nSteps)
S   = FLTARR(nSteps)
Z   = FLTARR(nSteps)
PSI = FLTARR(nSteps)
;
; Initialize them
;
P[0]   = P0
S[0]   = S0
Z[0]   = Z0
PSI[0] = PSI0
;
; Perform first call to LSODE
;
yold = FLTARR(4)
yold[0] = P[0]
yold[1] = S[0]
yold[2] = Z[0]
yold[3] = PSI[0]

status = 1
ynew = LSODE(yold, t[0], dt, "GKV_RHS", status)
print, 'status = ', status
status = 2
FOR i=1, nsteps-1 DO BEGIN
	IF(status LT 0) THEN GOTO, DONE
	P[i]   = ynew[0]
	S[i]   = ynew[1]
	Z[i]   = ynew[2]
	PSI[i] = ynew[3]
	yold = ynew
	ynew = LSODE(yold, t[i], dt, "GKV_RHS", status)
ENDFOR
DONE: 

tGrid = {GRID}
tGrid.mnemonic = 't'
tGrid.title = 't'
tGrid.values = PTR_NEW(t)
tGrid.boundary = "Open"
tGrid.uniform = 1b
tGrid.range = [0,tmax]
tGrid.irange = [0,nSteps-1]

Pstr = {GKVs1D}
Pstr.title = "P"
Pstr.mnemonic = "P"
Pstr.indices = PTR_NEW(['*'])
Pmin = MIN(P, MAX=Pmax)
Pstr.vrange = [Pmin, Pmax]
Pstr.values = PTR_NEW(P)
Pstr.Grid1 = tGrid

Pobj = OBJ_NEW("GKVs1D", Pstr)

Sobj = Pobj -> MAKECOPY(/novalues)
Smin = MIN(S, MAX=Smax)
Sobj -> Set, title='S', mnemonic='S', values=PTR_NEW(S), vrange=[Smin, Smax]

Zobj = Pobj -> MAKECOPY(/novalues)
Zmin = MIN(Z, MAX=Zmax)
Zobj -> Set, title='Z', mnemonic='Z', values=PTR_NEW(Z), vrange=[Zmin, Zmax]

PSIobj = Pobj -> MAKECOPY(/novalues)
sinPSI = SIN(PSI)
PSI = PSI MOD (2.*!PI)
PSI = PSI - 2.*!PI*(PSI GT !PI)
PSIobj -> Set, title='!4W!X', mnemonic='PSI', values=PTR_NEW(PSI), vrange=[-!PI,!PI]

SinPSIobj = PSIobj -> MakeCopy(/NoValues)
SinPSIOBj -> Set, title='Sin !4W!X', mnemonic='Sin_PSI', values = PTR_NEW(sinPSI), vrange=[-1,1]
;
; Display results
;
vMin = MIN([Pmin, Smin, Zmin, -!PI])
vMax = MAX([Pmax, Smax, Zmax,  !PI])
Pobj -> View, /pretty, charsize=1.4, vrange=[vmin, vmax], t=[0,tmax], Sobj, Zobj, SinPSIobj
PRINT, "BLACK = Pump"
PRINT, "RED   = Sideband"
PRINT, "BLUE  = Zonal Flow"
PRINT, "GREEN = Sin(PSI)"
;
; compute fixed-points of mapping
;
Zfix = SQRT( (delta^2 + gamma_d^2)/(2.*gamma_d) )
Pfix = SQRT(nu_z)*Zfix
Sfix = Pfix/SQRT(2.*gamma_d)
SIN_PSIfix = delta/SQRT(delta^2 + gamma_d^2)
PSIfix = ASIN(SIN_PSIfix)
;
; Form output structures
;
initialValues ={P0	:  P0		,	$
		S0	:  S0		,	$
		Z0	:  Z0		,	$
		PSI0	:  PSI0			}

inputs = {	tmax	:  tmax		,	$
		nSteps	:  nSteps	,	$
		delta	:  delta	,	$
		gamma_d	:  gamma_d	,	$
		nu_z	:  nu_z			}

fixedPoints = { P	:  Pfix		,	$		
		S	:  Sfix		,	$
		Z	:  Zfix		,	$
		PSI	:  PSIfix	,	$
		SIN_Psi	:  SIN_PSIfix		}


	
output = {	Name	:  "GKV_LSODE"	,	$
		inputs	:  inputs	,	$
	InitialValues	:  InitialValues,	$
	fixedPoints	:  fixedPoints	,	$
		P	:  Pobj		,	$
		S	:  Sobj		,	$
		Z	:  Zobj		,	$
		PSI	:  PSIobj	,	$
		sinPSI  :  SINPSIobj	,	$
		Vmin	:  vMin		,	$
		Vmax	:  vMax			}

RETURN, output

END ; ****** GKV_LSODE ****** ;
