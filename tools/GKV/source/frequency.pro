Function GKVs2D::QL_HeatFlux, Q_qlObj
;
; A quick patch to multiply Phi_k_sq_r by
; Q_ql to form quasilinear estimate of heat flux.
;
; Written by W.M. Nevins
;  7/29/07
;
Q_ql_2 = self -> MakeCopy(/NoValues, /NoErrorBars)
phiSq = *self.values
Q_ql  = *Q_qlObj.values
PhiSqInfo = SIZE(phiSq)
nt = PhiSqInfo[phiSqInfo[0]]
nky= PHiSqInfo[1]
values = FLTARR(nky,nt)
FOR it = 0, nt-1 DO values[*,it] = phiSq[*,it]*Q_ql
Q_ql_2.mnemonic = "Q_ql"
Q_ql_2.title = "!12Q!X!DQL!N"
Q_ql_2.units = ""
Q_ql_2.values = PTR_NEW(values)
vMin = MIN(values, MAX=vMax)
Q_ql_2.vrange = [vMin,vMax]
kString = self.Grid1.mnemonic
Q_Ql_1 = Q_QL_2 -> INT(kString, /SUM)
result = {	Name	:	"QL_HeatFlux",	$
		Q_ky	:	Q_QL_2,		$
		Chi_i	:	Q_QL_1		}
Return, result
 
END ; ****** QL_HeatFlux ****** ;

Function GKVs1D::Frequency, tRange=tRange_in
;
; Acts on amplitude (we assume phi) 
; in a k-space (vs. time) representation.
; Computes time-dependent frequency
; by taking the logarimitic derivative.
; This scheme only works if data is from
; a linear run.
; 
;
;  Keywords:
;
;	tRange	The time interval over which
;		the logarithmic dirivative should
;		be averaged to obtain an estimate
;		of the frequency.  Defaults to the 
;		last half of the data-set (optional).
;
; Written by W.M. Nevins
;   7/18/07
;
t=*(self.grid1.values)
values = *self.values
nt = N_ELEMENTS(t)
dt = t[1:nt-1] - t[0:nt-2]
dV = values[1:nt-1] - values[0:nt-2]
Vbar = 0.5*(values[1:nt-1] + values[0:nt-2])
omega = dV/(Vbar*dt)
gamma = FLOAT(omega)
omega = -IMAGINARY(omega)
tbar = 0.5*(t[1:nt-1] + t[0:nt-2])
tUnits = self.grid1.units
wUnits = '1/(' + tUnits + ')'
omegaObj = self -> MakeCopy(/NoValues, /NoErrorBars)
gammaObj = self -> MakeCopy(/NoValues, /NoErrorBars)
omegaObj.title = "!4x!X"
gammaObJ.title = "!4c!X"
omegaObj.mnemonic = "omega"
gammaObj.mnemonic = "gamma"
omegaObj.units = wUnits
gammaObj.units = wUnits
tmin = MIN(tbar)
tmax = MAX(tbar)
PTR_FREE, omegaObj.grid1.values
omegaObj.grid1.values = PTR_NEW(tbar)
omegaObj.grid1.irange = [0,nt-2]
omegaObj.grid1.range = [tmin, tmax]
PTR_FREE, gammaObj.grid1.values
gammaObj.grid1.values = PTR_NEW(tbar)
gammaObj.grid1.irange = [0,nt-2]
gammaObj.grid1.range = [tmin, tmax]
wmin = MIN(omega)
wmax = max(omega)
omegaObj.values = PTR_NEW(omega)
omegaObj.vrange = [wmin, wmax]
gmin = MIN(gamma)
gmax = MAX(gamma)
gammaObj.values = PTR_NEW(gamma)
gammaObj.vrange = [gmin, gmax]
;
; Compute "average" frequency, growth rate
; from average over last half of time interval
;
t0   = t[0]
tmax = t[nt-1]
tmin = 0.5*(t0+tmax)
tRange = [tmin, tmax]
IF(N_ELEMENTS(tRange_in) EQ 2) THEN tRange = tRange_in
temp1 = omegaObj -> Avg(t=tRange)
temp2 = gammaObj -> Avg(t=tRange)
omegaObj -> SignalWindow, t=[t0,tmax]
gammaObj -> SignalWindow, t=[t0,tmax]
omegaValue = temp1 -> GetValues()
omegaPM    = temp1 -> GetErrors()
gammaValue = temp2 -> GetValues()
gammaPM    = temp2 -> GetErrors()
temp1 -> Trash
temp2 -> Trash 

result = {	Name	:	"Frequency",	$	
		omegaAvg:	omegaValue,	$
		omegaPM	:	omegaPM,	$
		gammaAvg:	gammaValue,	$
		gammaPM	:	gammaPM,	$
		omegaObj:	omegaObj,	$
		gammaObj:	gammaObj	}
RETURN, result
END ; ****** Function GKVs1D::Frequency ****** ;


Function GKVs2D::Frequency, KEObj, tRange=tRange_in, deBug=deBug
;
; Acts on amplitude (we assume phi) 
; in a k-space (vs. time) representation.
; Computes time-dependent frequency
; by taking the logarimitic derivative.
; This scheme only works if data is from
; a linear run.
; 
;
; Argument (optional: The k-space (vs. time) representation of the 
; kinetic energy. If supplied, routine will assume that
; 'self' is the poential and routine will compute the
; quasi-linear response relating |phi(k)|^2 to the
; heat flux.
;
;  Keywords:
;
;	tRange	The time interval over which
;		the logarithmic dirivative should
;		be averaged to obtain an estimate
;		of the frequency.  Defaults to the 
;		last half of the data-set (optional).
;
; Written by W.M. Nevins
;   7/18/07
; Modified by W.M. Nevins
;   7/27/07
; to compute QL response functions.
;
t=*(self.grid2.values)
values = *self.values
nx = N_ELEMENTS(*self.grid1.values)
One = MAKE_ARRAY(nx, VALUE=1.)
nt = N_ELEMENTS(t)
dt = One#(t[1:nt-1] - t[0:nt-2])
dV = values[*,1:nt-1] - values[*,0:nt-2]
Vbar = 0.5*(values[*,1:nt-1] + values[*,0:nt-2])
omega = dV/(Vbar*dt)
gamma = FLOAT(omega)
omega = -IMAGINARY(omega)
tbar = 0.5*(t[1:nt-1] + t[0:nt-2])
tUnits = self.grid2.units
wUnits = '1/(' + tUnits + ')'
omegaObj = self -> MakeCopy(/NoValues, /NoErrorBars)
gammaObj = self -> MakeCopy(/NoValues, /NoErrorBars)
omegaObj.title = "!4x!X"
gammaObJ.title = "!4c!X"
omegaObj.mnemonic = "omega"
gammaObj.mnemonic = "gamma"
omegaObj.units = wUnits
gammaObj.units = wUnits
tmin = MIN(tbar)
tmax = MAX(tbar)
PTR_FREE, omegaObj.grid2.values
omegaObj.grid2.values = PTR_NEW(tbar)
omegaObj.grid2.irange = [0,nt-2]
omegaObj.grid2.range = [tmin, tmax]
PTR_FREE, gammaObj.grid2.values
gammaObj.grid2.values = PTR_NEW(tbar)
gammaObj.grid2.irange = [0,nt-2]
gammaObj.grid2.range = [tmin, tmax]
wmin = MIN(omega)
wmax = max(omega)
;
; Coherce omega, gamma =0 at k_y=0.
;
kyvalues = *(OmegaObj.Grid1.values)
epsilon = MIN(kyValues^2, iky0)
omega[iky0,*]=0.
gamma[iky0,*]=0.
;
omegaObj.values = PTR_NEW(omega)
omegaObj.vrange = [wmin, wmax]
gmin = MIN(gamma)
gmax = MAX(gamma)
gammaObj.values = PTR_NEW(gamma)
gammaObj.vrange = [gmin, gmax]
;
; Compute "average" frequency, growth rate
; from average over last half of time interval
;
t0   = t[0]
tmax = t[nt-1]
tmin = 0.5*(t0+tmax)
tRange = [tmin, tmax]
IF(N_ELEMENTS(tRange_in) EQ 2) THEN tRange = tRange_in
OmegaAvg = omegaObj -> Avg(t=tRange)
GAmmaAvg = gammaObj -> Avg(t=tRange)
omegaObj.grid2.iRange = [0,nt-2]
gammaObj.grid2.iRange = [0,nt-2]

result = {	Name	:	"Frequency",	$	
		omegaAvg:	omegaAvg,	$
		gammaAvg:	gammaAvg,	$
		omegaObj:	omegaObj,	$
		gammaObj:	gammaObj	}
N_Args = N_PARAMS()
IF(N_Args EQ 0) THEN RETURN, result
; 
; User supplied KeObj as second argument,
; so compute QL response function for
; heat flux.
;
; We will only use same part of time-
; record used to compute average frequency
;
phi_k = self -> MakeCopy()
phi_k -> SignalWindow, t=trange
phi_k -> Restrict
KE_k = KEobj -> MakeCopy()
KE_k -> SignalWindow, t=trange
KE_k -> Restrict
;
; Fourier transform both Phi_k and
; KE_k in time
;
Phi_k_w = phi_k -> FFT('t')
KE_k_w  = KE_k  -> FFT('t')
;
; Form |Phi|^2 and 
; cross-spectrum between
; dPhi/dy and kE
;
PhiSq_k_w = Phi_k_w -> AbsSq()
KEStar_k_w = KE_k_w -> EXECUTE("CONJ")
KE_k_w -> TRASH
MinusEye = COMPLEX(0.,-1.)
kyString = Phi_k_w.grid1.mnemonic
temp   = phi_k_w -> times(kyString)
dPhidy = temp -> Times(MinusEye)
dPhidySq = dPhidy -> AbsSq()
temp -> Trash
CrossSpect = dPhidy -> Times(KEStar_k_w)
;
; Form QL response function vs. (k, omega)
; and V_exb correlation time vs. (k, omega)
;
temp = CrossSpect -> Over(PhiSq_k_w)
Response_k_w = temp -> Execute("Float")
tauCorr = CrossSpect -> Over(dPhidySq)

IF(N_ELEMENTS(deBug) EQ 0) THEN temp -> Trash
;
; clean up
;
IF(N_ELEMENTS(debug) EQ 0) THEN BEGIN
	CrossSpect -> Trash
	PhiSq_k_w  -> Trash
	dPhidy     -> Trash
ENDIF
; 
; Now evaluate Response_k_w at linear
; frequency. Begin by finding index
; corresponding to omega=0
;
Omega  = *(Response_k_w.grid2.values)
epsilon = MIN(Omega^2, iOmega_0)
;
; Now compute index adjecent to
; linear frequency.
; 
dOmega = Omega[1] - Omega[0]
frequency = -(*(OmegaAvg.values))
nFrequency = frequency/dOmega
iFrequency = FIX(nFrequency)
iOmega = iFrequency + iOmega_0
nOmegaMax = N_ELEMENTS(Omega) - 1
iFrequency = iFrequency > 0
iFrequency = iFrequency < nOmegaMax
;
; Prepare for linear interpolation
;
delta = nFrequency - iFrequency
AbsDelta = ABS(delta)
sgn = delta/AbsDelta
jFrequency = iFrequency + sgn
jFrequency = jFrequency > 0
jFrequency = jFrequency < nOmegaMax
;
; and now evaluate response at function at linear frequency
;
Values = *(Response_k_w.values)
newValues = (1.-Absdelta)*Values(*,iFrequency) + Absdelta*Values(*,jFrequency)
tauValues = *(TauCorr.values)
newTauValues = (1.-Absdelta)*tauValues(*,iFrequency) + Absdelta*tauValues(*,jFrequency)
;
; Create Response_k object
;
Response_k = response_k_w -> slice(omega=0)
PTR_FREE, Response_k.values
response_k.values = PTR_NEW(newValues)
response_k.title = "!12Q!X!DQL!N"
response_k.mnemonic = "Q_ql"
response_k.Indices = PTR_NEW(["!4x!X=!4x!X!Dlin!N", "*"])
;
; Create tauCorr_k object
;
tauCorr_k = TauCorr -> Slice(omega=0)
PTR_FREE, tauCorr_k.values
tauCorr_k.values = PTR_NEW(newTauValues)
tauCorr_k.title = "!4s!Dv!X!N"
tauCorr_k.mnemonic = "TauCorr_Chi"
tauCorr_k.Indices = PTR_NEW(["!4x!X=!4x!X!Dlin!N", "*"])
;
; Clean up
;
IF(N_ELEMENTS(deBug) EQ 0) THEN BEGIN
	response_k_w -> Trash
	tauCorr -> Trash
ENDIF
;
; add Response_k, tauCorr_k objects to output structure and return
;
result = CREATE_STRUCT(result, "Q_ql", Response_k, "Tau_Chi", TauCorr_k)
RETURN, result
END ; ****** Function GKVs2D::Frequency ****** ;


Function GKVs3D::Frequency, KEObj, tRange=tRange_in, deBug=deBug
;
; Acts on amplitude (we assume phi) 
; in a k-space (vs. time) representation.
; Computes time-dependent frequency
; by taking the logarimitic derivative.
; This scheme only works if data is from
; a linear run.
; 
;
; Argument (optional: The k-space (vs. time) representation of the 
; kinetic energy. If supplied, routine will assume that
; 'self' is the poential and routine will compute the
; quasi-linear response relating |phi(k)|^2 to the
; heat flux.
;
;  Keywords:
;
;	tRange	The time interval over which
;		the logarithmic dirivative should
;		be averaged to obtain an estimate
;		of the frequency.  Defaults to the 
;		last half of the data-set (optional).
;
; Written by W.M. Nevins
;   7/18/07
; Modified by W.M. Nevins
;   7/27/07
; to compute QL response functions.
;
t=*(self.grid3.values)
values = *self.values
nx = N_ELEMENTS(*self.grid1.values)
ny = N_ELEMENTS(*self.grid2.values)
One = MAKE_ARRAY(ny, VALUE=1.)
nt = N_ELEMENTS(t)
dt = One#(t[1:nt-1] - t[0:nt-2])
dVdt = MAKE_ARRAY(nx,ny,nt-1, /COMPLEX)
Vbar = MAKE_ARRAY(nx,ny,nt-1, /COMPLEX)
FOR i=0,nx-1 DO dVdt[i,*,*] = (values[i,*,1:nt-1] - values[i,*,0:nt-2])/dt
FOR i=0,nx-1 DO Vbar[i,*,*] = 0.5*(values[i,*,1:nt-1] + values[i,*,0:nt-2])
omega = dVdt/(Vbar)
gamma = FLOAT(omega)
omega = -IMAGINARY(omega)
tbar = 0.5*(t[1:nt-1] + t[0:nt-2])
tUnits = self.grid3.units
wUnits = '1/(' + tUnits + ')'
omegaObj = self -> MakeCopy(/NoValues, /NoErrorBars)
gammaObj = self -> MakeCopy(/NoValues, /ErrorBars)
omegaObj.title = "!4x!X"
gammaObJ.title = "!4c!X"
omegaObj.mnemonic = "omega"
gammaObj.mnemonic = "gamma"
omegaObj.units = wUnits
gammaObj.units = wUnits
tmin = MIN(tbar)
tmax = MAX(tbar)
PTR_FREE, omegaObj.grid3.values
omegaObj.grid3.values = PTR_NEW(tbar)
omegaObj.grid3.irange = [0,nt-2]
omegaObj.grid3.range = [tmin, tmax]
PTR_FREE, gammaObj.grid3.values
gammaObj.grid3.values = PTR_NEW(tbar)
gammaObj.grid3.irange = [0,nt-2]
gammaObj.grid3.range = [tmin, tmax]
wmin = MIN(omega)
wmax = max(omega)
;
; Coherce omega, gamma =0 at k_y=0.
;
kyvalues = *(OmegaObj.Grid2.values)
epsilon = MIN(kyValues^2, iky0)
omega[*,iky0,*]=0.
gamma[*,iky0,*]=0.
;
omegaObj.values = PTR_NEW(omega)
omegaObj.vrange = [wmin, wmax]
gmin = MIN(gamma)
gmax = MAX(gamma)
gammaObj.values = PTR_NEW(gamma)
gammaObj.vrange = [gmin, gmax]
;s
; Compute "average" frequency, growth rate
; from average over last half of time interval
;
t0   = t[0]
tmax = t[nt-1]
tmin = 0.5*(t0+tmax)
tRange = [tmin, tmax]
IF(N_ELEMENTS(tRange_in) EQ 2) THEN tRange = tRange_in
OmegaAvg = omegaObj -> Avg(t=tRange)
GAmmaAvg = gammaObj -> Avg(t=tRange)
omegaObj.grid3.iRange = [0,nt-2]
gammaObj.grid3.iRange = [0,nt-2]

result = {	Name	:	"Frequency",	$	
		omegaAvg:	omegaAvg,	$
		gammaAvg:	gammaAvg,	$
		omegaObj:	omegaObj,	$
		gammaObj:	gammaObj	}
N_Args = N_PARAMS()
IF(N_Args EQ 0) THEN RETURN, result
; 
; User supplied KeObj as second argument,
; so compute QL response function for
; heat flux.
;
; We will only use same part of time-
; record used to compute average frequency
;
phi_k = self -> MakeCopy()
phi_k -> SignalWindow, t=trange
phi_k -> Restrict
KE_k = KEobj -> MakeCopy()
KE_k -> SignalWindow, t=trange
KE_k -> Restrict
;
; Fourier transform both Phi_k and
; KE_k in time
;
Phi_k_w = phi_k -> FFT('t')
KE_k_w  = KE_k  -> FFT('t')
;
; Form |Phi|^2 and 
; cross-spectrum between
; dPhi/dy and kE
;
PhiSq_k_w = Phi_k_w -> AbsSq()
KEStar_k_w = KE_k_w -> EXECUTE("CONJ")
KE_k_w -> TRASH
MinusEye = COMPLEX(0.,-1.)
kyString = Phi_k_w.grid2.mnemonic
temp   = phi_k_w -> times(kyString)
dPHidy = temp -> Times(MinusEye)
dPhidySq = dPhidy -> AbsSq()
temp -> Trash
CrossSpect = dPhidy -> Times(KEStar_k_w)
;
; Form QL response function vs. (k, omega)
;
temp = CrossSpect -> Over(PhiSq_k_w)
Response_k_w = temp -> Execute("Float")
;
; clean up
;
IF(N_ELEMENTS(debug) EQ 0) THEN BEGIN
	CrossSpect -> Trash
	PhiSq_k_w  -> Trash
	dPhidy     -> Trash
	temp       -> Trash
ENDIF
; 
; Now evaluate Response_k_w at linear
; frequency. Begin by finding index
; corresponding to omega=0
;
Omega  = *(Response_k_w.grid3.values)
epsilon = MIN(Omega^2, iOmega_0)
;
; Now compute index adjecent to
; linear frequency.
; 
dOmega = Omega[1] - Omega[0]
frequency = -(*(OmegaAvg.values))
nFrequency = frequency/dOmega
iFrequency = FIX(nFrequency)
iOmega = iFrequency + iOmega_0
nOmegaMax = N_ELEMENTS(Omega) - 1
iFrequency = iFrequency > 0
iFrequency = iFrequency < nOmegaMax
;
; Prepare for linear interpolation
;
delta = nFrequency - iFrequency
AbsDelta = ABS(delta)
sgn = delta/AbsDelta
jFrequency = iFrequency + sgn
jFrequency = jFrequency > 0
jFrequency = jFrequency < nOmegaMax
;
; and now evaluate response at function at linear frequency
;
Values = *(Response_k_w.values)
newValues = (1.-Absdelta)*Values(*,*,iFrequency) + Absdelta*Values(*,*,jFrequency)
;
; Create Response_k object
;
Response_k = response_k_w -> slice(omega=0)
IF(N_ELEMENTS(deBug) EQ 0) THEN PTR_FREE, Response_k.values
response_k.values = PTR_NEW(newValues)
response_k.title = "!12Q!X!DQL!N"
response_k.mnemonic = "Q_ql"
response_k.Indices = PTR_NEW(["!4x!X=!4x!X!Dlin!N", "*", "*"])
;
; Clean up
;
IF(N_ELEMENTS(deBug) EQ 0) THEN response_k_w -> Trash
;
; add Response_k object to output structure and return
;
result = CREATE_STRUCT(result, "Q_ql", Response_k)
RETURN, result
END ; ****** Function GKVs3D::Frequency ****** ;
