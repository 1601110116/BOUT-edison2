FUNCTION GKVs2D::ZonalAnalysis, dPhi, _Extra=Extra
;
; Purpose:
;
;	This function acts on the toroidally averaged 
;	potential vs. (r,t) and returns information on 
;	the initial development of the Zonal flows for 
;	comparisons with the 4-wave model [see Chen, Lin, 
;	and White, Phys. Plasmas <7>, 3129 (Aug. 200)].
;
; Keywords:
;
;	rhoStar		Value of rho_s/a (or, alternatively a/rho_s -- 
;			either will work).  If rhoStar is present, it will
;			be used to re-scale both the radial axis and phi
;
;	gamma_max 	Maximum linear growth rate (in units of a/c_s).
;			If present, gamma_max is used to rescale the 
;			time axis and the Zonal flow intensity.
; 
;  Written by W.M. Nevins
;	5/7/04
;
nArgs= N_PARAMS()
dPhiSq = 'None'
IF(nArgs EQ 1) THEN BEGIN
	dPhi_kr = dPhi -> FFT('1')
	dPhi_kZeta = dPHi
	temp1  = temp  -> avg(axis=1)
	dPhiSq = temp1 -> avg(axis=1)
	dPhiSq -> set, title='!12<!4du!X!U2!N!12>!X', mnemonic='dPHiSq', units='(T/e)!U2!N'
	temp  -> trash
	temp1 -> trash
ENDIF

phi = self -> MakeCopy()
;
; Check if the time coordinate needs to be re-scaled
;
gamma_max=1.
result = GetKeyWord('gamma_max', Extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN BEGIN
; keyword "gamma_max" was set, so assume that
; the time-base of 'self' needs to be rescaled
	gamma_max = result
	phi -> ScaleAxis, 2, const=gamma_max, units='1/!4c!X!Dmax!N'
	IF(N_ELEMENTS(dPHiSq) NE 0) THEN $
		dPhiSq -> ScaleAxis, 1, const=gamma_max, units='1/!4c!X!Dmax!N'
ENDIF
;
; Check if the radial coordinate needs to be re-scaled
;
rhoStar=1.
result = GetKeyWord('rhoStar', Extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN BEGIN
; keyword "rhoStar" was set, so assume that
; the radial axis of 'self' needs to be rescaled
	rhoStar = result
	IF(rhoStar LT 1.) THEN rhoStar = 1./rhoStar
	phi -> ScaleAxis, 1, const=rhostar, units='!4q!Ds!N!X'
	phi -> get, title=title, mnemonic=mnemonic, units=units
	phi = phi -> times(rhoStar)
	phi -> set, title=title, mnemonic=mnemonic, units= '(!4q!X/a)' + units
	IF(N_ELEMENTS(dPhiSq) NE 0) THEN BEGIN
		dPhiSq -> get, title=title, mnemonic=mnemonic, units=units
		dPHiSq = dPHiSq -> times(rhoStar*rhoStar)
		dPHiSq -> set, title=title, mnemonic=mnemonic, units= '(!4q!X/a)!U2!N' + units
	ENDIF
ENDIF
phi -> ScaleAxis, 1, /uniform
;
; Transform to k-space represntation
;
phi_k = phi -> fft(1)
;
; Form absSq
;
phi_ksq = phi_k -> AbsSq()
kRange = phi_ksq.grid1.range
phi_ksq -> signalwindow, axis=1, range=[0,kRange[1]]
phi_ksq -> restrict
;
; get k-values and broadcast them in time
;
k=*phi_ksq.grid1.values
dk = k[1]-k[0]
t=*self.grid2.values
nt = N_ELEMENTS(t)
ones = MAKE_ARRAY(nt, /FLOAT, VALUE=1.)
kk = k#ones
;
; Form k-moments
;
phiSq = *phi_ksq.values
intensity = TOTAL(phiSQ,1)
avg_k = TOTAL(phisq*kk,1)/intensity
avg_kSq = TOTAL(phiSQ*kk*kk,1)/intensity
delta_k = SQRT(avg_ksq - avg_k^2)
;
; Find value of k at which Phi_kSq is maximum
;
maxStr=phi_kSq -> slice(axis=1, /AtMax, /maxLocation)
kmaxObj = maxstr.maxlocation
maxIndex = kMaxObj -> Over(dk)
maxIndex -> set, title='I!Dmax!N', mnemonic='maxIndex', units=''
;
; Form <k> object objects
;
avg_kobj = phi_ksq -> slice(axis=1, index=0)
;
; reset object fields
;
ktitle = phi_ksq.grid1.title
kmnemonic = phi_ksq.grid1.mnemonic
avg_kobj -> set, title='!12<!X' + ktitle + '!12>!X
avg_kobj -> set, mnemonic = 'Avg_' + kmnemonic
avg_kobj -> set, units='1/!4q!X!Ds!N'
PTR_FREE, avg_kobj.indices
avg_kobj -> set, indices=PTR_NEW(['*'])
PTR_FREE, avg_kobj.values
avg_kobj -> set, values = PTR_NEW(avg_k)
avg_kobj -> set, errorBars = PTR_NEW(delta_k)
;
; form  delta k object
;
delta_kObj = avg_kobj -> MakeCopy(/NoValues, /NoErrorBars)
;
; reset object fields
;
delta_kObj -> set, title='!14d!X' + ktitle
delta_kObj -> set, mnemonic = 'delta_' + kmnemonic
delta_kObj -> set, values = PTR_NEW(delta_k)
;
; form intensity object
;
intensityObj = avg_kobj -> MakeCopy(/NoValues, /NoErrorBars)
;
; reset object fields
;
intensityObj -> set, title=phi_kSq.title
intensityObj -> set, mnemonic=phi_kSq.mnemonic
intensityObj -> set, units = phi_kSq.units
intensityObj -> set, values = PTR_NEW(intensity)
;
; prepare output structure
;
output = {	Name 	:  'ZonalAnalysis'	,	$
	     gamma_max	:  gamma_max		,	$
	       rhoStar	:  rhoStar		,	$
		dk	:  dk			,	$
		PhiAvg	:  Phi			,	$
		phi_k	:  phi_k		,	$
		phiSq_k	:  phi_kSq		,	$
	      Intensity	:  intensityObj		,	$
		dPhiSq	:  dPhiSq		,	$
		k	:  avg_kObj		,	$
		delta_k	:  delta_kObj		,	$
		kmax	:  kmaxObj		,	$
	      maxIndex	:  maxIndex			}
;
; Cleanup ?
;
maxstr.slice -> trash
;
; and we're done
;
RETURN, output
END  ;  ****** ZonalAnalysis ****** ;
