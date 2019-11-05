FUNCTION GKVs2D::GYRO_4WaveAvg, _Extra=Extra
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
		dk_r	:  dk			,	$
		x	:  Phi			,	$
		kr	:  phi_k		,	$
		kr_Sq	:  phi_kSq		,	$
	      Intensity	:  intensityObj		,	$
		avg_k_r	:  avg_kObj		,	$
	      delta_k_r	:  delta_kObj		,	$
		max_k_r	:  kmaxObj		,	$
	      maxIndex	:  maxIndex			}
;
; Cleanup ?
;
maxstr.slice -> trash
;
; and we're done
;
RETURN, output
END  ;  ****** GYRO_4WaveAvg ****** ;


FUNCTION GYRO_4Wave, 	phiStrIn, _Extra=Extra
;
; runtime keywords (which are searched for in 'Extra')
;
;			CaseID=CaseID,				$
;			rhoStar=rhoStar, gammaMax=gammaMax,	$
;			q=q, AspectRatio=AspectRatio
;
; Purpose:
;
;	This function acts on GYRO (flux tube, for the moment)
;	data, and performs analysis motivatged by the 4-wave
;	model [see Chen, Lin, and White, Phys. Plasmas <7>,
;	3129 (Aug. 2000)].j
;
; 
; Argument:
;
;	The argument of this function is the structure
;	output by the GKV function 'delta'.  If this argument
;	is not supplied, GYRO_4Wave will retrieve the required
;	data (using NETCDF_DATA) and produce same. (Optional)
;		 
;
; Keywords:
;
;	CaseID		A string variable containing a legal IDL 
;			file-name to be used to identify this 
;			case. (Optional)
;
;	rhoStar		The value of rho/a for this run.  If one
;			is not supplied, it will be assumed that
;			the data is already in gyrokinetic units,
;			and no scaling of either amplitude or
;			spatial variables will be performed.
;			(Optional)
;
;	gammaMax	The maximum linear growth rate.  If supplied,
;			gammaMax will be used to scale the time varaible.
;			(Optional)
;
;	q		The safety factor.  This is be used to scale
;			the 'zeta' axis to gyro-kinetic units (that is,
;			distance in the bi-normal). Defaults to 1.4
;			(CYCLONE base case value).
;			(Optional)
;
;	AspectRatio	The aspect ratio.  This is required 
;			(in in addition to q) to correctly
;			scale the 'zeta' axis to distance in the
;			bi-normal direction.  Defaults to 2.7775
;			(CYCLONE base case value). 
;			(Optional) 
;
;  Written by W.M. Nevins
;	5/8/04
;
; Gete CaseID
;
CaseID="GKV_4Wave"
result = GetKeyWord('Case_ID', Extra)
IF(TypeOf(result) EQ 7)    THEN BEGIN
IF(result NE 'undefined')  THEN BEGIN
	CaseID = result
	MESSAGE, CaseID, /INFORMATIONAL
ENDIF
ENDIF
;
; Save Inputs
;
Input = GKV_CopyStructure(Extra)
;
; Check for phistr.  If not there, prompt user for data
; (which is done by NETCDF_DATA) 
;
nArgs = N_PARAMS()
CASE nArgs OF
	0  :  BEGIN
	   	; User failed to supply a phistr, so 
	   	; get data using NETCDF_DATA
		dataStr = NETCDF_DATA()
		phi_n = dataStr.potential -> all_k('n')
		gkvdelete, dataStr
		phi = phi_n -> FFT('n', /INVERSE)
		phi_n -> Trash
		phistr = phi -> delta(axis='zeta')
		phi -> Trash
	     END
	1  :    PhiStr = GKV_CopyStructure(phistrIn)
ENDCASE
;
; Add 'nArgs' tag ot input structure
;
CASE TypeOf(Input) OF
	8	:  Input = CREATE_STRUCT('nArgs', nArgs, Input)
	Else	:  Input = { nArgs : nArgs }
ENDCASE
;
; Check for rhoStar and, if supplied by user,
; rescale radial and zeta axes.
;
result = GetKeyWord('rhoStar', Extra)
IF(TypeOf(result) NE 7) THEN BEGIN
;
; rhoStar was supplied by the user, so we need to re-scale dependent
; variable AND spatial co-ordinates.
;
	rhoStar = result
	IF(rhoStar GT 1) THEN rhoStar = 1./rhoStar
	dPhi = phiStr.delta -> OVER(rhoStar)
	phiAvg = phistr.avg -> OVER(rhostar)
	dPhi   -> Set, units='(!4q!X!Ds!N/a)(T/e)'
	phiAvg -> Set, units='(!4q!X!Ds!N/a)(T/e)'
	dPhi   -> ScaleAxis, 1, const=1./rhoStar, units='!4q!X!Ds!N'
	phiAvg -> ScaleAxis, 1, const=1./rhoStar, units='!4q!X!Ds!N'
;
; Check for q and AspectRatio
;
	q = 1.4
	result = GetKeyWord('q_safety', Extra)
	IF(TypeOf(result) NE 7) THEN q=result
	aspectRatio = 2.7775
	result = GetKeyWord('AspectRatio', Extra)
	IF(TypeOf(result) NE 7) THEN aspectRatio=result
;
; Form scale_factor to get from toroidal angle to bi-normal
;	dperp/zeta = (R/a)*(a/rho)/SQRT{1 + [(R/r)*q(r)]^2}
;
	scaleFactor = ((AspectRatio+0.5)/rhoStar)/SQRT( 1 + (2.*AspectRatio*q)^2 )
	dPhi -> ScaleAxis, 'zeta', const=scaleFactor, title='!9x!X!N', mnemonic='y', units='!4q!X!Ds!N'	
ENDIF ELSE BEGIN
	dPhi = phistr.delta -> MakeCopy()
	phiAvg = phistr.avg -> MakeCopy()
ENDELSE
GKVdelete, phistr
;
; Make sure that the radial axis has a uniform grid spacing
;
dPhi   -> ScaleAxis, 1, /Uniform
phiAvg -> ScaleAxis, 1, /Uniform 
;
; Now check for gammaMax. If it's there, scale time axis.
;
result = GetKeyWord('gammaMax', Extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	gammaMax = result
	dPhi   -> ScaleAxis, 3, const=gammaMax, title='t', mnemonic='t', units='1/!4c!X!Dmax!N'
	phiAvg -> ScaleAxis, 2, const=gammaMax, title='t', mnemonic='t', units='1/!4c!X!Dmax!N'
ENDIF
;
; Form intensity averaged over y vs. k_r and t.
; 
dPhi_k   = dPhi   -> FFT(1)
dPhi_kSq = dPhi_k -> AbsSq()
dPhi_k -> Trash
dPhi_krSq = dPhi_kSq -> Avg(2)
dPhi_kSq -> Trash
;
; get t-values, and number of time steps
;
dPhi_krSq -> Get, axis=2, gridValues = tPtr
nt = N_ELEMENTS(*tPtr)
ones = MAKE_ARRAY(nt, /FLOAT, VALUE=1.)
;
; get k_r-values, form <k_r> and standard deviation
;
dPhiSq = dPhi_krSq -> MakeCopy()
dPhiSq -> get, axis=1, gridValues=k_rPtr, gridTitle=kTitle, gridmnemonic=kmnemonic
k_r = *k_rPtr
dk_r = k_r[1] - k_r[0]
nks = N_ELEMENTS(k_r)
kmax = k_r[nks-1]
dPhiSq -> SignalWindow, k_r=[0,kmax]
dPhiSq -> Restrict
dPhiSq -> get, axis=1, gridValues=k_rPtr, gridTitle=kTitle, gridmnemonic=kmnemonic
k_r = *k_rPtr
nks = N_ELEMENTS(k_r)
kk = k_r#ones
dPhiSq -> Get, values=valuePtr
phiSq = *valuePtr
intensity = TOTAL(phiSQ,1)
avg_k = TOTAL(phisq*kk,1)/intensity
avg_kSq = TOTAL(phiSQ*kk*kk,1)/intensity
delta_k = SQRT(avg_ksq - avg_k^2)
;
; Form <k_r> object
;
avg_k_rObj = dPhi_krsq -> slice(axis=1, index=0)
;
; reset object fields
;
avg_k_rObj -> set, title='!12<!X' + ktitle + '!12>!X
avg_k_rObj -> set, mnemonic = 'Avg_' + kmnemonic
avg_k_rObj -> set, units='1/!4q!X!Ds!N'
avg_k_rObj -> get, indices=indices, values=valuePtr
PTR_FREE, indices
avg_k_rObj -> set, indices=PTR_NEW(['*'])
PTR_FREE, valuePtr
avg_k_rObj -> set, values = PTR_NEW(avg_k)
avg_k_rObj -> set, errorBars = PTR_NEW(delta_k)
;
; form  delta k object
;
delta_k_rObj = avg_k_rObj -> MakeCopy(/NoValues, /NoErrorBars)
;
; reset object fields
;
delta_k_rObj -> set, title='!14d!X' + ktitle
delta_k_rObj -> set, mnemonic = 'delta_' + kmnemonic
delta_k_rObj -> set, values = PTR_NEW(delta_k)
;
; Form intensity averaged over r vs. k_y and t.
;
dPhi_k   = dPhi   -> FFT(2)
dPhi_kSq = dPhi_k -> AbsSq()
dPhi_k -> Trash
dPhi_kySq = dPhi_kSq -> Avg(1)
dPHi_kSq -> Trash
;
; get k_y-values, form <k_y> and standard deviation
;
dPhiSq = dPhi_kySq -> MakeCopy()
dPhiSq -> get, axis=1, gridValues=k_yPtr, gridTitle=kTitle, gridmnemonic=kmnemonic
k_y = *k_yPtr
dk_y = k_y[1] - k_y[0]
nks = N_ELEMENTS(k_y)
kmax = k_y[nks-1]
dPhiSq -> SignalWindow, k_y=[0,kmax]
dPhiSq -> Restrict
dPhiSq -> get, axis=1, gridValues=k_yPtr, gridTitle=kTitle, gridmnemonic=kmnemonic
k_y = *k_yPtr
nks = N_ELEMENTS(k_y)
kk = k_y#ones
dPhiSq -> Get, values=valuePtr
phiSq = *valuePtr
intensity = TOTAL(phiSQ,1)
avg_k = TOTAL(phisq*kk,1)/intensity
avg_kSq = TOTAL(phiSQ*kk*kk,1)/intensity
delta_k = SQRT(avg_ksq - avg_k^2)
;
; Form <k_y> object
;
avg_k_yObj = dPhi_kysq -> slice(axis=1, index=0)
;
; reset object fields
;
avg_k_yObj -> set, title='!12<!X' + ktitle + '!12>!X
avg_k_yObj -> set, mnemonic = 'Avg_' + kmnemonic
avg_k_yObj -> set, units='1/!4q!X!Ds!N'
avg_k_yObj -> get, indices=indices, values=valuePtr
PTR_FREE, indices
avg_k_yObj -> set, indices=PTR_NEW(['*'])
PTR_FREE, valuePtr
avg_k_yObj -> set, values = PTR_NEW(avg_k)
avg_k_yObj -> set, errorBars = PTR_NEW(delta_k)
;
; form  delta k object
;
delta_k_yObj = avg_k_yObj -> MakeCopy(/NoValues, /NoErrorBars)
;
; reset object fields
;
delta_k_yObj -> set, title='!14d!X' + ktitle
delta_k_yObj -> set, mnemonic = 'delta_' + kmnemonic
delta_k_yObj -> set, values = PTR_NEW(delta_k)
;
; Form intensity vs. (r,t) and vs. t.
;
temp = dPhi -> AbsSq()
dPhiSq = temp -> Avg(2)
temp -> Trash
Intensity = dPhiSq -> Avg(1)
;
; Form dPhiStr
;
dPhiStr = {  Intensity	:	Intensity	,	$
		kr_Sq	:	dPhi_krSq	,	$
		avg_k_r	:	avg_k_rObj	,	$
	      delta_k_r :	delta_k_rObj	,	$
		ky_Sq	:	dPhi_kySq	,	$
		avg_k_y :	avg_k_yObj	,	$
	      delta_k_y :	delta_k_yObj		}
;
; Analyze avgPhi
;
PhiAvgStr = phiAvg -> GYRO_4WaveAvg()
;
; Form output stucture
;
output = {	Name	:	CaseId		,	$
		dk_r	:	dk_r		,	$
		dk_y	:	dk_y		,	$s
		Input	:	Input		,	$
		dPhi	:	dPhiStr		,	$
		AvgPhi	:	PhiAvgStr		}
;
; Clean-up Tasks 
;
dPhi -> Trash
;
RETURN, output
END  ;  ****** GYRO_4Wave ******  ;
