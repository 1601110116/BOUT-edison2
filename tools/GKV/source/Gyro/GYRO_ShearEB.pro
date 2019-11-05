FUNCTION GKVs2D::GYRO_ShearEB,	rhoStar=rho_star, rCorr=r_Corr, 	$
				tauCorr=tau_corr, tInterval=t_interval,		$
				rIntervals=r_intervals, PhiNormed=PhiNormed
;
;  This function expects "SELF" to be the toroidally averaged
;  Potential in units of (T/e) vs. the dependent variables
;  t (in units of a/c_s) and radius (in units of a).  
;  It returns a structure containing information about
;  the ExB shear.
;
;  Written by W.M. Nevins
;	10/25/03

rhoStar=rho_star
IF(rhoStar GT 1.) THEN rhoStar=1./rhoStar
;
; Form the ExB shear from the (toroidally averaged) potential
;
temp = self -> d2byd('r')
temp1 = temp -> Filter('t', Dt=tau_corr)
temp -> Trash
temp = temp1 -> Filter('r', DL=rhoStar*r_corr)
ShearEB = temp -> Times(rhoStar)
temp  -> Trash
temp1 -> Trash
IF(KEYWORD_SET(PhiNormed)) THEN BEGIN
	temp = ShearEB
	ShearEB = temp -> Times(rhoStar)
	temp -> Trash
ENDIF
ShearEB -> set, title='!9d!XV!DExB!N/!9d!Xr', mnemonic='Shear_ExB', units='c!Ds!N/a'
;
; Compute RMS shear over designated intervals in r and t
;
rIntervals=[0.1, 0.26, 0.42, 0.58, 0.74, 0.9]
IF(N_ELEMENTS(r_intervals) GT 1) THEN rIntervals=r_intervals
nIntervals = N_ELEMENTS(rIntervals) - 1
Shear_Traces = OBJARR(nIntervals)
temp = ShearEB -> AbsSq()
FOR i=0, nIntervals-1 DO BEGIN
	temp1 = temp -> Avg(r=rIntervals[i:i+1])
	Shear_Traces[i] = temp1 -> Execute('SQRT')
	temp1 -> Trash
	Shear_Traces[i] -> Set, title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X!S!Dr!R!U1/2!N', 	$
				mnemonic='avg_Shear_ExB', units='c!Ds!N/a'
ENDFOR
temp -> SignalWindow, axis=1, range=[.1,.9]

tInterval = [200,900]
IF(N_ELEMENTS(t_interval) EQ 2) THEN tInterval = t_interval
temp -> Signalwindow, t=tInterval
temp -> Restrict
AvgShearSq = temp -> Avg_t()
temp -> trash

temp = ShearEB -> MakeCopy()
temp -> SignalWindow, t=tInterval
temp -> Restrict
ShearSpect = OBJARR(nIntervals)
FOR i=0,nIntervals-1 DO BEGIN
	thisInterval = temp -> window(r=rIntervals[i:i+1], /NORM)
	thisInterval -> ScaleAxis, 1, const=1./rhoStar, units='!4q!X!Ds!N'
	thisSpect = thisInterval -> Xspect()
	ShearSpect[i] = thisSpect -> Int('omega')
	ShearSpect[i] -> Set, 	title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X', 	$
				mnemonic='ShearSpect', units='(c!Ds!N/a)!U2!N!4q!X!Ds!N'
	thisInterval -> Trash
	thisSpect -> Trash
ENDFOR
ShearEB -> SignalWindow, axis=1, range=[.1,.9]

temp -> Trash
AvgShearSq -> Set, title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X!Dt!N', 	$
		 mnemonic='avg_Shear_ExB', units='(c!Ds!N/a)!U2!N'
rmsShear = AvgShearSq -> SQRT(units='c!Ds!N/a', mnemonic='rmsShear_ExB')

V_ExB = self -> dbyd('r')
V_ExB -> set, title='V!DExB!N', mnemonic='V_ExB', units='(!4q!X!Ds!N/a)c!Ds!N'
temp = V_ExB -> MakeCopy()
temp -> SignalWindow, t=tInterval
temp -> Restrict
AvgV_ExB = temp -> avg_t()
temp -> trash
PhiAvg = self -> MakeCopy()

Result = { 	Name 		:	'GYRO_ShearEB',	$
		t		:	tInterval,	$
		r		:	rIntervals,	$
		PhiAvg		:	PhiAvg,		$
		V_ExB		:	V_ExB,		$
		avgV_ExB	:	avgV_ExB,	$
		Shear_ExB	:	ShearEB,	$
		ShearTraces	:	Shear_Traces,	$
		ShearSpects	:	ShearSpect,	$
		AvgShearSq	:	AvgShearSq,	$
		rmsShear_ExB	:	rmsShear	}
RETURN, result

END  ;  ********** GKVs2D::GYRO_ShearEB **********  ;
