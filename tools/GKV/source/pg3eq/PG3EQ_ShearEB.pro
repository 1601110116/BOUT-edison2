FUNCTION GKVs2D::PG3EQ_ShearEB,	rhoStar=rho_star, L_T = L__T, rCorr=r_Corr, 	$
				tauCorr=tau_corr, tInterval=t_interval
;
;  This function expects "SELF" to be the flux-surface averaged
;  Potential (i.e., 'xp_pot' from the *.x.nc file) in units of 
;  (rho_s/L_T)*(T/e) vs. the dependent variables
;  t (in units of L_T/c_s) and radius (in units of rho_s).  
;  It returns a structure containing information about
;  the ExB shear.
;
;  Written by W.M. Nevins
;	10/25/03
;
L_T=1.0
IF(N_ELEMENTS(L__T) EQ 1) THEN L_T=L__T
;
; rescale averaged potential to units of (rho_s/a)*(T/e)
; and time to units of a/c_s
;
AvgPhi = self -> Over(L_T)
AvgPhi -> set, title='!12<!4u!12>!X', mnemonic='Avg_Phi', units='(!4q!X!Ds!N/a)(T/e)'
AvgPhi -> ScaleAxis, 't', const=L_T, units='a/c!Ds!N'
AvgPhi -> Set, axis=1, GridTitle='r', GridMnemonic='r'
AvgPhi -> ScaleAxis, 'r', /Uniform
AvgPhi -> ScaleAxis, 't', /Uniform
;
; Form the ExB shear from the (toroidally averaged) potential
;
temp = AvgPhi -> d2byd('r')
temp1 = temp -> Filter('t', Dt=tau_corr)
temp -> Trash
ShearEB = temp1 -> Filter('r', DL=  r_corr)
ShearEB -> set, title='!9d!XV!DExB!N/!9d!Xr', mnemonic='Shear_ExB', units='c!Ds!N/a'
;
; Compute RMS shear averaged over r and t
;
temp = ShearEB -> AbsSq()
temp1 = temp -> Avg('r')
Shear_Trace = temp1 -> Execute('SQRT')
temp1 -> Trash
Shear_Trace -> Set, title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X!S!Dr!R!U1/2!N', 	$
				mnemonic='avg_Shear_ExB', units='c!Ds!N/a'
tInterval = [200,1600]
IF(N_ELEMENTS(t_interval) EQ 2) THEN tInterval = t_interval
temp -> Signalwindow, t=tInterval
temp -> Restrict
AvgShearSq = temp -> Avg_t()
temp -> trash

temp = ShearEB -> MakeCopy()
temp -> SignalWindow, t=tInterval
temp -> Restrict
thisSpect = temp -> Xspect()
ShearSpect = thisSpect -> Int('omega')
ShearSpect -> Set, 	title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X', 	$
				mnemonic='ShearSpect', units='(c!Ds!N/a)!U2!N!4q!X!Ds!N'
thisSpect -> Trash
temp -> Trash
AvgShearSq -> Set, title='!12<!X(!9d!XV!DExB!N/!9d!Xr!X)!U2!N!12>!X!Dt!N', 	$
		 mnemonic='avg_Shear_ExB', units='(c!Ds!N/a)!U2!N'
rmsShear = AvgShearSq -> SQRT(units='c!Ds!N/a', mnemonic='rmsShear_ExB')

V_ExB = AvgPhi -> dbyd('r')
V_ExB -> set, title='V!DExB!N', mnemonic='V_ExB', units='(!4q!X!Ds!N/a)c!Ds!N'
temp = V_ExB -> MakeCopy()
temp -> SignalWindow, t=tInterval
temp -> Restrict
AvgV_ExB = temp -> avg_t()
temp -> trash

rhoStar = 1.0
IF(N_ELEMENTS(rho_star) EQ 1) THEN BEGIN
	rhoStar=rho_star
	IF(rhoStar GT 1.) THEN rhoStar=1./rhoStar
	AvgPhi     -> ScaleAxis, 'r', const=rhoStar, units='a'
	V_ExB      -> ScaleAxis, 'r', const=rhoStar, units='a'
	avgV_ExB   -> ScaleAxis, 'r', const=rhoStar, units='a'
	ShearEB    -> ScaleAxis, 'r', const=rhoStar, units='a'
	AvgShearSq -> ScaleAxis, 'r', const=rhoStar, units='a'
	rmsShear   -> ScaleAxis, 'r', const=rhoStar, units='a'
ENDIF


Result = { 	Name 		:	'PG3EQ_ShearEB',$
		PhiAvg		:	AvgPhi,	$
		V_ExB		:	V_ExB,		$
		avgV_ExB	:	avgV_ExB,	$
		Shear_ExB	:	ShearEB,	$
		ShearTrace	:	Shear_Trace,	$
		ShearSpect	:	ShearSpect,	$
		AvgShearSq	:	AvgShearSq,	$
		rmsShear_ExB	:	rmsShear	}
RETURN, result

END  ;  ********** GKVs2D::PG3EQ_ShearEB **********  ;
