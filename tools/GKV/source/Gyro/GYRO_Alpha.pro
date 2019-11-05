PRO GYRO_AlphaPatch1, AlphaStr, rRange=rRangeIn, tRange=tRangeIn, 	$
				dChidPhiSq = slopeIn

rRange=[0.42,0.58]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

tRange=[150.,900.]
IF(N_ELEMENTS(tRangeIn) EQ 2) THEN tRange=tRangeIn

phisq_r = AlphaStr.phisq_r -> makeCopy()
chi_r = alphaStr.Chi_r -> MakeCopy()
phisq_t = AlphaStr.phisq_t -> makeCopy()
chi_t = alphaStr.Chi_t -> MakeCopy()
phisq_t -> signalwindow, axis=1, range=rRange
phiSq_r -> signalwindow, t=tRange
Chi_t -> signalwindow, axis=1, range=rRange
Chi_r -> signalwindow, t=trange
PhiSq_r -> restrict
Chi_r -> restrict
PhiSq_t -> restrict
Chi_t -> restrict
Chi_vs_dPhiSq_r = Chi_r -> Squash(phiSq_r)
Chi_vs_dPhiSq_t = Chi_t -> Squash(phiSq_t)
Chi_vs_dPhiSq_t -> set, axis=1, Range=[0.,100.]
Chi_vs_dPhiSq_t -> set, vrange=[0.,10.]

;Chi_vs_dPhiSq_t -> Find_Slope, slope=dChidPhiSq
;dChidPhiSq = 0.056
;IF(N_ELEMENTS(slopeIN) EQ 1) THEN dChidPhiSq = slopeIn

;temp_r = phiSq_r -> times(dChidPhiSq)
;offset_r = chi_r -> Minus(temp_r)
;temp_r -> Trash

;temp_t = phiSq_t -> times(dChidPhiSq)
;offset_t = chi_t -> Minus(temp_t)
;temp_t -> Trash

phisq_t -> Trash
Chi_t   -> trash
phisq_r -> Trash
Chi_r   -> trash

AlphaStr = Create_Struct(AlphaStr, 	'Chi_vs_dPhiSq_r', Chi_vs_dPhiSq_r,	$
					'Chi_vs_dPhiSq_t', Chi_vs_dPhiSq_t	)	
;					'offset_r'       , offset_r,		$
;					'offset_t'       , offset_t		)
RETURN
END


PRO GYRO_AlphaPatch, AlphaStr, rRange=rRangeIn, tRange=tRangeIn

rRange=[0.42,0.58]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

tRange=[0.,900.]
IF(N_ELEMENTS(tRangeIn) EQ 2) THEN tRange=tRangeIn

phisq = AlphaStr.phisq -> makeCopy()
chi_i = alphaStr.Chi_i -> MakeCopy()
phisq -> signalwindow, axis=1, range=rRange
phiSq -> signalwindow, t=tRange
Chi_i -> signalwindow, axis=1, range=rRange
Chi_i -> signalwindow, t=trange
PhiSq -> restrict
Chi_i -> restrict
phiSq_r = phiSq -> avg('r')
PhiSq -> trash
Chi_r   = Chi_i -> avg('r')
Chi_i -> trash
AlphaStr = Create_Struct(AlphaStr, 	'Chi_r', Chi_r,		$
					'PhiSq_r', PhiSq_r)


RETURN
END


PRO GYRO_AlphaPlus, AlphaStr, rRange=rRangeIn, tRange=tRangeIn

rRange=[0.42,0.58]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

tRange=[0.,900.]
IF(N_ELEMENTS(tRangeIn) EQ 2) THEN tRange=tRangeIn

phisq = AlphaStr.phisq -> makeCopy()
chi_i = alphaStr.Chi_i -> MakeCopy()
phisq -> signalwindow, axis=1, range=rRange
phiSq -> signalwindow, t=tRange
Chi_i -> signalwindow, axis=1, range=rRange
Chi_i -> signalwindow, t=trange
PhiSq -> restrict
Chi_i -> restrict
phiSq_r = phiSq -> avg('r')
PhiSq -> trash
Chi_r   = Chi_i -> avg('r')
Chi_i -> trash
alpha_r = chi_r -> over(phiSq_r)
alpha_r -> get, errorBars=errorPtr
PTR_FREE, errorPtr
alpha_r -> set, errorBars=PTR_NEW()
alphaStats = alpha_r -> stats()
alphaAvg = alphaStats.avg
alphaError = alphaStats.avgpm
AlphaStr = Create_Struct(AlphaStr, 	'Alpha_r', alpha_r, 	$
					'alpha_Avg', alphaAvg,	$
					'Chi_r', Chi_r,		$
					'PhiSq_r', PhiSq_r,	$
					'alpha_Error', alphaError)


RETURN
END


FUNCTION GYRO_Alpha,	GYROStr, path=path, rhoStar=rhoStar, 	$
			rRange=rRangeIn, tRange=tRangeIn,	$
			Fluxtube=fluxtube
;
; 

CD, CURRENT=currentWorkingDirectory
IF(TypeOF(path) EQ 7) THEN CD, path

;rRange=[0.3,.7]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

tRange=[150.,900.]
IF(N_ELEMENTS(tRangeIN) EQ 2) THEN tRange=tRangeIn

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar

nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN GyroStr = NetCDF_Data()
phi_n = GyroStr.potential -> times(rhoStar)
phi_k = phi_n -> all_k('n')
phi_n -> trash
phi = phi_k -> FFT('n', /INVERSE)
phi_k -> trash
phistr = phi -> delta(axis='zeta')
phi -> trash
temp = phistr.delta -> times (phistr.delta)
;gkvDelete, phistr
temp -> set, title='!4du!x!U2!N', mnemonic='phiSq'
phiSq = temp -> avg('zeta')
phiSq -> GET, axis=2, range=tRangeInit
phiSq -> signalwindow, t=trange
;phiSq -> restrict
phiSq_t = phiSq -> avg_t()
PhiSq -> SignalWindow, t=tRangeInit
IF(N_ELEMENTS(FluxTube) EQ 1) THEN BEGIN
	Chi_i = GYROStr.diff_t_1 -> MakeCopy()
ENDIF ELSE BEGIN
	Chi_i = GYROStr.diff_t -> MakeCopy()
ENDELSE
Chi_i -> GET, axis=2, range=tRangeInit
Chi_i -> signalwindow, t=trange
;Chi_i -> restrict
Chi_t = Chi_i -> avg_t()
Chi_i -> SignalWindow, t=tRangeInit
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN BEGIN
	Chi_i -> signalwindow, axis=1, range=rRange
	Chi_i -> restrict
ENDIF
alpha = Chi_i -> over(phiSq)
alpha -> get, errorbars=errorPtr
IF(PTR_VALID(errorPtr)) THEN PTR_FREE, errorPtr
alpha -> set, errorbars=PTR_NEW(), mnemonic='alpha'
alpha -> GET, axis=2, range=tRangeInit
alpha -> signalwindow, t=trange
alpha_t = alpha -> avg_t()
alpha -> SignalWindow, t=tRangeInit
result = {	Name	:	'GyroAlpha',	$
		phiStr	:	phiStr,		$
		phiSq	:	phiSq,		$
		phiSq_t	:	phiSq_t,	$
		Chi_i	:	Chi_i,		$
		Chi_t	:	Chi_t,		$
		alpha	:	alpha,		$
		alpha_t	:	alpha_t		}
IF(nArgs EQ 0) THEN gkvDelete, GYROStr
RETURN, result
END
	
