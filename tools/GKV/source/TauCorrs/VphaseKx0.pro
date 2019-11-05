PRO VphaseKx0, pg3eqDataX, t=tInterval, dT=dT, dL=dL, NoQperp = noQperp, fraction=ff;; This proceedure accepts the structure produced when reading; an *.x.nc file with pg3eq as input.  It computes the (radial); phase velocity associated with the cross correlation between; the flux surface averaged potential (the "zonal flows") and the; the fulx surface averaged heat flux.;;; Arguments:;;		An anonymous structure containing GKVs1D objects;		xp_pot and xp_qtper. (REQUIRED).;; Keywords:;;	t	Set this equal to the time interval over which you wish ;		to compute the cross correlations between potential and;		heat flux. Defaults to current signal window of the potential.;		(Optional).;;	dT	Full width at half maximum (in time) of low pass filter applied;		to correlation function to remove "noise".  Should be set to;		something like the correlation time of the ITG turbulence.;;	dL	Full width at half maximum (in space) of low pass fileter applied;		to correlation function to remove "noise".  Should be set to ;		something liike the radial correlation length of the ITG turbulence.;;;    NoQperp	Set this key word (i.e., put "/NoQperp" on command line) to;		compute phase velocity from autocorrelations of potential.;		(Optional).;;   fraction	Measures phase velocity over range in correlation function greater;		than 'fraction' of its maximum value.  Defaults to 0.5.;		(optional).;;  Written by W.M. Nevins;  12/17/00;pg3eqDataX.xp_pot -> get, axis='t', SignalWindow=otIntervalIF(KEYWORD_SET(tInterval)) THEN BEGIN	pg3eqDataX.xp_pot   -> signalwindow, t=tInterval	IF( NOT KEYWORD_SET(NoQperp) ) THEN pg3eqDataX.xp_qtper -> signalwindow, t=tIntervalENDIF;; Form cross correlation function;IF( KEYWORD_SET(NoQperp) ) THEN BEGIN	xPotCorrs1 = pg3eqDataX.xp_pot -> xcorr()ENDIF ELSE BEGIN	xPotCorrs1 = pg3eqDataX.xp_pot -> xcorr(ref=pg3eqDataX.xp_pot)ENDELSE;; Take Real part of cross correlation function;xPotCorrs2 = xPotCorrs1 -> Execute("FLOAT");; If KEYWORD dT is set, then filter in time to remove GAM oscillations and other noise;IF KEYWORD_SET(Dt) THEN BEGIN	xPotCorrs2a = xPotCorrs2 -> Filter("tau", dT=dT)ENDIF ELSE BEGIN	xPotCorrs2a = xPotCorrs2 -> MakeCopy()ENDELSE;; If KEYWORD dL is set, then filter in space to remove 'noise';IF KEYWORD_SET(dL) THEN BEGIN	xPotCorrs3 = xPotCorrs2a -> Filter("x", dL=dL)ENDIF ELSE BEGIN	xPotCorrs3 = xPotCorrs2a -> MakeCopy()ENDELSE;; Slice at maximum;potCorrsMax = xPotCorrs3 -> slice(axis=1, /max, /maxLocation);; Find half-width at half maximum in time;halfWidth = ( potCorrsMax.slice -> FullWidth(fraction=ff) )/2.; ; Find time variation of half-width at half maximum in space;xSlice1 = xPotCorrs3 -> slice(tau=0.)xSlice2 = xPotcorrs3 -> Slice(tau=halfWidth)dxWidth1 = 0.5*( xSlice1 -> FullWidth(fraction=ff) )dxWidth2 = 0.5*( xSlice2 -> FullWidth(fraction=ff) ) dxWidth = (dxWidth2-dxwidth1) > 0.;; Compute phase velocity (and associate error bars);potCorrsMax.maxlocation -> SignalWindow, tau=[-halfWidth, halfwidth]vPhase = potCorrsMax.maxlocation -> slope(error=error)vPhaseError1 = errorvPhaseError2 = dxWidth/halfWidthvPhaseError = SQRT( vPhaseError1^2 + vPhaseError2^2 );; Print result;PRINT, "V_phase = ", vPhase, " +/- ", vPhaseError, vPhaseError1, vPhaseError2;; Final duties;xPotCorrs1 -> TrashxPotCorrs2 -> TrashxPotcorrs2a-> TrashxPotCorrs3 -> TrashpotCorrsMax.slice -> TrashpotCorrsMax.maxlocation -> TrashxSlice1 -> TrashxSlice2 -> TrashIF(KEYWORD_SET(tInterval)) THEN BEGIN	pg3eqDataX.xp_pot   -> set, axis='t', SignalWindow=otInterval	pg3eqDataX.xp_qtper -> set, axis='t', SignalWindow=otIntervalENDIF;RETURNEND ; ****** VphaseKx0 ****** ;