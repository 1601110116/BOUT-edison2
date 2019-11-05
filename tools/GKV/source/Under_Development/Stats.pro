FUNCTION GKVs1D::Stats, lw=lagWindow
;
; Purpose:
;
;		This function returns structure containing
;		the average, RMS value, and standard deviation of 'self'.
;
valuePtr = self -> GetValues()
values = FLOAT(*valuePtr)
info = SIZE(values)
nDims = info[0]
nPoints = info[nDims+2]
;
; old coding
;
;sumVals = TOTAL(values)
;sumValsSq = TOTAL(values^2)
;avgVal = sumVals/nPoints
;avgValSq = sumValsSq/nPoints
;rmsVal = SQRT(avgValSq)
;std = SQRT(avgValSQ - avgVal^2)
;
; New coding using IDL's 'MOMENT' function
;
moments = MOMENT(values, MDEV=mDev)
avgVal = moments[0]
mvar = moments[1]
avgValSq = avgVal^2
rmsVal = SQRT(avgValSq + mvar)
errorBars = FLTARR(nPoints)
IF(PTR_VALID(self.errorBars)) THEN errorBars = *( self -> GetErrors() )
AvgErrSq = TOTAL(errorBars^2)/nPoints
stdError = SQRT(AvgErrSq + mvar)
std = SQRT(mvar)
skewness = moments[2]
kurtosis = moments[3]
;
; This algorithm for the expected error in our estimate of the
; average comes from Jenkins & Watts, "Spectral Analysis and its Applications"
; Holden-Day, San Francisco, 1968 (see pg. 180).  (I borrowed this text from
; Tom Casper ...)
;
avgPM = 0.
avgmaxPM = 0.
tauCorrmax = 0.
tauCorrInt = 0.
tauCorrOne = 0.
IF(nDims EQ 1) THEN BEGIN
	; 
	; compute data window
	;
	idw = nPoints/10 > 1
	dw = Window_Fcn(nPoints, idw)
	dwNorm1 = TOTAL(dw)
	;
	; apply data window (no need to correct rms value ...)
	;
	values = values*dw
	;
	; subtract average
	;
	avg = TOTAL(values)
	values = values - avg*dw/dwNorm1
	;
	dt = (*self.Grid1.values)[1] - (*self.Grid1.values)[0]
	ncorrs = 2L
	WHILE ncorrs LT nPoints DO ncorrs = ncorrs*2L
	ncorrs = ncorrs*2L
	temp = COMPLEXARR(ncorrs)
	temp[0:nPoints-1] = COMPLEX(values, 0.)
	temp1 = FFT(temp, -1)
	temp = FFT(temp1*conj(temp1), 1)
	temp = FLOAT(temp)*FLOAT(nCorrs)/FLOAT(nPoints)	; 'temp' now contains the auto-variance
	var = temp[0]
	corr = temp/var				; 'corr' contains the correlation function
	;
	; Compute lag window
	;
	IF( (Query_Real(lagWindow) + Query_Integer(LagWindow)) EQ 0) THEN BEGIN
		hCorr = HILBERT(corr)
		envelope = SQRT( corr*corr + hCorr*hCorr )
		delta = (envelope[0:(nCorrs/2)] - 0.5)^2
		eps = MIN(delta, iHalf)		; point of this call is to get 'iHalf'
		ilw = LONG(SQRT(2L*iHalf*nPoints))
	ENDIF
	IF Query_Real(lagWindow) THEN BEGIN
		ilw=LONG(lagWindow/dt)
	ENDIF
	IF Query_Integer(lagWindow) THEN BEGIN
		ilw=lagWindow
	ENDIF
	ilw = ilw < (nCorrs-1)/2
	ilw = 2L*LONG(ilw/2)
	corr = SHIFT(corr, ilw)
	corr = corr[0:2*ilw]*LagWindow_Fcn(2*ilw+1, ilw)
	corr = SHIFT(corr, -ilw)
	int = TOTAL(corr, /CUMULATIVE)
	tauInt = int[2*ilw]
	dAvgSq = (tauInt/nPoints)*var
	avgPM = SQRT(dAvgSq)
	tauCorrInt = tauInt*dt
	tauMax = 2.*MAX(int[0:ilw], iMax)
	dAvgMaxSq = (tauMax/nPoints)*var
	avgMaxPM = SQRT(dAvgMaxSq)
	tauCorrMax = tauMax*dt
	tauCorrOne = (iMax+0.5)*dt
	taulag = ilw*dt
	tauHalf = iHalf*dt
	output={avg:avgVal, avgPM:avgPM, avgMaxPM:avgMaxPM, tauCorrInt:tauCorrInt,  tauHalf:tauHalf,	$
		tauLag:tauLag, tauCorrmax:tauCorrMax, tauCorr1:tauCorrOne, rms:rmsVal, std:std,    	$
	 	stdError:stdError, var:mvar, skewness:skewness, kurtosis:kurtosis, meanAbsDev:mDev}
	 RETURN, output
ENDIF


output={avg:avgVal, rms:rmsVal, std:std, stdError:stdError,  	$
	 var:mvar, skewness:skewness, kurtosis:kurtosis, meanAbsDev:mDev}
RETURN, output
END ; ****** GKVs1D::Stats ****** ;

