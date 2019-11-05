PRO GKV_DeTrend, signal, t, order=iorder
;
;  Purpose:
;
;	This proceedure 'detrends' signal
;
nt = N_ELEMENTS(t)
order = 1
IF(N_ELEMENTS(iOrder) EQ 1) THEN order = iOrder
rSignal = FLOAT(signal)
IF Query_Complex(signal) THEN $
	iSignal = COMPLEX(signal)

coeffs = POLY_FIT(t, rsignal, order)

IF Query_Complex(signal) THEN BEGIN
	iSignal = IMAGINARY(signal)
	iCoeffs = POLY_FIT(t, iSignal, order)
	coeffs = COMPLEX(coeffs, iCoeffs)
ENDIF

trend = FLTARR(nt)
tt = REPLICATE(1., nt)
FOR i=0,order DO BEGIN
	trend = trend + coeffs[i]*tt
	tt = tt*t
ENDFOR
signal = signal - trend
RETURN
END ; ****** GKV_DeTrend ****** ;

FUNCTION GKVs2D::DeTrend, order=iorder
;
; Purpose:
;
;	Just a dummy to stop people from using DeTrend on data objects
;	of dimensionality greater than 1.
;
MESSAGE, "Only implimented for 1-D objects", /INFORMATIONAL
RETURN, 0
END

FUNCTION GKVs1D::DeTrend, order=iorder
;
; Purpose:
;
;	This proceedure removes avg, trend, etc. through
;	order 'iorder'
;
; Input Keywords:
;
;	iorder	The order to which you wish to detrend 'self'
;
; Written by W.M. Nevins
;	6/17/00
;
result = self -> Makecopy()
result -> Restrict
t = *result.grid1.values
values = *result.values
GKV_DeTrend, values, t, order=iorder
vmin = GKVsd_MIN(values, MAX=vmax)
result.values = PTR_NEW(values)
result.vrange = [vmin,vmax]
RETURN, result
END ; ****** GKVs1D::DeTrend ****** ;


FUNCTION GKVs1D::SubSample, axisID, _EXTRA=extra
;
; Purpose:
;
; 		This function returns a GKVsd object of the 
;		same dimensionality as 'self' in which
;		the data of 'self' has been 'subSampled' onto 
;		a uniform grid of 'nSteps' on the 'axisID' axis.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis	If no argument is provided, then this keyword may be 
;			used to identify independent variable to be subsampled. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be subsampled, and to reset the signal window on this axis (before subsampling).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before subsampling on the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over subsampling is performed
;
;	   nSteps	The desired number of gridpoint for the selected axis
;			Defaults to no change in number of grid points over
;			'self's SignalWindow. (Optional).
;
;	   DownBy	Decrease number of grid points by a factor of 'Downby'
;			Defaults to no change in number of grid points over
;			'self's SignalWindow. (Optional).
;
;		 Dt	Desired size of the uniform sampling interval in the output grid.
;			Defaults to no change in number of grid points over
;			'self's SignalWindow. (Optional).
;
; Written by W.M. Nevins
;	5/13/01
;
;
; Find axis identifier
;
temp = self -> MakeCopy()
CASE N_PARAMS() OF
	0	:	axis = temp -> AxisIrange(        _Extra=extra)
	1	:	axis = temp -> AxisIrange(axisID, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'SubSample called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get selected grid structure
;
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = temp.Grid' + axisString
ok = EXECUTE(commandString)
irange = Grid.irange
nStepsIn = 1 + irange[1] - irange[0]
gridRange = Grid.range
interval = gridRange[1] - gridRange[0]
;
; Get number of steps in the output grid
;
nSteps = nStepsIn
result = GetKeyWord('nSteps', extra)		; Check for keyword 'nSteps'
IF( Query_Integer(result) + Query_Real(result) ) THEN nSteps = FIX(result)
result = GetKeyWord('DownBy', extra)		; Check for keyword 'DownBy' 
IF( Query_Integer(result) + Query_Real(result) ) THEN BEGIN
	IF( result GE 1 ) THEN nSteps = FIX(nSteps/result) 
	IF( result LT 1 ) THEN nSteps = FIX(nSteps*result)
ENDIF
result = GetKeyWord('Dt', extra)		; Checee for keyword 'Dt'
IF( Query_Integer(result) + Query_Real(result) ) THEN BEGIN
	nSteps = ( FIX(interval/result) > 0 ) + 1
ENDIF
nSteps = nSteps < nStepsIn
nSteps = nSteps > 2
;
; Find desired sampling interval
;
delta = interval/(nSteps - 1)
;
; Compute new grid values
;
gridValues = gridRange[0] + delta*FINDGEN(nSteps)
;
; Create output grid structure
;
PTR_FREE, grid.values
grid.values = PTR_NEW(gridValues)
grid.irange = [0, nSteps-1]
;
; load corrected grid structure into 'temp'
;
commandString = 'temp.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
;
; Now, interpolate 'self' onto temp's grid
;
result = self -> Interpolate(temp)
;
; and we're done
;
RETURN, result
END ; ****** GKVs1D::SubSample ****** ;



PRO GKVS1D::RMSNorm, RMSValue=rmsValue
;
; Purspose:
;
;		This proceedure normalizes 'self' such that the
;		RMS value of the signal within the 'SignalWindow'
;		is equal to rmsValue.
;
; Arguments:
;
;		None
;
; Keywords:
;		
;	RMSValue	Set this keyword equal to the desired rms value of signal
;			Defaults to 1.0. (Optional)
;
; Written by W.M. Nevins
;	5/13/01
;
rms = 1.0
IF(N_ELEMENTS(rmsValue) EQ 1) THEN rms = rmsValue
valuesIn = *(self -> GetValues() )
IF Query_Complex(valuesIN) THEN BEGIN
	avgValuesSq = TOTAL( valuesIn*CONJ(valuesIn) )/N_ELEMENTS(valuesIn)
ENDIF ELSE BEGIN
	avgValuesSq = TOTAL(valuesIn^2)/N_ELEMENTS(valuesIn)
ENDELSE
valuesOut = rms*(*self.values)/SQRT(avgValuesSq)
PTR_FREE, self.values
self.values = PTR_NEW(valuesOut)
vmin = GKVsd_MIN(valuesOut, Max=vmax)
self.vrange= [vmin, vmax]
RETURN
END ; ****** GKVS1D::RMSNorm ****** ;



FUNCTION GKVs1D::BiSpect, arg1, arg2, Idw=idataWindow, ilw = iLagWindow, order = iorder
;
; Purpose:
;
;	This function computes the bispectrum of
;	'self' with the arguments 'arg1' and arg2'
;
;  Arguments:
;
;	arg1	A GKVs1D time series.  Defaults to 
;		'self'.  (Optional)
;
;	arg2	A GKVs1D time series.  Defaults to 
;		'self'.  (Optional)
;
;  Keywords
;
;	idw	The (integer) length of the data window.
;		Defaults to 1/10 of the length of the data series.
;
;	ilw	The (integer) length of the lag window.
;		Defaults to 1/2 of the length of the data series.
;
;	order	Order at which to detrend (0 to just remove average, 
;		1 to remove average and trend).  Defaults to 1.
;		(Optional)
;
;
;  Written by W.M. Nevins
;	5/1/01
;
nargs = N_PARAMS()
;
; Make temporary copy of 'self'
;
temp = OBJARR(3)
temp[0] = self -> MakeCopy()
temp[0] -> Restrict
temp[0] -> ScaleAxis, axis=1, /uniform
;
; now make temporary copies of arguments (or, self, if argument is missing)
;
CASE nargs OF
	0	:	BEGIN
				temp[1] = temp[0] -> MakeCopy()
				temp[2] = temp[0] -> MakeCopy()
			END
	1	:	BEGIN
				temp[1] =    arg1 -> INTERPOLATE(temp[0])
				temp[2] = temp[0] -> MakeCopy()
			END
	2	:	BEGIN
				temp[1] =  arg1 -> INTERPOLATE(temp[0])
				temp[2] =  arg2 -> INTERPOLATE(temp[0])
			END
	ELSE	:	BEGIN
				MESSAGE, "Called with too many argument -- returning", /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
;
; Normalize input objects such that each has an RMS value of 1.0
;
;FOR i=0,2 DO temp[i] -> RMSNorm
;
; Get info about number of time steps, etc
;
t = *((temp[0]).grid1.values)
nt = N_ELEMENTS(t)
dt = t[1] - t[0]
;
; Get values
;
vPtr = PTRARR(3)
values = FLTARR(nt,3)
FOR i=0,2 DO vPtr[i] = temp[i] -> GetValues()
v0 = *(vPtr[0])
v1 = *(vPtr[1])
v2 = *(vPtr[2])
FOR i=0,2 DO PTR_FREE, vPtr[i]
;
; remove average and trend
;
order = 1
IF(N_ELEMENTS(iOrder) EQ 1) THEN order = iOrder
GKV_DeTrend, v0, t, order=order
GKV_DeTrend, v1, t, order=order
GKV_DeTrend, v2, t, order=order
;
; Pack values into arrays with power-of-2 length
;
idw = (nt/10) > 1
IF(N_ELEMENTS(idataWindow) EQ 1) THEN idw = idataWindow
windowFcn = Window_Fcn(nt, idw)
windowNorm = SQRT( nt/TOTAL(windowFcn*windowFcn) )
windowFcn = windowNorm*windowFcn
j=0L
REPEAT j=j+1 UNTIL (2L^j GE nt) 	; Find n_corrs such that 2^n_corrs > nt
n_corrs = 2L^j
n_corrs = 2*n_corrs
v0arr = COMPLEXARR(n_corrs)
v1arr = COMPLEXARR(n_corrs)
v2arr = COMPLEXARR(n_corrs)
v0arr[0:nt-1] = v0*windowFcn
v1arr[0:nt-1] = v1*windowFcn
v2arr[0:nt-1] = v2*windowFcn
;
; Take complex conjugate of data from 'self'
;
v0arr = CONJ(v0Arr)
;
; Compute fourier transforms
;
ft0 = FFT(v0arr, /INVERSE)
ft1 = FFT(v1arr, /INVERSE)
ft2 = FFT(v2arr, /INVERSE)
;
; Form 'bare' bispectra
;
biSpectra = (ft1#ft2) 
ft0Backwards = SHIFT(REVERSE(ft0), 1)/nt  ; NORMALIZATION???
FOR i=0,n_Corrs-1 DO BEGIN
	biSpectra[*,i] = biSpectra[*,i]*SHIFT(ft0Backwards, -i)	
ENDFOR
;
; and 'bare' spectral densities
;
spect0 = ft0*CONJ(ft0)						
spect1 = ft1*CONJ(ft1)
spect2 = ft2*CONJ(ft2)
;
; Now form biCorrelation function
;
biCorrs = FFT(biSpectra)
;biCorrs = biCorrs			; Correct normalization of biCorrelation function			
;
; and auto correlation functions
;
corr0 = FFT(spect0)/FLOAT(nt)		; Correct normalization of auto correlation functions					
corr1 = FFT(spect1)/FLOAT(nt)
corr2 = FFT(spect2)/FLOAT(nt)
;
; Apply lag window to biCorrelation function
;
ilw = nt/2
IF(N_ELEMENTS(ilagWindow) EQ 1) THEN ilw=iLagWindow
lagWindow = LagWindow_Fcn(n_corrs, ilw)
lagWindow = SHIFT(lagWindow, n_corrs/2)
fullLagWindow = lagwindow#lagWindow
biCorrs = biCorrs*fullLagWindow
;
; and to auto correlation functions
;
corr0 = corr0*lagWindow
corr1 = corr1*lagWindow
corr2 = corr2*lagWindow
;
; Now compute (frequency-averaged) biSpectrum
;
biSpectra = FFT(biCorrs, /INVERSE)					
;
; and (frequency-averaged) spectral densities
;
spect0 =   FFT(corr0, /INVERSE)						
spect0 = FLOAT(spect0)

spect1 =   FFT(corr1, /INVERSE)
spect1 = FLOAT(spect1) 

spect2 =   FFT(corr2, /INVERSE)
spect2 = FLOAT(spect2) 
;
; Compute biCoherence
;
biCoherence = biSpectra*CONJ(biSpectra)				
biCoherence = FLOAT(biCoherence) 
spect0Backwards = SHIFT(REVERSE(spect0), 1)
norm = spect1#spect2						
FOR i=0,n_Corrs-1 DO BEGIN
	norm[*,i] = norm[*,i]*SHIFT(spect0Backwards, -i) 
ENDFOR

epsilon = float(ilw)/n_Corrs

normMax = MAX(norm)
norm = norm > epsilon*normMax
biCoherence = biCoherence/norm
;biCoherence = biCoherence/n_Corrs	; is this the right normalization???
;
; Shift arrays so that omega=0 will be at center of frame
;
spect0 = SHIFT(spect0, n_corrs/2)
spect1 = SHIFT(spect1, n_corrs/2)
spect2 = SHIFT(spect2, n_corrs/2)
biSpectra = SHIFT(biSpectra, n_corrs/2, n_corrs/2)
biCoherence = SHIFT(biCoherence, n_Corrs/2, n_Corrs/2)
;
; Prepare output structure for biSpectrum
;
resultStr = {GKVs2D}
resultStr.mnemonic = 'biS_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic

titleIndices = STRARR(3)
FOR i=0,2 DO BEGIN
	tIndices = temp[i] -> IndexRemove(1)
	IF(N_ELEMENTS(tIndices) EQ 1) THEN BEGIN
		IF(tIndices EQ '') THEN BEGIN
			titleIndices[i] = ''
		ENDIF ELSE BEGIN
			titleIndices[i] = "[" + tIndices + ']'
		ENDELSE
		tIndices = STRJOIN(tIndices, ',')
		titleIndices[i] = "[" + tIndices + ']'
	ENDIF
ENDFOR
resultStr.title = 'biS{'+ (temp[0]).title + titleIndices[0] + '|'	$
			+ (temp[1]).title + titleIndices[1] + ',' 	$
			+ (temp[2]).title + titleIndices[2] + '}'
indices = REPLICATE('*', 2)
resultStr.indices = PTR_NEW(indices)
resultStr.units = ''
resultStr.values = PTR_NEW(biSpectra)
vmin = GKVsd_MIN(biSpectra, max=vmax)
resultStr.vrange = [vmin, vmax]
resultStr.CodeName = (temp[0]).CodeName
resultStr.CodePI   = (temp[0]).CodePI
resultStr.RunID    = (temp[0]).RunID
resultStr.FileID   = (temp[0]).FileID
grid1 = {Grid}
grid1.mnemonic = 'omega1'
grid1.title    = '!4x!X!D1!N'
grid1.units    = '1/(' + (temp[0]).grid1.units + ')'
gridValues = 2.*!PI/(dt*n_corrs)*(FINDGEN(n_corrs) - n_corrs/2 + 1)
grid1.values = PTR_NEW(gridValues)
grid1.boundary = 'periodic (open)'
grid1.uniform = 1B
xmin = GKVsd_MIN(gridValues, max=xmax)
grid1.range = [xmin,xmax]
grid1.irange = [0, n_corrs-1]
resultStr.Grid1 = grid1
grid2 = GKVsd_GridCopy(grid1)
grid2.mnemonic = 'omega2'
grid2.title = '!4x!X!D2!N'
resultStr.Grid2 = grid2
;
; Create biSpectra object
;
biSpectraObj = OBJ_NEW('GKVs2D', resultStr)
;
; Prepare output structure for biCoherence
;
resultStr.mnemonic = 'biC_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic 
resultStr.title = 'biC{'+ (temp[0]).title + titleIndices[0] + '|'	$
			+ (temp[1]).title + titleIndices[1] + ',' 	$
			+ (temp[2]).title + titleIndices[2] + '}'
resultStr.values = PTR_NEW(biCoherence)
vmin = GKVsd_MIN(biCoherence, max=vmax)
vmin = vmin > 0.
resultStr.vrange = [vmin, vmax]
;
; Create biCoherence object
;
biCoherenceObj = OBJ_NEW('GKVs2D', resultStr)
;
; Create output structure
;
result = {Name:'biSpectra', biSpectra:biSpectraObj, biCoherence:biCoherenceObj}
;
; Clean up
;
FOR i=0,2 DO temp[i] -> trash
;
; and we're done ...
;
RETURN, result
END ; ****** GKVs1D::BiSpect ******  ;
