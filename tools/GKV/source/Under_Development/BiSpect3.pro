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



FUNCTION GKVs1D::BiSpect, arg1, arg2, 	dw=ddw, Idw=idataWindow, 	$
					lw=llw, ilw=iLagWindow, 	$
					order=iorder, epsilon=eps	
;
; Purpose:
;
;	This function a bispectral analysis of the data in
;	'self' (a GKVs1D object) with the data in arguments 
;	'arg1' and arg2'.  It returns a structure containing
;	the biSpectrum; the biCoherence; the biCorrelation function; 
;	the spectral densities of 'self', 'arg1', and 'arg2'; 
;	the correlation functions of 'self', 'arg1', and 'arg2';
;	and the length of the data and lag windows.
;
;	The algorithm used here is a REAL memory hog!  It requires
;	several arrays of (n_corrs x n_corrs), where n_corrs is a 
;	power of 2 which is MORE than twice the length of the data 
;	string.  I strongly recommend that you sub-sample the data 
;	(see 'GKVs1D::SubSample' above) to minimize the number of
;	data points.  
;
;	If after sub-sampling the data string is longer than about 
;	2000 points, you should consider using BiSpectL instead.
;	BispectL performs similar analysis using an algorithm which
;	does not require such large arrays.
;
;
;  Arguments:
;
;	arg1	A GKVs1D time series.  Defaults to 
;		'self'.  (Optional)
;
;	arg2	A GKVs1D time series.  Defaults to 
;		'self'.  (Optional)
;
;
;  Keywords:
;
;	dw	The length of the data window in proper time units.
;		(Optional).
;
;	idw	The (integer) length of the data window in time steps.
;		If 'idw' is set, its value over rides value given with 'dw'.
;		Defaults to 1/100 of the length of the data series.
;		(Optional).
;
;	lw	the length of the lagwindow in proper time units.
;		(Optional).
;
;	ilw	The (integer) length of the lag window.
;		if 'ilw' is set, its value over rides the value given with 'lw'. 
;		Defaults to 1/100 of the length of the data series.
;
;	order	Order at which to detrend (0 to just remove average, 
;		1 to remove average and trend).  Defaults to 1.
;		(Optional).
;
; 	epsilon	Fractional decrement allowed in norm used to compute
;		the bicoherence.  Such a limit is required because
;		alternative leads to spurious divides by small numbers
;		resulting from zeros in the fourier transform of the 
;		lag window.  If uncorrected, this would result in 
;		unwanted, large-amplitude noise in the estimate of the
;		bicoherence.  Defaults to 0.0085, the ratio of the first
;		subsidary maximum in the FT of a Hanning window (used
;		here as the lag window) to the maximum value.  
;		(Optional... and only experts should mess with this!).
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
idw = FIX(nt/100) > 1
IF(N_ELEMENTS(ddw) EQ 1) THEN idw = FIX(ddw/dt) > 1
IF(N_ELEMENTS(idataWindow) EQ 1) THEN idw = FIX(idataWindow)
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
ilw = FIX(nt/10) > 1 
IF(N_ELEMENTS(llw) GE 1) THEN ilw = FIX(llw/dt) > 1
IF(N_ELEMENTS(ilagWindow) GE 1) THEN ilw=FIX(iLagWindow)
IF(N_ELEMENTS(ilw) EQ 1) THEN ilw = [ilw, ilw, ilw]
IF(N_ELEMENTS(ilw) EQ 2) THEN ilw = [ilw[0], ilw[0], ilw[1]]
lagWindow0 = LagWindow_Fcn(n_corrs, ilw[0])
lagWindow0 = SHIFT(lagWindow0, n_corrs/2)
lagWindow1 = LagWindow_Fcn(n_corrs, ilw[1])
lagWindow1 = SHIFT(lagWindow1, n_corrs/2)
lagWindow2 = LagWindow_Fcn(n_corrs, ilw[2])
lagWindow2 = SHIFT(lagWindow2, n_corrs/2)
fullLagWindow = lagwindow1#lagWindow2
biCorrs = biCorrs*fullLagWindow
;
; and to auto correlation functions
;
corr0 = corr0*lagWindow0
corr1 = corr1*lagWindow1
corr2 = corr2*lagWindow2
;
; Now compute (frequency-averaged) biSpectrum
;
biSpectra = FFT(biCorrs, /INVERSE)					
;
; and (frequency-averaged) spectral densities
;
sspect0 =   FFT(corr0, /INVERSE)						
sspect0 = FLOAT(sspect0) > 0.

sspect1 =   FFT(corr1, /INVERSE)
sspect1 = FLOAT(sspect1) > 0.

sspect2 =   FFT(corr2, /INVERSE)
sspect2 = FLOAT(sspect2) > 0.
;
; Compute biCoherence
;
biCoherence = ABS(biSpectra)				
;
; Compute norm for biCoherence
;
spect0 = FLOAT(spect0)
spect1 = FLOAT(spect1)
spect2 = FLOAT(spect2)
 
spect0Backwards = SHIFT(REVERSE(spect0), 1)
norm = spect1#spect2						
FOR i=0,n_Corrs-1 DO BEGIN
	norm[*,i] = norm[*,i]*SHIFT(spect0Backwards, -i) 
ENDFOR

norm = SQRT(norm)/nt
;
; Norm is a 'proxy' for the biSpectra.  We now transform it to 
; tau1, tau1 space (where it is a proxy for the biCorrelation function)
; and apply the lag window (this insures norm has the SAME frequency
; filtering as the biSpectra).
;
norm = FFT(norm)
norm = norm*fullLagWindow
norm = FFT(norm, /INVERSE)
norm = FLOAT(norm)

epsilon = (0.0085)	; Choosen to remove structure in Bicoherence from sidebands in norm caused by
			; lagwindow.  Chosen as ratio of fLW(0) to value at first (subidiary maximum).
IF Query_REAL(eps) THEN epsilon=eps

normMax = MAX(norm)
norm = norm > epsilon*normMax
biCoherence = biCoherence/norm
;
; Renormalize biSpectra such that its integral over frequency is the power transfer.
;
biSpectra = biSpectra*( dt/(2.*!PI) )^2
;
; Renormalize sspect0, etc. such that their integral gives the variance of the 
; correspoding time series
;
sspect0 = sspect0*( dt/(2.*!PI) )
sspect1 = sspect1*( dt/(2.*!PI) )
sspect2 = sspect2*( dt/(2.*!PI) )
;
; Renormalize Corr0, etc. such that the value at zero lag is one.
;
corr0 = corr0/corr0[0]
corr1 = corr1/corr1[0]
corr2 = corr2/corr2[0]
;
; sspect0, corr0 contain spectral density, correlation function for the 
; complex conjugate of 'self'.  Switch this to those of 'self'
;
sspect0 = SHIFT(REVERSE(sspect0) ,1)
corr0   = CONJ(corr0)
;
; Shift arrays so that omega=0 will be at center of frame
;
sspect0 = SHIFT(sspect0, n_corrs/2)
sspect1 = SHIFT(sspect1, n_corrs/2)
sspect2 = SHIFT(sspect2, n_corrs/2)
corr0   = SHIFT(corr0  , n_Corrs/2)
corr1   = SHIFT(corr1  , n_Corrs/2)
corr2   = SHIFT(corr2  , n_Corrs/2)
biSpectra   = SHIFT(biSpectra,   n_corrs/2, n_corrs/2)
biCoherence = SHIFT(biCoherence, n_Corrs/2, n_Corrs/2)
;
; Prepare output structure for biSpectrum
;
resultStr = {GKVs2D}
resultStr.mnemonic = 'BiS_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic

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
		IF(tIndices NE '') THEN $
			titleIndices[i] = "[" + tIndices + ']'		
	ENDIF
ENDFOR
resultStr.title = '!18BiS[!X'+ (temp[0]).title + titleIndices[0] + '|'	$
			+ (temp[1]).title + titleIndices[1] + ',' 	$
			+ (temp[2]).title + titleIndices[2] + '!18]!X'
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
grid2 = GKVsd_GridCopy(Grid1)
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
resultStr.mnemonic = 'BiCoh_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic 
resultStr.title = '!18BiCoh[!X'+ (temp[0]).title + titleIndices[0] + '|'	$
			+ (temp[1]).title + titleIndices[1] + ',' 	$
			+ (temp[2]).title + titleIndices[2] + '!18]!X'
resultStr.indices = PTR_NEW(indices)
resultStr.values  = PTR_NEW(biCoherence)
;vmin = GKVsd_MIN(biCoherence, max=vmax)
;vmin = vmin > 0.
resultStr.vrange = [0.,1.]
newGrid1 = GKVsd_GridCopy(Grid1)
newGrid2 = GKVsd_GridCopy(Grid2)
resultStr.grid1 = newGrid1
resultStr.grid2 = newGrid2
;
; Create biCoherence object
;
biCoherenceObj = OBJ_NEW('GKVs2D', resultStr)
;
; prepare output stucture for biCorrelation function
;
resultStr.mnemonic = 'BiCorr_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic 
resultStr.title = '!18BiCorr[!X'+ (temp[0]).title + titleIndices[0] + '|'	$
			+ (temp[1]).title + titleIndices[1] + ',' 	$
			+ (temp[2]).title + titleIndices[2] + '!18]!X'
resultStr.indices = PTR_NEW(indices)
biCorrs = SHIFT(biCorrs, n_Corrs/2, n_Corrs/2)
keep = BiCorrs[(n_Corrs/2-ilw[1]):(n_Corrs/2+ilw[1]),(n_Corrs/2-ilw[2]):(n_Corrs/2+ilw[2])]
resultStr.values  = PTR_NEW(keep)
vmin = GKVsd_MIN(keep, max=vmax)
resultStr.vrange = [vmin, vmax]

tau1s = dt*(-ilw[1] + FINDGEN(2*ilw[1]+1))
tGrid1 = {Grid}
tGrid1.mnemonic = 'tau1'
tGrid1.title = '!4s!X!D1!N'
tGrid1.units = (temp[0]).grid1.units
tGrid1.values = PTR_NEW(tau1s)
tGrid1.boundary = 'periodic (closed)'
tGrid1.uniform = 1B
tGrid1.range = [-dt*ilw[1], dt*ilw[1]]
tGrid1.irange = [0, 2*ilw[1]]

tau2s = dt*(-ilw[2] + FINDGEN(2*ilw[2]+1))
tGrid2 = {Grid}
tGrid2.mnemonic = 'tau2'
tGrid2.title = '!4s!X!D2!N'
tGrid2.units = (temp[0]).grid1.units
tGrid2.values = PTR_NEW(tau2s)
tGrid2.boundary = 'periodic (closed)'
tGrid2.uniform = 1B
tGrid2.range = [-dt*ilw[2], dt*ilw[2]]
tGrid1.irange = [0, 2*ilw[2]]

resultStr.grid1=tGrid1
resultStr.grid2=tGrid2
;
; Create biCorrelation function object
;
biCorrelationObj = OBJ_NEW('GKVs2D', resultStr)
;
; Create output structures for spectral densities
;
spect0Str = {GKVs1D}
spect1Str = {GKVs1D}
spect2Str = {GKVs1D}
FOR i=0, N_TAGS({GKVsd}) DO BEGIN
	spect0Str.(i) = (temp[0]).(i)
	spect1Str.(i) = (temp[1]).(i)
	spect2Str.(i) = (temp[2]).(i)
ENDFOR
spect0Str.mnemonic='S_' + (temp[0]).mnemonic
spect1Str.mnemonic='S_' + (temp[1]).mnemonic
spect2Str.mnemonic='S_' + (temp[2]).mnemonic
spect0Str.title  = 'S{' + (temp[0]).title + '}'
spect1Str.title  = 'S{' + (temp[1]).title + '}'
spect2Str.title  = 'S{' + (temp[2]).title + '}'
indices1 = REPLICATE('*', 1)
spect0Str.Indices = PTR_NEW(indices1)
spect1Str.Indices = PTR_NEW(indices1)
spect2Str.Indices = PTR_NEW(indices1)
spect0Str.values = PTR_NEW(sspect0)
spect1Str.values = PTR_NEW(sspect1)
spect2Str.values = PTR_NEW(sspect2)
vmin = GKVsd_MIN(sspect0, MAX=vmax)
spect0Str.vrange = [vmin, vmax]
vmin = GKVsd_MIN(sspect1, MAX=vmax)
spect1Str.vrange = [vmin, vmax]
vmin = GKVsd_MIN(sspect2, MAX=vmax)
spect2Str.vrange = [vmin, vmax]
wGrid0 = GKVsd_Gridcopy(grid1)
wGrid0.mnemonic = 'omega0'
wGrid0.title = '!4x!X!D0!N'
spect0Str.Grid1 = wGrid0
spect1Str.Grid1 = GKVsd_GridCopy(grid1)
spect2Str.Grid1 = GKVsd_GridCopy(grid2)
;
; Create spectral density objects
;
spect0Obj = OBJ_NEW('GKVs1D', spect0Str)
spect1Obj = OBJ_NEW('GKVs1D', spect1Str)
spect2Obj = OBJ_NEW('GKVs1D', spect2Str)
;
; Create output structures for correlation functions
;
corr0Str = {GKVs1D}
corr1Str = {GKVs1D}
corr2Str = {GKVs1D}
FOR i=0, N_TAGS({GKVsd}) DO BEGIN
	corr0Str.(i) = (temp[0]).(i)
	corr1Str.(i) = (temp[1]).(i)
	corr2Str.(i) = (temp[2]).(i)
ENDFOR
corr0Str.mnemonic='C_' + (temp[0]).mnemonic
corr1Str.mnemonic='C_' + (temp[1]).mnemonic
corr2Str.mnemonic='C_' + (temp[2]).mnemonic
corr0Str.title  = 'C{' + (temp[0]).title + '}'
corr1Str.title  = 'C{' + (temp[1]).title + '}'
corr2Str.title  = 'C{' + (temp[2]).title + '}'
indices1 = REPLICATE('*', 1)
corr0Str.Indices = PTR_NEW(indices1)
corr1Str.Indices = PTR_NEW(indices1)
corr2Str.Indices = PTR_NEW(indices1)
corr0Values = corr0[(n_Corrs/2-ilw[0]):(n_Corrs/2+ilw[0])]
corr0Str.values = PTR_NEW(corr0Values)
corr1Values = corr1[(n_Corrs/2-ilw[1]):(n_Corrs/2+ilw[1])]
corr1Str.values = PTR_NEW(corr1Values)
corr2Values = corr2[(n_Corrs/2-ilw[2]):(n_Corrs/2+ilw[2])]
corr2Str.values = PTR_NEW(corr2Values)
vmin = GKVsd_MIN(corr0Values, MAX=vmax)
corr0Str.vrange = [vmin, vmax]
vmin = GKVsd_MIN(corr1Values, MAX=vmax)
corr1Str.vrange = [vmin, vmax]
vmin = GKVsd_MIN(corr2Values, MAX=vmax)
corr2Str.vrange = [vmin, vmax]

tau0s = dt*(-ilw[0] + FINDGEN(2*ilw[0]+1))
tGrid0 = {grid}
tGrid0.mnemonic = 'tau0'
tGrid0.title = '!4s!X!D0!N'
tGrid0.units = (temp[0]).grid1.units
tGrid0.values = PTR_NEW(tau0s)
tGrid0.boundary = 'periodic (closed)'
tGrid0.uniform = 1B
tGrid0.range = [-dt*ilw[0], dt*ilw[0]]
tGrid0.irange = [0, 2*ilw[0]]

corr0Str.Grid1 = tGrid0
corr1Str.Grid1 = GKVsd_GridCopy(tGrid1)
corr2Str.Grid1 = GKVsd_GridCopy(tGrid2)
;
; Create correlation function objects
;
corr0Obj = OBJ_NEW('GKVs1D', corr0Str)
corr1Obj = OBJ_NEW('GKVs1D', corr1Str)
corr2Obj = OBJ_NEW('GKVs1D', corr2Str)

;
; Create output structure
;
result = {	Name:'biSpectra', 		$
		tauLags:dt*ilw, 		$
		tauData:dt*idw,			$
		biSpectra:biSpectraObj, 	$
		biCoherence:biCoherenceObj, 	$
		biCorr:biCorrelationObj,	$
		spect0:spect0Obj,		$
		spect1:spect1Obj,		$
		spect2:spect2Obj,		$
		corr0:corr0Obj,			$
		corr1:corr1Obj,			$
		corr2:corr2Obj			}
;
; Clean up
;
FOR i=0,2 DO temp[i] -> trash
;
; and we're done ...
;
RETURN, result
END ; ****** GKVs1D::BiSpect ******  ;

FUNCTION GKVs2D::BiSpect, arg1, arg2, 	dw=ddw, Idw=idataWindow, 	$
					lw=llw, ilw=iLagWindow, order = iorder
;
; Purpose:
;
;	Dummy routine to short-stop attempts to 
;	call BiSpect on data with dimensionality 
;	greater than 1.
;
; Written by W.M. Nevins
;	7/6/01
;
MESSAGE, 'BiSpect only implimented for 1-D objects', /INFORMATIONAL
RETURN, 0
END ; ****** GKVs2D::BiSpect ****** ;