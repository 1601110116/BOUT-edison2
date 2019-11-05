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
IF Query_Complex(rsignal) THEN $
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
;	The algorithm used here is less of a memory hog than that
;	used by BiSpect, and should be employed on long data sets.
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
;		the bicoherence.  Such a limit is required because the
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
;	5/6/01
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
; Take complex conjugate of data from 'self'
;
v0 = CONJ(v0)
;
; Get detrending order from command line
;
order = 1
IF(N_ELEMENTS(iOrder) EQ 1) THEN order = iOrder
;
; Get length of the lag window from the command line
;
ilw = FIX(nt/10) > 2 
IF(N_ELEMENTS(llw) GE 1) THEN ilw = FIX(llw/dt) > 1
IF(N_ELEMENTS(ilagWindow) GE 1) THEN ilw=FIX(iLagWindow)
IF(N_ELEMENTS(ilw) EQ 1) THEN ilw = [ilw, ilw, ilw]
IF(N_ELEMENTS(ilw) EQ 2) THEN ilw = [ilw[0], ilw[0], ilw[1]]
ilws = ilw[1]
;
; Get length of the roll-off in the data window from the command line
;
idw = FIX(nt/100) > 1
IF(N_ELEMENTS(ddw) EQ 1) THEN idw = FIX(ddw/dt) > 1
IF(N_ELEMENTS(idataWindow) EQ 1) THEN idw = FIX(idataWindow)
IF(idw GT ilws/2) THEN idw=ilws/2
;
; Use lagwindow to determine the length of data subsets
;
j=0L
REPEAT j=j+1 UNTIL (2L^j GE (3*ilws+idw)) 	; Find n_corrs such that 2^n_corrs > 3*ilws+idw
n_Corrs = 2L^j
;
; Find number of subsets
;
nSubs = nt/(3*ilws)
nLeft = nt - nSubs*(3*ilws)
IF(nLeft GT (2*idw)) THEN nSubs=nSubs+1
;
; Set up index arrays to point at each data subset
;
startArr = LONARR(nSubs)
finishArr= LONARR(nSubs)
FOR i=0, nSubs-1 DO  startArr[i]=i*3*ilws
FOR i=0, nSubs-1 DO finishArr[i]=startArr[i]+3*ilws+idw-1
IF(nLeft GT (2*idw)) THEN finishArr[nSubs-1] = startArr[nSubs-1]+nLeft-1
ntEff=3*ilw*(nSubs-1)+nLeft*(nLeft GT (2*idw)) -idw
;
; Set up arrays to hold subSeries
;
v0arr = COMPLEXARR(n_corrs)
v1arr = COMPLEXARR(n_corrs)
v2arr = COMPLEXARR(n_corrs)
;
; and arrays in which to accumulate spectral densities, cross spectrum and normalization factor
;
spect0    = COMPLEXARR(n_Corrs)
spect1    = COMPLEXARR(n_Corrs)
spect2    = COMPLEXARR(n_Corrs)
biSpectra = COMPLEXARR(n_Corrs, n_Corrs)
norm      = COMPLEXARR(n_Corrs, n_Corrs)
;
; Start loop over subSeries
;
FOR isub=0, nSubs-1 DO BEGIN
	start = startArr[iSub]
	finish= finishArr[iSub]
	sv0 = v0[start:finish]
	sv1 = v1[start:finish]
	sv2 = v2[start:finish]
	st  =  t[start:finish]
	snt = finish -start + 1
	;
	; remove average and trend
	;
	GKV_DeTrend, sv0, st, order=order
	GKV_DeTrend, sv1, st, order=order
	GKV_DeTrend, sv2, st, order=order
	;
	; Pack values into arrays with power-of-2 length
	;
	windowFcn = Window_Fcn(snt, idw)
	windowNorm = SQRT( snt/TOTAL(windowFcn*windowFcn) ) ; ????****
	windowFcn = windowNorm*windowFcn


	v0arr[0:snt-1] = sv0*windowFcn
	v1arr[0:snt-1] = sv1*windowFcn
	v2arr[0:snt-1] = sv2*windowFcn
	;
	; Compute fourier transforms
	;
	ft0 = FFT(v0arr, /INVERSE)
	ft1 = FFT(v1arr, /INVERSE)
	ft2 = FFT(v2arr, /INVERSE)
	;
	; Form 'bare' bispectra
	;
	ttemp = (ft1#ft2) 
	ft0Backwards = SHIFT(REVERSE(ft0), 1)
	FOR i=0,n_Corrs-1 DO BEGIN
		biSpectra[*,i] = biSpectra[*,i]  + ttemp[*,i]*SHIFT(ft0Backwards, -i)	
	ENDFOR
	;
	; Compute 'bare' spectral densities
	;
	aft0Sq = FLOAT( ft0*CONJ(ft0) ) > 0.
	aft1Sq = FLOAT( ft1*CONJ(ft1) ) > 0.
	aft2Sq = FLOAT( ft2*CONJ(ft2) ) > 0.
	
	spect0 = spect0 + aft0Sq
	spect1 = spect1 + aft1Sq
	spect2 = spect2 + aft2Sq
	;
	; Compute normalization function for biCoherence
	;
	aft0 = ABS(ft0)
;	aft1 = abs(ft1)
;	aft2 = abs(ft2)
	aft0Backwards = SHIFT(REVERSE(aft0), 1)
	ttemp = ABS(ft1#ft2)			; See Kim and Powers, IEEE Trans. Plasma Sci. <PS-7> 2, 120 (June, 1979). 					
	FOR i=0,n_Corrs-1 DO BEGIN
		norm[*,i] = norm[*,i] + ttemp[*,i]*SHIFT(aft0Backwards, -i)
	ENDFOR
ENDFOR
;
; normalize spectral densities, biSpectra, norm
;
spect0 = spect0/nt
spect1 = spect1/nt
spect2 = spect2/nt
biSpectra = biSpectra/nt
norm = norm/nt
;
; Now form bare biCorrelation function
;
biCorrs = FFT(biSpectra)
;
; and bare auto correlation functions
;
corr0 = FFT(spect0)					
corr1 = FFT(spect1)
corr2 = FFT(spect2)
;
; Apply lag window to biCorrelation function
;
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
; Try filtering again?
;
;ttemp = FFT(biCoherence)
;ttemp = ttemp*fullLagWindow
;biCoherence = FFT(ttemp, /INVERSE)
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
		tIndices=tIndices[0]
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
END ; ****** GKVs2D::BiSpectL ****** ;