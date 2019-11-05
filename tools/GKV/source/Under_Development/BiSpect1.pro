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
GKV_DeTrend, v1, t, order=order
;
; Pack values into arrays with power-of-2 length
;
idw = (nt/10) > 1
IF(N_ELEMENTS(idataWindow) EQ 1) THEN idw = dataWindow
windowFcn = Window_Fcn(nt, idw)
j=0L
repeat j=j+1 until (2L^j gt nt) 	; Find n_corrs such that 2^n_corrs > nt
n_corrs = 2L^j
v0arr = COMPLEXARR(n_corrs)
v1arr = COMPLEXARR(n_corrs)
v2arr = COMPLEXARR(n_corrs)
v0arr[0:nt-1] = v0*windowFcn
v1arr[0:nt-1] = v1*windowFcn
v2arr[0:nt-1] = v2*windowFcn
;
; Compute fourier transforms
;
ft0 = FFT(v0arr)
ft1 = FFT(v1arr)
ft2 = FFT(v2arr)
;
; Form 'bare' bispectra
;
biSpectra = ft1#ft2
ft0star = CONJ(ft0)
FOR i=0,n_Corrs-1 DO BEGIN
	biSpectra[*,i] = biSpectra[*,i]*SHIFT(ft0star, -i)
ENDFOR
;
; Now form biCorrelation function
;
biCorrs = FFT(biSpectra, /INVERSE)
;
; Apply lag window
;
ilw = nt/2
IF(N_ELEMENTS(ilagWindow) EQ 1) THEN ilw=iLagWindow
lagWindow = LagWindow_Fcn(n_corrs, ilw)
lagWindow = SHIFT(lagWindow, n_corrs/2)
lagWindow = lagwindow#lagWindow
biCorrs = biCorrs*lagWindow
;
; Now compute averaged biSpectrum
;
biSpectra = FFT(biCorrs)
biSpectra = SHIFT(biSpectra, n_corrs/2, n_corrs/2)
;
; Prepare output structure
;
resultStr = {GKVs2D}
resultStr.mnemonic = 'B_' + (temp[0]).mnemonic + '_' + (temp[1]).mnemonic + '_' + (temp[2]).mnemonic

titleIndices = STRARR(3)
FOR i=0,2 DO BEGIN
	tIndices = temp[i] -> IndexRemove(1)
	tIndices = STRJOIN(tIndices, ',')
	titleIndices[i] = "[" + tIndices + ']'
ENDFOR
resultStr.title = 'B{'	+ (temp[0]).title + titleIndices[0] + '|'		$
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
gridValues = 2.*!PI/(dt*n_corrs)*(n_corrs/2 +1 - FINDGEN(n_corrs))
gridValues = SHIFT(GridValues, n_corrs/2)
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
; Create output object
;
result = OBJ_NEW('GKVs2D', resultStr)
;
; Clean up
;
FOR i=0,2 DO temp[i] -> trash
;
; and we're done ...
;
RETURN, result
END ; ****** GKVs1D::BiSpect ******  ;
