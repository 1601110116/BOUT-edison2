FUNCTION GKVs1D::MEM, FindSlope=FindSlope, growthRate=growth_rate, 	$
			tRange=trange, dTau=d_tau, nModes=n_modes, 	$
			Omegas=KeepOmegas, nOmegas=n_Omegas, 		$
			Norm=norm, debug=debug, _Extra=extra
;
; Purpose:
;
;	This function performs a "Maximum Entropy Method" spectral analysis
;	[see J.P. Burg, "A New Analysis Technique for Time Series Data",
;	Presented atthe nato Advanced Study Instituyte on Signal Processing 
;	with Emphasis on Underwater Acoustics, Aug. 12-23, 1968]
;	on the time series contained in 'Self'.
;
; Keywords
;
;	FindSlope	Set this keyword (i.e., put '/FindSlope' on the command line)
;			and this function will invoke the proto-widget 'Find_Slope'
;			to allow the user to select the relevant time interval and
;			and estimate any expotential growth rate in the signal.
;			Default is to use all data in 'self', and to assume
;			that there is no expotential growth. (Optional)
;
;	growthRate	Alternatively, the user may supply an estimate of the
;			maximum expotential growth in the signal.  Default
;			is to assume that there is no expotential growth
;			in the signal. (Optional)
;
;	Norm		Set this keyword (i.e., put '/Norm' on the command line)
;			to normalize the correlation function before performing
;			MEM spectral analysis (the integral over frequency of 
;			the spectral density should then be one, independent of
;			the amplitude of the input signal).  If a postive growthRate is
;			set (either through the keyword 'growthRate', or via the
;			Find_Slope proto-widget), then signal will be normed auto-
;			matically.  Otherwise, default is not to 'norm' the correlation
;			function. (Optional)
;
;	dTau		Interval for sampling the data.	 Default is to use
;			the (minimum) time step of 'Self'. (Optional)
;
;	nModes		Upper limit to the order of the MEM analysis.  
;			Actual order may be less than 'nModes' if a local
;			minimum in the 'FPE' is found [see H. Akaike,
;			"A New Look at the Statistical Model Identification"
;			IEEE Trans. Autom. Control, Vol. AC-19, pp. 716-723 
;			(Dec. 1974).  Defaults to 3.  (Optional)
;
;	Omegas		Set this keyword (i.e., put '/Omegas' on the command line)
;			to compute the complex normal mode frequencies of the
;			MEM forward prediction 'filter'.  Default is not to compute
;			and return these frequencies.  (Optional) 
;
;	nOmegas		Number of points in 'omega'-grid of the MEM spectral
;			density object which is returnedin the output structure.
;			Defaults to the number of points within 'irange' of
;			"self's" time grid. (Optional)
;			
;
;	
;
; Output
;
;		...
;
; Written by W.M. Nevins
;	12/30/01
; Modified by W.M. Nevins to improve efficiency
;	1/25/02
;
growthRate=0.
IF(N_ELEMENTS(growth_rate) EQ 1) THEN growthRate=growth_rate > 0.
PMgrowthRate=0.
;
; Select time interval, and check for 
; expotential growth of signal
;
IF KEYWORD_SET(FindSlope) THEN BEGIN
	self -> Find_Slope, slope=growthRate, Xout=t_range, PMslope=PMgrowthRate, /LinearAnalysis
	IF(N_ELEMENTS(trange) NE 2) THEN trange=t_range
ENDIF
;
; Make copy of 'self' (so self will be left unaltered by this proceedure)
;
thiscopy = self -> MakeCopy()
IF(KEYWORD_SET(norm)) THEN thisCopy -> NORM, /NoAvg
;
; Restrict data to time interval of interest
;
IF(N_ELEMENTS(trange) EQ 2) THEN BEGIN
	thiscopy -> signalwindow, axis=1, range=trange
	thiscopy -> restrict
ENDIF ELSE BEGIN
	trange = thisCopy.grid1.range
ENDELSE
;
; Remove any expotential growth
;
IF(growthRate GT 0) THEN BEGIN
	big = EXP( growthRate*(trange[1]-trange[0]) )
	big = big/thisCopy.vrange[1]
	temp = thisCopy -> Times(big)	
	newCopy = temp -> TimesExp(-growthRate)
	thiscopy -> trash
	temp -> trash
	thiscopy = newCopy
	thisCopy -> NORM, /NoAvg
ENDIF
;
; Make axis uniform (for computation of correlation function)
;
thiscopy -> ScaleAxis, 1, /Uniform
;
; Select order of MEM analysis
;
nModes = 3
IF(N_ELEMENTS(n_modes) EQ 1) THEN nModes=n_modes > 1
;
; Select sampling interval for MEM analysis
;
tValues = *thisCopy.Grid1.values
dt = tValues[1] - tValues[0]
dTau = dt
IF(N_ELEMENTS(d_tau) EQ 1) THEN dTau=d_tau > dTau
iTau = FIX(dTau/dt) > 1
dTau = iTau*dt
;
; Set up arrays for MEM analysis...
;
; Array for correlation function, sampled at intervals 'dTau'
;
corrs = COMPLEXARR(nModes+1)
theseValues = *thisCopy.values
info = size(theseValues)
nVals = info[1]
info[1] = nVals + nModes*iTau
values = MAKE_ARRAY(SIZE=info)
values[0:(nVals-1)] = theseValues
FOR i=0, nModes DO corrs[i] = TOTAL( values*CONJ(SHIFT(values,-i*iTau)) )/(nVals-i*iTau)
norm = corrs[0]
corrs = corrs/norm
;
; an array for Berg's lambda's (=exp(-i*omega_j*dTau)
; (we call this array 'a' to avoid carpal tunnel syndrome ...)
; 
a = COMPLEXARR(nModes+1)
a[0] = 1.
a[1] = -CONJ(corrs[1])/corrs[0]
aa = a
;
; an array for the residual power
;
p = FLTARR(nModes+1)
p[0] = 1.
p[1] = 1 - a[1]*CONJ(a[1])
;
; and an array for the FPE
;
irange = thisCopy.grid1.irange
nPoints = irange[1] - irange[0]
FPE = FLTARR(nModes+1)
FPE[0] = 1.*(nPoints + 1)/(nPoints - 1.)
FPE[1] = p[1]*(nPoints + 2.)/(nPoints - 2.)
;
; Now, use Berg's iteration to compute the 
; coefficients in each higher order
;
FOR m=2, nModes DO BEGIN
	delta = COMPLEX(0.,0.)
	FOR j=0,m-1 DO delta = delta + a[j]*CONJ(corrs[m-j])
	a[m] = - delta/p[m-1]
	p[m] = p[m-1]*( 1. - a[m]*CONJ(a[m]) )
	FPE[m] = p[m]*(nPoints + 1. + m)/(nPoints - 1. - m)
	IF(FPE[m] GE FPE[m-1]) THEN GOTO, DONE	; terminate iteration when FPE reaches a minimum
	IF(FPE[m] LE 0) THEN GOTO, DONE		;  or if FPE gets too strange ...
	FOR j=1,m-1 DO a[j]=aa[j] + a[m]*CONJ(aa[m-j])
	aa=a
ENDFOR
m = nModes + 1
DONE	:
nModesFound = m-1
;
; Create object to hold MEM version of spectral density
; and set fields
;
spectStr = {GKVs1D}
FOR i=0, N_TAGS({GKVsd})-1 DO spectSTR.(i) = self.(i)
spectStr.title = 'S!DMEM!N{' + thisCopy.title + '}'
spectStr.mnemonic = 'S_' + thisCopy.mnemonic
indices = *thisCopy.indices
spectStr.indices = PTR_New(indices)
;
; Create omega-grid
;
wGrid = {Grid}
tGrid = thisCopy.grid1
tvalues = *tGrid.values
omegaMax = 2.*!PI/dTau
IF((nPoints/2)*2 EQ nPoints) THEN nOmegas=nPoints+1 ELSE nOmegas=nPoints
IF(N_Elements(n_Omegas) EQ 1) THEN BEGIN
	IF((n_Omegas/2)*2 EQ n_Omegas) THEN nOmegas=n_Omegas+1 ELSE nOmegas=n_Omegas
ENDIF
dOmega = omegaMax/nOmegas
omega = dOmega*(FINDGEN(nOmegas) - nOmegas/2)
;
; Compute MEM estimate of spectral density
;
iImag = COMPLEX(0.,1.)
one   = COMPLEX(1.,0.)
x = exp(iImag*omega*dTau)
power = FLOAT(2.*p[nModesFound]*norm/omegaMax) > 0.
IF(power EQ 0.) THEN power=1.
dielec=REPLICATE(one, nOmegas)
FOR i=1, nModesFound DO dielec = dielec + aa[i]*x^i
spect = FLOAT( power/(dielec*CONJ(dielec)) ) > 0.
;
; Correct normalization
;
IntSpect = TOTAL(spect*dOmega)
spect = FLOAT(norm*spect/IntSpect)
;
; Store MEM spectral density in spectStr
;
spectStr.values = PTR_NEW(spect)
spectStr.ErrorBars = PTR_NEW()
vmin = MIN(spect, MAX=vmax)
spectStr.vrange=[vmin,vmax]
;
; Store wGrid in spectStr
;
wGrid.values = PTR_NEW(omega)
wGrid.mnemonic = 'omega'
wGrid.title = '!4x!X'
wGrid.units = '1/(' + tGrid.units + ')'
wGrid.boundary = 'periodic (closed)'
spectStr.grid1 = wGrid
spectStr.units = '(' + thisCopy.units + ')*(' + tGrid.units + ')'
IF KEYWORD_SET(FindSlope) THEN spectStr.units = ''
IF KEYWORD_SET(Norm)      THEN spectStr.units = ''
;
; Create GKVs1D object from spectStr
;
spectObj = OBJ_NEW('GKVs1D', spectStr)
;
; Create output structure
;

output = {	Name		:	'MEM_' + self.mnemonic,	$
		Spectrum	:	spectObj,		$
		nModes		:	nModesFound,		$
		gammaMax	:	growthRate,		$
		PMgammaMax	:	PMgrowthRate,		$
		P		:	p,			$
		FPE		:	FPE			}
;
; Compute roots of effective MEM dieletric
;
IF KEYWORD_SET(DeBug) THEN KeepOmegas = Debug
IF KEYWORD_SET(KeepOmegas) THEN BEGIN
	roots = FZ_ROOTS(aa[0:nModesFound])
	;
	; Now, solve for associated (complex) frequencies
	;
	omegas = COMPLEX(0.,-1.)*ALOG(roots)/dTau
	;
	; order frequencies such that those with smallest imaginary part
	; appear first (these make largest contributions to MEM spectral density ...)
	; Note that the 'gammas' are more about line width than growth vs. damping...
	;
	gammas  = IMAGINARY(omegas)
	absGammas  = ABS(gammas)
	indices = SORT(absGammas)
	gammas = gammas[indices]
	absGammas = absGammas[indices]
	;
	; Only keep real part of frequency ...
	; use any impressed growth rate for imaginary part
	;
	ReOmegas = FLOAT(omegas(indices))
	;
	; Print into to terminal if '/DeBug' is set
	;
	IF KEYWORD_SET(debug) THEN BEGIN
		Print, 'nModes found = ', nModesFound
		PRINT, 'FPE : ', fpe
		print, 'Omegas:'
		FOR i=0, nModesFound-1 DO PRINT, ReOmegas[i], ' +/- ', gammas[i]
		Print, 'Gamma:'
		Print, gammas
	ENDIF
	;
	; Add frequency into to output structure
	;
	output = CREATE_STRUCT(output, 'ReOmega', ReOmegas,'PMomega', absGammas)
ENDIF

;
; clean up
;
thisCopy -> Trash
;
; and, we're done ...
;
RETURN, output
END	; ****** GKVs1D::MEM ****** ;


FUNCTION GKVs2D::MEM, nOmega=n_Omega, _EXTRA=extra
;
; Purpose:
;
;	For each value of the first independent variable, an
;	MEM spectral analysis is performed over the second
;	independent variable.
;
; Written by W.M. Nevins
;	1/1/02
; Save contents of 'extra'
IF(N_ELEMENTS(extra) GT 0) THEN thisExtra = extra
;
; will build 'result' structure on copy of 'self'
;
result = self -> MakeCopy(/NoErrorBars)
result -> Restrict
;
; find number of omega bins for output
;
tGrid = result.grid2
irange_t = tGrid.irange
nt = irange_t[1] - irange_t[0] + 1
nOmega = nt
IF(N_ELEMENTS(n_Omega) EQ 1) THEN nOmega = n_Omega
IF( ((nOmega/2)*2) EQ nOmega) THEN nOmega = nOmega + 1
;
; find range of 1st independent variable
;
irange = result.grid1.irange
imin = irange[0]
imax = irange[1]
nx = imax - imin + 1
;
; get initial 'slice in 1st independent variable
;
firstSlice = result -> slice(axis=1, index=imin)
;
; Perform MEM analysis on first slice (and then trash it ...)
;
MEMstr = firstSlice -> MEM(nOmega=nOmega, _Extra=thisExtra)
firstSlice -> TRASH
;
; Extract needed info from this first spectrum
;
firstSpect = MEMstr.spectrum
wGrid = firstSpect.grid1
firstSpect.grid1 = {grid}	; so we can later trash firstSpect without loosing wGrid
result.title = firstSpect.title
result.mnemonic = firstSpect.mnemonic
result.units = firstSpect.units
theseValues = *firstSpect.values
;
; Set up arrays to hold values for the output objects
;
nOmegas = N_ELEMENTS(theseValues)
values = FLTARR(nx, nOmegas)
values[0,*] = theseValues
gammas = FLTARR(nx)
PMgammas = FLTARR(nx)
gammas[0]  = MEMstr.gammaMax
PMgammas[0]= MEMstr.PMgammaMax
FOR i=imin+1, imax DO BEGIN
	GKVdelete, MEMstr
	thisSlice = result -> slice(axis=1, index=i)
	IF(N_ELEMENTS(thisExtra) GT 0) THEN nextra = thisExtra
	MEMstr = thisSlice -> MEM(nOmega=nOmega, _Extra=nExtra)
	thisSpect = MEMstr.spectrum
	values[i,*] = *thisSpect.values
	gammas[i]  = MEMstr.gammaMax
	PMgammas[i]= MEMstr.PMgammaMax
	thisSlice -> Trash
ENDFOR
GKVdelete, MEMstr
vmin = MIN(values, MAX=vmax)
result.values = PTR_NEW(values)
result.vrange = [vmin,vmax]
result.grid2 = wGrid
;
; set up  structure for creation of gamma object
;
gammaStr = {GKVs1D}
FOR i=0, N_TAGS({GKVsd})-1 DO gammaStr.(i)=self.(i)
gammaStr.title='!4c!X!Dmax!N{' + self.title + '}'
gammaStr.mnemonic = 'Gamma_' + self.mnemonic
indices = self -> IndexRemove(2)
gammaStr.indices = PTR_NEW(indices)
gammaStr.units = '1/(' + tGrid.units + ')'
gammaStr.values = PTR_NEW(gammas)
vmin = MIN(gammas, MAX=vmax)
gammaStr.vrange = [vmin,vmax]
gammaStr.errorBars = PTR_NEW(PMgammas)
xGrid = GKVsd_Gridcopy(result.grid1)
gammaStr.grid1=xGrid
gammaObj = OBJ_NEW("GKVs1D", gammaStr)

output = {	Name		:	'MEM2D',	$
		spectrum	:	result,		$
		gammaMax	:	gammaObj	}

RETURN, output
END 	;  ****** GKVs2D::MEM ******  ;


FUNCTION GKVs3D::MEM, _Extra=extra
;
; Purpose:
;
;	Dummy to prevent calls to MEM acting on 
;	objects with dimensionality greater than 1.
;
;
; Written by W.M. Nevins
;	12/30/01
;
MESSAGE, "MEM is only defined for objects of dimensionality less than 3.  Will return 0", /INFORMATIONAL
RETURN, 0
END
