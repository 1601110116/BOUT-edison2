FUNCTION GKVs3D::TauCorrs,	arg, CorrFcns=CorrFcnArr, yCorrFcns=yCorrFcnArr, 	$
				rCorrFcns=rCorrFcnArr, kSpects=kSpectArr, _EXTRA=extra
;
; Purpose:
;
;		This function performs a correlation analysis of 'self'
;		(a GKVs3D object) as a function of the specified independent
;		variable (which is assumed to be  inhomogeneous). If the system 
;		is homogeneous, it is ***MUCH*** more efficient to simply 
;		compute the correlation function with Xcorr (with no reference 
;		object).
;
;
; Input Arguments:
;
;		GKVs3D::TauCorrs will accept any legal axis identifier
;		(that is, an integer between 1 and 3, or a valid
;		axis mnemonic) as an argument.  If no argument is 
;		provided, then GKVs3D::TauCorrs will expect the axis
;		corresponding to the inhomogeneous coordinate to be
;		identified by keywords (see below). (Optional)
;
;
;  Input Keywords:
;
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the inhomogeneous coordinate. Set axis 
;			equal to any legal axis identifier (see above).
;			(optional)
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the selected independent 
;			variable, and reset the signal window on this axis.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable (i.e., it is interpreted in the units
;			of the corresponding independent variable), NOT the integer 'irange'
;			(that is, NOT as an integer index into the grid.values array).
;			(optional)
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.  'irange'
;			defaults to the current signal window of 'self.
;			(optional)
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable.
;			(optional)
;
;	     skip	The sampling interval for the specified inhomogeneous independent
;			variable.  Defaults to 1 (Optional).
;
;	localNorm  	Set this keyword (i.e., put '/LocalNorm' on the command line) to normalize
;			the data from 'self' (but not 'self' itself, whose data remains unaltered)
;			such that the rms fluctuation at each value of the specified inhomogeneous
;			independent variable is the same BEFORE  computing the correlation time.
;			(Optional)
;
;	    debug	Set this keyword (i.e., put '/debug' on the command line) to 
;			print out intermediate information, allowing the user to keep track
;			of the progress (or lack thereof...) of this function. 
;			(optional)
;
;	 fraction	Normally TauCorrs returns the fullwidth at half maximum of the correlation
;			function vs. the specified axis.  If this keyword is set, RCORR will 
;			instead compute the full width at 'fraction' of maximum vs. time.
;			Defaults to 0.5 (Optional).
;			
;
;
; Output:  	GKVs3D::TauCorrs returns a a structure with the following tags.  
;		Additional tags can be included by setting output keywords 
;		(see "Output Keywords" below).
;		
;	TauCorr 	A GKVs1D object containing the correlation time
;			(full-width at half-maximum measured along the
;			maximum of the local correlation function) vs. 
;			the specified inhomogeneous independent variable.
;
;	vPhase		A GKVs1D object containing the phase velocity
;			(slope of the maximum in the local correlation 
;			function) vs. the inhomogenous independent variable. 
;
;	yCorr		A GKVs1D object containing the correlation length in
;			the remaining independent variable (full-width at
;			half maximum of the Hilbert-transform 'envelope'
;			of the oscillatory local correlation function evaluated
;			at zero time lag) vs. the inhomogeneous variable.;			
;
;	kAvg		A GKVs1D object containing the power-weighted avgerage  
;			wave number in the remaining (persumably homogeneous) 
;			independent variable vs. the inhomogeneous variable. 
;			The error bars on kAvg represent the power-weighted
;			standard deviation about this mean wave number.
;
;	kSpect		A GKVs1D object containing the average (over the inhomogeneous
;			variable) of the local frequency-integrated spectrum vs. the 
;			wave number associated with the remaining independent variable.
;			If the /LocalNorm keyword is set, then this is a volume average,
;			while if it is not set, then this is a power-weighted volume
;			average.
;
;	CorrFcn		A GKVs1D object containing the (volume or power-weighted volume as
;			with kSpect above) average over the inhomogeneous coordinate
;			of the maximum of the local correlation function vs. the time lag.
;
;	yCorrFcn	A GKVs1D object containing the (volume or power-weighted volume as
;			with kSpect above) average over the inhomogeneous coordinate of the
;			local corrlation function vs. the remaining independent variable
;			evaluated at zero time lag.
;
;
;  Output Keywords:
;
;	 CorrFcns	Set this keyword (i.e., put "/CorrFcns" on the command line)
;			to add the tag "TauCorrArr" to the output structure.  The
;			corresponding value is  an array of GKVs1D objects 
;			containing the local correlation function vs. tau for each location
;			at which the correlation time is computed. (Optional)
;
;	yCorrFcns	Set this keyword (i.e., put "/yCorrFcns" on the command line)
;			to add the tag "yCorrArr" to the output structure.  The
;			corresponding value is an array of GKVs1D objects
;			containing the local correlation function vs. the homogeneous
;			spatial coordinate at each radial location at which the correlation
;			function in computed. (Optional)
;
;	rCorrFcns	Set this keyword (i.e., put "rCorrFcns" on the command line)
;			to add the tag "rCorrArr" to the output sturcutre.  The 
;			corresponding value is an array of GKVs1D objects
;			containing the local correlation function vs. the 
;			inhomogeneous spatial coordinate at each radial location at
;			which the correlation function is computed.  (Optional)
;			
;	  kSpects	Set this keyword (i.e., put "/kSpects" on the command line)
;			to add the tag "kSpectArr" to the output structure.  The
;			corresponding value is an array of GKVs1D objects
;			containing the local k-spectrum (integrated over frequency)
;			vs. the homogeneous spatial coordinate at each radial location 
;			at which the correlation function in computed. (Optional)
;			
;
;
;
; Side Effects:
;
;		Resets the signal window of 'self' to specified 'range' or 'irange' if
;		one is provided with any of the keywords 'mnemonic', irange, or range
;	
;
; Written by W.M. Nevins
;	9/28/00
;
; Revised by W.M. Nevins
;	4/22/01
;
; Revised by W.M. Nevins
;	3/9/02
; to include radial correlation functions
FORWARD_FUNCTION GKVsd_MIN
;
; Get info from command line:
;
; First find identifier for axis corresponding to 
; the inhomogeneous independent variable
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	ELSE	:	BEGIN
				MESSAGE, 'tauCorrs called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
 ;
 ; Find (integer) sampling interval
 ; 
skip = 1
iskip = GetKeyWord('skip', extra)
IF(TypeOf(iskip) NE 7) THEN skip = FIX(iskip) > 1
;
; Form locally normalized (such that rms value is locally 1)
; copy of self
;
normedSelf = self
localNorm=0
locNorm = GetKeyWord('LocalNorm', extra)
IF (TypeOf(locNorm) NE 7) THEN localNorm = locNorm
IF KEYWORD_SET(localNorm) THEN BEGIN
	normedSelf = self -> LocalNorm(axis, AvgValSq=delSqObj)
ENDIF
;
; Check for 'fraction' keyword
;
fraction=0.5
fract = GetKeyWord('fraction', extra)
IF(TypeOF(fract) NE 7) THEN fraction=fract
;
; Check for 'debug' flag
;
debug=0
dbug = GetKeyWord('debug', extra)
IF(TypeOf(dbug) NE 7) THEN debug=dbug
;
; Create GKV structure to hold result,
; and copy fields from 'self' to 'result'
;
result = {GKVs1D}						
FOR i=0, N_TAGS(result) - 1 DO 	$
	result.(i) = self.(i)	
;
; Get grid structure corresponding to 
; the inhomogeneous independent variable
;
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'grid = normedSelf.Grid' + axisString
ok = EXECUTE(commandString)
IF NOT ok THEN BEGIN
	MESSAGE, "Could not extract inhomogeneous Grid Structure:  ", /INFORMATIONAL
	MESSAGE, "Command String = " + commandString, /INFORMATIONAL
	RETURN, 0
ENDIF
IF(N_ELEMENTS(arg) EQ 0) THEN arg=Grid.mnemonic	
oldGridValues = *Grid.values
;
; Get grid structure corresponding to the 
; homogeneous variable
;
CASE axis OF
	1 :	yGrid = normedSelf.grid2
	2 :	yGrid = normedSelf.grid1
	ELSE :	BEGIN
			MESSAGE, 'Could not extract homogeneous Grid structure', /INFORMATIONAL
			RETURN, 0
		END
ENDCASE
yGridValues = *yGrid.values
nypoints = N_ELEMENTS(yGridValues)
ly = yGridValues[nyPoints-1] - yGridValues[0]
dy = yGridValues[1] - yGridValues[0]
IF(yGrid.boundary EQ 'periodic (open)') THEN ly = ly+dy
;
; Set up loop over inhomogeneous independent variable
;
gridTitle = Grid.title
gridMnemonic = Grid.mnemonic
irange = Grid.irange
imax = irange[1]-irange[0] 
nPoints = (imax +1)/skip
values = FLTARR(nPoints)
newGridValues = FLTARR(nPoints)
values = FLTARR(nPoints)
errors = FLTARR(nPoints)
vPhase = FLTARR(nPoints)
vPhaseErrors = FLTARR(nPoints)
yValues = FLTARR(nPoints)
yErrors = FLTARR(nPoints)
yHalfValues = FLTARR(nPoints)
yHalfErrors = FLTARR(nPoints)
kValues = FLTARR(nPoints)
kErrors = FLTARR(nPoints)
rCorrValues = FLTARR(nPoints)
rCorrErrors = FLTARR(nPoints)
avgSqValues = FLTARR(nPoints)
CorrObjs = OBJARR(nPoints)
yCorrObjs = OBJARR(nPoints)
kSpectObjs = OBJARR(nPoints)
rCorrObjs = OBJARR(nPoints)
;
; Begin loop over inhomogeneous independent variable
;
FOR i=0, nPoints-1 DO BEGIN						; 'i' is index into 'result' arrays, while 'index' is the  
	index = irange[0]+i*skip					; index to the inhomogeneous independent variable of 'self'. 
	newGridValues[i] = oldGridValues[index]				; 	Set new grid values
	tempObj = normedSelf -> Slice(axis=axis, index=index)		;	Slice 'normedSelf' at current value of inhomogeneous variable		
	tempcorrs = tempObj -> XCorr()					;	Compute local auto-correlation function
	tempSpect = tempObj -> XSpect()					;	Compute local spectral density
	tempStats = tempObj -> Stats()					;	Compute local intensity
	avgSqValues[i] = tempStats.std^2
	; 
	; Slice vs. tau along the maximum of spatial separation in the remaining (homogeneous) spatial coordinate.
	;
	tempCorrsmax = tempCorrs -> Slice(axis=1, /AtMax, /maxlocation)
	;
	; Slice vs. homogeneous coordinate at tau=0.
	;
	tempCorrs1 = tempCorrs -> Slice(tau=0.)
	tempCorrs1 -> Get, title=yCorrTitle, mnemonic=yCorrMnemonic
	tempCorrs0 = tempCorrs1 -> Execute("FLOAT")
	tempCorrs0 -> Set, title=yCorrTitle, mnemonic=yCorrMnemonic
	tempCorrs2 = tempCorrs1 -> Envelope()
	tempCorrs3 = tempCorrs2 -> Execute("FLOAT")
	tempCorrs1 -> Trash
	tempCorrs2 -> Trash
	;
	; Compute full width at 1/2 maximum (and associate error bars)
	;
	values[i] = tempCorrsmax.slice -> FullWidth(debug=debug, Error=localError, Fraction=fraction)
	errors[i] = localError
	yValues[i] = tempCorrs3 -> FullWidth(debug=debug, Error=yError, Fraction=fraction)
	yErrors[i] = yError
	tempCorrs3 -> Trash
	yHalfValues[i] = tempCorrs0 -> FullWidth(debug=debug, Error=yHalfError, Fraction=fraction)
	yHalfErrors[i] = yHalfError
	;
	; Compute mean k-value and width
	;
	tempSpect -> get, title=spectTitle
	wMoments = tempSpect -> Moments(axis="omega")	; This delivers the integral of
	wMoments[1] -> Trash				; 'tempSpect' over omega.
	wMoments[2] -> Trash
	kSpect = wMoments[0]
	kSpect -> Get, axis=1, Range=kRange
	kSpect -> Set, title=spectTitle
	kRange[0] = 0.
	kStats = kSpect -> Moments(axis=1, Range=kRange, /Avg)
	kValues[i] = kStats[1] -> GetValues()
	kerrors[i] = kStats[1] -> GetErrors()
	FOR ixx=0,2 DO kStats[ixx] -> Trash
	;
	; Compute phase velocity (and associate error bars)
	;
	tempCorrsmax.maxlocation -> SignalWindow, tau=[-values[i]/2., values[i]/2.]
	tempCorrsMax.maxLocation -> Restrict
	tempCorrsMax.maxLocation -> RemoveShift, ly=ly
	vPhase[i] = tempcorrsmax.maxlocation -> slope(error=error)
	vPhaseErrors[i] = error
	; 
	; Provide some output if requested
	;
	;
	; Compute radial correlation function
	;
	rCorrObjs[i] = normedSelf -> XCORR0(axis=axis, ref=tempObj, /Norm)
	rCorrObjs[i] -> center, arg, index=index
	rCorrValues[i] = rCorrObjs[i] -> FullWidth(debug=debug, Error=rCorrError, fraction=fraction)
	rCorrErrors[i] = rCorrError

	IF(KEYWORD_SET(debug)) THEN BEGIN
		print, index, newGridValues[i], rCorrvalues[i], rCorrerrors[i], vPhase[i], vPhaseErrors[i]
		IF(i EQ 0) THEN tempCorrsMax.slice -> draw, /pretty ELSE tempCorrsMax.slice -> oplot
	ENDIF
	;
	; Clean up temporary objects 
	;
	tempObj    -> Trash
	tempcorrs  -> Trash
	tempcorrsmax.maxlocation -> Trash
	tempSpect -> Trash
	; 
	; Get Normalization for local corrlation functions
	;
	norm0 = tempCorrs0 -> slice(axis=1, value=0.)
	tempcorrs0 -> Get, title=title0, mnemonic=mnemonic0
	norm1 = tempCorrsMax.slice -> slice(axis=1, value=0.)
	tempCorrsMax.slice -> Get, title=title1, mnemonic=memonic1
	;
	; Save local correlation info
	;
	yCorrObjs[i]=tempCorrs0 -> over(norm0, title=title0, mnemonic=mnemonic0, units="")
	CorrObjs[i]=tempCorrsMax.slice -> Over(norm1, title=title1, mnemonic=mnemonic1, units="")
	kSpectObjs[i]=kSpect
	;
	; and clean up some more
	;
	tempcorrs0 -> Trash
	norm0 -> Trash
	GKVdelete, tempCorrsMax
	norm1 -> Trash
ENDFOR
;
; Clean up normedSelf if necessary
;
IF KEYWORD_SET(localNorm) THEN BEGIN
	normedSelf -> Trash
ENDIF
;
; Set fields of output Grid structure
;
Grid.values = PTR_NEW(newGridValues)
Grid.uniform = GKVsd_UniformGrid(newGridValues)
IF( (Grid.irange[0] NE 0) OR (Grid.irange[1] NE (N_ELEMENTS(oldGridValues)-1)) ) THEN Grid.boundary = 'open'
xmin = MIN(newGridValues, MAX=xmax)
Grid.range  = [xmin,xmax]
Grid.irange = [0,nPoints-1]
;
; Set fields of 'result' structure (containing taucorr vs. inhomogeneous coordinate)
;
axisMnemonic = Grid.mnemonic
axisTitle = Grid.title
result.mnemonic = 'tau_corr'
result.title = '!4s!X!Lc!N {' + self.title + '}'
CASE axis OF								;
	1:	indices = Self -> IndexRemove([2,3])			; Remove 'indices' corresponding to
	2:	indices = Self -> IndexRemove([1,3])			; the homogeneous independent variables
	3:	indices = Self -> IndexRemove([1,2])			; 
ENDCASE
result.indices = PTR_NEW(indices)
result.units = self.Grid3.units
result.values = PTR_NEW(values)
vmin = GKVsd_MIN(values, Max=vmax)
result.vrange = [0.0, vmax]
result.ErrorBars = PTR_NEW(errors)
result.Grid1 = Grid
;
; register 'result' as a GKVs1D object
;
tauCorrObj = OBJ_NEW('GKVs1D', result)
;
; Now make a vPhase object
;
vPhaseObj = taucorrObj -> MakeCopy(/noValues, /noErrorBars)
vPhaseObj.values = PTR_NEW(vPhase)
vmin = GKVsd_MIN(vPhase, Max=vmax)
vPhaseObj.vrange=[vmin, vmax]
vPhaseObj.title = 'V!I!4u!N!X'
vPhaseObj.mnemonic = 'v_phase'
CASE axis OF
	1:	units = self.Grid2.units + '/' + self.grid3.units
	2:	units = self.Grid1.units + '/' + self.grid3.units	
	3:	units = self.Grid1.units + '/' + self.grid2.units
ENDCASE
vPhaseObj.units = units
vPhaseObj.ErrorBars = PTR_NEW(vPhaseErrors)
;
; Now make a yCorr object
;
yCorrObj = taucorrObj -> MakeCopy(/noValues, /noErrorBars)
yCorrObj.values = PTR_NEW(yValues)
vmin = GKVsd_MIN(yValues, Max=vmax)
yCorrObj.vrange=[0.0, vmax]
CASE axis OF
	1:	BEGIN
			mnemonic = self.Grid2.mnemonic + "_corr"
			halfMnemonic = self.Grid2.mnemonic + "_1/2"
			rMnemonic = self.Grid1.mnemonic + "_corr"
			title = "!12l!X" + SubScript(self.Grid2.title) + "!N"
			halfTitle = "!12l!X!S" + SubScript(self.Grid2.title) + "!R!U(1/2)!N"
			rTitle = "!12l!X" + SubScript(self.Grid1.title) + "!N"
			units  = self.Grid2.units
			rUnits = self.Grid1.units
		END
	2:	BEGIN
			mnemonic = self.Grid1.mnemonic + "_corr"
			halfMnemonic = self.Grid1.mnemonic + "_1/2"
			rMnemonic = self.Grid2.mnemonic + "_corr"
			title = "!12l!X" + SubScript(self.Grid1.title) + "!N"
			halfTitle = "!12l!X!S" + SubScript(self.Grid1.title) + "!R!U(1/2)!N"
			rTitle = "!12l!X" + SubScript(self.Grid2.title) + "!N"
			units  = self.Grid1.units
			rUnits = self.Grid2.units
		END	
	3:	BEGIN
			mnemonic = self.Grid1.mnemonic + "_corr"
			halfMnemonic = self.Grid1.mnemonic + "_1/2"
			rMnemonic = self.Grid3.mnemonic + "_corr"
			title = "!12l!X" + SubScript(self.Grid1.title) + "!N"
			halfTitle = "!12l!X!S" + SubScript(self.Grid1.title) + "!R!U(1/2)!N"
			rTitle = "!12l!X" + SubScript(self.Grid3.title) + "!N"
			units  = self.Grid1.units
			rUnits = self.Grid3.units
		END
ENDCASE
yCorrObj.title = title
yCorrObj.mnemonic = mnemonic
yCorrObj.units = units
yCorrObj.ErrorBars = PTR_NEW(yErrors)
;
; Now make a yHalf object
;
yHalfObj = yCorrObj -> MakeCopy(/noValues, /noErrorBars)
yHalfObj.title = halfTitle
yHalfObj.mnemonic = halfMnemonic
yHalfObj.values = PTR_NEW(yHalfValues)
vmin = GKVsd_MIN(yHalfValues, Max=vmax)
vmin = vmin < 0.
IF(vmax EQ vmin) THEN vmax=vmin+1.
yhalfObj.vrange=[vmin,vmax]
yHalfObj.errorBars = PTR_NEW(yHalfErrors)
;
; Now make a kavg object
;
kAvgObj = taucorrObj -> MakeCopy(/noValues, /noErrorBars)
kAvgObj.values = PTR_NEW(kValues)
vmin = GKVsd_MIN(kValues, Max=vmax)
vmin = vmin < 0.
vmax = vmax > 1.
kAvgObj.vrange=[vmin, vmax]
CASE axis OF
	1:	BEGIN
			mnemonic = "k_" + self.Grid2.mnemonic
			title = "k" + SubScript(self.Grid2.title) + "!N"
			units = '1/(' + self.Grid2.units + ')'
		END
	2:	BEGIN
			mnemonic = "k_" + self.Grid1.mnemonic
			title = "k" + SubScript(self.Grid1.title) + "!N"
			units = '1/(' + self.Grid1.units + ')'
		END	
	3:	BEGIN
			mnemonic = "k_" + self.Grid1.mnemonic
			title = "k" + SubScript(self.Grid1.title) + "!N"
			units = '1/(' + self.Grid1.units + ')'
		END
ENDCASE
kAvgObj.title = "!12<!X" + title + "!12>!X"
kAvgObj.mnemonic = mnemonic
kAvgObj.units = units
kAvgObj.ErrorBars = PTR_NEW(kErrors)
;
; and AvgSQ object
;
avgSqObj = taucorrObj -> MakeCopy(/noValues, /NoErrorBars)
avgSqObj.values = PTR_NEW(avgSqValues)
vmin = GKVsd_MIN(avgSqValues, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin + 1.
avgSqObj.vrange = [vmin, vmax]
avgSqObj.title = "!12<!X(" + self.title + ")!U2!N!12>!X"
avgSqObj.units = ''
IF(self.units NE '') THEN avgSqObj.units = "(" + self.units + ")!U2!N"
avgSqobj.mnemonic = "Avg_" + self.mnemonic + "Sq"
;
; make an rCorr object
;
rCorrObj = tauCorrObj -> MakeCopy(/noValues, /NoErrorBars)
rCorrObj.title = rtitle
rCorrObj.mnemonic = rMnemonic
rCorrObj.units = rUnits
rCorrObj.values = PTR_NEW(rCorrValues)
rCorrObj.ErrorBars = PTR_NEW(rCorrErrors)
vmin = GKVsd_MIN(rCorrValues, MAX=vmax)
rcorrObj.vrange = [vmin, vmax]
;
; Now compute averages of local correlation info
;
kSpectObjs[0] -> get, title=kSpectTitle, mnemonic=kSpectMnemonic
kSpectTitle = '!12<!X' + kSpectTitle + '!12>!X' + SubScript(GridTitle) + "!N"
kSpectMnemonic = 'Avg' + kSpectMnemonic + '_' + gridMnemonic
kSpectObj =  kSpectObjs[0] -> ObjAvg(kSpectObjs, title=kSpectTitle, mnemonic=kSpectMnemonic)

corrObjs[0] -> get, title=corrFcnTitle, mnemonic=corrFcnMnemonic
corrFcnTitle = '!12<!X' + corrFcnTitle + '!12>!X' + SubScript(GridTitle) + "!N"
corrFcnMnemonic = 'Avg' + corrFcnMnemonic + '_' + gridMnemonic
temp =  corrObjs[0] -> ObjAvg(corrObjs)
;
; Normalize the correlation function such that its value at zero lag is 1.0
;
norm = temp -> slice(tau=0.)
corrFcnObj = temp -> Over(norm, title=corrFcnTitle, mnemonic=corrFcnMnemonic, units="")
temp -> Trash
norm -> Trash

yCorrObjs[0] -> get, title=yCorrFcnTitle, mnemonic=yCorrFcnMnemonic
yCorrFcnTitle = '!12<!X' + yCorrFcnTitle + '!12>!X' + SubScript(GridTitle) + "!N"
yCorrFcnMnemonic = 'Avg' + ycorrFcnMnemonic + '_' + gridMnemonic
temp =  yCorrObjs[0] -> ObjAvg(yCorrObjs)
;
; and normalize ...
;
norm = temp -> slice(axis=1, value=0.)
yCorrFcnObj = temp -> over(norm, title=yCorrFcnTitle, mnemonic=yCorrFcnMnemonic, units="")
temp -> Trash
norm -> Trash
;
; ****** Need to 'shift' the rCorrFcns BEFORE we average them! ******
;
rCorrObjs[0] -> get, title=rCorrFcnTitle, mnemonic=rCorrFcnMnemonic
rCorrFcnTitle = '!12<!X' + rCorrFcnTitle + '!12>!X' + SubScript(GridTitle) + "!N"
rCorrFcnMnemonic = 'Avg' + rcorrFcnMnemonic + '_' + gridMnemonic
temp =  rCorrObjs[0] -> ObjAvg(rCorrObjs)
;
; and normalize ...
;
norm = temp -> slice(axis=1, value=0.)
rCorrFcnObj = temp -> over(norm, title=yCorrFcnTitle, mnemonic=yCorrFcnMnemonic, units="")
temp -> Trash
norm -> Trash
;
; Construct output structure
;
output = {	Name	:	"TauCorrs", 	$
		TauCorr	:	taucorrObj, 	$
		vPhase	:	vPhaseObj, 	$
		yCorr	:	ycorrObj, 	$
		yHalf	:	yHalfObj,	$
		kAvg	:	kAvgObj, 	$
		rCorr	:	rCorrObj,	$
		avgSQ	:	avgSqObj,	$
		kSpect	:	kSpectObj,	$
		corrFcn	:	corrFcnObj,	$
		yCorrFcn:	yCorrFcnObj,	$
		rCorrFcn:	rCorrFcnObj	}
;
; Now, check for output keywords
;
IF KEYWORD_SET(LocalNorm) THEN BEGIN
	output = CREATE_STRUCT(output, "delSq", delSqObj)
ENDIF
IF KEYWORD_SET(CorrFcnArr) THEN BEGIN
	output = CREATE_STRUCT(output, "TauCorrArr", corrObjs)
ENDIF ELSE BEGIN
	GKVdelete, corrObjs
ENDELSE

IF(KEYWORD_SET(yCorrFcnArr)) THEN BEGIN
	output = CREATE_STRUCT(output, "yCorrArr", yCorrObjs)
ENDIF ELSE BEGIN
	GKVdelete, yCorrObjs
ENDELSE

IF(KEYWORD_SET(rCorrFcnArr)) THEN BEGIN
	output = CREATE_STRUCT(output, "rCorrArr", rCorrObjs)
ENDIF ELSE BEGIN
	GKVdelete, rCorrObjs
ENDELSE


IF(KEYWORD_SET(kSpectArr)) THEN BEGIN
	output = CREATE_STRUCT(output, "kSpectArr", kSpectObjs)
ENDIF ELSE BEGIN
	GKVdelete, kSpectObjs
ENDELSE
;
; and we're done!
;
RETURN, output
;
END ; ****** GKVs3D::TauCorrs ****** ;
