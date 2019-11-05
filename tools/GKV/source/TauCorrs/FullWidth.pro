
PRO GKVsd::Norm, MaxValue=maxVal, AvgValueSq=avgValuesSq, NoAvg=NoAvg
;
; Purpose:
;
;		Normalize input object such that it's RMS value, 
;		over interval selected by the 'irange' of each
;		independent variable, is maxVal.
;
; Arguments:
;
;		None
;
; Input Keywords:
;
;		MaxValue	On return 'self' is normalized such that the RMS value is
;				equal to 'MaxValue'. Defaults to 1 (Optional).
;
;		NoAvg		Set this keyword (i.e., put '/NoAvg' on the command line)
;				to remove the average BEFORE norming the data in 'self'
;
; Output Keywords:
;
;		AvgValueSq	On return, AvgValueSq contains the average of the squares
;				of 'self.values' (that is, the initial values of 'self.values')
;				over the interval selected by the 'irange' for each
;				independent variable
;
; Written By W.M. Nevins
; 	5/7/00
;
; Modified by W.M. Nevins
;	5/23/01
; to restrict average to interval selected by 'irange' and include  the output keyword 'AvgValueSq'
;
FORWARD_FUNCTION GKVsd_MIN
maxvalue=1.
IF(KEYWORD_SET(maxval)) THEN maxvalue=maxval
valuePtr = self -> GetValues()
values = *valuePtr
IF(KEYWORD_SET(NoAvg)) THEN BEGIN
	avgValue = TOTAL(values)/N_ELEMENTS(values)
	values = values - avgValue
ENDIF
avgValuesSq = TOTAL(values*CONJ(values))/N_ELEMENTS(values)
avgValuesSq = FLOAT(avgValuesSq)
oldVAlues = *self.values
values = maxValue*oldValues/SQRT(avgValuesSq)
PTR_FREE, self.values
PTR_FREE, valuePtr
self.values = PTR_NEW(values)
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errors = *self.ErrorBars*maxValue/SQRT(avgValuesSq)
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(errors)
ENDIF
vmin=GKVsd_Min(values, Max=vmax)
self.vrange=[vmin,vmax]
END  ; ****** PRO GKVsd::Norm ****** ;

FUNCTION GKVsd::Norm, _Extra=extra
result = self -> MakeCopy()
result -> Restrict
result -> NORM, _Extra=extra
RETURN, result
END  ; ****** FUNCTION GKVsd::Norm ****** ;


FUNCTION GKVs3D::LocalNorm, arg, AvgValSq=intensityObj, _Extra=extra
;
; Purpose:
;
;		This function returns a copy of 'self' in which
;		the values have been normalized at each value of 
;		the independent variable identified on the command line
;		such that the local RMS value is 'MaxValue' (defaults to 1).
;		this is useful for removing an overall bias due to variations
;		in the fluctuation level when computing correlation lengths.
;
; Output:
;
;		GKVs3D::LocalNorm returns a GKV object of the same 
;		dimensionality as self containing values which have
;		been scaled at each value of the specified independent
;		variable such that the rms value is locally equal to 
;		rmsValue (defaults to 1.0).
;
;
; Input Arguments:
;
;		GKVs3D::LocalNorm will accept any legal axis identifier
;		(that is, an integer between 1 and 3, or a valid
;		axis mnemonic) as an argument.  If no argument is 
;		provided, then GKVs3D::LocalNorm will expect the axis
;		corresponding to the inhomogeneous coordinate to be
;		identified by keywords (see below). Optional
;
;
; Output Arguments:
;
;		None
;
;
;	Input Keywords:
;
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the selected independent variable. Set axis 
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
;			is interpreted as an index into the grid.values array.
;			(optional)
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable.
;			(optional)
;
;	     skip	The sampling interval for the specified inhomogeneous independent
;			variable.  Defaults to 1 (Optional).
;
;	    debug	Set this keyword (i.e., put '/debug' on the command line) to 
;			print out intermediate information, allowing the user to keep track
;			of the progress (or lack thereof...) of this function. 
;			(optional)
;
;
;	Output Keywords:
;
;			None
;
; Written by W.M. Nevins
;	9/25/00
;
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
				MESSAGE, 'LocalNorm called with too many arguments', /INFORMATIONAL
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
IF(typeOf(iskip) NE 7) THEN skip = FIX(iskip) > 1
;
; Find rmsValue
;
rmsValue=1.
rmsv = GetKeyWord('rmsValue', extra)
IF(typeOf(rmsv) NE 7) THEN rmsValue = ABS(rmsv)
IF(rmsValue EQ 0) THEN rmsValue = 1.
;
; Create GKV object to hold result,
;
result = self -> MakeCopy()
result -> RESTRICT
;
; Create another GKV structure to hold intensities
;
intensityStr = {GKVs1D}
FOR i=0, N_TAGS({GKVsd})-1 DO intensityStr.(i) = Result.(i)
intensityStr.indices = PTR_NEW(["*"])
intensityStr.values = PTR_NEW()
intensityStr.ErrorBars = PTR_NEW()						
;
; Get grid structure corresponding to 
; the inhomogeneous independent variable
;
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'grid = result.Grid' + axisString
ok = EXECUTE(commandString)
IF NOT ok THEN BEGIN
	MESSAGE, "Could not extract inhomogeneous Grid Structure:  ", /INFORMATIONAL
	MESSAGE, "Command String = " + commandString, /INFORMATIONAL
	RETURN, 0
ENDIF	
oldGridValues = *Grid.values
;
; Set up loop over inhomogeneous independent variable
;
irange = Grid.irange
imax = irange[1]-irange[0] 
nPoints = (imax +1)/skip
newGridValues = FLTARR(nPoints)
intensities = FLTARR(npoints)
values = *(result -> GetValues())
info = SIZE(values)
CASE axis OF
	1:	newValues = FLTARR(nPoints, info[2], info[3])
	2:	newValues = FLTARR(info[1], nPoints, info[3])
	3:	newValues = FLTARR(info[1], info[2], nPoints)	
	ELSE:	BEGIN
			messageString =  'Invalid axis identifier:  ' + STRING(axis, FORMAT='(I1)')
			MESSAGE, messageString, /INFORMATIONAL
			RETURN, 0
		END
ENDCASE
CASE axis OF
	1:	nLocal = info[2]*info[3]
	2:	nLocal = info[1]*info[3]
	3:	nLocal = info[1]*info[2]	
ENDCASE
;
; Begin loop over inhomogeneous independent variable
;
FOR i=0, nPoints-1 DO BEGIN				; 'i' is index into 'result' arrays, while 'index' is the  
	index = irange[0]+i*skip			; index to the inhomogeneous independent variable of 'self'. 
	newGridValues[i] = oldGridValues[index]		; Set new grid values
	CASE axis OF
		1:	BEGIN
				localValues=values[index,*,*]
				avg = TOTAL(localValues)/nLocal
				avgSq = TOTAL(localValues^2)/nLocal
				stdSq = avgSq - avg^2
				intensities[i] = stdSq
				newvalues[i,*,*] = (rmsValue/SQRT(stdSq))*(values[index,*,*] - avg)
			END
		2:	BEGIN
				localValues=values[*,index,*]
				avg = TOTAL(localValues)/nLocal
				avgSq = TOTAL(localValues^2)/nLocal
				stdSq = avgSq - avg^2
				intensities[i] = stdSq
				newvalues[*,i,*] = (rmsValue/SQRT(stdSq))*(values[*,index,*] - avg)
			END
		3:	BEGIN
				localValues=values[*,*,index]
				avg = TOTAL(localValues)/nLocal
				avgSq = TOTAL(localValues^2)/nLocal
				stdSq = avgSq - avg^2
				intensities[i] = stdSq
				newvalues[*,*,i] = (rmsValue/SQRT(stdSq))*(values[*,*,index] - avg)
			END
	ENDCASE
ENDFOR
;
; Set fields of output Grid structure
;
PTR_FREE, grid.values
Grid.values = PTR_NEW(newGridValues)
Grid.uniform = GKVsd_UniformGrid(newGridValues)
IF( (Grid.irange[0] NE 0) OR (Grid.irange[1] NE (N_ELEMENTS(oldGridValues)-1)) ) THEN Grid.boundary = 'open'
xmin = MIN(newGridValues, MAX=xmax)
Grid.range  = [xmin,xmax]
Grid.irange = [0,nPoints-1]
;
; Set fields of 'result' structure
;
result.mnemonic = result.mnemonic + '_Normed'
result.title = result.title + '!X!N|!LNormed!N'
result.values = PTR_NEW(newValues)
result.units = ''
vmin = GKVsd_MIN(newValues, Max=vmax)
result.vrange = [vmin, vmax]
CASE axis OF
	1 : result.Grid1 = Grid
	2 : result.Grid2 = Grid
	3 : result.Grid3 = Grid
ENDCASE
;
; Set fields of 'intensity' structure
;
intensityStr.title = "!12<(!X" + intensityStr.title + ")!U2!N!12>!X"
intensityStr.units = "(" + intensityStr.units + ")!U2!N"
intensityStr.values = PTR_NEW(intensities)
vmin = GKVsd_MIN(intensities, Max=vmax)
intensityStr.vrange = [vmin, vmax]
IntensityStr.Grid1 = GKVsd_GridCopy(Grid)
;
; register 'intensityStr' as a GKVs1D object
;
intensityObj = OBJ_NEW('GKVs1D', intensityStr)
;
; and we're done!
;
RETURN, result

END ; ****** GKVs3D::LocalNorm ****** ;



FUNCTION GKVs1D::Slope, error=error
;
; Returns a least-squares-fit to the the (scalar) slope of a
; GKVs1D object
;
; Written by W.M. Nevins
;	5/7/00
;
imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
yvalues = (*self.values)[imin:imax]
xvalues = (*self.Grid1.values)[imin:imax]
Trend_Line, xvalues, yvalues, slope=AA, PMSlope=error
RETURN, AA
END ; ****** GKVs1D::Slope ****** ; 


FUNCTION GKVs1D::FullWidth, Debug=debug, Error=error, Fraction=fract
;
; Purpose:
;
;		Returns the full width at 'fraction' of maximum of the input GKVs1D object's 
;		data.  An estimate of the error in the full width is returned via the "Error"
;		keyword.
;
;  Input KeyWords:
;
;		Fraction	The fraction of the maximum value at which to evaluate the 
;				width. Defaults to 0.5 (thereby returning the full width at
;				hald maximum).  Optional.
;
;		Debug		Set this keyword (i.e., put '/Debug' on the command line) to
;				cause FullWidth to print out several intermediate results.
;
;  Output Keywords:
;
;		Error		On return 'Error' contains an estimate of the error in the halfwidth.
;				(Optional)
;
; Compute full width of 'self'
; using quadratic interpolation
; Interpolation algorithm does NOT assume a uniform grid 
; ... but recall that algorithm for computing 2-time correlation function 
; (Xcorr) does.
;
;
; Written by W.M. Nevins
;	4/15/00
;
; Fixed error in quadratic interpolation at xback
;	W.M. Nevins 2/9/04
;
fraction=0.5
IF KEYWORD_SET(fract) THEN fraction=fract
imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
values=FLOAT((*self.values)[imin:imax])
x=(*self.Grid1.values)[imin:imax]
jmax=N_ELEMENTS(values) - 1
xmin = x[0]
xmax = x[jmax]
maxVal=MAX(values, imaxVal)
epsilon=maxVal*1.e-6
xMaxVal=x[imaxVal]
noiseLevel = SQRT( TOTAL(values*values)/(imax-imin+1.) )
;
; Find grid point nearest half-value for indices greater than jmax
;
temp=FLOAT(values - fraction*maxVal)
FOR jUp=imaxVal,jmax DO BEGIN
	vUp = temp[jUp]
	IF(vUp LE 0) THEN GOTO, GOTjUp
ENDFOR
jUp=jmax
GOTjUp:
IF(jUP EQ jmax) THEN BEGIN
	dxUp = 0.
	dxtotal = 2.*(x[jUp] - x[jUp-1])
	errorUp = 0
	GOTO, DoneUp
ENDIF
IF(ABS(vUp) LT epsilon) THEN BEGIN
	dxUP = 0.
	dxtotal= (x[jUp+1] - x[jUp-1])
	dxminus= x[jUp]   - x[jUp-1]
	dVminusdx= (temp[jUp]   - temp[jUp-1])/dxminus
	errorUp = ABS(noiseLevel/dVminusdx)
	GOTO, DoneUP
ENDIF
;
; Now, compute coefficients for quadratic interpolation about jUp
;
dxplus = x[jUp+1] - x[jUp]
dxminus= x[jUp]   - x[jUp-1]
dxtotal= x[jUp+1] - x[jUp-1]
dVplusdx = (temp[jUp+1] - temp[jUp  ])/dxplus
dVminusdx= (temp[jUp]   - temp[jUp-1])/dxminus
c = temp[jUp]
b = (dxplus*dvminusdx + dxminus*dvplusdx)/dxtotal
a = (dVplusdx - dVminusdx)/dxtotal
discriminant = 1. - 4.*a*c/(b^2)
IF(discriminant GE 0) THEN BEGIN
	dxUp = -2.*(c/b)/(1.0 + SQRT(discriminant) )		; This is the 'nearby' root to the quadratic
ENDIF ELSE BEGIN
	dxUp = -c/dVminusdx					; Use linearinterpolation if discriminant is negative...								
ENDELSE
errorUp = ABS(noiseLevel/dVminusdx)
DoneUP : IF(ABS(dxUp) GT dxtotal) THEN dxup=0.			; avoid doing something stupid...
xUp = (x[jUp] + dxUp) < xmax
; 
; Same for quadratic interpolation about jback
;
FOR jback=imaxVal,0,-1 DO BEGIN
	vBack = temp[jback]
	IF(vBack LE 0) THEN GOTO, GOTjBack
ENDFOR
jback=0
GOTjBack:
IF(jback EQ 0) THEN BEGIN
	dxBack = 0.
	dxtotal= 2.*(x[jBack+1] - x[jBack])
	errorBack = 0
	GOTO, DoneBack
ENDIF
IF(ABS(vBack) LT epsilon) THEN BEGIN
	dxBack = 0.
	dxtotal= (x[jBack+1] - x[jBack-1])
	dxplus = x[jBack+1] - x[jBack]
	dVplusdx = (temp[jBack+1] - temp[jBack  ])/dxplus
	errorBack = ABS(noiseLevel/dVplusdx )
	GOTO, DoneBack
ENDIF
dxplus = x[jBack+1] - x[jBack]
dxminus= x[jBack]   - x[jBack-1]
dxtotal= x[jBack+1] - x[jBack-1]
dVplusdx = (temp[jBack+1] - temp[jBack  ])/dxplus
dVminusdx= (temp[jBack]   - temp[jBack-1])/dxminus
c = temp[jBack]
b = (dxplus*dvminusdx + dxminus*dvplusdx)/dxtotal
a = (dVplusdx - dVminusdx)/dxtotal
discriminant = 1. - 4.*a*c/(b^2)
IF(discriminant GE 0) THEN BEGIN
	dxBack = -2.*(c/b)/( 1. + SQRT(discriminant) )
ENDIF ELSE BEGIN
	dxBack = -c/dVplusdx							; Use linearinterpolation if discriminant is negative...								
ENDELSE
errorBack = ABS(noiseLevel/dVplusdx )
DoneBack : IF(ABS(dxBack) GT dxtotal) THEN dxBack=0.
xBack = (x[jBack] + dxBack) > xmin
fullWidth = ABS(xUp-xBack)
error = SQRT(errorUp^2 + errorBack^2)
IF(KEYWORD_SET(debug)) THEN PRINT, xUp, dxUp, xBack, dxBack, fullwidth, error
RETURN, fullWidth
END ; ****** GKVs1D::FullWidth ****** ;
