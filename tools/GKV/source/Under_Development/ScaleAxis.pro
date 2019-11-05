PRO GKVs1D::ScaleAxis, arg, _extra=extra
;
; Purpose:
;
;	This proceedure rescales selected grid, allowing user
;	to change units.  If no valid scaling information is 
;	supplied, then the grid values are left unchanged, but
;	grid units, grid mnemonic, and grid title may still be 
;	changed. 
;
;
; Arguments:
;
;	arg 		(optional). If present, this argument is a legal axis identifier, either:
;			an integer between 1 and the dimensionality of this object, or
;			a string containing an axis mnemonic
;
; Keywords:
;
;	Axis		Set equal to a legal axis identifier, either:
;			an integer between 1 and the dimensionality of this object, or
;			a string containing an axis mnemonic
;			Optional
;
;	Const		Constant multiplier for rescaling axis.
;			Optional
;
;	Uniform		Set this keyword (i.e., put "/uniform" on the command line)
;			to coherced the grid spacing for the specified grid to be uniform.
;			(Optional)
;
;	OffSet		Offset for rescaled axis (relative to original axis)
;			Optional
;
;	Delta		Constant separation for uniformly spaced rescaled axis
;			(only valid if original axis was uniformly scaled).  ; or
;			Optional
;
;    'Mnemonic'		A runtime keyword consisting of the mnemonic for one of the
;			object's axis.  Set equal to an array containing the values
;			you desire at each grid point for this independent variable.
;			(optional). 
;
;   'dMnemonic'		A runtime keyword consisting of 'd' concatenated with
;			the mnemonic of the selected axis.  Has the effect as 
;			the keyword 'delta' above.
;			Optional
;
;	Start		Initial value for rescaled uniform grid formed with either
;			the keyword 'delta', or the runtime keyword 'dMnemonic'.
;			Optional
;
; 'Mnemonic_0'		A runtime keyword consisting of the mnemonic of the selected
;			axis concatnated with '_0'.  Has an effect identical to the
;			keyword 'start' above.
;			Optional
;
;	Units		A string constant containing the units of rescaled grid
;			(use Hershey Vector Font for Greek or special characters)
;			Optional
;
;	Title		A string constant containing the title of the rescaled grid
;			(use Hershey Vector Font for Greek or special characters)
;			Optional
;
;	Mnemonic	A string constant containing the 'mnemonic' of the rescaled grid.
;			The Mnemonic should contain no embedded blanks, or any character which
;			IDL will interpret as an operator.
;
; Written by W.M. Nevins
;	7/13/00
;
nDims = self -> NumDims()
IF(N_ELEMENTS(arg) NE 0) THEN BEGIN		; Get axis ID from 'arg' if one is present
	iaxis = arg					;	if 'arg' is a string, try to match with a mnemonic
	IF(TypeOF(arg) EQ 7) THEN iaxis = self -> AxisNumber(arg)
								; Otherwise, get axis ID from keyword "Axis"
	IF(NOT Query_Integer(iaxis)) THEN GOTO, tryagain	
ENDIF ELSE BEGIN
tryagain:
	iaxis=0
	result = GetKeyWord('axis', extra)
	IF(Query_Integer(result)) THEN iaxis = result				; command line of form axis = axisnumber
	IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result)	; command line of the form axis = 'mnemnic'
ENDELSE
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN
	ENDIF
	axisInfo = self -> GetAxis(extra)
	FOR i=0, nDims-1 DO BEGIN
		iaxis = i+1
		IF(typeOf(axisInfo.(i)) NE 7) THEN BEGIN
			gridValues = axisInfo.(i)
			GOTO, DONE1
		ENDIF
	ENDFOR
;
; an axis mnemonic did not appear on the command line, so print error message and return
;
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN
ENDIF
DONE1:
;
; Get grid structure corresponding to selected axis
;
axisStr = STRING(iaxis, FORMAT='(i1)')
commandStr = 'grid = self.Grid' + axisStr
ok = EXECUTE(commandStr)
;
; If values have been provided via 'mnemonic' = [...], then use these
;
IF(N_ELEMENTS(gridValues) NE 0) THEN BEGIN
	IF( N_ELEMENTS(gridValues) NE N_ELEMENTS(*grid.values) ) THEN BEGIN
		PRINT, 	grid.mnemonic, "-axis has ", N_ELEMENTS(*grid.values), "grid points, while ", $
				N_ELEMENTS(gridValues), " were supplied."
		RETURN
	ENDIF
	gridValues = FLOAT(gridValues)
	PTR_FREE, grid.values
	grid.values = PTR_NEW(gridValues)
	gridMin = MIN(gridValues, max=gridMax)
	Grid.range = [gridMin, gridMax]
	GOTO, setUniform
ENDIF
;
; Get information for scaling selected axis
;	first look for 'const=xx' on command line
;
result = GetKeyWord('const', extra)
IF( TypeOf(result) NE 7) THEN BEGIN
	const = result
	;
	; Check for 'OffSet' keyword
	;
	offSet = 0.
	result = GetKeyWord('OffSet', extra)
	IF( TypeOf(result) NE 7) THEN offSet = result
	gridValues = *Grid.values
	gridValues = offSet + const*gridValues
	PTR_FREE, Grid.values
	Grid.values = PTR_NEW(gridValues)
	gridMin = MIN(gridValues, max=gridMax)
	Grid.range = [gridMin, gridMax]
	GOTO, setUniform
ENDIF
IF(NOT Grid.uniform) THEN GOTO, setUniform
;
; Now check for either 'delta' or 'dMnemonic'
;
dMnemonic = 'd' + Grid.mnemonic
result = GetKeyWord(dMnemonic, extra)
IF(TypeOF(result) EQ 7) THEN result = GetKeyWord('delta', extra)
IF(TypeOF(result) EQ 7) THEN BEGIN
	MESSAGE, 'No valid scaling information supplied', /INFORMATIONAL
	GOTO, setUniform
ENDIF
delta = result
;
; Check for either 'start' or 'Mnemonic_0'
;
start=0.
Mnemonic_0 = Grid.mnemonic + '_0'
result = GetKeyWord(Mnemonic_0, extra)
IF(TypeOF(result) EQ 7) THEN result = GetkeyWord('start', extra)
IF(TypeOF(result) NE 7) THEN start = result
gridValues = *Grid.values
npoints = N_ELEMENTS(gridValues)
PTR_FREE, Grid.values
gridValues = start + delta*FINDGEN(npoints)
Grid.values = PTR_NEW(gridValues)
gridMin = MIN(gridValues, MAX=gridMax)
Grid.range = [gridMin, gridMax]
GOTO, setUniform

setUniform	:
;
; Before setting units, etc. check if selected grid needs to be made uniform
;
result = GetKeyWord('Uniform', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
IF(result NE 1) THEN GOTO, SetUnits
;
; First reset grid so that "/Uniform" will work with other commands
;
	commandStr = 'self.Grid' + axisStr + ' = Grid'
	ok = EXECUTE(commandStr)
;
; Now get current grid info
;
	oldGridValuesPtr = grid.Values
	gridValues = *grid.values
	irange = Grid.irange
	imin = irange[0]
	imax = irange[1]
	dValues = gridValues[(imin+1):imax] - gridValues[imin:(imax-1)]
	dMin = MIN(dValues)
	vMin = gridValues[imin]
	vMax = gridValues[imax]
;
; Compute values of the independpendent variable on an appropriate uniform grid
;
	nPoints = ROUND((vMax - vMin)/dMin + 0.4999999) + 1
	dvalue = (vmax - vmin)/(nPoints - 1)
	newValues = vmin + dValue*FINDGEN(nPoints)
;
; Reset values, etc. in the structure "Grid" 
;
	grid.Values = PTR_NEW(newValues)
	grid.irange =[0, npoints-1] 
	grid.uniform = 1b
;
; Make dummy object 
;
	target = self -> MakeCopy(/NoValues, /NoErrorBars)
;
; Replace appropriate grid with uniform grid just created
;
	targetGrid = GKVsd_Gridcopy(grid)
	commandStr = 'target.Grid' + axisStr + ' = targetGrid'
	ok = EXECUTE(commandStr)
;
; Interpolate 'self' onto new grid
;
	uniformSelf = self -> Interpolate(target)
;
; Reload interpolated data into 'self'
;
	PTR_FREE, self.values
	self.values = uniformSelf.values
;
; Now, load new grid structure into 'self'
;
	PTR_FREE, oldGridValuesPtr
	commandStr = 'self.Grid' + axisStr + ' = Grid'
	ok = EXECUTE(commandStr)
;
; and clean up 
;
	target -> trash
	uniformSelf.values = PTR_NEW()
	uniformSelf -> trash
ENDIF
;
;  Now Check for new units, etc.
;
setUnits	:
result = GetKeyWord('units', extra)
IF(TypeOF(result) EQ 7) THEN BEGIN
	IF(result NE 'undefined')  THEN Grid.units = result
ENDIF
;
; Check for new 'title'
;
result = GetKeyWord('title', extra)
IF(TypeOF(result) EQ 7) THEN BEGIN
	IF(result NE 'undefined') THEN Grid.title = result
ENDIF
;
; Check for new 'mnemonic'
;
result = GetKeyWord('mnemonic', extra)
IF(TypeOF(result) EQ 7) THEN BEGIN
	IF(result NE 'undefined') THEN Grid.mnemonic = STRCOMPRESS(result, /REMOVE_ALL)
ENDIF
;
; Now replace rescaled grid structure
;
commandStr = 'self.Grid' + axisStr + ' = Grid'
ok = EXECUTE(commandStr)

RETURN
END ; ****** GKVs1D::ScaleAxis ****** ;