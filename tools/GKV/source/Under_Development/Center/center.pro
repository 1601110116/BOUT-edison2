PRO GKVs1D::Center, arg, _Extra=Extra
;
; Purpose:
;
;		This proceedure shifts data in 'self'' such that it
;		is centered about the specified location of the specified
;		axis.
;
; Argument:	Any legal axis identifier.  This specifies the axis over which
;		the data in 'self; is to be shifted.
;
; Keywords:
;
;	axis		Set 'axis' equal to any legal axis identifier to specify the axis
;			over which the data in 'self' is to be shifted.
;
;	index		Set 'index' equal to the index into the specified independent variable 
;			array about which you wish to center the data.
;
;	value		Set 'value' equal to the value of the specified independent varialbe about
;			which you wish to center the data.  Defaults to 0. (Optional)
;
;	'mnemonic'	Set 'mnemnonic' (where 'mnemonic' is a legal axis mnemonic) equal
;			to the value of this independent variable about which you wish to 
;			center the data to simultaneously specify the axis and center value.
;
;	
; Written by W.M. Nevins
;	4/29/02
;
; Set defaults:
;
gridValue = 0.
;
; Find axis identifier
;
nDims = self -> NumDims()
iaxis = -1
IF(N_PARAMS() EQ 1) THEN BEGIN			; Take argument as axis identifier if present
	iaxis = arg
ENDIF ELSE BEGIN				; Else, look for the keyword 'Axis'
	result = GetKeyWord('Axis', extra)
	IF(Query_Integer(result)) THEN iaxis = result
	IF(TypeOf(result) EQ 7) THEN BEGIN
		IF(result NE 'undefined') THEN iaxis = result
	ENDIF
ENDELSE
;
; If 'iaxis' is a mnemonic, then convert it to an integer axis identifier
;
IF(TypeOF(iaxis) EQ 7) THEN iaxis = self -> AxisNumber(iaxis)
	IF(NOT Query_Integer(iaxis)) THEN BEGIN
		MESSAGE, "Illegal axis identifier, unable to center data", /INFORMATIONAL
		RETURN
	ENDIF
;
; Check if we have a valid (integer) axis identifier
;
IF(iaxis NE -1) THEN BEGIN
	IF((iaxis LT 1) OR (iaxis GT nDims)) THEN BEGIN
		MESSAGE, "Illegal axis number, unable to center data", /INFORMATIONAL
		RETURN	
	ENDIF
;
; or (if we haven't found an axis identifer yet)
; look for 'mnemonic = ...'
;
ENDIF ELSE BEGIN
	axisInfo = self -> GetAxis(extra)
	FOR iaxis = 1, nDIms DO BEGIN
		axisValue = axisInfo.(iaxis-1)
		IF(TypeOf(axisValue) NE 7) THEN BEGIN
			gridValue = axisValue
			GOTO, GotAxis
		ENDIF
	ENDFOR
	;
	; Failed to find a valid axis mnemonic 
	;
	MESSAGE, 'No valid axis identifier was supplied, unable to center data', /INFORMATIONAL
	RETURN
ENDELSE
;
; You should only get here if you already have an axis ID
;
GotAxis : axisString = STRING(iaxis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
objValues = *self.values
info = SIZE(objValues)
nPoints = info[iaxis]
IF STRCMP(Grid.boundary, 'periodic (closed)', /FOLD_CASE) THEN nPoints = nPoints-1
gridValues = (*Grid.values)[0:nPoints-1]
;
; Now check command line for 'value':
;
result = GetKeyWord('value', Extra)
IF(TypeOf(result) NE 7) THEN gridValue = result
;
; Compute index of grid point closest to 'value'
;
temp = (gridValues - gridValue)^2
epsilon = MIN(temp, index)
;
; Check command line for index = ...
;
result = GetKeyWord('index', Extra)
IF Query_Integer(result) THEN index=result
;
; Fold 'index' back into first zone, using periodic boundary conditions.
; 
index = index MOD nPoints
IF(index LT 0) then index = index + nPoints
;
; correct 'gridValue' to value of nearest grid point
;
gridValue = gridValues[index]
;
; Now, shift the 'objValues' array:
;
odd = nPoints - 2*(npoints/2)
indexShift = -index + nPoints/2
newValues = MAKE_ARRAY(SIZE=info)
CASE nDims OF
	1:	newValues[0:(nPoints-1)] = SHIFT(objValues[0:(nPoints-1)], indexShift)	
	2:	CASE iaxis OF
			1:	newValues[0:(nPoints-1), *] = SHIFT(objvalues[0:(nPoints-1), *], indexShift, 0)
			2:	newValues[*, 0:(nPoints-1)] = SHIFT(objvalues[*, 0:(nPoints-1)], 0, indexShift)
		ENDCASE
	3:	CASE iaxis OF
			1:	newValues[0:(nPoints-1), *, *] = SHIFT(objvalues[0:(nPoints-1), *, *], indexShift, 0, 0)
			2:	newValues[*, 0:(nPoints-1), *] = SHIFT(objvalues[*, 0:(nPoints-1), *], 0, indexShift, 0)
			3:	newValues[*, *, 0:(nPoints-1)] = SHIFT(objvalues[*, *, 0:(nPoints-1)], 0, 0, indexShift)
		ENDCASE
	4:	CASE iaxis OF
			1:	newValues[0:(nPoints-1), *, *, *] = SHIFT(objvalues[0:(nPoints-1), *, *, *], indexShift, 0, 0, 0)
			2:	newValues[*, 0:(nPoints-1), *, *] = SHIFT(objvalues[*, 0:(nPoints-1), *, *], 0, indexShift, 0, 0)
			3:	newValues[*, *, 0:(nPoints-1), *] = SHIFT(objvalues[*, *, 0:(nPoints-1), *], 0, 0, indexShift, 0)
			4:	newValues[*, *, *, 0:(nPoints-1)] = SHIFT(objvalues[*, *, *, 0:(nPoints-1)], 0, 0, 0, indexShift)
		ENDCASE
ENDCASE
;
; Check for ErrorBars
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	oldErrors = *self.ErrorBars
	errorInfo = SIZE(oldErrors)
	newErrors = MAKE_ARRAY(SIZE=errorInfo)
	CASE nDims OF
		1:	newErrors[0:(nPoints-1)] = SHIFT(oldErrors[0:(nPoints-1)], indexShift)	
		2:	CASE iaxis OF
				1:	newErrors[0:(nPoints-1), *] = SHIFT(oldErrors[0:(nPoints-1), *], indexShift, 0)
				2:	newErrors[*, 0:(nPoints-1)] = SHIFT(oldErrors[*, 0:(nPoints-1)], 0, indexShift)
			ENDCASE
		3:	CASE iaxis OF
				1:	newErrors[0:(nPoints-1), *, *] = SHIFT(oldErrors[0:(nPoints-1), *, *], indexShift, 0, 0)
				2:	newErrors[*, 0:(nPoints-1), *] = SHIFT(oldErrors[*, 0:(nPoints-1), *], 0, indexShift, 0)
				3:	newErrors[*, *, 0:(nPoints-1)] = SHIFT(oldErrors[*, *, 0:(nPoints-1)], 0, 0, indexShift)
			ENDCASE
		4:	CASE iaxis OF
				1:	newErrors[0:(nPoints-1), *, *, *] = SHIFT(oldErrors[0:(nPoints-1), *, *, *], indexShift, 0, 0, 0)
				2:	newErrors[*, 0:(nPoints-1), *, *] = SHIFT(oldErrors[*, 0:(nPoints-1), *, *], 0, indexShift, 0, 0)
				3:	newErrors[*, *, 0:(nPoints-1), *] = SHIFT(oldErrors[*, *, 0:(nPoints-1), *], 0, 0, indexShift, 0)
				4:	newErrors[*, *, *, 0:(nPoints-1)] = SHIFT(oldErrors[*, *, *, 0:(nPoints-1)], 0, 0, 0, indexShift)
			ENDCASE
	ENDCASE
	IF STRCMP(Grid.boundary, 'periodic (closed)', /FOLD_CASE) THEN BEGIN
		CASE nDims OF
			1:	newErrors[nPoints] = newErrors[0]
			2:	CASE iaxis OF
					1:	newErrors[nPoints, *] = newErrors[0, *]
					2:	newErrors[*, nPoints] = newErrors[*, 0]
				ENDCASE
			3:	CASE iaxis OF
					1:	newErrors[nPoints, *, *] = newErrors[0, *, *]
					2:	newErrors[*, nPoints, *] = newErrors[*, 0, *]
					3:	newErrors[*, *, nPoints] = newErrors[*, *, 0]
				ENDCASE
			4:	CASE iaxis OF
					1:	newErrors[nPoints, *, *, *] = newErrors[0, *, *, *]
					2:	newErrors[*, nPoints, *, *] = newErrors[*, 0, *, *]
					3:	newErrors[*, *, nPoints, *] = newErrors[*, *, 0, *]
					4:	newErrors[*, *, *, nPoints] = newErrors[*, *, *, 0]
				ENDCASE
			ENDCASE
	ENDIF
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(newErrors)
ENDIF
;
; Shift gridValues array
;
lx = (gridValues[nPoints-1] - gridValues[0])
dx = (gridValues[1] - gridValues[0])
lx = lx + dx
gridValues = Shift(gridValues, indexShift) - gridValue
;
; Put 'indexShift' into first zone
;
indexShift = indexShift MOD nPoints
IF(indexShift LT 0) THEN indexShift = indexShift + nPoints
;
; Remove dis-continuities in gridvalues if necessary
;
IF(indexShift GT 0) THEN BEGIN
	IF(gridValues[(indexShift)-1] GT  lx/2.) THEN gridValues[0:(indexShift-1)] = gridValues[0:(indexShift-1)] - lx
	IF(gridValues[(indexShift)  ] LT -lx/2.) THEN gridValues[(indexShift):(npoints-1)] = gridValues[(indexShift):(npoints-1)] + lx
ENDIF
irange = [0,nPoints-1]
IF STRCMP(Grid.boundary, 'periodic (closed)', /FOLD_CASE) THEN BEGIN
	gridValues = [gridValues, gridValues[0] + lx]
	irange = [0, nPoints]
	CASE nDims OF
		1:	newValues[nPoints] = newValues[0]
		2:	CASE iaxis OF
				1:	newValues[nPoints, *] = newValues[0, *]
				2:	newValues[*, nPoints] = newValues[*, 0]
			ENDCASE
		3:	CASE iaxis OF
				1:	newValues[nPoints, *, *] = newValues[0, *, *]
				2:	newValues[*, nPoints, *] = newValues[*, 0, *]
				3:	newValues[*, *, nPoints] = newValues[*, *, 0]
			ENDCASE
		4:	CASE iaxis OF
				1:	newValues[nPoints, *, *, *] = newValues[0, *, *, *]
				2:	newValues[*, nPoints, *, *] = newValues[*, 0, *, *]
				3:	newValues[*, *, nPoints, *] = newValues[*, *, 0, *]
				4:	newValues[*, *, *, nPoints] = newValues[*, *, *, 0]
			ENDCASE
	ENDCASE
ENDIF
;
; Construct new Grid structure
;
PTR_FREE, Grid.values
Grid.Values = PTR_NEW(gridValues)
IF( NOT STRCMP(grid.title, '!4D!X', 1) ) THEN grid.title    = '!4D!X' + grid.title
IF( NOT STRCMP(grid.mnemonic, 'd' , 1) ) THEN grid.mnemonic = 'd' + grid.mnemonic
grid.irange = irange
xmin = MIN(gridValues, MAX=xmax)
grid.range = [xmin, xmax]
;
; Replace grid structure
;
commandString = 'self.Grid' + axisString + ' = Grid'
ok = EXECUTE(commandString)

;
; Replace values in 'self'
;
PTR_FREE, self.values
self.Values = PTR_NEW(newValues)
RETURN
END	; ******  GKVs1D::Center  ****** ;
	 