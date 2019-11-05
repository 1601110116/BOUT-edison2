Pro GKVs1D::Smooth, _EXTRA=Extra
;
; Apply '1:2:1' smoothing to 1-D data
;
temp = *self.values
values = 0.25*(SHIFT(temp, -1) + 2.*temp + SHIFT(temp, 1))
self.values = PTR_NEW(values)
self -> Set, _EXTRA=Extra
RETURN
END ; ****** GKVs1D::Smooth ****** ;



Function GKVs1D::Avg,arg, _Extra=extra
;
; Purpose:
;
;		This function returns the average of 'self' over the selected
;		independent variable.
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable to be averaged 
;			over. Set 'axis' equal to any legal axis identifier 
;			(see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be averaged over, and to set the range in this variable
;			over which this average is to be taken.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the average over the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the averate 
;			is to be taken.
;
;	IgnoreErrors	'Setting' this keyword causes 'Avg' to ignore the 'ErrorBars' in
;			'self' when estimating the error in the average. 
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
; Written by W.M. Nevins
;	6/18/00
;
; Revised by W.M. Nevins
;	3/20/01
;
; Further revised by W.M. Nevins
;	4/30/02
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
resultStr = {GKVsd}
FOR i=0, N_TAGS(resultStr)-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title

resultStr.Title = '!12<!X' + self.title + '!12>!X!D' + axisTitle  + '!N'
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Compute Average and standard deviation
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
lx = ( (*grid.values)[irange[1]] - (*grid.values)[irange[0]] )
values = 0.5*((*self.values)[(irange[0]+1):irange[1]] + (*self.values)[irange[0]:(irange[1]-1)]) 
avgValue   = TOTAL(values*dx)/lx
avgValueSq = TOTAL(ABS(values)^2*dx)/lx
var = avgValueSq - ABs(avgValue)^2
;
; Set up array of errors
;
result = GetKeyWord('IgnoreErrors', Extra)
	IF Query_Integer(IgnoreErrors) THEN BEGIN
		errorBars = SQRT(var)
	ENDIF ELSE BEGIN
	IF PTR_VALID(self.ErrorBars) THEN BEGIN
		errors = 0.5*((*self.ErrorBars)[(irange[0]+1):irange[1]] + (*self.ErrorBars)[irange[0]:(irange[1]-1)])
		avgErrorSq = TOTAL(errors*errors*dx)/lx
		errorBars = SQRT(avgErrorSq + var)
	ENDIF ELSE BEGIN
		errorBars = SQRT(var)
	ENDELSE
ENDELSE
;
; Load results into 'result' structure
;
resultStr.values = PTR_NEW(avgValue)
resultStr.vrange = [0, 0]
resultStr.ErrorBars = PTR_NEW(errorBars)
result = OBJ_NEW("GKVsd", resultStr)
IF Query_Integer(Extra) THEN RETURN, result
result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs1D::Avg ****** ;


Function GKVs2D::Avg, arg, _Extra=extra
;
; Purpose:
;
;		This function returns the average of 'self' over the selected
;		independent variable.
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable to be averaged 
;			over. Set 'axis' equal to any legal axis identifier 
;			(see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be averaged over, and to set the range in this variable
;			over which this average is to be taken.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the average over the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the averate 
;			is to be taken.
;
;	IgnoreErrors	'Setting' this keyword causes 'Avg' to ignore the 'ErrorBars' in
;			'self' when estimating the error in the average. 
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
; Written by W.M. Nevins
;	6/18/00
;
; Revised by W.M. Nevins
;	3/20/01
;
; Further revised by W.M. Nevins
;	4/30/02
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
resultStr = {GKVs1D}

FOR i=0, N_TAGS({GKVsd})-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
irange=irange - irange[0]
axisTitle = grid.title

resultStr.Title = '!12<!X' + self.title + '!12>!X!D' + axisTitle  + '!N'
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create 2-D 'dx' array
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
CASE axis OF
	1:	BEGIN
			ones = MAKE_ARRAY(info[2], value=1)
			dx = dx#ones
			values = 0.5*( (*valuePtr)[(irange[0]+1):irange[1], *] + (*valuePtr)[irange[0]:(irange[1]-1), *]) 
		END
	2:	BEGIN
			ones = MAKE_ARRAY(info[1], value=1)
			dx = ones#dx
			values = 0.5*( (*valuePtr)[*, (irange[0]+1):irange[1]] + (*valuePtr)[*, irange[0]:(irange[1]-1)])
		END
ENDCASE
lx = ( (*grid.values)[irange[1]] - (*grid.values)[irange[0]] )
;
; Compute Average and standard deviation
;
avgValue   = TOTAL(values*dx,        axis)/lx
avgValueSq = TOTAL(ABS(values)^2*dx, axis)/lx
var = avgValueSq - ABS(avgValue)^2
;
; Set up array of errors
;
avgErrorSq = 0.
result = GetKeyWord('IgnoreErrors', Extra)
IF Query_Integer(IgnoreErrors) THEN BEGIN
	errorBars = SQRT(var)
ENDIF ELSE BEGIN
	IF PTR_VALID(self.ErrorBars) THEN BEGIN
		errorPtr = self -> GetErrors()
		CASE axis OF
			1:	errors = 0.5*((*errorPtr)[(irange[0]+1):irange[1], *] + (*errorPtr)[irange[0]:(irange[1]-1), *])
			2:	errors = 0.5*((*errorPtr)[*, (irange[0]+1):irange[1]] + (*errorPtr)[*, irange[0]:(irange[1]-1)])
		ENDCASE
		PTR_FREE, errorPtr
		avgErrorSq = TOTAL(errors*errors*dx, axis)/lx
	ENDIF
	errorBars = SQRT(var + avgErrorSq)
ENDELSE
;
; Get output grid structure
;
CASE axis OF
	1:	GridOut = GKVsd_GridCopy(self.grid2)
	2:	GridOut = GKVsd_GridCopy(self.grid1)
ENDCASE
;
; Load results into 'result' structure
;
resultStr.values = PTR_NEW(avgValue)
resultStr.vrange = [0, 0]
resultStr.ErrorBars = PTR_NEW(errorBars)
resultStr.Grid1 = GridOut
;
; Clean up
;
PTR_FREE, ValuePTR

result = OBJ_NEW("GKVs1D", resultStr)
result -> GridRestrict, 1

;IF Query_Integer(Extra) THEN RETURN, result
;result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs2D::Avg ****** ;


Function GKVs3D::Avg,arg, _Extra=extra
;
; Purpose:
;
;		This function returns the average of 'self' over the selected
;		independent variable.
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable to be averaged 
;			over. Set 'axis' equal to any legal axis identifier 
;			(see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be averaged over, and to set the range in this variable
;			over which this average is to be taken.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the average over the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the averate 
;			is to be taken.
;
;	IgnoreErrors	'Setting' this keyword causes 'Avg' to ignore the 'ErrorBars' in
;			'self' when estimating the error in the average. 
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
; Written by W.M. Nevins
;	6/18/00
;
; Revised by W.M. Nevins
;	3/20/01
;
; Further revised by W.M. Nevins
;	4/30/02
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'AVG called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
resultStr = {GKVs2D}

FOR i=0, N_TAGS({GKVsd})-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title

resultStr.Title = '!12<!X' + self.title + '!12>!X!D' + axisTitle  + '!N'
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create dx-weighted values, valuesSQ arrays
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values   = MAKE_ARRAY(SIZE=Info)
valuesSq = MAKE_ARRAY(SIZE=Info) 
CASE axis OF
	1:	BEGIN
			FOR i=0, info[1]-2 DO   values[i, *, *] =     0.5*( (*valuePtr)[i, *, *] + (*valuePtr)[i+1, *, *])*dx[i]
			FOR i=0, info[1]-2 DO valuesSq[i, *, *] = ABS(0.5*( (*valuePtr)[i, *, *] + (*valuePtr)[i+1, *, *]))^2*dx[i]
		END
	2:	BEGIN
			FOR i=0, info[2]-2 DO   values[*, i, *] =     0.5*( (*valuePtr)[*, i, *] + (*valuePtr)[*, i+1, *])*dx[i]
			FOR i=0, info[2]-2 DO valuesSq[*, i, *] = ABS(0.5*( (*valuePtr)[*, i, *] + (*valuePtr)[*, i+1, *]))^2*dx[i]
		END
	3:	BEGIN
			FOR i=0, info[3]-2 DO   values[*, *, i] =     0.5*( (*valuePtr)[*, *, i] + (*valuePtr)[*, *, i+1])*dx[i]
			FOR i=0, info[3]-2 DO valuesSq[*, *, i] = ABS(0.5*( (*valuePtr)[*, *, i] + (*valuePtr)[*, *, i+1]))^2*dx[i]
		END
ENDCASE
lx = ( (*grid.values)[irange[1]] - (*grid.values)[irange[0]] )
;
; Compute Average and standard deviation
;
avgValue   = TOTAL(values,   axis)/lx
avgValueSq = TOTAL(valuesSq, axis)/lx
var =avgValueSq - ABS(avgValue)^2
;
; Set up array of errors
;
avgErrorSq = 0.
result = GetKeyWord('IgnoreErrors', Extra)
IF Query_Integer(IgnoreErrors) THEN BEGIN
	errorBars = SQRT(var)
ENDIF ELSE BEGIN
	IF PTR_VALID(self.ErrorBars) THEN BEGIN
		errorPtr = self -> GetErrors()
		ErrorSq = MAKE_ARRAY(SIZE=Info) 
			CASE axis OF
			1:	FOR i=0, info[1]-2 DO   ErrorSq[i, *, *] = ( 0.5*(*errorPtr)[i, *, *] + (*errorPtr)[i+1, *, *] )^2*dx[i]
			2:	FOR i=0, info[2]-2 DO   ErrorSq[*, i, *] = ( 0.5*(*errorPtr)[*, i, *] + (*errorPtr)[*, i+1, *] )^2*dx[i]
			3:	FOR i=0, info[3]-2 DO   ErrorSq[*, *, i] = ( 0.5*(*errorPtr)[*, *, i] + (*errorPtr)[*, *, i+1] )^2*dx[i]
		ENDCASE
		avgErrorSq = TOTAL(errorSq, axis)/lx
		PTR_FREE, errorPtr
	ENDIF
	errorBars = SQRT(var + avgErrorSq)
ENDELSE
;
; Get output grid structures
;
CASE axis OF
	1:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid2)
			Grid2 = GKVsd_GridCopy(self.grid3)
		END
	2:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid1)
			Grid2 = GKVsd_GridCopy(self.grid3)
		END
	3:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid1)
			Grid2 = GKVsd_GridCopy(self.grid2)
		END
ENDCASE
;
; Load results into 'result' structure
;
resultStr.values = PTR_NEW(avgValue)
resultStr.vrange = [0, 0]
resultStr.ErrorBars = PTR_NEW(errorBars)
resultStr.Grid1 = Grid1
resultStr.Grid2 = Grid2
;
; Clean up
;
PTR_FREE, ValuePTR

result = OBJ_NEW("GKVs2D", resultStr)
result -> GridRestrict, 1
result -> GridRestrict, 2

;IF Query_Integer(Extra) THEN RETURN, result
;result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs3D::Avg ****** ;



Function GKVs4D::Avg,arg, _Extra=extra
;
; Purpose:
;
;		This function returns the average of 'self' over the selected
;		independent variable.
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable to be averaged 
;			over. Set 'axis' equal to any legal axis identifier 
;			(see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be averaged over, and to set the range in this variable
;			over which this average is to be taken.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the average over the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the averate 
;			is to be taken.
;
;	IgnoreErrors	'Setting' this keyword causes 'Avg' to ignore the 'ErrorBars' in
;			'self' when estimating the error in the average. 
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
; Written by W.M. Nevins
;	6/18/00
;
; Revised by W.M. Nevins
;	3/20/01
;
; Further revised by W.M. Nevins
;	4/30/02
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
resultStr = {GKVs3D}

FOR i=0, N_TAGS({GKVsd})-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title

resultStr.Title = '!12<!X' + self.title + '!12>!X!D' + axisTitle  + '!N'
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create dx-weighted values, valuesSQ arrays
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values   = MAKE_ARRAY(SIZE=Info)
valuesSq = FLTARR(Info[1], info[2], info[3], info[4]) 
CASE axis OF
	1:	BEGIN
			FOR i=1, info[1]-2 DO   values[i, *, *, *] =     0.5*( (*valuePtr)[i, *, *, *] + (*valuePtr)[i+1, *, *, *])*dx[i]
			FOR i=1, info[1]-2 DO valuesSq[i, *, *, *] = ABS(0.5*( (*valuePtr)[i, *, *, *] + (*valuePtr)[i+1, *, *, *]))^2*dx[i]
		END
	2:	BEGIN
			FOR i=1, info[2]-2 DO   values[*, i, *, *] =     0.5*( (*valuePtr)[*, i, *, *] + (*valuePtr)[*, i+1, *, *])*dx[i]
			FOR i=1, info[2]-2 DO valuesSq[*, i, *, *] = ABS(0.5*( (*valuePtr)[*, i, *, *] + (*valuePtr)[*, i+1, *, *]))^2*dx[i]
		END
	3:	BEGIN
			FOR i=1, info[3]-2 DO   values[*, *, i, *] =     0.5*( (*valuePtr)[*, *, i, *] + (*valuePtr)[*, *, i+1, *])*dx[i]
			FOR i=1, info[3]-2 DO valuesSq[*, *, i, *] = ABS(0.5*( (*valuePtr)[*, *, i, *] + (*valuePtr)[*, *, i+1, *]))^2*dx[i]
		END
	4:	BEGIN
			FOR i=1, info[4]-2 DO   values[*, *, *, i] =     0.5*( (*valuePtr)[*, *, *, i] + (*valuePtr)[*, *, *, i+1])*dx[i]
			FOR i=1, info[4]-2 DO valuesSq[*, *, *, i] = ABS(0.5*( (*valuePtr)[*, *, *, i] + (*valuePtr)[*, *, *, i+1]))^2*dx[i]
		END
ENDCASE
lx = ( (*grid.values)[irange[1]] - (*grid.values)[irange[0]] )
;
; Compute Average and standard deviation
;
avgValue   = TOTAL(values,   axis)/lx
avgValueSq = TOTAL(valuesSq, axis)/lx
var = avgValueSq - ABS(avgValue)^2
;
; Set up array of errors
;
avgErrorSq = 0.
result = GetKeyWord('IgnoreErrors', Extra)
IF Query_Integer(IgnoreErrors) THEN BEGIN
	errorBars = SQRT(var)
ENDIF ELSE BEGIN
	IF PTR_VALID(self.ErrorBars) THEN BEGIN
		errorPtr = self -> GetErrors()
		ErrorSq = MAKE_ARRAY(SIZE=Info) 
		CASE axis OF
			1:	FOR i=1, info[1]-1 DO   ErrorSq[i, *, *, *] = ( 0.5*(*errorPtr)[i, *, *, *] + (*errorPtr)[i+1, *, *, *] )^2*dx[i]
			2:	FOR i=1, info[2]-1 DO   ErrorSq[*, i, *, *] = ( 0.5*(*errorPtr)[*, i, *, *] + (*errorPtr)[*, i+1, *, *] )^2*dx[i]
			3:	FOR i=1, info[3]-1 DO   ErrorSq[*, *, i, *] = ( 0.5*(*errorPtr)[*, *, i, *] + (*errorPtr)[*, *, i+1, *] )^2*dx[i]
			4:	FOR i=1, info[4]-1 DO   ErrorSq[*, *, *, i] = ( 0.5*(*errorPtr)[*, *, *, i] + (*errorPtr)[*, *, *, i+1] )^2*dx[i]
		ENDCASE
		avgErrorSq = TOTAL(errorSq, axis)/lx
		PTR_FREE, errorPtr
	ENDIF
	errorBars = SQRT(var + avgErrorSq)
ENDELSE

;
; Get output grid structures
;
CASE axis OF
	1:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid2)
			Grid2 = GKVsd_GridCopy(self.grid3)
			Grid3 = GKVsd_GridCopy(self.grid4)
		END
	2:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid1)
			Grid2 = GKVsd_GridCopy(self.grid3)
			Grid3 = GKVsd_GridCopy(self.grid4)
		END
	3:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid1)
			Grid2 = GKVsd_GridCopy(self.grid2)
			Grid3 = GKVsd_GridCopy(self.grid4)
		END
	4:	BEGIN
			Grid1 = GKVsd_GridCopy(self.grid1)
			Grid2 = GKVsd_GridCopy(self.grid2)
			Grid3 = GKVsd_GridCopy(self.grid3)
		END
ENDCASE
;
; Load results into 'result' structure
;
resultStr.values = PTR_NEW(avgValue)
resultStr.vrange = [0, 0]
resultStr.ErrorBars = PTR_NEW(errorBars)
resultStr.Grid1 = Grid1
resultStr.Grid2 = Grid2
resultStr.Grid3 = Grid3
;
; Clean up
;
PTR_FREE, ValuePTR

result = OBJ_NEW("GKVs3D", resultStr)
result -> GridRestrict, 1
result -> GridRestrict, 2
result -> GridRestrict, 3

;IF Query_Integer(Extra) THEN RETURN, result
;result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs4D::Avg ****** ;
