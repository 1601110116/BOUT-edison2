Function GKVs1D::Int, arg, _Extra=extra
;
; Purpose:
;
;		This function returns the integral of 'self' over the selected
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
;			used to identify independent variable to be integrated 
;			over. Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be integrated over, and to set the range of integration.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to to set the
;			range of integration through indices into the corresponding
;			grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of integration.
;
;	      sum	Set this keyword (i.e., put "/Sum" on the command line)
;			to perform a sum (rather than integral) over the selected
;			range.  Default is sum=0 (that is, to perform an integral).
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
				MESSAGE, 'INT called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF

Sum = 0b
result = GetKeyWord("Sum", Extra)
IF( Query_Integer(result)) THEN Sum=result

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
axisUnits = grid.units

resultStr.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
resultStr.units = '(' + self.units + ')*(' + axisUnits +')'
IF(Sum) THEN BEGIN
	resultStr.Title = '!7R!X!D' + axisTitle  + '!N' + self.title 
	resultStr.Units = self.units
ENDIF 
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Compute integral 
;
IF(Sum) THEN BEGIN
	dx = 1.+ 0.*( (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)] )
ENDIF ELSE BEGIN
	dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
ENDELSE
values = 0.5*((*self.values)[(irange[0]+1):irange[1]] + (*self.values)[irange[0]:(irange[1]-1)]) 
intValue   = TOTAL(values*dx)
;
; Set up array of errors
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errors = 0.5*((*self.ErrorBars)[(irange[0]+1):irange[1]] + (*self.ErrorBars)[irange[0]:(irange[1]-1)])
	avgErrorSq = TOTAL(errors*errors*dx)/lx
	errorBars = SQRT(avgErrorSq)
	resultStr.ErrorBars = PTR_NEW(errorBars)
ENDIF
;
; Load results into 'result' structure
;
resultStr.values = PTR_NEW(intValue)
resultStr.vrange = [0, 0]
result = OBJ_NEW("GKVsd", resultStr)

IF(TypeOF(Extra) EQ 8) THEN	$
	result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs1D::INT ****** ;


Function GKVs2D::INT, arg, _Extra=extra
;
; Purpose:
;
;		This function returns the integral of 'self' over the selected
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
;			used to identify independent variable to be integrated 
;			over. Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be integrated over, and to set the range of integration.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to to set the
;			range of integration through indices into the corresponding
;			grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of integration.
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
				MESSAGE, 'INT called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF

Sum = 0b
result = GetKeyWord("Sum", Extra)
IF( Query_Integer(result)) THEN Sum=result

resultStr = {GKVs1D}

FOR i=0, N_TAGS({GKVsd})-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title
axisUnits = grid.units

resultStr.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
resultStr.units = '(' + self.units + ')*(' + axisUnits +')'
IF(Sum) THEN BEGIN
	resultStr.Title = '!7R!X!D' + axisTitle  + '!N' + self.title 
	resultStr.Units = self.units
ENDIF 
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create 2-D 'dx' array
;
;
IF(Sum) THEN BEGIN
	dx = 1.+ 0.*( (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)] )
ENDIF ELSE BEGIN
	dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
ENDELSE

CASE axis OF
	1:	BEGIN
			ones = MAKE_ARRAY(info[2], value=1)
			dx = dx#ones
			values = 0.5*( (*valuePtr)[1:(irange[1]-irange[0]), *] + (*valuePtr)[0:(irange[1]-irange[0]-1), *]) 
		END
	2:	BEGIN
			ones = MAKE_ARRAY(info[1], value=1)
			dx = ones#dx
			values = 0.5*( (*valuePtr)[*, 1:(irange[1]-irange[0])] + (*valuePtr)[*, 0:(irange[1]-irange[0]-1)])
		END
ENDCASE
;
; Compute Average and standard deviation
;
intValue   = TOTAL(values*dx, axis)
;
; Set up array of errors
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errorPtr = self -> GetErrors()
	CASE axis OF
		1:	errors = 0.5*( (*errorPtr)[1:(irange[1]-irange[0]),  *] + (*errorPtr)[0:(irange[1]-irange[0]-1), *]  )
		2:	errors = 0.5*( (*errorPtr)[*, 1:(irange[1]-irange[0]) ] + (*errorPtr)[*, 0:(irange[1]-irange[0]-1)]  )
	ENDCASE
	PTR_FREE, errorPtr
	avgErrorSq = TOTAL(errors*errors*dx, axis)
	errorBars = SQRT(avgErrorSq)
	PTR_FREE, resultStr.ErrorBars
	resultStr.ErrorBars = PTR_NEW(errorBars)
ENDIF
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
resultStr.values = PTR_NEW(intValue)
resultStr.vrange = [0, 0]
resultStr.Grid1 = GridOut
;
; Clean up
;
PTR_FREE, ValuePTR

result = OBJ_NEW("GKVs1D", resultStr)
result -> GridRestrict, 1

IF(TypeOF(Extra) EQ 8) THEN	$
	result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs2D::INT ****** ;


Function GKVs3D::INT, arg, _Extra=extra
;
; Purpose:
;
;		This function returns the integral of 'self' over the selected
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
;			used to identify independent variable to be integrated 
;			over. Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be integrated over, and to set the range of integration.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to to set the
;			range of integration through indices into the corresponding
;			grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of integration.
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
				MESSAGE, 'INT called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF

Sum = 0b
result = GetKeyWord("Sum", Extra)
IF( Query_Integer(result)) THEN Sum=result

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
axisUnits = grid.units

resultStr.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
resultStr.units = '(' + self.units + ')*(' + axisUnits +')'
IF(Sum) THEN BEGIN
	resultStr.Title = '!7R!X!D' + axisTitle  + '!N' + self.title 
	resultStr.Units = self.units
ENDIF 
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create dx-weighted values
;
IF(Sum) THEN BEGIN
	dx = 1. + 0.*( (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)] )
ENDIF ELSE BEGIN
	dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
ENDELSE

values   = MAKE_ARRAY(SIZE=Info)
CASE axis OF
	1:	FOR i=0, info[1]-2 DO   values[i, *, *] =  0.5*( (*valuePtr)[i, *, *] + (*valuePtr)[i+1, *, *] )*dx[i]
	2:	FOR i=0, info[2]-2 DO   values[*, i, *] =  0.5*( (*valuePtr)[*, i, *] + (*valuePtr)[*, i+1, *] )*dx[i]
	3:	FOR i=0, info[3]-2 DO   values[*, *, i] =  0.5*( (*valuePtr)[*, *, i] + (*valuePtr)[*, *, i+1] )*dx[i]
ENDCASE
;
; Compute Average and standard deviation
;
intValue   = TOTAL(values, axis)
;
; Set up array of errors
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errorPtr = self -> GetErrors()
	;ErrorSq = MAKE_ARRAY(SIZE=Info) 
CASE axis OF
	1:	FOR i=0, info[1]-2 DO   ErrorSq[i, *, *] = ( 0.5*(*errorPtr)[i, *, *] + (*errorPtr)[i+1, *, *] )^2*dx[i]
	2:	FOR i=0, info[2]-2 DO   ErrorSq[*, i, *] = ( 0.5*(*errorPtr)[*, i, *] + (*errorPtr)[*, i+1, *] )^2*dx[i]
	3:	FOR i=0, info[3]-2 DO   ErrorSq[*, *, i] = ( 0.5*(*errorPtr)[*, *, i] + (*errorPtr)[*, *, i+1] )^2*dx[i]
ENDCASE
	avgErrorSq = TOTAL(errorSq, axis)
	errorBars = SQRT(avgErrorSq)
	PTR_FREE, errorPtr
	PTR_FREE, resultStr.ErrorBars
	resultStr.ErrorBars = PTR_NEW(errorBars)
ENDIF
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
resultStr.values = PTR_NEW(intValue)
resultStr.vrange = [0, 0]
resultStr.Grid1 = Grid1
resultStr.Grid2 = Grid2
;
; Clean up
;
PTR_FREE, ValuePTR

result = OBJ_NEW("GKVs2D", resultStr)
result -> GridRestrict, 1
result -> GridRestrict, 2

IF(TypeOF(Extra) EQ 8) THEN	$
	result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs3D::Int ****** ;



Function GKVs4D::Int,arg, _Extra=extra
;
; Purpose:
;
;		This function returns the integral of 'self' over the selected
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
;			used to identify independent variable to be integrated 
;			over. Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			to be integrated over, and to set the range of integration.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to to set the
;			range of integration through indices into the corresponding
;			grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of integration.
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
				MESSAGE, 'INT called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF

Sum = 0b
result = GetKeyWord("Sum", Extra)
IF( Query_Integer(result)) THEN Sum=result

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
axisUnits = grid.units

resultStr.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
resultStr.units = '(' + self.units + ')*(' + axisUnits +')'
IF(Sum) THEN BEGIN
	resultStr.Title = '!7R!X!D' + axisTitle  + '!N' + self.title 
	resultStr.Units = self.units
ENDIF 
newIndices = self -> IndexRemove(axis)
resultStr.indices = PTR_NEW(newIndices)
;
; Create dx-weighted values
IF(Sum) THEN BEGIN
	dx = 1.+ 0.*( (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)] )
ENDIF ELSE BEGIN
	dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
ENDELSE

values   = MAKE_ARRAY(SIZE=Info)
CASE axis OF
	1:	FOR i=0, info[1]-2 DO   values[i, *, *, *] =  0.5*( (*valuePtr)[i, *, *, *] + (*valuePtr)[i+1, *, *, *])*dx[i]
	2:	FOR i=0, info[2]-2 DO   values[*, i, *, *] =  0.5*( (*valuePtr)[*, i, *, *] + (*valuePtr)[*, i+1, *, *])*dx[i]
	3:	FOR i=0, info[3]-2 DO   values[*, *, i, *] =  0.5*( (*valuePtr)[*, *, i, *] + (*valuePtr)[*, *, i+1, *])*dx[i]
	4:	FOR i=0, info[4]-2 DO   values[*, *, *, i] =  0.5*( (*valuePtr)[*, *, *, i] + (*valuePtr)[*, *, *, i+1])*dx[i]
ENDCASE
;
; Compute Average and standard deviation
;
intValue = TOTAL(values, axis)
;
; Set up array of errors
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errorPtr = self -> GetErrors()
	ErrorSq = MAKE_ARRAY(SIZE=Info) 
CASE axis OF
	1:	FOR i=1, info[1]-1 DO   ErrorSq[i, *, *, *] = ( 0.5*(*errorPtr)[i, *, *, *] + (*errorPtr)[i+1, *, *, *] )^2*dx[i]
	2:	FOR i=1, info[2]-1 DO   ErrorSq[*, i, *, *] = ( 0.5*(*errorPtr)[*, i, *, *] + (*errorPtr)[*, i+1, *, *] )^2*dx[i]
	3:	FOR i=1, info[3]-1 DO   ErrorSq[*, *, i, *] = ( 0.5*(*errorPtr)[*, *, i, *] + (*errorPtr)[*, *, i+1, *] )^2*dx[i]
	4:	FOR i=1, info[4]-1 DO   ErrorSq[*, *, *, i] = ( 0.5*(*errorPtr)[*, *, *, i] + (*errorPtr)[*, *, *, i+1] )^2*dx[i]
ENDCASE
	avgErrorSq = TOTAL(errorSq, axis)
	errorBars = SQRT(avgErrorSq)
	PTR_FREE, errorPtr
	PTR_FREE, resultStr.ErrorBars
	resultStr.ErrorBars = PTR_NEW(errorBars)
ENDIF
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
resultStr.values = PTR_NEW(intValue)
resultStr.vrange = [0, 0]
resultStr.Grid1 = Grid1
resultStr.Grid2 = Grid2
resultStr.Grid3 = Grid3
;
; Clean up
;
PTR_FREE, ValuePtr

result = OBJ_NEW("GKVs3D", resultStr)
result -> GridRestrict, 1
result -> GridRestrict, 2
result -> GridRestrict, 3

IF(TypeOF(Extra) EQ 8) THEN	$
	result -> Set, _EXTRA=Extra

RETURN, result

END ; ****** GKVs4D::Int ****** ;
