Function GKVs1D::CINT, arg, _Extra=extra
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
;	    start	Set 'start' to the (floating point) value at which the cumulative
;			integral should be zero.  Defaults to first element of 'range'.
;			(Optional).
;
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
; Further revised by W.M. Nevins
;	11/5/03
; Added 'start' keyword to GKVs1D::CInt and GKVs2D::CInt
;
;
FORWARD_FUNCTION GKVsd_MIN
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
result = self -> MakeCopy()
result -> restrict
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = result -> GetValues()
info = SIZE(*valuePtr)
nPoints = info[axis]
values = MAKE_ARRAY(SIZE=info)

irange = grid.irange
axisTitle = grid.title
axisUnits = grid.units

result.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
result.units = '(' + self.units + ')*(' + axisUnits +')'
;
; Compute integral 
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values[1:(nPoints-1)] = 0.5*((*self.values)[(irange[0]+1):irange[1]] + (*self.values)[irange[0]:(irange[1]-1)])*dx
values[0] = 0.
intValue  = TOTAL(values, /CUMULATIVE)
;
; Chect for "start" keyword
;
startValue = GetKeyWord('start', Extra)
IF(TypeOf(startValue) NE 7) THEN BEGIN
	test = (*Grid.values - startValue)^2
	eps = MIN(test, iStart)
	intValue = intValue - intValue[iStart]
ENDIF
;
; Set up array of errors
;
IF PTR_VALID(result.ErrorBars) THEN BEGIN
	errors = 0.5*((*result.ErrorBars)[(irange[0]+1):irange[1]] + (*result.ErrorBars)[irange[0]:(irange[1]-1)])
	errorSQ = FLTARR(nPoints)
	errorSq[1:(nPoints-1)] = FLOAT(errors*errors*dx)
	ErrorSq = TOTAL(errorSq, /CUMULATIVE)
	errorBars = SQRT(ErrorSq)
	PTR_FREE, result.errorbars
	result.ErrorBars = PTR_NEW(errorBars)
ENDIF
;
; Load results into 'result' structure
;
PTR_FREE, valuePtr
PTR_FREE, result.values
result.values = PTR_NEW(intValue)
vmin = GKVsd_Min(intValue, MAX=vmax)
result.vrange = [vmin, vmax]

IF(TypeOF(Extra) EQ 8) THEN result -> Set, _EXTRA=Extra


RETURN, result
END ; ****** GKVs1D::INT ****** ;


Function GKVs2D::CINT, arg, _Extra=extra
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
;	    start	Set 'start' to the (floating point) value at which the cumulative
;			integral should be zero.  Defaults to first element of 'range'.
;			(Optional).
;
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
; Further revised by W.M. Nevins
;	11/5/03
; Added 'start' keyword to GKVs1D::CInt and GKVs2D::CInt
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
result = self -> MakeCopy()
result -> restrict

axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]

irange = grid.irange
axisTitle = grid.title
axisUnits = grid.units

result.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
result.units = '(' + self.units + ')*(' + axisUnits +')'
;
; Create 2-D 'dx' array
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values = MAKE_ARRAY(SIZE=info)
CASE axis OF
	1:	BEGIN
			ones = MAKE_ARRAY(info[2], value=1)
			dx = dx#ones
			values[1:(nPoints-1), *] = 0.5*( (*valuePtr)[(irange[0]+1):irange[1], *] + (*valuePtr)[irange[0]:(irange[1]-1), *])*dx
		END
	2:	BEGIN
			ones = MAKE_ARRAY(info[1], value=1)
			dx = ones#dx
			values[ *, 1:(nPoints-1)] = 0.5*( (*valuePtr)[*, (irange[0]+1):irange[1]] + (*valuePtr)[*, irange[0]:(irange[1]-1)])*dx
		END
ENDCASE
;
; Compute integral
;
intValue   = TOTAL(values, axis, /CUMULATIVE)
;
; Chect for "start" keyword
;
startValue = GetKeyWord('start', Extra)
IF(TypeOf(StartValue) NE 7) THEN BEGIN
	test = (*Grid.values - startValue)^2
	eps = MIN(test, iStart)
	CASE axis OF
		1:	intValue0 = MAKE_ARRAY(info[1], VALUE=1.)#REFORM(intValue[iStart, *])
		2:	intValue0 = REFORM(intValue[*, iStart])#MAKE_ARRAY(info[2], VALUE=1.)
	ENDCASE
	intValue = intValue - intValue0
ENDIF
;
; Set up array of errors
;
IF PTR_VALID(result.ErrorBars) THEN BEGIN
	errorPtr = result -> GetErrors()
	errorSq = FLTARR(info[1], info[2])
	CASE axis OF
		1:	errorSq[1:(nPoints-1), *] = (0.5*((*errorPtr)[(irange[0]+1):irange[1], *] + (*errorPtr)[irange[0]:(irange[1]-1), *]))^2*dx
		2:	errorSq[*, 1:(nPoints-1)] = (0.5*((*errorPtr)[*, (irange[0]+1):irange[1]] + (*errorPtr)[*, irange[0]:(irange[1]-1)]))^2*dx
	ENDCASE
	PTR_FREE, errorPtr
	ErrorSq = TOTAL(errorSq, axis, /CUMULATIVE)
	errorBars = SQRT(ErrorSq)
	PTR_FREE, result.ErrorBars
	result.ErrorBars = PTR_NEW(errorBars)
ENDIF
;
; Load results into 'result' structure
;
PTR_FREE, result.values
result.values = PTR_NEW(intValue)
vmin = GKVsd_MIN(intValue, MAX=vmax)
result.vrange = [vmin, vmax]
;
; Clean up
;
PTR_FREE, ValuePtr

IF(TypeOF(Extra) EQ 8) THEN result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs2D::CINT ****** ;


Function GKVs3D::CINT, arg, _Extra=extra
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
;	    start	Set 'start' to the (floating point) value at which the cumulative
;			integral should be zero.  Defaults to first element of 'range'.
;			(Optional).  	*** NOT YET IMPLIMENTED IN 3-D ***
;
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
result = self -> MakeCopy()
result -> restrict

axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title
axisUnits = grid.units

result.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
result.units = '(' + self.units + ')*(' + axisUnits +')'
;
; Create dx-weighted values
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values   = MAKE_ARRAY(SIZE=Info)
CASE axis OF
	1:	FOR i=1, info[1]-1 DO   values[i, *, *] =  0.5*( (*valuePtr)[i, *, *] + (*valuePtr)[i+1, *, *])*dx[i]
	2:	FOR i=1, info[2]-1 DO   values[*, i, *] =  0.5*( (*valuePtr)[*, i, *] + (*valuePtr)[*, i+1, *])*dx[i]
	3:	FOR i=1, info[3]-1 DO   values[*, *, i] =  0.5*( (*valuePtr)[*, *, i] + (*valuePtr)[*, *, i+1])*dx[i]
ENDCASE
;
; Compute Average and standard deviation
;
intValue   = TOTAL(values, axis, /CUMULATIVE)
;
; Set up array of errors
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	errorPtr = result -> GetErrors()
	ErrorSq = FLTARR(info[1], info[2], info[3]) 
CASE axis OF
	1:	FOR i=1, info[1]-1 DO   ErrorSq[i, *, *] = (0.5*( *errorPtr)[i, *, *] + (*errorPtr)[i+1, *, *])^2*dx[i]
	2:	FOR i=1, info[2]-1 DO   ErrorSq[*, i, *] = (0.5*( *errorPtr)[*, i, *] + (*errorPtr)[*, i+1, *])^2*dx[i]
	3:	FOR i=1, info[3]-1 DO   ErrorSq[*, *, i] = (0.5*( *errorPtr)[*, *, i] + (*errorPtr)[*, *, i+1])^2*dx[i]
ENDCASE
	gErrorSq = TOTAL(errorSq, axis, /CUMULATIVE)
	errorBars = SQRT(gErrorSq)
	PTR_FREE, errorPtr
	PTR_FREE, result.ErrorBars
	result.ErrorBars = PTR_NEW(errorBars)
ENDIF
;
; Load results into 'result' 
;
PTR_FREE, result.values
result.values = PTR_NEW(intValue)
vmin = GKVsd_MIN(intValue, MAX=vmax)
result.vrange = [vmin, vmax]
;
; Clean up
;
PTR_FREE, ValuePTR

IF(TypeOF(Extra) EQ 8) THEN result -> Set, _EXTRA=Extra


RETURN, result
END ; ****** GKVs3D::CInt ****** ;



Function GKVs4D::CInt,arg, _Extra=extra
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
result = self -> MakeCopy()
result -> restrict

axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
valuePtr = self -> GetValues()

info = SIZE(*valuePtr)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title
axisUnits = grid.units

result.Title = '!9i!X' + self.title + ' d' + axisTitle  + '!N'
result.units = '(' + self.units + ')*(' + axisUnits +')'
;
; Create dx-weighted values
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values   = MAKE_ARRAY(SIZE=Info)
CASE axis OF
	1:	FOR i=1, info[1]-1 DO   values[i, *, *, *] =  0.5*( (*valuePtr)[i, *, *, *] + (*valuePtr)[i+1, *, *, *] )*dx[i]
	2:	FOR i=1, info[2]-1 DO   values[*, i, *, *] =  0.5*( (*valuePtr)[*, i, *, *] + (*valuePtr)[*, i+1, *, *] )*dx[i]
	3:	FOR i=1, info[3]-1 DO   values[*, *, i, *] =  0.5*( (*valuePtr)[*, *, i, *] + (*valuePtr)[*, *, i+1, *] )*dx[i]
	4:	FOR i=1, info[4]-1 DO   values[*, *, *, i] =  0.5*( (*valuePtr)[*, *, *, i] + (*valuePtr)[*, *, *, i+1] )*dx[i]
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
	ErrorSq = FLTARR(info[1], info[2], info[3], info[4])
CASE axis OF
	1:	FOR i=1, info[1]-1 DO   ErrorSq[i, *, *, *] = (  0.5*(*errorPtr)[i, *, *, *] + (*errorPtr)[i+1, *, *, *] )^2*dx[i]
	2:	FOR i=1, info[2]-1 DO   ErrorSq[*, i, *, *] = (  0.5*(*errorPtr)[*, i, *, *] + (*errorPtr)[*, i+1, *, *] )^2*dx[i]
	3:	FOR i=1, info[3]-1 DO   ErrorSq[*, *, i, *] = (  0.5*(*errorPtr)[*, *, i, *] + (*errorPtr)[*, *, i+1, *] )^2*dx[i]
	4:	FOR i=1, info[4]-1 DO   ErrorSq[*, *, *, i] = (  0.5*(*errorPtr)[*, *, *, i] + (*errorPtr)[*, *, *, i+1] )^2*dx[i]
ENDCASE
	ErrorSq = TOTAL(errorSq, axis, /CUMULATIVE)
	errorBars = SQRT(ErrorSq)
	PTR_FREE, errorPtr
	PTR_FREE, result.ErrorBars
	result.ErrorBars = PTR_NEW(errorBars)
ENDIF
;
; Load results into 'result' structure
;
PTR_FREE, result.values
result.values = PTR_NEW(intValue)
vmin = GKVsd_MIN(intValues, MAX=vmax)
result.vrange = [vmin, vmax]
;
; Clean up
;
PTR_FREE, ValuePTR

IF(TypeOF(Extra) EQ 8) THEN result -> Set, _EXTRA=Extra


RETURN, result
END ; ****** GKVs4D::CInt ****** ;
