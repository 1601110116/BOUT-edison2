FORWARD_FUNCTION GKVsd_MIN
FUNCTION GKVs1D::Cat, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10
;
; Purpose:
;
;	This function concatenates its argument(s)
;	onto 'self' and returns the result.  
;	'Self' is left unaltered. 
;
; Output:
;
;	The output is a GKVsd object of the same dimensionality as 'self'
;	in which the data and grid values from 'arg1' & Co. have been concatenated
;	behind those of 'self'.
;
; Input Arguments:
;
;	Cat must be called with at least one argument, which must be
;	a GKVsd object of the same dimensionality as "self".  If called
;	with no argument, or with inappropriate arguments, then CAT returns
;	a (deep) copy of 'self'.
;
; Output Arguments:
;
;	NONE
;
;
;
; Written by W.M. Nevins
;	2/15/01
;
result = self -> MakeCopy()			; Make deep copy of 'self' for return if arg1 is inappropriate
;
; Check that argument is a GKVs1D object
;
IF(TypeOf(arg1) EQ 0)  THEN RETURN, result	; 'self' concatenated with nothing is self
IF(TypeOf(arg1) NE 11) THEN BEGIN		; 'arg1' is not an object
	MESSAGE, "Argument is not an Object", /INFORMATIONAL
	RETURN, result
ENDIF
CASE OBJ_ISA(arg1, "GKVs1D") OF
	0	:	BEGIN
				MESSAGE, "Argument is not a GKVs1D object", /INFORMATIONAL
				RETURN, result
			END
	1	:	; arg1 is a GKVs1D object. Check for like title, mnemonic, units, etc. here?
ENDCASE
nElements1 = N_ELEMENTS(arg1) 
IF(nElements1 GT 1) THEN BEGIN		; arg1 is an array
	previousResult = result		; so must call CAT recursively
	result = previousResult -> Cat(arg1[0], arg1[1:nElements1-1], arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
	previousResult -> Trash
	RETURN, result
ENDIF 
;
; Get 'time' grid of 'self' and first argument
;
selfGrid = self.grid1
argGrid  = arg1.grid1
sIrange = selfGrid.irange
aIrange =  argGrid.irange
selfGridValues = *selfGrid.values
 argGridValues =  *argGrid.values
;
; Check that grid values don't overlap.  If they do, find appropriate
; starting index in 'arg'
;
IF(selfGridValues[sIrange[1]] GE argGridValues[aIrange[0]]) THEN BEGIN
	MESSAGE, "Time grids overlap.  Will disregard beginning of 'arg'", /INFORMATIONAL
	temp = (argGridValues - selfGridVAlues[sIrange[1]])^2
	smallStuff = MIN(temp, aStart)
	IF(aStart GT aIrange[1]) THEN RETURN, result
	IF(selfGridValues[sIrange[1]] GE argGridValues[aStart]) THEN aStart = aStart + 1
	IF(aStart GT aIrange[1]) THEN RETURN, result
	aIrangeSave=aIrange
	aIrange[0] = aStart
	arg1 -> signalWindow, axis=1, Irange=aIrange	
ENDIF
;
; Get pointers to required values from 'self' and 'arg'
;
selfValues = self -> GetValues()
argValues  = arg1 -> GetValues()
nTself = sIrange[1] - sIrange[0] + 1
nTarg  = airange[1] - aIrange[0] + 1
nT = nTself + nTarg
;
; set up "resultInfo" array for creating array to hold output
;
selfInfo = SIZE(*selfValues)
 argInfo = SIZE( *argValues)
resultInfo = selfInfo
resultInfo[selfInfo[0]] = selfInfo[selfInfo[0]] + argInfo[argInfo[0]]
resultInfo[selfInfo[0]+2] = selfInfo[selfInfo[0]+2] + argInfo[argInfo[0]+2]
;
; Create array to hold output
;
resultValues = MAKE_ARRAY(SIZE=resultInfo)
;
; Create array to hold error bars if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	errorInfo = resultInfo
	errorInfo[selfInfo[0]+1] = 4 ; set 'type' of errorbars to floating point
	resultErrors = MAKE_ARRAY(SIZE=errorInfo)
ENDIF
;
; Copy information from 'self' and 'arg' into output array
;
resultValues[ 0:(nTself-1)] = *selfValues
resultValues[ntSelf:(nT-1)] =  *argValues
;
; Copy error bars from 'self' and 'arg' into output array if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	selfErrors = self -> GetErrors()
	argErrors  = arg1 -> GetErrors()
	resultErrors[0:(nTself-1)] = *selfErrors
	resultErrors[ntSelf:(nT-1)] =  *argErrors
	PTR_FREE, selfErrors, argErrors, result.errorBars
	result.ErrorBars = PTR_NEW(resultErrors)
ENDIF
;
; Set up output grid structure
;
resultGridValues = FLTARR(nT)
resultGridValues[ 0:(nTself-1)] = selfGridValues[sIrange[0]:sIrange[1]]
resultGridValues[ntSelf:(nT-1)] =  argGridValues[aIrange[0]:aIrange[1]]
;
; Load new values (and new grid values) into 'result'
;
PTR_FREE, result.values, selfValues, argValues
result.values = PTR_NEW(resultValues)
vmin = GKVsd_MIN(resultValues, MAX=vmax)
result.vrange = [vmin, vmax]
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(resultGridValues)
result.Grid1.irange = [0,nT-1]
tmin = GKVsd_Min(resultGridValues, Max=tmax)
result.Grid1.range = [tmin, tmax]
;
; Restore Irange field of arg
;
IF(N_ELEMENTS(aIrangeSave) EQ 2) THEN	$
	arg1 -> Signalwindow, axis=1, Irange=aIrangeSave
;
; Call 'cat' recursively to concatenate remaining arguments
;
IF(N_PARAMS() GT 1) THEN BEGIN
	previousResult = result
	result = previousResult -> Cat(arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
	previousResult -> trash
ENDIF
;
; and we're done
;
RETURN, result
END ; ****** GKVs1D::Cat ****** ;


FUNCTION GKVs2D::Cat, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10
;
; Purpose:
;
;	This function concatenates its argument(s)
;	onto 'self' and returns the result.  
;	'Self' is left unaltered. 
;
; Output:
;
;	The output is a GKVsd object of the same dimensionality as 'self'
;	in which the data and grid values from 'arg1' & Co. have been concatenated
;	behind those of 'self'.
;
; Input Arguments:
;
;	Cat must be called with at least one argument, which must be
;	a GKVsd object of the same dimensionality as "self".  If called
;	with no argument, or with inappropriate arguments, then CAT returns
;	a (deep) copy of 'self'.
;
; Output Arguments:
;
;	NONE
;
;
;
; Written by W.M. Nevins
;	2/15/01
;
result = self -> MakeCopy()			; Make deep copy of 'self' for return if arg1 is inappropriate
;
; Check that argument is a GKVs2D object
;
IF(TypeOf(arg1) EQ 0)  THEN RETURN, result	; 'self' concatenated with nothing is self
IF(TypeOf(arg1) NE 11) THEN BEGIN		; 'arg1' is not an object
	MESSAGE, "Argument is not an Object", /INFORMATIONAL
	RETURN, result
ENDIF
CASE OBJ_ISA(arg1, "GKVs2D") OF
	0	:	BEGIN
				MESSAGE, "Argument is not a GKVs2D object", /INFORMATIONAL
				RETURN, result
			END
	1	:	; arg1 is a GKVs2D object. Check for like title, mnemonic, units, etc. here?
ENDCASE 
nElements1 = N_ELEMENTS(arg1) 
IF(nElements1 GT 1) THEN BEGIN		; arg1 is an array
	previousResult = result		; so must call CAT recursively
	result = previousResult -> Cat(arg1[0], arg1[1:nElements1-1], arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
	previousResult -> Trash
	RETURN, result
ENDIF 
;
; Check that Grid1 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid1, arg1.Grid1) ) THEN BEGIN
	MESSAGE, "Self.Grid1 does not match arg.Grid1", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Get 'time' grid of 'self' and first argument
;
selfGrid = self.grid2
argGrid  = arg1.grid2
sIrange = selfGrid.irange
aIrange =  argGrid.irange
selfGridValues = *selfGrid.values
 argGridValues =  *argGrid.values
;
; Check that grid values don't overlap.  If they do, find appropriate
; starting index in 'arg'
;
IF(selfGridValues[sIrange[1]] GE argGridValues[aIrange[0]]) THEN BEGIN
	MESSAGE, "Time grids overlap.  Will disregard beginning of 'arg'", /INFORMATIONAL
	temp = (argGridValues - selfGridVAlues[sIrange[1]])^2
	smallStuff = MIN(temp, aStart)
	IF(aStart GT aIrange[1]) THEN RETURN, result
	IF(selfGridValues[sIrange[1]] GE argGridValues[aStart]) THEN aStart = aStart + 1
	IF(aStart GT aIrange[1]) THEN RETURN, result
	aIrangeSave=aIrange
	aIrange[0] = aStart
	arg1 -> signalWindow, axis=2, Irange=aIrange	
ENDIF
;
; Get pointers to required values from 'self' and 'arg'
;
selfValues = self -> GetValues()
argValues  = arg1 -> GetValues()
nTself = sIrange[1] - sIrange[0] + 1
nTarg  = airange[1] - aIrange[0] + 1
nT = nTself + nTarg
;
; set up "resultInfo" array for creating array to hold output
;
selfInfo = SIZE(*selfValues)
 argInfo = SIZE( *argValues)
resultInfo = selfInfo
resultInfo[selfInfo[0]] = selfInfo[selfInfo[0]] + argInfo[argInfo[0]]
resultInfo[selfInfo[0]+2] = selfInfo[selfInfo[0]+2] + argInfo[argInfo[0]+2]
;
; Create array to hold output
;
resultValues = MAKE_ARRAY(SIZE=resultInfo)
;
; Create array to hold error bars if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	errorInfo = resultInfo
	errorInfo[selfInfo[0]+1] = 4 ; set 'type' of errorbars to floating point
	resultErrors = MAKE_ARRAY(SIZE=errorInfo)
ENDIF
;
; Copy data from 'self' and 'arg' into output array
;
resultValues[*, 0:(nTself-1)] = *selfValues
resultValues[*,ntSelf:(nT-1)] =  *argValues
PTR_FREE, selfValues, argValues
;
; Copy error bars from 'self' and 'arg' into output array if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	selfErrors = self -> GetErrors()
	argErrors  = arg1 -> GetErrors()
	resultErrors[*, 0:(nTself-1)] = *selfErrors
	resultErrors[*,ntSelf:(nT-1)] =  *argErrors
	PTR_FREE, result.errorBars, selfErrors, argErrors
	result.ErrorBars = PTR_NEW(resultErrors)
ENDIF
;
; Set up output grid structure
;
resultGridValues = FLTARR(nT)
resultGridValues[ 0:(nTself-1)] = selfGridValues[sIrange[0]:sIrange[1]]
resultGridValues[ntSelf:(nT-1)] =  argGridValues[aIrange[0]:aIrange[1]]
;
; Load new values (and new grid values) into 'result'
;
PTR_FREE, result.values
result.values = PTR_NEW(resultValues)
vmin = GKVsd_MIN(resultValues, MAX=vmax)
result.vrange = [vmin, vmax]
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(resultGridValues)
result.Grid2.irange = [0,nT-1]
tmin = GKVsd_Min(resultGridValues, Max=tmax)
result.Grid2.range = [tmin, tmax]
;
; Restore Irange field of arg
;
IF(N_ELEMENTS(aIrangeSave) EQ 2) THEN	$
	arg1 -> Signalwindow, axis=2, Irange=aIrangeSave
;
; Call 'cat' recursively to concatenate remaining arguments
;
IF(N_PARAMS() GT 1) THEN BEGIN
	previousResult = result
	result = previousResult -> Cat(arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
	previousResult -> trash
ENDIF
;
; and we're done
;
RETURN, result
END ; ****** GKVs2D::Cat ****** ;


FUNCTION GKVs3D::Cat, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10
;
; Purpose:
;
;	This function concatenates its argument(s)
;	onto 'self' and returns the result.  
;	'Self' is left unaltered. 
;
; Output:
;
;	The output is a GKVsd object of the same dimensionality as 'self'
;	in which the data and grid values from 'arg1' & Co. have been concatenated
;	behind those of 'self'.
;
; Input Arguments:
;
;	Cat must be called with at least one argument, which must be
;	a GKVsd object of the same dimensionality as "self".  If called
;	with no argument, or with inappropriate arguments, then CAT returns
;	a (deep) copy of 'self'.
;
; Output Arguments:
;
;	NONE
;
;
;
; Written by W.M. Nevins
;	2/15/01
;
result = self -> MakeCopy()			; Make deep copy of 'self' for return if arg1 is inappropriate
;
; Check that argument is a GKVs3D object
;
IF(TypeOf(arg1) EQ 0)  THEN RETURN, result	; 'self' concatenated with nothing is self
IF(TypeOf(arg1) NE 11) THEN BEGIN		; 'arg1' is not an object
	MESSAGE, "Argument is not an Object", /INFORMATIONAL
	RETURN, result
ENDIF
CASE OBJ_ISA(arg1, "GKVs3D") OF
	0	:	BEGIN
				MESSAGE, "Argument is not a GKVs3D object", /INFORMATIONAL
				RETURN, result
			END
	1	:	; arg1 is a GKVs3D object. Check for like title, mnemonic, units, etc. here?
ENDCASE 
nElements1 = N_ELEMENTS(arg1) 
IF(nElements1 GT 1) THEN BEGIN		; arg1 is an array
	previousResult = result		; so must call CAT recursively
	result = previousResult -> Cat(arg1[0], arg1[1:nElements1-1], arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
	previousResult -> Trash
	RETURN, result
ENDIF 
;
; Check that Grid1 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid1, arg1.Grid1) ) THEN BEGIN
	MESSAGE, "Self.Grid1 does not match arg.Grid1", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Check that Grid2 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid2, arg1.Grid2) ) THEN BEGIN
	MESSAGE, "Self.Grid2 does not match arg.Grid2", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Get 'time' grid of 'self' and first argument
;
selfGrid = self.grid3
argGrid  = arg1.grid3
sIrange = selfGrid.irange
aIrange =  argGrid.irange
selfGridValues = *selfGrid.values
 argGridValues =  *argGrid.values
;
; Check that grid values don't overlap.  If they do, find appropriate
; starting index in 'arg'
;
IF(selfGridValues[sIrange[1]] GE argGridValues[aIrange[0]]) THEN BEGIN
	MESSAGE, "Time grids overlap.  Will disregard beginning of 'arg'", /INFORMATIONAL
	temp = (argGridValues - selfGridVAlues[sIrange[1]])^2
	smallStuff = MIN(temp, aStart)
	IF(aStart GT aIrange[1]) THEN RETURN, result
	IF(selfGridValues[sIrange[1]] GE argGridValues[aStart]) THEN aStart = aStart + 1
	IF(aStart GT aIrange[1]) THEN RETURN, result
	aIrangeSave=aIrange
	aIrange[0] = aStart
	arg1 -> signalWindow, axis=3, Irange=aIrange	
ENDIF
;
; Get pointers to required values from 'self' and 'arg'
;
selfValues = self -> GetValues()
argValues  = arg1 -> GetValues()
nTself = sIrange[1] - sIrange[0] + 1
nTarg  = airange[1] - aIrange[0] + 1
nT = nTself + nTarg
;
; set up "resultInfo" array for creating array to hold output
;
selfInfo = SIZE(*selfValues)
 argInfo = SIZE( *argValues)
resultInfo = selfInfo
resultInfo[selfInfo[0]] = selfInfo[selfInfo[0]] + argInfo[argInfo[0]]
resultInfo[selfInfo[0]+2] = selfInfo[selfInfo[0]+2] + argInfo[argInfo[0]+2]
;
; Create array to hold output
;
resultValues = MAKE_ARRAY(SIZE=resultInfo)
;
; Create array to hold error bars if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	errorInfo = resultInfo
	errorInfo[selfInfo[0]+1] = 4 ; set 'type' of errorbars to floating point
	resultErrors = MAKE_ARRAY(SIZE=errorInfo)
ENDIF
;
; Copy data from 'self' and 'arg' into output array
;
resultValues[*, *, 0:(nTself-1)] = *selfValues
resultValues[*, *,ntSelf:(nT-1)] =  *argValues
;
; Copy error bars from 'self' and 'arg' into output array if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	selfErrors = self -> GetErrors()
	argErrors  = arg1 -> GetErrors()
	resultErrors[*, *, 0:(nTself-1)] = *selfErrors
	resultErrors[*, *,ntSelf:(nT-1)] =  *argErrors
	PTR_FREE, result.errorBars, selfErrors, argErrors
	result.ErrorBars = PTR_NEW(resultErrors)
ENDIF
;
; Set up output grid structure
;
resultGridValues = FLTARR(nT)
resultGridValues[ 0:(nTself-1)] = selfGridValues[sIrange[0]:sIrange[1]]
resultGridValues[ntSelf:(nT-1)] =  argGridValues[aIrange[0]:aIrange[1]]
;
; Load new values (and new grid values) into 'result'
;
PTR_FREE, result.values, selfValues, argValues
result.values = PTR_NEW(resultValues)
vmin = GKVsd_MIN(resultValues, MAX=vmax)
result.vrange = [vmin, vmax]
PTR_FREE, result.Grid3.values
result.Grid3.values = PTR_NEW(resultGridValues)
result.Grid3.irange = [0,nT-1]
tmin = GKVsd_Min(resultGridValues, Max=tmax)
result.Grid3.range = [tmin, tmax]
;
; Restore Irange field of arg
;
IF(N_ELEMENTS(aIrangeSave) EQ 2) THEN	$
	arg1 -> Signalwindow, axis=3, Irange=aIrangeSave
;
; Call 'cat' recursively to concatenate remaining arguments
;
IF(N_PARAMS() GT 1) THEN BEGIN
	previousResult = result
	result = previousResult -> Cat(arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
	previousResult -> trash
ENDIF
;
; and we're done
;
RETURN, result
END ; ****** GKVs3D::Cat ****** ;


FUNCTION GKVs4D::Cat, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10
;
; Purpose:
;
;	This function concatenates its argument(s)
;	onto 'self' and returns the result.  
;	'Self' is left unaltered. 
;
; Output:
;
;	The output is a GKVsd object of the same dimensionality as 'self'
;	in which the data and grid values from 'arg1' & Co. have been concatenated
;	behind those of 'self'.
;
; Input Arguments:
;
;	Cat must be called with at least one argument, which must be
;	a GKVsd object of the same dimensionality as "self".  If called
;	with no argument, or with inappropriate arguments, then CAT returns
;	a (deep) copy of 'self'.
;
; Output Arguments:
;
;	NONE
;
;
;
; Written by W.M. Nevins
;	2/15/01
; Modified by W.M. Nevins
; to concatenate 4D objects
;	2/11/09
;
result = self -> MakeCopy()			; Make deep copy of 'self' for return if arg1 is inappropriate
;
; Check that argument is a GKVs4D object
;
IF(TypeOf(arg1) EQ 0)  THEN RETURN, result	; 'self' concatenated with nothing is self
IF(TypeOf(arg1) NE 11) THEN BEGIN		; 'arg1' is not an object
	MESSAGE, "Argument is not an Object", /INFORMATIONAL
	RETURN, result
ENDIF
CASE OBJ_ISA(arg1, "GKVs4D") OF
	0	:	BEGIN
				MESSAGE, "Argument is not a GKVs4D object", /INFORMATIONAL
				RETURN, result
			END
	1	:	; arg1 is a GKVs3D object. Check for like title, mnemonic, units, etc. here?
ENDCASE 
nElements1 = N_ELEMENTS(arg1) 
IF(nElements1 GT 1) THEN BEGIN		; arg1 is an array
	previousResult = result		; so must call CAT recursively
	result = previousResult -> Cat(arg1[0], arg1[1:nElements1-1], arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
	previousResult -> Trash
	RETURN, result
ENDIF 
;
; Check that Grid1 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid1, arg1.Grid1) ) THEN BEGIN
	MESSAGE, "Self.Grid1 does not match arg.Grid1", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Check that Grid2 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid2, arg1.Grid2) ) THEN BEGIN
	MESSAGE, "Self.Grid2 does not match arg.Grid2", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Check that Grid3 of 'self' and 'arg1' match
;
IF( NOT GKV_GridSame(self.Grid3, arg1.Grid3) ) THEN BEGIN
	MESSAGE, "Self.Grid3 does not match arg.Grid3", /INFORMATIONAL
	RETURN, result
ENDIF
;
; Get 'time' grid of 'self' and first argument
;
selfGrid = self.grid4
argGrid  = arg1.grid4
sIrange = selfGrid.irange
aIrange =  argGrid.irange
selfGridValues = *selfGrid.values
 argGridValues =  *argGrid.values
;
; Check that grid values don't overlap.  If they do, find appropriate
; starting index in 'arg'
;
IF(selfGridValues[sIrange[1]] GE argGridValues[aIrange[0]]) THEN BEGIN
	MESSAGE, "Time grids overlap.  Will disregard beginning of 'arg'", /INFORMATIONAL
	temp = (argGridValues - selfGridVAlues[sIrange[1]])^2
	smallStuff = MIN(temp, aStart)
	IF(aStart GT aIrange[1]) THEN RETURN, result
	IF(selfGridValues[sIrange[1]] GE argGridValues[aStart]) THEN aStart = aStart + 1
	IF(aStart GT aIrange[1]) THEN RETURN, result
	aIrangeSave=aIrange
	aIrange[0] = aStart
	arg1 -> signalWindow, axis=4, Irange=aIrange	
ENDIF
;
; Get pointers to required values from 'self' and 'arg'
;
selfValues = self -> GetValues()
argValues  = arg1 -> GetValues()
nTself = sIrange[1] - sIrange[0] + 1
nTarg  = airange[1] - aIrange[0] + 1
nT = nTself + nTarg
;
; set up "resultInfo" array for creating array to hold output
;
selfInfo = SIZE(*selfValues)
 argInfo = SIZE( *argValues)
resultInfo = selfInfo
resultInfo[selfInfo[0]] = selfInfo[selfInfo[0]] + argInfo[argInfo[0]]
resultInfo[selfInfo[0]+2] = selfInfo[selfInfo[0]+2] + argInfo[argInfo[0]+2]
;
; Create array to hold output
;
resultValues = MAKE_ARRAY(SIZE=resultInfo)
;
; Create array to hold error bars if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	errorInfo = resultInfo
	errorInfo[selfInfo[0]+1] = 4 ; set 'type' of errorbars to floating point
	resultErrors = MAKE_ARRAY(SIZE=errorInfo)
ENDIF
;
; Copy data from 'self' and 'arg' into output array
;
resultValues[*, *, *,  0:(nTself-1)] = *selfValues
resultValues[*, *, *, ntSelf:(nT-1)] =  *argValues
;
; Copy error bars from 'self' and 'arg' into output array if required
;
IF( PTR_VALID(self.errorbars) ) THEN BEGIN
	selfErrors = self -> GetErrors()
	argErrors  = arg1 -> GetErrors()
	resultErrors[*, *, *, 0:(nTself-1)] = *selfErrors
	resultErrors[*, *, *,ntSelf:(nT-1)] =  *argErrors
	PTR_FREE, result.errorBars, selfErrors, argErrors
	result.ErrorBars = PTR_NEW(resultErrors)
ENDIF
;
; Set up output grid structure
;
resultGridValues = FLTARR(nT)
resultGridValues[ 0:(nTself-1)] = selfGridValues[sIrange[0]:sIrange[1]]
resultGridValues[ntSelf:(nT-1)] =  argGridValues[aIrange[0]:aIrange[1]]
;
; Load new values (and new grid values) into 'result'
;
PTR_FREE, result.values, selfValues, argValues
result.values = PTR_NEW(resultValues)
vmin = GKVsd_MIN(resultValues, MAX=vmax)
result.vrange = [vmin, vmax]
PTR_FREE, result.Grid4.values
result.Grid4.values = PTR_NEW(resultGridValues)
result.Grid4.irange = [0,nT-1]
tmin = GKVsd_Min(resultGridValues, Max=tmax)
result.Grid4.range = [tmin, tmax]
;
; Restore Irange field of arg
;
IF(N_ELEMENTS(aIrangeSave) EQ 2) THEN	$
	arg1 -> Signalwindow, axis=4, Irange=aIrangeSave
;
; Call 'cat' recursively to concatenate remaining arguments
;
IF(N_PARAMS() GT 1) THEN BEGIN
	previousResult = result
	result = previousResult -> Cat(arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
	previousResult -> trash
ENDIF
;
; and we're done
;
RETURN, result
END ; ****** GKVs4D::Cat ****** ;
