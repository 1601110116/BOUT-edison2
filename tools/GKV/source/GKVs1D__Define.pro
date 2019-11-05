;
; *****************************************************************************************************************
; ******************************************     Copyright Notice     *********************************************
; *                                                                                                               *
; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *
; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *
; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *
; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *
; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *
; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * 
; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *
; *                                                                                                               *
; *****************************************************************************************************************
;
; *****************************************************************************************************************
; **********************************************     DISCLAIMER     ***********************************************
; *                                                                                                               *
; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *
; *  Neither the United States Government nor the University of California nor any of their employees, makes      *
; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *
; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *
; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *
; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *
; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *
; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *
; *  States Government or the University of California, and shall not be used for advertising or product          *
; *  endorsement purposes.                                                                                        *
; *                                                                                                               *
; *****************************************************************************************************************
;

;
; Defining routine for 1-dimensional signal objects
;
; Written by W.M. Nevins
;	1/15/00
;
FUNCTION GKVs1D::AxisNumber, stringIn, Debug=d
;
; Returns number of axis whose mnemonic matches 'stringIn'
; If there is no match, returns 0
; 
; Written by W.M. Nevins
;	2/20/00
;
debug=0
result=0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
IF(N_ELEMENTS(stringIn) EQ 0) THEN RETURN, result
argType = typeOF(stringIN)								; Check for proper argument type
IF(argType NE 7) THEN BEGIN
	MESSAGE, 'Argument must be a string', Informational=d
	RETURN, result
ENDIF
stringIn = STRTRIM(stringIn, 2)							; Remove both leading and trailing blanks
axis1Mnemonic = STRTRIM(self.Grid1.mnemonic, 2)
IF( STRCMP(stringIn, axis1Mnemonic, /FOLD_CASE) ) THEN result = 1
RETURN, result
END ; ****** GKVs1D::AxisNumber ****** ;


FUNCTION GKVs1D::GetAxis, structure, Debug=d
;
; Searches 'structure' for a tag which is the same as the mnemonic of axis1
; Returns a structure with tag 'axis1'.  
; The associated value is whatever value was associated with the mnemonic of axis1.
;
; On return, 'structure' has tag (and value) corresponding to the mnemonic of axis1
; removed.
;
; If the mnemonic of axis1 is the only tag in structure, then structure=-1 on return
; (null structures are not allowed in IDL 5.3)
;
; NOTE:  the argument structure may be undefined on entry.  In this case it is important
; not to do any thing with 'structure (beyond calling N_ELEMENTS, TypeOf, or SIZE), as
; this would result in an error.  Just return with a null result.
; 
; Written by W.M. Nevins
;	2/20/00
;
debug=1
result={axis1:'no match'}
otherTags=0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
IF(N_ELEMENTS(structure) EQ 0) THEN RETURN, result
arg1Type = typeOF(structure)							; Check for proper argument type
IF(arg1Type NE 8) THEN BEGIN
	structure = -1
	RETURN, result
ENDIF
nTags = N_TAGS(structure)
IF(nTags EQ 0) THEN RETURN, result
tagNames = TAG_NAMES(structure)
tagNames = STRTRIM(tagNames, 2)							; Remove both leading and trailing blanks
command_str = 'structure = {'
axisMnemonic = STRTRIM(self.Grid1.mnemonic)				; Remove both leading and trailing blanks
FOR i=0, ntags-1 DO BEGIN								; Search tags of 'structure'
	IF( STRCMP(axisMnemonic, tagNames[i], /FOLD_CASE) ) THEN BEGIN
		result = {axis1:structure.(i)}					; IF mulitple occurances of axis Mnemonic in
	ENDIF ELSE BEGIN								;	'structure', only the LAST occurance is significant
		i_str = STRING(i, FORMAT='(I3)')				; Construct command string for 'Nstructure'
		i_str = STRTRIM(i_str, 2)						; Remove both leading and trailing blanks
		IF(otherTags NE 0) THEN command_str = command_str + ', '
		command_str = command_str + tagNames[i] + ':' + 'structure.(' + i_str + ')'
		otherTags = otherTags + 1
	ENDELSE
ENDFOR
command_str = command_str + '}'
IF(otherTags NE 0) THEN BEGIN
	ok = EXECUTE(command_str)
ENDIF ELSE BEGIN
	structure = -1
ENDELSE
RETURN, result
END ; ****** GKVs1D::GetAxis ****** ;

FUNCTION GKVs1D::AxisIndex, axisNumber, value, Debug=d
;
; We should only get here if axisNumber=1.  
; Function returns index go Grid1 closest to value
;
debug=0
IF(N_ELEMENTS(d) NE 0) then debug=d
IF(axisNumber NE 1) THEN BEGIN
	MESSAGE, 'axisNumber not equal to 1', INFORMATIONAL=d
	RETURN, 0
ENDIF
temp=(Grid1.values - values)^2
eps=MIN(temp,index)
RETURN, index
END ; ****** GKVs1D::AxisIndex ****** ;


FUNCTION GKVs1D::Interpolate, arg
;
; Interpolate values from self onto the signal window of arg's grid.
;
IF (OBJ_ISA(arg, 'GKVs1D') NE 1) THEN BEGIN
	MESSAGE, "GKVs1D::Interpolate:  Argument is not a valid GKVs1D object", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Check for common units, independent variable, interval
;
IF (self.Grid1.units NE arg.Grid1.units) THEN BEGIN
	MESSAGE, "GKVs1D::Interpolate:  Incompatible units (independent variable)", /INFORMATIONAL
	RETURN, 0
ENDIF
IF (self.Grid1.title NE arg.Grid1.title) THEN BEGIN
	MESSAGE, "GKVs1D::Interpolate:  Incompatible independent variables", /INFORMATIONAL
	RETURN, 0
ENDIF
IF(self -> SameGRID(arg)) THEN BEGIN
	MESSAGE, "Objects have common grid, no interpolation necessary", /INFORMATIONAL
	result = self -> MakeCopy()
	result -> Restrict
	RETURN, result
ENDIF

x1 = *self.Grid1.values					; Get grid values of Obj to be interpolated
x2 =  *arg.Grid1.values					; Get values of target grid
values1 = *self.values					; Get values of to be interpolated
imin = arg.Grid1.irange[0]					; Get range of indices on target grid
imax = arg.Grid1.irange[1]
x = x2[imin:imax]						; For target grid
values = INTERPOL(values1, x1, x)			; Interpolate input values

result = self -> MakeCopy(/NoValues, /NoErrorBars)
;PTR_FREE, result.values					; Load pointer to interpolated values
result.values = PTR_NEW(values)				;	into 'result'
vmin = GKVsd_MIN(values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1
result.vrange = [vmin, vmax]
PTR_FREE, result.Grid1.values			
result.Grid1.values = PTR_NEW(x)			
result.Grid1.range=arg.Grid1.range			; Set signal window of 'result'
result.Grid1.irange = [0,imax-imin]

IF PTR_VALID(self.ErrorBars) THEN BEGIN		; Interpolate error bars...
	error1 = *self.ErrorBars				;	(probably we could think
	errors = INTERPOL(error1, x1, x)		;	 of something a bit more
	;PTR_FREE, result.ErrorBars			;	 rigorous?)
	result.ErrorBars = PTR_NEW(errors)
ENDIF

RETURN, result

END ; ***** GKVs1D::Interpolate ***** ;


FUNCTION GKVs1D::MakeCopy, _Extra=extra
;
; Make "deep" copy of self
;
result = self -> GKVsd::MakeCopy( _Extra=extra)	; Creates new GKVsxD objects, and
									; 	copies elements of GKVsd class
										
result.Grid1 = GKVsd_GridCopy(self.Grid1)		; Make 'deep' copy of the Grid-class structure Grid1,			
RETURN, result

END ; ***** GKVs1D::MakeCopy ***** ;
		

FUNCTION GKVs1D::Plus, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  addition
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg						; Make proxy for argument
try_again:
argInfo = SIZE(arg)					; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN								; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		IF(iaxis NE 1) THEN BEGIN
			MESSAGE, "Bad Axis ID -- Returning"
			RETURN, 0
		ENDIF
		result = self -> MakeCopy(/noValues)
		argValues = *self.Grid1.values	
		selfValues = *self.values
		resultValues = selfValues + argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "+" + self.Grid1.title + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again						;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs1D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs1D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "GKVs1D::Plus:  Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			selfValues = (*self.values)[self.Grid1.irange[0]:self.Grid1.irange[1]]
			values = argValues+selfValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			argErrors = 0.
			gotErrors = 0
			IF(PTR_VALID(result.ErrorBars)) THEN BEGIN
				argErrors = *result.ErrorBars
				gotErrors = 1
			ENDIF
			selfErrors = 0.
			IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
				selfErrors = (*self.ErrorBars)[self.Grid1.irange[0]:self.Grid1.irange[1]]
				gotErrors = 1
			ENDIF
			errors = SQRT( (selfErrors)^2 + (argErrors)^2 )
			IF( gotErrors EQ 1 ) THEN BEGIN
				PTR_FREE, result.ErrorBars
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + "+" + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
			grid1=result.grid1
			gridvalues = *grid1.values
			imin=grid1.irange[0]
			imax=grid1.irange[1]
			newValues = gridvalues[imin:imax]
			PTR_FREE, grid1.values
			grid1.values = PTR_NEW(newValues)
			grid1.irange=[0, imax-imin]
			result.grid1 = grid1
		ENDIF ELSE BEGIN
			IF OBJ_ISA(arg, 'GKVsd') THEN BEGIN
				IF (self.units NE arg.units) THEN BEGIN
					MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
					RETURN, 0
				ENDIF
				result = self -> MakeCopy(/NoValues, /NoErrorBars)
				argValues  = *arg.values
				selfValues = *self.values
				values = selfValues + argValues
				result.values = PTR_NEW(values)
				selfErrors = 0
				IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
				argErrors = 0
				IF PTR_VALID( arg.ErrorBars) THEN  argErrors =  *arg.ErrorBars
				errors = argErrors^2 + selfErrors*2
				errors = SQRT(errors)
				IF(N_ELEMENTS(errors) EQ 1) THEN BEGIN
					nVals = N_ELEMENTS(values)
					errors = REPLICATE(errors, nVals)
				ENDIF
				result.ErrorBars = PTR_NEW(errors)
				vmin = GKVsd_MIN(values, MAX=vmax)
				IF(vmax EQ vmin) THEN vmax = vmin+1
				result.vrange = [vmin, vmax]
				IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
				result.Title = self.Title + "+" + arg.Title
				IF(TypeOf(title) EQ 7) THEN result.title=title
				IF(TypeOf(units) EQ 7) THEN result.units=units
			ENDIF ELSE BEGIN
				MESSAGE, "Argument is unrecognized object", /INFORMATIONAL
				RETURN, 0
			ENDELSE
		ENDELSE 
			
	END
ELSE:		result = self -> GKVsd::Plus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs1D::Plus ***** ;


FUNCTION GKVs1D::Minus, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  Subtraction
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN								; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		IF(iaxis NE 1) THEN BEGIN
			MESSAGE, "Bad Axis ID -- Returning"
			RETURN, 0
		ENDIF
		result = self -> MakeCopy(/noValues)
		argValues = *self.Grid1.values	
		selfValues = *self.values
		resultValues = selfValues - argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "-" + self.Grid1.title + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs1D') EQ 1) THEN BEGIN
		;
		; Take difference between two GKVs1D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "GKVs1D::Minus:  Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			selfValues = (*self.values)[self.Grid1.irange[0]:self.Grid1.irange[1]]
			values = selfValues - argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			argErrors = 0.
			gotErrors = 0
			IF(PTR_VALID(result.ErrorBars)) THEN BEGIN
				argErrors = *result.ErrorBars
				gotErrors = 1
			ENDIF
			selfErrors = 0.
			IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
				selfErrors = (*self.ErrorBars)[self.Grid1.irange[0]:self.Grid1.irange[1]]
				gotErrors = 1
			ENDIF
			errors = SQRT( (selfErrors)^2 + (argErrors)^2 )
			IF( gotErrors EQ 1 ) THEN BEGIN
				PTR_FREE, result.ErrorBars
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
			result.vrange = [vmin, vmax]
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + "-" + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
			grid1=result.grid1
			gridvalues = *grid1.values
			imin=grid1.irange[0]
			imax=grid1.irange[1]
			newValues = gridvalues[imin:imax]
			PTR_FREE, grid1.values
			grid1.values = PTR_NEW(newValues)
			grid1.irange=[0, imax-imin]
			result.grid1 = grid1
		ENDIF ELSE BEGIN
			IF OBJ_ISA(arg, 'GKVsd') THEN BEGIN
				IF (self.units NE arg.units) THEN BEGIN
					MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
					RETURN, 0
				ENDIF
				result = self -> MakeCopy(/NoValues, /NoErrorBars)
				argValues  = *arg.values
				selfValues = *self.values
				values = selfValues - argValues
				result.values = PTR_NEW(values)
				selfErrors = 0
				IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
				argErrors = 0
				IF PTR_VALID( arg.ErrorBars) THEN  argErrors =  *arg.ErrorBars
				errors = argErrors^2 + selfErrors*2
				IF(N_ELEMENTS(errors) EQ 1) THEN BEGIN
					nVals = N_ELEMENTS(values)
					errors = REPLICATE(errors, nVals)
				ENDIF
				errors = SQRT(errors)
				result.ErrorBars = PTR_NEW(errors)
				vmin = GKVsd_MIN(values, MAX=vmax)
				IF(vmax EQ vmin) THEN vmax = vmin+1
				result.vrange = [vmin, vmax]
				IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
				result.Title = self.Title + "-" + arg.Title
				IF(TypeOf(title) EQ 7) THEN result.title=title
				IF(TypeOf(units) EQ 7) THEN result.units=units
			ENDIF ELSE BEGIN
				MESSAGE, "Argument is unrecognized object", /INFORMATIONAL
				RETURN, 0
			ENDELSE
		ENDELSE 
			
	END
ELSE:		result = self -> GKVsd::Minus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs1D::Minus ***** ;


FUNCTION GKVs1D::Times, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  Multiplication
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN								; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		IF(iaxis NE 1) THEN BEGIN
			MESSAGE, "Bad Axis ID -- Returning"
			RETURN, 0
		ENDIF
		result = self -> MakeCopy(/noValues, /NoErrorBars)
		argValues = *self.Grid1.values	
		selfValues = *self.values
		resultValues = selfValues*argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF PTR_VALID(self.ErrorBars) THEN BEGIN
			errorBars = *self.ErrorBars*argValues
			result.ErrorBars = PTR_NEW(errorBars)
		ENDIF
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "*" + self.Grid1.title + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs1D') EQ 1) THEN BEGIN
			
		;
		; Multiply values of two GKVs1D objects
		;
			result = arg -> Interpolate(self)
			; argcopy -> Trash
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "GKVs1D::Times:  Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			selfValues = (*self.values)[self.Grid1.irange[0]:self.Grid1.irange[1]]
			values = argValues*selfValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			argErrors = 0.
			gotErrors = 0
			IF(PTR_VALID(result.ErrorBars)) THEN BEGIN
				argErrors = *result.ErrorBars
				gotErrors = 1
			ENDIF
			selfErrors = 0.
			IF(PTR_VALID(self.ErrorBars)) THEN BEGIN
				selfErrors = (*self.ErrorBars)[self.Grid1.irange[0]:self.Grid1.irange[1]]
				gotErrors = 1
			ENDIF
			errors = SQRT( (selfErrors*argValues)^2 + (argErrors*selfValues)^2 )
			IF(gotErrors EQ 1) THEN BEGIN
				PTR_FREE, result.ErrorBars
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + "*" + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			result.units = "(" + self.units + ")*(" + arg.units + ")"
			IF(TypeOf(units) EQ 7) THEN result.units=units
			grid1=result.grid1
			gridvalues = *grid1.values
			imin=grid1.irange[0]
			imax=grid1.irange[1]
			newValues = gridvalues[imin:imax]
			PTR_FREE, grid1.values
			grid1.values = PTR_NEW(newValues)
			grid1.irange=[0, imax-imin]
			result.grid1 = grid1
		ENDIF ELSE BEGIN
			IF OBJ_ISA(arg, 'GKVsd') THEN BEGIN
				result = self -> MakeCopy(/NoValues, /NoErrorBars)
				argValues  = *arg.values
				selfValues = *self.values
				values = selfValues*argValues
				result.values = PTR_NEW(values)
				selfErrors = 0
				IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
				argErrors = 0
				IF PTR_VALID( arg.ErrorBars) THEN  argErrors =  *arg.ErrorBars
				errors = (selfValues*argErrors)^2 + (argValues*selfErrors)*2
				errors = SQRT(errors)
				result.ErrorBars = PTR_NEW(errors)
				vmin = GKVsd_MIN(values, MAX=vmax)
				IF(vmax EQ vmin) THEN vmax = vmin+1
				result.vrange = [vmin, vmax]
				IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
				result.Title = self.Title + "*" + arg.Title
				IF(TypeOf(title) EQ 7) THEN result.title=title
				result.units = "(" + self.units + ")*(" + arg.units + ")"
				IF(TypeOf(units) EQ 7) THEN result.units=units
			ENDIF ELSE BEGIN
				MESSAGE, "Argument is unrecognized object", /INFORMATIONAL
				RETURN, 0
			ENDELSE
		ENDELSE 
			
	END
ELSE:		result = self -> GKVsd::Times(arg)

ENDCASE

RETURN, result

END ; ***** GKVs1D::Times ***** ;


FUNCTION GKVs1D::Over, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  Division
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN								; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		IF(iaxis NE 1) THEN BEGIN
			MESSAGE, "Bad Axis ID -- Returning"
			RETURN, 0
		ENDIF
		result = self -> MakeCopy(/noValues, /NoErrorBars)
		argValues = *self.Grid1.values	
		selfValues = *self.values
		resultValues = selfValues/argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF PTR_VALID(self.ErrorBars) THEN BEGIN
			errorBars = *self.ErrorBars/argValues
			result.ErrorBars = PTR_NEW(errorBars)
		ENDIF
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "/" + self.Grid1.title + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again						;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs1D') EQ 1) THEN BEGIN
		;
		; Divide values of self (a GKVs1D object) by values of arg 
		;	(another GKVs1D object)
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "GKVs1D::Over:  Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			selfValues = (*self.values)[self.Grid1.irange[0]:self.Grid1.irange[1]]
			values = selfValues/argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			argErrors = 0.
			IF(PTR_VALID(result.ErrorBars)) THEN argErrors = *result.ErrorBars
			selfErrors = 0.
			IF(PTR_VALID(self.ErrorBars)) THEN selfErrors = (*self.ErrorBars)[self.Grid1.irange[0]:self.Grid1.irange[1]]
			errors = values*SQRT( (selfErrors/selfValues)^2 + (argErrors/argValues)^2 )
			IF( N_ELEMENTS(errors) EQ N_ELEMENTS(values) ) THEN BEGIN 
				PTR_FREE, result.ErrorBars
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + "/" + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			result.units = "(" + self.units + ")/(" + arg.units + ")"
			IF(TypeOf(units) EQ 7) THEN result.units=units
			grid1=result.grid1
			gridvalues = *grid1.values
			imin=grid1.irange[0]
			imax=grid1.irange[1]
			newValues = gridvalues[imin:imax]
			PTR_FREE, grid1.values
			grid1.values = PTR_NEW(newValues)
			grid1.irange=[0, imax-imin]
			result.grid1 = grid1
		ENDIF ELSE BEGIN
			IF OBJ_ISA(arg, 'GKVsd') THEN BEGIN
				result = self -> MakeCopy(/NoValues, /NoErrorBars)
				argValues  = *arg.values
				selfValues = *self.values
				values = selfValues/argValues
				result.values = PTR_NEW(values)
				selfErrors = 0
				IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
				argErrors = 0
				IF PTR_VALID( arg.ErrorBars) THEN  argErrors =  *arg.ErrorBars
				errors = (selfErrors/argValues)^2 + (values*argErrors/argValues)^2
				errors = SQRT(errors)
				result.ErrorBars = PTR_NEW(errors)
				vmin = GKVsd_MIN(values, MAX=vmax)
				IF(vmax EQ vmin) THEN vmax = vmin+1
				result.vrange = [vmin, vmax]
				IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
				result.Title = self.Title + "/" + arg.Title
				IF(TypeOf(title) EQ 7) THEN result.title=title
				result.units = "(" + self.units + ")/(" + arg.units + ")"
				IF(TypeOf(units) EQ 7) THEN result.units=units
			ENDIF ELSE BEGIN
				MESSAGE, "Argument is unrecognized object", /INFORMATIONAL
				RETURN, 0
			ENDELSE
		ENDELSE 
			
	END
ELSE:		result = self -> GKVsd::Over(arg)

ENDCASE

RETURN, result

END ; ***** GKVs1D::Over ***** ;


FUNCTION GKVs1D::Slice,_Extra=extra 
;
; Virtual compile-time keywords:	Axis=axisNum, Value=value, Index=index, AtMax=max 
;
; Generic slice routine for GKV data objects.
; Reduces diminsionality by 1 by selecting only
; data with selected value from the selected axis.
;
; If 'AtMax' keyword is set 
; (or if user enters AxisMnemonic= 'max')
; then the maximumvalue over selected axis is retained in the result
;
; Written to be Polymorphic:
;
;	INPUT		OUTPUT
;
;	GKVs1D ->  GKVsd
;	GKVs2D ->  GKVs1D
;	GKVs3D ->  GKVs2D
;	GKVs4D ->  GKVs3D
;
; The _Extra keyword inheritance feature of IDL is used to
; allow axes to be identified by mnemonics
;
;	Written 2/1/00
;	by W.M. Nevins
;
;	Modified 2/12/00 to include Max keyword
;
GKVsClass = STRARR(5)
GKVsClass = ['GKVSD', 'GKVS1D', 'GKVS2D', 'GKVS3D', 'GKVS4D']
;
; Parse 'extra' for compile-time keywords
;
result = GetKeyWord('Axis', extra)				; Equivalent to having
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN axisNum = result	;	Axis = axisNum
IF(TypeOf(result) NE 7) THEN axisNum = result 
result = GetKeyWord('Value', extra)				
IF(TypeOf(result) NE 7) THEN value = result		;	Value = value 
result = GetKeyWord('Index', extra)				
IF(TypeOf(result) NE 7) THEN index = result		;	Index = index 
result = GetKeyWord('AtMax', extra)				
IF(TypeOf(result) NE 7) THEN max = result				;	Max = max
;
;										; on the command line


ndims = self -> NumDims()					; Get number of dimensions
IF (N_ELEMENTS(axisNum) EQ 0) THEN BEGIN			; Try to get valid axis mnemonic
	AxisInfo = self -> GetAxis(extra, /debug) 		;	from _Extra
	FOR iaxis = 1, ndims DO BEGIN				
		axisValue = AxisInfo.(iaxis-1)			; Get value associated with the "axis'i'" tag
		IF(TypeOf(axisValue) NE 7) THEN BEGIN
			value = axisValue			; Get value from corresponding element of 'extra'
			GOTO, GOT_IT
		ENDIF ELSE BEGIN		
			axisValue = STRTRIM(axisValue, 2)	; Strip out leading and trailing blanks
			IF( STRCMP(axisValue, 'MAX', /FOLD_CASE) ) THEN BEGIN
				max=1
				GOTO, GOT_IT
			ENDIF
		ENDELSE
	ENDFOR								; Couldn't match axis... 
	iaxis = ndims							; 	Default to final (probably "time") axis
GOT_IT : axisNum = iaxis
ENDIF
IF(TypeOf(axisNum) EQ 7) THEN axisNum = self -> AxisNumber(axisNum, /debug)
str_axis = STRING(axisNum, FORMAT ='(I1)')	
command_str = 'selectedGrid = self.Grid' + str_axis
ok = EXECUTE(command_str)
value_ptr = selectedGrid.values
values = *value_ptr						; Dereference pointer to array of Grid values
IF(N_ELEMENTS(index) EQ 0) THEN BEGIN				; Index keyword was not set... 
	IF (N_ELEMENTS(value) EQ 0) THEN value = values[0]	; Default value...
	temp = (values - value)^2
	vmin = MIN(temp, index)					; Find index closest to 'value'
ENDIF
IF (index LT 0) THEN INDEX = 0
IF (index GE N_ELEMENTS(values)) THEN index = N_ELEMENTS(values)-1

command_str = "result = {" + GKVsClass[ndims-1] + "}"
ok = EXECUTE(command_str)					; Create Class Structure of appropriate type

aGKVsd = {GKVsd}
nTagsGKVsd = N_TAGS(aGKVsd)
PTR_FREE, aGKVsd.values
PTR_FREE, aGKVsd.ErrorBars
FOR i=0, nTAgsGKVsd-1 DO 		$
	result.(i) = self.(i)					; Populate GKVsd field of result

iResultAxis = 1
iSelfAxis = 1 				
WHILE iResultAxis LE ndims-1 DO BEGIN			; Populate Grid structures of result
	IF(axisNum EQ iResultAxis) THEN	$		; Skip self.Grid'axisNum'
		iSelfAxis=iSelfAxis + 1
	str_Result	= STRING(iResultAxis, FORMAT ='(I1)')	
	str_Self 	= STRING(  iSelfAxis, FORMAT ='(I1)')	
	command_str = "result.Grid" + str_Result + " = self.Grid" + str_Self
	ok = EXECUTE(command_str)
	command_str = "result.Grid" + str_Result + ".values = PTR_NEW(*self.Grid" + str_Self + ".values)"
	ok = EXECUTE(command_str)
	iResultAxis = iResultAxis + 1
	ISelfAxis = iSelfAxis + 1
ENDWHILE
info = SIZE(*self.values)
;
; check for error bars
;
errorBars = 0
IF(PTR_VALID(self.errorBars)) THEN errorBars = 1
;
; Check for 'AtMax' keyword
;
IF(KEYWORD_SET(Max)) THEN BEGIN				; MAX keyword was set, so must search for maximum
irange = selectedGrid.irange				; at each element of selected axis within irange
IF(ndims GT 3) THEN BEGIN					; (which may not correspond to current signal window...)
	MESSAGE, 'AtMax keyword not yet implimented for ndims > 3', /informational
	RETURN, 0
ENDIF
imin=irange[0]
imax=irange[1]
xvalues = (*selectedGrid.values)[imin:imax]
CASE ndims OF							; Create 'temp', xtemp and (if appropriarte) erroor arrays of appropriate size
	1 :	BEGIN
			temp = 0.0
			xtemp= 0.0
			values = (*self.values)[imin:imax]
			IF(errorBars) THEN errors = (*self.ErrorBars)[imin:imax]
		END
				
	2 :	CASE axisNum OF
		1 :	BEGIN
				temp = FLTARR(info[2])
				xtemp= FLTARR(info[2])
				values = (*self.values)[imin:imax,*]
				etemp= FLTARR(info[2]) 
				IF(errorBars) THEN errors = (*self.ErrorBars)[imin:imax,*]
			END
		2 :	BEGIN
				temp = FLTARR(info[1])
				xtemp= FLTARR(info[1])
				values = (*self.values)[*,imin:imax]
				IF(errorBars) THEN errors = (*self.ErrorBars)[*,imin:imax]
				etemp= FLTARR(info[1]) 
			END
		ENDCASE
	3 :	CASE axisNum OF
		1 :	BEGIN
				temp = FLTARR(info[2], info[3])
				xtemp= FLTARR(info[2], info[3])
				values = (*self.values)[imin:imax,*,*]
				IF(errorBars) THEN errors = (*self.ErrorBars)[imin:imax,*,*]
				etemp= FLTARR(info[2], info[3]) 
			END
		2 :	BEGIN
				temp = FLTARR(info[1], info[3])
				xtemp= FLTARR(info[1], info[3])
				values = (*self.values)[*,imin:imax,*]
				IF(errorBars) THEN errors = (*self.ErrorBars)[*,imin:imax,*]
				etemp= FLTARR(info[1], info[3]) 
			END
		3 :	BEGIN
				temp = FLTARR(info[1], info[2])
				xtemp= FLTARR(info[1], info[2])
				values = (*self.values)[*,*,imin:imax]
				IF(errorBars) THEN errors = (*self.ErrorBars)[*,*,imin:imax]
				etemp= FLTARR(info[1], info[2]) 
			END
		ENDCASE
ENDCASE
CASE ndims OF							; set values in sliced object & value of independent variable where max occurs
	1 :	BEGIN
			temp = MAX(values, maxIndex)
			xtemp= xvalues[maxIndex]
			IF(errorBars) THEN eTemp = errors[maxIndex]
		END
				
	2 :	CASE axisNum OF
		1 :	FOR j=0, info[2]-1 DO BEGIN
				temp[j] = Max(values[*,j], maxIndex)
				xtemp[j]= xvalues[maxIndex]
				IF(errorBars) THEN eTemp[j] = errors[maxIndex,j]
			ENDFOR
		2 :	FOR j=0, info[1]-1 DO BEGIN
				temp[j] = Max(values[j,*], maxIndex)
				xtemp[j]= xvalues[maxIndex]
				IF(errorBars) THEN eTemp[j] = errors[j,maxIndex]
			ENDFOR
		ENDCASE
	3 :	CASE axisNum OF
		1 :	FOR i=0, info[2]-1 DO BEGIN
				FOR j=0, info[3]-1 DO BEGIN
					temp[i,j] = Max(values[*,i,j], maxIndex)
					xtemp[i,j]= xvalues[maxIndex]
					IF(errorBars) THEN eTemp[i,j] = errors[maxIndex,i,j]
				ENDFOR
			ENDFOR
		2 :	FOR i=0, info[1]-1 DO BEGIN
				FOR j=0, info[3]-1 DO BEGIN
					temp[i,j] = Max(values[i,*,j], maxIndex)
					xtemp[i,j]= xvalues[maxIndex]
					IF(errorBars) THEN eTemp[i,j] = errors[i,maxIndex,j]
				ENDFOR
			ENDFOR
		3 :	FOR i=0, info[1]-1 DO BEGIN
				FOR j=0, info[2]-1 DO BEGIN
					temp[i,j] = Max(values[i,j,*], maxIndex)
					xtemp[i,j]= xvalues[maxIndex]
					IF(errorBars) THEN eTemp[i,j] = errors[i,j,maxIndex]
				ENDFOR
			ENDFOR
		ENDCASE
ENDCASE

ENDIF ELSE BEGIN						; MAX keyword not set, so just select using 'index'

CASE ndims OF							; set values in sliced object
	1 :		temp = (*self.values)[index]
				
	2 :	CASE axisNum OF
		1 :	temp = REFORM((*self.values)[index, *], info[2])
		2 :	temp = REFORM((*self.values)[*, index], info[1])
		ENDCASE
	3 :	CASE axisNum OF
		1 :	temp = REFORM((*self.values)[index, *, *], info[2], info[3])
		2 :	temp = REFORM((*self.values)[*, index, *], info[1], info[3])
		3 :	temp = REFORM((*self.values)[*, *, index], info[1], info[2])
		ENDCASE
	4 :	CASE axisNum OF
		1 :	temp = REFORM((*self.values)[index, *, *, *], info[2], info[3], info[4])
		2 :	temp = REFORM((*self.values)[*, index, *, *], info[1], info[3], info[4])
		3 :	temp = REFORM((*self.values)[*, *, index, *], info[1], info[2], info[4])
		4 :	temp = REFORM((*self.values)[*, *, *, index], info[1], info[2], info[3])
		ENDCASE
ENDCASE
;
; get errorbars is they are there
;
IF(errorBars) THEN BEGIN
CASE ndims OF							; set values in sliced object
	1 :		errors = (*self.errorBars)[index]
				
	2 :	CASE axisNum OF
		1 :	errors = REFORM((*self.errorBars)[index, *], info[2])
		2 :	errors = REFORM((*self.errorBars)[*, index], info[1])
		ENDCASE
	3 :	CASE axisNum OF
		1 :	errors = REFORM((*self.errorBars)[index, *, *], info[2], info[3])
		2 :	errors = REFORM((*self.errorBars)[*, index, *], info[1], info[3])
		3 :	errors = REFORM((*self.errorBars)[*, *, index], info[1], info[2])
		ENDCASE
	4 :	CASE axisNum OF
		1 :	errors = REFORM((*self.errorBars)[index, *, *, *], info[2], info[3], info[4])
		2 :	errors = REFORM((*self.errorBars)[*, index, *, *], info[1], info[3], info[4])
		3 :	errors = REFORM((*self.errorBars)[*, *, index, *], info[1], info[2], info[4])
		4 :	errors = REFORM((*self.errorBars)[*, *, *, index], info[1], info[2], info[3])
		ENDCASE
ENDCASE
ENDIF

ENDELSE
result.values = PTR_NEW(temp)
IF(errorbars) THEN result.errorBars = PTR_NEW(errors)
value = values[index]
IF(KEYWORD_SET(Max))	THEN valueStr = "Max"	$
				ELSE valueStr = STRING(value, FORMAT='(G10.3)')

result.Title=self.Title

plotAxis = INTARR(ndims+1)
plotAxis[axisNum] = -index	; ***NOTE*** if INDEX=0, then plotAxis won't do anything, and INDICES of result will have an extra '*'
indices = Self -> IndexString(plotAxis, /pretty)
;
IF(index EQ 0) THEN BEGIN  ; an inelegant patch to fix the INDEX=0 problem
	indexStr = selectedGrid.title + '=' + valueStr
	indexStr = STRCOMPRESS(indexStr, /REMOVE_ALL)
	indices = *self.indices
	info = SIZE(indices)
	Nindices = info[1]
	iaxis=0
	FOR i=0, nindices-1 DO BEGIN
	IF(indices[i] EQ '*') THEN BEGIN
		iaxis=iaxis+1
		IF(iaxis EQ axisNum) THEN indices[i] = indexStr
	ENDIF
	ENDFOR
ENDIF					; end of inelegant patch...
;
IF(KEYWORD_SET(Max)) THEN indices = Self -> IndexMax(axisNum)
result.Indices = PTR_NEW(Indices)

result.mnemonic =self.mnemonic

vmin = MIN(*result.values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1
result.vrange = [vmin, vmax]

command_str = "sliced = OBJ_NEW('" + GKVsClass[ndims-1] + "', result)"
ok = EXECUTE(command_str)

IF(KEYWORD_SET(MAX)) THEN BEGIN
	result=0
	IF(TypeOf(extra) EQ 8) THEN result=GetKeyWord('MaxLocation', extra)
	IF(KEYWORD_SET(result)) THEN BEGIN
		maxloc = sliced -> MakeCopy(/noValues)
		maxloc.values = PTR_NEW(xtemp)
		xmax = MAX(xtemp, MIN=xmin)
		maxloc.vrange = [xmin, xmax]
		maxloc.title = selectedGrid.title + '!Imax!N'
		maxloc.mnemonic = selectedGrid.mnemonic + '_max'
		maxloc.units = selectedGrid.units
		outputStr = {Name:"slice", slice:sliced, maxLocation:maxloc}
		RETURN, outputStr
	ENDIF
ENDIF

RETURN, sliced
END ; ***** GKVs1D::Slice ***** ;
	

FUNCTION GKVs1D::Delta, _Extra=extra 
;
; Virtual compile-time keywords:	Axis=axisNum, Range=Range, Irange=irange
;
; Generic routine to compute average and deviations of GKV objects
; over selected axis.  Results are returned in a structure with tags
; "Avg" (for the average) and "Delta" (for the deviations).
;
; Written to be Polymorphic:
;
;	INPUT		OUTPUT
;
; The _Extra keyword inheritance feature of IDL is used to
; allow axes to be identified by mnemonics
;
;	Written 2/11/00
;	by W.M. Nevins
;
; First, parse 'extra' for keywords which we wish to intercept 
; (this is the moral equivalent of defining them on the PRO line...
;  but prevents IDL from enforcing keyword abreviation, which would
;  interfer with use of mnemonics as 'run-time' keywords)
;
result = GetKeyWord('Axis', extra)			; Equivalent to having
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN axisNum = result	;	Axis = axisNum 
IF(TypeOf(result) NE 7) THEN axisNum = result
result = GetKeyWord('Range', extra)				
IF(TypeOf(result) NE 7) THEN range = result		;	Range = range 
result = GetKeyWord('Irange', extra)				
IF(TypeOf(result) NE 7) THEN irange = result		;	Irange = irange 
;
;							; on the command line


GKVsClass = STRARR(5)
GKVsClass = ['GKVSD', 'GKVS1D', 'GKVS2D', 'GKVS3D', 'GKVS4D']
ndims = self -> NumDims()
iaxis=axisNum
IF ( ndims LE 0 ) THEN BEGIN
	MESSAGE, 'Avg undefined for GKVsd Class', /Informational
	RETURN, 0
ENDIF
IF(TypeOF(axisNum) EQ 7) THEN BEGIN				; axisNum is a "string".  Assume it is an axis mnemonic
	iaxis = self -> AxisNumber(axisNum,/debug)		; 	try to match grid mnemonics
	IF(iaxis NE 0) THEN GOTO, GOT_IT
	MESSAGE, 'Invalid axis mnemonic', /Informational
	RETURN, 0
ENDIF
IF (N_ELEMENTS(axisNum) EQ 0) THEN BEGIN			; Try to get valid axis mnemonic
	AxisInfo = self -> GetAxis(extra, /debug) 		;	from _Extra
	FOR iaxis = 1, ndims DO BEGIN				
		axisValue = AxisInfo.(iaxis-1)			; Get value associated with the "axis'i'" tag
		IF(TypeOf(axisValue) NE 7) THEN BEGIN
			range = axisValue			;  Assume value is desired range
			GOTO, GOT_IT
		ENDIF 
	ENDFOR							; Couldn't match axis... 
	iaxis = 1
	IF(ndims EQ 1) THEN GOTO, GOT_IT			; If only one dimension, default to single available axis,
	MESSAGE, 'No valid axis ID', /Informational	; otherwise send error message and return 0.
	RETURN, 0
ENDIF
GOT_IT : axisNum = iaxis
IF(axisNum GT ndims) THEN BEGIN				; Invalid axis number...
	MESSAGE, 'Invalid axis ID', /Informational	;	send error message and return 0
	RETURN, 0
ENDIF
str_axis = STRING(axisNum, FORMAT ='(I1)')	
command_str = 'Grid = self.Grid' + str_axis
ok = EXECUTE(command_str)				; Get selected Grid structure
gridValues_ptr = Grid.values
GridValues = *gridValues_ptr				; Dereference pointer to array of Grid values
nGridPoints = N_ELEMENTS(GridValues)
boundary = Grid.boundary
uniform = Grid.uniform
iGridRange = Grid.irange
IF(N_ELEMENTS(Irange) NE 2) THEN BEGIN		; Irange keyword was not set... 
Irange = INTARR(2)
CASE N_ELEMENTS(Range) OF
	0:	Irange = iGridRange				; Use signal window for default range
	2:	BEGIN						; take specified values of 'Range' as upper and lower limits
			FOR i=0,1 DO BEGIN
				temp = (GridValues - Range[i])^2
				rmin = MIN(temp, index)
				Irange[i] = index
			ENDFOR
		IF(Irange[0] LT 0) THEN Irange[0]=0
		IF(Irange[1] GE nGridPoints) THEN Irange[1] = nGridPoints - 1
		END
	ELSE:	BEGIN
			Message, 'Syntax:  result = Obj -> DELTA( ..., Range =[lower, upper])', /Informational
			RETURN, 0
	ENDELSE
ENDCASE
ENDIF
values = *self.values
imin = irange[0]
imax = irange[1]
range = [GridValues[imin], GridValues[imax]]
gridLength = range[1] - range[0]
npoints = imax - imin + 1
info = SIZE(values)
;
; Check for periodic grid
;
gridInfo = SIZE(GridValues)
maxPoints= gridInfo[1]
Check = maxPoints EQ npoints				; Are all grid points included?
Check = Check AND uniform					; Is Grid uniform?
boundary = STRUPCASE(BOUNDARY)
Check = Check AND (boundary EQ 'PERIODIC')		; Are Grid boundary conditions periodic?
IF(Check) THEN BEGIN						; Irange covers entire periodic grid, so
dGrid = GridValues[imin+1] - GridValues[imin]	; 	compute endpoint contributions assuming periodic BC's
gridLength = gridLength + dGrid
CASE ndims OF							

	1:	avgValues = (values[imin] + values[imax])*dGrid
	2:	CASE axisNum OF	
		1 :	avgValues = (values[imin, *] + values[imax, *])*dGrid
						
		2 :	avgValues = (values[*, imin] + values[*, imax])*dGrid
		ENDCASE
	3 :	CASE axisNum OF
		1 :	avgValues	= (values[imin, *, *] + values[imax, *, *])*dGrid
						
		2 :	avgValues	= (values[*, imin, *] + values[*, imax, *])*dGrid
						
		3 :	avgValues	= (values[*, *, imin] + values[*, *, imax])*dGrid
		ENDCASE
	4 :	CASE axisNum OF
		1 :	avgValues 	= (values[imin, *, *, *] + values[imax, *, *, *])*dGrid
						
		2 :	avgValues	= (values[*, imin, *, *] + values[*, imax, *, *])*dGrid
						
		3 :	avgValues	= (values[*, *, imin, *] + values[*, *, imax, *])*dGrid
						
		4 :	avgValues	= (values[*, *, *, imin] + values[*, *, *, imax])*dGrid
		ENDCASE
ENDCASE
ENDIF ELSE BEGIN
CASE ndims OF							; Compute endpoint contributions, and store in avgValues array

	1:	avgValues =	0.5*	(	values[imin]*(GridValues[imin+1]-GridValues[imin])	$
						+	values[imax]*(GridValues[imax]-GridValues[imax-1])	)	
	2:	CASE axisNum OF	
		1 :	avgValues =0.5*	( values[imin, *]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[imax, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		2 :	avgValues =0.5*	( values[*, imin]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[*, imax]*(GridValues[imax]-GridValues[imax-1])	)
		ENDCASE
	3 :	CASE axisNum OF
		1 :	avgValues	=0.5*	( values[imin, *, *]*(GridValues[imin+1]-GridValues[imin]) $
						+ values[imax, *, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		2 :	avgValues	=0.5*	( values[*, imin, *]*(GridValues[imin+1]-GridValues[imin]) $
						+ values[*, imax, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		3 :	avgValues	=0.5*	( values[*, *, imin]*(GridValues[imin+1]-GridValues[imin]) $
						+ values[*, *, imax]*(GridValues[imax]-GridValues[imax-1])	)
		ENDCASE
	4 :	CASE axisNum OF
		1 :	avgValues=0.5*	( values[imin, *, *, *]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[imax, *, *, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		2 :	avgValues =0.5*	( values[*, imin, *, *]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[*, imax, *, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		3 :	avgValues =0.5*	( values[*, *, imin, *]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[*, *, imax, *]*(GridValues[imax]-GridValues[imax-1])	)
						
		4 :	avgValues =0.5*	( values[*, *, *, imin]*(GridValues[imin+1]-GridValues[imin]) 	$
						+ values[*, *, *, imax]*(GridValues[imax]-GridValues[imax-1])	)
		ENDCASE
ENDCASE
ENDELSE
IF uniform THEN BEGIN
dGrid = GridValues[imin+1] - GridValues[imin]	 
CASE ndims OF							; Compute interior contributions using Simpson's rule

	1:		avgValues = avgValues + TOTAL(values[imin+1:imax-1])*dGrid
	
	2:	CASE axisNum OF	
		1 :	avgValues = avgValues + TOTAL(values[imin+1:imax-1,*], 1)*dGrid
						
		2 :	avgValues = avgValues + TOTAL(values[*,imin+1:imax-1], 2)*dGrid
		ENDCASE
	3:	CASE axisNum OF
		1 :	avgValues = avgValues + TOTAL(values[imin+1:imax-1,*,*], 1)*dGrid
						
		2 :	avgValues = avgValues + TOTAL(values[*,imin+1:imax-1,*], 2)*dGrid
						
		3 :	avgValues = avgValues + TOTAL(values[*,*,imin+1:imax-1], 3)*dGrid
		ENDCASE
	4:	CASE axisNum OF
		1 :	avgValues = avgValues + TOTAL(values[imin+1:imax-1,*,*,*], 1)*dGrid
						
		2 :	avgValues = avgValues + TOTAL(values[*,imin+1:imax-1,*,*], 2)*dGrid
						
		3 :	avgValues = avgValues + TOTAL(values[*,*,imin+1:imax-1,*], 3)*dGrid
						
		4 :	avgValues = avgValues + TOTAL(values[*,*,*,imin+1:imax-1], 4)*dGrid
		ENDCASE
ENDCASE
ENDIF ELSE BEGIN
FOR i=imin+1, imax-1 DO BEGIN
CASE ndims OF							; Compute interior contributions using Simpson's rule

	1:		avgValues = avgValues + 0.5*values[i]*(GridValues[i+1]-GridValues[i-1]) 
	
	2:	CASE axisNum OF	
		1 :	avgValues = avgValues + 0.5*values[i,*]*(GridValues[i+1]-GridValues[i-1])
						
		2 :	avgValues = avgValues + 0.5*values[*,i]*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
	3:	CASE axisNum OF
		1 :	avgValues = avgValues + 0.5*values[i,*,*]*(GridValues[I+1]-GridValues[i-1])
						
		2 :	avgValues = avgValues + 0.5*values[*,i,*]*(GridValues[i+1]-GridValues[i-1])
						
		3 :	avgValues = avgValues + 0.5*values[*,*,i]*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
	4:	CASE axisNum OF
		1 :	avgValues = avgValues + 0.5*values[i,*,*,*]*(GridValues[i+1]-GridValues[i-1])
						
		2 :	avgValues = avgValues + 0.5*values[*,i,*,*]*(GridValues[i+1]-GridValues[i-1])
						
		3 :	avgValues = avgValues + 0.5*values[*,*,i,*]*(GridValues[i+1]-GridValues[i-1])
						
		4 :	avgValues = avgValues + 0.5*values[*,*,*,i]*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
ENDCASE
ENDFOR
ENDELSE
;
CASE ndims OF							; Reform avgValues array
	1:								; Average is scalar
	2:	CASE axisNum OF
		1 :	avgValues = REFORM(avgValues, info[2], /OVERWRITE)	
		2 :	avgValues = REFORM(avgValues, info[1], /OVERWRITE)
		ENDCASE
	3 :	CASE axisNum OF
		1 :	avgValues = REFORM(avgValues, info[2], info[3], /OVERWRITE)
		2 :	avgValues = REFORM(avgValues, info[1], info[3], /OVERWRITE)
		3 :	avgValues = REFORM(avgValues, info[1], info[2], /OVERWRITE)
		ENDCASE
	4 :	CASE axisNum OF
		1 :	avgValues = REFORM(avgValues, info[2], info[3], info[4], /OVERWRITE)
		2 :	avgValues = REFORM(avgValues, info[1], info[3], info[4], /OVERWRITE)
		3 :	avgValues = REFORM(avgValues, info[1], info[2], info[4], /OVERWRITE)
		4 :	avgValues = REFORM(avgValues, info[1], info[2], info[3], /OVERWRITE)
		ENDCASE
ENDCASE
avgValues = avgValues/gridLength
; 
; Compute deviations from average
;
CASE ndims OF							; Form array to hold deviations from average over selected axis
	1:		delta = FLTARR(npoints)		
	2:	CASE axisNum OF
		1 :	delta = FLTARR(npoints, info[2])	
		2 :	delta = FLTARR(info[1], npoints)
		ENDCASE
	3 :	CASE axisNum OF
		1 :	delta = FLTARR(npoints, info[2], info[3])			
		2 :	delta = FLTARR(info[1], npoints, info[3])
		3 :	delta = FLTARR(info[1], info[2], npoints)
		ENDCASE
	4 :	CASE axisNum OF
		1 :	delta = FLTARR(npoints, info[2], info[3], info[4])
		2 :	delta = FLTARR(info[1], npoints, info[3], info[4])
		3 :	delta = FLTARR(info[1], info[2], npoints, info[4])
		4 :	delta = FLTARR(info[1], info[2], info[3], npoints)
		ENDCASE
ENDCASE
j=0
FOR i=imin, imax DO BEGIN
CASE ndims OF							
	1 :		delta(j) = values(i) - avgValues
	2:	CASE axisNum OF
		1 :	delta(j,*) = values(i,*) - avgValues
		2 :	delta(*,j) = values(*,i) - avgValues
		ENDCASE
	3:	CASE axisNum OF
		1 :	delta(j,*,*) = values(i,*,*) - avgValues
		2 : 	delta(*,j,*) = values(*,i,*) - avgValues
		3 :	delta(*,*,j) = values(*,*,i) - avgValues
		ENDCASE
	4:	CASE axisNum OF
		1 :	delta(j,*,*,*) = values(i,*,*,*) - avgValues
		2 : 	delta(*,j,*,*) = values(*,i,*,*) - avgValues
		3 :	delta(*,*,j,*) = values(*,*,i,*) - avgValues
		4 :	delta(*,*,*,j) = values(*,*,*,i) - avgValues
		ENDCASE
ENDCASE
j=j+1
ENDFOR
;
; Store avgValues, delta in appropriate GKVs... structures
; along with associated info
;
command_str =   "avgClass = {" + GKVsClass[ndims-1] + "}"
ok = EXECUTE(command_str)					; Create Class Structure for AVG
command_str = "deltaClass = {" + GKVsClass[ndims]   + "}"
ok = EXECUTE(command_str)					; Create Class Strucdture for DELTA 
aGKVsd = {GKVsd}
nTagsGKVsd = N_TAGS(aGKVsd)
PTR_FREE, aGKVsd.values
PTR_FREE, aGKVsd.ErrorBars
FOR i=0, nTAgsGKVsd-1 DO BEGIN
	  avgClass.(i) = self.(i)				; Populate GKVsd field of avgClass
	deltaClass.(i) = self.(i)				; Populate GKVsd field of deltaClass
ENDFOR
deltaIndices = *self.Indices
deltaClass.indices = PTR_NEW(deltaIndices)
deltaClass.values = PTR_NEW(delta)
vmin = GKVsd_MIN(delta, MAX=vmax)
deltaClass.vrange = [vmin,vmax]
str_axis = STRING(axisNum, FORMAT='(I1)')
command_str = "axistitle = self.Grid" + str_axis + ".title"
ok = EXECUTE(command_str)					; Get axis title 
deltaClass.title = "!4d!X" + self.title 

avgClass.values = PTR_NEW(avgValues)
vmin = GKVsd_MIN(avgValues,MAX=vamx)
avgClass.vrange = [vmin, vmax]
avgClass.title = "!12<!X" + self.title + "!12>!X!I" + axistitle + "!N!X"
indices = self -> IndexRemove(axisNum)
avgClass.Indices = PTR_NEW(indices)
 
iavgAxis = 1
FOR iSelfAxis = 1, ndims DO BEGIN			; Populate Grid structures of avgClass and deltaClass
	str_Self 	= STRING(  iSelfAxis, FORMAT ='(I1)')	
	command_str = "Grid = self.Grid" + str_Self
	ok = EXECUTE(command_str)				; Get iSelfAxisth Grid structure from self
	Gvalues = *Grid.values				; make new copy of grid values
	IF(axisNum EQ iSelfAxis) THEN BEGIN		; Must correct Grid.values, .irange, .range for axis 'axisNum'
		Gvalues=GridValues[imin:imax]
		Grid.irange = [0, imax-imin]
		Grid.range  = range
	ENDIF
	Grid.values = PTR_NEW(Gvalues)			; Create a new pointer to Grid values.
	command_str = "deltaClass.Grid" + str_self + " = Grid"
	ok = EXECUTE(command_str)		
	IF(axisNum NE iselfAxis) THEN BEGIN		; Skip self.Grid'axisNum'
		aGridValues = *Grid.values
		Grid.values = PTR_NEW(aGridValues)
		str_avg	= STRING(iavgAxis, FORMAT ='(I1)')
		command_str = "avgClass.Grid" + str_avg + " = Grid"
		ok = EXECUTE(command_str)
		iavgAxis = iavgAxis + 1
	ENDIF
ENDFOR
Delta_Obj = OBJ_NEW(GKVsClass[ndims], deltaClass)
Avg_Obj   = OBJ_NEW(GKVsClass[ndims-1], avgClass)
RETURN, {Avg:Avg_Obj, Delta:Delta_Obj}
END ; ****** GKVs1D::Delta ****** ;


FUNCTION GKVs1D::Moments, Axis=axisNum, Range=Range, Irange=irange, Avg=Avg, _Extra=extra
;
; Generic routine to compute moments of GKV objects
; by integrating over selected axis, thereby
; Reducing the  diminsionality by 1.
;
; Written to be Polymorphic:
;
;	INPUT		OUTPUT
;;
; The _Extra keyword inheritance feature of IDL is used to
; allow axes to be identified by mnemonics
;
;	Written 2/9/00
;	by W.M. Nevins
;
GKVsClass = STRARR(5)
GKVsClass = ['GKVSD', 'GKVS1D', 'GKVS2D', 'GKVS3D', 'GKVS4D']
ndims = self -> NumDims()
IF ( ndims LE 0 ) THEN BEGIN
	MESSAGE, 'Moments undefined for GKVsd Class', /Informational
	RETURN, 0
ENDIF
IF(TypeOF(axisNum) EQ 7) THEN BEGIN			; axisNum is a "string".  Assume it is an axis mnemonic
	axisName = STRUPCASE(axisNum)
	FOR iaxis = 1, ndims DO BEGIN				
		str_axis = STRING(iaxis, FORMAT ='(I1)')
		command_str = 'mnemonic = self.Grid' + str_axis + '.mnemonic'
		ok = EXECUTE(command_str)			; Get mnemonic associated with the iaxis^th axis
		mnemonic = STRUPCASE(mnemonic)		; Force all characters to be upper case
		IF(mnemonic EQ axisName) THEN GOTO, GOT_IT	; Found a match
	ENDFOR							; Couldn't match axis... 
ENDIF
IF (N_ELEMENTS(axisNum) EQ 0) THEN BEGIN		; Try to get valid axis mnemonic
	numTags = N_TAGS(extra)				;	from _Extra
	IF(TypeOf(extra) EQ 8) THEN tagName = TAG_NAMES(extra)
	IF (numTags EQ 0) THEN BEGIN
		iaxis = 1
		IF(ndims EQ 1) THEN GOTO, GOT_IT	; If only one dimension, default to single available axis,
		MESSAGE, 'No valid axis ID', /Informational	; otherwise send error message and return 0.
		RETURN, 0
	ENDIF
	FOR iaxis = 1, ndims DO BEGIN				
		str_axis = STRING(iaxis, FORMAT ='(I1)')
		command_str = 'mnemonic = self.Grid' + str_axis + '.mnemonic'
		ok = EXECUTE(command_str)			; Get mnemonic associated with the iaxis^th axis
		mnemonic = STRUPCASE(mnemonic)		; Force all characters to be upper case
		FOR iTag = 0, numTags-1 DO BEGIN
			tname = STRUPCASE(tagName(itag))
			IF ( tname EQ mnemonic ) THEN BEGIN
				Range = extra.(itag)		; Get Range of values from corresponding element of 'extra'
				GOTO, GOT_IT
			ENDIF
		ENDFOR
	ENDFOR							; Couldn't match axis... 
	MESSAGE, 'No valid axis ID',/Informational	; 	send error message and return 0.
	RETURN, 0
GOT_IT : axisNum = iaxis
ENDIF
IF(axisNum GT ndims) THEN BEGIN				; Invalid axis number...
	MESSAGE, 'Invalid axis ID', /Informational	;	send error message and return 0
	RETURN, 0
ENDIF
str_axis = STRING(axisNum, FORMAT ='(I1)')	
command_str = 'gridValues_ptr = self.Grid' + str_axis + '.values'
ok = EXECUTE(command_str)
GridValues = *gridValues_ptr					; Dereference pointer to array of Grid values
IF(N_ELEMENTS(Irange) NE 2) THEN BEGIN		; Irange keyword was not set... 
Irange = INTARR(2)
CASE N_ELEMENTS(Range) OF
	0:	BEGIN							; Use signal window as default range
			command_str = 'Irange = self.Grid' + str_axis + '.irange'
			ok = EXECUTE(command_str)
		END
	1:	BEGIN							; take specified value of 'Range' as one limit, default other limit.
			command_str = 'Irange = self.Grid' + str_axis + '.irange'
			ok = EXECUTE(command_str)
			temp = (GridValues - Range)^2
			rmin = MIN(temp, index)
			IF(index LE Irange[1]) THEN BEGIN
				Irange[0] = index 
			ENDIF ELSE BEGIN
				Irange[1] = index
			ENDELSE
		END
	2:	BEGIN							; take specified values of 'Range' as upper and lower limits
			FOR i=0,1 DO BEGIN
				temp = (GridValues - Range[i])^2
				rmin = MIN(temp, index)
				Irange[i] = index
			ENDFOR
		END
	ELSE:	BEGIN
			Message, 'Syntax:  Result = Obj -> Moments( ..., Range =[lower, upper])', /Informational
			RETURN, 0
	ENDELSE
ENDCASE
ENDIF
values = *self.values
imin = irange[0]
imax = irange[1]
range = [GridValues[imin], GridValues[imax]]
npoints = imax - imin + 1
info = SIZE(values)
CASE ndims OF							; Form array to hold first three moments over selected axis
	1:	moments = FLTARR(3)
	2:	CASE axisNum OF
		1 :	moments = FLTARR(3, info[2])
		2 :	moments = FLTARR(3, info[1])
		ENDCASE
	3 :	CASE axisNum OF
		1 :	moments = FLTARR(3, info[2], info[3])
		2 :	moments = FLTARR(3, info[1], info[3])
		3 :	moments = FLTARR(3, info[1], info[2])
		ENDCASE
	4 :	CASE axisNum OF
		1 :	moments = FLTARR(3, info[2], info[3], info[4])
		2 :	moments = FLTARR(3, info[1], info[3], info[4])
		3 :	moments = FLTARR(3, info[1], info[2], info[4])
		4 :	moments = FLTARR(3, info[1], info[2], info[3])
		ENDCASE
ENDCASE
FOR m=0,2 DO BEGIN
CASE ndims OF							; Compute endpoint contributions, and store in moments array

	1:	moments[m] 		=0.5*	( values[imin]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin])			$
							+ values[imax]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])			)
	2:	CASE axisNum OF	
		1 :	moments[m,*] 	=0.5*	( values[imin, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 		$
							+ values[imax, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])		)
						
		2 :	moments[m,*] 	=0.5*	( values[*, imin]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 		$
							+ values[*, imax]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])		)
		ENDCASE
	3 :	CASE axisNum OF
		1 :	moments[m,*,*]	=0.5*	( values[imin, *, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[imax, *, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])	)
						
		2 :	moments[m,*,*]	=0.5*	( values[*, imin, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[*, imax, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])		)
						
		3 :	moments[m,*,*]	=0.5*	( values[*, *, imin]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[*, *, imax]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])		)
		ENDCASE
	4 :	CASE axisNum OF
		1 :	moments[m,*,*,*]=0.5*	( values[imin, *, *, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[imax, *, *, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])	)
						
		2 :	moments[m,*,*,*]=0.5*	( values[*, imin, *, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[*, imax, *, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])	)
						
		3 :	moments[m,*,*,*]=0.5*	( values[*, *, imin, *]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[*, *, imax, *]*(GridValues[imax]^m)*(GridValues[imax]-GridValues[imax-1])	)
						
		4 :	moments[m,*,*,*]=0.5*	( values[*, *, *, imin]*(GridValues[imin]^m)*(GridValues[imin+1]-GridValues[imin]) 	$
							+ values[*, *, *, imax]*(GridValues[imin]^m)*(GridValues[imax]-GridValues[imax-1])	)
		ENDCASE
ENDCASE
FOR i=imin+1, imax-1 DO BEGIN
CASE ndims OF							; Compute interior contributions using Simpson's rule
	1:	moments[m] = moments[m] +	$
					0.5*values[i]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
	2:	CASE axisNum OF	
		1 :	moments[m,*] = moments[m,*] +	$
 					0.5*values[i,*]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
						
		2 :	moments[m,*] = 	moments[m,*] +	$
					0.5*values[*,i]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
	3:	CASE axisNum OF
		1 :	moments[m,*,*] = moments[m,*,*] +	$
					0.5*values[i,*,*]*(GridValues[i]^m)*(GridValues[I+1]-GridValues[i-1])
						
		2 :	moments[m,*,*] = moments[m,*,*] +	$
					0.5*values[*,i,*]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
						
		3 :	moments[m,*,*] = moments[m,*,*] +	$
					0.5* values[*,*,i]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
	4:	CASE axisNum OF
		1 :	moments[m,*,*,*] = moments[m,*,*,*] +		$
					0.5*values[i,*,*,*]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
						
		2 :	moments[m,*,*,*] = moments[m,*,*,*] +		$
					0.5*values[*,i,*,*]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
						
		3 :	moments[m,*,*,*] = moments[m,*,*,*] +		$
					0.5*values[*,*,i,*]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
						
		4 :	moments[m,*,*,*] = moments[m,*,*,*] +		$
					0.5*values[*,*,*,i]*(GridValues[i]^m)*(GridValues[i+1]-GridValues[i-1])
		ENDCASE
ENDCASE
ENDFOR
CASE ndims OF
	1:	moments[2] = moments[2] - moments[1]^2/moments[0]
	2:	moments[2,*] = moments[2,*] - moments[1,*]^2/moments[0,*]	
	3:	moments[2,*,*] = moments[2,*,*] - moments[1,*,*]^2/moments[0,*,*]
	4:	moments[2,*,*,*] = moments[2,*,*,*] - moments[1,*,*,*]^2/moments[0,*,*,*]
ENDCASE
ENDFOR
IF KEYWORD_SET(Avg) THEN BEGIN
CASE ndims OF
	1:	BEGIN
		moments[2] = moments[2]/moments[0]
		moments[1] = moments[1]/moments[0]
		moments[0] = moments[0]/(GridValues[imax]-GridValues[imin])
		END
	2:	BEGIN
		moments[2,*] = moments[2,*]/moments[0,*]
		moments[1,*] = moments[1,*]/moments[0,*]
		moments[0,*] = moments[0,*]/(GridValues[imax]-GridValues[imin])
		END
	3:	BEGIN
		moments[2,*,*] = moments[2,*,*]/moments[0,*,*]
		moments[1,*,*] = moments[1,*,*]/moments[0,*,*]
		moments[0,*,*] = moments[0,*,*]/(GridValues[imax]-GridValues[imin])
		END
	4:	BEGIN
		moments[2,*,*,*] = moments[2,*,*,*]/moments[0,*,*,*]
		moments[1,*,*,*] = moments[1,*,*,*]/moments[0,*,*,*]
		moments[0,*,*,*] = moments[0,*,*,*]/(GridValues[imax]-GridValues[imin])
		END
ENDCASE
ENDIF
moment_ptrs = PTRARR(3)
FOR m=0,2 DO BEGIN
CASE ndims OF
	1:	moment_ptrs[m] = PTR_NEW(moments[m])
	2:	moment_ptrs[m] = PTR_NEW(REFORM(moments[m,*]))
	3:	moment_ptrs[m] = PTR_NEW(REFORM(moments[m,*,*]))
	4:	moment_ptrs[m] = PTR_NEW(REFORM(moments[m,*,*,*]))
ENDCASE
ENDFOR
command_str = "result = {" + GKVsClass[ndims-1] + "}"
ok = EXECUTE(command_str)					; Create Class Structure of appropriate type

aGKVsd = {GKVsd}
nTagsGKVsd = N_TAGS(aGKVsd)
PTR_FREE, aGKVsd.values
PTR_FREE, aGKVsd.ErrorBars
FOR i=0, nTAgsGKVsd-1 DO 		$
	result.(i) = self.(i)					; Populate GKVsd field of result

iResultAxis = 1
iSelfAxis = 1 				
WHILE iResultAxis LE ndims-1 DO BEGIN			; Populate Grid structures of result
	IF(axisNum EQ iResultAxis) THEN	$		; Skip self.Grid'axisNum'
		iSelfAxis=iSelfAxis + 1
	str_Result	= STRING(iResultAxis, FORMAT ='(I1)')	
	str_Self 	= STRING(  iSelfAxis, FORMAT ='(I1)')	
	command_str = "result.Grid" + str_Result + " = self.Grid" + str_Self
	ok = EXECUTE(command_str)
	command_str = "result.Grid" + str_Result + ".values = PTR_NEW(*self.Grid" + str_Self + ".values)"
	ok = EXECUTE(command_str)
	iResultAxis = iResultAxis + 1
	ISelfAxis = iSelfAxis + 1
ENDWHILE
command_str = 'axisTitle = self.Grid' + str_axis + '.title'
ok = EXECUTE(command_str)
min_str = STRING(range[0], FORMAT='(G10.3)')
max_str = STRING(range[1], FORMAT='(G10.3)')
min_str = STRTRIM(min_str, 2)
max_str = STRTRIM(max_str, 2)
Titles 	= STRARR(3)
Titles[0]	= '!MI!S!A!E!8' + max_str + '!R!B!I' + min_str + '!N!X ' + self.title + ' d' + axisTitle
Titles[1] 	= '!MI!S!A!E!8' + max_str + '!R!B!I' + min_str + '!N!X ' + axisTitle + ' ' + self.title + ' d' + axisTitle
Titles[2] 	= '!MI!S!A!E!8' + max_str + '!R!B!I' + min_str + '!N!X (' + axisTitle + '!e2!N - !12<!X' + axisTitle + '!12>!X!e2!N) '	$
			+ self.title + ' d' + axisTitle
IF KEYWORD_SET(Avg) THEN BEGIN
	Titles[0]	= '!12<!X' + self.title + '!12>!X!I' + axisTitle  + '!N'
	Titles[1] 	= '!12<!X' + axistitle  + '!12>!X!I' + self.title + '!N'
	Titles[2] 	= '!12<!X' + axistitle + '!E2!N - !12<!X' + axistitle + '!12>!X!E2!N!12>!X' + Subscript(self.title) + '!N'
ENDIF


command_str = 'gridUnits = self.Grid' + str_axis + '.units'
ok = EXECUTE(command_str)
Units		= STRARR(3)
units[0] = self.units + '*(' + gridUnits + ')'
units[1] = self.units + '*(' + gridUnits + '!e2!N)'
units[2] = self.units + '*(' + gridUnits + '!e3!N)'
IF KEYWORD_SET(Avg) THEN BEGIN
	Units[0]	= self.units
	Units[1] 	= gridUnits
	Units[2] 	= gridUnits + '!E2!N'
ENDIF
indices = self -> IndexRemove(axisNum)
output = OBJARR(3)
result.title = titles[0]
result.indices = PTR_NEW(indices)
result.units = units[0]
result.values = moment_ptrs[0]
result.vrange = [0.,0.]
output[0] = OBJ_NEW(GKVsClass[ndims-1], result)
FOR m=1,2 DO BEGIN
	result = output[0] -> MakeCopy(/NoValues, /NoErrorBars)
	result.title = titles[m]
	result.indices = PTR_NEW(indices)
	result.units = units[m]
	result.values = moment_ptrs[m]
	vmin = GKVsd_MIN(*result.values, MAX=vmax)
	result.vrange = [vmin, vmax]
	output[m] = result
ENDFOR
IF KEYWORD_SET(Avg) THEN BEGIN
	CASE ndims OF
		1:	errors = SQRT( moments[2] )
		2:	errors = SQRT( moments[2,*] )
		3:	errors = SQRT( moments[2,*,*] )
		4:	errors = SQRT( moments[2,*,*,*] )
	ENDCASE
	output[1] -> SET, errorBars = PTR_NEW(errors)
ENDIF
RETURN, output
END ; ***** GKVs1D::Moments ***** ;
	

PRO GKVs1D::SignalWindow, Axis=axisNum, range=rangeIn, irange=irangeIn, _Extra=extra
;
; Generic SignalWindow routine for GKV data objects.
; Sets signal window to values specified
; either by giving axis numer and range 
; (index or values); or by giving axis mnemonic
; and range in the form:	
;
;		mnemonic = [min, max]
;
; where [min,max] is a two-element array.  
; SignalWindow finds the indices to self.Gridx.values
; whose values are closest to these array values.   
;
; Written to be Polymorphic:
;
; The _Extra keyword inheritance feature of IDL is used to
; allow axes to be identified by mnemonics
;
;	Written 2/3/00
;	by W.M. Nevins
;
IF(N_ELEMENTS( rangeIN) GT 0) THEN  range=rangeIn
IF(N_ELEMENTS(irangeIN) GT 0) THEN irange=irangeIn
aGKVsd = {GKVsd}
nTagsGKVsd = N_TAGS(aGKVsd)
PTR_FREE, aGKVsd.values
PTR_FREE, aGKVsd.ErrorBars
GKVsClass = STRARR(5)
GKVsClass = ['GKVSD', 'GKVS1D', 'GKVS2D', 'GKVS3D', 'GKVS4D']
Class = OBJ_CLASS(self)					; Get class of 'self'
Class = STRUPCASE(Class)					; Make sure it's uppercase
CASE Class OF							; Get number of dimensions
	"GKVSD"  :	ndims = 0
	"GKVS1D" :	ndims = 1
	"GKVS2D" :	ndims = 2
	"GKVS3D" :	ndims = 3
	"GKVS4D" :	ndims = 4
ELSE : BEGIN							; self isn't a GKV object???
	messageTxt = 'SignalWindow method undefined for ' + Class + ' Class'
	MESSAGE, messageTxt, /Informational
	RETURN
ENDELSE
ENDCASE
IF ( ndims EQ 0 ) THEN BEGIN
	MESSAGE, 'SignalWindow method undefined for GKVsd Class', /Informational
	RETURN
ENDIF
IF (N_ELEMENTS(axisNum) EQ 0) THEN BEGIN		; Try to get valid axis mnemonic
	numTags = N_TAGS(extra)				;	from _Extra
	IF (numTags EQ 0) THEN BEGIN
		MESSAGE, 'No axis identifier', /INFORMATIONAL
		RETURN
	ENDIF
	tagName = TAG_NAMES(extra)					
	FOR iaxis = 1, ndims DO BEGIN				
		str_axis = STRING(iaxis, FORMAT ='(I1)')
		command_str = 'mnemonic = self.Grid' + str_axis + '.mnemonic'
		ok = EXECUTE(command_str)			; Get mnemonic associated with the iaxis^th axis
		mnemonic = STRUPCASE(mnemonic)		; Force all characters to be upper case
		FOR iTag = 0, numTags-1 DO BEGIN
			tname = STRUPCASE(tagName(itag))
			IF ( tname EQ mnemonic ) THEN BEGIN
				mrange = extra.(itag)		; Get value from corresponding element of 'extra'
				GOTO, GOT_IT
			ENDIF
		ENDFOR
	ENDFOR							; Couldn't match axis... 
	MESSAGE, 'No valid axis mnemonic', /INFORMATIONAL
	RETURN
GOT_IT : axisNum = iaxis
ENDIF
str_axis = STRING(axisNum, FORMAT ='(I1)')	

IF(N_ELEMENTS(mrange) EQ 2) THEN BEGIN
	mtype = TypeOf(mrange)
	CASE mtype OF
		0:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
		1:	range = mrange
		2:	range = mrange
		3:	range = mrange
		4:	range = mrange
		5:	range = mrange
		6:	range = FLOAT(mrange)
		7:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
		8:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
		9:	 range = FLOAT(mrange)
		10:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
		11:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
		ELSE:	BEGIN
				MESSAGE, 'Invalid range specifier', /Informational
				RETURN
			END
	ENDCASE
ENDIF
	
IF(N_ELEMENTS(range) EQ 2) THEN BEGIN
	vmin = range[0]
	vmax = range[1]
	command_str = 'value_ptr = self.Grid' + str_axis + '.values'
	ok = EXECUTE(command_str)
	values = *value_ptr						; Dereference pointer to array of Grid values
	temp = (values - vmin)^2
	vvmin = MIN(temp, index)					; Find index closest to 'vmin=range[0]'
	IF (index LT 0) THEN index = 0
	IF (index GE N_ELEMENTS(values)) THEN index = N_ELEMENTS(values)-1
	irange = LONARR(2)
	irange[0] = index
	temp = (values - vmax)^2
	vvmax = MIN(temp, iindex)					; Find index closest to 'vmax=range[1]'
	IF (iindex LE index) THEN iindex = index+1
	IF (iindex GE N_ELEMENTS(values)) THEN iindex = N_ELEMENTS(values)-1
	irange[1] = iindex
ENDIF
IF(N_ELEMENTS(irange) EQ 2) THEN BEGIN
;
; Reset 'irange' field of grid structure
;
	command_str = 'grid = self.Grid' + str_axis
	ok = EXECUTE(command_str)
	grid.irange = irange
;
; Reset 'uniform' field of grid structure to reflect uniformity
; of grid over current irange
;
	gridValues = *grid.values
	grid.uniform = GKVsd_UniformGrid( gridValues[irange[0]:irange[1]] )
	command_str = 'self.Grid' + str_axis + ' = grid'
	ok = EXECUTE(command_str)
;
; Now, reset vrange to max/min over values within reset irange
;
   IF( PTR_VALID(self.values) ) THEN BEGIN
	values = *self.values
	imin1 = (self.Grid1.irange)[0]
	imax1 = (self.Grid1.irange)[1]
	IF(ndims GT 1) THEN BEGIN
		imin2 = (self.Grid2.irange)[0]
		imax2 = (self.Grid2.irange)[1]
		IF(ndims GT 2) THEN BEGIN
			imin3 = (self.Grid3.irange)[0]
			imax3 = (self.Grid3.irange)[1]
			IF(ndims GT 3) THEN BEGIN
				imin4 = (self.Grid4.irange)[0]
				imax4 = (self.Grid4.irange)[1]
			ENDIF
		ENDIF
	ENDIF
	CASE ndims OF
		1:	values = values[imin1:imax1]
		2:	values = values[imin1:imax1, imin2:imax2]
		3:	values = values[imin1:imax1, imin2:imax2, imin3:imax3]
		4:	values = values[imin1:imax1, imin2:imax2, imin3:imax3, imin4:imax4]
		ELSE:	MESSAGE, 'SignalWindow only supports GKV objects with no more than 4 dimensions', /INFORMATIONAL
	ENDCASE
	vmin = GKVsd_MIN(values, MAX=vmax)
	IF(vmin EQ vmax) THEN vmax = vmin + 1
	self.vrange = [vmin, vmax]
   ENDIF
   RETURN	; THIS is the normal return
ENDIF
MESSAGE, "Couldn't successfully parse keywords", /Informational
RETURN ; ****** GKVs1D::SignalWindow ****** ;
END


FUNCTION GKVs1D::GetValues, All=all, Open=open
;
; Returns pointer to an array containing values which fall within the 
; signal window of 'self'
;
;	Written by W.M. Nevins
;	2/6/00
;
iimax = N_ELEMENTS(*self.Grid1.values) - 1
IF(KEYWORD_SET(open)) THEN BEGIN
	IF(self.grid1.boundary EQ "periodic (closed)") THEN iimax=iimax-1 
ENDIF
imin = self.Grid1.irange[0]
imax = self.Grid1.irange[1] < iimax
IF N_ELEMENTS(all) THEN BEGIN
	imin = 0
	imax = iimax
ENDIF
values = (*self.values)[imin:imax]
RETURN, PTR_NEW(values)
END ; ****** GKVs1d::GetValues ******


FUNCTION GKVs1D::GetErrors, All=all
;
; Returns pointer to an array containing errors which fall within the 
; signal window of 'self'
;
;	Written by W.M. Nevins
;	10/10/00
;
imin = self.Grid1.irange[0]
imax = self.Grid1.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	imin = 0
	imax = N_ELEMENTS(self.Grid1.values) - 1
ENDIF
errors = (*self.ErrorBars)[imin:imax]
RETURN, PTR_NEW(errors)
END ; ****** GKVs1d::GetErrors ******


FUNCTION GKVs1D::SameGrid, argObj, All=all
;
; Returns 1 if 'self' and argObj have same grids, 0 otherwise
; Only need to check Grid1 -- other grids are checked by
; GKVs2D::SameGrid, GKVs3D::SameGrid, ...
;
;	Written 2/6/00
;	by W.M. Nevins
;
selfDims = self -> NumDims()
argDims  = argObj-> NumDims()
IF(argDims NE selfDims) THEN RETURN, 0
imin = self.Grid1.irange[0]
imax = self.Grid1.irange[1]
IF KEYWORD_SET(all) THEN BEGIN
	imin = 0
	imax = N_ELEMENTS(*self.Grid1.values) - 1
ENDIF
iargmin = argObj.Grid1.irange[0]
iargmax = argObj.Grid1.irange[1]
IF KEYWORD_SET(all) THEN BEGIN
	iargmin = 0
	iargmax = N_ELEMENTS(*argObj.Grid1.values) - 1
ENDIF
IF((imax-imin) NE (iargmax-iargmin)) THEN RETURN, 0
sValues = (*self.Grid1.values)[imin:imax]
argValues = (*argObj.Grid1.values)[iargmin:iargmax]
Err = TOTAL((sValues - argValues)^2)/(1+imax-imin)
L = sValues[imax-imin] - sValues[0]
dL = L/(imax-imin)
IF(ERR GT 1.e-3*(dL^2)) THEN RETURN, 0
RETURN, 1
END ; ****** GKVs1D::SameGrid ****** ;


FUNCTION GKVs1D::XCORR, Ref = RefObj, No_Avg=noAverage, Dw=dataWindow, _Extra=extra
;
; Returns a GKVsd object containing the cross correlation function between the
; data in 'self' and the data in RefOjb (which BOTH must be GKVsd objects)
;
; If no Ref_Obj is specified, then returns auto-correlation function of data in 'self'.
; 
; Remaining keywords are passed to XCORR (see definition of XCORR to understand their use)
;
;	Written 2/6/00 by W.M. Nevins
;
;*** NEED TO CHECK FOR UNIFORM GRIDS!!! *****
rflag = 0
No_Avg = 1
IF(KEYWORD_SET(noAverage)) THEN No_Avg = noAverage
ndims = self -> NumDims()						; Get dimensionaliity of data in 'self'
IF((ndims LT 1) OR (ndims GT 4)) THEN BEGIN
	MESSAGE, 'Bad dimensionality, cannot form cross correlations function', /Informational
	RETURN, 0
ENDIF
spacing = FLTARR(ndims+1)
ngrids  = LONARR(ndims+1)
FOR i=1, ndims DO BEGIN						; Check for uniform grids
	dimStr = STRING(i, FORMAT='(I1)')
	command_str = 'Grid = self.Grid' + dimStr
	ok = EXECUTE(command_str)
	imin = Grid.irange[0]
	imax = Grid.irange[1]
	uniform = Grid.uniform					 ; was GKVsd_UniformGrid((*Grid.values)[imin:imax])
	IF(NOT uniform) THEN BEGIN
		MESSAGE, 'WARNING: Grid(s) may not be uniform', /Informational
	ENDIF
	spacing(i)	= (*Grid.values)[imin+1] - (*Grid.values)[imin]
	ngrids(i)	= imax-imin+1
ENDFOR
dt = spacing(ndims)
selfValuePtr = self -> GetValues()	; makes a pointer to a COPY of self.values
IF (N_ELEMENTS(RefObj) NE 0) THEN BEGIN			; Check for Reference object
	rflag = 1
	refType = typeOf(RefObj)
	IF (refType NE 11) THEN BEGIN				; RefObj is NOT an object
		MESSAGE, 'Ref is not an Object', /Informational
		RETURN, 0
	ENDIF
	sameGrid = self -> SameGrid(refObj)
	IF(NOT sameGrid) THEN BEGIN
		MESSAGE, 'Reference object grid is incompatible with self grid', /Informational
		RETURN, 0
	ENDIF
	refValuesPtr = RefObj -> GetValues()	; makes a pointer to a COPY of refObj.values
ENDIF
corrValuesPtr = XCorr(selfValuePtr, Ref=refValuesPtr, /Ptr, No_Avg=No_Avg, DT=dt,Dw=dataWindow, _Extra=extra)
PTR_FREE, selfValuePtr			; Clear pointer to copy of self values
IF(rflag) THEN PTR_FREE, refValuesPtr	; Clear pointer to copy of reference values
corrInfo = SIZE(*corrValuesPtr)
ngrids(ndims) = corrInfo[ndims]

; Make GKVsd object to contain correlation function
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
;
; Start filling tags of result
;
result.mnemonic = 'XCorr_' + self.mnemonic			; Set Mnemonic
IF(rflag) THEN result.mnemonic = result.mnemonic + '_' + RefObj.mnemonic
result.title = 'C{' + self.title + '}' 		; Set Title
IF(rflag) THEN result.title = 'C{' + self.title + ', ' + RefObj.title + '}'
;result.title = result.title + '!N!X'
;result.Indices = PTR_NEW(*self.Indices)
result.units = '(' + result.units + ')!U2!N'

PTR_FREE, result.values						; Set Values pointer
result.values = corrValuesPtr
vmin = GKVsd_MIN(*result.values, MAX=vmax)			; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange
PTR_FREE, result.ErrorBars					; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on correlation functions?
iZero = LONARR(ndims+1)						; Array to hold index of value at zero lag
FOR i=1, ndims DO BEGIN						; Loop over ndims instances of the Grid Class
	dimStr = STRING(i, FORMAT='(I1)')			
	GridStr = 'result.Grid' + dimStr
	command_str = 'Grid = ' + gridStr			; Copy ith instance of the Grid Class
	ok = EXECUTE(command_str)				;	 into Grid
	Mnemonic = Grid.mnemonic				; Set Grid Mnemoinc
	title = Grid.title					; Set Grid Title					
	IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
		Grid.mnemonic = 'tau'
		Grid.title = '!4s!X'
	ENDIF ELSE BEGIN
		Grid.title = '!4D!X' + title
	ENDELSE							; Grid Units can be left as is ...
	PTR_FREE, Grid.values					; Clear old grid.values pointer						
	n = ngrids[i]
	gridvals = spacing(i)*(FINDGEN(n) -n/2)			; Generate array of grid values
	GridPtr = PTR_NEW(gridvals)
	iZero[i] = n/2	
	Grid.values = GridPtr					; Set Grid Values pointer
	Grid.boundary = 'Periodic'				; Correlation function BC's are always Periodic?		
	Grid.uniform = 1					; Correlation functions computed on uniform grid
	min = MIN(gridvals, Max=max)				; Set Grid plot Range
	range = [min, max]
	Grid.range = range
	irange = [0, n-1]					; Set Grid signal window range, irange
	Grid.irange = irange
	command_str = GridStr + ' = Grid'			; Copy ith instance of Grid Class
	ok = EXECUTE(command_str)				;	back into result
ENDFOR
	norm = GetKeyWord('Norm', extra)			; Check for keyword 'norm'
	IF Query_Integer(norm) THEN BEGIN			; If 'norm' is set, then normalize
		IF(norm EQ 1) THEN BEGIN			; correlation function such that the value
			result.units = ''			; Normed correlation functions are unitless.
			values = *result.values			; at zero lag is 1.0.
			PTR_FREE, result.values
			CASE ndims OF
				1:	values = values/values[iZero[1]]
				2:	values = values/values[iZero[1], iZero[2]]
				3:	values = values/values[iZero[1], iZero[2], iZero[3]]
				4:	values = values/values[iZero[1], iZero[2], iZero[3], iZero[4]]
			ENDCASE
			result.values = PTR_NEW(values)
			result.vrange = [-1,1]
		ENDIF
	ENDIF

RETURN, result	
END ; ****** GKVs1D::XCORR ****** ;


FUNCTION GKVs1D::XSPECT, 	Ref = RefObj, No_Avg=noAverage, Lw=lagWindow,	iLw=iLagWindow, $
				Dw=dataWindow, iDw=iDataWindow,  _Extra=extra, force=force
;
; Returns a GKVsd object containing the cross spectrum between the
; data in 'self' and the data in RefOjb (which BOTH must be GKVsd objects)
;
; If no Ref_Obj is specified, then returns spectral density of data in 'self'.
;
;
; 
; Remaining keywords are passed to XCORR (see definition of XCORR to understand their use)
;
;	Written 2/6/00 by W.M. Nevins
;

rflag = 0
No_Avg = 1
IF(KEYWORD_SET(noAverage)) THEN No_Avg = noAverage
ndims = self -> NumDims()						; Get dimensionaliity of data in 'self'
IF((ndims LT 1) OR (ndims GT 4)) THEN BEGIN
	MESSAGE, 'Bad dimensionality, cannot form cross correlations function', /Informational
	RETURN, 0
ENDIF
spacing = FLTARR(ndims+1)
ngrids  = LONARR(ndims+1)
FOR i=1, ndims DO BEGIN						; Check for uniform grids
	dimStr = STRING(i, FORMAT='(I1)')
	command_str = 'Grid = self.Grid' + dimStr
	ok = EXECUTE(command_str)
	imin = Grid.irange[0]
	imax = Grid.irange[1]
	uniform = GKVsd_UniformGrid((*Grid.values)[imin:imax])
	IF(NOT uniform) THEN BEGIN
		MESSAGE, 'WARNING: Grid ' + dimstr + ' is not uniform.  Result may be invalid!', /Informational
		; RETURN, 0
	ENDIF
	spacing(i)	= (*Grid.values)[imin+1] - (*Grid.values)[imin]
	ngrids(i)	= imax-imin+1
ENDFOR
dt = spacing(ndims)
selfValuePtr = self -> GetValues()
IF N_ELEMENTS(RefObj) THEN BEGIN				; Check for Reference object
	rflag = 1
	refType = typeOf(RefObj)
	IF (refType NE 11) THEN BEGIN				; RefObj is NOT an object
		MESSAGE, 'Ref is not an Object', /Informational
		RETURN, 0
	ENDIF
	sameGrid = self -> SameGrid(refObj) 
	IF KEYWORD_SET(force) THEN sameGrid=1
	IF(NOT sameGrid) THEN BEGIN
		MESSAGE, 'Reference object grid is incompatible with self grid', /Informational
		RETURN, 0
	ENDIF
	refValuesPtr = RefObj -> GetValues()	
ENDIF
;
; Set default lag window and data windows
;
iLw = ngrids[nDims]/2 > 1
iDw = ngrids[nDims]/10 < SQRT(ngrids[nDims])
iDw = iDw > 1
IF(N_ELEMENTS( lagWindow) GT 0) THEN BEGIN
	iLw=FIX( lagWindow/dt)
	iDw = FIX(iLw/10)
ENDIF
IF(N_ELEMENTS(dataWindow) GT 0) THEN iDw=FIX(dataWindow/dt)

spectValuesPtr = XSpect(selfValuePtr, Ref=refValuesPtr, /Ptr, No_Avg=No_Avg,  $
				  Dt=dt, Lw=iLw, Dw=iDw, _Extra=extra)
spectInfo = SIZE(*spectValuesPtr)
ngrids[ndims] = spectInfo[ndims]					; Correct lasts element of ngrids for ncorrs > nt
;
; clear pointers
;
PTR_FREE, selfValuePtr
IF(rflag) THEN PTR_FREE, refValuesPtr
;
; Make GKVsd object to contain correlation function
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
;
; Start filling tags of result
;
result.mnemonic = 'XSpect_' + self.mnemonic		; Set Mnemonic
IF(rflag) THEN result.mnemonic = result.mnemonic + '_' + RefObj.mnemonic
result.title = 'S{' + self.title				; Set Title
IF(rflag) THEN result.title = result.title + ', ' + RefObj.title
result.title = result.title + '}'
IF(rflag) THEN BEGIN						; Set units
	result.units = '(' + result.units + ')*(' + refObj.units + ')'
ENDIF ELSE BEGIN						;   for cross spectra
	result.units = '(' + result.units + ')!U2!N'		;   and for auto spectra
ENDELSE 
PTR_FREE, result.ErrorBars					; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on correlation functions?
spectNorm = 1.
FOR i=1, ndims DO BEGIN						; Loop over ndims instances of the Grid Class
	dimStr = STRING(i, FORMAT='(I1)')			
	GridStr = 'result.Grid' + dimStr
	command_str = 'Grid = ' + gridStr			; Copy ith instance of the Grid Class
	ok = EXECUTE(command_str)				;	 into Grid
	result.units = result.units + '*(' + Grid.units + ')'
								; Correct units for current wavenumber/frequency
	Mnemonic = Grid.mnemonic				; Set Grid Mnemoinc
	IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
		Grid.mnemonic = 'omega'
	ENDIF ELSE BEGIN
		Grid.mnemonic = 'k_' + Grid.mnemonic
	ENDELSE
	title = Grid.title					; Set Grid Title					
	IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
		Grid.title = '!4x!X'
	ENDIF ELSE BEGIN
		Grid.title = 'k!I' + title + '!N!X'
	ENDELSE
	Grid.units = '1/(' + Grid.units	+ ')'			; Set Grid Units
	PTR_FREE, Grid.values						
	n = ngrids[i]
	dgrid = 2*!PI/(spacing[i]*nGrids[i])
	gridvals = dgrid*(FINDGEN(n) -n/2)	; Generate array of grid values
	spectNorm = spectNorm/dgrid
	;
	; Patch to compensate for the "Reverse"
	;
	GridPtr = PTR_NEW(gridvals)	
	Grid.values = GridPtr						; Set Grid Values pointer
	Grid.boundary = 'periodic (open)'				; Correlation function BC is always 'periodic (open)'?
	Grid.uniform = 1						; Correlation functions computed on uniform grid
	min = MIN(gridvals, Max=max)					; Set Grid plot Range
	range = [min, max]
	Grid.range = range
	irange = [0, n-1]						; Set Grid signal window range, irange
	Grid.irange = irange
	command_str = GridStr + ' = Grid'				; Copy ith instance of Grid Class
	ok = EXECUTE(command_str)					;	back into result
ENDFOR
PTR_FREE, result.values						; Set Values pointer
result.values = PTR_NEW(spectNorm*(*spectValuesPtr))
PTR_FREE, spectValuesPtr
vmin = GKVsd_MIN(*result.values, MAX=vmax)			; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange


result -> set, _EXTRA=Extra

RETURN, result	

END ; ****** GKVs1D::XSPECT ****** ;


FUNCTION GKVs1D_Gen, Nt=n1, Amplitude=a, omega=k1, 	$
			dt=dt,	dOmega=bandwidth
;
; Set up sample GKVs2D signal object
;
nx=64L
IF ( N_ELEMENTS(n1) NE 0 ) THEN nx=n1
dx = 2*!PI/(nx-1)
IF(N_ELEMENTS(dt) NE 0) THEN dx=dt
k0 = 2.*!PI/(nx*dx)

amplitude=1.
IF KEYWORD_SET(a)  then amplitude=a

kx = 2.
Del_k = 2.
IF N_ELEMENTS(k1) then kx=k1
IF N_ELEMENTS(bandwidth) then del_k=bandwidth
kx =FIX(kx/k0)
del_k = FIX(del_k/k0)

xgrid = dx*indgen(nx)
i=COMPLEX(0.,1.)
argx= -0.5d*((FINDGEN(nx) - nx/2)^2/Del_k^2 > (-25.0d) )	; Put max values in center of array
argx = SHIFT(argx, kx-nx/2)					; Now shift these values to desired location
rarray=RANDOMN(seed, nx)	
iarray=RANDOMN(nseed, nx)
array=COMPLEX(rarray, iarray)*( exp(argx) )
array=FFT(array, 1, /OVERWRITE)
msarray= TOTAL(array*CONJ(array))/( LONG(nx) ) 
values=amplitude*array/sqrt(msarray)
values = CONJ(values)
vmin = GKVsd_MIN(values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1
;errors = 0.-5*amplitude*REPLICATE(1.,nx)

Indices = REPLICATE('*', 1)

signal={GKVs1D}
;
; GKVsd tags
;
	signal.mnemonic		= "sgen1D"
	signal.Title		= "!13Title!N!X"
	signal.Indices		= PTR_NEW(Indices)
	signal.units		= "units"
	PTR_FREE, signal.values
	signal.values		= PTR_NEW(values)
IF (N_ELEMENTS(ErrorBars) NE 0) THEN	BEGIN
	PTR_FREE, signal.ErrorBars
	signal.ErrorBars	= PTR_NEW(errors)
ENDIF
	signal.codeName		= "GKV selftest"
	signal.codePI		= "W.M. Nevins"	
	signal.RunID		= "Null run"
	signal.FileID		= "Self-Generated"
;
; remaining GKVs1D tags
;
	signal.Grid1.mnemonic	= "t"
	signal.Grid1.title 	= "t"
	signal.Grid1.Units 	= ""
	PTR_FREE, signal.Grid1.values
	signal.Grid1.values	= PTR_NEW(xgrid)
	signal.Grid1.Boundary	= "open"
	signal.Grid1.uniform	= 1
	signal.Grid1.range	= [xgrid[0], xgrid[nx-1]]
	signal.Grid1.irange	= [0, nx-1]
;
; Create signal object 
;
result = OBJ_NEW("GKVs1D", signal)
 
RETURN, result

END ; ***** GKVs1D_Gen ***** ;


Pro GKVs1D::Draw, _Extra=extra
;
; 'Virtual compile-time keywords:
;			  xrange=x_range, xstyle=x_style, title=T_title, xtitle=x_title,	$
;			  ytitle=y_title, yrange=y_range, Pretty=pretty, Log=log, YLog=ylog
;
; Default Plotting routine for GKV objects (2-D or more).
; Allows any legal (to PLOT) graphics keyword, but gets defaults from
; object definition if not supplied on command line.
;
; In additon, allows user to set range and indices using axis mnemonics
;
; First, parse 'extra' for keywords which we wish to intercept 
; (this is the moral equivalent of defining them on the PRO line...
;  but prevents IDL from enforcing keyword abreviation, which would
;  interfer with use of mnemonics as 'run-time' keywords)
;
result = GetKeyWord('xrange', extra)				; Equivalent to having
IF(TypeOf(result) NE 7) THEN x_range = result		;	xrange = x_range 
result = GetKeyWord('xtitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN x_title=result	;	xtitle = x_title
result = GetKeyWord('title', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN T_title=result	;	title = T_title
result = GetKeyWord('ytitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN y_title=result	;	ytitle = y_title
result = GetKeyWord('yrange', extra)
IF(TypeOf(result) NE 7) THEN y_range=result		;	yrange = y_range
result = GetKeyWord('vRange', extra)
IF(typeOf(result) NE 7) THEN vrange = result		;	vrange = vrange
result = GetKeyWord('Pretty', extra)
IF(TypeOf(result) NE 7) THEN pretty=result			;	Pretty = pretty
result = GetKeyWord('Log', extra)
IF(TypeOf(result) NE 7) THEN log=result			;	Log = log
result = GetKeyWord('YLog', extra)
IF(TypeOf(result) NE 7) THEN ylog=result			;	YLog = ylog
result = GetKeyWord('RealOnly', extra)
IF(TypeOf(result) NE 7) THEN realOnly=result		;	RealOnly = realOnly
result = GetKeyWord('Thick', extra)				;	Thick = thick
IF(TypeOf(result) NE 7) THEN thick = result
ShowErrorBars = 1
result = GetKeyWord('NoErrorBars', extra)
IF(TypeOf(result) NE 7) THEN ShowErrorBars = 1 - KEYWORD_SET(result)
;										;   NoErrorBars = ...
;
;										; on the command line

imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
x = (*self.Grid1.Values)[imin:imax]
y = (*self.values)[imin:imax]
info=SIZE(y)
index=info[0]+1
ytype=info[index]
;
; Check command line for axis mnemonics
;
axisInfo = self -> GetAxis(extra)
axisValue=axisInfo.(0)
IF(TypeOF(axisValue) NE 7) THEN BEGIN	; Mnemonic for axis is used on command line
	IF(N_ELEMENTS(AxisValue) EQ 2) THEN x_range=AxisValue
ENDIF
;
; Use default plot keyword values from Object definition if they
; are not over ridden on command line.
;
xrange=self.Grid1.range
if KEYWORD_SET(x_range) THEN xrange=x_range
extendPlot = 0
IF( (xrange[0] LT x[0]) OR (xrange[1] GT x[imax-imin]) ) THEN extendPlot = 1

indices = self -> IndexString([0,1], pretty=pretty)
indexStr = '[' + STRJOIN(indices, ', ') + ']'
IF(KEYWORD_SET(pretty)) THEN BEGIN
	title=self.title + indexStr + " (" + self.units + ")"
ENDIF ELSE BEGIN
	title=self.mnemonic + indexStr + " (" + self.units + ")"
ENDELSE
IF KEYWORD_SET(T_title) THEN title=T_title

xstyle=1
x_style = GetKeyWord('xstyle', extra)
IF(TypeOf(x_style) NE 7) THEN xstyle=x_style

x_title = self.Grid1.mnemonic + " (" + self.Grid1.units + ")"
IF N_ELEMENTS(pretty)		THEN x_title = self.Grid1.title + " (" + self.Grid1.units + ")"
if KEYWORD_SET(x_title) 	THEN xtitle=x_title

ytitle = "(" + self.units + ")"
if KEYWORD_SET(y_title) THEN ytitle=y_title
ymin = GKVsd_MIN(y, MAX=ymax)				; Get MIN and MAX of dependent variable

yrange = self.vrange
IF KEYWORD_SET(vrange)  THEN yrange = vrange
IF KEYWORD_SET(y_range) THEN yrange=y_range

IF(KEYWORD_SET(log)) THEN BEGIN
	ylog=log
	yOrig = y
	y = ABS(y)
ENDIF

lineThickness = 2
IF KEYWORD_SET(thick) THEN lineThickness = thick > .1

IF(N_ELEMENTS(extra) NE 0) THEN BEGIN
	IF(TypeOF(extra) EQ 8) THEN nextra = extra
ENDIF
;
; Choose an appropriate character size
; (scaled to size of currently open graphics window)
;
newSizes = GKV_CharSize(Old_Sizes = oldSizes)
;
; Reset character size
;
DEVICE, SET_CHARACTER_SIZE = newSizes
;
; Set line thickness (If you do this with "thick", PS plots won't catch the change ...)
;
!P.THICK = !P.THICK*lineThickness

PLOT,	x, y, xrange=xrange, xstyle=xstyle, 				$
	xtitle=xtitle,  yrange=yrange, 					$
	position=[0.15, 0.15, 0.90, 0.90], YLog=ylog, _Extra=nextra 
;
; Restore line  thickness
;
!P.THICK = !P.THICK/linethickness > 1
	
IF( PTR_VALID(self.ErrorBars) AND ShowErrorBars ) THEN BEGIN
	err=(*self.ErrorBars)[imin:imax]
	
	ERRPLOT, x, y-err, y+err
ENDIF

IF ((ytype EQ 6) OR (ytype EQ 9)) THEN BEGIN	; Plot imaginary part of signal as a dashed line.
	IF(KEYWORD_SET(log)) THEN BEGIN
		!P.THICK = !P.THICK*lineThickness	; Set line thickness
		oplot, x, ABS(FLOAT(yOrig))
		!P.THICK = !P.THICK/linethickness > 1	; Restore line  thickness
	ENDIF
	IF(KEYWORD_SET(realOnly)) THEN GOTO, doneReal
	!P.THICK = !P.THICK*lineThickness	; Set line thickness
	IF(KEYWORD_SET(log)) THEN BEGIN
		oplot, x, ABS(IMAGINARY(yOrig)), Linestyle=2
	ENDIF ELSE BEGIN
		OPLOT, x, IMAGINARY(y), Linestyle=2
	ENDELSE
	!P.THICK = !P.THICK/linethickness > 1	; Restore line  thickness
	IF( PTR_VALID(self.ErrorBars) AND ShowErrorBars ) THEN	$
		ERRPLOT, x, IMAGINARY(y)-err, IMAGINARY(y)+err
ENDIF
;
; Extend plot to include values outside the signalwindow?
;
doneReal:

;IF extendPlot THEN BEGIN
;	allx = *self.grid1.values
;	ally = *self.values
;	OPLOT, allx, ally, thick=lineThickness/2
;	IF ((ytype EQ 6) OR (ytype EQ 9)) THEN BEGIN	; Plot imaginary part of signal as a dashed line.
;		IF(KEYWORD_SET(realOnly)) THEN GOTO, finalWords
;		OPLOT, allx, IMAGINARY(ally), Linestyle=2, thick=linethickness/2
;	ENDIF
;ENDIF

finalwords :
;
; Put title at top of Graphics window
;

result = GetKeyWord("NoTitle", Extra)
IF(Query_Integer(result)) THEN GOTO, NoTitle
charSize = GKV_TitleSize(title, Width = (0.90-0.15), yPosition=ywrite)
xwrite = (0.15 + 0.90)/2.*!D.X_SIZE
XYOUTS, xwrite, ywrite, title, ALIGNMENT=0.5, CHARSIZE=charSize, /DEVICE 
NoTitle:

;
; Now write code and run info in lower corners
;
result = GetKeyWord("NoFooter", Extra)
IF(Query_Integer(result)) THEN GOTO, NoFooter 

maxChars =ROUND(0.4*!D.X_SIZE/!D.X_CH_SIZE)		; Compute maximum allowed characters in a line
xwrite=!D.x_ch_size					; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower left-hand corner of plot window.
codeName = STRMID(self.CodeName, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, CodeName, /Device		; Write CodeName to lower left-hand corner.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
codePI = STRMID(self.CodePI, 0, maxChars)		; Truncate 'CodePI' if necessary
XYOUTS, xwrite, ywrite, CodePI, /Device			; Write CodePI below CodeName.
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window.
runID = STRMID(self.RunID, 0, maxChars)			; Truncate 'runID' if necessary
XYOUTS, xwrite, ywrite, RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
;fileID = STRMID(self.FileID, 0, maxChars)		; Truncate 'CodeName' if necessary
fileId = self.FileID
XYOUTS, xwrite, ywrite, FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID
;
; Return character size to input values
;
NoFooter: DEVICE, SET_CHARACTER_SIZE = oldSizes

END ; ***** GKVs1d::Draw ***** ;

PRO GKVs1D::oPlot, RealOnly=realOnly, NoErrorBars=noErrorBars, _Extra=extra

imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
x = (*self.Grid1.Values)[imin:imax]
y = (*self.values)[imin:imax]
info=SIZE(y)
index=info[0]+1
ytype=info[index]
;
; Check for 'color' in _Extra to pass to ERRPLOT
;
localExtra=-1
IF(TypeOf(extra) NE 0) THEN localExtra=extra
result = GetKeyWord('color', localExtra)
IF Query_Integer(result) THEN BEGIN
	color=result
ENDIF ELSE BEGIN
	color=1
ENDELSE

oplot, x, y, _Extra=extra

IF( PTR_VALID(self.ErrorBars) AND (TypeOf(noErrorBars) EQ 0) ) THEN BEGIN
	err=(*self.ErrorBars)[imin:imax]
	PCOLOR=!P.COLOR
	!P.COLOR=color
	ERRPLOT, x, y-err, y+err
	!P.COLOR=PCOLOR
ENDIF
IF( KEYWORD_SET(realOnly) ) THEN RETURN
IF ((ytype EQ 6) OR (ytype EQ 9)) THEN BEGIN	; Plot imaginary part of signal as a dashed line.
	OPLOT, x, IMAGINARY(y), Linestyle=2, _Extra=extra
	IF( PTR_VALID(self.ErrorBars) AND (TypeOf(noErrorBars) EQ 0) ) THEN	BEGIN
		PCOLOR=!P.COLOR
		!P.COLOR=color
		ERRPLOT, x, IMAGINARY(y)-err, IMAGINARY(y)+err
		!P.COLOR=PCOLOR
	ENDIF
ENDIF
RETURN
END ; ***** GKVs1D::oPlot ***** ;


FUNCTION GKVs1D::GA_DATA
;
; Returns a GA_DATA object which points to the **SAME*** data as "self"
;	(so any data manipulation done on the GA_DATA object will also 
;	 change the underlying GKVs1D object...)
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	result = OBJ_NEW('GA_DATA', self.Grid1.values, self.values, 			$
				XName = self.Grid1.title + ' (' + self.Grid1.units + ')',	$
				YName = self.title + ' (' + self.units + ')',			$
				ErrorBar=*self.ErrorBars)
ENDIF ELSE BEGIN
	result = OBJ_NEW('GA_DATA', self.Grid1.values, self.values, 			$
				XName = self.Grid1.title + ' (' + self.Grid1.units + ')',	$
				YName = self.title + ' (' + self.units + ')' )
ENDELSE

RETURN, result

END ; ***** GKVs1D::GA_DATA ***** ;


PRO GKVs1D::GET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices,		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Get values of elements of the realization of the GRID class, Grid1;
; and then call GKVsd::GET to get values of elements of the GKVsd Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs1D,
; or any or its subclasses
;
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 1) THEN 	$
			GKVsd_GetGrid, self.Grid1,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	2:	IF(axisID EQ 1) THEN 	$
			GKVsd_GetGrid, self.Grid1,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	3:	IF(axisID EQ 1) THEN 	$
			GKVsd_GetGrid, self.Grid1,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid1.mnemonic) THEN 	$
			GKVsd_GetGrid, self.Grid1,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self -> GKVsd::GET, 	mnemonic=mnemonic, Title=Title, Indices=indices,$		
				units=units, values=values, vrange=vrange,		$
				ErrorBars=ErrorBars, CodeName=CodeName,		$
				CodePI=CodePI, RunID=RunID, FileID=FIleID	
RETURN
END ; ****** GKVs1D::GET ****** ;


PRO GKVs1D::SET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices,		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
;
; Get values of elements of the realization of the GRID class, Grid1;
; and then call GKVsd::GET to get values of elements of the GKVsd Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs1D,
; or any or its subclasses
;
arg = self.Grid1
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 1) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	2:	IF(axisID EQ 1) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	3:	IF(axisID EQ 1) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid1.mnemonic) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self.Grid1 = arg
self -> GKVsd::SET, 	mnemonic=mnemonic, Title=Title, Indices=indices,$
				units=units, values=values, vrange=vrange,		$
				ErrorBars=ErrorBars, CodeName=CodeName,		$
				CodePI=CodePI, RunID=RunID, FileID=FIleID	
RETURN
END ; ****** GKVs1D::SET ****** ;


PRO GKVs1D::Info
;
; Prints information about contents of GKVs1D objects
;
self -> GKVsd::Info
PRINT, 'Grid1:'
GKVsd_PrintGrid, self.Grid1
RETURN
END ; ***** GKVs1D::Info ***** ;


PRO GKVs1D::CleanUp

Ptr_Free, self.Grid1.Values
self -> GKVsd::CleanUp

RETURN

END ; ***** GKVs1D::CleanUp ***** ;	


FUNCTION GKVs1D::INIT, signal
;
; Preliminary init for testing object methods ...
;
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF 
;
; Set tags for GKVsd class
;
ok = self-> GKVsd::INIT(signal)
;
; Set Tags for Grid1
;
NumTags=N_TAGS(self.Grid1)
FOR itag=0, NumTags-1 DO $
	self.Grid1.(itag) = signal.Grid1.(itag)
;
; remove embedded blanks in grid title, mnemonic, and units
;
self.grid1.title    = STRCOMPRESS(self.grid1.title,    /REMOVE_ALL)
self.grid1.mnemonic = STRCOMPRESS(self.grid1.mnemonic, /REMOVE_ALL)
self.grid1.units    = STRCOMPRESS(self.grid1.units,    /REMOVE_ALL)


IF(PTR_VALID(self.grid1.values)) THEN BEGIN
	;
	; Check for uniform grid
	;
	self.Grid1.uniform = GKVsd_UniformGrid(self.Grid1.values)
	;
	; Check if Grid1.range is set ... and set if necessary
	;
	IF((self.Grid1.range[0] EQ 0.) AND (self.Grid1.range[1] EQ 0.)) THEN $
		self.Grid1.range = [MIN(*self.Grid1.Values, Max=max), max]
	;
	; Check if Grid1.irange is set ... and set if necessary
	;
	IF((self.Grid1.irange[0] EQ 0) AND (self.Grid1.irange[1] EQ 0)) THEN $
		self.Grid1.irange = [0, N_ELEMENTS(*self.Grid1.values)-1]
ENDIF
;
; Return on successful completion
;
RETURN, ok
END ; ***** GKVs1D::INIT ***** ;


PRO GKVs1D__Define
struct = {	GKVs1D,				$	; "GK Visualization signal (1-D)"
		INHERITS GKVsd,			$	; GKVs1D is a subclass of GKVsd
		Grid1:{Grid}			}	; Include one 'Grid' class
								; 
END ; ***** GKVs1D__Define ***** ;
