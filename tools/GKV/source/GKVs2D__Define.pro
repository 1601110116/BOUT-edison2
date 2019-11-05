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
; Defining routine for 2-dimensional signal objects
;
; Written by W.M. Nevins
;	1/17/00
;
FUNCTION GKVs2D::AxisNumber, stringIn, Debug=d
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
result = self -> GKVs1D::AxisNumber(stringIn, debug=debug)
axis2Mnemonic = STRTRIM(self.Grid2.mnemonic, 2)
IF( STRCMP(stringIn, axis2Mnemonic, /FOLD_CASE) ) THEN result = 2
RETURN, result
END ; ****** GKVs2D::AxisNumber ****** ;


FUNCTION GKVs2D::GetAxis, structure, Debug=d
;
; Searches 'structure' for tags which are the same as the mnemonics of axis1 and axis2
; Returns a structure with tags 'axis1' and 'axis2'.  
; The associated values are whatever values were associated with the corresponding mnemonic.
; 
; Written by W.M. Nevins
;	2/20/00
;
debug=0
otherTags = 0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
;
; Check for axis1 mnemonic by calling GKVs1D::GetAxis
;
result1 = self -> GKVs1D::GetAxis(structure, Debug=debug)
;
; Structure now contains input structure with tag to mnemonic of axis1 (if found) deleted
;
result = {axis1:result1.axis1, axis2:'no match'}
IF(N_ELEMENTS(structure) EQ 0) THEN RETURN, result			; Structure is undefined
IF(TypeOf(structure) EQ 2) THEN RETURN, result				; No tags were left to return
nTags = N_TAGS(structure)
IF(nTags EQ 0) THEN RETURN, result
tagNames = TAG_NAMES(structure)
tagNames = STRTRIM(tagNames, 2)							; Remove both leading and trailing blanks
command_str = 'structure = {'
axisMnemonic = STRTRIM(self.Grid2.mnemonic)				; Remove both leading and trailing blanks
FOR i=0, ntags-1 DO BEGIN								; Search tags of 'structure'
	IF( STRCMP(axisMnemonic, tagNames[i], /FOLD_CASE) ) THEN BEGIN
		result = {axis1:result1.axis1, axis2:structure.(i) }	; If multiple occurances of axis3 mnemonic in 'structure'
	ENDIF ELSE BEGIN								;	only the LAST occurance is signifacant
		i_str = STRING(i, FORMAT='(I3)')				; Construct command string for 'Nstructure'
		i_str = STRTRIM(i_str, 2)						; Remove both leading and trailing blanks
		IF (otherTags NE 0) THEN command_str = command_str + ', '
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
END ; ****** GKVs2D::GetAxis ****** ;

FUNCTION GKVs2D::AxisIndex, axisNumber, value, Debug=d
;
; We should only get here if axisNumber is � 2.  
; Function returns index of axisNumber^th  Grid closest to 'value'
;
debug=0
IF(N_ELEMENTS(d) NE 0) then debug=d
IF(axisNumber GT 2) THEN BEGIN
	MESSAGE, 'axisNumber � 2', INFORMATIONAL=d
	RETURN, 0
ENDIF
axisStr = STRING(axisNumber,  FORMAT='(I1)')
commandStr = 'temp = (*self.Grid' + axisStr + '.values - value)^2'
ok = EXECUTE(commandStr)
eps=MIN(temp,index)
RETURN, index
END ; ****** GKVs2D::AxisIndex ****** ;


FUNCTION GKVs2D::GetValues, All=all, open=open
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

jjmax = N_ELEMENTS(*self.Grid2.values) - 1
IF(KEYWORD_SET(open)) THEN BEGIN
	IF(self.grid2.boundary EQ "periodic (closed)") THEN jjmax=jjmax-1 
ENDIF
jmin = self.Grid2.irange[0]
jmax = self.Grid2.irange[1] < jjmax
IF N_ELEMENTS(all) THEN BEGIN
	jmin = 0
	jmax = jjmax
ENDIF
values = (*self.values)[imin:imax, jmin:jmax]
RETURN, PTR_NEW(values)
END ; ****** GKVs2d::GetValues ******


FUNCTION GKVs2D::GetErrors, All=all, open=open
;
; Returns pointer to an array containing ErrorBars which fall within the 
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

jjmax = N_ELEMENTS(*self.Grid2.values) - 1
IF(KEYWORD_SET(open)) THEN BEGIN
	IF(self.grid2.boundary EQ "periodic (closed)") THEN jjmax=jjmax-1 
ENDIF
jmin = self.Grid2.irange[0]
jmax = self.Grid2.irange[1] < jjmax
IF N_ELEMENTS(all) THEN BEGIN
	jmin = 0
	jmax = jjmax
ENDIF
errors = (*self.ErrorBars)[imin:imax, jmin:jmax]
RETURN, PTR_NEW(errors)
END ; ****** GKVs2d::GetErrors ******

FUNCTION GKVs2D::SameGrid, argObj, All=all
;
; Returns 1 if 'self' and argObj have same grid, 0 otherwise
;
;	Written 2/6/00
;	by W.M. Nevins
;
selfDims = self  -> NumDims()
argDims  = argObj-> NumDims()
IF(argDims NE selfDims) THEN RETURN, 0
imin = self.Grid2.irange[0]
imax = self.Grid2.irange[1]
IF  KEYWORD_SET(all) THEN BEGIN
	imin = 0
	imax = N_ELEMENTS(*self.Grid2.values) - 1
ENDIF
sGrid = (*self.Grid2.values)[imin:imax]
iargmin = argObj.Grid2.irange[0]
iargmax = argObj.Grid2.irange[1]
IF KEYWORD_SET(all) THEN BEGIN
	iargmin = 0
	iargmax = N_ELEMENTS(*argObj.Grid2.values) - 1
ENDIF
IF((imax-imin) NE (iargmax-iargmin)) THEN RETURN, 0
argGrid = (*argObj.Grid2.values)[iargmin:iargmax]
Err = TOTAL((sGrid - argGrid)^2)/(1+imax-imin)
L1 = sGrid[imax-imin] - sGrid[0]
d1 = L1/(1+imax-imin)
IF(ERR LT 1.e-3*(d1^2)) THEN RETURN, self -> GKVs1D::SameGrid(argObj)
RETURN, 0
END ; ****** GKVs2D::SameGrid ****** ;


FUNCTION GKVs2D::Interpolate, arg
;
; Interpolate values from self onto the signal window of arg's grid.
;
IF (OBJ_ISA(arg, 'GKVs2D') NE 1) THEN BEGIN
	MESSAGE, "Argument is not a valid GKVs2D object", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Check for common units, independent variable, interval
;
;IF ((self.Grid1.units NE arg.Grid1.units) OR (self.Grid2.units NE arg.Grid2.units)) THEN BEGIN
;	MESSAGE, "Incompatible units (independent variables)"
;	RETURN, 0
;ENDIF
;IF ((self.Grid1.mnemonic NE arg.Grid1.mnemonic) OR (self.Grid2.mnemonic NE arg.Grid2.mnemonic)) THEN BEGIN
;	MESSAGE, "Interpolate:  Incompatible independent variables"
;	RETURN, 0
;ENDIF
IF(self -> SameGRID(arg)) THEN BEGIN
	MESSAGE, "Objects have common grid, no interpolation necessary", /Informational
	result = self -> MakeCopy()
	result -> Restrict
	RETURN, result
ENDIF

x1 = *self.Grid1.values					; *1 -- grid values of self (Obj to be interpolated)
x2 =  *arg.Grid1.values					; *2 -- grid values of arg  (target grid)
y1 = *self.Grid2.values
y2 =  *arg.Grid2.values
values1 = *self.values					; Get values of to be interpolated
imin = arg.Grid1.irange[0]					; Get range of indices on target grid
imax = arg.Grid1.irange[1]
jmin = arg.Grid2.irange[0]
jmax = arg.Grid2.irange[1]
x = x2[imin:imax]						; axis1 for target grid
y = y2[jmin:jmax]						; axis2 for target grid
;
; Get index into old grid arrays (that is, those of 'self') for each element of the new grid arrays (that is, those of 'arg').
;
jjx = VALUE_LOCATE(x1,x)
jjy = VALUE_LOCATE(y1,y)
;
; jjx, jjy will be = -1 if new grid point does not lie within range of old grid points.
; First create jx, jy arrays (with 
jx = jjx > 0
jy = jjy > 0
;
; Compute old grid spacing
;
info = SIZE(values1)
nx1 = info[1] - 1
ny1 = info[2] - 1
dx1 = x1[1:nx1] - x1[0:(nx1-1)]
dy1 = y1[1:ny1] - y1[0:(ny1-1)]
;
; Compute fractional "address" of new x and y grid points within old grid
;
sx = jx + ((x - x1[jx])/dx1[jx])*(jjx GT 0)		; the factor (jjx GT 0) has the effect of zeroing out
sy = jy + ((y  -y1[jy])/dy1[jy])*(jjy GT 0)		;  this term for points outside of target grid.

values = INTERPOLATE(values1, sx, sy, /GRID, cubic=-0.5)	; Interpolate input values to target grid

result = self -> MakeCopy(/NoValues, /NoErrorBars)
PTR_FREE, result.values					; Load pointer to interpolated values
result.values = PTR_NEW(values)				;	into 'result'
vmin = GKVsd_MIN(values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1
result.vrange = [vmin, vmax]

PTR_FREE, result.Grid1.values			
result.Grid1.values = PTR_NEW(x)			; Load Grid1			
result.Grid1.uniform = GKVsd_UniformGrid(x)	; 	...
result.Grid1.range=arg.Grid1.range			
result.Grid1.irange = [0, imax-imin]

PTR_FREE, result.Grid2.values				; Load Grid2
result.Grid2.values = PTR_NEW(y)			;	...
result.Grid2.uniform = GKVsd_UniformGrid(y)
result.Grid2.range = arg.Grid2.range
result.Grid2.irange = [0, jmax-jmin]


IF PTR_VALID(self.ErrorBars) THEN BEGIN			; Interpolate error bars...
	error1 = *self.ErrorBars			;	(probably we could think
	errors = INTERPOLATE(error1,sx,sy, /GRID)	;	 of something a bit more
	PTR_FREE, result.ErrorBars			;	 rigorous?)
	result.ErrorBars = PTR_NEW(errors)
ENDIF
RETURN, result
END ; ***** GKVs2D::Interpolate ***** ;


FUNCTION GKVs2D::MakeCopy, _Extra=extra
;
; Make "deep" copy of self
;
result = self -> GKVs1D::MakeCopy( _Extra=extra)	; Creates new GKVsxD objects, and
										; 	copies elements of GKVs1D class
										
result.Grid2 = GKVsd_GridCopy(self.Grid2)			; Make 'deep' copy of the Grid-class structure Grid2,			
RETURN, result
END ; ***** GKVs2D::MakeCopy ***** ;
		

FUNCTION GKVs2D::Plus, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  addition
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
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						nVals = N_ELEMENTS(*self.Grid2.values)
						argValues = gridValues#Replicate(1,nVals)
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						nVals = N_ELEMENTS(*self.Grid1.values)
						argValues = Replicate(1,nVals)#gridValues
						argTitle = self.Grid2.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORAMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues)
		selfValues = *self.values
		resultValues = selfValues + argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "+" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs2D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)"
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs2D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables"
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax]
			values = argValues + selfValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + " + " + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Plus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs2D::Plus ***** ;


FUNCTION GKVs2D::Minus, argg, title=title, units=units, mnemonic=mnemonic
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
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						nVals = N_ELEMENTS(*self.Grid2.values)
						argValues = gridValues#Replicate(1,nVals)
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						nVals = N_ELEMENTS(*self.Grid1.values)
						argValues = Replicate(1,nVals)#gridValues
						argTitle = self.Grid2.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORAMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues)
		selfValues = *self.values
		resultValues = selfValues - argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "-" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs2D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "WARNING:  Incompatible units (dependent variable)", /INFORMATIONAL
				;RETURN, 0
			ENDIF
		;
		; Difference of data from two GKVs2D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables"
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax]
			values = selfValues - argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + " - " + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Minus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs2D::Minus ***** ;


FUNCTION GKVs2D::Times, argg, title=title, units=units, mnemonic=mnemonic
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
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						nVals = N_ELEMENTS(*self.Grid2.values)
						argValues = gridValues#Replicate(1,nVals)
						argTitle = self.Grid1.title
						argUnits = self.Grid1.units
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						nVals = N_ELEMENTS(*self.Grid1.values)
						argValues = Replicate(1,nVals)#gridValues
						argTitle = self.Grid2.title
						argUnits = self.Grid2.units
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORAMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues, /NoErrorBars)
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
		result.Title = '(' + self.Title + "*" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		result.units = '(' + self.units + ')*(' + argUnits + ')'
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs2D') EQ 1) THEN BEGIN
		;
		; Multiply data from two GKVs2D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables"
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax]
			values = selfValues*argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + " * " + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Times(arg)

ENDCASE

RETURN, result

END ; ***** GKVs2D::Times ***** ;


FUNCTION GKVs2D::Over, argg, title=title, units=units, mnemonic=mnemonic
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
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						nVals = N_ELEMENTS(*self.Grid2.values)
						argValues = gridValues#Replicate(1,nVals)
						argTitle = self.Grid1.title
						argUnits = self.Grid1.units
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						nVals = N_ELEMENTS(*self.Grid1.values)
						argValues = Replicate(1,nVals)#gridValues
						argTitle = self.Grid2.title
						argUnits = self.Grid2.units
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORAMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues, /NoErrorBars)
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
		result.Title = '(' + self.Title + "/" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		result.units = '(' + self.units + ')/(' + argUnits + ')'
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again						;	and try again.
	END
11:	BEGIN								; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs2D') EQ 1) THEN BEGIN
		;
		; Divide data from two GKVs2D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables"
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax]
			values = selfValues/argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = self.Title + " / " + arg.Title
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Over(arg)

ENDCASE

RETURN, result

END ; ***** GKVs2D::Over ***** ;


FUNCTION GKVs2D::GA_DATA
;
; Returns a GA_DATA3 object which points to the **SAME*** data as "self"
;	(so any data manipulation done on the GA_DATA object will also 
;	 change the underlying GKVs1D object...)
;
IF PTR_VALID(self.ErrorBars) THEN BEGIN
	result = OBJ_NEW('GA_DATA3', self.values, self.Grid1.values, self.Grid2.values,	$
				XName = self.Grid1.title + ' (' + self.Grid1.units + ')',		$
				YName = self.Grid2.title + ' (' + self.Grid2.units + ')',		$
				ErrorBar=*self.ErrorBars)
ENDIF ELSE BEGIN
	result = OBJ_NEW('GA_DATA3', self.values, self.Grid1.values, self.Grid2.values,	$
				XName = self.Grid1.title + ' (' + self.Grid1.units + ')',		$
				YName = self.Grid2.title + ' (' + self.Grid2.units + ')')		
ENDELSE

RETURN, result

END ; ***** GKVs2D::GA_DATA ***** ;
	

FUNCTION GKVs2D_Gen, Nx=n1, Ny=n2, Amplitude=a, kx=k1, ky=k2, 	$
				Del_k=bandwidth
;
; Set up sample GKVs2D signal object
;
nx=64L
ny=64L
IF ( N_ELEMENTS(n1) NE 0 ) THEN nx=n1
IF ( N_ELEMENTS(n2) NE 0 ) THEN ny=n2
dx = 2*!PI/(nx-1)
dy = 2*!PI/(ny-1)

amplitude=1.
IF KEYWORD_SET(a)  then amplitude=a

kx = 2.
ky = 2.
Del_k = 2.
IF N_ELEMENTS(k1) then kx=k1
IF N_ELEMENTS(k2) then ky=k2
IF N_ELEMENTS(bandwidth) then del_k=bandwidth

xgrid = dx*indgen(nx)
ygrid = dy*indgen(ny)
i=COMPLEX(0.,1.)
argx= -0.5d*((INDGEN(nx) - kx)^2/Del_k > (-30.0d) )
argy= -0.5d*((INDGEN(ny) - ky)^2/Del_k > (-30.0d))
rarray=RANDOMN(seed, nx, ny)
iarray=RANDOMN(seed, nx, ny)
array=COMPLEX(rarray, iarray)*( exp(argx)#exp(argy) )
array=FFT(array, 1, /OVERWRITE)
msarray= TOTAL(array*CONJ(array))/( LONG(nx)*LONG(ny) ) 
values=amplitude*array/sqrt(msarray)
vmin = GKVsd_MIN(values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1

Indices = REPLICATE('*', 2)

signal={GKVs2D}
;
; GKVsd tags
;
	signal.mnemonic	= "sgen2D"
	signal.Title	= "!13Title!N!X"
	signal.Indices	= PTR_NEW(Indices)
	signal.units	= "units"
	signal.values	= ptr_new(values)
	signal.vrange	= [vmin, vmax]
	signal.codeName	= "GKV selftest"
	signal.codePI	= "W.M. Nevins"	
	signal.RunID	= "Null run"
	signal.FileID	= "Self-Generated"
;
; Grid1 tags
;
	signal.Grid1.mnemonic	= 'x'
	signal.Grid1.title 	= "x-title"
	signal.Grid1.units 	= "x-units"
	signal.Grid1.values  = ptr_new(xgrid)
	signal.Grid1.boundary	= "open"
	signal.Grid1.range 	= [xgrid[0], xgrid[nx-1]]
	signal.Grid1.irange 	= [0, nx-1]
;
; Grid2 tags
;
	signal.Grid2.mnemonic	= 'y'
	signal.Grid2.title 	= "y-title"
	signal.Grid2.units 	= "y-units"
	signal.Grid2.values 	= PTR_NEW(ygrid)
	signal.Grid2.boundary	= "open"
	signal.Grid2.range	= [ygrid[0], ygrid[ny-1]]
	signal.Grid2.irange	= [0, ny-1]
;
; Create signal object
;
	signal_obj=OBJ_NEW("GKVs2D", signal) 
RETURN, signal_obj
END

FUNCTION GKVs2D::XCORR, Ref = RefObj, All=all, No_Avg=noAverage, Dw=dataWindow, _Extra=extra
;
; Returns a GKVsd object containing the cross correlation function between the
; data in 'self' and the data in RefOjb (which BOTH must be GKVsd objects)
;
; If no RefObj is specified, then returns auto-correlation function of data in 'self'.
; 
; Remaining keywords are passed to XCORR (see definition of XCORR to understand their use)
;
;	Written 2/6/00 by W.M. Nevins
;
;*** NEED TO CHECK FOR UNIFORM GRIDS!!! *****
rflag = 0
No_Avg = 1
IF(N_ELEMENTS(noAverage)) THEN No_Avg = noAverage
ndims = self -> NumDims()						; Get dimensionaliity of data in 'self'
IF((ndims LT 1) OR (ndims GT 4)) THEN BEGIN
	MESSAGE, 'Bad dimensionality, cannot form cross correlations function', /Informational
	RETURN, 0
ENDIF
IF(N_ELEMENTS(RefObj) EQ 0) THEN BEGIN			; No reference object... so pass on to GKVs1D
	result = self -> GKVs1D::XCORR(No_Avg=noAverage, Dw=dataWindow, _Extra=extra)
	RETURN, result
ENDIF
refType = TypeOf(RefObj)
IF (refType NE 11) THEN BEGIN				; RefObj is NOT an object
	MESSAGE, 'Ref is not an Object', /Informational
	RETURN, 0
ENDIF
refdims = RefObj -> NumDims()					; Get dimensionality of data in 'RefObj'
IF(refdims EQ ndims) THEN BEGIN					; Dimensionality is same... so pass on to GKVs1D
	result = self -> GKVs1D::XCORR(Ref = RefObj, No_Avg=noAverage, Dw=dataWindow, _Extra=extra)
	RETURN, result
ENDIF
IF(refdims NE (ndims-1) ) THEN BEGIN
	MESSAGE, "Expect dimensionality of Ref = dimensionality 'self', or self - 1" , /informational
	RETURN, 0
ENDIF
spacing = FLTARR(ndims+1)
ngrids  = LONARR(ndims+1)
iref = 1
iaxis = 0
FOR iself=1, ndims DO BEGIN						; 
	iselfStr = STRING(iself, FORMAT='(I1)')
	command_str = 'selfGrid = self.Grid' + iselfStr		; Get iself^th Grid structure of 'self'
	ok = EXECUTE(command_str)
	irefstr = STRING(iref, FORMAT='(I1)')
	command_str = 'refGrid = RefObj.Grid' + irefstr
	ok = EXECUTE(command_str)						; Get iref^th Grid structure of 'RefObj'
	ismin = selfGrid.irange[0]
	ismax = selfGrid.irange[1]
	IF(GKV_GridSame(selfGrid,RefGrid,All=all, /force)) THEN BEGIN	; Check to see if the Grid structures match (force match if possible)
		IF(KEYWORD_SET(All)) THEN BEGIN
			ismin = 0
			ismax = N_ELEMENTS(*selfGrid.values) - 1
		ENDIF
		uniform = selfGrid.uniform							; was GKVsd_UniformGrid((*selfGrid.values)[ismin:ismax])
		IF(NOT uniform) THEN BEGIN					; Check for uniform grids
			MESSAGE, 'WARNING:  Grid ' + iSelfStr + ' is nonuniform.  Result may be invalid!', /Informational
			;RETURN, 0
		ENDIF
		spacing(iself)	= (*selfGrid.values)[ismin+1] - (*selfGrid.values)[ismin]
		ngrids(iself)	= ismax-ismin+1
		iref = iref + 1
	ENDIF ELSE BEGIN							; This selfGrid does not match corresponding refGrid
		IF(iaxis NE 0) THEN BEGIN					;	this should only happen once!
			MESSAGE, 'Bad match between reference grids and self grids', /Informational
			RETURN, 0
		ENDIF
		iaxis= iself							; Prepare for loop which calls XCORR
		imin = 0
		imax = ismax - ismin
	ENDELSE
ENDFOR
IF(iaxis EQ 0) THEN BEGIN
	MESSAGE, "Couldn't find cross correlation axis", /Informational
	RETURN, 0
ENDIF

dt = spacing(ndims)
selfValuePtr =  self -> GetValues() 						; Returns pointer to values within signal window
refValuePtr = RefObj -> GetValues() 
refValues = *refValuePtr
refInfo = SIZE(refValues)
result = GetKeyWord('Norm', extra)
localNorm=0
IF(Query_Integer(result)) THEN localNorm=result
IF(localNorm EQ 1) THEN BEGIN
	nRefDims = refInfo[0]
	refType = refInfo[nRefdims+1]
	nRefElements = refInfo[nRefDims+2]
	IF(Query_Complex(refValues)) THEN BEGIN
		refNormSq = TOTAL(refValues*CONJ(refValues))/nRefElements
	ENDIF ELSE BEGIN 
		refNormSq = TOTAL(refValues*refValues)/nRefElements
	ENDELSE
ENDIF

FOR i = imin, imax DO BEGIN
CASE ndims OF
	2:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *], refInfo[1])
			2:	selfValues = REFORM((*selfValuePtr)[*, i], refInfo[1])
		ENDCASE
	3:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *], refInfo[1], refInfo[2])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *], refInfo[1], refinfo[2])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i], refInfo[1], refInfo[2])
		ENDCASE
	4:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *, *], refInfo[1], refInfo[2], refInfo[3])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *, *], refInfo[1], refInfo[2], refInfo[3])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i, *], refInfo[1], refInfo[2], refInfo[3])
			4:	selfValues = REFORM((*selfValuePtr)[*, *, *, i], refInfo[1], refInfo[2], refInfo[3])
		ENDCASE
ENDCASE
IF(TypeOf(extra) EQ 8) THEN BEGIN
	CorrValues_i = XCorr(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, Dw=dataWindow, _Extra=extra)
ENDIF ELSE BEGIN
	CorrValues_i = XCorr(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, Dw=dataWindow)
ENDELSE

IF(i EQ imin) THEN BEGIN							; Initialization duties
	selfInfo = SIZE(*selfValuePtr)
	nSelfValueElements = N_ELEMENTS(selfValues)
	corrInfo = SIZE(CorrValues_i)
	corrDims = corrInfo[0]
	ncorrs = corrinfo[corrDims]
	corrType = 4
	IF( (Query_Complex(selfValues) + Query_Complex(RefValues)) GT 0 ) THEN corrType=6
	CASE ndims OF
		2:	CorrValues = MAKE_ARRAY(selfInfo[1], ncorrs, TYPE=corrType)
		3:	CorrValues = MAKE_ARRAY(selfinfo[1], selfInfo[2], ncorrs, TYPE=corrType)
		4:	CorrValues = MAKE_ARRAY(selfinfo[1], selfInfo[2], selfInfo[3], ncorrs, TYPE=corrType)
	ENDCASE
ENDIF
;
; Apply local norm if required
;
IF(localNorm EQ 1) THEN BEGIN
	IF(Query_Complex(selfValues)) THEN BEGIN
		selfNormSq_i = TOTAL(selfValues*CONJ(selfValues))/nSelfValueElements
	ENDIF ELSE BEGIN 
		selfNormSq_i = TOTAL(selfValues*selfValues)/nSelfValueElements
	ENDELSE
	CorrValues_i = CorrValues_i/SQRT(refNormSq*selfNormSq_i)
ENDIF

CASE ndims OF
	2:	CorrValues[i, *] = CorrValues_i
	3:	CASE iaxis OF
		1:	CorrValues[i,*,*] = CorrValues_i
		2:	CorrValues[*,i,*] = CorrValues_i
		ENDCASE
	4:	CASE iaxis OF
		1:	CorrValues[i,*,*,*] = CorrValues_i
		2:	CorrValues[*,i,*,*] = CorrValues_i
		3:	CorrValues[*,*,i,*] = CorrValues_i
		ENDCASE
ENDCASE
ENDFOR
corrValuesPtr = PTR_NEW(CorrValues)
corrInfo = SIZE(corrValues)
ngrids(ndims) = corrInfo[ndims]
;
; Clean up ... free up selfValuePtr and refValuePtr
;
PTR_FREE, selfValuePtr
PTR_FREE, refValuePtr
;
; Make GKVsd object to contain cross correlation function
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
result -> GridRestrict, iaxis
;
; Start filling tags of result
;
result.mnemonic = 'XCorr_' + self.mnemonic	 + '_' + RefObj.mnemonic	; Set Mnemonic
result.title = 'C{' + self.title + ', ' + RefObj.title + '}'			; Set Title
result.units = '(' + self.units + ')*(' + refObj.units + ')'			; set units
IF(localNorm EQ 1) THEN result.units = ''
PTR_FREE, result.values								; Set Values pointer
result.values = corrValuesPtr
vmin = GKVsd_MIN(CorrValues, MAX=vmax)						; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange
PTR_FREE, result.ErrorBars					; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on correlation functions?
FOR i=1, ndims DO BEGIN						; Loop over ndims instances of the Grid Class
	IF(i NE iaxis) THEN BEGIN				; Change grid array for all axis EXCEPT 'iaxis'
		dimStr = STRING(i, FORMAT='(I1)')			
		GridStr = 'result.Grid' + dimStr
		command_str = 'Grid = ' + gridStr		; Copy ith instance of the Grid Class
		ok = EXECUTE(command_str)			;	 into Grid
		Mnemonic = Grid.mnemonic			; Set Grid Mnemoinc
		title = Grid.title				; Set Grid Title					
		IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
			Grid.mnemonic = 'tau'
			Grid.title = '!4s!X'
		ENDIF ELSE BEGIN
			Grid.title = '!4D!X' + title
		ENDELSE
		PTR_FREE, Grid.values						
		n = ngrids[i]
		gridvals = spacing(i)*(FINDGEN(n) -n/2)		; Generate array of grid values
		GridPtr = PTR_NEW(gridvals)	
		Grid.values = GridPtr				; Set Grid Values pointer
		Grid.boundary = 'periodic (open)'		; Correlation function BC's are always Periodic?		
		Grid.uniform = 1				; Correlation functions computed on uniform grid
		min = MIN(gridvals, Max=max)			; Set Grid plot Range
		range = [min, max]
		Grid.range = range
		irange = [0, n-1]				; Set Grid signal window range, irange
		Grid.irange = irange
		command_str = GridStr + ' = Grid'		; Copy ith instance of Grid Class
		ok = EXECUTE(command_str)			;	back into result
	ENDIF 
ENDFOR
RETURN, result	
END ; ****** GKVs2D::XCORR ****** ;


FUNCTION GKVs2D::XSPECT, 	Ref = RefObj, All=all, No_Avg=noAverage,  	$
				Lw=lagWindow,	iLw=iLagWindow, 		$
				Dw=dataWindow, iDw=iDataWindow, _Extra=extra
;
; Returns a GKVsd object containing the cross spectrum between the
; data in 'self' and the data in RefOjb (which BOTH must be GKVsd objects)
;
; If no RefObj is specified, then returns auto-spectrum of data in 'self'.
; 
; Remaining keywords are passed to XSPECT (see definition of XSPECT to understand their use)
;
;	Written 2/6/00 by W.M. Nevins
;
;	Revised 4/18/01 by W.M. Nevins
;	to fix bug relating to scaling of independent varialbe grids
;
rflag = 0
No_Avg = 1
IF(N_ELEMENTS(noAverage)) THEN No_Avg = noAverage
ndims = self -> NumDims()						; Get dimensionaliity of data in 'self'
IF((ndims LT 1) OR (ndims GT 4)) THEN BEGIN
	MESSAGE, 'Bad dimensionality, cannot form cross correlations function', /Informational
	RETURN, 0
ENDIF
IF(N_ELEMENTS(RefObj) EQ 0) THEN BEGIN			; No reference object... so pass on to GKVs1D
	result = self -> GKVs1D::XSPECT(No_Avg=noAverage, Lw=lagWindow, iLw=iLagWindow, Dw=dataWindow, iDw=iDataWindow, _Extra=extra)
	RETURN, result
ENDIF
refType = TypeOf(RefObj)
IF (refType NE 11) THEN BEGIN				; RefObj is NOT an object
	MESSAGE, 'Ref is not an Object', /Informational
	RETURN, 0
ENDIF
refdims = RefObj -> NumDims()					; Get dimensionality of data in 'RefObj'
IF(refdims EQ ndims) THEN BEGIN					; Dimensionality is same... so pass on to GKVs1D
	result = self -> GKVs1D::XSPECT(Ref = RefObj, No_Avg=noAverage, Dw=dataWindow, _Extra=extra)
	RETURN, result
ENDIF
IF(refdims NE (ndims-1) ) THEN BEGIN
	MESSAGE, "Expect dimensionality of Ref = dimensionality 'self', or self - 1" , /informational
	RETURN, 0
ENDIF
spacing = FLTARR(ndims+1)
ngrids  = LONARR(ndims+1)
iref = 1
iaxis = 0
FOR iself=1, ndims DO BEGIN						; 
	iselfStr = STRING(iself, FORMAT='(I1)')
	command_str = 'selfGrid = self.Grid' + iselfStr		; Get iself^th Grid structure of 'self'
	ok = EXECUTE(command_str)
	irefstr = STRING(iref, FORMAT='(I1)')
	command_str = 'refGrid = RefObj.Grid' + irefstr
	ok = EXECUTE(command_str)						; Get iref^th Grid structure of 'RefObj'
	ismin = selfGrid.irange[0]
	ismax = selfGrid.irange[1]
	IF(GKV_GridSame(selfGrid,RefGrid,All=all, /force)) THEN BEGIN	; Check to see if the Grid structures match (force match if possible)
		IF(KEYWORD_SET(All)) THEN BEGIN
			ismin = 0
			ismax = N_ELEMENTS(*selfGrid.values) - 1
		ENDIF
		uniform = GKVsd_UniformGrid((*selfGrid.values)[ismin:ismax])
		IF(NOT uniform) THEN BEGIN					; Check for uniform grids
			MESSAGE, 'XCORR only defined for data on uniform grid(s)', /Informational
			RETURN, 0
		ENDIF
		spacing(iself)	= 2.0*!PI/( (*selfGrid.values)[ismax] - (*selfGrid.values)[ismin] )
		ngrids(iself)	= ismax-ismin+1
		iref = iref + 1
	ENDIF ELSE BEGIN							; This selfGrid does not match corresponding refGrid
		IF(iaxis NE 0) THEN BEGIN						;	this should only happen once!
			MESSAGE, 'Bad match between reference grids and self grids', /Informational
			RETURN, 0
		ENDIF
		iaxis= iself							; Prepare for loop which calls XCORR
		imin = 0
		imax = ismax - ismin
	ENDELSE
	IF(iself EQ nDims) THEN dt = (*selfGrid.values)[ismin+1] - (*selfGrid.values)[ismin]
ENDFOR
IF(iaxis EQ 0) THEN BEGIN
	MESSAGE, "Couldn't find cross correlation axis", /Informational
	RETURN, 0
ENDIF

dt = spacing(ndims)
;
; Set default lag window and data windows
;
iLw = ngrids[nDims]/2 > 1
iDw = ngrids[nDims]/10 < SQRT(ngrids[nDims])
iDw = iDw > 1
IF(N_ELEMENTS( lagWindow) GT 0) THEN iLw=FIX( lagWindow/dt)
IF(N_ELEMENTS(dataWindow) GT 0) THEN iDw=FIX(dataWindow/dt)

selfValuePtr =  self -> GetValues(All=all)				; Returns pointer to values within signal window
refValuePtr = RefObj -> GetValues(All=all)
refValues = *refValuePtr
refInfo = SIZE(refValues)
FOR i = imin, imax DO BEGIN
CASE ndims OF
	2:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *], refInfo[1])
			2:	selfValues = REFORM((*selfValuePtr)[*, i], refInfo[1])
		ENDCASE
	3:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *], refInfo[1], refInfo[2])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *], refInfo[1], refinfo[2])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i], refInfo[1], refInfo[2])
		ENDCASE
	4:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *, *], refInfo[1], refInfo[2], refInfo[3])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *, *], refInfo[1], refInfo[2], refInfo[3])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i, *], refInfo[1], refInfo[2], refInfo[3])
			4:	selfValues = REFORM((*selfValuePtr)[*, *, *, i], refInfo[1], refInfo[2], refInfo[3])
		ENDCASE
ENDCASE
SpectValues_i = XSpect(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, _Extra=extra, Lw=iLw, Dw=iDw)
IF(i EQ imin) THEN BEGIN							; Initialization duties
	selfInfo = SIZE(*selfValuePtr)
	spectInfo = SIZE(SpectValues_i)
	spectDims = spectInfo[0]
	nOmegas = spectinfo[spectDims]
	CASE ndims OF
		2:	spectValues = COMPLEXARR(selfInfo[1], nOmegas)
		3:	spectValues = COMPLEXARR(selfinfo[1], selfInfo[2], nOmegas)
		4:	spectValues = COMPLEXARR(selfinfo[1], selfInfo[2], selfInfo[3], nOmegas)
	ENDCASE
ENDIF
CASE ndims OF
	2:	spectValues[i, *] = spectValues_i
	3:	CASE iaxis OF
		1:	spectValues[i,*,*] = spectValues_i
		2:	spectValues[*,i,*] = spectValues_i
		ENDCASE
	4:	CASE iaxis OF
		1:	spectValues[i,*,*,*] = spectValues_i
		2:	spectValues[*,i,*,*] = spectValues_i
		3:	spectValues[*,*,i,*] = spectValues_i
		ENDCASE
ENDCASE
ENDFOR
PTR_FREE, selfValuePtr							; Free up pointers
PTR_FREE, refValuePtr
spectValuesPtr = PTR_NEW(spectValues)
;spectInfo = SIZE(spectValues)
;ngrids[ndims] = spectInfo[ndims]
;
; Make GKVsd object to contain cross spectrum
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
result -> GridRestrict, iaxis
;
; Start filling tags of result
;
result.mnemonic = 'XSpect_' + self.mnemonic	 + '_' + RefObj.mnemonic	; Set Mnemonic
result.title = 'S{' + self.title + ', ' + RefObj.title + '}'			; Set Title
result.units = '(' + result.units + ')*(' + refObj.units + ')'			; Set units (but see addtions in ndim loop below)
PTR_FREE, result.ErrorBars					; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on spectral density?
spectNorm = 1.
FOR i=1, ndims DO BEGIN						; Loop over ndims instances of the Grid Class
	IF(i NE iaxis) THEN BEGIN				; Change grid array for all axis EXCEPT 'iaxis'
		dimStr = STRING(i, FORMAT='(I1)')			
		GridStr = 'result.Grid' + dimStr
		command_str = 'Grid = ' + gridStr		; Copy ith instance of the Grid Class
		ok = EXECUTE(command_str)			;	 into Grid
		result.units = result.units + '*(' + Grid.units + ')'
		Mnemonic = Grid.mnemonic			; Get Grid Mnemonic
		norm = 1.
		IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
			Grid.mnemonic = 'omega'
			norm = FLOAT(ngrids[i]-1)/nOmegas	; Constant to correct dOmega
			Grid.title = '!4x!X'
		ENDIF ELSE BEGIN
			Grid.mnemonic = 'k_' + Grid.mnemonic
			Grid.title = 'k!I' + grid.title + '!N!X'
		ENDELSE
		Grid.units = '1/' + Grid.units			; Update grid units
		PTR_FREE, Grid.values						
		n = ngrids[i]
		IF(i EQ ndims) THEN n=nOmegas
		dGrid = norm*spacing(i)
		spectNorm = spectNorm/dGrid
		gridvals = dGrid*(FINDGEN(n) -n/2)		; Generate array of grid values
		GridPtr = PTR_NEW(gridvals)	
		Grid.values = GridPtr				; Set Grid Values pointer
		Grid.boundary = 'Periodic'			; Spectral density BC's are always Periodic?		
		Grid.uniform = 1				; Spectral density computed on uniform grid
		min = MIN(gridvals, Max=max)			; Set Grid plot Range
		range = [min, max]
		Grid.range = range
		irange = [0, n-1]				; Set Grid signal window range, irange
		Grid.irange = irange
		command_str = GridStr + ' = Grid'		; Copy ith instance of Grid Class
		ok = EXECUTE(command_str)			;	back into result
	ENDIF
ENDFOR

PTR_FREE, result.values								; Set Values pointer
result.values = PTR_NEW(spectNorm*(*spectValuesPtr))
PTR_FREE, spectValuesPtr
vmin = GKVsd_MIN(spectValues, MAX=vmax)						; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange


RETURN, result	
END ; ****** GKVs2D::XSPECT ****** ;



Pro GKVs2D::Draw, _Extra=extra
;
;  'Virtual' keywords:
;
;			  xrange=x_range,  title=T_title, xtitle=x_title, ytitle=y_title, 		$
;			  yrange=y_range, Grid1=Grid_1, Grid2=Grid_2, indx1=indx_1, indx2=indx_2,	$
;			  Pretty=pretty, Polar=polar, Log=log, Vrange=vrange, Bottom
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
result = GetKeyWord('Grid1', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_1=result	;	Grid1 = Grid_1
IF(TypeOf(result) NE 7) THEN Grid_1=result
result = GetKeyWord('Grid2', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_2=result	; 	Grid2 = Grid_2
IF(TypeOf(result) NE 7) THEN Grid_2=result
result = GetKeyWord('indx1', extra)
IF(TypeOf(result) NE 7) THEN indx_1=result			; 	indx1 = indx_1
result = GetKeyWord('indx2', extra)
IF(TypeOf(result) NE 7) THEN indx_2=result			; 	indx2	 = indx_2
result = GetKeyWord('Pretty', extra)
IF(TypeOf(result) NE 7) THEN pretty=result			;	Pretty = pretty
result = GetKeyWord('Polar', extra)
IF(TypeOf(result) NE 7) THEN polar=result			;	Polar = polar
result = GetKeyWord('Log', extra)
IF(TypeOF(result) NE 7) THEN Log=result			;	Log = log
result = GetKeyWord('Vrange', extra)
IF(TypeOF(result) NE 7) THEN vrange=result			;	Vrange=vrange
result = GetKeyWord(self.mnemonic, extra)			;	or, set vrange
IF(TypeOF(result) NE 7) THEN vrange=result			;	using 'mnemonic' = [vmin, vmax]  on the command line
bottom = !COLOR_SETUP_NCOLORS + 1
result = GetKeyWord('bottom', extra)
IF(TypeOF(result) NE 7) THEN bottom = result		;	Bottom = bottom
ncolors = !D.TABLE_SIZE - bottom
result = GetKeyWord('ncolors', extra)				;	ncolors = ncolors
IF(TypeOF(result) NE 7) THEN ncolors = result		;
noAvg=0
result = GetKeyWord('noavg', extra)   ; no avg
IF(TypeOF(result) NE 7) THEN noAvg=result

Grid1default = 1
Grid2defalut = 2
;
; Find Grid structure for first axis
;
IF(N_ELEMENTS(Grid_1) EQ 0) THEN Grid_1=1			; default value
IF(TypeOF(Grid_1) EQ 7) THEN Grid_1=(self -> AxisNumber(Grid_1)) > 1 
str_1 = STRING(Grid_1, FORMAT='(I1)')
command_str = 'Grid1=self.Grid' + str_1
ok = EXECUTE(command_str)
;
; Find Grid structure for second axis
;
IF(N_ELEMENTS(Grid_2) EQ 0) THEN Grid_2=Grid_1+1	; default value
IF(TypeOF(Grid_2) EQ 7) THEN Grid_2=(self -> AxisNumber(Grid_2)) > (Grid_1+1) 
str_2 = STRING(Grid_2, FORMAT='(I1)')
command_str = 'Grid2=self.Grid' + str_2
ok = EXECUTE(command_str)
ndims = self -> NumDims()			; Get number of dimensions
plotAxis=INTARR(ndims+1)
plotAxis[Grid_1]=1
plotAxis[Grid_2]=2
indx=-1
FOR iaxis=1, ndims DO BEGIN
	IF(plotAxis[iaxis] GT 0) THEN GOTO, DONE01
	plotAxis[iaxis] = indx
	indx=indx-1
DONE01 :  ; moral equivalent of 'continue'? 
ENDFOR
	
;
; Check command line for axis mnemonics
;
axisInfo = self -> GetAxis(extra)
;
; 
FOR iaxis=1, ndims DO BEGIN
	axisValue=axisInfo.(iaxis-1)
	IF(TypeOF(axisValue) NE 7) THEN BEGIN	; Mnemonic for iaxis is used on command line
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -1)) THEN indx_1= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -2)) THEN indx_2= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 1)) THEN x_range=AxisValue
		IF((N_Elements(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 2)) THEN y_range=AxisValue
	ENDIF
ENDFOR
;
; put 'indx_1', 'indx_2' into plot axis if these are set.
;
IF(N_ELEMENTS(indx_1) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -1) THEN plotAxis[iaxis]=-indx_1
	ENDFOR
ENDIF
IF(N_ELEMENTS(indx_2) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -2) THEN plotAxis[iaxis]=-indx_2
	ENDFOR
ENDIF
; 
; Get limits to signal window for each axis
;
imin=Grid1.irange[0]
imax=Grid1.irange[1]
jmin=Grid2.irange[0]
jmax=Grid2.irange[1]
;
; use default plot keyword values from Object definition if they
; are not over ridden on command line.
;
xrange=Grid1.range
IF(N_ELEMENTS(x_range) EQ 2) THEN xrange = x_range

xmin=xrange[0]
xmax=xrange[1]
nx = N_ELEMENTS(*Grid1.values)
temp=(*Grid1.values - xmin)^2			; Finds iimin such that Grid1.values(iimin)
smallest = MIN(temp, iimin)			;	takes value closest to xmin
IF(iimin eq imax) THEN iimin=imin
imin = imin > iimin > 0
xmin=(*Grid1.values)[imin]

temp=(*Grid1.values - xmax)^2			; sets iimax such that Grid1.values(iimax)
smallest = MIN(temp, iimax)			;	takes value closest to xmax
IF(iimax eq imin) THEN iimax=imax
imax = imax < iimax < (nx-1)
xmax=(*Grid1.values)[imax]
xrange = [xmin, xmax]

indices = self -> IndexString(plotAxis, pretty=pretty)
indexStr = '[' + STRJOIN(indices, ', ') + ']'
IF(KEYWORD_SET(pretty)) THEN BEGIN
	title=self.title + indexStr + " (" + self.units + ")"
ENDIF ELSE BEGIN
	title=self.mnemonic + indexStr + " (" + self.units + ")"
ENDELSE
IF KEYWORD_SET(T_title) THEN title=T_title

x_title=Grid1.mnemonic + " (" + Grid1.units + ")"
IF KEYWORD_SET(Pretty)  THEN x_title=Grid1.title + " (" + Grid1.units + ")"
IF KEYWORD_SET(x_title) THEN xtitle=x_title

ytitle=Grid2.mnemonic + " (" + Grid2.units + ")"
IF KEYWORD_SET(Pretty)  THEN ytitle=Grid2.title + " (" + Grid2.units + ")"
IF KEYWORD_SET(y_title) THEN ytitle=y_title

yrange=Grid2.range
IF (N_ELEMENTS(y_range) EQ 2) THEN yrange = y_range

ymin=yrange[0]
ymax=yrange[1]
ny = N_ELEMENTS(*Grid2.values)
temp=(*Grid2.values - ymin)^2			; Finds jjmin such that Grid2.values(iimin)
smallest = MIN(temp, jjmin)			;	takes value closest to ymin
IF(jjmin eq jmax) THEN jjmin=jmin
jmin = jmin > jjmin > 0
ymin=(*Grid2.values)[jmin]

temp=(*Grid2.values - ymax)^2			; sets jjmax such that Grid2.values(iimax)
smallest = MIN(temp, jjmax)			;	takes value closest to ymax
IF(jjmax eq jmin) THEN jjmax=jmax
jmax = jmax < jjmax < (ny-1)
ymax=(*Grid2.values)[jmax]
yrange = [ymin, ymax]

x = (*Grid1.values )[imin:imax]
y = (*Grid2.values )[jmin:jmax]
;
; POLYMORPHISM!  Make this work for subclasses of GKVs2D
;
indx1=0
indx2=0
IF(N_ELEMENTS(indx_1)) THEN indx1=indx_1
IF(N_ELEMENTS(indx_2)) THEN indx2=indx_2
CASE ndims OF
	2:	values = (*self.values)[imin:imax, jmin:jmax]
	3:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1]	
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=3', /Informational
				ENDCASE
				END
				
			2:	values = (*self.values)[indx1, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, "Bad value for Grid_1, ndims=3", /Informational
		ENDCASE
		END 
	4:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1, indx2] 
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax, indx2]
					4	:	values = (*self.values)[imin:imax, indx1, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=4', /Informational
				ENDCASE
				END
			2:	BEGIN
				CASE Grid_2 OF
					3	:	values = (*self.values)[indx1, imin:imax, jmin:jmax, indx2] 
					4	:	values = (*self.values)[indx1, imin:imax, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=2, ndims=4', /Informational
				ENDCASE
				END
			3:	 values = (*self.values)[indx1, indx2, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, 'Bad value for Grid_1, ndims=4', /Informational
		ENDCASE
		END
	ELSE	:	MESSAGE, 'Bad value for ndims', /Informational
ENDCASE
values = FLOAT(REFORM(values, imax-imin+1, jmax-jmin+1, /OVERWRITE))
info=SIZE(values)
index=info[0]+1
vtype=info[index]

IF(KEYWORD_SET(Polar)) THEN BEGIN
	Nrad = N_ELEMENTS(x)
	Rmax = MAX(x)
	dx = Rmax/Nrad
	spacing=[dx,dx]
	bounds=[-Rmax, -Rmax, Rmax, Rmax]
	values = POLAR_SURFACE(values, x, y, /GRID, SPACING=spacing, BOUNDS=bounds)
	info = SIZE(values)
	index=info[0]+1
	vtype=info[index]
	Nx = info[1]
	x = -Rmax + dx*FINDGEN(NX)
	y = x
	xrange = [-Rmax, Rmax]
	yrange = xrange
	xtitle = 'R-R!Io!N (a)'
	ytitle = 'Z (a)'
ENDIF

vmin = self.vrange[0]
vmax = self.vrange[1]
IF(N_ELEMENTS(vrange) EQ 2) THEN BEGIN
	vmin=vrange[0]
	vmax=vrange[1]
ENDIF
vvmin = vmin
vvmax = vmax
IF(KEYWORD_SET(Log)) THEN BEGIN
	values = ALOG10(values)
	vvmax = ALOG10(vmax)
	vvmin = vvmax-12
	IF(vmin GT 0) THEN vvmin = ALOG10(vmin)
	vmin = exp(vvmin*ALOG(10.))
ENDIF

IF(N_ELEMENTS(extra) NE 0) THEN BEGIN
	IF(TypeOf(extra) EQ 8) THEN nextra = extra
ENDIF
;
; Set up 'image' array
;
ncolors = BYTE(ncolors)
bottom  = BYTE(bottom)
IF(ncolors GE (256 - bottom)) THEN ncolors = 256 - bottom
image = BYTSCL(values, Min=vvmin, Max=vvmax, Top=ncolors-1, /NaN) + bottom
;
; Display image
;
TVImage, image, position=[0.15, 0.15, 0.9, 0.8], /Erase, /Minus_One, _Extra=nextra 
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
; Draw axis, etc. over image
;

IF(vmin EQ vmax) THEN values[0,0]=values[0,0] +MAX(values) + 1.	
IF( NOT FINITE(vmax)) THEN BEGIN
	values = FLOAT(values) < 100.*vmin
	values[0,0]=values[0,1]+1.
ENDIF
;
; These 'IF's' are there to insure that the following 'CONTOUR' (which only draws the axises)
;	does not fail when there is no variation in 'values'.
; 
CONTOUR, 	values, x, y, position=[0.15, 0.15, 0.9, 0.8], 	$
		xrange=xrange, yrange=yrange,					$
		xtitle=xtitle, ytitle=ytitle, XStyle=1, YStyle=1, 	$
		Ticklen=0.02, /NODATA, /NOErase, _Extra=nextra
;
; Add color bar
;
ColorBar, Ncolors=ncolors, bottom = bottom, range=[vmin,vmax], 	$
		Format='(G10.2)', position=[0.145, 0.85, 0.9, 0.9], Xlog=log
;
; Mark average value on color bar
;
IF (NOT KEYWORD_SET(Log))AND (noAvg EQ 0) THEN BEGIN
	avg=TOTAL(values)/N_elements(values)
	ARROW,  avg, 0.45*(vmax-vmin)+vmin,  avg, vmin, Thick=2., /Data
	XYOUTS, avg, 0.55*(vmax-vmin)+vmin, 'AVG', Alignment=0.5, /Data
ENDIF
;
; Put title at top of Graphics window
;

result = GetKeyWord("NoTitle", Extra)
IF(Query_Integer(result)) THEN GOTO, NoTitle
charSize = GKV_TitleSize(title, Width = (0.9-0.15), yPosition=ywrite)
xwrite = (0.15 + 0.9)/2.*!D.X_SIZE
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
XYOUTS, xwrite, ywrite, CodeName, /Device		; Write CodeName to lower left-hand corne.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
codePI = STRMID(self.CodePI, 0, maxChars)		; Truncate 'CodePI' if necessary
XYOUTS, xwrite, ywrite, CodePI, /Device			; Write CodePI below CodeName.
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window.
runID = STRMID(self.RunID, 0, maxChars)			; Truncate 'runID' if necessary
XYOUTS, xwrite, ywrite, RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
fileID = STRMID(self.FileID, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID
;
; Return character size to input values
;
NoFooter:  DEVICE, SET_CHARACTER_SIZE = oldSizes

RETURN
END ; ***** GKVs2D::Draw ***** ;


Pro GKVs2D::Shade_Surf, _Extra=extra
;
;  'Virtual' keywords:
;
;			  xrange=x_range,  title=T_title, xtitle=x_title, ytitle=y_title, 		$
;			  yrange=y_range, Grid1=Grid_1, Grid2=Grid_2, indx1=indx_1, indx2=indx_2,	$
;			  zrange=z_range, ztitle=z_title, position=pstn, charsize=chsz,	$
;			  Pretty=pretty, Polar=polar, Log=log, Vrange=vrange
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

result = GetKeyWord('ztitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN z_title=result	;	ytitle = y_title
result = GetKeyWord('zrange', extra)
IF(TypeOf(result) NE 7) THEN z_range=result		;	yrange = y_range

result = GetKeyWord('position', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN pstn=result		;	position = pstn
result = GetKeyWord('charsize', extra)
IF(TypeOf(result) NE 7) THEN chsz=result			;	charsize = chsz


result = GetKeyWord('Grid1', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_1=result	;	Grid1 = Grid_1
IF(TypeOf(result) NE 7) THEN Grid_1=result
result = GetKeyWord('Grid2', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_2=result	; 	Grid2 = Grid_2
IF(TypeOf(result) NE 7) THEN Grid_2=result
result = GetKeyWord('indx1', extra)
IF(TypeOf(result) NE 7) THEN indx_1=result			; 	indx1 = indx_1
result = GetKeyWord('indx2', extra)
IF(TypeOf(result) NE 7) THEN indx_2=result			; 	indx2	 = indx_2
result = GetKeyWord('Pretty', extra)
IF(TypeOf(result) NE 7) THEN pretty=result			;	Pretty = pretty
result = GetKeyWord('Polar', extra)
IF(TypeOf(result) NE 7) THEN polar=result			;	Polar = polar
result = GetKeyWord('Log', extra)
IF(TypeOF(result) NE 7) THEN Log=result			;	Log = log
result = GetKeyWord('Vrange', extra)
IF(TypeOF(result) NE 7) THEN vrange=result			;	Vrange=vrange
result = GetKeyWord(self.mnemonic, extra)			;	or, set vrange
IF(TypeOF(result) NE 7) THEN vrange=result			;	using 'mnemonic' = [vmin, vmax] on the command line
bottom = !COLOR_SETUP_NCOLORS + 1
result = GetKeyWord('bottom', extra)				;   Bottom = bottom
IF(TypeOF(result) NE 7) THEN bottom = result		; 
ncolors = !D.TABLE_SIZE - bottom
result = GetKeyWOrd('ncolors', extra)				;	ncolors = ncolors
IF(TypeOF(result) NE 7) THEN ncolors = result		;
phaseShading=1
result = GetKeyWord('PhaseShading', extra)			;	PhaseShading = phaseShading
IF(TypeOf(result) NE 7) THEN phaseShading=result	;
result = GetKeyWord('Shades', extra)
IF(TypeOF(result) NE 7) THEN shades = result		;	Shades = shades

Grid1default = 1
Grid2defalut = 2
;
; Find Grid structure for first axis
;
IF(N_ELEMENTS(Grid_1) EQ 0) THEN Grid_1=1			; default value
IF(TypeOF(Grid_1) EQ 7) THEN Grid_1=(self -> AxisNumber(Grid_1)) > 1 
str_1 = STRING(Grid_1, FORMAT='(I1)')
command_str = 'Grid1=self.Grid' + str_1
ok = EXECUTE(command_str)
;
; Find Grid structure for second axis
;
IF(N_ELEMENTS(Grid_2) EQ 0) THEN Grid_2=Grid_1+1	; default value
IF(TypeOF(Grid_2) EQ 7) THEN Grid_2=(self -> AxisNumber(Grid_2)) > (Grid_1+1) 
str_2 = STRING(Grid_2, FORMAT='(I1)')
command_str = 'Grid2=self.Grid' + str_2
ok = EXECUTE(command_str)
ndims = self -> NumDims()			; Get number of dimensions
plotAxis=INTARR(ndims+1)
plotAxis[Grid_1]=1
plotAxis[Grid_2]=2
indx=-1
FOR iaxis=1, ndims DO BEGIN
	IF(plotAxis[iaxis] GT 0) THEN GOTO, DONE01
	plotAxis[iaxis] = indx
	indx=indx-1
DONE01 :  ; moral equivalent of 'continue'? 
ENDFOR
	
;
; Check command line for axis mnemonics
;
axisInfo = self -> GetAxis(extra)
;
; 
FOR iaxis=1, ndims DO BEGIN
	axisValue=axisInfo.(iaxis-1)
	IF(TypeOF(axisValue) NE 7) THEN BEGIN	; Mnemonic for iaxis is used on command line
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -1)) THEN indx_1= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -2)) THEN indx_2= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 1)) THEN x_range=AxisValue
		IF((N_Elements(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 2)) THEN y_range=AxisValue
	ENDIF
ENDFOR
;
; put 'indx_1', 'indx_2' into plot axis if these are set.
;
IF(N_ELEMENTS(indx_1) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -1) THEN plotAxis[iaxis]=-indx_1
	ENDFOR
ENDIF
IF(N_ELEMENTS(indx_2) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -2) THEN plotAxis[iaxis]=-indx_2
	ENDFOR
ENDIF
; 
; Get limits to signal window for each axis
;
imin=Grid1.irange[0]
imax=Grid1.irange[1]
jmin=Grid2.irange[0]
jmax=Grid2.irange[1]
;
; use default plot keyword values from Object definition if they
; are not over ridden on command line.
;
xrange=Grid1.range
IF(N_ELEMENTS(x_range) EQ 2) THEN xrange = x_range

xmin=xrange[0]
xmax=xrange[1]
nx = N_ELEMENTS(*Grid1.values)
temp=(*Grid1.values - xmin)^2			; Finds iimin such that Grid1.values(iimin)
smallest = MIN(temp, iimin)			;	takes value closest to xmin
IF(iimin eq imax) THEN iimin=imin
imin = imin > iimin > 0
xmin=(*Grid1.values)[imin]

temp=(*Grid1.values - xmax)^2			; sets iimax such that Grid1.values(iimax)
smallest = MIN(temp, iimax)			;	takes value closest to xmax
IF(iimax eq imin) THEN iimax=imax
imax = imax < iimax < (nx-1)
xmax=(*Grid1.values)[imax]
xrange = [xmin, xmax]

indices = self -> IndexString(plotAxis, Pretty=pretty)
indexStr = '[' + STRJOIN(indices, ', ') + ']'
IF(KEYWORD_SET(pretty)) THEN BEGIN
	title=self.title + indexStr + " (" + self.units + ")"
ENDIF ELSE BEGIN
	title=self.mnemonic + indexStr + " (" + self.units + ")"
ENDELSE
IF KEYWORD_SET(T_title) THEN title=T_title

x_title=Grid1.mnemonic + " (" + Grid1.units + ")"
IF KEYWORD_SET(Pretty)  THEN x_title=Grid1.title + " (" + Grid1.units + ")"
IF KEYWORD_SET(x_title) THEN xtitle=x_title

ytitle=Grid2.mnemonic + " (" + Grid2.units + ")"
IF KEYWORD_SET(Pretty)  THEN ytitle=Grid2.title + " (" + Grid2.units + ")"
IF KEYWORD_SET(y_title) THEN ytitle=y_title

yrange=Grid2.range
IF (N_ELEMENTS(y_range) EQ 2) THEN yrange = y_range

ymin=yrange[0]
ymax=yrange[1]
ny = N_ELEMENTS(*Grid2.values)
temp=(*Grid2.values - ymin)^2			; Finds jjmin such that Grid2.values(iimin)
smallest = MIN(temp, jjmin)			;	takes value closest to ymin
IF(jjmin eq jmax) THEN jjmin=jmin
jmin = jmin > jjmin > 0
ymin=(*Grid2.values)[jmin]

temp=(*Grid2.values - ymax)^2			; sets jjmax such that Grid2.values(iimax)
smallest = MIN(temp, jjmax)			;	takes value closest to ymax
IF(jjmax eq jmin) THEN jjmax=jmax
jmax = jmax < jjmax < (ny-1)
ymax=(*Grid2.values)[jmax]
yrange = [ymin, ymax]

x = (*Grid1.values )[imin:imax]
y = (*Grid2.values )[jmin:jmax]
;
; POLYMORPHISM!  Make this work for subclasses of GKVs2D
;
indx1=0
indx2=0
IF(N_ELEMENTS(indx_1)) THEN indx1=indx_1
IF(N_ELEMENTS(indx_2)) THEN indx2=indx_2
CASE ndims OF
	2:	values = (*self.values)[imin:imax, jmin:jmax]
	3:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1]	
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=3', /Informational
				ENDCASE
				END
				
			2:	values = (*self.values)[indx1, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, "Bad value for Grid_1, ndims=3", /Informational
		ENDCASE
		END 
	4:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1, indx2] 
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax, indx2]
					4	:	values = (*self.values)[imin:imax, indx1, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=4', /Informational
				ENDCASE
				END
			2:	BEGIN
				CASE Grid_2 OF
					3	:	values = (*self.values)[indx1, imin:imax, jmin:jmax, indx2] 
					4	:	values = (*self.values)[indx1, imin:imax, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=2, ndims=4', /Informational
				ENDCASE
				END
			3:	 values = (*self.values)[indx1, indx2, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, 'Bad value for Grid_1, ndims=4', /Informational
		ENDCASE
		END
	ELSE	:	MESSAGE, 'Bad value for ndmis', /Informational
ENDCASE
values = REFORM(values, imax-imin+1, jmax-jmin+1, /OVERWRITE)
info=SIZE(values)
index=info[0]+1
vtype=info[index]

IF(KEYWORD_SET(Polar)) THEN BEGIN
	Nrad = N_ELEMENTS(x)
	Rmax = MAX(x)
	dx = Rmax/Nrad
	spacing=[dx,dx]
	bounds=[-Rmax, -Rmax, Rmax, Rmax]
	values = POLAR_SURFACE(values, x, y, /GRID, SPACING=spacing, BOUNDS=bounds)
	info = SIZE(values)
	index=info[0]+1
	vtype=info[index]
	Nx = info[1]
	x = -Rmax + dx*FINDGEN(NX)
	y = x
	xrange = [-Rmax, Rmax]
	yrange = xrange
	xtitle = 'R-R!Io!N (a)'
	ytitle = 'Z (a)'
ENDIF

ncolors=220
vmin = self.vrange[0]
vmax = self.vrange[1]
vvmin = vmin
vvmax = vmax
IF(N_ELEMENTS(vrange) EQ 2) THEN BEGIN
	vmin=vrange[0]
	vmax=vrange[1]
ENDIF
IF(KEYWORD_SET(Log)) THEN BEGIN
;	values = ALOG10(values)
	vvmax = ALOG10(vmax)
	vvmin = vvmax-12
	IF(vmin GT 0) THEN vvmin = ALOG10(vmin)
	vmin = exp(vvmin*ALOG(10.))
ENDIF

zrange=[vmin, vmax]
IF(N_ELEMENTS(z_range) EQ 2) THEN zrange=z_range

ztitle=self.mnemonic + " (" + self.units + ")"
IF KEYWORD_SET(Pretty)  THEN ztitle=self.title + " (" + self.units + ")"
IF KEYWORD_SET(z_title) THEN ztitle=z_title

position=[0.2, 0.15, 0.95, 0.9] 
IF KEYWORD_SET(pstn) THEN position=pstn
charsize=2.
IF KEYWORD_SET(chsz) THEN charsize=chsz

IF(N_ELEMENTS(extra) NE 0) THEN BEGIN
	IF(TypeOf(extra) EQ 8) THEN nextra = extra
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

IF( phaseShading AND Query_Complex(values) ) THEN BEGIN
	cosinePhase = FLOAT(values)
	values = ABS(values)
	IF(vmin LT 0) THEN vmin = 0
	; vmax = 2*vmax
	IF(N_ELEMENTS(z_range) NE 2) THEN zrange = [vmin, vmax]
	cosinePhase = cosinePhase/values
	;
	; Set up 'phaseShading' array
	;
	ncolors = BYTE(ncolors)
	bottom  = BYTE(bottom)
	IF(ncolors GE (256 - bottom)) THEN ncolors = 256 - bottom
	phaseShading = BYTSCL(cosinePhase, Min=-1, Max=1, Top=ncolors-1, _Extra=nextra, /NaN) + bottom
	SHADE_SURF, values, x, y, xrange=xrange, yrange=yrange, max_value=vmax, min_value=vmin, zrange=zrange, $
		 	xtitle=xtitle, ytitle=ytitle, zlog=log, $
		 	Shades=phaseShading, position=position, charsize=charsize, xstyle=1, ystyle=1, _Extra=nextra 
	GOTO, DONE_SHADE_SURF
ENDIF 
IF(N_ELEMENTS(shades) NE 0) THEN BEGIN
	IF(TypeOF(shades) EQ 11) THEN BEGIN			; 'shades' is an object, so must extract values
		shadeDims = shades -> NumDims()
		IF(Numdims NE 2) THEN BEGIN
			MESSAGE, "Can't shade with objects of dimensionality greater than 2-D", /INFORMATIONAL
			GOTO, noShades
		ENDIF
		nshades = shades -> INTERPOLATE(self)
		IF(nshades EQ 0) THEN BEGIN
			MESSAGE, "Can't interpolate object onto self's Grid", /INFORMATIONAL
			GOTO, noShades 
		ENDIF
		shades = *(nshades -> GetValues())
		shadeRange = nshades.vrange
		nshades -> trash		
	ENDIF
	shades = FIX(shades)
	IF(N_ELEMENTS(shadeRange) NE 2) THEN BEGIN
		shadeMin = MIN(shades, MAX=shadeMax)
		shadeRange = [shadeMin, shadeMax]
	ENDIF
	bottom = !COLOR_SETUP_NCOLORS + 1
	ncolors = !D.TABLE_SIZE - bottom
	shades = shades < shadeRange[1]
	shades = shades > shadeRange[0]
	shades = LONG(bottom + ( (shades-shadeRange[0])/(shadeRange[1]-shadeRange[0]) )*ncolors)
	SHADE_SURF, values, x, y, xrange=xrange, yrange=yrange, max_value=vmax, min_value=vmin, zrange=zrange, $
		 	xtitle=xtitle, ytitle=ytitle, zlog=log, shades=shades, $
		 	position=position, charsize=charsize, xstyle=1, ystyle=1, _Extra=nextra 
ENDIF ELSE BEGIN 
noShades:
	SHADE_SURF, values, x, y, xrange=xrange, yrange=yrange, max_value=vmax, min_value=vmin, zrange=zrange, $
		 	xtitle=xtitle, ytitle=ytitle, zlog=log, $
		 	position=position, charsize=charsize, xstyle=1, ystyle=1, _Extra=nextra 
ENDELSE
DONE_SHADE_SURF:
;
; Put title at top of Graphics window
;
charSize = GKV_TitleSize(title, Width = (0.9-0.15), yPosition=ywrite)
xwrite = (0.15 + 0.9)/2.*!D.X_SIZE
XYOUTS, xwrite, ywrite, title, ALIGNMENT=0.5, CHARSIZE=charSize, /DEVICE 
;
; Now write code and run info in lower corners
;
xwrite=!D.x_ch_size						; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower left-hand corner of plot window
XYOUTS, xwrite, ywrite, self.CodeName, /Device	; Write CodeName to lower left-hand corner
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines
XYOUTS, xwrite, ywrite, self.CodePI, /Device	; Write CodePI below CodeName
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window
XYOUTS, xwrite, ywrite, self.RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines
XYOUTS, xwrite, ywrite, self.FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID
;
; Return character size to input values
;
DEVICE, SET_CHARACTER_SIZE = oldSizes

RETURN
END ; ****** GKVs2D::Shade_Surf ****** ;


PRO GKVs2D::GET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices,		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Get values of elements of the realization of the GRID class, Grid2;
; and then call GKVs1D::GET to get values of elements of the GKVs1D Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs2D,
; or any or its subclasses
;
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 2) THEN 	$
			GKVsd_GetGrid, self.Grid2,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	2:	IF(axisID EQ 2) THEN 	$
			GKVsd_GetGrid, self.Grid2,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	3:	IF(axisID EQ 2) THEN 	$
			GKVsd_GetGrid, self.Grid2,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid2.mnemonic) THEN 	$
			GKVsd_GetGrid, self.Grid2,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self -> GKVs1D::GET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
RETURN
END ; ****** GKVs2D::GET ****** ;


PRO GKVs2D::SET,axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices,		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
;
; Get values of elements of the realization of the GRID class, Grid2;
; and then call GKVs1D::GET to get values of elements of the GKVs1D Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs2D,
; or any or its subclasses
;
arg = self.Grid2
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 2) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	2:	IF(axisID EQ 2) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	3:	IF(axisID EQ 2) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid2.mnemonic) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self.Grid2 = arg
self -> GKVs1D::SET, 	axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 			$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
RETURN
END ; ****** GKVs2D::SET ****** ;


PRO GKVs2D::Info
;
; Prints information about contents of GKVs2D objects
;
self -> GKVs1D::Info
PRINT, 'Grid2:'
GKVsd_PrintGrid, self.Grid2
RETURN
END ; ***** GKVs2D::Info ***** ;


Function GKVs2D::INIT, signal
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
; Set tags for GKVs1D class
;
ok = self -> GKVs1D::INIT(signal)
;
; Set Tags for Grid2
;
NumTags=N_TAGS(self.Grid2)
FOR itag=0, NumTags-1 DO $
	self.Grid2.(itag) = signal.Grid2.(itag)

IF(PTR_VALID(self.grid2.values)) THEN	 BEGIN
	;
	; Check for uniform grid
	;
	self.Grid2.uniform = GKVsd_UniformGrid(self.Grid2.values)
	;
	; Check if Grid2.range is set ... and set if necessary
	;
	IF((self.Grid2.range[0] EQ 0.) AND (self.Grid2.range[1] EQ 0.)) THEN $
		self.Grid2.range = [MIN(*self.Grid2.Values, Max=max), max]
	;
	; Check if Grid2.irange is set ... and set if necessary
	;
	IF((self.Grid2.irange[0] EQ 0) AND (self.Grid2.irange[1] EQ 0)) THEN $
		self.Grid2.irange = [0, N_ELEMENTS(*self.Grid2.values)-1]
ENDIF
;
; Return on successful completion
;
RETURN, ok
END ; ***** GKVs2D::INIT ***** ;


PRO GKVs2D::CleanUp

PTR_FREE, self.Grid2.values
self -> GKVs1D::CleanUp

RETURN

END ; ***** GKVs2D::CleanUp ***** ;	

PRO GKVs2D__Define
struct = {	GKVs2D,				$	; "GK Visualization signal (2-D)"
		INHERITS GKVs1D,		$	; GKVs2D is a subclass of GKVs1D
		GRID2:{GRID}			}	; Include a second 'Grid' class
	
END ; ***** GKVs2D__Define ***** ;
