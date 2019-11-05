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
; Define Base Class GKVsd (Gyrokinetic Visualization signal data)
; 
; This class contains information about the code run and the data
; to be analyzed.
;
; Written by W.M. Nevins
;	1/15/00.  
;	Modified 1/29/00 to include error bars
;

PRO GKVsd::debug
;
; Purpose:
;
;		This proceedure intentionally creates an error, causing
;		IDL to halt with in this proceedure.  Since 'debug' is
;		a method for all GKVsd objects, you are now able to examine
;		the contents of 'self' using IDL's variable window, 'print',
;		etc.  

y=x		; An intentional error (x is undefined).
END ; ****** GKVsd::debug ****** ;

FUNCTION GKVsd_MIN, argg, MAX=mymax
;
; Local version of MIN function which returns correct scales 
; when argument is complex
;
arg=argg
tryagain:
info = SIZE(arg)
indx = info[0]+1
argType = info[indx]
IF (argType EQ 10) THEN BEGIN			; Argument is a pointer, 
	arg = *arg				; 	so dereferece, and
	GOTO, tryagain				;	try again.
ENDIF
IF ( (argType NE 6) AND (argType NE 9) ) THEN BEGIN
	rmin = MIN(arg, MAX=mymax, /NAN)
	RETURN, rmin
ENDIF
;
; argument is complex
;
IF(info[0] EQ 0) THEN BEGIN	; argument is complex scalar
	rarg = FLOAT(arg)
	rmin = rarg
	rmax = rarg
	iarg = IMAGINARY(arg)
	imin = iarg
	imax = iarg
	mymin = rmin < imin
	mymax = rmax > imax
	RETURN, mymin
ENDIF
; 
; Argument is a complex array.  
; 
arg = REFORM(arg)
info1 = SIZE(arg)
ndims = info1[0]
rarg = FLOAT(arg)
rmin = MIN(rarg, MAX=mymax, /NAN)
rmax=mymax
iarg = IMAGINARY(arg)
imin = MIN(iarg, MAX=imax, /NAN)
mymin = rmin < imin
mymax = rmax > imax
;
; for complex arrays with dimension greater than one
; return maximum absolute value (rather than max over
; real and imaginary parts).
;
IF(ndims GT 1) THEN mymax = MAX(ABS(arg))
RETURN, mymin
END ; ***** GKVsd_MIN ***** ;

FUNCTION GKVsd_UniformGrid, argg
; 
; Tests if argument is a uniform grid
; (or a pointer to a uniform grid)
;
arg=argg
tryagain:
info = SIZE(arg)
indx = info[0]+1
argType = info(indx)
IF (argType EQ 10) THEN BEGIN			; Argument is a pointer, 
	IF(PTR_VALID(arg)) THEN BEGIN		
		arg = *arg					; 	Dereference, and
		GOTO, tryagain				;	try again.
	ENDIF ELSE BEGIN
		ok = Error_Message("Invalid pointer", /Trace)
		RETURN, 0
	ENDELSE
ENDIF
IF (argType EQ 8) THEN BEGIN			; Argument is a structure
	IF( TAG_NAMES(arg, /STRUCTURE) NE 'GRID' ) THEN return, 0	; return if not a GRID structure
	irange = arg.irange
	imin=irange[0]
	imax=irange[1]
	arg = (*arg.values)[imin:imax]
	GOTO, tryagain
ENDIF
IF (info[0] NE 1) THEN BEGIN			; Argument is not a 1-D array
	ok = Error_Message("Bad Grid",/Trace)	; 	(or pointer to same)
	RETURN, 0
ENDIF
n = info[1]					; Number of elements in Grid
IF(n EQ 1) THEN RETURN, 1			; single-element grids are taken to be uniform
dx = arg[1:n-1] - arg[0:n-2]			; An array of grid spacings...
dmin = MIN(dx, MAX=dmax)				; Find Min and Max grid spacings
IF( dmin*dmax LE 0. ) THEN RETURN, 0		; Spacing changes sign (or vanishes...)
eps = (dmax - dmin)/dmax
IF ( eps LT 0. ) THEN eps = - eps
IF ( eps LT 1.e-2 ) THEN RETURN, 1		; 1 part in 1000 is good enough.
RETURN, 0
END ; ***** GKVsd_UniformGrid ***** ;


FUNCTION GKVsd::NumDims, ndims=ndims
Class = OBJ_CLASS(self)					; Get class of 'self'
Class = STRUPCASE(Class)					; Make sure it's uppercase
CASE Class OF							; Get number of dimensions
	"GKVSD"  :	ndims = 0
	"GKVS1D" :	ndims = 1
	"GKVS2D" :	ndims = 2
	"GKVS3D" :	ndims = 3
	"GKVS4D" :	ndims = 4
ELSE : BEGIN							; self isn't a GKV object???
	messageTxt = 'Dimensionality undefined for ' + Class + ' Class'
	MESSAGE, messageTxt, /Informational
	RETURN, -1
ENDELSE
ENDCASE
RETURN, ndims
END ; ****** GKVsd::NumDims ****** ;
	


FUNCTION GKVsd::MakeCopy, NoValues=novalues, NoErrorBars=noerrorbars, _Extra=extra
;
; Make "deep" copy of self
;
class = OBJ_CLASS(self)				; Remember POLYMORPHISM! 
								;	(self may be a subclass of GKVsd)
ok = EXECUTE('copy = { ' + class + ' }')	; Create 'copy', an instance of 
								;	the 'class' class-structure	
copy.mnemonic	= self.mnemonic		
copy.Title		= self.Title			; Copy fields from "self" into "copy"
Indices		= *self.indices			; Get Indices
copy.indices	= PTR_NEW(Indices)		; 	and make a new pointer for 'copy'
copy.units		= self. units
IF ( (NOT KEYWORD_SET(novalues)) AND PTR_VALID(self.values) ) THEN BEGIN	
								; Keywords NoValues and NoErrorBars
								;	suppress inclusion of value and
	values = *self.values				;	ErrorBar fields (which can take
	copy.values = PTR_NEW(values)		;	substantial memory) in the
ENDIF								;	"deep" copy.
copy.vrange	= self.vrange
IF ( (NOT KEYWORD_SET(noErrorBars)) AND PTR_VALID(self.ErrorBars) ) THEN BEGIN
	ErrorBars = *self.errorBars
	copy.ErrorBars = PTR_NEW(ErrorBars)
ENDIF	
copy.CodeName	= self.codename			
copy.CodePI	= self.CodePI
copy.RunID		= self.RunID
copy.FileID	= self.FileID

exestring = "result=OBJ_NEW('" + class + "', copy)"
ok = EXECUTE(exestring)				; Register object
RETURN, result
END ; ***** GKVsd::MakeCopy ***** ;


PRO GKVSD::Trash
OBJ_DESTROY, self
;heap_gc, /verbose
END ; ****** GKVSD::Trash ****** ;


FUNCTION GKVsd::Plus, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  addition
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0		; Undefined argument
arg=argg					; Make proxy for argument
try_again:
argInfo = SIZE(arg)				; Check proxy type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]
;
; First eliminate illegal argument types
;
CASE argType OF

0:	BEGIN							; Argument is undefined
		MESSAGE, "GKVsd::Plus:  Argument is undefined"
		RETURN, 0
	END
7:	BEGIN							; Argument is a STRING
		MESSAGE, "GKVsd::Plus:  Strings are not valid arguments"
		RETURN, 0
	END
8:	BEGIN							; Argument is a STRING
		MESSAGE, "GKVsd::Plus:  Structures are not valid arguments"
		RETURN, 0
	END
10:	BEGIN							; Argument is a POINTER
		arg = *arg					; Dereference pointer
		GOTO, try_again				;	and try again.
	END
11:	BEGIN							; Argument is an OBJECT reference
		IF (OBJ_ISA(arg, 'GKVsd') NE 1) THEN BEGIN
			MESSAGE, "GKVsd::Plus:  Argument is not a valid GKVsd object"
			RETURN, 0
		ENDIF						; Check units
		IF (self.units NE arg.units) THEN BEGIN
			MESSAGE, "GKVsd::Plus:  Incompatible units"
			RETURN, 0
		ENDIF
	; 
	; Have valid GKVsd object as argument, so proceed
	;
		result = self -> MakeCopy(/NoValues, /NoErrorBars)
		selfValues = *self.values
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic=STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.title = self.title + " + " + arg.title
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOF(units) EQ 7) THEN result.units=units
		result.values = PTR_NEW(selfValues + *arg.values)
		vmin = GKVsd_MIN(*result.values, MAX=vmax)
		result.vrange = [vmin, vmax]
		selfErrors = 0.
		IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
		argErrors = 0.
		IF PTR_VALID(arg.ErrorBars) THEN argErrors = *arg.ErrorBars
		ErrorBars = SQRT(selfErrors^2 + argErrors^2)
		IF (ErrorBars NE 0.) THEN result.ErrorBars = PTR_NEW(ErrorBars)
		RETURN, result
	END
ELSE: 
ENDCASE
;
; argument is a number or an array (Remember POLYMORPHISM!) 
;
result = self -> MakeCopy(/NoValues)
selfValues = *self.values
selfInfo = SIZE(selfValues)
selfDims = selfInfo[0]
CASE argDims OF
0	:	BEGIN					; Argument is a scalar
			values = selfValues + arg
		END
selfDims:	BEGIN					; Argument is a 1-D array
			valueString = "selfValues = (*self.values)["
			errorString = "    errors = (*self.ErrorBars)["
			FOR ndim = 1, selfDims DO BEGIN
				next="*"
				IF (selfInfo[ndim] NE argInfo[ndim]) THEN BEGIN
					dimnum = STRING(FORMAT='(I1)',ndim)
					iminstring = 'imin=self.Grid' + dimnum + '.imin'
					imaxstring = 'imin=self.Grid' + dimnum + '.imax'
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
					IF ( arginfo[ndim] NE (imax-imin+1) ) THEN BEGIN
						MESSAGE, "GKVsd::Plus:  Incompatible dimensions for array argument"
						OBJ_DESTROY, result
						RETURN, 0
					ENDIF
					next=STRING(imin) + ":" + STRING(imax)
					imax=imax-imin
					imin=0
					iminstring = 'result.Grid' + dimnum + '.imin =' + STRING(imin)
					imaxstring = 'result.Grid' + dimnum + '.imax =' + STRING(imax)
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
				ENDIF
				IF (ndim LT selfDims) THEN next=next + ","
				valueString = valueString + next
				errorString = errorString + next

			ENDFOR
			valueString = valueString + ']'
			errorString = errorString + ']'
			ok = EXECUTE(valueString)
			values = selfValues + arg
			IF PTR_VALID(self.ErrorBars) THEN BEGIN
				ok = EXECUTE(errorString)
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
ELSE:	  	BEGIN
			MESSAGE, "GKVsd::Plus:  Incompatible dimensions for array argument"
			OBJ_DESTROY, result
			RETURN, 0
		END
ENDCASE
result.values = PTR_NEW(values)
vmin = GKVsd_MIN(values, MAX=vmax)
result.vrange = [vmin, vmax]
IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic	= STRCOMPRESS(mnemonic, /REMOVE_ALL)
IF(TypeOf(title)    EQ 7) THEN result.Title 	= title
IF(TypeOf(units)    EQ 7) THEN result.units	= units
RETURN, result

END ; ***** GKVsd::Plus ***** :


FUNCTION GKVsd::Minus, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  subtraction
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0				; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Check proxy type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]
;
; First eliminate illegal argument types
;
CASE argType OF

0:	BEGIN								; Argument is undefined
		MESSAGE, "GKVsd::Minus:  Argument is undefined"
		RETURN, 0
	END
7:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Minus:  Strings are not valid arguments"
		RETURN, 0
	END
8:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Minus:  Structures are not valid arguments"
		RETURN, 0
	END
10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an OBJECT reference
		IF (OBJ_ISA(arg, 'GKVsd') NE 1) THEN BEGIN
			MESSAGE, "GKVsd::Minus:  Argument is not a valid GKVsd object"
			RETURN, 0
		ENDIF							; Check units
		IF (self.units NE arg.units) THEN BEGIN
			MESSAGE, "GKVsd::Minus:  Incompatible units"
			RETURN, 0
		ENDIF
	; 
	; Have valid GKVsd object as argument, so proceed
	;
		result = self -> MakeCopy(/NoValues, /NoErrorBars)
		selfValues = *self.values
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.title = self.title + " - " + arg.title
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOF(units) EQ 7) THEN result.units=units
		result.values = PTR_NEW(selfValues - *arg.values)
		vmin = GKVsd_MIN(*result.values, MAX=vmax)
		result.vrange = [vmin, vmax]
		selfErrors = 0.
		IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
		argErrors = 0.
		IF PTR_VALID(arg.ErrorBars) THEN argErrors = *arg.ErrorBars
		ErrorBars = SQRT(selfErrors^2 + argErrors^2)
		IF (ErrorBars NE 0.) THEN result.ErrorBars = PTR_NEW(ErrorBars)
		RETURN, result
	END
ELSE: 
ENDCASE
;
; argument is a number or an array (Remember POLYMORPHISM!) 
;
result = self -> MakeCopy(/NoValues)
selfValues = *self.values
selfInfo = SIZE(selfValues)
selfDims = selfInfo[0]
CASE argDims OF
0	:	BEGIN					; Argument is a scalar
			values = selfValues - arg
		END
selfDims:	BEGIN					; Argument is a 1-D array
			valueString = "selfValues = (*self.values)["
			errorString = "    errors = (*self.ErrorBars)["
			FOR ndim = 1, selfDims DO BEGIN
				next="*"
				IF (selfInfo[ndim] NE argInfo[ndim]) THEN BEGIN
					dimnum = STRING(FORMAT='(I1)',ndim)
					iminstring = 'imin=self.Grid' + dimnum + '.imin'
					imaxstring = 'imin=self.Grid' + dimnum + '.imax'
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
					IF ( arginfo[ndim] NE (imax-imin+1) ) THEN BEGIN
						MESSAGE, "GKVsd::Minus:  Incompatible dimensions for array argument"
						OBJ_DESTROY, result
						RETURN, 0
					ENDIF
					next=STRING(imin) + ":" + STRING(imax)
					imax=imax-imin
					imin=0
					iminstring = 'result.Grid' + dimnum + '.imin =' + STRING(imin)
					imaxstring = 'result.Grid' + dimnum + '.imax =' + STRING(imax)
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
				ENDIF
				IF (ndim LT selfDims) THEN next=next + ","
				valueString = valueString + next
				errorString = errorString + next

			ENDFOR
			valueString = valueString + ']'
			errorString = errorString + ']'
			ok = EXECUTE(valueString)
			values = selfValues - arg
			IF PTR_VALID(self.ErrorBars) THEN BEGIN
				ok = EXECUTE(errorString)
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
ELSE:	  	BEGIN
			MESSAGE, "GKVsd::Minus:  Incompatible dimensions for array argument"
			OBJ_DESTROY, result
			RETURN, 0
		END
ENDCASE
result.values = PTR_NEW(values)
vmin = GKVsd_MIN(values, MAX=vmax)
result.vrange = [vmin, vmax]
IF(TypeOF(mnemonic) EQ 7) THEN result.mnemonic 	= STRCOMPRESS(mnemonic, /REMOVE_ALL)
IF(TypeOf(title)    EQ 7) THEN result.Title 	= title
IF(TypeOf(units)    EQ 7) THEN result.units	= units
RETURN, result

END ; ***** GKVsd::Minus ***** :


FUNCTION GKVsd::Times, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  Multiplication
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Check proxy type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]
;
; First eliminate illegal argument types
;
CASE argType OF

0:	BEGIN								; Argument is undefined
		MESSAGE, "GKVsd::Times:  Argument is undefined"
		RETURN, 0
	END
7:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Times:  Strings are not valid arguments"
		RETURN, 0
	END
8:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Times:  Structures are not valid arguments"
		RETURN, 0
	END
10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an OBJECT reference
		IF (OBJ_ISA(arg, 'GKVsd') NE 1) THEN BEGIN
			MESSAGE, "GKVsd::Times:  Argument is not a valid GKVsd object"
			RETURN, 0
		ENDIF	
	; 
	; Have valid GKVsd object as argument, so proceed
	;
		result = self -> MakeCopy(/NoValues, /NoErrorBars)
		selfValues = *self.values
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.title = self.title + " * " + arg.title
		IF(TypeOf(title) EQ 7) THEN result.title=title
		result.units = "(" + self.units + ")*(" + arg.units + ")"
		IF(TypeOf(units) EQ 7) THEN result.units=units
		result.values = PTR_NEW( selfValues * (*arg.values) )
		vmin = GKVsd_MIN(*result.values, MAX=vmax)
		result.vrange = [vmin, vmax]
		selfErrors = 0.
		IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
		argErrors = 0.
		IF PTR_VALID(arg.ErrorBars) THEN argErrors = *arg.ErrorBars
		ErrorBars = SQRT(selfErrors^2 + argErrors^2)
		IF (ErrorBars NE 0.) THEN result.ErrorBars = PTR_NEW(ErrorBars)
		RETURN, result
	END
ELSE: 
ENDCASE
;
; argument is a number or an array (Remember POLYMORPHISM!) 
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
selfValues = *self.values
selfInfo = SIZE(selfValues)
selfDims = selfInfo[0]
CASE argDims OF
0	:	BEGIN					; Argument is a scalar
			values = selfValues*arg
			IF PTR_VALID(self.errorBars) THEN BEGIN
				errors = *self.ErrorBars*arg
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
selfDims:	BEGIN					; Argument is a 1-D array
			valueString = "selfValues = (*self.values)["
			errorString = "    errors = (*self.ErrorBars)["
			FOR ndim = 1, selfDims DO BEGIN
				next="*"
				IF (selfInfo[ndim] NE argInfo[ndim]) THEN BEGIN
					dimnum = STRING(FORMAT='(I1)',ndim)
					iminstring = 'imin=self.Grid' + dimnum + '.imin'
					imaxstring = 'imin=self.Grid' + dimnum + '.imax'
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
					IF ( arginfo[ndim] NE (imax-imin+1) ) THEN BEGIN
						MESSAGE, "GKVsd::Times:  Incompatible dimensions for array argument"
						OBJ_DESTROY, result
						RETURN, 0
					ENDIF
					next=STRING(imin) + ":" + STRING(imax)
					imax=imax-imin
					imin=0
					iminstring = 'result.Grid' + dimnum + '.imin =' + STRING(imin)
					imaxstring = 'result.Grid' + dimnum + '.imax =' + STRING(imax)
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
				ENDIF
				IF (ndim LT selfDims) THEN next=next + ","
				valueString = valueString + next
				errorString = errorString + next

			ENDFOR
			valueString = valueString + ']'
			errorString = errorString + ']'
			ok = EXECUTE(valueString)
			values = selfValues * arg
			IF PTR_VALID(self.ErrorBars) THEN BEGIN
				ok = EXECUTE(errorString)
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
ELSE:	  	BEGIN
			MESSAGE, "GKVsd::Times:  Incompatible dimensions for array argument"
			OBJ_DESTROY, result
			RETURN, 0
		END
ENDCASE
result.values = PTR_NEW(values)
vmin = GKVsd_MIN(values, MAX=vmax)
result.vrange = [vmin, vmax]
IF(TypeOF(mnemonic) EQ 7) THEN result.mnemonic 	= STRCOMPRESS(mnemonic, /REMOVE_ALL)
IF(TypeOf(title)    EQ 7) THEN result.Title 	= title
result.units = ''
IF(TypeOf(units)    EQ 7) THEN result.units	= units
RETURN, result

END ; ***** GKVsd::Times ***** :


FUNCTION GKVsd::Over, argg, title=title, units=units, mnemonic=mnemonic
;
; Define basic arithmetic operation:  Division
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0			; Undefined argument
arg=argg								; Make proxy for argument
try_again:
argInfo = SIZE(arg)						; Check proxy type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]
;
; First eliminate illegal argument types
;
CASE argType OF

0:	BEGIN								; Argument is undefined
		MESSAGE, "GKVsd::Over:  Argument is undefined"
		RETURN, 0
	END
7:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Over:  Strings are not valid arguments"
		RETURN, 0
	END
8:	BEGIN								; Argument is a STRING
		MESSAGE, "GKVsd::Over:  Structures are not valid arguments"
		RETURN, 0
	END
10:	BEGIN								; Argument is a POINTER
		arg = *arg						; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN								; Argument is an OBJECT reference
		IF (OBJ_ISA(arg, 'GKVsd') NE 1) THEN BEGIN
			MESSAGE, "GKVsd::Over:  Argument is not a valid GKVs1D object"
			RETURN, 0
		ENDIF	
	; 
	; Have valid GKVsd object as argument, so proceed
	;
		result = self -> MakeCopy(/NoValues, /NoErrorBars)
		selfValues = *self.values
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.title = self.title + " / " + arg.title
		IF(TypeOf(title) EQ 7) THEN result.title=title
		result.units = "(" + self.units + ")/(" + arg.units + ")"
		IF(TypeOf(units) EQ 7) THEN result.units=units
		argValues = *arg.values
		result.values = PTR_NEW(selfValues / argValues)
		vmin = GKVsd_MIN(*result.values, MAX=vmax)
		result.vrange = [vmin, vmax]
		selfErrors = 0.
		IF PTR_VALID(self.ErrorBars) THEN selfErrors = *self.ErrorBars
		argErrors = 0.
		IF PTR_VALID(arg.ErrorBars) THEN argErrors = *arg.ErrorBars
		errorSq = (selfErrors/argValues)^2 + (selfValues*argErrors/argValues)^2
		errorBars = SQRT(errorSq)
		IF (MAX(ErrorBars) NE 0.) THEN result.ErrorBars = PTR_NEW(ErrorBars)
		RETURN, result
	END
ELSE: 
ENDCASE
;
; argument is a number or an array (Remember POLYMORPHISM!) 
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
selfValues = *self.values
selfInfo = SIZE(selfValues)
selfDims = selfInfo[0]
CASE argDims OF
0	:	BEGIN					; Argument is a scalar
			values = selfValues/arg
			IF PTR_VALID(self.errorBars) THEN BEGIN
				errors = *self.ErrorBars/arg
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
selfDims:	BEGIN					; Argument is a 1-D array
			valueString = "selfValues = (*self.values)["
			errorString = "    errors = (*self.ErrorBars)["
			FOR ndim = 1, selfDims DO BEGIN
				next="*"
				IF (selfInfo[ndim] NE argInfo[ndim]) THEN BEGIN
					dimnum = STRING(FORMAT='(I1)',ndim)
					iminstring = 'imin=self.Grid' + dimnum + '.imin'
					imaxstring = 'imin=self.Grid' + dimnum + '.imax'
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
					IF ( arginfo[ndim] NE (imax-imin+1) ) THEN BEGIN
						MESSAGE, "GKVsd::Over:  Incompatible dimensions for array argument"
						OBJ_DESTROY, result
						RETURN, 0
					ENDIF
					next=STRING(imin) + ":" + STRING(imax)
					imax=imax-imin
					imin=0
					iminstring = 'result.Grid' + dimnum + '.imin =' + STRING(imin)
					imaxstring = 'result.Grid' + dimnum + '.imax =' + STRING(imax)
					ok = EXECUTE(iminstring)
					ok = EXECUTE(imaxstring)
				ENDIF
				IF (ndim LT selfDims) THEN next=next + ","
				valueString = valueString + next
				errorString = errorString + next

			ENDFOR
			valueString = valueString + ']'
			errorString = errorString + ']'
			ok = EXECUTE(valueString)
			values = selfValues / arg
			IF PTR_VALID(self.ErrorBars) THEN BEGIN
				ok = EXECUTE(errorString)
				result.ErrorBars = PTR_NEW(errors)
			ENDIF
		END
ELSE:	  	BEGIN
			MESSAGE, "GKVsd::Over:  Incompatible dimensions for array argument"
			OBJ_DESTROY, result
			RETURN, 0
		END
ENDCASE
result.values = PTR_NEW(values)
vmin = GKVsd_MIN(values, MAX=vmax)
result.vrange = [vmin, vmax]
IF(TypeOF(mnemonic) EQ 7) THEN result.mnemonic 	= STRCOMPRESS(mnemonic, /REMOVE_ALL)
IF(TypeOf(title)    EQ 7) THEN result.Title 	= title
result.units = ''
IF(TypeOf(units)    EQ 7) THEN result.units	= units
RETURN, result

END ; ***** GKVsd::Over ***** :


FUNCTION GKVsd::IndexString, plotAxis, Pretty=pretty
;
; Returns a string array containing index elements for a GKVsd object.
;
; The argument 'plotAxis' is a 1-D array whose length as equal to the 
; dimensionality of the GKVsd object (GKVs1D -> 1, GKVs2D -> 2, etc).
; If an element of 'plotAxis' is positive, the value corresponds to the
; number of one of the axis being plotted:
;
;	1 for x (in the sense of IDL... that is, horizontal), 
;	2 for y (in the sense of IDL... that is, vertical), 
;	etc.
;
; Negative elements of 'plotAxis' correspond to axis that are being 'sliced' in
; this plot, and the absolute value of this elements is the index at which the 
; data is being sliced.
;
; Zero values of 'plotAxis' signify that the corresponding element of 
; 'Indices' should be left alone.
;
indices = *self.indices				; Get GKVsd object's Indices
info = SIZE(indices)
Nindices = info[1]					; Find length of 'indices' array
IF(info[0] EQ 0) THEN nIndices=1
;
; Loop to replace *'s in 'indices' array with corresponding axis title
; and (for negative values of 'plotAxis') grid values.
;
iaxis = 0							; index into 'plotAxis' array (elements 1 thru ndims)
FOR i=0, nIndices-1 DO BEGIN			; i indexes 'Indices' array, 
	IF(indices[i] EQ '*') THEN BEGIN	; which can have more than ndims elements
		iaxis=iaxis+1
		IF(plotAxis[iaxis] EQ 0) THEN GOTO, DONE
		axis_str = STRING(iaxis, FORMAT='(i1)')
		command_str = 'grid = self.grid' + axis_str
		ok = EXECUTE(command_str)
		IF(plotAxis[iaxis] GT 0) THEN BEGIN
			indices[i] = grid.mnemonic
			IF KEYWORD_SET(pretty) THEN indices[i] = grid.title
		ENDIF ELSE BEGIN
			Value = (*grid.values)[-plotAxis[iaxis]]
			valueStr = STRING(value, FORMAT='(G10.3)')
			indexStr = grid.mnemonic + '=' + valueStr
			IF KEYWORD_SET(pretty) THEN indexStr = grid.title + '=' + valueStr
			indexStr = STRCOMPRESS(indexStr, /REMOVE_ALL)
			indices[i] = indexStr
		ENDELSE
	ENDIF
	DONE	:	; Don't do anything
ENDFOR
RETURN, indices
END ; ****** GKVsd::IndexString ****** ;


FUNCTION GKVsd::IndexRemove, axisNum
;
; returns string array of 'indices' with element 
; corresponding to 'axisNum', or (if 'axisNum' is an array) 
; each element of 'axisNum', removed.
;
; Written by W.M. Nevins
;	4/6/00
;
; Modified by W.M. Nevibs
;	9/25/00
; to allow for multiple elements to axisNum
;
indices = *self.indices
info=SIZE(indices)
Nindices=info[1]
nAxes = N_ELEMENTS(axisNum)
IF((Nindices - nAxes) LE 0) THEN RETURN, ''
result = STRARR(Nindices-nAxes)
;
iaxis=0
j=0
FOR i=0, Nindices-1 DO BEGIN
	IF(indices[i] EQ "*") THEN BEGIN
		iaxis=iaxis+1
		FOR k=0,nAxes-1 DO IF(iaxis EQ axisNum[k]) THEN GOTO, DONEK
		result[j] = indices[i]
		j=j+1
		DONEK:
	ENDIF ELSE BEGIN
		result[j] = indices[i]
		j=j+1
	ENDELSE
ENDFOR
RETURN, result
END ; ****** GKVsd::IndexRemove ****** ;

FUNCTION GKVsd::IndexMax, axisNum
;
; Returns string array of 'indices' with element
; corresponding to 'axisNum' set equal to "'Grid.title' = '=Max'"
;
; Written by W.M. Nevins
;	4/12/00
;
indices = *self.indices
info=SIZE(indices)
Nindices=info[1]
;
iaxis=0
FOR i=0, Nindices-1 DO BEGIN
	IF(indices[i] EQ "*") THEN BEGIN
		iaxis=iaxis+1
		IF(iaxis EQ axisNum) THEN BEGIN
			axis_str = STRING(iaxis, FORMAT='(i1)')
			command_str = 'grid = self.grid' + axis_str
			ok = EXECUTE(command_str)
			indexStr = grid.title + '=Max'
			indexStr = STRCOMPRESS(indexStr, /REMOVE_ALL)
			indices[i] = indexStr
			GOTO, DONE
		ENDIF
	ENDIF
ENDFOR
DONE	:  ; finished...
RETURN, indices
END ; ****** GKVsd::IndexMax ****** ;



Pro GKVsd::Save, GKVsdOjb, FileName=file_name, Path=path, DeBug=d
;
; Purpose:
;
; 		This proceedure SAVEs a GKVsd object to a disk
;		for later retreval with GKV_Restore (see below).
;
; KEYWORDS:
;
;	FileName		Name of file.  Must include either full path, or else the path from the 
;				current working directory.  If FileName is not provided, DIALOG_PICKFILE will 
;				be called to allow the user to select an appropriate file name.
;
;	Path			Sets path to initial directory selected by DIALOG_PICKFILE (which will allow user
;				to change to other directories...).  If no path is provided, DIALOG_PICKFILE will
;				default to current working directory. 
; 
IF(KEYWORD_SET(file_name)) THEN BEGIN
	subStrings=STRSPLIT(file_name, '.', /EXTRACT)
	nSubStrings = N_ELEMENTS(subStrings)
	extension = subStrings[nSubStrings-1]
	filename=file_name
	IF(extension NE '.gkv') THEN filename=file_name + '.gkv'
ENDIF
IF(N_ELEMENTS(fileName) EQ 0) THEN BEGIN
	filename = self.mnemonic + '.gkv'
	saveFile = DIALOG_PICKFILE(PATH=path, FILE=filename, GET_PATH=newPath, /write)
	CD, newPath 
ENDIF ELSE BEGIN
	saveFile = fileName
ENDELSE

GKVsdClass = OBJ_CLASS(self)
;
SAVE, GKVsdClass, self, FIleName=saveFile, /COMPRESS, VERBOSE=d
RETURN
END ; ****** GKVsd::Save ****** ;


FUNCTION GKV_RESTORE, FileName=fileName, PATH=path, DeBug=D
;
; Proceedure to Restore the GKVsd object previously saved 
; (using GKVsd::Save) in the file 'file_name'
;
IF(NOT KEYWORD_SET(filename)) THEN BEGIN
	fileName = DIALOG_PICKFILE(Path=path, Filter='*.gkv', GET_PATH=newPath)
	CD, newPath
ENDIF ELSE BEGIN
	ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
	IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... so try  DIALOG_PICKFILE
		fileName = DIALOG_PICKFILE(Path=path, Filter='*.gkv', GET_PATH=newPath)
		CD, newPath
	ENDIF		
ENDELSE
;
IF(filename[0] EQ '') THEN RETURN, 0
RESTORE, filename, /RELAXED_STRUCTURE_ASSIGNMENT, VERBOSE=d
result = OBJ_NEW(GKVsdClass, self)
RETURN, result
END ; ****** GKVsd_RESTORE ****** ;


PRO GKVsd::Info
;
; Prints information about GKV Data objects
;
GKVsClass = OBJ_CLASS(self)
GKVSCLASS = STRUPCASE(GKVsClass)
IF(PTR_VALID(self.values)) THEN BEGIN
	ValInfo = SIZE(*self.values)
ENDIF ELSE BEGIN
	ValInfo = "Null Pointer"
ENDELSE
ValueText = "ValueInfo"
ErrorInfo = "no error bars"
IF(PTR_VALID(self.ErrorBars)) THEN	$
	ErrorInfo = SIZE(*self.ErrorBars)
ErrorText = "ErrorInfo"
IF (GKVSCLASS EQ 'GKVSD') THEN BEGIN
	valInfo   = *self.values
	ValueText ="   values"
	IF(PTR_VALID(self.ErrorBars)) THEN BEGIN
		ErrorInfo = *self.ErrorBars
		ErrorText = "ErrorBars"
	ENDIF
ENDIF
Print, GKVsClass
PRINT, " mnemonic", " = ", self.mnemonic
PRINT, "    Title", " = ", self.Title
IF(PTR_VALID(self.Indices)) THEN BEGIN
	PRINT, "  Indices", " = [", STRJOIN(*self.Indices, ", "), "]"
ENDIF ELSE BEGIN
	PRINT, "  Indices = Null Pointer"
ENDELSE
PRINT, "    units", " = ", self.units
PRINT,  ValueText , " = ", ValInfo
PRINT, "   vrange", " = ", self.vrange
PRINT,  ErrorText , " = ", ErrorInfo
PRINT, " CodeName", " = ", self.CodeName
PRINT, "   CodePI", " = ", self.CodePI
PRINT, "    RunID", " = ", self.RunID
PRINT, "   FileID", " = ", self.FileID

RETURN
END ; ****** GKVsd::Info ****** ;


FUNCTION GKVsd::GetValues
;
; Returns value of a GKVsd object
;
;	Written by W.M. Nevins
;	6/20/00

value = *self.values
RETURN, value
END ; ****** GKVsd::GetValues ******


FUNCTION GKVsd::GetErrors
;
; Returns value of the errorbars for a GKVsd object
;
;	Written by W.M. Nevins
;	10/10/00

error = *self.ErrorBars
RETURN, error
END ; ****** GKVsd::GetErrors ******


PRO GKVsd::GET,	mnemonic=mnemonic, Title=Title, Indices=indices, $
			units=units,values=values, vrange=vrange,		$
			ErrorBars=ErrorBars, CodeName=CodeName,		$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Get values of elements of GKVsd Class
;
mnemonic	= self.mnemonic
Title		= self.Title
Indices	= self.Indices
units		= self.units
values	= self.values
vrange	= self.vrange
ErrorBars	= self.ErrorBars
CodeName	= self.CodeName
CodePI	= self.CodePI
RunID		= self.RunID
FileID	= self.FIleID
END ; ***** GKVsd::GET ***** ;

	
PRO GKVsd::SET,	mnemonic=mnemonic, Title=Title, Indices=indices,$
			units=units,values=values, vrange=vrange,		$
			ErrorBars=ErrorBars, CodeName=CodeName,		$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Set values of elements of GKVsd Class
;
type = TypeOf(mnemonic)
IF(type EQ 7) THEN	$
	self.mnemonic = mnemonic

type = TypeOf(Title)
IF(type EQ 7) THEN	$
	self.Title	= Title
	
type = TypeOf(Indices)
IF(type EQ 10) THEN BEGIN
	PTR_FREE, self.Indices
	self.Indices = Indices
ENDIF
IF(Type EQ 7) THEN BEGIN
	PTR_FREE, self.Indices
	self.Indices = PTR_NEW(Indices)
ENDIF

type = TypeOf(units)
IF(type EQ 7) THEN	$
	self.units	= units

type = TypeOf(values)
IF(type EQ 10) THEN BEGIN
	PTR_FREE, self.values
	self.values = values
	IF(TypeOf(vrange) EQ 0) THEN BEGIN
		vmin = GKVsd_Min(*values, Max=vmax)
	ENDIF ELSE BEGIN
		vmin=vrange[0]
		vmax=vrange[1]
	ENDELSE
	IF(vmax EQ vmin) THEN vmax = vmin + 1
	vrange = [vmin, vmax]
ENDIF
IF(((type GT 0) AND (type LT 7)) OR (type EQ 9)) THEN BEGIN
	PTR_FREE, self.values
	self.values = PTR_NEW(values)
	IF(TypeOf(vmin) EQ 0) THEN vmin = GKVsd_Min(values, Max=vmax)
	IF(vmax EQ vmin) THEN vmax = vmin + 1
	vrange = [vmin, vmax]
ENDIF

type = TypeOf(vrange)
IF(((type GT 0) AND (type LT 7)) OR (type EQ 9)) THEN			$
	 self.vrange = vrange

type = TypeOf(ErrorBars)
IF(type EQ 10) THEN BEGIN
	PTR_FREE, self.ErrorBars
	self.ErrorBars	= ErrorBars
ENDIF
IF(((type GT 0) AND (type LT 7)) OR (type EQ 9)) THEN BEGIN
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(ErrorBars)
ENDIF

type = TypeOf(CodeName)
IF(type EQ 7) THEN	$
	self.CodeName = CodeName

type = TypeOf(CodePI)
IF(type EQ 7) THEN	$
	self.CodePI = CodePI

type = TypeOf(RunID)
IF(type EQ 7) THEN	$
	self.RunID	= RunID

type = TypeOf(FileID)
IF(type EQ 7) THEN	$
	self.FIleID = FileID

END ; ***** GKVsd::SET ***** ;


PRO GKVsD::CleanUp

PTR_FREE, self.values
PTR_FREE, self.ErrorBars
PTR_FREE, self.indices

RETURN
END ; ***** GKVsD::CleanUp ***** ;


FUNCTION GKVsd::INIT, signal
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
; Setting self=signal DOESN'T work, so ...
;
; Set tags for GKVsd class
;
	self.mnemonic	= STRCOMPRESS(signal.mnemonic, /REMOVE_ALL)
	self.Title	= STRCOMPRESS(signal.Title,    /REMOVE_ALL)
	self.Indices	= signal.Indices
	self.units	= STRCOMPRESS(signal.units,    /REMOVE_ALL)
	self.values	= signal.values
	IF( (signal.vrange[0] EQ 0) AND (signal.vrange[1] EQ 0) ) THEN BEGIN
		IF(PTR_VALID(signal.values)) THEN BEGIN
			vmin = GKVsd_MIN(*signal.values, MAX=vmax)
		ENDIF ELSE BEGIN
			vmin = 0.
			vmax = 0.
		ENDELSE
		IF(vmax EQ vmin) THEN vmax = vmin + 1
		self.vrange	= [vmin, vmax]
	ENDIF ELSE BEGIN
		self.vrange = signal.vrange
	ENDELSE
	self.ErrorBars	= signal.ErrorBars
	self.CodeName	= signal.CodeName
	self.CodePI	= signal.CodePI
	self.RunID	= signal.RunID
	self.FileID	= signal.FileID
RETURN, 1
END

PRO GKVsd__Define
struct = {	GKVsd,				$	; Define "Gyrokinetic Visualization signal data"
							; 	class
		mnemonic	: "",		$	; Alpha/numeric mnemonic (as used in code)
		Title		: "",		$	; "Pretty" name	(use IDL Vector Font)
		Indices		: PTR_NEW(),	$	; Pointer to string array containing values of indices
		units		: "",		$	; Units		(use IDL Vector Font)
		values		: PTR_NEW(),	$	; pointer to value(s)
		vrange		: FLTARR(2),	$	; Range of values for displays
		ErrorBars	: PTR_NEW(),	$	; Pointer to error bar(s)
		CodeName	: "",		$	; Name of code which produced data
		CodePI		: "",		$	; Name of PI (or team) responsible for code runs
		RunID		: "",		$	; Run identifier
		FileID		: ""		}	; Name of file from which data was extracted
								
END ; ***** GKVsd__Define ***** ;

								
		
