FUNCTION GKVs2D::Rotate, arg, _EXTRA=extra
;
; Purpose:
;
;	This proceedure rotates the data in 'self' through the angle
;	specified by the argument 'arg'.
;
;
; Argument:
;
;	The argument to this proceedure is the desired (clockwise) angle of 
;	rotation in degrees
;
;
; Input Keywords:
;
;	x0		Value of first independent variable rotate about 
;			('center of rotation').  Defaults to 0. (Optional)
;
;	'mnemonic1'_0 	A synonym for x0 (see above), where 'mnemonic1' is the 
;			mnemonic for the first independent variable. (Optional)
;
;	y0		Value of second independent variable rotate about 
;			('center of rotation').  Defaults to 0. (Optional)
;
;	'mnemonic2'_0 	A synonym for y0 (see above), where 'mnemonic2' is the 
;			mnemonic for the first independent variable. (Optional)
;
;
;	title1		Set to an ascii string containing the title associated 
;			with the first independent variable.  Defaults to 'x!D1!N'.
;			(Optional)
;
;	mnemonic1	Set to an ascii string containing the mnemonic associated
;			with the first independentend variable.  Defaults to 'x1'.
;			(Optional)
;
;	units1		Set to an ascii string containing the units assciated with
;			the first independent variable.  Defaults to input value of
;			these units.  (Optional)
;
;	title2		Set to an ascii string containing the title associated 
;			with the second independent variable.  Defaults to 'x!D2!N'.
;			(Optional)
;
;	mnemonic2	Set to an ascii string containing the mnemonic associated
;			with the second independentend variable.  Defaults to 'x2'.
;			(Optional)
;
;	units2		Set to an ascii string containing the units assciated with
;			the second independent variable.  Defaults to input value of
;			these units.  (Optional)
;
;
;  Output KeyWords:	None
;
;
; Written by W.M. Nevins
;	6/6/01
;
; First make sure that the independent variables are defined on uniform grids
;
self -> ScaleAxis, 1, /Uniform
self -> ScaleAxis, 2, /Uniform
;
; Now, determine which grid has finest spacing
;
grid1Values = *self.grid1.values
grid2Values = *self.grid2.values
dx1 = grid1Values[1] - grid1Values[0]
dx2 = grid2Values[1] - grid2Values[0]
;
; Interpolate data to spacing of the finer grid
;
target = self -> MakeCopy(/noValues, /noErrorBars)
IF(dx1 GT dx2) THEN BEGIN
	target.grid1.values = self.grid2.values
	target.grid1.irange = self.grid2.irange
	target.grid1.range  = self.grid2.range
ENDIF
IF(dx2 GT dx1) THEN BEGIN
	target.grid2.values = self.grid1.values
	target.grid2.irange = self.grid1.irange
	target.grid2.range  = self.grid1.range
ENDIF
result = self -> Interpolate(target)
;
; Find indices to center of rotation
;
x0=0.
xx0 = GetKeyWord('x0', extra)
IF(TypeOf(xx0) NE 7) THEN x0=xx0
xx0 = GetKeyWord(result.grid1.mnemonic + '_0', extra)
IF(TypeOf(xx0) NE 7) ThEN x0=xx0
y0=0.
yy0 = GetKeyWord('y0', extra)
IF(TypeOf(yy0) NE 7) THEN y0=yy0
yy0 = GetKeyWord(result.grid2.mnemonic + '_0', extra)
IF(TypeOf(yy0) NE 7) ThEN y0=yy0

x1Values = *result.grid1.values
temp1 = (x1Values - x0)^2
eps1 = MIN(temp1, i)	; the point here is to get subscript of minimum ('i'), NOT to find the minimum ('eps1')

x2Values = *result.grid2.values
temp2 = (x2Values - y0)^2
eps2 = MIN(temp2, j)	; the point here is to get subscript of minimum ('j'), NOT to find the minimum ('eps2')
;
; Now, get values from 'result' and rotate them through the angle 'arg'
;
values = *result.values
newValues = ROT(values, arg, 1.,i, j, /INTERP, MISSING=0., /PIVOT) 
PTR_FREE, result.values
result.values = PTR_NEW(newValues)
IF PTR_VALID(result.ErrorBars) THEN BEGIN
	errorBars = *result.errorBars
	PTR_FREE, result.errorBars
	errorBars = ROT(errorBars, arg, x0=x0, y0=y0, /INTERP, MISSING=0.) 
	result.errorBars = PTR_NEW(errorBars)
ENDIF
;
; Update titles, mnemonics, etc.
;
title1='x!D1!N'
nTitle1 = GetKeyWord('title1', extra)
IF(TypeOf(ntitle1) EQ 7) THEN BEGIN
	IF(ntitle1 NE 'undefined') THEN title1=ntitle1
ENDIF 
result.grid1.title = title1

mnemonic1 = 'x1'
nMnemonic1 = GetKeyWord('mnemonic1', extra)
IF(TypeOf(nMnemonic1) EQ 7) THEN BEGIN
	IF(nMnemonic1 NE 'undefined') THEN mnemonic1=nMnemonic1
ENDIF
result.grid1.mnemonic = mnemonic1

nUnits1 = GetKeyWord('units1', extra)
IF(TypeOF(nUnits1) EQ 7) THEN BEGIN
	IF(nUnits1 NE 'undefined') THEN result.grid1.units=nUnits1
ENDIF

title2='x!D2!N'
nTitle2 = GetKeyWord('title2', extra)
IF(TypeOf(ntitle2) EQ 7) THEN BEGIN
	IF(ntitle2 NE 'undefined') THEN title2=ntitle2
ENDIF 
result.grid2.title = title2

mnemonic2 = 'x2'
nMnemonic2 = GetKeyWord('mnemonic2', extra)
IF(TypeOf(nMnemonic2) EQ 7) THEN BEGIN
	IF(nMnemonic2 NE 'undefined') THEN mnemonic2=nMnemonic2
ENDIF
result.grid2.mnemonic = mnemonic2

nUnits2 = GetKeyWord('units2', extra)
IF(TypeOF(nUnits2) EQ 7) THEN BEGIN
	IF(nUnits2 NE 'undefined') THEN result.grid2.units=nUnits2
ENDIF

RETURN, result

END ; ****** GKVs2D::Rotate ****** ;




FUNCTION GKVs3D::Rotate, arg, _extra=extra
;
; Purpose:
;
;	This proceedure rotates the data in 'self' through the angle
;	specified by the argument 'arg' about axis #3.
;
;
; Argument:
;
;	The argument to this proceedure is the desired (clockwise) angle of 
;	rotation in degrees
;
;
; Input Keywords:
;
;	x0		Value of first independent variable rotate about 
;			('center of rotation').  Defaults to 0. (Optional)
;
;	'mnemonic1'_0 	A synonym for x0 (see above), where 'mnemonic1' is the 
;			mnemonic for the first independent variable. (Optional)
;
;	y0		Value of second independent variable rotate about 
;			('center of rotation').  Defaults to 0. (Optional)
;
;	'mnemonic2'_0 	A synonym for y0 (see above), where 'mnemonic2' is the 
;			mnemonic for the first independent variable. (Optional)
;
;
;	title1		Set to an ascii string containing the title associated 
;			with the first independent variable.  Defaults to 'x!D1!N'.
;			(Optional)
;
;	mnemonic1	Set to an ascii string containing the mnemonic associated
;			with the first independentend variable.  Defaults to 'x1'.
;			(Optional)
;
;	units1		Set to an ascii string containing the units assciated with
;			the first independent variable.  Defaults to input value of
;			these units.  (Optional)
;
;	title2		Set to an ascii string containing the title associated 
;			with the second independent variable.  Defaults to 'x!D2!N'.
;			(Optional)
;
;	mnemonic2	Set to an ascii string containing the mnemonic associated
;			with the second independentend variable.  Defaults to 'x2'.
;			(Optional)
;
;	units2		Set to an ascii string containing the units assciated with
;			the second independent variable.  Defaults to input value of
;			these units.  (Optional)
;
;
;  Output KeyWords:	None
;
;
; Written by W.M. Nevins
;	6/6/01
;
; First make sure that the independent variables are defined on uniform grids
;
self -> ScaleAxis, 1, /Uniform
self -> ScaleAxis, 2, /Uniform
;
; Now, determine which grid has finest spacing
;
grid1Values = *self.grid1.values
grid2Values = *self.grid2.values
dx1 = grid1Values[1] - grid1Values[0]
dx2 = grid2Values[1] - grid2Values[0]
;
; Interpolate data to spacing of the finer grid
;
target = self -> MakeCopy(/noValues, /noErrorBars)
IF(dx1 GT dx2) THEN BEGIN
	target.grid1.values = self.grid2.values
	target.grid1.irange = self.grid2.irange
	target.grid1.range  = self.grid2.range
ENDIF
IF(dx2 GT dx1) THEN BEGIN
	target.grid2.values = self.grid1.values
	target.grid2.irange = self.grid1.irange
	target.grid2.range  = self.grid1.range
ENDIF
result = self -> Interpolate(target)
result -> Restrict
;
; Find indices to center of rotation
;
x0=0.
xx0 = GetKeyWord('x0', extra)
IF(TypeOf(xx0) NE 7) THEN x0=xx0
xx0 = GetKeyWord(result.grid1.mnemonic + '_0', extra)
IF(TypeOf(xx0) NE 7) ThEN x0=xx0
y0=0.
yy0 = GetKeyWord('y0', extra)
IF(TypeOf(yy0) NE 7) THEN y0=yy0
yy0 = GetKeyWord(result.grid2.mnemonic + '_0', extra)
IF(TypeOf(yy0) NE 7) ThEN y0=yy0

x1Values = *result.grid1.values
temp1 = (x1Values - x0)^2
eps1 = MIN(temp1, i)	; the point here is to get subscript of minimum ('i'), NOT to find the minimum ('eps1')

x2Values = *result.grid2.values
temp2 = (x2Values - y0)^2
eps2 = MIN(temp2, j)	; the point here is to get subscript of minimum ('j'), NOT to find the minimum ('eps2')
;
; For each value on grid3, get values from 'result' and rotate them about grid3 through the angle 'arg'
;
irange3=result.grid3.irange	
n3=irange3[1] - irange3[0]
values = *result.values
valuesInfo = SIZE(values)
newValues = MAKE_ARRAY(SIZE=valuesInfo)
FOR k=0,n3 DO BEGIN
	temp = values[*,*,k]
	newValues[*,*,k] = ROT(temp, arg, 1.,i, j, /INTERP, MISSING=0., /PIVOT)
ENDFOR 
PTR_FREE, result.values
result.values = PTR_NEW(newValues)
IF PTR_VALID(result.ErrorBars) THEN BEGIN
	errorBars = *result.errorBars
	errorBarsInfo = SIZE(errorBars)
	newErrorBars = MAKE_ARRAY(SIZE=errorBarsInfo)
	FOR k=0,n3 DO BEGIN
		temp = errorBars[*,*,k]
		newErrorBars[*,*,k] = ROT(temp, arg, 1.,i, j, /INTERP, MISSING=0., /PIVOT)
	ENDFOR 
	PTR_FREE, result.errorBars
	result.errorBars = PTR_NEW(newErrorBars)
ENDIF
;
; Update titles, mnemonics, etc.
;
title1='x!D1!N'
nTitle1 = GetKeyWord('title1', extra)
IF(TypeOf(ntitle1) EQ 7) THEN BEGIN
	IF(ntitle1 NE 'undefined') THEN title1=ntitle1
ENDIF 
result.grid1.title = title1

mnemonic1 = 'x1'
nMnemonic1 = GetKeyWord('mnemonic1', extra)
IF(TypeOf(nMnemonic1) EQ 7) THEN BEGIN
	IF(nMnemonic1 NE 'undefined') THEN mnemonic1=nMnemonic1
ENDIF
result.grid1.mnemonic = mnemonic1

nUnits1 = GetKeyWord('units1', extra)
IF(TypeOF(nUnits1) EQ 7) THEN BEGIN
	IF(nUnits1 NE 'undefined') THEN result.grid1.units=nUnits1
ENDIF

title2='x!D2!N'
nTitle2 = GetKeyWord('title2', extra)
IF(TypeOf(ntitle2) EQ 7) THEN BEGIN
	IF(ntitle2 NE 'undefined') THEN title2=ntitle2
ENDIF 
result.grid2.title = title2

mnemonic2 = 'x2'
nMnemonic2 = GetKeyWord('mnemonic2', extra)
IF(TypeOf(nMnemonic2) EQ 7) THEN BEGIN
	IF(nMnemonic2 NE 'undefined') THEN mnemonic2=nMnemonic2
ENDIF
result.grid2.mnemonic = mnemonic2

nUnits2 = GetKeyWord('units2', extra)
IF(TypeOF(nUnits2) EQ 7) THEN BEGIN
	IF(nUnits2 NE 'undefined') THEN result.grid2.units=nUnits2
ENDIF

RETURN, result

END ; ****** GKVs3D::Rotate ****** ;




