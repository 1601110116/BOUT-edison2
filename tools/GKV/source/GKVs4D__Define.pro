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
; Defining routine for 4-dimensional signal objects
;
; Written by W.M. Nevins
;	1/31/00
;
FUNCTION GKVs4D::AxisNumber, stringIn, Debug=d
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
result = self -> GKVs3D::AxisNumber(stringIn, debug=debug)
axis4Mnemonic = STRTRIM(self.Grid4.mnemonic, 2)
IF( STRCMP(stringIn, axis4Mnemonic, /FOLD_CASE) ) THEN result = 4
RETURN, result
END ; ****** GKVs4D::AxisNumber ****** ;


FUNCTION GKVs4D::GetAxis, structure, Nstructure, Debug=d
;
; Searches 'structure' for tags which are the same as the mnemonics of axis1, axis2, axis3, and axis4
; Returns a structure with tags 'axis1', 'axis2', 'axis3', and 'axis4'.  
; The associated values are whatever values were associated with the corresponding mnemonic.
; 
; Written by W.M. Nevins
;	2/20/00
;
debug=0
otherTags = 0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
;
; Check for axis1, axis2, and axis3 mnemonics by calling GKVs3D::GetAxis
;
result3 = self -> GKVs3D::GetAxis(structure, Debug=debug)
;
; Structure now contains input structure with tag to mnemonic of axis1, axis2 and axis3 (if found) deleted
;
result = {axis1:result3.axis1, axis2:result3.axis2, axis3:result3.axis3, axis4:'no match'}
IF(N_ELEMENTS(structure) EQ 0) THEN RETURN, result			; Structure is undefined
IF(TypeOf(structure) EQ 2) THEN RETURN, result				; No tags were left to return
nTags = N_TAGS(structure)
IF(nTags EQ 0) THEN RETURN, result
tagNames = TAG_NAMES(structure)
tagNames = STRTRIM(tagNames, 2)							; Remove both leading and trailing blanks
command_str = 'structure = {'
axisMnemonic = STRTRIM(self.Grid4.mnemonic)				; Remove both leading and trailing blanks
FOR i=0, ntags-1 DO BEGIN								; Search tags of 'structure'
	IF( STRCMP(axisMnemonic, tagNames[i], /FOLD_CASE) ) THEN BEGIN
		result = {axis1:result3.axis1, axis2:result3.axis2, axis3:result3.axis3, axis4:structure.(i) }	
												; If multiple occurances of axis4 mnemonic in 'structure'
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
END ; ****** GKVs4D::GetAxis ****** ;

FUNCTION GKVs4D::AxisIndex, axisNumber, value, Debug=d
;
; We should only get here if axisNumber is � 4.  
; Function returns index of axisNumber^th  Grid closest to 'value'
;
debug=0
IF(N_ELEMENTS(d) NE 0) then debug=d
IF(axisNumber GT 4) THEN BEGIN
	MESSAGE, 'axisNumber � 4', INFORMATIONAL=d
	RETURN, 0
ENDIF
axisStr = STRING(axisNumber,  FORMAT='(I1)')
commandStr = 'temp = (self.Grid' + axisStr + '.values - values)^2'
ok = EXECUTE(commandStr)
eps=MIN(temp,index)
RETURN, index
END ; ****** GKVs4D::AxisIndex ****** ;


FUNCTION GKVs4D::GetValues, All=all, open=open
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

kkmax = N_ELEMENTS(*self.Grid3.values) - 1
IF(KEYWORD_SET(open)) THEN BEGIN
	IF(self.grid3.boundary EQ "periodic (closed)") THEN kkmax=kkmax-1 
ENDIF
kmin = self.Grid3.irange[0]
kmax = self.Grid3.irange[1] < kkmax
IF N_ELEMENTS(all) THEN BEGIN
	kmin = 0
	kmax = kkmax
ENDIF

llmax = N_ELEMENTS(*self.Grid4.values) - 1
IF(KEYWORD_SET(open)) THEN BEGIN
	IF(self.grid4.boundary EQ "periodic (closed)") THEN llmax=llmax-1 
ENDIF
lmin = self.Grid4.irange[0]
lmax = self.Grid4.irange[1] < llmax
IF N_ELEMENTS(all) THEN BEGIN
	lmin = 0
	lmax = llmax
ENDIF
values = (*self.values)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
RETURN, PTR_NEW(values)
END ; ****** GKVs4D::GetValues ******


FUNCTION GKVs4D::GetErrors, All=all
;
; Returns pointer to an array containing errorBars which fall within the 
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
jmin = self.Grid2.irange[0]
jmax = self.Grid2.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	jmin = 0
	jmax = N_ELEMENTS(self.Grid2.values) - 1
ENDIF
kmin = self.Grid3.irange[0]
kmax = self.Grid3.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	kmin = 0
	kmax = N_ELEMENTS(self.Grid3.values) - 1
ENDIF
lmin = self.Grid4.irange[0]
lmax = self.Grid4.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	lmin = 0
	lmax = N_ELEMENTS(self.Grid4.values) - 1
ENDIF
errors = (*self.ErrorBars)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
RETURN, PTR_NEW(errors)
END ; ****** GKVs4D::GetErrors ******


FUNCTION GKVs4D::SameGrid, argObj, All=all
;
; Returns 1 if 'self' and argObj have same grid, 0 otherwise
;
;	Written 2/6/00
;	by W.M. Nevins
;
selfDims = self  -> NumDims()
argDims  = argObj-> NumDims()
IF(argDims NE selfDims) THEN RETURN, 0
imin = self.Grid4.irange[0]
imax = self.Grid4.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	imin = 0
	imax = N_ELEMENTS(*self.Grid4.values) - 1
ENDIF
sGrid = (*self.Grid4.values)[imin:imax]
iargmin = argObj.Grid4.irange[0]
iargmax = argObj.Grid4.irange[1]
IF N_ELEMENTS(all) THEN BEGIN
	iargmin = 0
	iargmax = N_ELEMENTS(*argObj.Grid4.values) - 1
ENDIF
IF((imax-imin) NE (iargmax-iargmin)) THEN RETURN, 0
argGrid = (*argObj.Grid4.values)[iargmin:iargmax]
Err = TOTAL((sGrid - argGrid)^2)/(1+imax-imin)
L1 = sGrid[imax-imin] - sGrid[0]
d1 = L1/(1+imax-imin)
IF(ERR LT 1.e-3*(d1^2)) THEN RETURN, self -> GKVs3D::SameGrid(argObj)
RETURN, 0
END ; ****** GKVs4D::SameGrid ****** ;



FUNCTION GKVs4D::Interpolate, arg
;
; Interpolate values from self onto the signal window of arg's grid.
;
; Written by W.M. Nevins
;     1/13/09
IF (OBJ_ISA(arg, 'GKVs4D') NE 1) THEN BEGIN
	MESSAGE, "Argument is not a valid GKVs4D object"
	RETURN, 0
ENDIF
;
; Check for common units, independent variable, interval
;
IF (	(self.Grid1.units NE arg.Grid1.units) OR 	$
	(self.Grid2.units NE arg.Grid2.units) OR	$
	(self.Grid3.units NE arg.Grid3.units) OR	$
	(self.Grid4.units NE arg.Grid4.units)		) 	THEN BEGIN
	MESSAGE, "Incompatible units (independent variables)"
	RETURN, 0
ENDIF
IF (	(self.Grid1.title NE arg.Grid1.title) OR 	$
	(self.Grid2.title NE arg.Grid2.title) OR	$
	(self.Grid3.title NE arg.Grid3.title) OR	$
	(self.Grid4.units NE arg.Grid4.units)		) 	THEN BEGIN
	MESSAGE, "Interpolate:  Incompatible independent variables"
	RETURN, 0
ENDIF
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
z1 = *self.Grid3.values
z2 =  *arg.Grid3.values
t1 = *self.Grid4.values
t2 =  *arg.Grid4.values
values1 = *self.values					; Get values of to be interpolated
imin = arg.Grid1.irange[0]				; Get range of indices on target grid
imax = arg.Grid1.irange[1]
jmin = arg.Grid2.irange[0]
jmax = arg.Grid2.irange[1]
kmin = arg.Grid3.irange[0]
kmax = arg.Grid3.irange[1]
lmin = arg.Grid4.irange[0]
lmax = arg.Grid4.irange[1]
x = x2[imin:imax]					; axis 1 for target grid
y = y2[jmin:jmax]					; axis 2 for target grid
z = z2[kmin:kmax]					; axis 3 for target grid
t = t2[lmin:lmax]					; axis 4 for target grid
;
; Get index into old grid arrays (that is, those of 'self') for each element of the new grid arrays (that is, those of 'arg').
;
jjx = VALUE_LOCATE(x1,x)
jjy = VALUE_LOCATE(y1,y)
jjz = VALUE_LOCATE(z1,z)
jjt = VALUE_LOCATE(t1,t)
;
; jjx, jjy, jjz will be = -1 if new grid point does not lie within range of old grid points.
; First create jx, jy, jz arrays (with 
jx = jjx > 0
jy = jjy > 0
jz = jjz > 0
jt = jjt > 0
;
; Compute old grid spacing
;
info = SIZE(values1)
nx1 = info[1] - 1
ny1 = info[2] - 1
nz1 = info[3] - 1
nt1 = info[4] - 1
dx1 = x1[1:nx1] - x1[0:(nx1-1)]
dy1 = y1[1:ny1] - y1[0:(ny1-1)]
dz1 = z1[1:nz1] - z1[0:(nz1-1)]
dt1 = t1[1:nt1] - t1[0:(nt1-1)]
;
; Compute fractional "address" of new x and y grid points within old grid
;
sx = jx + ((x - x1[jx])/dx1[jx])*(jjx GT 0)		; the factor (jjx GT 0) has the effect of zeroing out
sy = jy + ((y - y1[jy])/dy1[jy])*(jjy GT 0)		;  this term for points outside of target grid.
sz = jz + ((z - z1[jz])/dz1[jz])*(jjz GT 0)
st =      ((t - t1[jt])/dt1[jt])*(jjt GT 0)             ; just include residual for 4th variable
;
; We will do the interpolation in two stages:
;     First use native IDL function INTERPOLATE to interpolate 1st three dimensions
;     Then  use linear interpolation in 4th dimension
nx = N_ELEMENTS(x)
ny = N_ELEMENTS(y)
nz = N_ELEMENTS(z)
nt = N_ELEMENTS(t)
jtp1 = jt < (nt-1)
type = typeOf(values1)
values = MAKE_ARRAY(nx,ny,nz,nt, TYPE=type)
;
; Perform 3D interpolation at neighboring values of 4th independent variable,
; then linearly interpolate these results in 4th independent variable.
;
FOR i=0,nt-1 DO $
	values[*,*,*,i] = INTERPOLATE(values1[*,*,*,  jt[i]],sx,sy,sz, /GRID)*(1-st[i]) $
			+ INTERPOLATE(values1[*,*,*,jtp1[i]],sx,sy,sz, /GRID)*st[i]

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
result.Grid2.irange = [0, jmax -jmin]

PTR_FREE, result.Grid3.values				; Load Grid3
result.Grid3.values = PTR_NEW(z)			;	...
result.Grid3.uniform = GKVsd_UniformGrid(z)
result.Grid3.range = arg.Grid3.range
result.Grid3.irange = [0, kmax-kmin]

PTR_FREE, result.Grid4.values				; Load Grid4
result.Grid4.values = PTR_NEW(t)			;	...
result.Grid4.uniform = GKVsd_UniformGrid(t)
result.Grid4.range = arg.Grid4.range
result.Grid4.irange = [0, lmax-lmin]
;
; we will make no attempt to interpolate error bars ...
;

RETURN, result
END ; ***** GKVs4D::Interpolate ***** ;


FUNCTION GKVs4D::Plus, argg, title=title, mnemonic=mnemonic, units=units
;
; Define basic arithmetic operation:  addition
;
; Adds argument to self. 
; 
; If argument is of STRING type, then this string
; is compared to the four axis mnemonics. If there
; is a match, then at each value of the corresponding
; dependent variable, the value of 'self' at the remaining
; three dependent variables is incremented by the corresponding
; value of this axis.
;
; If argument is a GKVs4D object, then argument values are
; interpolated onto 'self''s grid and summed.
;
; Keywords:
;
;	title		Set to desired title for the sum. 
;			(Optional)
;
;	mnemonic	Set to desired mnemonic for the sum.
;			(Optional)
;
;	units		Set to desired units for the sum.
;			(Optional)
;
; Written by W.M. Nevins
;	2/11/09
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0		; Undefined argument
arg=argg					; Make proxy for argument
try_again:
argInfo = SIZE(arg)				; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN					; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						kVals = N_ELEMENTS(gridValues)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = gridValues#Replicate(1,lVals)
						ENDFOR
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(gridValues)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = Replicate(1,kVals)#gridValues
						ENDFOR
						argTitle = self.Grid2.title
					END	
			3	:	BEGIN
						gridValues = *self.Grid3.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(gridValues)
						nVals = N_ELEMENTS(*slef.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = gridValues#REPLICATE(1, nVals)
						ENDFOR
						argTitle = self.Grid3.title
					END	
			4	:	BEGIN
						gridValues = *self.Grid4.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(gridValues)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN 
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = Replicate(1,mVals)#gridValues
						ENDFOR
						argTitle = self.Grid4.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORMATIONAL
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
		result.mnemonic = self.mnemonic + 'p' + arg
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "+" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN							; Argument is a POINTER
		arg = *arg					; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN							; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs4D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs4D objects
		;
			result = arg -> Interpolate(self)
	
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			kmin = self.Grid3.irange[0]
			kmax = self.Grid3.irange[1]
			lmin = self.Grid4.irange[0]
			lmax = self.Grid4.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
			values = argValues + selfValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			result.mnemonic = self.mnemonic + "p" + arg.mnemonic
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = '(' + self.Title + "+" + arg.Title + ')'
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units	
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Plus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs4D::Plus ***** ;


FUNCTION GKVs4D::Minus, argg, title=title, mnemonic=mnemonic, units=units
;
; Define basic arithmetic operation:  subtraction
;
; Subtracts argument from self. 
; 
; If argument is of STRING type, then this string
; is compared to the four axis mnemonics. If there
; is a match, then at each value of the corresponding
; dependent variable, the value of 'self' at the remaining
; three dependent variables is decremented by the corresponding
; value of this axis.
;
; If argument is a GKVs4D object, then argument values are
; interpolated onto 'self''s grid and subtracted.
;
; Keywords:
;
;	title		Set to desired title for the difference. 
;			(Optional)
;
;	mnemonic	Set to desired mnemonic for the difference.
;			(Optional)
;
;	units		Set to desired units for the difference.
;			(Optional)
;
; Written by W.M. Nevins
;	2/11/09
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0		; Undefined argument
arg=argg					; Make proxy for argument
try_again:
argInfo = SIZE(arg)				; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN					; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						kVals = N_ELEMENTS(gridValues)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = gridValues#Replicate(1,lVals)
						ENDFOR
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(gridValues)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = Replicate(1,kVals)#gridValues
						ENDFOR
						argTitle = self.Grid2.title
					END	
			3	:	BEGIN
						gridValues = *self.Grid3.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(gridValues)
						nVals = N_ELEMENTS(*slef.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = gridValues#REPLICATE(1, nVals)
						ENDFOR
						argTitle = self.Grid3.title
					END	
			4	:	BEGIN
						gridValues = *self.Grid4.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(gridValues)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN 
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = Replicate(1,mVals)#gridValues
						ENDFOR
						argTitle = self.Grid4.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORMATIONAL
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
		result.mnemonic = self.mnemonic + 'm' + arg
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "-" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN							; Argument is a POINTER
		arg = *arg					; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN							; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs4D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs4D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			kmin = self.Grid3.irange[0]
			kmax = self.Grid3.irange[1]
			lmin = self.Grid4.irange[0]
			lmax = self.Grid4.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
			values = selfValues - argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			result.mnemonic = self.mnemonic + 'm' + arg.mnemonic
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = '(' + self.Title + "-" + arg.Title + ')'
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units	
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Minus(arg)

ENDCASE

RETURN, result

END ; ***** GKVs4D::Minus ***** ;


FUNCTION GKVs4D::Times, argg, title=title, mnemonic=mnemonic, units=units
;
; Define basic arithmetic operation:  Multiplication
;
; Multiplies self by argument. 
; 
; If argument is of STRING type, then this string
; is compared to the four axis mnemonics. If there
; is a match, then at each value of the corresponding
; dependent variable, the value of 'self' at the remaining
; three dependent variables is multiplied by the corresponding
; value of this axis.
;
; If argument is a GKVs4D object, then argument values are
; interpolated onto 'self''s grid and then multiplied.
;
; Keywords:
;
;	title		Set to desired title for the product. 
;			(Optional)
;
;	mnemonic	Set to desired mnemonic for the product.
;			(Optional)
;
;	units		Set to desired units for the product.
;			(Optional)
;
; Written by W.M. Nevins
;	2/11/09
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0		; Undefined argument
arg=argg					; Make proxy for argument
try_again:
argInfo = SIZE(arg)				; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN					; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						kVals = N_ELEMENTS(gridValues)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = gridValues#Replicate(1,lVals)
						ENDFOR
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(gridValues)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = Replicate(1,kVals)#gridValues
						ENDFOR
						argTitle = self.Grid2.title
					END	
			3	:	BEGIN
						gridValues = *self.Grid3.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(gridValues)
						nVals = N_ELEMENTS(*slef.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = gridValues#REPLICATE(1, nVals)
						ENDFOR
						argTitle = self.Grid3.title
					END	
			4	:	BEGIN
						gridValues = *self.Grid4.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(gridValues)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN 
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = Replicate(1,mVals)#gridValues
						ENDFOR
						argTitle = self.Grid4.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues)
		selfValues = *self.values
		resultValues = selfValues*argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		result.mnemonic = self.mnemonic + 't' + arg
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "*" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN							; Argument is a POINTER
		arg = *arg					; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN							; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs4D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs4D objects
		;
			result = arg -> Interpolate(self)
			
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			kmin = self.Grid3.irange[0]
			kmax = self.Grid3.irange[1]
			lmin = self.Grid4.irange[0]
			lmax = self.Grid4.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
			values = selfValues*argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			result.mnemonic = self.mnemonic + 't' + arg.mnemonic
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = '(' + self.Title + "*" + arg.Title + ')'
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units	
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Times(arg)

ENDCASE

RETURN, result

END ; ***** GKVs4D::Times ***** ;



FUNCTION GKVs4D::Over, argg, title=title, mnemonic=mnemonic, units=units
;
; Define basic arithmetic operation:  Division
;
; Divides self by argument. 
; 
; If argument is of STRING type, then this string
; is compared to the four axis mnemonics. If there
; is a match, then at each value of the corresponding
; dependent variable, the value of 'self' at the remaining
; three dependent variables is divided by the corresponding
; value of this axis.
;
; If argument is a GKVs4D object, then argument values are
; interpolated onto 'self''s grid before dividing.
;
; Keywords:
;
;	title		Set to desired title for the ratio. 
;			(Optional)
;
;	mnemonic	Set to desired mnemonic for the ratio.
;			(Optional)
;
;	units		Set to desired units for the ratio.
;			(Optional)
;
; Written by W.M. Nevins
;	2/11/09
;
arggInfo=SIZE(argg)
typeIndex = arggInfo[0] + 1
arggType = arggInfo(typeIndex)
IF(arggType EQ 0) THEN RETURN, 0		; Undefined argument
arg=argg					; Make proxy for argument
try_again:
argInfo = SIZE(arg)				; Get argument variable type
argDims = argInfo[0]
typeIndex = argDims + 1
argType = argInfo[typeIndex]

CASE argType OF

7:	BEGIN					; Argument is a String.  Check if it is a valid axis mnemonic
		iaxis = self -> AxisIrange(arg)
		CASE iaxis OF
			1	:	BEGIN
						gridValues = *self.Grid1.values
						kVals = N_ELEMENTS(gridValues)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = gridValues#Replicate(1,lVals)
						ENDFOR
						argTitle = self.Grid1.title
					END
			2	:	BEGIN
						gridValues = *self.Grid2.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(gridValues)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(*self.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR m=0, mVals-1 DO BEGIN
							FOR n=0, nVals-1 DO argValues[*,*,m,n] = Replicate(1,kVals)#gridValues
						ENDFOR
						argTitle = self.Grid2.title
					END	
			3	:	BEGIN
						gridValues = *self.Grid3.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(gridValues)
						nVals = N_ELEMENTS(*slef.Grid4.values)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = gridValues#REPLICATE(1, nVals)
						ENDFOR
						argTitle = self.Grid3.title
					END	
			4	:	BEGIN
						gridValues = *self.Grid4.values
						kVals = N_ELEMENTS(*self.Grid1.values)
						lVals = N_ELEMENTS(*self.Grid2.values)
						mVals = N_ELEMENTS(*self.Grid3.values)
						nVals = N_ELEMENTS(gridValues)
						argValues = FLTARR(kVals, lVals, mVals, nVals)
						FOR k=0, kVals-1 DO BEGIN 
							FOR l=0, lVals-1 DO argValues[k,l,*,*] = Replicate(1,mVals)#gridValues
						ENDFOR
						argTitle = self.Grid4.title
					END	
			ELSE	:	BEGIN
						MESSAGE, "Illegal argument -- Returning", /INFORMATIONAL
						RETURN, 0
					END
		ENDCASE
		result = self -> MakeCopy(/noValues)
		selfValues = *self.values
		resultValues = selfValues/argValues
		result.values =  PTR_NEW(resultValues)
		vmin = GKVsd_MIN(resultValues, MAX=vmax)
		IF(vmax EQ vmin) THEN vmax = vmin+1
		result.vrange = [vmin, vmax]
		result.mnemonic = self.mnemonic + 'o' + arg
		IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
		result.Title = '(' + self.Title + "/" + argTitle + ')'
		IF(TypeOf(title) EQ 7) THEN result.title=title
		IF(TypeOf(units) EQ 7) THEN result.units=units	
	END

10:	BEGIN							; Argument is a POINTER
		arg = *arg					; Dereference pointer
		GOTO, try_again					;	and try again.
	END
11:	BEGIN							; Argument is an Object Reference
		IF (OBJ_ISA(arg, 'GKVs4D') EQ 1) THEN BEGIN
			IF (self.units NE arg.units) THEN BEGIN
				MESSAGE, "Incompatible units (dependent variable)", /INFORMATIONAL
				RETURN, 0
			ENDIF
		;
		; Sum data from two GKVs4D objects
		;
			result = arg -> Interpolate(self)
			IF ( NOT OBJ_VALID(result) ) THEN BEGIN
				MESSAGE, "Can't form common grid for independent variables", /INFORMATIONAL
				RETURN, 0
			ENDIF
			argValues = *result.values
			imin = self.Grid1.irange[0]
			imax = self.Grid1.irange[1]
			jmin = self.Grid2.irange[0]
			jmax = self.Grid2.irange[1]
			kmin = self.Grid3.irange[0]
			kmax = self.Grid3.irange[1]
			lmin = self.Grid4.irange[0]
			lmax = self.Grid4.irange[1]
			selfValues = (*self.values)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
			values = selfValues/argValues
			PTR_FREE, result.values
			result.values = PTR_NEW(values)
			vmin = GKVsd_MIN(values, MAX=vmax)
			IF(vmax EQ vmin) THEN vmax = vmin+1
			result.vrange = [vmin, vmax]
			result.mnemonic = self.mnemonic + 'o' + arg.mnemonic
			IF(TypeOf(mnemonic) EQ 7) THEN result.mnemonic = STRCOMPRESS(mnemonic, /REMOVE_ALL)
			result.Title = '(' + self.Title + "/" + arg.Title + ')'
			IF(TypeOf(title) EQ 7) THEN result.title=title
			IF(TypeOf(units) EQ 7) THEN result.units=units	
		ENDIF 
			
	END
ELSE:		result = self -> GKVsd::Over(arg)

ENDCASE

RETURN, result

END ; ***** GKVs4D::Over ***** ;


FUNCTION GKVs4D::XCORR, Ref = RefObj, All=allin, 		  	$
			Lw=lagWindow, Dw=dataWindow, _Extra=extra
;
; Returns a GKVsd object containing the cross correlation function
; between the data in 'self' and the data in RefObj (both must be
; GKV objects).  This routine handles only the case where
; 'self' is a GKVs4D object and RefObj is a GKVs2D object.
; Other cases are handled by GKVs2D::XCORR.
;
; Written by W.M. Nevins
;     1/16/09
;
; Check for presence of RefObj
;
IF(TypeOF(RefObj) EQ 0) THEN BEGIN
	result =  self -> GKVs2D::XCORR(_EXTRA=extra)
	RETURN, result
ENDIF
;
; Check if RefObj is an object
;
IF(TypeOF(RefObj) NE 11) THEN BEGIN
	MESSAGE, 'Ref is not an Object', /Informational
	RETURN, 0
ENDIF
;
; Check if RefObj is a GKV object
;
IF(NOT OBJ_ISA(RefObj, "GKVsd")) THEN BEGIN
	MESSAGE, 'Ref is not a GKV Object', /Informational
	RETURN, 0
ENDIF
;
; Check dimensionality of RefObj
;
refDims = RefObj -> NumDims()
IF(refDims EQ 3) THEN BEGIN
	result = self -> GKVs2D::XCORR(ref=RefObj, _Extra=extra)
	RETURN, result
ENDIF
IF(refDims LT 2) THEN BEGIN
	MESSAGE, 'ref must be GKVs2D, GKVs3D or GKVs4D object', /Informational
	RETURN, 0
ENDIF
nDims = self -> NumDims()
IF(nDims NE 4) THEN BEGIN
	MESSAGE, 'XCORR not implimented for more than 4 dimensions', /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Now, deal with the case in which 'self' is 4D and 'ref' is 2D
;
all = [0,0]
IF(Query_Integer(allin)) THEN all = allin
IF(N_ELEMENTS(all) EQ 1) THEN all = [all,all]
;
; Find grid spacing and number of elements 
; within the "signalwindow" for each dimension
;
spacing = FLTARR(nDims+1)
nGrids  = LONARR(nDims+1)
FOR i=1, nDims DO BEGIN
	dimStr = STRING(i, FORMAT='(I1)')			
	GridStr = 'self.Grid' + dimStr
	command_str = 'Grid = ' + gridStr	; Copy ith instance of the Grid Class
	ok = EXECUTE(command_str)		; into Grid
	gridValues = *grid.values
	spacing[i] = gridValues[1] - gridValues[0]
	nGrids[i]  = grid.irange[1] - grid.irange[0] + 1
ENDFOR
;
; Find "ignorable" independent variables. That is, the 
; Grid structures from 'self' which match the two Grid structures 
; in RefObj
;
;
; Find 1st "ignorable" independent variable.
;
refGrid1 = RefObj.Grid1
refGrid2 = RefObj.Grid2
FOR iVar1=1,nDims DO BEGIN
	iVar1str = STRING(iVar1, FORMAT='(I1)')
	command = 'selfGrid1 = self.Grid' + iVar1str
	ok = EXECUTE(command)
	iVar1min = selfGrid1.irange[0]
	iVar1max = selfGrid1.irange[1]
	IF(GKV_GridSame(selfGrid1,RefGrid1,All=all[0], /force)) THEN GOTO, Done1
ENDFOR
;
; only get here if ignorable coordinates don't match
;
MESSAGE, 'No match between self grids and ref.Grid1, returning', /INFORMATIONAL
RETURN, 0
;
; Have match between refGrid 1 and selfgrid ivar1
;
Done1	:
IF(KEYWORD_SET(all[0])) THEN BEGIN
	iVar1min=0
	iVar1max = N_ELEMENTS(*selfGrid1.values)
ENDIF
uniform = selfGrid1.uniform
IF(NOT uniform) THEN BEGIN
	MESSAGE, 'WARNING:  Grid ' + iVar1Str + ' is nonuniform.  Result may be invalid!', /INFORMATIONAL
ENDIF
spacing[iVar1] = (*selfGrid1.values)[iVar1min+1] - (*selfGrid1.values)[iVar1min]
nGrids[iVar1] = iVar1max - iVar1min + 1
;
; Find 2nd "ignorable" independent variable
;
FOR iVar2=iVar1+1,nDims DO BEGIN
	iVar2str = STRING(iVar2, FORMAT='(I1)')
	command = 'selfGrid2 = self.Grid' + iVar2str
	ok = EXECUTE(command)
	iVar2min = selfGrid2.irange[0]
	iVAr2max = selfGrid2.irange[1]
	IF(GKV_GridSame(selfGrid2,RefGrid2,All=all[1], /force)) THEN GOTO, Done2
ENDFOR
;
; only get here if ignorable coordinates don't match
;
MESSAGE, 'No match between self grids and ref.Grid2, returning', /INFORMATIONAL
RETURN, 0
;
; Have match between refGrid 2 and selfgrid ivar2
;
Done2	:
IF(KEYWORD_SET(all[1])) THEN BEGIN
	iVar2min=0
	iVar2max=N_ELEMENTS(*selfGrid2.values)
ENDIF
uniform = selfGrid2.uniform
IF(NOT uniform) THEN BEGIN
	MESSAGE, 'WARNING:  Grid ' + iVar2Str + ' is nonuniform.  Result may be invalid!', /Informational
ENDIF
spacing[iVar2] = (*selfGrid2.values)[iVar2min+1] - (*selfGrid2.values)[iVar2min]
nGrids[iVar2] = iVar2max - iVar2min + 1
;
;  Find non-ignorable coordinates
;
CASE iVar1 OF
	1:	BEGIN
	CASE iVar2 OF
		2:	BEGIN
			ifix1 = 3
			ifix2 = 4
			irange = self.grid3.irange
			jrange = self.grid4.irange	
			END
		3:	BEGIN
			ifix1 = 2
			ifix2 = 4
			irange = self.grid2.irange
			jrange = self.grid4.irange
			END
		4:	BEGIN
			ifix1 = 2
			ifix2 = 3
			irange = self.grid2.irange
			jrange = self.grid3.irange
			END
	ENDCASE
		END
	2:	BEGIN
	CASE iVar2 OF
		3:	BEGIN
			ifix1 = 1
			ifix2 = 4
			irange = self.grid1.irange
			jrange = self.grid4.irange
			END
		4:	BEGIN
			ifix1 = 1
			ifix2 = 3
			irange = self.grid1.irange
			jrange = self.grid3.irange
			END
	ENDCASE
		END
	3:	BEGIN
			ifix1 = 1
			ifix2 = 2
			irange = self.grid1.irange
			jrange = self.grid2.jrange
		END
ENDCASE
dt = spacing[iVar2]
ilw = ngrids[iVar2]/2  > 1
idw = ngrids[iVar2]/10 < SQRT(ngrids[iVar2])
idw = idw > 1
IF(N_ELEMENTS( lagWindow) GT 0) THEN ilw = FIX(lagWindow/dt)
IF(N_ELEMENTS(dataWindow) GT 0) THEN idw = FIX(dataWindow/dt)
;
; get pointer to self's values and actual ref values
;
selfValuePtr =   self -> GetValues()
refValuePtr = refObj -> GetValues()
refValues    = *refValuePtr
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

FOR i = irange[0], irange[1] DO BEGIN
	FOR j = jrange[0], jrange[1] DO BEGIN

CASE iVar1 OF
	1:	BEGIN
	CASE iVar2 OF
		2:	BEGIN
			selfValues = REFORM((*selfValuePtr)[*,*,i,j], refInfo[1], RefInfo[2])
			END
		3:	BEGIN
			selfValues = REFORM((*selfValuePtr)[*,i,*,j], refInfo[1], RefInfo[2])
			END
		4:	BEGIN
			selfValues = REFORM((*selfValuePtr)[*,i,j,*], refInfo[1], RefInfo[2])
			END
	ENDCASE
		END
	2:	BEGIN
	CASE iVar2 OF
		3:	BEGIN
			selfValues = REFORM((*selfValuePtr)[i,*,*,j], refInfo[1], RefInfo[2])
			END
		4:	BEGIN
			selfValues = REFORM((*selfValuePtr)[i,*,j,*], refInfo[1], RefInfo[2])
			END
	ENDCASE
		END
	3:	BEGIN
			selfValues = REFORM((*selfValuePtr)[i,j,*,*], refInfo[1], RefInfo[2])
		END
ENDCASE
IF(TypeOf(extra) EQ 8) THEN BEGIN
	CorrValues_ij = XCorr(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, Dw=dataWindow, _Extra=extra)
ENDIF ELSE BEGIN
	CorrValues_ij = XCorr(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, Dw=dataWindow)
ENDELSE
	
IF((i EQ irange[0]) AND (j EQ jrange[0])) THEN BEGIN							; Initialization duties
	selfInfo = SIZE(*selfValuePtr)
	nSelfValueElements = N_ELEMENTS(selfValues)
	corrInfo = SIZE(CorrValues_ij)
	corrDims = corrInfo[0]
	ncorrs = corrinfo[corrDims]
	corrType = 4
	IF( (Query_Complex(selfValues) + Query_Complex(RefValues)) GT 0 ) THEN corrType=6
CASE iVar1 OF
	1:	BEGIN
	CASE iVar2 OF
		2:	BEGIN
			corrValues = MAKE_ARRAY(corrInfo[1], corrInfo[2], selfInfo[3], selfInfo[4], TYPE=corrType)
			END
		3:	BEGIN
			corrValues = MAKE_ARRAY(corrInfo[1], selfInfo[2], corrInfo[2], selfInfo[4], TYPE=corrType)
			END
		4:	BEGIN
			corrValues = MAKE_ARRAY(corrInfo[1], selfInfo[2], selfInfo[3], corrInfo[2], TYPE=corrType)
			END
	ENDCASE
		END
	2:	BEGIN
	CASE iVar2 OF
		3:	BEGIN
			corrValues = MAKE_ARRAY(selfInfo[1], corrInfo[1], corrInfo[2], selfInfo[4], TYPE=corrType)
			END
		4:	BEGIN
			corrValues = MAKE_ARRAY(selfInfo[1], corrInfo[1], selfInfo[3], corrInfo[2], TYPE=corrType)
			END
	ENDCASE
		END
	3:	BEGIN
			corrValues = MAKE_ARRAY(selfInfo[1], selfInfo[2], corrInfo[1], corrInfo[2], TYPE=corrType)
		END
ENDCASE
ENDIF
;
; Apply local norm if required
;
IF(localNorm EQ 1) THEN BEGIN
	IF(Query_Complex(selfValues)) THEN BEGIN
		selfNormSq_ij = TOTAL(selfValues*CONJ(selfValues))/nSelfValueElements
	ENDIF ELSE BEGIN 
		selfNormSq_ij = TOTAL(selfValues*selfValues)/nSelfValueElements
	ENDELSE
	CorrValues_ij = CorrValues_ij/SQRT(refNormSq*selfNormSq_ij)
ENDIF

CASE iVar1 OF
	1:	BEGIN
	CASE iVar2 OF
		2:	BEGIN
			corrValues[*,*,i,j] = corrValues_ij
			END
		3:	BEGIN
			corrValues[*,i,*,j] = corrValues_ij
			END
		4:	BEGIN
			corrValues[*,i,j,*] = corrValues_ij
			END
	ENDCASE
		END
	2:	BEGIN
	CASE iVar2 OF
		3:	BEGIN
			corrValues[i,*,*j] = corrValues_ij
			END
		4:	BEGIN
			corrValues[i,*,j,*] = corrValues_ij
			END
	ENDCASE
		END
	3:	BEGIN
			corrValues[i,j,*,*] = corrValues_ij
		END
ENDCASE
	ENDFOR
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
result -> GridRestrict, iVar1
result -> GridRestrict, iVar2
;
; Start filling tags of result
;
result.mnemonic = 'XCorr_' + self.mnemonic	 + '_' + RefObj.mnemonic	; Set Mnemonic
refIndices = *refObj.Indices
RefObjTitle = refObj.title + '[' + refIndices[iFix1-1] + ',' + refIndices[iFix2-1] + ']'
result.title = 'C{' + self.title + ', ' + RefObjTitle + '}'
result.units = '(' + self.units + ')*(' + refObj.units + ')'			; set units
IF(localNorm EQ 1) THEN result.units = ''
PTR_FREE, result.values								; Set Values pointer
result.values = corrValuesPtr
vmin = GKVsd_MIN(CorrValues, MAX=vmax)						; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange
PTR_FREE, result.ErrorBars					; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on correlation functions?
;
; correct grid structures
; corresponding to ignorable
; dependent variables
;
iVarstr = STRING(iVar1, FORMAT='(I1)')
command  = "varGrid1 = result.grid" + iVar1Str
OK = EXECUTE(command)
mnemonic = vargrid1.mnemonic
title = vargrid1.title
IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
	varGrid1.mnemonic = 'tau'
	varGrid1.title = '!4s!X'
ENDIF ELSE BEGIN
	varGrid1.title = '!4D!X' + title
	varGrid1.mnemonic = 'd' + mnemonic
ENDELSE
PTR_FREE, varGrid1.values						
n = ngrids[ivar1]
gridvals = spacing(ivar1)*(FINDGEN(n) -n/2)	; Generate array of grid values
GridPtr = PTR_NEW(gridvals)	
varGrid1.values = GridPtr			; Set Grid Values pointer
varGrid1.boundary = 'periodic (open)'		; Correlation function BC's are always Periodic?		
varGrid1.uniform = 1				; Correlation functions computed on uniform grid
min = MIN(gridvals, Max=max)			; Set Grid plot Range
range = [min, max]
varGrid1.range = range
irange = [0, n-1]				; Set Grid signal window range, irange
varGrid1.irange = irange
command = "result.Grid" + ivar1str + "=varGrid1"; Copy corrected Grid
ok = EXECUTE(command)				; back into result


ivar2str = STRING(ivar2, FORMAT='(I1)')
command  = "varGrid2 = result.grid" + ivar2Str
OK = EXECUTE(command)
mnemonic = vargrid2.mnemonic
title = vargrid2.title
IF((Mnemonic EQ 't') OR STRCMP(Mnemonic, 'time', 4, /FOLD_CASE)) THEN BEGIN
	varGrid2.mnemonic = 'tau'
	varGrid2.title = '!4s!X'
ENDIF ELSE BEGIN
	varGrid2.title = '!4D!X' + title
	varGrid2.mnemonic = 'd' + mnemonic
ENDELSE
PTR_FREE, varGrid2.values						
n = ngrids[ivar2]
gridvals = spacing(ivar2)*(FINDGEN(n) -n/2)	; Generate array of grid values
GridPtr = PTR_NEW(gridvals)	
varGrid2.values = GridPtr			; Set Grid Values pointer
varGrid2.boundary = 'periodic (open)'		; Correlation function BC's are always Periodic?		
varGrid2.uniform = 1				; Correlation functions computed on uniform grid
min = MIN(gridvals, Max=max)			; Set Grid plot Range
range = [min, max]
varGrid2.range = range
irange = [0, n-1]				; Set Grid signal window range, irange
varGrid2.irange = irange
command = "result.Grid" + ivar2str + "=varGrid2"; Copy corrected Grid
ok = EXECUTE(command)				; back into result

; 
; and we're done ...
;
RETURN, result

END ; ***** FUNCTION GKVs4D::XCORR ****** ;


FUNCTION GKVs4D::MakeCopy, _Extra=extra
;
; Make "deep" copy of self
;
result = self -> GKVs3D::MakeCopy( _Extra=extra)	; Creates new GKVsxD objects, and
							; 	copies elements of GKVs3D class
result.Grid4 = GKVsd_GridCopy(self.Grid4)		; Make 'deep' copy the Grid-class structure Grid4,			
RETURN, result

END ; ***** GKVs4D::MakeCopy ***** ;


FUNCTION GKVs4D_Gen, Nx=n1, Ny=n2, Nz=n3, Nt=n4, Amplitude=a,  	$
				kx=k1, ky=k2, kz=k3, omega=k4, Del_k=bandwidth
;
; Set up sample GKVs2D signal object
;
nx=32L
ny=32L
nz=32L
nt=32L
IF ( N_ELEMENTS(n1) NE 0 ) THEN nx=n1
IF ( N_ELEMENTS(n2) NE 0 ) THEN ny=n2
IF ( N_ELEMENTS(n3) NE 0 ) THEN nz=n3
IF ( N_ELEMENTS(n4) NE 0 ) THEN nt=n4
dx = 2*!PI/(nx-1)
dy = 2*!PI/(ny-1)
dz = 2*!PI/(nz-1)
dt = 2*!PI/(nt-1)

amplitude=1.
IF KEYWORD_SET(a)  then amplitude=a

kx = 2.
ky = 2.
kz = 2.
omega = 2.
Del_k = 2.
IF N_ELEMENTS(k1) THEN kx=k1
IF N_ELEMENTS(k2) THEN ky=k2
IF N_ELEMENTS(k3) THEN kz=k3
IF N_ELEMENTS(k4) THEN omega=k4
IF N_ELEMENTS(bandwidth) then del_k=bandwidth

xgrid = dx*indgen(nx)
ygrid = dy*indgen(ny)
zgrid = dz*indgen(nz)
tgrid = dt*indgen(nt)
i=COMPLEX(0.,1.)
argx= -0.5d*((INDGEN(nx) - kx)^2/Del_k > (-30.0d))
argy= -0.5d*((INDGEN(ny) - ky)^2/Del_k > (-30.0d))
argz= -0.5d*((INDGEN(nz) - kz)^2/Del_k > (-30.0d))
argt= -0.5d*((INDGEN(nt) - omega)^2/Del_k > (-30.0d))
rarray=RANDOMN(seed, nx, ny, nz, nt)
iarray=RANDOMN(seed, nx, ny, nz, nt)
array=COMPLEX(rarray, iarray)
FOR i=0, nz-1 DO	$
	FOR j=0, nt-1 DO $
		array(*,*,i,j)=array(*,*,i,j)*( exp(argx)#exp(argy) )*( exp(argz(i))*exp(argt(j)) )
array=FFT(array, 1, /OVERWRITE)
msarray= TOTAL(array*CONJ(array))/( LONG(nx)*LONG(ny)*LONG(nz)*LONG(nt) ) 
values=amplitude*array/sqrt(msarray)
vmin = GKVsd_MIN(values, MAX=vmax)
IF(vmax EQ vmin) THEN vmax = vmin+1

Indices = REPLICATE('*', 4)

signal={GKVs4D}
;
; GKVsd tags
;
	signal.mnemonic	= "sgen4D"
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
; Grid3 tags
;
	signal.Grid3.mnemonic	= 'z'
	signal.Grid3.title 	= "z-title"
	signal.Grid3.units 	= "z-units"
	signal.Grid3.values 	= PTR_NEW(zgrid)
	signal.Grid3.boundary	= "open"
	signal.Grid3.range	= [zgrid[0], zgrid[nz-1]]
	signal.Grid3.irange	= [0, nz-1]
;
; Grid4 tags
;
	signal.Grid4.mnemonic	= 't'
	signal.Grid4.title 	= "t-title"
	signal.Grid4.units 	= "t-units"
	signal.Grid4.values 	= PTR_NEW(tgrid)
	signal.Grid4.boundary	= "open"
	signal.Grid4.range	= [tgrid[0], tgrid[nt-1]]
	signal.Grid4.irange	= [0, nt-1]
;
; Create signal object
;
	signal_obj=OBJ_NEW("GKVs4D", signal) 
RETURN, signal_obj
END ; ***** GKVs4D_Gen ***** ;


PRO GKVs4D::GET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Get values of elements of the realization of the GRID class, Grid4;
; and then call GKVs3D::GET to get values of elements of the GKVs3D Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs4D,
; or any or its subclasses
;
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 4) THEN 	$
			GKVsd_GetGrid, self.Grid4,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	2:	IF(axisID EQ 4) THEN 	$
			GKVsd_GetGrid, self.Grid4,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	3:	IF(axisID EQ 4) THEN 	$
			GKVsd_GetGrid, self.Grid4,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid4.mnemonic) THEN 	$
			GKVsd_GetGrid, self.Grid4,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
								GridUnits=GridUnits,	GridValues=GridValues, 			$
								boundary=boundary, uniform=uniform,range=range, 		$
								irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self -> GKVs3D::GET, axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
RETURN
END ; ****** GKVs4D::GET ****** ;


PRO GKVs4D::SET,axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
;
; Get values of elements of the realization of the GRID class, Grid4;
; and then call GKVs3D::GET to get values of elements of the GKVs3D Class.
;
; This call will be used polymorphically, so 'self' may be of class GKVs4D,
; or any or its subclasses
;
arg = self.Grid4
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
	1:	IF(axisID EQ 4) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	2:	IF(axisID EQ 4) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	3:	IF(axisID EQ 4) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
	7:	IF(axisID EQ self.Grid4.mnemonic) THEN 	$
			GKVsd_SetGrid, arg,	GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
							GridUnits=GridUnits,	GridValues=GridValues, 			$
							boundary=boundary, uniform=uniform,range=range, 		$
							irange=irange, _Extra=extra
ELSE	 :	; just continue on else 
ENDCASE
self.Grid4 = arg
self -> GKVs3D::SET, 	axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle, 				$
			GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,		$
			range=range, irange=irange, mnemonic=mnemonic, Title=Title, Indices=indices, 		$
			units=units, values=values, vrange=vrange, ErrorBars=ErrorBars, CodeName=CodeName,	$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra
RETURN
END ; ****** GKVs4D::SET ****** ;


PRO GKVs4D::Info
;
; Prints information about contents of GKVs4D objects
;
self -> GKVs3D::Info
PRINT, 'Grid4:'
GKVsd_PrintGrid, self.Grid4
RETURN
END ; ***** GKVs4D::Info ***** ;


Function GKVs4D::INIT, signal
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
; Set tags for GKVs3D class
;
ok = self -> GKVs3D::INIT(signal)
;
; Set Tags for Grid4
;
NumTags=N_TAGS(self.Grid4)
FOR itag=0, NumTags-1 DO $
	self.Grid4.(itag) = signal.Grid4.(itag)
;
; remove embedded blanks in grid title, mnemonic, and units
;
self.grid4.title    = STRCOMPRESS(self.grid4.title,    /REMOVE_ALL)
self.grid4.mnemonic = STRCOMPRESS(self.grid4.mnemonic, /REMOVE_ALL)
self.grid4.units    = STRCOMPRESS(self.grid4.units,    /REMOVE_ALL)


IF(PTR_VALID(self.Grid4.values)) THEN BEGIN
	;
	; Check for uniform grid
	;
	self.Grid4.uniform = GKVsd_UniformGrid(self.Grid4.values)
	;
	; Check if Grid4.range is set ... and set if necessary
	;
	IF((self.Grid4.range[0] EQ 0.) AND (self.Grid4.range[1] EQ 0.)) THEN $
		self.Grid4.range = [MIN(*self.Grid4.Values, Max=max), max]
	;
	; Check if Grid4.irange is set ... and set if necessary
	;
	IF((self.Grid4.irange[0] EQ 0) AND (self.Grid4.irange[1] EQ 0)) THEN $
		self.Grid4.irange = [0, N_ELEMENTS(*self.Grid4.values)-1]
ENDIF
;
; Return on successful completion
;
RETURN, ok
END ; ***** GKVs4D::INIT ***** ;


PRO GKVs4D::CleanUp


PTR_FREE, self.Grid4.values
self -> GKVs3D::CleanUp

RETURN

END ; ***** GKVs4D::CleanUp ***** ;	

PRO GKVs4D__Define
struct = {	GKVs4D,				$	; "GK Visualization signal (4-D)"
		INHERITS GKVs3D,		$	; GKVs4D is a subclass of GKVs3D
		GRID4:{GRID}			}	; Include a fourth 'Grid' class
	
END ; ***** GKVs4D__Define ***** ;
