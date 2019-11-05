

FUNCTION GKVs2D::Interpolate, arg
;
; Interpolate values from self onto the signal window of arg's grid.
;
IF (OBJ_ISA(arg, 'GKVs2D') NE 1) THEN BEGIN
	MESSAGE, "Argument is not a valid GKVs2D object"
	RETURN, 0
ENDIF
;
; Check for common units, independent variable, interval
;
IF ((self.Grid1.units NE arg.Grid1.units) OR (self.Grid2.units NE arg.Grid2.units)) THEN BEGIN
	MESSAGE, "Incompatible units (independent variables)"
	RETURN, 0
ENDIF
IF ((self.Grid1.mnemonic NE arg.Grid1.mnemonic) OR (self.Grid2.mnemonic NE arg.Grid2.mnemonic)) THEN BEGIN
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


IF PTR_VALID(self.ErrorBars) THEN BEGIN		; Interpolate error bars...
	error1 = self.ErrorBars				;	(probably we could think
	errors = INTERPOLATE(error1,sx,sy, /GRID)	;	 of something a bit more
	PTR_FREE, result.ErrorBars				;	 rigorous?)
	result.ErrorBars = PTR_NEW(errors)
ENDIF
RETURN, result
END ; ***** GKVs2D::Interpolate ***** ;


FUNCTION GKVs3D::Interpolate, arg
;
; Interpolate values from self onto the signal window of arg's grid.
;
IF (OBJ_ISA(arg, 'GKVs3D') NE 1) THEN BEGIN
	MESSAGE, "Argument is not a valid GKVs3D object"
	RETURN, 0
ENDIF
;
; Check for common units, independent variable, interval
;
IF (	(self.Grid1.units NE arg.Grid1.units) OR 	$
	(self.Grid2.units NE arg.Grid2.units) OR	$
	(self.Grid3.units NE arg.Grid3.units) ) 	THEN BEGIN
	MESSAGE, "Incompatible units (independent variables)"
	RETURN, 0
ENDIF
IF (	(self.Grid1.mnemonic NE arg.Grid1.mnemonic) OR 	$
	(self.Grid2.mnemonic NE arg.Grid2.mnemonic) OR	$
	(self.Grid3.mnemonic NE arg.Grid3.mnemonic)) 	THEN BEGIN
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
values1 = *self.values					; Get values of to be interpolated
imin = arg.Grid1.irange[0]				; Get range of indices on target grid
imax = arg.Grid1.irange[1]
jmin = arg.Grid2.irange[0]
jmax = arg.Grid2.irange[1]
kmin = arg.Grid3.irange[0]
kmax = arg.Grid3.irange[1]
x = x2[imin:imax]					; axis 1 for target grid
y = y2[jmin:jmax]					; axis 2 for target grid
z = z2[kmin:kmax]					; axis 3 for target grid
;
; Get index into old grid arrays (that is, those of 'self') for each element of the new grid arrays (that is, those of 'arg').
;
jjx = VALUE_LOCATE(x1,x)
jjy = VALUE_LOCATE(y1,y)
jjz = VALUE_LOCATE(z1,z)
;
; jjx, jjy, jjz will be = -1 if new grid point does not lie within range of old grid points.
; First create jx, jy, jz arrays (with 
jx = jjx > 0
jy = jjy > 0
jz = jjz > 0
;
; Compute old grid spacing
;
info = SIZE(values1)
nx1 = info[1] - 1
ny1 = info[2] - 1
nz1 = info[3] - 1
dx1 = x1[1:nx1] - x1[0:(nx1-1)]
dy1 = y1[1:ny1] - y1[0:(ny1-1)]
dz1 = z1[1:nz1] - z1[0:(nz1-1)]
;
; Compute fractional "address" of new x and y grid points within old grid
;
sx = jx + ((x - x1[jx])/dx1[jx])*(jjx GT 0)		; the factor (jjx GT 0) has the effect of zeroing out
sy = jy + ((y  -y1[jy])/dy1[jy])*(jjy GT 0)		;  this term for points outside of target grid.
sz = jz + ((z  -z1[jz])/dz1[jz])*(jjz GT 0)

values = INTERPOLATE(values1,sx,sy,sz, /GRID)	; Interpolate input values to target grid

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


IF PTR_VALID(self.ErrorBars) THEN BEGIN		; Interpolate error bars...
	error1 = self.ErrorBars				;	(probably we could think
	errors = INTERPOLATE(error1,sx,sy,sz, /GRID);	 of something a bit more
	PTR_FREE, result.ErrorBars				;	 rigorous?)
	result.ErrorBars = PTR_NEW(errors)
ENDIF
RETURN, result
END ; ***** GKVs3D::Interpolate ***** ;

