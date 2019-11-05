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

FUNCTION GKVs2D::XtoK
;
; Transforms object from configuration space representation
; to fourier space representation.
;


; **** Add check for uniform grid! ****
FORWARD_FUNCTION GKVsd_MIN
valuePtr = self -> GetValues()			; Gets pointer to values within 
values=*valuePtr				;	signal window.
info=SIZE(values)
IF(info[0] ne 2) then begin
	MESSAGE, 'only 2-D for now', /Informational
	RETURN, 0
ENDIF
nx  = info[1]
ikx = nx
ny  = info[2]
iky = ny
ikyplus = (iky+1)/2
;
; Make array to hold k-space representation of data
; (implicitly using reality condition such that only
;  positive values of ky need to be saved)
;
kSpaceValues = COMPLEXARR(ikx, ikyplus)
;
; Fourier transform
;
temp = FFT(values, -1)
;
; load positive values of ky from 'temp' into 
; 'kSpaceValues' (the output array)
;
kSpaceValues = temp[*, 0:(ikyplus-1)]
;
; Shift 'kSpaceValues' such that negative
; values of kx preceed positive values of kx.
;
kxShift = (ikx/2) + 1
kSpaceValues = SHIFT(kSpaceValues,  -kxShift, 0)
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
xGrid = *self.Grid1.Values
imin = self.Grid1.irange[0]
imax = self.Grid1.irange[1]
lx = xGrid[imax] - xGrid[imin]
dx = xGrid[1] - xgrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dx larger
;
IF(STRCMP(self.grid1.boundary,"periodic (open)", /FOLD_CASE)) THEN lx = lx + dx
;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx/2))
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
xMnemonic = STRTRIM(self.Grid1.mnemonic, 2)
result.Grid1.mnemonic = 'k' + xMnemonic
;
; Make Grid1 title equal to 'k' sub 'xMnemonic'
;
result.Grid1.title = 'k!I' + xMnemonic + '!N'
;
;
; Make Grid1 units equal to 1/xUnits
;
result.Grid1.Units = "1/" + self.Grid1.units
result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]
yGrid = *self.Grid2.values
jmin = self.Grid2.irange[0]
jmax = self.Grid2.irange[1]
ly = yGrid[jmax] - yGrid[jmin]
dy = yGrid[1] - yGrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dx larger
;
IF(STRCMP(self.grid2.boundary,"periodic (open)", /FOLD_CASE)) THEN ly = ly + dy
dky = 2*!PI/ly
kyGrid = dky*FINDGEN(ikyplus)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(kyGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
yMnemonic = STRTRIM(self.Grid2.mnemonic, 2)
result.Grid2.mnemonic = 'k' + yMnemonic
;
; Make Grid2 title equal to 'k' sub 'yMnemonic'
;
result.Grid2.title = 'k!I' + yMnemonic + '!N'
;
;
; Make Grid2 units equal to 1/yUnits
;
result.Grid2.Units = "1/" + self.Grid2.units
result.Grid2.range = [0., kyGrid[ikyplus-1]]
result.Grid2.irange = [0,ikyplus-1]
RETURN, result
END ; ****** GKVs2D::XtoK ****** ;


FUNCTION GKVs2D::XtToKt
;
; Transforms object from configuration space representation
; to fourier space representation.  Assumes that the first 
; dimension correspond to an indepent varialbes in configuration
; space (e.g., 'x'), while the second dimension
; is time (which we will not transform over).
;
ndims = self -> Numdims()
IF(ndims ne 2) then begin
	MESSAGE, 'Works only in 2-D at present', /Informational
	RETURN, 0
ENDIF
boundary=self.grid1.boundary
irangeIn = self.grid1.irange
nxIn = N_ELEMENTS(*self.grid1.values)
CASE boundary OF
	'periodic (open)'	:	irange=[0,nxIn-1]
	'periodic (closed)'	:	irange=[0,nxIn-2]
	'periodic'		:	irange=[0,nxIn-1]
ELSE:	irange=irangeIn
ENDCASE

self -> Set, axis=1, irange=irange
values=*(self -> GetValues())
self -> Set, axis=1, irange=irangeIn
info=SIZE(values)
nx  = info[1]
nt = info[2]
;
; Make array to hold k-space representation of data
;
nxEven = 2*(nx/2)-nx+1	; will be =1 if nx is even, =0 if nx is odd
ikx = nx
kSpaceValues = COMPLEXARR(ikx, nt)
;
; Begin loop over time-slices
;
FOR it = 0, nt-1 DO BEGIN
;
; Fourier transform
;
	temp = FFT(values[*,it], -1)
;
; load positive values of ky from 'temp' into 
; 'kSpaceValues' (the output array)
;
	kSpaceValues[*,it] = temp
;
ENDFOR
;
; Shift 'kSpaceValues' such that negative
; values of kx preceed positive values of kx.
;
kxShift = (ikx-1)/2 + 1
kSpaceValues = SHIFT(kSpaceValues,  -kxShift, 0)
IF(nxEven) THEN BEGIN
	temp = COMPLEXARR(ikx+nxEven, nt)
	temp[0:(nx-1),*] = kSpaceValues
	temp[nx,*] = temp[0,*]
	ikx=ikx+nxEven
ENDIF
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
xGrid = *self.Grid1.Values
lx = xGrid[nx-1] - xGrid[0]
dx = xGrid[1] - xgrid[0]
CASE boundary OF
	'periodic (open)'	:	lx=lx+dx
	'periodic (closed)'	:	lx=lx
	'periodic'		:	lx=lx+dx
ELSE:	lx=lx
ENDCASE
;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx-1)/2)
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
result.Grid1.irange = [0,ikx-1]
;
; Create k-space mnemonic by prefacing the
; configuration-space mnemonic with  a 'k_'
; after removing any embedded blanks.
; 
xMnemonic = STRTRIM(self.Grid1.mnemonic, 2)
result.Grid1.mnemonic = 'k_' + xMnemonic
;
; Make Grid1 title equal to 'k' sub 'xTitle'
;
xtitle=self.grid1.title
result.Grid1.title = 'k!D' + xTitle + '!N'
;
;
; Make Grid1 units equal to 1/xUnits
;
result.Grid1.Units = "1/" + self.Grid1.units

result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]

RETURN, result
END ; ****** GKVs2D::XtToKt ****** ;




FUNCTION GKVs3D::XtoK
;
; Transforms object from configuration space representation
; to fourier space representation.  Assumes that the first 
; two dimensions correspond to indepent varialbes in configuration
; space (e.g., 'x' and 'y'), while the third dimension
; is time (which we will not transform over).
;
values=*self.values
info=SIZE(values)
IF(info[0] ne 3) then begin
	MESSAGE, 'only up to 3-D for now', /Informational
	RETURN, 0
ENDIF
nx  = info[1]
ikx = nx
ny  = info[2]
iky = ny
ikyplus = (iky+1)/2
nt = info[3]
;
; Make array to hold k-space representation of data
; (implicitly using reality condition such that only
;  positive values of ky need to be saved)
;
kSpaceValues = COMPLEXARR(ikx, ikyplus, nt)
;
; Begin loop over time-slices
;
FOR it = 0, nt-1 DO BEGIN
;
; Fourier transform
;
	temp = FFT(values[*,*,it], -1)
;
; load positive values of ky from 'temp' into 
; 'kSpaceValues' (the output array)
;
	kSpaceValues[*,*,it] = temp[*, 0:(ikyplus-1)]
;
ENDFOR
;
; Shift 'kSpaceValues' such that negative
; values of kx preceed positive values of kx.
;
kxShift = (ikx/2) + 1
kSpaceValues = SHIFT(kSpaceValues,  -kxShift, 0, 0)
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
xGrid = *self.Grid1.Values
lx = xGrid[nx-1] - xGrid[0]
dx = xGrid[1] - xgrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dx larger
;
IF(STRCMP(self.grid1.boundary,"periodic (open)", /FOLD_CASE)) THEN lx = lx + dx
;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx/2))
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
xMnemonic = STRTRIM(self.Grid1.mnemonic, 2)
result.Grid1.mnemonic = 'k' + xMnemonic
;
; Make Grid1 title equal to 'k' sub 'xMnemonic'
;
result.Grid1.title = 'k!I' + xMnemonic + '!N'
;
; Make Grid1 units equal to 1/ xUnits
;
result.Grid1.Units = '1/' + self.Grid1.units
result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]
;
yGrid = *self.Grid2.values
ly = yGrid[ny-1] - yGrid[0]
dy = yGrid[1] - yGrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dx larger
;
IF(STRCMP(self.grid2.boundary,"periodic (open)", /FOLD_CASE)) THEN ly = ly + dy
dky = 2*!PI/ly
kyGrid = dky*FINDGEN(ikyplus)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(kyGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
yMnemonic = STRTRIM(self.Grid2.mnemonic, 2)
result.Grid2.mnemonic = 'k' + yMnemonic
;
; Make Grid2 title equal to 'k' sub 'yMnemonic'
;
result.Grid2.title = 'k!I' + yMnemonic + '!N'
;
; Make Grid2 units equal to 1/ yUnits
;
result.Grid2.Units = '1/' + self.Grid2.units
result.Grid2.range = [0., kyGrid[ikyplus-1]]
result.Grid2.irange = [0,ikyplus-1]
RETURN, result
END ; ****** GKVs3D::XtoK ****** ;


FUNCTION GKVs2D::KtoX
;
; Transforms object from fourier space representation to 
; configuration space representation.
;
values=*self.values
info=SIZE(values)
IF(info[0] ne 2) then begin
	MESSAGE, 'only 2-D for now', /Informational
	RETURN, 0
ENDIF
ikx=info[1]
ikyplus=info[2]
iky=2*info[2]-1
;
; Find index of kx=0
;
kxGrid=*self.grid1.values
eps=MIN(kxGrid^2, ikxshift)
;
; Shift 'values' such that index corresponding to kx=0 is first index
;
values = SHIFT(values, -ikxshift, 0)
kxGrid = SHIFT(kxGrid, -ikxshift)
;
; Make array to hold full k-space representation of data
; and populate 1st half with data from values
;
kSpaceValues = COMPLEXARR(ikx, iky)
;
; load positive values of ky, then construct negative
; ky elements using reality condition.
;
kSpaceValues[0:ikx-1, 0:ikyplus-1] = values
kSpaceValues[1:ikx-1, ikyplus:(iky-1)] = CONJ(REVERSE(REVERSE(values[1:ikx-1, 1:ikyplus-1], 1), 2))
kSpaceValues[0,ikyplus:(iky-1)] = CONJ(REVERSE(values[0,1:ikyplus-1], 2))
;
; Fourier transform
;
xSpaceValues = FLOAT(FFT(kSpaceValues, 1))
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(xSpaceValues)
vmin = GKVsd_MIN(xSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
lx = 2*!PI/kxGrid[1]
nx=ikx
dx = lx/(nx)
xGrid = dx*FINDGEN(nx)
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(xGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid1.mnemonic, 'kN_') 
xMnemonic = STRMID(self.Grid1.mnemonic, firstChars[0]) 
result.Grid1.mnemonic = STRTRIM(xMnemonic, 2)
result.Grid1.title = xMnemonic
result.Grid1.boundary = 'periodic (open)'
result.Grid1.range = [0, lx]
result.Grid1.irange = [0,nx-1]
kyGrid = *self.Grid2.values
ny=iky
ly = 2*!PI/kyGrid[1]
dy = ly/(ny)
yGrid = dy*FINDGEN(ny)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(yGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid2.mnemonic, 'kN_') 
yMnemonic = STRMID(self.Grid2.mnemonic, firstChars[0]) 
result.Grid2.mnemonic = STRTRIM(yMnemonic, 2)
result.Grid2.title = yMnemonic
result.Grid2.boundary = 'periodic (open)'
result.Grid2.range = [0., ly]
result.Grid2.irange = [0,ny-1]
RETURN, result
END ; ****** GKVs2D::KtoX ****** ;


FUNCTION GKVs3D::KtoX
;
; Transforms object from fourier space representation to configuration space representation
; Assumes that the first two dimensions correspond to k-values, while the third dimension
; is time (which we will not transform over.
;
values=*self.values
info=SIZE(values)
IF(info[0] ne 3) then begin
	MESSAGE, 'only up to 3-D for now', /Informational
	RETURN, 0
ENDIF
ikx=info[1]
ikyplus=info[2]
iky=2*info[2]-1
nt = info[3]
;
; Find index of kx=0
;
kxGrid=*self.grid1.values
eps=MIN(kxGrid^2, ikxshift)
;
; Shift 'values' such that index corresponding to kx=0 is first index
;
values = SHIFT(values, -ikxshift, 0, 0)
kxGrid = SHIFT(kxGrid, -ikxshift)
;
; Make array to hold full k-space representation of data
; at one time slice, and array to hold resulting xSpaceValues
; at all times available.
;
kSpaceValues = COMPLEXARR(ikx, iky)
xSpaceValues = FLTARR(ikx, iky, nt)
;
; Begin loop over time-slices
;
FOR it = 0, nt-1 DO BEGIN
;
; load positive values of ky, then construct negative
; ky elements using reality condition.
;
	kSpaceValues[0:ikx-1, 0:ikyplus-1] = values[*,*,it]
	kSpaceValues[1:ikx-1, ikyplus:(iky-1)] = CONJ(REVERSE(REVERSE(values[1:ikx-1, 1:ikyplus-1, it], 1), 2))
	kSpaceValues[0,ikyplus:(iky-1)] = CONJ(REVERSE(values[0,1:ikyplus-1, it], 2))
;
; Fourier transform this time-slice from k-space to x-space
;
	xSpaceValues[*,*,it] = FLOAT(FFT(kSpaceValues, 1))
;
ENDFOR
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(xSpaceValues)
vmin = GKVsd_MIN(xSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
; Load up x-grid
;
lx = 2*!PI/kxGrid[1]
nx=ikx
dx = lx/(nx)
xGrid = dx*FINDGEN(nx)
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(xGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid1.mnemonic, 'kN_') 
xMnemonic = STRMID(self.Grid1.mnemonic, firstChars[0]) 
result.Grid1.mnemonic = STRTRIM(xMnemonic, 2)
result.Grid1.title = xMnemonic
result.Grid1.uniform = 1b
result.Grid1.boundary = 'periodic (open)'
result.Grid1.range = [0, lx]
result.Grid1.irange = [0,nx-1]
;
; Extract units by removing leading "1/"
;
lastChars = STRSPLIT(self.Grid1.units, '/', /EXTRACT)
nChars = N_ELEMENTS(lastChars)
xUnits = lastChars[nChars-1]
result.Grid1.units = STRTRIM(xUnits, 2)
;
; Now, load up y-grid
;
kyGrid = *self.Grid2.values
ny=iky
ly = 2*!PI/kyGrid[1]
dy = ly/(ny)
yGrid = dy*FINDGEN(ny)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(yGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid2.mnemonic, 'kN_') 
yMnemonic = STRMID(self.Grid2.mnemonic, firstChars[0]) 
result.Grid2.mnemonic = STRTRIM(yMnemonic, 2)
result.Grid2.title = yMnemonic
result.Grid1.uniform = 1b
result.Grid2.boundary = 'periodic (open)'
result.Grid2.range = [0., ly]
result.Grid2.irange = [0,ny-1]
;
; Extract units by removing leading "1/" 
;
lastChars = STRSPLIT(self.Grid2.units, '/', /EXTRACT)
nChars = N_ELEMENTS(lastChars)
yUnits = lastChars[nChars-1]
result.Grid2.units = STRTRIM(yUnits, 2)
RETURN, result
END ; ****** GKVs3D::KtoX ****** ;


FUNCTION GKVs4D::XtoK
;
; Transforms object from configuration space representation
; to fourier space representation.  Assumes that the first 
; two dimensions correspond to indepent varialbes in configuration
; space (e.g., 'x' and 'y'), while the third dimension
; is time (which we will not transform over).
;
values=*self.values
info=SIZE(values)
IF(info[0] ne 4) then begin
	MESSAGE, 'only up to 4-D for now', /Informational
	RETURN, 0
ENDIF
nx  = info[1]
ikx = nx
ny  = info[2]
iky = ny
ikyplus = (iky+1)/2
nz = info[3]
nt = info[4]
;
; Make array to hold k-space representation of data
; (implicitly using reality condition such that only
;  positive values of ky need to be saved)
;
kSpaceValues = COMPLEXARR(ikx, ikyplus, nz, nt)
;
; Begin loop over z-slices
;
FOR iz = 0, nz-1 DO BEGIN
;
; Begin loop over time-slices
;
FOR it = 0, nt-1 DO BEGIN
;
; Fourier transform
;
	temp = FFT(values[*,*,iz,it], -1)
;
; load positive values of ky from 'temp' into 
; 'kSpaceValues' (the output array)
;
	kSpaceValues[*,*,iz,it] = temp[*, 0:(ikyplus-1)]
;
ENDFOR  ; end loop over z-slices
ENDFOR  ; end loop over t-slices
;
; Shift 'kSpaceValues' such that negative
; values of kx preceed positive values of kx.
;
kxShift = (ikx/2) + 1
kSpaceValues = SHIFT(kSpaceValues,  -kxShift, 0, 0, 0)
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
xGrid = *self.Grid1.Values
lx = xGrid[nx-1] - xGrid[0]
dx = xGrid[1] - xgrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dx larger
;
IF(STRCMP(self.grid1.boundary,"periodic (open)", /FOLD_CASE)) THEN lx = lx + dx

;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx/2))
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
xMnemonic = STRTRIM(self.Grid1.mnemonic, 2)
result.Grid1.mnemonic = 'k' + xMnemonic
;
; Make Grid1 title equal to 'k' sub 'xMnemonic'
;
result.Grid1.title = 'k!I' + xMnemonic + '!N'
;
; Make Grid1 units equal to 1/ xUnits
;
result.Grid1.Units = '1/' + self.Grid1.units
result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]
;
yGrid = *self.Grid2.values
ly = yGrid[ny-1] - yGrid[0]
dy = yGrid[1] - yGrid[0]
;
; If boundary condition is "periodic (open),
; then the system length is dy larger
;
IF(STRCMP(self.grid2.boundary,"periodic (open)", /FOLD_CASE)) THEN ly = ly + dy
dky = 2*!PI/ly
kyGrid = dky*FINDGEN(ikyplus)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(kyGrid)
;
; Create k-space mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
yMnemonic = STRTRIM(self.Grid2.mnemonic, 2)
result.Grid2.mnemonic = 'k' + yMnemonic
;
; Make Grid2 title equal to 'k' sub 'yMnemonic'
;
result.Grid2.title = 'k!I' + yMnemonic + '!N'
;
; Make Grid2 units equal to 1/ yUnits
;
result.Grid2.Units = '1/' + self.Grid2.units
result.Grid2.range = [0., kyGrid[ikyplus-1]]
result.Grid2.irange = [0,ikyplus-1]
RETURN, result
END ; ****** GKVs4D::XtoK ****** ;


FUNCTION GKVs4D::KtoX
;
; Transforms object from fourier space representation to configuration space representation
; Assumes that the first two dimensions correspond to transverse k-values, while the third 
; and fourth dimensions are distance-along-B (theta) and time (which we will not transform over).
;
values=*self.values
info=SIZE(values)
IF(info[0] ne 4) then begin
	MESSAGE, 'only up to 4-D for now', /Informational
	RETURN, 0
ENDIF
ikx=info[1]
ikyplus=info[2]
iky=2*info[2]-1
nz = info[3]
nt = info[4]
;
; Find index of kx=0
;
kxGrid=*self.grid1.values
eps=MIN(kxGrid^2, ikxshift)
;
; Shift 'values' such that index corresponding to kx=0 is first index
;
values = SHIFT(values, -ikxshift, 0, 0, 0)
kxGrid = SHIFT(kxGrid, -ikxshift)
;
; Make array to hold full k-space representation of data
; at one (z,t) slice, and array to hold resulting xSpaceValues
; at all (z,t) available.
;
kSpaceValues = COMPLEXARR(ikx, iky)
xSpaceValues = FLTARR(ikx, iky, nz, nt)
; begin loop over z-slices
;
FOR iz = 0, nz-1 DO BEGIN
;
; Begin loop over time-slices
;
FOR it = 0, nt-1 DO BEGIN
;
; load positive values of ky, then construct negative
; ky elements using reality condition.
;
	kSpaceValues[0:ikx-1, 0:ikyplus-1] = values[*,*,iz,it]
	kSpaceValues[1:ikx-1, ikyplus:(iky-1)] = CONJ(REVERSE(REVERSE(values[1:ikx-1, 1:ikyplus-1, iz, it], 1), 2))
	kSpaceValues[0,ikyplus:(iky-1)] = CONJ(REVERSE(values[0,1:ikyplus-1, iz, it], 2))
;
; Fourier transform this time-slice from k-space to x-space
;
	xSpaceValues[*,*,iz,it] = FLOAT(FFT(kSpaceValues, 1))
;
ENDFOR  ; end loop over z-slices
ENDFOR  ; end loop over t-slices
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues)
result.values = PTR_NEW(xSpaceValues)
vmin = GKVsd_MIN(xSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
; Load up x-grid
;
lx = 2*!PI/kxGrid[1]
nx=ikx
dx = lx/(nx)
xGrid = dx*FINDGEN(nx)
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(xGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid1.mnemonic, 'kN_') 
xMnemonic = STRMID(self.Grid1.mnemonic, firstChars[0]) 
result.Grid1.mnemonic = STRTRIM(xMnemonic, 2)
result.Grid1.title = xMnemonic
result.Grid1.uniform = 1b
result.Grid1.boundary = 'periodic (open)'
result.Grid1.range = [0, lx]
result.Grid1.irange = [0,nx-1]
;
; Extract units by removing leading "1/"
;
lastChars = STRSPLIT(self.Grid1.units, '/', /EXTRACT)
nChars = N_ELEMENTS(lastChars)
xUnits = lastChars[nChars-1]
result.Grid1.units = STRTRIM(xUnits, 2)
;
; Now, load up y-grid
;
kyGrid = *self.Grid2.values
ny=iky
ly = 2*!PI/kyGrid[1]
dy = ly/(ny)
yGrid = dy*FINDGEN(ny)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(yGrid)
;
; Extract configuration-space mnemonic from
; k-space mnemonic by removing any leading 
; 'k', 'N', and/or '_''s.  Then remove any
; embedded blanks.
; 
firstChars = STRSPLIT(self.Grid2.mnemonic, 'kN_') 
yMnemonic = STRMID(self.Grid2.mnemonic, firstChars[0]) 
result.Grid2.mnemonic = STRTRIM(yMnemonic, 2)
result.Grid2.title = yMnemonic
result.Grid1.uniform = 1b
result.Grid2.boundary = 'periodic (open)'
result.Grid2.range = [0., ly]
result.Grid2.irange = [0,ny-1]
;
; Extract units by removing leading "1/" 
;
lastChars = STRSPLIT(self.Grid2.units, '/', /EXTRACT)
nChars = N_ELEMENTS(lastChars)
yUnits = lastChars[nChars-1]
result.Grid2.units = STRTRIM(yUnits, 2)
RETURN, result
END ; ****** GKVs4D::KtoX ****** ;
