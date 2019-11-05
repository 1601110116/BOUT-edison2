FUNCTION GKVs2D::XtoK
;
; Transforms object from configuration space representation
; to fourier space representation.
;
values=*self.values
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
lx = xGrid[nx-1] - xGrid[0]
dx = xGrid[1] - xgrid[0]
;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx/2))
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
;
; Create k-sace mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
xMnemonic = STRTRIM(self.Grid1.mnemonic, 2)
result.Grid1.mnemonic = 'k' + xMnemonic
;
; Make Grid1 title equal to 'k' sub 'xMnemonic'
;
result.Grid1.title = 'k!I' + xMnemonic + '!N'
result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]
yGrid = *self.Grid2.values
ly = yGrid[ny-1] - yGrid[0]
dy = yGrid[1] - yGrid[0]
dky = 2*!PI/ly
kyGrid = dky*FINDGEN(ikyplus)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(kyGrid)
;
; Create k-sace mnemonic by appending a 'k'
; to the configuration-space mnemonic after
; removing any embedded blanks.
; 
yMnemonic = STRTRIM(self.Grid2.mnemonic, 2)
result.Grid1.mnemonic = 'k' + yMnemonic
;
; Make Grid2 title equal to 'k' sub 'yMnemonic'
;
result.Grid1.title = 'k!I' + yMnemonic + '!N'
result.Grid2.range = [0., kyGrid[ikyplus-1]]
result.Grid2.irange = [0,ikyplus-1]
RETURN, result
END ; ****** GKVs2D::KtoX ****** ;


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
;kxMax = !PI/dx
dkx = 2*!PI/lx
kxGrid = dkx*(FINDGEN(ikx) - (ikx/2))
PTR_FREE, result.Grid1.values
result.Grid1.values = PTR_NEW(kxGrid)
result.Grid1.mnemonic = 'kx'
result.Grid1.title = 'k!Ix!N'
result.Grid1.range = [kxGrid[0], kxGrid[ikx-1]]
result.Grid1.irange = [0,ikx-1]
;
yGrid = *self.Grid2.values
ly = yGrid[ny-1] - yGrid[0]
dy = yGrid[1] - yGrid[0]
dky = 2*!PI/ly
kyGrid = dky*FINDGEN(ikyplus)
PTR_FREE, result.Grid2.values
result.Grid2.values = PTR_NEW(kyGrid)
result.Grid2.mnemonic = 'ky'
result.Grid2.title = 'k!Iy!N'
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
dx = lx/(nx-1)
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
result.Grid1.range = [0, lx]
result.Grid1.irange = [0,nx-1]
kyGrid = *self.Grid2.values
ny=iky
ly = 2*!PI/kyGrid[1]
dy = ly/(ny-1)
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
dx = lx/(nx-1)
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
result.Grid1.range = [0, lx]
result.Grid1.irange = [0,nx-1]
;
; Now, load up y-grid
;
kyGrid = *self.Grid2.values
ny=iky
ly = 2*!PI/kyGrid[1]
dy = ly/(ny-1)
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
result.Grid2.range = [0., ly]
result.Grid2.irange = [0,ny-1]
RETURN, result
END ; ****** GKVs3D::KtoX ****** ;

