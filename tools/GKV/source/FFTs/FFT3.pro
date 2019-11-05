
FUNCTION GKVs3D::FFT3, Inverse=inverse
;
; Transforms object from configuration space representation
; in third independent variable to fourier space representation
; in third independent variable.
;
; Keywords:
;
;	inverse		Set this keyword (i.e., put "/Inverse" on command line)
;			to perform inverse transform. (Optional)
;
;	offset		Set equal to desired minimum value of dependent
;			variable on inverse transform.  Defaults to 
;			zero (optional).
;
; Written by W.M. Nevins
;  11/24/01
;
nDims = Self -> NumDims()
IF(nDims EQ 4) THEN BEGIN
	result = self -> GKVs4D::FFT3(Inverse=inverse)
	RETURN, result
END
direction=-1
inverse = KEYWORD_SET(inverse)
IF inverse THEN direction=1
uniform = self.Grid3.uniform
IF(NOT uniform) THEN self -> ScaleAxis, 3, /Uniform
boundary=self.Grid3.boundary
irangeIn = self.Grid3.irange
nzIn = N_ELEMENTS(*self.Grid3.values)
CASE boundary OF
	'periodic (open)'	:	irange=[0,nzIn-1]
	'periodic (closed)'	:	irange=[0,nzIn-2]
	'periodic'		:	irange=[0,nzIn-2]
ELSE:	irange=irangeIn
ENDCASE

self -> Set, axis=3, irange=irange
;values=*(self -> GetValues())
valueptr = self -> getvalues()
values = *valueptr
PTR_FREE, valueptr
self -> Set, axis=3, irange=irangeIn
info=SIZE(values)
nx = info[1]
ny = info[2]
nz = info[3]
nzEven = 2*(nz/2)-nz+1	; will be =1 if nz is even, =0 if nz is odd
nzOdd = 1 - nzEven
ikz = nz
;
; Shift k-values if this is an inverse transform
;
kzShift = (ikz-1)/2 + 1
IF inverse THEN values = SHIFT(values, 0, 0, kzShift)
;
; Make array to hold k-space representation of data
;
kSpaceValues = COMPLEXARR(nx, ny, ikz)
;
; Begin loop over y-slices and time-slices
;
FOR ix = 0, nx-1 DO BEGIN	; Fourier transform
	FOR iy = 0, ny-1 DO BEGIN
		kSpaceValues[ix,iy,*] = FFT(values[ix,iy,*], direction)
	ENDFOR
ENDFOR
;
; Shift 'kSpaceValues' such that negative
; values of kz preceed positive values of kz
; if this is a foward transform.
;
IF( NOT INVERSE) THEN kSpaceValues = SHIFT(kSpaceValues, 0, 0, -kzShift)
;
; Add boundary layer to 'close' transform data
; if even number of data points
;
IF(nzEven) THEN BEGIN
	temp = COMPLEXARR(nx, ny, ikz+nzEven)
	temp[*,*,0:(nz-1)] = kSpaceValues
	temp[*,*,nz] = temp[*,*,0]
	ikz=ikz+nzEven
	kSpaceValues = temp
ENDIF
;
; Correct normalization
;
;CASE inverse OF
;	0	:	kSpaceValues = kSpaceValues/iky
;	1	:	kSpaceValues = kSpaceValues*iky
;ENDCASE
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues, /NoErrorBars)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
zGrid = *self.Grid3.Values
lz = zGrid[nx-1] - zGrid[0]
dz = zGrid[1] - zGrid[0]
CASE boundary OF
	'periodic (open)'	:	lz=lz + dz
	'periodic (closed)'	:	lz=lz + dz
	'periodic'		:	lz=lz + dz
ELSE:	lz=lz
ENDCASE
dkz = 2.*!PI/lz
zMnemonic = STRTRIM(self.Grid3.mnemonic, 2)
zTitle=self.grid3.title
zUnits=self.grid3.units
IF(STRMATCH(zUnits, 'dimensionless')) THEN zUnits=''
new_zUnits = GKV_InvertString(zUnits)
result.Grid3.units = new_zUnits
IF inverse THEN BEGIN
	kzGrid = dkz*FINDGEN(ikz)
	;
	; Create x-space mnemonic by deleting the
	; 'k_' which should preface the configuration-space 
	; mnemonic (and remove any embedded blanks).
	; 
CASE STRLOWCASE(zMnemonic) OF
		'omega':	BEGIN
			result.grid3.mnemonic = 't'
			result.grid3.title = 't'
				END
		'n'	:	BEGIN
			result.grid3.mnemonic = 'zeta'
			result.grid3.title = '!4f!X'
				END
		ELSE	: 	BEGIN
			mnemonicSubStrings = STRSPLIT(xMnemonic, 'k_', /EXTRACT)
			result.Grid3.mnemonic = STRJOIN(mnemonicSubStrings, '_')
			titleChars = STRLEN(zTitle)
			newLength = titleChars - 5
			result.Grid3.title = STRMID(zTitle, 3, newLength)
		END
	ENDCASE	
ENDIF ELSE BEGIN
	kzGrid = dkz*(FINDGEN(ikz) - (ikz-1)/2)
	;
	; Create k-space mnemonic by prefacing the
	; configuration-space mnemonic with  a 'k_'
	; after removing any embedded blanks.
	; 
	IF(zMnemonic EQ 't') THEN BEGIN
		result.grid3.mnemonic = 'omega'
		result.grid3.title = '!4x!X'
	ENDIF ELSE BEGIN
		result.Grid3.mnemonic = 'k_' + zMnemonic
		result.Grid3.title = 'k!D' + zTitle + '!N'
	ENDELSE
ENDELSE
PTR_FREE, result.Grid3.values
result.Grid3.values = PTR_NEW(kzGrid)
result.Grid3.irange = [0,ikz-1]
result.Grid3.range = [kzGrid[0], kzGrid[ikz-1]]
CASE nzEven OF
	0	:	result.grid3.boundary = 'periodic (open)'
	1	:	result.grid3.boundary = 'periodic (closed)'
ENDCASE

IF(KEYWORD_SET(offset)) THEN result -> ScaleAxis, 3, const=1., offset=offset

RETURN, result
END ; ****** GKVs3D::FFT3 ****** ;


FUNCTION GKVs4D::FFT3, Inverse=inverse, offset=offset
;
; Transforms object from configuration space representation
; in second independent variable to fourier space representation
; in second independent variable.
;
; Keywords:
;
;	inverse		Set this keyword (i.e., put "/Inverse" on command line)
;			to perform inverse transform. (Optional)
;
;	offset		Set equal to desired minimum value of dependent
;			variable on inverse transform.  Defaults to 
;			zero (optional).
;
; Written by W.M. Nevins
;  11/24/01
;
direction=-1
inverse = KEYWORD_SET(inverse)
IF inverse THEN direction=1
uniform = self.Grid3.uniform
IF(NOT uniform) THEN self -> ScaleAxis, 3, /Uniform
boundary=self.Grid3.boundary
irangeIn = self.Grid3.irange
nzIn = N_ELEMENTS(*self.Grid3.values)
CASE boundary OF
	'periodic (open)'	:	irange=[0,nzIn-1]
	'periodic (closed)'	:	irange=[0,nzIn-2]
	'periodic'		:	irange=[0,nzIn-2]
ELSE:	irange=irangeIn
ENDCASE
self -> Set, axis=3, irange=irange
;values=*(self -> GetValues())
valueptr = self -> getvalues()
values = *valueptr
PTR_FREE, valueptr
self -> Set, axis=3, irange=irangeIn
info=SIZE(values)
nx = info[1]
ny = info[2]
nz = info[3]
nt = info[4]
nzEven = 2*(nz/2)-nz+1	; will be =1 if ny is even, =0 if ny is odd
nzOdd = 1 - nzEven
ikz = nz
;
; Shift k-values if this is an inverse transform
;
kzShift = (ikz-1)/2 + 1
IF inverse THEN values = SHIFT(values, 0, 0,kzShift, 0)
;
; Make array to hold k-space representation of data
;
kSpaceValues = COMPLEXARR(nx, ny,ikz, nt)
;
; Begin loop over y-slices and time-slices
;
FOR ix = 0, nx-1 DO BEGIN	; Fourier transform
    FOR iy = 0, ny-1 DO BEGIN
        FOR it = 0, nt-1 DO BEGIN
            kSpaceValues[ix,iy,*,it] = FFT(values[ix,iy,*,it], direction)
        ENDFOR
    ENDFOR
ENDFOR
;
; Shift 'kSpaceValues' such that negative
; values of ky preceed positive values of ky
; if this is a foward transform.
;
IF( NOT INVERSE) THEN kSpaceValues = SHIFT(kSpaceValues, 0, 0, -kzShift, 0)
;
; Add boundary layer to 'close' transform data
; if even number of data points
;
IF(nzEven) THEN BEGIN
	temp = COMPLEXARR(nx, ny, ikz+nzEven, nt)
	temp[*,*,0:(nz-1),*] = kSpaceValues
	temp[*,*,nz,*] = temp[*,*,0,*]
	ikz=ikz+nzEven
	kSpaceValues = temp
    ENDIF
;
; Correct normalization
;
;CASE inverse OF
;	0	:	kSpaceValues = kSpaceValues/ikz
;	1	:	kSpaceValues = kSpaceValues*ikz
;ENDCASE
;
; Make copy of 'self' to store results in
;
Result = self -> MakeCopy(/NoValues, /NoErrorBars)
result.values = PTR_NEW(kSpaceValues)
vmin = GKVsd_MIN(kSpaceValues, MAX=vmax)
result.vrange = [vmin,vmax]
;
zGrid = *self.Grid3.Values
lz = zGrid[nz-1] - zGrid[0]
dz = zGrid[1] - zGrid[0]
CASE boundary OF
	'periodic (open)'	:	lz=lz + dz
	'periodic (closed)'	:	lz=lz + dz
	'periodic'		:	lz=lz + dz
ELSE:	lz=lz
ENDCASE
dkz = 2.*!PI/lz
zMnemonic = STRTRIM(self.Grid3.mnemonic, 2)
zTitle=self.grid3.title
zUnits=self.grid3.units
IF(STRMATCH(zUnits, 'dimensionless')) THEN zUnits=''
new_zUnits = GKV_InvertString(zUnits)
result.Grid3.units = new_zUnits
IF inverse THEN BEGIN
	kzGrid = dkz*FINDGEN(ikz)
	;
	; Create x-space mnemonic by deleting the
	; 'k_' which should preface the configuration-space 
	; mnemonic (and remove any embedded blanks).
	; 
	CASE STRLOWCASE(zMnemonic) OF
		'omega':	BEGIN
			result.grid3.mnemonic = 't'
			result.grid3.title = 't'
				END
		'n'	:	BEGIN
			result.grid3.mnemonic = 'zeta'
			result.grid3.title = '!4f!X'
				END
		ELSE	: 	BEGIN
			mnemonicSubStrings = STRSPLIT(zMnemonic, 'k_', /EXTRACT)
			result.Grid3.mnemonic = STRJOIN(mnemonicSubStrings, '_')
			titleChars = STRLEN(zTitle)
			newLength = titleChars - 5
			result.Grid3.title = STRMID(zTitle, 3, newLength)
		END
	ENDCASE	
ENDIF ELSE BEGIN
	kzGrid = dkz*(FINDGEN(ikz) - (ikz-1)/2)
	;
	; Create k-space mnemonic by prefacing the
	; configuration-space mnemonic with  a 'k_'
	; after removing any embedded blanks.
	; 
	IF(zMnemonic EQ 't') THEN BEGIN
		result.grid3.mnemonic = 'omega'
		result.grid3.title = '!4x!X'
	ENDIF ELSE BEGIN
		result.Grid3.mnemonic = 'k_' + zMnemonic  ; ymemonic replaced with zmenmonic
		result.Grid3.title = 'k!D' + zTitle + '!N'
	ENDELSE
ENDELSE
PTR_FREE, result.Grid3.values
result.Grid3.values = PTR_NEW(kzGrid)
result.Grid3.irange = [0,ikz-1]
result.Grid3.range = [kzGrid[0], kzGrid[ikz-1]]
CASE nzEven OF
	0	:	result.grid3.boundary = 'periodic (open)'
	1	:	result.grid3.boundary = 'periodic (closed)'
ENDCASE

IF(KEYWORD_SET(offset)) THEN result -> ScaleAxis, 3, const=1., offset=offset

RETURN, result
END ; ****** GKVs4D::FFT3 ****** ;
