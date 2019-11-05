
FUNCTION GKVs4D::FFT4, Inverse=inverse
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
; Written by E. Wang
;  5-6-09
; Only change from FFT4 to FFT3 is all axis calls are to Grid 4 etc.
;
direction=-1
inverse = KEYWORD_SET(inverse)
IF inverse THEN direction=1
uniform = self.Grid4.uniform
IF(NOT uniform) THEN self -> ScaleAxis, 4, /Uniform
boundary=self.Grid4.boundary
irangeIn = self.Grid4.irange
nzIn = N_ELEMENTS(*self.Grid4.values)
CASE boundary OF
	'periodic (open)'	:	irange=[0,ntIn-1]
	'periodic (closed)'	:	irange=[0,ntIn-2]
	'periodic'		:	irange=[0,ntIn-2]
ELSE:	irange=irangeIn
ENDCASE

self -> Set, axis=4, irange=irange
;values=*(self -> GetValues())
valueptr = self -> getvalues()
values = *valueptr
PTR_FREE, valueptr
self -> Set, axis=4, irange=irangeIn
info=SIZE(values)
nx = info[1]
ny = info[2]
nz = info[3]
nt = info[4]
ntEven = 2*(nt/2)-nt+1	; will be =1 if nz is even, =0 if nz is odd
ntOdd = 1 - ntEven
ikt = nt
;
; Shift k-values if this is an inverse transform
;
ktShift = (ikt-1)/2 + 1
IF inverse THEN values = SHIFT(values, 0, 0, 0, ktShift)
;
; Make array to hold k-space representation of data
;
kSpaceValues = COMPLEXARR(nx, ny, nz, ikt)
;
; Begin loop over y-slices and time-slices
;
FOR ix = 0, nx-1 DO BEGIN	; Fourier transform
	FOR iy = 0, ny-1 DO BEGIN
            FOR iz = 0, nz-1 DO BEGIN
                kSpaceValues[ix,iy,iz,*] = FFT(values[ix,iy,iz,*], direction)
            ENDFOR
        ENDFOR
ENDFOR
;
; Shift 'kSpaceValues' such that negative
; values of kz preceed positive values of kz
; if this is a foward transform.
;
IF( NOT INVERSE) THEN kSpaceValues = SHIFT(kSpaceValues, 0, 0, 0,-ktShift)
;
; Add boundary layer to 'close' transform data
; if even number of data points
;
IF(ntEven) THEN BEGIN
	temp = COMPLEXARR(nx, ny, nz, ikt+ntEven)
	temp[*,*,*,0:(nt-1)] = kSpaceValues
	temp[*,*,*,nt] = temp[*,*,*,0]
	ikt=ikt+ntEven
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
tGrid = *self.Grid4.Values
l_t = tGrid[nt-1] - tGrid[0]
dt = tGrid[1] - tGrid[0]
CASE boundary OF
	'periodic (open)'	:	l_t=l_t + dt
	'periodic (closed)'	:	l_t=l_t + dt
	'periodic'		:	l_t=l_t + dt
ELSE:	l_t=l_t
ENDCASE
dkt = 2.*!PI/l_t
tMnemonic = STRTRIM(self.Grid4.mnemonic, 2)
tTitle=self.grid4.title
tUnits=self.grid4.units
IF(STRMATCH(tUnits, 'dimensionless')) THEN tUnits=''
new_tUnits = GKV_InvertString(tUnits)
result.Grid4.units = new_tUnits
IF inverse THEN BEGIN
	ktGrid = dkt*FINDGEN(ikt)
	;
	; Create x-space mnemonic by deleting the
	; 'k_' which should preface the configuration-space 
	; mnemonic (and remove any embedded blanks).
	; 
CASE STRLOWCASE(tMnemonic) OF
		'omega':	BEGIN
			result.grid4.mnemonic = 't'
			result.grid4.title = 't'
				END
		'n'	:	BEGIN
			result.grid4.mnemonic = 'zeta'
			result.grid4.title = '!4f!X'
				END
		ELSE	: 	BEGIN
			mnemonicSubStrings = STRSPLIT(tMnemonic, 'k_', /EXTRACT)
			result.Grid4.mnemonic = STRJOIN(mnemonicSubStrings, '_')
			titleChars = STRLEN(tTitle)
			newLength = titleChars - 5
			result.Grid4.title = STRMID(tTitle, 3, newLength)  
		END
	ENDCASE	
ENDIF ELSE BEGIN
	ktGrid = dkt*(FINDGEN(ikt) - (ikt-1)/2)
	;
	; Create k-space mnemonic by prefacing the
	; configuration-space mnemonic with  a 'k_'
	; after removing any embedded blanks.
	; 
	IF(tMnemonic EQ 't') THEN BEGIN
		result.grid4.mnemonic = 'omega'
		result.grid4.title = '!4x!X'
	ENDIF ELSE BEGIN
		result.Grid4.mnemonic = 'k_' + tMnemonic
		result.Grid4.title = 'k!D' + tTitle + '!N'
	ENDELSE
ENDELSE
PTR_FREE, result.Grid4.values
result.Grid4.values = PTR_NEW(ktGrid)
result.Grid4.irange = [0,ikt-1]
result.Grid4.range = [ktGrid[0], ktGrid[ikt-1]]
CASE ntEven OF
	0	:	result.grid4.boundary = 'periodic (open)'
	1	:	result.grid4.boundary = 'periodic (closed)'
ENDCASE

IF(KEYWORD_SET(offset)) THEN result -> ScaleAxis, 4, const=1., offset=offset

RETURN, result
END ; ****** GKVs4D::FFT4 ****** ;


