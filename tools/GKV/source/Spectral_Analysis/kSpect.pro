FUNCTION GKVs3D::kSpect, _Extra=Extra
;
; Compute k-spectrum for 3-D data (assuming that 
; third independent variable is time) by sequence of
; 2-D FFT's.  This uses less memory, and so may be
; better with large data sets.
;
; Contentx of "Extra" are passed to xSpect.
;
; Written by W.M. Nevins
;  12/8/03
;
tGrid = self.grid3
irange=self.grid3.irange
imin = irange[0]
imax = irange[1]
tValues = (*(tGrid.values))[imin:imax]
nt = N_ELEMENTS(tValues)

keyWords = 0
IF(typeOF(Extra) EQ 8) THEN BEGIN
	nExtra = Extra
	keyWords = 1
ENDIF

firstSlice = self -> slice(axis=3, index=0)
result = firstSlice -> xSpect(_Extra=Extra)
firstSlice -> Trash

IF(nt EQ 1) THEN RETURN, result

values = *result.values
FOR i=1,nt-1 DO BEGIN
	thisSlice = self -> slice(axis=3, index=i)
	IF(keyWords) THEN Extra = nExtra
	thisSpect = thisSlice -> xSpect(_Extra=Extra)
	thisSlice -> Trash
	values = values + *(thisSpect.values)
	thisSpect -> Trash
ENDFOR
values = FLOAT(values/nt)

PTR_FREE, result.values
result.values = PTR_NEW(values)
vMin = Min(values, Max=vMax)
result.vrange = [vMin, vMax]
result.indices = PTR_NEW(['*','*'])
RETURN, result
END  ;  ****** GKVs3D::kSpect ******  ;
