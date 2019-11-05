FUNCTION GKVs3D::k_max, debug=debug
;
; Assumes 'self' is a GKVs3D object (2 space plus time) for which the
; last two independent variables are ignorable (e.g., toroidal angle plus time);
; while the first independent variable (e.g., radius) is non-ignorable.
;
; For each value of first independent variable, computes mean wave number of
; corresponding to the second independent variable.  
;
; Returns a GKVs1D object whose values are the average wavenumber corresponding
; to the second (e.g., toroidal) independent variable vs. the first (e.g., radial)
; independent variable.  Error bars are the RMS width of the k-spectrum about the mean.
;
; Written by W.M. Nevins
;	6/20/00
;
result = {GKVs1D}						; Create GKVs1D structure
FOR i=0, N_TAGS(result) - 1 DO 	$
	result.(i) = self.(i)					; Copy fields from 'self' to 'result'
result.mnemonic = 'k_max'
result.title = 'k!Imax!N{' + self.title + '}'
result.units = '1/' + self.Grid2.Units
IF( STRTRIM(result.units, 2) EQ '1/') THEN result.units=' '
indices = ['*']							
result.indices = PTR_NEW(indices)
irange = self.Grid1.irange
;
; Check boundary condition
;
bc = STRTRIM(self.Grid1.boundary, 2)			; Get bc with both leading and trailing blanks removed
IF( NOT STRCMP(bc, 'periodic', 8, /FOLD_CASE) ) THEN BEGIN
	IF(irange[0] EQ 0) THEN irange[0]=1		; Avoid end-points if bc is not 'periodic' as spectral density here
	irange[1] = irange[1]-1				; tends to be dominated by boundary conditions
ENDIF
imax = irange[1]-irange[0] 
result.Grid1.values = PTR_NEW((*self.Grid1.values)[irange[0]:irange[1]])
result.Grid1.irange = [0,imax]
vmin = (*self.Grid1.values)[irange[0]]
vmax = (*self.Grid1.values)[irange[1]]
result.Grid1.range = [vmin, vmax]
k_max = FLTARR(imax+1)
k_rms = FLTARR(imax+1)
FOR i=0, imax DO BEGIN						; loop over 1st (radial) index
	index = irange[0]+i
	tempObj = self -> Slice(axis=1, index=index)	; Take slice at constant radius
	tempObj -> Norm							; Normalize 
	tempSpect = tempObj -> xspect()				; Form spectral density
	tempSpectK = tempSpect -> Avg(axis=2)			; Average over frequency
	tempSpectK -> Get, axis=1, range = krange
	krange[0] = 0.
	k_Moments = tempSpectK -> Moments(axis=1, Range=krange, /Avg)
	k_max[i] = k_Moments[1] -> GetValues()
	k_rms[i] = SQRT( k_Moments[2] -> GetValues() )
;
; Final Duties, and then clean-up
;
	IF(KEYWORD_SET(debug)) THEN print, index, (*self.Grid1.values)[index], k_max[i], k_rms
	tempObj -> Trash
	tempSpect -> Trash
	tempSpectK -> Trash
	HEAP_GC, Verbose=debug
ENDFOR
vmin = MIN(k_max - k_rms) > 0.0
vmax = MAX(k_max + k_rms)
result.vrange = [vmin, vmax]
result.values = PTR_NEW(k_max)
result.errorBars = PTR_NEW(k_rms)
resultObj = OBJ_NEW('GKVs1D', result)
return, resultObj
END ; ****** k_max ****** ;
 