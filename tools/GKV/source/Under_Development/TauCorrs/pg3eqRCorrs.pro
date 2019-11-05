FUNCTION GKVs3D::pg3eqRCorrs, skip = iskip, debug=debug
;
; Takes a GKVs3D object as an argument.
; Assumes that axises 1 and 3 correspond to
; ignorable coordinates.
;
; For each values of axis2, computes full width at
; half maximum of cross correlation function function.
;
; Returns GKVs1D object containing Rcorr (that is,
; the full width at half maximum of correlation function
; vs. axis2).
;
; Written by W.M. Nevins
;	9/24/00
;
skip = 1
IF(N_ELEMENTS(iskip) EQ 1) THEN skip = iskip > 1
result3D = {GKVs3D}						; Create GKV structure to hold result
FOR i=0, N_TAGS(result3D) - 1 DO 	$
	result3D.(i) = self.(i)				; Copy fields from 'self' to 'result2D'
axisMnemonic = self.Grid2.mnemonic
axisTitle = self.Grid2.title
result3D.mnemonic = axisMnemonic + '_corr'
result3D.title = '!12l!X' + Subscript(axisTitle) + '!X!N{' + self.title + '}'
indices = Self -> IndexRemove(3)			; Return 'indices' of Self with time-field removed.
result3D.indices = PTR_NEW(indices)
; **** go over use of result.grid1, ...
irange = self.Grid2.irange
imax = irange[1]-irange[0] 
result3D.Grid1.values = PTR_NEW((*self.Grid2.values)[irange[0]:irange[1]])
result3D.Grid1.irange = [0,imax]
vmin = (*self.Grid2.values)[irange[0]]
vmax = (*self.Grid2.values)[irange[1]]
result3D.Grid1.range = [vmin, vmax]
values = FLTARR(irange[1]-irange[0]+1)
FOR i=0, imax, skip DO BEGIN
	index = irange[0]+i
	tempObj = self -> Slice(axis=2, index=index)
	tempObj -> Norm
	tempcorrs = self -> xcorr(ref=tempObj)
	tempcorrsmax = tempcorrs -> Slice(axis=3, /max)
	tcmaxmax  = tempcorrsmax -> Slice(axis=1, /max)
	values[i] = tcmaxmax -> FullWidth(debug=debug)
	IF(KEYWORD_SET(debug)) THEN print, index, (*self.Grid2.values)[index], values[i]
	tempObj -> Trash
	tempcorrs -> Trash
	tempcorrsmax -> Trash
	tcmaxmax -> Trash
	HEAP_GC, Verbose=debug
ENDFOR
IF(skip GT 1) THEN BEGIN					; Interpolate to get skipped values
	FOR i=0, imax, skip DO BEGIN
		index = irange[0]+i
		nextIndex = index+skip
		dValues=values[nextIndex]-values[index]
		dx = (*self.Grid2.values)[nextIndex] - (*self.Grid2.values)[Index]
		FOR j=1, iskip-1 DO BEGIN
			delx = (*self.Grid2.values)[Index+j] - (*self.Grid2.values)[Index]
			values[index+j] = dValues/dx*delx
		ENDFOR
	ENDFOR 
ENDIF
vmin = GKVsd_MIN(values, Max=vmax)
result3D.vrange = [vmin, vmax]
result3D.values = PTR_NEW(values)
threeDObj = OBJ_NEW('GKVs3D', result3D)
indices = threeDObj -> IndexMax(3)
PTR_FREE, result3D.indices
result3D.indices = PTR_NEW(indices)
result2D = {GKVs2D}
FOR i=0,N_TAGS(result2D)-1 DO		$
	result2D.(i) = result3D.(i)
twoDObj = OBJ_NEW('GKVs2D', result2D)
indices = twoDObj -> IndexMax(1)
PTR_FREE, result2D.indices
result2D.indices = PTR_NEW(indices)
result1D = {GKVs1D}
FOR i=0, N_TAGS(result1D) - 1 DO 	$
	result1D.(i) = result2D.(i)				; Copy fields from 'result2D' to 'result1D'
	result1D.Grid1 = result2D.Grid2
RCorrObj = OBJ_NEW('GKVs1D', result1D)
;GKVdelete, result2D							; Get rid of result2D
return, RcorrObj
END ; ****** GKVs3D::pg3eqTauCorrs ****** ;
