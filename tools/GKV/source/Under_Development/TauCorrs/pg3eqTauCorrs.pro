FUNCTION GKVs3D::pg3eqTauCorrs, debug=debug
;
; Takes a GKVs3D object as an argument.
; Assumes that axises 1 and 3 correspond to
; ignorable coordinates.
;
; For each values of axis2, computes full width at
; half maximum of auto correlation function function.
;
; Returns GKVs1D object containing Taucorr (that is,
; the full width at half maximum of correlation function
; vs. axis2).
;
; Written by W.M. Nevins
;	4/18/00
;
result2D = {GKVs2D}					; Create GKVs1D structure to hold result
FOR i=0, N_TAGS(result2D) - 1 DO 	$
	result2D.(i) = self.(i)				; Copy fields from 'self' to 'result2D'
result2D.mnemonic = 'tau_corr'
result2D.title = '!4t!N!Ic!N{' + self.title + '}'
indices = Self -> IndexRemove(3)		; Return 'indices' of Self with time-field removed.
result2D.indices = PTR_NEW(indices)
irange = self.Grid2.irange
IF(irange[0] EQ 0) THEN irange[0]=1
irange[1] = irange[1]-1
imax = irange[1]-irange[0] 
result2D.Grid1.values = PTR_NEW((*self.Grid2.values)[irange[0]:irange[1]])
result2D.Grid1.irange = [0,imax]
vmin = (*self.Grid2.values)[irange[0]]
vmax = (*self.Grid2.values)[irange[1]]
result2D.Grid1.range = [vmin, vmax]
values = FLTARR(irange[1]-irange[0]+1)
vPhase = values
vPhaseErrors = values
FOR i=0, imax DO BEGIN
	index = irange[0]+i
	tempObj = self -> Slice(axis=2, index=index)
	tempObj -> Norm
	tempcorrs = tempObj -> xcorr()
	tempcorrsmax = tempcorrs -> Slice(axis=1, /max, /maxlocation)
	values[i] = tempcorrsmax.slice -> FullWidth(debug=debug)
	tempcorrsmax.maxlocation -> SignalWindow, tau=[-values[i]/2., values[i]/2.]
	vPhase[i] = tempcorrsmax.maxlocation -> slope(error=error)
	vPhaseErrors[i] = error
	IF(KEYWORD_SET(debug)) THEN print, index, (*self.Grid2.values)[index], values[i]
	tempObj -> Trash
	tempcorrs -> Trash
	tempcorrsmax.slice -> Trash
	tempcorrsmax.maxlocation -> Trash
	HEAP_GC, Verbose=debug
ENDFOR
vmin = GKVsd_MIN(values, Max=vmax)
result2D.vrange = [vmin, vmax]
result2D.values = PTR_NEW(values)
twoDObj = OBJ_NEW('GKVs2D', result2D)
indices = twoDObj -> IndexMax(1)
PTR_FREE, result2D.indices
result2D.indices = PTR_NEW(indices)
result1D = {GKVs1D}
FOR i=0, N_TAGS(result1D) - 1 DO 	$
	result1D.(i) = result2D.(i)				; Copy fields from 'result2D' to 'result1D'
tauCorrObj = OBJ_NEW('GKVs1D', result1D)
vPhaseObj = taucorrObj -> MakeCopy(/noValues)
vPhaseObj.values = PTR_NEW(vPhase)
vmin = GKVsd_MIN(vPhase, Max=vmax)
vPhaseObj.vrange=[vmin, vmax]
vPhaseObj.title = 'V!I!4u!N!X'
vPhaseObj.mnemonic = 'v_phase'
vPhaseObj.ErrorBars = PTR_NEW(vPhaseErrors)
;GKVdelete, result2D							; Get rid of result2D
return, {TauCorr:taucorrObj, vPhase:vPhaseObj}
END ; ****** GKVs3D::pg3eqTauCorrs ****** ;
