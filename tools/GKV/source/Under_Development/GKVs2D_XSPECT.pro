FUNCTION GKVs2D::XSPECT, Ref = RefObj, All=all, No_Avg=noAverage, Dw=dataWindow, _Extra=extra
;
; Returns a GKVsd object containing the cross spectrum between the
; data in 'self' and the data in RefOjb (which BOTH must be GKVsd objects)
;
; If no RefObj is specified, then returns auto-spectrum of data in 'self'.
; 
; Remaining keywords are passed to XSPECT (see definition of XSPECT to understand their use)
;
;	Written 2/6/00 by W.M. Nevins
;
;	Revised 4/18/01 by W.M. Nevins
;	to fix bug relating to scaling of independent varialbe grids
;
rflag = 0
No_Avg = 1
IF(N_ELEMENTS(noAverage)) THEN No_Avg = noAverage
ndims = self -> NumDims()						; Get dimensionaliity of data in 'self'
IF((ndims LT 1) OR (ndims GT 4)) THEN BEGIN
	MESSAGE, 'Bad dimensionality, cannot form cross correlations function', /Informational
	RETURN, 0
ENDIF
IF(N_ELEMENTS(RefObj) EQ 0) THEN BEGIN			; No reference object... so pass on to GKVs1D
	result = self -> GKVs1D::XSPECT(No_Avg=noAverage, Dw=dataWindow, _Extra=extra)
	RETURN, result
ENDIF
refType = TypeOf(RefObj)
IF (refType NE 11) THEN BEGIN				; RefObj is NOT an object
	MESSAGE, 'Ref is not an Object', /Informational
	RETURN, 0
ENDIF
refdims = RefObj -> NumDims()					; Get dimensionality of data in 'RefObj'
IF(refdims EQ ndims) THEN BEGIN					; Dimensionality is same... so pass on to GKVs1D
	result = self -> GKVs1D::XSPECT(Ref = RefObj, No_Avg=noAverage, Dw=dataWindow, _Extra=extra)
	RETURN, result
ENDIF
IF(refdims NE (ndims-1) ) THEN BEGIN
	MESSAGE, "Expect dimensionality of Ref = dimensionality 'self', or self - 1" , /informational
	RETURN, 0
ENDIF
spacing = FLTARR(ndims+1)
ngrids  = LONARR(ndims+1)
iref = 1
iaxis = 0
FOR iself=1, ndims DO BEGIN						; 
	iselfStr = STRING(iself, FORMAT='(I1)')
	command_str = 'selfGrid = self.Grid' + iselfStr		; Get iself^th Grid structure of 'self'
	ok = EXECUTE(command_str)
	irefstr = STRING(iref, FORMAT='(I1)')
	command_str = 'refGrid = RefObj.Grid' + irefstr
	ok = EXECUTE(command_str)						; Get iref^th Grid structure of 'RefObj'
	ismin = selfGrid.irange[0]
	ismax = selfGrid.irange[1]
	IF(GKV_GridSame(selfGrid,RefGrid,All=all, /force)) THEN BEGIN	; Check to see if the Grid structures match (force match if possible)
		IF(KEYWORD_SET(All)) THEN BEGIN
			ismin = 0
			ismax = N_ELEMENTS(*selfGrid.values) - 1
		ENDIF
		uniform = GKVsd_UniformGrid((*selfGrid.values)[ismin:ismax])
		IF(NOT uniform) THEN BEGIN					; Check for uniform grids
			MESSAGE, 'XCORR only defined for data on uniform grid(s)', /Informational
			RETURN, 0
		ENDIF
		spacing(iself)	= 2.0*!PI/( (*selfGrid.values)[ismax] - (*selfGrid.values)[ismin] )
		ngrids(iself)	= ismax-ismin+1
		iref = iref + 1
	ENDIF ELSE BEGIN							; This selfGrid does not match corresponding refGrid
		IF(iaxis) THEN BEGIN						;	this should only happen once!
			MESSAGE, 'Bad match between reference grids and self grids', /Informational
			RETURN, 0
		ENDIF
		iaxis= iself							; Prepare for loop which calls XCORR
		imin = ismin
		imax = ismax
	ENDELSE
	IF(iself EQ nDims) THEN dt = (*selfGrid.values)[ismin+1] - (*selfGrid.values)[ismin]
ENDFOR
IF(iaxis EQ 0) THEN BEGIN
	MESSAGE, "Couldn't find cross correlation axis", /Informational
	RETURN, 0
ENDIF

;dt = spacing(ndims)
selfValuePtr =  self -> GetValues(All=all)				; Returns pointer to values within signal window
refValuePtr = RefObj -> GetValues(All=all)
refValues = *refValuePtr
refInfo = SIZE(refValues)
FOR i = imin, imax DO BEGIN
CASE ndims OF
	2:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *], refInfo[1])
			2:	selfValues = REFORM((*selfValuePtr)[*, i], refInfo[1])
		ENDCASE
	3:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *], refInfo[1], refInfo[2])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *], refInfo[1], refinfo[2])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i], refInfo[1], refInfo[2])
		ENDCASE
	4:	CASE iaxis OF
			1:	selfValues = REFORM((*selfValuePtr)[i, *, *, *], refInfo[1], refInfo[2], refInfo[3])
			2:	selfValues = REFORM((*selfValuePtr)[*, i, *, *], refInfo[1], refInfo[2], refInfo[3])
			3:	selfValues = REFORM((*selfValuePtr)[*, *, i, *], refInfo[1], refInfo[2], refInfo[3])
			4:	selfValues = REFORM((*selfValuePtr)[*, *, *, i], refInfo[1], refInfo[2], refInfo[3])
		ENDCASE
ENDCASE
SpectValues_i = XSpect(selfValues, Reference = RefValues, No_Avg = No_Avg, Dt=dt, _Extra=extra)
IF(i EQ imin) THEN BEGIN							; Initialization duties
	selfInfo = SIZE(*selfValuePtr)
	spectInfo = SIZE(SpectValues_i)
	spectDims = spectInfo[0]
	nOmegas = spectinfo[spectDims]
	CASE ndims OF
		2:	spectValues = FLTARR(selfInfo[1], nOmegas)
		3:	spectValues = FLTARR(selfinfo[1], selfInfo[2], nOmegas)
		4:	spectValues = FLTARR(selfinfo[1], selfInfo[2], selfInfo[3], nOmegas)
	ENDCASE
ENDIF
CASE ndims OF
	2:	spectValues[i, *] = spectValues_i
	3:	CASE iaxis OF
		1:	spectValues[i,*,*] = spectValues_i
		2:	spectValues[*,i,*] = spectValues_i
		ENDCASE
	4:	CASE iaxis OF
		1:	spectValues[i,*,*,*] = spectValues_i
		2:	spectValues[*,i,*,*] = spectValues_i
		3:	spectValues[*,*,i,*] = spectValues_i
		ENDCASE
ENDCASE
ENDFOR
PTR_FREE, selfValuePtr							; Free up pointers
PTR_FREE, refValuePtr
spectValuesPtr = PTR_NEW(spectValues)
;spectInfo = SIZE(spectValues)
;ngrids[ndims] = spectInfo[ndims]
;
; Make GKVsd object to contain cross spectrum
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
;
; Start filling tags of result
;
result.mnemonic = 'XSpect_' + self.mnemonic	 + '_' + RefObj.mnemonic	; Set Mnemonic
result.title = 'S{' + self.title + ', ' + RefObj.title + '}'			; Set Title
PTR_FREE, result.values								; Set Values pointer
result.values = spectValuesPtr
vmin = GKVsd_MIN(spectValues, MAX=vmax)						; Set Vrange
vrange = [vmin, vmax]
result.vrange = vrange
PTR_FREE, result.ErrorBars						; Set ErrorBars pointer ...
result.ErrorBars = PTR_NEW()					; No error bars on spectral density?
FOR i=1, ndims DO BEGIN						; Loop over ndims instances of the Grid Class
	dimStr = STRING(i, FORMAT='(I1)')			
	GridStr = 'result.Grid' + dimStr
	command_str = 'Grid = ' + gridStr			; Copy ith instance of the Grid Class
	ok = EXECUTE(command_str)					;	 into Grid
	Grid.Mnemonic = 'k_' + Grid.mnemonic 			; Set Grid Mnemonic
	norm = 1.
	IF(Grid.Mnemonic EQ 'k_t') THEN BEGIN
		Grid.mnemonic = 'omega'
		norm = FLOAT(ngrids[i]-1)/nOmegas			; Constant to correct dOmega
	ENDIF
	title = Grid.title						; Set Grid Title					
	IF(title EQ 't') THEN BEGIN
		Grid.title = '!4x!X'
	ENDIF ELSE BEGIN
		Grid.title = '!8k!X!I' + title + '!N!X'
	ENDELSE
	IF(i NE iaxis) THEN BEGIN					; Change grid array for all axis EXCEPT 'iaxis'
		Grid.units = '1/' + Grid.units			; Update grid units
		PTR_FREE, Grid.values						
		n = ngrids[i]
		IF(i EQ ndims) THEN n=nOmegas
		gridvals = norm*spacing(i)*(FINDGEN(n) -n/2)	; Generate array of grid values
		GridPtr = PTR_NEW(gridvals)	
		Grid.values = GridPtr					; Set Grid Values pointer
		Grid.boundary = 'Periodic'				; Spectral density BC's are always Periodic?		
		Grid.uniform = 1					; Spectral density computed on uniform grid
		min = MIN(gridvals, Max=max)			; Set Grid plot Range
		range = [min, max]
		Grid.range = range
		irange = [0, n-1]					; Set Grid signal window range, irange
		Grid.irange = irange
		command_str = GridStr + ' = Grid'		; Copy ith instance of Grid Class
		ok = EXECUTE(command_str)				;	back into result
	ENDIF
ENDFOR
RETURN, result	
END ; ****** GKVs2D::XSPECT ****** ;

