FUNCTION GKVsd::IndexString, plotAxis
;
; Returns a string array containing index elements for a GKVsd object.
;
; The argument 'plotAxis' is a 1-D array whose length as equal to the 
; dimensionality of the GKVsd object (GKVs1D -> 1, GKVs2D -> 2, etc).
; If an element of 'plotAxis' is positive, the value corresponds to the
; number of one of the axis being plotted:
;
;	1 for x (in the sense of IDL... that is, horizontal), 
;	2 for y (in the sense of IDL... that i , vertical), 
;	etc.
;
; Negative elements of 'plotAxis' correspond to axis that are being 'sliced' in
; this plot, and the absolute value of this elements is the index at which the 
; data is being sliced.
;
; Zero values of 'plotAxis' signify that the corresponding element of 
; 'Indices' should be left alone.
;
indices = *self.indices				; Get GKVsd object's Indices
info = SIZE(indices)
Nindices = info[1]					; Find length of 'indices' array
;
; Loop to replace *'s in 'indices' array with corresponding axis title
; and (for negative values of 'plotAxis') grid values.
;
iaxis = 0							; index into 'plotAxis' array (elements 1 thru ndims)
FOR i=0, nIndices-1 DO BEGIN			; i indexes 'Indices' array, 
	IF(indices[i] EQ '*') THEN BEGIN	; which can have more than ndims elements
		iaxis=iaxis+1
		IF(plotAxis[iaxis] EQ 0) THEN GOTO, DONE
		axis_str = STRING(iaxis, FORMAT='(i1)')
		command_str = 'grid = self.grid' + axis_str
		ok = EXECUTE(command_str)
		IF(plotAxis[iaxis] GT 0) THEN BEGIN
			indices[i] = grid.title
		ENDIF ELSE BEGIN
			Value = (*grid.values)[-plotAxis[iaxis]]
			valueStr = STRING(value, FORMAT='(G10.3)')
			indexStr = grid.title + '=' + valueStr
			indexStr = STRCOMPRESS(indexStr, /REMOVE_ALL)
			indices[i] = indexStr
		ENDELSE
	ENDIF
	DONE	:	; Don't do anything
ENDFOR
RETURN, indices
END ; ****** GKVsd::IndexString ****** ;