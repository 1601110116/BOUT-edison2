FUNCTION GKVsd::IndexRemove, axisNum
;
; returns string array of 'indices' with element 
; corresponding to 'axisNum' removed.
;
; Written by W.M. Nevins
;	4/6/00
;
indices = *self.indices
info=SIZE(indices)
Nindices=info[1]
result = STRARR(Nindices-1)
;
iaxis=0
j=0
FOR i=0, Nindices-1 DO BEGIN
	IF(indices[i] EQ "*") THEN BEGIN
		iaxis=iaxis+1
		IF(iaxis NE axisNum) THEN BEGIN
		result[j] = indices[i]
		j=j+1
		ENDIF
	ENDIF ELSE BEGIN
		result[j] = indices[i]
		j=j+1
	ENDELSE
ENDFOR
RETURN, result
END ; ****** GKVsd::IndexRemove ****** ;