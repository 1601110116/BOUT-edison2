PRO GKVsd::repairValues
;
;
values = *self.values
badIndices = WHERE(FINITE(values) EQ 0)
IF(badIndices[0] NE -1) THEN BEGIN
	values(badIndices) = 0.
	ptr_free, self.values
	self.values = PTR_NEW(values)
	print, N_ELEMENTS(badIndices), " bad values out of ", N_ELEMENTS(values)
	print, "****** REPAIRED ******"
ENDIF ELSE BEGIN
	print, "No bad values found"
ENDELSE
RETURN
END
