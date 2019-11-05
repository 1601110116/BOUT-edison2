PRO CompareStructures, arg1, arg2
;
; Compares two structures
;
; Written by W.M. Nevins
;	3/12/2009
;
type1 = TypeOf(arg1)
IF(type1 NE 8) THEN BEGIN
	MESSAGE, "Arg1 is not a structure", /INFORMATIONAL
	RETURN
ENDIF
type2 = TypeOf(arg2)
IF(type2 NE 8) THEN BEGIN
	MESSAGE, "Arg2 is not a structure", /INFORMATIONAL
	RETURN
ENDIF

nTags1 = N_TAGS(arg1)
nTags2 = N_TAGS(arg2)
tags1 = TAG_NAMES(arg1)
tags2 = TAG_NAMES(arg2)

FOR i=0, nTags1-1 DO BEGIN
	type1 = TypeOf(arg1.(i))
	type2 = TypeOf(arg2.(i))
	IF(type1 NE type2) THEN $
	PRINT, tags1[i], TypeOf(arg1.(i)), tags2[i], TypeOf(arg2.(i)), FORMAT = "(5X, A20, I3, 5X, A20, I3)"
	IF(TypeOf(arg1.(i)) EQ 8) THEN BEGIN
		PRINT, tags1[i], TypeOf(arg1.(i)), tags2[i], TypeOf(arg2.(i)), FORMAT = "(A20, I3, 5X, A20, I3)"		
		CompareStructures, arg1.(i), arg2.(i)
	ENDIF
ENDFOR
RETURN
END  ; ****** CompareStructures ****** ;

