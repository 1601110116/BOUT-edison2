PRO GKVs1D::GridScale, arg, _extra=extra
;
; Purpose:
;
;	This proceedure rescales selected grid, allowing user
;	to change units.
;
;
; Arguments:
;
;	arg (optional). If present, this argument is a legal axis identifier, either:
;		an integer between 1 and the dimensionality of this object, or
;		a string containing an axis mnemonic
;
; Keywords:
;
;	axis	Set equal to a legal axis identifier, either:
;		an integer between 1 and the dimensionality of this object, or
;		a string containing an axis mnemonic
;
; Written by W.M. Nevins
;	7/13/00
;
IF(N_ELEMENTS(arg) NE 0) THEN BEGIN		; Get axis ID from 'arg' if one is present
	iaxis = arg					;	if 'arg' is a string, try to match with a mnemonic
	IF(TypeOF(arg) EQ 7) THEN iaxis = self -> AxisNumber(arg)
								; Otherwise, get axis ID from keyword "Axis"
	IF(NOT Query_Integer(iaxis)) THEN GOTO, tryagain	
ENDIF ELSE BEGIN
tryagain:
	iaxis=0
	result = GetKeyWord('axis', extra)
	IF(Query_Integer(result)) THEN iaxis = result					; command line of form axis = axisnumber
	IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result)	; command line of the form axis = 'mnemnic'
ENDELSE
IF(NOT Query_Integer(iaxis)) THEN BEGIN
	MESSAGE, 'No legal axis identifier was supplied', /INFORMATIONAL
	RETURN
ENDIF
;
; Get grid structure corresponding to selected axis
;
axisStr = STRING(iaxis, FORMAT='(i1)')
commandStr = 'grid = self.Grid' + axisStr
ok = EXECUTE(commandStr)
;
;