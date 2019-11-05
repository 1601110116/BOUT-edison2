FUNCTION GKVs3D::DeltaSq, arg, _Extra=Extra
;
; Purpose:
;
;	Decomposes 'self' into average over selected axis
;	and deviations from this average using GKV Delta 
;	routine. In addition, DeltaSq also computes
;	(and returns) the intensity of the deviations
;	from the average -- both as a function of the 
;	remaining indepedent variables and as a function
;	of the third dependent variable (usually time).
;
; Argument:
;
;	Any valid axis identifier. Used to identify the axis to
;	be averaged over.
;
; Keywords:
;
;	The axis to be averaged over (and the range over 
;	which this average is performed) can also be identified 
;	by using keywords (see routine AxisIrange).
;
;
; Written by W.M. Nevins
;	7/25/08
;
; Use AxisIrange to parse input line	
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DeltaSq called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Compute average over selected axis and deviations from this average
;
; Check boundary contidions, don't include last point in average
; if "Periodic (closed)"
;
self -> get, axis=axis, boundary=boundary, irange=iRange, gridValues=yvaluePtr
yValues = *yValuePtr
IF( STRCMP(boundary, "Periodic (closed)", /FOLD_CASE) ) THEN iRange[1] = N_ELEMENTS(yValues) - 2
output = self -> Delta(Axis=axis, irange=iRange)
IF( STRCMP(boundary, "Periodic (closed)", /FOLD_CASE) ) THEN   $
  output.delta -> Set, Axis=axis, boundary="Periodic (open)"
output = CREATE_STRUCT("Name", "DeltaSq", output)
;
; Compute intensities
;
Title="!12I!X{" + self.title + "}"
commandString = "Mnemonic = 'I_' + self.Grid" + STRTRIM(STRING(axis),2) + ".Mnemonic"
ok = EXECUTE(commandString)
temp = output.delta -> AbsSq()
I_x = temp -> AVG(axis, irange=Irange)
temp -> Trash
I_x -> Set, Title=title, Mnemonic=mnemonic
;
; Check boundary contidions, don't include last point in average
; if "Periodic (closed)"
;
I_x -> get, axis=1, boundary=boundary, irange=iRange, gridValues=xvaluePtr
xValues = *xValuePtr
IF( STRCMP(boundary, "Periodic (closed)", /FOLD_CASE) )  THEN BEGIN
	IF(iRange[1] EQ ( N_Elements(xValues) - 1 ) ) THEN iRange[1] = iRange[1] - 1
ENDIF
I   = I_x  -> AVG(axis=1, irange=iRange)
I   -> Set, Title=title, Mnemonic="I"
;
; Compute intensity of axis-averaged quantity
;
temp=output.avg -> AbsSq()
temp -> get, errorbars=errors
ptr_free, errors
Iavg= temp -> AVG(axis=1, irange=iRange)
temp -> Trash
;
; add intensities to output structure and return 
;
output = CREATE_STRUCT(output, "I_x", I_x, "I", I, "Iavg", Iavg)
RETURN, output
END ; ****** GKVs3D::DeltaSq ****** ;
