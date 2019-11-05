PRO GKVs1D::Repair, arg, _Extra=extra
;
; Purpose:
;
;	This proceedure 'repairs' objects by interpolating 
;	from neighboring grid points to replace missing data
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.  The repair
;			will be performed by interpolating over this axis.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable which is to be
;			interpolate over in making the repair, together with 
;			a range of values which enclose the bad data (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for  which is to be interpolate over in making the repair, together with 
;			a range of values which enclose the bad data (see above).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array  
;			corresponding to indices into the independent variable's 'value'
;			array which enclose the bad data (see above).
;			This two-element array is interpreted as the desired (integer) 
;			IRANGE in the independent variable, NOT the (real) 'RANGE'
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable which enclose the bad data (see above).
;			This two-element array is interpreted as the desired (real) 
;			RANGE in the independent variable, NOT the (integer) 'IRANGE'
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
; Written by W.M. Nevins
;	1/29/02
;
; Revised by W.M. Nevins
;	3/20/01
;

CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Avg called with too many arguments', /INFORMATIONAL
				RETURN, self
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, self
ENDIF
resultStr = {GKVsd}
FOR i=0, N_TAGS(resultStr)-1 DO resultStr.(i) = self.(i)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = grid.irange
axisTitle = grid.title

resultStr.Title = '!12<!X' + self.title + '!12>!X!D' + axisTitle  + '!N'
newIndices = self -> IndexRemove(iaxis)
resultStr.indices = PTR_NEW(newIndices)
;
; Compute Average
;
dx = (*grid.values)[(irange[0]+1):irange[1]] - (*grid.values)[irange[0]:(irange[1]-1)]
values = 0.5*((*self.values)[(irange[0]+1):irange[1]] + (*self.values)[irange[0]:(irange[1]-1)])*dx
values = values/( (*grid.values)[irange[1]] - (*grid.values)[irange[0]] )
values = TOTAL(values)
resultStr.values = PTR_NEW(values)
result = OBJ_NEW("GKVsd", resultStr)

RETURN, result
