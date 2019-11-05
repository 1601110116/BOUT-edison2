
FUNCTION GKVs1D::AxisIrange, arg, _Extra=extra
;
; Purpose:
;
;	This function parses input line, and returns the number of the selected axis. 
;	Information on 'irange' for this axis can be supplied on the command line.
;	If it is, then the 'SignalWIndow' will be invoked to modify 'irange' in self.
;
;	Output:	Returns the (integer) value of the selected axis.  If
;			no axis can be identified, -1 is returned.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Input Keywords:
;
;	     Axis		If no argument is provided, then this keyword may be 
;			used to identify Grid to be 'Restricted'. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the selected independent 
;			variable, and reset the signal window on this axis.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable (i.e., it is interpreted in the units
;			of the corresponding independent variable), NOT the integer 'irange'
;			(that is, NOT as an integer index into the grid.values array).
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable.
;
;	Output Keywords:
;
;	AxisValue	The range of values for the selected axis.  If only one axis value
;			is input, then returns this value.
;
;      AxisIrange	The Irange for the selected axis.  If only one axis value is input,
;			then returns the index closest to the selected value.
;
; Written by W.M. Nevins
;	8/17/00
; Modified by W.M. Nevins
;	11/6/2008
; for consistent use of the AxisValue keyword, and to add the AxisIrange keyword.
;
; First find axis identifier
;
nDims = self -> NumDims()
iaxis = -1
IF(N_PARAMS() EQ 1) THEN BEGIN			; Take argument as axis identifier if present
	iaxis = arg
ENDIF ELSE BEGIN				; Else, look for the keyword 'Axis'
	result = GetKeyWord('Axis', extra)
	IF(Query_Integer(result)) THEN iaxis = result
	IF(TypeOf(result) EQ 7) THEN BEGIN
		IF(result NE 'undefined') THEN iaxis = result
	ENDIF
ENDELSE
;
; If 'iaxis' is a mnemonic, then convert it to an integer axis identifier
;
IF(TypeOF(iaxis) EQ 7) THEN iaxis = self -> AxisNumber(iaxis)
IF(iaxis EQ 0) THEN BEGIN 
	IF(NOT Query_Integer(iaxis)) THEN BEGIN
		MESSAGE, "Illegal axis identifier", /INFORMATIONAL
		RETURN, -1
	ENDIF
ENDIF
;
; Check if we have a valid (integer) axis identifier
;
IF(iaxis NE -1) THEN BEGIN
	IF((iaxis LT 1) OR (iaxis GT nDims)) THEN BEGIN
		MESSAGE, "Illegal axis number", /INFORMATIONAL
		RETURN, -1	
	ENDIF
;
; or (if we haven't found an axis identifer yet)
; look for 'mnemonic = ...'
;
ENDIF ELSE BEGIN
	axisInfo = self -> GetAxis(extra)
	FOR iaxis = 1, nDIms DO BEGIN
		axisValue = axisInfo.(iaxis-1)
		IF(TypeOf(axisValue) NE 7) THEN BEGIN
			thisAxisValue = axisValue
			thisRange = axisValue
			IF(N_ELEMENTS(axisValue) EQ 1) THEN thisRange=[axisvalue, axisValue]
			self -> SignalWindow, axis=iaxis, range = thisRange
			CASE N_ELEMENTS(axisValue) OF
				1:	BEGIN
					self -> Get, axis=iaxis, irange=irange
					thisAxisIrange = irange[0]
					END
				2:	self -> Get, axis=iaxis, irange=thisAxisIrange
			     ELSE:
			ENDCASE
			RETURN, iaxis
		ENDIF
	ENDFOR
	;
	; Failed to find a valid axis mnemonic 
	;
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, -1
ENDELSE
;
; You should only get here if you already have an axis ID
;
; Check for 'range' on command line 
; if it is not already set (by 'mnemonic = ...')
;
result = GetKeyWord('range', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	self -> SignalWindow, axis=iaxis, range = result
	CASE N_ELEMENTS(result) OF
		1:	BEGIN
			self -> Get, axis=iaxis, irange=irange
			thisAxisIrange = irange[0]
			END
		2:	self -> Get, axis=iaxis, irange=thisAxisIrange
		ELSE:
	ENDCASE
ENDIF
;
; Check for 'irange' on the command line
;
result = GetKeyWord('irange', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	thisAxisIrange=result	
	CASE N_ELEMENTS(result) OF
		1:	self -> SignalWindow, axis=iaxis, irange = [result, result]
		2:	self -> SignalWindow, axis=iaxis, irange = result
		ELSE:
	ENDCASE
ENDIF
;
; and we're done ...
;
RETURN, iaxis

END ; ****** GKVs1D::AxisIrange ****** ;
