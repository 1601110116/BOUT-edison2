
PRO GKVs1D::Restrict
;
; Purpose:
;
;	This proceedure restricts 'self', removing
;	values outside the current signalwindow.
;
;	Arguments:
;
;			None
;
;	KeyWords:
;
;			None
;
; Written by W.M. Nevins
;	8/21/00
;
; Get current 'values' array
;
oldValues = *self.values
; 
; get current signal window
;
irange = self.grid1.irange
imin = irange[0]
imax = irange[1]
;
; form new values array
;
newValues = oldValues[imin:imax]
;
; Load new values array into 'self'
;
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
;
; Check for errorbars
;
IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
	oldErrors = self.ErrorBars
	newErrors = (*oldErrors)[imin:imax]
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(newErrors)
ENDIF
;
; Restrict Grid
;
self -> GridRestrict, 1
;
; ... and we're done
;
RETURN
END ; ****** GKVs1D::Restrict ****** ;


PRO GKVs2D::Restrict
;
; Purpose:
;
;	This proceedure restricts 'self', removing
;	values outside the current signalwindow.
;
;	Arguments:
;
;			None
;
;	KeyWords:
;
;			None
;
; Written by W.M. Nevins
;	8/21/00
;
; Get current 'values' array
;
oldValues = *self.values
; 
; get current signal window
;
irange = self.grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.grid2.irange
jmin = jrange[0]
jmax = jrange[1]
;
; form new values array
;
newValues = oldValues[imin:imax, jmin:jmax]
;
; Load new values array into 'self'
;
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
;
; Check for errorbars
;
IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
	oldErrors = self.ErrorBars
	newErrors = (*oldErrors)[imin:imax, jmin:jmax]
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(newErrors)
ENDIF
;
; Restrict Grids
;
FOR i=1,2 DO self -> GridRestrict, i
;
; ... and we're done
;
RETURN
END ; ****** GKVs2D::Restrict ****** ;

PRO GKVs3D::Restrict
;
; Purpose:
;
;	This proceedure restricts 'self', removing
;	values outside the current signalwindow.
;
;	Arguments:
;
;			None
;
;	KeyWords:
;
;			None
;
; Written by W.M. Nevins
;	8/21/00
;
; Get current 'values' array
;
oldValues = *self.values
; 
; get current signal window
;
irange = self.grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.grid2.irange
jmin = jrange[0]
jmax = jrange[1]
krange = self.grid3.irange
kmin = krange[0]
kmax = krange[1]
;
; form new values array
;
newValues = oldValues[imin:imax, jmin:jmax, kmin:kmax]
;
; Load new values array into 'self'
;
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
;
; Check for errorbars
;
IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
	oldErrors = self.ErrorBars
	newErrors = (*oldErrors)[imin:imax, jmin:jmax, kmin:kmax]
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(newErrors)
ENDIF
;
; Restrict Grids
;
FOR i=1,3 DO self -> GridRestrict, i
;
; ... and we're done
;
RETURN
END ; ****** GKVs3D::Restrict ****** ;

PRO GKVs4D::Restrict
;
; Purpose:
;
;	This proceedure restricts 'self', removing
;	values outside the current signalwindow.
;
;	Arguments:
;
;			None
;
;	KeyWords:
;
;			None
;
; Written by W.M. Nevins
;	8/21/00
;
; Get current 'values' array
;
oldValues = *self.values
; 
; get current signal window
;
irange = self.grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.grid2.irange
jmin = jrange[0]
jmax = jrange[1]
krange = self.grid3.irange
kmin = krange[0]
kmax = krange[1]
lrange = self.grid4.irange
lmin = lrange[0]
lmax = lrange[1]
;
; form new values array
;
newValues = oldValues[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
;
; Load new values array into 'self'
;
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
;
; Check for errorbars
;
IF( PTR_VALID(self.ErrorBars) ) THEN BEGIN
	oldErrors = self.ErrorBars
	newErrors = (*oldErrors)[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(newErrors)
ENDIF
;
; Restrict Grids
;
FOR i=1,4 DO self -> GridRestrict, i
;
; ... and we're done
;
RETURN
END ; ****** GKVs4D::Restrict ****** ;

PRO GKVs1D::GridRestrict, arg, _EXTRA=extra
;
; Purpose:
;
;	This proceedure reforms the Grid identified by either
;	the argument ('arg'), the (runtime) keyword 'axis', or any 
;	valid grid mnemonic to restrict the grid to only those  
;	values within its signal window. 
;	
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	  Axis	If no argument is provided, then this keyword may be 
;			used to identify Grid to be 'Restricted'. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the Grid to be restricted 
;			and reset the signal window on this axis (before restriction).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before restricting the selected Grid.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable.
;
;	Side Effects:
;
;			Generally, restricting the grid can leave 'self' in an invalid
;			state (the grid dimension will not match the corresponding dimension
;			of 'values').  Hence, this method should not be called directly by 
;			users, but is rather part of the 'private' interface for use by
;			developers of GKV.  The 'values' array must be reformed by the 
;			developer in conjunction with the use of GridRestrict.
;
; Written by W.M. Nevins
;	8/14/00

;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	iaxis = self -> AxisIrange(     _Extra=extra)
	1	:	iaxis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'GridRestrict called with too many arguments', /INFORMATIONAL
				RETURN
			END
ENDCASE
IF(iaxis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN
ENDIF
axisString = STRING(iaxis, FORMAT='(i1)')
;
; Get grid structure
;
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *Grid.values
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
npoints = imax - imin + 1
IF (N_ELEMENTS(gridValues) EQ npoints) THEN RETURN
;
; Form new grid values
;
newValues = gridValues[imin:imax]
PTR_FREE, Grid.values
Grid.values = PTR_NEW(newValues)
Grid.irange = [0, npoints-1]
grid.uniform = GKVsd_UniformGrid(newValues)
;
; Restrict 'range' so that it falls within
; the restricted grid if necessary
;
range = Grid.range
xmin = range[0]
xmax = range[1]
xmin = xmin > newValues[0]
xmax = xmax < newValues[npoints-1]
Grid.range = [xmin, xmax]
;
; put result back into 'self'
;
commandString = 'self.Grid' + axisString + ' = Grid'
ok = EXECUTE(commandString)
; 
; and we're done ...
;
RETURN
END ; ****** GKVs1D::GridRestrict ****** ;


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
;	Keywords:
;
;		Axis	If no argument is provided, then this keyword may be 
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
; Written by W.M. Nevins
;	8/17/00

;
; First find axis identifier
;
nDims = self -> NumDims()
iaxis = -1
IF(N_PARAMS() EQ 1) THEN BEGIN			; Take argument as axis identifier if present
	iaxis = arg
ENDIF ELSE BEGIN					; Else, look for the keyword 'Axis'
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
	IF(NOT Query_Integer(iaxis)) THEN BEGIN
		MESSAGE, "Illegal axis identifier", /INFORMATIONAL
		RETURN, -1
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
			self -> SignalWindow, axis=iaxis, range = axisValue
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
IF(TypeOf(result) NE 7) THEN 			$
	IF(N_ELEMENTS(result) EQ 2) THEN	$
		self -> SignalWindow, axis=iaxis, range = result
;
; Check for 'irange' on the command line
;
result = GetKeyWord('irange', extra)
	IF(TypeOf(result) NE 7) THEN 		$
		IF(N_ELEMENTS(result) EQ 2) then	$
			self -> SignalWindow, axis=iaxis, irange = LONG(result)
;
; and we're done ...
;
RETURN, iaxis

END ; ****** GKVs1D::AxisIrange ****** ;

