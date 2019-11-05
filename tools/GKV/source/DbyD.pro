;
; *****************************************************************************************************************
; ******************************************     Copyright Notice     *********************************************
; *                                                                                                               *
; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *
; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *
; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *
; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *
; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *
; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * 
; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *
; *                                                                                                               *
; *****************************************************************************************************************
;
; *****************************************************************************************************************
; **********************************************     DISCLAIMER     ***********************************************
; *                                                                                                               *
; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *
; *  Neither the United States Government nor the University of California nor any of their employees, makes      *
; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *
; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *
; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *
; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *
; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *
; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *
; *  States Government or the University of California, and shall not be used for advertising or product          *
; *  endorsement purposes.                                                                                        *
; *                                                                                                               *
; *****************************************************************************************************************
;

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
			IF(N_ELEMENTS(axisValue) EQ 2) THEN BEGIN
				self -> SignalWindow, axis=iaxis, range = axisValue
				RETURN, iaxis
			ENDIF
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


FUNCTION GKVs1D::DbyD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the partial derivative
;			will be returned for BOTH uniform and nonuniform grids.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis	If no argument is provided, then this keyword may be 
;			used to identify independent variable with respect to 
;			which the partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the (partial) derivative 
;			is to be taken.
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
;	5/8/00
; Revised by W.M. Nevins
;	8/18/00
;
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
result = self -> MakeCopy(/noValues)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
dValues = SHIFT(values, -1) - SHIFT(values, 1)
xold = *Grid.values
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.	; centered values in case grid is not uniform
dx = SHIFT(xold, -1) - SHIFT(xold, 1)
derivative = values						; make array of correct size and type
derivative[imin:imax] = dValues[imin:imax]/dx[imin:imax]
;
; Now pick up end points
;
boundary = STRLOWCASE(Grid.boundary)
IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
IF(imin GT irange[0]) THEN BEGIN
	xnew[0] = xold[0]
	CASE boundary OF
		'periodic (open)'	:	derivative[0] = (values[1] - values[npoints-1])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
		'periodic (closed)'	:	derivative[0] = (values[1] - values[npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
		ELSE				:	derivative[0] = 2.0*(values[1] - values[0])/(xold[1]-xold[0]) - derivative[1]
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	derivative[npoints-1] = (values[0] - values[npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
		'periodic (closed)'	:	derivative[npoints-1] = (values[1] - values[npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
		ELSE				:	derivative[npoints-1] = 2.0*(values[npoints-1] - values[npoints-2])/(xold[npoints-1]-xold[npoints-2]) - derivative[npoints-2]
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

;
; Now load values, etc. into 'result'
; 
imin=irange[0]
imax=irange[1]
values = derivative[imin:imax]
result.values = PTR_NEW(values)
result.title = "!9d!X" + self.title + '/!9d!X' + Grid.title
result.mnemonic = 'd' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')'
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
Gridvalues = xnew[imin:imax]
Grid.values = PTR_NEW(GridValues)
Grid.irange = [0,imax-imin]
Grid.range  = GridValues[[0,imax-imin]]
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)

IF Query_Integer(Extra) THEN RETURN, result
result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs1D::DbyD ****** ;


FUNCTION GKVs2D::DbyD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the partial derivative
;			will be returned for BOTH uniform and nonuniform grids.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis	If no argument is provided, then this keyword may be 
;			used to identify independent variable with respect to 
;			which the partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the (partial) derivative 
;			is to be taken.
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
;
; Written by W.M. Nevins
;	5/8/00
; Revised by W.M. Nevins
;	8/18/00
;
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
result = self -> MakeCopy(/noValues)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
dx = SHIFT(xold, -1) - SHIFT(xold, 1)
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
derivative = values
CASE axis OF
	1	:	Begin
		dValues = SHIFT(values, -1, 0) - SHIFT(values, 1, 0)
		FOR i=imin, imax DO derivative[i,*] = dvalues[i, *]/dx[i]
			END
	2	:	Begin
		dValues = SHIFT(values, 0, -1) - SHIFT(values, 0, 1)
		FOR i=imin, imax DO derivative[*,i] = dvalues[*, i]/dx[i]
			END
ENDCASE
;
; Now pick up end points
;
boundary = STRLOWCASE(Grid.boundary)
IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
IF(imin GT irange[0]) THEN BEGIN
	xnew[0] = xold[0]
	IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*] = (values[1,*] - values[npoints-1,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0] = (values[*,1] - values[*,npoints-1])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*] = (values[1,*] - values[npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0] = (values[*,1] - values[*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[0,*] = 2.0*(values[1,*] - values[0,*])/(xold[1]-xold[0]) - derivative[1,*]
				2	:	derivative[*,0] = 2.0*(values[*,1] - values[*,0])/(xold[1]-xold[0]) - derivative[*,1]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*] = (values[0,*] - values[npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1] = (values[*,0] - values[*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*] = (values[1,*] - values[npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1] = (values[*,1] - values[*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*] = 2.0*(values[npoints-1,*] - values[npoints-2,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[npoints-2,*]
				2	:	derivative[*,npoints-1] = 2.0*(values[*,npoints-1] - values[*,npoints-2])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,npoints-2]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'
;
; Now load values, etc. into 'result'
; 
;
; Note that the meaning of 'irange' now 
; changes.  Above here 'irange' refers to 
; the axis over which we are taking the derivative,
; while below here it refers to axis1.
;
; Get 'irange' for each independent variable
;
irange = self.Grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.Grid2.irange
jmin = jrange[0]
jmax = jrange[1]
;
; Put derivative values on this restricted range into 'result'
;
values = derivative[imin:imax, jmin:jmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X" + self.title + '/!9d!X' + Grid.title
result.mnemonic = 'd' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')'
;
; Reload grid corresponding to axis for which partial derivative was taken
;	(and meaning of 'irange' changes back to refering to the axis     )
;	(corresponding to variable with which we have taken the derivative)
;
irange = Grid.irange
imin = irange[0]
imax = irange[1]
gridValues = xnew[imin:imax]
Grid.values = PTR_NEW(gridValues)
nlast = imax-imin
Grid.irange = [0, nlast]
Grid.range  = gridValues[[0, nlast]]
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
;
; Restrict all grids to present signal windows
;
nDims = self -> NumDims()
FOR i=1, nDims DO result -> GridRestrict, i

IF Query_Integer(Extra) THEN RETURN, result
result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs2D::DbyD ****** ;


FUNCTION GKVs3D::DbyD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the partial derivative
;			will be returned for BOTH uniform and nonuniform grids.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis	If no argument is provided, then this keyword may be 
;			used to identify independent variable with respect to 
;			which the partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the (partial) derivative 
;			is to be taken.
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
;
; Written by W.M. Nevins
;	5/8/00
; Revised by W.M. Nevins
;	8/18/00
;
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
result = self -> MakeCopy(/noValues)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
dx = SHIFT(xold, -1) - SHIFT(xold, 1)
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
derivative = values
CASE axis OF
	1	:	BEGIN
		dValues = SHIFT(values, -1, 0, 0) - SHIFT(values, 1, 0, 0)
		FOR i=imin, imax DO derivative[i,*,*] = dvalues[i,*,*]/dx[i]
			END
	2	:	BEGIN
		dValues = SHIFT(values, 0, -1, 0) - SHIFT(values, 0, 1, 0)
		FOR i=imin, imax DO derivative[*,i,*] = dvalues[*,i,*]/dx[i]
			END
	3	:	BEGIN
		dValues = SHIFT(values, 0, 0, -1) - SHIFT(values, 0, 0, 1)
		FOR i=imin, imax DO derivative[*,*,i] = dvalues[*,*,i]/dx[i]
			END
ENDCASE
;
; Now pick up end points
;
boundary = STRLOWCASE(Grid.boundary)
IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
IF(imin GT irange[0]) THEN BEGIN
	xnew[0] = xold[0]
	IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*] = (values[1,*,*] - values[npoints-1,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0,*] = (values[*,1,*] - values[*,npoints-1,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,0] = (values[*,*,1] - values[*,*,npoints-1])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*] = (values[1,*,*] - values[npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0,*] = (values[*,1,*] - values[*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,0] = (values[*,*,1] - values[*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*] = 2.0*(values[1,*,*] - values[0,*,*])/(xold[1]-xold[0]) - derivative[1,*,*]
				2	:	derivative[*,0,*] = 2.0*(values[*,1,*] - values[*,0,*])/(xold[1]-xold[0]) - derivative[*,1,*]
				3	:	derivative[*,*,0] = 2.0*(values[*,*,1] - values[*,*,0])/(xold[1]-xold[0]) - derivative[*,*,1]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*] = (values[0,*,*] - values[npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1,*] = (values[*,0,*] - values[*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,npoints-1] = (values[*,*,0] - values[*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*] = (values[1,*,*] - values[npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1,*] = (values[*,1,*] - values[*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,npoints-1] = (values[*,*,1] - values[*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*] = 2.0*(values[npoints-1,*,*] - values[npoints-2,*,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[npoints-2,*,*]
				2	:	derivative[*,npoints-1,*] = 2.0*(values[*,npoints-1,*] - values[*,npoints-2,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,npoints-2,*]
				3	:	derivative[*,*,npoints-1] = 2.0*(values[*,*,npoints-1] - values[*,*,npoints-2])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,*,npoints-2]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

;
; Now load values, etc. into 'result'
; 
;
; Note that the meaning of 'irange' now 
; changes.  Above here 'irange' refers to 
; the axis over which we are taking the derivative,
; while below here it refers to axis1.
;
; Get 'irange' for each independent variable
;
irange = self.Grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.Grid2.irange
jmin = jrange[0]
jmax = jrange[1]
krange = self.Grid3.irange
kmin = krange[0]
kmax = krange[1]
;
; Put derivative values on this restricted range into 'result'
;
values = derivative[imin:imax, jmin:jmax, kmin:kmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X" + self.title + '/!9d!X' + Grid.title
result.mnemonic = 'd' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')'
;
; Reload grid corresponding to axis for which partial derivative was taken
;	(and meaning of 'irange' changes back to refering to the axis     )
;	(corresponding to variable with which we have taken the derivative)
;
irange = Grid.irange
imin = irange[0]
imax = irange[1]
gridValues = xnew[imin:imax]
Grid.values = PTR_NEW(gridValues)
nlast = imax-imin
Grid.irange = [0, nlast]
Grid.range  = gridValues[[0, nlast]]
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
;
; Restrict all grids to present signal windows
;
nDims = self -> NumDims()
FOR i=1, nDims DO result -> GridRestrict, i

IF Query_Integer(Extra) THEN RETURN, result
result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs3D::DbyD ****** ;


FUNCTION GKVs4D::DbyD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the partial derivative
;			will be returned for BOTH uniform and nonuniform grids.
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
;	Keywords:
;
;	   Axis		If no argument is provided, then this keyword may be 
;			used to identify independent variable with respect to 
;			which the partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the (partial) derivative 
;			is to be taken.
;
;	Side Effects:
;
;			If a 'range' or 'irange' is specified on the command line 
;			(either directly, or via 'mnemonic' = ...) then the 
;			SignalWindow method will be invoked on 'self' and,
;			on return, the signal window of the selected independent
;			variable will have been modified.
;
;
; Written by W.M. Nevins
;	5/8/00
; Revised by W.M. Nevins
;	8/18/00
;
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'DbyD called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
result = self -> MakeCopy(/noValues)
axisString = STRING(axis, FORMAT='(I1)')
commandString = 'Grid = self.Grid' + axisString
ok = EXECUTE(commandString)
values = *self.values
info = SIZE(values)
npoints = info[axis]
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
dx = SHIFT(xold, -1) - SHIFT(xold, 1)
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
derivative = values
CASE axis OF
	1	:	BEGIN
		dValues = SHIFT(values, -1, 0, 0, 0) - SHIFT(values, 1, 0, 0, 0)
		FOR i=imin, imax DO derivative[i,*,*,*] = dvalues[i,*,*,*]/dx[i]
			END
	2	:	BEGIN
		dValues = SHIFT(values, 0, -1, 0, 0) - SHIFT(values, 0, 1, 0, 0)
		FOR i=imin, imax DO derivative[*,i,*,*] = dvalues[*,i,*,*]/dx[i]
			END
	3	:	BEGIN
		dValues = SHIFT(values, 0, 0, -1, 0) - SHIFT(values, 0, 0, 1, 0)
		FOR i=imin, imax DO derivative[*,*,i,*] = dvalues[*,*,i,*]/dx[i]
			END
	4	:	BEGIN
		dValues = SHIFT(values, 0, 0, 0, -1) - SHIFT(values, 0, 0, 0, 1)
		FOR i=imin, imax DO derivative[*,*,*,i] = dvalues[*,*,*,i]/dx[i]
			END
ENDCASE
;
; Now pick up end points
;
boundary = STRLOWCASE(Grid.boundary)
IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
IF(imin GT irange[0]) THEN BEGIN
	xnew[0] = xold[0]
	IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*,*] = (values[1,*,*,*] - values[npoints-1,*,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0,*,*] = (values[*,1,*,*] - values[*,npoints-1,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,0,*] = (values[*,*,1,*] - values[*,*,npoints-1,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				4	:	derivative[*,*,*,0] = (values[*,*,*,1] - values[*,*,*,npoints-1])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*,*] = (values[1,*,*,*] - values[npoints-2,*,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,0,*,*] = (values[*,1,*,*] - values[*,npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,0,*] = (values[*,*,1,*] - values[*,*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				4	:	derivative[*,*,*,0] = (values[*,*,*,1] - values[*,*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[0,*,*,*] = 2.0*(values[1,*,*,*] - values[0,*,*,*])/(xold[1]-xold[0]) - derivative[1,*,*,*]
				2	:	derivative[*,0,*,*] = 2.0*(values[*,1,*,*] - values[*,0,*,*])/(xold[1]-xold[0]) - derivative[*,1,*,*]
				3	:	derivative[*,*,0,*] = 2.0*(values[*,*,1,*] - values[*,*,0,*])/(xold[1]-xold[0]) - derivative[*,*,1,*]
				4	:	derivative[*,*,*,0] = 2.0*(values[*,*,*,1] - values[*,*,*,0])/(xold[1]-xold[0]) - derivative[*,*,*,1]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	Grid.boundary = 'open'
	CASE boundary OF
		'periodic (open)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*,*] = (values[0,*,*,*] - values[npoints-2,*,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1,*,*] = (values[*,0,*,*] - values[*,npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,npoints-1,*] = (values[*,*,0,*] - values[*,*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				4	:	derivative[*,*,*,npoints-1] = (values[*,*,*,0] - values[*,*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*,*] = (values[1,*,*,*] - values[npoints-2,*,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				2	:	derivative[*,npoints-1,*,*] = (values[*,1,*,*] - values[*,npoints-2,*,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				3	:	derivative[*,*,npoints-1,*] = (values[*,*,1,*] - values[*,*,npoints-2,*])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
				4	:	derivative[*,*,*,npoints-1] = (values[*,*,*,1] - values[*,*,*,npoints-2])/((xold[1]-xold[0]) + (xold[npoints-1] - xold[npoints-2]))
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	derivative[npoints-1,*,*,*] = 2.0*(values[npoints-1,*,*,*] - values[npoints-2,*,*,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[npoints-2,*,*,*]
				2	:	derivative[*,npoints-1,*,*] = 2.0*(values[*,npoints-1,*,*] - values[*,npoints-2,*,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,npoints-2,*,*]
				3	:	derivative[*,*,npoints-1,*] = 2.0*(values[*,*,npoints-1,*] - values[*,*,npoints-2,*])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,*,npoints-2,*]
				4	:	derivative[*,*,*,npoints-1] = 2.0*(values[*,*,*,npoints-1] - values[*,*,*,npoints-2])/(xold[npoints-1]-xold[npoints-2]) - derivative[*,*,*,npoints-2]
			ENDCASE
							END
	ENDCASE
ENDIF
;
; Now load values, etc. into 'result'
; 
;
; Note that the meaning of 'irange' now 
; changes.  Above here 'irange' refers to 
; the axis over which we are taking the derivative,
; while below here it refers to axis1.
;
; Get 'irange' for each independent variable
;
irange = self.Grid1.irange
imin = irange[0]
imax = irange[1]
jrange = self.Grid2.irange
jmin = jrange[0]
jmax = jrange[1]
krange = self.Grid3.irange
kmin = krange[0]
kmax = krange[1]
lrange = self.Grid4.irange
lmin = lrange[0]
lmax = lrange[1]
;
; Put derivative values on this restricted range into 'result'
;
values = derivative[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X" + self.title + '/!9d!X' + Grid.title
result.mnemonic = 'd' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')'
;
; Reload grid corresponding to axis for which partial derivative was taken
;	(and meaning of 'irange' changes back to refering to the axis     )
;	(corresponding to variable with which we have taken the derivative)
;
irange = Grid.irange
imin = irange[0]
imax = irange[1]
gridValues = xnew[imin:imax]
Grid.values = PTR_NEW(gridValues)
nlast = imax-imin
Grid.irange = [0, nlast]
Grid.range  = gridValues[[0, nlast]]
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
;
; Restrict all grids to present signal windows
nDims = self -> NumDims()
FOR i=1, nDims DO result -> GridRestrict, i

IF Query_Integer(Extra) THEN RETURN, result
result -> Set, _EXTRA=Extra

RETURN, result
END ; ****** GKVs4D::DbyD ****** ;
