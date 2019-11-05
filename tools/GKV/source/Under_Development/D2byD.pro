FUNCTION GKVs1D::D2byD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the 2nd partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the 2nd partial derivative
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
;			which the 2nd partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the 2nd derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the 2nd (partial) derivative 
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
d2Values = SHIFT(values, -1) - 2.*values + SHIFT(values, 1)
xold = *Grid.values
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.	; centered values in case grid is not uniform
dxPlus = SHIFT(xold, -1) - xold
dxMinus= xold - SHIFT(xold, 1)
dxSQ = dxPlus*dxMinus						; dxSq and dDxdx are introduced so that the algorithm works for a non-uniform grid.
dDxdx= (dxPlus - dxMinus)/(dxPlus + dxMinus)	; If grid is uniform, then dxSq = dx^2, and dDxdx=0.  This algorithm produces second-order
d2 = values							; accurate 1st and second derivatives for both uniform and non-uniform grids.
d2[imin:imax] = d2Values[imin:imax]/dxSq[imin:imax] - (dvalues[imin:imax]/dxSq[imin:imax])*dDxdx[imin:imax]
;
; Now pick up end points
;
boundary = STRLOWCASE(Grid.boundary)
IF(boundary EQ 'periodic') THEN boundary = 'periodic (open)'
IF(imin GT irange[0]) THEN BEGIN
	xnew[0] = xold[0]
	CASE boundary OF
		'periodic (open)'	:	BEGIN														; Can't really know separation between separation between last value
			dxSqLeft = ((xold[1]-xold[0]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			d2[0] = (values[1] - 2.*values[0] + values[npoints-1])/dxSqLeft							; BC's, so just doing our best by assuming that Ddxdx is smooth.
							END
		'periodic (closed)'	:	BEGIN
			dxSqLeft  =  (xold[1] - xold[0]) * (xold[npoints-1] - xold[npoints-2])
			dDxdxLeft = ((xold[1] - xold[0]) - (xold[npoints-1] - xold[npoints-2]))/((xold[1] - xold[0]) + (xold[npoints-1] + xold[npoints-2]))
			d2[0] = (values[1] - 2.*values[0] + values[npoints-2])/dxSqLeft	- (values[1] - values[npoints-2])/dxSqLeft * dDxdxLeft
							END
		ELSE				:	d2[0] = d2[1]
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative object is 'open'if we reached here because the new grid does not extend to the grid boundary of 'self'.

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN																		; Can't really know separation between separation between last value
			dxSqRight = ((xold[npoints-1] - xold[npoints-2]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			d2[npoints-1] = (values[0] - 2.*values[npoints-1] + values[npoints-2])/dxSqRight							; BC's, so just doing our best by assuming that Ddxdx is smooth.
							END
		'periodic (closed)'	:	d2[npoints-1] = d2[0]
		ELSE				:	d2[npoints-1] = d2[npoints-2]
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative object is 'open'if we reached here because the new grid does not extend to the grid boundary of 'self'.

;
; Now load values, etc. into 'result'
; 
imin=irange[0]
imax=irange[1]
values = d2[imin:imax]
result.values = PTR_NEW(values)
result.title = "!9d!X!U2!N" + self.title + '/!9d!X' + Grid.title + '!U2!N'
result.mnemonic = 'd2' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')!U2!N'
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
Gridvalues = xnew[imin:imax]
Grid.values = PTR_NEW(GridValues)
Grid.irange = [0,imax-imin]
Grid.range  = GridValues[[0,imax-imin]]
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
RETURN, result
END ; ****** GKVs1D::D2byD ****** ;


FUNCTION GKVs2D::D2byD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the 2nd partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the 2nd partial derivative
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
;			which the 2nd partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the 2nd derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the 2nd (partial) derivative 
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
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
dxPlus = SHIFT(xold, -1) - xold
dxMinus= xold - SHIFT(xold, 1)
dxSQ = dxPlus*dxMinus						; dxSq and dDxdx are introduced so that the algorithm works for a non-uniform grid.
dDxdx= (dxPlus - dxMinus)/(dxPlus + dxMinus)	; If grid is uniform, then dxSq = dx^2, and dDxdx=0.  This algorithm produces second-order
d2 = values							; accurate 1st and second derivatives for both uniform and non-uniform grids.
CASE axis OF
	1	:	Begin
		dValues = SHIFT(values, -1, 0) - SHIFT(values, 1, 0)
		d2Values= SHIFT(values, -1, 0) - 2.*values + SHIFT(values, 1, 0) 
		FOR i=imin, imax DO d2[i,*] = d2Values[i,*]/dxSQ[i] - (dvalues[i,*]/dxSQ[i])*dDxdx[i]
			END
	2	:	Begin
		dValues = SHIFT(values, 0, -1) - SHIFT(values, 0, 1)
		d2Values= SHIFT(values, 0, -1) - 2.*values + SHIFT(values, 0, 1) 
		FOR i=imin, imax DO d2[*,i] = d2Values[*,i]/dxSQ[i] - (dvalues[*,i]/dxSQ[i])*dDxdx[i]
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
		'periodic (open)'	:	BEGIN														; Can't really know separation between separation between last value
			dxSqLeft = ((xold[1]-xold[0]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[0,*] = (values[1,*] - 2.*values[0,*] + values[npoints-1,*])/dxSqLeft
				2	:	d2[*,0] = (values[*,1] - 2.*values[*,0] + values[*,npoints-1])/dxSqLeft
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			dxSqLeft  =  (xold[1] - xold[0]) * (xold[npoints-1] - xold[npoints-2])
			dDxdxLeft = ((xold[1] - xold[0]) - (xold[npoints-1] - xold[npoints-2]))/((xold[1] - xold[0]) + (xold[npoints-1] + xold[npoints-2]))
			CASE axis OF
				1	:	d2[0,*] = (values[1,*] - 2.*values[0,*] + values[npoints-2,*])/dxSqLeft - (values[1,*] - values[npoints-2,*])/dxSqLeft * dDxdxLeft				
				2	:	d2[*,0] = (values[*,1] - 2.*values[*,0] + values[*,npoints-2])/dxSqLeft - (values[*,1] - values[*,npoints-2])/dxSqLeft * dDxdxLeft	
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[0,*] = d2[1,*]
				2	:	d2[*,0] = d2[*,1]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative is 'open'if we reached here because the new grid does not extend to the old boundary.


IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN																		; Can't really know separation between separation between last value
			dxSqRight = ((xold[npoints-1] - xold[npoints-2]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																				; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[npoints-1,*] = (values[0,*] - 2.*values[npoints-1,*] + values[npoints-2,*])/dxSqRight
				2	:	d2[*,npoints-1] = (values[*,0] - 2.*values[*,npoints-1] + values[*,npoints-2])/dxSqRight
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*] = d2[0,*]
				2	:	d2[*,npoints-1] = d2[*,0]
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*] = d2[npoints-2,*]
				2	:	d2[*,npoints-1] = d2[*,npoints-2]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative is 'open'if we reached here because the new grid does not extend to the old boundary.
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
values = d2[imin:imax, jmin:jmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X!U2!N" + self.title + '/!9d!X' + Grid.title + '!U2!N'
result.mnemonic = 'd2' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')!U2!N'
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
RETURN, result
END ; ****** GKVs2D::D2byD ****** ;


FUNCTION GKVs3D::D2byD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the 2nd partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the 2nd partial derivative
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
;			which the 2nd partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the 2nd derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the 2nd (partial) derivative 
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
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
dxPlus = SHIFT(xold, -1) - xold
dxMinus= xold - SHIFT(xold, 1)
dxSQ = dxPlus*dxMinus						; dxSq and dDxdx are introduced so that the algorithm works for a non-uniform grid.
dDxdx= (dxPlus - dxMinus)/(dxPlus + dxMinus)	; If grid is uniform, then dxSq = dx^2, and dDxdx=0.  This algorithm produces second-order
d2 = values							; accurate 1st and second derivatives for both uniform and non-uniform grids.
CASE axis OF
	1	:	BEGIN
		dValues = SHIFT(values, -1, 0, 0) - SHIFT(values, 1, 0, 0)
		d2Values= SHIFT(values, -1, 0, 0) - 2.*values + SHIFT(values, 1, 0, 0) 
		FOR i=imin, imax DO d2[i,*,*] = d2Values[i,*,*]/dxSq[i] - (dvalues[i,*,*]/dxSq[i])*dDxdx[i]
			END
	2	:	BEGIN
		dValues = SHIFT(values, 0, -1, 0) - SHIFT(values, 0, 1, 0)
		d2Values= SHIFT(values, 0, -1, 0) - 2.*values + SHIFT(values, 0, 1, 0) 
		FOR i=imin, imax DO d2[*,i,*] = d2Values[*,i,*]/dxSq[i] - (dvalues[*,i,*]/dxSq[i])*dDxdx[i]
			END
	3	:	BEGIN
		dValues = SHIFT(values, 0, 0, -1) - SHIFT(values, 0, 0, 1)
		d2Values= SHIFT(values, 0, 0, -1) - 2.*values + SHIFT(values, 0, 0, 1) 
		FOR i=imin, imax DO d2[*,*,i] = d2Values[*,*,i]/dxSq[i] - (dvalues[*,*,i]/dxSq[i])*dDxdx[i]
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
		'periodic (open)'	:	BEGIN														; Can't really know separation between separation between last value
			dxSqLeft = ((xold[1]-xold[0]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[0,*,*] = (values[1,*,*] - 2.*values[0,*,*] + values[npoints-1,*,*])/dxSqLeft
				2	:	d2[*,0,*] = (values[*,1,*] - 2.*values[*,0,*] + values[*,npoints-1,*])/dxSqLeft
				3	:	d2[*,*,0] = (values[*,*,1] - 2.*values[*,*,0] + values[*,*,npoints-1])/dxSqLeft
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			dxSqLeft  =  (xold[1] - xold[0]) * (xold[npoints-1] - xold[npoints-2])
			dDxdxLeft = ((xold[1] - xold[0]) - (xold[npoints-1] - xold[npoints-2]))/((xold[1] - xold[0]) + (xold[npoints-1] + xold[npoints-2]))
			CASE axis OF
				1	:	d2[0,*,*] = (values[1,*,*] - 2.*values[0,*,*] + values[npoints-2,*,*])/dxSqLeft - (values[1,*,*] - values[npoints-2,*,*])/dxSqLeft * dDxdxLeft				
				2	:	d2[*,0,*] = (values[*,1,*] - 2.*values[*,0,*] + values[*,npoints-2,*])/dxSqLeft - (values[*,1,*] - values[*,npoints-2,*])/dxSqLeft * dDxdxLeft
				3	:	d2[*,*,0] = (values[*,*,1] - 2.*values[*,*,0] + values[*,*,npoints-2])/dxSqLeft - (values[*,*,1] - values[*,*,npoints-2])/dxSqLeft * dDxdxLeft
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[0,*,*] = d2[1,*,*]			; At this level of accuracy the 2nd derivative is constant within a grid cell 
				2	:	d2[*,0,*] = d2[*,1,*]			; (a local parabolic approximation is the basis of algorithm used for interior points)
				3	:	d2[*,*,0] = d2[*,*,1]			; So this is the best we can do without any boundary information.
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative object is 'open'if we reached here because the new grid does not extend to the grid boundary of 'self'.

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN																		; Can't really know separation between separation between last value
			dxSqRight = ((xold[npoints-1] - xold[npoints-2]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																				; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[npoints-1,*,*] = (values[0,*,*] - 2.*values[npoints-1,*,*] + values[npoints-2,*,*])/dxSqRight
				2	:	d2[*,npoints-1,*] = (values[*,0,*] - 2.*values[*,npoints-1,*] + values[*,npoints-2,*])/dxSqRight
				3	:	d2[*,*,npoints-1] = (values[*,*,0] - 2.*values[*,*,npoints-1] + values[*,*,npoints-2])/dxSqRight
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*,*] = d2[0,*,*]
				2	:	d2[*,npoints-1,*] = d2[*,0,*]
				3	:	d2[*,*,npoints-1] = d2[*,*,0]
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*,*] = d2[npoints-2,*,*]
				2	:	d2[*,npoints-1,*] = d2[*,npoints-2,*]
				3	:	d2[*,*,npoints-1] = d2[*,*,npoints-2]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative is 'open'if we reached here because the new grid does not extend to the old boundary.

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
values = d2[imin:imax, jmin:jmax, kmin:kmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X!U2!N" + self.title + '/!9d!X' + Grid.title + '!U2!N'
result.mnemonic = 'd2' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')!U2!N'
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
RETURN, result
END ; ****** GKVs3D::2DbyD ****** ;


FUNCTION GKVs4D::D2byD, arg, _EXTRA=extra
;
;	Purpose:
;
; 			Returns (as a GKVsd object of same dimensionality)
; 			the 2nd partial derivative of input object with respect to
; 			the indicated axis.  
;
;			A reasonable approximation to the 2nd partial derivative
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
;			which the 2nd partial derivative is to be taken. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the independent variable
;			for which the (partial) derivative is to be taken, and to 
;			reset the signal window on this axis (before taking the 2nd derivative).
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable, NOT the integer 'irange'
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window before taking the derivative w.r.t. the selected independent variable.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable over which the 2nd (partial) derivative 
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
irange = Grid.irange
;
; First do interior points
;
imin = irange[0] > 1
imax = irange[1] < (npoints - 2)
xold = *Grid.values
xnew = (SHIFT(xold, 1) + SHIFT(xold, -1))/2.
dxPlus = SHIFT(xold, -1) - xold
dxMinus= xold - SHIFT(xold, 1)
dxSQ = dxPlus*dxMinus						; dxSq and dDxdx are introduced so that the algorithm works for a non-uniform grid.
dDxdx= (dxPlus - dxMinus)/(dxPlus + dxMinus)	; If grid is uniform, then dxSq = dx^2, and dDxdx=0.  This algorithm produces second-order
d2 = values							; accurate 1st and second derivatives for both uniform and non-uniform grids.
CASE axis OF
	1	:	BEGIN
		dValues = SHIFT(values, -1, 0, 0, 0) - SHIFT(values, 1, 0, 0, 0)
		d2Values= SHIFT(values, -1, 0, 0, 0) - 2.*values + SHIFT(values, 1, 0, 0, 0) 
		FOR i=imin, imax DO d2[i,*,*,*] = d2Values[i,*,*,*]/dxSq[i] - (dvalues[i,*,*,*]/dxSq[i])*dDxdx[i]
			END
	2	:	BEGIN
		dValues = SHIFT(values, 0, -1, 0, 0) - SHIFT(values, 0, 1, 0, 0)
		d2Values= SHIFT(values, 0, -1, 0, 0) - 2.*values + SHIFT(values, 0, 1, 0, 0) 
		FOR i=imin, imax DO d2[*,i,*,*] = d2Values[*,i,*,*]/dxSq[i] - (dvalues[*,i,*,*]/dxSq[i])*dDxdx[i]
			END
	3	:	BEGIN
		dValues = SHIFT(values, 0, 0, -1, 0) - SHIFT(values, 0, 0, 1, 0)
		d2Values= SHIFT(values, 0, 0, -1, 0) - 2.*values + SHIFT(values, 0, 0, 1, 0) 
		FOR i=imin, imax DO d2[*,*,i,*] = d2Values[*,*,i,*]/dxSq[i] - (dvalues[*,*,i,*]/dxSq[i])*dDxdx[i]
			END
	4	:	BEGIN
		dValues = SHIFT(values, 0, 0, 0, -1) - SHIFT(values, 0, 0, 0, 1)
		d2Values= SHIFT(values, 0, 0, 0, -1) - 2.*values + SHIFT(values, 0, 0, 0, 1) 
		FOR i=imin, imax DO d2[*,*,*,i] = d2Values[*,*,*,i]/dxSq[i] - (dvalues[*,*,*,i]/dxSq[i])*dDxdx[i]
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
		'periodic (open)'	:	BEGIN														; Can't really know separation between separation between last value
			dxSqLeft = ((xold[1]-xold[0]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[0,*,*,*] = (values[1,*,*,*] - 2.*values[0,*,*,*] + values[npoints-1,*,*,*])/dxSqLeft
				2	:	d2[*,0,*,*] = (values[*,1,*,*] - 2.*values[*,0,*,*] + values[*,npoints-1,*,*])/dxSqLeft
				3	:	d2[*,*,0,*] = (values[*,*,1,*] - 2.*values[*,*,0,*] + values[*,*,npoints-1,*])/dxSqLeft
				4	:	d2[*,*,*,0] = (values[*,*,*,1] - 2.*values[*,*,*,0] + values[*,*,*,npoints-1])/dxSqLeft
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			dxSqLeft  =  (xold[1] - xold[0]) * (xold[npoints-1] - xold[npoints-2])
			dDxdxLeft = ((xold[1] - xold[0]) - (xold[npoints-1] - xold[npoints-2]))/((xold[1] - xold[0]) + (xold[npoints-1] + xold[npoints-2]))
			CASE axis OF
				1	:	d2[0,*,*,*] = (values[1,*,*,*] - 2.*values[0,*,*,*] + values[npoints-2,*,*,*])/dxSqLeft - (values[1,*,*,*] - values[npoints-2,*,*,*])/dxSqLeft * dDxdxLeft				
				2	:	d2[*,0,*,*] = (values[*,1,*,*] - 2.*values[*,0,*,*] + values[*,npoints-2,*,*])/dxSqLeft - (values[*,1,*,*] - values[*,npoints-2,*,*])/dxSqLeft * dDxdxLeft
				3	:	d2[*,*,0,*] = (values[*,*,1,*] - 2.*values[*,*,0,*] + values[*,*,npoints-2,*])/dxSqLeft - (values[*,*,1,*] - values[*,*,npoints-2,*])/dxSqLeft * dDxdxLeft
				4	:	d2[*,*,*,0] = (values[*,*,*,1] - 2.*values[*,*,*,0] + values[*,*,*,npoints-2])/dxSqLeft - (values[*,*,*,1] - values[*,*,*,npoints-2])/dxSqLeft * dDxdxLeft
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[0,*,*,*] = d2[1,*,*,*]			; At this level of accuracy the 2nd derivative is constant within a grid cell 
				2	:	d2[*,0,*,*] = d2[*,1,*,*]			; (a local parabolic approximation is the basis of algorithm used for interior points)
				3	:	d2[*,*,0,*] = d2[*,*,1,*]			; So this is the best we can do without any boundary information.
				4	:	d2[*,*,*,0] = d2[*,*,1,*]			
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative object is 'open'if we reached here because the new grid does not extend to the grid boundary of 'self'.

IF(imax lt irange[1]) THEN BEGIN
	xnew[npoints-1] = xold[npoints-1]
	CASE boundary OF
		'periodic (open)'	:	BEGIN																		; Can't really know separation between separation between last value
			dxSqRight = ((xold[npoints-1] - xold[npoints-2]) * ((xold[npoints-1] - xold[npoints-2]) + (xold[1]-xold[0]))/2.)	; and first value (modulo the peridicity length) with 'open' periodic
			CASE axis OF																				; BC's, so just doing our best by assuming that Ddxdx is smooth.
				1	:	d2[npoints-1,*,*,*] = (values[0,*,*,*] - 2.*values[npoints-1,*,*,*] + values[npoints-2,*,*,*])/dxSqRight
				2	:	d2[*,npoints-1,*,*] = (values[*,0,*,*] - 2.*values[*,npoints-1,*,*] + values[*,npoints-2,*,*])/dxSqRight
				3	:	d2[*,*,npoints-1,*] = (values[*,*,0,*] - 2.*values[*,*,npoints-1,*] + values[*,*,npoints-2,*])/dxSqRight
				4	:	d2[*,*,*,npoints-1] = (values[*,*,*,0] - 2.*values[*,*,*,npoints-1] + values[*,*,*,npoints-2])/dxSqRight
			ENDCASE
							END
		'periodic (closed)'	:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*,*,*] = d2[0,*,*,*]
				2	:	d2[*,npoints-1,*,*] = d2[*,0,*,*]
				3	:	d2[*,*,npoints-1,*] = d2[*,*,0,*]
				4	:	d2[*,*,*,npoints-1] = d2[*,*,*,0]
			ENDCASE
							END
		ELSE				:	BEGIN
			CASE axis OF
				1	:	d2[npoints-1,*,*,*] = d2[npoints-2,*,*,*]
				2	:	d2[*,npoints-1,*,*] = d2[*,npoints-2,*,*]
				3	:	d2[*,*,npoints-1,*] = d2[*,*,npoints-2,*]
				4	:	d2[*,*,*,npoints-1] = d2[*,*,npoints-2,*]
			ENDCASE
							END
	ENDCASE
ENDIF ELSE Grid.boundary = 'open'	; Boundary of the second derivative is 'open'if we reached here because the new grid does not extend to the old boundary.
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
values = d2[imin:imax, jmin:jmax, kmin:kmax, lmin:lmax]
result.values = PTR_NEW(values)
;
; Set appropriate 'vrange'
;
vmin = GKVsd_Min(values, Max=vmax)
result.vrange = [vmin, vmax]
;
; set title, etc.
;
result.title = "!9d!X!U2!N" + self.title + '/!9d!X' + Grid.title + '!U2!N'
result.mnemonic = 'd2' + self.mnemonic + 'd' + Grid.mnemonic
result.units = '('+ self.units + ')/(' + Grid.units + ')!U2!N'
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
RETURN, result
END ; ****** GKVs4D::D2byD ****** ;
