;
; Purpose:
;
; Computes weighting function for "PROJECT" function.
; Result is returned as a GKV object.  Weighting function
; has the form:
;
;	w = f_0 EXP{ - |(x-center[0])/width[0]|^power 
;		     - |(y-center[1])/width[1]|^power }
;
; and is normalized such that the integral of w =1.
;
; Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
; Input Keywords:
;
;	iaxis	Axis which you are projecting along.
;		must be an a positive integer
;		less than 4.  Defaluts to 3.
;		(Optional)
;
;	Center	Center of weighting function. This
;		can be a scalar or a 2-element array.
;		Defaults to center of grid. (Optional)
;
;	width	1/e half-width of weighting function. 
;		This can be a scalar or a 2-element 
;		array. Defaults to geometric mean of
;		grid spacing and grid range. (Optional)
;
;	power	Exponent controlling rate of 
;		fall-off of weighting function.
;		This can be a scalar or a 2-element 
;		array. Defaults to 2s.0. (Optional)
;
; Written by W.M. Nevins
;	4/16/04


FUNCTION GKVs3D::ProjectWeight1, arg, _Extra=extra
;
;
; Purpose:
;
;	This FUNCTION method computes ProjectWeight 
;	for projections to iaxis = 1
;
; Written by W.M. Nevins
;	4/16/04
;
;
; Get "template" -- a 2-D object at fixed value of 
; first independent variable.
;
iMin = self.Grid1.irange[0]
iMax = self.Grid1.irange[1]
template = self -> slice(axis=1, index=iMin)
;
; Now use ProjectWeight3 to compute weight object
;
weight = template -> ProjectWeight3(_Extra=extra)
template -> Trash
RETURN, weight
END  ;  ****** GKVs3D::ProjectWeight1 ******  ;


FUNCTION GKVs3D::ProjectWeight2, arg, _Extra=extra
;
;
; Purpose:
;
;	This FUNCTION method computes ProjectWeight 
;	for projections to iaxis = 2
;
; Written by W.M. Nevins
;	4/16/04
;
;
; Get "template" -- a 2-D object at fixed value of 
; second independent variable.
;
iMin = self.Grid2.irange[0]
iMax = self.Grid2.irange[1]
template = self -> slice(axis=2, index=iMin)
;
; Now use ProjectWeight3 to compute weight object
;
weight = template -> ProjectWeight3(_Extra=extra)
template -> Trash
RETURN, weight
END  ;  ****** GKVs3D::ProjectWeight2 ******  ;


FUNCTION GKVs2D::ProjectWeight3, arg, _Extra=extra
;
; Purpose:
;
;	This FUNCTION method computes ProjectWeight 
;	for projections to iaxis = 3
;
; Written by W.M. Nevins
;	4/16/04
;
; Get "template" -- a 2-D object at fixed value of 
; third independent variable.
;
nDims = self -> NumDims()
CASE nDims OF
	2   :	BEGIN
		template = self -> MakeCopy()
		END
	3   :	BEGIN
		iMin = self.Grid3.irange[0]
		iMax = self.Grid3.irange[1]
		template = self -> slice(axis=3, index=iMax)
		END
	else :	BEGIN
		MESSAGE, "Only implimented for 2- and 3-D Objects", /INFORMATIONAL
		RETURN, 0
		END
ENDCASE
;
; Get independent variables and compute center position
; on their grids.
;
x = *template.grid1.values
y = *template.grid2.values
;
; Compute default center
;
xRange = template.grid1.range
xCenter = TOTAL(xRange)/2.
yRange = template.grid2.range
yCenter = TOTAL(yrange)/2.
;
; Check command line for user-supplied "center"
;
result = GetKeyWord('center', extra)
IF( Query_Real(result) OR Query_Integer(result) ) THEN BEGIN
	center = FLOAT(result)
	CASE N_ELEMENTS(center) OF
		1   :	BEGIN
			xrange = center
			yrange = center
			END
		2   :	BEGIN
			xrange = center[0]
			yrange = center[2]
			END
		else :	BEGIN
			MESSAGE, "User-Supplied CENTER has more than 2 elements, " +	$
				 "Will center weights on center of grid", /INFORMATIONAL
			END
	ENDCASE
ENDIF
;
; Compute default widths
;
dx = x[1] - x[0]
Lx = xRange[1] - xRange[0]
xWidth = SQRT(dx*Lx)

dy = y[1] - y[0]
Ly = yRange[1] - yRange[0]
yWidth = SQRT(dy*Ly)
;
; Check command line for user-supplied "width"
;
result = GetKeyWord('width', extra)
IF( Query_Real(result) OR Query_Integer(result) ) THEN BEGIN
	width = FLOAT(result)
	CASE N_ELEMENTS(width) OF
		1   :	BEGIN
			xWidth = width
			yWidth = width
			END
		2   :	BEGIN
			xWidth = width[0]
			yWidth = width[1]
			END
		else :	BEGIN
			MESSAGE, "User-Supplied WIDTH has more than 2 elements, " +	$
				 "Will use default values", /INFORMATIONAL
			END
	ENDCASE
ENDIF
;
; Set default values of 'Power'
;
xPower = 2.
yPower = 2.
;
; Check command line for user-supplied "power"
;
result = GetKeyWord('power', extra)
IF( Query_Real(result) OR Query_Integer(result) ) THEN BEGIN
	power = FLOAT(result)
	CASE N_ELEMENTS(power) OF
		1   :	BEGIN
			xPower = power
			yPower = power
			END
		2   :	BEGIN
			xPower = power[0]
			yPower = power[1]
			END
		else :	BEGIN
			MESSAGE, "User-Supplied POWER has more than 2 elements, " +	$
				 "Will use default value of 2.0", /INFORMATIONAL
			END
	ENDCASE
ENDIF
;
; Compute weighting function
;
xArg = ABS((x-xcenter)/xWidth)^xPower
yArg = ABS((y-ycenter)/yWidth)^yPower
weights = EXP(-xArg)#EXP(-yArg)
;
; Convert 'template' into output object
;
PTR_FREE, template.values
output = template ; new name makes coding more transparant
output.values = PTR_NEW(weights)
;
; Normalize output
;
temp = output -> INT(1)
norm = temp   -> INT(1)
normV = norm -> GetValues()
output = output -> OVER(normV)
;
; clean up
;
template -> trash
temp     -> trash
norm     -> trash
;
; Set output fields for title, etc.
;
output.mnemonic = 'Weight'
output.title = 'W'
PTR_FREE, output.indices
output.Indices = PTR_NEW(['*','*'])
output.units = '1/(' + output.grid1.units + ')(' + output.grid2.units + ')'
vMax = MAX(*output.values)
output.vrange = [0,vMax]

RETURN, output
END  ;  ****** GKVs3D::ProjectWeight3 ******  ;


FUNCTION GKVs3D::ProjectWeight, arg, _Extra=extra
;
; Purpose:
;
; Computes weighting function for "PROJECT" function.
; Result is returned as a GKV object.  Weighting function
; has the form:
;
;	w = f_0 EXP{ - |(x-center[0])/width[0]|^power 
;		     - |(y-center[1])/width[1]|^power }
;
; and is normalized such that the integral of w =1.
;
; Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
; Input Keywords:
;
;	iaxis	Axis which you are projecting along.
;		must be an a positive integer
;		less than 4.  Defaluts to 3.
;		(Optional)
;
;	Center	Center of weighting function. This
;		can be a scalar or a 2-element array.
;		Defaults to center of grid. (Optional)
;
;	width	1/e half-width of weighting function. 
;		This can be a scalar or a 2-element 
;		array. Defaults to 1.0. (Optional)
;
;	power	Exponent controlling rate of 
;		fall-off of weighting function.
;		This can be a scalar or a 2-element 
;		array. Defaults to 2s.0. (Optional)
;
; Written by W.M. Nevins
;	4/16/04
;
; Get axis we are projecting about
;
iaxis = 3
nArgs = N_PARAMS()
IF(NArgs EQ 1) THEN iaxis = self -> AxisIrange(arg)
;
; Call ProjectWeight routine appropriate to projection axis
;
CASE iaxis OF
	1   :	result = self -> ProjectWeight1(_Extra=extra)
	2   :	result = self -> ProjectWeight2(_Extra=extra)
	3   :   result = self -> ProjectWeight3(_Extra=extra)
	else:	BEGIN
		MESSAGE, 'Inappropriate choice if projection axis, ' + iaxis, /INFORMATIONAL
		result = 0
		END
ENDCASE

return, result
END	; ****** GKVs3D::ProjectWeight ****** ;

