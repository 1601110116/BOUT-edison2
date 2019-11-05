PRO GKVs1D::Window, arg, _Extra=Extra
;
; This proceedure returns "self" restricted
; to the range indicated, with hanning window 
; applied in the indicated variable.
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
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the independent varialbe to be windowed.
;			Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the selected independent 
;			variable, and to specify the size of window on this axis.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable (i.e., it is interpreted in the units
;			of the corresponding independent variable), NOT the integer 'irange'
;			(that is, NOT as an integer index into the grid.values array).
;
;	   irange	Set 'irange' to a two-element (integer) array to select the
;			range of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of the window in the specifiec independent variable.
;
;	      idw	The (integer) width of the roll-off at each end of the data
;			window.
;
;	       dw	The (real) width of the roll-off at each end of the dat widow.
;
;	     Norm	Set this keyword (i.e., put '/Norm' on the command line) to
;			renorm the data to (approximately) preserve the RMS value.
;
;
; Written by W.M. Nevins
;	8/17/00
;
FORWARD_FUNCTION GKVsd_MIN
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'WINDOW called with too many arguments', /INFORMATIONAL
				RETURN
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier, returning', /INFORMATIONAL
	RETURN
ENDIF
self -> Restrict

grid = self.grid1

gridValues = *(grid.values)
nVals = N_ELEMENTS(gridValues)
idw = nVals/10
result = GetKeyWord('idw', Extra)
IF(TypeOf(result) NE 7) THEN idw = result
result = GetKeyWord('dw', Extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	dx = gridValues[1] - gridValues[0]
	idw = LONG(dw/dx)
ENDIF

windowFcn = Window_Fcn(nVals, idw)

norm = 1.
result = GetKeyWord('Norm', Extra)
IF(TypeOF(result) NE 7) THEN norm = nVals/TOTAL(windowFcn)

values = norm*(*self.values)*windowFcn
PTR_FREE, self.values
self.values = PTR_NEW(values)
vmin = GKVsD_MIN(values, MAX=vmax)
self.vrange=[vmin, vmax]

RETURN
END  ; ********** PRO GKVs1D::Window **********  ;


FUNCTION GKVs1D::WINDOW, arg, _Extra=Extra
;
; function version of previous
;
result = self -> MakeCopy()
nArgs = N_PARAMS()
CASE nArgs of
	0	:	result -> Window, _Extra=Extra
	1	:	result -> Window, arg, _Extra=Extra
ENDCASE
RETURN, result
END  ;  ********** FUNCTION GKVs1D::WINDOW **********  ;

	


PRO GKVs2D::Window, arg, _Extra=Extra
;
; This proceedure returns "self" restricted
; to the range indicated, with hanning window 
; applied in the indicated variable.
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
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the independent varialbe to be windowed.
;			Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the selected independent 
;			variable, and to specify the size of window on this axis.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable (i.e., it is interpreted in the units
;			of the corresponding independent variable), NOT the integer 'irange'
;			(that is, NOT as an integer index into the grid.values array).
;
;	   irange	Set 'irange' to a two-element (integer) array to select the
;			range of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of the window in the specifiec independent variable.
;
;	      idw	The (integer) width of the roll-off at each end of the data
;			window.
;
;	       dw	The (real) width of the roll-off at each end of the dat widow.
;
;	     Norm	Set this keyword (i.e., put '/Norm' on the command line) to
;			renorm the data to (approximately) preserve the RMS value.
;
;
; Written by W.M. Nevins
;	8/17/00
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'WINDOW called with too many arguments', /INFORMATIONAL
				RETURN
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN
ENDIF
self -> Restrict

CASE axis OF
	1	:	grid = self.grid1
	2	:	grid = self.grid2
ENDCASE

gridValues = *(grid.values)
nVals = N_ELEMENTS(gridValues)
idw = nVals/10
result = GetKeyWord('idw', Extra)
IF(TypeOf(result) NE 7) THEN idw = result
result = GetKeyWord('dw', Extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	dx = gridValues[1] - gridValues[0]
	idw = LONG(dw/dx)
ENDIF

windowFcn = Window_Fcn(nVals, idw)

CASE axis OF
	1	:	nOther = N_ELEMENTS(*self.grid2.values)
	2	:	nOther = N_ELEMENTS(*self.grid1.values)
ENDCASE

CASE axis OF
	1	: 	windowFcn = windowFcn#MAKE_ARRAY(nOther, /FLOAT, VALUE=1.)
	2	:	windowFcn = MAKE_ARRAY(nOther, /FLOAT, VALUE=1.)#windowFcn
ENDCASE

norm = 1.
result = GetKeyWord('Norm', Extra)
IF(TypeOF(result) NE 7) THEN norm = nVals*nOther/TOTAL(windowFcn)

values = norm*(*self.values)*windowFcn
PTR_FREE, self.values
self.values = PTR_NEW(values)
vmin = GKVsD_MIN(values, MAX=vmax)
self.vrange=[vmin, vmax]

RETURN
END  ; ********** PRO GKVs2D::Window **********  ;


FUNCTION GKVs2D::WINDOW, arg, _Extra=Extra
;
; function version of previous
;
result = self -> MakeCopy()
nArgs = N_PARAMS()
CASE nArgs of
	0	:	result -> Window, _Extra=Extra
	1	:	result -> Window, arg, _Extra=Extra
ENDCASE
RETURN, result
END  ;  ********** FUNCTION GKVs2D::WINDOW **********  ;

	


PRO GKVs3D::Window, arg, _Extra=Extra
;
; This proceedure returns "self" restricted
; to the range indicated, with hanning window 
; applied in the indicated variable.
;
;
; 	Argument:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims-1, or
;			a STRING containing an axis mnemonic. (This proceedure
;			DOES NOT SUPPORT putting a "window" on the final
;			(generally time) dependent varialbe.
;
;	Keywords:
;
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the independent varialbe to be windowed.
;			Set axis equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two-element 
;			array, [min, max], to both identify the selected independent 
;			variable, and to specify the size of window on this axis.
;			This two-element array is interpreted as the desired RANGE in
;			the independent variable (i.e., it is interpreted in the units
;			of the corresponding independent variable), NOT the integer 'irange'
;			(that is, NOT as an integer index into the grid.values array).
;
;	   irange	Set 'irange' to a two-element (integer) array to select the
;			range of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range of the window in the specifiec independent variable.
;
;	      idw	The (integer) width of the roll-off at each end of the data
;			window.
;
;	       dw	The (real) width of the roll-off at each end of the dat widow.
;
;	     Norm	Set this keyword (i.e., put '/Norm' on the command line) to
;			renorm the data to (approximately) preserve the RMS value.
;
;
; Written by W.M. Nevins
;	8/17/00
;
; Find axis identifier
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'WINDOW called with too many arguments', /INFORMATIONAL
				RETURN
			END
ENDCASE
IF(axis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN
ENDIF
self -> Restrict

CASE axis OF
	1	:	grid = self.grid1
	2	:	grid = self.grid2
	3	:	BEGIN
				MESSAGE, 'Can not put WINDOW on third (time) independent variable', /INFORMATIONAL
				RETURN
			END
ENDCASE

gridValues = *(grid.values)
nVals = N_ELEMENTS(gridValues)
idw = nVals/10
result = GetKeyWord('idw', Extra)
IF(TypeOf(result) NE 7) THEN idw = result
result = GetKeyWord('dw', Extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	dx = gridValues[1] - gridValues[0]
	idw = LONG(dw/dx)
ENDIF

windowFcn = Window_Fcn(nVals, idw)

CASE axis OF
	1	:	nOther = N_ELEMENTS(*self.grid2.values)
	2	:	nOther = N_ELEMENTS(*self.grid1.values)
ENDCASE

CASE axis OF
	1	: 	windowFcn = windowFcn#MAKE_ARRAY(nOther, /FLOAT, VALUE=1.)
	2	:	windowFcn = MAKE_ARRAY(nOther, /FLOAT, VALUE=1.)#windowFcn
ENDCASE

norm = 1.
result = GetKeyWord('Norm', Extra)
IF(TypeOF(result) NE 7) THEN norm = nVals*nOther/TOTAL(windowFcn)

tGrid = self.grid3
tValues = *(tGrid.values)
nt = N_ELEMENTS(tValues)
FOR i=0, nt-1 DO BEGIN
	(*(self.values))[*,*,i] = (*(self.values))[*,*,i]*windowFcn
ENDFOR

RETURN
END  ; ********** PRO GKVs3D::Window **********  ;
