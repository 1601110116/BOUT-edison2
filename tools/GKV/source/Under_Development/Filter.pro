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

FUNCTION GKVs1D::Filter, arg, _Extra=extra
;
; Purpose:
;
;		This function returns a GKVs1D object containing
;		data from 'self' passed through a band-pass filter.
;
;
; Arguments:
;
;			Any legal axis identifier. Defaults to the final axis.
;			The independent variable may also be identified using an
;			Axis mnemonic. (Optional)
;
;
; Keywords:
;
;   'mnemonic'	Where 'mnemonic' is the mnemonic for the independent variable
;			over which the filter is to be applied.  'Mnemonic' should be 
;			set equal to the range in this variable over which you wish to 
;			have filtered data returned.  Defaults to the final axis and
;			current irange. (Optional)
;
;	 Omega_0	Central frequency of the digital filter to be applied. Defaults
;			to zero. (Optional)
;
;	    k_0	Central wavenumber of the digital filter to be applied.  Defaults
;			to zero.  This keyword is really a synonym for 'Omega_0', and is
;			included to maintain an intuitive interface when the axis in 
;			question corresponds to a spatial (rather than time-like) variable.
;			(Optional).
;
;	dOmega	Width of the filter in frequency.  Defaults to 10/T, where 'T' is
;			the length of the time interval over which data is available from
;			'self'. (Optional)
;
;	dk		A synonym for dOmega--the width of the filter in wavenumber space
;			which defaults to 10/L, where 'L is the length of the interval over
;			which data is available from 'self'.  (Optional)
;
;	dT		An alternative means of specifying the filter width--if set, dT
;			gives the width of the support of the digital filter in the time-domain,
;			corresponding to a frequency width of dOmega=2¼/dT.
;
;	dL		A snyonym for dT--the width of the support of the digital filter in the 
;			spatial-domain, corresponding to a frequency width of dk=2¼/dT.	
;
;
; Written by W.M. Nevins
;	10/3/00
;
;
; Find axis identifier from command line 
; (this also resets the signal window if appropriate
;  keywords have been included on the command line)
;
;
; Make copy of 'self'
;
output = self -> MakeCopy()
CASE N_PARAMS() OF
	0	:	iaxis = output -> AxisIrange(     _Extra=extra)
	1	:	iaxis = output -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Filter called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LT 0) THEN iaxis=1
output -> Restrict
;
; Get grid structure
;
axisString = STRING(iaxis, FORMAT='(i1)')
commandString = 'Grid = output.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *Grid.values
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
;
; Check that grid is uniform over specified interval
;
IF(NOT Grid.uniform) THEN BEGIN
	MESSAGE, "Grid is not uniform over specified interval.", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get time step (deltaT) and total length of signal (Ttot)
; 
deltaT = gridValues[1] - gridValues[0]
Ttot = gridValues[imax] - gridvalues[imin]
;
; set up 'temp' array to process data from 'self' (which is now in output...)
;
periodic = 0b
IF(grid.boundary EQ "periodic (open)") THEN BEGIN
	Ttot = Ttot + deltaT
	periodic = 1b
ENDIF
nPoints = imax - imin + 1
IF(grid.boundary EQ "periodic (closed)") THEN BEGIN
	nPoints = nPoints-1
	periodic = 1b
ENDIF
IF(periodic) THEN BEGIN
	ntemp = npoints
	GOTO, DONE1
ENDIF
FOR i=0,24 DO BEGIN
	ntemp = 2^i
	IF(nTemp GE 2*nPoints) THEN GOTO, DONE1
ENDFOR
MESSAGE, "too many elements in input array", /INFORMATIONAL
RETURN, 0

DONE1:
valuePtr = output -> GetValues()
values = *valuePtr
valueType = TypeOf(values)
temp = MAKE_ARRAY(ntemp, TYPE=valueType)
temp[0:(npoints-1)] = values[0:(npoints-1)]
;
; Get central frequency and width
;
Omega_0 = 0
dT=Ttot/10.
sgn = -1.
result = GetKeyWord("Omega_0", extra)
IF(TypeOf(result) NE 7) THEN Omega_0 = result
result = GetKeyWord("k_0", extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	Omega_0 = result
	sgn = 1.
ENDIF
result = GetKeyWord("dOmega", extra)
IF(TypeOf(result) NE 7) THEN dT = 2*!Pi/result
result = GetKeyWord("dk", extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	dT = 2*!Pi/result
	sgn = 1.
ENDIF
result = GetKeyWord("dT", extra)
IF(TypeOf(result) NE 7) THEN dT = result
result = GetKeyWord("dL", extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	dT = result
	sgn = 1.
ENDIF
;
; Create digital filter
;
filterLength = FIX(dT/deltaT+0.5) < nPoints
filterLength = filterLength > 2
filter = COMPLEXARR(nTemp)
cosineWave = EXP( COMPLEX(0., sgn)*Omega_0*deltaT*(FINDGEN(filterLength) - (filterLength/2)) )
filterEnvelope = MyHanning(filterLength)
filterNorm = 1./TOTAL(filterEnvelope)
filter[0:filterLength-1] = filterEnvelope*cosineWave
filter = filterNorm*SHIFT( filter, -(filterLength+1)/2 )
;
; Perform a convolution between data in 'temp' and 'filter' using FFT's
;
temp = FFT(temp, 1, /OVERWRITE)
filter = FFT(filter, 1, /OVERWRITE)
temp = temp*filter
temp = FFT(temp, -1, /OVERWRITE)
IF(Grid.boundary EQ "periodic (closed)") THEN temp = [temp, temp[0]]
valuesOut = temp[0:(npoints-1)]
IF(NOT Query_Complex(values) ) THEN valuesOut = FLOAT(valuesOut)
;
; Prepare output object
;
PTR_FREE, output.values
output.values = PTR_NEW(valuesOut)
vmin = GKVsd_MIN(valuesOut, Max=vmax)
output.vrange = [vmin, vmax]
;
; The filtered array has a different time/space grid
; (same spacing, but possibly different length)
;
IF(NOT periodic) THEN BEGIN
	newGridValues = gridValues[0] + deltaT*FINDGEN(nTemp)
	PTR_FREE, Grid.values
	Grid.values = PTR_NEW(newGridValues)
	commandString = 'output.Grid' + axisString + ' = Grid'
	ok = EXECUTE(commandString)
	output -> signalwindow, t=[gridValues[imin], gridValues[imax]]
ENDIF
output -> restrict
;
; and we're done
;
RETURN, output
END ; ****** GKVs1D::Filter ****** ;


FUNCTION GKVs2D::Filter, arg, _Extra=extra
;
; Purpose:
;
;		This function returns a GKVs2D object containing
;		data from 'self' passed through a band-pass filter.
;
;
; Arguments:
;
;			Any legal axis identifier. Defaults to the final (time?) axis.
;			The independent variable over which the filter will be applied
;			may also be identified using an Axis mnemonic (see Keywords below). 
;			(Optional)
;
;
; Keywords:
;
;   'mnemonic'	Where 'mnemonic' is the mnemonic for the independent variable
;			over which the filter is to be applied.  'Mnemonic' should be 
;			set equal to the range in this variable over which you wish to 
;			have filtered data returned.  Defaults to the final axis and
;			current irange. (Optional)
;
;	 Omega_0	Central frequency of the digital filter to be applied. Defaults
;			to zero. (Optional)
;
;	    k_0	Central wavenumber of the digital filter to be applied.  Defaults
;			to zero.  This keyword is really a synonym for 'Omega_0', and is
;			included to maintain an intuitive interface when the axis in 
;			question corresponds to a spatial (rather than time-like) variable.
;			(Optional).
;
;	dOmega	Width of the filter in frequency.  Defaults to 10/T, where 'T' is
;			the length of the time interval over which data is available from
;			'self'. (Optional)
;
;	dk		A synonym for dOmega--the width of the filter in wavenumber space
;			which defaults to 10/L, where 'L is the length of the interval over
;			which data is available from 'self'.  (Optional)
;
;	dT		An alternative means of specifying the filter width--if set, dT
;			gives the width of the support of the digital filter in the time-domain,
;			corresponding to a frequency width of dOmega=2¼/dT.
;
;	dL		A snyonym for dT--the width of the support of the digital filter in the 
;			spatial-domain, corresponding to a frequency width of dk=2¼/dT.	
;
;
; Written by W.M. Nevins
;	10/3/00
;
; Make copy of self
;
result = self -> MakeCopy()
;
; Find axis identifier from command line 
; (this also resets the signal window if appropriate
;  keywords have been included on the command line)
;
CASE N_PARAMS() OF
	0	:	iaxis = result -> AxisIrange(     _Extra=extra)
	1	:	iaxis = result -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Filter called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LT 0) THEN iaxis=2
IF(iaxis GT 2) THEN BEGIN
	MESSAGE, "Couldn't successfully parse command line for axis identifer", /INFORMATIONAL
	RETURN, 0
ENDIF
result -> Restrict
;
; Get grid structure
;
axisString = STRING(iaxis, FORMAT='(i1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *Grid.values
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
;
; Check that grid is uniform over specified interval
;
IF(NOT Grid.uniform) THEN BEGIN
	MESSAGE, "Grid is not uniform over specified interval.", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get time step (deltaT) and total length of signal (Ttot)
; 
deltaT = gridValues[1] - gridValues[0]
Ttot = gridValues[imax] - gridvalues[imin]
;
; set up 'temp' array to process data from 'self' (which has been copied into result...)
;
periodic = 0b
IF(grid.boundary EQ "periodic (open)") THEN BEGIN
	Ttot = Ttot + deltaT
	periodic = 1b
ENDIF
nPoints = imax - imin + 1
IF(grid.boundary EQ "periodic (closed)") THEN BEGIN
	nPoints = nPoints-1
	periodic = 1b
ENDIF
IF(periodic) THEN BEGIN
	ntemp = npoints
	GOTO, DONE1
ENDIF
FOR i=0,24 DO BEGIN
	ntemp = 2^i
	IF(nTemp GE 2*nPoints) THEN GOTO, DONE1
ENDFOR
MESSAGE, "too many elements in input array", /INFORMATIONAL
RETURN, 0

DONE1:
valuePtr = result -> GetValues()
values = *valuePtr
info = SIZE(values)
valueType = TypeOf(values)
CASE iaxis OF
	1:	BEGIN
			temp = MAKE_ARRAY(ntemp, info[2], TYPE=valueType)
			temp[0:(npoints-1), *] = values[0:(npoints-1), *]
		END
	2:	BEGIN
			temp = MAKE_ARRAY(info[1], ntemp, TYPE=valueType)
			temp[ *, 0:(npoints-1)] = values[ *, 0:(npoints-1)]
		END
ENDCASE
;
; Get central frequency and width
;
Omega_0 = 0
dT=Ttot/10.
sgn = -1.
commandLineInfo = GetKeyWord("Omega_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN Omega_0 = commandLineInfo
commandLineInfo = GetKeyWord("k_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	Omega_0 = commandLineInfo
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dOmega", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = (2*!Pi/commandLineInfo) > deltaT
commandLineInfo = GetKeyWord("dk", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = (2*!Pi/commandLineInfo) > deltaT
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dT", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = commandLineInfo > deltaT
commandLineInfo = GetKeyWord("dL", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = commandLineInfo > deltaT
	sgn = 1.
ENDIF
;
; Create digital filter
;
filterLength = FIX(dT/deltaT+0.5) < nPoints
filterLength = filterLength > 2
filter = COMPLEXARR(nTemp)
cosineWave = EXP( COMPLEX(0., sgn)*Omega_0*deltaT*(FINDGEN(filterLength) - (filterLength/2)) )
filterEnvelope = MyHanning(filterLength)
filterNorm = 1./TOTAL(filterEnvelope)
filter[0:filterLength-1] = filterEnvelope*cosineWave
filter = filterNorm*SHIFT( filter, -(filterLength+1)/2 )
;
; Perform a convolution between data in 'temp' and 'filter' using FFT's
;
temp = FFT(temp, 1, /OVERWRITE)
filter = FFT(filter, 1, /OVERWRITE)
CASE iaxis OF
	1:	BEGIN
			fill = REPLICATE(1., info[2])
			filter = filter#fill
		END
	2:	BEGIN
			fill = REPLICATE(1., info[1])
			filter = fill#filter
		END
ENDCASE
temp = temp*filter
;
; Invert transform to obtain filter data
;
temp = FFT(temp, -1, /OVERWRITE)
;
; If necessary, deal with boundary conditions 
;
IF(NOT Query_Complex(values) ) THEN temp = FLOAT(temp)
IF(Grid.boundary EQ "periodic (closed)") THEN BEGIN
	valuesOut = MAKE_ARRAY(info[1], info[2], TYPE=valueType)
	CASE iaxis OF
		1:	BEGIN
				valuesOut[0:(npoints-1), *] = temp[0:(npoints-1), *]
				valuesOut[      npoints, *] = temp[0, *]
			END 
		2:	BEGIN
				valuesOut[*, 0:(npoints-1)] = temp[*, 0:(npoints-1)]
				valuesOut[*, npoints      ] = temp[*, 0]
			END 
	ENDCASE
ENDIF ELSE BEGIN
	CASE iaxis OF
		1:	BEGIN
				valuesOut = temp[0:(npoints-1), *]
			END
		2:	BEGIN
				valuesOut = temp[*, 0:(npoints-1)]
			END
	ENDCASE
ENDELSE
;
; Prepare output object
;
PTR_FREE, result.values
result.values = PTR_NEW(valuesOut)
vmin = GKVsd_MIN(valuesOut, Max=vmax)
result.vrange = [vmin, vmax]
;
; The filtered array has a different time/space grid
; (same spacing, but possibly different length)
;
IF(NOT periodic) THEN BEGIN
	newGridValues = gridValues[0] + deltaT*FINDGEN(nTemp)
	PTR_FREE, Grid.values
	Grid.values = PTR_NEW(newGridValues)
	commandString = 'result.Grid' + axisString + ' = Grid'
	ok = EXECUTE(commandString)
	result -> signalwindow, t=[gridValues[imin], gridValues[imax]]
ENDIF
result -> Restrict
;
; and we're done
;
RETURN, result
END ; ****** GKVs2D::Filter ****** ;


FUNCTION GKVs3D::Filter, arg, _Extra=extra
;
; Purpose:
;
;		This function returns a GKVs3D object containing
;		data from 'self' passed through a band-pass filter.
;
;
; Arguments:
;
;			Any legal axis identifier. Defaults to the final (time?) axis.
;			The independent variable over which the filter will be applied
;			may also be identified using an Axis mnemonic (see Keywords below). 
;			(Optional)
;
;
; Keywords:
;
;   'mnemonic'	Where 'mnemonic' is the mnemonic for the independent variable
;			over which the filter is to be applied.  'Mnemonic' should be 
;			set equal to the range in this variable over which you wish to 
;			have filtered data returned.  Defaults to the final axis and
;			current irange. (Optional)
;
;	 Omega_0	Central frequency of the digital filter to be applied. Defaults
;			to zero. (Optional)
;
;	    k_0	Central wavenumber of the digital filter to be applied.  Defaults
;			to zero.  This keyword is really a synonym for 'Omega_0', and is
;			included to maintain an intuitive interface when the axis in 
;			question corresponds to a spatial (rather than time-like) variable.
;			(Optional).
;
;	dOmega	Width of the filter in frequency.  Defaults to 10/T, where 'T' is
;			the length of the time interval over which data is available from
;			'self'. (Optional)
;
;	dk		A synonym for dOmega--the width of the filter in wavenumber space
;			which defaults to 10/L, where 'L is the length of the interval over
;			which data is available from 'self'.  (Optional)
;
;	dT		An alternative means of specifying the filter width--if set, dT
;			gives the width of the support of the digital filter in the time-domain,
;			corresponding to a frequency width of dOmega=2¼/dT.
;
;	dL		A snyonym for dT--the width of the support of the digital filter in the 
;			spatial-domain, corresponding to a frequency width of dk=2¼/dT.	
;
;
; Written by W.M. Nevins
;	10/3/00
;
; Make copy of self
;
result = self -> MakeCopy()
;
; Find axis identifier from command line 
; (this also resets the signal window if appropriate
;  keywords have been included on the command line)
;
CASE N_PARAMS() OF
	0	:	iaxis = result -> AxisIrange(     _Extra=extra)
	1	:	iaxis = result -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Filter called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LT 0) THEN iaxis=3
IF(iaxis GT 3) THEN BEGIN
	MESSAGE, "Couldn't successfully parse command line for axis identifer", /INFORMATIONAL
	RETURN, 0
ENDIF
result -> Restrict
;
; Get grid structure
;
axisString = STRING(iaxis, FORMAT='(i1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *Grid.values
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
;
; Check that grid is uniform over specified interval
;
IF(NOT Grid.uniform) THEN BEGIN
	MESSAGE, "Grid is not uniform over specified interval.", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get time step (deltaT) and total length of signal (Ttot)
; 
deltaT = gridValues[1] - gridValues[0]
Ttot = gridValues[imax] - gridvalues[imin]
;
; set up 'temp' array to process data from 'self' (which has been copied into result...)
;
periodic = 0b
IF(grid.boundary EQ "periodic (open)") THEN BEGIN
	Ttot = Ttot + deltaT
	periodic = 1b
ENDIF
nPoints = imax - imin + 1
IF(grid.boundary EQ "periodic (closed)") THEN BEGIN
	nPoints = nPoints-1
	periodic = 1b
ENDIF
IF(periodic) THEN BEGIN
	ntemp = npoints
	GOTO, DONE1
ENDIF
FOR i=0,24 DO BEGIN
	ntemp = 2^i
	IF(nTemp GE 2*nPoints) THEN GOTO, DONE1
ENDFOR
MESSAGE, "too many elements in input array", /INFORMATIONAL
RETURN, 0

DONE1:
valuePtr = result -> GetValues()
values = *valuePtr
info = SIZE(values)
valueType = TypeOf(values)
CASE iaxis OF
	1:	BEGIN
			temp = MAKE_ARRAY(ntemp, info[2], info[3], TYPE=valueType)
			temp[ 0:(npoints-1), *, *] = values[ 0:(npoints-1), *, *]
		END
	2:	BEGIN
			temp = MAKE_ARRAY(info[1], ntemp, info[3], TYPE=valueType)
			temp[ *, 0:(npoints-1), *] = values[ *, 0:(npoints-1), *]
		END
	3:	BEGIN
			temp = MAKE_ARRAY(info[1],info[2],  ntemp, TYPE=valueType)
			temp[ *, *, 0:(npoints-1)] = values[ *, *, 0:(npoints-1)]
		END
ENDCASE
;
; Get central frequency and width
;
Omega_0 = 0
dT=Ttot/10.
sgn = -1.
commandLineInfo = GetKeyWord("Omega_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN Omega_0 = commandLineInfo
commandLineInfo = GetKeyWord("k_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	Omega_0 = commandLineInfo
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dOmega", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = (2*!Pi/commandLineInfo) > deltaT
commandLineInfo = GetKeyWord("dk", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = (2*!Pi/commandLineInfo) > deltaT
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dT", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = commandLineInfo > deltaT
commandLineInfo = GetKeyWord("dL", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = commandLineInfo > deltaT
	sgn = 1.
ENDIF
;
; Create digital filter
;
filterLength = FIX(dT/deltaT+0.5) < nPoints
filterLength = filterLength > 2
filter = COMPLEXARR(nTemp)
cosineWave = EXP( COMPLEX(0., sgn)*Omega_0*deltaT*(FINDGEN(filterLength) - (filterLength/2)) )
filterEnvelope = MyHanning(filterLength)
filterNorm = 1./TOTAL(filterEnvelope)
filter[0:filterLength-1] = filterEnvelope*cosineWave
filter = filterNorm*SHIFT( filter, -(filterLength+1)/2 )
;
; Perform a convolution between data in 'temp' and 'filter' using FFT's
;
temp = FFT(temp, 1, /OVERWRITE)
filter = FFT(filter, 1, /OVERWRITE)
CASE iaxis OF
	1:	BEGIN
			fill = REPLICATE(1., info[2])
			filter = filter#fill
		END
	2:	BEGIN
			fill = REPLICATE(1., info[1])
			filter = fill#filter
		END
	3:	BEGIN
			fill = REPLICATE(1., info[2])
			filter = fill#filter
		END
ENDCASE
CASE iaxis OF
	1:	FOR i=0, info[3]-1 DO temp[*,*,i] = temp[*,*,i]*filter
	2:	FOR i=0, info[3]-1 DO temp[*,*,i] = temp[*,*,i]*filter 
	3:	FOR i=0, info[1]-1 DO temp[i,*,*] = temp[i,*,*]*filter
ENDCASE
;
; Invert transform to obtain filter data
;
temp = FFT(temp, -1, /OVERWRITE)
;
; If necessary, deal with boundary conditions 
;
IF(NOT Query_Complex(temp) ) THEN temp = FLOAT(temp)
IF(Grid.boundary EQ "periodic (closed)") THEN BEGIN
	valuesOut = MAKE_ARRAY(info[1], info[2], info[3], TYPE=valueType)
	CASE iaxis OF
		1:	BEGIN
				valuesOut[0:(npoints-1), *, *] = temp[0:(npoints-1), *, *]
				valuesOut[      npoints, *, *] = temp[0, *, *]
			END 
		2:	BEGIN
				valuesOut[*, 0:(npoints-1), *] = temp[*, 0:(npoints-1), *]
				valuesOut[*, npoints      , *] = temp[*, 0, *]
			END 
		3:	BEGIN
				valuesOut[*, *, 0:(npoints-1)] = temp[*, *, 0:(npoints-1)]
				valuesOut[*, *, npoints      ] = temp[*, *, 0]
			END 
	ENDCASE
ENDIF ELSE BEGIN
	CASE iaxis OF
		1:	BEGIN
				valuesOut = temp[0:(npoints-1), *, *]
			END
		2:	BEGIN
				valuesOut = temp[*, 0:(npoints-1), *]
			END
		3:	BEGIN
				valuesOut = temp[*, *, 0:(npoints-1)]
			END
	ENDCASE
ENDELSE
;
; Prepare output object
;
PTR_FREE, result.values
result.values = PTR_NEW(valuesOut)
vmin = GKVsd_MIN(valuesOut, Max=vmax)
result.vrange = [vmin, vmax]
;
; The filtered array has a different time/space grid
; (same spacing, but possibly different length)
;
IF(NOT periodic) THEN BEGIN
	newGridValues = gridValues[0] + deltaT*FINDGEN(nTemp)
	PTR_FREE, Grid.values
	Grid.values = PTR_NEW(newGridValues)
	commandString = 'result.Grid' + axisString + ' = Grid'
	ok = EXECUTE(commandString)
	result -> signalwindow, t=[gridValues[imin], gridValues[imax]]
ENDIF
result -> Restrict
;
; and we're done
;
RETURN, result
END ; ****** GKVs3D::Filter ****** ;


FUNCTION GKVs4D::Filter, arg, _Extra=extra
;
; Purpose:
;
;		This function returns a GKVs4D object containing
;		data from 'self' passed through a band-pass filter.
;
;
; Arguments:
;
;			Any legal axis identifier. Defaults to the final (time?) axis.
;			The independent variable over which the filter will be applied
;			may also be identified using an Axis mnemonic (see Keywords below). 
;			(Optional)
;
;
; Keywords:
;
;   'mnemonic'	Where 'mnemonic' is the mnemonic for the independent variable
;			over which the filter is to be applied.  'Mnemonic' should be 
;			set equal to the range in this variable over which you wish to 
;			have filtered data returned.  Defaults to the final axis and
;			current irange. (Optional)
;
;	 Omega_0	Central frequency of the digital filter to be applied. Defaults
;			to zero. (Optional)
;
;	    k_0	Central wavenumber of the digital filter to be applied.  Defaults
;			to zero.  This keyword is really a synonym for 'Omega_0', and is
;			included to maintain an intuitive interface when the axis in 
;			question corresponds to a spatial (rather than time-like) variable.
;			(Optional).
;
;	dOmega	Width of the filter in frequency.  Defaults to 10/T, where 'T' is
;			the length of the time interval over which data is available from
;			'self'. (Optional)
;
;	dk		A synonym for dOmega--the width of the filter in wavenumber space
;			which defaults to 10/L, where 'L is the length of the interval over
;			which data is available from 'self'.  (Optional)
;
;	dT		An alternative means of specifying the filter width--if set, dT
;			gives the width of the support of the digital filter in the time-domain,
;			corresponding to a frequency width of dOmega=2¼/dT.
;
;	dL		A snyonym for dT--the width of the support of the digital filter in the 
;			spatial-domain, corresponding to a frequency width of dk=2¼/dT.	
;
;
; Written by W.M. Nevins
;	10/3/00
;
; Make copy of self
;
result = self -> MakeCopy()
;
; Find axis identifier from command line 
; (this also resets the signal window if appropriate
;  keywords have been included on the command line)
;
CASE N_PARAMS() OF
	0	:	iaxis = result -> AxisIrange(     _Extra=extra)
	1	:	iaxis = result -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Filter called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LT 0) THEN iaxis=4
IF(iaxis GT 4) THEN BEGIN
	MESSAGE, "Couldn't successfully parse command line for axis identifer", /INFORMATIONAL
	RETURN, 0
ENDIF
result -> Restrict
;
; Get grid structure
;
axisString = STRING(iaxis, FORMAT='(i1)')
commandString = 'Grid = result.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *Grid.values
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
;
; Check that grid is uniform over specified interval
;
IF(NOT Grid.uniform) THEN BEGIN
	MESSAGE, "Grid is not uniform over specified interval.", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get time step (deltaT) and total length of signal (Ttot)
; 
deltaT = gridValues[1] - gridValues[0]
Ttot = gridValues[imax] - gridvalues[imin]
;
; set up 'temp' array to process data from 'self' (which has been copied into result...)
;
periodic = 0b
IF(grid.boundary EQ "periodic (open)") THEN BEGIN
	Ttot = Ttot + deltaT
	periodic = 1b
ENDIF
nPoints = imax - imin + 1
IF(grid.boundary EQ "periodic (closed)") THEN BEGIN
	nPoints = nPoints-1
	periodic = 1b
ENDIF
IF(periodic) THEN BEGIN
	ntemp = npoints
	GOTO, DONE1
ENDIF
FOR i=0,24 DO BEGIN
	ntemp = 2^i
	IF(nTemp GE 2*nPoints) THEN GOTO, DONE1
ENDFOR
MESSAGE, "too many elements in input array", /INFORMATIONAL
RETURN, 0

DONE1:
valuePtr = result -> GetValues()
values = *valuePtr
info = SIZE(values)
valueType = TypeOf(values)
CASE iaxis OF
	1:	BEGIN
			temp = MAKE_ARRAY(ntemp, info[2], info[3], info[4], TYPE=valueType)
			temp[ 0:(npoints-1), *, *] = values[ 0:(npoints-1), *, *, *]
		END
	2:	BEGIN
			temp = MAKE_ARRAY(info[1], ntemp, info[3], info[4], TYPE=valueType)
			temp[ *, 0:(npoints-1), *] = values[ *, 0:(npoints-1), *, *]
		END
	3:	BEGIN
			temp = MAKE_ARRAY(info[1],info[2],  ntemp, info[4], TYPE=valueType)
			temp[ *, *, 0:(npoints-1)] = values[ *, *, 0:(npoints-1), *]
		END
	4:	BEGIN
			temp = MAKE_ARRAY(info[1],info[2],  info[3], ntemp, TYPE=valueType)
			temp[ *, *, 0:(npoints-1)] = values[ *, *, *, 0:(npoints-1)]
		END
ENDCASE
;
; Get central frequency and width
;
Omega_0 = 0
dT=Ttot/10.
sgn = -1.
commandLineInfo = GetKeyWord("Omega_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN Omega_0 = commandLineInfo
commandLineInfo = GetKeyWord("k_0", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	Omega_0 = commandLineInfo
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dOmega", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = (2*!Pi/commandLineInfo) > deltaT
commandLineInfo = GetKeyWord("dk", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = (2*!Pi/commandLineInfo) > deltaT
	sgn = 1.
ENDIF
commandLineInfo = GetKeyWord("dT", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN dT = commandLineInfo > deltaT
commandLineInfo = GetKeyWord("dL", extra)
IF(TypeOf(commandLineInfo) NE 7) THEN BEGIN
	dT = commandLineInfo > deltaT
	sgn = 1.
ENDIF
;
; Create digital filter
;
filterLength = FIX(dT/deltaT+0.5) < nPoints
filterLength = filterLength > 2
filter = COMPLEXARR(nTemp)
cosineWave = EXP( COMPLEX(0., sgn)*Omega_0*deltaT*(FINDGEN(filterLength) - (filterLength/2)) )
filterEnvelope = MyHanning(filterLength)
filterNorm = 1./TOTAL(filterEnvelope)
filter[0:filterLength-1] = filterEnvelope*cosineWave
filter = filterNorm*SHIFT( filter,  -(filterLength+1)/2 )
;
; Perform a convolution between data in 'temp' and 'filter' using FFT's
;
temp = FFT(temp, 1, /OVERWRITE)
filter = FFT(filter, 1, /OVERWRITE)
CASE iaxis OF
	1:	BEGIN
			fill = REPLICATE(1., info[2])
			filter = filter#fill
		END
	2:	BEGIN
			fill = REPLICATE(1., info[1])
			filter = fill#filter
		END
	3:	BEGIN
			fill = REPLICATE(1., info[2])
			filter = fill#filter
		END
	4:	BEGIN
			fill = REPLICATE(1., info[3])
			filter = fill#filter
		END
ENDCASE
CASE iaxis OF
	1:	FOR i=0, info[3]-1 DO 		$
			FOR j=0, info[4]-1 DO		$
				temp[*,*,i,j] = temp[*,*,i,j]*filter
	2:	FOR i=0, info[3]-1 DO 		$
			FOR j=0, info[4]-1 DO		$
				temp[*,*,i,j] = temp[*,*,i,j]*filter 
	3:	FOR i=0, info[1]-1 DO 		$
			FOR j=0, info[4] DO		$
				temp[i,*,*,j] = temp[i,*,*,j]*filter
	4:	FOR i=0, info[1]-1 DO 		$
			FOR j=0, info[2] DO		$
				temp[i,j,*,*] = temp[i,j,*,*]*filter
ENDCASE
;
; Invert transform to obtain filter data
;
temp = FFT(temp, -1, /OVERWRITE)
;
; If necessary, deal with boundary conditions 
;
IF(NOT Query_Complex(temp) ) THEN temp = FLOAT(temp)
IF(Grid.boundary EQ "periodic (closed)") THEN BEGIN
	valuesOut = MAKE_ARRAY(info[1], info[2], info[3], info[4], TYPE=valueType)
	CASE iaxis OF
		1:	BEGIN
				valuesOut[0:(npoints-1), *, *, *] = temp[0:(npoints-1), *, *, *]
				valuesOut[      npoints, *, *] = temp[0, *, *]
			END 
		2:	BEGIN
				valuesOut[*, 0:(npoints-1), *, *] = temp[*, 0:(npoints-1), *, *]
				valuesOut[*, npoints      , *, *] = temp[*, 0, *, *]
			END 
		3:	BEGIN
				valuesOut[*, *, 0:(npoints-1), *] = temp[*, *, 0:(npoints-1), *]
				valuesOut[*, *, npoints      , *] = temp[*, *, 0, *]
			END 
		4:	BEGIN
				valuesOut[*, *, *, 0:(npoints-1)] = temp[*, *, *, 0:(npoints-1)]
				valuesOut[*, *, *, npoints      ] = temp[*, *, *, 0]
			END 
	ENDCASE
ENDIF ELSE BEGIN
	CASE iaxis OF
		1:	BEGIN
				valuesOut = temp[0:(npoints-1), *, *, *]
			END
		2:	BEGIN
				valuesOut = temp[*, 0:(npoints-1), *, *]
			END
		3:	BEGIN
				valuesOut = temp[*, *, 0:(npoints-1), *]
			END
		4:	BEGIN
				valuesOut = temp[*, *, *, 0:(npoints-1)]
			END
	ENDCASE
ENDELSE
;
; Prepare output object
;
PTR_FREE, result.values
result.values = PTR_NEW(valuesOut)
vmin = GKVsd_MIN(valuesOut, Max=vmax)
result.vrange = [vmin, vmax]
;
; Prepare output object
;
;
; The filtered array has a different time/space grid
; (same spacing, but possibly different length)
;
IF(NOT periodic) THEN BEGIN
	newGridValues = gridValues[0] + deltaT*FINDGEN(nTemp)
	PTR_FREE, Grid.values
	Grid.values = PTR_NEW(newGridValues)
	commandString = 'result.Grid' + axisString + ' = Grid'
	ok = EXECUTE(commandString)
	result -> signalwindow, t=[gridValues[imin], gridValues[imax]]
ENDIF
result -> Restrict
;
; and we're done
;
RETURN, result
END ; ****** GKVs4D::Filter ****** ;






