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

FUNCTION GTC_FileName, File_Path
;
; Cracks argument 'File_Path', and returns the name of the file
; 	(with extensions stripped off)
;
separator="/"							; Path separator for unix devices
IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
IF (!D.NAME EQ 'WIN') THEN separator="\"			; or windows systems
substrings=STRSPLIT(File_Path, separator, /Extract)		; Break file_name into substrings
n_strings=N_ELEMENTS(substrings)
FileName=substrings[n_strings-1]				; Strip leading directory names of filename
separator="."
subsubstrings=STRSPLIT(filename, separator, /Extract)		; Break filename at dots
Runname=subsubstrings[0]					; Run_name is leading piece of "filename"
RETURN, RunName
END ; ****** GTC_FileName ****** ;


FUNCTION GTC_PoloidalAvg, ncdFileID, ncdVarID, Rundiminfo
;
; Function returns a structure containing:
;	
;		ValuePtr		Pointer to an array containing values of data averaged
;					over poloidal angle (-� < theta � �) vs. radius and 
;					toroidal angle (a field line label, so average is a 
;					field-line average) 
;
;		ValueSqPtr		Pointer to an array containing values of the data
;					squared, and then averaged over toroidal angle 
;					vs. radius and poloidal angle.
;
;		avgThetaPtr	Pointer to an array containing values of the poloidal
;					angle (as per usual toroidal coordinates--NOT field-
;					line label) averaged over field line.
;
;		stdevThetaPtr	Pointer to an array containing values of the standard
;					deviation in the poloidal angle from previous field-
;					line average.
;
;		radialPtr		Pointer to an array containing values of the radial 
;					coordinate corresponding to the first index of the   
;					'value' array.  The radial range defaults to 
;					0.1 � r/a < 0.9.
;
;		zetaPtr		Pointer to an array containing values of the toroidal angles
;					correspoing to the second index of the 'value' array.
;					Theta is on the interval -� � zeta � �.
;
;	Written by W.M. Nevins
;	3/8/00
;
; NOTE:  This algorithm ***DOES NOT WORK PROPERLY*** when |zeta0*iota| > pi.  Since
; output is on -pi � zeta0 � pi, this only causes a problem in regions where iota >1.
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Just set RadialRange for the moment...
; 
RadialRange = [0.1, 0.9]
;
; Read ALL data from ncdFile
;
NCDF_VARGET, ncdFileID, ncdVarId, data
;
; Now, interpolate this data onto a rectangular grid in radius, poloidal angle space.
;	This rectangular grid is fit to the number of poloidal angles on the final radial
;	surface (nThetas0
;
;
; Compute number of zeta0 (toroidal angle viewed as a field-line label) grid points needed for last flux surface in radius
;
nZeta0s = FIX(iota[nRadialGrids-1]*(poloidalGrids[nRadialGrids-1]-1)) + 1
IF(2*(nZeta0s/2) EQ nZeta0s) THEN nZeta0s=nZeta0s+1		; Coherce nZetas to be an ODD number.
values     = FLTARR(nRadialGrids, nZeta0s)				; Create arrays to hold field line
ValuesSq   = FLTARR(nRadialGrids, nZeta0s)				;	average of data, average of data^2,
avgTheta   = FLTARR(nRadialGrids, nZeta0s)				;	data^2 weighted average of the poloidal angle,
stdevTheta = FLTARR(nRadialGrids, nZeta0s) 			;	and data^2 weighted standard deviation of the poloidal angle.
dZetaOne = 1./(nZeta0s-1)							
zetaOne  = dZetaOne*FINDGEN(nZeta0s)					; Toloidal angle normalized so that one cycle covers 0 � zeta � 1
 FOR i=0, nRadialGrids-1 DO BEGIN
	nPGrids = poloidalGrids[i]						; Get number of theta grids at this radius	
	theta_0 = 2.*!PI*FINDGEN(nPGrids-1)/(nPGrids-1)		; Create array of the field-line label, theta_0
	theta_0 = theta_0 - 2.*!PI*(theta_0 GT !PI)		; Use periodicity in theta_0 to wrap theta_0's back onto the interval -� < theta_0 � �
	nTGrids = FLOAT(q[i]*(nToroidalGrids-1))			; Number of toroidal grids in one full poloidal circuit
	deltaTheta = 2.*!PI*iota[i]/(nToroidalGrids-1.)		; Rotational transform across ONE toroidal grid
	thetas = FLTARR(nPGrids-1, nToroidalGrids-1)		; Make  'thetas' array which contains local value of poloidal angle at each grid point
	FOR j=0, nToroidalGrids-2 DO thetas[*,j]=theta_0+j*deltaTheta
	;				  						
 	; Average data over field line from theta= -� to �
 	;
 	temp   = FLTARR(nPGrids-1)
 	tempSq = FLTARR(nPGrids-1)
 	temp1  = FLTARR(nPGrids-1)
 	temp2  = FLTARR(nPGRids-1)
	nCycles = LONG(q[i]) + 1 ; (*** '+ 3' didn't help***)		; Number of toroidal cycles required to go a full poloidal cycle (rounded up to next full integer)
 	FOR iCycle = 0, nCycles DO BEGIN				; Sum over first NCycles toroidally
 		;
 		; Create forward and backward masks (=1 if grid point should be included in average, 0 otherwise)
 		;
 		 forwardMask = ((Thetas + 2.*!PI*iota[i]*iCycle)     GT -!PI) AND ((Thetas + 2.*!PI*iota[i]*iCycle)     LE !PI)
 		backwardMask = ((Thetas - 2.*!PI*iota[i]*(iCycle+1)) GT -!PI) AND ((Thetas - 2.*!PI*iota[i]*(iCycle+1)) LE !PI)
		;
		; Compute shift (in theta_0) of this cycle relative to fundamental cycle
		;
		ithetaForward 	= (    iCycle*IndexShift[i]) MOD (nPGrids - 1)
		ithetaBack		= ((iCycle+1)*IndexShift[i]) MOD (nPGrids - 1)
		;
		; Shift 'data' (over first--poloidal--index, theta0) to obtain 'data' values in 'forward' and 'backward' zones in zeta
		;
		forward  = SHIFT(data[start[i]:finish[i]-1, 0:nToroidalGrids-2],  -iThetaForward, 0)	; Shift 'data' back from 'iThetaForward' zone.
		backward = SHIFT(data[start[i]:finish[i]-1, 0:nToroidalGrids-2],    +iThetaBack, 0)	; Shift 'data' forward from 'iThetaBack' zone.
		;
		; Apply masks (so data is zero when field line passes out of first poloidal circuit) 
		;
		forward  =  forward*forwardMask
		backward =	backward*backwardMask	
		;
		; Sum data along field lines (zeta index), computing the contribution 
		;	of these zones to the 0th, 1st, and 2nd moments (with respect to poloidal angle)  
		;
		temp  = temp  	+ TOTAL(forward , 2) 									$
					+ TOTAL(backward, 2)
		tempSq= tempSq	+ TOTAL( forward^2, 2)								$
					+ TOTAL(backward^2, 2)
		temp1 = temp1 	+ TOTAL( forward^2*(thetas + 2.*!PI*iota[i]*iCycle),       2) 	$
					+ TOTAL(backward^2*(thetas - 2.*!PI*iota[i]*(iCycle+1)),   2)	
		temp2 = temp2 	+ TOTAL( forward^2*(thetas + 2.*!PI*iota[i]*iCycle)^2,     2) 	$
					+ TOTAL(backward^2*(Thetas - 2.*!PI*iota[i]*(iCycle+1))^2, 2)
	ENDFOR
	;
	; Compute <theta>, standard deviation = SQRT(<theta^2> - <theta>^2)
	;
	temp1 = temp1/tempSq
	temp2 = temp2/tempSq
	temp2 = SQRT(temp2 - temp1^2)	
	;
	; Divide by number of toroidal grid points that were summed over
	;
	temp   =   temp/nTGrids
	tempSq = tempSq/nTgrids
	; 
	; Interpolate from local (in radius) theta0 grid [0 through (nPgrids)-1] onto global zeta0 grid
	;	to do this, we will need 
	;
	temp   = [  temp[nPGrids-3:nPGrids-2],  temp,  temp[0:2]]	; Add some boundary points at front and back
	tempSq = [tempSq[nPGrids-3:nPGrids-2],tempSq,tempSq[0:2]]
	temp1  = [ temp1[nPGrids-3:nPGrids-2], temp1, temp1[0:2]]
	temp2  = [ temp2[nPGrids-3:nPGrids-2], temp2, temp2[0:2]]
	;
	; Set up an array, containing the theta0's corresponding to the 
	; 	(global) values of zeta0.
	;
	x=-iota[i]*ZetaOne
	break = -(0.5 - 2*iota[i]*dZetaOne)
	y=(1+x)*(x GT break) + (iota[i]+x)*(x LE break)
	x = (nPGrids-1)*y + 2		; 2 to account for 2 extra values at front) 
	
	temp   = INTERPOLATE(  temp, x, Cubic=-.5)				; and interpolate temp (tempSq, temp1, temp2) onto this grid
	tempSq = INTERPOLATE(tempSq, x, Cubic=-.5)				; This should approximate band-limited interpolation.
	temp1  = INTERPOLATE( temp1, x, Cubic=-.5)				;	(See IDL Ref. Manual,INTERPOLATE function, 
	temp2  = INTERPOLATE( temp2, x, Cubic=-.5)				;	 keyword 'Cubic')
	values[i,*]      = temp							; Load result into output array(s)
	ValuesSq[i, *]   = tempSq - temp^2
	avgTheta[i,*]    = temp1
	stdevTheta[i, *] = temp2
ENDFOR
;
; Set up radial Grid arrays
;
radius = RadialRange[0] + (RadialRange[1] - RadialRange[0])*FINDGEN(nRadialGrids)/(nRadialGrids - 1)
Zeta0 = 2.*!PI*ZetaOne
;
;
; Find elements of 'Zeta_0' array closest to �
;
smallStuff = MIN((Zeta0-!PI)^2, iZetaShift)
IF(Zeta0[iZetaShift] GT !PI) THEN iZetaShift = iZetaShift - 1 
;
; Shift values of 'Zeta0' back into zone between -� and �
;
Zeta0[iZetashift:nZeta0s-1] = Zeta0[iZetashift:nZeta0s-1] - 2.*!PI
;
; Now shift 'Zeta0', and 'values' arrays so that elements corresponding to Zeta -� come first
;
Zeta0[        0:nZeta0s-2] = SHIFT(     Zeta0[   0:nZeta0s-2],    -iZetaShift)
values[    *, 0:nZeta0s-2] = SHIFT(    values[*, 0:nZeta0s-2], 0, -iZetaShift)
valuesSq[  *, 0:nZeta0s-2] = SHIFT(  valuesSq[*, 0:nZeta0s-2], 0, -iZetaShift) 
avgTheta[  *, 0:nZeta0s-2] = SHIFT(  avgTheta[*, 0:nZeta0s-2], 0, -iZetaShift)
stdevTheta[*, 0:nZeta0s-2] = SHIFT(stdevTheta[*, 0:nZeta0s-2], 0, -iZetaShift) 
;
; Now set endpoints using periodicity
;
     Zeta0[   nZeta0s-1] =      Zeta0[   0] + 2.*!PI
     
     
    values[*, nZeta0s-1] =     values[*, 0]
  valuesSq[*, nZeta0s-1] =   valuesSq[*, 0]
  avgTheta[*, nZeta0s-1] =   avgTheta[*, 0]
stdevTheta[*, nZeta0s-1] = stdevTheta[*, 0]
;
; Form output structure, and return result
;
result = {	valuePtr:PTR_NEW(values), 		valueSqPtr:PTR_NEW(valuesSq),		$
		avgThetaPtr:PTR_NEW(avgTheta), 	stdevThetaPtr:PTR_NEW(stdevTheta), 		$
		radialPtr:PTR_NEW(radius), 	zetaPtr:PTR_NEW(Zeta0) }
RETURN, result

END ; ****** GTC_PoloidalAvg ****** ;


Function GTC_ToroidalAvg, ncdFileID, ncdVarID, Rundiminfo
;
; Function returns a structure containing:
;	
;		ValuePtr		Pointer to an array containing values of data averaged
;					over toroidal angle (-� < zeta � �) vs. radius and 
;					poloidal angle (a field line label, so average is a 
;					field-line average) 
;
;		ValueSqPtr		Pointer to an array containing values of the data
;					squared, and then averaged over toroidal angle 
;					vs. radius and poloidal angle.
;
;		avgThetaPtr	Pointer to an array containing values of the poloidal
;					angle (as per usual toroidal coordinates--NOT field-
;					line label) averaged over field line.
;
;		stdevThetaPtr	Pointer to an array containing values of the standard
;					deviation in the poloidal angle from previous field-
;					line average.
;
;		radialPtr		Pointer to an array containing values of the radial 
;					coordinate corresponding to the first index of the   
;					'value' array.  The radial range defaults to 
;					0.1 � r/a < 0.9.
;
;		ThetaPtr		Pointer to an array containing values of the poloidal angles
;					correspoing to the second index of the 'value' array.
;					Theta is on the interval -� � theta � �.
;
;	Written by W.M. Nevins
;	3/8/00
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Just set RadialRange for the moment...
; 
RadialRange = [0.1, 0.9]
;
; Read ALL data from ncdFile
;
NCDF_VARGET, ncdFileID, ncdVarId, data
;
; Now, interpolate this data onto a rectangular grid in radius, poloidal angle space.
;	This rectangular grid is fit to the number of poloidal angles on the final radial
;	surface (nThetas0
;
nThetas    = poloidalGrids[nRadialGrids-1]
IF(2*(nThetas/2) EQ nThetas) THEN nThetas=nThetas+1		; Coherce nThetas to be an ODD number.
values     = FLTARR(nRadialGrids, nThetas)
ValuesSq   = FLTARR(nRadialGrids, nThetas)
avgTheta   = FLTARR(nRadialGrids, nThetas)
stdevTheta = FLTARR(nRadialGrids, nThetas) 
dThetaOne = 1./(nThetas-1)							
thetaOne  = dThetaOne*FINDGEN(nThetas)				; Poloidal angle normalized so that one cycle covers 0 � theta � 1
 FOR i=0, nRadialGrids-1 DO BEGIN
	nPGrids = poloidalGrids[i]				; Get number of theta grids at this radius	
	theta_0 = 2.*!PI*FINDGEN(nPGrids-1)/(nPGrids-1)		; Create array of the field-line label, theta_0
	theta_0 = theta_0 - 2.*!PI*(theta_0 GT !PI)		; Use periodicity in theta_0 to wrap theta_0's back onto the interval -� < theta_0 � �
	nTGrids = FLOAT(q[i]*(nToroidalGrids-1))		; Number of toroidal grids in one full poloidal circuit
	deltaTheta = 2.*!PI*iota[i]/(nToroidalGrids-1.)		; Rotational transform across ONE toroidal grid
	thetas = FLTARR(nPGrids-1, nToroidalGrids-1)		; Make  'thetas' array which contains local value of poloidal angle at each grid point
	FOR j=0, nToroidalGrids-2 DO thetas[*,j]=theta_0+j*deltaTheta
	;										
 	; Average data over field line from zeta= -� to �
 	;
 	temp   = FLTARR(nPGrids-1)
 	tempSq = FLTARR(nPGrids-1)
 	temp1  = FLTARR(nPGrids-1)
 	temp2  = FLTARR(nPGRids-1)
	nCycles = LONG(q[i]) + 1					; Number of toroidal cycles required to go a full poloidal cycle (rounded up to next full integer)
 	FOR iCycle = 0, nCycles DO BEGIN				; Sum over first NCycles toroidally
 		;
 		; Create forward and backward masks (=1 if grid point should be included in average, 0 otherwise)
 		;
 		 forwardMask = ((Thetas + 2.*!PI*iota[i]*iCycle)     GT -!PI) AND ((Thetas + 2.*!PI*iota[i]*iCycle)     LE !PI)
 		backwardMask = ((Thetas - 2.*!PI*iota[i]*(iCycle+1)) GT -!PI) AND ((Thetas - 2.*!PI*iota[i]*(iCycle+1)) LE !PI)
		;
		; Compute shift (in theta_0) of this cycle relative to fundamental cycle
		;
		ithetaForward 	= (    iCycle*IndexShift[i]) MOD (nPGrids - 1)
		ithetaBack		= ((iCycle+1)*IndexShift[i]) MOD (nPGrids - 1)
		;
		; Shift 'data' (over first--poloidal--index, theta0) to obtain 'data' values in 'forward' and 'backward' zones in zeta
		;
		forward  = SHIFT(data[start[i]:finish[i]-1, 0:nToroidalGrids-2], -iThetaForward, 0)	; Shift 'data' back from 'iThetaForward' zone.
		backward = SHIFT(data[start[i]:finish[i]-1, 0:nToroidalGrids-2],     iThetaBack, 0)	; Shift 'data' forward from 'iThetaBack' zone.
		;
		; Apply masks (so data is zero when field line passes out of first poloidal circuit) 
		;
		forward  =  forward*forwardMask
		backward =	backward*backwardMask	
		;
		; Sum data along field lines (zeta index), computing the contribution 
		;	of these zones to the 0th, 1st, and 2nd moments (with respect to poloidal angle)  
		;
		temp  = temp  	+ TOTAL(forward , 2) 									$
					+ TOTAL(backward, 2)
		tempSq= tempSq	+ TOTAL( forward^2, 2)								$
					+ TOTAL(backward^2, 2)
		temp1 = temp1 	+ TOTAL( forward^2*(thetas + 2.*!PI*iota[i]*iCycle),       2) 	$
					+ TOTAL(backward^2*(thetas - 2.*!PI*iota[i]*(iCycle+1)),   2)	
		temp2 = temp2 	+ TOTAL( forward^2*(thetas + 2.*!PI*iota[i]*iCycle)^2,     2) 	$
					+ TOTAL(backward^2*(Thetas - 2.*!PI*iota[i]*(iCycle+1))^2, 2)
	ENDFOR
	;
	; Compute <theta>, standard deviation = SQRT(<theta^2> - <theta>^2)
	;
	temp1 = temp1/tempSq
	temp2 = temp2/tempSq
	temp2 = SQRT(temp2 - temp1^2)	
	;
	; Divide by number of toroidal grid points that were summed over
	;
	temp   =   temp/nTGrids
	tempSq = tempSq/nTgrids
	; 
	; Interpolate from local (in radius) theta0 grid [0 through (nPgrids)-1] onto global Theta0 grid
	;
	temp   = [  temp[nPGrids-3:nPGrids-2],  temp,  temp[0:2]]	; Add some boundary points at front and back
	tempSq = [tempSq[nPGrids-3:nPGrids-2],tempSq,tempSq[0:2]]
	temp1  = [ temp1[nPGrids-3:nPGrids-2], temp1, temp1[0:2]]
	temp2  = [ temp2[nPGrids-3:nPGrids-2], temp2, temp2[0:2]]
	x = (nPGrids-1)*thetaOne + 2						; Set up array of theta values (+2 to account for 2 extra values at front) 
	temp   = INTERPOLATE(  temp, x, Cubic=-.5)				; and interpolate temp (tempSq, temp1, temp2) onto this grid
	tempSq = INTERPOLATE(tempSq, x, Cubic=-.5)				; This should approximate band-limited interpolation.
	temp1  = INTERPOLATE( temp1, x, Cubic=-.5)				;	(See IDL Ref. Manual,INTERPOLATE function, 
	temp2  = INTERPOLATE( temp2, x, Cubic=-.5)				;	 keyword 'Cubic')
	values[i,*]      = temp							; Load result into output array(s)
	ValuesSq[i, *]   = tempSq - temp^2
	avgTheta[i,*]    = temp1
	stdevTheta[i, *] = temp2
ENDFOR
;
; Set up radial Grid arrays
;
radius = RadialRange[0] + (RadialRange[1] - RadialRange[0])*FINDGEN(nRadialGrids)/(nRadialGrids - 1)
theta = 2.*!PI*ThetaOne
;
; Find elements of 'theta_0' array closest to �
;
smallStuff = MIN((theta-!PI)^2, iThetaShift)
IF(theta[iThetaShift] GT !PI) THEN iThetaShift = iThetaShift - 1 
;
; Shift values of 'theta' back into zone between -� and �
;
theta[iThetashift:nThetas-1] = theta[iThetashift:nThetas-1] - 2.*!PI
;
; Now shift 'theta', and 'values' arrays so that elements corresponding to theta -� come first
;
theta[        0:nThetas-2] = SHIFT(     theta[   0:nThetas-2],    -iThetaShift)
values[    *, 0:nThetas-2] = SHIFT(    values[*, 0:nThetas-2], 0, -iThetaShift)
valuesSq[  *, 0:nThetas-2] = SHIFT(  valuesSq[*, 0:nThetas-2], 0, -iThetaShift) 
avgTheta[  *, 0:nThetas-2] = SHIFT(  avgTheta[*, 0:nThetas-2], 0, -iThetaShift)
stdevTheta[*, 0:nThetas-2] = SHIFT(stdevTheta[*, 0:nThetas-2], 0, -iThetaShift) 
;
; Now set endpoints using periodicity
;
     theta[   nThetas-1] =      theta[   0] + 2.*!PI
     
     
    values[*, nThetas-1] =     values[*, 0]
  valuesSq[*, nThetas-1] =   valuesSq[*, 0]
  avgTheta[*, nThetas-1] =   avgTheta[*, 0]
stdevTheta[*, nthetas-1] = stdevTheta[*, 0]
;
; Form output structure, and return result
;
result = {	valuePtr:PTR_NEW(values), 		valueSqPtr:PTR_NEW(valuesSq),			$
		avgThetaPtr:PTR_NEW(avgTheta), 	stdevThetaPtr:PTR_NEW(stdevTheta), 		$
		radialPtr:PTR_NEW(radius), 	thetaPtr:PTR_NEW(theta) }
RETURN, result

END  ; ****** GTC_ToroidalAvg ****** ;


FUNCTION GTC_PoloidalCut, ncdFileID, ncdVarID, theta, Rundiminfo
;
; Function returns a structure containing:
;	
;		ValuePtr		Pointer to an array containing values of data at 
;					fixed poloidal angle (-� < theta � �) vs. radius and 
;					toroidal angle (in the usual sense--that is the toroidal
;					angle of the returned data is NOT to be intepreted as either
;					a field-line label nor as distance along a field line). 
;
;		radialPtr		Pointer to an array containing values of the radial 
;					coordinate corresponding to the first index of the   
;					'value' array.  The radial range defaults to 
;					0.1 � r/a < 0.9.
;
;		zetaPtr			Pointer to an array containing values of the toroidal angles
;					correspoing to the second index of the 'value' array.
;					Theta is on the interval -� � zeta � �.
;
;	Written by W.M. Nevins
;	12/21/00
;
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Just set RadialRange for the moment...
; 
RadialRange = [0.1, 0.9]
;
; Read ALL data from ncdFile
;
NCDF_VARGET, ncdFileID, ncdVarId, data
;
;	Now, slice this data at fixed (laboratory) poloidal angle, and store result
;	on a regular rectangular grid in (radius, toroidal angle)-space.
;
oldDzeta = 2.*!PI/(nToroidalGrids - 1.)
oldZetas = oldDzeta*FINDGEN(nToroidalGrids)				; Create array of toroidal angles at input resolution
nZetas = LONG( MAX(iota*poloidalGrids) ) 				; Compute required resolution in zeta
newDzeta = 2.*!PI/(nZetas-1.)
newZetas = newDzeta*FINDGEN(nZetas)					; Create array of toriodal angles at this resolution
iNewZetas = LONG(newZetas/oldDzeta)
dNewZetas = newZetas/oldDzeta - iNewZetas

values = FLTARR(nRadialGrids, nZetas)						; Create array to hold "sliced" data
;				  						
; Interpolate data from field line coordinates to fixed theta value
;
 
FOR i=0, nRadialGrids-1 DO BEGIN
	nPGrids = poloidalGrids[i]						; Get number of theta grids at this radius
	dTheta = 2.*!PI/(nPGrids-1)						; Compute spacing in poloidal angle at this radius
	oldData = FLTARR(nPGrids, nToroidalGrids+1)				; Create array to hold netCDF file data at this radius (including one boundary cell in zeta)
	oldData[*,0:(nToroidalGrids-1)] = data[start[i]:finish[i], *]		;	and load array with data from netCDF file.
	;
	; Try using toroidal + shift periodicity to "patch" the problem with input data values
	; 	at toroidal index=0
	;
	indices = (INDGEN(poloidalGrids[i]) - Indexshift[i]) MOD (poloidalGrids[i] - 1)
	indices = indices + (indices LT 0)*(poloidalGrids[i] - 1)
	oldData[*, 0] = oldData[indices, NtoroidalGrids-1]
	;
	; interpolate data onto higher resolution (in zeta) magnetic coordinate grid
	;
	newData = FLTARR(nPGrids+1, nZetas)					; Create array to hold (higher-resolution in zeta) data in magentic coordinates at this radius
	FOR j=0, nPGrids-1 DO	newData[j,*] = (1-dNewZetas)*oldData[ j, iNewZetas] + dNewZetas*oldData[ j, iNewZetas + 1]	
	;
	; Now, interpolate values from high-resolution magnetic coordinate grid to fixed poloidal angle
	;
	deltaTheta = newDzeta*iota[i]						; Compute rotational transform across ONE toroidal grid	(at the higher resolution) at this radius 
 	magTheta = ( theta - deltaTheta*FINDGEN(nZetas) )/dTheta		; Compute array of values of the desired poloidal angles vs. toroidal index 
	magTheta = magTheta MOD (nPGrids-1)					;	(in units of the current poloidal grid spacing) in magentic coordinates
	magTheta = magTheta + (nPGrids-1)*(magTheta LT 0.)			;	and use periodicity to map magTheta into the interval 0 � magTheta < nPGrids.
 	iTheta = LONG(magTheta)							; Compute poloidal grid index (iTheta) at each
 	delta = magTheta-iTheta							; 	toroidal angle, and (normalized to grid spacing) displacement from that poloidal grid point (delta). 
 	star = FINDGEN(nZetas)
	values[i,star] = (1-delta)*newData[iTheta, star] + delta*newData[(iTheta+1), star]
 ENDFOR
;
; Set up arrays of radial and toroidal Grid values
;
radius = RadialRange[0] + (RadialRange[1] - RadialRange[0])*FINDGEN(nRadialGrids)/(nRadialGrids - 1)
;
;
; Find elements of 'Zeta' array closest to �
;
;smallStuff = MIN((Zeta-!PI)^2, iZetaShift)
;IF(Zeta[iZetaShift] LE !PI) THEN iZetaShift = iZetaShift + 1 
;
; Shift values of 'Zeta' back into zone between -� and �
;
;Zeta[iZetashift:(nToroidalGrids-1)] = Zeta[iZetashift:nToroidalGrids-1] - 2.*!PI
;
; Now shift 'Zeta', and 'values' arrays so that elements corresponding to Zeta -� come first
;
;Zeta   = SHIFT(  Zeta,    -iZetaShift)
;values = SHIFT(values, 0, -iZetaShift)
;
; Form output structure, and return result
;
result = { valuePtr:PTR_NEW(values), radialPtr:PTR_NEW(radius), zetaPtr:PTR_NEW(newZetas) }
RETURN, result

END ; ****** GTC_PoloidalCut ****** ;



Function GTC_ToroidalCut, ncdFIleID, ncdVarID, zeta, RundimInfo
;
; Function returns a structure containing:
;
;		ValuePtr	Pointer to an array containing values of data at fixed
;				toloidal angle vs. radius and poloidal angle.    
;
;		radialPtr	Pointer to an array containing values of the radial 
;				coordinate corresponding to the first index of the   
;				'value' array.  The radial range defaults to 
;				 0.1 � r/a < 0.9.
;
;		ThetaPtr	Pointer to an array containing values of the poloidal angles
;				correspoing to the second index of the 'value' array.
;				Theta is on the interval -� � zeta � �.
;
;	Written by W.M. Nevins
;	3/7/00
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Just set RadialRange for the moment...
; 
RadialRange = [0.1, 0.9]
;
; Find toroidal index required, and offset of requested 
;	value of zeta from the zeta grid
;
dZeta = 2.*!PI/(nToroidalGrids - 1)
iZeta =  LONG(zeta/dZeta)
deltaZeta = zeta - iZeta*dZeta
IF(nToroidalGrids EQ 1) THEN BEGIN 
	deltaZeta = 0.
	iZeta=0
ENDIF
;
; Compute number of data points to be read
;
npoints = LONG(TOTAL(poloidalGrids))
;
; Read data from ncdFile
;
IF(nToroidalGrids EQ 1) THEN BEGIN
	NCDF_VARGET, ncdFileID, ncdVarId, data, Count=[npoints, 1], Offset=[0, iZeta]
ENDIF ELSE BEGIN
	NCDF_VARGET, ncdFileID, ncdVarId, data, Count=[npoints, 2], Offset=[0, iZeta]
ENDELSE
;
; Now, interpolate this data onto a rectangular grid in radius, poloidal angle space.
;	This rectangular grid is fit to the number of poloidal angles on the final radial
;	surface (nThetas0
;
nThetas = poloidalGrids[nRadialGrids-1]
values = FLTARR(nRadialGrids, nThetas)
dTheta = 1./(nThetas-1)							; Normalize poloidal angle to a cycle of 0 � theta � 1
theta = dTheta*FINDGEN(nThetas)						;	over one full cycle for now
 FOR i=0, nRadialGrids-1 DO BEGIN					
 	IF(deltaZeta EQ 0.) THEN BEGIN
 		temp = data[start[i]:finish[i], 0]
 	ENDIF ELSE BEGIN					
		temp = (1.	- deltaZeta)*data[start[i]:finish[i], 0] $
			 + deltaZeta*data[start[i]:finish[i], 1]	; Linear interpolation to desired zeta value
	ENDELSE
	nPGrids = poloidalGrids[i]
	dTheta_i = 2.*!PI/(nPGrids-1)				; Spacing in theta_0 of gridpoints on ith flux surface
	deltaTheta = iota[i]*zeta				; Now, we must 'rotate' in theta_0 through the angle
	iDeltaTheta = LONG(deltaTheta/dTheta_i)			;	over which field lines have spiraled in the interval
	dDeltaTheta = deltaTheta - iDeltaTheta*dTheta_i		;	0 to zeta (remember that 'theta+0' is a field-line-label!).
	temp1 = SHIFT(temp[0:nPGrids-2], iDeltaTheta)		; Shift in Theta
	temp1 = [temp1,temp1[0]]				;	inforce periodicity
	temp2 = Shift(temp[0:nPGrids-2], (iDeltaTheta+1))	; Same thing, but shifted one more
	temp2 = [temp2,temp2[0]]				;	... still inforcing periodicity
	temp = (1. - dDeltaTheta)*temp1 + dDeltaTheta*temp2	; Interpolate to correct value of back shift in theta_0
	temp1 = [temp[nPGrids-3:nPGRids-2],temp,temp[1:2]]	; Add some boundary points at front and back
	x = (nPGrids-1)*theta + 2						; Set up array of theta values and interpolate temp onto this grid
	temp2 = INTERPOLATE(temp1, x, Cubic=-.5)			; This should approximate band-limited interpolation--see IDL Ref. Manual
	values[i,*] = temp2
ENDFOR
;
; Form grid arrays
;
radius = RadialRange[0] + (RadialRange[1]-RadialRange[0])*FINDGEN(nRadialGrids)/(nRadialGrids-1)
theta = 2.*!PI*theta
;
; Find elements of 'theta' array closest to �
;
smallStuff = MIN((theta-!PI)^2, iThetaShift)
IF(theta[iThetaShift] GT !PI) THEN iThetaShift = iThetaShift - 1
;
; Roll 'theta' back into zone between -� and �
;
theta[iThetashift:nThetas-1] = theta[iThetashift:nThetas-1] - 2.*!PI
;
; Now shift 'theta', and 'values' arrays so that elements corresponding to theta -� come first
;
 theta[   0:nThetas-2] = SHIFT( theta[   0:nThetas-2],    -iThetaShift)
values[*, 0:nThetas-2] = SHIFT(values[*, 0:nThetas-2], 0, -iThetaShift) 
;
; Use periodicity to set final element(s) of theta, values
;
theta[    nThetas-1] =  theta[   0] + 2.*!PI 
values[*, nthetas-1] = values[*, 0]
;
; form structure to return
;
result = { valuePtr:PTR_NEW(values), radialPtr:PTR_NEW(radius), thetaPtr:PTR_NEW(theta) } 
;
; and we're done
;
RETURN, result
END	; ****** GTC_ToroidalCut ****** ;


FUNCTION GTC_Surface_FieldLine, ncdFIleID, ncdVarID, iflux, RundimInfo
;
; Function returns a structure containing:
;
;		ValuePtr	Pointer to an array containing values of 'phi' on the iflux 
;				flux-surface vs. poloidal and toroidal angle in field line 
;				coordinates. That is, the poloidal angle, theta, measures 
;				distance along the field line, and toroidal angle, zeta, 
;				(the location in angle at which the field line crosses the 
;				midplane) is a field line label.   
;
;		ThetaPtr	Pointer to an array containing values of the poloidal angles
;				corresponding to the first index of the 'value' array.  
;				Theta is on the interval -� � theta < �
;
;		ZetaPtr	Pointer to an array containing values of the toroidal angles
;				correspoing to the second index of the 'value' array.
;				Zeta is on the interval -� � zeta � �.
;
;	Written by W.M. Nevins
;	3/5/00
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Read Potential (or whatever) from ncdFile
;
NCDF_VARGET, ncdFileID, ncdVarId, phi, Count=[poloidalGrids[iflux], NtoroidalGrids], Offset=[start[iflux], 0]
;
; As read from ncdFIle, array is in field-line co-ordinates with TOROIDAL angle labeling position
; along the field line.  We must extend the range in by a factor of q[iflux] (the safety factor)  
; so that the extended range covers poloidal angles from -pi to pi.
;
nThetaGrids = FIX(q[iflux]*NtoroidalGrids)
IF(2*(nThetaGrids/2) EQ nThetaGrids) THEN nThetaGrids = nThetaGrids + 1	; Coherce nThetaGrids to an odd number
nZetaGrids  = FIX(iota[iflux]*(poloidalGrids[iflux]-1)) + 1
;
; Try using toroidal + shift periodicity to "patch" the problem with phi values
; 	at toroidal index=0
;
indices = (INDGEN(poloidalGrids[iflux]) - Indexshift[iflux]) MOD (poloidalGrids[iflux] - 1)
indices = indices + (indices LT 0)*(poloidalGrids[iflux]-1)
phi[*, 0] = phi[indices, NtoroidalGrids-1]
;
; Create arrays for output values, grid values
;
values = FLTARR(nThetaGrids, NzetaGrids)
zeta   = FLTARR(nZetaGrids)
theta  = FINDGEN(nThetaGrids)
dToroidalAngle = 2.*!PI/(NtoroidalGrids - 1)			; Toroidal angle increment in ncdFile
dTheta = iota[iflux]*dToroidalAngle				; Corresponding increment in poloidal angle
iThetaZero = (nThetaGrids+1)/2 - 1				; this is the INDEX corresponding to the theta=0 element 
theta = dTheta*(theta - iThetaZero)				;	of the 'values' and 'theta' arrays
;
; Transform from coordinates in which poloidal angle is the field line label
; to coordinates in which the toroidal angle is the field line label as we 
; load Potential (or whatever) into 'values' array
;
iThetaMax = NthetaGrids - 1						; Highest valid theta index to 'values' array
thetaMax = theta[iThetaMax]						; Maximum value of theta in 'values' array
thetaMin = theta[0]
iupper = NtoroidalGrids - 1						; Highest valid zeta  index input from ncdFile
ilast = iupper - 1								; 	(but don't generally read first data point as it seems to be corrupted ...)
dPoloidalAngle = 2.*!PI/(poloidalGrids[iflux] - 1)		; Poloidal angle increment in ncdFile
max1 = NzetaGrids ; < poloidalGrids[iflux]
;
; Begin loop over field line labels, zeta
;
FOR i=0, max1-1 DO BEGIN
	ii = (max1-1) - i							; Need index to load 'zeta', 'values' grids from back to front
	zeta[ii] = - i*dPoloidalAngle*q[iflux]			; Toroidal angle label for the ith field line
	IF(zeta[ii] LT -!pi) THEN zeta[ii] = zeta[ii]+2.*!PI	; 	(shift into -� � zeta � � zone)
	thetaNCD0 = i*dPoloidalAngle + dTheta				; Value of theta at first zeta grid point from NCD input, fundamental zone
	delta = thetaNCD0 MOD dTheta					; Displacement of Value grid relative to temp grid 
	minZone = 0								; Find 'minZone' -- the number of zeta-zones from the
	thetaNCDmin = thetaNCD0						;	NCDF input data which we must shift through to 
	WHILE (thetaNCDmin GT thetaMin) DO BEGIN			;	cover Theta values down to -� at this value of zeta
		minZone = minZone - 1
		thetaNCDmin = thetaNCD0 + minZone*2.*!PI*iota[iflux] - dTheta
	ENDWHILE
	maxZone = minZone								; Find 'maxZone- -- the number of zeta-zones from the
	thetaNCDmax = thetaNCD0 + (maxZone+1)*2.*!PI*iota[iflux] 	;	NCDF input data which we must shift through to
	WHILE (thetaNCDmax LE thetaMax) DO BEGIN				;	cover theta values up to � at this zeta value
		maxZone = maxZone + 1
		thetaNCDmax = thetaNCD0 + (maxZone+1)*2.*!PI*iota[iflux]
	ENDWHILE
;
; Get next further back value of potential (in case it is needed for interpolation, and also for parallelism ...
;

	
	iback = (i+(minzone-1)*indexShift[iflux]) MOD (poloidalGrids[iflux] - 1)	; Set theta index to access negative values of zeta in input
	IF(iback LT 0) THEN iback = iback + poloidalGrids[iflux] - 1			; shift to first positive zone if necessary (generally is ...)
	endValue = phi[iback, iupper]
	
	iThetaMin = 0L								; Theta index of next element of 'values' array to be written
	
	; 
	; Begin loop over zeta-zones in NCD input data
	;	
	FOR zoneShift = minZone, maxZone DO BEGIN
		iback = (i+zoneShift*indexShift[iflux]) MOD (poloidalGrids[iflux]-1)	; Set theta index to access current zeta-zone
		IF(iback LT 0) THEN iback = iback + poloidalGrids[iflux] - 1		; shift to first positive modulus if necessary (generally is ...)
		temp = phi[iback, 1:iupper]								; get potential from current zeta zone
		temp = REFORM(temp, iupper)
		temp = [endValue, temp]									; Append last value from previous zeta-cycle at begining ...
		thetaTemp = thetaNCD0 + 2.*!PI*zoneShift*iota[iflux] + dTheta*(FINDGEN(iupper+1) -1) 	; Form array of theta values for present 'temp' array  
		FOR imin = 1L, iupper DO IF(thetaTemp[imin] GT theta[iThetaMin]) THEN GOTO, GOTMIN
			MESSAGE, "we shouldn't be here!", /informational
			imin = iupper
GOTMIN:	imin = imin - 1
		imax = iupper					
		IF(thetaTemp[ilast] GE thetaMax) THEN BEGIN
			FOR imax = imin, ilast DO IF(thetaTemp[imax] GE thetaMax) THEN GOTO, GOTMAX
			MESSSAGE, "we shouldn't be here either!", /informational
			imax = ilast
		ENDIF
GOTMAX:	imax = imax-1
		iThetaMax = iThetaMin + (imax - imin)
		; ********** IF statement is perhaps a cluge??? was patch to fix problem when
		; thetaTemp[1] = thetaMax (which resulted in imax=0 and imin=1...)
		IF(imax GE imin) THEN values[ithetaMin:iThetaMax, ii] = (1-delta)*temp[imin:imax] + delta*temp(imin+1:imax+1)
		iThetaMin = iThetaMax + 1
		endValue = temp[imax+1]
	ENDFOR
	
ENDFOR
;
; Shift 'zeta' and 'values' arrays in zeta subscript to put -� (instead of 0) as first element
;	(thereby insuring that 'zeta' increases monotonically)
;
zetamax = MAX(zeta, ishift)					; Use IDL's MAX function to get index of last positive value of zeta
ishift = ishift + 1							; 	increment by 1 so that the next value (-�) will be first element
zeta = SHIFT(zeta[0:nZetaGrids-2], -ishift)		;	of zeta array after SHIFT
values = SHIFT(values[*,0:nZetaGrids-2], 0, -ishift)	; Shift 'values' array by same amount
;
; Prepare output structure
;
result={valuePtr: PTR_NEW(values[*,0:nZetaGrids-2]), thetaPtr: PTR_NEW(theta), zetaPtr: PTR_NEW(zeta[0:nZetaGrids-2])}
;
; 
RETURN, result
END ; ****** GTC_Surface_FieldLine ****** ;



FUNCTION GTC_Surface_Normal, ncdFIleID, ncdVarID, iflux, RundimInfo
;
; Function returns a structure containing:
;
;		ValuePtr	Pointer to an array containing values of 'phi' vs. poloidal 
; 				and toroidal angle (that is, poloidal and toroidal angle as 
;				per a standard toroidal coordinate system, NOT a field-line-
;				following coordinate system) on the iflux flux-surface.
;
;		ThetaPtr	Pointer to an array containing values of the poloidal angles
;				corresponding to the first index of the 'value' array.  
;				Theta is on the interval -� � theta < �
;
;		ZetaPtr	Pointer to an array containing values of the toroidal angles
;				correspoing to the second index of the 'value' array.
;				Zeta is on the interval 0 � zeta � 2�.
;
;	Written by W.M. Nevins
;	3/4/00
;
; First extract information from Rundiminfo structure
;
poloidalGrids	= *RundimInfo.poloidalGrids
indexShift		= *RundimInfo.indexshift
nToroidalGrids	=  RundimInfo.nToroidalGrids
nRadialGrids	=  RundimInfo.nRadialGrids
start			= *RundimInfo.start
finish		= *RundimInfo.finish
iota			= *RundimInfo.iota
q			= *RundimInfo.q
;
; Read Potential (or whatever) from ncdFile
;
NCDF_VARGET, ncdFileID, ncdVarId, data, Count=[poloidalGrids[iflux], NtoroidalGrids], Offset=[start[iflux], 0]

nPGrids = poloidalGrids[iflux]						; Get number of theta grids at this radius
dTheta = 2.*!PI/(nPGrids-1)						; Compute spacing in poloidal angle at this radius
oldData = FLTARR(nPGrids, nToroidalGrids+1)				; Create array to hold netCDF file data at this radius (including one boundary cell in zeta)
oldData[*,0:(nToroidalGrids-1)] = data[*, *]				;	and load array with data from netCDF file.
;
; Try using toroidal + shift periodicity to "patch" the problem with phi values
; 	at toroidal index=0
;
indices = (INDGEN(poloidalGrids[iflux]) - Indexshift[iflux]) MOD (poloidalGrids[iflux] - 1)
indices = indices + (indices LT 0)*(poloidalGrids[iflux] - 1)
oldData[*, 0] = oldData[indices, NtoroidalGrids-1]
;
; Compute number of toroidal grid points required to achieve 
; the same resolution on 'normal' grid as was achieved on the
; field line following grid of GTC
;
oldDzeta = 2.*!PI/(nToroidalGrids - 1.)
oldZetas = oldDzeta*FINDGEN(nToroidalGrids)				; Create array of toroidal angles at input resolution
nZetas = LONG( iota[iflux]*poloidalGrids[iflux] ) 			; Compute required number of toroidal grid points
newDzeta = 2.*!PI/(nZetas-1.)						;	and associated grid spacing.
newZetas = newDzeta*FINDGEN(nZetas)					; Create array of toriodal angles at this resolution
iNewZetas = LONG(newZetas/oldDzeta)					; Compute array containing indices to old toroidal grid
dNewZetas = newZetas/oldDzeta - iNewZetas				; and (upward) displacement from this old grid point
;
; interpolate data onto higher resolution (in zeta) magnetic coordinate grid
;
newData = FLTARR(nPGrids, nZetas)					; Create array to hold (higher-resolution in zeta) data in magentic coordinates at this radius
FOR j=0, nPGrids-1 DO	$
	newData[j,*] = (1-dNewZetas)*oldData[ j, iNewZetas] + dNewZetas*oldData[ j, iNewZetas + 1]	


;
; Create array for output values
;
values = FLTARR(poloidalGrids[iflux], nZetas)
;
; Transform from field-line-following coordinates as we Load Potential 
; (or whatever) into 'output' array
;
thetaShift = 0.
dTheta = iota[iflux]*FLOAT(nPGrids-1)/FLOAT(nZetas-1)
FOR i=0, nZetas-1 DO BEGIN
	ishift = FIX(thetaShift)
	delta = thetaShift - ishift
	slice = REFORM(newData[*,i] ,nPGrids)
	values[*,i] = (1-delta)*SHIFT(slice, ishift) + delta*SHIFT(slice, ishift+1)
	thetaShift=thetaShift + dTheta
ENDFOR
;
; Set up coordinate arrays for theta (the poloidal angle) and zeta (the toroidal angle)
;
deltaTheta = 2*!PI/(nPGrids)
theta = deltaTheta*FINDGEN(poloidalGrids[iflux])
zeta  = newDzeta*FINDGEN(nZetas)
;
; Find elements of 'theta' and 'zeta' arrays closest to �
;
 smallStuff = MIN((theta-!PI)^2, iThetaShift)
osmallStuff = MIN(( zeta-!PI)^2,  iZetaShift)
IF(theta[iThetaShift] GT !PI) THEN iThetaShift = iThetaShift - 1
IF( zeta[ izetaShift] GT !PI) THEN  iZetaShift =  iZetaShift - 1
iThetaShift = iThetaShift + 1
 iZetaShift =  iZetaShift + 1
;
; Roll 'theta' (but not  'zeta') back into zone between -� and �
;
theta[iThetashift:poloidalGrids[iflux]-1] = theta[iThetashift:poloidalGrids[iflux]-1] - 2.*!PI
;zeta[ iZetashift:NtoroidalGRids-1]       =  zeta[ iZetashift:NtoroidalGRids-1]       - 2.*!PI ; This line is for the (disabled) attempt to roll zeta back into zone between -� and �
;
; Now shift 'theta', and 'values' arrays so that elements corresponding to theta -� come first
;	DON'T shift zeta, as this puts problem at zeta = 0 RIGHT in the center of the plots
;
theta = SHIFT( theta, -iThetaShift)
;zeta = SHIFT(  zeta,  -iZetaShift) ;
values= SHIFT(values, -iThetaShift, 0) ; -iZetaShift
;
; prepare structure to be returned
;
 result = {valuePtr:PTR_NEW(values), thetaPtr:PTR_NEW(theta), zetaPtr:PTR_NEW(zeta)} 

RETURN, result
END ; ****** GTC_Surface_Normal ****** ;


FUNCTION GTC_Data,	IsurfaceNormal=IsurfaceNormal, IsurfaceFieldLine=IsurfaceFieldLine,		$
			IradialNormal= IradialNormal,  PoloidalAvg=PoloidalAvg,				$
			IfluxLimits=IfluxLimits, TimeSteps=TImeSteps, Info=Info,			$
			Qprofile=Qprofile, IotaProfile=IotaProfile, Q_0=q0, Iota_0=Iota0,		$
			Variable=variableName, Dt=delta_t, R_range=rrange, zetaCut=zetaCut,		$
			thetacut=thetaCut, ToroidalAvg=ToroidalAvg, ncd=fileSuffix, Debug=Debug, path=path
;
; Function reads data from the RUNdimenID, and Pxxxx output files of GTC
; and returns a structure containing requested data within GKVsd objects
; 
; This version will only allow ONE type of data from a given call to GTC_Data.
;
;	Keywords:
;
;		Info				If this keyword is set (e.g, /Info), GTC_Data prints 
;						information about the grid structure and the RETURNS.  
;						Selected information is PRINTed to the IDL output 
;						log, and more information within the structure returned 
;						on this call to GTC_Data. NO GKVsd data objects are 
;						returned.
;
;		IsurfaceNormal		Set to an integer array containing indices to the
;						flux surfaces on which data is desired vs. usual
;						poloidal (theta) and toroidal (zeta) angle.
;
;		IsurfaceFieldLine	Set to an integer array containing indices to the
;						flux surfaces on which data is desired vs. 
;						poloidal angle (theta) and the field-line label
;						zeta_0 (the toroidal angle at which field line in
;						question passes through the outer midplane).
;
;		IradialNormal		Set to a floating point array containing the values
;						of the poloidal angle at which slices of the data
;						vs. radial coordinate (r?) and toroidal angle (zeta)
;						are desired.
;
;		PoloidalAvg		Set this keyword (e.g., /PoloidalAvg) to get
;						data averaged over fieldlines (over poloidal angles
;						-� to �) vs. radial coordinate (r?) and field-line
;						label zeta_o
;
;		ZetaCut			Set this keyword to a floating point array of zeta
;						values at which cuts of the data vs. radius and
;						poloidal angle are desired.
;
;		Thetacut			Set this keyword to a floating point array of theta 
;						values at which cuts of the data vs. radius and
;						toroidal angle are desired.
;
;		ToroidalAvg		Set this keyword (e.g., /ToroidalAvg) to get data
;						averaged along field line (over one cycle of poloidal
;						angle from -� to �) vs. radius and poloidal angle
;						theta_0 (a field-line label).
;
;		IfluxLimits		Set to a two-element integer array containing the
;						indices of the first and last flux surfaces to be 
;						included in forming IradialNormal or IradailFieldLIne
;						data.
;
;		TimeSteps			Set this keyword to the number of timesteps which you
;						desire in the output file.  If not set, defaults to
;						timeslice from selected Pxxxx file.  If it is set,
;						then you can only get ONE GKVsd object from this call
;						to GTC_Data.  If timeSteps is set to an integer array (and the
;						second value of this array is less than the first element)
;						then the first element determines the number of files to 
;						be read (in total) while the second element dedtermines
;						the interval between files (defaults to 1, or every file)
;						
;		Dt				Set this keyword to the time interval between between
;						GTC data file files.
;
;		R_range			Set this keyword to a two-element floating-point array
;						containing the range of 'r/q' values corresponding
;						to the full range of 'iflux'.  Default is [0.1, 0.9].
;
;		Qprofile			Set this keyword (e.g., /Qprofile) to obtain the
;						q-profile as a GKVsd object.
;
;		IotaProfile		Set this keyword (e.g., /IotaProfile) to obtain the
;						iota-profile as a GKVsd object.
;
;		Q_0				Set this keyword to the value of q on axis. 
;						(defaults to the range 0.5 < Q_0 � 1.0)
;
;		Iota_0			Alternatively, you can set this keyword to the
;						value of iota on axis (defaults to the range
;						1 � iota < 2).
;
;		Variable			Set this keyword to a STRING containing the name of
;						the variable within the Pxxxx file desired.  This
;						defaults to 'potential'.  If 'variable' is set to
;						a non-zero integer, then will attempt to extract
;						data corresponding to the (Variable+1)th variable
;						within the Pxxxx file (set Variable = 0 to get first,
;						and probably only variable in the Pxxxx file).
;
;		ncd				Set this keyword to the suffix of the sequence of 
;						netcdf files.  Defaults to "ncd" (optional).
;
;		Path				Path to begin search for the GTC netcdf files.
;						Defaults to current directory. (Optional)
;
;	Written by W.M. Nevins
;	3/4/00
;
IF(N_ELEMENTS(DeBug) EQ 0) THEN BEGIN				; Set Keyword 'DeBug' to avoid error trap
Catch, error 								; Set up error trap
IF error NE 0 then Begin 
   Catch, /Cancel             					; Cancel error trap 
   ok=Error_Message(/Trace)  					; Print out a trace-back in the output log.
   IF(RUNdimenID GE 0) THEN NCDF_CLOSE, RUNdimenID	; Close any open NCDF files
   IF(PxxxxID    GE 0) THEN NCDF_CLOSE, PxxxxID
   IF(N_ELEMENTS(Result) NE 0) 	THEN	RETURN, result	; Return SOMETHING useful
   IF(N_ELEMENTS(RundimInfo) NE 0) THEN	$		;	(if possible ...)
   	RETURN, RundimInfo
   RETURN, 0                  					; Return 0 if all else fails. 
ENDIF 
ENDIF
;
; Initialize some flags...
;
RUNdimenID = -1
PxxxxID = -1
;
; FIrst find and open 'Rundimen.ncd' file
;
PRINT, 'Select RUNdimen.ncd File'
RUNdimenPath = DIALOG_PICKFILE(Filter='*.ncd', GET_PATH=savePath, PATH=path)
RUNdimenID = NCDF_OPEN(RUNdimenPath, /NoWrite)
;
; Find 'poloidal-grids' and 'index-shift' arrays
;
pgID = NCDF_VARID(RUNdimenID, 'poloidal-grids')
isID = NCDF_VARID(RUNdimenID, 'index-shift')
;
; Read 'poloidal-grids' and 'index-shift' arrays
;
NCDF_VARGET, RUNdimenID, pgID, poloidalGrids
NCDF_VARGET, RUNdimenID, isID, indexShift
;
; Find dimension 'flux-surfaces'
;
FluxSurID = NCDF_DIMID(RUNdimenID, 'flux-surfaces')
NCDF_DIMINQ, RUNdimenID, FluxSurID, FluxSurfaces, nFluxSurfaces
;
; Close Rundimen.ncd
;
NCDF_CLOSE, RUNdimenID
RUNdimenID = -1
;
; Compute poloidal grid index to start and finish of poloidal data for each flux surface
;
finish = LONARR(nFluxSurfaces)
FOR i=0, nFluxSurfaces-1 DO finish[i] = TOTAL(poloidalGrids[0:i])
start =[0,finish[0:nFluxSurfaces-2]]
finish=finish-1
;
; Compute iota-profile [assume that iota(0) is slightly greater than one...]
;
iota_shift=1
IF(N_ELEMENTS(Q0) 	EQ 1) THEN iota_shift = FLOAT(1/Q0)
IF(N_ELEMENTS(Iota0) 	EQ 1) THEN iota_shift = FLOAT(iota0)
iota=FLTARR(nFluxSurfaces)
iota[0] = iota_shift + FLOAT(indexShift[0])/FLOAT(poloidalGrids[0]-1)
FOR i=1, nFluxSurfaces-1 DO BEGIN
	IF((indexShift[i]-indexshift[i-1]) GT poloidalGrids[i]/2) THEN iota_shift=iota_shift-1
	iota[i]=iota_shift + FLOAT(indexShift[i])/FLOAT(poloidalGrids[i]-1)
ENDFOR
q = 1./iota
;
; Get structure to hold information regarding poloidal structure of data in Pxxxx files
;
RundimInfo = {GTC_Str}
;
; Populate with poloidal structure info developed so far
;
RundimInfo.poloidalGrids 	= PTR_NEW(poloidalGrids)
RundimInfo.indexShift		= PTR_NEW(indexShift)
RundimInfo.nRadialGrids	= nFluxSurfaces
RundimInfo.start		= PTR_NEW(start)
RundimInfo.finish		= PTR_NEW(finish)
RundimInfo.iota			= PTR_NEW(iota)
RundimInfo.q			= PTR_NEW(q)
;
; Find and open the 'Pxxxx.ncd' file
;
ncd="ncd"
IF(KEYWORD_SET(fileSuffix)) THEN ncd=fileSuffix
ncd = STRCOMPRESS(ncd, /REMOVE_ALL)
PRINT, "Select 'Pxxxx.ncd' file"
PxxxxPath = DIALOG_PICKFILE(Filter='*.'+ncd, PATH=savePath)
IF(STRCMP(PxxxxPath,'',/FOLD_CASE)) THEN GOTO, NoPxxxxFile
PxxxxID = NCDF_OPEN(PxxxxPath, /NoWrite)
;
; Find dimensions 'poloidal_grids', 'toroidal_grids'
;
PoloidalGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
IF(PoloidalGridID LT 0) THEN PoloidalGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
ToroidalGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
IF(ToroidalGridID LT 0) THEN ToroidalGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
NCDF_DIMINQ, PxxxxID, PoloidalGridID, PoloidalGrid, nPoloidalGrids
NCDF_DIMINQ, PxxxxID, ToroidalGridID, ToroidalGrid, nToroidalGrids
;
; Check that poloidal grid dimension from Pxxxx file is consistant with info from
; RUNdimen file
;
IF(nPoloidalGrids NE TOTAL(poloidalGRids)) THEN BEGIN
	MESSAGE, 'Structure of Pxxxx.ncd file is inconsistent with RUNdimen.ncd', /informational
	NCDF_CLOSE, PxxxxID
	PxxxxID = -1
	GOTO, NoPxxxxFile 
ENDIF
RundimInfo.nToroidalGrids = nToroidalGrids
;
; Find name(s) of variable(s) in Pxxxx file
;
InquireStr = NCDF_INQUIRE(PxxxxID)
nVars = InquireStr.nvars
varNames = STRARR(nVars)
FOR ivarindex=0, nVars-1 DO BEGIN
	varInqStr = NCDF_VARINQ(PxxxxID, ivarindex)
	varNames[ivarindex] = varInqStr.name
ENDFOR
RundimInfo.varNames = PTR_NEW(varNames)
;
; Check for input variable names
;
Variable='Potential'
IF Query_Integer(variableName) THEN Variable = varNames[variableName]
IF  Query_String(variableName) THEN Variable = VariableName
;
; Process /Info requests
;
IF KEYWORD_SET(Info) THEN GOTO, NoPxxxxFile 
;
; Find VariableID for selected variable
;
variableID = NCDF_VARID(PxxxxID, variable)
IF(variableID EQ -1) THEN BEGIN
	MESSAGE, "Can't locate " + Variable, /Informational
	NCDF_CLOSE, PxxxxID
	PxxxxID = -1
	GOTO, NoPxxxxFile
ENDIF
;
; Create GKVsd structure to hold general information about
; 	this data
;
	dataInfo = {GKVsd}					; Get a GKVsd structure to hold data info, and
	dataInfo.Mnemonic	= variable			; Populate GKVsd structure with common values
	IF(variable EQ 'Potential') THEN $
	dataInfo.Title		= '!4u!X'			; a greek 'phi' if variable is 'Potential'
	IF(variable EQ 'Potential') THEN $
	dataInfo.Units		= 'T/e'			; Units for GTC potential as per conversatin with Z. Lin on 3/6/00
	dataInfo.CodeName	= 'GTC'			; Assume data is from GTC code
	dataInfo.CodePi		= 'Z. Lin'			;	and that Zhihong Lin is the code PI
	dataInfo.FileID		= GTC_FileName(PxxxxPath)	
;
; Set R_range
;
R_range = [0.1, 0.9]
IF(N_ELEMENTS(rrange) EQ 2) THEN R_range = rrange
rGrid = R_range[0] + (R_range[1] - R_range[0])*FINDGEN(nFluxSurfaces)/nFluxSurfaces
;
; Process requests for q-profiles
;
IF(KEYWORD_SET(Qprofile) OR KEYWORD_SET(IotaProfile)) THEN BEGIN
	result = { q: OBJ_NEW(), iota: OBJ_NEW() }
	q_str = {GKVs1D}							; Get a GKVs1D structure to hold info on q-profile
	FOR i=0, N_TAGS(dataInfo)-1 DO q_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	q_str.Mnemonic		= 'q'						; Reset mnemonic
	q_str.Title		= 'q'						;	and title
	q_str.Indices		= PTR_NEW(['*'])
	q_str.Values		= PTR_NEW(q)				; Set values field
	q_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	q_str.Grid1.title	= 'r'						;	...
	q_str.Grid1.Units	= 'a'
	q_str.Grid1.Values 	= PTR_NEW(rGrid)
	q_str.Grid1.boundary	= 'Open'
	q_str.Grid1.uniform	= 1B
	q_str.Grid1.range	= R_range
	q_str.Grid1.irange	= [0, nFluxSurfaces-1]
	result.q=OBJ_NEW('GKVs1D', q_str)				; Create GKVs1D object
	result.iota = result.q -> MakeCopy(/NoValues)		; Make deep copy (without values) in iota field
	result.iota -> SET, values=PTR_NEW(iota)			; Set iota values
	result.iota -> SET, Mnemonic='iota', Title='!4i!X'	;	also mnemonic and title
	RETURN, result								; and we're done
ENDIF
;
; Get time step size
;
Dt=1.
IF N_ELEMENTS(Delta_t) THEN Dt = Delta_t
;
; Check if this is request for a time-series
;
IF(N_ELEMENTS(TimeSteps) GE 1) THEN GOTO, TimeSeries
;
; Crack 'PxxxxPath', to get the name of the file
; 	(with extensions stripped off)
;
separator="/"								; Path separator for unix devices
IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
n_strings=N_ELEMENTS(substrings)
FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
subSeparator="."
subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
TimeString = STRING(t0, FORMAT='(g10.3)')
TimeString = STRTRIM(TimeString, 2)
Indices = REPLICATE('*',4)
Indices[3] = 't=' + TimeString
;
; Process IsurfaceNormal requests (data on flux surface vs. two toroidal angles)
;
Nsurfaces = N_ELEMENTS(IsurfaceNormal)
IF(Nsurfaces GT 0) THEN BEGIN
	data_str = {GKVs2D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'theta'					; Set fields of Grid1
	data_str.Grid1.title		= '!4h!X'					;	(a greek 'theta')
	data_str.Grid1.Units		= ''
	data_str.Grid1.boundary	= 'Periodic (closed)'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X'					;	(a greek 'zeta')
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	CASE Nsurfaces OF
		1	:	CommandString = 'result = { ' + Variable + ': OBJ_NEW() }'
		else	:	CommandString = 'result = { ' + Variable + ': OBJARR(Nsurfaces) }'
	ENDCASE
	ok = EXECUTE(CommandString)
	FOR isurface = 0, Nsurfaces - 1 DO BEGIN
		iflux = IsurfaceNormal(isurface)
		IF((iflux LT 0) OR (iflux GE nFluxSurfaces)) THEN BEGIN
			PRINT, 'Flux surface out of range, iflux = ', iflux, ' Range is 0 � iflux � ', nFluxSurfaces
			commandString = 'result.' + Variable + '[i] = OBJ_NEW()
			ok = EXECUTE(commandString)
			GOTO, DONE1
		ENDIF
		tempIndices = Indices
		rvalue = rGrid[iflux]
		iFluxString = STRING(rvalue, FORMAT='(g10.3)')
		iFluxString = STRTRIM(iFluxString, 2)
		tempIndices[0] = 'r/a=' + iFluxString
		data_str.Indices = PTR_NEW(tempIndices)
		surfNorm_str = GTC_Surface_Normal(PxxxxID, variableID, iflux, RundimInfo)
		data_str.Values = surfNorm_str.valuePtr
		data_str.Grid1.Values = surfNorm_str.thetaPtr
		data_str.Grid2.Values = surfNorm_str.zetaPtr
		commandString = "result." + Variable + "[isurface] = OBJ_NEW('GKVs2D', data_str)"
		ok = EXECUTE(commandString)
	DONE1:
	ENDFOR
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process IsurfaceFieldLine requests (data on flux surface vs. two toroidal angles)
;
Nsurfaces = N_ELEMENTS(IsurfaceFieldLine)
IF(Nsurfaces GT 0) THEN BEGIN
	data_str = {GKVs2D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'theta'					; Set fields of Grid1
	data_str.Grid1.title		= '!4h!X'					;	(a greek 'theta')
	data_str.Grid1.Units		= ''
	data_str.Grid1.boundary	= 'Periodic (closed)'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta_0'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X!D0!N'				;	(a greek 'zeta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	CASE Nsurfaces OF
		1	:	CommandString = 'result = { ' + Variable + ': OBJ_NEW() }'
		else	:	CommandString = 'result = { ' + Variable + ': OBJARR(Nsurfaces) }'
	ENDCASE
	ok = EXECUTE(CommandString)
	FOR isurface = 0, Nsurfaces - 1 DO BEGIN
		iflux = IsurfaceFieldLine(isurface)
		IF((iflux LT 0) OR (iflux GE nFluxSurfaces)) THEN BEGIN
			PRINT, 'Flux surface out of range, iflux = ', iflux, ' Range is 0 � iflux � ', nFluxSurfaces
			commandString = 'result.' + Variable + '[i] = OBJ_NEW()
			ok = EXECUTE(commandString)
			GOTO, DONE2
		ENDIF
		tempIndices = Indices
		rvalue = rGrid[iflux]
		iFluxString = STRING(rvalue, FORMAT='(g10.3)')
		iFluxString = STRTRIM(iFluxString, 2)
		tempIndices[0] = 'r/a=' + iFluxString
		data_str.Indices = PTR_NEW(tempIndices)
		surfFieldLine_str = GTC_Surface_FieldLine(PxxxxID, variableID, iflux, RundimInfo)
		data_str.Values = surfFieldLine_str.valuePtr
		data_str.Grid1.Values = surfFieldLine_str.thetaPtr
		data_str.Grid2.Values = surfFieldLine_str.zetaPtr
		commandString = "result." + Variable + "[isurface] = OBJ_NEW('GKVs2D', data_str)"
		ok = EXECUTE(commandString)
	DONE2:
	ENDFOR
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process zetaCut requests
;	(data vs. radius and poloidal angle at fixed toroidal angle)
;
Nzetas = N_ELEMENTS(ZetaCut)
IF(Nzetas GT 0) THEN BEGIN
	data_str = {GKVs2D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;	(radial coordinate)
	data_str.Grid1.Units		= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'theta_0'					; Set fields of Grid2
	data_str.Grid2.title		= '!4h!X!D0!N'				;	(a greek 'theta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic'
	data_str.Grid2.uniform	= 1B
	CASE Nzetas OF
		1	:	CommandString = 'result = { ' + Variable + ': OBJ_NEW() }'
		else	:	CommandString = 'result = { ' + Variable + ': OBJARR(Nzetas) }'
	ENDCASE
	ok = EXECUTE(CommandString)
	FOR izeta = 0, Nzetas - 1 DO BEGIN
		zeta = zetaCut(izeta)
		IF((zeta LT 0) OR (zeta GT 2.*!PI)) THEN BEGIN
			PRINT, 'zeta is out of range, !4f!X = ', zeta, ' Range is -!4p !9i !4f !9i !4p!X', zeta
			commandString = 'result.' + Variable + '[i] = OBJ_NEW()
			ok = EXECUTE(commandString)
			GOTO, DONE3
		ENDIF
		tempIndices = Indices
		izetaString = STRING(zeta, FORMAT='(g10.3)')
		izetaString = STRTRIM(izetaString, 2)
		tempIndices[2] = '!4f!X=' + izetaString
		data_str.Indices = PTR_NEW(tempIndices)
		zetaCutStr = GTC_ToroidalCut(PxxxxID, variableID, zeta, RundimInfo)
		data_str.Values = zetaCutstr.valuePtr
		data_str.Grid1.Values = zetaCutstr.radialPtr
		data_str.Grid2.Values = zetaCutstr.thetaPtr
		commandString = "result." + Variable + "[izeta] = OBJ_NEW('GKVs2D', data_str)"
		ok = EXECUTE(commandString)
	DONE3:
	ENDFOR
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process thetaCut requests
;	(data vs. radius and toroidal angle at fixed poloidal angle)
;
nThetas = N_ELEMENTS(thetaCut)
IF(nThetas GT 0) THEN BEGIN
	data_str = {GKVs2D}						; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'r'					; Set fields of Grid1
	data_str.Grid1.title	= 'r'				;	(radial coordinate)
	data_str.Grid1.Units	= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta'				; Set fields of Grid2
	data_str.Grid2.title	= '!4f!X'				;	(a greek 'zeta')
	data_str.Grid2.Units	= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	CASE nThetas OF
		1	:	CommandString = 'result = { ' + Variable + ': OBJ_NEW() }'
		else	:	CommandString = 'result = { ' + Variable + ': OBJARR(nThetas) }'
	ENDCASE
	ok = EXECUTE(CommandString)
	FOR iTheta = 0, nThetas - 1 DO BEGIN
		theta = thetaCut(itheta)
		IF((theta LT -2.*!PI) OR (theta GT 2.*!PI)) THEN BEGIN
			PRINT, 'theta is out of range, !4h!X = ', theta, ' Range is -!4p !9i !4f !9i !4p!X'
			commandString = 'result.' + Variable + '[i] = OBJ_NEW()
			ok = EXECUTE(commandString)
			GOTO, DONE3a
		ENDIF
		tempIndices = Indices
		iThetaString = STRING(theta, FORMAT='(g10.3)')
		iThetaString = STRTRIM(iThetaString, 2)
		tempIndices[2] = '!4h!X=' + iThetaString
		data_str.Indices = PTR_NEW(tempIndices)
		thetaCutStr = GTC_PoloidalCut(PxxxxID, variableID, theta, RundimInfo)
		data_str.Values = thetaCutstr.valuePtr
		data_str.Grid1.Values = thetaCutstr.radialPtr
		data_str.Grid2.Values = thetaCutstr.zetaPtr
		commandString = "result." + Variable + "[itheta] = OBJ_NEW('GKVs2D', data_str)"
		ok = EXECUTE(commandString)
	DONE3a:
	ENDFOR
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process ToroidalAvg requests
;	(data vs. radius and poloidal angle--as field-line label-- averaged over toroidal angle
;	 ... that is, averaged along the field line.)
;
IF KEYWORD_SET(ToroidalAvg) THEN BEGIN
	data_str = {GKVs2D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Title			= '<' + data_str.Title + '>!I!4f!N!X'	; Set title to <title>_zeta
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;	(radial coordinate)
	data_str.Grid1.Units		= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'theta_0'				; Set fields of Grid2
	data_str.Grid2.title		= '!4h!X!L0!N'				;	(a greek 'theta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'			; '(closed)' means that endpoints are 
	data_str.Grid2.uniform	= 1B
	CommandString = 'result = { ' + Variable +   ': OBJ_NEW(), '	$
	                              + Variable + 'Sq: OBJ_NEW(), '	$
	                              + 'avgTheta: OBJ_NEW(), stdevTheta: OBJ_NEW() }'
	ok = EXECUTE(CommandString)						; Create structure to hold output
	ToroidalAvgStr = GTC_ToroidalAvg(PxxxxID, variableID, RundimInfo)
	data_str.Values = ToroidalAvgStr.valuePtr
	data_str.Grid1.Values = ToroidalAvgStr.radialPtr
	data_str.Grid2.Values = ToroidalAvgStr.thetaPtr
	tempIndices = Indices[1:3]
	data_str.Indices = PTR_NEW(tempIndices)
	valueObj = OBJ_NEW('GKVs2D', data_str)				; Create GKVs2D object containing fieldline average of 'variable'
	commandString = "result." + Variable + " = valueObj"		;
	ok = EXECUTE(commandString)						; 	and store in 'result'
	;
	; Modify 'valueObj' to describe <variable^2>_zeta
	;
	valueSqObj = valueObj -> MakeCopy(/noValues)			; Make a 'deep' copy of valueObj (without copying 'values')
	valueSqObj -> Set,	mnemonic	= Variable + 'Sq',	$	; 	modify mnemonic,
					Title		= '<' + dataInfo.Title + '!E2!N>!D!4f!N!X',	$
												; 	set title to <title^2>_zeta,
					values	= ToroidalAvgStr.valueSqPtr,				$
												;	set 'values',
					units		='(T/e)!E2!N'			;	and set units.
	commandString = "result." + Variable + "Sq = valueSqObj"
	ok = EXECUTE(commandString)						; Store 'variableSq' object in 'result'
	;
	; Modify 'valueOb'j to describe <theta>_variableSQ
	;
	thetaObj = valueObj -> MakeCopy(/noValues)				; Make a 'deep' copy of valueObj (without copying 'values') 
	thetaObj -> Set,	mnemonic	= 'avgTheta',	$		; 	modify mnemonic,
					Title		= '<!4h!X>',	$		; 	set title to <theta>,
					values	= ToroidalAvgStr.avgThetaPtr,			$
												;	set 'values',
					units		= ''					;	and units
	result.avgTheta = thetaObj							; Store thetaObj in 'result'
	;
	; Modify 'valueObj' to describe standard deviations in theta
	;
	stdevThetaObj = valueObj -> MakeCopy(/noValues)			; Make a 'deep' copy of valueObj (without copying 'values') 
	stdevThetaObj -> Set,	mnemonic	= 'StDevTheta',	$		; 	modify mnemonic,
					Title		= '(<!4h!X!E2!N> - <!4h!X>!E2!N)!E1/2!N',	$	
												; 	set title to <theta^2> - <theta>^2,
					values	= ToroidalAvgStr.stdevThetaPtr,			$
												;	set 'values'
					units		= ''					;	and units.
					
	result.stdevTheta = stdevThetaObj					; Store stdevThetaObj in 'result'
	

	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF


;
; Process IradialNormal requests 
;	(data on cut at fixed poloidal angle from magnetic axis to surface)
;
Nangles = N_ELEMENTS(IradialNormal)
IF(Nangles GT 0) THEN BEGIN
	MESSAGE, 'IradialNormal not yet implimented', /Informational
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, 0
ENDIF
;
; Process PoloidalAvg requests
;	(data vs. radius and toroidal angle--as field-line label-- averaged over poloidal angle
;	 ... that is, averaged along the field line.)
;
IF KEYWORD_SET(PoloidalAvg) THEN BEGIN
	data_str = {GKVs2D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Title			= '<' + data_str.Title + '>!I!4h!N!X'	; Set title to <title>_theta
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;	(radial coordinate)
	data_str.Grid1.Units		= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta_0'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X!L0!N'				;	(a greek 'Zeta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'			; '(closed)' means that endpoints are 
	data_str.Grid2.uniform	= 1B
	CommandString = 'result = { ' + Variable +   ': OBJ_NEW(), '	$
	                              + Variable + 'Sq: OBJ_NEW(), '	$
	                              + 'avgTheta: OBJ_NEW(), stdevTheta: OBJ_NEW() }'
	ok = EXECUTE(CommandString)						; Create structure to hold output
	PoloidalAvgStr = GTC_PoloidalAvg(PxxxxID, variableID, RundimInfo)
	data_str.Values = PoloidalAvgStr.valuePtr
	data_str.Grid1.Values = PoloidalAvgStr.radialPtr
	data_str.Grid2.Values = PoloidalAvgStr.ZetaPtr
	tempIndices = Indices[1:3]
	data_str.Indices = PTR_NEW(tempIndices)
	valueObj = OBJ_NEW('GKVs2D', data_str)				; Create GKVs2D object containing fieldline average of 'variable'
	commandString = "result." + Variable + " = valueObj"		;
	ok = EXECUTE(commandString)						; 	and store in 'result'
	;
	; Modify 'valueObj' to describe <variable^2>_zeta
	;
	valueSqObj = valueObj -> MakeCopy(/noValues)			; Make a 'deep' copy of valueObj (without copying 'values')
	valueSqObj -> Set,	mnemonic	= Variable + 'Sq',	$	; 	modify mnemonic,
					Title		= '<' + dataInfo.Title + '!E2!N>!D!4f!N!X',	$
												; 	set title to <title^2>_zeta,
					values	= PoloidalAvgStr.valueSqPtr,				$
												;	set 'values',
					units		='(T/e)!E2!N'			;	and set units.
	commandString = "result." + Variable + "Sq = valueSqObj"
	ok = EXECUTE(commandString)						; Store 'variableSq' object in 'result'
	;
	; Modify 'valueOb'j to describe <theta>_variableSQ
	;
	thetaObj = valueObj -> MakeCopy(/noValues)				; Make a 'deep' copy of valueObj (without copying 'values') 
	thetaObj -> Set,	mnemonic	= 'avgTheta',	$		; 	modify mnemonic,
					Title		= '<!4h!X>',	$		; 	set title to <theta>,
					values	= PoloidalAvgStr.avgThetaPtr,			$
												;	set 'values',
					units		= ''					;	and units
	result.avgTheta = thetaObj							; Store thetaObj in 'result'
	;
	; Modify 'valueObj' to describe standard deviations in theta
	;
	stdevThetaObj = valueObj -> MakeCopy(/noValues)			; Make a 'deep' copy of valueObj (without copying 'values') 
	stdevThetaObj -> Set,	mnemonic	= 'avgThetaSq',	$		; 	modify mnemonic,
					Title		= '(<!4h!X!E2!N> - <!4h!X>!E2!N)!E1/2!N',	$	
												; 	set title to <theta^2> - <theta>^2,
					values	= PoloidalAvgStr.stdevThetaPtr,			$
												;	set 'values'
					units		= ''					;	and units.
					
	result.stdevTheta = stdevThetaObj					; Store stdevThetaObj in 'result'
	

	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF


RETURN, result ; THIS WILL BECOME THE NORMAL RETURN STATEMENT****************************


TimeSeries:
Indices = REPLICATE('*',4)
;
; Process IsurfaceNormal requests (data on flux surface vs. two toroidal angles and time)
;
Nsurfaces = N_ELEMENTS(IsurfaceNormal)
IF(Nsurfaces GT 0) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs3D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'theta'					; Set fields of Grid1
	data_str.Grid1.title		= '!4h!X'					;	(a greek 'theta')
	data_str.Grid1.Units		= ''						; 
	data_str.Grid1.boundary	= 'Periodic (closed)'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X'					;	(a greek 'zeta')
	data_str.Grid2.Units		= ''						; 
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B
	iflux = IsurfaceNormal[0]
	IF((iflux LT 0) OR (iflux GE nFluxSurfaces)) THEN BEGIN
		PRINT, 'Flux surface out of range, iflux = ', iflux, ' Range is 0 � iflux � ', nFluxSurfaces
		InSurObjs(0)=OBJ_NEW()
		Return, 0
	ENDIF
	rvalue = rGrid[iflux]
	iFluxString = STRING(rvalue, FORMAT='(g10.3)')
	iFluxString = STRTRIM(iFluxString, 2)
	tempIndices = Indices
	tempIndices[0] = 'r/a=' + iFluxString
	data_str.Indices = PTR_NEW(tempIndices)
	surfNorm_str = GTC_Surface_Normal(PxxxxID, variableID, iflux, RundimInfo)
	data_str.Grid1.Values = surfNorm_str.thetaPtr
	data_str.Grid2.Values = surfNorm_str.zetaPtr
	surfNorm_str.thetaPtr = PTR_NEW()					; Reload null pointers so grid pointers don't
	surfNorm_str.zetaPtr = PTR_NEW()					; 	get trashed when these get cleaned up
	;
	; Get theta, zeta dimensions
	;
	nTheta = N_ELEMENTS(*data_str.Grid1.Values)
	nZeta  = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index, nFiles
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nTheta, nZeta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"								; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *surfNorm_str.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append a '/' to the front of 'File_Path'
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE10				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, DONE10
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, DONE10
		ENDIF
		;
		; Free up pointers from previous call to GTC_Surface_Normal
		;
		FOR isnElements = 0, N_TAGS(SurfNorm_str)-1  DO PTR_FREE, surfNorm_str.(isnElements)
		;
		; Get next set of data values
		;
		surfNorm_str = GTC_Surface_Normal(PxxxxID, variableID, iflux, RundimInfo)
		Values[*,*,it] = *surfNorm_str.valuePtr		; Append data values at new time step to 'Values' array	
	ENDFOR
DONE10:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process IsurfaceFieldLine requests (data on flux surface vs. zeta_0, an ignorable field-line label,
;	and theta, which measures position along a field line)
;
Nsurfaces = N_ELEMENTS(IsurfaceFieldLine)
IF(Nsurfaces GT 0) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs3D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Grid1.Mnemonic	= 'theta'					; Set fields of Grid1
	data_str.Grid1.title		= '!4h!X'					;	(a greek 'theta')
	data_str.Grid1.Units		= ''						; 
	data_str.Grid1.boundary	= 'Periodic (closed)'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta_0'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X!D0!N'				;	(a greek 'zeta' sub 0)
	data_str.Grid2.Units		= ''						; 
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B
	iflux = IsurfaceFieldLine[0]
	IF((iflux LT 0) OR (iflux GE nFluxSurfaces)) THEN BEGIN
		PRINT, 'Flux surface out of range, iflux = ', iflux, ' Range is 0 � iflux � ', nFluxSurfaces
		InSurObjs(0)=OBJ_NEW()
		Return, 0
	ENDIF
	rvalue = rGrid[iflux]
	iFluxString = STRING(rvalue, FORMAT='(g10.3)')
	iFluxString = STRTRIM(iFluxString, 2)
	tempIndices = Indices
	tempIndices[0] = 'r/a=' + iFluxString
	data_str.Indices = PTR_NEW(tempIndices)
	surfFIeldLIne_str = GTC_Surface_FieldLine(PxxxxID, variableID, iflux, RundimInfo)
	data_str.Grid1.Values = surfFIeldLIne_str.thetaPtr
	data_str.Grid2.Values = surfFIeldLIne_str.zetaPtr
	surfFieldLIne_str.thetaPtr = PTR_NEW()				; Reload null pointers so grid pointers don't
	surfFieldLIne_str.zetaPtr = PTR_NEW()					; 	get trashed when these get cleaned up
	;
	; Get theta, zeta dimensions
	;
	nTheta = N_ELEMENTS(*data_str.Grid1.Values)
	nZeta  = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index (nFiles)
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nTheta, nZeta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"								; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *surfFIeldLIne_str.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append '/' to front of File_Path
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE11				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, DONE11
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, DONE11
		ENDIF
		;
		; Free up pointers from previous call to GTC_Surface_Normal
		;
		FOR isnElements = 0, N_TAGS(surfFIeldLIne_str)-1  DO PTR_FREE, surfFIeldLIne_str.(isnElements)
		;
		; Get next set of data values
		;
		surfFIeldLIne_str = GTC_Surface_FieldLine(PxxxxID, variableID, iflux, RundimInfo)
		Values[*,*,it] = *surfFIeldLIne_str.valuePtr		; Append data values at new time step to 'Values' array	
	ENDFOR
DONE11:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process ZetaCut requests (data at fixed toroidal angle vs. radius and poloidal angle)
;
;
Ncuts = N_ELEMENTS(ZetaCut)
IF(Ncuts GT 0) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs3D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into data_str
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;
	data_str.Grid1.Units		= 'a'						; 
	data_str.Grid1.boundary	= 'open'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'theta'					; Set fields of Grid2
	data_str.Grid2.title		= '!4h!X'					;	(a greek 'theta')
	data_str.Grid2.Units		= ''						; 
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B
	zeta = ZetaCut[0]
	IF((zeta LT 0) OR (zeta GT 2.*!PI)) THEN BEGIN
		PRINT, 'zeta is out of range, !4f!X = ', zeta, ' Range is -!4p !9i !4f !9i !4p!X', zeta
		Return, 0
	ENDIF
		izetaString = STRING(zeta, FORMAT='(g10.3)')
		izetaString = STRTRIM(izetaString, 2)
		tempIndices = Indices
		tempIndices[2] = '!4f!X=' + izetaString
		data_str.Indices = PTR_NEW(tempIndices)
		zetaCutStr = GTC_ToroidalCut(PxxxxID, variableID, zeta, RundimInfo)
		data_str.Values = zetaCutstr.valuePtr
		data_str.Grid1.Values = zetaCutstr.radialPtr
		data_str.Grid2.Values = zetaCutstr.thetaPtr
		zetaCutstr.radialPtr = PTR_NEW()
		zetaCutstr.thetaPtr  = PTR_NEW()
	;
	; Get radial, theta dimensions
	;
	nRadius = N_ELEMENTS(*data_str.Grid1.Values)
	nTheta  = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index (nFiles)
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nRadius, nTheta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"								; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *zetaCutStr.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append '/' to front of File_Path
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE12				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, DONE12
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, DONE12
		ENDIF
		;
		; Free up pointers from previous call to GTC_ToroidalCut
		;
		FOR isnElements = 0, N_TAGS(zetaCutStr)-1  DO PTR_FREE, zetaCutStr.(isnElements)
		;
		; Get next set of data values
		;
		zetaCutStr = GTC_ToroidalCut(PxxxxID, variableID, zeta, RundimInfo)
		Values[*,*,it] = *zetaCutStr.valuePtr		; Append data values at new time step to 'Values' array	
	ENDFOR
DONE12:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process thetaCut requests (data at fixed poloidal angle vs. radius and toroidal angle)
;
;
Ncuts = N_ELEMENTS(thetaCut)
IF(Ncuts GT 0) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs3D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into data_str
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;
	data_str.Grid1.Units		= 'a'						; 
	data_str.Grid1.boundary	= 'open'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X'					;	(a greek 'zeta')
	data_str.Grid2.Units		= ''						; 
	data_str.Grid2.boundary	= 'Periodic (closed)'
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B
	theta = thetaCut[0]
	IF((theta LT -2.*!PI) OR (theta GT 2.*!PI)) THEN BEGIN
		PRINT, 'theta is out of range, theta = ', theta, ' Range is -pi < theta < pi'
		Return, 0
	ENDIF
		ithetaString = STRING(theta, FORMAT='(g10.3)')
		ithetaString = STRTRIM(ithetaString, 2)
		tempIndices = Indices
		tempIndices[2] = '!4h!X=' + ithetaString
		data_str.Indices = PTR_NEW(tempIndices)
		thetaCutStr = GTC_PoloidalCut(PxxxxID, variableID, theta, RundimInfo)
		data_str.Values = thetaCutstr.valuePtr
		data_str.Grid1.Values = thetaCutstr.radialPtr
		data_str.Grid2.Values = thetaCutstr.zetaPtr
		thetaCutstr.radialPtr = PTR_NEW()
		thetaCutstr.zetaPtr  = PTR_NEW()
	;
	; Get radial, zeta dimensions
	;
	nRadius = N_ELEMENTS(*data_str.Grid1.Values)
	nzeta   = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index (nFiles)
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nRadius, nzeta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"							; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)		; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]					; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex			; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)					; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *thetaCutStr.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append '/' to front of File_Path
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE12				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, DONE12
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, DONE12a
		ENDIF
		;
		; Free up pointers from previous call to GTC_ToroidalCut
		;
		FOR isnElements = 0, N_TAGS(zetaCutStr)-1  DO PTR_FREE, zetaCutStr.(isnElements)
		;
		; Get next set of data values
		;
		thetaCutStr = GTC_PoloidalCut(PxxxxID, variableID, theta, RundimInfo)
		Values[*,*,it] = *thetaCutStr.valuePtr			; Append data values at new time step to 'Values' array	
	ENDFOR
DONE12a:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process ToroidalAvg requests
;	(data vs. radius and poloidal angle--as field-line label-- averaged over toroidal angle
;	 ... that is, averaged along the field line.)
;
IF KEYWORD_SET(ToroidalAvg) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Title			= '<' + data_str.Title + '>!I!4f!N!X'	; Set title to <title>_zeta
	tempIndices = Indices[1:3]
	data_str.Indices = PTR_NEW(tempIndices)
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;	(radial coordinate)
	data_str.Grid1.Units		= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'theta_0'				; Set fields of Grid2
	data_str.Grid2.title		= '!4h!X!L0!N'				;	(a greek 'theta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'			; '(closed)' means that endpoints are 
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B

	ToroidalAvgStr = GTC_ToroidalAvg(PxxxxID, variableID, RundimInfo)
	data_str.Values = ToroidalAvgStr.valuePtr
	data_str.Grid1.Values = ToroidalAvgStr.radialPtr
	data_str.Grid2.Values = ToroidalAvgStr.zetaPtr
	ToroidalAvgStr.radialPtr = PTR_NEW()
	ToroidalAvgStr.zetaPtr  = PTR_NEW()
	;
	; Get radial, theta dimensions
	;
	nRadius = N_ELEMENTS(*data_str.Grid1.Values)
	nTheta  = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index (nFiles)
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nRadius, nTheta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"								; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *ToroidalAvgStr.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append '/' to front of File_Path
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE11				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, DONE13
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, DONE13
		ENDIF
		;
		; Free up pointers from previous call to GTC_ToroidalAvg
		;
		FOR isnElements = 0, N_TAGS(zetaCutStr)-1  DO PTR_FREE, ToroidalAvgStr.(isnElements)
		;
		; Get next set of data values
		;
		ToroidalAvgStr = GTC_ToroidalAvg(PxxxxID, variableID, RundimInfo)
		Values[*,*,it] = *ToroidalAvgStr.valuePtr		; Append data values at new time step to 'Values' array	
	ENDFOR
DONE13:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF
;
; Process PoloidalAvg requests
;	(data vs. radius and toroidal angle--as field-line label-- averaged over poloidal angle
;	 ... that is, averaged along the field line.)
;
IF KEYWORD_SET(PoloidalAvg) THEN BEGIN
	data_str = {GKVs3D}								; Get a GKVs2D structure to hold data
	FOR i=0, N_TAGS(dataInfo)-1 DO data_str.(i)=dataInfo.(i)	; Copy dataInfo into q_str
	data_str.Title			= '<' + data_str.Title + '>!I!4h!N!X'	; Set title to <title>_theta
	tempIndices = Indices[1:3]
	data_str.Indices = PTR_NEW(tempIndices)
	data_str.Grid1.Mnemonic	= 'r'						; Set fields of Grid1
	data_str.Grid1.title		= 'r'						;	(radial coordinate)
	data_str.Grid1.Units		= 'a'
	data_str.Grid1.boundary	= 'Neumann'
	data_str.Grid1.uniform	= 1B
	data_str.Grid2.Mnemonic	= 'zeta_0'					; Set fields of Grid2
	data_str.Grid2.title		= '!4f!X!L0!N'				;	(a greek 'Zeta' sub 0)
	data_str.Grid2.Units		= ''
	data_str.Grid2.boundary	= 'Periodic (closed)'			; '(closed)' means that endpoints are 
	data_str.Grid2.uniform	= 1B
	data_str.Grid3.Mnemonic	= 't'						; Set fields of Grid3
	data_str.Grid3.title		= 't'
	data_str.Grid3.Units		= ''						; ***** NEED UNITS! *****
	data_str.Grid3.boundary	= 'Open'
	data_str.Grid3.uniform	= 1B
	
	PoloidalAvgStr = GTC_PoloidalAvg(PxxxxID, variableID, RundimInfo)
	data_str.Values = PoloidalAvgStr.valuePtr
	data_str.Grid1.Values = PoloidalAvgStr.radialPtr
	data_str.Grid2.Values = PoloidalAvgStr.ZetaPtr
	PoloidalAvgStr.radialPtr = PTR_NEW()
	PoloidalAvgStr.zetaPtr  = PTR_NEW()
	;
	; Get radial, zeta dimensions
	;
	nRadius = N_ELEMENTS(*data_str.Grid1.Values)
	nZeta  = N_ELEMENTS(*data_str.Grid2.Values)
	;
	; Work out number of time steps (nT) and number total increment to the file index (nFiles)
	;
	iskip = 1
	nT = 10
	CASE N_ELEMENTS(TimeSteps) OF
		0	:	nFIles = nT-1
		1	:	BEGIN
				nFiles = TimeSteps[0]-1
				nT = TimeSteps[0]
				END
		else	:	BEGIN
				nFiles = TimeSteps[0]
				nT = Timesteps[0]
				IF(TimeSteps[1] LT TimeSteps[0]) THEN BEGIN
					iskip= TimeSteps[1]
					nFiles = iSkip*(TimeSteps[0]-1)
				ENDIF
				END
	ENDCASE
	;
	; Create array to hold data values
	;
	Values = FLTARR(nRadius, nZeta, nT)
	TimeSeries = FLTARR(nT)
	;
	; Crack 'PxxxxPath', to get the name of the file
	; 	(with extensions stripped off)
	;
	separator="/"								; Path separator for unix devices
	IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
	substrings=STRSPLIT(PxxxxPath, separator, /Extract)	; Break file_name into substrings
	n_strings=N_ELEMENTS(substrings)
	FileNameNCD=substrings[n_strings-1]				; Strip leading directory names of filename
	subSeparator="."
	subSubStrings=STRSPLIT(fileNameNCD, subSeparator, /Extract)	; Break filenameNCD at dots
	FIleName=subsubstrings[0]						; Run_name is leading piece of "filename"
	FileNameLen = STRLEN(FileName)					; length of FileName (in characters)
	FileHeadLen = FileNameLen - 4					; Length of ASCII string preceeding index number in file name
	READS, STRMID(FileName, FileHeadLen), FileIndex		; This SHOULD get the 'xxxx' out of the FileName, and put it in FIleIndex
	FileIndex = LONG(FileIndex)					; Force 'FileIndex' to be an integer
	it = 0
	t0 = DT*FLOAT(fileIndex)						; Use 'FileIndex' to guess value of 1st time step
	Dt = Dt*iskip
	timeSeries = t0 + Dt*FINDGEN(nT)				; Generate time grid		
	Values[*,*,it] = *PoloidalAvgStr.valuePtr
	FOR iFileIndex=FileIndex + iskip, FileIndex + nfiles, iskip DO BEGIN 
		it = it+1								; Increment time index
		FileIndexStr = STRING(iFileIndex, FORMAT='(I4)'); Turn FileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)		; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number
		subSubStrings[0] = FileName
		FileNameNCD = STRJOIN(subSubStrings, subSeparator)
		subStrings[n_strings-1] = FileNameNCD
		File_Path = STRJOIN(subStrings, separator)		; Create new file path
		IF(!D.NAME EQ 'X') THEN 	$				; Append '/' to front of File_Path
			File_Path = '/' + File_Path			;	if this is a UNIX device
		NCDF_CLOSE, PxxxxID						; Close previous file
		PxxxxID = -1
		ok = FINDFILE(File_Path)					; Check if next NCDF file exists
		IF(ok[0] EQ '') THEN GOTO, DONE11				; 	(no file found which matches File_Path)		
		PxxxxID = NCDF_OPEN(File_Path, /NoWrite)		; Open new one
		variableID = NCDF_VARID(PxxxxID, variable)
		IF(variableID EQ -1) THEN BEGIN				; Check that it contains 'variable'
			MESSAGE, "Can't locate " + Variable + 'in ' + File_Path, /Informational
			GOTO, Done14
		ENDIF
		PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grids')
		TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grids')
		IF(PGridID LT 0) THEN PGridID = NCDF_DIMID(PxxxxID, 'poloidal_grid')
		IF(TGridID LT 0) THEN TGridID = NCDF_DIMID(PxxxxID, 'toroidal_grid')
		NCDF_DIMINQ, PxxxxID, PGridID, PGrid, nPGrids
		NCDF_DIMINQ, PxxxxID, TGridID, TGrid, nTGrids
		;
		; Check that poloidal & toroidal grid dimensions from this Pxxxx file 
		; is consistant with info from previous file
		;
		IF((nPoloidalGrids NE npGrids) OR (nToroidalGrids NE nTGrids)) THEN BEGIN
			MESSAGE, 'Structure of' + File_Path + 'is inconsistent with previous files', /informational
			GOTO, Done14
		ENDIF
		;
		; Free up pointers from previous call to GTC_PoloidalAvg
		;
		FOR isnElements = 0, N_TAGS(zetaCutStr)-1  DO PTR_FREE, ToroidalAvgStr.(isnElements)
		;
		; Get next set of data values
		;
		PoloidalAvgStr = GTC_PoloidalAvg(PxxxxID, variableID, RundimInfo)
		Values[*,*,it] = *PoloidalAvgStr.valuePtr		; Append data values at new time step to 'Values' array	
	ENDFOR
Done14:	;
		; Now, make pointers, and create the GKVs3D object
		;
	ValuePtr = PTR_NEW(Values)
	tGridPtr = PTR_NEW(timeSeries)
	data_str.Values = ValuePTR
	data_str.Grid3.Values = tGridPTR
	result = OBJ_NEW('GKVs3D', data_str)
	IF(PxxxxID GE 0) THEN NCDF_CLOSE, PxxxxID
	RETURN, result
ENDIF




MESSAGE, 'This time Series request not yet implimented'
RETURN, 0





noPxxxxFIle:
;
; Process abnormal return--no valid Pxxxxfile found, 
; 	requested variable couldn't be found, or 
;	/Info keyword was set
;
IF KEYWORD_SET(Info) THEN BEGIN
	PRINT, 'nToroidalGrids = ', nToroidalGrids
	PRINT, 'nPoloidalGrids = ', TOTAL(poloidalGrids)
	PRINT, ' nFluxSurfaces = ', nFluxSurfaces
	PRINT, '      Variable = ', Variable
	PRINT, 'Variables in file: ', varNames
ENDIF
IF(PxxxxID    GE 0) 	THEN NCDF_CLOSE, PxxxxID
IF(RundimenID GE 0)	THEN NCDF_CLOSE, RundimenID 
RETURN, Rundiminfo
END ; ****** GTC_DATA ****** ;

