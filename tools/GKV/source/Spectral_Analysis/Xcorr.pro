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

FUNCTION Query_Integer, ARG
;
; Returns 1 if ARG is an integer, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 0 ) then return, 0
if (arg_type lt 4 ) then return, 1
if (arg_type lt 12) then return, 0
if (arg_type lt 15) then return, 1
return, 0
END


Function Query_Real, ARG
;
; Returns 1 if ARG is a real number, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 4 ) then return, 1
if (arg_type eq 5 ) then return, 1
return, 0
END


Function Query_Complex, ARG
;
; Returns 1 if ARG is a Complex number, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 6 ) then return, 1
if (arg_type eq 9 ) then return, 1
return, 0
END


Function Query_String, ARG
;
; Returns 1 if ARG is a string, 0 otherwise
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 7 ) then return, 1
return, 0
END


Function Window_Fcn, N_VALS, IDW, DEBUG=d, DOUBLE=Dbl
;
; Compute window function for convolution with data.
;	N_VALS is (integer) length of window
;	IDW is (ingeter) length of roll-off at each end
;	DEBUG is flag to turn on/of diagnostics
;	if DOUBLE is set, then returns double-precision values
;
; Error traps
;
;if (N_PARAMS(0) ne 2) then begin
;	PRINT, "ERROR in Window_Fcn:  called with too few/too many arguments"
;	return, R ; will return undefined result
;endif
;if(Query_Integer(n_vals) ne 1) then begin
;	PRINT, "ERROR in Window_Fcn: N_VALS is not an integer"
;	return, R ; will return undefined result
;endif
;if(Query_Integer(idw) ne 1) then begin
;	PRINT, "ERROR in Window_Fcn: IDW is not an integer"
;	return, R ; will return undefined result
;endif

if(idw gt n_vals/2 ) then idw = n_vals/2
;
; Check return array packed with '1' if idw=0
;
IF(idw EQ 0) THEN RETURN, REPLICATE(1.0,n_Vals)
;
;  Compute window function
;

j = LINDGEN(n_vals)
dw = (j LE idw)*0.5*(1. - COS(!PI*j/idw)) + (J GT idw)	$
    - (j GT (n_vals-1L-idw))*0.5*(1. + COS(!PI*(n_vals -1L-j)/idw))

;if KEYWORD_SET(Dbl) then begin
;	dw= DOUBLE(dw)
;	for j = 0L, n_vals-1 do begin
;		if ( j lt idw ) then dw(j) = 0.5d - 0.5d*cos(!dpi*j/idw)
;		if ( j gt n_vals-1L-idw ) then dw(j) = 0.5d - 0.5d*cos(!dpi*(n_vals-1L-j)/idw)
;	endfor
;endif
;for j = 0L, n_vals-1 do begin
;	if ( j lt idw ) then dw(j) = 0.5 - 0.5d*cos(!pi*j/idw)
;	if ( j gt n_vals-1L-idw ) then dw(j) = 0.5 - 0.5*cos(!pi*(n_vals-1L-j)/idw)
;endfor

;if KEYWORD_SET(d) then plot, dw
return, dw
END


Function LagWindow_Fcn, N_VALSin, ILWin, DEBUG=d, DOUBLE=Dbl
;
; Computes lag  window for convolution with correlation function.
; N_VALS is (integer) length of array to be returned
; ILW is (integer) length of of the lagWindow
; DEBUG is flag to turn on/of diagnostics
; if DOUBLE is set, then returns double-precision values
;
; Error traps
; (disabled 6/20/01 to speed this function up!)
;
;if (N_PARAMS(0) ne 2) then begin
;	PRINT, "ERROR in LagWindow_Fcn:  called with too few/too many arguments"
;	return, R ; will return undefined result
;endif
;if(Query_Integer(n_vals) ne 1) then begin
;	PRINT, "ERROR in LagWindow_Fcn: N_VALS is not an integer"
;	return, R ; will return undefined result
;endif
;if(Query_Integer(ilw) ne 1) then begin
;	PRINT, "ERROR in LagWindow_Fcn: ILW is not an integer"
;	return, R ; will return undefined result
;endif
ilw=ilwIn
n_Vals = n_ValsIn
IF(ilw EQ 0) THEN BEGIN
	lagWindow = REPLICATE(1.0, n_vals)
	RETURN, LagWindow
ENDIF
if(ilw gt n_vals/2 ) then ilw = n_vals/2
;
;  Compute lag window
;
centerPoint = n_vals/2
;IF KEYWORD_SET(dbl) THEN BEGIN
;	j = DINDGEN(n_vals) - centerPoint
;	lagWindow = 0.5d*( 1.0d + COS(!DPI*j/ilw) )*( ABS(j) LE ilw )
;ENDIF ELSE BEGIN
	j = FINDGEN(n_vals) - centerPoint
	lagWindow = 0.5*( 1.0 + COS(!PI*j/ilw) )*( ABS(j) LE ilw )
;ENDELSE

;if KEYWORD_SET(d) then plot, lagWindow
return, lagWindow
END



Function MyHanning, N_Vals, Weight, DOUBLE=Dbl
;
; A Hanning window on [0, N_Vals-1] which is Symmetric about the element indexed (N_Vals/2 - 1)
; 	i.e., Window(0) = Window(N_Vals-2) ­ 0 (but close to 0) while  = Window(N_Vals-1) = 0 (exactly).
;
; Weight is weighting function allowing change to HAMMING window (with same symmetry)
;
; if the keyword DOUBLE is set, then a double-precision window is returned.
;
;	Check arguments
;
n_args = N_PARAMS()
if( n_args eq 0) then begin
	PRINT, "ERROR in MyHanning:  No arguemnts passed"
	return, error
endif

if( (Query_Integer(N_vals) + Query_Real(N_vals)) ne 1) then begin
	PRINT, "ERROR in MyHanning:  N_VALS neither a Real nor an Integer"
	return, error
endif

j = LINDGEN(N_vals)

if KEYWORD_SET(Dbl) then begin
	W = 0.5d
	if( (Query_Integer(Weight) + Query_Real(Weight)) ne 0) then W=Double(Weight)
	Window = DBLARR(N_vals) + W
	Window = Window - (1.0d - W )*Cos(2.0d*!dpi*j/N_Vals)
	return, Window
endif

W = 0.5
if( (Query_Integer(Weight) + Query_Real(Weight)) ne 0) then W = Weight
Window = FLTARR(N_vals) + W
Window = Window - (1.0 - W )*Cos(2.0*!pi*j/N_Vals)
return, Window

END


Function XSpect, Input_Signal, DT=Delta_T, Reference = ref,			$
	No_Avg=No_Average, Hamming=Ham, Hanning=Han,  DW=Data_Window, 		$
	LW=lagWindow, DEBUG=DeBug, Double=Dbl, Ptr=Ptr, _Extra=e
;
; Computes cross spectrum between input_signal and 'ref'.  If no 'ref'
; is supplied, then computes spectral density of input_signal
;
;	Written by W.M. Nevins
;		2/7/00
;
; First compute cross-correlation function
;
FORWARD_FUNCTION XCorr
DT=1.
IF(KEYWORD_SET(Delta_T)) THEN DT = Delta_T 
;
; Check validity of input_signal
;
info=SIZE(Input_Signal)
ndims=info[0]
signal_type = info[ndims+1]
IF(signal_type EQ 10) THEN BEGIN		; input_signal is a pointer
	temp = *input_signal			; Dereference pointer
	info=SIZE(temp)
	ndims = info[0]
	signal_type = info[ndims+1]
	ptrflag=1
ENDIF
nt_in = info(nDims)
;
; Check if lag window is set.  If not, the use default of nt_in/2
;
IF(N_ELEMENTS(lw) EQ 0) THEN lw=LONG(nt_in/2L)
;
; Compute correlation function
;
temp = XCorr(Input_Signal, DT=Delta_T, Reference = ref,		$
	No_Avg=No_Average, Hamming=Ham, Hanning=Han,  DW=Data_Window, 		$
	LW=lagWindow, DEBUG=DeBug, Double=Dbl, _Extra=e)

rflag=N_ELEMENTS(ref)
;
; Check validity of cross correlation function
;
info=SIZE(temp)
ndims=info[0]
signal_type = info[ndims+1]
IF((ndims GT 4) OR (ndims LT 1)) THEN BEGIN
	MESSAGE, '1, 2, 3, or 4-D array expected as input', Informational = DeBug
	RETURN, 0
ENDIF
IF(ndims gt 1) THEN n_x = info[ndims-1]
IF(ndims gt 2) THEN n_y = info[ndims-2]
IF(ndims gt 3) THEN n_z = info[ndims-3]
n_t=info[ndims]
n_points=info[ndims+2]
;
; Check for 'Norm' keyword
;
norm = GetKeyWord('Norm', e)
IF Query_Integer(norm) THEN BEGIN
	IF(NORM EQ 1)  THEN BEGIN
		CASE ndims OF
			1:	tNorm = temp[n_t/2]
			2:	tNorm = temp[n_x/2, n_t/2]
			3:	tNorm = temp[n_x/2, n_y/2, n_t/2]
			4:	tNorm = temp[n_x/2, n_y/2, n_z/2, n_t/2]
		ENDCASE
		temp = temp/tNorm
	ENDIF
ENDIF
;
; Cross spectrum is returned from XCorr with values for 
; negative frequencies (and negative wavenumbers) at front of the "temp" array 
; We now repack cross correlation function into 'spect' such that values 
; for negative lags and dispalcements come at the end of the array;
;
CASE ndims OF
	1:	spect = SHIFT(temp, -n_t/2)
	2:	spect = SHIFT(temp, -n_x/2, -n_t/2)
	3:	spect = SHIFT(temp, -n_y/2, -n_x/2, -n_t/2)
	4:	spect = SHIFT(temp, -n_z/2, -n_y/2, -n_x/2, -n_t/2)
ENDCASE
;
; Transform back to frequency space to form the spectral density
;
temp = FFT(spect, -1, DOUBLE=Dbl)
;
; Cross spectrum is returned from FFT with values for 
; negative frequencies (and negative wavenumbers) at end of the "temp" array 
; (see desciption of FFT function in IDL function reference document)
; We now repack cross spectrum into 'spect' such that values 
; for negative lags and dispalcements come before postive lags;
;
CASE ndims OF
	1:	spect = SHIFT(temp, n_t/2)
	2:	spect = SHIFT(temp, n_x/2, n_t/2)
	3:	spect = SHIFT(temp, n_y/2, n_x/2, n_t/2)
	4:	spect = SHIFT(temp, n_z/2, n_y/2, n_x/2, n_t/2)
ENDCASE
;
; free up some space...
;
temp = 0.
IF(rflag EQ 0) THEN spect = FLOAT(spect)	; Make cross spectrum function real...
;
; The IDL FFT routine takes both space and time transforms as exp(-i...), and
; inverts to exp(+i...).  The general convention of physicists is to take the
; sign of 'i' in the time transforms opposite that of the sign of 'i' in the 
; space transforms.  We fix that her with a call to "REVERSE" (assuming that
; the final independent variable corresponds to time/frequency).
;
temp = REVERSE(spect, ndims)
CASE ndims OF
	1:	spect = SHIFT(temp, 1)
	2:	spect = SHIFT(temp, 0, 1)
	3:	spect = SHIFT(temp, 0, 0, 1)
	4:	spect = SHIFT(temp, 0, 0, 0, 1)
ENDCASE


IF(KEYWORD_SET(Ptr)) THEN BEGIN
	CrossSpectPtr = PTR_NEW(spect)
	RETURN, CrossSpectPtr
ENDIF
RETURN, spect


IF(KEYWORD_SET(Ptr)) THEN RETURN, PTR_NEW(spect)
return, spect
END ; ****** Xspect ****** ;


Function XCORR, Input_Signal, DT=Delta_T, Reference = ref,				$
	No_Avg=No_Average, Hamming=Ham, Hanning=Han,  DW=Data_Window, 		$
	LW=lagWindow, DEBUG=DeBug, Double=Dbl, Ptr=Ptr, _Extra=e
;---------------------------------------------------------------------------
; Written and debugged by WM Nevins, 4/30/99.
;
; Modified to allow 0, 1, 2 or 3 spatial dimensions (plus time)
; by WM Nevins, 8/29/99
;
; Modified 2/6/00 by WM Nevins to improve readability and make self=contained
;
; Computes cross-correlation function of INPUT_SIGNAL (which is left unaltered by call to XCORR),
; Stores result in CrossCorr_Fcn, with lag-values in Tau and spatial displacements in DELTA*****
; Assumes data is uniformally spaced in time, with time interval
;
; 	DT (defaults to 1.0). 
;
; Assumes rectangular grid 
; On return, the array TAU contains the time lags (from -tau_max to tau_max-DT)
; where 2*taumax is DT times the length of the array into which data was packed
; (generally, first power of 2 larger than number of time steps in input data); 
; while DELTA contains spatial displacements [from -(NR/2)*DR to (NR/2 - 1)*DR].
; Periodic boundary conditions are assumed on spatial variables.
;
; If you set the keyword:
;
;	Ptr,		function returns a pointer to the cross correlation function
;			(instead of returning the values of the cross correlation function)
;
;	No_Avg, 	then average and trend of data are removed
; 	Hamming, 	then Hamming window (in time) is applied to data
;	Hanning, 	then Hanning window (in time) is applied to data
;	DW, 		then cosine roll-off of length DW (if DW is real)
;			or of length DT*DW (if DW is an integer) in time is applied to data.
;	DOUBLE	then double precision arithematic is used
;	DEBUG, 	then function quits on decting bad inputs, etc. 
;
;---------------------------------------------------------------------------------------
one = COMPLEX(1.0, 0.0)
IF KEYWORD_SET(Dbl) THEN one = dCOMPLEX(1.0, 0.0)
zero = COMPLEX(0.0, 0.0)
IF KEYWORD_SET(Dbl) THEN zero = DCOMPLEX(0.0, 0.0)
DT=1.
IF(KEYWORD_SET(Delta_T)) THEN DT = Delta_T 
rflag=N_ELEMENTS(ref)
IF(rflag EQ 0) THEN refReal = 1
ptrflag=0
;
; Check validity of input_signal
;
info=SIZE(Input_Signal)
ndims=info[0]
signal_type = info[ndims+1]
IF(signal_type EQ 10) THEN BEGIN		; input_signal is a pointer
	temp = one*(*input_signal)		; Dereference pointer
	info=SIZE(temp)
	ndims = info[0]
	signal_type = info[ndims+1]
	ptrflag=1
ENDIF
IF((ndims GT 4) OR (ndims LT 1)) THEN BEGIN
	MESSAGE, '1, 2, 3, or 4-D array expected as input', Informational = DeBug
	RETURN, 0
ENDIF
IF(ndims gt 1) THEN n_x = info[ndims-1]
IF(ndims gt 2) THEN n_y = info[ndims-2]
IF(ndims gt 3) THEN n_z = info[ndims-3]
n_t=info[ndims]
nt_in = info[nDims]
n_points=info[ndims+2]
;
; check TYPE of input signal
;
IF(ptrflag) THEN BEGIN					; Check for real/complex input signal
	IF((Query_Real(temp) + Query_Complex(temp)) NE 1) THEN BEGIN
		MESSAGE, 'Input signal must be real or complex valued', Informational = DeBug
		RETURN, 0
	ENDIF
	inReal = Query_Real(temp) 
ENDIF ELSE BEGIN
	IF ((Query_Real(Input_Signal) + Query_Complex(Input_Signal)) NE 1) THEN BEGIN
		MESSAGE, 'Input signal must be real or complex valued', Informational = DeBug
		RETURN, 0
	ENDIF
	inReal = Query_Real(Input_Signal)
ENDELSE

IF(N_ELEMENTS(temp) EQ 0) THEN	$		; Pack Input_signal in to complex array,
	temp = one*Input_Signal			; Coherce temp to double precision if necessary.
;
; form data window
;
IF KEYWORD_SET(Han) THEN BEGIN		; Use Hanning window
	dw = Hanning(n_t)
	GOTO, Norm
ENDIF
IF KEYWORD_SET(Ham) THEN BEGIN		; Use Hamming window
	dw= Hanning(n_t, ALPHA=0.54)
	GOTO, Norm 
ENDIF
idw = LONG(n_t/10L)
IF ( Query_Real(Data_Window) + Query_Integer(Data_Window) ) THEN BEGIN
	; Data_Window is either a real or an integer.  Hence, it has been set.
	; Note that KEYWORD_SET returns FALSE if Data_Window is "set" to zero!
	if ( Query_Real(Data_Window)    eq 1 ) then idw=LONG(Data_Window/dt) 
	if ( Query_Integer(Data_Window) eq 1 ) then idw=LONG(Data_Window)
ENDIF
dw = Window_Fcn(n_t, idw, DOUBLE=Dbl, Debug=d)	; Default is 0.5*(1-cos ...) roll-off of width idw
Norm: $
window_norm = SQRT(n_t/TOTAL(dw*dw))
dw = dw*window_norm				; Normalize data window

IF(ndims GT 1) THEN	$
	dw = REPLICATE(1,n_x)#dw		; Broadcast data window into n_x by n_t array
newWindowNorm = TOTAL(dw)

CASE ndims OF					; Apply data window to input signal
	1:	temp = temp*dw
	2:	temp = temp*dw
	3:	FOR j=0L, n_y-1 DO temp(j,*,*) = temp(j,*,*)*dw
	4:	BEGIN					; two sequential 1-D loops instead of a nested 2-D loop
			ddw=replicate(1.0, n_y, n_x, n_t)
			FOR j=0L,n_y-1 DO ddw(j,*,*) = ddw(j,*,*)*dw
			FOR k=0L,n_z-1 DO temp(k,*,*,*) = temp(k,*,*,*)*ddw
		END
ENDCASE


IF KEYWORD_SET(No_Average) THEN BEGIN
	avg=TOTAL(temp)						; compute average
	CASE ndims OF						; and then remove average
		1:	temp = temp - avg*dw/newWindowNorm	; using data window to insure that average
		2:	temp = temp - avg*dw/(newWindowNorm)	; is REALLY zero AND that signal vanishes outside of data window
		3:	FOR j=0L,n_y-1 DO temp[j,*,*] = temp[j,*,*] - avg*dw/(newWindowNorm*n_y)
		4:	FOR k=0L,n_z-1 DO temp[k,*,*,*] = temp[k,*,*,*] - avg*ddw/(newWindowNorm*n_y*n_z)
	ENDCASE
ENDIF
;
; Check for reference signal
;
IF(rflag NE 0) THEN BEGIN
	rinfo = SIZE(ref)
	rdims = rinfo[0]
	ref_type = rinfo[rdims+1]
	IF(ref_type EQ 10) THEN BEGIN			; reference signal is a pointer
		rtemp = *ref				; Dereference pointer
		rinfo = SIZE(rtemp)
		rdims = rinfo[0]
		ref_type = rinfo[rdims+1]
	ENDIF ELSE BEGIN
		rtemp = ref
	ENDELSE
	FOR i=0,rdims DO BEGIN
		IF(rinfo[i] NE info[i]) THEN BEGIN
			MESSAGE, 'Dimensions of reference signal do not match input signal', Informational = DeBug
			RETURN, 0
		ENDIF
	ENDFOR
	IF ( (ref_type LE 3) OR ((ref_type GE 7) AND (ref_type NE 9))) THEN BEGIN
; Reference signal is neither real nor complex
		MESSAGE, "Reference signal must be real or complex valued"
		RETURN, 0
	ENDIF
	IF(N_ELEMENTS(rtemp) EQ 0) THEN	$	; Pack reference signal into complex array,
		rtemp = one*ref				; Coherce rtemp to double precision if necessary.
	refReal = Query_Real(rtemp)
ENDIF


IF(rflag NE 0) THEN BEGIN				; Apply data window to reference signal
	CASE ndims OF					
		1:	rtemp = rtemp*dw
		2:	rtemp = rtemp*dw
		3:	FOR j=0L,n_y-1 DO rtemp(j,*,*) = rtemp(j,*,*)*dw
		4:	BEGIN				; two sequential 1-D loops instead of a nested 2-D loop
				ddw=replicate(1.0, n_y, n_x, n_t)
				FOR j=0L,n_y-1 DO ddw(j,*,*) = ddw(j,*,*)*dw
				FOR k=0L,n_z-1 DO rtemp(k,*,*,*) = rtemp(k,*,*,*)*ddw
				ddw=0				; Release ddw
			END
	ENDCASE

	IF(KEYWORD_SET(No_Average)) THEN BEGIN	
		ravg = TOTAL(rtemp)/newWindowNorm
		CASE ndims OF					
			1:	rtemp = rtemp - ravg*dw
			2:	rtemp = rtemp - ravg*dw
			3:	FOR j=0L,n_y-1 DO rtemp(j,*,*) = rtemp(j,*,*) - ravg*dw
			4:	BEGIN				; two sequential 1-D loops instead of a nested 2-D loop
				ddw=replicate(1.0, n_y, n_x, n_t)
				FOR j=0L,n_y-1 DO ddw(j,*,*) = ddw(j,*,*)*dw
				FOR k=0L,n_z-1 DO rtemp(k,*,*,*) = rtemp(k,*,*,*) - ravg*ddw
				ddw=0				; Release ddw
				END
		ENDCASE
	ENDIF
ENDIF



j=0L
repeat j=j+1 until (2L^j gt n_t) 	; Find n_corrs such that 2^n_corrs > n_t
n_corrs = 2L^j
pack = GetKeyWord('pack', e)
IF Query_Integer(pack) THEN BEGIN	; Add an extra power of 2 if requested (so there is no wrapping of correlations at large lag)
	IF(pack EQ 1)  THEN n_corrs = n_corrs*2
ENDIF
;
; pack input signal into complex array of length n_corrs
;
CASE ndims OF		; Form array with time dimension = n_corrs
	1:	BEGIN
			CrossCorr_fcn = COMPLEXARR(n_corrs)
			IF KEYWORD_SET(Dbl) THEN CrossCorr_fcn = DCOMPLEXARR(n_corrs)
		END
	2:	BEGIN
			CrossCorr_fcn = COMPLEXARR(n_x, n_corrs)
			IF KEYWORD_SET(Dbl) THEN CrossCorr_fcn = DCOMPLEXARR(n_x, n_corrs)
		END
	3:	BEGIN
			CrossCorr_fcn = COMPLEXARR(n_y, n_x, n_corrs)
			if KEYWORD_SET(Dbl) then CrossCorr_fcn = DCOMPLEXARR(n_y, n_x, n_corrs)
		END
	4:	BEGIN
			CrossCorr_fcn = COMPLEXARR(n_z, n_y, n_x, n_corrs)
			IF KEYWORD_SET(Dbl) THEN CrossCorr_fcn = DCOMPLEXARR(n_z, n_y, n_x, n_corrs)
		END
ENDCASE
CASE ndims OF		; Pack data into CrossCorr_fcn
	1:	CrossCorr_fcn(0:n_t-1) = temp(0:n_t-1)
	2:	CrossCorr_fcn(*,0:n_t-1) = temp(*,0:n_t-1)
	3:	CrossCorr_fcn(*,*,0:n_t-1) = temp(*,*,0:n_t-1)
	4:	CrossCorr_fcn(*,*,*,0:n_t-1) = temp(*,*,*,0:n_t-1) 
ENDCASE
temp = FFT(CrossCorr_fcn, -1, DOUBLE=Dbl)	; Form (forward) fourier transform in all dimension
IF(rflag NE 0) THEN BEGIN					; Pack reference signal into CrossCorr_fcn
	CASE ndims OF					; 	(which is being used as a tempory storage at the moment)
		1:	CrossCorr_fcn(0:n_t-1) = rtemp(0:n_t-1)
		2:	CrossCorr_fcn(*,0:n_t-1) = rtemp(*,0:n_t-1)
		3:	CrossCorr_fcn(*,*,0:n_t-1) = rtemp(*,*,0:n_t-1)
		4:	CrossCorr_fcn(*,*,*,0:n_t-1) = rtemp(*,*,*,0:n_t-1) 
	ENDCASE
	rtemp = FFT(CrossCorr_fcn, -1, DOUBLE=Dbl)	; Form (forward) fourier transform in all dimension
ENDIF
;
; Compute the correlation function via fourier transform method 
; 
IF (rflag NE 0)	THEN	BEGIN
	CrossCorr_fcn = temp*CONJ(rtemp)
ENDIF 	ELSE	BEGIN 
	CrossCorr_fcn = temp*CONJ(temp)		; Take absolute value squared
ENDELSE
IF(rflag NE 0) THEN rtemp = 0.				; release rtemp
;
temp = FFT(CrossCorr_fcn,  1, DOUBLE=Dbl)	; Invert transform to form two-time, two-point correlation function
IF(inReal*refReal) THEN BEGIN			; Convert to real if both inputs were real
	IF(KEYWORD_SET(Dbl)) THEN temp = DOUBLE(temp) ELSE temp = FLOAT(temp)
ENDIF
;
; Correlation function is returned from FFT with values for 
; negative lags (and negative displacements) at end of the "temp" array 
; (see desciption of FFT function in IDL function reference document)
; We now repack correlation function into corr_fcn such that values 
; for negative lags and dispalcements come before postive lags;
; and also correct the normalization.
;
CASE ndims OF
	1:	CrossCorr_fcn = SHIFT(temp, n_corrs/2)
	2:	CrossCorr_fcn = SHIFT(temp, n_x/2, n_corrs/2)
	3:	CrossCorr_fcn = SHIFT(temp, n_y/2, n_x/2, n_corrs/2)
	4:	CrossCorr_fcn = SHIFT(temp, n_z/2, n_y/2, n_x/2, n_corrs/2)
ENDCASE
;
; form lag window
;
IF KEYWORD_SET(Han) THEN BEGIN		; Use Hanning window
	lw = Hanning(n_t)
	GOTO, GotLW
ENDIF
IF KEYWORD_SET(Ham) THEN BEGIN		; Use Hamming window
	lw= Hanning(n_t, ALPHA=0.54)
	GOTO, GotLW 
ENDIF
ilw = 0		; If user deos not set the lag window (that is, keyword 'lw')  then none is applied
IF ( Query_Real(lagWindow) + Query_Integer(lagWindow) ) THEN BEGIN
	; lagWindow is either a real or an integer.  Hence, it has been set.
	; Note that KEYWORD_SET returns FALSE if lagWindow is "set" to zero!
	if ( Query_Real(lagWindow)    eq 1 ) then ilw=LONG(lagWindow/dt) 
	if ( Query_Integer(lagWindow) eq 1 ) then ilw=LONG(lagWindow)
ENDIF
lw = LagWindow_Fcn(n_corrs, ilw, DOUBLE=Dbl, Debug=d)	; Default is 0.5*(1-cos ...) roll-off of width ilw
GotLW: $
IF(ndims GT 1) THEN	$
	lw = REPLICATE(1,n_x)#lw		; Broadcast lag window into n_x by n_t array

CASE ndims OF					; Apply lag window to input signal
	1:	temp = CrossCorr_fcn*lw
	2:	temp = CrossCorr_fcn*lw
	3:	FOR j=0L, n_y-1 DO temp(j,*,*) = CrossCorr_fcn(j,*,*)*lw
	4:	BEGIN					; two sequential 1-D loops instead of a nested 2-D loop
			llw=replicate(1.0, n_y, n_x, n_corrs)
			FOR j=0L,n_y-1 DO llw(j,*,*) = llw(j,*,*)*lw
			FOR k=0L,n_z-1 DO temp(k,*,*,*) = CrossCorr_fcn(k,*,*,*)*llw
			llw=0				; Release llw
		END
ENDCASE

CrossCorr_fcn = (FLOAT(n_corrs)/n_t)*temp		; Correct normalization of CrossCorr_fcn
IF(KEYWORD_SET(Ptr)) THEN BEGIN
	CrossCorrPtr = PTR_NEW(CrossCorr_fcn)
	RETURN, CrossCorrPtr
ENDIF


RETURN, CrossCorr_fcn
END ; ****** XCORR ****** ;
