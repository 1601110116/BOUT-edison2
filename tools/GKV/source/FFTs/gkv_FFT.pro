FUNCTION GKVs1D::FFT, arg, _EXTRA=extra
;
;  Purpose:

; 	Transforms data from/to a FOURIER representation
; 	vs. the specified independent variable
;
; Arguments:
;
;			Any legal axis identifier. Defaults to the final axis.
;			The independent variable may also be identified using an
;			Axis mnemonic. (Optional)
;
;
; Input Keywords:
;
;   'mnemonic'		Where 'mnemonic' is the mnemonic for the independent variable
;			over which the fourier transform is to be applied.  'Mnemonic' should be 
;			set equal to the range in this variable over which you wish to 
;			transform the data.  Defaults to the final axis and
;			current irange. (Optional)
;
;   Inverse		'Set' this keyword (i.e., put '/Inverse') to perform an inverse
;			Fourier transform vs. the indicated variable.  Defaults to forward
;			transform.  (Optional)
;
;    Offset		Set this keyword to the (floating point) offset value which the specified
;			independent variable will assume on completion of the inverse transform.
;			Defaults to zero.  (Optional)
;
;  Written by W.M. Nevins
;	11/25/01
;
; Get axis ID
;
CASE N_PARAMS() OF
	0	:	iaxis = self -> AxisIrange(     _Extra=extra)
	1	:	iaxis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'FFT called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LE 0) THEN BEGIN
	MESSAGE, 'FFT called with axis identifier LE 0', /INFORMATIONAL
	RETURN, 0
ENDIF
CASE iaxis OF
	1	:	result = self -> FFT1(_EXTRA=extra)
	2	:	result = self -> FFT2(_EXTRA=extra)
	3	:	result = self -> FFT3(_EXTRA=extra)
	4	:	result = self -> FFT4(_Extra=extra)
	ELSE	:	BEGIN
				MESSAGE, 'FFT called with axis identifier GT 4', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE

RETURN, result
END ; ****** GKVs1D::FFT ****** ;
