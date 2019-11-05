FUNCTION GKVs2D::PROJECT, arg, _Extra=extra
;
; Purpose:
;
;	Projects the GKV object onto selected axis
;	using specified weighting funciton.
;
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
;			project the data.  Defaults to the final axis and
;			current irange. (Optional)
;
;      Weight		Set equal to a GKV object containing the desired weighting
;			function.  Defaults to 1.0.
;
;  Written by W.M. Nevins
;	4/9/04
;
;
; Get axis ID
;
CASE N_PARAMS() OF
	0	:	iaxis = self -> AxisIrange(     _Extra=extra)
	1	:	iaxis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'PROJECTION called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE

nDims = self -> NumDims()
IF(iaxis LT 0) THEN iaxis=nDims
IF( (iaxis GT nDims) OR (iaxis EQ 0) ) THEN BEGIN
	MESSAGE, 'PROJECTION called with illegal axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Call "PROJECT" method of appropriate dimensionality
;
CASE iaxis OF
	1	:	result = self -> PROJECT1(_EXTRA=extra)
	2	:	result = self -> PROJECT2(_EXTRA=extra)
	3	:	result = self -> PROJECT3(_EXTRA=extra)
	ELSE	:	BEGIN
				MESSAGE, 'PROJECT called with axis identifier GT 3', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE


RETURN, result
END ; ****** GKVs2D::PROJECTs ****** ;
