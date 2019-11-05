FUNCTION GKVs1D::Hilbert, arg, Rotation=d, _Extra=extra
;
; Purpose:
;
;		This function returns the Hilbert transform of 'self'
;		(see discussion of the Hilbert transform in the IDL
;		reference manual).  
;
;	Arguments:
;
;		A legal axis identifier. Default to Grid1.
;		(Optional)
;
;	Keywords:
;
;	Rotation	The direction of the 90 degree phase rotation.
;			Defaults to 1.  Alternative is -1. (Optional)
;
;	{axisId}	The range in the independent variable 'self'
;			over which the envelop is computed can be
;			set using any legal axis ID.  Defaults to
;			current signal window. (Optional)

; Written by W.M. Nevins
;	9/30/00
;
ndims = self -> NumDims()
IF(ndims NE 1) THEN BEGIN
	MESSAGE, 'Hilbert only implimented for GKVs1D objects', /INFORMATIONAL
	RETURN, 0
ENDIF
rotation=1
IF KEYWORD_SET(d) THEN BEGIN
	IF(d EQ -1) THEN rotation=-1
ENDIF
;
copy = self -> MakeCopy()
;
; Use AxisIrange to reset signal window of 'self' if necessary
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Hilbert called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
copy -> Restrict
valuePtr = copy -> GetValues()
values = *valuePtr
newValues = HILBERT(values, rotation)
PTR_FREE, copy.values
copy.values = PTR_NEW(newValues)
copy.title = 'H{' + copy.title + '}'
vmin = GKVsd_MIN(newvalues, Max=vmax)
copy.vrange=[vmin,vmax]
RETURN, copy
END ; ****** GKVs1D::Hilbert ****** ;

Function GKVs1D::Envelope, arg, _Extra=extra
;
;  Purpose:
;
;		This function returns an 'Envelope' for the
;		(persumably) oscillatory signal in 'self'.
;		The envelope is produced using a Hilbert 
;		transform.
;
;  Arguments:
;
;			A legal axis identifier. Defaults to Grid1.
;			(Optional)
;
;  Keywords:
;
;	{axisId}	The range in the independent variable 'self'
;			over which the envelop is computed can be
;			set using any legal axis ID.  Defaults to
;			current signal window. (Optional)
;
;  Written by W.M. Nevins
;	9/30/00
ndims = self -> NumDims()
IF(ndims NE 1) THEN BEGIN
	MESSAGE, 'Envelope only implimented for GKVs1D objects', /INFORMATIONAL
	RETURN, 0
ENDIF
;
copy = self -> MakeCopy()
;
; Use AxisIrange to reset signal window of 'self' if necessary
;
CASE N_PARAMS() OF
	0	:	axis = self -> AxisIrange(     _Extra=extra)
	1	:	axis = self -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'Envelope called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
copy -> Restrict
valuePtr = copy -> GetValues()
values = *valuePtr
hValues = HILBERT(values)
newValues = values + COMPLEX(0.,1.)*hValues
newValues = ABS(newValues)
PTR_FREE, copy.values
copy.values = PTR_NEW(newValues)
copy.title = 'E!DH!N{' + copy.title + '}'
vmin = GKVsd_MIN(newvalues, Max=vmax)
copy.vrange=[vmin,vmax]
RETURN, copy
END ; ******  GKVs1D::Envelope ****** ;

FUNCTION GKVs1D_FullWidths, ObjArr, arg, _Extra=extra
;
;  Purpose:
;
;			This function accepts an array of GKVs1D objects
;			as the argument.  It returns the full width
;			at half maximum (but, see KeyWord 'Fraction' below) of
;			the two-time correlation function for each of these
;			objects.  Optionally, the correlation functions are 
;			displayed.
;
; Arguments:
;
;	ObjArr	An array of GKVs1D objects whose correlation time/length
;			is to be computed. REQUIRED
;
;			A legal axis identifier. Defaults to Grid1.
;			(Optional)
;
;  Keywords:
;
;	{axisId}	The range in the independent variable 'self'
;			over which the envelop is computed can be
;			set using any legal axis ID.  Defaults to
;			current signal window. (Optional)
;
;	Indices	Set this keyword to a two element integer
;			array to specify the starting and final
;			elements of ObjArr to examine.  Defaults
;			to all elements of ObjArr. (Optional)
;
;	Fraction	If this keyword is set, then GKVs1D_FullWidths computes
;			the full width at 'fraction' of the maximum.
;			Defaults to 1/2. (Optional)
;
;	View		Set this keyword (that is, put '/view' on the 
;			command line) to display correlation and envelope
;			of each object in an individual window.  
;			Default is not to display intermediates.
;			(Optional)
;
;			If 'View' is set, then any additional keywords
;			provided on the command line will be forwarded
;			to the plotting routine.
;
;  Written by W.M. Nevins
;	9/30/00
nObjs = N_ELEMENTS(objArr)
IF(nObjs EQ 0) THEN BEGIN
	MESSAGE, "No arguments provided", /INFORMATIONAL
	RETURN, 0
ENDIF
IF(TypeOf(Objarr) NE 11) THEN BEGIN
	MESSAGE, "First argument must be an array of GKVs1D objects", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Set up array to hold the computed fullwidths
;
fullWidths = FLTARR(nObjs)
;
; Check for 'View' keyword
;
view=0
result = GetKeyWord('View', extra)
IF(TypeOf(result) NE 7) THEN view=KEYWORD_SET(result)
;
; Check for 'Fraction' keyword
;
fraction = 0.5
result = GetKeyWord('Fraction', extra)
IF(TypeOf(result) NE 7) THEN fraction = result < 0.995
;
; Check for 'Rotation' keyword
;
rotation = 1
result = GetKeyWord('Rotation', extra)
IF(TypeOf(result) NE 7) THEN rotation=result
;
; Check of Indices keyword
;
indices=[0, nObjs-1]
result = GetKeyWord('Indices', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	CASE N_ELEMENTS(result) OF
		0:	
		1:	indices[0]= FIX(result)
		2:	indices   = FIX(result)
	  ELSE:
	  ENDCASE
ENDIF 
;
; Use AxisIrange to reset signal window of 'self' if necessary
;
FOR i=indices[0], indices[1] DO BEGIN
	copy = ObjArr[i] -> MakeCopy()
	IF( NOT Query_Integer(extra) ) THEN cExtra = extra
	CASE N_PARAMS() OF
		1	:	axis = copy -> AxisIrange(     _Extra=cExtra)
		2	:	axis = copy -> AxisIrange(arg, _Extra=cExtra)
		else	:	BEGIN
					MESSAGE, 'GKVs1D_FullWidths called with too many arguments', /INFORMATIONAL
					RETURN, 0
			END
	ENDCASE
	copy -> Restrict
	copy -> get, axis=1, uniform=gridUniform
	IF(NOT gridUniform) THEN BEGIN	
		copy -> set, axis=1, uniform=1B
		PRINT, 'Warning, grid is not uniform.  Will force correlations'
	ENDIF
	corrs = copy -> xcorr()
	Ecorrs = corrs -> Envelope(rotation=rotation)
	FullWidths[i] = Ecorrs -> FullWidth(fraction=fraction)
	PRINT, i, FullWidths[i]
	IF KEYWORD_SET(view) THEN BEGIN
		corrs -> view, _Extra=cExtra
		Ecorrs -> oPlot
	ENDIF ELSE BEGIN
	corrs  -> Trash
	Ecorrs -> Trash 
	ENDELSE
	copy   -> Trash
ENDFOR
RETURN, fullWidths
END ; ****** GKVs1D_FullWidths ****** ;
