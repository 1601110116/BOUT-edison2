FUNCTION GKVs1D::CHI_MIX, tauCorr, rho_star=rhoStar, NoErrorBars=noErrorBars
;
; Purpose:
;
;	This function computes the mixing length Chi (as a function of 
;	the inhomogeneous coordinate, 'r').  The radial correlation length
;	is assumed to be in 'self', while the argument is the correlation
;	time.  RhoStar used to correct the normalization.
;
; Written by W.M. Nevins
;	3/5/02
;
const=1.0
IF(N_ELEMENTS(rhoStar) EQ 1) THEN const=1./(rhoStar*rhoStar)
temp1=self -> times(self)
temp2=temp1 -> over(tauCorr)
result=temp2 -> times(const)

temp1 -> trash
temp2 -> trash

result.title='!4v!X!Dmix!N'
result.mnemonic = 'CHi_mix'
result.units = '!4q!X!S!Ds!U!R!U2!Nc!Ds!N/a'

IF(KEYWORD_SET(NoErrorBars)) THEN BEGIN
	PTR_FREE, result.ErrorBars
	result.ErrorBars = PTR_NEW()
ENDIF

return, result

END  ; ****** GKVs1D::CHI_MIX ******  ;


FUNCTION GKV_CopyStructure, arg
;
; This function returns a 'deep' copy of its argument, which
; is assumed to be a structure
;
; Written by W.M. Nevins
;  9/30/03
;
IF(TypeOf(arg) NE 8) THEN BEGIN
	MESSAGE, "Argument is not a structure, returning", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Create output structure
;
Output = arg
;
; Get some basic info
;
nTags = N_TAGS(arg)
tagNames = TAG_NAMES(arg)
;
; process tags, putting copies of all GKVsd objects in the OUTPUT structure
;
FOR i=0, nTags-1 DO BEGIN
	value = arg.(i)
	IF( TypeOf(value) EQ 11 ) THEN BEGIN
		n = N_ELEMENTS(value)
		IF(n EQ 1) THEN BEGIN
			value = value -> MakeCopy()
		ENDIF ELSE BEGIN
			array=value
			value = OBJARR(n)
			FOR j=0,n-1 DO value[j] = array[j] -> MakeCopy()
		ENDELSE
	Output.(i) = value
	ENDIF
ENDFOR
RETURN, Output
END  ; ****** FUNCTION GKV_CopyStructure ****** ;


FUNCTION GYRO_GKNorm, 	TauCorrStr, RhoStar=rho_star, q=q_safety, AspectRatio=A_mid, 	$
			PhiAvg=Phi_Avg, FluxTube=FluxTube, PhiNormed=PhiNormed
;
; Just the same as the proceedure with same name (see below),
; except that it leaves input structure unaltered, while
; returning an output structure in which dependent and
; dependent variables have been rescaled to appropriate
; gyro kinetic normalizations.
;
Output = GKV_CopyStructure(TauCorrStr)
GYRO_GKNorm, Output, RhoStar=rho_star, q=q_safety, AspectRatio=A_mid, PhiAvg=Phi_Avg, FluxTube=FluxTube, PhiNormed=PhiNormed
RETURN, Output
END  ;  ****** FUNCTION GYRO_GKNorm ****** ;



PRO GYRO_GKNorm, 	TauCorrStr, RhoStar=rho_star, q=q_safety, AspectRatio=A_mid, 	$
			PhiAvg=Phi_Avg, FluxTube=FluxTube, PhiNormed=PhiNormed
;
; This proceedure accepts the output structucture from
; TauCorrs in the usual GYRO variables.  It re-scales
; each object such that both independent and dependent
; variables have an appropriate gyro-kinetic normalization
;
;	Arguments:
;
;		TauCorrstr	The output structure produced by the proceedure TauCorrs
;				On return both the dependent and independent variables
;				of the GKV objects within this structure have a gyro-
;				kinetic normalization.
;
;	Input Keywords: 
;
;		RhoStar		The value of RhoStar for this GYRO run.  If a number greater
;				than one is supplied it is assumed to be the reciprocal of rhoStar.
;
;		q		The safety factor vs. radius (in units of r/a)
;
;	AspectRatio		The value of the aspect ratio (the full aspect ratio, R/a)
;
;		
;
; Written by W.M. Nevins
;  9/29/03
;
; First check for valid arguments and keywords
;
nArgs = N_PARAMS()
IF(nArgs NE 1) THEN BEGIN
	MESSAGE, 'GYRO_GKNorm called with wrong number of arguments -- returning', /INFORMATIONAL
	RETURN
ENDIF
IF(TypeOF(q_safety) EQ 0) THEN BEGIN
	MESSAGE, 'GYRO_GKNorm called without specifying q -- returning', /INFORMATIONAL
ENDIF

IF(TypeOf(q_safety) NE 11) THEN BEGIN
	MESSAGE, 'q is not an object -- returning', /INFORMATIONAL
	RETURN	
ENDIF

IF(NOT OBJ_ISA(q_safety, 'GKVs1D')) THEN BEGIN
	MESSAGE, 'q is not a GKVs1D object -- returning', /INFORMATIONAL
	RETURN
ENDIF

IF(N_ELEMENTS(rho_star) NE 1) THEN BEGIN
	MESSAGE, 'GYRO_GKNorm called without specifying RhoStar -- returning', /INFORMATIONAL
	RETURN
ENDIF
rhoStar = rho_star
IF(rho_star GT 1.) THEN rhoStar=1./rho_star

IF(N_ELEMENTS(A_mid) NE 1) THEN BEGIN
	MESSAGE, 'GYRO_GKNorm called without specifying Amid -- returning', /INFORMATIONAL
	RETURN
ENDIF
;
; Now, compute (r-dependent) scaling factor
; which takes you from toroidal angle to distance 
; perpendicular to B:
;
;	dperp/rho = (R/a)*(a/rho)/SQRT{1 + [(R/r)*q(r)]^2}
;
;
q_safety -> Get, values=qptr			; get pointer to q-values
q_safety -> Get, axis=1, gridValues=rptr	; get pointer to r-values (in units of 1/a)
q=*qptr						; dereference pointers
r=*rptr
IF(KEYWORD_SET(FluxTube)) THEN r=0.5+0.0*r
Amid=A_mid					; Amid is really the FULL aspect ratio, R/a
BoverBt = SQRT(1 + ( (r/Amid)*(1/q) )^2 )
BoverBtPtr = PTR_NEW(BoverBT)
BoverBt = q_safety -> MakeCopy(/Novalues)
BoverBt -> set, values=BoverBtPtr					
scale_factor = ((Amid+r)/rhoStar)/SQRT(1 + ( (Amid/r)*q )^2 )
scalePtr = PTR_NEW(scale_Factor)			; make pointer to scale factor values
scaleFactor = q_safety -> MakeCopy(/NoValues)		; and the make new object containing
scaleFactor -> set, values=scalePtr		; the scale factor vs. r.
;
; See what is in TauCorrStr
;
nTags = N_TAGS(TauCorrStr)
tagNames = TAG_NAMES(TauCorrStr)
IF( TOTAL(STRMATCH(tagNames, 'NAME', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.Name = 'Normed' + TauCorrStr.Name
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'TAUCORR', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	temp = BoverBt
	BoverBt = BoverBt -> Interpolate(TauCorrStr.TAUCORR)
	temp -> trash
	rScaleFactor = ScaleFactor -> Interpolate(TauCorrStr.TAUCORR)
	rScaleFactor -> get, values=rScalePtr
	rScale_factor = *rScalePtr
	n_rValues = N_ELEMENTS(rScale_Factor)
ENDIF ELSE BEGIN
	MESSAGE, 'No tauCorr object in this structure, returning', /INFORMATIONAL
	RETURN
ENDELSE

IF( TOTAL(STRMATCH(tagNames, 'VPHASE', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.vphase -> Get, title=title, mnemonic=mnemonic
	temp = TauCorrStr.vphase
	TauCorrStr.vphase = temp -> times(scaleFactor)
	temp -> trash
	TauCorrStr.vphase -> set, title=title, mnemonic=mnemonic, units='(!4q!X!Ds!N/a)*c!Ds!N)'
ENDIF 

IF( TOTAL(STRMATCH(tagNames, 'YCORR', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.YCORR -> Get, title=title, mnemonic=mnemonic
	temp = TauCorrStr.YCORR
	TauCorrStr.YCORR = temp -> times(scaleFactor)
	temp -> trash
	TauCorrStr.YCORR -> set, title='!12l!X!D!9x!X!N', mnemonic='perp_corr', units='!4q!X!Ds!N'
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'YHALF', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.YHALF -> Get, title=title, mnemonic=mnemonic
	temp = TauCorrStr.YHALF
	TauCorrStr.YHALF = temp -> times(scaleFactor)
	temp -> trash
	TauCorrStr.YHALF -> set, title='!12l!X!S!D!9x!X!R!U(1/2)!N', mnemonic='perp_1/2', units='!4q!X!Ds!N'
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'KAVG', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.KAVG -> Get, title=title, mnemonic=mnemonic
	temp = TauCorrStr.KAVG
	TauCorrStr.KAVG = temp -> Over(scaleFactor)
	temp -> trash
	TauCorrStr.KAVG -> set,	title='!12<!Xk!d!9x!X!N!12>!X', 	$
				mnemonic='k_perp', 			$
				units='1/!4q!X!Ds!N',			$
				vrange=[0.,1.]
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'RCORR', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.RCORR -> Get, title=title, mnemonic=mnemonic
	temp = TauCorrStr.RCORR
	TauCorrStr.RCORR = temp -> Over(rhoStar)
	temp -> trash
	TauCorrStr.RCORR -> set, title=title, mnemonic=mnemonic, units='!4q!X!Ds!N'
	ChiMix = tauCorrStr.RCORR -> Chi_Mix(tauCorrStr.taucorr)
	TauCorrStr = CREATE_STRUCT(tauCorrStr, "ChiMix", ChiMix)	
ENDIF

IF( KEYWORD_SET(PhiNormed) NE 1 )THEN BEGIN
IF( TOTAL(STRMATCH(tagNames, 'AVGSQ', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.AVGSQ -> Get, title=title, mnemonic=mnemonic, units=units
	temp = TauCorrStr.AVGSQ
	TauCorrStr.AVGSQ = temp -> Over(rhoStar*rhoStar)
	temp -> trash
	TauCorrStr.AVGSQ -> set, title=title, mnemonic=mnemonic, units='(!4q!X!Ds!N/a)!U2!N*' + units
ENDIF
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'KSPECT', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	temp = TauCorrStr.KSPECT
	i=n_rValues/2
	temp -> scaleAxis, 1, const=1/rScale_factor[i], title='k!D!9x!X!N', mnemonic='k_perp', units='1/!4q!X!Ds!N'
	IF(KEYWORD_SET(PhiNormed)) THEN BEGIN
		TauCorrStr.KSPECT = temp -> times(rScale_factor[i], units='((!4q!X!Ds!N/a)*(T/e))!U2!N')
	ENDIF ELSE BEGIN
		TauCorrStr.KSPECT = temp -> times(rScale_factor[i]/rhoStar^2, units='((!4q!X!Ds!N/a)*(T/e))!U2!N')
	ENDELSE
	temp -> trash
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'YCORRFCN', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	i=n_rValues/2
	TauCorrStr.YCORRFCN -> scaleAxis, 1, const=rScale_factor[i], title='!12l!D!9x!X!N', mnemonic='d_perp', units='!4q!X!Ds!N'
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'RCORRFCN', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	TauCorrStr.RCORRFCN -> scaleAxis, 1, const=1/rhoStar, units='!4q!X!Ds!N'
ENDIF


IF( TOTAL(STRMATCH(tagNames, 'RCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	nPoints = N_ELEMENTS(TauCorrStr.RCorrArr)
	FOR i=0,nPoints-1 DO TauCorrStr.RCorrArr[i] -> scaleAxis, 1, const=1/RhoStar, units='!4q!X!Ds!N'
ENDIF



IF( TOTAL(STRMATCH(tagNames, 'KSPECTARR', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	nPoints = N_ELEMENTS(TauCorrStr.KSPECTARR)
	IF(nPoints EQ n_rValues) THEN BEGIN
		FOR i=0,nPoints-1 DO BEGIN
			temp = TauCorrStr.KSPECTARR[i]
			temp -> scaleAxis, 1, const=1/rScale_factor[i], title='k!D!9x!X!N', mnemonic='k_perp', units='1/!4q!X!Ds!N'
			IF(KEYWORD_SET(PhiNormed)) THEN BEGIN
				TauCorrStr.KSPECTARR[i] = temp -> times(rScale_factor[i], units='((!4q!X!Ds!N/a)*(T/e))!U2!N')
			ENDIF ELSE BEGIN
				TauCorrStr.KSPECTARR[i] = temp -> times(rScale_factor[i]/rhoStar^2, units='((!4q!X!Ds!N/a)*(T/e))!U2!N')
			ENDELSE
			TauCorrStr.KSPECTARR[i] -> Set, units='((!4q!X!Ds!N/a)*(T/e))!U2!N'
			temp -> trash
		ENDFOR
	ENDIF
	kspectAvg = TauCorrStr.KSPECTARR[0] -> ObjAvg(TauCorrStr.KSPECTARR)
	TauCorrStr = CREATE_STRUCT(tauCorrStr, "kSpectAvg", kSpectAvg)
ENDIF


IF( TOTAL(STRMATCH(tagNames, 'YCORRARR', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	nPoints = N_ELEMENTS(TauCorrStr.YCORRARR)
	IF(nPoints EQ n_rValues) THEN BEGIN
		FOR i=0,nPoints-1 DO BEGIN
			TauCorrStr.YCORRARR[i] -> scaleAxis, 1, const=rScale_factor[i], title='!12l!D!9x!X!N', mnemonic='d_perp', units='!4q!X!Ds!N'
		ENDFOR
	ENDIF
	yCorrAvg = TauCorrStr.YCORRARR[0] -> ObjAvg(TauCorrStr.YCORRARR)
	TauCorrStr = CREATE_STRUCT(tauCorrStr, "yCorrAvg", yCorrAvg)
ENDIF

IF( TypeOf(Phi_Avg) EQ 11) THEN BEGIN		; Compute ExB velocity
	PhiAvgDims = Phi_Avg -> NumDims()
	CASE PhiAvgDims OF
		1	: 	PhiAvg = Phi_Avg -> MakeCopy()
		2	:	PhiAvg = Phi_Avg -> Avg(axis=2)
		ELSE	:	GOTO, DONE
	ENDCASE
	PhiAvg -> restrict
	temp = PhiAvg -> dbyd('r')
	V_ExB = temp -> over(BoverBt)
	temp -> trash
	IF(KEYWORD_SET(PhiNormed)) THEN BEGIN
		temp = V_ExB
		V_ExB = temp -> times(rhoStar)
		temp -> trash
	ENDIF
	V_ExB -> Set, title='V!DExB!N', mnemonic='V_ExB', units='(!4q!X!Ds!N/a)*c!Ds!N'
	PhiAvg -> trash
	TauCorrStr = CREATE_STRUCT(tauCorrStr, "V_ExB", V_ExB)	
ENDIF
DONE	:

test =	  TOTAL(STRMATCH(tagNames, 'RCORR', /FOLD_CASE) )	$
	+ TOTAL(STRMATCH(tagNames, 'YHALF', /FOLD_CASE) )	$
	+ TOTAL(STRMATCH(tagNames, 'avgSq', /FOLD_CASE) )
IF(test EQ 3) THEN BEGIN	; Compute eddy turn over time
	taucorrstr.AVGSQ -> get, values=phiSqPtr
	taucorrstr.RCORR -> get, values=rCorrPtr
	taucorrstr.yHalf -> get, values=yHalfPtr
	BoverBt -> get, values=Bptr
	phirms = SQRT(*phiSqPtr)
	tauEddy = !PI*(*yHalfPtr)*(*rCorrPtr)*(*Bptr)/phiRMS
	tauEddyPtr = PTR_NEW(tauEddy)
	tauEddy = taucorrstr.avgsq -> MakeCopy(/NoValues, /NoErrorbars)
	tauEddy -> Set, values=tauEddyPtr, title='!4s!X!DEddy!N', mnemonic='TauEddy', units='a/c!Ds!N'
	TauCorrStr = CREATE_STRUCT(tauCorrStr, "TauEddy", TauEddy)	
ENDIF	
;
; clean up
;
ScaleFactor -> Trash
BoverBt -> Trash

RETURN 
END  ;  ****** PRO GYRO_GKNorm ****** ;
