FUNCTION REPAIR_NormedTauCorr, tauCorrStr, _Extra=Extra
IF(TypeOF(tauCorrStr) NE 8) THEN BEGIN
	MESSAGE, "Argument is not a structure, returning", /INFORMATIONAL
	RETURN, 0
ENDIF
result = GKV_CopyStructure(tauCorrStr)
REPAIR_NormedTauCorr, result, _Extra=Extra
RETURN, result
END

PRO REPAIR_NormedTauCorr, tauCorrStr, rhoStar=rho_Star, q=q_safety, AspectRatio=Aspect_Ratio
;
; This proceedure accepts a "normed" taucorrs output structure
; as input, and recomputes:
;
;	taucorr
;	rcorr
;	yHalf
;	ycorr
;	ChiMix
;	tauEddy
;
; using the corrected version of fullwidth (bug patch of 2/9/04)
;
; Written by W.M. Nevins
;    2/9/04
;
IF(TypeOF(tauCorrStr) NE 8) THEN BEGIN
	MESSAGE, "Argument is not a structure, returning", /INFORMATIONAL
	RETURN
ENDIF
tagNames = TAG_NAMES(tauCorrStr)

newChiMix=0
IF( TOTAL(STRMATCH(tagNames, 'ChiMix', /FOLD_CASE)) EQ 1 ) THEN newChiMix=1
newTauEddy=0
message = ""
;
; Check for rhoStar keyword
;
rhoStar=1.
IF(N_ELEMENTS(rho_Star) NE 0) THEN rhoStar=rho_star
IF(rhoStar GT 1.) THEN rhoStar=1./rhoStar
;
; Construct BoverBt if Aspect_Ratio and q are provided
;
BoverBt=1.
IF( (TypeOf(Aspect_Ratio)*TypeOf(q_safety)) NE 0) THEN BEGIN
	IF(Aspect_Ratio LT 1.) THEN Aspect_Ratio = 1/Aspect_ratio
	IF(TypeOf(q_safety) EQ 11) THEN BEGIN
		q_safety -> Get, values=qptr			; get pointer to q-values
		q_safety -> Get, axis=1, gridValues=rptr	; get pointer to r-values (in units of 1/a)
		q=*qptr						; dereference pointers
		r=*rptr
		AspectRatio = Aspect_Ratio			; AspectRatio is the full aspect ratio, R/a
		BoverBt = SQRT(1 + ( (r/AspectRatio)*(1/q) )^2 )
		BoverBtPtr = PTR_NEW(BoverBT)
		BoverBt = q_safety -> MakeCopy(/Novalues, /NoErrorBars)
		BoverBt -> set, values=BoverBtPtr
		temp = BoverBt
		BoverBt = BoverBt -> Interpolate(TauCorrStr.TAUCORR)
		temp -> trash					
	ENDIF ELSE BEGIN
		BoverBt = SQRT(1 + ( (0.5/Aspect_Ratio)*(1/q_safety) )^2 )
	ENDELSE
ENDIF
;
; Check for tauCorr function array
;
IF( TOTAL(STRMATCH(tagNames, 'TauCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	tauCorr = tauCorrStr.taucorr
	tauCorr -> Get, values=tauCorrPtr, ErrorBars=ErrorPtr
	tauValues = *tauCorrPtr
	nPoints = N_ELEMENTS(tauValues)
	tauErrors = FLTARR(nPoints) 
	FOR i=0, nPoints-1 DO 	$
		tauValues[i] = tauCorrStr.TauCorrArr[i] -> FullWidth(Error=tauErrors[i])
	PTR_FREE, tauCorrPtr
	tauCorr -> Set, values = PTR_NEW(tauValues), ErrorBars=PTR_NEW(tauErrors)
	IF(PTR_VALID(ErrorPtr)) THEN PTR_FREE, ErrorPtr
	message = message + " updated tauCorr"
ENDIF
;
; Check for rCorr function array
;
IF( TOTAL(STRMATCH(tagNames, 'rCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	rCorr = tauCorrStr.rcorr
	rCorr -> Get, values=rCorrPtr, ErrorBars=ErrorPtr
	rValues = *rCorrPtr
	nPoints = N_ELEMENTS(rValues)
	rErrors = FLTARR(nPoints) 
	FOR i=0, nPoints-1 DO 	$
		rValues[i] = tauCorrStr.rCorrArr[i] -> FullWidth(Error=rErrors[i])
	tauCorrStr.rCorrArr[0] -> get, axis=1, Gridunits=rUnits
	IF(STRCMP(rUnits, "a", 1)) THEN BEGIN
		rValues = rValues/rhoStar
		rErrors = rErrors/rhoStar
	ENDIF
	PTR_FREE, rCorrPtr
	rCorr -> Set, values = PTR_NEW(rValues), ErrorBars=PTR_NEW(rErrors), units='!4q!X!Ds!N'
	IF(PTR_VALID(ErrorPtr)) THEN PTR_FREE, ErrorPtr
	message = message + ", rCorr"
	newChiMix=1
	newTauEddy = 1
ENDIF
;
; Check for yCorr function array
;
IF( TOTAL(STRMATCH(tagNames, 'yCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	yHalf = tauCorrStr.yHalf
	yHalf -> Get, values=yHalfPtr, ErrorBars=yHalfErrorPtr
	yValues = *yHalfPtr
	nPoints = N_ELEMENTS(yValues)
	yHalfErrors = FLTARR(nPoints) 
	FOR i=0, nPoints-1 DO 	$
		yValues[i] = tauCorrStr.yCorrArr[i] -> FullWidth(Error=yHalfErrors[i])
	PTR_FREE, yHalfPtr
	yHalf -> Set, values = PTR_NEW(yValues), ErrorBars=PTR_NEW(yHalfErrors)
	IF(PTR_VALID(yHalfErrorPtr)) THEN PTR_FREE, yHalfErrorPtr
	message = message + ", yHalf"
;
; Now, on to re-computing yCorr
;
	yCorr = tauCorrStr.yCorr
	yCorr -> Get, values=yCorrPtr, ErrorBars=yCorrErrorPtr
	yCorrErrors = FLTARR(nPoints)
	FOR i=0, nPoints-1 DO BEGIN
		envelope = tauCorrStr.yCorrArr[i] -> Envelope()
		yValues[i] = envelope -> FullWidth(error=yCorrErrors[i])
	ENDFOR
	PTR_FREE, yCorrPtr
	yCorr -> Set, values = PTR_NEW(yValues), ErrorBars=PTR_NEW(yCorrErrors)
	IF(PTR_VALID(yCorrErrorPtr)) THEN PTR_FREE, yCorrErrorPtr
	message = message + ", yCorr"
	NewTauEddy = 1
ENDIF
;
; Check if a new ChiMix needs to be computed
;
IF(newChiMix) THEN BEGIN
	ChiMix = tauCorrStr.rCorr -> Chi_Mix(tauCorrStr.tauCorr)
	IF( TOTAL(STRMATCH(tagNames, 'yCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
		tauCorrStr.ChiMix -> Trash
		tauCorrStr.ChiMix = ChiMix
	ENDIF ELSE BEGIN
		tauCorrStr = MAKE_STRUCT(tauCorrStr, "ChiMix", ChiMix)
	ENDELSE
	message = message + ", ChiMix"
ENDIF
;
; Check if a new tauEddy needs to be computed
;
IF(newTauEddy) THEN BEGIN
	taucorrstr.AVGSQ -> get, values=phiSqPtr
	taucorrstr.RCORR -> get, values=rCorrPtr
	taucorrstr.yHalf -> get, values=yHalfPtr
	phirms = SQRT(*phiSqPtr)
	tauEddy = !PI*(*yHalfPtr)*(*rCorrPtr)/phiRMS
	IF(TypeOF(BoverBt) NE 11) THEN tauEddy = tauEddy*BoverBt
	tauEddyPtr = PTR_NEW(tauEddy)
	tauEddy = taucorrstr.avgsq -> MakeCopy(/NoValues, /NoErrorbars)
	tauEddy -> Set, values=tauEddyPtr
	IF(TypeOf(BoverBt) EQ 11) THEN BEGIN
		temp = tauEddy
		tauEddy = tauEddy -> times(BoverBt)
		temp -> Trash
	ENDIF
	tauEddy -> Set, title='!4s!X!DEddy!N', mnemonic='TauEddy', units='a/c!Ds!N'
	IF( TOTAL(STRMATCH(tagNames, 'tauEddy', /FOLD_CASE)) EQ 1 ) THEN BEGIN
		tauCorrStr.tauEddy -> Trash
		tauCorrStr.tauEddy = tauEddy
	ENDIF ELSE BEGIN
		TauCorrStr = CREATE_STRUCT(tauCorrStr, "TauEddy", TauEddy)
	ENDELSE
		message = message + ", tauEddy"
ENDIF

name = tauCorrStr.name
taucorrStr.name = 'r' + name

IF( TOTAL(STRMATCH(tagNames, 'MESSAGE', /FOLD_CASE)) EQ 1 ) THEN BEGIN
		tauCorrStr.message = tauCorrStr.message + " " + message
	ENDIF ELSE BEGIN
		tauCorrStr = CREATE_STRUCT(tauCorrStr, "message", message)
	ENDELSE
RETURN
END
