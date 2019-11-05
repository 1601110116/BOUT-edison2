FUNCTION CORRAVG, nTauCorrStr, rBins = r_bins, RhoStar=rho_star
;
; This function accepts a normalized correlation structure
; (as produced by TauCorrs, and normalized by GYRO_GKNorm)
; and produces an array of GKVs1D objects containing the 
; correlation functions averaged over the indicated radial
; bins.
;  
;  Argumemnt
;
;	nTauCorrStr	An structure containing correlation information
;			as produced by TauCorrs and normalized by 
;			GYRO_GKNorm.  This argument is required.
;
;  Keywords 
;
;	nBins		An array of bin-boundaries.  Defaults to 
;			[0.3, 0.42, 0.58, 0.7]
;			(optional)
;
;	RhoStar		The value of a/rho (or of rho/a), which is
;			used to rescale the axis of correlation functions
;
;  Written by W.M. Nevins
; 	10/9/03
;
IF (TypeOF(nTauCorrStr) NE 8) THEN BEGIN
	MESSAGE, "CorrAvg called with illegal argument, Returning", /INFORMATIONAL
	RETURN, 0
ENDIF

rScale = 0
IF(N_ELEMENTS(rho_star) EQ 1) THEN BEGIN
	RhoStar=rho_star > 1./rho_star
	rScale = 1
ENDIF

rBINS = [0.3, 0.42, 0.58, 0.7]
IF(N_ELEMENTS(r_bins) GT 1) THEN rBins = r_bins
PRINT, "Bin Boundaries:  ", rBins

nBins = N_ELEMENTS(rBins) - 1

nTauCorrStr.taucorr -> Get, axis=1, gridValues=rPtr
rValues = *rPtr

iBins = LONARR(nBins+1)
FOR i=0,nBins DO BEGIN
	test=(rValues-rBins[i])^2
	epsilon = MIN(test, index)
	iBins[i] = index
ENDFOR
PRINT, "Bin Indices:    ", iBins

Result = {	Name		:	"CorrAvg"	,		$
		rBins		:	rBins		,		$
		iBins		:	iBins		}

tagNames = TAG_NAMES(nTauCorrStr)

IF( TOTAL(STRMATCH(tagNames, 'TauCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	aTauCorrFcns = OBJARR(nBins)
	FOR i=0, nBins-1 DO	$
		aTauCorrFcns[i] = nTaucorrStr.TauCorrArr[iBins[i]] -> ObjAvg( nTaucorrStr.TauCorrArr[ iBins[i]:iBins[i+1] ] )
	Result = CREATE_STRUCT(Result, 'TauCorrFcns', aTauCorrFcns)
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'RCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	aRCorrFcns = OBJARR(nBins)
	FOR i=0, nBins-1 DO BEGIN
		aRCorrFcns[i]   =   nTaucorrStr.RCorrArr[iBins[i]] -> ObjAvg(   nTaucorrStr.RCorrArr[ iBins[i]:iBins[i+1] ] )
		IF(rScale) THEN aRCorrFcns[i] -> ScaleAxis, 1, const= RhoStar, units='!4q!X!Ds!N'
	ENDFOR
	Result = CREATE_STRUCT(Result, 'RCorrFcns', aRCorrFcns)
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'YCorrArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	aYCorrFcns = OBJARR(nBins)
	FOR i=0, nBins-1 DO	$
		aYCorrFcns[i]   =   nTaucorrStr.YCorrArr[iBins[i]] -> ObjAvg(   nTaucorrStr.YCorrArr[ iBins[i]:iBins[i+1] ] )
	Result = CREATE_STRUCT(Result, 'YCorrFcns', aYCorrFcns)
ENDIF

IF( TOTAL(STRMATCH(tagNames, 'kSpectArr', /FOLD_CASE)) EQ 1 ) THEN BEGIN
	aKspect = OBJARR(nBins)
	FOR i=0, nBins-1 DO	$
		aKSpect[i]   =   nTaucorrStr.kSpectArr[iBins[i]] -> ObjAvg(   nTaucorrStr.kSpectArr[ iBins[i]:iBins[i+1] ] )
	Result = CREATE_STRUCT(Result, 'kspectFcns', aKspect)
ENDIF
		
RETURN, result
END   ; ****** CorrAvg ****** :
