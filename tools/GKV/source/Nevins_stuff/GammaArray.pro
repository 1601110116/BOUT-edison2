FUNCTION GKVs1D::GammaArray, nGammas=nGammasIn, dGamma=dGammaIn, 	$
		Gamma_0=Gamma_0In, Norm=NormIn
;
; returns an object array containing nGammas 
; expotentials with growth (decay) rates separated by dGamma
;
; Written by W.M. Nevins
;	4/18/2007
;
; Process command line
;
nGammas=10
IF(N_ELEMENTS(nGammaIn) NE 0) THEN nGammas=nGammasIn
dGamma=-0.1
IF(N_ELEMENTS(dGammaIn) NE 0) THEN dGamma=dGammaIn
Gamma_0=0.
IF(N_ELEMENTS(Gamma_0In) NE 0) THEN Gamma_0=Gamma_0In
Norm=1.0
IF(N_ELEMENTS(NormIn) NE 0) THEN Norm=NormIn
;
; Make output object array
;
output = OBJARR(nGammas+1)
FOR i=0,nGammas DO output[i] = self -> MakeCopy(/NoValues, /NoErrorBars)
FOR i=0,nGammas DO output[i] -> Set, Title="Exp(-gamma*t)", mnemonic="Exp", units=""
;
; get time base
tValues = *self.grid1.values
Gamma = Gamma_0
FOR i=0, nGammas DO BEGIN
	values = Norm*EXP(gamma*tValues)
	gamma = gamma + dGamma
	gammaText = STRCOMPRESS(STRING(gamma))
	output[i] -> Set, values=PTR_NEW(values), RunID="!4c!X="+GammaText
ENDFOR
RETURN, output
END  ;  ****** GammaArray ******  ;


