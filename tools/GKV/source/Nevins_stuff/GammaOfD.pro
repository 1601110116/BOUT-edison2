FUNCTION GKVs1D::GammaOfD, _Extra=extra
;
; Acts on D_noise(t) and returns an object array
; containing the linear growth rates, corrected
; to include the effect of D_noise, as a function
; of time
;
;	Keywords:
;
;	Gamma_0		A 1-D array containing the linear growth
;			rates in the absence of D_noise.
;			Defaults to [0.012, 0.029, 0.037, 0.031, 0.016]
;			(optional)
;
;	k_y		A 1-D array containing the corresponding values k_y. 
;			Defaults to [0.1, 0.2, 0.3, 0.4, 0.5]
;			(optional)
;
;	alpha		Parameter in computation of gamma(D):
;				gamma(D) = gamma_0 - alpha*k_y^2*D
;			Defaults to 2.4
;			(optional)
;
; Written by W.M. Nevins
;	6/30/05
;
Gamma_0 = [	0.0121,		$
		0.0295,		$
		0.0366,		$
		0.0312,		$
		0.0163		]
result = GetKeyWord("Gamma_0", extra)
IF( Query_Real(result) ) THEN Gamma_0 = result
nGammas = N_ELEMENTS(Gamma_0)

k_y = [		0.1,		$
		0.2,		$
		0.3,		$
		0.4,		$
		0.5		]
result = GetKeyWord("k_y", extra)
IF( Query_Real(result) ) THEN k_y=result
nKys = N_ELEMENTS(k_y)
IF(nKys GT nGammas) THEN k_y=k_y[0:nGammas-1]
IF(nKys LT nGammas) THEN Gamma_0 = Gamma_0[0:nKys-1]
nGammas = N_ELEMENTS(Gamma_0)
nKys = N_ELEMENTS(k_y)

alpha = 2.4
result = GetKeyWord("alpha", extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN alpha = result

D = *self.values
nt = N_ELEMENTS(D)
gammas = FLTARR(nGammas+1, nt)
result = OBJARR(nGammas+1)
 
result[0] = self -> MakeCopy(/NoValues)
result[0] -> Set, values=PTR_NEW(FLTARR(nt))
result[0] -> Set, title="!4c!X", mnemonic="gamma", units="v!Dt!N/L!DT!N"
result[0] -> Set, indices = PTR_NEW(["*"])
vMin0 = 0.
vMax0 = 0.
FOR i=1, nGammas DO BEGIN
	result[i] = result[0] -> MakeCopy(/NoValues)
	values = gamma_0[i-1] - alpha*k_y[i-1]^2*D
	indices = ["k!Dy!N=" + STRTRIM(STRING(k_y[i-1], FORMAT='(F4.2)'),2), "*"]
	vMax = MAX(values, MIN=vMin)
	vMin0 = vMin0 < vMin
	vMax0 = vMax0 > vMax
	result[i] -> Set,	values=PTR_NEW(values), 	$
				Indices=PTR_NEW(indices),	$
				vrange = [vmin, vmax]
ENDFOR
vMin0 = vMin0 > (-vMax0)
result[0] -> set, vrange=[vMin0, vMax0]
Return, result
END	
