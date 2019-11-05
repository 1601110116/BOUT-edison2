FUNCTION GKV_DnoiseFcn, Dnoise
COMMON DnoiseInfo, ThisConst, coeff, NormedNuStar
COMMON moreDnoiseInfo, thisDobs, NuTurbConst
Dturb = (thisDobs - Dnoise) > 0.0
result = Dnoise - ThisConst*ALOG(1.0+coeff/(Dnoise+NormedNuStar + NuTurbConst*Dturb))
RETURN, result
END ; ****** GKV_DnoiseFcn ****** ;

FUNCTION GKVs1D::Dnoise, guessOut, guess1Out, guess2out, _Extra=extra
;
; Computes Dnoise from <(ePhi/t)^2>_noise
; using Greg Hammett's formulea
;
; Argument:
;
;	Argument is a GKVs1D object containing
;	<(ePhi/t)^2>_noise in units of [(rho/L_T)(T/e)]^2
;
; Keywords:
;
;	ChiObs		The Chi observed in the simulation
;			(as a GKVs1D object).  If ChiObs is
;			provided, then Dnoise will compute
;			the turbulence decorrelation rate, and
;			include it in the computation of 
;			Dnoise. (Optional)
;
;	a_0		Greg Hammett's a_0 parameter.
;			Defaults to 1.0
;
;	GradT		R/L_T.  Defaults to 6.9
;
;	c_0		Greg Hammett's C_0 parameter.
;			For the moment, this defaults 
;			to 1.0 (get right number from Greg!)
;
;	bmax		Maximum value of (k_perp*rho)^2 passed
;			by filters used in the field-solve.
;			Defaults to 1.0
;
;	kpMax		Maximum value of k_parallel in units of
;			1./(parallel grid spacing). 
;			Defaults to 1.0
;
;	nparallel	Number of grid points along field line.
;			Defaults to 32.
;
;	q		Safety factor.
;			Defaults to 1.4
;
;	eps		Inverse Aspect ratio at flux-surface in question.
;			Defaults to 0.18
;
;	tol		Tolerated error in root finder.
;			Defaults to 1.e-4.
;
;	Newton		Set this keyword (i.e., put "/NEWTON"
;			on command line to use vector Newton's
;			method for root solve (Default is
;			something called FX_SOLVE).
;
;	NuTurbMult	Multiplier on the turbulent decorrelation rate
;			allowing you to access the sensitivity of
;			Dnoise to the turbulent decorrelation rate.
;			Defaults to 1.0 (Optional).
;
;	
;
;  Written by W.M. Nevins
;	5/5/2005
;
COMMON DnoiseInfo, ThisConst, coeff, NormedNuStar
COMMON moreDnoiseInfo, thisDobs, nuTurbConst

a_0 = 1.0
result = GetKeyWord("a_0", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN a_0=result

GradT=6.92
result = GetKeyWord("GradT", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN GradT=result

c_0 = 0.55 ; from Greg's note of 6/16/05
result = GetKeyWord("c_0", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN c_0=result

bmax = 1.1813438
result = GetKeyWord("bmax", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN bmax=result

kpMax = 0.935668
result = GetKeyWord("kpmax", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN kpMax=result

kSqHammett = 0.349341
result = GetKeyWord("kSqHammett", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN kSqHammett=result

kSqVolume = 79.4446
result = GetKeyWord("kSqVolume", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN kSqVolume=result

VolumeHammett = 139.165
result = GetKeyWord("VolumeHammett", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN VolumeHammett=result

nparallel = 32.
result = GetKeyWord("nparallel", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN nparallel=result

q=1.4
result = GetKeyWord("q", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN q=result

sHat=0.78
result = GetKeyWord("sHat", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN shat=result

eps=0.18
result = GetKeyWord("eps", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN eps=result

tol = 1.e-4
result = GetKeyWord("tol", extra)
IF(Query_Real(result) + Query_Integer(result)) THEN tol=result

Newton = 0b
result = GetKeyWord("Newton", extra)
IF(Query_Integer(result)) THEN Newton=result
;
; compute coefficient second term in argument to Logarithm in expression for Dnoise
;
dzOverR = (2.*!PI*q/nParallel)*SQRT(1.+(eps/q)^2)
LTOverdz = 1./(GradT*dzOverR)
phiSq = *self.values
;coeff = c_0*SQRT(2./!PI)*( (2.0+bmax)/(1.0+bmax) )*(kpMax/dzOverR)/GradT/bmax
coeff = (SQRT(2.)/3.05)*kpMax*LTOverdz/kSqHammett
;
; Compute constant out front
;
;bmaxSq = bmax*bmax
;h_Hammett = 4.*( ALOG(1+bmax) - bmax/(1.+bmax) )/ALOG(1.+2.*bmax)
;
;const = a_0*(SQRT(!PI/2.)/8.)*(dzOverR)*(GradT/kpMax)*h_Hammett

const = (a_0/12.)*(VolumeHammett/kSqVolume)/coeff
;
; Compute NormedNuStar
;
NuStarBshear = (sHat/q)*(SQRT(kSqHammett)/GradT)
NuStarDrift  = SQRT(kSqHammett/2)/GradT
NuStar = NuStarBshear + NustarDrift

nt = N_ELEMENTS(phiSq)
D_obs = FLTARR(nt)
kr_0 = 0.044
ky_0=0.15
kPerp_0Sq = kr_0^2 + ky_0^2
kPerp_0 = SQRT(kPerp_0Sq)
sqrt_2 = SQRT(2.)
NuTurbConst= kr_0/SQRT( kperp_0*SQRT(kSqHammett) ) 
NuTurbMult = 1.0
result = GetKeyWord("NuTurbMult", Extra)
IF(Query_Integer(result)+Query_Real(result)) THEN NuTurbMult = result
NuTurbConst = NuTurbMult*NuTurbConst
result = GetKeyWord("ChiObs", extra)
IF(TypeOf(result) EQ 11) THEN BEGIN
	self -> GET, axis=1, GridUnits=GridUnits
	result -> Set, axis=1, GridUnits=GridUnits
	ChiObs = result -> Interpolate(self)
	ChiObs -> Get, values=ChiPtr
	D_obs = *ChiPtr/1.5
ENDIF

result = GetKeyWord("NuStar", extra)
IF( Query_real(result) + Query_Integer(result) ) THEN NuStar=result

NormedNuStar = NuStar/kSqHammett
print, NuStarBshear, NuStarDrift, NuStarBshear + NuStarDrift
;
; Guess of Dnoise (large Dnoise limit)
;
; guess1 follows on expanding Alog(1+coeff/phiSq) ~= coeff/phiSq.
; This expansion is valid in limit that phiSq >> coeff/const
;
guess1 = SQRT(const*phisq*coeff)
;
; guess2 is obtained by functional iteration.  This seems to
; work well in limit that phiSq << coeff/const
;
guess2 = const*ALOG(1.+coeff/(const*phiSq))
FOR i=1,10 DO guess2 = const*PhiSQ*ALOG(1.+coeff/guess2)
;
; Use vector logic to select proper guess for each element
; of phiSq
;
guess = guess1*(phiSQ GT coeff/const) + guess2*(phiSQ LE coeff/const)
Dnoise=guess
;
; use rootsolver to compute Dnoise
;
;IF(0b) THEN BEGIN ; NEWTON is a REALLY BAD idea!
;	Dnoise = NEWTON(guess,"GKV_DnoiseFcn", /DOUBLE)
;ENDIF ELSE BEGIN
;
; Use FX_ROOT solve to compute Dnoise
;
FOR i=0L,nt-1 DO	BEGIN
	ThisConst = phiSq[i]*const
	ThisDobs = D_obs[i]
	IF( ABS(GKV_DnoiseFcn(Dnoise[i])) GT tol) THEN	$
		Dnoise[i]=FX_ROOT([guess[i],1.e-20,100.],"GKV_DnoiseFcn", /DOUBLE, TOL=tol)
;		print, GKV_DnoiseFcn(Dnoise[i]), GKV_DnoiseFcn(guess[i])
;	ENDIF
ENDFOR
;ENDELSE

result = self -> MakeCopy(/NoValues, /NoErrorBars)

result.title="D!Dnoise!N"
result.mnemonic="Dnoise"
result.units = "(!4q!X/L!DT!N)!4q!Xv!Dt!N"
result.values = PTR_NEW(Dnoise)
vmax = MAX(Dnoise)
result.vrange=[0,vmax]

IF(N_PARAMS() GT 0) THEN BEGIN
	guessOut = result -> MakeCopy(/NoValues)
	guessOut.title="D!Dguess!N"
	guessOut.mnemonic = "Dguess"
	guessOut.values=PTR_NEW(guess)
ENDIF

IF(N_PARAMS() GT 1) THEN BEGIN
	guess1Out = result -> MakeCopy(/NoValues)
	guess1Out.title="D!Dguess1!N"
	guess1Out.mnemonic = "Dguess1"
	guess1Out.values=PTR_NEW(guess1)
ENDIF

IF(N_PARAMS() GT 2) THEN BEGIN
	guess2Out = result -> MakeCopy(/NoValues)
	guess2Out.title="D!Dguess2!N"
	guess2Out.mnemonic = "Dguess2"
	guess2Out.values=PTR_NEW(guess2)
ENDIF

RETURN, result

END ; ****** GKVs1D::Dnoise ****** ;
 			 
