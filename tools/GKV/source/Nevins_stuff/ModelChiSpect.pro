FUNCTION GKVs1D::ModelChiSpect, omega_max=omega_in, OmegaRange=OmegaRange_in, LW=LW_in
;
;  This function acts on a (steady-state) time
;  series.  It computes the spectral density
;  using standard (Blackman Tukey) techniques.
;  It then estimates the spectral density at omega=0
;  assuming that S(omega) ~= S_0*exp(-omega*tau)
;  at small frequencies.  It returns a structure
;  containing the intercept (that is, S_0), the
;  slope (that is, tau), and a GKVs1D function 
;  containing the spectral density of "self", 
;  and a second GKVs1D function containing 
;  the resulting fit to the spectral density.
;
;  Written by W.M. Nevins
;	2/14/2007
;
avgObj = self -> avg('t')
avg = *avgObj.values

corr = self -> XCorr(/NORM)
tauCorr = corr -> FullWidth()
omega_max = !PI/tauCorr
IF(N_ELEMENTS(omega_in) EQ 1) THEN omega_max=omega_in

tInterval = self.grid1.range[1] - self.grid1.range[0]
LW = SQRT(tInterval*tauCorr)
IF(N_ELEMENTS(LW_in) EQ 1) THEN LW=LW_in
spect = self -> XSpect(LW=LW)

temp = spect -> MakeCopy()
temp -> SignalWindow, omega=[0,omega_max]
temp -> Restrict
values = ALOG(*temp.values)
omegas = *temp.grid1.values
coeffs = POLY_FIT(omegas, values, 1, vFit, vError,sigma, CorrM)
PTR_FREE, temp.values
temp.values = PTR_NEW(EXP(vFit))
IF(PTR_VALID(temp.errorbars)) THEN PTR_FREE, temp.errorbars
temp.errorbars = PTR_NEW((*temp.values)*vError/vFit)
;
; Brew our own least squares fit so that we can estimate error
; how we want to do it!
;
nPoints = N_ELEMENTS(omegas)
xAvg = TOTAL(omegas)/nPoints
yAVg = TOTAL(values)/npoints
xSqAvg = TOTAL(omegas*omegas)/nPoints
ySqAvg = TOTAL(values*values)/nPoints
xyAvg  = TOTAL(omegas*values)/nPoints
SigmaX = xSqAvg - xAvg*xAvg
sigmaXY= xyAvg - xAvg*yAvg
a = sigmaXY/sigmaX
b = yAvg - a*xAvg
PRINT, "CONST", b, coeffs[0]
PRINT, "slope", a, coeffs[1]
yFit = a*omegas + b
delta_0 = TOTAL((values - yFit)^2)/nPoints
daSq = delta_0/SigmaX
da = -SQRT(daSq)
db = -da*xAvg 
PRINT, da, db

longFit = spect -> MakeCopy(/NoValues)
omegas = *longFit.grid1.values
fitValues = exp(coeffs[0] + omegas*coeffs[1])
longFit.values=PTR_NEW(fitValues)
longFit -> Get, axis=1, range=range
longFit -> SignalWindow, omega=[0,range[1]]
longFit -> Restrict
omegas = *longFit.grid1.values
fitValues = exp(b+db + omegas*(a+da))
;
; Now make linear fit to high-frequency tail
;
temp1 = spect -> makeCopy()
OmegaRange = [.5,1.5]
IF(N_ELEMENTS(OmegaRange_in) EQ 2) THEN OmegaRange=OmegaRange_in
temp1 -> signalwindow, omega=OmegaRange
temp1 -> restrict
nValues = ALOG(*temp1.values > 1.e-12)
coeffs1 = POLY_FIT(*temp1.grid1.values, nValues, 1, vFit1, vError1) 
PTR_FREE, temp1.values
temp1.values = PTR_NEW(EXP(vFit1))
IF(PTR_VALID(temp1.errorbars)) THEN PTR_FREE, temp1.errorbars
temp1.errorbars = PTR_NEW((*temp1.values)*vError1/vFit1)
;fitValues = fitValues + exp(coeffs1[0] + omegas*coeffs1[1])
PTR_FREE, longFit.values
longFit.values=PTR_NEW(fitValues)
;
; Make estimates from our model spectral density
;
S_0    = EXP(coeffs[0])
S_0max = EXP(b+db)
tau_0   = -1.*coeffs[1]
S_1    = EXP(coeffs1[0])
tau_1  = -1.*coeffs1[1]
estVar0= 2.*(S_0/tau_0)
estVar1= 2.*(S_1/tau_1)
estVar = estVar0 + estVar1
estTauInt = !PI*(estVar0*tau_0 + estVar1*tau_1)/estVar
estError = SQRT(estVar*estTauInt/tInterval)
estError0= SQRT(2.*!PI*S_0/tInterval)
maxError = SQRT(2.*!PI*S_0max/tInterval)

estCorr = corr -> MakeCopy(/NoValues)
estCorrValues = estVar0/(1.+(*estCorr.grid1.values/tau_0)^2)
estCorrValues = estCorrValues + estVar1/(1.+(*estCorr.grid1.values/tau_1)^2)
estCorrValues = estCorrValues/estVar
estCorr.values = PTR_NEW(estCorrValues)

output = {	name	:	"ModelChiSpect",	$
		avg	:	avg,			$
		tauCorr	:	tauCorr,		$
		omega_max:	omega_max,		$
		S_0	:	S_0,			$
		S_0max	:	S_0max,			$
		tau_0	:	tau_0,			$
		tau_1	:	tau_1,			$
		estVar	:	estVar,			$
		Var	:	VARIANCE(*self.values),	$
		estTauInt:	estTauInt,		$
		estError:	estError,		$
		estError0:	estError0,		$
		maxError:	maxError,		$
		Corr	:	Corr,			$
		estCorr	:	estCorr,		$
		Spect	:	spect,			$
		shortFit:	temp,			$
		fit	:	longFit			}
RETURN, output
END ;***** FUNCTION GKVs1D::ModelChiSpect *****;	



