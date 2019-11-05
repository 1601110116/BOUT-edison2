FUNCTION ChiValues, dt=dtIn, nt=ntIn, Power=powerIn, TauInt=tauIntIn
;
; Create a realizatin of chi values selected from
; ensemble with known statistics
;
dt=.1
IF(N_ELEMENTS(dtIn) EQ 1) THEN dt=dtIn
nt=10001
IF(N_ELEMENTS(ntIn) EQ 1) THEN nt=ntIN
tauInt = 100.
IF(N_ELEMENTS(tauIntIn) EQ 1) THEN tauInt = tauIntIn
power=0
IF(N_ELEMENTS(powerIN) EQ 1) THEN power=powerIN
IF(Power LT 0) THEN power = - power
;
; create time base
;
t = dt*FINDGEN(nt)
tMin=0
tMax=dt*(nt-1)
;
; create frequency base
;
dOmega=2.*!PI/tMax
OmegaMax = !PI/dt
Omega = -omegaMax + dOmega*FINDGEN(nt)
I_w0 = OmegaMax/dOmega
;
; Exponential spectrum
;
OmegaInt = !PI/tauInt
chiSpect = EXP(-ABS(Omega)/OmegaInt)
;
; Power law spectrum
;
chiSpect_1 = ( (Power/EXP(1.))*(OmegaInt/Omega) )^(power)
chiSpect_1[I_w0] = 0.
test = ABS(omega) GT Power*OmegaInt
IF(power NE 0) Then ChiSpect = ChiSpect*(1-test) + ChiSpect_1*test
RETURN, ChiSpect
END ; ****** ModelCHi ***** ;
