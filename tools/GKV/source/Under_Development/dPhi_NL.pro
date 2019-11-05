

FUNCTION dPhi_NL,HStr =Hstr, Xstr=Xstr, Cstr=Cstr, Trange=Trange, Vprime=Vprime
;
; Purpose:
;
;		This proceedure provides standard analysis of
;		turbulence data contained in structures specified 
;		on input line.
;
; Written by W.M. Nevins
;	10/6/00
;
;
; Repair data files 
;
FOR i=0, N_TAGS(Hstr)-1 DO Hstr.(i) -> Repair
FOR i=0, N_TAGS(Xstr)-1 DO Xstr.(i) -> Repair
FOR i=0, N_TAGS(Cstr)-1 DO Cstr.(i) -> Repair
;
; Set Signal Window if 'Trange' has been set
;
IF(N_ELEMENTS(Trange) EQ 2) THEN BEGIN
	FOR i=0, N_TAGS(Hstr)-1 DO Hstr.(i) -> signalwindow, t=Trange
;	FOR i=0, N_TAGS(Xstr)-1 DO Xstr.(i) -> signalwindow, t=Trange
	FOR i=0, N_TAGS(Cstr)-1 DO Cstr.(i) -> signalwindow, t=Trange
;
; and then 'restrict' objects to relevant time interval
;
;	FOR i=0, N_TAGS(Hstr)-1 DO Hstr.(i) -> Restrict
;	FOR i=0, N_TAGS(Xstr)-1 DO Xstr.(i) -> Restrict
;	FOR i=0, N_TAGS(Cstr)-1 DO Cstr.(i) -> Restrict
ENDIF
;
; Compute Avg of 'efluxi'
;
efluxiObj = Hstr.efluxi -> MakeCopy()
efluxiObj -> Restrict
efluxiStats = efluxiObj -> Stats()
efluxi_avg = efluxiStats.avg
efluxi_std = efluxiStats.std
Print, '<efluxi> = ', efluxi_avg, ' ± ', efluxi_std
efluxiObj -> Trash
;
; Compute <dPhi^2>
;
Pot_mid = Cstr.pot_mid -> MakeCopy()
Pot_mid -> Restrict
;
; Unwrap coordinate system if 'Vprime' has been set
;
IF(N_ELEMENTS(Vprime) EQ 1) THEN BEGIN
	Pot_mid  -> EddyBasis, Vprime=Vprime
ENDIF
potStr = pot_mid -> delta(axis='y')
dPhi = potStr.delta
dPhiSq = dPhi -> Times(dPhi)
dPhiSq_y = dPhiSq -> Avg(axis='y')
dPhiSq_yt = dPhiSq_y -> Avg(axis='t')
dPhiSq_ytStats = dPhiSq_yt -> Stats()
dPhiSq_avg = dPhiSq_ytStats.avg
dPhiSq_std = dPhiSq_ytStats.std
Print, '<dPhiSq> = ', dPhiSq_avg, ' ± ', dPhiSq_std
Pot_mid -> Trash
;
; Compute correlation time
;
tauCorrStr = dPhi -> TauCorrs('x')
tauCorr = tauCorrStr.taucorr
;
; Compute Avg of 'tauCorr'
;
tauCorrStats = tauCorr -> Stats()
tauCorr_avg = tauCorrStats.avg
tauCorr_std = tauCorrStats.std
Print, '<tauCorr> = ', tauCorr_avg, ' ± ', tauCorr_std
;
; Compute radial correlength
;
corrs = dPhi -> xcorr()
corrs_tau = corrs -> Slice(tau='max')
corrs_tau_y = corrs_tau -> Slice(y='max')
rCorr =corrs_tau_y ->  FullWidth( error=drcorrs)
PRINT, 'Rcorr = ', rCorr, ' ± ', drcorrs
;
; Compute 'toroidal' correlation length
;
yCorrObj = tauCorrStr.yCorr
yCorrStats = yCorrObj -> Stats()
yCorrAvg = yCorrStats.avg
yCorrSTD = yCorrStats.std
PRINT, 'yCorr = ', yCorrAvg, ' ± ', yCorrSTD
;
; Compute ExB velocity
;
V_ExB_bare = potStr.avg -> dbyd('x')
temp2 = V_ExB_bare -> Avg(axis='t')
temp1= temp2 -> FILTER('x', dL = rCorr)
temp2 -> Trash
IF(N_ELEMENTS(Vprime) EQ 1) THEN BEGIN
	VprimeString = STRING(Vprime)
	temp2 = temp1 -> Execute(VprimeString + "*(*self.Grid1.values) + 1.0*")
	temp1 -> Trash
	temp1 = temp2
ENDIF
V_ExB = temp1 
;
; Now, compute difference between Vphase and V_ExB
;
dV = tauCorrStr.vPhase -> Minus(V_ExB)
dVstats =dV -> Stats()
dV_avg = dVstats.avg
dV_std = dVstats.std
Print, '<dV> = ', dV_avg, ' ± ', dV_std
;
; Put up graphic showing V_phase and V_ExB
;
tauCorrStr.vPhase -> View, V_ExB, /pretty
;
; Now Compute effective ExB shearing rate
;
Vprime_bare = V_ExB_bare -> dbyd('x')
V_ExB_bare -> Trash
temp1 = Vprime_bare -> Filter('x', dL=rCorr)
Vprime_bare -> Trash
temp2 = temp1 -> FILTER('t', dT=tauCorr_avg)
temp1 -> Trash
IF KEYWORD_SET(vPrime) THEN BEGIN
	VprimeObj = temp2 -> Plus(vPrime)
	temp2 -> Trash
ENDIF ELSE BEGIN
	VprimeObj = temp2
ENDELSE
vPrimeStats = VprimeObj -> Stats()
vPrimeAvg = vPrimeStats.avg
vPrimeRMS = vPrimeStats.rms
vPrimeSTD = vPrimeStats.std
PRINT, "<V'> = ", vPrimeAvg, " ± ", vPrimeSTD, "     <V'>_RMS = ", vPrimeRMS 
;
; Put up graphic showing Vprime in full 2-D)
;
VprimeObj -> Set, mnemonic='Vprime', title="V!S`!R!LExB!N", units="(!4q!X/L!LT!N)c!Ls!N"
VprimeObj -> View, /pretty

;
; Prepare output structure
;
output = {	efluxi_avg:efluxi_avg, 	efluxi_std:efluxi_std, 	$
		dPhiSq_avg:dPhiSq_avg, 	dPhiSq_std:dPhiSq_std	,	$
		tauCorr_avg:tauCorr_avg,	tauCorr_std:tauCorr_std,	$
		rCorr:rCorr,			drcorrs:drcorrs,		$
		yCorrLength:yCorrAvg,		dyCorrLength:yCorrSTD,	$
		dV_avg:dV_avg,			dV_std:dV_std,			$
		dV:dV,				Vphase:tauCorrStr.vPhase,	$
		V_ExB:V_ExB,			VprimeAVG:vPrimeAvg,		$
		VprimeRMS:vPrimeRMS,		VprimeSTD:vPrimeSTD,		$
		Vprime:VprimeObj,		yCorr:yCorrObj			}
;
; clean up misc. objects...
;
dPhi -> Trash
dPhiSq -> Trash
dPhiSq_y -> Trash
dPhiSQ_yt -> Trash
potStr.avg -> Trash
corrs -> Trash
corrs_tau -> Trash
corrs_tau_y -> Trash
tauCorrStr.taucorr -> trash
; tauCorrStr.vPhase -> Trash		; leave these arround so that so that them   
; V_ExB -> Trash				; (through 'view' widget) if desired.
; dV -> Trash
;
; and we're done...
;
RETURN, output
END ; ****** dPhi_N ****** ; 
