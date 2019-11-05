FUNCTION PG3EQ_alpha, pg3eq_hStr, pg3eq_cStr, L_over_a=Lnorm, trange=t_range
;
; 
IF(N_ELEMENTS(Lnorm) EQ 0) THEN BEGIN
	MESSAGE, "Must supply value for L_over_a", /INFORMAIONAL
	RETURN, 0
ENDIF

trange=[300.,1700.]
IF(N_ELEMENTS(t_range) EQ 2) THEN trange=t_range

chi=pg3eq_hStr.efluxi -> times(1.5/Lnorm)
chi -> set, title="!4v!X!Di!N", mnemonic="Chi_i", units="(!4q!X!Ds!N/a)!4q!X!Ds!Nc!Ds!N"
chi -> ScaleAxis, 't', const=Lnorm, units="a/c!Ds!N"

phi=pg3eq_cStr.pot_mid -> OVER(Lnorm)
phi -> set, title="!4u!X", mnemonic="Phi", units="(!4q!X!Ds!N/a)(T/e)"
phi -> ScaleAxis, 't', const=Lnorm, units="a/c!Ds!N"
phi -> RepairValues

phiStr = phi -> delta(axis='y')
temp = phiStr.delta -> times(phiStr.delta)
tempSq = temp -> avg(axis='y')
temp -> trash
dPHiSq = tempSq -> avg(axis='x')
tempSq -> trash
dPhiSq -> set,	title='!12<!4du!X!U2!N!12>!X!Dx,y!N', 	$
		mnemonic='dPhiSq', 			$
		units='(!4q!X!Ds!N/a)!U2!N(T/e)!U2!N'

alpha = chi -> over(dPHiSq)
alpha -> Set,	title="!4a", 				$
		mnemonic="alpha",			$
		units="!4q!X!Ds!Nc!Ds!N/(T/e)"	

alpha_nl = alpha -> MakeCopy()
alpha_nl -> signalWindow, t=trange
alpha_nl -> restrict

chi_nl = chi -> MakeCopy()
chi_nl -> signalwindow, t=trange
chi_nl -> restrict

dPhiSq_nl = dPHiSq -> MakeCopy()
dPhiSq_nl -> signalwindow, t=trange
dPhiSq_nl -> restrict

Chi_vs_dPhiSq = chi_nl -> Squash(dPhiSq_nl)

alphaStats = alpha_nl -> Stats()
ChiStats = chi_nl -> Stats()
dPHiSqStats = dPhiSq_nl -> Stats()

result = {	Name	:	'pg3eqAlpha',		$
		phiStr	:	phiStr,			$
		dPhiSq	:	dPhiSq,			$
		dPhiSq_t:	dPhiSqStats.avg,	$
		dPhiSq_error:	dPhiSqStats.avgpm,	$
		chi	:	chi,			$
		chi_avg	:	chiStats.avg,		$
		chi_error:	chiStats.avgpm,		$
		alpha	:	alpha,			$
		alpha_avg:	alphaStats.avg,		$
		alpha_error:	alphaStats.avgpm,	$
		Chi_vs_dPhiSq:	Chi_vs_dPhiSq,		$
		trange	:	trange,			$
		L_over_a:	Lnorm			}
RETURN, result
END  ; ******  PG3EQ_Alpha  ****** ;
