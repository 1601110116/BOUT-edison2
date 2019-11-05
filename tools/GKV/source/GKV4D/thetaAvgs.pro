FUNCTION thetaAvgs, thisSet
;
; Computes various useful things ...
;
thisSet -> CloseTheta, ky0=2.*!PI/126.281, Lx=211.225*256./255., sHat=0.786, nx=256, nky=16
ystr = thisSet -> deltaSq(axis=2)
thetaAvg = ystr.delta -> avg(axis='theta')
thetaAvg -> Get, errorbars=errorptr
PTR_FREE, errorptr
thetaStr = thetaAvg -> deltaSq(axis=2)
I_xtot =ystr.I_x -> avg('theta')
avg = ystr.avg -> avg('theta')
avg -> Get, errorBars=errorptr
PTR_FREE, errorPtr
thetaAvg_k = thetaAvg -> FFT(axis=2)
thetaAvg_k -> get, axis=2, gridValues=kyptr
kyValues = *kyptr
dky = kyValues[1]-kyValues[0]
I_xk = OBJARR(11)
I_xk[0] = thetaStr.I_x
FOR i=1,10 DO BEGIN
	temp = thetaAvg_k -> slice(axis=2, value=i*dky)
	I_xk[i] = temp -> AbsSq()
	temp -> Trash
ENDFOR



ref = thetastr.delta -> slice(axis=1, value=26)
Corr = ystr.delta -> xcorr(ref=ref)
;
cat = {	I_xtot	:	I_xtot,		$
	I_xk	:	I_xk,		$
	avg	:	avg,		$
	dvar	:	thetaStr.delta	}
;
avg = {	Corr	:	corr		}
;
; cleanup
;
gkvdelete, ystr, thetastr.I, thetaStr.Iavg, ref
gkvdelete, thetaAvg, thetaAvg_k, thetaStr.Avg
;
; make output structure
;
output = { Cat: cat, Avg: avg }
RETURN, output
END ; thetaAvgs ;
