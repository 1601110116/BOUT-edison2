PRO StatArrays, ChiArr, Avg=avg, error=error
;
; compute estimated mean and error in 
; this estimate of the mean for array
; of Chi's
;
; Written by W.M. Nevins
;	5/19/04
;
n=N_ELEMENTS(ChiArr)
avg = FLTARR(n)
error = FLTARR(n)

FOR i=0,n-1 DO BEGIN
	ChiStats = ChiArr[i] -> Alt1Stats()
	avg[i] = ChiStats.avg
	error[i] = ChiStats.avgpm
	Print, avg[i], ' +/- ', error[i]
ENDFOR
RETURN
END
