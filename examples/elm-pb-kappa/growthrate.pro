p = collect(path='data', var='P')
moment_xyzt, p, rms=rms
plot, deriv(alog(rms(327,32,*)))
;n0 = collect(path='data',var='n0')
print, deriv(alog(rms(327,32,*))); * sqrt(n0(317,32))

