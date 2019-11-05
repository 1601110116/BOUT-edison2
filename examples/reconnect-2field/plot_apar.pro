xmid=32

Apar2 = collect(var='Apar',path='data_eta-2')
dim2  = size(Apar2)
tstep2=500.
time2 = collect(var='t_array',path='data_eta-2')*tstep2
arec2 = fltarr(dim2(4))
for t=0,dim2[4]-1 do arec2[t]=max(apar2[xmid,*,*,t])
amax2=max(arec2,imax2)
tmax2=imax2*tstep2
ddtarec2 = deriv(time2,arec2)
ddtamax2=max(ddtarec2,imax2)
ddtmax2=imax2*tstep2
ddtazero2 = min(ddtarec2(imax2+1:dim2[4]-1),izero2,/abs)
ddtzero2 = izero2*tstep2

Apar3 = collect(var='Apar',path='data_eta-3')
dim3  = size(Apar3)
tstep3=500.
time3 = collect(var='t_array',path='data_eta-3')*tstep3
arec3 = fltarr(dim3(4))
for t=0,dim3[4]-1 do arec3[t]=max(apar3[xmid,*,*,t])
amax3=max(arec3,imax3)
tmax3=imax3*tstep3
ddtarec3 = deriv(time3,arec3)
ddtamax3=max(ddtarec3,imax3)
ddtmax3=imax3*tstep3
ddtazero3 = min(ddtarec3(imax3+1:dim3[4]-1),izero3,/abs)
ddtzero3 = izero3*tstep3

Apar4 = collect(var='Apar',path='data_eta-4')
dim4  = size(Apar4)
tstep4=500.
time4 = collect(var='t_array',path='data_eta-4')*tstep4
arec4 = fltarr(dim4(4))
for t=0,dim4[4]-1 do arec4[t]=max(apar4[xmid,*,*,t])
amax4=max(arec4,imax4)
tmax4=imax4*tstep4
ddtarec4 = deriv(time4,arec4)
ddtamax4 = max(ddtarec4,imax4)
ddtmax4  = imax4*tstep4
ddtazero4 = min(ddtarec4(imax4+1:dim4[4]-1),izero4,/abs)
ddtzero4 = izero4*tstep4

Apar5 = collect(var='Apar',path='data_eta-5')
dim5  = size(Apar5)
tstep5=500.
time5 = collect(var='t_array',path='data_eta-5')*tstep5
arec5 = fltarr(dim5(4))
for t=0,dim5[4]-1 do arec5[t]=max(apar5[xmid,*,*,t])
amax5=max(arec5,imax5)
tmax5=imax5*tstep5
ddtarec5 = deriv(time5,arec5)
ddtamax5 = max(ddtarec5,imax5)
ddtmax5  = imax5*tstep5
ddtazero5 = min(ddtarec5(imax5+1:dim5[4]-1),izero5,/abs)
ddtzero5 = izero5*tstep5

Apar6 = collect(var='Apar',path='data_eta-6')
dim6  = size(Apar6)
tstep6=500.
time6 = collect(var='t_array',path='data_eta-6')*tstep6
arec6 = fltarr(dim6(4))
for t=0,dim6[4]-1 do arec6[t]=max(apar6[xmid,*,*,t])
amax6=max(arec6,imax6)
tmax6=imax6*tstep6
ddtarec6 = deriv(time6,arec6)
ddtamax6 = max(ddtarec6,imax6)
ddtmax6  = imax6*tstep6
ddtazero6 = min(ddtarec6(imax6+1:dim6[4]-1),izero6,/abs)
ddtzero6 = izero6*tstep6

chars=2

;window, 1
 plot, time5, arec5, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars
oplot, time2, arec2
oplot, time3, arec3
oplot, time4, arec4
;oplot, time5, arec5


;window, 2
 plot, time2, ddtarec2, title='Reconnection rate', ytitle='dArec/dt', xtitle='time', chars=chars, /ylog, yrange=[1e-12,1e-6]
;oplot, time2, ddtarec2 
oplot, time3, ddtarec3 
oplot, time4, ddtarec4 
oplot, time5, ddtarec5


;window, 3
eta = [1e-5,1e-4,1e-3,1e-2]
tmax =[tmax5,tmax4,tmax3,tmax2]
plot, eta, tmax, /xlog, /ylog, psym=1, title='Time of peak flux', xtitle='eta', ytitle='time', chars=chars
oplot, eta, tmax
 oplot, eta, (eta/1e-5)^(-2./3.)*1e5, color=100

;window, 4
eta = [1e-5,1e-4,1e-3]
ddtmax =[ddtmax5,ddtmax4,ddtmax3]
plot, eta, ddtmax, /xlog, /ylog, psym=1, title='Time of peak reconnection rate', xtitle='dArec/dt', ytitle='time', chars=chars
 oplot, eta, ddtmax
 oplot, eta, (eta/1e-5)^(-2./3.)*1e4, color=100


;window, 5
eta = [1e-5,1e-4,1e-3,1e-2]
ddtamax =[ddtamax5,ddtamax4,ddtamax3,ddtamax2]
plot, eta, ddtamax, /xlog, /ylog, psym=1, title='Peak reconnection rate', xtitle='eta', ytitle='time', chars=chars
oplot, eta, ddtamax
;oplot, eta, eta^(1./3.)
 oplot, eta, (eta/1e-2)^(2./3.)*5e-8, color=100
