device, decomposed=0
loadct, 39
xmid=32

;CASE etab0x-1 takes too many internal steps, have not completed reconnection cycle
path='etab0x-1'
Apar1 = collect(var='Apar',path=path)
dim1  = size(Apar1)
tstep1 = collect(var='TIMESTEP',path=path) 
time1 = collect(var='t_array',path=path)*tstep1
arec1 = fltarr(dim1(4))
for t=0,dim1[4]-1 do arec1[t]=max(apar1[xmid,*,*,t])
;amax1=max(arec1,imax1)
;tmax1=imax1*tstep1
;ddtarec1 = deriv(time1,arec1)
;ddtamax1 = max(ddtarec1,imax1)
;ddtmax1 = imax1*tstep1
;ddtazero1 = min(ddtarec1(floor(dim1[4]/2.):dim1[4]-1),izero1,/abs)
;ddtzero1 = (floor(dim1[4]/2.)+izero1)*tstep1


path='etab0x-2'
Apar2 = collect(var='Apar',path=path)
dim2  = size(Apar2)
tstep2 = collect(var='TIMESTEP',path=path) 
time2 = collect(var='t_array',path=path)*tstep2
arec2 = fltarr(dim2(4))
for t=0,dim2[4]-1 do arec2[t]=max(apar2[xmid,*,*,t])
amax2=max(arec2,imax2)
tmax2=imax2*tstep2
ddtarec2 = deriv(time2,arec2)
ddtamax2=max(ddtarec2,imax2)
ddtmax2=imax2*tstep2
ddtazero2 = min(ddtarec2(floor(dim2[4]/2.):dim2[4]-1),izero2,/abs)
ddtzero2 = (floor(dim2[4]/2.)+izero2)*tstep2

path='etab0x-3'
Apar3 = collect(var='Apar',path=path)
dim3  = size(Apar3)
tstep3 = collect(var='TIMESTEP',path=path) 
time3 = collect(var='t_array',path=path)*tstep3
arec3 = fltarr(dim3(4))
for t=0,dim3[4]-1 do arec3[t]=max(apar3[xmid,*,*,t])
amax3=max(arec3,imax3)
tmax3=imax3*tstep3
ddtarec3 = deriv(time3,arec3)
ddtamax3=max(ddtarec3,imax3)
ddtmax3=imax3*tstep3
ddtazero3 = min(ddtarec3(floor(dim3[4]/2.):dim3[4]-1),izero3,/abs)
ddtzero3 = (floor(dim3[4]/2.)+izero3)*tstep3

path='etab0x-4'
Apar4 = collect(var='Apar',path=path)
dim4  = size(Apar4)
tstep4 = collect(var='TIMESTEP',path=path) 
time4 = collect(var='t_array',path=path)*tstep4
arec4 = fltarr(dim4(4))
for t=0,dim4[4]-1 do arec4[t]=max(apar4[xmid,*,*,t])
amax4=max(arec4,imax4)
tmax4=imax4*tstep4
ddtarec4 = deriv(time4,arec4)
ddtamax4 = max(ddtarec4,imax4)
ddtmax4  = imax4*tstep4
ddtazero4 = min(ddtarec4(floor(dim4[4]/2.):dim4[4]-1),izero4,/abs)
ddtzero4 = (floor(dim4[4]/2.)+izero4)*tstep4

 
path='etab0x-5'
Apar5 = collect(var='Apar',path=path)
time5 = collect(var='t_array',path=path)
tstep5 = collect(var='TIMESTEP',path=path) 
time5 = time5*tstep5
dim5  = size(Apar5)
arec5 = fltarr(dim5(4))
for t=0,dim5[4]-1 do arec5[t]=max(apar5[xmid,*,*,t])
amax5=max(arec5,imax5)
tmax5=imax5*tstep5
ddtarec5 = deriv(time5,smooth(arec5,10,/edge))
ddtamax5 = max(ddtarec5,imax5)
ddtmax5  = imax5*tstep5
ddtazero5 = min(ddtarec5(floor(dim5[4]/2.):dim5[4]-1),izero5,/abs)
ddtzero5 = (floor(dim5[4]/2.)+izero5)*tstep5

;CASE etab0x-6 has not reached max flux yet
path = 'etab0x-6'
Apar6 = collect(var='Apar',path=path)
dim6  = size(Apar6)
print, dim6
tstep6 = collect(var='TIMESTEP',path=path) 
time6 = collect(var='t_array',path=path)*tstep6
arec6 = fltarr(dim6(4))
for t=0,dim6[4]-1 do arec6[t]=max(apar6[xmid,*,*,t])
amax6=max(arec6,imax6)
tmax6=imax6*tstep6
ddtarec6 = deriv(time6,arec6)
ddtamax6 = max(ddtarec6,imax6)
ddtmax6  = imax6*tstep6
ddtazero6 = min(ddtarec6(floor(dim6[4]/2.):dim6[4]-1),izero6,/abs)
ddtzero6 = (floor(dim6[4]/2.)+izero6)*tstep6

chars=2

;window, 1
 plot, time5, arec5, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /xlog, xrange=[1e2,5e8]
oplot, time2, arec2
oplot, time3, arec3
oplot, time4, arec4
;oplot, time5, arec5


;window, 2
 plot, time2, ddtarec2, title='Reconnection rate', ytitle='dArec/dt', xtitle='time', chars=chars, /ylog, yrange=[1e-12,1e-4], /xlog, xrange=[1e2,5e8]
;oplot, time2, ddtarec2 
oplot, time3, ddtarec3 
oplot, time4, ddtarec4 
oplot, time5, ddtarec5


;window, 3
eta = [1e-5,1e-4,1e-3,1e-2]
imax =[imax5,imax4,imax3,imax2]
tmax =[tmax5,tmax4,tmax3,tmax2]
plot, eta, tmax, /xlog, /ylog, psym=1, title='Time of peak flux', xtitle='eta', ytitle='time', chars=chars
oplot, eta, tmax
;;oplot, eta, tmax(1)*(eta/eta(1))^(-2./3.), color=100
 oplot, eta, tmax(1)*(eta/eta(1))^(-0.6), color=100


;window, 4
eta = [1e-5,1e-4,1e-3,1e-2]
ddtmax =[ddtmax5,ddtmax4,ddtmax3,ddtmax2]
plot, eta, ddtmax, /xlog, /ylog, psym=1, title='Time of peak reconnection rate', xtitle='dArec/dt', ytitle='time', chars=chars
 oplot, eta, ddtmax
 oplot, eta, ddtmax(1)*(eta/eta(1))^(-0.8), color=100


;window, 5
eta = [1e-5,1e-4,1e-3,1e-2]
ddtamax =[ddtamax5,ddtamax4,ddtamax3,ddtamax2]
plot, eta, ddtamax, /xlog, /ylog, psym=1, title='Peak reconnection rate', xtitle='eta', ytitle='max(dA/dt)', chars=chars
oplot, eta, ddtamax
;oplot, eta, eta^(1./3.)
 oplot, eta, ddtamax(1)*(eta/eta(1))^(4./3.), color=100

;window, 6
eta = [1e-5,1e-4,1e-3,1e-2]
izero =[izero5,izero4,izero3,izero2]
tzero =[ddtzero5,ddtzero4,ddtzero3,ddtzero2]
plot, eta, tzero, /xlog, /ylog, psym=1, title='Time of peak flux', xtitle='eta', ytitle='time', chars=chars
oplot, eta, tzero
 oplot, eta, tzero(1)*(eta/eta(1))^(-0.6), color=100
