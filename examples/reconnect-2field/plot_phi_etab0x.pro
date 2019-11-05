device, decomposed=0
loadct, 39
xmid=32

path='phi+0_etab0x-4'
Apar0 = collect(var='Apar',path=path)
dim0  = size(Apar0)
tstep0 = collect(var='TIMESTEP',path=path) 
time0 = collect(var='t_array',path=path)*tstep0
arec0 = fltarr(dim0(4))
for t=0,dim0[4]-1 do arec0[t]=max(apar0[xmid,*,*,t])
amax0=max(arec0,imax0)
tmax0=imax0*tstep0
ddtarec0 = deriv(time0,arec0)
ddtamax0 = max(ddtarec0,imax0)
ddtmax0 = imax0*tstep0
ddtazero0 = min(ddtarec0(floor(dim0[4]/2.):dim0[4]-1),izero0,/abs)
ddtzero0 = (floor(dim0[4]/2.)+izero0)*tstep0

path='phi5+0_etab0x-4'
Apar5e0 = collect(var='Apar',path=path)
dim5e0  = size(Apar5e0)
tstep5e0 = collect(var='TIMESTEP',path=path) 
time5e0 = collect(var='t_array',path=path)*tstep5e0
arec5e0 = fltarr(dim5e0(4))
for t=0,dim5e0[4]-1 do arec5e0[t]=max(apar5e0[xmid,*,*,t])
amax5e0=max(arec5e0,imax5e0)
tmax5e0=imax5e0*tstep5e0
ddtarec5e0 = deriv(time5e0,arec5e0)
ddtamax5e0 = max(ddtarec5e0,imax5e0)
ddtmax5e0 = imax5e0*tstep5e0
ddtazero5e0 = min(ddtarec5e0(floor(dim5e0[4]/2.):dim5e0[4]-1),izero5e0,/abs)
ddtzero5e0 = (floor(dim5e0[4]/2.)+izero5e0)*tstep5e0

path='phi+1_etab0x-4'
Apar1 = collect(var='Apar',path=path)
dim1  = size(Apar1)
tstep1 = collect(var='TIMESTEP',path=path) 
time1 = collect(var='t_array',path=path)*tstep1
arec1 = fltarr(dim1(4))
for t=0,dim1[4]-1 do arec1[t]=max(apar1[xmid,*,*,t])
amax1=max(arec1,imax1)
tmax1=imax1*tstep1
ddtarec1 = deriv(time1,arec1)
ddtamax1 = max(ddtarec1,imax1)
ddtmax1 = imax1*tstep1
ddtazero1 = min(ddtarec1(floor(dim1[4]/2.):dim1[4]-1),izero1,/abs)
ddtzero1 = (floor(dim1[4]/2.)+izero1)*tstep1

path='phi5+1_etab0x-4'
Apar5e1 = collect(var='Apar',path=path)
dim5e1  = size(Apar5e1)
tstep5e1 = collect(var='TIMESTEP',path=path) 
time5e1 = collect(var='t_array',path=path)*tstep5e1
arec5e1 = fltarr(dim5e1(4))
for t=0,dim5e1[4]-1 do arec5e1[t]=max(apar5e1[xmid,*,*,t])
amax5e1=max(arec5e1,imax5e1)
tmax5e1=imax5e1*tstep5e1
ddtarec5e1 = deriv(time5e1,arec5e1)
ddtamax5e1 = max(ddtarec5e1,imax5e1)
ddtmax5e1 = imax5e1*tstep5e1
ddtazero5e1 = min(ddtarec5e1(floor(dim5e1[4]/2.):dim5e1[4]-1),izero5e1,/abs)
ddtzero5e1 = (floor(dim5e1[4]/2.)+izero5e1)*tstep5e1

path='phi+2_etab0x-4'
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

path='phi5+2_etab0x-4'
Apar5e2 = collect(var='Apar',path=path)
dim5e2  = size(Apar5e2)
tstep5e2 = collect(var='TIMESTEP',path=path) 
time5e2 = collect(var='t_array',path=path)*tstep5e2
arec5e2 = fltarr(dim5e2(4))
for t=0,dim5e2[4]-1 do arec5e2[t]=max(apar5e2[xmid,*,*,t])
amax5e2=max(arec5e2,imax5e2)
tmax5e2=imax5e2*tstep5e2
ddtarec5e2 = deriv(time5e2,arec5e2)
ddtamax5e2 = max(ddtarec5e2,imax5e2)
ddtmax5e2 = imax5e2*tstep5e2
ddtazero5e2 = min(ddtarec5e2(floor(dim5e2[4]/2.):dim5e2[4]-1),izero5e2,/abs)
ddtzero5e2 = (floor(dim5e2[4]/2.)+izero5e2)*tstep5e2

path='phi+4_etab0x-4'
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

chars=2
stop
;window, 1
 plot, time0, arec0, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2

;window, 2
phi=[1,5,10,50,100,500]
arec = [ arec0[dim0[4]-1],arec5e0[dim5e0[4]-1],arec1[dim1[4]-1],arec5e1[dim5e1[4]-1],arec2[dim2[4]-1],arec5e2[dim5e2[4]-1] ]

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
