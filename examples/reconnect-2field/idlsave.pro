; IDL Version 8.0 (linux x86_64 m64)
; Journal File for joseph5@smaug
; Working directory: /afs/fepcluster/home/joseph5/ws/bout-nersc-grendel/examples/reconnect-2field
; Date: Wed Aug 15 14:28:14 2012
 
@plot_phi_etab0x
;Reading from phi+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+1_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+1_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+2_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+2_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi5+0_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+4_etab0x-4/BOUT.dmp.0.nc
;Reading from phi+4_etab0x-4/BOUT.dmp.0.nc
; % Prematurely closing batch file: /afs/fepcluster/home/joseph5/ws/bout-nersc-grendel/examples/reconnect-2field/plot_phi_etab0x.pro.
 plot, time0, arec0, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /xlog, xrange=[1e2,5e8]
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
 plot, time0, arec0, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
phi=[1,5,10,50,100,500]
arec = [ arec0[dim0[4]],arec5e0[dim5e0[4]],arec1[dim1[4]],arec5e1[dim5e1[4]],arec2[dim2[4]],arec5e2[dim5e2[4]] ]
; % Attempt to subscript AREC0 with <LONG     (         101)> is out of range.
retall
phi=[1,5,10,50,100,500]
arec = [ arec0[dim0[4]-1],arec5e0[dim5e0[4]-1],arec1[dim1[4]-1],arec5e1[dim5e1[4]-1],arec2[dim2[4]-1],arec5e2[dim5e2[4]-1] ]
plot, phi, arec
plot, arec2
plot, arec5e2
path='phi5+1_etab0x-4'
Apar5e1 = collect(var='Apar',path=path)
dim5e1  = size(Apar5e1)
tstep5e1 = collect(var='TIMESTEP',path=path) 
;Reading from phi5+1_etab0x-4/BOUT.dmp.0.nc
time5e1 = collect(var='t_array',path=path)*tstep5e1
;Reading from phi5+1_etab0x-4/BOUT.dmp.0.nc
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
;Reading from phi+2_etab0x-4/BOUT.dmp.0.nc
time2 = collect(var='t_array',path=path)*tstep2
;Reading from phi+2_etab0x-4/BOUT.dmp.0.nc
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
;Reading from phi5+2_etab0x-4/BOUT.dmp.0.nc
time5e2 = collect(var='t_array',path=path)*tstep5e2
;Reading from phi5+2_etab0x-4/BOUT.dmp.0.nc
arec5e2 = fltarr(dim5e2(4))
for t=0,dim5e2[4]-1 do arec5e2[t]=max(apar5e2[xmid,*,*,t])
amax5e2=max(arec5e2,imax5e2)
tmax5e2=imax5e2*tstep5e2
ddtarec5e2 = deriv(time5e2,arec5e2)
ddtamax5e2 = max(ddtarec5e2,imax5e2)
ddtmax5e2 = imax5e2*tstep5e2
ddtazero5e2 = min(ddtarec5e2(floor(dim5e2[4]/2.):dim5e2[4]-1),izero5e2,/abs)
ddtzero5e2 = (floor(dim5e2[4]/2.)+izero5e2)*tstep5e2
phi=[1,5,10,50,100,500]
arec = [ arec0[dim0[4]-1],arec5e0[dim5e0[4]-1],arec1[dim1[4]-1],arec5e1[dim5e1[4]-1],arec2[dim2[4]-1],arec5e2[dim5e2[4]-1] ]
plot, phi, arec
plot, phi, arec, /ylog
plot, phi, arec, /ylog, /clog
; % Keyword CLOG not allowed in call to: PLOT
plot, phi, arec, /ylog, /xlog
plot, phi, 1/phi, /ylog, /xlog
; % PLOT: Warning: Infinite plot range.
print, phi
;       1       5      10      50     100     500
plot, phi, 1./phi, /ylog, /xlog
plot, phi, 1./phi^1.5, /ylog, /xlog
plot, phi, arec, /ylog, /clog
; % Keyword CLOG not allowed in call to: PLOT
plot, phi, arec, /ylog, /xlog
plot, phi, arec, /ylog, /xlog, chars=3
surface, apar5e1[xmid,*,*,100]
surface, apar5e1[xmid,*,*,100], chars=2
aext  = collect(var='Apar_ext',path=path)
surface, aext[xmid,*,*,100], chars=2
; % Attempt to subscript AEXT with <INT      (     100)> is out of range.
surface, aext[xmid,*,*], chars=2
surface, aext[31,*,*], chars=2
surface, aext[32,*,*], chars=2
surface, aext[33,*,*], chars=2
surface, aext[31,*,*], chars=2
plot, total( aext[31,*,*],2), chars=2
plot, total( aext[1,*,*],2), chars=2
ak = fft(total(aext[1,*,*],2)
; % Syntax error.
ak = fft(total(aext[1,*,*],2))
print, ak
;(      0.00000,      0.00000)
;(   0.00412821,   -0.0866215)
;(      0.00000,      0.00000)
;( -4.22542e-09,  2.51430e-09)
;(      0.00000,      0.00000)
;(  2.36277e-09, -1.21099e-09)
;(      0.00000,     -0.00000)
;(  1.43010e-09, -4.43545e-09)
;(      0.00000,     -0.00000)
;(  1.43010e-09,  4.43545e-09)
;(      0.00000,      0.00000)
;(  2.36277e-09,  1.21099e-09)
;(      0.00000,     -0.00000)
;( -4.22542e-09, -2.51430e-09)
;(      0.00000,     -0.00000)
;(   0.00412821,    0.0866215)
help aext
; % Syntax error.
help, aext
print, atan2(1,0)
; % Variable is undefined: ATAN2.
print, atan(1,0)
;      1.57080
print, atan(0,1)
;      0.00000
print, atan(0+j)
; % Variable is undefined: J.
j=imaginary(1)
print, j
;      0.00000
j=imaginary(1,0)
; % IMAGINARY: Incorrect number of arguments.
j=imaginary((1,0))
; % Syntax error.
j=complex(0,1)
print, j
;(      0.00000,      1.00000)
print, atan(j)
;(      0.00000,          Inf)
; % Program caused arithmetic error: Floating divide by 0
print, atan(real_part(j),imaginary(j))
;      0.00000
print, atan(imaginary(j),real_part(j))
;      1.57080
print, atan(real_part(ak[1]),imaginary(ak[1]))
;      3.09397
print, atan(real_part(ak[1]),imaginary(ak[1]))*180/!pi
;      177.271
help, ak
print, ak[1]
;(   0.00412821,   -0.0866215)
print, atan(imaginary(ak[1]),real_part(ak(1)))*180/!pi
;     -87.2715
print, atan(imaginary(ak[1]),real_part(ak[1]))*180/!pi
;     -87.2715
print, atan(imaginary(ak[1]),real_part(ak[1]))
;     -1.52317
Aext = collect(var='Apar_ext',path=path)
Aextk = fft(total(Aext[1,*,*],2))
Aext_phase =  atan(imaginary(Aextk[1]),real_part(Aextk[1]))
print, Aext_phase
;     -1.52317
Apar0k = fft(total(Apar0[1,*,*],2))
Apar0_phase =  atan(imaginary(Apar0k[1]),real_part(Apar0k[1]))
print, Apar0_phase
;      0.00000
Apar0k = fft(total(Apar0[xmid,*,*],2))
Apar0_phase =  atan(imaginary(Apar0k[1]),real_part(Apar0k[1]))
print, Apar0_phase
;      0.00000
plot, total(Apar0[xmid,*,*],2)
surface, Apar0[xmid,*,*]
surface, Apar0[xmid,*,*,100]
print, total(Apar0[xmid,*,*,dim0[4]-1],2)
;    0.0588740
;     0.193757
;     0.299143
;     0.358985
;     0.364177
;     0.313925
;     0.215881
;    0.0849720
;   -0.0588747
;    -0.193757
;    -0.299142
;    -0.358986
;    -0.364177
;    -0.313926
;    -0.215882
;   -0.0849713
plott, total(Apar0[xmid,*,*,dim0[4]-1],2)
; % Attempt to call undefined procedure/function: 'PLOTT'.
plot, total(Apar0[xmid,*,*,dim0[4]-1],2)
Apar0k = fft(total(Apar0[xmid,*,*,dim0[4]-1],2))
Apar0_phase =  atan(imaginary(Apar0k[1]),real_part(Apar0k[1]))
print, Apar0_phase
;     -1.41052
print, Apar0_phase-Aext_phas
; % Variable is undefined: AEXT_PHAS.
print, Apar0_phase-Aext_phase
;     0.112656
Apar5e0k = fft(total(Apar5e0[xmid,*,*,dim5e0[4]-1],2))
Apar5e0_phase =  atan(imaginary(Apar5e0k[1]),real_part(Apar5e0k[1]))
print, Apar5e0_phase-Aext_phase
;     0.996152
print, (Apar5e0_phase-Aext_phase)*180/!pi
;      57.0753
Apar1k = fft(total(Apar5e0[xmid,*,*,dim5e0[4]-1],2))
Apar1_phase =  atan(imaginary(Apar1k[1]),real_part(Apar1k[1]))
print, (Apar1_phase-Aext_phase)*180/!pi
;      57.0753
Apar1k = fft(total(Apar1[xmid,*,*,dim1[4]-1],2))
Apar1_phase =  atan(imaginary(Apar1k[1]),real_part(Apar1k[1]))
print, (Apar1_phase-Aext_phase)*180/!pi
;      87.8199
Apar2k = fft(total(Apar2[xmid,*,*,dim2[4]-1],2))
Apar2_phase =  atan(imaginary(Apar2k[1]),real_part(Apar2k[1]))
print, (Apar2_phase-Aext_phase)*180/!pi
;     -58.0664
Apar2_phase =  atan(imaginary(Apar2k[1,*]),real_part(Apar2k[1,*]))
; % Attempt to subscript APAR2K with <INT      (       1)> is out of range.
Apar2k = fft(total(Apar2[xmid,*,*,*],2))
surface, real_part(apar2k)
surface, real_part(apar2k), chars=3
Apar2k = fft(total(Apar2[xmid,*,*,*],2),dim=1)
help, apar2k
Apar2k = fft(reform(total(reform(Apar2[xmid,*,*,*],2)),dim=1)
; % Syntax error.
Apar2k = fft(reform(total(reform(Apar2[xmid,*,*,*]),2)),dim=1)
help, apar2k
surface, real_part(apar2k), chars=3
surface, imaginary(apar2k), chars=3
Apar2k = fft(reform(total(Apar2[xmid,*,*,*],2)),dim=1)
help, apar2k
surface, imaginary(apar2k), chars=3
surface, real_part(apar2k), chars=3
surface, atan(imaginary(apar2k),real_part(apar2k)), chars=3
Apar2k = fft(reform(total(Apar2[xmid,*,*,*],2)),dim=1)
Apar2_phase =  atan(imaginary(Apar2k[1,*]),real_part(Apar2k[1,*]))
plot, apar2_phase
Apar5e1k = fft(reform(total(Apar5e1[xmid,*,*,*],2)),dim=1)
Apar5e1_phase =  atan(imaginary(Apar5e1k[1,*]),real_part(Apar5e1k[1,*]))
plot, apar5e1_phase
plot, apar5e1_phase-aext_phase
plot, (apar5e1_phase-aext_phase)*180/!pi
Apar1k = fft(reform(total(Apar1[xmid,*,*,*],2)),dim=1)
Apar1_phase =  atan(imaginary(Apar1k[1,*]),real_part(Apar1k[1,*]))
plot, Apar1_phase
Apar1k = fft(total(Apar1[xmid,*,*,dim1[4]-1],2))
Apar1_phase =  atan(imaginary(Apar1k[1]),real_part(Apar1k[1]))
print, apar1_phase
;   0.00957165
Apar5e0k = fft(reform(total(Apar5e0[xmid,*,*,*],2)),dim=1)
Apar5e0_phase =  atan(imaginary(Apar5e0k[1,*]),real_part(Apar5e0k[1,*]))
plot, apar5e0_phase
Apar0k = fft(reform(total(Apar0[xmid,*,*,*],2)),dim=1)
Apar0_phase =  atan(imaginary(Apar0k[1,*]),real_part(Apar0k[1,*]))
plot, apar0_phase
Apar5e2k = fft(reform(total(Apar5e2[xmid,*,*,*],2)),dim=1)
Apar5e2_abs = abs(Apar5e2k[1,*])
Apar5e2_phase =  atan(imaginary(Apar5e2k[1,*]),real_part(Apar5e2k[1,*]))
plot, apar5e2_abs
plot, apar5e2_phase
Apar0k = fft(fft(Apar0,dim=2),dim=3)
surface, abs(apar0k[xmid,*,*,100])
surface, abs(apar0k[xmid,*,*,100]), chars=3
surface, abs(apar0k[xmid,0:10,0:5,100]), chars=3
surface, abs(apar0k[xmid,*,1,*]), chars=3
; % SURFACE: Array must have 2 dimensions: <FLOAT     Array[1, 32, 1, 101]>.
surface, abs(reform(apar0k[xmid,*,1,*])), chars=3
surface, abs(reform(apar0k[xmid,*,1,*])), chars=3, /ylog
; % SURFACE: Warning: Infinite plot range.
surface, abs(reform(apar0k[xmid,*,1,*])), chars=3, /zlog
; % SURFACE: Warning: Infinite plot range.
showdata, apar0[xmid,*,*,*]
;chars=       2
showdata, apar0[*,*,4,*]
;chars=       2
showdata, apar1[*,*,4,*]
;chars=       2
showdata, apar2[*,*,4,*]
;chars=       2
showdata, apar5e2[*,*,4,*]
;chars=       2
showdata, apar5e2[*,*,0,*]
;chars=       2
showdata, apar5e1[*,*,0,*]
;chars=       2
showdata, apar2[*,*,0,*]
;chars=       2
surface, apar2[*,*,4,100]
surface, apar2[*,*,4,101]
; % Attempt to subscript APAR2 with <INT      (     101)> is out of range.
retall
surface, apar2[*,*,4,100], chars=3
surface, apar0[*,*,4,100], chars=3
surface, apar1[*,*,4,100], chars=3
surface, apar2[*,*,4,100], chars=3
surface, apar0[*,*,0,100], chars=3
surface, apar1[*,*,0,100], chars=3
surface, apar2[*,*,0,100], chars=3
surface, apar2[*,*,4,100], chars=3
surface, apar2[*,*,0,100], chars=3
 plot, time0, arec0, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
 plot, time0, arec0, yrange=[0,0.0125], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
; % PLOT: Warning: Infinite plot range.
oplot, time2, arec2
 plot, time0, arec0, yrange=[1e-6,1e-2], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
 plot, time0, arec0, yrange=[1e-6,1e-1], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
oplot, time1, arec1
oplot, time2, arec2
oplot, timewindow, 1
; % OPLOT: Expression must be an array in this context: TIMEWINDOW.
 plot, time0, arec0, yrange=[1e-6,1e-1], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
window, 1
 plot, time0, arec0, yrange=[1e-6,1e-1], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
 plot, time0, arec0, yrange=[1e-6,1e-1], title='Reconnected flux', ytitle='Arec', xtitle='time', chars=chars, /ylog
oplot, time5e0, arec5e0
oplot, time1, arec1
oplot, time5e1, arec5e1
oplot, time2, arec2
oplot, time5e2, arec5e2
device, /close & set_plot, md
;window, 2
phi=[1.,5.,10.,50.,100.,500.]
arec = [ arec0[dim0[4]-1],arec5e0[dim5e0[4]-1],arec1[dim1[4]-1],arec5e1[dim5e1[4]-1],arec2[dim2[4]-1],arec5e2[dim5e2[4]-1] ]
plot, phi, arec, /ylog, /xlog, chars=3
plot, phi, arec, /ylog, /xlog, chars=3, title='Reconnected flux', ytitle='Arec', xtitle='E/B'
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
plot, phi, arec, /ylog, /xlog, chars=3, title='Reconnected flux', ytitle='Arec', xtitle='Phi0'
device, /close & set_plot, md
$ps2pdf idl.ps plots/arec_v_phi0.pdf
surface, apar0[32,*,*,100]
surface, apar0[32,*,*,100], chars=chars
surface, apar0[*,*,0,100], chars=chars
surface, apar0[*,*,4,100], chars=chars
surface, apar0[*,*,4,100], chars=chars, title='Apar for phi=x'
surface, apar0[*,*,4,100], chars=chars, ztitle='Apar for phi=x'
surface, apar0[*,*,4,100], chars=chars, ytitle='Apar', ztitle='phi=x'
surface, apar0[*,*,4,100], chars=chars, xtitle='x',ytitle='y', ztitle='phi=x'
surface, apar0[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='phi=x'
surface, apar0[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar0[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi1e0.pdf
surface, apar1[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar1[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi1e1.pdf
surface, apar12[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
; % Variable is undefined: APAR12.
surface, apar2[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar1[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
device, /close & set_plot, md
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar2[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi1e2.pdf
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar5e0[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi5e0.pdf
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar5e1[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=50*x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi5e1.pdf
md=!d.name & set_plot,'ps' & device,/col,bits=8,/landscape,file='idl.ps'
surface, apar5e2[*,*,4,100], chars=chars, xtitle='IX',ytitle='IY', ztitle='Apar', title='Phi0=500*x'
device, /close & set_plot, md
$ps2pdf idl.ps plots/apar_phi5e2.pdf
