
;path='data256_Nov29C0_steadystate'

path='data'

;;; input data collect

g=file_import("data/circle.grd.hl2a.nc")

print,'collecting normalized parameters'
Ni_x=collect(path=path,var="Ni_x")
Te_x=collect(path=path,var="Te_x")
Lbar=collect(path=path,var="Lbar")
tbar=collect(path=path,var="tbar")
t=collect(path=path,var="t_array")

;print,'collecting profiles';

;ni=collect(path=path,var="Ni")
;ti=collect(path=path,var="Ti")
;te=collect(path=path,var="Te")
;jpar=collect(path=path,var="jpar_BS0")


mu = 4.e-7*3.1415926

t0=349+100
density_unit = 1.e19                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 
P0 = ni*(te+ti)                               ; normalized with Ni_x Te_x
P0 = reform(P0[*,*,0,t0])
P0 = P0*Ni_x*density_unit*ee*Te_x             ; SI unit Pascals

J0 = jpar                                     ; normalized with g.bxy/lbar/mu
J0=reform(J0[*,*,0,t0])
J0 = J0*g.bxy/lbar/mu                         ; SI unit Ampere/m^2 





;;gfile re-write with profiles TNT

;filename= "grd_smbi_elmpb.nc"
filename= "circle.grd.hl2a.nc"


handle = file_open(filename, /write)
;handle = file_open(filename, /write, /create)

s=file_write(handle,'t_P0',t[t0])
s=file_write(handle,'nt_P0',t0)
print,'Pressure profile at time ', t[t0], 'nt ', t0

s = file_write(handle, 'pressure',P0)
print,s,'pressure min',min(P0)
s = file_write(handle, 'Jpar0',J0)
print,s,'Jpar0 min',min(J0)

print,'Please double check whether there is some missing data especially when any min() value is ZERO!'

window,1
surface,ni[*,*,0,t0],az=-20,chars=3,title='Ni[No]'
window,2
surface,te_x*ti[*,*,0,t0],az=-20,chars=3,title='Ti[eV]'
window,3
surface,te_x*te[*,*,0,t0],az=-20,chars=3,title='Te[eV]'

;s = file_export(filename,g)

;print,s,'g'

 file_close, handle
