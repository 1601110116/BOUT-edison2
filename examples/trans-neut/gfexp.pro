


Ni_exp = nne   ; unit 1x10^19/m^3=N0
Ti_exp = ti    ; eV
Te_exp = te    ; eV

Ni_exp = double(Ni_exp)
Ti_exp = double(Ti_exp)
Te_exp = double(Te_exp)


Te_x= 10.

g=file_import('circle.grd.hl2a16246_500ms.nc')

mu = 4.e-7*3.1415926

t0=349+100
density_unit = 1.e19                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 
;P0 = ni*(te+ti)                               ; normalized with Ni_x Te_x
;P0 = reform(P0[*,*,0,t0])
;P0 = P0*Ni_x*density_unit*ee*Te_x             ; SI unit Pascals

;J0 = jpar                                     ; normalized with g.bxy/lbar/mu
;J0=reform(J0[*,*,0,t0])
;J0 = J0*g.bxy/lbar/mu                         ; SI unit Ampere/m^2 

lbar=g.rmag

unit_psi=g.rmag*g.rmag*g.bmag
jy=30

rxy=g.Rxy[*,jy]/lbar
bpxy=g.Bpxy[*,jy]/g.bmag

x=(g.psixy[*,jy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
xreal=g.psixy[*,jy]/unit_psi

NX=g.NX

;;gfile re-write with profiles TNT

;filename= "grd_smbi_elmpb.nc"
filename= "circle.grd.hl2a16246_500ms.nc"


handle = file_open(filename, /write)
;handle = file_open(filename, /write, /create)

s = file_write(handle, 'Ni_exp',Ni_exp)
print,s,'Ni experiment min',min(Ni_exp)

s = file_write(handle, 'Te_exp',Te_exp)
print,s,'Te experiment min',min(Te_exp)

s = file_write(handle, 'Ti_exp',Ti_exp)
print,s,'Ti experiment min',min(Ti_exp)

print,'Please double check whether there is some missing data especially when any min() value is ZERO!'


dnidx=deriv(xreal,ni_exp[*,jy])
dtidx=deriv(xreal,ti_exp[*,jy])
dtedx=deriv(xreal,te_exp[*,jy])

print,'Gradient', '    Core',    '      Edge'
print,'ni : ', dnidx[0],'---', dnidx[NX-1]

print,'ti : ', dtidx[0],'---', dtidx[NX-1]
print,'te : ', dtedx[0],'---', dtedx[NX-1]
print, '!!Normalized Te gradient at core should divided by Te_x=',Te_x
window,1
surface,ni_exp,az=-20,chars=3,title='Ni[No]'
window,2
surface,ti,az=-20,chars=3,title='Ti[eV]'
window,3
surface,te,az=-20,chars=3,title='Te[eV]'

;s = file_export(filename,g)

;print,s,'g'

 file_close, handle
