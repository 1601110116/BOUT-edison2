
g=file_import('CMod_1100223012_1150_260x64y_0.9psi_v1_diff.bout.nc')
filename= "CMod_1100223012_1150_260x64y_0.9psi_v1_diff_source.bout.nc"
mu = 4.e-7*3.1415926
Ni_x=1.
Te_x= 10.
Ti_x= 10.
density_unit = 1.e20                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 

lbar=g.rmag
unit_psi=g.rmag*g.rmag*g.bmag
jy=38
rxy=g.Rxy[*,jy]/lbar
bpxy=g.Bpxy[*,jy]/g.bmag
x=(g.psixy[*,jy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
xreal=g.psixy[*,jy]/unit_psi

NX=g.NX
NY=g.NY
jx=g.IXSEPS1

ni=g.Niexp/Ni_x
ti=g.Tiexp/Ti_x
te=g.Teexp/Te_x

dnidx=deriv(xreal,ni[*,jy])
dtidx=deriv(xreal,ti[*,jy])
dtedx=deriv(xreal,te[*,jy])

print,'Gradient', '    Core','     separatrix','       Edge'
print,'ni : ', dnidx[0],'  ---',dnidx[jx],'  ---',dnidx[NX-1]
print,'ti_no_source : ', dtidx[0],'  ---',dtidx[jx],'  ---',dtidx[NX-1]
print,'te_no_source : ', dtedx[0],'  ---',dtedx[jx],'  ---',dtedx[NX-1]

diff0 = dnidx[0]/dnidx
chi0=dtidx[0]/dtidx
che0=dtedx[0]/dtedx

V_conv = diff0/ni[*,jy]*dnidx
window,1
plot,V_conv,chars=2,title='velorcity V_conv'

term_conv_e = te[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtedx
term_conv_i = ti[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtidx
Si=fltarr(NX)
Se=fltarr(NX)
for i=0,NX-1 do Si[i] = 0.0
for i=0,NX-1 do Se[i] = 0.0

for i=1,NX-1 do Si[i] = Si[i-1]-(term_conv_i[i]+term_conv_i[i-1])*(xreal[i]-xreal[i-1])/2.
for i=1,NX-1 do Se[i] = Se[i-1]-(term_conv_e[i]+term_conv_e[i-1])*(xreal[i]-xreal[i-1])/2.

window,2
plot,Si,chars=2,title='Source'
oplot,Se,color=2

chi = dtidx[0]/dtidx + Si/dtidx
che = dtedx[0]/dtedx + Se/dtedx

diff=fltarr(NX,NY)
chi_i=fltarr(NX,NY)
chi_e=fltarr(NX,NY)

for i= 0,63 do diff[*,i]=diff0
for i= 0,63 do chi_i[*,i]=chi
for i= 0,63 do chi_e[*,i]=che

print,'coefficient', '  separatrix'
print,'ni :  ', diff0[jx]
print,'ti :  ', chi[jx]
print,'te :  ', che[jx]

;window,1
;surface,diff,az=-20,chars=3,title='$D_{\perp}$'
;window,2
;surface,chi_i,az=-20,chars=3,title='$\chi_{i\perp}$'
;window,3
;surface,chi_e,az=-20,chars=3,title='$\chi_{e\perp}$'

window,4
plot,diff0,chars=2,title='perpendicular transport coefficient D'
window,5
plot,chi0,chars=2,title='chi_i'
oplot,chi,color=2
window,6
plot,che0,chars=2,title='chi_e'
oplot,che,color=2


handle = file_open(filename, /write)
s = file_write(handle, 'diff',diff)
s = file_write(handle, 'chi_i',chi_i)
s = file_write(handle, 'chi_e',chi_e)
file_close, handle
