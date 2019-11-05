

PRO fitexp,yne_sim,yte_sim,yti_sim

;; columns means Lie
;; rows means Hang

;s=size(var)
g=file_import('circle.grd.hl2a16246_500ms.nc')
lbar=g.rmag
bbar=g.bmag

safe_colors

read_data,nne,'ne16246.txt',header,columns=16

read_data,te,'Te16246.txt',header,columns=12

unit_density = 1.e13    ;1/cm^3
minor_a = 0.42          ; m
unit_Te = 1000.         ; eV

Te_x = 10               ;eV

ee = 1.6022e-19                               ; elementary charge 
Ni_x = 1

ne_edge = 0.05      ; 1e19/m^3
te_edge = 1.         ;eV


unit_pressure = Ni_x*unit_density*1e6*ee             ; SI unit Pascals

NX = g.nx
NY = g.ny
psi = (g.psixy-g.psi_axis)/(g.psi_bndry-g.psi_axis)
psinl = g.psixy/(lbar*lbar*bbar)

psi = double(psi)
psinl=double(psinl)

yte_sim =Make_Array(NX,NY,/Double)
yti_sim =Make_Array(NX,NY,/Double)
yne_sim =Make_Array(NX,NY,/Double)

jy0 =32

nx0=40      ; core, r/a=0
nt_ne=2

nt_te=5000

order_fitne = 6
order_fitte = 5


;;;; loading input variables 

xne=reform(nne[0,nx0:*])
yne=reform(nne[nt_ne,nx0:*])/unit_density
xne = double(xne)
yne = double(yne)

s=size(xne)


xne1=Make_Array(s[1]+1,/Double)
yne1=Make_Array(s[1]+1,/Double)

FOR j=0,s[1]-1 DO BEGIN

 xne1[j]=xne[j]
 yne1[j]=yne[j]

ENDFOR

 xne1[s[1]]=psi[NX-1,jy0]
 yne1[s[1]]=ne_edge

xte = - reform(te[1:*,0])/minor_a           ; unit a
yte = reform(te[1:*,nt_te])*unit_Te         ; eV


s=size(yte)
yte[s[1]-1]=yte[s[1]-2]*0.95

xte=double(xte)
yte=double(yte)

nxfak=100
xte1=Make_Array(s[1]+nxfak,/Double)
yte1=Make_Array(s[1]+nxfak,/Double)

FOR j=0,s[1]-1 DO BEGIN

 xte1[j]=xte[j]
 yte1[j]=yte[j]

ENDFOR

xte1[s[1]+nxfak-1]=psi[NX-1,jy0]
yte1[s[1]+nxfak-1]=te_edge
Lyte = (yte1[s[1]-1]-yte1[s[1]+nxfak-1])
Lxte = -(xte1[s[1]-1]-xte1[s[1]+nxfak-1])
dydxfak = Lyte/Lxte
dxfak = Lxte/(nxfak-1)

FOR j=0,nxfak-1 DO BEGIN

  xte1[s[1]-1+j] = xte1[s[1]-1] + double(j*dxfak)
   
  yte1[s[1]+j-1] = yte1[s[1]-1] - double(j*dxfak*dydxfak)

ENDFOR


;;;; NE calculation and fitting

x=xne1
y=yne1


coef=POLY_FIT(x,y,order_fitne,y_fit)
coef= double(coef)
y_fit = double(y_fit)

yne_sim =0.
FOR j=0,order_fitne DO BEGIN

    yne_sim += coef[j]*psi^j

ENDFOR

yne_sim = double(yne_sim)

;;y_fit[i]=coef[0]+coef[1]*x[i]+coef[2]*x[i]^2 +...

dnedx=deriv(psinl[*,jy0],yne_sim[*,jy0])
dnedx= double(dnedx)
d2nedx2=deriv(psinl[*,jy0],dnedx)
d2nedx2=double(d2nedx2)

window,0
plot,x,y, thick=2, title='Ne_exp[N0]',xrange=[0,1.1]
oplot,x,y_fit,color=2
oplot,psi[*,jy0],yne_sim[*,jy0],color=5,thick=2
window,5
plot,psi[*,jy0],dnedx,thick=2,title='deriv(Ne)'
window,6
plot,psi[*,jy0],d2nedx2,color=2,title='D2DX2(Ne)'

;;;;;
;;; TE calculation and fitting

xte1=double(xte1)
yte1=double(yte1)

x=xte1
y=yte1

coef=POLY_FIT(x,y,order_fitte,y_fit)
coef = double(coef)
y_fit = double(y_fit)

;print,'Te','Coef_fit',coef

yte_sim = 0.

FOR j=0,order_fitte DO BEGIN

    yte_sim += coef[j]*psi^j

ENDFOR 


FOR i=0,NY-1 DO BEGIN

   jx_botom = where(deriv(yte_sim[*,i])eq max(deriv(yte_sim[*,i])))
   jx_botom = min(jx_botom) - NX/16
   Lyte = (yte_sim[jx_botom,i]-te_edge)
   Lxte = -(psi[jx_botom,i]-psi[NX-1,i])
   dydxfak = Lyte/Lxte
   dxfak = Lxte/(NX-1-jx_botom)

   FOR j=jx_botom,NX-1 DO BEGIN
       yte_sim[j,i] = yte_sim[jx_botom,i] - double((j-jx_botom)*dxfak*dydxfak)

   ENDFOR
ENDFOR

yte_sim = double(yte_sim)

;;y_fit[i]=coef[0]+coef[1]*x[i]+coef[2]*x[i]^2 +...
dtedx = deriv(psinl[*,jy0],yte_sim[*,jy0])
d2tedx2 = deriv(psinl[*,jy0],dtedx)

window,1
plot,x,y, thick=2, title='Te_exp[]',xrange=[0,1.1]
oplot,x,y_fit,color=2
oplot,psi[*,jy0],yte_sim[*,jy0],color=5,thick=2

window,2
plot,psi[*,jy0],dtedx,thick=2,title='deriv(Te)'
window,3
plot,psi[*,jy0],d2tedx2,thick=2,title='D2DX2(Te)',color=4


yti_sim = yte_sim

P=yne_sim*(yti_sim+yte_sim)*unit_pressure

window,4
plot,psi[*,jy0],p[*,jy0],thick=2,title='P[pa]'


END 
	
