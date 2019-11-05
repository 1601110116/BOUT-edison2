;pro result 
;; analytical solution for initial f(x)=x, Boundary condition f(x=0)=0,f(x=p)=0



safe_colors, /first
path="data"
g=file_import("data/circle.grd.nc")
;g=file_import("data/uedge.grd_Up_Ni_Tei_2d.nc")
;ni=collect(path="data",var="Ni")
;nn=collect(path="data",var="NN")
t=collect(path="data",var="t_array")
Te_x=collect(path=path,var="Te_x")
Lbar=collect(path=path,var="Lbar")
tbar=collect(path=path,var="tbar")
Tn_x=0.1
Vix=Lbar/tbar
NX=g.NX
dx=findgen(NX)
Pei=findgen(NX)*0.
dN_dx=findgen(NX)*0.

jy0=g.NY/2
dx=g.DX[*,jy0]

;t=findgen(200)/200

x=g.Rxy[*,jy0]
y=g.Rxy[90,*]
t0=0
dt=2

;y=findgen(nxmax)*0

;      t1=t[dt] 

    ; IF nx  EQ nxmax-2  THEN BEGIN

     ;   print,' hahah '
      
    ; ENDIF

; print,nx, x[nx],xx,y[nx]

pei=ni*(Te+Ti)*Te_x
pn=nn*Tn_x*Te_x

;SET_PLOT, 'PS'
;DEVICE, file="NTP.ps", /color, /landscape


;!P.MULTI=[0,2,2,0,0]

!x.margin=[10,3]
!y.margin=[3,3]

window,0
plot,x,ni[*,32,0,t0],thick=2.0,title='Ni',chars=2, yrange=[0.1,200],/xst,/ylog
oplot,x,ni[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,ni[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,ni[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,ni[*,32,0,t0+4*dt],thick=2.0,color=2

window,1
plot,x,nn[*,32,0,t0],thick=2.0,title='Na',chars=2 ,yrange=[0.,4],/xst
oplot,x,nn[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,nn[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,nn[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,nn[*,32,0,t0+4*dt],thick=2.0,color=2

window,2

plot,x,Te_x*te[*,32,0,t0],thick=2.0,title='Te',chars=2,/xst
oplot,x,Te_x*te[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,Te_x*te[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,Te_x*te[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,Te_x*te[*,32,0,t0+4*dt],thick=2.0,color=2

window,3

plot,x,-Vix*vmx[*,32,0,0],thick=2.0,title='-Vmx',chars=2 ,yrange=[-600.,600],/xst
oplot,x,-Vix*vmx[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,-Vix*vmx[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,-Vix*vmx[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,-Vix*vmx[*,32,0,t0+4*dt],thick=2.0,color=2

window,4

plot,x,pei[*,32,0,t0],thick=2.0,title='P',yrange=[0,15000],chars=2,/xst
oplot,x,pei[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,pei[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,pei[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,pei[*,32,0,t0+4*dt],thick=2.0,color=2


window,5

plot,x,nm[*,32,0,t0],thick=2.0,title='Nm',yrange=[0,2],chars=2,/xst
oplot,x,nm[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,nm[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,nm[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,nm[*,32,0,t0+4*dt],thick=2.0,color=2

window,6
f=Te_x*ti
plot,x,f[*,32,0,t0],thick=2.0,title='Ti',chars=2,/xst
oplot,x,f[*,32,0,t0+dt],thick=2.0,color=4
oplot,x,f[*,32,0,t0+2*dt],thick=2.0,color=6
oplot,x,f[*,32,0,t0+3*dt],thick=2.0,color=5
oplot,x,f[*,32,0,t0+4*dt],thick=2.0,color=2

 ; DEVICE, /close
 ; SET_PLOT, 'X'


;!P.MULTI=0


;save,y,filename='.sav'
;end



    
  
