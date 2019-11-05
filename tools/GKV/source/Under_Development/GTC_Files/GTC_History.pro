FUNCTION GTC_History


output={null:0}

openr, plotidh, 'history.out'

	nmode=5
        nradial=5
	nother=18
	nrun=0
	ntime=1
	nbegin=1
	naverage=0
	delt=1.0
; nplot=# of plots, nradial=# of radial bin, ntime=# of time step
; nlast=int. points of flux surface

        readf,plotidh, nrun,nother,nradial,nmode,ntime,delt

	print,'total time step=',ntime
        print,'read # of last time step'
	new=''
        read,new
	if new ne '' then ntime=0+new
;	nbegin=1
;	print,"read # of first time step"
;	new=''
;        read,new
;	if new ne '' then nbegin=0+new

	f=fltarr(4*nmode+nradial+nother,ntime)
	readf, plotidh, f
	close, plotidh

	rmajor=1.0/0.00142
	xtime=indgen(ntime)

	np=(ntime-nbegin)/100
	nm=(ntime-nbegin)/4
;	L_T for a=400
	glength=162.8 
;	L_T for a=320
	glength=130.24 
;	L_T for a=240
	glength=97.68 
;	L_T for a=160
	glength=65.12 
;	L_T for a=120
	glength=48.84
;	L_T for a=80
	glength=32.56 
	glength=1.0

hname=strarr(45)
hname=["Exit","???","orbit poloidal", "orbit flux", "weight",$
	"V_para","energy","angular momentum", "energy flux","potential",$
	"marker", "density", "radial mode","field energy","kinetic energy",$
	"entropy","flow","delta-f flow","particle flux","momentum flux",$
	"local fluxes","shear rate1","shear rate2","shear rate3",$
	"shear rate4","shear rate5","shear rate6","shear rate7",$
	"shear rate8","(m,n) mode1 ","(m,n) mode2","(m,n) mode3","(m,n) mode4",$
	"(m,n) mode5","(m,n) mode6","(m,n) mode7","(m,n) mode8","peak chi_i",$
	"nbegin,naverage","frequency","X window2","X_size","PS file"]

xmenu,hname,BASE=hbase,SPACE=5,TITLE='history',column=2,xpad=5,ypad=5
widget_control,hbase,/realize
xmanager,"history",hbase,/NO_BLOCK

end
;**************************************************************************
