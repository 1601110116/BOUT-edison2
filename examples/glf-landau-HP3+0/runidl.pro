;loadct,39
tek_color
 
nvar=collect(path='./data', var='nvar')
uvar=collect(path='./data', var='uvar')
tvar=collect(path='./data', var='tvar')
qvar_l=collect(path='./data', var='qvar_l')
qvar_nl=collect(path='./data', var='qvar_nl')


g=file_import('./slab.grd.nc')
rxy=g.rxy
zxy=g.zxy


sz=SIZE(nvar)

nx=sz[1]
ny=sz[2]
nz=sz[3]
nt=sz[4]

;;-pick a particular field line
ix=nx/2
iz=nz/2
it=0
itMax=nt-1


;;-vertical coordinate
zz=REFORM(g.zxy[ix,*])

;;-parallel coordinate
b0=MEDIAN(g.bxy)
bp0=MEDIAN(g.bpxy)
spar=(B0/Bp0)*zz


tek_color
!p.multi=[0,3,4,0,1]

it0=0
tit1=' (t='+STRTRIM(STRING(it0),2)+')'
plot, spar, nvar[ix,*,iz,it0], tit='n'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, uvar[ix,*,iz,it0], tit='u'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, tvar[ix,*,iz,it0], tit='T'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, 1e4*qvar_nl[ix,*,iz,it0], tit='q'+tit1, chars=2, yr=[-1,1],/yst,/xst
oplot, spar, 1e4*qvar_l[ix,*,iz,it0], col=2, lin=2


it0=itMax/2
tit1=' (t='+STRTRIM(STRING(it0),2)+')'
plot, spar, nvar[ix,*,iz,it0], tit='n'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, uvar[ix,*,iz,it0], tit='u'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, tvar[ix,*,iz,it0], tit='T'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, 1e4*qvar_nl[ix,*,iz,it0], tit='q'+tit1, chars=2, yr=[-1,1],/yst,/xst
oplot, spar, 1e4*qvar_l[ix,*,iz,it0], col=2, lin=2


it0=itMax
tit1=' (t='+STRTRIM(STRING(it0),2)+')'
plot, spar, nvar[ix,*,iz,it0], tit='n'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, uvar[ix,*,iz,it0], tit='u'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, tvar[ix,*,iz,it0], tit='T'+tit1, chars=2, yr=[-1,1],/yst,/xst
plot, spar, 1e4*qvar_nl[ix,*,iz,it0], tit='q'+tit1, chars=2, yr=[-1,1],/yst,/xst
oplot, spar, 1e4*qvar_l[ix,*,iz,it0], col=2, lin=2
!p.multi=0

;================================================================================;
