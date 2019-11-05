;g=file_import("data/g056129.05550_2_x260y64_psi090to106_fitTe_decresePF_gaussianY_RF_thermal.nc")
g=file_import("data/g056129.05550_2_x260y64_psi090to106_thermalSBC.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi070to106_fitTe_decreasePF_gaussianY.nc")

r1=collect(path='data',var='r1')
z1=collect(path='data',var='z1')
energy1=collect(path='data',var='energy1')
Pzeta1=collect(path='data',var='pzeta1')


g.rxy=g.rxy/g.rmag
g.zxy=g.zxy/g.rmag

nstep=100
print,energy1(1:nstep)/energy1(1)
print,pzeta1(1:nstep)/pzeta1(1)
;!p.charsize=2
;plot,g.rxy(259,*),g.zxy(259,*),psym=4
;oplot,g.rxy(0,*),g.zxy(0,*),psym=4
;oplot,g.rxy(*,0),g.zxy(*,0),psym=4
;oplot,g.rxy(*,63),g.zxy(*,63),psym=4
;oplot,g.rxy(216,8:55),g.zxy(216,8:55),psym=3
;oplot,r1(1:nstep-1),z1(1:nstep-1),thick=3


;plot,energy1(1:nstep-1),yrange=[-0.1,1.3],thick=3,charsize=2.0
;oplot,pzeta1(1:nstep-1),thick=3
