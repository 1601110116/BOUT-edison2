pro plotgc,file,npoint,n1,n2

  ;filename='6639966'
  l_b=3.497172
;filename=string(file)
  openr,lun,file,/get_lun
  rz = dblarr(9,npoint)
  readf,lun,rz,format='(2f10.5,7e15.5)'
  free_lun,lun

  ;entry_device=!d.name
  ;set_plot,'ps'

  ;device,filename=filename + '.ps'
  !p.multi =[0,2,2,0,0]
  !p.font=0

  ;plot,rz(0,index),rz(1,index),psym=3,position=[0.25,0.15,0.75,0.85],$
  ;     xtitle='R(m)',ytitle='Z(m)',xrange=[1.8,2.2],yrange=[-0.2,0.2]
  plot,rz(0,n1:n2-1)*l_b,rz(1,n1:n2-1)*l_b,psym=3,xtitle='R',ytitle='Z'
  plot,rz(6,n1:n2-1),psym=3,xtitle='steps',ytitle='mu'
  plot,rz(7,n1:n2-1),psym=3,xtitle='steps',ytitle='E'
  plot,rz(8,n1:n2-1),psym=3,xtitle='steps',ytitle='P'
  
  !p.multi =[0,1,1,0,0]



  ;device,/close_file
  ;set_plot,entry_device

end
