pro cal_er
  gfile=findfile('../data/cbm*.nc') & gfile=reform(gfile)
  restore,gfile+'.dat'
  ;g=file_import(gfile)
  ;phi0=collect(var='phi0',path='../nosbc_t1/data/')
  ;save,phi0,filename='phi0.dat'
  restore,'../phi0.dat'
  phi0=phi0*n.phibar
  phi1=g.phi0*n.phibar
  er=phi0
  s=size(er,/dim)
  rxy=sqrt((g.rxy-g.rmag)^2+g.zxy^2)

  for j=0,s[1]-1 do begin
     er[*,j]=deriv(g.psixy[*,j],phi1[*,j]) - deriv(g.psixy[*,j],phi0[*,j])
     ;er[*,j]=-deriv(g.psixy[*,j],phi1[*,j]-phi0[*,j])
  endfor
  save,er,filename='er.dat'
  surface,er,chars=4
  
  handle = file_open(gfile, /write)
  s = file_write(handle, 'Dphi0', er)
  file_close, handle
  
end
