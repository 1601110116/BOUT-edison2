pro test_phi0
  normal
  restore,'nosbc_phi0.dat'
  filename=findfile('../data/cbm*.nc') & filename=filename[0]
  restore,filename+'.dat'
  g=file_import(filename)
  phi0=phi0*n.phibar
  phi1=g.phi0*n.phibar
  er=phi0
  s=size(er,/dim)
  for j=0,s[1]-1 do begin
     er[*,j]=deriv(g.psixy[*,j],phi1[*,j])-deriv(g.psixy[*,j],phi0[*,j])
  endfor
  save,er,filename='er.dat'
  plot,er[*,32]
  handle = file_open(filename, /write)
  s = file_write(handle, 'Dphi0', er)
  file_close, handle

end
