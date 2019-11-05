pro phi_fit_3
  
  filename=findfile('../data/cbm*.nc') & filename=filename[0]
  g=file_import(filename)
  x=(reform(g.psixy[*,32])-g.psi_axis)/(g.psi_bndry-g.psi_axis)
  
  restore,'nosbc_phi0.dat'
  phi0=phi0*g.bxy*1.0e4
  p032=reform(phi0[*,32])

  restore,'phi.dat'
  phi1=phi0*g.bxy*1.0e4
  p32=reform(phi1[*,32])
  
  s=size(phi1,/dim)
  dis=x[s[0]-1]-1.0
  ind=0
  for i=0,s[0]-1 do begin
     if (x[i]-1.0 gt 0.0) and (x[i]-1.0 lt dis) then begin
        dis=x[i]-1.0
        ind=i
     endif
  endfor
  print, dis, ind
  dis=p32[ind]-p32[ind-1]
  dis1=p32[ind-1]-p32[ind-2]
  p32[0:ind-1]=p32[0:ind-1]+dis-dis1
  p032=p032+dis-dis1

  win,0
  plot,x,p032,xstyle=1,thick=2
  oplot,x,p32,thick=2,color=2
  win,1
  ;plot,x,-deriv(x,p032),xstyle=1,thick=2,yrange=[min([min(-deriv(x,p032)),min(-deriv(x,p32))]),max([max(-deriv(x,p032)),max(-deriv(x,p32))])]
  plot,x,-deriv(x,p032),xstyle=1,thick=2,yrange=[-20,70],ystyle=1
  oplot,x,-deriv(x,p32),thick=2,color=2

  xcurve,x1,y1,/data
  restore,'points.dat'
  ;oplot,smooth(xa,20),smooth(ya,20),color=4
  er0=-deriv(x,p32)
  tmp=xa
  for i=0,n_elements(xa)-1 do tmp[i]=xa[n_elements(xa)-1-i]
  xa=smooth(tmp,40)
  tmp=ya
  for i=0,n_elements(ya)-1 do tmp[i]=ya[n_elements(ya)-1-i]
  ya=smooth(tmp,40)

  ind=min(where(x ge min(xa)))
  ind1=min(where(x ge max(xa)))
  er1=[er0[0:ind-1],ya,er0[ind1:515]]
  x1=[x[0:ind-1],xa,x[ind1:515]]
  help,x1,er1
  ;oplot,x1,er1,thick=2,color=3
  efit=er0
  for i=2,n_elements(x)-3 do begin
     ind=where(x1 ge x[i-2] and x1 le x[i+2])
     efit[i]=mean(er1[ind])
  endfor
  help,x,efit
  oplot,x,efit,color=2,thick=2
  save,x,efit,filename='efit.dat'
  phifit=-1.0*integral(x,efit,/accumulate)+min(p32)
  oplot,x,-deriv(x,phifit),color=4

  win,0
  plot,x,-deriv(x,p032),xstyle=1,thick=2,yrange=[-20,70]
  oplot,x,-deriv(x,p32),thick=2,color=2
  oplot,x,-deriv(x,phifit),thick=2,color=4

  win,2
  plot,x,p032,xstyle=1,thick=2,yrange=[-10,8]
  oplot,x,p32,thick=2,color=2
  oplot,x,phifit,thick=2,color=4

  phifit0=phi0
  for j=0,s[1]-1 do phifit0[*,j]=phifit
  phifit0=phifit0/(g.bxy*1.0e4)
  
  handle = file_open(filename, /write)
  s = file_write(handle, 'Phi0', phifit0)
  file_close, handle
end
