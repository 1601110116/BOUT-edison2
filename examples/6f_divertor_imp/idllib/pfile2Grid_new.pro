pro pfile2Grid_new, pfile, gridname, spline=spline ;ext=ext,

  ON_ERROR, 2

  if ( size(pfile,/type) ne 7 ) then begin
     print, "Wrong input pfile."
     return
  endif

  if ( size(gridname,/type) ne 7 ) then begin
     print, "Wrong input BOUT++ grid file."
     return
  endif

  g=file_import(gridname)
  tmp=file_import(pfile)

  print,max(g.rxy[g.nx-1,*], ypeak), ypeak
  ixsep=g.ixseps1
;  if keyword_set(ext) then begin
;     ixsep = FIX(ixsep*1.05)
;  endif

  if keyword_set(spline) then begin
     print,'Using tension spline function to smooth.'
     xpls = intarr(9)
     xpls[0] = 0    ; set default x index
     xpls[8] = g.nx-1
     xpls[1] = FIX(g.nx*0.39)
     xpls[2] = FIX(g.nx*0.59)
     xpls[3] = FIX(g.nx*0.64)
     xpls[4] = FIX(g.nx*0.70)
     xpls[5] = FIX(g.nx*0.74)
     xpls[6] = FIX(g.nx*0.82)
     xpls[7] = FIX(g.nx*0.91) 
;     xpls = [0,200,300,330,360,380,420,470,515]
     tens = 5
     nodes = n_elements(xpls)
;     help,xpls,tens,nodes
  endif

  if (ixsep lt 0. ) then begin
     ixsep = g.nx-1
  endif else begin
     if (ixsep ge g.nx) then begin
        ixsep = g.nx-1
     endif
  endelse 

  psn=(g.psixy[*,ypeak]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
  psn2=psn[0:ixsep]

;;Interpolate ne
  n_ne=81
  print,'Number of ne is: ', n_ne
  psn_ne=fltarr(n_ne)
  ne0=fltarr(n_ne)
  psn_ne=tmp.psin
  ne0=tmp._ne
  
  nefit=interpol(ne0,psn_ne,psn2,/spline)
  ;ssetn=bspline_iterfit(psn2,ne_tmp,nord=5,bkspace=10)
  ;nefit=bspline_valu(psn2,ssetn)
  nee=fltarr(g.nx)
;dnlast=nefit[ixsep]-nefit[ixsep-1]
;print,dnlast

  for i=0,ixsep do nee[i]=nefit[i]
  for i=ixsep+1,g.nx-1 do nee[i]=nee[i-1]/1.01
  ;nee=smooth(nee,10)
;print,ixsep
;print,g.nx
window,0
plot,nee
  
;;interpolate Te
  n_te=81
  print,'Number of Te is: ', n_te
  psn_te=fltarr(n_te)
  te0=fltarr(n_te)
  psn_te=tmp.psin
  te0=tmp.te/1.0e3

  tefit=interpol(te0,psn_te,psn2,/spline)
  ;ssete=bspline_iterfit(psn2,te_tmp,nord=9,bkspace=10)
  ;tefit=bspline_valu(psn2,ssete)
  te=fltarr(g.nx)

  for i=0,ixsep do te[i]=tefit[i]*1e3
  for i=ixsep+1,g.nx-1 do te[i]=te(i-1)/1.01
  ;te=smooth(te,10)
window,1
plot,te

;;Interpolate ni
 
  n_ni=81
  print,'Number of ni is: ', n_ni
  psn_ni=fltarr(n_ni)
  ni0=fltarr(n_ni)
  psn_ni=tmp.psin
  ni0=tmp.ni
  
  nifit=interpol(ni0,psn_ni,psn2,/spline)
  ;ssetn=bspline_iterfit(psn2,ni_tmp,nord=5,bkspace=10)
  ;nifit=bspline_valu(psn2,ssetn)
  ni=fltarr(g.nx)

  for i=0,ixsep do ni[i]=nifit[i]
  for i=ixsep+1,g.nx-1 do ni[i]=ni(i-1)/1.01
  ;ni=smooth(ni,10)
window,2
plot,ni

;;Interpolate Ti
  n_ti= 81
  print,'Number of Ti is: ', n_ti
  psn_ti=fltarr(n_ti)
  ti0=fltarr(n_ti)
  psn_ti=tmp.psin
  ti0=tmp.ti/1.0e3

  tifit=interpol(ti0,psn_ti,psn2,/spline)
  ;sseti=bspline_iterfit(psn2,ti_tmp,nord=9,bkspace=10)
  ;tifit=bspline_valu(psn2,sseti)
  ti=fltarr(g.nx)

  for i=0,ixsep do ti[i]=tifit[i]*1e3
  for i=ixsep+1,g.nx-1 do ti[i]=ti[i-1]/1.01
  ;ti=smooth(ti,10)
window,3
plot,ti

;;Interpolate Er  

  n_er=81
  print,'Number of Er is: ', n_er
  psn_er=fltarr(n_er)
  er0=fltarr(n_er)
  psn_er=tmp.psin
  er0=tmp.er
  er0[n_er-1]=0.

  erfit=interpol(er0,psn_er,psn2,/spline)
  ;sseti=bspline_iterfit(psn2,er_tmp,nord=9,bkspace=10)
  ;erfit=bspline_valu(psn2,sseti)
  er=fltarr(g.nx)

  for i=0,ixsep do er[i]=erfit[i]
  for i=ixsep+1,g.nx-1 do er[i]=er[i-1]


;;Interpolate omgeb  
  n_om=81
  print,'Number of Omgeb is: ', n_om
  psn_om=fltarr(n_om)
  om0=fltarr(n_om)
  psn_om=tmp.psin
  om0=tmp.er
  om0[n_om-1]=0.

  omfit=interpol(om0,psn_om,psn2,/spline)
 ; sseti=bspline_iterfit(psn2,om_tmp,nord=9,bkspace=10)
 ; omfit=bspline_valu(psn2,sseti)
  om=fltarr(g.nx)

  for i=0,ixsep do om[i]=omfit[i]
  for i=ixsep+1,g.nx-1 do om[i]=om[i-1]

  p0 = g.pressure[*,ypeak]
 plt=plot(psn,(ni*ti+nee*te)*1.6e1,'b2');

 plt=plot(psn,p0,'r2',/overplot);
 plt=plot(psn,(ni*ti+nee*te)*1.6e1/p0,'b2');
;;;;;;;; smooth data
  if keyword_set(spline) then begin
     help,xpls,tens,nodes
;     print,xpls
     window,0
     pf=smspline(p0,tens,nodes,xpls,/plot)
     window,1
     nef=smspline(nee,tens,nodes,xpls,/plot)
     window,2
     nif=smspline(ni,tens,nodes,xpls,/plot)
     window,3
     tif=smspline(ti,tens,nodes,xpls,/plot)
     window,4
     tef=smspline(te,tens,nodes,xpls,/plot)
     window,5
     erf=smspline(er,tens,nodes,xpls,/plot)
     window,6
     omf=smspline(om,tens,nodes,xpls,/plot)
     opt = get_yesno("Is this smooth OK?")
     if opt ne 1 then begin
        repeat begin
           print,"Using input x index for spline-nodes:" 
           nodes = get_integer("Please enter the total number of input x index:")
           xpls = intarr(nodes)
           xpls[0] = 0
           xpls[-1] = g.nx-1
           print,"The x index for position 1 is set to be: ",xpls[0]
           for i =1, nodes-2 do begin
              print,"Please enter x index for position",i+1," :"
              xpls[i]=get_integer()
           endfor
           print,"The x index for the last position is set to be: ",xpls[-1]
           window,0
           pf=smspline(p0,tens,nodes,xpls,/plot)
           window,1
           nef=smspline(nee,tens,nodes,xpls,/plot)
           window,2
           nif=smspline(ni,tens,nodes,xpls,/plot)
           window,3
           tif=smspline(ti,tens,nodes,xpls,/plot)
           window,4
           tef=smspline(te,tens,nodes,xpls,/plot)
           window,5
           erf=smspline(er,tens,nodes,xpls,/plot)
           window,6
           omf=smspline(om,tens,nodes,xpls,/plot)
           opt = get_yesno("Is this smooth OK?")
        endrep until opt eq 1
     endif


     p0 = pf
     nee = nef
     ni = nif
     ti = tif
     te = tef
     er = erf
     om = omf
  endif    
;;;;;;;; Save data

  save,nee,nefit,f='ne.sav'
  save,te,tefit,f='te.sav'
  save,ni,nifit,f='ni.sav'
  save,ti,tifit,f='ti.sav'
  save,er,erfit,f='er.sav'
  save,om,omfit,f='omgeb.sav'
  save,p0,f='p0.sav'
  
  save,psn_ni,psn_ne,psn_ti,psn_te,psn_er,psn_om,ni0,ne0,ti0,te0,er0,om0,f='prof_pfile.sav'
  
;;;;;;;;;; Output profiles
  
  set_plot,'ps'
  device,file='prof_fit.ps',/color,bits=8

  !p.multi=[0,3,2]
  
  plot,psn,ni,xtitle='!7W!3!iN!n',ytitle='n!ii!n (10!e20!n m!e-3!n)$',xrange=[min(psn),1.0],thick=3
  oplot,psn_ni,ni0,col=2,thick=3

  plot,psn,nee,xtitle='!7W!3!iN!n',ytitle='n!ii!n (10!e20!n m!e-3!n)$',xrange=[min(psn),1.0],thick=3
  oplot,psn_ne,ne0,col=2,thick=3

  plot,psn,ti,xtitle='!7W!3!iN!n',ytitle='T!ii!n (eV)',thick=3,xrange=[min(psn),1.0],yrange=[0,1500]
  oplot,psn_ti,ti0*1000,col=2,thick=3

  plot,psn,te,xtitle='!7W!3!iN!n',ytitle='T!ie!n (eV)',thick=3,xrange=[min(psn),1.0],yrange=[0,1500]
  oplot,psn_te,te0*1000,col=2,thick=3

  plot,psn,er,xtitle='!7W!3!iN!n',ytitle='E!ir!n (kV/m)',thick=3,xrange=[min(psn),1.0] ;,yrange=[0,500]
  oplot,psn_er,er0,col=2,thick=3

  !p.multi=[0,1,1]
  
  device,/close
  set_plot,'x'

end

;l1=plot(psn,ni,xtitle='$\Psi_N$',ytitle='$n_i (10^{20} m^{-3})$')
;l2=plot(psn0,ne0*1e-14,xrange=[0.85,1.0],/overplot,color='red')

;l3=plot(psn,ti,xtitle='$\Psi_N$',ytitle='$T_i (eV)$')
;l4=plot(psn0,ti0*1000,xrange=[0.85,1.0],/overplot,color='red',yrange=[0,500])

;l3=plot(psn,te,xtitle='$\Psi_N$',ytitle='$T_e (eV)$')
;l4=plot(psn0,te0*1000,xrange=[0.85,1.0],/overplot,color='red',yrange=[0,500])


