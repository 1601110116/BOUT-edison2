FUNCTION val_interp, var,valnx,valny

;; interpolate values linely in the middle of initial grid points
;; make arrays be bouble size expanded, i.e., nxmax=nxmax*2 and nymax=nymax*2

nxmax=valnx
nymax=valny
test2_slice=fltarr(2*valnx,2*valny)
test2D=var
  for jy=0,nymax-2 do begin
        for jx=0,nxmax-2 do begin
        test2_slice[2*jx,2*jy]=test2D[jx,jy]
        test2_slice[2*jx,2*jy+1]=0.5*(test2D[jx,jy]+test2D[jx,jy+1])
        test2_slice[2*jx+1,2*jy]=0.5*(test2D[jx,jy]+test2D[jx+1,jy])
        test2_slice[2*jx+1,2*jy+1]=0.25*(test2D[jx,jy]+test2D[jx+1,jy+1]+test2D[jx+1,jy]+test2D[jx,jy+1])
        endfor
    endfor

     for jx=0,nxmax-2 do begin
       test2_slice[2*jx,2*nymax-2]=test2D[jx,nymax-1]
       test2_slice[2*jx+1,2*nymax-2]=0.5*(test2D[jx,nymax-1]+test2D[jx+1,nymax-1])
   endfor

     for jy=0,nymax-2 do begin
       test2_slice[2*nxmax-2,2*jy]=test2D[nxmax-1,jy]
       test2_slice[2*nxmax-2,2*jy+1]=0.5*(test2D[nxmax-1,jy]+test2D[nymax-1,jy+1])
   endfor

    test2_slice[2*nxmax-2,2*nymax-2]=test2D[nxmax-1,nymax-1]

     for jx=0,2*nxmax-2 do begin

      ;test2_slice[jx,2*nymax-1]=2.*test2_slice[jx,2*nymax-2]-test2_slice[jx,2*nymax-3]
      test2_slice[jx,2*nymax-1]=test2_slice[jx,2*nymax-2]
      endfor

     for jy=0,2*nymax-1 do begin

      ;test2_slice[2*nxmax-1,jy]=2.*test2_slice[2*nxmax-2,jy]-test2_slice[2*nxmax-3,jy]
      test2_slice[2*nxmax-1,jy]=test2_slice[2*nxmax-2,jy]
  endfor

  RETURN, test2_slice
END


pro lineaverz,varin,grid=g,ntmax=ntmax,R=R_lineave_position,Z=Z_lineave_position
;; calculate line average of Ni

safe_colors, /first
;g=file_import("data/circle.grd.hl2a.nc")
;g=file_import("data/circle.grd.nc")
print,'Circular Geometry g file without x point applied, NB: out mid-plane at theta=!PI/2\n'
;print,'Circular Geometry g file with x point applied, NB: out mid-plane not at theta=!PI/2\n'

path="data"
;ni=collect(path=path,var="Ni")
;ti=collect(path=path,var="Ti")
;te=collect(path=path,var="Te")
t=collect(path=path,var="t_array")
tbar=collect(path=path,var="tbar")
lbar=collect(path=path,var="Lbar")
tex=collect(path=path,var="Te_x")
vix=lbar/tbar

nxmax=g.nx
nymax=g.ny
Rxy=g.rxy
;ntmax=70
RTOL=4.e-3
Interplate=1
Interplate_twice=1
 IF NOT KEYWORD_SET(R_lineave_position) THEN R_lineave_position = g.rxy[10,32] ; default
 IF NOT KEYWORD_SET(Z_lineave_position) THEN Z_lineave_position = g.zxy[5,48] ; default


;R_lineave_position=g.rxy[10,32]
print,'R_lineave_position',R_lineave_position
;Z_lineave_position=g.zxy[5,16+32]
print,'Z_lineave_position',Z_lineave_position

;var0=reform(ni[*,*,0,*])
;var0=reform(tex*ti[*,*,0,*])
var0=reform(varin[*,*,0,*])
val_background=0.1
;t=t*tbar

sum_ni=fltarr(ntmax)
sum_ni_R=fltarr(ntmax)
val_lineave=fltarr(ntmax)
val_lineave_R=fltarr(ntmax)
ni_inline=fltarr(nxmax,nymax)
ni_inline_R=fltarr(nxmax,nymax)
t_lineave=fltarr(ntmax)
test=fltarr(nxmax,nymax,ntmax)

for nt=0,ntmax-1 do begin
    for jy=0,nymax-1 do begin
        for jx=0,nxmax-1 do begin
            test[jx,jy,nt]=1.
        endfor
    endfor
 endfor
 ;var0=test

  grxy=g.rxy
  gzxy=g.zxy
  if Interplate gt 0. then begin
    fill_grxy=val_interp(g.rxy,nxmax,nymax)
    fill_gzxy=val_interp(g.zxy,nxmax,nymax)
     if Interplate_twice gt 0. then begin
       fill_grxy_twice=val_interp(fill_grxy,2*nxmax,2*nymax)
       fill_gzxy_twice=val_interp(fill_gzxy,2*nxmax,2*nymax)
     endif
  endif

  nxmax0=nxmax
  nymax0=nymax
  for nt=0,ntmax-1 do begin

    var1=reform(var0[*,*,nt])

    ;;interpolate values
    if Interplate gt 0. then begin
       var1=val_interp(var1,nxmax0,nymax0)
       grxy=fill_grxy
       gzxy=fill_gzxy
       nxmax=2*nxmax0
       nymax=2*nymax0
       if Interplate_twice gt 0. then begin
         var1=val_interp(var1,2*nxmax0,2*nymax0)
         grxy=fill_grxy_twice
         gzxy=fill_gzxy_twice 
         nxmax=2*nxmax
         nymax=2*nymax
       endif
    endif
    ;; end of interpolate 
    num_z=0
    firstz=0.
    secondz=0.
    num_R=0
    firstR=0.
    secondR=0.
    for jy=0,nymax-1 do begin
        for jx=0,nxmax-1 do begin 
         if abs(grxy[jx,jy]-R_lineave_position) le RTOL then begin  
             if num_z eq 0 then begin
              firstz=gzxy[jx,jy]
           endif
           if num_z eq 1 then begin
              secondz=gzxy[jx,jy]
           endif
          num_z=num_z+1
         ;print,'z0',firstz,'z1',secondz
        endif

         if abs(gzxy[jx,jy]-Z_lineave_position) le RTOL then begin  
             if num_R eq 0 then begin
              firstR=grxy[jx,jy]
           endif
           if num_R eq 1 then begin
              secondR=grxy[jx,jy]
           endif
          num_R=num_R+1

         endif
        endfor
     endfor   ;; n loop

     if abs(secondz) lt 1.e-10 then begin
        print, 'WARNING: only one grid point is selected, RTOL should be increase!','   secondz',secondz
     endif

     if abs(secondR) lt 1.e-10 then begin
        print, 'WARNING: only one grid point is selected, RTOL should be increase!','   secondR',secondR
     endif

    sum_ni[nt]=0.
    n_sum=0  
    lengthz=0.
    ni_inline=0.*var1+val_background

    sum_ni_R[nt]=0.
    n_sum_R=0
    lengthR=0.
    ni_inline_R=0.*var1+val_background

    for jy=0,nymax-1 do begin
        for jx=0,nxmax-1 do begin

         if abs(grxy[jx,jy]-R_lineave_position) le RTOL then begin  

           if n_sum eq 0 then begin 
              dz0=abs(secondz-firstz)
              sum_ni[nt]=sum_ni[nt]+var1[jx,jy]*dz0
              lengthz=lengthz+dz0
           endif else begin
              dz=abs(z_previous-gzxy[jx,jy])
              sum_ni[nt]=sum_ni[nt]+var1[jx,jy]*dz
              ni_inline[jx,jy]=var1[jx,jy]
              lengthz=lengthz+dz
           endelse
           z_previous=gzxy[jx,jy]
           n_sum=n_sum+1
        endif

       if abs(gzxy[jx,jy]-Z_lineave_position) le RTOL then begin  

           if n_sum_R eq 0 then begin 
              dR0=abs(secondR-firstR)
              sum_ni_R[nt]=sum_ni_R[nt]+var1[jx,jy]*dR0
              lengthR=lengthR+dR0
           endif else begin
              dR=abs(R_previous-grxy[jx,jy])
              sum_ni_R[nt]=sum_ni_R[nt]+var1[jx,jy]*dR
              ni_inline_R[jx,jy]=var1[jx,jy]
              lengthR=lengthR+dR
           endelse
           R_previous=grxy[jx,jy]
           n_sum_R=n_sum_R+1
        endif

        endfor ;;jx loop
   endfor   ;; jy loop
   val_lineave[nt]=sum_ni[nt]/lengthz
   val_lineave_R[nt]=sum_ni_R[nt]/lengthR

   t_lineave[nt]=t[nt]
   ;print,'nt',nt,'  sum',sum_ni[nt],'  line length',lengthz,'  ave_line',val_lineave[nt]
   print,'nt',nt,'   ave_line in Z:',val_lineave_R[nt],'   ave_line in R:',val_lineave[nt]

 endfor ;nt loop

;A = 2 
;B = 4 
 
;IF (A EQ 2) AND (B EQ 3) THEN BEGIN 
;   PRINT, 'A = ', A 
  
;ENDIF ELSE BEGIN 
  ; IF A NE 2 THEN PRINT, 'A <> 2' ELSE PRINT, 'B <> 3' 
;ENDELSE 


;ni_inline=ni_inline*tex

print,'R_lineave_position',R_lineave_position
print,'Z_lineave_position',Z_lineave_position

window,0
;contour,ni_inline,grxy,gzxy,/fill,nlevel=100,chars=2,/iso,xtitle='R[m]',ytitle='Z[m]',title='const=1 test'
;contour,ni_inline,grxy,gzxy,/fill,nlevel=100,chars=2,/iso,xtitle='R[m]',ytitle='Z[m]',title='Ti on the line'
contour,ni_inline,grxy,gzxy,/fill,nlevel=100,chars=2,/iso,xtitle='R[m]',ytitle='Z[m]',title='Ni on a const R line'

window,2
;contour,ni_inline,/fill,nlevel=100,chars=2,/iso,xtitle='nx',ytitle='ny',/xst,/yst,title='const=1 test'
;contour,ni_inline,/fill,nlevel=100,chars=2,/iso,xtitle='nx',ytitle='ny',/xst,/yst,title='Ti on the line'
contour,ni_inline,/fill,nlevel=100,chars=2,/iso,xtitle='nx',ytitle='ny',/xst,/yst,title='Ni on a const R line'

;val_lineave=val_lineave*tex
window,3
;plot,t_lineave,val_lineave,chars=2,xtitle='t[s]',ytitle='[a.u.]',title='const=1 test',yrange=[0,1.2]
;plot,t_lineave,val_lineave,chars=2,xtitle='t[s]',ytitle='Ti[eV]',title='line averaged Ti'
plot,t_lineave,val_lineave,chars=2,xtitle='t[t0]',ytitle='Ni[No]',title='Ni averaged along a const R Line',yrange=[0.5,1.5]

window,10
contour,ni_inline_R,grxy,gzxy,/fill,nlevel=100,chars=2,/iso,xtitle='R[m]',ytitle='Z[m]',title='Ni on a const Z line '


window,12
contour,ni_inline_R,/fill,nlevel=100,chars=2,/iso,xtitle='nx',ytitle='ny',/xst,/yst,title='Ni on a const Z line '

window,13
plot,t_lineave,val_lineave_R,chars=2,xtitle='t[t0]',ytitle='Ni[No]',title=' Ni averaged along a const Z Line ',yrange=[0.5,1.5]

;window,4
;surface,fill_grxy,title='g.rxy interpolated',chars=3

;window,14
;surface,fill_grxy_twice,title='g.rxy interpolated twice',chars=3

;window,5
;surface,g.rxy,title='g.rxy',chars=3

;window,6
;surface,var1,title='Ni interpolated',chars=3

;window,7
;surface,var0[*,*,ntmax-1],title='Ni',chars=3


end
