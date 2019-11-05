 pro laststepdata,var,filename=filename,time=time

; NB: filename should be better the same as the var name
; i.e., if var=ni, then filename='ni'

;if ~keyword_set(g) then 


;if ~keyword_set(path) then path='data'
if ~keyword_set(filename) then begin
    filename='varname'
    print, 'filename should be better the same as the var name'
    print,'i.e., if var=ni, then filename=ni'
 endif

if ~keyword_set(time) then begin
   time = 0
   print, 'The time of the last step is needed to input via: time= ..., set to be 0 now!'
endif


lsvar = reform(var[*,*,*,time])
print,'at last time step min(var)  : ', min(lsvar)

print,'Please double check whether there is some missing data once any min() value is ZERO!'

safe_colors


;s=size(lsvar,/dim)
s=size(lsvar)

IF s[0] EQ 2 THEN BEGIN

  PRINT,"variables are 2D"
  nx=s[1]
  ny=s[2]
  nz=0

window,0
surface,lsvar,chars=3,title=filename,az=5

ENDIF 

IF s[0] EQ 3 THEN BEGIN

nx=s[1]
ny=s[2]
nz=s[3]

window,1
surface,lsvar[*,*,32],chars=3,title=filename,az=5

ENDIF


;openw,lun,'./dp0smbi'+fname_t1+'.txt',/get_lun
openw,lun,'./data/lstime_'+filename+'.txt',/get_lun


for i=0,nx-1 do begin
   for j=0,ny-1 do begin
      for k=0,nz-1 do begin
      printf,lun,lsvar[i,j,k]
      endfor
   endfor
endfor
free_lun,lun

print,filename,'file written done!'

end
