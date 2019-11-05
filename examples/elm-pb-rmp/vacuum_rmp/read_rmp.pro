pro read_rmp,path=path,period=period,gfile=gfile,output=output

if ~keyword_set(period) then period=3
if ~keyword_set(path) then path='data'
if ~keyword_set(gfile) then gfile=findfile('./'+path+'/'+'cbm*.nc')
gfile=reform(gfile)
g=file_import(gfile)


rmp=collect(path=path,var='rmp_psi')

surface,rmp[*,*,0],chars=3,xtitle='X',ytitle='Y',title='Rmp_Psi'
window,1
plotpolslice,rmp,g,period=period
print,max(rmp)
s=size(rmp,/dim)

openw,lun,'./'+path+'/'+output+'.txt',/get_lun

for i=0,s[0]-1 do begin
   for j=0,s[1]-1 do begin
      for k=0,s[2]-1 do begin
      printf,lun,rmp[i,j,k]
      endfor
   endfor
endfor
free_lun,lun

end
