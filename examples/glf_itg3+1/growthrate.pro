;init=10 ; for the case w/ LD, w/o TC
;nt=15 ; for the case w/ LD, w/o TC
init=50 
nt=70
dt=1
n=20
factor=1.0/0.8197 ; [Vti/a] --> [Vti/Ln]

fint=nt-1
timedim='time('+strcompress(string(dt,format='(f3.1)'),/remove_all)+'x a/Cs0)'
x1=nt*0.2
x2=nt*0.4

gfile='./data/cyclone_516x64.nc'
xind2=258
yind2=32
g=file_import(gfile)

fpath ='./data'
p=collect(var = "Phi",path=fpath,xind=xind2,tind=[0,fint])
psi=p
t_array = collect(var="t_array", path=fpath)
moment_xyzt,psi,rms=rms
rms=reform(total(rms,2))

pgam=p[0,0,0,init:fint]
rms_avergrow=deriv(alog(rms[*]))/dt
mrms=factor*mean(rms_avergrow[init:fint])

window,1, xsize=600,ysize=600,xpos=100,ypos=100
!p.charsize=1.5
plot,rms_avergrow,yr=[-1,1],$
xtitle=timedim,ytitle='growth rate (Vti/Ln)'
xyouts,x2,0.8,'mean value='+strcompress(string(mrms),/remove_all)
xyouts,x1,0.8,'n='+strcompress(string(n),/remove_all)+':'
print,'Averaged growth rate for n=',n,' is  ',mrms

end

