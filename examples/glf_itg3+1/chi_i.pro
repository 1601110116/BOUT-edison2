fint=300
fpath='./data'
time_conv=collect(var="time_conv",path=fpath,tind=fint)
time = collect(var="t_array", path=fpath)
time = time*time_conv[87,0]

Chi_i_conv=collect(var="Chi_i_conv",path=fpath,tind=fint)
Chi_i=collect(var="Chi_i",path=fpath,tind=[0,fint],xind=87,yind=0)
Chi_i=Chi_i*Chi_i_conv[87,0]
@spenv
window,/free,xsize=600,ysize=600,xpos=100,ypos=100
plot,time,Chi_i[0,0,*],thick=2,xstyle=1,charsize=2,xtitle='time',ytitle='Thermal diffusivity'

STOP
end

