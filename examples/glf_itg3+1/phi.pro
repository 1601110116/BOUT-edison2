zperiod=3
fint=300
fpath='./data'
g=file_import("./data/cyclone_196x32_nl.nc")

phitilde=collect(var="Phi",path=fpath,tind=fint)
@spenv
window,/free,xsize=600,ysize=600,xpos=100,ypos=100
plotpolslice2,reform(phitilde[*,*,*,0]),g,nlev=260,period=zperiod;, charsize=2.

STOP
end

