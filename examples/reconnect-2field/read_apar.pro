path='etab2x-3'

;READ INPUT PARAMETERS
etab0x = collect(var='etab0x',path=path)
etab2x = collect(var='etab2x',path=path)
muix2x = collect(var='muix2x',path=path)
;muix4b = collect(var='muix4b',path=path)
;nout = collect(var='NOUT',path=path)
;timestep = collect(var='TIMESTEP',path=path)

;length_norm = collect(var='length_norm',path=path)
;speed_norm = collect(var='speed_norm',path=path)
;time_norm = collect(var='time_norm',path=path)
;timestep_norm = collect(var='timestep_norm',path=path)


Apar = collect(var='Apar',path=path)
;Jpar = collect(var='Jpar',path=path)

dim = size(Apar) & print,dim

Atot=Apar
;Apar_coil = collect(var='Apar_coil',path=path)
;for t = 0,dim[4]-1 do Atot[*,*,*,t]+=Apar_coil
Apar_ext = collect(var='Apar_ext',path=path)
for t = 0,dim[4]-1 do Atot[*,*,*,t]+=Apar_ext

;showdata, atot[*,*,0,*], chars=3


amax=fltarr(dim[4])
for t=0,dim[4]-1 do amax[t]=max(apar[32,*,*,t])
plot, amax, title='Max Apar[32,*,*,t]'

aparkz=fft(apar,dim=3)
aparxykz=aparkz

g=file_import('slab_68x32.nc')
y=total(g.dy,2,/cum)
yphase = y
for ix=0,g.nx-1 do yphase[ix,*] = -yphase[ix,*]*g.shiftangle[ix]/(2*!pi)
for kz=0,dim(3)/2-1 do aparxykz[*,*,kz]=exp(complex(0,1)*kz*yphase)*aparkz[*,*,kz]
for kz=dim(3)/2,dim(3)-1 do aparxykz[*,*,kz]=exp(complex(0,1)*(kz-dim(3))*yphase)*aparkz[*,*,kz]

amaxkz=fltarr(dim[4])
for t=0,dim[4]-1 do amaxkz[t]=max(aparkz[32,*,0,t])

