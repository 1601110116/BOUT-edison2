FUNCTION GKVs1D::DelSq
; Return (as a GKVsd object of same dimensionality)
; the second derivative of the input object of input object. 
;
; Written by W.M. Nevins
;	6/23/00
;
result = self -> MakeCopy(/noValues)
Grid = self.Grid1
values = *self.values
irange = Grid.irange
imin = irange[0] + 1
imax = irange[1] - 1
dValues = SHIFT(values, 1) - values
xold = *Grid.values
dx = SHIFT(xold, 1) - xold
xnew = 0.25*(SHIFT(xold, -1) + 2.*xold + SHIFT(xold, 1))
derivative = dValues/dx
ddx = 0.5*(dx + SHIFT(dx,-1))
secondDerivative = (derivative - SHIFT(derivative, -1))/ddx
secondDerivative = secondDerivative[imin:imax]
xnew = xnew[imin:imax]
;
; 
result.values = PTR_NEW(secondDerivative)
result.title = "!9d!X!E2!N" + self.title + '/!9d!X' + Grid.title + '!E2!N'
result.mnemonic = 'd2' + self.mnemonic + 'd' + Grid.mnemonic + 2
result.units = '('+ self.units + ')/(' + Grid.units + ')!E2!N'
vmin = GKVsd_Min(secondDerivative, Max=vmax)
result.vrange = [vmin, vmax]
Grid.irange = [0, imax-imin]
Grid.values = PTR_NEW(xnew)
commandString = 'result.Grid' + axisString + '= Grid'
ok = EXECUTE(commandString)
RETURN, result
END ; ****** GKVs1D::DbyD ****** ;


FUNCTION GKVs2D::DelSq
;
; Return (as a GKVsd object of same dimensionality)
; the laplacian of the input object taken over the first
; two independent variables.
;
; Written by W.M. Nevins
;	6/23/00
;
result = self -> MakeCopy(/noValues)

grid1=self.grid1
values = *self.values
info = SIZE(values)
nx = info[1]
ny = info[2]
irange1 = Grid1.irange
imin1 = irange1[0] + 1
imax1 = irange1[1] - 1
isize1= imax1 - imin1 + 1

grid2=self.grid2
irange2 = Grid2.irange
imin2 = irange2[0] + 1
imax2 = irange2[1] - 1
isize2= imax2 - imin2 + 1

xold = *Grid1.values
dx = SHIFT(xold, 1) - xold
xnew = 0.25*(SHIFT(xold, -1) + 2.*xold + SHIFT(xold, 1))
ddx = 0.5*(dx + SHIFT(dx,-1))
dValues = SHIFT(values, 1, 0) - values
xderivative = FLTARR(nx,ny)
xsecondDerivative = FLTARR(nx,ny)
FOR i=1, nx-2 DO BEGIN
	xderivative[i,*] = dValues[i,*]/dx[i]
	xsecondDerivative[i,*] = (xderivative[i,*] -xderivative[i-1,*])/ddx[i]
ENDFOR
xsecondDerivative = xsecondDerivative[imin1:imax1, imin2:imax2]


yold = *Grid2.values
dy = SHIFT(yold, 1) - yold
ynew = 0.25*(SHIFT(yold, -1) + 2.*yold + SHIFT(yold, 1))
ddy = 0.5*(dy + SHIFT(dy,-1))
dValues = SHIFT(values, 0, 1) - values
yderivative = FLTARR(nx,ny)
ysecondDerivative = FLTARR(nx,ny)
FOR j=1, ny-2 DO BEGIN
	yderivative[*,j] = dValues[*,j]/dy[j]
	ysecondDerivative[*,j] = (yderivative[*,j] - yderivative[*,j-1])/ddy[i]
ENDFOR
ysecondDerivative = ysecondDerivative[imin1:imax1, imin2:imax2]

Delsq = xsecondDerivative + ysecondDerivative
xnew = xnew[imin1:imax1]
grid1.values = PTR_NEW(xnew)
grid1.irange = [0, isize1-1]
grid1.range = [xnew[0], xnew[isize1-1]]

ynew = ynew[imin2:imax2]
grid2.values = PTR_NEW(ynew)
grid2.irange = [0, isize2-1]
grid2.range = [ynew[0], ynew[isize2-1]]
;
; 
result.values = PTR_NEW(DelSq)
result.title = "!9G!X!E2!N" + self.title
result.mnemonic = 'DelSa' + self.mnemonic
result.units = '('+ self.units + ')/(' + Grid1.units + ')!E2!N'
vmin = GKVsd_Min(delSq, Max=vmax)
result.vrange = [vmin, vmax]
result.Grid1=grid1
result.Grid2=grid2
RETURN, result
END ; ****** GKVs2D::DbyD ****** ;


FUNCTION GKVs3D::DelSq
;
; Return (as a GKVsd object of same dimensionality)
; the laplacian of the input object taken over the first
; two independent variables.
;
; Written by W.M. Nevins
;	6/23/00
;
result = self -> MakeCopy(/noValues)

grid1=self.grid1
values = *self.values
info = SIZE(values)
nx = info[1]
ny = info[2]
nt = info[3]
irange1 = Grid1.irange
imin1 = irange1[0] + 1
imax1 = irange1[1] - 1
isize1= imax1 - imin1 + 1

grid2=self.grid2
irange2 = Grid2.irange
imin2 = irange2[0] + 1
imax2 = irange2[1] - 1
isize2= imax2 - imin2 + 1

xold = *Grid1.values
dx = SHIFT(xold, 1) - xold
xnew = 0.25*(SHIFT(xold, -1) + 2.*xold + SHIFT(xold, 1))
ddx = 0.5*(dx + SHIFT(dx,-1))
dValues = SHIFT(values, 1, 0, 0) - values
xderivative = FLTARR(nx, ny, nt)
xsecondDerivative = FLTARR(nx, ny, nt)
FOR i=10, nx-2 DO BEGIN
	xderivative[i,*,*] = dValues[i,*,*]/dx[i]
	xsecondDerivative[i,*,*] = (xderivative[i,*,*] -xderivative[i-1,*,*])/ddx[i]
ENDFOR
xsecondDerivative = xsecondDerivative[imin1:imax1, imin2:imax2,*]


yold = *Grid2.values
dy = SHIFT(yold, 1) - yold
ynew = 0.25*(SHIFT(yold, -1) + 2.*yold + SHIFT(yold, 1))
ddy = 0.5*(dy + SHIFT(dy,-1))
dValues = SHIFT(values, 0, 1, 0) - values
yderivative = FLTARR(nx, ny, nt)
ysecondDerivative = FLTARR(nx, ny, nt)
FOR j=1, ny-2 DO BEGIN
	yderivative[*,j,*] = dValues[*,j,*]/dy[j]
	ysecondDerivative[*,j,*] = (yderivative[*,j,*] - yderivative[*,j-1,*])/ddy[j]
ENDFOR
ysecondDerivative = ysecondDerivative[imin1:imax1, imin2:imax2,*]

Delsq = xsecondDerivative + ysecondDerivative

xnew = xnew[imin1:imax1]
grid1.values = PTR_NEW(xnew)
grid1.irange = [0, isize1-1]
grid1.range = [xnew[0], xnew[isize1-1]]

ynew = ynew[imin2:imax2]
grid2.values = PTR_NEW(ynew)
grid2.irange = [0, isize2-1]
grid2.range = [ynew[0], ynew[isize2-1]]
;
; 
result.values = PTR_NEW(DelSq)
result.title = "!9G!X!E2!N" + self.title
result.mnemonic = 'DelSa' + self.mnemonic
result.units = '('+ self.units + ')/(' + Grid1.units + ')!E2!N'
vmin = GKVsd_Min(delSq, Max=vmax)
result.vrange = [vmin, vmax]
result.Grid1=grid1
result.Grid2=grid2
RETURN, result
END ; ****** GKVs3D::DbyD ****** ;

