FUNCTION GKVs1D::SubSample, n
;
; 
result = self -> makecopy()
result -> subsample, n
RETURN, result
END


PRO GKVs1D::SubSample, n
;
;
values=*self.values
info=SIZE(values)
tIndex=info[0]
nt_in=LONG(info[tIndex])
IF(N_ELEMENTS(n) EQ 0) THEN n=2
n=n > 1
n=LONG(n)
nt_out=nt_in/n
info[tIndex] = nt_out
newValues = MAKE_ARRAY(SIZE=info)

tValues=*(self.grid1.values)
tInfo=SIZE(tvalues)
tInfo[1]=nt_out
newTvalues = MAKE_ARRAY(SIZE=tinfo)

FOR i=0L, nt_out-1 DO BEGIN
	newValues[i]=Values[n*i]
	newTvalues[i] = tValues[n*i]
ENDFOR

tGrid=self.grid1
tMin = MIN(newTValues, MAX=tMax)
tGrid.irange=[0,nt_out-1]
tGrid.range=[tmin, tmax]
PTR_FREE, tGrid.values
tGrid.values = PTR_NEW(newTvalues)
self.grid1 = tGrid

PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END


PRO GKVs2D::SubSample, n
;
;
values=*self.values
info=SIZE(values)
tIndex=info[0]
nt_in=LONG(info[tIndex])
IF(N_ELEMENTS(n) EQ 0) THEN n=2
n=n > 1
n=LONG(n)
nt_out=nt_in/n
info[tIndex] = nt_out
newValues = MAKE_ARRAY(SIZE=info)

tValues=*(self.grid2.values)
tInfo=SIZE(tvalues)
tInfo[1]=nt_out
newTvalues = MAKE_ARRAY(SIZE=tinfo)

FOR i=0L, nt_out-1 DO BEGIN
	newValues[*,i]=Values[*,n*i]
	newTvalues[i] = tValues[n*i]
ENDFOR

tGrid=self.grid2
tMin = MIN(newTValues, MAX=tMax)
tGrid.irange=[0,nt_out-1]
tGrid.range=[tmin, tmax]
PTR_FREE, tGrid.values
tGrid.values = PTR_NEW(newTvalues)
self.grid2 = tGrid

PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END


PRO GKVs3D::SubSample, n
;
;
values=*self.values
info=SIZE(values)
tIndex=info[0]
nt_in=LONG(info[tIndex])
IF(N_ELEMENTS(n) EQ 0) THEN n=2
n=n > 1
n=LONG(n)
nt_out=nt_in/n
info[tIndex] = nt_out
newValues = MAKE_ARRAY(SIZE=info)

tValues=*(self.grid3.values)
tInfo=SIZE(tvalues)
tInfo[1]=nt_out
newTvalues = MAKE_ARRAY(SIZE=tinfo)

FOR i=0L, nt_out-1 DO BEGIN
	newValues[*,*,i]=Values[*,*,n*i]
	newTvalues[i] = tValues[n*i]
ENDFOR

tGrid=self.grid3
tMin = MIN(newTValues, MAX=tMax)
tGrid.irange=[0,nt_out-1]
tGrid.range=[tmin, tmax]
PTR_FREE, tGrid.values
tGrid.values = PTR_NEW(newTvalues)
self.grid3 = tGrid

PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END

PRO GKVs4D::SubSample, n
;
;
values=*self.values
info=SIZE(values)
tIndex=info[0]
nt_in = LONG(info[tIndex])
IF(N_ELEMENTS(n) EQ 0) THEN n=2
n=n>1
n=Long(n)
nt_out=nt_in/n
info[tIndex] = nt_out
newValues = MAKE_ARRAY(SIZE=info)
tValues=*(self.grid4.values)
tInfo=SIZE(tvalues)
tInfo[1]=nt_out
newTvalues= MAKE_ARRAY(SIZE=tinfo)

For i=0L, nt_out-1 DO BEGIN
    newValues[*,*,*,i]=Values[*,*,*,n*i]
    newTvalues[i] = tvalues[n*i]
ENDFOR
tGrid=self.grid4
tMin = MIN(newTValues, MAX=tMax)
tGrid.irange=[0,nt_out-1]
tGrid.range=[tmin,tmax]
PTR_FREE, tGrid.values
tGrid.values = PTR_NEW(newTvalues)
self.grid4 = tGrid

PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END
