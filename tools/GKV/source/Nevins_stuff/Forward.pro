FUNCTION GKVs2D::Forward
;
; Purpose:
;
;	To separate a GKVs2D signal into
;	a piece that is propagating "forward"
;	(that is, in the positive sense vs. axis 1)
;	and a piece that is propagating "backward"
;	(that is, in a negative sense vs. axis 1).
;	
; Written by W.M. Nevins
;	8/14/08
;
; Fourier transform in space and time
;
temp1 = self -> FFT(1)
temp = temp1 -> FFT(2)
temp1 -> TRASH
;
; Find index in k (axis 1) where k=0
;
kValues = *temp.grid1.values
nK = N_ELEMENTS(kValues)
kSq = kValues*kValues
eps = MIN(kSq, k0)
dk = kValues[k0+1]-kValues[k0]
zerok = (eps LT dk*dk/10.)
;
; and find the index in omega (axis2) where omega=0
;
wValues = *temp.grid2.values
nW = N_ELEMENTS(wValues)
wSq = wValues*wValues
eps = MIN(wSq, w0)
dw = wValues[w0+1]-wValues[w0]
zerow = (eps LT dw*dw/10.)
;
; Consruct forward and backward filters
;
One    = MAKE_ARRAY(k0+1,w0+1, VALUE=1.)
altOne = MAKE_ARRAY(nk-k0-1, nw-w0-1, VALUE=1.)
BackwardFilter  = FLTARR(nk,nw)
BackwardFilter[0:k0,0:w0] = ONE
BackwardFilter[k0+1:nK-1,w0+1:nw-1] = altOne
IF( zeroK) THEN BackwardFilter[k0, *] = 0.5
IF(zeroW ) THEN BAckwardFilter[*, w0] = 0.5
One = MAKE_ARRAY(nk, nw, VALUE=1, /INT)
ForwardFilter = One - BackwardFilter
;
; Apply Filters to data
;
fdata = temp -> Times(ForwardFilter)
bdata = temp -> Times(backwardFilter)
temp -> TRASH
;
; Invert transforms
;
temp1   = fData -> FFT(1, /Inverse)
forward = temp1 -> FFT(2, /Inverse)
temp1 -> trash
temp1    = bdata -> FFT(1, /Inverse)
backward = temp1 -> FFT(2, /Inverse)
;
; fuss with units, etc.
;
self -> get, title=title, units=units
self -> get, axis=1, gridunits=grid1units
self -> get, axis=2, gridunits=grid2units
forward  -> set, title=title + "!D+!N"
forward  -> set, axis=1, gridunits=grid1units
forward  -> set, axis=2, gridunits=grid2units
backward -> set, title=title + "!D-!N"
backward -> set, axis=1, gridunits=grid1units
backward -> set, axis=2, gridunits=grid2units
;
; Create output structure
;
output = 	{	Name	:	"Forward",	$
			Forward	:	Forward,	$
			Backward:	Backward	}
RETURN, output
END ; ****** GKVs2D::Forward ****** ;	

FUNCTION GKVs3D::Forward
;
; Purpose:
;
;	To separate a GKVs2D signal into
;	a piece that is propagating "forward"
;	(that is, in the positive sense vs. axis 2)
;	and a piece that is propagating "backward"
;	(that is, in a negative sense vs. axis 2).
;	
; Written by W.M. Nevins
;	8/14/08
;
; Fourier transform in space and time
;
temp1 = self -> FFT(2)
temp = temp1 -> FFT(3)
temp1 -> TRASH
;
; Find index in k (axis 2) where k=0
;
kValues = *temp.grid2.values
nK = N_ELEMENTS(kValues)
kSq = kValues*kValues
eps = MIN(kSq, k0)
dk = kValues[k0+1]-kValues[k0]
zerok = (eps LT dk*dk/10.)
;
; and find the index in omega (axis3) where omega=0
;
wValues = *temp.grid3.values
nW = N_ELEMENTS(wValues)
wSq = wValues*wValues
eps = MIN(wSq, w0)
dw = wValues[w0+1]-wValues[w0]
zerow = (eps LT dw*dw/10.)
;
; Consruct forward and backward filters
;
One    = MAKE_ARRAY(k0+1,w0+1, VALUE=1.)
altOne = MAKE_ARRAY(nk-k0-1, nw-w0-1, VALUE=1.)
BackwardFilter  = FLTARR(nk,nw)
BackwardFilter[0:k0,0:w0] = ONE
BackwardFilter[k0+1:nK-1,w0+1:nw-1] = altOne
IF( zeroK) THEN BackwardFilter[k0, *] = 0.5
IF(zeroW ) THEN BAckwardFilter[*, w0] = 0.5
One = MAKE_ARRAY(nk, nw, VALUE=1, /INT)
ForwardFilter = One - BackwardFilter
;
; Apply Filters to data
;
nx = N_ELEMENTS( *(temp.grid1.values) )
filter = FLTARR(nx, nk, nw)
FOR i=0,nx-1 DO filter[i,*,*] = ForwardFilter
fdata = temp -> Times(Filter)
FOR i=0,nx-1 DO filter[i,*,*] = BackwardFilter
bdata = temp -> Times(Filter)
temp -> TRASH
;
; Invert transforms
;
temp1   = fData -> FFT(2, /Inverse)
forward = temp1 -> FFT(3, /Inverse)
temp1 -> trash
temp1    = bdata -> FFT(2, /Inverse)
backward = temp1 -> FFT(3, /Inverse)
;
; fuss with units, etc.
;
self -> get, title=title, units=units
self -> get, axis=2, gridunits=grid2units
self -> get, axis=3, gridunits=grid3units
forward  -> set, title=title + "!D+!N"
forward  -> set, axis=2, gridunits=grid2units
forward  -> set, axis=3, gridunits=grid3units
backward -> set, title=title + "!D-!N"
backward -> set, axis=2, gridunits=grid2units
backward -> set, axis=3, gridunits=grid3units
;
; Create output structure
;
output = 	{	Name	:	"Forward",	$
			Forward	:	Forward,	$
			Backward:	Backward	}
RETURN, output
END ; ****** GKVs3D::Forward ****** ;	

