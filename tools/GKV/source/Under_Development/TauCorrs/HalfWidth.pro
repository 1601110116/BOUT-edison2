PRO GKVsd::Norm, Maxvalue=maxval
maxvalue=1.
IF(KEYWORD_SET(maxval)) THEN maxvalue=maxval
maximum = MAX(*self.values)
values = maxvalue*(*self.values)/maximum
PTR_FREE, self.values
self.values = PTR_NEW(values)
vmin=GKVsd_Min(values, Max=vmax)
self.vrange=[vmin,vmax]
END


FUNCTION GKVs1D::FullWidth, debug=debug
;
; Compute full Width of 'self'
; using quadratic interpolation
; Interpolation algorithm does NOT assume a uniform grid (in time)
; ... but algorithm for computing 2-time correlation function does.
;
; Returns correlation time (full-width at half-maximum)
;
; Written by W.M. Nevins
;	4/15/00
;
imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
values=(*self.values)[imin:imax]
x=(*self.Grid1.values)[imin:imax]
maxVal=MAX(values, imaxVal)
epsilon=maxVal*1.e-6
xmax=values[imaxVal]
jmax=N_ELEMENTS(values) - 1
;
; Find grid point nearest half-value for indices greater than jmax
;
temp=values - (maxVal/2.)
temp2=temp^2
vUp  = MIN(temp2[imaxVal:jmax], jUp )
jUp = imaxVal + jUp
IF(jUP EQ jmax) THEN BEGIN
	dxUp = 0.
	dxtotal = 2.*(x[jUp] - x[jUp-1])
	GOTO, DoneUp
ENDIF
IF(ABS(vUp) LT epsilon) THEN BEGIN
	dxUP = 0.
	dxtotal= x[jUp+1] - x[jUp-1]
	GOTO, DoneUP
ENDIF
;
; Now, compute coefficients for quadratic interpolation about jUp
;
dxplus = x[jUp+1] - x[jUp]
dxminus= x[jUp]   - x[jUp-1]
dxtotal= x[jUp+1] - x[jUp-1]
dVplusdx = (temp[jUp+1] - temp[jUp  ])/dxplus
dVminusdx= (temp[jUp]   - temp[jUp-1])/dxminus
c = temp[jUp]
b = (dxplus*dvminusdx + dxminus*dvplusdx)/dxtotal
a = (dVplusdx - dVminusdx)/dxtotal
discriminant = 1 - 4.*a*c/(b^2)
IF(discriminant GE 0) THEN BEGIN
	dxUp = -2.*(c/b)/(1.0 + SQRT(discriminant) )		; This is the 'nearby' root to the quadratic
ENDIF ELSE BEGIN
	dxUp = 0									; Just 'punt' if discriminant is negative...								
ENDELSE
DoneUP : IF(ABS(dxUp) GT dxtotal) THEN dxup=0.
xUp = x[jUp] + dxUp
; 
; Same for quadratic interpolation about jback
;
vBack = MIN(temp2[   0:imaxVal], jback)
IF(jback EQ 0) THEN BEGIN
	dxBack = 0.
	dxtotal= 2.*(x[jBack+1] - x[jBack])
	GOTO, DoneBack
ENDIF
IF(ABS(vBack) LT epsilon) THEN BEGIN
	dxBack = 0.
	dxtotal= x[jBack+1] - x[jBack-1]
	GOTO, DoneBack
ENDIF
dxplus = x[jBack+1] - x[jBack]
dxminus= x[jBack]   - x[jBack-1]
dxtotal= x[jBack+1] - x[jBack-1]
dVplusdx = (temp[jBack+1] - temp[jBack  ])/dxplus
dVminusdx= (temp[jBack]   - temp[jBack-1])/dxminus
c = temp[jBack]
b = (dxplus*dvminusdx + dxminus*dvplusdx)/dxtotal
a = (dVplusdx - dVminusdx)/dxtotal
discriminant = 1 - 4.*a*c/(b^2)
IF(discriminant GE 0) THEN BEGIN
	dxBack = -2.*c*dxtotal/(b + SQRT(discriminant) )
ENDIF ELSE BEGIN
	dxBack = 0
ENDELSE
DoneBack : IF(ABS(dxBack) GT dxtotal) THEN dxBack=0.
xBack = x[jBack] + dxBack
fullWidth = ABS(xUp-xBack)
IF(KEYWORD_SET(debug)) THEN PRINT, xUp, dxUp, xBack, dxBack, fullwidth
RETURN, fullWidth
END ; ****** GKVs1D::FullfWidth ****** ;
