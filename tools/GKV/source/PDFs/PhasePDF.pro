FUNCTION GKVs3D::PhasePDFs, arg, kMax=kMax
;
; Acts on a GKVs3D object containing
; data vs. space and time to
; computes 2-D array of phase PDF's 
; for each kx, ky pair with kx, ky < kMax.
;
; Arguments:
;
;	arg	A GKVs3D object on same grid as "self".
;		Phases of self are computed (Fourier mode
;		by Fourier mode) relative to 
;		the phases of arg.
;
; Keywords:
;
;	kMax	Maximum value of kx, ky for
;		which PhasePDFs will be computed.
;		defaults to MIN(kxMax, kyMax)
;		where kxMax = pi/dx, etc.
;
; Written by W.M. Nevins
;	April 17, 2008
;
; Transform 'self' into kx, ky representation.
;
temp   = self -> FFT(1)
self_k = temp -> FFT(2)
temp -> TRASH
temp = arg -> FFT(1)
arg_k = temp -> FFT(2)
temp -> TRASH
k1Range = self_k.grid1.range
k1Max = k1range[1]
k2range = self_k.grid1.range
k2max = k2Range[1]
k1Values = *(self_k.grid1.values)
k2Values  = *(self_k.grid1.values)
dk1 = k1Values[1] - k1Values[0]
dk2 = k2Values[1] - k2Values[0]
IF(N_ELEMENTS(kMax) NE 0) THEN BEGIN
	k1Max = kMax
	k2Max = kMax
ENDIF
nk1 = 2.*k1Max/dk1 + 1
nk2 = k2Max/dk2 + 1
k1Min = -k1Max
k2Min = 0.
output1 = OBJARR(nk1, nk2)
output2 = OBJARR(nk1, nk2)
FOR i=0,nk1-1 DO BEGIN
	FOR j=0,nk2-1 DO BEGIN
		k1 = k1Min + i*dk1
		k2 = k2Min + j*dk2
		sk1 = self_k -> slice(axis=1, value=k1)
		sk  = sk1    -> slice(axis=1, value=k2)
		ak1 = arg_k  -> slice(axis=1, value=k1)
		ak  = ak1    -> slice(axis=1, value=k2)
		output1[i,j] = sk -> phasePDF(ak, /Weight)
		output2[i,j] = sk -> xSpect(ref=ak)
		sk1 -> TRASH
		ak1 -> TRASH
		sk  -> TRASH
		ak  -> TRASH
	ENDFOR
ENDFOR
self_k -> TRASH
arg_k  -> TRASH
output	=	{	Name      :	"PhsePDFs",	$
			phasePDFs :	output1,	$
			  xSpects :	output2		}
RETURN, output
END  ;  ****** GKVs3D::PhasePDFs ******  ;






FUNCTION GKVs1D::Phase
;
; Assumes that "self" is a complex number, 
; comptes the phase of "self"
;
; Arguments:	none
;
; Written by W.M. Nevins
;	4/14/2008
;
x=FLOAT(*self.values)
y=IMAGINARY(*self.values)
phase=ATAN(y,x)
result = self -> MakeCopy(/noValues, /NoErrorbars)
result.title = "!12Arg!X(" + self.title + ")"
result.mnemonic="Arg_" + self.mnemonic
result.units = ""
result.values = PTR_NEW(phase)
result.vrange=[-!PI,!PI]
Return, result
END  ; ****** FUNCTION GKVs1D::Phase ****** ;

FUNCTION GKVs1D::phasePDF, obj, Weight=weight
;
; Computes pdf of the relative phase of "self" 
; and obj
;
; Arguments:
;
;	obj	A GKV object. Should be defined on same spatial/temporal
;		domain as "self".
;
; Keywords:
;
;	Weight	This keyword allows you to compute a "weighted" pdf.
;		if "set" (i.e., "/Weights" or "Weights=1" on command line)
;		the PDF is weighted by self -> absSq(). If you wish  
;		to use a different weighting function, set Weights equal
;		to an object containing the desired weighting function.
;
;
; Written by W.M. Nevins
;	4/14/2008
ratio = self -> over(Obj)
phase = ratio -> phase()

CASE TypeOf(weight) OF
	0:	phasePDF = phase -> PDF(xMin=-!PI, xMax=!PI)
	11:	phasePDF = phase -> PDF(xMin=-!PI, xMax=!PI, Weight=weight)
	ELSE:	BEGIN
		weight = self -> absSq()
		phasePDF = phase -> PDF(xMin=-!PI, xMax=!PI, Weight=weight)
		weight -> Trash
		END
ENDCASE
;
;Clean up
;
ratio -> TRASH
phase -> TRASH
;
RETURN, phasePDF
END ; ****** FUNCTION GKVs1D::phasePDF ****** ;
	
