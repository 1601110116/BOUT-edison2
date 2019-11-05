FUNCTION GKV_bMaxFcn, bMax
;
COMMON bMaxInfo, VshieldConst
result = ALOG(1.+2.*bMax) - VshieldConst
RETURN, result
END  ;  ******  GKV_bMaxFcn  ******  ;


Function NoiseCoefficients, _Extra=extra
;
; Purpose:
;
;	This function computes coefficients needed to
;	evaluate the fluctuation intentity and fluctuation
;	energy associated with discrete particle noise.
;
;	These coefficients can be viewed as definite integrals
;	over a 3-D k-space.
;
; Keywords:
;
;	k_max=kmax	Noise is computed over range -k_max < k_perp < k_max.
;			defaults to pi/dx.
;
;	kxcut		Value of k_x at which the filter "cuts off". 
;			Defaults to 0.7.
;
;	kycut		Value of k_y at which the filter "cuts off". 
;			Defaults to 0.7.
;
;	G = Gin		Parameter in Hammett's computation of ion self-screening.
;			Defaults to 1.0.
;
;	d_nz = d_nzIn	Parameter in Hammett's computation of electron shielding.
;			Acts only on noise computation at k_y=0.  Set to 1 to
;			allow Debye shielding of zonal flows, and set to 0 to
;			prevent Debyeshielding of zonal flows.  Defaults to 1.
;
;	dx = dxIn	Transverse (to B) grid size in units of rho.
;			Defaults to 1.
;
;	ell = ellIn or	(synomyms).  Describes the filter.  In 'x'  
;	  l = ellIn	direction, the filter is 1/(1.+b^ell[0])
;			in the 'y' direction the filter
;			is exp(-(k_y*a_y)^(2*ell[1])), and 
;			in the z-direction the filter is
;			1/(1.+b^ell[2]).  If a scalar ellIn is a
;			scalar, then this scalar is used in all three
;			directions.  If ellIn is of length 2, then
;			ellIn[2] is set equal to ellIn[0].
;
;	a_x = a_xIn	Used in defining the filter in the x-direction.
;			Defaults to 1.
;
;	a_y = a_yIn	Used in defining the filter in the y-direction. 
;			Defaults to 1.
;
;	a_z = a_zIn	Used in defining the filter in the y-direction. 
;			Defaults to 1.
;
;	DifSq = DifSq	"Set" this keyword (i.e., put /DifSq on command
;			line, or set DifSq=1) to include NGP particle 
;			shapefunctions in noise Estimates.
;
;	TOL = tol	Tolerated error in root finder.  
;			Defaults to 1.e-4.
;
;
; Written by W.M. Nevins
;  3/22/05
;
PRINT, "********************************************"
PRINT, "****** YOU HAVE INVOKED THE GEM Noise ******"
PRINT, "******      Coefficient routine       ******"
PRINT, "********************************************"
;
; Parse Input line (which is now found in "Extra")
;
COMMON bMaxInfo, VshieldConst
DifSq = 0b
result = GetKeyWord('DifSq', extra)
IF( Query_Integer(result) ) THEN DifSq = result

dx = 1.
result = GetKeyWord('dx', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN dx=FLOAT(result)

k_max = !PI/dx
result = GetKeyWord('k_max', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN k_max=FLOAT(result)

kxcut = 0.7
result = GetKeyWord('kxcut', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN kxcut=FLOAT(result)

kycut = 0.7
result = GetKeyWord('kycut', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN kycut=FLOAT(result)

g = 1.
result = GetKeyWord('g', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN g=FLOAT(result)

d_nz = 1.
result = GetKeyWord('d_nz', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN d_nz=FLOAT(result)

tol = 1.e-4
result = GetKeyword("TOL", extra)
IF( Query_Real(result)) THEN tol = result

ell = 3
result = GetKeyWord('ell', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN ell=result
result = GetKeyWord('l', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN ell=result
CASE N_ELEMENTS(ell) OF 
	1	:	ell = ell*MAKE_ARRAY(3, /LONG, VALUE=1.)
	2	:	ell = [ell[0], ell[1], ell[0]]
	3	:
ENDCASE

a_x = 1.
result = GetKeyWord('a_x', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN a_x=FLOAT(result)

a_y = 1.
result = GetKeyWord('a_y', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN a_y=FLOAT(result)

a_z = 1.
result = GetKeyWord('a_z', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN a_z=FLOAT(result)
;
; Set up arrays of k-values
;
k_perp = -k_max + (2.*k_max/100.)*FINDGEN(101)
k_z = -!PI + (2.*!PI/100.)*FINDGEN(101)
;
; Compute interesting intermediates
;
b = k_perp*k_perp
b_num = 2.*(1.-cos(k_perp*dx))/(dx*dx)
b_z = 2.*(1.-cos(k_z))	; dz is taken to be one...
;
;sGamma_0t = BESELI(b,0)*exp(-b)	; Theoretical version of Gamma_0
;sGamma_0 = 1./(1.+b)        	; (version of Gammma_0 actually used in PG3EQ)
;sGamma_0x= 1./(1.+b_num)    	; (but finite-differenced in x-direction)
;
; A note on aliasing and difSq
;
;	The potential, and so also <phi(k)^2>, must be periodic in the
;	grid vector.  This is automatically true of all operators which
;	are applied on the grid (like b_num, ...).  Not obviously true of
;	particle weighting functions (like DifSq, etc), where periodicity
;	in the grid vector comes as a consequence of the sum over grid
;	aliases.  
;
;	For nearest-grid-point weighting (what is employed in PG3EQ),
;	there is a remarkable identity,
;
;	Sum_p{DifSq[(k+pk_g)*dx/2]} = 1 (!)
;
;	This means that we can allow for the grid-aliasing by
;	replacing the square of the particle weighting factor
;	(in k-space) everywhere it appears by 1.0.  
;
;
; Compute k-space representation of parallel derivative ...
;
d_parallel = Sin(k_z)/k_z
d_parallel[50] = 0.
;
; Compute k_perp-space representation of Gamma_0
;
one = MAKE_ARRAY(101, /FLOAT, value=1.0)
b_perp =  one#b + b_num#one
Gamma_0 = 1./(1. + b_perp)
;altGamma_0 = 1./(1. + 0.5*b_perp)
;
; Compute k-space representation of energy operator
;
Eop = one#b_num + b_num#one
;
; Compute Difs
;
Dif_x = ( SIN(k_perp*dx/2.)/(k_perp*dx/2.) )
Dif_x[50] = 1.
Dif_x[50] = 1.
Dif_y = ( SIN(k_perp*dx/2.)/(k_perp*dx/2.) )
Dif_y[50] = 1.
Dif_z = ( SIN(k_z/2.)/(k_z/2.) )
DIf_z[50] = 1.
KappaSq_x = (SIN(k_perp*dx)/dx)^2
KappaSq_y = (SIN(k_perp*dx)/dx)^2
KappaSq = KappaSq_x#one + one#KappaSq_y
kxSQ = k_perp*k_perp
kySq = k_perp*k_perp
barekSq = kxSq#One + one#kySq
;
; Compute k-space representation of filters applied on the grid
; (as used in GEM)
;
filter_y = FLOAT(ABS(k_perp) LE kyCut)
IF(DifSq) THEN filter_y = filter_y*Dif_y
filter_x= FLOAT(ABS(k_perp) LE kyCut)
IF(DifSq) THEN filter_x = filter_x*Dif_x
filter_z= ( 0.5*(1.+COS(k_z))*1.5*(1.-(1./3.)*COS(k_z)) )^ell[2]
IF(DifSq) THEN filter_z = filter_z*Dif_z
;
; And form k_perp-space representation;
;
filter_perp = filter_x#filter_y
;
; What we're reduced to is a DO-loop on k_z to compute 
; needed coefficients (MUST BE A BETTER WAY!)
;
quotientHammett = 0.
quotientNevins = 0.
energyHammett = 0.
energyNevins = 0.
kSqHammett = 0.
kSqBare = 0.
kSqVolume = 0.
qHammett_z = FLTARR(101,101)
 qNevins_z = FLTARR(101,101)
Noise_HxAvg = 0.*k_perp
Noise_HyAvg = 0.*k_perp
Noise_NxAvg = 0.*k_perp
Noise_NyAvg = 0.*k_perp
FOR i=0,99 DO BEGIN
	qHammett = filter_z[i]^2*filter_perp^2*Gamma_0/( (2.-Gamma_0)*(2.-(1.-g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0) )
	qNevins  = filter_z[i]^2*filter_perp^2*Gamma_0/( (2.-Gamma_0)*(2.-Gamma_0) )
	kSq_Hammett = 	3.0*KappaSQ*filter_z[i]^2*filter_perp^2*Gamma_0^2/	$
			( (2.-(1.-g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0)*(2.-(1.-0.5*g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0) )
	kSq_bare = 	3.0*barekSQ*filter_z[i]^2*filter_perp^2*Gamma_0^2/	$
			( (2.-(1.-g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0)*(2.-(1.-0.5*g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0) )
	
	kSq_Volume = 	3.0*filter_z[i]^2*filter_perp^2*Gamma_0^2/	$
			( (2.-(1.-g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0)*(2.-(1.-0.5*g*d_parallel[i]*filter_z[i]*filter_perp)*gamma_0) )
	quotientHammett = quotientHammett + TOTAL(qHammett[0:99,0:99])
	quotientNevins  =  quotientNevins + TOTAL(qNevins[0:99,0:99])
	energyHammett = energyHammett + TOTAL((Eop*qHammett)[0:99,0:99])
	energyNevins  =  energyNevins + TOTAL((Eop*qNevins )[0:99,0:99])
	kSqHammett = kSqHammett + TOTAL(kSq_Hammett[0:99,0:99])
	kSqVolume = kSqVolume + TOTAL(kSQ_Volume[0:99,0:99])
	kSqBare = kSqBare + TOTAL(kSq_bare[0:99,0:99])
	qHammett_z  = qHammett_z + qHammett
	qNevins_z   =  qNevins_z + qNevins 
ENDFOR
quotientHammett = quotientHammett*(2.*k_max/100.)*(2.*k_max/100.)*(2.*!PI/100.)
quotientNevins  =  quotientNevins*(2.*k_max/100.)*(2.*k_max/100.)*(2.*!PI/100.)
VolumeHammett = (2.*!PI)^3/quotientHammett
VolumeNevins  = (2.*!PI)^3/quotientNevins
energyHammett   =   energyHammett*(2.*k_max/100.)*(2.*k_max/100.)*(2.*!PI/100.)
energyNevins    =    energyNevins*(2.*k_max/100.)*(2.*k_max/100.)*(2.*!PI/100.)
KappaSqHammett = energyHammett/quotientHammett
KappaSqNevins  =  energyNevins/quotientNevins
kSqHammett = kSqHammett/kSqVolume
kSqBare = kSqBare/kSqVolume
kSqVolume = (2.*!PI)^3/( (2.*k_max/100.)*(2.*k_max/100.)*(2.*!PI/100.)*kSqVolume )
Noise_HxAvg = TOTAL(qHammett_z,2)
Noise_HyAvg = TOTAL(qHammett_z,1)
Noise_NxAvg = TOTAL( qNevins_z,2)
Noise_NyAvg = TOTAL( qNevins_z,1)
Noise_HxAvg     =     Noise_HxAvg*(2.*k_max/100.)*(2.*!PI/100.)
Noise_HyAvg     =     Noise_HyAvg*(2.*k_max/100.)*(2.*!PI/100.)
Noise_NxAvg     =     Noise_NxAvg*(2.*k_max/100.)*(2.*!PI/100.)
Noise_NyAvg     =     Noise_NyAvg*(2.*k_max/100.)*(2.*!PI/100.)
;
; Compute constants needed for Hammett's computation of D_noise
;
kpplMax = (!PI/100.)*TOTAL(filter_z^2)
VshieldConst = (2./!PI)*quotientHammett/(kpplMax)
bMax_guess = 0.5*VshieldConst
bMax = FX_ROOT([bMax_guess, 1.e-20, !PI], "GKV_bMaxFcn", /DOUBLE, TOL=tol)

;
output = {	Name		:	"NoiseCoefficients",		$
		quotientHammett	:	quotientHammett,		$
		quotientNevins	:	quotientNevins,			$
		energyHammett	:	energyHammett,			$
		energyNevins	:	energyNevins,			$
		Noise_HxAvg	:	Noise_HxAvg,			$
		Noise_HyAvg	:	Noise_HyAvg,			$
		Noise_NxAvg	:	Noise_NxAvg,			$
		Noise_NyAvg	:	Noise_NyAvg,			$
		kpplMax		:	kpplMax,			$
		bMax		:	bMax,				$
		VolumeHammett	:	VolumeHammett,			$
		VolumeNevins	:	VolumeNevins,			$
		KappaSqHammett	:	KappaSqHammett,			$
		KappaSqNevins	:	KappaSqNevins,			$
		kSqHammett	:	kSqHammett,			$
		kSqBare		:	kSqBare,			$
		kSQVolume	:	kSqVolume			}
;
RETURN, output
END ;  ****** Function NoiseCoefficients ******  ;	
