FUNCTION GKVs3D::Eop
;
;  Purpose:
;
;	Computes energy density as a function of time 
;	from fluctuating potential in plane perpendicular
;	to B.
;
;	This function assumes that "self" contains the
;	fluctuating potential, and that system has 
;	"open" periodic boundary conditions in the
;	perpendiculat (to B) plane.
;
;  Written by W.M. Nevins
;	3/22/05
;
FORWARD_FUNCTION NoiseCoefficients
grid1 = self.grid1
grid1 = self.grid2
dx = (*self.grid1.values)[1] - (*self.grid1.values)[0]
dy = (*self.grid2.values)[1] - (*self.grid2.values)[0]

Nx = N_ELEMENTS(*self.grid1.values)
Ny = N_ELEMENTS(*self.grid2.values)
Lx = (*self.grid1.values)[Nx-1] - (*self.grid1.values)[0] + dx
Ly = (*self.grid2.values)[Ny-1] - (*self.grid2.values)[0] + dy

phi = *self.values
rho = 	  ( SHIFT(phi,1,0,0) - 2.*phi + SHIFT(phi,-1,0,0) )/(dx*dx) 	$
	+ ( SHIFT(phi,0,1,0) - 2.*phi + SHIFT(phi,0,-1,0) )/(dy*dy)
temp = TOTAL(-1.*rho*phi,1)
E = TOTAL(temp,1)/(Lx*Ly)
Estr = {GKVs1D}
nTags = N_TAGS(Estr)
FOR i=0,nTags-1 DO Estr.(i) = self.(i)
Estr.mnemonic = "E_phi"
Estr.title = "!12E!D!4u!X!N"
Estr.Indices = PTR_NEW(["*"])
Estr.Units = "(!4q!X/L!DT!N)!U2!NnT"
Estr.values = PTR_NEW(E)
vMax = MAX(E)
Estr.Vrange = [0.,vMax]
Estr.Grid1 = GKVsd_GridCopy(self.Grid3)
Eobj = OBJ_NEW("GKVs1D", Estr)

IObj = Eobj -> MakeCopy(/NoValues)
temp = TOTAL(phi*phi,1)
Intensity = TOTAL(temp,1)/(Lx*Ly)
IObj -> Set, 	title="!8I!D!4u!X!N", mnemonic = "I_phi", 	$
			values = PTR_NEW(Intensity), 			$
			Units ="(!4q!X/L!DT!N)!U2!N(T/e)!U2!N" 
Output = {	Name	:	"Eop",		$
		Ephi	:	Eobj,		$
		Iphi	:	Iobj		}
RETURN, output
END  ;  ****** FUNCTION GKVs3D::Eop ******  ;
 

FUNCTION GKVs1D::WhiteNoise, _Extra=extra
;
; Acts on Wsq(t) and computes the spectral density
; of dPhi associated with white noise at the given
; amplitude of wsq(t)
;
; Effective keywords:
;
;	k_max=kmax	Noise is computed over range -k_max < k_perp < k_max.
;			defaults to pi/dx.
;
;	G = Gin		Parameter in Hammett's computation of ion self-screening.
;			Defaults to 1.0.
;
;	d_nz = d_nzIn	Parameter in Hammett's computation of electron shielding.
;			Acts only on noise computation at k_y=0.  Set to 1 to
;			allow Debye shielding of zonal flows, and set to 0 to
;			prevent Debyeshielding of zonal flows.  Defaults to 1.
;
;	Norm=NormIn	Set to 1/(effective number of particles)
;
;	Nx		Number of grid points in x direction
;			(defaults to 128.)
;
;	Ny		Number of grid points in y direction
;			(defaluts to 128.)
;
;	Npx = NpxIn		Alternatively, compute by setting Npx to 
;				number of particles in x (Defaults to 128).
;
;	Npy = NyIn		Npy to number of particles in y (Defaults to 128),
;
;	Ncell_z = NcellzIn	and Ncell_z to number of particles/cell in z.
;				(Defaults to 1).
;
;	dx = dxIn	Transverse (to B) grid size in units of rho.
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
;
; Written by W.M. Nevins
;  1/6/05
; Revised 2/27/05 to reflect Hammett's improved calculation of 
; the noise spectrum.
;
;
result = GetKeyWord('Norm', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN BEGIN
	norm=result
ENDIF ELSE BEGIN
	Npx=128.
	result = GetKeyWord('Npx', extra)
	IF( Query_Real(result) + Query_Integer(result) ) THEN Npx=FLOAT(result)
	Npy = 128.
	result = GetKeyWord('Npy', extra)
	IF( Query_Real(result) + Query_Integer(result) ) THEN Npy=FLOAT(result)
	Ncell_z = 1.
	result = GetKeyWord('NCell_z', extra)
	IF( Query_Real(result) + Query_Integer(result) ) THEN Ncell_z=FLOAT(result)
	Norm = 1./(Npx*Npy*Ncell_z)
ENDELSE

Nx=128
result = GetKeyWord('Nx', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN Nx=FIX(result)

Ny=128
result = GetKeyWord('Ny', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN Ny=FIX(result)

DifSq = 0b
result = GetKeyWord('DifSq', extra)
IF( Query_Integer(result) ) THEN DifSq = result

dx = 0.9817
result = GetKeyWord('dx', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN dx=FLOAT(result)

k_max = !PI/dx
result = GetKeyWord('k_max', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN k_max=FLOAT(result)

g = 1.
result = GetKeyWord('g', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN g=FLOAT(result)

d_nz = 1.
result = GetKeyWord('d_nz', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN d_nz=FLOAT(result)

ell = 2
result = GetKeyWord('ell', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN ell=result
result = GetKeyWord('l', extra)
IF( Query_Real(result) + Query_Integer(result) ) THEN ell=result
CASE N_ELEMENTS(ell) OF 
	1	:	ell = ell*MAKE_ARRAY(3, /LONG, VALUE=1.)
	2	:	ell = [ell[0], ell[1], ell[0]]
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

k_perp = -k_max + (2.*k_max/100.)*FINDGEN(101)
k_z = -!PI + (2.*!PI/100.)*FINDGEN(101)
kGrid = {Grid}
kGrid.mnemonic = 'k_perp'
kGrid.title = 'k!D!9x!X!N!4q!X'
kGrid.values = PTR_NEW(k_perp)
kGrid.boundary = 'Open'
kGrid.uniform = 1b
kGrid.range= [-k_max, k_max]
kGrid.irange = [0,100]

b = k_perp*k_perp
b_num = 2.*(1.-cos(k_perp*dx))/(dx*dx)
b_z = 2.*(1.-cos(k_z))
sGamma_0t = BESELI(b,0)*exp(-b)
sGamma_0 = 1./(1.+b)        ; (version of Gammma_0 actually used in PG3EQ)
sGamma_0x= 1./(1.+b_num)    ; (but finite-differenced in x-direction)
;
; A note on aliasing and difSq
;
;	The potential, and so also <phi(k)^2>, must be periodic in the
;	grid vector.  This is automatically true of all operators which
;	are applied on the grid (like b_num, ...).  Not obviously true of
;	particle weighting functins (like DifSq, etc), where periodicity
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
;	This allowed me to remove all of the coding below focused on 
;	summing over grid aliases.
;	
;
sdif = SIN(k_perp*dx/2.)/(k_perp*dx/2.)
sdif[50] = 1.0
;sdifSq = sdif*sdif
;
; add poor-man's correction for aliasing in x
;
;kg = 2*!PI/dx
one_x = 0.*k_perp + 1.
;p = FINDGEN(100) + 1.
;one_2 = 0.*p + 1.
;kp = kg*p
;term = b#one_2*( 1./(k_perp#one_2 + one_x#kp)^2 + 1./(k_perp#one_2 - one_x#kp)^2 )
;correction = one_x + TOTAL(term, 2)
csdifSq = one_x
;
zDif = sin(k_z/2.)/(k_z/2.)
zDif(50) = 1.0
zDifSq = zDif*zDif
;
; add poor-man's correction for aliasing in z
;
;kg = 2*!PI
one_z = 0.*k_z + 1.
;one_2 = 0.*p + 1.
;kp = kg*p
;term = b#one_2*( 1./(k_z#one_2 + one_z#kp)^2 + 1./(k_z#one_2 - one_z#kp)^2 )
;term = one_z#one_2*( 1./(k_z#one_2 + one_z#kp)^2 + 1./(k_z#one_2 - one_z#kp)^2 )
;correction = one_x + TOTAL(term, 2)
;czdifSq = zDifSq*correction
;czdifSq = one_z
;
; Compute k-space representatin of parallel derivative ...
;
zd_parallel = Sin(k_z)/k_z
zd_parallel[50] = 0.
;
; add aliasing to this?
;
;kg = 2*!PI
;p = FINDGEN(100) + 1.
;one_2 = 0.*p + 1.
;kp = kg*p
;term = (k_z#one_2)*( 1./(k_z#one_2 + one_z#kp) + 1./(k_z#one_2 - one_z#kp) )
;correction = one_z + TOTAL(term, 2)
;zd_parallel = zd_Parallel*correction
;
; Compute k-space representation of filters applied on the grid
;
power=2*ell[1]
filter_y = exp(-(k_perp*a_y*dx)^power)
IF(DifSq) THEN filter_y=filter_y*sdif
filter_x= 1./(1. + (a_x^2*b_num*dx^2)^ell[0])
IF(DifSq) THEN filter_x=filter_x*sdif
filter_z= 1./(1. + (a_z^2*b_z  )^ell[2])
IF(DifSq) THEN filter_z=filter_z*zdif
one = MAKE_ARRAY(101, /FLOAT, value=1.0)
filter = filter_y#filter_z
filterx= filter_x#filter_z
gamma_0 = sgamma_0#one
gamma_0x= sGamma_0x#one
;dif = sdif#zDif
; difsq = sdifsq#zDifsq
; difsq = csdifsq#czDifsq
d_parallel = one#zd_parallel

quotientHammett = filter^2*Gamma_0/( (2.-Gamma_0)*(2.-(1.-g*d_parallel*filter)*gamma_0) )
quotientNevins  = filter^2*Gamma_0/( (2.-Gamma_0)*(2.-gamma_0) )
quotientHammett = TOTAL(quotientHammett, 2)/100.
quotientNevins  = TOTAL(quotientNevins , 2)/100.

quotientHammettx = filterx^2*Gamma_0x/( (1.+d_nz-Gamma_0x)*(1.+d_nz-(1.-g*d_parallel*filterx)*gamma_0x) )
quotientNevinsx  = filterx^2*Gamma_0x/( (1.+d_nz-Gamma_0x)*(1.+d_nz-gamma_0x) )
quotientHammettx = TOTAL(quotientHammettx, 2)/100.
quotientNevinsx  = TOTAL(quotientNevinsx , 2)/100.
;
; NoiseIntensity estimate based on integral over 2-D k_perp space.
; However, recall that we are in the "sum" norm, so k_perp needs to
; be put into units where dk_perp*rho_i = 2*!PI/100.
; so 2*!PI*k_perp*dk_perp -> 2.*!PI*(k_perp*rho_i/d_kperp*rho_i)*(d_kperp*rho_i/d_kperp*rho_i) = 100.*(k_perp*rho_i)
;
IntensityCoefficient = (Nx/100.)*(Ny/100.)*TOTAL( 100.*k_perp[50:100]*quotientHammett[50:100] )

wSq = *self.values
nT = N_ELEMENTS(wSq)
values = FLTARR(101, nt)
FOR i=0L,nt-1 DO values[*,i] = norm*wSq[i]*quotientHammett

Intensity = FLTARR(nt)
FOR i=0L,nt-1 DO Intensity[i] = norm*wSq[i]*IntensityCoefficient
IntensityStr = {GKVs1D}
IntensityStr.mnemonic = "Intensity"
IntensityStr.title    = "!8I!3!Dnoise!N"
IntensityStr.indices  = PTR_NEW(["*"])
IntensityStr.units    = "(!4q!X/L!DT!N)!U2!N(T/e)!U2!N"	; Need to revisit the units!
IntensityStr.values   = PTR_NEW(Intensity)
 vmax = MAX(Intensity)
IntensityStr.vrange   = [0.,vmax]
IntensityStr.codename = self.codename
IntensityStr.CodePi   = self.codePI
IntensityStr.RunId    = self.runID
IntensityStr.FileId   = self.FileID
IntensityStr.Grid1    = GKVsd_GridCopy(self.grid1)
IntensityObj = OBj_NEW("GKVs1D", IntensityStr)
;
; A more elaborate Intensity (and energy) calculation
;
zSize = TOTAL(filter_z^2)/100.
eFilter = filter_x#filter_y
;ed_parallel = TOTAL(zd_parallel*filter_z)/100.
ed_Parallel = 1.
eb = b_num#one_x + one_x#b
eeb = b_num#one_x + one_x#b_num
eGamma_0 = 1./(1.+eb)
eQuotient = eFilter^2*zSize*eGamma_0/( (2. - eGamma_0)*(2. - (1.-g*ed_parallel*eFilter)*eGamma_0) )
altIntensityCoefficient = (Nx/100.)*(Ny/100.)*TOTAL(eQuotient)
altEnergyCoefficient    = (Nx/100.)*(Ny/100.)*TOTAL(eeB*eQuotient)
altIntensity = FLTARR(nt)
altEnergy    = FLTARR(nt)
FOR i=0L,nt-1 DO BEGIN
	altIntensity[i] = norm*wSq[i]*altIntensityCoefficient
	   altEnergy[i] = norm*wSq[i]*altEnergyCoefficient
ENDFOR
altIntensityObj = IntensityObj -> MakeCopy(/noValues)
vMax = MAX(altIntensity)
altIntensityObj -> set, values = PTR_NEW(altIntensity), vrange=[0.,vmax]
altEnergyObj = IntensityObj -> MakeCopy(/NoValues)
vMax = MAX(altEnergy)
altEnergyObj -> set,	values = PTR_NEW(altEnergy), Title="!12E!3!DNoise!N", 	$
			mnemonic="Energy_noise", Units="(!4q!X/L!DT!N)!U2!NnT",	$
			vrange=[0.,vMax]

;
; Yet ANOTHER effort at intensity and energy associated with the noise
;
Ncoeff = NoiseCoefficients(	k_max=k_max, g=g, d_nz=d_nz, dx=dx, 	$
				ell=ell, a_x=a_x, a_y=a_y, a_z=a_z, DifSq=DifSq)
nIntensity  = FLTARR(nt)
nEnergy     = FLTARR(nt)
Noise_HxAvg = FLTARR(101, nt)
Noise_HyAvg = FLTARR(101, nt)
Noise_NxAvg = FLTARR(101, nt)
Noise_NyAvg = FLTARR(101, nt)
FOR i=0L,nt-1 DO BEGIN
	nIntensity[i] = norm*wSq[i]*(0.5*Nx/k_max)*(0.5*Ny/k_max)*(1./(2.*!PI))*Ncoeff.quotientHammett
	   nEnergy[i] = norm*wSq[i]*(0.5*Nx/k_max)*(0.5*Ny/k_max)*(1./(2.*!PI))*Ncoeff.energyHammett
     Noise_HxAvg[*,i] = norm*wSq[i]*(0.5*Ny/k_max)*(1./(2.*!PI))*Ncoeff.Noise_HxAvg
     Noise_HyAvg[*,i] = norm*wSq[i]*(0.5*Nx/k_max)*(1./(2.*!PI))*Ncoeff.Noise_HyAvg
     Noise_NxAvg[*,i] = norm*wSq[i]*(0.5*Ny/k_max)*(1./(2.*!PI))*Ncoeff.Noise_NxAvg
     Noise_NyAvg[*,i] = norm*wSq[i]*(0.5*Nx/k_max)*(1./(2.*!PI))*Ncoeff.Noise_NyAvg
ENDFOR
nIntensityObj = IntensityObj -> MakeCopy(/noValues)
vMax = MAX(nIntensity)
nIntensityObj -> set, values = PTR_NEW(nIntensity), vrange=[0.,vMax]
nEnergyObj = IntensityObj -> MakeCopy(/NoValues)
vMax = MAX(nEnergy)
nEnergyObj -> set,	values = PTR_NEW(nEnergy), Title="!12E!3!DNoise!N", 	$
			mnemonic="Energy_noise", Units="(!4q!X/L!DT!N)!U2!NnT",	$
			vrange=[0.,vMax]


Hammett = {GKVs2D}
Hammett.mnemonic = "noise_H"
Hammett.title    = "!13<!4du!X!U2!N!13>!X!U(H)!N"
Hammett.indices  = PTR_NEW(["k!Dx!N=0", "*","*"])
Hammett.units    = self.units + '(T/e)!U2!N'
Hammett.values   = PTR_NEW(values)
 vmax = MAX(values)
Hammett.vrange   = [0.,vmax]
Hammett.codename = self.codename
Hammett.CodePi   = self.codePI
Hammett.RunId    = self.runID
Hammett.FileId   = self.FileID
Hammett.Grid1    = kGrid
Hammett.Grid2    = GKVsd_GridCopy(self.grid1)
Noise_H = OBj_NEW("GKVs2D", Hammett)
Noise_H -> Set, axis=1, GridMnemonic="k_y", GridTitle="k!Dy!N"

Noise_HyAvgObj = Noise_H -> MakeCopy(/noValues)
vMax = MAX(noise_HyAvg)
Noise_HyAvgObj -> Set,	title="!13<!4du!X!U2!N!13>!S!X!Dx!R!U(H)!N",		$
			mnemonic="Noise_HyAvg", values=PTR_NEW(Noise_HyAvg),	$
			vrange=[0,vMax], Indices=PTR_NEW(["*","*"])

Noise_Hx = Noise_H -> MakeCopy(/noValues)
FOR i=0L,nt-1 DO values(*,i) = norm*wSq(i)*quotientHammettx
vMax = MAX(values)
Noise_Hx -> set, title="!13<!4du!X!U2!N!13>!X!U(H)!N",	$
		Mnemonic = "noise_Hx",			$
		values = PTR_NEW(values),		$
		vrange = [0, vmax],			$
		Indices = PTR_NEW(["*","k!Dy!N=0","*"])
Noise_Hx -> Set, axis=1, GridMnemonic="k_x", GridTitle="k!Dx!N"

Noise_HxAvgObj = Noise_Hx -> MakeCopy(/NoValues)
vMax = MAX(Noise_HxAvg)
Noise_HxAvgObj-> Set,	title="!13<!4du!X!U2!N!13>!S!X!Dy!R!U(H)!N",		$
			mnemonic="Noise_HxAvg", values=PTR_NEW(Noise_HxAvg),	$
			vrange=[0,vMax], Indices=PTR_NEW(["*","*"])

Noise_N = Noise_H -> MakeCopy(/noValues)
FOR i=0L,nt-1 DO values(*,i) = norm*wSq(i)*quotientNevins
vMax = MAX(values)
Noise_N -> set, title="!13<!4du!X!U2!N!13>!X!U(N)!N",	$
		Mnemonic = "noise_N",				$
		values = PTR_NEW(values),			$
		vrange = [0, vmax]

Noise_NyAvgObj = Noise_N -> MakeCopy(/NoValues)
vMax = MAX(Noise_NyAvg)
Noise_NyAvgObj-> Set,	title="!13<!4du!X!U2!N!13>!S!X!Dy!R!U(N)!N",		$
			mnemonic="Noise_NyAvg", values=PTR_NEW(Noise_NyAvg),	$
			vrange=[0,vMax], Indices=PTR_NEW(["*","*"])

Noise_Nx = Noise_Hx -> MakeCopy(/noValues)
FOR i=0L,nt-1 DO values(*,i) = norm*wSq(i)*quotientNevinsx
vMax = MAX(values)
Noise_Nx -> set, title="!13<!4du!X!U2!N!13>!X!U(N)!N",	$
		Mnemonic = "noise_Nx",				$
		values = PTR_NEW(values),			$
		vrange = [0, vmax]

Noise_NxAvgObj = Noise_Nx -> MakeCopy(/NoValues)
vMax = MAX(Noise_NxAvg)
Noise_NxAvgObj-> Set,	title="!13<!4du!X!U2!N!13>!S!X!Dy!R!U(N)!N",		$
			mnemonic="Noise_HxAvg", values=PTR_NEW(Noise_NxAvg),	$
			vrange=[0,vMax], Indices=PTR_NEW(["*","*"])

filterObj = Noise_N -> slice(axis=2, index=0)
filterObj -> get,	values=OldValues, Indices=OldIndices
PTR_FREE, OldValues
PTR_FREE, OldIndices
filterObj -> set,	title="!12F!X",		$
			mnemonic="Filter",	$
			units = "",		$
			values=PTR_NEW(filter_y),	$
			vrange=[0, 1.],		$
			Indices=PTR_NEW(["*"])
;
;DifObj = FilterObj -> MakeCopy(/NoValues)
;DifObj -> set,	title = "Dif",		$
;		Mnemonic = "Dif",	$
;		values = PTR_NEW(sdif)
;
Gamma_Obj = FilterObj -> MakeCopy(/NoValues)
Gamma_Obj -> set,	Title="!4C!X!S!I0!R!E(num)!N",	$
			mnemonic = "Gamma_0num",	$
			values = PTR_NEW(sGamma_0)

GammaTheoryObj = Gamma_Obj -> MakeCopy(/NoValues)
GammaTheoryObj -> Set,	Title="!4C!X!S!I0!R!E(th)!N",	$
			mnemonic = "Gamma_0th",		$
			values = PTR_NEW(sGamma_0t)

Nevins = Gamma_Obj -> MakeCopy(/NoValues)
Nevins -> Set,	Title="!12N!X",		$
		mnemonic="Nevins",	$
		values = PTR_NEW(QuotientNevins)

Hammett = Gamma_Obj -> MakeCopy(/NoValues)
Hammett -> Set,	Title="!12H!X",		$
		mnemonic="Hammett",	$
		values = PTR_NEW(QuotientHammett)




result = {	Name	:	"WhiteNoise",		$
		Noise_H	:	Noise_H,		$
		Noise_N	:	Noise_N,		$
		Noise_Hx:	Noise_Hx,		$
		Noise_Nx:	Noise_Nx,		$
	    Noise_HxAvg	:	Noise_HxAvgObj,		$
	    Noise_HyAvg	:	Noise_HyAvgObj,		$
	    Noise_NxAvg	:	Noise_NxAvgObj,		$
	    Noise_NyAvg	:	Noise_NyAvgObj,		$
		Filter	:	FilterObj,		$
;		Dif	:	DifObj,			$
		Gamma_0num :	Gamma_Obj,		$
		Gamma_0th  :	GammaTheoryObj,		$
		Nevins	:	Nevins,			$
		Hammett	:	Hammett,		$
	      Intensity	:	IntensityObj,		$
	   altIntensity :	altIntensityObj,	$
		Energy	:	altEnergyObj,		$
	     nIntensity	:	nIntensityObj,		$
	        nEnergy	:	nEnergyObj		}

return, result
END 
