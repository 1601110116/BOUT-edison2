PRO GKVS1D::MaxNorm
values=*self.values
maxVal=max(values)
values=values/maxVal
self.values = PTR_NEW(values)
self.units=""
RETURN
END



FUNCTION nModes, dataStr, nVals=N, nMax=n_max
nVals=8
IF(N_ELEMENTS(n) EQ 1) THEN nVals=N

nMax=4088
IF(N_ELEMENTS(n_max) EQ 1) THEN nMax=4088

nModes = OBJARR(nVals)
tempArr = OBJARR(nVals)
FOR i=0,nVals-1 DO  nModes[i] = dataStr.potential_r -> slice(n=i*nMax/nVals)
FOR i=0,nVals-1 DO tempArr[i] = nModes[i] -> AbsSq()
FOR i=0,nVals-1 DO  nModes[i] -> Trash
FOR i=0,nVals-1 DO  nModes[i] = tempArr[i] -> avg('r')
FOR i=0,nVals-1 DO tempArr[i] -> Trash
RETURN, nModes
END



FUNCTION alt_nModes, dataStr, nVals=N, nMax=n_max
nVals=8
IF(N_ELEMENTS(n) EQ 1) THEN nVals=N

nMax=4088
IF(N_ELEMENTS(n_max) EQ 1) THEN nMax=4088

nModes = OBJARR(nVals)
tempArr = OBJARR(nVals)
FOR i=0,nVals-1 DO  nModes[i] = dataStr.potential -> slice(n=i*nMax/nVals)
FOR i=0,nVals-1 DO tempArr[i] = nModes[i] -> AbsSq()
FOR i=0,nVals-1 DO  nModes[i] -> Trash
FOR i=0,nVals-1 DO  nModes[i] = tempArr[i] -> avg('r')
FOR i=0,nVals-1 DO tempArr[i] -> Trash
RETURN, nModes
END


FUNCTION GYRO_ETG, GyroData, _Extra=Extra
;
; Data analysis protocol for GYRO ETG benchmark and 
; related magnetic shear scan.
;
; Argument:
;
;		GyroData is the structure which NetCDF_DATA 
;		returns after reading the GYRO NetCDF data
;		file. (Required)
;
; Keywords:
;
;	RunID		Text string to appear in the lower
;			right-hand corner (upper line) of plots.
;			Defaults to no alteration (Optional).
;
;
;	FileID		Text string to appear in the lower
;			right-hand corner (lower line) of plots.
;			Defaults to no alteration (Optional).
;
;	AspectRatio	The aspect ratio from the GYRO input file
;			(that is, R/a, NOT R/r).
;			Defaults to 2.75 (optional)
;
;	RhoStar		That is, rho_e/a (if a numbe greater 
;			than 1 is entered, will assume that 
;			you have entered a/rho_e).
;			Defaults to 6e-5. (Optional).
;
;	q_safety	The safety factor.
;			Defaults to 1.4 (Optional)
;
;	L_T		The temperature gradient scale length
;			in units of 'a' (as per the GYRO
;			input file).
;			Defaults to 0.3986 (Optional).
;
;	tRange		Time interval over which correlation functions
;			are to be computed.  
;			Defaults to entire interval (Optional).
;
;	nModes		Number of toroidal modes to be retained
;			in mode plots.  
;			Defaults to 5 (Optional).
;
; Written by W.M. Nevins
;	10/4/05
;
output = {Name:"Gyro_ETG", RunDate:GyroData.RunDate}

result = GetKeyWord("RunID", Extra)
IF(Query_String(result)) THEN BEGIN
	RunID=Result
	nTags = N_TAGS(GyroData)
	FOR i=4, nTags-1 DO GyroData.(i) -> set, RunID=RunID
	output = CREATE_STRUCT(output, "RunID", RunID)
ENDIF
;
result = GetKeyWord("FileID", Extra)
IF(Query_String(result)) THEN BEGIN
	FileID=Result
	nTags = N_TAGS(GyroData)
	FOR i=4, nTags-1 DO GyroData.(i) -> set, FileID=FileID
	output = CREATE_STRUCT(output, "FileID", FileID)
END

AspectRatio = 2.75
result = GetKeyWord("AspectRatio", Extra)
IF(Query_Real(result)+Query_Integer(result)) THEN AspectRatio = result > 1./result
output = CREATE_STRUCT(output, "AspectRatio", AspectRatio)

RhoStar = 6.e-5
result = GetKeyWord("RhoStar", Extra)
IF(Query_Real(result)+Query_Integer(result)) THEN RhoStar = result < 1./result
output = CREATE_STRUCT(output, "RhoStar", RhoStar)

q_safety = 1.4
result = GetKeyWord("q_safety", Extra)
IF(Query_Real(result)+Query_Integer(result)) THEN q_safety=result
output = CREATE_STRUCT(output, "q_safety", q_safety)

L_T = 0.3986
result = GetKeyWord("L_T", Extra)
IF(Query_Real(result)+Query_Integer(result)) THEN L_T=result
output = CREATE_STRUCT(output, "L_T", L_T)

result = GetKeyWord("tRange", Extra)
IF(N_ELEMENTS(result) EQ 2) THEN BEGIN
	tRange=result
	output = CREATE_STRUCT(output, "tRange", tRange)
ENDIF

nModes = 5
result = GetKeyWord("nModes", Extra)
IF(Query_Real(result)+Query_Integer(result)) THEN nModes = result
output = CREATE_STRUCT(output, "nModes", nModes)
;
; Compute radial average of Chi_e
;
temp = GyroData.Diff_t -> Avg('r')
Chi_e = temp -> times(L_T)
temp -> trash
Chi_e -> Set, mnemonic="Chi_e", units="(!4q!X!De!N/L!DT!N)!4q!3!De!Nv!Dte!N"
Chi_e -> ScaleAxis, 't', const=1./L_T, units="L!DT!N/v!Dte!N"
output = CREATE_STRUCT(output, "Chi_e", Chi_e)
;
; Get potential vs (r, N, t).
; Scale r to units of rho_e
; Scale N to units of k_perp*rho_e  
; Scale amplitude to units of (rho_e/L_T)(T/e) 
;
Phi_r = GyroData.Potential_r -> Over(RhoStar/L_T)
Phi_r -> Set, mnemonic="Phi", units="(!4q!X!De!N/L!DT!N)(T/e)"
phi_r -> Set, axis=1, boundary="periodic (open)"

Bratio = SQRT(1.+(2.*AspectRatio*q_safety)^2)
Phi_r -> Get, axis=1, range=range
Phi_r -> ScaleAxis, 'r', const=1./RhoStar, offset=-range[0]/rhoStar, units="!4q!X!De!N"
Phi_r -> ScaleAxis, "N", const=Bratio*rhoStar/AspectRatio, title="k!D!9x!X!N", units='1/!4q!X!De!N'
Phi_r -> ScaleAxis, "t", const=1./L_T, units="L!DT!N/v!Dte!N"
Phi_r -> Get, axis=2, irange=irange
nMax = irange[1]
nMax = nMax < nModes
Intensity = OBJARR(nMax+1)
rCorrs    = OBJARR(nMax+1)
rCorrsMax = OBJARR(nMax+1)
FOR i=0, nMax DO BEGIN
	phi_kperp = Phi_r -> slice(axis=2, index=i)
	temp = phi_kperp -> AbsSQ()
	Intensity[i] = temp -> Avg('r')
	temp -> trash
	IF(N_ELEMENTS(tRange) EQ 2) THEN BEGIN
		phi_kperp -> ScaleAxis, 2, range=trange
		phi_kperp -> Restrict
	ENDIF
	temp = Phi_Kperp -> Xcorr(/Norm)
	rCorrs[i]    = temp -> slice(tau=0)
	rCorrsMax[i] = temp -> slice(tau="max")
	rCorrs[i] -> Get, vRange=vRange
	vRange[0] = 0 < vRange[0]
	rCorrs[i] -> Set, vRange=vRange
	rCorrsMax[i] -> Get, vRange=vRange
	vRange[0] = 0 < vRange[0]
	rCorrsMax[i] -> Set, vRange=vRange
ENDFOR
output = CREATE_STRUCT(output, 	"IntensityArr", Intensity, 	$
				"rCorrArr", rCorrs, 		$
				"rCorrMaxArr", rCorrsMax	)

;
; Get potential vs (theta, N, t).
; Scale N to units of k_perp*rho_e  
; Scale amplitude to units of (rho_e/L_T)(T/e) 
;
Phi_theta = GyroData.Potential_theta -> Over(RhoStar/L_T)
Phi_theta -> Set, title="!4u!X", mnemonic="Phi", units="(!4q!X!De!N/L!DT!N)(T/e)"
phi_theta -> Set, axis=1, gridtitle="!4h!X", gridUnits="", boundary="open"

Phi_theta -> ScaleAxis, "N", const=Bratio*rhoStar/AspectRatio, title="k!D!9x!X!N", units='1/!4q!X!De!N'
Phi_theta -> ScaleAxis, "t", const=1./L_T, units="L!DT!N/v!Dte!N"
Phi_theta -> Get, axis=2, irange=irange
nMax = irange[1]
nMax = nMax < nModes
I_vs_theta = OBJARR(nMax+1)
thetaCorrs    = OBJARR(nMax+1)
thetaCorrsMax = OBJARR(nMax+1)
FOR i=0, nMax DO BEGIN
	phi_kperp = Phi_theta -> slice(axis=2, index=i)
	IF(N_ELEMENTS(tRange) EQ 2) THEN BEGIN
		phi_kperp -> ScaleAxis, 2, range=trange
		phi_kperp -> Restrict
	ENDIF
	temp = phi_kperp -> AbsSQ()
	I_vs_theta[i] = temp -> Avg('t')
	temp -> trash
	I_vs_theta[i] -> MaxNorm
	I_vs_theta[i] -> Sqrt 
	ref = phi_kperp -> slice(theta=0.)
	temp = Phi_Kperp -> Xcorr(Ref=ref)
	temp -> MaxNorm
	thetaCorrs[i]    = temp -> slice(tau=0)
	thetaCorrsMax[i] = temp -> slice(tau="max")
	thetaCorrs[i] -> Get, vRange=vRange
	vRange[0] = 0 < vRange[0]
	thetaCorrs[i] -> Set, vRange=vRange
	thetaCorrsMax[i] -> Get, vRange=vRange
	vRange[0] = 0 < vRange[0]
	thetaCorrsMax[i] -> Set, vRange=vRange
ENDFOR
output = CREATE_STRUCT(output, 	"I_vs_thetaArr", I_vs_theta, 	$
				"thetaCorrArr", thetaCorrs, 	$
				"thetaCorrMaxArr", thetaCorrsMax)


RETURN, output
END  ; ****** GYRO_ETG ****** ; 

FUNCTION PG3EQ_H_ETG, pg3eq_Hstr, _Extra=Extra
;
; Fixes normalization on Chi, wSq.  Computes 
; expected noise level for chi, intensity, etc.
;
; Argument:
;
;	Argument is the data structure returned by 
;	PG3EQ_Data upon reading the *.h.nc file.
;	You MUST select ELS, ETOT, and EFLUXI
;	when you read this file.
;
; Keywords:
;
;	RunID		Text string to appear in the lower
;			right-hand corner (upper line) of plots.
;			Defaults to no alteration (Optional).
;
;
;	FileID		Text string to appear in the lower
;			right-hand corner (lower line) of plots.
;			Defaults to no alteration (Optional).
;
; Keywords passed to GKVs1D::WhiteNoise
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
;	GradT = GradT	The temperature gradient, R/L_T.  Used in computation
;			of the noise diffusion coefficient.  
;			Defaults to 6.92.
;
;
; 
; Written by W.M. Nevins
;	10/6/05
;
result = GetKeyWord("RunID", Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE "undefined") THEN BEGIN
		RunID=result
		pg3eq_hstr.efluxi -> set, runID=runID
		pg3eq_hstr.etot   -> set, runID=runID
		pg3eq_hstr.els    -> set, runID=runID
	ENDIF
ENDIF

result = GetKeyWord("FileID", Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE "undefined") THEN BEGIN
		fileID=result
		pg3eq_hstr.efluxi -> set, fileID=fileID
		pg3eq_hstr.etot   -> set, fileID=fileID
		pg3eq_hstr.els    -> set, fileID=fileID
	ENDIF
ENDIF
;
; Get Chi_e
;
Chi_e = pg3eq_hstr.efluxi -> times(1.5)
Chi_e -> set, 	title="!4v!X!De!N", 	$
		mnemonic="Chi_e", 	$
		units="(!4q!X!De!N/L!DT!N)!4q!X!De!Nv!Dte!N"
;
; Get wSq
;
wSq = pg3eq_hstr.etot -> times(2.)
wSq -> set,	title="!13<!Xw!U2!N!13>!X", 	$
		mnemonic="Wsq",			$
		units=""
;
; Get fluctuation intensity
;
Iphi = pg3eq_hstr.els
Iphi -> set,	title="!13<!4u!X!U2!N!13>!X",	$
		mnemonic="Iphi",		$
		units="(((!4q!X!De!N/L!DT!N)(T/e))!U2!N"
;
; Prepare output structure
;
output = {	Name	:	"pg3eq_h_etg",	$
		Chi_e	:	Chi_e,		$
		wSq	:	wSq,		$
		Iphi	:	Iphi		}
;
; Compute estimates of discrete particle noise 
; in this pg3eq simulation
;
NoiseStr = wsq -> whiteNoise(_Extra=Extra)
;
; Update output structure
;
output = CREATE_STRUCT(output, "NoiseStr", NoiseStr)


RETURN, output
END ; ****** FUNCTION PG3EQ_H_ETG ****** ;
