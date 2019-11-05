FUNCTION GKVs3D::GYRO_ScalePhi_theta, q=qin, aspectRatio=aspectRatioIn, 		$
			rhoStar=rhoStarIn, L_perp = Lperp,	$
			perpScale = PerpScaleIn, Electrons=electrons

;
; This function accepts Phi(theta,n,t) as read in from a netCDF file
; in the usual GYRO format (no negative values of n).  It scales the
; magnitude of the potential by rho/L_perp so that the result can be 
; directly compared to output from a fluxtube code.
;
;
; Argument of function is Phi(r,n,t) as returned from NETCDF_DATA
;
; Keywords:
;
;	q		Value of the safety factor at the center of the
;			flux tube.  Defaults to 1.4.
;
;	AspectRatio	Aspect ratio of the tokamak being simulated.
;			Defaults to 2.7775
;
;	rhoStar		Ratio of gyro radius to tokamak monor radius
;			(or inverse of same).  Defaults to 4.e-3
;			for ion simulations and 1./16935.8 for
;			electron simulations
;
;	dr		Alternatively, rhoStar may be specifiec by setting
;			dr to the radial grid spacing in units of gyro-radii.
;
;	r		Value of minor radius at center of flux tube
;			in units of a.  Defaults to 0.5.
;
;	L_perp		Perpendicular scale length (in units of the
;			minor radius, a) to be used for 
;			scaling the magnitude of Phi and the time
;			scale.  Defaults to 1.
;
;	PerpScale	ASCI string identifying (in HERSHEY vector
;			font coding) the perpendicular scale.
;			Defaults to 'a'.
;
;	Electrons	Set this keyword (i.e., put /ELECTRONS on the
;			command line) if GYRO simulation is for a
;			kinetic electron species (and rho, etc. is
;			electron gyroradius).  
;
; Written by W.M. Nevins
;	5/4/2005
;
rhoStar=4.e-3
IF(KEYWORD_SET(electrons)) THEN rhoStar = 16935.8
IF(Query_Real(rhoStarIn) + Query_Integer(rhoStarIn) ) THEN rhoStar=rhoStarIn

IF(Query_Real(drIn) + Query_Integer(drIn) ) THEN BEGIN
	r=*self.grid1.values
	dr_a=r[1]-r[0]
	rhoStar = drIn/dr_a
ENDIF

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar
print, "rhoStar = ", rhoStar ; 
;
; Note: the convention with the code that follows is that 
; the variable "rhoStar" contains a/rho.
;
r=0.5
IF(Query_Real(rin)) THEN r=rin

q=1.4
IF(Query_Real(qin)) THEN q=qin

aspectRatio=2.775
IF(N_ELEMENTS(aspectRatioIn) EQ 1) THEN aspectRatio=aspectRatioIn

L_perp = 0.4
IF(Query_Real(Lperp)+Query_Integer(Lperp)) THEN L_perp = Lperp

perpScale = '!3L!DT!N'
IF(Query_String(perpScaleIn)) THEN perpScale=perpScaleIn

dPerpdRho = (aspectRatio+r)*(rhoStar)/SQRT(1. + ((aspectRatio/r)*q)^2)

;
; Scale potential to gyrokinetic units
;
phi = self -> times(rhoStar*L_perp)
;
; set units, etc.
;
self -> Get, units=units
phiUnits = '(!4q!X/' + perpScale + ')(' + units + ')'
phi -> set, title='!4u!X', mnemonic='phi', units=phiUnits
;
; Scale time axis
;
vScale = 'c!Ds!N'
IF(KEYWORD_SET(electrons)) THEN vScale = 'v!Dte!N'
tUnits = perpScale + '/' + vScale
phi -> scaleAxis, 't', const=1./L_perp, title='t', 		$
			mnemonic='t', units=tUnits
;
; scale n axis to k_perp*rho
;
phi -> ScaleAxis, 2, const=1./dPerpdRho, title='k!9!Dx!3!N',	$
			mnemonic='k_perp', units='1/!4q!3'
;
return, phi

END ; ****** FUNCTION GYRO_ScalePhi_Theta ****** ;

FUNCTION GKVs3D::GYRO_ScalePhi, q=qin, aspectRatio=aspectRatioIn, 		$
			rhoStar=rhoStarIn, dr=drIn, r=rin, L_perp = Lperp,	$
			perpScale = PerpScaleIn, Electrons=electrons

;
; This function accepts Phi(r,n,t) as read in from a netCDF file
; in the usual GYRO format (no negative values of n).  It inverts
; the fourier transform to get phi(r,zeta,t).  It then scales the
; r-axis to units of rho and scales the zeta axis to units of rho 
; (in the bi-normal direction).  The magnitude of the potential
; is scaled by rho/L_perp so that the result can be directly compared to 
; output from other fluxtube codes.
;
;
; Argument of function is Phi(r,n,t) as returned from NETCDF_DATA
;
; Keywords:
;
;	q		Value of the safety factor at the center of the
;			flux tube.  Defaults to 1.4.
;
;	AspectRatio	Aspect ratio of the tokamak being simulated.
;			Defaults to 2.7775
;
;	rhoStar		Ratio of gyro radius to tokamak monor radius
;			(or inverse of same).  Defaults to 4.e-3
;			for ion simulations and 1./16935.8 for
;			electron simulations
;
;	dr		Alternatively, rhoStar may be specifiec by setting
;			dr to the radial grid spacing in units of gyro-radii.
;
;	r		Value of minor radius at center of flux tube
;			in units of a.  Defaults to 0.5.
;
;	L_perp		Perpendicular scale length (in units of the
;			minor radius, a) to be used for 
;			scaling the magnitude of Phi and the time
;			scale.  Defaults to 1.
;
;	PerpScale	ASCI string identifying (in HERSHEY vector
;			font coding) the perpendicular scale.
;			Defaults to 'a'.
;
;	Electrons	Set this keyword (i.e., put /ELECTRONS on the
;			command line) if GYRO simulation is for a
;			kinetic electron species (and rho, etc. is
;			electron gyroradius).  
;
; Written by W.M. Nevins
;	5/4/2005
;
rhoStar=4.e-3
IF(KEYWORD_SET(electrons)) THEN rhoStar = 16935.8
IF(Query_Real(rhoStarIn) + Query_Integer(rhoStarIn) ) THEN rhoStar=rhoStarIn

IF(Query_Real(drIn) + Query_Integer(drIn) ) THEN BEGIN
	r=*self.grid1.values
	dr_a=r[1]-r[0]
	rhoStar = drIn/dr_a
ENDIF

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar
print, "rhoStar = ", rhoStar ; , "dr = ", dr_a
;
; Note: the convention with the code that follows is that 
; the variable "rhoStar" contains a/rho.
;
r=0.5
IF(Query_Real(rin)) THEN r=rin

q=1.4
IF(Query_Real(qin)) THEN q=qin

aspectRatio=2.775
IF(N_ELEMENTS(aspectRatioIn) EQ 1) THEN aspectRatio=aspectRatioIn

L_perp = 0.4
IF(Query_Real(Lperp)+Query_Integer(Lperp)) THEN L_perp = Lperp

perpScale = '!3L!DT!N'
IF(Query_String(perpScaleIn)) THEN perpScale=perpScaleIn

dPerpdRho = (aspectRatio)*(rhoStar)/SQRT(1. + ((aspectRatio/r)*q)^2)

;
; Transform potential to (r,zeta,t)-space

temp = self -> all_k('n')
temp1 = temp -> FFT('n', /INVERSE)
temp -> trash
;
; Scale potential to gyrokinetic units
;
phi = temp1 -> times(rhoStar*L_perp)
;
; set units, etc.
;
self -> Get, units=units
phiUnits = '(!4q!X/' + perpScale + ')(' + units + ')'
phi -> set, title='!4u!X', mnemonic='phi', units=phiUnits
;
; Scale time axis
;
vScale = 'c!Ds!N'
IF(KEYWORD_SET(electrons)) THEN vScale = 'v!DT!N'
tUnits = perpScale + '/' + vScale
phi -> scaleAxis, 't', const=1./L_perp, title='t', 		$
			mnemonic='t', units=tUnits
;
; scale zeta axis to rho in the bi-normal direction
;
phi -> ScaleAxis, 2, const=dPerpdRho, title='d!9!Dx!3!N',	$
			mnemonic='d_perp', units='!4q!3'

rRange=self.grid1.range
phi -> ScaleAxis, 1, const=rhoStar, offset=-rhoStar*rRange[0],	$
			title='r', mnemonic='r', units='!4q!X'
;
return, phi

END ; ****** FUNCTION GYRO_ScalePhi ****** ;



FUNCTION GKVsd_Gyro_ITG_Data, GyroStr, q=qin, aspectRatio=aspectRatioIn, 		$
		rhoStar=rhoStarIn, dr=drIn, r=rin, L_perp = Lperp,	$
		perpScale = PerpScaleIn, Mu=muIn, _Extra=Extra
;
; Routine for reading gyro flux tube data and returing info
; in gyrokinetic units.
;
; Argument:
;
;	gyrostr		Structure returned by NetCDF_Data when it reads 
;			gyro .ncd file.  If no such structure is provided,
;			then this routine will invoke NetCDF_Data.
;
; Keywords:
;
;
;	q		Value of the safety factor at the center of the
;			flux tube.  Defaults to 1.4.
;
;	AspectRatio	Aspect ratio of the tokamak being simulated.
;			Defaults to 2.775
;
;	rhoStar		Ratio of gyro radius to tokamak monor radius
;			(or inverse of same).  Defaults to 1./0.0004
;
;	r		Value of minor radius at center of flux tube
;			in units of a.  Defaults to 0.5.
;
;	L_perp		Perpendicular scale length (in units of the
;			minor radius, a) to be used for 
;			scaling the magnitude of Phi and the time
;			scale.  Defaults to 0.4.
;
;	Mu		Square root of ion/electron mass ratio.
;			Only used for transforming kinetic ion runs
;			Defaults to 1.0 (correct for adiabatic ions).
;			Should probably set to 20 for kinetic ions 
;			(but see variable MU_2 in GYRO input file).
;			
;
;	PerpScale	ASCI string identifying (in HERSHEY vector
;			font coding) the perpendicular scale.
;			Defaults to 'L!DT!N'.
;
;	path		path to directory containing GYRO data.  
;			Defaults to current working directory.
;
;	ncdFile		name of NetCDF file
;
;
; written by W.M. Nevins
;	12/17/2005
;
; Modified by W.M.Nevins
;	 7/19/2006
; to properly convert kinetic ion data by
; including the "Mu" keyword.
;
; Modified by W.M. Nevins
;	6/29/07
; To read ITG data with adiabatic electrons,
; Including Moment_n and Moment_E fields
;

rhoStar = 0.0004
IF(Query_Real(rhoStarIn) + Query_Integer(rhoStarIn) ) THEN rhoStar=rhoStarIn
;
; Note: the convention with the code that follows is that 
; the variable "rhoStar" contains a/rho.
;

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar
;
;Check for "MU"
;
mu=1.0
IF(Query_Real(MuIn) + Query_Integer(MuIn) ) THEN Mu=MuIn

r=0.5
IF(Query_Real(rin)) THEN r=rin

q=1.4
IF(Query_Real(qin)) THEN q=qin

aspectRatio=2.7775
IF(N_ELEMENTS(aspectRatioIn) EQ 1) THEN aspectRatio=aspectRatioIn

L_perp = 0.4
IF(Query_Real(Lperp)+Query_Integer(Lperp)) THEN L_perp = Lperp

perpScale = 'L!DT!N'
IF(Query_String(perpScaleIn)) THEN perpScale=perpScaleIn

nArgs = N_PARAMS(0)
;
; Get gyro data if necessary
;
IF(nArgs EQ 0) THEN GyroStr = NetCDF_Data(_Extra=Extra)
;
; Get Diff_t, Diff_p
;
Gyro_tags = Tag_Names(GyroStr)
nTags = N_Tags(GyroStr)
FOR i=0, nTags-1 DO BEGIN
	IF( STRCMP(Gyro_tags(i), "DIFF_T"  ,     /FOLD_CASE) ) THEN DIFF_T = GyroStr.(i) -> Times(Mu*L_perp)
	IF( STRCMP(Gyro_tags(i), "DIFF_T_1"  ,   /FOLD_CASE) ) THEN DIFF_T = GyroStr.(i) -> Times(Mu*L_perp)
	IF( STRCMP(Gyro_tags(i), "DIFF_P"  ,     /FOLD_CASE) ) THEN DIFF_P = GyroStr.(i) -> Times(Mu*L_perp)
	IF( STRCMP(Gyro_tags(i), "DIFF_P_1"  ,   /FOLD_CASE) ) THEN DIFF_P = GyroStr.(i) -> Times(Mu*L_perp)
ENDFOR
;
; Form properly normalized Chi
;
Diff_T -> ScaleAxis, 't', const=Mu/L_perp, title="t", mnemonic="t", units=perpScale + "/v!Dti!N"
Chi = Diff_t -> Avg(axis='r')
Chi -> set, title="!4v!X!Di!N", mnemonic="Chi_i", units="(!4q!3!Di!N/" + PerpScale + ")!4q!3!Di!Nv!Dti!N"
;
; Form properly normalized momentum diffusivity, D_p 
;
Diff_P -> ScaleAxis, 't', const=Mu/L_perp, title="t", mnemonic="t", units=perpScale + "/v!Dti!N"
D_p = Diff_P -> Avg(axis='r')
D_p -> Set, title="!13D!X!Dp!N", mnemonic="D_p", units="(!4q!3!Di!N/" + PerpScale + ")!4q!3!Di!Nv!Dti!N"
;
; Get potential, density moment, and energy moment vs. bi-normal and radius
;
Gyro_tags = Tag_Names(GyroStr)
nTags = N_Tags(GyroStr)
FOR i=0, nTags-1 DO BEGIN
	IF( STRCMP(Gyro_tags(i), "POTENTIAL"  ,     /FOLD_CASE) ) THEN potential = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "POTENTIAL_R"  ,   /FOLD_CASE) ) THEN potential = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "POTENTIAL_THETA", /FOLD_CASE) ) THEN pot_theta = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "MOMENT_N", /FOLD_CASE) ) THEN aDensity = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "MOMENT_E", /FOLD_CASE) ) THEN aPressure = GyroStr.(i) -> Times(2./3.)
ENDFOR
IF(TypeOf(potential) EQ 11) THEN BEGIN
	phi_r = potential -> GYRO_ScalePhi(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	phi_r -> Set, axis=3, gridunits = PerpScale + "/v!Dti!N"
	phi_r -> ScaleAxis, 't', const=mu
;	potential -> Trash
ENDIF 
IF(TypeOF(pot_theta) EQ 11) THEN BEGIN
	phi_theta = pot_theta -> GYRO_ScalePhi_Theta(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	phi_theta -> Set, axis=3, gridunits = PerpScale + "/v!Dti!N"
	phi_Theta -> ScaleAxis, 't', const=mu
;	pot_theta -> Trash	
ENDIF 
IF(TypeOF(aDensity) EQ 11) THEN BEGIN
	Density = aDensity -> GYRO_ScalePhi(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	Density -> Set, Title="!12N!X", mnemonic="n" ;, units="(!4q!X!Di!N/L!DT!N)n"
	Density -> Set, axis=3, gridunits = PerpScale + "/v!Dti!N"
	Density -> ScaleAxis, 't', const=mu
;	aDensity -> Trash	
ENDIF 
IF(TypeOF(aPressure) EQ 11) THEN BEGIN
	Pressure = aPressure -> GYRO_ScalePhi(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	Pressure -> Set, Title="!13P!X", mnemonic="P" ;, units="(!4q!X!Di!N/L!DT!N)nT"
	Pressure -> Set, axis=3, gridunits = PerpScale + "/v!Dti!N"
	Pressure -> ScaleAxis, 't', const=mu
	aPressure -> Trash	
ENDIF 
;
IF(nArgs EQ 0) THEN GKVdelete, GyroStr
output = {	Name		:	"GyroData",	$
		Chi		:	Chi,		$
		Diff_T		:	Diff_T,		$		
		D_p		:	D_p,		$
		Diff_p		:	Diff_p		}
IF(TypeOF(phi_r) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "phi_r", phi_r)
IF(TypeOF(phi_theta) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "phi_theta", Phi_theta)
IF(TypeOF(Density) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "Density", density)
IF(TypeOF(Pressure) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "Pressure", pressure)

RETURN, output
END ; ****** GKVsd_Gyro_ITG_Data ****** ;
