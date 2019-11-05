FUNCTION GKVsd_GyroKiData, GyroStr, q=qin, aspectRatio=aspectRatioIn, 		$
		rhoStar=rhoStarIn, dr=drIn, r=rin, L_perp = Lperp,	$
		perpScale = PerpScaleIn, Mu=muIn, _Extra=Extra
;
; Routine for reading gyro flux tube data and returing info
; in units used for our ETG benchmark.
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
;			(or inverse of same).  Defaults to 1./16935.8
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

rhoStar = 16935.8
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

aspectRatio=2.775
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
; Get Diff_T_1
;
Diff_T_1 = GyroStr.Diff_T_1 -> Times(Mu*L_perp)
Diff_T_1 -> ScaleAxis, 't', const=Mu/L_perp, title="t", mnemonic="t", units=perpScale + "/v!Dte!N"
Chi_i = Diff_T_1 -> Avg('r')
Chi_i -> set, title="!4v!X!Di!N", mnemonic="Chi_i", units="(!4q!3!De!N/" + PerpScale + ")!4q!3!De!Nv!Dte!N"
;
; Get Diff_T_2
;
Diff_T_2 = GyroStr.Diff_T_2 -> Times(Mu*L_perp)
Diff_T_2 -> ScaleAxis, 't', const=Mu/L_perp, title="t", mnemonic="t", units=perpScale + "/v!Dte!N"
Chi_e = Diff_T_2 -> Avg('r')
Chi_e -> set, title="!4v!X!De!N", mnemonic="Chi_e", units="(!4q!3!De!N/" + PerpScale + ")!4q!3!De!Nv!Dte!N"
;
; Get potential vs. bi-normal and radius
;
Gyro_tags = Tag_Names(GyroStr)
nTags = N_Tags(GyroStr)
FOR i=0, nTags-1 DO BEGIN
	IF( STRCMP(Gyro_tags(i), "POTENTIAL"  ,     /FOLD_CASE) ) THEN potential = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "POTENTIAL_R"  ,   /FOLD_CASE) ) THEN potential = GyroStr.(i)
	IF( STRCMP(Gyro_tags(i), "POTENTIAL_THETA", /FOLD_CASE) ) THEN pot_theta = GyroStr.(i)

ENDFOR
IF(TypeOf(potential) EQ 11) THEN BEGIN
	phi_r = potential -> GYRO_ScalePhi(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	phi_r -> Set, axis=3, gridunits = PerpScale + "/v!Dte!N"
	phi_r -> ScaleAxis, 't', const=mu
ENDIF 
IF(TypeOF(pot_theta) EQ 11) THEN BEGIN
	phi_theta = pot_theta -> GYRO_ScalePhi_Theta(q=q, AspectRatio=AspectRatio, rhoStar=rhoStar, 	$
					L_Perp=L_Perp, perpScale=perpScale)
	phi_theta -> Set, axis=3, gridunits = PerpScale + "/v!Dte!N"
	phi_Theta -> ScaleAxis, 't', const=mu	
ENDIF 
;
IF(nArgs EQ 0) THEN GKVdelete, GyroStr
output = {	Name		:	"GyroData",	$
		Chi_i		:	Chi_i,		$
		Chi_e		:	Chi_e		}
IF(TypeOF(phi_r) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "phi_r", phi_r)
IF(TypeOF(phi_theta) EQ 11) THEN	$
	output = CREATE_STRUCT(output, "phi_theta", Phi_theta)

RETURN, output
END ; ****** GKVsd_GyroData ****** ;
