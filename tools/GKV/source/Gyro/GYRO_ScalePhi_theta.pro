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
phiUnits = '!4q!X/' + perpScale + '(T/e)'
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



 
