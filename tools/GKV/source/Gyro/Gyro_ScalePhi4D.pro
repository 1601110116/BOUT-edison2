FUNCTION GKVs4D::GYRO_ScalePhi4D, q=qin, aspectRatio=aspectRatioIn, 		$
			rhoStar=rhoStarIn, dr=drIn, r=rin, L_perp = Lperp,	$
			perpScale = PerpScaleIn, Electrons=electrons,           $
                        AparIn=AparIn

;
; This function accepts Phi(r,n,t) as read in from a netCDF file
; in the usual GYRO format (no negative values of n).  It inverts
; the fourier transform to get phi(r,zeta,t).  It then scales the
; r-axis to units of rho and scales the zeta axis to units of rho 
; (in the bi-normal direction).  The magnitude of the potential
; is scaled by rho/L_perp so that therResult can be directly compared to 
; output from a fluxtube code.
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
;			minor radius, a) to be used 
;for 
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
;       AparIn          Set this keyword to 1 if scaling GYRO
;                       a_parallel data instead of phi data.
;
;
; Written by E. Wang
;	5/6/2009
;

rhoStar=4.e-3
IF(KEYWORD_SET(electrons)) THEN rhoStar = 16935.8
IF(Query_Real(rhoStarIn) + Query_Integer(rhoStarIn) ) THEN rhoStar=rhoStarIn

IF(Query_Real(drIn) + Query_Integer(drIn) ) THEN BEGIN
	r=*self.grid2.values
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

dPerpdRho = (aspectRatio+r)*(rhoStar)/SQRT(1. + ((aspectRatio/r)*q)^2)
;
; Transform potential to (r,zeta,t)-space
;
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
phiUnits = '!4q!X/' + perpScale + '(T/e)'
phi -> set, title='!4u!X', mnemonic='phi', units=phiUnits
IF(AparIn EQ 1) THEN BEGIN  ; a_parallel
phiUnits = '!4q!X/' + perpScale + '(T c/e c_s)'
phi -> set, title='!4u!X', mnemonic='a_parallel', units=phiUnits
ENDIF
IF(AparIn EQ 2) THEN BEGIN  ; ion density
phiUnits = '!4q!X/' + perpScale + '(1/n_e)'
phi -> set, title='!4u!X', mnemonic='delta n_i', units=phiUnits
ENDIF
IF(AparIn EQ 3) THEN BEGIN  ; electron density
phiUnits = '!4q!X/' + perpScale + '(n_e)'
phi -> set, title='!4u!X', mnemonic='delta n_e', units=phiUnits
ENDIF
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
phi -> ScaleAxis, 'zeta', const=dPerpdRho, title='d!9!Dx!3!N',	$
			mnemonic='d_perp', units='!4q!3'

rRange=self.grid2.range
phi -> ScaleAxis, 'r', const=rhoStar, offset=-rhoStar*rRange[0],	$
			title='r', mnemonic='r', units='!4q!X'
;
; Rearrange grids and values so that order is r,d_perp,theta,t
;

g1 = GKVsd_Gridcopy(phi.grid1)
g2 = GKVsd_Gridcopy(phi.grid2)
g3 = GKVsd_Gridcopy(phi.grid3)


temp = TRANSPOSE(*phi.values, [1,2,0,3])
*phi.values = temp
;PTR_FREE, temp

;phi.grid1 -> TRASH
;phi.grid2 -> TRASH
;phi.grid3 -> TRASH

phi.grid1 = g2
phi.grid2 = g3
phi.grid3 = g1

return, phi

END ; ****** FUNCTION GYRO_ScalePhi4D ****** ;



 
