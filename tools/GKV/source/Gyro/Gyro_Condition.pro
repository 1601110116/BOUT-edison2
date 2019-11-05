FUNCTION GKVs1D::Gyro_Condition, params=params, _Extra=extra

;
; Routine to condition gyro data from the original netCDF(GYRO) format into
; coordinates of (r,d_perp,theta,t).  
; 
; Argument:
;
;            rhoStar,l_perp,mu,r,q,aspectRatio-  manually set any of
;                                               these variables.
;                                               Overwrites defaults as
;                                               well as values found
;                                               in run.out
;
;
; Keywords:
;   
;            None
;
;
;
; written by E. Wang
;      5/14/09
;      10/5/10 Edited abs_of_lx_over_rho_s for gyro 9.3  Needs update !
;
;



result = self -> MakeCopy()
rhoStar = 2500.

L_perp = 2.7775
mu = 1
perpScale = 'R!D0!N'
r=0.5
q=1.4
aspectRatio=2.7775

IF (N_ELEMENTS(params) EQ 0) THEN params = self -> Gyro_Params()


IF (N_ELEMENTS(params)) EQ 0 THEN BEGIN
    PRINT, "No parameters given or found, using default or user supplied values"
ENDIF ELSE BEGIN
PRINT, "using parameters found in run.out"
L_perp = params.local_parameters.aspect_ratio + 0.
mu = params.local_parameters.mu_electron + 0.
r = params.local_parameters.radius + 0. 
q = params.local_parameters.safety_factor + 0.
aspectRatio = L_perp 
rhoStar = params.local_parameters.rho_star + 0.
rhoStar = 1/rhoStar
ENDELSE




;
; Parse _EXTRA to change physical parameters
;
test = GetKeyWord('rhoStar', extra)
IF (TypeOf(test) NE 7) THEN rhoStar=test
test = GetKeyWord('L_perp', extra)
IF (TypeOf(test) NE 7) THEN L_perp=test
test = GetKeyWord('mu', extra)
IF (TypeOf(test) NE 7) THEN mu=test
test = GetKeyWord('r', extra)
IF (TypeOf(test) NE 7) THEN r=test
test = GetKeyWord('q', extra)
IF (TypeOf(test) NE 7) THEN q=test
test = GetKeyWord('aspectRatio', extra)
IF (TypeOf(test) NE 7) THEN aspectRatio=test

dPerpdRho = (aspectRatio)*(rhoStar)/SQRT(1. + ((aspectRatio/r)*q)^2)
IF (L_perp EQ 1) THEN perpScale = 'a'


test = GetKeyWord('rhoStar', extra)
IF (TypeOf(test) NE 7) THEN rhoStar=test
;PRINT, "using  rhoStar = ", rhoStar

;
; set units 
;
myName = self.mnemonic
CASE myName OF
    'potential'         : nameSet = 1
    'potential_r'       : nameSet = 1
    'a_parallel'        : nameSet = 2
    'potential_full'    : nameSet = 1
    'a_parallel_full'   : nameSet = 2
    'moment_n'          : nameSet = 3
    'moment_n_full'     : nameSet = 3
    'moment_e'          : nameSet = 4
    'moment_e_full'     : nameSet = 4
    'diff_t_1'          : nameset = 5
    'diff_t_em_1'       : nameset = 6
    'diff_t_em_2'       : nameset = 7
    ELSE                :  MESSAGE, 'Unrecognized data set for conditioning', /INFORMATIONAL
ENDCASE


IF (nameSet EQ 1) THEN BEGIN
resultUnits = '!4q!X/' + perpScale + '(T/e)'
result -> set, title='!4u!X', mnemonic='phi', units=resultUnits
ENDIF

IF (nameSet EQ 2) THEN BEGIN
resultUnits = '!4q!X/' + perpScale + '(T c/e c_s)'
result -> set, title='A!D!9#!X!N', mnemonic='a_parallel', units=resultUnits
ENDIF 

IF (nameSet EQ 3) THEN BEGIN
resultUnits = '!4q!X/' + perpScale + '(1/n_e)'
result -> set, title='!4n!X', mnemonic='delta n_i', units=resultUnits
ENDIF

IF (nameSet EQ 4) THEN BEGIN
resultUnits = '!4q!X/' + perpScale + '(n_e)'
result -> set, title='!4t!X', mnemonic='delta t', units=resultUnits
ENDIF 

IF (nameSet EQ 5) THEN BEGIN
result -> set, title="!4v!X!Di!N", mnemonic="Chi", units="(!4q!3!Di!N/" + PerpScale + ")!4q!3!Di!Nv!Dti!N"
ENDIF

IF (nameSet EQ 6) THEN BEGIN
result -> set, title="!4v!X!DEMi!N", mnemonic="Chi_em_i", units="(!4q!3!Di!N/" + PerpScale + ")!4q!3!Di!Nv!Dti!N"
ENDIF

IF (nameSet EQ 7) THEN BEGIN
result -> set, title="!4v!X!DEMe!N", mnemonic="Chi_em_e", units="(!4q!3!Di!N/" + PerpScale + ")!4q!3!Di!Nv!Dti!N"
ENDIF
;
; adjust GRID 
; 

r_axis = self -> AxisNumber('r')
t_axis = self -> AxisNumber('t')
n_axis = self -> AxisNumber('n')

IF (t_axis NE 0) THEN BEGIN
vScale = 'c!Ds!N'
tUnits = perpScale + '/' + vScale
result -> scaleAxis, 't', const=1./L_perp, title='t', $
                           mnemonic='t', Units=tUnits
result -> Set, axis=taxis, gridunits = PerpScale + "/v!Dte!N"
;result -> ScaleAxis, 't', const=mu
ENDIF

IF (r_axis NE 0 ) THEN BEGIN
    CASE r_axis OF
        1 : rRange= self.grid1.range
        2 : rRange= self.grid2.range
        3 : rRange= self.grid3.range
        4 : rRange= self.grid4.range
    ENDCASE
    result -> ScaleAxis, 'r', const=rhoStar, offset=-rhoStar*rRange[0],	$
      title='r', mnemonic='r', units='!4q!X'
ENDIF

IF (n_axis NE 0) THEN BEGIN
temp = result -> all_k('n')
result -> trash
temp1 = temp -> FFT('n', /INVERSE)
temp -> trash
result = temp1 -> times(rhoStar*L_perp)
temp1 -> trash
result -> ScaleAxis, 'zeta', const=dPerpdRho, title='d!9!Dx!3!N',	$
			mnemonic='d_perp', units='!4q!3'
ENDIF

;
;  For 4D objects, rearrange grids and values so order is r, d_perp, theta, t
;

nDims = self -> NumDims()
IF (nDims EQ 4) THEN BEGIN
g1 = GKVsd_Gridcopy(result.grid1)
g2 = GKVsd_Gridcopy(result.grid2)
g3 = GKVsd_Gridcopy(result.grid3)

PTR_FREE, result.grid1.values
PTR_FREE, result.grid2.values
PTR_FREE, result.grid3.values

temp = TRANSPOSE(*result.values, [1,2,0,3])
*result.values = temp

result.grid1 = g2
result.grid2 = g3
result.grid3 = g1

;
; Close theta grid if argument /CloseTheta is given
;

test = GetKeyWord('CloseTheta', extra)

IF (TypeOf(test) NE 7) THEN BEGIN
IF (N_ELEMENTS(params)) EQ 0 THEN BEGIN
nky=16  ;  dummy parameters if no run.out found
nx=120
shat=.8
lx=59.52381
ky0=.084
ENDIF ELSE BEGIN  ; reading parameters, stored as strings so FIX or FLOAT for type conversion
nky = FIX(params.grid_dimensions.n_n)
nx = FIX(params.grid_dimensions.n_x)
shat = FLOAT(params.local_parameters.shear)
lx = FLOAT(params.central_box_size.abs_of_lx_over_rho_s)   ; EDITED FOR GYRO 9.3, need conditional ideally
ly = FLOAT(params.central_box_size.abs_of_ly_over_rho_s) 
ky0 = 2*!PI/ly
ENDELSE
PRINT, "closeTheta argument given, closing theta grid"

result -> closeTheta, ky0=ky0, lx=lx, shat=shat, nx=nx,nky=nky
ENDIF  ;close theta IF statement
ENDIF  ; 4D object IF statement

RETURN, result
END  ; **** Gyro_Condition.pro ****
