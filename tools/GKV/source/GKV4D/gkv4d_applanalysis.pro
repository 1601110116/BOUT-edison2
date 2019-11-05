FUNCTION GKV4D_ApplAnalysis, thisSet
;
; Analysis proceedure for A-parallel data in (x, y, theta, t) format,
; with particular attention to issue of magnetic stochasticity
; and island formation.
; 
;
; Pick poloidal angles at top (up) and bottom (down) of torus
;
theta_up   =  !PI/2.
theta_down = -!PI/2.
;
; pick radial locations of the rational surfaces
; for the fundamental mode in the bi-normal.
; Assuming that first rational surface is at x=0,
; the remaining rational surfaces can be computed 
; from the formulea Delta x = 1/(s*k_y)
; 
q    = 1.4
Shat = 0.786
k_y0 = 0.05
eps  = 0.18
;
x_rat = 1./(sHat*k_y0)*FINDGEN(4)
dx_rat2 = 0.5/(sHat*k_y0)
x_rat2 = x_rat + dx_rat2
;
; Form flux surface average
; 
Appl_y = thisSet -> Avg('y')
appl_y -> Set, axis="theta", boundary="periodic (open)"
ApplAvg= appl_y -> Avg('theta')
Appl_y -> Trash
;
; Compute modifications to shat
;
alpha = SQRT(1.+(eps/q)^2)
temp  = ApplAvg -> d2byd('x')
applAvg -> Get, ErrorBars=ErrorPtr
IF(PTR_VALID(ErrorPtr)) THEN PTR_FREE, ErrorPtr
dShat = temp -> Times(alpha*q)
temp -> Trash
dShat -> Set, title="!4d!Xs", mnemonic="dShat"
dShat -> Get, ErrorBars=ErrorPtr
IF(PTR_VALID(ErrorPtr)) THEN PTR_FREE, ErrorPtr
; 
; Estimate the width of the 
; islands, etc. at rational surfaces
;
iMax = 3
Iwidth = OBJARR(4)
Irot   = OBJARR(4)
Wavg   = OBJARR(4)
Wavg2  = OBJARR(4)
;
FOR i=0,imax DO BEGIN
  ref   = thisSet   -> slice(x=x_rat[i])
  ref -> Set, axis="theta", boundary="periodic (open)"
  dAppl_ky = ref -> FFT('y')
  ref -> Trash
  dAppl_ky0 = dAppl_ky -> Slice(k_y=k_y0)
  dAppl_ky -> Trash
  AbsdAppl_0 = dAppl_ky0 -> Abs()
  temp = dAppl_ky0 -> Avg(axis ='theta')
  AvgdAppl_0 = temp -> Abs()
  temp -> Trash
  dAppl_ky0 -> Trash
  temp = AbsdAppl_0 -> Times(16.*q/Shat/alpha)
  Iwidth[i] = temp -> SQRT()
  temp -> Trash
  Iwidth[i] -> Set, title="!4D!X!DIsland!N", mnemonic="Iwidth", units="!4q!X!Ds!N"
  temp = AbsdAppl_0 -> Times(0.5*q/alpha*k_y0^2)
  Irot[i] = temp -> CInt('theta')
  temp -> Trash
  Irot[i] -> Set, title="!4D!9P!X", mnemonic="rotation", units=""
  temp = AvgdAppl_0 -> Times(16.*q/Shat/alpha)
  AbsdAppl_0 -> Trash
  AvgdAppl_0 -> Trash
  Wavg[i] = temp -> SQRT()
  temp -> Trash
  Wavg[i] -> Set, title="W!Davg!N", mnemonic="Iwidth", units="!4q!X!Ds!N"
ENDFOR
;
; Now compute Wavg for secondary rational surfaces
; 
FOR i=0,imax DO BEGIN
  ref   = thisSet   -> slice(x=x_rat2[i])
  ref -> Set, axis="theta", boundary="periodic (open)"
  dAppl_ky = ref -> FFT('y')
  ref -> Trash
  dAppl_ky2 = dAppl_ky -> Slice(k_y=2.*k_y0)
  dAppl_ky -> Trash
  temp = dAppl_ky2 -> Avg(axis ='theta')
  AvgdAppl_2 = temp -> Abs()
  temp -> Trash
  dAppl_ky2 -> Trash
  temp = AvgdAppl_2 -> Times(16.*q/Shat/alpha)
  AvgdAppl_2 -> Trash
  Wavg2[i] = temp -> SQRT()
  temp -> Trash
  Wavg2[i] -> Set, title="W2!Davg!N", mnemonic="I2width", units="!4q!X!Ds!N"
ENDFOR

;
; create lower and upper 
; 2-D reference signals
; at x_rat[3] rational surface
;
A_up     = thisSet -> slice(theta=theta_up)
ref_up   = A_up    -> slice(x=x_rat[3])
A_down   = thisSet -> slice(theta=theta_down)
ref_down = A_down  -> slice(x=x_rat[3])
;
; Compute correlation function about x=x_rat[3] (top and bottom)
;
Corr_up   = A_up   -> Xcorr(ref=ref_up,   /Norm)
Corr_down = A_down -> Xcorr(ref=ref_down, /Norm)
A_up     -> Trash
ref_up   -> Trash
A_down   -> Trash
ref_down -> Trash
;
; create output structures
;
CAT = { Applavg  :   ApplAvg,  $
        dsHat    :   dshat,    $
        Iwidth   :   Iwidth,   $
        Irot     :   Irot,     $
        Wavg     :   Wavg,     $
        W2avg    :   Wavg2     }
  
AVG = { Corr_up    : Corr_up,  $
        Corr_down  : Corr_down  }
;
; Clean up
; 
;
;
; and, we're done ...
;
result = {  Name :  "ApplAnalysis",  $ 
            CAT  :  CAT,             $
            AVG : AVG                }
RETURN, Result
END  ;  ****** GKV4D_ApplAnalysis  ******  ;
