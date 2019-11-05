FUNCTION GKVs3D::EnergyInfo, p
;
; This function acts on a GKVs3D object containing
; the electrostatic potential (normalized to 
; gyrokinetic units in an x,y representation)
; and returns a structure containing information 
; about the wave energy.
;
; Argument: A GKVs3D object containing 
; "pressure" vs. x,y (really 2/3 of the 
; skinetic energy moment).  (Optional)
;
; Written by W.M. Nevins
;    7/17/07
;
nArgs = N_PARAMS()
;
; Form Fourier representation of spatial data
;
temp = self -> FFT(1)
phi_k = temp -> FFT(2)
temp -> Trash
IF(nArgs EQ 1) THEN BEGIN
	ke = p -> Times(1.5)
	temp = ke -> FFT(1)
	ke_k = temp -> FFT(2)
	ke_k -> Set, title="!13K!X!Dk!N", mnemonic="K_k"
	temp -> Trash
	ke -> Trash
ENDIF
phi_k -> GET, axis=1, gridValues=k1ptr
phi_k -> GET, axis=2, gridValues=k2ptr
phi_k -> GET, axis=3, gridValues=tPtr
k1Sq = (*k1ptr)^2
k2Sq = (*k2ptr)^2
epsilon = MIN(k2SQ, k20)
n1 = N_ELEMENTS(k1Sq)
n2 = N_ELEMENTS(k2Sq)
nt = N_ELEMENTS(*tPtr)
t1 = MAKE_ARRAY(n1, VALUE=1.)
t2 = MAKE_ARRAY(n2, VALUE=1.)
shortkSq = (k1Sq#t2 + t1#k2Sq)
ksq = MAKE_ARRAY(n1,n2,nt)
FOR i=0,nt-1 DO ksq[*,*,i] = shortkSq[*,*]
phi_k_Sq = phi_k -> AbsSq()

Phi_ky0_Sq = phi_k_sq -> Slice(axis=2, value=0.)
GradPhi_k_Sq = phi_k_Sq -> Times(kSq)
GradPHi_k_Sq -> Set, title="!9!!G!4du!X!Dk!N!9!!!X!U2!N", mnemonic="Grad_Phi_sq"
phiSqValues = *(phi_k_Sq.values)
phiSqValues[*,k20,*] = 0.
PTR_FREE, phi_k_Sq.values
Phi_k_Sq -> Set, values=PTR_NEW(phiSqValues)
Phi_k_Sq -> Set, units=""
phi_ky_Sq = Phi_k_Sq -> INT(1, /SUM)
I_phi = Phi_ky_Sq -> Int(1, /SUM)
I_ky0 = Phi_ky0_Sq -> INT(1, /SUM)
GradPhi_k_Sq -> Set, units=""
GradPhi_ky_Sq = GradPHi_k_Sq -> INT(1, /SUM)
I_GradPhi = GradPHi_ky_sq -> Int(1, /SUM)
E_phi_k = Phi_k_sq -> Plus(GradPhi_k_sq)
E_phi_k -> Set, title="!12E!d!4du!X!N", mnemonic="E_phi"
E_phi_ky = E_phi_k -> INT(1, /SUM)
E_phi = E_phi_ky -> INT(1, /SUM)

result = {	Name		:	"EnergyINfo",	$
		Phi_k		:	Phi_k,		$
		Phi_k_Sq	:	Phi_k_sq,	$
		Phi_ky_Sq	:	Phi_ky_Sq,	$
		GradPhi_k_sq	:	GradPhi_k_sq,	$
		GradPhi_ky_sq	:	GradPhi_ky_sq,	$
		E_phi_k		:	E_phi_k,	$
		E_phi_ky	:	E_phi_ky,	$
		I_phi		:	I_phi,		$
		I_GradPhi	:	I_GradPhi,	$
		E_phi		:	E_phi,		$
		Phi_ky0_sq	:	Phi_ky0_Sq,	$
		I_ky0		:	I_ky0		}
IF(nArgs EQ 1) THEN BEGIN
	E_total_k  = E_phi_k -> Plus(ke_k)
	E_total_k -> Set, title="!12E!d!Xtotal!N", mnemonic="E_total"
	E_Total_ky = E_total_k -> INT(1, /SUM)
	E_total = E_total_ky -> INT(1, /SUM)
	result = CREATE_STRUCT(result,	"ke_k", ke_k, 	$
				"E_Total_k", E_total_k,	$
				"E_Total_ky", E_Total_ky,$
				"E_Total", E_Total)
ENDIF

RETURN, result
END ; ****** GKVs3D::EnergyInfo ****** ;

 
