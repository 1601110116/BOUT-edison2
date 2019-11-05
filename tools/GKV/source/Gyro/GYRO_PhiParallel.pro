Function GYRO_PhiParallel, 	gyroStr, _Extra=Extra
;
;  Protocol for analyzing the parallel structure of phi
;
; Written by W.M. Nevins
;   3/2/06
;
nArgs = N_PARAMS(0)
FORWARD_FUNCTION GKVsd_GyroData
CASE nArgs OF
	0:	thisData=GKVsd_GyroData(         _Extra=Extra)
	1:	thisData=GKVsd_GyroData(GyroStr, _Extra=Extra)
ENDCASE
phi_theta = thisData.phi_theta
phi_theta -> get, axis=2, Irange=irange
n_theta = 1 + irange[1] - irange[0]
Intensity = OBJARR(n_Theta)
Intensity_t=OBJARR(n_Theta)
I_even    = OBJARR(n_Theta)
I_odd     = OBJARR(n_Theta)
I_even_t  = OBJARR(n_Theta)
I_odd_t   = OBJARR(n_Theta)
CorrEven  = OBJARR(n_Theta)
CorrOdd   = OBJARR(n_Theta)
CorrEven_tau  = OBJARR(n_Theta)
CorrOdd_tau   = OBJARR(n_Theta)
SpectEven = OBJARR(n_Theta)
SpectOdd  = OBJARR(n_Theta)

FOR i=irange[0], irange[1] DO BEGIN
	thisPhi = phi_theta -> slice(axis=2, index=i)
	temp = thisPhi -> AbsSq()
	Intensity[i] = temp -> avg('t')
	Intensity_t[i] = temp -> avg('theta')
	temp -> trash
	str = thisPhi ->  EvenOdd(theta=0.)
	phi_0 = str.even -> slice(theta=0.)
	phi_1 = str.odd  -> slice(theta=1.)
	temp  = str.even -> xcorr(ref=phi_0)
	temp1 = temp -> slice(tau=0.)
	temp2 = temp1-> slice(theta=0)
	temp1  -> Get, title=title
	CorrEven[i] = temp1 -> over(temp2)
	CorrEven[i] -> Set, title=title, units=''
	temp1 -> Trash
	temp1 = temp -> slice(theta=0.)
	temp1  -> Get, title=title
	CorrEven_tau[i] = temp1 -> Over(temp2)
	CorrEven_tau[i] -> Set, title=title, units=''
	temp  -> trash
	temp1 -> Trash
	temp2 -> Trash
	temp = Str.odd -> xcorr(ref=phi_1)
	temp1= temp -> slice(tau=0.)
	temp2= temp1-> slice(theta=1.)
	temp1  -> Get, title=title
	CorrOdd[i] = temp1 -> Over(temp2)
	CorrOdd[i] -> Set, title=title, units=''
	temp1 -> Trash
	temp1 = temp -> slice(theta=1.)
	temp1  -> Get, title=title
	CorrOdd_tau[i] = temp1 -> Over(temp2)
	CorrOdd_tau[i] -> Set, title=title, units=''
	temp  -> Trash
	temp1 -> Trash
	temp2 -> Trash
	SpectEven[i] = phi_0 -> xspect()
	SpectOdd[i]  = phi_1 -> xspect()
	temp = str.even -> AbsSq()
	I_even[i]   = temp -> avg('t')
	I_even_t[i] = temp -> avg('theta')
	temp  -> Trash
	temp = str.odd  -> AbsSq()
	I_odd[i]   = temp -> Avg('t')
	I_odd_t[i] = temp -> Avg('theta')
	temp  -> Trash
	phi_0 -> Trash
	phi_1 -> Trash
	GKVdelete, str
ENDFOR

output = {	Name		:	"GYRO_PhiParallel",		$
		GyroData	:	thisData,			$
		Intensity	:	intensity,			$
		Intensity_t	:	intensity_t,			$
		I_even		:	I_even,				$
		I_odd		:	I_odd,				$
		I_even_t	:	I_even_t,			$
		I_odd_t		:	I_odd_t,			$
		CorrEven	:	CorrEven,			$
		CorrOdd		:	CorrOdd,			$
		CorrEven_tau	:	CorrEven_tau,			$
		CorrOdd_tau	:	CorrOdd_tau,			$
		SpectEven	:	SpectEven,			$
		SpectOdd	:	SpectOdd			}

RETURN, output
END ; ****** GYRO_PhiParallel ****** ;
