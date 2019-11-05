FUNCTION GKVs2D::EvenOdd
;
; returns even and odd (vs. theta) parts of correlation functions of self
;
FORWARD_FUNCTION GKVsd_MIN
phi_0 = self -> slice(theta=0.)
phi_1 = self -> slice(theta=!PI/2.)
phi_0Corrs = self -> xcorr(ref=phi_0)
phi_1Corrs = self -> xcorr(ref=phi_1)
;
; force phi_0Corrs to be an even function
;
phi_0CorrValues = *phi_0Corrs.values
even = phi_0CorrValues
thetaRange=self.grid1.irange
jMax=thetaRange[1]
jMid=jMax/2
even[0,*]=0.
even[jMid,*]=phi_0CorrValues[jMid,*]
For i=1,jMid DO BEGIN
	ii=jMax-i+1
	even[i ,*] = 0.5*(phi_0CorrValues[i ,*]+phi_0CorrValues[ii,*])
	even[ii,*] = 0.5*(phi_0CorrValues[ii,*]+phi_0CorrValues[ i,*])
ENDFOR
eMin = GKVsd_MIN(even, MAX=eMax)
vRange = [eMIn, eMax]
evenPtr=PTR_NEW(even)
evenValues = even
even = self -> MakeCopy(/NoValues)
even.vrange=vrange
even.values = evenPtr
even1 = even -> slice(theta=!PI/2.)
even0 = even -> slice(theta=0)
even00 = even0 -> slice(tau=0)
even10=even1 -> slice(tau=0.)
norm = *even10.values/*even00.values
oddValues = *phi_1Corrs.values - norm*evenValues
odd = even -> MakeCopy(/NoValues)
odd.values=PTR_NEW(oddValues)
;
; force "odd" to be odd
;
nOddValues = 0.*oddValues
oddValues[0,*] = 0.
FOR i=1,jMid DO BEGIN
	ii=jMax-i+1
	nOddValues[ i,*] = 0.5*(oddValues[ i,*] - oddValues[ii,*])
	nOddValues[ii,*] = 0.5*(oddValues[ii,*] - oddValues[ i,*])
ENDFOR
even1 -> trash
PTR_FREE, odd.values
odd.values = PTR_NEW(nOddValues)
oMin = GKVsd_MIN(nOddValues, MAX=oMax)
odd.vrange=[oMIn,oMax]
even -> set, vrange=[oMin,eMax]
odd0 = odd -> slice(tau=0.)
even0=even -> slice(tau=0.)
result = {	Name:	"EvenOdd",	$
		Even:	even,		$
		 Odd:	 odd,		$		
		Even0:	even0,		$
		 Odd0:	 Odd0,		$
		Phi0Corrs:phi_0Corrs,	$
		Phi1Corrs:phi_1corrs	}	
return, result
END
