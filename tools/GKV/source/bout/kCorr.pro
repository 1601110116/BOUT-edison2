FUNCTION GKVs3D::kCorr
;
;
; assumes 'self' is phi(theta, zeta, t)
;
phi = self
phi -> set, axis=2, boundary='periodic (closed)'
phi_n = phi -> fft('zeta')
phi_n -> get, axis=2, range=nRange, irange=irange
nMax = nRange[1]
nPoints = irange[1]/2-1
phi_n_arr = objarr(nPoints)
for i=0,nPoints-1 do phi_n_arr[i]=phi_n -> slice(k_zeta=(i+1)*nMax/nPoints)
phi_n -> trash
ref_n_arr = objarr(nPoints)
for i=0,nPoints-1 do ref_n_arr[i] = phi_n_arr[i] -> slice(theta=38.)
kCorr_n_arr = objarr(nPoints)
for i=0, nPoints-1 do kCorr_n_arr[i] = phi_n_arr[i] -> xcorr(ref=ref_n_arr[i])
gkvdelete, ref_n_arr
gkvdelete, phi_n_arr
temp = objarr(nPoints)
for i=0, nPoints-1 DO temp[i] = kCorr_n_arr[i] -> slice(tau='max')
refarr = objarr(nPoints)
for i=0, nPoints-1 do refarr[i] = temp[i] -> slice(theta=38.)
result = objarr(nPoints)
for i=0, nPoints-1 do result[i] = temp[i] -> over(refarr[i])
gkvdelete, temp
gkvdelete, refarr
return, result
end

