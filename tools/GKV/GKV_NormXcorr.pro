Function GKV_NormXCorr, object1, object2, R0=R0
;
;-Calculate normed cross-correlation for two objects
; at given R
;--------------------------------------------------

slice1 = object1 -> slice(r=r0)
slice2 = object2 -> slice(r=r0)
slice2 -> RMSnorm
objectCorrs = object1 -> xcorr(ref=slice2)
stats1 = slice1 -> stats()
result = objectCorrs -> over(stats1.std)
objectcorrs -> trash
slice1 -> trash
slice2 -> trash
return, result
end
