FUNCTION GKV_NormedCorrs, obj, index
;
; script for computing normed cross-correlation function
;
objstr=obj -> delta(axis='t')
refObj = objstr.delta -> slice(theta=index)
dobjCorrs = objstr.delta -> xcorr(ref=refObj)
refStats = refObj -> stats()
norm = refStats.std^2
result = dObjCorrs -> over(norm)
;
; clean up
;
dObjCorrs -> trash
refobj -> trash
gkvdelete, objstr
;
; and we're done
;
RETURN, result
END

