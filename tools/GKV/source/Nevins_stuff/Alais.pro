FUNCTION GKVs1D::Alais, Extra=_Extra
;
; Computes alaising correction to "gridless" estimate of the particle noise
;
dx = 0.9817
result = GetKeyWord('dx', Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN dx=result
kx = *self.grid1.values
kxSq = kx*kx
values = *self.values
n=FINDGEN(100)+1.
one_1 = kx*0. + 1.
one_2 = n*0. + 1.
kp = (2.*!PI/dx)*n

term = kxSq#one_2*( 1./(kx#one_2 + one_1#kp)^2 + 1./(kx#one_2 - one_1#kp)^2 )
correction = one_1 + TOTAL(term, 2)
newValues = values*correction

result = self -> MakeCopy(/novalues)

result.values = PTR_NEW(newValues)
vMin = min(newValues, MAX=vMax)
result.vrange = [vMin, vMax]

RETURN, result

END     ;  ** GKVs1D::Alais **  ;

