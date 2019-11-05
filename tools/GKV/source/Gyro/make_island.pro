function make_island, obj
;
; Takes apar data from gyro, returns 1d object of island width in units of rho_s
;

a = obj -> gyro_condition(l_perp=1,/closetheta)
a_t = a -> avg('theta')
a_tkr = a_t -> fft('r')
a_tkrd = a_tkr -> fft('d_perp')
a_sq = a_tkrd -> abssq()
a_sq -> set, errorbars=PTR_NEW()
a_sq_r = a_sq -> int('k_r', /sum)
a_sq_rt = a_sq_r -> avg('t')

a_1 = a_sq_rt -> sqrt()
a_sqrt = a_1 -> sqrt()

result = a_sqrt -> times(4*sqrt(1.4/.786*2.8))


RETURN, result
END

FUNCTION printQ, obj

q1 = obj.flux_t_1 -> avg('r')
q1 -> signalwindow, t=[150,500]
q1 -> restrict
q1stats = q1 ->stats()

q2 = obj.flux_t_2 -> avg('r')
q2 -> signalwindow, t=[150,500]
q2 -> restrict
q2stats = q2 -> stats()

q1e = obj.flux_t_em_1 -> avg('r')
q1e -> signalwindow, t=[150,500]
q1e -> restrict
q1estats = q1e -> stats()

q2e = obj.flux_t_em_2 -> avg('r')
q2e -> signalwindow, t=[150,500]
q2e -> restrict
q2estats = q2e -> stats()

PRINT, q1stats.avg, q1stats.avgpm 
PRINT, q1estats.avg, q1estats.avgpm
PRINT, q2stats.avg, q2stats.avgpm
PRINT, q2estats.avg, q2estats.avgpm


result=0
RETURN, result
END

FUNCTION printflux, obj

arg = obj -> makecopy()
arg -> signalwindow, t=[150,500]
arg -> restrict
argstats = arg -> stats()
print, argstats.avg, argstats.avgpm

result=0
RETURN, result
END


FUNCTION makemag, obj

for i = 0,99 do obj[i,0] -> poincare_norm, 59.823
a1 = REFORM(obj[*,0],10,10)
result = dmag(a1)

RETURN, result
END


