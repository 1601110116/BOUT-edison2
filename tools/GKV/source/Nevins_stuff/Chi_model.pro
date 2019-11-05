FUNCTION Chi_Model, QLconst=QL_const, beta=betaIn
;
; Take follwing data, and make a GKVs1D object out of it
;
QLconst = 0.09668
IF(N_ELEMENTS(QL_const) EQ 1) THEN QLconts = QL_const

beta = 0.011401
IF(N_ELEMENTS(betaIn) EQ 1) THEN beta = betaIn


dPhiSq = FINDGEN(201)
chi=QLconst*dPhiSQ/(1.+beta*dPhiSq)

ObjStr = {GKVs1D}
ObjStr.mnemonic = 'Chi_model'
ObjStr.Title = '!4v!X!Dmodel!N'
indices = REPLICATE('*', 1)
ObjStr.Indices = PTR_NEW(indices)
ObjStr.units = '(!4q!X!Ds!N/a)!4q!X!Ds!Nc!Ds!N'
ObjStr.values = PTR_NEW(chi)
Emin = MIN(chi, max=Emax)
ObjStr.vrange = [Emin, Emax]
;
; Define grid for independent variable, 'dPhiSq'
;
grid1 = {Grid}
grid1.Mnemonic = 'dPhiSq'
grid1.title = '!4du!X!U2!N'
grid1.units = '((a/!4q!X!Ds!N)(T/e))!U2!N'
grid1.values = PTR_NEW(dPhiSq)
grid1.boundary = 'open'
grid1.uniform = 1b

grid1.range = [0.,200.]
imin=0
imax = 200
grid1.irange = [imin, imax]
;
;
ObjStr.Grid1 = grid1
;
result = OBJ_NEW('GKVs1D', ObjStr)
;
RETURN, result
END

