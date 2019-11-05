FUNCTION GKV_Random, nPoints=nn, sigma=ss, dt=ddt
;
; Purpose:
;
;	This function returns a GKVs1D object containing (uncorrelated?) random numbers.
;
;
; Input Keywords:
;
;	nPoints		Length of the signal.  Defaults to 1000.
;			(Optional)
;
;
;	Sigma		Standard Deviation of the signal.  Defaults to 1.0.
;			(Optional)
;
;	dt		Step size for the independent variabal
;
;
;  Written by W.M. Nevins
;	6/16/00
;
; First set defaults and check for keywords
;
nPoints = 1000
IF Query_Integer(nn) THEN	$
	IF(nn GT 0) THEN nPoints = nn
sigma=1.0
IF Query_Real(ss) THEN		$
	IF(ss GT 0.) THEN sigma = ss
dt=1.0
IF Query_Real(ddt) THEN		$
	IF(ddt GT 0.) THEN dt=ddt
;
; Get GKVs1D structure, and initialize fields
;
out = {GKVs1D}
out.title = "Random Numbers"
out.mnemonic = "Random"
values = RANDOMN(seed, nPoints)		; 'seed', which is RETURNED by 'RANDOMN' reflects the state of the random number generator
out.values = PTR_NEW(values)
out.vrange = [-4,4]
indices = ['*']
out.Indices = PTR_NEW(indices)
CodeName = 'GKV'
CodePI = 'W.M. Nevins'
FileID = 'Random Numbers'
grid1 = {Grid}
grid1.mnemonic = 't'
grid1.title = 't'
gridValues = dt*FINDGEN(nPoints)
grid1.values = PTR_NEW(gridValues)
tmin = MIN(gridValues, max=tmax)
grid1.boundary = 'open'
grid1.uniform = 1b
grid1.range = [tmin, tmax]
grid1.irange = [0,nPoints-1]
out.grid1 = grid1
output = OBJ_NEW('GKVs1D', out)
RETURN, output
END ; ****** Random ****** ;


	