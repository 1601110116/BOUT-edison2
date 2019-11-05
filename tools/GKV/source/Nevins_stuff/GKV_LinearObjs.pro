FUNCTION GKV_LinearObjs, Structure, _Extra=Extra
;
; Creates GKV objects containing the real frequency and growth
; rates vs. D from the data in "Structure"
;
; Written by W.M. Nevins
;	6/25/05
;
ky=0.1
result = GetKeyWord('k_y', Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN ky=result

GradT=6.92
result = GetKeyWord('GradT', Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN GradT=result

GKVstr = {GKVs1D}
GKVstr.mnemonic	= "gamma"
GKVstr.title	= "!4c!X"
Indices = STRARR(2)
Indices[0] = "k!Dy!N!4q!X!De!N=" + STRTRIM(STRING(ky, FORMAT='(F4.2)'),2)
Indices[1] = "*"
GKVstr.Indices	= PTR_NEW(Indices)
GKVstr.units	= "v!Dte!N/L!DT!N"
values = Structure.gamma/GradT
GKVstr.values	=  PTR_NEW(values)
vMin = MIN(values, MAX=vMax)
vrange		= [0, vMax]
CodeName	= "GS2"
CodePI		= "Bill Dorland"
RunID		= "Linear Solve"
FileID		= "ETG w/ Diffusion"

dGrid = {Grid}
dGrid.mnemonic	= "D"
dGrid.title	= "D"
dGrid.units	= "(!4q!X!De!N/L!DT!N)(!4q!X!De!Nv!Dte!N"
dValues = Structure.D/GradT
dGrid.values	= PTR_NEW(dValues)
dMin = MIN(dValues, MAX=dMax)
dGrid.range	= [dMin, dMax]
dGrid.irange	= [0, N_ELEMENTS(dValues)-1]

GKVstr.grid1 = dGrid

gammaObj = OBJ_NEW("GKVs1D", GKVstr)

wValues = Structure.omega/GradT
wMin=MIN(wValues, MAX=wMax)
omegaObj = gammaObj -> MakeCopy(/NoValues)
omegaObj -> Set,	title = "!4x!X",		$
			mnemonic = "omega",		$
			values = PTR_NEW(wValues),	$
			vrange = [0., wMax]

result = {	NAME	:	"LinearDispersion " + Indices[0],	$
		gamma	:	gammaObj,				$
		omega	:	omegaObj				}

return, result
END
				
