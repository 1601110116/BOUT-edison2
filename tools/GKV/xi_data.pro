FUNCTION Xi_Data
;
; This function reads Xi's data files and returns GKV objects
;
; Written by W.M. Nevins
;   11/28/2012
;
CD, CURRENT=CurrentWorkingDirectory
XiPath = DIALOG_PICKFILE(/DIRECTORY)
CD, XiPath
;
; Get file names
; 
XiFiles = FILE_SEARCH(XiPath, '*.dat',COUNT=nFiles)
;
; Restore Xi's data
;
IF(nFiles GT 0) THEN BEGIN
FOR i=0, nFIles-1 DO RESTORE, XiFiles[i], /VERBOSE ; second file has some problem ... must fix this!!!
ENDIF
;
; Construct 'x' and 'y' grids
;
xValues = REFORM(TOTAL(G.dx[*,0], /CUMULATIVE))
yValues = REFORM(TOTAL(G.dy[0,*], /CUMULATIVE))
Nx = N_ELEMENTS(xValues)
Ny = N_ELEMENTS(yValues)
xGrid = {grid}
xGrid.title = 'x'
xGrid.mnemonic = 'x'
xGrid.units = '' ; need units!!!
xGrid.boundary = "open"
xGrid.irange = [0,Nx-1]
xMax = MAX(xValues, MIN=xMin)
xGrid.range = [xMin, xMax]
xGrid.values = PTR_NEW(xValues)

yGrid = {grid}
yGrid.title = 'y'
yGrid.mnemonic = 'y'
yGrid.units = '' ; need units!!!
yGrid.boundary = "periodic (open)"
yGrid.irange = [0,Ny-1]
yMax = MAX(yValues, MIN=yMin)
yGrid.range = [yMin, yMax]
yGrid.values = PTR_NEW(yValues)

info=size(P_T500_2F_XU3)
Nz = info[3]
Nt = info[4]

zGrid = {grid}
zGrid.title='z'
zGrid.mnemonic='z'
zGrid.boundary = 'periodic (open)'
zGrid.irange = [0,Nz-1]
zGrid.range = [0., 2.*!PI/5.]
zValues = (2.*!PI/5.)*FINDGEN(Nz)/Nz
zGrid.values = PTR_NEW(zValues)

tGrid = {grid}
tGrid.title='t'
tGrid.mnemonic='t'
tGrid.boundary = 'open'
tGrid.irange = [0,Nt-1]
tGrid.range = [0., Nt-1.]
tValues = FINDGEN(Nt)
tGrid.values = PTR_NEW(tValues)


;
; Make a few objects ...
;
ni0Values = G.ni0[*,0]
ni0Max = MAX(ni0Values)
ni0_str = {GKVs1D}
ni0_str.title = "n!Dio!N"
ni0_str.mnemonic = "ni0"
ni0_str.CodeName = "BOUT++"
ni0_str.codePi = "Xi"
ni0_str.vRange = [0,ni0Max]
ni0_str.values= PTR_NEW(ni0Values)
ni0_str.indices = PTR_NEW(["*"])
ni0_str.Grid1 = GKVsd_GridCopy(xGrid)
ni0 = OBJ_NEW('GKVs1D', ni0_str)

Ti0Values = G.Ti0[*,0]
Ti0Max = MAX(Ti0Values)
Ti0_str = {GKVs1D}
Ti0_str.title = "T!Dio!N"
Ti0_str.mnemonic = "Ti0"
Ti0_str.codePi = "Xi"
Ti0_str.CodeName = "BOUT++"
Ti0_str.vRange = [0,Ti0Max]
Ti0_str.values = PTR_NEW(Ti0Values)
Ti0_str.indices = PTR_NEW(["*"])
Ti0_str.Grid1 = GKVsd_GridCopy(xGrid)
Ti0 = OBJ_NEW('GKVs1D', Ti0_str)


P0Values = G.pressure[*,0]
P0Max = MAX(P0Values)
P0_str = {GKVs1D}
P0_str.title = "P!D0!N"
P0_str.mnemonic = "P_0"
P0_str.codePi = "Xu Xi"
P0_str.CodeName = "BOUT++"
P0_str.vRange = [0,P0Max]
P0_str.values = PTR_NEW(P0Values)
P0_str.indices = PTR_NEW(["*"])
P0_str.Grid1 = GKVsd_GridCopy(xGrid)
P0 = OBJ_NEW('GKVs1D', P0_str)

P_values = P_T500_2F_XU3
P_Max = MAX(P_values, MIN=P_min)
P_str = {GKVs4D}
P_str.title = "P"
P_str.mnemonic = "P"
P_str.codePi = "Xu Xi"
P_str.CodeName = "BOUT++"
P_str.vRange = [P_min,P_Max]
P_str.values = PTR_NEW(P_values)
P_str.indices = PTR_NEW(["*","*","*","*"])
P_str.Grid1 = GKVsd_GridCopy(xGrid)
P_str.Grid2= GKVSD_GridCopy(yGrid)
P_str.Grid3 = GKVsd_Gridcopy(zGrid)
P_str.Grid4 = GKVSd_GridCopy(tGrid)
P = OBJ_NEW('GKVs4D', P_str)


result = {name:"xiData", ni0:ni0, Ti0:Ti0, P0:P0, P:P}
CD, CurrentWorkingDirectory
RETURN, result
END
