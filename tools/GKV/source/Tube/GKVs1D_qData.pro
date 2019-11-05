FUNCTION GKVs1D_qData, path=path, _extra=extra, userInput=userInput
;
; Purpose:
;
;	This routine reads ascii data file, and produces a
;	GKVs1D object containing q-profile.
;
;
FORWARD_FUNCTION GKVsd_MIN
CD, current=currentWorkingDirectory
IF(Query_String(path)) THEN CD, path
;
; Checke for input file name
;
inFile = "q_vs_x.txt"
;result=GetKeyWord('inFile', extra)
;IF(TypeOF(result) EQ 7) AND THEN infile=result
;
; Open input file
;
GET_LUN, lun
OPENR, lun, inFile
header = STRARR(1)
READF, lun, header
q = FLTARR(50)
READF, lun, q, FORMAT='(5E13.4)'
READF, lun, header
x = FLTARR(50)
READF, lun, x, FORMAT='(5E13.4)'
FREE_LUN, lun

objStr={GKVs1D}
xGrid={grid}
nPoints=N_ELEMENTS(x)

objStr.mnemonic = 'q'
objStr.title = 'q'
indices=REPLICATE('*',1)
objStr.indices=PTR_NEW(indices)
units = ''
objStr.Units = units
objStr.values = PTR_NEW(q)
vmin=GKVsd_MIN(q, MAX=vmax)
objStr.vrange=[vmin,vmax]

codename='GEM'
objStr.CodeName = codeName

codePI = 'S. Parker & Y. Chen'
objStr.CodePi = codePi

runID = ''
result=GetKeyWord('runID', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, runID,	PROMPT='Enter Run ID: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, runID,	PROMPT='Enter Run ID: '
	ENDIF ELSE BEGIN
		runID = result
	ENDELSE
ENDELSE
objStr.RunID = runID

fileID = ''
result=GetKeyWord('fileID', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, fileID,	PROMPT='Enter File ID:'
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, fileID,	PROMPT='Enter File ID:'
	ENDIF ELSE BEGIN
		fileID = result
	ENDELSE
ENDELSE
objStr.FileID = fileID
;

gridMnemonic = 'x'
xGrid.Mnemonic = gridMnemonic

gridTitle = 'x'
xGrid.title = gridTitle

xGrid.values = PTR_NEW(x)

gridUnits = 'r/a'
xGrid.units = gridUnits

gridBoundary = ''
xGrid.boundary = gridBoundary
 
xmin=GKVsd_MIN(*(xGrid.values), MAX=xmax)

xGrid.range = [xmin, xmax]
xGrid.irange= [0, nPoints-1]
objStr.Grid1 = xGrid

obj = OBJ_NEW('GKVs1D', objStr)
cd, currentWorkingDirectory
RETURN, obj
END
