FUNCTION GKVs1D_MakeObj, path=path, _extra=extra, userInput=userInput
;
; Purpose:
;
;	This routine reads ascii data file, and produces a
;	GKVs1D object from its contents.
;
;
FORWARD_FUNCTION GKVsd_MIN
CD, current=currentWorkingDirectory
IF(Query_String(path)) THEN CD, path
IF( KEYWORD_SET(userInput) ) THEN BEGIN
	inputvalues  = GetKeyWord('values', extra)
	inputxvalues = GetKeyWord('xValues', extra)
	datastr={ xValues:inputXvalues, values:inputValues}
	errorBars = GetKeyWord('errorBars', extra)
	IF(TypeOf(errorBars) NE 7) THEN dataStr = CREATE_STRUCT(dataStr, 'ErrorBars', errorBars)
ENDIF ELSE BEGIN
	GetGridValues=0b
	result=GetKeyWord('GetGridValues', extra)
	IF(TypeOf(result) NE 0) THEN GetGridValues = result
	dataStr=READ_ASCII(dataFile, template=ASCII_TEMPLATE(dataFile))
ENDELSE
objStr={GKVs1D}
xGrid={grid}
nPoints=N_ELEMENTS(dataStr.(0))

mnemonic='mnemonic'
result=GetKeyWord('mnemonic', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, mnemonic, PROMPT='Enter Mnemonic: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, mnemonic, PROMPT='Enter Mnemonic: '
	ENDIF ELSE BEGIN
		mnemonic = result
	ENDELSE
ENDELSE
objStr.mnemonic = mnemonic

title = ''
result=GetKeyWord('title', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, title,	PROMPT='Enter Title: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, title,	PROMPT='Enter Title: '
	ENDIF ELSE BEGIN
		title = result
	ENDELSE
ENDELSE
objStr.title = title

indices=REPLICATE('*',1)
objStr.indices=PTR_NEW(indices)

units = ''
result=GetKeyWord('units', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, units,	PROMPT='Enter Units: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, units,	PROMPT='Enter Units: '
	ENDIF ELSE BEGIN
		units = result
	ENDELSE
ENDELSE
objStr.Units = units

IF(GetGridValues EQ 0) THEN BEGIN
	objStr.values = PTR_NEW(dataStr.(1))
	vmin=GKVsd_MIN(dataStr.(1), MAX=vmax)
	objStr.vrange=[vmin,vmax]
ENDIF ELSE BEGIN
	objStr.values = PTR_NEW(dataStr.(0))
	vmin=GKVsd_MIN(dataStr.(0), MAX=vmax)
	objStr.vrange=[vmin,vmax]
ENDELSE


codename='BOUT'
result=GetKeyWord('CodeName', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, codeName,	PROMPT='Enter Code Name: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, codeName,	PROMPT='Enter Code Name: '
	ENDIF ELSE BEGIN
		IF(result NE '') THEN CodeName = result
	ENDELSE
ENDELSE
objStr.CodeName = codeName

codePI = 'X.Q. Xu'
result=GetKeyWord('CodePI', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, codePI,	PROMPT='Enter Code PI: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, codePI,	PROMPT='Enter Code PI: '
	ENDIF ELSE BEGIN
		IF(result NE '') THEN CodePI = result
	ENDELSE
ENDELSE
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

gridMnemonic = ''
result=GetKeyWord('gridMnemonic', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, gridMnemonic,	PROMPT='Enter Grid Mnemonic: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, gridMnemonic,	PROMPT='Enter Grid Mnemonic: '
	ENDIF ELSE BEGIN
		gridMnemonic = result
	ENDELSE
ENDELSE
xGrid.Mnemonic = gridMnemonic

gridTitle = ''
result=GetKeyWord('gridTitle', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, gridTitle,	PROMPT='Enter Grid Title: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, gridTitle,	PROMPT='Enter Grid Title: '
	ENDIF ELSE BEGIN
		gridtitle = result
	ENDELSE
ENDELSE
xGrid.title = gridTitle

IF(GetGridValues EQ 1)THEN BEGIN
	gridValuesStr=READ_ASCII(gridFile, template=ASCII_TEMPLATE(gridFile))
	xGrid.values = PTR_NEW(gridValuesStr.(0))
ENDIF ELSE BEGIN
	xGrid.values = PTR_NEW(dataStr.(0))
ENDELSE

gridUnits = ''
result=GetKeyWord('gridUnits', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ, gridUnits,	PROMPT='Enter Grid Units: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ, gridUnits,	PROMPT='Enter Grid Units: '
	ENDIF ELSE BEGIN
		gridUnits = result
	ENDELSE
ENDELSE
xGrid.units = gridUnits

gridBoundary = ''
result=GetKeyWord('gridBoundary', extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	READ,gridBoundary,	PROMPT='Enter Grid Boundary condition: '
ENDIF ELSE BEGIN
	IF(result EQ 'undefined') THEN BEGIN
		READ,gridBoundary,	PROMPT='Enter Grid Boundary condition: '
	ENDIF ELSE BEGIN
		gridBoundary = result
	ENDELSE
ENDELSE
xGrid.boundary = gridBoundary
 
xmin=GKVsd_MIN(*(xGrid.values), MAX=xmax)

xGrid.range = [xmin, xmax]
xGrid.irange= [0, nPoints-1]
objStr.Grid1 = xGrid

IF(N_TAGS(dataStr) EQ 3) THEN BEGIN
	objStr.ErrorBars = PTR_NEW(dataStr.(2))
ENDIF

obj = OBJ_NEW('GKVs1D', objStr)
cd, currentWorkingDirectory
RETURN, obj
END
