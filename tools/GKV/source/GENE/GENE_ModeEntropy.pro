FUNCTION GENE_ModeEntropy, path=Path, CodePI=CodePIin, RunID=RunIDin, FileID=fileIDin,	$
				EigenFunctions=ef
;
; Reads the EigenVectors.dat and corresonding Parameters.dat
; files for a linear GENE run, and returns the (normalized)
; entropy vs. mode number.
;
; Keywords:
;
;	Path=Path	Path to the directory containing
;			the EigenVectors.dat file and 
;			corresponding Parameters.dat
;			file. Default is to look for
;			these files in the current
;			working directory. (Optional)
;
;	Eigenfunctions	Set this keyword to include 
;			a GKV object array containing 
;			the eigenfunctions. Only works
;			(for now) if there is only one
;			kx. (Optional)
;
;  Written by W.M. Nevin
;	1/28/2009
;
CD, CURRENT = CurrentWorkingDirectory
IF(Query_String(path)) THEN CD, path
;
; Read parameters.dat file
;
Params = GENE_ReadParameters()
;
; get dimensions of each linear eigenvector
;
n_spec = Params.box.n_spec
nx = Params.box.nx0
nz = Params.box.nz0
nv = Params.box.nv0
nw = Params.box.nw0
;
; Compute expected number of linear modes,
; and create an array to hold (normalized)
; entropy of each mode
;
nModes = n_spec*nx*nz*nv*nw
entropy = FLTARR(nModes)
modeNumber = FLTARR(nModes)
EigenFunctions = 0
IF(Query_Integer(ef)) THEN EigenFunctions = (ef GT 0)
IF(nx NE 1) THEN EigenFunctions = 0
IF(n_spec NE 1) THEN EigenFunctions = 0
IF(EigenFunctions) THEN BEGIN
	efStr = {GKVs3D}
	efStr.title = '!9!!!4d!Xf!9!!!X!U2!N'
	efStr.mnemonic = 'dfSq'
	indices = ['*','*','*']
	efStr.Indices = PTR_NEW(indices)
	efStr.units = ''
	efStr.CodeName = "GENE"
	efStr.CodePI = 'F. Jenko'
	efStr.RunID = 'Linear Run'
	efStr.FileID = ''
	IF(Query_String(CodePIin)) THEN efStr.CodePI = CodePIin
	IF(Query_String( RunIDin)) THEN efStr.RunID  = RunIDin
	IF(Query_STring(FileIDin)) THEN efStr.FileID = FileIDin
	;
	; Create grid structures
	;
	zGrid = {Grid}
	zGrid.mnemonic  = 'theta'
	zGrid.title     = '!4h!X'
	zGrid.units     = ''
	Lz		= 2.*!PI
	dz		= Lz/nz
	zValues		= -Lz/2. + dz*FINDGEN(nz)
	zGrid.values    = PTR_NEW(zValues)
	zGrid.boundary  = "open"
	zGrid.range     = [-Lz/2., (Lz/2.-dz)]
	zGrid.irange    = [0, nz-1]
	efStr.grid1 = zGrid

	vGrid = {Grid}
	vGrid.mnemonic  = 'Vppl'
	vGrid.title     = 'v!D!9#!X!N'
	vGrid.units     = 'v!Dth!N'
	lv = Params.box.Lv
	dV		= 2.*LV/(nV-1)
	Vmax		= Lv
	vValues		= -Vmax + dV*(FINDGEN(nV))
	vGrid.values    = PTR_NEW(vValues)
	vGrid.boundary  = "open"
	vGrid.range     = [-Vmax, Vmax]
	vGrid.irange    = [0, nV-1]
	efStr.grid2 = vGrid

	muGrid = {Grid}
	muGrid.mnemonic = 'Mu'
	muGrid.title    = '!4l!X'
	muGrid.units    = 'T!Dref!N/B!Dref!N'
	lMu = Params.box.Lw
;	dMu		= Lmu/(nMu-1)
;	muValues	= dMu*FINDGEN(nMu)
	muValues	= GetMuKNots(nW, Lmu)
	muGrid.values    = PTR_NEW(muValues)
	muGrid.boundary  = "closed/open"
	muGrid.range     = [0, Lmu]
	muGrid.irange    = [0, nW-1]
	efStr.grid3 = muGrid

	efObj = OBJ_NEW("GKVs3D", efStr)
	efObjs = OBJARR(nModes)
ENDIF
;
; Check Precision
;
precision = Params.info.precision
;
; Create array to hold one (unnormalized) eigenvector
; As per David Hatch e-mail of 1/28/09,
; "The order of the indices in the eigenvector
;  is [kx, z, v||, mu, species].
;
IF(STRCMP(precision, 'Double', /FOLD_CASE)) THEN BEGIN
	thisVector = DCOMPLEXARR(nx*nz*nv*nw*n_spec)
ENDIF ELSE BEGIN
	thisVector = COMPLEXARR(nx*nz*nv*nw*n_spec)
ENDELSE
;
; Open the eigenvectors.dat file
;
fileName = 'eigenvectors.dat'
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening eigenvectors.dat file:
	PRINT, "Error Message = ", !ERR_STRING
	RETURN, thisVector
ENDIF
;
; Begin reading the eivenvectors.dat file
;
numModes = 0L
thisMode = 0L
FOR i=nModes, 1, -1 DO BEGIN
	READF, ioUnit, thisMode
	READF, ioUnit, thisVector
	absThisVector = ABS(thisVector)
	absSqThisVector = FLOAT(thisVector*CONJ(thisVector))
	entropy[i-1] = TOTAL(absSqThisVector)/(TOTAL(absThisVector)^2)
	modeNumber[i-1] = thisMode
	numModes = numModes + 1L
	IF(EigenFunctions) THEN BEGIN
		efObjs[i-1] = efObj -> MakeCopy(/NoValue)
		values = REFORM(absThisVector,nz,nV,nW)
;		values = TOTAL(values,1)
		ModeIndex = STRCOMPRESS(STRING(thisMode))
		thisIndex = [modeIndex, indices]
		efObjs[i-1] -> Set, values=PTR_NEW(values), indices=PTR_NEW(thisIndex)		
	ENDIF
	IF(EOF(ioUnit)) THEN GOTO, DONE
ENDFOR
DONE	: FREE_LUN, ioUnit
;
; Create a GKV object to hold the (normalized) mode entropy
;
result = {GKVs1D}
result.mnemonic = 'entropy'
result.title = '!8S!X!DN!N'
result.Indices = PTR_NEW(['*'])
result.units = ''
result.values = PTR_NEW(entropy)
result.CodeName = "GENE"
result.CodePI = 'F. Jenko'
result.RunID = 'Linear Run'
Result.FileID = ''
IF(Query_String(CodePIin)) THEN result.CodePI = CodePIin
IF(Query_String( RunIDin)) THEN result.RunID  = RunIDin
IF(Query_STring(FileIDin)) THEN result.FileID = FileIDin
;
; Create a grid object
;
grid1 = {Grid}
grid1.mnemonic = "modeNumber"
grid1.title = "Mode Number"
grid1.units = ''
grid1.values = PTR_NEW(modeNumber)
grid1.range = [0,nModes]
grid1.irange= [nModes-numModes,nModes-1]
result.grid1 = grid1
;
; creat output object
;
entropyObj = OBJ_NEW('GKVs1D', result)
support = 1./entropy
eMax = MAX(entropy)
sMax = MAX(support)
nBins= FIX(SQRT(nModes))
vMin = MIN(support, MAX=vMax)
supportObj = entropyObj -> MakeCopy(/NoValues)
supportObj -> set, values = PTR_NEW(support), vrange=[vmin, vMax]
supportObj -> set, title = '!8A!X', mnemonic="A"
entropyPDF = entropyObj -> pdf(xMin=0., xMax=eMax, nBins=nBins)
supportPDF = supportObj -> pdf(xMIn=0, xMax=sMax, nBins=nBins, Reverse_Indices=r)
output = {	Name 	:	"GENE_ModeEntropy",	$
		Entropy	:	entropyObj,		$
		Support	:	supportObj,		$
		EntropyPDF :	EntropyPDF,		$
		SupportPDF :	SupportPDF,		$		
		nBins	:	nBins,			$
		Reverse_Indices : r			}
IF(eigenFunctions) THEN output = CREATE_STRUCT(output, "EigenFunctions", efObjs)
CD, CurrentWorkingDirectory

RETURN, output
END ; ****** GENE_ModeEntropy ****** ;


