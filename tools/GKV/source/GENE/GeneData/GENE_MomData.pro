

FUNCTION GENE_MomData, 	fileName, ParamStr=ParamStr, 	$
			CorrectMoments=CorrectMoments,	$
			path=path, _Extra=Extra
;
; This function reads the mom_xxx.dat files output
; by GENE, and returns these fields as GeneField objects.
;
; Arguments:
;
;	fileName	Name of the mom_xxx.dat file to be analyzed.
;			Defaults to "mom_ions.dat". (optional).
;
; Keywords:
;
;	Path		Path to structure containing the 
;			desired mom_xxx.dat file. Defaults
;			to current working directory.
;			(Optional)
;
;	ParamStr	Structure created by GENE_ReadParameters
;			which contains various parameters (from
;			the parameters.dat file) required
;			to read the moment files. 
;
;	Endian		The "Endianness" of the mom_xxx.dat file.
;			Defaults to value in the 'parameters.dat' 
;			file (which, in my limited experience, is
;			usually 'little'). (Optional)
;
;	Start		Length (in bytes) of header in binary files
;			Defaults to 4. (Optional)
;
;	EOR		Length (in bytes) of "end of record" mark
;			in binary files. Defaults to 8. (Optional)
;
;	L_ref		Set this to a Hershey vector font string identifying the the (macroscopic) 
;			reference length.  Defaults to "R!D0!N".  (Optional)
;
;	CodeName	Set to a string constant containing the desired 'CodeName' (ascii field appearing at the
;			penultimate line at the lower left hand corner of the plot).
;			Defaults to 'GENE'. (Optional)
;
;	CodePI		Set to a string constant containing the desired 'CodePI' (ascii field appearing at the
;			last line at the lower left hand corner of the plot).  
;			Defaults to "F. Jenko". (Optional)
;
;
;	FileID		Set to a string constant containing the desired 'FileID' (ascii field appearing at the 
;			last line in the lower right hand corner of the plot).  This value will over ride 
;			the GENE .dat filename. (Optional)
;
;	RunID		Set to a string constant containing the desired 'RunID' (ascii field appearing in the 
;			penultimate line at the lower right hand corner of the plot).  (Optional)
;
;  Written by W.M. Nevins
;	10/22/2008
; Modified From GENE_READMOM by W.M. Nevins
;	11/8/2008
; to return GeneField objects 
; (instead of GKV objects)
;
nArgs = N_PARAMS()
If(nArgs EQ 0) THEN BEGIN
	fileName="mom_ions.dat"
	MESSAGE, "Using default filename, 'mom_ions.dat'", /INFORMATIONAL
ENDIF
CD, CURRENT=CurrentWorkingDirectory
IF(Query_String(path)) THEN CD, path
ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find Moment file", /INFORMATIONAL
	CD, CurrentWorkingDirectory
	RETURN, 0
ENDIF
CD, CURRENT=momPath
;
; Get FileID, etc.
;
L_ref = "R!D0!N"
result = GetKeyWord('L_ref', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN L_ref=result
ENDIF

CodeName =  "GENE"
result = GetKeyWord('CodeName', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN CodeName=result
ENDIF

CodePI   = "F. Jenko"
result = GetKeyWord('CodePI', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN CodePI=result
ENDIF

FileID = fileName
result = GetKeyWord('FileID', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN FileID=result
ENDIF

RunID = ""
result = GetKeyWord('RunID', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN RunID=result
ENDIF
;
; Get a free I/O unit, and open mom_xxx.dat file
;
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening MOMENT file: ", fileName
	PRINT, "Error Message = ", !ERR_STRING
	CD, CurrentWorkingDirectory
	RETURN, 0
ENDIF
;
; Compute size of Moment data using data from 
; parameter structure.
;
IF(TypeOf(ParamStr) EQ 8) THEN BEGIN
	nMoments = ParamStr.INFO.N_MOMS
	nkx      = ParamStr.BOX.NX0
	nky	 = ParamStr.BOX.NKY0
	nZ	 = ParamStr.BOX.NZ0
	nt	 = ParamStr.INFO.ITIME/ParamStr.IN_OUT.ISTEP_MOM
	Endian	 = ParamStr.INFO.ENDIANNESS
	Endian	 = STRLOWCASE(Endian)
	Endian	 = STRCOMPRESS(Endian, /REMOVE_ALL)
	Precision=ParamStr.INFO.Precision
	Precision=STRLOWCASE(Precision)
	Precision=STRCOMPRESS(Precision, /REMOVE_ALL)
ENDIF ELSE BEGIN
	MESSAGE, "Invalid Parameter Structure, Returning 0", /INFORMATIONAL
	RETURN, 0
ENDELSE
;
; Some key constants garnared from close
; inspection of a sample field.dat file:
;
start = 4L	; The mom_xxx.dat file begins with a 4 byte field
EOR = 8L	; Each record is separated by an 8 byte field
;
; These may be operating system dependent,
; (but, as per note from Florian Merz
;  of 10/22/08 they are not) 
; so have allow the user to change 
; "start" and "EOR" from the command line
;
result = GetKeyWord("start", Extra)
IF(Query_Integer(result)) THEN start = result

result = GetKeyWord("EOR", Extra)
IF(Query_Integer(result)) THEN EOR = result
;
; Size of each Moment block (in bytes)
;
WordSize = 8
dataType = 5		; recall that we only read t-field
IF( STRCMP(Precision, "single") ) THEN BEGIN
	WordSize=4
	dataType=4	; recall that we only read t-field
ENDIF
nMomentBytes = LONG64(2L*WordSize*nkx*nky*nz)
;
; Compute number of bytes between time values
;
iSkip = 8 + EOR + nMoments*(nMomentBytes + EOR)
;
; Read values of t at which fields are written
;
t=DBLARR(nt)
FOR i=0,nt-1 DO	$
	t[i] = READ_BINARY(ioUnit, DATA_START=start+i*iSkip, DATA_TYPE=dataType, DATA_DIMS=1, ENDIAN=Endian)
;
; Get species name by "cracking" fileName
;
speciesName = "?"
subStrings = STRSPLIT(fileName, ".", /EXTRACT, COUNT=nSubStrings)
IF(nSubStrings GT 0) THEN	$
	subStrings = STRSPLIT(subStrings[0], "_", /EXTRACT, COUNT=nSubStrings)
IF(nSubStrings GT 1) THEN BEGIN
	subStrings  = STRCOMPRESS(subStrings[1:nSubstrings-1], /REMOVE_ALL)
	speciesName = STRJOIN(substrings, "_")
ENDIF
;
; Create GKV Object
; 
ObjStr = {GKVs4D}
ObjStr.mnemonic = "dn_"+ speciesName
ObjStr.title    = "!4d!Xn!D" + speciesName + "!N"
ObjStr.indices  = PTR_NEW(["*","*","*","*"])
ObjStr.units    = "(!4q!X!Ds!N/" + L_ref + ")n!D0!N"
ObjStr.values   = PTR_NEW()
ObjStr.vrange   = [0.,1.]
ObjStr.CodeName = CodeName
ObjStr.CodePI   = CodePI
ObjStr.RunID    = RunID
ObjStr.fileID   = FileID
;
; eigenmode mode runs don't seem to define Lx in the parameters.dat file (sigh!).
; Need a work-arround
;
boxTags = TAG_NAMES(ParamStr.BOX)
boxTags = STRCOMPRESS(boxTags, /REMOVE_ALL)
nTags = N_TAGS(ParamStr.BOX)
Lx=1.0
thisTag = STRCOMPRESS("LX", /REMOVE_ALL)
FOR i=0,nTags-1 DO BEGIN
	IF( STRCMP(boxTags[i], thisTag, /FOLD_CASE) ) THEN Lx=ParamStr.BOX.LX
ENDFOR

xGrid = {Grid}
	xGrid.mnemonic  = 'k_x'
	xGrid.title     = 'k!Dx!N'
	xGrid.units     = '1/!4q!X!Ds!N'
	dkx		= 2.*!PI/Lx
	kxMax		= dkx*nkx/2.
	kxValues	= -kxMax + dkx*FINDGEN(nkx)
	xGrid.values    = PTR_NEW(kxValues)
	xGrid.boundary  = "periodic (open)"
	xGrid.range     = [-kxMax, kxMax]
	xGrid.irange    = [0, nkx-1]
ObjStr.grid1 = xGrid

yGrid = {Grid}
	yGrid.mnemonic  = 'k_y'
	yGrid.title     = 'k!Dy!N'
	yGrid.units     = '1/!4q!X!Ds!N'
	dky		= ParamStr.BOX.KYMIN
	kyMax		= (nky-1)*dky
	kyValues	= dky*FINDGEN(nky)
	yGrid.values    = PTR_NEW(kyValues)
	yGrid.boundary  = "periodic (open)"
	yGrid.range     = [0, kyMax]
	yGrid.irange    = [0, nky-1]
ObjStr.grid2 = yGrid

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
ObjStr.grid3 = zGrid

tGrid = {Grid}
	tGrid.mnemonic  = 't'
	tGrid.title     = 't'
	tGrid.units     = L_ref + '/c!Ds!N'
	tGrid.values    = PTR_NEW(t)
	tGrid.boundary  = "open"
	tGrid.range     = [t[0], t[nt-1]]
	tGrid.irange    = [0,nt-1]
ObjStr.grid4 = tGrid

nObj = OBJ_NEW("GKVs4D", ObjStr)
nDataStr = {GeneField}
nDataStr.GeneDataDir	= momPath
nDataStr.ParamStr	= PTR_NEW(ParamStr)
nDataStr.GkvObj		= nObj
nDataStr.FieldFile	= FileName
nDataStr.nField		= nMoments
nDataStr.iField		= 1
nDataObj = OBJ_NEW("GeneField", nDataStr)
;
; Prepare output structure
;
output = {	Name	:	"GENE_MomData",		$
		n	:	nDataObj		}
;
; Make T_|| object 
;
TpplObj = nObj -> MakeCopy(/NoValues)
TpplObj -> Set,	title="!4d!XT!D!9#!X" + speciesName + "!N"
TpplObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")T!D0!N"
TpplObj -> Set,	mnemonic="Tppl_" + SpeciesName	
TpplDataObj = nDataObj -> MakeCopy()
TpplDataObj -> Set, iField=2
TpplDataObj -> Set, gkvObj = TpplObj 
output  = CREATE_STRUCT(output, "Tppl", TpplDataObj)
;
; Make T_perp object 
;
TperpObj = TpplObj -> MakeCopy(/NoValues)
TperpObj -> Set, title="!4d!XT!D!9x!X" + speciesName + "!N"
TperpObj -> Set, mnemonic="Tperp_" + speciesName	
TperpDataObj = nDataObj -> MakeCopy()
TperpDataObj -> Set, iField=3
TperpDataObj -> Set, gkvObj = TperpObj 
output  = CREATE_STRUCT(output, "Tperp", TperpDataObj)
;
; Make a Q_|| object 
;
QpplObj =  nObj -> MakeCopy(/NoValues)
QpplObj -> Set,	title="Q!D!9#!X" + speciesName + "!N"
QpplObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")c!Ds!NT!D0!N"
QpplObj -> Set,	mnemonic="Qppl_" + speciesName	
QpplDataObj = nDataObj -> MakeCopy()
QpplDataObj -> Set, iField=4
QpplDataObj -> Set, gkvObj = QpplObj 
output  = CREATE_STRUCT(output, "Qppl", QpplDataObj)
;
; Make a Q_perp object 
;
QperpObj = QpplObj -> MakeCopy(/NoValues)
QperpObj -> Set, title="Q!D!9x!X" + speciesName + "!N"
QperpObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")c!Ds!NT!D0!N"
QperpObj -> Set, mnemonic="Qperp_" + speciesName
QperpDataObj = nDataObj -> MakeCopy()
QperpDataObj -> Set, iField=5
QperpDataObj -> Set, gkvObj = QperpObj 
output  = CREATE_STRUCT(output, "Qperp", QperpDataObj)
;
; Make a u_|| object 
;
UpplObj =  nObj -> MakeCopy(/NoValues)
UpplObj -> Set,	title="U!D!9#!X" + speciesName + "!N"
UpplObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")c!Ds!N"
UpplObj -> Set,	mnemonic="Uppl_" + speciesName	
UpplDataObj = nDataObj -> MakeCopy()
UpplDataObj -> Set, iField=6
UpplDataObj -> Set, gkvObj = UpplObj 
output  = CREATE_STRUCT(output, "Uppl", UpplDataObj)
;
; ane we're done!
;
FREE_LUN, ioUnit
CD, CurrentWorkingDirectory
RETURN, output

END ; ****** GENE_ReadMom ****** ;
