
FUNCTION GENE_FieldData, 	fileName, ParamStr=ParamStr,	$
				path=path, nFields=nFieldsIn, 	$
				nkx=nkxIn, nky=nkyIn, nz=nzIn, 	$
				nt_FIELD=nt_fieldIN,		$
				 _Extra=Extra
;
; This function reads the FIELDS.dat file output
; by GENE, and returns these fields as GeneData objects.
;
; Arguments:
;
;	fileName	Name of the FIELD.dat file to be analyzed.
;			Defaults to "field.dat". (optional).
;
; Keywords:
;
;	ParamStr	Structure created by GENE_ReadParameters
;			which contains various parameters required
;			to read the field files.	
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
; Written (as GENE_ReadFields) by W.M. Nevins
;	10/19/2008
; Modified (to GENE_FieldData) by W.M. Nevins
;	11/5/2008
; to return GeneField objects 
; (instead of GKV objects)
;
CD, CURRENT=CurrentWorkingDirectory
IF(N_ELEMENTS(path) NE 0) THEN BEGIN
	CD, path
ENDIF ELSE BEGIN
	path=CurrentWorkingDirectory
ENDELSE
nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN fileName="field.dat"
ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find FIELD file", /INFORMATIONAL
	RETURN, 0
ENDIF
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

FileID = "field.dat"
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
; Get a free I/O unit, and open FIELD.dat file
;
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening FIELD file: ", fileName
	PRINT, "Error Message = ", !ERR_STRING
	RETURN, 0
ENDIF
;
; Make a parameter structure if none exists
;
IF( TypeOF(ParamStr) NE 8 ) THEN ParamStr = GENE_ReadParameters()
;
; Compute size of FIELD data using data from 
; parameter structure.
;
IF(TypeOf(ParamStr) EQ 8) THEN BEGIN
	nFields = ParamStr.INFO.N_FIELDS
	nkx     = ParamStr.BOX.NX0
	nky	= ParamStr.BOX.NKY0
	nZ	= ParamStr.BOX.NZ0
	nt	= ParamStr.INFO.ITIME/ParamStr.IN_OUT.ISTEP_FIELD
	Endian	= ParamStr.INFO.ENDIANNESS
	Endian	= STRLOWCASE(Endian)
	Endian	= STRCOMPRESS(Endian, /REMOVE_ALL)
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
start = LONG64(4)	; The FIELD.dat file begins with a 4 byte field
EOR   = LONG64(8)	; Each record is separated by an 8 byte field
;
; These may be operating system dependent, 
; (but, as per note from Florian Merz
;  of 10/22/08 they are not) 
; so must allow the user to change 
; "start" and "EOR" from the command line
;
result = GetKeyWord("start", Extra)
IF(Query_Integer(result)) THEN start = result

result = GetKeyWord("EOR", Extra)
IF(Query_Integer(result)) THEN EOR = result
;
; Size of each field block (in bytes)
;
WordSize = LONG64(8)
dataType = LONG64(5)
IF( STRCMP(Precision, "single") ) THEN BEGIN
	WordSize=LONG64(4)
	dataType=LONG64(4)
;	EOR     = LONG64(4)
ENDIF
nFieldBytes = LONG64(2L*WordSize*nkx*nky*nz)
;
; Compute number of bytes between time values
;
nFields = LONG64(nFields)
iSkip = WordSize + EOR + nFields*(nFieldBytes + EOR)
;
; Read values of t at which fields are written
;
t=MAKE_ARRAY(nt, TYPE=dataType)
FOR i=0L, nt-1 DO BEGIN
;	print, i, start+i*iSkip
	t[i] = READ_BINARY(ioUnit, DATA_START=start+i*iSkip, DATA_TYPE=dataType, DATA_DIMS=1, ENDIAN=Endian)
ENDFOR
;
; Create GKV Object
;
ObjStr = {GKVs4D}
ObjStr.mnemonic = "phi"
ObjStr.title    = "!4u!X"
ObjStr.indices  = PTR_NEW(["*","*","*","*"])
ObjStr.units    = "(!4q!X!Ds!N/" + L_ref + ")(T/e)"
ObjStr.values   = PTR_NEW()
ObjStr.vrange   = [0, 1]
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
	xGrid.range     = [-kxMax, kxMax-dkx]
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

PhiObj = OBJ_NEW("GKVs4D", ObjStr)
phiDataStr = {GeneField}
phiDataStr.GeneDataDir	= path
phiDataStr.ParamStr	= PTR_NEW(ParamStr)
phiDataStr.GkvObj	= PhiObj
phiDataStr.FieldFile	= FileName
phiDataStr.nField	= nFields
phiDataStr.iField	= 1
PhiDataObj = OBJ_NEW("GeneField", phiDataStr)
;
; Prepare output structure
;
output = {	Name	:	"GENE_FieldData",	$
		phi	:	PhiDataObj		}
;
; Make A_|| object if the A_|| field is present
;
IF(nFields GT 1) THEN BEGIN
	ApplObj = PhiObj -> MakeCopy(/NoValues)
	ApplObj -> Set,	title="A!D!9#!X!N"
	ApplObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")!4q!X!Ds!NB!D0!N"
	ApplObj -> Set,	mnemonic="Appl"	
	ApplDataObj = phiDataObj -> MakeCopy()
	ApplDataObj -> Set, iField = 2
	ApplDataObj -> Set, gkvObj = ApplObj
	output  = CREATE_STRUCT(output, "Appl", ApplDataObj)
ENDIF
;
; Make a B_|| object if the B_|| field is present
;
IF(nFields GT 2) THEN BEGIN
	BpplObj = PhiObj -> MakeCopy(/NoValues)
	BpplObj -> Set,	title="B!D!9#!X!N"
	BpplObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")B!D0!N"
	BpplObj -> Set,	mnemonic="Bppl"
	BpplDataObj = phiDataObj -> MakeCopy()
	BpplDataObj -> Set, iField = 3
	BpplDataObj -> Set, gkvObj = BpplObj	
	output  = CREATE_STRUCT(output, "Bppl", BpplDataObj)
ENDIF
;
; ane we're done!
;
FREE_LUN, ioUnit
CD, CurrentWorkingDirectory
RETURN, output

END ; ****** GENE_ReadFields ****** ;
