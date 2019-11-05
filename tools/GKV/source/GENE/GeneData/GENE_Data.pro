;
; *****************************************************************************************************************
; ******************************************     Copyright Notice     *********************************************
; *                                                                                                               *
; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *
; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *
; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *
; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *
; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *
; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * 
; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *
; *                                                                                                               *
; *****************************************************************************************************************
;
; *****************************************************************************************************************
; **********************************************     DISCLAIMER     ***********************************************
; *                                                                                                               *
; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *
; *  Neither the United States Government nor the University of California nor any of their employees, makes      *
; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *
; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *
; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *
; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *
; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *
; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *
; *  States Government or the University of California, and shall not be used for advertising or product          *
; *  endorsement purposes.                                                                                        *
; *                                                                                                               *
; *****************************************************************************************************************
; 
; 
FUNCTION GENE_ReadParameters, fileName, StartGroup=groupstart, EndGroup=groupend
;
; Reads the GENE file "parameters.dat"
; which is in namelist format. Routine
; should work for pretty much any 
; NAMELIST file with one (scalar) variable
; per line.  Data may be (optionally)
; divided into groups (stored as a substructure
; in the output structure returned by this 
; routine) by delimiting elements of the
; group with a line of the form
; 
;&groupName
;
; at the start of a group, and a line
; of the form
;
; /
;
; at the end of the group
;
; ARGUMENT:
;	(ASCII) name of the GENE "parameters.dat" file to be read.
;	Defaluts to "parameters.dat"
;
; Keywords:
;
;	StartGroup	Delimiter indicating the start of a new
;			group of variables.  Defaults to "&"
;
;	EndGroup	Delimiter indicating the end of a
;			group of variables. Defaluts to "/"
;			
;
; Written by W.M. Nevins
;	10/18/2008
;
; Get a free I/O unit, and open PARAMETER.dat file
;
nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN fileName="parameters.dat"
StartGroup = "&"
IF(N_ELEMENTS(groupstart) EQ 1) THEN StartGroup=groupstart
EndGroup = "/"
IF(N_ELEMENTS(groupend) EQ 1)   THEN EndGroup=groupend

T = "True"		; set logical variables T (true) 
F = "False"		; and F (false), as symbols T and F are
;			; used in this context in PARAMETERS.dat files
DOUBLE = "Double"	; Similarly for DOUBLE (double precisons) 
LITTLE = "Little"	; and LITTLE (Little endian).

iSpecies=0
output = 0
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening PARAMETER.dat file:"
	PRINT, "Error Message = ", !ERR_STRING
	RETURN, output
ENDIF
output = { Name : "GENE_ReadParametes" }

thisLine = STRARR(1)
NewGroup : 	blankCount=0
ReadLine :	IF( EOF(ioUnit) ) THEN GOTO, Done
		READF, ioUnit, thisLine

		shortString = STRCOMPRESS(thisLine)
		IF( STRCMP(shortString, "") ) THEN BEGIN ; line is blank
			blankCount=blankCount+1
			GOTO, ReadLine
		ENDIF 
		
		locStart = STRSPLIT(thisLine, startGroup)
		IF(locStart EQ 1) THEN BEGIN		; this line starts a new group
			splitLine = STRSPLIT(thisLine, startGroup, /EXTRACT, Count=nStrings)
			groupName = splitLine[0]
			commandLine = "thisStr = { myName : '" + groupName + "'}"
			OK = EXECUTE(commandLine) 
		ENDIF

		locEnd = STRCMP(STRCOMPRESS(thisLine, /REMOVE_ALL), endGroup)
		IF( locEnd EQ 1) THEN BEGIN ; line ends group
			IF(STRCMP(groupName, "species", /FOLD_CASE)) THEN BEGIN
				groupName = groupName + "_" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
				iSpecies = iSpecies + 1
			ENDIF
			commandLine = "output = CREATE_STRUCT(output, '" + groupName + "', thisStr)"
			OK = EXECUTE(commandLine)
			GOTO, NewGroup
		ENDIF

		splitLIne = STRSPLIT(thisLine, "=", /EXTRACT, Count=nStrings)
		IF(nStrings EQ 2) THEN BEGIN
			; Check for multiple, space-delimited values
			pieces = STRSPLIT(splitLine[1], COUNT=nPieces, /EXTRACT)
			IF(nPieces GT 1) THEN splitLine[1] = '[' + STRJOIN(pieces, ',') + ']'
			OK = EXECUTE("value ="+ splitLine[1])
			IF(OK EQ 0) THEN BEGIN
				commandLine = "value = '"+ splitLine[1] + "'"
				OK = EXECUTE(commandLine)
			ENDIF
			varName = STRCOMPRESS(splitLine[0], /REMOVE_ALL) 
			commandLine = "thisStr = CREATE_STRUCT(thisStr, '" + varName + "', value)"
			OK = EXECUTE(commandLine)
			IF(STRCMP(varName, "n_spec", /FOLD_CASE)) THEN nSpecies=value
		ENDIF
		GOTO, ReadLine
;
; you should only get here if you reached EOF
;
Done:	
FREE_LUN, ioUnit
RETURN, output
END ; ****** GENE_ReadParameters ****** ;



FUNCTION GetMuKnots, totalnw, upperbound
;
; From VSP_LIB.pro as recieved from Florian Merz 10/23/08.
;
; velocity space diags need mu knots as implemented in GENE

  muknots = DBLARR(totalnw)
  
  ; ===================================================
  ; =========== Gauss-Legendre quadrature =============
  ; ===================================================
  ; Definition of knots and weights
  ; according to  Abramowitz, Stegun, Table 25.4, S.916    
  
  CASE totalnw OF
    4  : glknots = [0.339981043584856,0.861136311594053]
    6  : glknots = [0.238619186083197,0.661209386466265,0.932469514203152]
    8  : glknots = [0.183434642495650,0.525532409916329,0.796666477413627,$
                    0.960289856497536]
    10 : glknots = [0.148874338981631,0.433395394129247,0.679409568299024,$
                    0.865063366688985,0.973906528517172]
    16 : glknots = [0.095012509837637440185,0.281603550779258913230,$
           0.458016777657227386342,0.617876244402643748447,0.755404408355003033895,$
           0.865631202387831743880,0.944575023073232576078,0.989400934991649932596]
    24 : glknots = [0.064056892862605626085,0.191118867473616309159,$
            0.315042679696163374387,0.433793507626045138487,$
            0.545421471388839535658,0.648093651936975569252,0.740124191578554364244,$
            0.820001985973902921954,0.886415527004401034213,0.938274552002732758524,$
            0.974728555971309498198,0.995187219997021360180]
    32 : glknots = [0.048307665687738316235,0.144471961582796493485,$
            0.239287362252137074545,0.331868602282127649780,0.421351276130635345364,$
            0.506899908932229390024,0.587715757240762329041,0.663044266930215200975,$
            0.732182118740289680387,0.794483795967942406963,0.849367613732569970134,$
            0.896321155766052123965,0.934906075937739689171,0.964762255587506430774,$
            0.985611511545268335400,0.997263861849481563545]
  ELSE : BEGIN
         PRINT, 'WARNING: Using equidistant mu points'
         glknots=(1.0+2.0*INDGEN(totalnw/2))/totalnw
     END
  ENDCASE
  
  muknots[totalnw/2:totalnw-1] = 0.5*upperBound * glknots + 0.5*upperBound
  FOR i=0,totalnw/2-1 DO $
     muknots[i] = 0.5*upperBound * (-glknots[totalnw/2-1-i]) + 0.5*upperBound
  RETURN, muknots  
  
END ; --- GetMuKnots

;##########################################################################


FUNCTION GENE_ReadVSP, 	fileName, ParamStr=ParamStr,	$
			path=path, _Extra=Extra
;
; This function reads the vsp.dat file output
; by GENE, and returns these fields as GKVs4D objects
; (vs. Vppl, Mu, theta, and t).
;
; Arguments:
;
;	fileName	Name of the vsp.dat file to be analyzed.
;			Defaults to "vsp.dat". (optional).
;
; Keywords:
;
;	Path		Path to structure containing the 
;			desired vsp.dat file. Defaults
;			to current working directory.
;			(Optional)
;
;	ParamStr	Structure created by GENE_ReadParameters
;			which contains various parameters (from
;			the parameters.dat file) normmally required
;			to read the vsp file.  (but see keywords
;			below to work-around absence of a
;			ParamStr.	
;
;	nz		Number of grid points in z.
;			Default is to compute from data in
;			'parameters.dat' file.
;			(Optional)
;
;	nVppl		Number of grid points in v_parallel.
;			Defaults to value in 'parameters.dat'
;			file. (Optional)
;
;	nMu		Number of grid points in Mu.
;			Defaluts to value in 'parameters.dat'
;			file. (Optional)
;
;	nSpecies	Number of kinetic species in the
;			GENE simulation. Defaults to value
;			in the 'parameters.dat' file.
;			(Optional)
;
;	Lv		Extent of the simulation box in v_parallel
;			in units of the thermal velocity, SQRT(2T/M)
;			(note SQRT(2) in this definition).  Simulation
;			domain extends from -Lvppl*vth to +Lvppl*vth.
;			Defaults to value in 'parameters.dat' file.
;			(Optional)
;
;	Lmu		Extent of the simulation box in magnetid 
;			moment in units of the T_j0/B_ref.  Simulation
;			domain extends from 0 to +Lmu*T_0j/B_ref.
;			Defaults to value in 'parameters.dat' file.
;			(Optional)
;
;	nt_vsp		Set to number of time-steps in the
;			vsp.dat file. Default is to compute  
;			from data in 'parameters.dat' file.
;			(Optional)
;
;	Endian		The "Endianness" of the vsp.dat file.
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
; Written by W.M. Nevins
;	10/22/2008
;
nArgs = N_PARAMS()
If(nArgs EQ 0) THEN BEGIN
	fileName="vsp.dat"
	MESSAGE, "Using default filename, 'vsp.dat'", /INFORMATIONAL
ENDIF
CD, CURRENT=CurrentWorkingDirectory
IF(Query_String(path)) THEN CD, path
ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find vsp file, " + fileName, /INFORMATIONAL
	CD, CurrentWorkingDirectory
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

FileID = "VSP File"
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
; Get a free I/O unit, and open VSP.dat file
;
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening FIELD file:
	PRINT, "Error Message = ", !ERR_STRING
	CD, CurrentWorkingDirectory
	RETURN, 0
ENDIF
;
; Compute size of vsp data using data from 
; parameter structure.
;
IF(TypeOf(ParamStr) EQ 8) THEN BEGIN
	nMoments = 5L
	nZ	 = ParamStr.BOX.NZ0
	nVppl	 = ParamStr.BOX.NV0
	nMu	 = ParamStr.BOX.NW0
	nSpecies = ParamStr.BOX.N_SPEC
	Lv	 = ParamStr.BOX.LV
	Lmu	 = ParamStr.BOX.LW
	nt	 = ParamStr.INFO.ITIME/ParamStr.IN_OUT.ISTEP_VSP + 1L
	Endian	 = ParamStr.INFO.ENDIANNESS
	Endian	 = STRLOWCASE(Endian)
	Endian	 = STRCOMPRESS(Endian, /REMOVE_ALL)
ENDIF
;
; Check if this information is available from the command line.
; 

result = GetKeyWord('nz', EXTRA)
IF(Query_Integer(result)) THEN nz=result

result = GetKeyWord('nVppl', EXTRA)
IF(Query_Integer(result)) THEN nVppl=result

result = GetKeyWord('nMu', EXTRA)
IF(Query_Integer(result)) THEN nMu=result

result = GetKeyWord('nSpecies', EXTRA)
IF(Query_Integer(result)) THEN nSpecies=result

result = GetKeyWord('Lv', EXTRA)
IF(Query_Integer(result)+Query_Real(result)) THEN Lv=result

result = GetKeyWord('LMu', EXTRA)
IF(Query_Integer(result)+Query_Real(result)) THEN LMu=result

result = GetKeyWord('nt_vsp', EXTRA)
IF(Query_Integer(result)) THEN nt=result

result = GetKeyWord('Endian', Extra)
IF(Query_String(result)) THEN BEGIN
	IF(result NE 'undefined') THEN BEGIN
		Endian=STRLOWCASE(result)
		Endian=STRCOMPRESS(Endian, /REMOVE_ALL)
	ENDIF
ENDIF
;
; Create array of species names
;
SpeciesNames = STRARR(nSpecies)
FOR i=0, nSpecies-1 DO SpeciesName[i] = STRCOMPRESS(STRING(i), /REMOVE_ALL)
IF(TypeOf(paramStr) EQ 8) THEN BEGIN
	FOR i=0,nSpecies-1 DO BEGIN
		commandLine = "thisName = ParamStr.species_" + speciesName[i] + ".Name"
		OK = EXECUTE(commandLine)
		IF(OK) THEN speciesName[i] = thisName
	ENDFOR
ENDIF
;
; Some key constants garnared from close
; inspection of a sample field.dat file:
;
start = 4L	; The FIELD.dat file begins with a 4 byte field
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
; Compute size of one block of VSP data
;
nMomentBytes = 8L*nVppl*nMu*nZ*nSpecies
;
; Compute number of bytes between time values
; (recalling that data is REAL)
;
iSkip = 8 + EOR + nMoments*nMomentBytes + EOR 
;
; Read values of t at which fields are written
;
t=DBLARR(nt)
FOR i=0,nt-1 DO	$
	t[i] = READ_BINARY(ioUnit, DATA_START=start+i*iSkip, DATA_TYPE=5, DATA_DIMS=1, ENDIAN=Endian)
;
; Read Gamma_ES values
;
thisStart = start + 8L + EOR
Gamma_ES = DCOMPLEXARR(nVppl, nMu, nz, nSpecies, nt)
FOR i=0,nt-1 DO	BEGIN
	Gamma_ES[*,*,*,*,i]=READ_BINARY(ioUnit, DATA_START=thisStart+i*iSkip, DATA_TYPE=5, DATA_DIMS=[nVppl,nMu,nz,nSpecies], ENDIAN=Endian)
ENDFOR
;
; Read Gamma_EM values 
;
thisStart = start + 8L + EOR + nMomentBytes
Gamma_EM = DCOMPLEXARR(nVppl, nMu, nz, nSpecies, nt)
FOR i=0,nt-1 DO	BEGIN
	Gamma_EM[*,*,*,*,i]=READ_BINARY(ioUnit, DATA_START=thisStart+i*iSkip, DATA_TYPE=5, DATA_DIMS=[nVppl,nMu,nz,nSpecies], ENDIAN=Endian)
ENDFOR
;
; Read Q_ES values 
;
thisStart = start + 8L + EOR + 2L*nMomentBytes
Q_ES = DCOMPLEXARR(nVppl, nMu, nz, nSpecies, nt)
FOR i=0,nt-1 DO	BEGIN
	Q_ES[*,*,*,*,i]=READ_BINARY(ioUnit, DATA_START=thisStart+i*iSkip, DATA_TYPE=5, DATA_DIMS=[nVppl,nMu,nz,nSpecies], ENDIAN=Endian)
ENDFOR
;
; Read Q_EM values 
;
thisStart = start + 8L + EOR + 3L*nMomentBytes
Q_EM = DCOMPLEXARR(nVppl, nMu, nz, nSpecies, nt)
FOR i=0,nt-1 DO	BEGIN
	Q_EM[*,*,*,*,i]=READ_BINARY(ioUnit, DATA_START=thisStart+i*iSkip, DATA_TYPE=5, DATA_DIMS=[nVppl,nMu,nz,nSpecies], ENDIAN=Endian)
ENDFOR
;
; Read df_rms values 
;
thisStart = start + 8L + EOR + 4L*nMomentBytes
dF_rms = DCOMPLEXARR(nVppl, nMu, nz, nSpecies, nt)
FOR i=0,nt-1 DO	BEGIN
	dF_rms[*,*,*,*,i]=READ_BINARY(ioUnit, DATA_START=thisStart+i*iSkip, DATA_TYPE=5, DATA_DIMS=[nVppl,nMu,nz,nSpecies], ENDIAN=Endian)
ENDFOR
;
; Create GKV Object to hold Gamma_ES
;
iSpecies=0
ObjStr = {GKVs4D}
mnemonic = "Gamma_ES" 
IF(nSpecies GT 1) THEN mnemonic = Objstr.Mnemonic + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
ObjStr.mnemonic = mnemonic
ObjStr.title    = "!4C!X!DES!N"
IF(nSpecies GT 1) THEN ObjStr.title = ObjStr.title + "!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
ObjStr.indices  = PTR_NEW(["*","*","*","*"])
ObjStr.units    = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!Nc!Ds!N"
values = GAMMA_ES[*,*,*,iSpecies,*]
values = REFORM(values, nVppl,nMu,nz,nt)
ObjStr.values   = PTR_NEW(values)
	vMin    = GKVsd_MIN(values, MAX=vMax)
ObjStr.vrange   = [vMin, vMax]
ObjStr.CodeName = CodeName
ObjStr.CodePI   = CodePI
ObjStr.RunID    = RunID
ObjStr.fileID   = FileID

vGrid = {Grid}
	vGrid.mnemonic  = 'Vppl'
	vGrid.title     = 'v!D!9#!X!N'
	vGrid.units     = 'v!Dth!N'
	dV		= 2.*LV/(nVppl+1)
	Vmax		= Lv
	vValues		= -Vmax + dV*(FINDGEN(nVppl)+1)
	vGrid.values    = PTR_NEW(vValues)
	vGrid.boundary  = "open"
	vGrid.range     = [-Vmax, Vmax]
	vGrid.irange    = [0, nVppl-1]
ObjStr.grid1 = vGrid

muGrid = {Grid}
	muGrid.mnemonic = 'Mu'
	muGrid.title    = '!4l!X'
	muGrid.units    = 'T!Dref!N/B!Dref!N'
;	dMu		= Lmu/(nMu-1)
;	muValues	= dMu*FINDGEN(nMu)
	muValues	= GetMuKNots(nMu, Lmu)
	muGrid.values    = PTR_NEW(muValues)
	muGrid.boundary  = "closed/open"
	muGrid.range     = [0, Lmu]
	muGrid.irange    = [0, nMu-1]
ObjStr.grid2 = muGrid

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

Gamma_ESObj = OBJ_NEW("GKVs4D", ObjStr)
;
; Prepare output structure
;
CommandLine = "output = {Name:'GENE_ReadVSP'," + mnemonic + ":Gamma_ESObj}"
ok = EXECUTE(CommandLine)

IF(nSpecies GT 1) THEN BEGIN
	FOR iSpecies=1,nSpecies-1 DO BEGIN
		mnemonic = "Gamma_ES" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
		title    = "!4C!X!DES" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
		values	 = GAMMA_ES[*,*,*,iSpecies,*]
		values = REFORM(values, nVppl,nMu,nz,nt)
		vMin     = GKVsd_MIN(values, MAX=vMax)
		vrange   = [vMin, vMax]
		thisObj	 = GAMMA_ESObj -> MakeCopy(/NoValues)
		thisObj -> Set, Mnemonic=mnemonic, Title=title, Values=PTR_NEW(values), Vrange=vrange
		output = CREATE_STRUCT(output, mnemonic, thisObj)
	ENDFOR
ENDIF
;
; Make Gamma_EM object(s) 
;
iSpecies = 0
Gamma_EMObj = Gamma_ESObj -> MakeCopy(/NoValues)
mnemonic="Gamma_EM"
IF(nSpecies GT 1) THEN mnemonic=mnemonic + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
title = "!4C!X!DEM!N"
IF(nSpecies GT 1) THEN title = title + "!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
Gamma_EMObj -> Set, title=title
Gamma_EMObj -> Set, mnemonic=mnemonic
values = Gamma_EM[*,*,*,iSpecies,*]	
values = REFORM(values, nVppl,nMu,nz,nt)
Gamma_EMObj -> Set, Values=values
vMin    = GKVsd_MIN(values, MAX=vMax)
Gamma_EMObj -> Set, vRange=[vMin, vMax]
output  = CREATE_STRUCT(output, mnemonic, Gamma_EMObj)

IF(nSpecies GT 1) THEN BEGIN
	FOR iSpecies=1,nSpecies-1 DO BEGIN
		mnemonic = "Gamma_EM" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
		title    = "!4C!X!DEM" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
		values	 = GAMMA_EM[*,*,*,iSpecies,*]
		values = REFORM(values, nVppl,nMu,nz,nt)
		vMin     = GKVsd_MIN(values, MAX=vMax)
		vrange   = [vMin, vMax]
		thisObj	 = GAMMA_EMObj -> MakeCopy(/NoValues)
		thisObj -> Set, Mnemonic=mnemonic, Title=title, Values=PTR_NEW(values), Vrange=vrange
		output = CREATE_STRUCT(output, mnemonic, thisObj)
	ENDFOR
ENDIF
;
; Make Q_ES object(s)
;
iSpecies=0
Q_ESObj = Gamma_ESObj -> MakeCopy(/NoValues)
mnemonic="Q_ES"
IF(nSpecies GT 1) THEN mnemonic=mnemonic + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
title="Q!DES!N"
IF(nSpecies GT 1) THEN title = title + "!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
units="(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!NT!D0!Nc!Ds!!N"
Q_ESObj -> Set, Title=title
Q_ESObj -> Set, Mnemonic=Mnemonic
Q_ESObj -> Set, Units=units
values = Q_ES[*,*,*,iSpecies,*]	
values = REFORM(values, nVppl,nMu,nz,nt)
Q_ESObj -> Set, values = PTR_NEW(values)
vMin    = GKVsd_MIN(values, MAX=vMax)
Q_ESObj -> Set, vRange=[vMin, vMax]
output  = CREATE_STRUCT(output, mnemonic, Q_ESObj)

IF(nSpecies GT 1) THEN BEGIN
	FOR iSpecies=1,nSpecies-1 DO BEGIN
		mnemonic = "Q_ES" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
		title    = "Q!DDES" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
		values	 = Q_ES[*,*,*,iSpecies,*]
		values = REFORM(values, nVppl,nMu,nz,nt)
		vMin     = GKVsd_MIN(values, MAX=vMax)
		vrange   = [vMin, vMax]
		thisObj	 = GAMMA_EMObj -> MakeCopy(/NoValues)
		thisObj -> Set, Mnemonic=mnemonic, Title=title, Values=PTR_NEW(values), Vrange=vrange
		output = CREATE_STRUCT(output, mnemonic, thisObj)
	ENDFOR
ENDIF
;
; Make a Q_EM object(s) 
;
iSpecies=0
Q_EMObj =  Q_ESObj -> MakeCopy(/NoValues)
mnemonic="Q_EM"	
IF(nSpecies GT 1) THEN mnemonic=mnemonic + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
title="Q!DEM!N"
IF(nSpecies GT 1) THEN title = title + "!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
Q_EMObj -> Set,	Title=title
Q_EMObj -> Set,	Mnemonic=mnemonic
values = Q_EM[*,*,*,iSpecies,*]
values = REFORM(values, nVppl,nMu,nz,nt)
Q_EMObj -> Set,	values = PTR_NEW(values)
vMin    = GKVsd_MIN(values, MAX=vMax)
Q_EMObj -> Set, vRange=[vMin, vMax]
output  = CREATE_STRUCT(output, mnemonic, Q_EMObj)

IF(nSpecies GT 1) THEN BEGIN
	FOR iSpecies=1,nSpecies-1 DO BEGIN
		mnemonic = "Q_EM" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
		title    = "Q!DDEM" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N"
		values	 = Q_EM[*,*,*,iSpecies,*]
		values = REFORM(values, nVppl,nMu,nz,nt)
		vMin     = GKVsd_MIN(values, MAX=vMax)
		vrange   = [vMin, vMax]
		thisObj	 = GAMMA_EMObj -> MakeCopy(/NoValues)
		thisObj -> Set, Mnemonic=mnemonic, Title=title, Values=PTR_NEW(values), Vrange=vrange
		output = CREATE_STRUCT(output, mnemonic, thisObj)
	ENDFOR
ENDIF
;
; Make a dF_rms object(s)
;
iSpecies=0
mnemonic="dF_rms"
IF(nSpecies GT 1) THEN mnemonic=mnemonic + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
title="!12<!4d!Xf"
IF(nSpecies GT 1) THEN title = title + "!S!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N!R"
title=title + "!U2!N!12>!X!U1/2!N"
dF_rmsObj = Q_EMObj -> MakeCopy(/NoValues)
dF_rmsObj -> Set, Title=title
dF_rmsObj -> Set, Units="(!4q!X!Ds!N/" + L_ref + ")(n!D0!N/v!Dth!U3!N"
dF_rmsObj -> Set, Mnemonic=mnemonic
values = dF_rms[*,*,*,iSpecies,*]	
values = REFORM(values, nVppl,nMu,nz,nt)
dF_rmsObj -> Set, Values=values
vMin    = GKVsd_MIN(values, MAX=vMax)
dF_rmsObj -> Set, vRange=[vMin, vMax]
output  = CREATE_STRUCT(output, mnemonic, dF_rmsObj)

IF(nSpecies GT 1) THEN BEGIN
	FOR iSpecies=1,nSpecies-1 DO BEGIN
		mnemonic = "dF_rms" + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
		title ="!12<!4d!Xf"
		title = title + "!S!D"+ STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL) + "!N!R"
		title =title + "!U2!N!12>!X!U1/2!N"
		values	 = Q_EM[*,*,*,iSpecies,*]
		values = REFORM(values, nVppl,nMu,nz,nt)
		vMin     = GKVsd_MIN(values, MAX=vMax)
		vrange   = [vMin, vMax]
		thisObj	 = GAMMA_EMObj -> MakeCopy(/NoValues)
		thisObj -> Set, Mnemonic=mnemonic, Title=title, Values=PTR_NEW(values), Vrange=vrange
		output = CREATE_STRUCT(output, mnemonic, thisObj)
	ENDFOR
ENDIF
;
; ane we're done!
;
FREE_LUN, ioUnit
CD, CurrentWorkingDirectory
RETURN, output

END ; ****** GENE_ReadVSP ****** ;


FUNCTION GENE_NRG_Data, fileName, _Extra=extra
;
; This function routine reads the NRG.dat file from 
; GENE, and returns the intensity of the density,
; parallel velocity, perpendicular and parallel 
; temperture fluctuations, and the particle and
; heat fluxes of each species.
;
; Arguments:
;
;	fileName	Name of the NRG.dat file to be analyzed.
;			Defaults to "nrg.dat". (optional).
;
; Keywords:
;
;	paramStr	The structure returned by GENE_ReadParameters
;			which contains information from the namelist
;			input file, parameters.dat
;
;
;	nt_NRG		Set to number of time-steps in the
;			nrg file to over-ride computation 
;			of nt from data in 'parameter.dat' file, 
;			or if no 'parameter.dat' file exists.
;			(Optional)
;
;	nSpecies	Set to number of species in GENE
;			simulation if no "parameter.dat" file is
;			available. Defaults to 1 in absence
;			of 'parameter.dat' file. 
;			(Optional)
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
; Written by W.M. Nevins NRG_Data (a stand-alone function)
;	6/16/05
;
; Modified by W.M. Nevins
;	10/18/2008
; To work as a module within GENE_Data
;
;
ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find NRG file", /INFORMATIONAL
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

FileID = "NRG"
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
; Get a free I/O unit, and open nrg.dat file
;
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening NRG file:
	PRINT, "Error Message = ", !ERR_STRING
	RETURN, 0
ENDIF

result = GetKeyWord('nt_NRG', Extra)
IF(Query_Integer(result)) THEN nt=result
IF(N_ELEMENTS(nt) EQ 0) THEN BEGIN
	PRINT, "Number of timesteps in NRG FILE"
	PRINT, "should be specified."
	PRINT, "will use nt_NRG=10000"
	nt=10000
ENDIF
IF(nt LT 1) THEN nt=10000

nSpecies = 1
result = GetKeyWord('nSpecies', Extra)
IF(Query_Integer(result)) THEN nSpecies=result		
;
; Create array of species names
;
result = GetKeyWord("Paramstr", Extra)
IF(TypeOf(result) EQ 8) THEN paramstr=result
SpeciesName = STRARR(nSpecies)
FOR i=0, nSpecies-1 DO SpeciesName[i] = STRCOMPRESS(STRING(i), /REMOVE_ALL)
IF(TypeOf(paramStr) EQ 8) THEN BEGIN
	FOR i=0,nSpecies-1 DO BEGIN
		commandLine = "thisName = ParamStr.species_" + speciesName[i] + ".Name"
		commandLine = STRCOMPRESS(commandLine, /REMOVE_ALL)
		OK = EXECUTE(commandLine)
		IF(OK) THEN speciesName[i] = STRCOMPRESS(thisName, /REMOVE_ALL)
	ENDFOR
ENDIF

time    = FLTARR(nt)

nSq     = FLTARR(nspecies, nt)
uparSq  = FLTARR(nspecies, nt)
TparSq  = FLTARR(nspecies, nt) 
TperSq  = FLTARR(nspecies, nt)
GammaES = FLTARR(nspecies, nt)
GammaEM = FLTARR(nspecies, nt)
Q_es    = FLTARR(nspecies, nt)
Q_em    = FLTARR(nspecies, nt)

nSqArr     = OBJARR(nSpecies)
uparSqArr  = OBJARR(nSpecies)
TparSqArr  = OBJARR(nSpecies)
TperSqArr  = OBJARR(nSpecies)
GammaESArr = OBJARR(nSpecies)
GammaEMArr = OBJARR(nSpecies)
Q_esArr    = OBJARR(nSpecies)
Q_emArr    = OBJARR(nSpecies)

line = FLTARR(8)

FOR j=0L,nt-1 DO BEGIN
	IF(EOF(ioUnit)) THEN GOTO, DONE
	READF, ioUnit, t
	time[j]=t
	FOR i=0, nSpecies-1 DO BEGIN
		READF, ioUnit, line
		nSQ[i,j]     = line[0]
		uparSq[i,j]  = line[1]
		TparSq[i,j]  = line[2]
		TperSq[i,j]  = line[3]
		GammaES[i,j] = line[4]
		GammaEM[i,j] = line[5]
		Q_es[i,j]    = line[6]
		Q_em[i,j]    = line[7]
	ENDFOR
ENDFOR
DONE : nnt=j
FREE_LUN, ioUnit
;
; Create GKVs1D objects to contain n^2 for first species
;
nSqObj = {GKVs1D}
nSqObj.mnemonic = "nSq_" + speciesName[0]
nSqObj.title    = "!12<!Xn!U2!N!12>!X!D" + speciesName[0] + "!N"
nSqObj.indices  = PTR_NEW(["*"])
nSqObj.units    = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!U2!N"
       values   = REFORM(nSq[0,0:nnt-1])
nSqObj.values   = PTR_NEW(values)
nSqObj.vrange   = [0,MAX(values)]
nSqObj.CodeName = CodeName
nSqObj.CodePI   = CodePI
nSqObj.RunID    = RunID
nSqObj.fileID   = FileID

tGrid =  {Grid}
	tGrid.mnemonic  = 't'
	tGrid.title     = 't'
	tGrid.units     = L_ref + "/c!Ds!N"
	nTime           = time[0:nnt-1]
	tGrid.values    = PTR_NEW(nTime)
	tGrid.boundary  = "open"
	tGrid.range     = [time[0], time[nnt-1]]
	tgrid.irange    = [0,nnt-1]
nSqObj.grid1 = tGrid
nSqObj = OBJ_NEW("GKVs1D", nSqObj)

nSqArr[0] = nSqObj
IF(nSpecies GT 1) THEN BEGIN
	FOR i=1,nSpecies-1 DO BEGIN
		thisSpecies = speciesName[i]
		nSqArr[i] = nSqObj -> MakeCopy(/NoValues)
		values    = REFORM(nSq[i,0:nnt-1])
		nSqArr[i] -> Set, values=PTR_NEW(values),	$
				  vrange=[0,MAX(values)],	$
				  title ="!12<!Xn!U2!N!12>!X!D"+thisSpecies+"!!N",	$
				  mnemonic = "nSq_" + thisSpecies
	ENDFOR
ENDIF

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	uparSqArr[i] = nSqObj -> MakeCopy(/NoValues)
	values       = REFORM(uparSq[i,0:nnt-1])
	uparSqArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [0, MAX(values)],		$
				mnemonic="uparSq_" + thisSpecies, 	$
				title= "!12<!Xu!S!L!9#!X!R!U2!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nc!Ds!U2!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	TparSqArr[i] = nSqObj -> MakeCopy(/NoValues)
	values       = REFORM(TparSq[i,0:nnt-1])
	TparSqArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [0, MAX(values)],		$
				mnemonic="TparSq_" + thisSpecies, 	$
				title="!12<!XT!S!L!9#!X!R!U2!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!NT!D0!N!U2!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	TperSqArr[i] = nSqObj -> MakeCopy(/NoValues)
	values       = REFORM(TperSq[i,0:nnt-1])
	TperSqArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [0, MAX(values)],		$
				mnemonic="TperSq_" + thisSpecies, 	$
				title="!12<!XT!S!L!9x!X!R!U2!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!NT!D0!N!U2!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	GammaESArr[i] = nSqObj -> MakeCopy(/NoValues)
	values        = REFORM(GammaES[i,0:nnt-1])
	vmax          = MAX(values, MIN=vmin)
	GAmmaESArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [vmin, vmax],			$
				mnemonic="Gamma_ES_" + thisSpecies, 	$
				title="!12<!4C!X!DES!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!Nc!Ds!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	GammaEMArr[i] = nSqObj -> MakeCopy(/NoValues)
	values        = REFORM(GammaEM[i,0:nnt-1])
	vmax          = MAX(values, MIN=vmin)
	GAmmaEMArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [vmin, vmax],			$
				mnemonic="Gamma_EM_" + thisSpecies, 	$
				title="!12<!4C!X!DEM!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!Nc!Ds!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	Q_esArr[i] = nSqObj -> MakeCopy(/NoValues)
	values        = REFORM(Q_es[i,0:nnt-1])
	vmax          = MAX(values, MIN=vmin)
	Q_esArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [vmin, vmax],			$
				mnemonic="Q_ES_" + thisSpecies, 	$
				title="!12<!XQ!DES!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!NTc!Ds!N" 
ENDFOR

FOR i=0,nSpecies-1 DO BEGIN
	thisSpecies = speciesName[i]
	Q_emArr[i] = nSqObj -> MakeCopy(/NoValues)
	values        = REFORM(Q_em[i,0:nnt-1])
	vmax          = MAX(values, MIN=vmin)
	Q_emArr[i] -> Set, 	values= PTR_NEW(values),		$
				vrange= [vmin, vmax],			$
				mnemonic="Q_EM_" + thisSpecies, 	$
				title="!12<!XQ!DEM!N!12>!X!D"+thisSpecies+"!N", 	$
				units = "(!4q!X!Ds!N/" + L_ref + ")!U2!Nn!D0!NTc!Ds!N" 
ENDFOR

result = {	Name		:	"GENE_NRG_Data",			$
		nSqArr		:	nSqArr,				$
		uparSqArr	:	uparSqArr,			$
		TparSqArr	:	TparSqArr,			$
		TperSqArr	:	TperSqArr,			$
		Gamma_ESArr	:	GammaESArr,			$
		Gamma_EMArr	:	GammaEMArr,			$
		Q_ESArr		:	Q_esArr,			$
		Q_EMArr		:	Q_emArr				}

RETURN, result

END ;  ****** GENE_NRG_Data.pro ******  ;


FUNCTION GENE_Data, 	FileName=fileName, Path=path, Search=search,	$
			nt = nt, nSpecies=nSpecies,			$
			Debug=d, Stick=stick, _Extra=Extra
;
; 
;
; Reads "*.dat" files produced by GENE and returns data as GKV objects.
;
; KEYWORDS:
;
;	FileName		Name of GENE .dat file(s).  Must include either full path, or else the path from the 
;				current working directory.  If FileName is not provided, all "*.dat" files in the 
;				selected directory (see "path" keyword below) will be read.  DIALOG_PICKFILE can 
;				be called to allow the user to select the directory containing the appropriate 
;				GENE .dat file(s).
;				(Optional).
;
;	Path			Sets path to directory containing the GENE .dat files. If the SEARCH keyword is set
;				(see below) then DIALOG_PICKFILE will start it's search in this directory.
;				(Optional)
;
;	Search			Set this keyword (i.e., put "/Search" on the command line) to use DIALOG_PICKFILE
;				to search for the directory containing the GENE .dat files. The "Path" keyword (see above)
;				can be used to specify the initial directory to be examined by DIALOG_PICKFILE 
;				(which will allow user to change to other directories...).  If no path is provided, 
;				DIALOG_PICKFILE will default to current working directory. 
;				(Optional).
;
;	nt_NRG			Number of time steps in the NRG file. This is required in the absence of a 'parameters.dat' file.
;				If enetered on the GENE_Data command line, this value will override contents of
;				the 'parameters.dat' file. (Optional)
;
;	nSpecies		Number of species. This is required in the absence of a 'parameters.dat' file.
;				If enetered on the GENE_Data command line, this value will override contents of
;				the 'parameters.dat' file. (Optional)
;
;	ArrayStyle		Set to 'c' if netCDF file was written by code (like C or PASCAL) which stores array
;				data in 'column major' format.  Default is 'fortran', appropriate when the netCDF file
;				was written by code (like FORTRAN or IDL) which stores array data in 'row major' format.
;				This issue is discussed in Chapter 5 of 'Building IDL Applications'.
;				Defaults to 'fortran' (Optional).
;
;	L_ref			Set this to a Hershey vector font string identifying the the (macroscopic) 
;				reference length.  Defaults to "R!D0!N".  (Optional)
;
;	CodeName		Set to a string constant containing the desired 'CodeName' (ascii field appearing at the
;				penultimate line at the lower left hand corner of the plot).
;				Defaults to 'GENE'. (Optional)
;
;	CodePI			Set to a string constant containing the desired 'CodePI' (ascii field appearing at the
;				last line at the lower left hand corner of the plot).  
;				Defaults to "F. Jenko". (Optional)
;
;
;	FileID			Set to a string constant containing the desired 'FileID' (ascii field appearing at the 
;				last line in the lower right hand corner of the plot).  This value will over ride 
;				the GENE .dat filename. (Optional)
;
;	RunID			Set to a string constant containing the desired 'RunID' (ascii field appearing in the 
;				penultimate line at the lower right hand corner of the plot).  (Optional)
;
;	Stick			Set this keyword (i.e., put '/Stick' on the command line) to cause GENE_Data to
;				exit with the current working directory set to the directory containing the GENE.dat
;				files.  Default is to leave current working directory unchanged. (Optional)
;
;
;	Debug			Set (i.e., put '/Debug' within the argument string) to turn off error trapping.  Turning
;				off error trapping is useful if you are debugging GENE_Data.	
;
;
; Written by W.M. Nevins
;	10/18/2008
deBug=0
IF(N_ELEMENTS(d) NE 0) THEN deBug=d
IF(deBug EQ 0) THEN BEGIN					; Set Keyword 'DeBug' to avoid error trap
	Catch, error 						; Set up error trap
	IF error NE 0 then Begin 
		Catch, /Cancel             			; Cancel error trap 
		ok=Error_Message(/Trace)  			; Print out a trace-back in the output log.
;******		IF(cdfID GE 0) THEN *****, cdfID		; Close any open files
		IF(N_ELEMENTS(output) NE 0) THEN $
   			RETURN, output				; Return SOMETHING useful if possible;
   		RETURN, 0                  			;	return 0 if all else fails. 
	ENDIF
ENDIF
CD, CURRENT=CurrentWorkingDirectory
IF(N_ELEMENTS(path) EQ 1) THEN CD, path
;
; Use DIALOG_PICKFILE to find directory containing GENE .dat files if needed
;
IF( KEYWORD_SET(Search) ) THEN BEGIN
	GeneDir=DIALOG_PICKFILE(/DIRECTORY) 
	CD, GeneDir
ENDIF
nFiles=N_ELEMENTS(fileName)
IF(nFiles EQ 0) THEN BEGIN
	datFiles=FILE_SEARCH("*.dat")
	nDatFiles = N_ELEMENTS(datFiles)
ENDIF ELSE BEGIN
		datFiles=STRARR(nFiles)
		checkFiles=fileName
		nDatFiles=0
	FOR i=0,nFiles-1 DO BEGIN
		result = STRSPLIT(fileName[i], ".", /EXTRACT, COUNT=nSubstrings)
		IF(STRCMP(result[nSubStrings-1], "dat", /FOLDCASE)) THEN BEGIN
			checkFiles[i]=fileName[i]
		ENDIF ELSE BEGIN
			checkFiles[i]=fileName[i] + ".dat"
		ENDELSE
		datFiles[nDatfiles] = FILE_SEARCH(checkFiles[i])
		IF(STRCMP(datFiles[nDatFiles], "") EQ 1) THEN nDatFiles=nDatFiles+1
	ENDFOR
	IF(nDatFiles EQ 0) THEN BEGIN
		MESSAGE, "Couldn't find any *.dat files, RETURNING", /INFORMATIONAL
		RETURN, 0
	ENDIF
	datFiles=datFiles[0:nDatFiles-1]
ENDELSE


IF (nDatfiles EQ 0) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find *.dat file", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Create output structure
;
output = { NAME : "GENE_Data" }
;
; Find "parameters.dat" file, and read it
;
noParams = 1b
FOR i=0,nDatFiles-1 DO BEGIN
	thisFile = STRSPLIT(datFiles[i], ".", /EXTRACT, COUNT=nSubstrings)
	IF( STRCMP(thisFile[0], "parameters", /FOLD_CASE) ) THEN BEGIN
		params = GENE_ReadParameters("parameters.dat")
		output = CREATE_STRUCT(output, "Parameters", params)
		IF(N_ELEMENTS(nt_NGR) EQ 0)   THEN nt_NRG=params.info.iTime/params.in_out.istep_nrg
		IF(N_ELEMENTS(nSpecies) EQ 0) THEN nSpecies = params.box.n_spec
	ENDIF
	noParams = 0b
ENDFOR
;
; Find "nrg.dat" file, and read it
;
IF(TypeOF(extra) EQ 8) THEN nrgExtra=Extra
FOR i=0,nDatFiles-1 DO BEGIN
	thisFile = STRSPLIT(datFiles[i], ".", /EXTRACT, COUNT=nSubstrings)
	IF( STRCMP(thisFile[0], "nrg", /FOLD_CASE) ) THEN BEGIN
		nrgStr = GENE_NRG_DATA(datFiles[i], nt_NRG=nt_NRG, nSpecies=nSpecies, ParamStr=params, _Extra=nrgExtra)
		output = CREATE_STRUCT(output, "NRG", nrgStr)
	ENDIF
ENDFOR
;
; Can't proceed without info from "parameters.dat" file.
; In extremis, you and run GENE_READFIELDS directly,
; entering values for nFields, nkx, nky, nz, and nt_field.
; (see comments at top of GENE_READFIELDS).
;
IF(noParams) THEN RETURN, output
;
; Find 'field.dat" file, and read it
;
IF(TypeOF(extra) EQ 8) THEN fieldExtra=Extra
FOR i=0,nDatFiles-1 DO BEGIN
	thisFile = STRSPLIT(datFiles[i], ".", /EXTRACT, COUNT=nSubstrings)
	IF( STRCMP(thisFile[0], "field", /FOLD_CASE) ) THEN BEGIN
		fieldStr = GENE_FieldData(datFiles[i], paramStr=params, _Extra=fieldExtras)
		output = CREATE_STRUCT(output, "FIELD", fieldStr)
	ENDIF
ENDFOR
;
; Find mom_xxx.dat files, and read them
;
speciesName = STRARR(nSpecies)
FOR iSpecies=0,nSpecies-1 DO BEGIN
	commandLine = "speciesName[iSpecies] = params.species_" 
	commandLIne = commandLine + STRCOMPRESS(STRING(iSpecies), /REMOVE_ALL)
	commandLine = commandLine + ".name"
	ok = EXECUTE(commandLine)
	IF(NOT ok) THEN PRINT, commandLine
ENDFOR

FOR iSpecies=0,nSpecies-1 DO BEGIN
	IF(TypeOF(Extra) EQ 8) THEN momExtra=EXTRA
	thisSpecies = STRCOMPRESS(speciesName[iSpecies], /REMOVE_ALL) 
	thisTag = "mom_" + thisSpecies
	momFile = thistag + ".dat"
	momFile = STRCOMPRESS(momFile, /REMOVE_ALL)
	thisFile = FILE_SEARCH(momFile, count=nMomFile)
	IF(nMomFile EQ 1) THEN BEGIN
		momStr = GENE_MomData(momFile, ParamStr=params, _EXTRA=momExtra)
		output = CREATE_STRUCT(output, thisTag, momStr)
	ENDIF
ENDFOR

CD, CurrentWorkingDirectory
RETURN, output
END  ; ****** FUNCTION GENE_Data ****** ;		
