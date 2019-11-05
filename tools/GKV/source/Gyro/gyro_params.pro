FUNCTION Gyro_Params, filename

CD, CURRENT=current_working_directory

separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'

IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out', Get_Path=NetCDF_DIR) ELSE fileIn=fileName
ok = FINDFILE(fileIn, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ; 	we didn't find a file... 
    MESSAGE, "Couldn't find run.out file", /INFORMATIONAL
    RETURN, 0
ENDIF

found:
output = 0

GET_LUN, ioUnit
OPENR, ioUnit, fileIn, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening run.out file:"
    PRINT, "Error Message = ", !ERR_STRING
    RETURN, output
ENDIF

output = { Name : "Gyro_ReadParameters" }

thisLine = STRARR(1)
READF, ioUnit, thisLine  ;=======
READF, ioUnit, thisLine  ; GYRO FRANKLIN

READF, ioUnit, thisLine  ;(C) 
READF, ioUnit, thisLine  ; =====
READF, ioUnit, thisLine  ; Platform
READF, ioUnit, thisLine  ; Machine ID
READF, ioUnit, thisLine  ; [GYRO MPI tasks
READF, ioUnit, thisLine  ; [Parsing data
blankCount = 0


groupName = "Methods"
commandLine = "thisStr = { myName: '" + groupName +"'}"
OK = EXECUTE(commandLine)

NewGroup : 	groupCount=3
ReadLine :	IF( EOF(ioUnit) ) THEN GOTO, Done
READF, ioUnit, thisLine


shortString = STRCOMPRESS(thisLine)
IF( STRCMP(shortString, "") ) THEN BEGIN ; line is blank
    blankCount=blankCount+1
    GOTO, ReadLine
ENDIF      
          
INFO = STRSPLIT(shortString, "INFO") ; info line
IF INFO[0] EQ 4 THEN BEGIN
    GOTO, ReadLine
ENDIF

Note = STRSPLIT(shortString, "Note", count=nStrings) ;note line
IF nStrings GT 8  THEN BEGIN
    ;PRINT, "got it"
    GOTO, ReadLine
ENDIF
;print, shortString, " and ", nStrings

out = STRSPLIT(thisLine, " out")
IF out[0] EQ 4 THEN BEGIN
    GOTO,ReadLine
ENDIF

explicit = STRSPLIT(thisLine, " exp") ;skip line FOR NOW
IF explicit[0] EQ 4 THEN BEGIN
 ;   PRINT, "got it"
    GOTO, ReadLine
ENDIF

explicit2 = STRSPLIT(thisLine, " EXP") ;skip line FOR NOW
IF explicit2[0] EQ 4 THEN BEGIN
 ;   PRINT, "GOT IT"
    GOTO, ReadLine
ENDIF

explicit3 = STRSPLIT(thisLine, " off") ;skip line FOR NOW
IF explicit3[0] EQ 4 THEN BEGIN
  ;  PRINT, "Got It"
    GOTO, ReadLine
ENDIF

nonzero1 = STRSPLIT(thisLine, " nonzeros") ;skip line FOR NOW
IF nonzero1[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero2 = STRSPLIT(thisLine, " values") ;skip line FOR NOW
IF nonzero2[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero3 = STRSPLIT(thisLine, " indices") ;skip line FOR NOW
IF nonzero3[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero4 = STRSPLIT(thisLine, " iterations") ;skip line FOR NOW
IF nonzero4[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

quasi = STRSPLIT(thisLine, " n_i") ; skip quasineutrality line
IF quasi[0] EQ 4 THEN GOTO, ReadLine

please = STRSPLIT(thisLIne, " PLE") ; end of file, with Please see comment
IF please[0] EQ 4 THEN GOTO, DONE


;
; check for line full of -------------------- only
;
ignore = STRSPLIT(thisLine, "-", /EXTRACT, Count=nStrings)
ignore2 = STRSPLIT(ignore[0], " pro")
IF (ignore2[0] EQ 4) THEN BEGIN 
    GOTO, Skip 
ENDIF
IF( STRCMP(ignore[0], "") ) THEN BEGIN ; line is only -------
    GOTO, ReadLine
ENDIF 

Skip: ;skip ignore ------ command for s-alpha boundary

;
; Start group processing
;
separate = STRSPLIT(shortString, "--", Count=gr, /EXTRACT)
IF(gr GT 1) THEN BEGIN    ;First type of group  ----- name -----
    ;PRINT, "FOUND GROUP ", separate
    splitLine = STRSPLIT(thisLine, ":", /EXTRACT, Count=nStrings)
    IF(nStrings EQ 2) THEN GOTO, notgroup
    IF (groupCount GT 1) THEN BEGIN ;not the first group, so close the group
        commandLine = "output = CREATE_STRUCT(output, '" + groupName + "',thisStr)"
        OK = EXECUTE(commandLine)
    ENDIF
    gSegs = n_elements(separate)
    groupName = separate[1]
    IF gSegs EQ 3 THEN groupName = groupName + " " + separate[2]
    gLength = STRLEN(groupName)
    groupName = STRMID(groupName,1,gLength-2)
    ;
    ; handle parenthesis in group name
    ; 
    groupParens = STRSPLIT(groupName, "(", /EXTRACT)
    parenLength = n_elements(groupParens)
    IF parenLength NE 1 THEN BEGIN 
        groupName = groupParens[0] +groupParens[1]
        groupParens = STRSPLIT(groupName, ")", /EXTRACT)        
        groupName = groupParens[0] 
    ENDIF
    ;
    ; handle commas in group name
    ;
    comma = STRSPLIT(groupName, ",", /EXTRACT)
    commaLength = N_elements(comma)
    IF commaLength EQ 3 THEN     groupName = comma[0] + " " + comma[1] + " " + comma[2]  
    commandLine = "thisStr = { myName: '" + groupName +"'}"
    OK = EXECUTE(commandLine)
    groupCount = groupCount + 1
    GOTO, ReadLine
ENDIF
;
; Input a value with name.  
;
splitLine = STRSPLIT(thisLine, ":", /EXTRACT, Count=nStrings)
notgroup: ; Should only get here if there is a value to be stored.
IF(nStrings EQ 2) THEN BEGIN 
                                ; Check for multiple, space-delimited values
    pieces = STRSPLIT(splitLine[1], COUNT=nPieces, /EXTRACT)
    IF(nPieces GT 1) THEN splitLine[1] = STRJOIN(pieces, ' ')
    ;OK = EXECUTE("value ="+ splitLine[1])
    ;IF(OK EQ 0) THEN BEGIN
    ;    commandLine = "value = '"+ splitLine[1] + "'"
    ;    OK = EXECUTE(commandLine)
    ;ENDIF
    value = splitLine[1]
    period = STRSPLIT(thisLine," .", Count=nStrings)

    varName = STRCOMPRESS(splitLine[0], /REMOVE_ALL) 
    ;
    ; check for period at beginning of name
    ;
    IF (period[0] EQ 2) THEN BEGIN
        varname = Strmid(varName, 1)        
    ENDIF
    ;
    ; check for period in middle of name
    ;
    periodMid = STRSPLIT(varName, ".",Count = nStrings, /EXTRACT)  
    IF (nStrings GT 1) THEN varName = periodMid[0] + "_" + periodMid[1]
    ;
    ; check for period at end
    ;    
    periodEnd = STRPOS(varName, ".") ; check for period at end
    varL = STRLEN(periodEnd)
    IF (periodEnd GT 0) THEN varName = STRMID(varName, 0, varL-1)
    ;
    ; check for =>
    ;
    arrow = STRSPLIT(varName, "=>", Count=nStrings)
    IF (arrow[0] GT 1) THEN varname = STRMID(varname,2)    
    ;
    ; replace "/" with "over"
    ;
    over = STRSPLIT(varName, "/", /EXTRACT,Count=nStrings)
    IF (nStrings GT 1) THEN varName = over[0] + " over " + over[1]
    ;
    ; replace "asdf(fdsa)" with "asdf_of_fdsa"
    ;
    parens = STRSPLIT(varName, "(", /EXTRACT)
    parenLength = n_elements(parens)
    IF parenLength NE 1 THEN BEGIN 
        varName = parens[0] + " of " + parens[1]
        varL = STRLEN(varName)
        varName = STRMID(varName, 0, varL-1)
    ENDIF
    ;
    ; replace "(asdf)"with "asdf"
    ;
    fullParens = STRSPLIT(varName, "(")
    varLen= STRLEN(varName)
    IF fullParens[0] EQ 1 THEN varName = STRMID(varName,1,varLen-2)
    ;
    ; replace "a^n" with "a_exp_n"
    ;
    exponent = STRSPLIT(varName, "^", count= nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = exponent[0] + "exp" + exponent[1]
    ;
    ; replace "asdf-fdsa" with "asdf_fdsa"
    ;
    dash = STRSPLIT(varName, "-", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = dash[0] + "_" + dash[1]
    ;
    ; replace "asdf*fdsa" with "asdf_times_fdsa"
    ; 
    mult = STRSPLIT(varName, "*", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = mult[0] + "_times_" + mult[1]
    ;
    ; handle "asdf #1" situations to "asdf number 1"
    ;
    sharp = STRSPLIT(varName, "#", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varName = sharp[0] + "_num_" + sharp[1]
    ;
    ; For "# Ion1" ignore
    ;
    sharp2 = STRSPLIT(varName, '#')
    IF (sharp2[0] EQ 1) THEN GOTO, ReadLine    
    ;
    ; done editing varName.  ready to save value!
    ;
    ;PRINT, typeOf(value), " ", value, "command line ", commandline
    IF (TypeOf(value) EQ 7) THEN BEGIN 
        valueShort =  STRCOMPRESS(value)
        valueCount = N_ELEMENTS(valueShort)
        value = valueShort[0]
        FOR i = 1,valueCount-1 DO BEGIN 
            value = value + " " + valueShort[ir]        
            ;PRINT, "value = ", value
        ENDFOR
        commandLine = "thisStr = CREATE_STRUCT(thisStr, '" + varName + "', '" + value +"')"  
    ENDIF ELSE commandLine = "thisStr = CREATE_STRUCT(thisStr, '" + varName + "', value)"
    
    OK = EXECUTE(commandLine)
    IF(STRCMP(varName, "n_spec", /FOLD_CASE)) THEN nSpecies=value
ENDIF
GOTO, ReadLine

;
; EOF reached
;
DONE:

;  Need to close last group
commandLine = "output = CREATE_STRUCT(output, '" + groupName + "',thisStr)"
OK = EXECUTE(commandLine)

FREE_LUN, ioUnit                

CD, current_working_directory


RETURN, output
END

FUNCTION GKVs1D::Gyro_Params, filename

;
; Routine to read in parameters found in GYRO run.out file
;
; Argument:
;           filename-  string text containing full path to run.out
;                     file, just name of run.out if user is in current directory.
;
; 
; Written by E.Wang
;      9-17-09



;
; Without an input file, look for self.fileid".run.out file in current directory 
;
CD, CURRENT=current_working_directory

separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'
subNames = STRSPLIT(self.fileID, separator, /EXTRACT)		; Crack 'fileID' on separators
nSubNames = N_ELEMENTS(subNames)
defaultName = subNames[nSubNames-1]
defaultDir = ''
IF nSubNames GT 1 THEN BEGIN
    FOR i = 0,nSubNames-2 DO BEGIN
        defaultDir = defaultDir + separator + subNames[i]
    ENDFOR
CD, defaultDir
ENDIF

idLen = STRLEN(defaultName)
runOut = STRMID(defaultName, 0, idLen-3) + ".run.out"
commandLine = "ok = FINDFILE('" + runOut + "', count=nfiles)"
ok = EXECUTE(commandLine)
;ok =FINDFILE(runOut, count=nfiles)
IF (nfiles EQ 1) THEN BEGIN 
filename = runOut
GOTO,found
ENDIF


IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out', Get_Path=NetCDF_DIR) ELSE fileIn=fileName
ok = FINDFILE(fileIn, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ; 	we didn't find a file... 
    MESSAGE, "Couldn't find run.out file", /INFORMATIONAL
    RETURN, 0
ENDIF

found:
output = 0

GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening run.out file:"
    PRINT, "Error Message = ", !ERR_STRING
    RETURN, output
ENDIF

output = { Name : "Gyro_ReadParameters" }

thisLine = STRARR(1)
READF, ioUnit, thisLine  ;=======
READF, ioUnit, thisLine  ; GYRO FRANKLIN

READF, ioUnit, thisLine  ;(C) 
READF, ioUnit, thisLine  ; =====
READF, ioUnit, thisLine  ; Platform
READF, ioUnit, thisLine  ; Machine ID
READF, ioUnit, thisLine  ; [GYRO MPI tasks
READF, ioUnit, thisLine  ; [Parsing data
blankCount = 0


groupName = "Methods"
commandLine = "thisStr = { myName: '" + groupName +"'}"
OK = EXECUTE(commandLine)

NewGroup : 	groupCount=3
ReadLine :	IF( EOF(ioUnit) ) THEN GOTO, Done
READF, ioUnit, thisLine


shortString = STRCOMPRESS(thisLine)
IF( STRCMP(shortString, "") ) THEN BEGIN ; line is blank
    blankCount=blankCount+1
    GOTO, ReadLine
ENDIF      
          
INFO = STRSPLIT(shortString, "INFO") ; info line
IF INFO[0] EQ 4 THEN BEGIN
    GOTO, ReadLine
ENDIF

Note = STRSPLIT(shortString, "Note", count=nStrings) ;note line
IF nStrings GT 8  THEN BEGIN
    ;PRINT, "got it"
    GOTO, ReadLine
ENDIF
;print, shortString, " and ", nStrings

out = STRSPLIT(thisLine, " out")
IF out[0] EQ 4 THEN BEGIN
    GOTO,ReadLine
ENDIF

explicit = STRSPLIT(thisLine, " exp") ;skip line FOR NOW
IF explicit[0] EQ 4 THEN BEGIN
 ;   PRINT, "got it"
    GOTO, ReadLine
ENDIF

explicit2 = STRSPLIT(thisLine, " EXP") ;skip line FOR NOW
IF explicit2[0] EQ 4 THEN BEGIN
 ;   PRINT, "GOT IT"
    GOTO, ReadLine
ENDIF

explicit3 = STRSPLIT(thisLine, " off") ;skip line FOR NOW
IF explicit3[0] EQ 4 THEN BEGIN
  ;  PRINT, "Got It"
    GOTO, ReadLine
ENDIF

nonzero1 = STRSPLIT(thisLine, " nonzeros") ;skip line FOR NOW
IF nonzero1[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero2 = STRSPLIT(thisLine, " values") ;skip line FOR NOW
IF nonzero2[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero3 = STRSPLIT(thisLine, " indices") ;skip line FOR NOW
IF nonzero3[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF

nonzero4 = STRSPLIT(thisLine, " iterations") ;skip line FOR NOW
IF nonzero4[0] EQ 17 THEN BEGIN 
    GOTO, ReadLine
ENDIF


quasi = STRSPLIT(thisLine, " n_i") ; skip quasineutrality line
IF quasi[0] EQ 4 THEN GOTO, ReadLine

please = STRSPLIT(thisLIne, " PLE") ; end of file, with Please see comment
IF please[0] EQ 4 THEN GOTO, DONE


;
; check for line full of -------------------- only
;
ignore = STRSPLIT(thisLine, "-", /EXTRACT, Count=nStrings)
ignore2 = STRSPLIT(ignore[0], " pro")
IF (ignore2[0] EQ 4) THEN BEGIN 
    GOTO, Skip 
ENDIF
IF( STRCMP(ignore[0], "") ) THEN BEGIN ; line is only -------
    GOTO, ReadLine
ENDIF 

Skip: ;skip ignore ------ command for s-alpha boundary

;
; Start group processing
;
separate = STRSPLIT(shortString, "--", Count=gr, /EXTRACT)
IF(gr GT 1) THEN BEGIN    ;First type of group  ----- name -----
    ;PRINT, "FOUND GROUP ", separate
    splitLine = STRSPLIT(thisLine, ":", /EXTRACT, Count=nStrings)
    IF(nStrings EQ 2) THEN GOTO, notgroup
    IF (groupCount GT 1) THEN BEGIN ;not the first group, so close the group
        commandLine = "output = CREATE_STRUCT(output, '" + groupName + "',thisStr)"
        OK = EXECUTE(commandLine)
    ENDIF
    gSegs = n_elements(separate)
    groupName = separate[1]
    IF gSegs EQ 3 THEN groupName = groupName + " " + separate[2]
    gLength = STRLEN(groupName)
    groupName = STRMID(groupName,1,gLength-2)
    ;
    ; handle parenthesis in group name
    ; 
    groupParens = STRSPLIT(groupName, "(", /EXTRACT)
    parenLength = n_elements(groupParens)
    IF parenLength NE 1 THEN BEGIN 
        groupName = groupParens[0] +groupParens[1]
        groupParens = STRSPLIT(groupName, ")", /EXTRACT)        
        groupName = groupParens[0] 
    ENDIF
    ;
    ; handle commas in group name
    ;
    comma = STRSPLIT(groupName, ",", /EXTRACT)
    commaLength = N_elements(comma)
    IF commaLength EQ 3 THEN     groupName = comma[0] + " " + comma[1] + " " + comma[2]  
    commandLine = "thisStr = { myName: '" + groupName +"'}"
    OK = EXECUTE(commandLine)
    groupCount = groupCount + 1
    GOTO, ReadLine
ENDIF
;
; Input a value with name.  
;
splitLine = STRSPLIT(thisLine, ":", /EXTRACT, Count=nStrings)
notgroup: ; Should only get here if there is a value to be stored.
IF(nStrings EQ 2) THEN BEGIN 
                                ; Check for multiple, space-delimited values
    pieces = STRSPLIT(splitLine[1], COUNT=nPieces, /EXTRACT)
    IF(nPieces GT 1) THEN splitLine[1] = STRJOIN(pieces, ' ')
    ;OK = EXECUTE("value ="+ splitLine[1])
    ;IF(OK EQ 0) THEN BEGIN
    ;    commandLine = "value = '"+ splitLine[1] + "'"
    ;    OK = EXECUTE(commandLine)
    ;ENDIF
    value = splitLine[1]
    period = STRSPLIT(thisLine," .", Count=nStrings)

    varName = STRCOMPRESS(splitLine[0], /REMOVE_ALL) 
    ;
    ; check for period at beginning of name
    ;
    IF (period[0] EQ 2) THEN BEGIN
        varname = Strmid(varName, 1)        
    ENDIF
    ;
    ; check for period in middle of name
    ;
    periodMid = STRSPLIT(varName, ".",Count = nStrings, /EXTRACT)  
    IF (nStrings GT 1) THEN varName = periodMid[0] + "_" + periodMid[1]
    ;
    ; check for period at end
    ;    
    periodEnd = STRPOS(varName, ".") ; check for period at end
    varL = STRLEN(periodEnd)
    IF (periodEnd GT 0) THEN varName = STRMID(varName, 0, varL-1)
    ;
    ; check for =>
    ;
    arrow = STRSPLIT(varName, "=>", Count=nStrings)
    IF (arrow[0] GT 1) THEN varname = STRMID(varname,2)    
    ;
    ; replace "/" with "over"
    ;
    over = STRSPLIT(varName, "/", /EXTRACT,Count=nStrings)
    IF (nStrings GT 1) THEN varName = over[0] + " over " + over[1]
    ;
    ; replace "asdf(fdsa)" with "asdf_of_fdsa"
    ;
    parens = STRSPLIT(varName, "(", /EXTRACT)
    parenLength = n_elements(parens)
    IF parenLength NE 1 THEN BEGIN 
        varName = parens[0] + " of " + parens[1]
        varL = STRLEN(varName)
        varName = STRMID(varName, 0, varL-1)
    ENDIF
    ;
    ; replace "(asdf)"with "asdf"
    ;
    fullParens = STRSPLIT(varName, "(")
    varLen= STRLEN(varName)
    IF fullParens[0] EQ 1 THEN varName = STRMID(varName,1,varLen-2)
    ;
    ; replace "a^n" with "a_exp_n"
    ;
    exponent = STRSPLIT(varName, "^", count= nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = exponent[0] + "exp" + exponent[1]
    ;
    ; replace "asdf-fdsa" with "asdf_fdsa"
    ;
    dash = STRSPLIT(varName, "-", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = dash[0] + "_" + dash[1]
    ;
    ; replace "asdf*fdsa" with "asdf_times_fdsa"
    ; 
    mult = STRSPLIT(varName, "*", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varname = mult[0] + "_times_" + mult[1]
    ;
    ; handle "asdf #1" situations to "asdf number 1"
    ;
    sharp = STRSPLIT(varName, "#", count=nStrings, /EXTRACT)
    IF (nStrings GT 1) THEN varName = sharp[0] + "_num_" + sharp[1]
    ;
    ; For "# Ion1" ignore
    ;
    sharp2 = STRSPLIT(varName, '#')
    IF (sharp2[0] EQ 1) THEN GOTO, ReadLine    
    ;
    ; done editing varName.  ready to save value!
    ;
    ;PRINT, typeOf(value), " ", value, "command line ", commandline
    IF (TypeOf(value) EQ 7) THEN BEGIN 
        valueShort =  STRCOMPRESS(value)
        valueCount = N_ELEMENTS(valueShort)
        value = valueShort[0]
        FOR i = 1,valueCount-1 DO BEGIN 
            value = value + " " + valueShort[ir]        
            ;PRINT, "value = ", value
        ENDFOR
        commandLine = "thisStr = CREATE_STRUCT(thisStr, '" + varName + "', '" + value +"')"  
    ENDIF ELSE commandLine = "thisStr = CREATE_STRUCT(thisStr, '" + varName + "', value)"
    
    OK = EXECUTE(commandLine)
    IF(STRCMP(varName, "n_spec", /FOLD_CASE)) THEN nSpecies=value    
ENDIF

GOTO, ReadLine

;
; EOF reached
;
DONE:

;  Need to close last group
commandLine = "output = CREATE_STRUCT(output, '" + groupName + "',thisStr)"
OK = EXECUTE(commandLine)

FREE_LUN, ioUnit                

CD, current_working_directory

RETURN, output

END                             ;  ****** GKVs1D::Gyro_Params ******
