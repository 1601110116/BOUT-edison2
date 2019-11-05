FUNCTION poincare_movie, _Extra=Extra

;
; Designed to read in directories of poincare surface of section plots
; Arguments:  
;     start-     initial time point to start at.  Default=1 
;
;     irange-     Max range of time to be sampled
;
;     iskip-      integer for how many time slices skip between plots.  
;                default = 1
;
;     sourcedir- directory structure should be "prefix#", with # being
;                the time slice of surface of section plot.  In such a
;                case, sourcedir should be set = "prefix"
;


CD, CURRENT=current_working_directory

; check _Extra for 'range' argument
irange = 10
result = GETKEYWORD('irange',Extra)
IF(QUERY_INTEGER(result)) THEN irange = result

; check _Extra for 'sourcedir' argument.  String should read 'G#t'
predir = "G8t"
result=GetKeyWord('sourcedir', Extra)
IF(QUERY_STRING(result) ) THEN predir = result

iskip = 1
result=GetKeyWord('iskip', Extra)
IF(QUERY_INTEGER(result)) THEN iskip = result
 
start = 1
result=GetKeyWord('start',Extra)
IF(QUERY_INTEGER(result)) THEN start = result 

range= irange/iskip


outputarr = fltarr(range+1)
outputerror = fltarr(range+1)
skiparr = fltarr(range+1)
surf = "/surface.out"



count = 0
skipdir = 0
output = fltarr(100,2,3000, range+1)

FOR i_time = 1,range+1 DO BEGIN

dir = predir + STRCOMPRESS(STRING(start-1+i_time*iskip), /REMOVE_ALL)

PRINT, "reading directory ", dir
surf = "/surface.out"

surf2 = dir + surf

ok = FINDFILE(surf2, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ;   we didn't find a file... 
    MESSAGE, "Couldn't find surf.out file", /INFORMATIONAL
    skiparr[i_time-1] = 1
    skipdir = skipdir + 1
    GOTO, noFile
ENDIF

fileName=surf2

found:
;output = 0



GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening surface.out file:"
    PRINT, "Error Message = ", !ERR_STRING
     skiparr[i_time-1] = 1
     skipdir = skipdir + 1
    GOTO, NoFile
ENDIF
  

;skiparr[count]=i_time
count = count + 1

IF( EOF(ioUnit) ) THEN GOTO, noFile
READF, ioUnit, thisLine

thisTheta = 0
iTheta = 0
i =  0
n_theta = thisLine + 0.
n_iterations = 3000

thisLine = STRARR(1)


;===============================   reading file now =======================
ReadLine :  IF( EOF(ioUnit) ) THEN GOTO, DoneReading
READF, ioUnit, thisLine

EqualSplit = STRSPLIT(thisLine, "=", /Extract, count=nStrings)

IF (nStrings EQ 2) THEN BEGIN
    thisTheta = EqualSplit[1] + 0.
    iTheta = iTheta + 1
    i=0
    GOTO, ReadLine
ENDIF


shortString = STRCOMPRESS(thisLine)
IF( STRCMP(shortString, "") ) THEN BEGIN ; line is blank
    blankCount=blankCount+1
    GOTO, ReadLine
ENDIF      

thisR = strmid(thisLine, 1,12) + 0.
thisPhi = strmid(thisLine, 14,12) + 0.
IF (thisR LT -.5) THEN BEGIN
    thisR = thisR + 1
ENDIF
IF (thisR GT .5) THEN BEGIN
    thisR = thisR - 1
ENDIF
output[iTheta,0,i,i_time-1] = thisR
output[iTheta,1,i,i_time-1] = thisPhi



i = i+1
;itheta = itheta +1
IF (i EQ n_iterations) THEN BEGIN
itheta = itheta+1
i = 0
ENDIF

GOTO, ReadLine

DoneReading:  

FREE_LUN, ioUnit                

CD, current_working_directory

noFile:
ENDFOR ; i_time loop




outputarr = objarr(100,range+1-skipdir)
i_insert=0
for i_time = 0,range DO BEGIN 
  IF (skiparr[i_time] EQ 1) THEN GOTO, noObject
  FOR  i = 0,99 DO BEGIN

  gridStructure = {Grid}
  gridStructure.Mnemonic = "r"
  gridStructure.Title = "r"
  gridStructure.units = "!4q!x"
  gridValues= output[i,0,*,i_time]
  gridValues = REFORM(gridvalues)
  gridStructure.Values= PTR_NEW(gridValues)
  nPoints = N_ELEMENTS(gridValues)
  gridStructure.Range=[-.5, .5]
  gridStructure.irange=[0,nPoints-1]


  objStructure = {GKVs1D}
  objStructure.Grid1 = gridStructure
  objStructure.mnemonic = STRMID(predir,0,STRLEN(predir)-1)+'_Poincare_Fieldmap'
  objStructure.title = 'Poincare Fieldmap'
  objStructure.indices = PTR_NEW('*')
  objStructure.units = ''
  objStructure.codename = 'Gyro'
  objStructure.CodePI = 'J. Candy'
  objValues = output[i,1,*,i_time]
  objValues = REFORM(objValues)
  objStructure.values = PTR_NEW(objValues)
  objStructure.vrange = [MIN(objValues), MAX(objValues)]
  obj = obj_new("gkvs1d", objstructure)
  
  outputarr[i,i_insert] = obj
  

  ENDFOR  ;i  
i_insert = i_insert+1
noObject:  
  
ENDFOR  ;i_time
RETURN, outputarr
END
