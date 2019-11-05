FUNCTION Gyro_Poincare, filename, _Extra=extra

;
; Routine to read in delta B_r from surface.dat file
;
; Argument:
;           filename- surface.out file location.  Without argument
;                     will invoke DIALOG_PICKFILE()
;           nturns- number of turns per fieldline.  Default = 3000          
; 
; 
; Written by E.Wang
;      9-3-09
;      10-5-09  Updated to handle new surf.out format



;
; Without an input file, look for self.fileid".run.out file in current directory 
;
CD, CURRENT=current_working_directory

IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out', Get_Path=NetCDF_DIR) ELSE fileIn=fileName
ok = FINDFILE(fileIn, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ; 	we didn't find a file... 
    MESSAGE, "Couldn't find surf.out file", /INFORMATIONAL
    RETURN, 0
ENDIF


separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'
;subNames = STRSPLIT(self.fileID, separator, /EXTRACT)		; Crack 'fileID' on separators
;nSubNames = N_ELEMENTS(subNames)
;defaultName = subNames[nSubNames-1]
;defaultDir = ''
;IF nSubNames GT 1 THEN BEGIN
;    FOR i = 0,nSubNames-2 DO BEGIN
;        defaultDir = defaultDir + separator + subNames[i]
;    ENDFOR
;CD, defaultDir
;ENDIF

;idLen = STRLEN(defaultName)
;runOut = STRMID(defaultName, 0, idLen-3) + ".surf.dat"
;commandLine = "ok = FINDFILE('" + runOut + "', count=nfiles)"
;ok = EXECUTE(commandLine)
;ok =FINDFILE(runOut, count=nfiles)
;IF (nfiles EQ 1) THEN BEGIN 
;filename = runOut
;GOTO,found
;ENDIF



fileName=fileIN

found:
output = 0

GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening run.out file:"
    PRINT, "Error Message = ", !ERR_STRING
    RETURN, output
ENDIF
  
output = { Name : "Delta_b" }
IF( EOF(ioUnit) ) THEN RETURN,0
READF, ioUnit, thisLine

thisTheta = 0
iTheta = 0
i =  0
n_theta = thisLine + 0.
n_iterations = 3000
output = fltarr(n_theta,2,n_iterations)
output2 = fltarr(n_theta,n_iterations)
output3 = fltarr(n_theta,n_iterations)
thisLine = STRARR(1)
mycount=  0
test = GetKeyWord('nturns', extra)
IF (TypeOf(test) NE 7) THEN n_iterations=test
;===============================   reading file now =======================
ReadLine :	IF( EOF(ioUnit) ) THEN GOTO, Done
READF, ioUnit, thisLine
;print, mycount, thisLine
mycount = mycount+1
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
;thisB = strmid(thisLine, 28,12) + 0.
thisB = strmid(thisLine,40,12) + 0.   ;MODIFYING FOR THETA, used to be 40
thisB_T = strmid(thisLine, 53,12) + 0.
thisR2 = thisR + 0.5
thisR3 = thisR2 mod 1 + (thisR2 LT 0)
thisRfinal = thisR3 - 0.5

IF (thisR LT -.5) THEN BEGIN
    thisR = thisR + 1
ENDIF
IF (thisR GT .5) THEN BEGIN
    thisR = thisR - 1
ENDIF
output[iTheta,0,i] = thisRfinal
output[iTheta,1,i] = thisPhi
output2[iTheta,i] = thisB
output3[iTheta,i] = thisB_T
;output[iTheta,2,i] = thisPhi
;output[iTheta,3,i] = thisB
;output[iTheta,4,i] = thisTheta



i = i+1
;itheta = itheta +1
IF (i EQ n_iterations) THEN BEGIN
itheta = itheta+1
i = 0
ENDIF

GOTO, ReadLine

;
; EOF reached
;
DONE:

FREE_LUN, ioUnit                

CD, current_working_directory

;
; Compute correlation function by fourier transform method
;

;arr = fltarr(n_iterations)

;FOR i=0,n_theta-1 DO BEGIN
;arr = output[i,*,*]
;output2[i,*] = FFT(arr,-1)
;ENDFOR



n_steps = n_iterations/n_theta

objarr = OBJARR(n_theta,3)

;fullGridValues = REFORM(output[0,0,*])
;gridOutput = REFORM(fullGridValues, n_theta, n_steps)

;fullObjValues = REFORM(output[0,1,*])
;objOutput = REFORM(fullObjValues,n_theta, n_steps)


FOR i = 0,n_theta -1 DO BEGIN

gridStructure = {Grid}
gridStructure.Mnemonic = "r"
gridStructure.Title = "r"
gridStructure.units = "!4q!x"
gridValues= output[i,0,*]
gridValues = REFORM(gridvalues)
gridStructure.Values= PTR_NEW(gridValues)
nPoints = N_ELEMENTS(gridValues)
gridStructure.Range=[-.5, .5]
gridStructure.irange=[0,nPoints-1]

objStructure = {GKVs1D}
objStructure.Grid1 = gridStructure
objStructure.mnemonic = 'Poincare Fieldmap'
objStructure.title = 'Poincare Fieldmap'
objStructure.indices = PTR_NEW('*')
objStructure.units = ''
objStructure.codename = 'Gyro 8.1'
objStructure.CodePI = 'J. Candy'
;objStructure.fileid = "!4b!x!de!n=0.1%"
;objStructure.runid = "cyclone-like !4b!x scan"
objValues = output[i,1,*]
objValues = REFORM(objValues)
objStructure.values = PTR_NEW(objValues)
objStructure.vrange = [MIN(objValues), MAX(objValues)]
obj = obj_new("gkvs1d", objstructure)
objarr[i,0] = obj

gridStructure2 = {Grid}
gridStructure2.Mnemonic = "Turns"
gridStructure2.Title = "Turns"
gridStructure2.Values = PTR_NEW(INDGEN(n_iterations))
nPoints = n_iterations
gridStructure2.Range= [0,nPoints-1]
gridStructure2.irange=[0,nPoints-1]

objStructure2 = {GKVs1D}
objStructure2.Grid1 = gridStructure2
objStructure2.mnemonic = 'B_y'
objStructure2.title = 'B_y'
objStructure2.indices = PTR_NEW('*')
objStructure2.units = 'put in units'
objStructure2.values = PTR_NEW(output2[i,*])
objStructure2.vrange = [MIN(output2[i,*]), MAX(output2[i,*])]
obj2 = obj_new("gkvs1d", objstructure2)
objarr[i,1] = obj2


objStructure3 = {GKVs1D}
objStructure3.Grid1 = gridStructure2
objStructure3.mnemonic = 'B_r'
objStructure3.title = 'B_r'
objStructure3.indices = PTR_NEW('*')
objStructure3.units = 'put in units'
objStructure3.values = PTR_NEW(output3[i,*])
objStructure3.vrange = [MIN(output3[i,*]), MAX(output3[i,*])]
obj3 = obj_new("gkvs1d", objstructure3)
objarr[i,2] = obj3

ENDFOR


RETURN, objarr

END                             ;  ****** Gyro_Poincare ******


PRO GKVs1D::Poincare_norm, xLength, _Extra=Extra
;
; Maps radial grid of poincare GKV object from [-.5,.5] to [0,xLength]
;
; yLength-  Optional argument to scale y axis
;



xRange = self.grid1.range
yRange = self.vrange
IF (xrange[0] NE -.5) OR (xrange[1] NE .5) THEN BEGIN
    PRINT, "radial range is not correct"
RETURN
ENDIF

self -> set, axis=1, range=[0,xlength]
xPtr = *(self.grid1.values)
xpts = N_ELEMENTS(xptr)

FOR i = 0, xpts-1 DO BEGIN
xptr[i] = (xptr[i] + .5) * xLength
ENDFOR



yLength=0
result=GetKeyWord('yLength', Extra)
IF(TYPEOF(result) NE 7 ) THEN yLength=result


IF (yLength EQ 0) THEN BEGIN
self -> set, vrange=[-!pi,!pi]
valuePtr = *(self.values)
nVals = N_ELEMENTS(valuePtr)

FOR i = 0, nVals-1 DO BEGIN
valuePtr[i] = (valuePtr[i]) * 2* !pi
ENDFOR
ENDIF ELSE BEGIN
self -> set, vrange=[0,yLength]
valuePtr=*(self.values)
nvals = N_ELEMENTS(valuePtr)

for i = 0,nVals-1 DO BEGIN
valuePtr[i] = (valuePtr[i] +.5)*yLength
ENDFOR
ENDELSE


PTR_FREE, self.grid1.values
PTR_FREE, self.values
self.grid1.values = PTR_NEW(xptr)
self.values = PTR_NEW(valuePtr)

RETURN
END ;  Poincare_Norm

FUNCTION Gyro_Poincare2, filename, _Extra=extra

;
; Routine to read in delta B_r from surface.dat file
;
; Argument:
;           filename- surface.out file location.  Without argument
;                     will invoke DIALOG_PICKFILE()
;           nturns- number of turns per fieldline.  Default = 3000   
;           niter- number of iterations for 1 turn       
; 
; 
; Written by E.Wang
;      9-3-09
;      10-5-09  Updated to handle new surf.out format



;
; Without an input file, look for self.fileid".run.out file in current directory 
;
CD, CURRENT=current_working_directory

IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out', Get_Path=NetCDF_DIR) ELSE fileIn=fileName
ok = FINDFILE(fileIn, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ;   we didn't find a file... 
    MESSAGE, "Couldn't find surf.out file", /INFORMATIONAL
    RETURN, 0
ENDIF


separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'
;subNames = STRSPLIT(self.fileID, separator, /EXTRACT)    ; Crack 'fileID' on separators
;nSubNames = N_ELEMENTS(subNames)
;defaultName = subNames[nSubNames-1]
;defaultDir = ''
;IF nSubNames GT 1 THEN BEGIN
;    FOR i = 0,nSubNames-2 DO BEGIN
;        defaultDir = defaultDir + separator + subNames[i]
;    ENDFOR
;CD, defaultDir
;ENDIF

;idLen = STRLEN(defaultName)
;runOut = STRMID(defaultName, 0, idLen-3) + ".surf.dat"
;commandLine = "ok = FINDFILE('" + runOut + "', count=nfiles)"
;ok = EXECUTE(commandLine)
;ok =FINDFILE(runOut, count=nfiles)
;IF (nfiles EQ 1) THEN BEGIN 
;filename = runOut
;GOTO,found
;ENDIF



fileName=fileIN

found:
output = 0

GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening run.out file:"
    PRINT, "Error Message = ", !ERR_STRING
    RETURN, output
ENDIF
  
output = { Name : "Delta_b" }
IF( EOF(ioUnit) ) THEN RETURN,0
READF, ioUnit, thisLine

thisTheta = 0
iTheta = 0
i =  0
n_theta = thisLine + 0.
n_iterations = 3000
thisLine = STRARR(1)
mycount=  0
n_i=1
test = GetKeyWord('nturns', extra)
IF (TypeOf(test) NE 7) THEN n_iterations=test
test = GetKeyWord('niter', extra)
IF (TypeOf(test) NE 7) THEN n_i=test
output = fltarr(n_theta,2,n_iterations/n_i)
output2 = fltarr(n_theta,n_iterations)
output3 = fltarr(n_theta,n_iterations)
;===============================   reading file now =======================
ReadLine :  IF( EOF(ioUnit) ) THEN GOTO, Done
READF, ioUnit, thisLine
;print, mycount, thisLine
mycount = mycount+1
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
;thisB = strmid(thisLine, 28,12) + 0.
thisB = strmid(thisLine,40,12) + 0.   ;MODIFYING FOR THETA, used to be 40
thisB_T = strmid(thisLine, 53,12) + 0.
thisR2 = thisR + 0.5
thisR3 = thisR2 mod 1 + (thisR2 LT 0)
thisRfinal = thisR3 - 0.5

IF (thisR LT -.5) THEN BEGIN
    thisR = thisR + 1
ENDIF
IF (thisR GT .5) THEN BEGIN
    thisR = thisR - 1
ENDIF

IF ( i mod n_i EQ n_i-1) THEN BEGIN
output[iTheta,0,i/n_i] = thisRfinal
output[iTheta,1,i/n_i] = thisPhi
ENDIF


output2[iTheta,i] = thisB
output3[iTheta,i] = thisB_T
;output[iTheta,2,i] = thisPhi
;output[iTheta,3,i] = thisB
;output[iTheta,4,i] = thisTheta



i = i+1
;itheta = itheta +1
IF (i EQ n_iterations) THEN BEGIN
itheta = itheta+1
i = 0
ENDIF

GOTO, ReadLine

;
; EOF reached
;
DONE:

FREE_LUN, ioUnit                

CD, current_working_directory

;
; Compute correlation function by fourier transform method
;

;arr = fltarr(n_iterations)

;FOR i=0,n_theta-1 DO BEGIN
;arr = output[i,*,*]
;output2[i,*] = FFT(arr,-1)
;ENDFOR



n_steps = n_iterations/n_theta

objarr = OBJARR(n_theta,3)

;fullGridValues = REFORM(output[0,0,*])
;gridOutput = REFORM(fullGridValues, n_theta, n_steps)

;fullObjValues = REFORM(output[0,1,*])
;objOutput = REFORM(fullObjValues,n_theta, n_steps)


FOR i = 0,n_theta -1 DO BEGIN

gridStructure = {Grid}
gridStructure.Mnemonic = "r"
gridStructure.Title = "r"
gridStructure.units = "!4q!x"
gridValues= output[i,0,*]
gridValues = REFORM(gridvalues)
gridStructure.Values= PTR_NEW(gridValues)
nPoints = N_ELEMENTS(gridValues)
gridStructure.Range=[-.5, .5]
gridStructure.irange=[0,nPoints-1]

objStructure = {GKVs1D}
objStructure.Grid1 = gridStructure
objStructure.mnemonic = 'Poincare Fieldmap'
objStructure.title = 'Poincare Fieldmap'
objStructure.indices = PTR_NEW('*')
objStructure.units = ''
objStructure.codename = 'Gyro Fieldmap'
objStructure.CodePI = 'J. Candy'
;objStructure.fileid = "!4b!x!de!n=0.1%"
;objStructure.runid = "cyclone-like !4b!x scan"
objValues = output[i,1,*]
objValues = REFORM(objValues)
objStructure.values = PTR_NEW(objValues)
objStructure.vrange = [MIN(objValues), MAX(objValues)]
obj = obj_new("gkvs1d", objstructure)
objarr[i,0] = obj

gridStructure2 = {Grid}
gridStructure2.Mnemonic = "Turns"
gridStructure2.Title = "Turns"
gridStructure2.Values = PTR_NEW(INDGEN(n_iterations))
nPoints = n_iterations
gridStructure2.Range= [0,nPoints-1]
gridStructure2.irange=[0,nPoints-1]

objStructure2 = {GKVs1D}
objStructure2.Grid1 = gridStructure2
objStructure2.mnemonic = 'B_y'
objStructure2.title = 'B_y'
objStructure2.indices = PTR_NEW('*')
objStructure2.units = 'put in units'
objStructure2.values = PTR_NEW(output2[i,*])
objStructure2.vrange = [MIN(output2[i,*]), MAX(output2[i,*])]
obj2 = obj_new("gkvs1d", objstructure2)
objarr[i,1] = obj2


objStructure3 = {GKVs1D}
objStructure3.Grid1 = gridStructure2
objStructure3.mnemonic = 'B_r'
objStructure3.title = 'B_r'
objStructure3.indices = PTR_NEW('*')
objStructure3.units = 'put in units'
objStructure3.values = PTR_NEW(output3[i,*])
objStructure3.vrange = [MIN(output3[i,*]), MAX(output3[i,*])]
obj3 = obj_new("gkvs1d", objstructure3)
objarr[i,2] = obj3

ENDFOR


RETURN, objarr

END                             ;  ****** Gyro_Poincare ******

FUNCTION GKVs1D::GENE_flip, _Extra=extra
;
; Takes poincare plot produced from gene data and 
;   1) Change values (x,y) to (x,-y) then shift yvalues by half box
;   2) OPTIONAL:  Adjust y axis range  (-L_y/2,L_y/2)
;   3) OPTIONAL:  Adjust x axis range (-L_x/2,L_y/2)
;
;   Arguments:  xlength, ylength- 
;
;  5-18-11 EW
;




xLength=1
result=GetKeyWord('xLength', Extra)
IF(TYPEOF(result) NE 7 ) THEN xLength=result

yLength=1
result=GetKeyWord('yLength', Extra)
IF(TYPEOF(result) NE 7 ) THEN yLength=result




result = self -> makecopy()

old_y_vals = *self.values

old_x_vals = *self.grid1.values

;new_x_vals = old_x_vals*(-1)
;PTR_FREE, (result.grid1.values)
;result.grid1.values=PTR_NEW(new_x_vals)

new_y_vals = old_y_vals*(-1)   ;flip y axis 
new_y_vals = new_y_vals +.5 - 1*ROUND(new_y_vals+.5) ; adjust box periodically by l_y/2
PTR_FREE, result.values
result.values=PTR_NEW(new_y_vals)

result -> scaleaxis, axis=1, const=xlength
result2 = result -> times(ylength)


RETURN, result2

END


