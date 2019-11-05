
FUNCTION GKV_ReadTUBEr, grdin, nm, Title=title, Units=units, 	$
			Yang=yang, iSkip=iSkip_in, 		$
			aScale=aScaleIn, rhoStar=rhoStar
;
;
; Get data from 'grdin'r.out
;
rOutFile = 'a' + grdin+'r.out'
rOutFile = STRCOMPRESS(rOutFile,/REMOVE_ALL)
result = FINDFILE(rOutFile, COUNT=ok)
IF(NOT ok) THEN BEGIN
	MESSAGE, " Couldn't find " + rOutFile, /INFORMATIONAL
	RETURN, 0
ENDIF

ind = FLTARR(4)

OPENR, lun, rOutFile, /GET_LUN
READF, lun, ind


im = LONG(ind[0])
;tm = LONG(ind[1])
tm = nm/ind[2]+1
;nplot=LONG(ind[2])
;tm=LONG(nm/nplot)+1
; ****** PATCH FOR YANG'S EM_BECHMARK OF  6/08 *********
IF KEYWORD_SET(Yang) THEN BEGIN
	km=LONG(ind[2])
	nm=LONG(ind[3])
	nPlot=LONG(ind[4])
	tm = LONG(nm/nplot)+1L
ENDIF
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
;	IF KEYWORD_SET(Yang) THEN tm=ind[2] 
; ****************************************************** 
iStart = 0L
iEnd = tm - 1L
CASE N_ELEMENTS(iSkip_in) OF
	0 :	iSkip = 1L
	1 :	iSkip = LONG(iSkip_in)
	2 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = 1L
		END
	3 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = LONG(iSkip_in[2])
		END
	ELSE :	BEGIN
			MESSAGE, "too many elements in iSkip", /INFORMATIONAL
			RETURN, 0
		END
ENDCASE

tLength = (iEnd+1L-iStart)/iSkip + 1
dt=ind[2]*iskip
timeStamp = STRARR(1)
holdgrd=FLTARR(im+1L)
aphir=FLTARR(im+1L,tLength)
timeStep = LONARR(tm)
FOR k=iStart, iEnd DO BEGIN
	READF, lun, timeStamp
	crackedTimeStamp = STR_SEP(timeStamp, "=")
	READS, crackedTimeStamp[1], fTimeStep
	timeStep[(k-iStart)/iSkip]=LONG(fTimeStep)
	READF, format='(e12.5)',lun, holdgrd
	IF( iSkip*((k-iStart)/iSkip) EQ (k-iStart) ) THEN	$
		aphir(*,(k-iStart)/iSkip) = holdgrd
	IF(EOF(lun))THEN GOTO, MakeObj
ENDFOR
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip
;
; Make a GKV object out of this data
;
ObjStr = {GKVs2D}				; Get a GKVs2D object structure
ObjStr.CodeName = 'TUBE'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker'		;	data.
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
	IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
; ****************************************************** 
;
; Set up the grid structure for the
; independent variable, time ('t').
;
;timex = dt*iStart + dt*iskip*FINDGEN(tLast)
timex = timeStep[0:tLast]
tGrid = {Grid}
tGrid.mnemonic = 't'
tGrid.title = 't'
tGrid.units = '1/!4X!X!Dci!N'
tGrid.values = PTR_NEW(timex)
tGrid.boundary = 'Open'
tGrid.uniform = 1B
tmin = MIN(timex, MAX=tmax)
tGrid.range = [tmin, tmax]
tGrid.irange = [0, tLast]
;
; set up the grid structure for the 
; independent variable, x
;
xValues = FINDGEN(im+1)
xGrid = {Grid}
xGrid.mnemonic = 'x'
xGrid.title = 'x'
xGrid.units = '!4q!X!Ds!N'
xGrid.values = PTR_NEW(xValues)
xGrid.boundary = 'Periodic'
IF((im MOD 2) EQ 0) THEN xGrid.boundary = 'Periodic (closed)'
xGrid.uniform = 1B
xMin = GKVsd_MIN(xvalues, MAX=xMax)
xGrid.range = [xMin, xMax]
xGrid.irange = [0, im]
;
;
ObjStr.mnemonic = grdin + '_avg'
objtitle = '!13<!X' + grdin + '!13>!4!Dw!X!N'
IF(N_ELEMENTS(title) NE 0) THEN objtitle='!13<!X' + title + '!13>!4!Dw!X!N'
ObjStr.title = objtitle
indices = REPLICATE('*', 2)
ObjStr.Indices = PTR_NEW(indices)
objUnits = 'T/e'
IF(N_ELEMENTS(units) NE 0) THEN objUnits = units
ObjStr.units = objUnits
ObjStr.values = PTR_NEW(aphir[*,0:tLast])
vMin= GKVsd_MIN(aphir[*,0:tLast], MAX=vMax)
ObjStr.vRange = [vmin, vmax]
;
ObjStr.Grid1 = xGrid
ObjStr.Grid2 = tGrid
;
result = OBJ_NEW('GKVs2D', ObjStr)
IF(N_ELEMENTS(rhoStar) EQ 1) THEN BEGIN
	rhoStar= rhostar < 1./rhoStar
	aScale = "L!DT!N"
	IF(N_ELEMENTS(aScaleIn) EQ 1) THEN aScale=aScaleIn
	result -> ScaleAxis, 't', const=rhoStar, units=aScale+"/c!Ds!N"
	temp = result -> over(rhoStar)
	result -> trash
	result=temp
	result -> get, units=units
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBEr ****** ;


FUNCTION GKV_ReadTUBExz, grdin, nm, Title=title, Units=units, 	$
			 Yang=yang, iSkip=iSkip_in, 		$
			 aScale=aScaleIn, rhoStar=rhoStar
;
; Get data from 'grdin'xz.out
;
xzOutFile = grdin+'xz.out'
xzOutFile = STRCOMPRESS(xzOutFile,/REMOVE_ALL)
result = FINDFILE(xzOutFile, COUNT=ok)
IF(NOT ok) THEN BEGIN
	MESSAGE, "Couldn't find " + xzOutFile, /INFORMATIONAL
	RETURN, 0
ENDIF

ind = FLTARR(5)

OPENR, lun, xzOutFile, /GET_LUN
READF, lun, ind


im = FIX(ind[0])
jm = FIX(ind[1])
;tm = FIX(ind[2])
nplot=FIX(ind[3])
tm = FIX(nm/nplot)+1
; ****** PATCH FOR YANG'S EM_BECHMARK OF  6/08 *********
IF (N_ELEMENTS(Yang) EQ 4) THEN BEGIN
	im = LONG(Yang[0])
	jm = LONG(Yang[2])
	nm = LONG(Yang[3])
	nPlot=LONG(ind[4])
	tm = LONG(nm/nplot)+1L
ENDIF
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
;	IF KEYWORD_SET(Yang) THEN tm=ind[2]
; ****************************************************** 
iStart = 0L
iEnd = tm - 1L
CASE N_ELEMENTS(iSkip_in) OF
	0 :	iSkip = 1L
	1 :	iSkip = LONG(iSkip_in)
	2 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = 1L
		END
	3 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = LONG(iSkip_in[2])
		END
	ELSE :	BEGIN
			MESSAGE, "too many elements in iSkip", /INFORMATIONAL
			RETURN, 0
		END
ENDCASE

tLength = (iEnd+1L-iStart)/iSkip + 1
dt=ind[3]*iskip
timeStamp = STRARR(1)
holdgrd=FLTARR(im+1L, jm+1)
grdq = FLTARR(im+1, jm+1, tm)
timeStep = LONARR(tm)
FOR k=iStart, iEnd DO BEGIN
	READF, lun, timeStamp
	crackedTimeStamp = STR_SEP(timeStamp, "=")
	READS, crackedTimeStamp[1], fTimeStep
	timeStep[(k-iStart)/iSkip]=LONG(fTimeStep)
	READF, format='(e12.5)',lun, holdgrd
	IF( iSkip*((k-iStart)/iSkip) EQ (k-iStart) ) THEN	$
		grdq(*,*,(k-iStart)/iSkip) = holdgrd
	IF(EOF(lun))THEN GOTO, MakeObj
ENDFOR
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip

;
; Make a GKV object out of this data
;
ObjStr = {GKVs3D}				; Get a GKVs3D object structure
ObjStr.CodeName = 'TUBE'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker'		;	data.
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
	IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
; ****************************************************** 
;
; Set up the grid structure for the
; independent variable, time ('t').
;
;timex = dt*iStart + dt*iskip*FINDGEN(tLast+1)
timex = timeStep[0:tLast]
tGrid = {Grid}
tGrid.mnemonic = 't'
tGrid.title = 't'
tGrid.units = '1/!4X!X!Dci!N'
tGrid.values = PTR_NEW(timex)
tGrid.boundary = 'Open'
tGrid.uniform = 1B
tmin = MIN(timex, MAX=tmax)
tGrid.range = [tmin, tmax]
tGrid.irange = [0, tLast]
;
; set up the grid structure for the 
; independent variable, x
;
xValues = FINDGEN(im+1L)
xGrid = {Grid}
xGrid.mnemonic = 'x'
xGrid.title = 'x'
xGrid.units = '!4q!X!Ds!N'
xGrid.values = PTR_NEW(xValues)
xGrid.boundary = 'Periodic'
IF((im MOD 2) EQ 0) THEN xGrid.boundary = 'Periodic (closed)'
xGrid.uniform = 1B
xMin = GKVsd_MIN(xvalues, MAX=xMax)
xGrid.range = [xMin, xMax]
xGrid.irange = [0, im]
;
; set up the grid structure for the 
; independent variable, z
;
zValues = FINDGEN(jm+1)/FLOAT(jm)-0.5
zGrid = {Grid}
zGrid.mnemonic = 'z'
zGrid.title = 'z'
zGrid.units = '2!4p!XqR'
zGrid.values = PTR_NEW(zValues)
zGrid.boundary = 'Periodic'
IF((jm MOD 2) EQ 0) THEN xGrid.boundary = 'Periodic (closed)'
zGrid.uniform = 1B
zMin = GKVsd_MIN(zvalues, MAX=zMax)
zGrid.range = [zMin, zMax]
zGrid.irange = [0, jm]
;
;
ObjStr.mnemonic = grdin + 'xz'
objtitle = grdin + 'xz'
IF(N_ELEMENTS(title) NE 0) THEN objtitle=title
ObjStr.title = objtitle
indices = REPLICATE('*', 3)
ObjStr.Indices = PTR_NEW(indices)
objUnits = 'T/e'
IF(N_ELEMENTS(units) NE 0) THEN objUnits = units
ObjStr.units = objUnits
ObjStr.values = PTR_NEW(grdq[*,*,0:tLast])
vMin= GKVsd_MIN(grdq, MAX=vMax)
ObjStr.vRange = [vmin, vmax]
;
ObjStr.Grid1 = xGrid
ObjStr.Grid2 = zGrid
ObjStr.Grid3 = tGrid
;
result = OBJ_NEW('GKVs3D', ObjStr)
IF(N_ELEMENTS(rhoStar) EQ 1) THEN BEGIN
	rhoStar= rhostar < 1./rhoStar
	aScale = "L!DT!N"
	IF(N_ELEMENTS(aScaleIn) EQ 1) THEN aScale=aScaleIn
	result -> ScaleAxis, 't', const=rhoStar, units=aScale+"/c!Ds!N"
	temp = result -> over(rhoStar)
	result -> trash
	result=temp
	result -> get, units=units
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBExz ****** ;


FUNCTION GKV_ReadTUBExy, grdin, nm, Title=title, Units=units, 	$
			 Yang=yang, iSkip=iSkip_in, 		$
			 aScale=aScaleIn, rhoStar=rhoStar
;
; Get data from 'grdin'xy.out
;
xyOutFile = grdin+'xy.out'
xyOutFile = STRCOMPRESS(xyOutFile,/REMOVE_ALL)
result = FINDFILE(xyOutFile, COUNT=ok)
IF(NOT ok) THEN BEGIN
	MESSAGE, "Couldn't find " + xyOutFile, /INFORMATIONAL
	RETURN, 0
ENDIF

ind = FLTARR(5)
timeStamp = STRARR(1)

OPENR, lun, xyOutFile, /GET_LUN
READF, lun, ind


im=LONG(ind[0])
jm=LONG(ind[1])
;tm=LONG(ind[2])
nplot=LONG(ind[3])
tm = LONG(nm/nplot)+1L
; ****** PATCH FOR YANG'S EM_BECHMARK OF  6/08 *********
IF(N_ELEMENTS(Yang) EQ 4) THEN BEGIN
	im = LONG(Yang[0])
	jm = LONG(Yang[1])	
	nm = LONG(YANG[3])
	nPlot=LONG(ind[4])
	tm = LONG(nm/nplot)+1L
ENDIF
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
;	IF KEYWORD_SET(Yang) THEN tm=ind[2]
; ****************************************************** 
iStart = 0L
iEnd = tm - 1L
CASE N_ELEMENTS(iSkip_in) OF
	0 :	iSkip = 1L
	1 :	iSkip = LONG(iSkip_in)
	2 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = 1L
		END
	3 :	BEGIN
			iStart = LONG(iSkip_in[0])
			iEnd   = LONG(iSkip_in[1])
			iSkip  = LONG(iSkip_in[2])
		END
	ELSE :	BEGIN
			MESSAGE, "too many elements in iSkip", /INFORMATIONAL
			RETURN, 0
		END
ENDCASE

tLength = (iEnd+1L-iStart)/iSkip + 1
dt = ind[3]*iskip
timeStamp = STRARR(1)
holdgrd=FLTARR(im+1L, jm+1)
grdq = FLTARR(im+1L, jm+1L, tm)
timeStep = LONARR(tm)
FOR k=iStart, iEnd DO BEGIN
	IF(EOF(lun))THEN GOTO, MakeObj
	READF, lun, timeStamp
	crackedTimeStamp = STR_SEP(timeStamp, "=")
	READS, crackedTimeStamp[1], fTimeStep
	timeStep[(k-iStart)/iSkip]=LONG(fTimeStep)
	READF, format='(e12.5)',lun, holdgrd
	IF( iSkip*((k-iStart)/iSkip) EQ (k-iStart) ) THEN	$
		grdq(*,*,(k-iStart)/iSkip) = holdgrd
ENDFOR
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip 
;
; Make a GKV object out of this data
;
ObjStr = {GKVs3D}				; Get a GKVs3D object structure
ObjStr.CodeName = 'TUBE'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker'			;	data.
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
	IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
; ****************************************************** 

;
; Set up the grid structure for the
; independent variable, time ('t').
;
;timex = dt*iStart + dt*iskip*FINDGEN(tLast+1)
timex = timeStep[0:tLast]
tGrid = {Grid}
tGrid.mnemonic = 't'
tGrid.title = 't'
tGrid.units = '1/!4X!X!Dci!N'
tGrid.values = PTR_NEW(timex)
tGrid.boundary = 'Open'
tGrid.uniform = 1B
tmin = MIN(timex, MAX=tmax)
tGrid.range = [tmin, tmax]
tGrid.irange = [0, tLast]
;
; set up the grid structure for the 
; independent variable, x
;
xValues = FINDGEN(im+1L)
xGrid = {Grid}
xGrid.mnemonic = 'x'
xGrid.title = 'x'
xGrid.units = '!4q!X!Ds!N'
xGrid.values = PTR_NEW(xValues)
xGrid.boundary = 'Periodic'
IF((im MOD 2) EQ 0) THEN xGrid.boundary = 'Periodic (closed)'
xGrid.uniform = 1B
xMin = GKVsd_MIN(xvalues, MAX=xMax)
xGrid.range = [xMin, xMax]
xGrid.irange = [0, im]
;
; set up the grid structure for the 
; independent variable, x
;
yValues = FINDGEN(jm+1)
yGrid = {Grid}
yGrid.mnemonic = 'y'
yGrid.title = 'y'
yGrid.units = '!4q!X!Ds!N'
yGrid.values = PTR_NEW(yValues)
yGrid.boundary = 'Periodic'
IF((jm MOD 2) EQ 0) THEN yGrid.boundary = 'Periodic (closed)'
yGrid.uniform = 1B
yMin = GKVsd_MIN(yvalues, MAX=yMax)
yGrid.range = [ymin, ymax]
yGrid.irange = [0, jm]
;
;
ObjStr.mnemonic = grdin + 'xy'
objtitle = grdin + 'xy'
IF(N_ELEMENTS(title) NE 0) THEN objtitle=title
ObjStr.title = objtitle
indices = REPLICATE('*', 3)
ObjStr.Indices = PTR_NEW(indices)
objUnits = 'T/e'
IF(N_ELEMENTS(units) NE 0) THEN objUnits = units
ObjStr.units = objUnits
ObjStr.values = PTR_NEW(grdq[*,*,0:tLast])
vMin= GKVsd_MIN(grdq, MAX=vMax)
ObjStr.vRange = [vmin, vmax]
;
ObjStr.Grid1 = xGrid
ObjStr.Grid2 = yGrid
ObjStr.Grid3 = tGrid
;
result = OBJ_NEW('GKVs3D', ObjStr)
IF(N_ELEMENTS(rhoStar) EQ 1) THEN BEGIN
	rhoStar= rhostar < 1./rhoStar
	aScale = "L!DT!N"
	IF( N_ELEMENTS(aScaleIn) EQ 1 ) THEN aScale=aScaleIn
	result -> ScaleAxis, 't', const=rhoStar, units=aScale+"/c!Ds!N"
	temp = result -> over(rhoStar)
	result -> trash
	result=temp
	result -> get, units=units
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBExy ****** ;




FUNCTION Tube_Data,	Path=startPath, TubePath=tubePath, nt=nm, Yang=yang,		$ 
			CodeName=Codename, CodePI=CodePi, RunId=RunID, 			$
			FileId=FIleId, Debug=debug, iSkip=iSkip_in, 			$
			aScale=aScaleIN, rhoStar=rhoStar
;
; Purpose:
;
;		This function reads output files from the 'TUBE' code
;		(S. Parker et al, U. of Co), creates GKV objects of 
;		appropriate dimensionality from the data contained 
;		within them, and returns these objects within an
;		anonymous structure.
;
; Arguments:
;
;		None
;
; Keywords:
;
;	Path		The path at which DIALOGUE_PICKFILE begins its search 
;			for the directory containing the TUBE output files.  
;			Defaults to the current working directory.
;			(Optional)
;
;	TubePath	The full and correct path to the directory containing
;			the TUBE output files.  If this keyword is set, then 
;			DIALOGUE_PICKFILE will not be used.
;			(Optional)
;
;	nt		Number of time steps the history file.  This is only
;			used if there is NOT a 'hist.out' file in the selected
;			directory.  
;
;	Yang		Set this keyword (i.e., put "/Yang" on the command line)
;			so that the third integer in the first record of the 
;			phi**.out, etc. files is interpreted as the number of time
;			slices in the file (instead of the interval between
;			timeslices in these files relative to the history
;			files).
;
;	Codename	Set this keyword to the name of the code which generated
;			the output files. Defaults to "TUBE".  (Optional)
;
;	CodePi		Set this this keyword to the name of the PI responsible
;			for this simulation data.  Defaults to "S.E. Parker"
;			(unless the keyword 'Yang' is set, in which case 'CodePi'
;			defaults to "Y. Chen").  (Optional)
;
;	RunId		Set this keyword to a terse phrase describing this run.
;			(Hershey vector font conventions can be employed).
;			No default, so that the RunID field on output will be blank
;			if this keyword is not set.  (Optional) 
;
;	FileId		Set this keyword to a terse phrase further describing this run.
;			(Hershey vector font conventions can be employed).
;			No default, so that the FileID field on output will be blank
;			if this keyword is not set.  (Optional)
;
;	iSkip		Controls which timeslices are read from the phixy.out, etc.
;			files:
;				iSkip = i		read every ith timeSlice
;				iSkip = [iStart, iEnd]	start at iStart, finish at iEND
;				iSkip = [iStart, iEnd, iSkip]
;							as above, except only read
;							every iSkipth timeSlice
;			Defaluts to iSkip = 1 (optional). 
;
;	aScale		String variable containing Hershey vector font representation
;			of the symbol to be used for the macroscopic length scale
;			in the gyrokinetic units used in labeling plots.  
;			Defaults to "L!DT!N" (that is, the tempergradient length).
;			Optional.
;
;	rhoStar		Value of rhoStar. If a value of rhoStar is supplied, it is
;			used to put both timebase and values of phi, etc into gyro
;			kinetic units. (Optional) 
;
;
;	Debug		Set this keyword to obtain relatively verbose behavior
;			which may be useful in debuggin this routine.
;			(Optional)
;
;
; Side Effects:
;
;			On return from Tube_Data the current working directory will be set
;			the directory containing the TUBE output files.
;
; Written by W.M. Nevins
;	7/25/00
; Revised by W.M. Nevins
;	1/18/06
; Revised by W.M. Nevins
;	2/10/07
;   Modified to use info from Hist.out file 
;   to estimate lenght of data in *xy.out, 
;   *xz.out, etc. files.
;
; Find directory containing TUBE output files
;
IF(N_ELEMENTS(tubePath) EQ 0) THEN BEGIN
	;
	; Use DIALOGUE_PICK_FILE to locate directory containing TUBE output files
	;
	CD, CURRENT=Current_Working_Directory
	path=Current_Working_Directory
	IF(N_ELEMENTS(startPath) NE 0) THEN path = startPath
	result = DIALOG_PICKFILE(/DIRECTORY, TITLE='GKV:  Select TUBE Output Directory', PATH=path, GET_PATH = tubePath)
ENDIF

;
; Set 'current working directory' to tubePath
;
CD, tubePath
;
; Get data from hist.out
;
histFile = 'hist.out'
result = FINDFILE(histFile, COUNT=ok)
IF(NOT ok) THEN BEGIN
	MESSAGE, "Couldn't find hist.out", /INFORMATIONAL
	output = {FirstTag:0}
	IF(TypeOf(nm) EQ 0) THEN BEGIN		; nm is undefined.  Must get value
		PRINT, "**** Enter number of time steps ******"
		READ, nm, PROMPT="nt="
	ENDIF
	ObjStr = {GKVs1D}				; Get a GKVs1D object structure
	ObjStr.CodeName = 'TUBE'			; 	and load in the common
	ObjStr.CodePI = 'S.E. Parker'		;	data.
	; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
		IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
	; ****************************************************** 
	IF KEYWORD_SET(CodeName) THEN Objstr.CodeName=CodeName
	IF KEYWORD_SET(CodePi  ) THEN Objstr.CodePI  =CodePi
	IF KEYWORD_SET(RunID   ) THEN Objstr.RunID   =RunId
	IF KEYWORD_SET(FileId  ) THEN Objstr.FileID  =FIleID
	GOTO, ReadXYout
ENDIF

ind=lonarr(3)	
; 
;	ind[0] = number of species
;	ind[1] = number of time steps
;	ind[2] = number of 'modes' in this history file
; 

OPENR, hist_unit, histFile, /GET_LUN
READF, hist_unit, ind, dt
modes = INTARR(3,ind[2])
READF, hist_unit, modes
IF(KEYWORD_SET(debug)) THEN PRINT, modes
nm=ind[1]
IF(KEYWORD_SET(debug)) THEN PRINT, "nm=", nm
pmodehis = COMPLEXARR(ind[2],ind[1])
te = FLTARR(ind[1])
ke = FLTARR(ind[0],ind[1])
pfl= FLTARR(ind[0],ind[1])
efl= FLTARR(ind[0],ind[1])
fe = FLTARR(ind[1])
rmsphi = FLTARR(ind[1])
avgWSq = FLTARR(ind[1])
timex  = FINDGEN(ind[1])*dt

holdpmode = COMPLEXARR(ind[2]*ind[1])
holdke  = FLTARR(ind[0],ind[1])
holdpfl = FLTARR(ind[0],ind[1])
holdefl = FLTARR(ind[0],ind[1])

READF, format='(2e12.5)', hist_unit, holdpmode
FOR mm=0L,ind[2]-1 DO BEGIN
	FOR n=0L,ind[1]-1 DO BEGIN
		pmodehis[mm,n] = holdpmode[n+mm*ind[1]]
  	ENDFOR
ENDFOR

READF, format='(e12.5)', hist_unit, holdke
READF, format='(e12.5)', hist_unit, holdpfl
READF, format='(e12.5)', hist_unit, holdefl
FOR m=0L,ind[0]-1 DO BEGIN
	FOR n=0L,ind[1]-1 DO BEGIN
		ke[m,n]  = holdke[n+m*ind[0]]
		pfl[m,n] = holdpfl[n+m*ind[0]]
		efl[m,n] = holdefl[n+m*ind[0]]
		te[n] = te[n] + ke[m,n]
	ENDFOR
ENDFOR

READF, format='(e12.5)', hist_unit, fe
READF, format='(e12.5)', hist_unit, rmsphi
;READF, format='(e12.5)', hist_unit, avgWsq

FOR n=0L,ind[1]-1 DO BEGIN
	te[n] = te[n]+fe[n]
ENDFOR

FREE_LUN, hist_unit

;
; Make GKV objects out of all this data
;
ObjStr = {GKVs1D}				; Get a GKVs1D object structure
ObjStr.CodeName = 'TUBE'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker'		;	data.
; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
	IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
; ****************************************************** 
IF KEYWORD_SET(CodeName) THEN Objstr.CodeName=CodeName
IF KEYWORD_SET(CodePi  ) THEN Objstr.CodePI  =CodePi
IF KEYWORD_SET(RunID   ) THEN Objstr.RunID   =RunId
IF KEYWORD_SET(FileId  ) THEN Objstr.FileID  =FIleID
;
; Set up the grid structure for the
; independent variable, time ('t').
;
tGrid = {Grid}
tGrid.mnemonic = 't'
tGrid.title = 't'
tGrid.units = '1/!4X!X!Dci!N'
tGrid.values = PTR_NEW(timex)
tGrid.boundary = 'Open'
tGrid.uniform = 1B
tmin = MIN(timex, MAX=tmax)
tGrid.range = [tmin, tmax]
tGrid.irange = [0, ind[1]-1]
;
; Start with total energy object
;
teStr = ObjStr
teStr.Mnemonic = 'TE'
teStr.Title = '!13E!X'
teStr.Indices = PTR_NEW(REPLICATE('*',1))
teStr.units = 'nT'				; ******GET HELP WITH THE UNITS *******
teStr.values = PTR_NEW(ke[0,*])
teMin = MIN(te[0, *], MAX=teMAX)
teStr.vrange = [teMin, teMax]
teStr.Grid1 = GKVsd_GridCopy(tGrid)
teObj = OBJ_NEW('GKVs1D', teStr)
;
histOutStr = {te:teObj}
;
; Make the kinetic energy object(s)
;
keOutArr = OBJARR(ind[0]+1)
FOR i=0, ind[0]-1 DO BEGIN  
	keStr = ObjStr
	iString = STRCOMPRESS( STRING(i+1, FORMAT='(i2)'), /REMOVE_ALL )
	keStr.mnemonic = 'KE_' + iString
	keStr.Title = '!12KE!X!D' + iString + '!N'
	keStr.Indices = PTR_NEW(REPLICATE('*',1))
	keStr.units = 'nT'				; ******GET HELP WITH THE UNITS *******
	keStr.values = PTR_NEW(ke[i,*])
	keMin = MIN(ke[i,*], MAX=keMAX)
	keStr.vrange = [keMin, keMax]
	keStr.Grid1 = GKVsd_GridCopy(tGrid)
	keObj = OBJ_NEW('GKVs1D', keStr)
	keOutArr[i+1] = keObj
	IF(i EQ 0) THEN BEGIN
		keSum = keObj -> MakeCopy()
	ENDIF ELSE BEGIN
		keSum = keSum -> plus(keObj)
	ENDELSE
ENDFOR
keSum -> Set, title='!12KE!X!Dtotal!N', mnemonic = 'KE_total'
keOutArr[0] = keSum
histOutStr = CREATE_STRUCT(histOutStr, 'KE', keOutArr)
;
; Make Field Energy object
;
feStr = ObjStr
feStr.Mnemonic = 'FE'
feStr.Title = '!12FE!X'
feStr.Indices = PTR_NEW(REPLICATE('*',1))
feStr.units = 'nT'				; ******GET HELP WITH THE UNITS *******
feStr.values = PTR_NEW(fe)
feMin = MIN(fe, MAX=feMAX)
feStr.vrange = [feMin, feMax]
feStr.Grid1 = GKVsd_GridCopy(tGrid)
feObj = OBJ_NEW('GKVs1D', feStr)
histOutStr = CREATE_STRUCT(histOutStr, 'FE', feObj)
;
; Make 'root-mean-square of potential' ojbect
;
rmsphiStr = ObjStr
rmsphiStr.Mnemonic = 'rmsPhi'
rmsphiStr.Title = '!4u!X!Drms!N'
rmsphiStr.Indices = PTR_NEW(REPLICATE('*',1))
rmsphiStr.units = 'T/e'				; ******GET HELP WITH THE UNITS *******
rmsphiStr.values = PTR_NEW(rmsphi)
rmsphiMin = MIN(rmsphi, MAX=rmsphiMAX)
rmsphiStr.vrange = [rmsphiMin, rmsphiMax]
rmsphiStr.Grid1 = GKVsd_GridCopy(tGrid)
rmsphiObj = OBJ_NEW('GKVs1D', rmsphiStr)
histOutStr = CREATE_STRUCT(histOutStr, 'rmsPhi', rmsphiObj)
;
; Make 'mean-squared-weight' ojbect
;
avgWsqStr = ObjStr
avgWsqStr.Mnemonic = 'Wsq'
avgWsqStr.Title = '!13<!XW!U2!N!13>!X'
avgWsqStr.Indices = PTR_NEW(REPLICATE('*',1))
avgWsqStr.units = ''				; ******GET HELP WITH THE UNITS *******
avgWsqStr.values = PTR_NEW(avgWsq)
avgWsqMin = MIN(avgWsq, MAX=avgWsqMAX)
avgWsqStr.vrange = [avgWsqMin, avgWsqMax]
avgWsqStr.Grid1 = GKVsd_GridCopy(tGrid)
avgWsqObj = OBJ_NEW('GKVs1D', avgWsqStr)
histOutStr = CREATE_STRUCT(histOutStr, 'avgWsq', avgWsqObj)
;
; Make particle flux object(s)
;
pflOutArr = OBJARR(ind[0]+1)
FOR i=0, ind[0]-1 DO BEGIN  
	pflStr = ObjStr
	iString = STRCOMPRESS( STRING(i+1, FORMAT='(i2)'), /REMOVE_ALL )
	pflStr.mnemonic = 'Pflux_' + iString
	pflStr.Title = '!4C!X!D' + iString + '!N'
	pflStr.Indices = PTR_NEW(REPLICATE('*',1))
	pflStr.units = 'n!De0!Nc!Ds!N'				; ******GET HELP WITH THE UNITS *******
	pflStr.values = PTR_NEW(pfl[i,*])
	pflMin = MIN(pfl[i,*], MAX=pflMAX)
	pflStr.vrange = [pflMin, pflMax]
	pflStr.Grid1 = GKVsd_GridCopy(tGrid)
	pflObj = OBJ_NEW('GKVs1D', pflStr)
	pflOutArr[i+1] = pflObj
	IF(i EQ 0) THEN BEGIN
		pflSum = pflObj -> MakeCopy()
	ENDIF ELSE BEGIN
		pflSum = pflSum -> plus(pflObj)
	ENDELSE
ENDFOR
pflsum -> Set, title='!4C!X!Dtotal!N', mnemonic = 'pflux_total'
pflOutArr[0] = pflSum
histOutStr = CREATE_STRUCT(histOutStr, 'pFlux', pflOutArr)
;
; Make energy flux object(s)
;
eflOutArr = OBJARR(ind[0]+1)
FOR i=0, ind[0]-1 DO BEGIN  
	eflStr = ObjStr
	iString = STRCOMPRESS( STRING(i+1, FORMAT='(i1)'), /REMOVE_ALL )
	eflStr.mnemonic = 'Eflux_' + iString
	eflStr.Title = 'Q!D' + iString + '!N'
	eflStr.Indices = PTR_NEW(REPLICATE('*',1))
	eflStr.units = 'n!De0!NTc!Ds!N'				; ******GET HELP WITH THE UNITS *******
	eflStr.values = PTR_NEW(efl[i,*])
	eflMin = MIN(efl[i,*], MAX=eflMAX)
	eflStr.vrange = [eflMin, eflMax]
	eflStr.Grid1 = GKVsd_GridCopy(tGrid)
	eflObj = OBJ_NEW('GKVs1D', eflStr)
	eflOutArr[i+1] = eflObj
	IF(i EQ 0) THEN BEGIN
		eflSum = eflObj -> MakeCopy()
	ENDIF ELSE BEGIN
		eflSum = eflSum -> plus(eflObj)
	ENDELSE
ENDFOR
eflSum -> Set, title='Q!Dtotal!N', mnemonic = 'Eflux_total'
eflOutArr[0] = eflSum
histOutStr = CREATE_STRUCT(histOutStr, 'eFlux', eflOutArr)
;
; Make 'phi mode' objects
;
modeOutArr = OBJARR(ind[2])
indices = STRARR(4)
indices[3] = '*'
FOR m = 0, ind[2]-1 DO BEGIN
	mString = STRCOMPRESS( STRING(m, FORMAT='(i2)'), /REMOVE_ALL)
	pModeHisStr = ObjStr
	pModeHisStr.mnemonic = 'Phi_' + mString
	pModeHisStr.Title = '!4u!X'
	FOR i=0,2 DO indices[i] = STRCOMPRESS( STRING(modes[i,m], FORMAT='(i2)'), /REMOVE_ALL ) 
	pModeHisStr.Indices = PTR_NEW(indices)
	pModeHisStr.units = 'T/e'				; ******GET HELP WITH THE UNITS *******
	pModeHisStr.values = PTR_NEW(pModeHis[m,*])
	pModeHisMin = GKVsd_MIN(pModeHis[m,*], MAX=pModeHisMAX)
	pModeHisStr.vrange = [pModeHisMin, pModeHisMax]
	pModeHisStr.Grid1 = GKVsd_GridCopy(tGrid)
	pModeHisObj = OBJ_NEW('GKVs1D', pModeHisStr)
	modeOutArr[m] =  pModeHisObj
ENDFOR
histOutStr = CREATE_STRUCT(histOutStr, 'Mode', modeOutArr)	

output = CREATE_STRUCT('hist', histOutStr)

ReadXYout: 
;
; Read 'phi' 
;
	phi = {firstTag:0}

	PhixyObj = GKV_ReadTUBExy('phi', nm, Title = '!4u!X', Units='T/e', 	$
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIn, rhoStar=rhoStar	)
	IF(TypeOf(PhixyObj) EQ 11) THEN BEGIN
		PhixyObj -> Set, CodeName=ObjStr.Codename
		PhixyObj -> Set, CodePI  =ObjStr.CodePI
		PhixyObj -> Set, RunID   =ObjStr.RUNID
		PhixyObj -> Set, FileID  =ObjStr.FileID
		PhixyObj -> ScaleAxis, 't', const=dt
		phi = CREATE_STRUCT(phi, 'xy', PhixyObj)
	ENDIF
	
	PhixzObj = GKV_ReadTUBExz('phi', nm, Title = '!4u!X', Units='T/e', 	$
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIN, rhoStar=rhoStar	)
	IF(TypeOf(PhixzObj) EQ 11) THEN BEGIN
		PhixzObj -> Set, CodeName=ObjStr.Codename
		PhixzObj -> Set, CodePI  =ObjStr.CodePI
		PhixzObj -> Set, RunID   =ObjStr.RUNID
		PhixzObj -> Set, FileID  =ObjStr.FileID
		PhixzObj -> ScaleAxis, 't', const=dt
		phi = CREATE_STRUCT(phi, 'xz', PhixzObj)
	ENDIF
	
	PhirObj  = GKV_ReadTUBEr( 'phi', nm, Title = '!4u!X', Units='T/e', 	$
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIN, rhoStar=rhoStar	)
	IF(TypeOf(PhirObj)  EQ 11) THEN BEGIN
		PhirObj -> Set, CodeName=ObjStr.Codename
		PhirObj -> Set, CodePI  =ObjStr.CodePI
		PhirObj -> Set, RunID   =ObjStr.RUNID
		PhirObj -> Set, FileID  =ObjStr.FileID
		PhirObj -> ScaleAxis, 't', const=dt
		phi = CREATE_STRUCT(phi, 'r' , PhirObj )
	ENDIF
	
	output = CREATE_STRUCT(output, 'phi', phi)
;
; Read 'rho'
;
	rho = {firstTag:0}
	
	RhoxyObj = GKV_ReadTUBExy('rho', nm, Title = '!4q!X', Units='en!De0!N', $
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIN, rhoStar=rhoStar	) 
	IF(TypeOf(RhoxyObj) EQ 11) THEN BEGIN
		RhoxyObj -> Set, CodeName=ObjStr.Codename
		RhoxyObj -> Set, CodePI  =ObjStr.CodePI
		RhoxyObj -> Set, RunID   =ObjStr.RUNID
		RhoxyObj -> Set, FileID  =ObjStr.FileID
		RhoxyObj -> ScaleAxis, 't', const=dt
		rho = CREATE_STRUCT(rho, 'xy', RhoxyObj)
	ENDIF
	
	RhoxzObj = GKV_ReadTUBExz('rho', nm, Title = '!4q!X', Units='en!De0!N', $
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIN, rhoStar=rhoStar	)
	IF(TypeOf(RhoxzObj) EQ 11) THEN BEGIN
		RhoxzObj -> Set, CodeName=ObjStr.Codename
		RhoxzObj -> Set, CodePI  =ObjStr.CodePI
		RhoxzObj -> Set, RunID   =ObjStr.RUNID
		RhoxzObj -> Set, FileID  =ObjStr.FileID
		RhoxzObj -> ScaleAxis, 't', const=dt
		rho = CREATE_STRUCT(rho, 'xz', RhoxzObj)
	ENDIF
	
	RhorObj  = GKV_ReadTUBEr( 'rho', nm, Title = '!4q!X', Units='en!De0!N', $
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIN, rhoStar=rhoStar	) 
	IF(TypeOf(RhorObj)  EQ 11) THEN BEGIN
		RhorObj -> Set, CodeName=ObjStr.Codename
		RhorObj -> Set, CodePI  =ObjStr.CodePI
		RhorObj -> Set, RunID   =ObjStr.RUNID
		RhorObj -> Set, FileID  =ObjStr.FileID
		RhorObj -> ScaleAxis, 't', const=dt
		rho = CREATE_STRUCT(rho, 'r' , RhorObj )
	ENDIF
	
	output = CREATE_STRUCT(output, 'rho', rho)
;
; Read 'dNe'
;
	dNe = {firstTag:0}
	
	dNexyObj = GKV_ReadTUBExy('dne', nm, Title = '!4d!Xn!de!N', Units='n!De0!N',	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		) 
	IF(TypeOf(dNexyObj) EQ 11) THEN BEGIN
		dNexyObj -> Set, CodeName=ObjStr.Codename
		dNexyObj -> Set, CodePI  =ObjStr.CodePI
		dNexyObj -> Set, RunID   =ObjStr.RUNID
		dNexyObj -> Set, FileID  =ObjStr.FileID
		dNexyObj -> ScaleAxis, 't', const=dt
		dNe = CREATE_STRUCT(dNe, 'xy', dNexyObj)
	ENDIF
	
	dNexzObj = GKV_ReadTUBExz('dne', nm, Title = '!4d!Xn!de!N', Units='n!De0!N', 	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(dNexzObj) EQ 11) THEN BEGIN
		dNexzObj -> Set, CodeName=ObjStr.Codename
		dNexzObj -> Set, CodePI  =ObjStr.CodePI
		dNexzObj -> Set, RunID   =ObjStr.RUNID
		dNexzObj -> Set, FileID  =ObjStr.FileID
		dNexzObj -> ScaleAxis, 't', const=dt
		dNe = CREATE_STRUCT(dNe, 'xz', dNexzObj)
	ENDIF
	
	dNerObj  = GKV_ReadTUBEr( 'dne', nm, Title = '!4d!Xn!de!N', Units='n!De0!N', 	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(dNerObj)  EQ 11) THEN BEGIN
		dNerObj -> Set, CodeName=ObjStr.Codename
		dNerObj -> Set, CodePI  =ObjStr.CodePI
		dNerObj -> Set, RunID   =ObjStr.RUNID
		dNerObj -> Set, FileID  =ObjStr.FileID
		dNerObj -> ScaleAxis, 't', const=dt
		dNe = CREATE_STRUCT(dNe, 'r' , dNerObj )
	ENDIF
	
	output = CREATE_STRUCT(output, 'dNe', dNe)
;
; Read 'dNi'
;
	dNi = {firstTag:0}
	
	dNixyObj = GKV_ReadTUBExy('dni', nm, Title = '!4d!Xn!di!N', Units='n!De0!N', 	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(dNixyObj) EQ 11) THEN BEGIN
		dNixyObj -> Set, CodeName=ObjStr.Codename
		dNixyObj -> Set, CodePI  =ObjStr.CodePI
		dNixyObj -> Set, RunID   =ObjStr.RUNID
		dNixyObj -> Set, FileID  =ObjStr.FileID
		dNixyObj -> ScaleAxis, 't', const=dt
		dNi = CREATE_STRUCT(dNi, 'xy', dNixyObj)
	ENDIF
	
	dNixzObj = GKV_ReadTUBExz('dni', nm, Title = '!4d!Xn!di!N', Units='n!De0!N', 	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(dNixzObj) EQ 11) THEN BEGIN
		dNixzObj -> Set, CodeName=ObjStr.Codename
		dNixzObj -> Set, CodePI  =ObjStr.CodePI
		dNixzObj -> Set, RunID   =ObjStr.RUNID
		dNixzObj -> Set, FileID  =ObjStr.FileID
		dNixzObj -> ScaleAxis, 't', const=dt
		dNi = CREATE_STRUCT(dNi, 'xz', dNixzObj)
	ENDIF
	
	dNirObj  = GKV_ReadTUBEr( 'dni', nm, Title = '!4d!Xn!di!N', Units='n!De0!N', 	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(dNirObj)  EQ 11) THEN BEGIN
		dNirObj -> Set, CodeName=ObjStr.Codename
		dNirObj -> Set, CodePI  =ObjStr.CodePI
		dNirObj -> Set, RunID   =ObjStr.RUNID
		dNirObj -> Set, FileID  =ObjStr.FileID
		dNirObj -> ScaleAxis, 't', const=dt
		dNi = CREATE_STRUCT(dNi, 'r' , dNirObj )
	ENDIF
	
	output = CREATE_STRUCT(output, 'dNi', dNi)
;
; Read 'apa'
;
	apa = {firstTag:0}
	
	apaxyObj = GKV_ReadTUBExy('apa', nm, Title = 'A!D!9#!X!N', units="(1/!4q!X!Ds!N)(T/e)",	$
					Yang=yang, iSkip=iSkip_in,				$
					aScale=aScaleIN, rhoStar=rhoStar			)
	IF(TypeOf(apaxyObj) EQ 11) THEN BEGIN
		apaxyObj -> Set, CodeName=ObjStr.Codename
		apaxyObj -> Set, CodePI  =ObjStr.CodePI
		apaxyObj -> Set, RunID   =ObjStr.RUNID
		apaxyObj -> Set, FileID  =ObjStr.FileID
		apaxyObj -> ScaleAxis, 't', const=dt
		apa = CREATE_STRUCT(apa, 'xy', apaxyObj)
	ENDIF
	
	apaxzObj = GKV_ReadTUBExz('apa', nm, Title = 'A!D!9#!X!N', units="(1/!4q!X!Ds!N)(T/e)",	$
					Yang=yang, iSkip=iSkip_in,				$
					aScale=aScaleIN, rhoStar=rhoStar			)
	IF(TypeOf(apaxzObj) EQ 11) THEN BEGIN 
		apaxzObj -> Set, CodeName=ObjStr.Codename
		apaxzObj -> Set, CodePI  =ObjStr.CodePI
		apaxzObj -> Set, RunID   =ObjStr.RUNID
		apaxzObj -> Set, FileID  =ObjStr.FileID
		apaxzObj -> ScaleAxis, 't', const=dt
		apa = CREATE_STRUCT(apa, 'xz', apaxzObj)
	ENDIF
	
	aparObj  = GKV_ReadTUBEr( 'apa', nm, Title = 'A!D!9#!X!N', units="(1/!4q!X!Ds!N)(T/e)",	$
					Yang=yang, iSkip=iSkip_in,				$
					aScale=aScaleIN, rhoStar=rhoStar			)
	IF(TypeOf(aparObj)  EQ 11) THEN BEGIN
		aparObj -> Set, CodeName=ObjStr.Codename
		aparObj -> Set, CodePI  =ObjStr.CodePI
		aparObj -> Set, RunID   =ObjStr.RUNID
		aparObj -> Set, FileID  =ObjStr.FileID
		aparObj -> ScaleAxis, 't', const=dt
		apa = CREATE_STRUCT(apa, 'r',  aparObj )
	ENDIF
	
	output = CREATE_STRUCT(output, 'apa', apa)
;
; Read 'upa'
;
	upa = {firstTag:0}
	
	upaxyObj = GKV_ReadTUBExy('uap', nm, Title = 'u!D!9#!X!N', units="c!Ds!N",	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		) 
	IF(TypeOf(upaxyObj) EQ 11) THEN BEGIN
		upaxyObj -> Set, CodeName=ObjStr.Codename
		upaxyObj -> Set, CodePI  =ObjStr.CodePI
		upaxyObj -> Set, RunID   =ObjStr.RUNID
		upaxyObj -> Set, FileID  =ObjStr.FileID
		upaxyObj -> ScaleAxis, 't', const=dt
		upa = CREATE_STRUCT(upa, 'xy', upaxyObj)
	ENDIF

	upaxzObj = GKV_ReadTUBExz('uap', nm, Title = 'u!D!9#!X!N', units="c!Ds!N",	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(upaxzObj) EQ 11) THEN BEGIN
		upaxzObj -> Set, CodeName=ObjStr.Codename
		upaxzObj -> Set, CodePI  =ObjStr.CodePI
		upaxzObj -> Set, RunID   =ObjStr.RUNID
		upaxzObj -> Set, FileID  =ObjStr.FileID
		upaxzObj -> ScaleAxis, 't', const=dt
		upa = CREATE_STRUCT(upa, 'xz', upaxzObj)
	ENDIF

	uparObj  = GKV_ReadTUBEr( 'uap', nm, Title = 'u!D!9#!X!N', units="c!Ds!N",	$
					Yang=yang, iSkip=iSkip_in,			$
					aScale=aScaleIN, rhoStar=rhoStar		)
	IF(TypeOf(uparObj)  EQ 11) THEN BEGIN
		uparObj -> Set, CodeName=ObjStr.Codename
		uparObj -> Set, CodePI  =ObjStr.CodePI
		uparObj -> Set, RunID   =ObjStr.RUNID
		uparObj -> Set, FileID  =ObjStr.FileID
		uparObj -> ScaleAxis, 't', const=dt
		upa = CREATE_STRUCT(upa, 'r' , uparObj)
	ENDIF
	
	output = CREATE_STRUCT(output, 'upa', upa)

RETURN, output
END	; ****** Tube_Data ****** ;
