FUNCTION GKV_ReadPEFluxes, 	nSpecies, nZones, nSteps, 	$
				FluxFile = FluxFile, 		$
				rRange = rRangeIn,		$
				RunID = RunID, FileID = FileID			
;
; Purpose:
;
;	Gets energy and particle flux data from 
;       pefluxes.out file (if it exists)
;	and returns this data as arrays of GKV objects
;
; Arguments:
;
;	nSpecies	Number of particle species
;	
;	nZones		Number of radial zones
;
;	nSteps		Number of time slices
;
; KeyWords:
;
;	FluxFile	Set to name of GEM .out file 
;			containing flux data. 
;			Defaults to "pefluxes.out"
;			(Optional).
;
;	RunID		Set to (ascii) run identifier,
;			which will appear as upper field
;			in the lower left-hand corner of
;			plots created from output flux
;			objects. (Optional) 
;
;	FileID		Set to File Identifier,
;			which will appear as upper field
;			in the lower left-hand corner of
;			plots created from output flux
;			objects. (Optional)
;
;	rRange		Set to a two-element floating-point 
;			array containing the minimum and
;			maximum radii of the simulation 
;			volume. Default is to obtain this
;			info from thte PEFluxses.out file.
;			(Optional)
;			 
;
;
;  Written by W.M. Nevins
;	7/21/08
;
FORWARD_FUNCTION GKVsd_MIN
nArgs = N_PARAMS()
IF(nArgs NE 3) THEN BEGIN
	MESSAGE, "Routine called with wrong number of arguments", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Check for pefluxes.out file
;
OutFile = "pefluxes.out"
IF(TypeOF(FluxFile) EQ 7) THEN OutFIle = FluxFIle
result = FINDFILE(OutFile, COUNT=ok)
IF(NOT ok) THEN BEGIN
	MESSAGE, "Couldn't find pefluxes.out", /INFORMATIONAL
	RETURN, 0
ENDIF
MESSAGE, "Reading " + OutFile, /INFORMATIONAL

OPENR, lun, OutFile, /GET_LUN
;
; Read header lines.
; (to get past them ...)
;
HeaderLines = STRARR(4)
READF, lun, HeaderLines
;
; Now read associated values
; (For the mement, we don't actually use this info ...)
;
READF, lun, nnSpecies, nnSteps, eUnit, refDensity, LengthUnit,	$
	FORMAT='(I5, I10, 5X, 3E20.6)'
READF, lun, frequencyUnit, velocityUnit, rMin, rMax, FORMAT='(4E20.6)'
;
; Set radail range
;
rRange=[rmin, rmax]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange = rRangeIn
;
; Now loop over particle species
;
SpeciesHeader = STRARR(2)
FluxHeader    = STRARR(3)
pFluxData = FLTARR(nZones+1, nSteps+1, nSpecies)
eFluxData = FLTARR(nZones+1, nSteps+1, nSpecies)
thisArray = FLTARR(nZones+1, nSteps)
thisFormat= '(I5, ' + STRTRIM(STRING(nZones),2) + 'E13.5)'
FOR iSpecies = 0, nSpecies-1 DO BEGIN
	READF, lun, SpeciesHeader     	; Get past species header
	READF, lun, FluxHeader		; Get past particle flux header
	READF, lun, thisArray, FORMAT=thisFormat
	pFluxData[*, 1:nSteps, iSpecies] = thisArray
	READF, lun, FluxHeader		; Get past energy flux header
	READF, lun, thisArray, FORMAT=thisFormat
	eFluxData[*, 1:nSteps, iSpecies] = thisArray
ENDFOR
;
; Create Arrays to hold output objects
;
pFlux = OBJARR(nSpecies)
eFlux = OBJARR(nSpecies)
;
; Create master object
;
objStr = {GKVs2D}
;
; Store mnemonics it string arrays
;
pMnemonic = STRARR(nSpecies)
eMnemonic = STRARR(nSpecies)
;
pMnemonic[0] = "Gamma_e"
eMnemonic[0] = "Q_e"
FOR iSpecies = 1, nSpecies-1 DO BEGIN
	SpeciesString = STRTRIM(STRING(iSpecies), 2)
	pMnemonic[iSpecies] = "Gamma_ion" + SpeciesString
	eMnemonic[iSpecies] = "Q_ion"     + SpeciesString
ENDFOR
;
; Store titles in string arrays
;
pTitle = STRARR(nSpecies)
eTitle = STRARR(nSpecies)
pTitle[0] = "!4C!X!De!N"
eTitle[0] = "Q!De!N"
FOR iSpecies = 1, nSpecies-1 DO BEGIN
	SpeciesString = STRTRIM(STRING(iSpecies), 2) + "!N"
	pTitle[iSpecies] = "!4C!X!Dion=" + SpeciesString
	eTitle[iSpecies] = "Q!Dion="     + SpeciesString
ENDFOR
objStr.Indices = PTR_NEW(["*", "*"])
; pUnits = ?
; eUnits = ?
objStr.CodeName = "GEM"
objStr.CodePI   = "S. Parker & Y. Chen"
IF(TypeOF(RunID)  EQ 7) THEN objStr.RunID  = RunID
IF(TypeOF(FileID) EQ 7) THEN objstr.FileID = FileID
;
; Create radial grid structure
;
Grid1 = {Grid}
Grid1.mnemonic = 'r'
Grid1.title    = 'r'
Grid1.units    = 'r/a'
dr = FLOAT(rRange[1] - rRange[0])/(nZones-1)
rValues = rRange[0] + dr*FINDGEN(nZones)
Grid1.Values   = PTR_NEW(rValues)
Grid1.boundary = "open"
Grid1.uniform  = 1b
Grid1.range    = rRange
Grid1.irange   = [0, (nZones-1)]
;
; Create time grid structure
;
Grid2 = {Grid}
Grid2.mnemonic = 't'
Grid2.title    = 't'
Grid2.units    = "timeSteps"
tValues = FINDGEN(nSteps + 1)
Grid2.values   = PTR_NEW(tValues)
Grid2.range    = FLOAT([0, nSteps])
Grid2.irange   = [0, nSteps]
;
objStr.Grid1 = Grid1
objStr.Grid2 = Grid2
;
masterObj = OBJ_NEW("GKVs2D", ObjStr)
;
; Create particle and energy flux objects
;
FOR iSpecies=0, nSpecies-1 DO BEGIN
	;
	; Create particle flux objects
	;
	pFlux[iSpecies] = masterObj -> MakeCopy()
	pFlux[iSpecies] -> Set, Mnemonic=pMnemonic[iSpecies]
	pFlux[iSpecies] -> Set, Title=pTitle[iSpecies]
	pFlux[iSPecies] -> Set, Values=PTR_NEW(pFluxData[1:nZones,*,iSpecies])
	;
	; and energy flux objects
	;
	eFlux[iSpecies] = masterObj -> MakeCopy()
	eFlux[iSpecies] -> Set, Mnemonic=eMnemonic[iSpecies]
	eFlux[iSpecies] -> Set, Title=eTitle[iSpecies]
	eFlux[iSPecies] -> Set, Values=PTR_NEW(eFluxData[1:nzones,*,iSpecies])
ENDFOR
;
; Create output structure
;
result = {	Name	:	"GEM_FluxData",	$
		Particle:	pFlux,		$
		Energy	:	eFlux		}
		 


FREE_LUN, lun
RETURN, result
END	; ****** GKV_ReadPEFluxes ****** ;



FUNCTION GKV_ReadTUBEr, grdin, nm, Title=title, Units=units, 	$
			Yang=yang, iSkip=iSkip_in, 		$
			aScale=aScaleIn, rhoStar=rhoStar, 	$
			rRange=rRange, rTitle=rTitle, 		$
			rMnemonic=rMnemonic, rUnits=rUnits
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
MESSAGE, "Reading "+ rOutFile, /INFORMATIONAL

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

tLength = (iEnd-iStart)/iSkip + 1
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
k=k-1
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip
;
; Make a GKV object out of this data
;
ObjStr = {GKVs2D}				; Get a GKVs2D object structure
ObjStr.CodeName = 'GEM'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker & Y. Chen'		;	data.
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
IF(TypeOF(rMnemonic) EQ 7) THEN xGrid.mnemonic = rMnemonic
xGrid.title = 'x'
IF(TypeOF(rTitle) EQ 7) THEN xGrid.title = rTitle
xGrid.units = '!4q!X!Ds!N'
IF(TypeOF(rUnits) EQ 7) THEN xGrid.Units = rUnits
;
IF(N_ELEMENTS(rRange) EQ 2) THEN BEGIN
	Dx = im
	Dr = rRange[1] - rRange[0]
	xValues = rRange[0] + Dr*xValues/im
ENDIF
;
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
	result -> get, units=units
	result -> trash
	result=temp
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBEr ****** ;


FUNCTION GKV_ReadTUBExz, grdin, nm, Title=title, Units=units, 	$
			 Yang=yang, iSkip=iSkip_in, 		$
			 aScale=aScaleIn, rhoStar=rhoStar,	$
			rRange=rRange, rTitle=rTitle, 		$
			rMnemonic=rMnemonic, rUnits=rUnits
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
MESSAGE, "Reading " + xzOutFile, /INFORMATIONAL

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

tLength = (iEnd-iStart)/iSkip + 1
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
k=k-1
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip

;
; Make a GKV object out of this data
;
ObjStr = {GKVs3D}				; Get a GKVs3D object structure
ObjStr.CodeName = 'GEM'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker & Y. Chen'		;	data.
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
xValues = FINDGEN(im+1)
xGrid = {Grid}
xGrid.mnemonic = 'x'
IF(TypeOF(rMnemonic) EQ 7) THEN xGrid.mnemonic = rMnemonic
xGrid.title = 'x'
IF(TypeOF(rTitle) EQ 7) THEN xGrid.title = rTitle
xGrid.units = '!4q!X!Ds!N'
IF(TypeOF(rUnits) EQ 7) THEN xGrid.Units = rUnits
;
IF(N_ELEMENTS(rRange) EQ 2) THEN BEGIN
	Dx = im
	Dr = rRange[1] - rRange[0]
	xValues = rRange[0] + Dr*xValues/im
ENDIF
;
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
	result -> get, units=units
	result -> trash
	result=temp
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBExz ****** ;


FUNCTION GKV_ReadTUBExy, grdin, nm, Title=title, Units=units, 	$
			 Yang=yang, iSkip=iSkip_in, 		$
			 aScale=aScaleIn, rhoStar=rhoStar,	$
			rRange=rRange, rTitle=rTitle, 		$
			rMnemonic=rMnemonic, rUnits=rUnits

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
MESSAGE, "Reading " + xyOutFile, /INFORMATIONAL

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

tLength = (iEnd-iStart)/iSkip + 1
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
k=k-1
MakeObj:	FREE_LUN, lun
tLast = (k-iStart)/iSkip 
;
; Make a GKV object out of this data
;
ObjStr = {GKVs3D}			; Get a GKVs3D object structure
ObjStr.CodeName = 'GEM'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker & Y. Chen'	;	data.
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
xValues = FINDGEN(im+1)
xGrid = {Grid}
xGrid.mnemonic = 'x'
IF(TypeOF(rMnemonic) EQ 7) THEN xGrid.mnemonic = rMnemonic
xGrid.title = 'x'
IF(TypeOF(rTitle) EQ 7) THEN xGrid.title = rTitle
xGrid.units = '!4q!X!Ds!N'
IF(TypeOF(rUnits) EQ 7) THEN xGrid.Units = rUnits
;
IF(N_ELEMENTS(rRange) EQ 2) THEN BEGIN
	Dx = im
	Dr = rRange[1] - rRange[0]
	xValues = rRange[0] + Dr*xValues/im
ENDIF
;
xGrid.values = PTR_NEW(xValues)
xGrid.boundary = 'Periodic'
IF((im MOD 2) EQ 0) THEN xGrid.boundary = 'Periodic (closed)'
xGrid.uniform = 1B
xMin = GKVsd_MIN(xvalues, MAX=xMax)
xGrid.range = [xMin, xMax]
xGrid.irange = [0, im]
;
; set up the grid structure for the 
; independent variable, y
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
	result -> get, units=units
	result -> trash
	result=temp
	result -> set, units="(!4q!X!Ds!N/" + aScale + ")" + units
ENDIF
RETURN, result
END	; ****** GKV_ReadTUBExy ****** ;




FUNCTION GEM_Data,	Path=startPath, TubePath=tubePath, nt=nt_in, Yang=yang,		$ 
			CodeName=Codename, CodePI=CodePi, RunId=RunID, 			$
			FileId=FIleId, Debug=debug, iSkip=iSkip_in, 			$
			aScale=aScaleIN, rhoStar=rhoStar,	$
			rRange=rRange, rTitle=rTitle, 		$
			rMnemonic=rMnemonic, rUnits=rUnits,	$
			nZones = nZones_in, 			$
			nSpecies = nSpecies_in, Electrons=electrons

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
;	nt		Number of time steps the history file.  Defaults to
;			IND(1) from the 'hist.out' file. (Optional)
;
;  	nSpecies	Number of species in the PEFluxes.out file. Defaults to
;			IND(0) +1 from the hist.out file.  [IND(0) is the number 
;			of ION species, add one for electrons]. (Optional)
;
;	nZones		Number of radial zones in the PEFluxes.out file.  Don't
;			have a logical default. Will use 8 since this is number of
;			zones in first PEFluxes.out file obtained (7/3/08).
;			(Optional)
;
;	Yang		Set this keyword (i.e., put "/Yang" on the command line)
;			so that the third integer in the first record of the 
;			phi**.out, etc. files is interpreted as the number of time
;			slices in the file (instead of the interval between
;			timeslices in these files relative to the history
;			files).
;
;	Codename	Set this keyword to the name of the code which generated
;			the output files. Defaults to "GEM".  (Optional)
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
;	rRange		Set to a two element array array containing the minimum and
;			and maximum values of the radial variable for use in 
;			setting the radial scale. Defaults to labeling radial
;			scale by gridpoint number.
;
;	FluxFile	Name of PEFLuxes.out file. Defaults to "PEFluxes.out".
;			(Optional)
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
;   to estimate length of data in *xy.out, 
;   *xz.out, etc. files.
;
; Revised by W.M. Nevins
;	Sept. 2008
;   added capability to read peflux.out files
;   added capability to scale hist.out data
;   to gyrokinetic units
;   Fixed units field for scaled field data
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
	ObjStr = {GKVs1D}			; Get a GKVs1D object structure
	ObjStr.CodeName = 'GEM'			; 	and load in the common
	ObjStr.CodePI = 'S.E. Parker & Y. Chen'	;	data.
	; ****** PATCH FOR YANG'S KINETIC e's OF 12/00 *********
		IF KEYWORD_SET(Yang) THEN ObjStr.CodePI = 'Y. Chen'
	; ****************************************************** 
	IF KEYWORD_SET(CodeName) THEN Objstr.CodeName=CodeName
	IF KEYWORD_SET(CodePi  ) THEN Objstr.CodePI  =CodePi
	IF KEYWORD_SET(RunID   ) THEN Objstr.RunID   =RunId
	IF KEYWORD_SET(FileId  ) THEN Objstr.FileID  =FIleID
	GOTO, ReadXYout
ENDIF
;
; check rhoStar, and compute gyrokinetic scaling factos
;
eScale = 1.
fScale = 1.
tScale = 1.
IF(N_ELEMENTS(rhoStar) EQ 1) THEN BEGIN
	rhoStar = rhoStar < 1./rhoStar
	eScale = 1./rhoStar^2
	fScale = 1./rhoStar
	tScale  = rhoStar
ENDIF

aScale = 'L!DT!N'
IF(TypeOf(aScaleIn) EQ 7) THEN aScale=aScaleIn 

ind=lonarr(3)	
; 
;	ind[0] = number of ion species
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
;
; patch for versions of GEM with kinetic electrons
;
IF KEYWORD_SET(electrons) THEN ind[0] = ind[0] + 1
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
ObjStr = {GKVs1D}			; Get a GKVs1D object structure
ObjStr.CodeName = 'GEM'			; 	and load in the common
ObjStr.CodePI = 'S.E. Parker & Y. Chen'	;	data.
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
IF(N_ELEMENTS(rhoStar) EQ 1) THEN tGrid.units = '!4q!X!Ds!N/' + aScale
values = timex*tScale
tGrid.values = PTR_NEW(values)
tGrid.boundary = 'Open'
tGrid.uniform = 1B
tmin = MIN(values, MAX=tmax)
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
IF(N_ELEMENTS(rhoStar) EQ 1) THEN teStr.units = '(!4q!X!Ds!N/' + aScale + ')!U2!NnT'
values = ke[0,*]*eScale
teStr.values = PTR_NEW(values)
teMin = MIN(values, MAX=teMAX)
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
	IF(N_ELEMENTS(rhoStar) EQ 1) THEN keStr.units = '(!4q!X!Ds!N/' + aScale + ')!U2!NnT'
	values = ke[i,*]*eScale
	keStr.values = PTR_NEW(values)
	keMin = MIN(values, MAX=keMAX)
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
IF(N_ELEMENTS(rhoStar) EQ 1) THEN feStr.units = '(!4q!X!Ds!N/' + aScale + ')!U2!NnT'
values = fe*eScale
feStr.values = PTR_NEW(values)
feMin = MIN(values, MAX=feMAX)
feStr.vrange = [feMin, feMax]
feStr.Grid1 = GKVsd_GridCopy(tGrid)
feObj = OBJ_NEW('GKVs1D', feStr)
histOutStr = CREATE_STRUCT(histOutStr, 'FE', feObj)
;
; Make 'root-mean-square of potential' object
;
rmsphiStr = ObjStr
rmsphiStr.Mnemonic = 'rmsPhi'
rmsphiStr.Title = '!4u!X!Drms!N'
rmsphiStr.Indices = PTR_NEW(REPLICATE('*',1))
rmsphiStr.units = 'T/e'				; ******GET HELP WITH THE UNITS *******
IF(N_ELEMENTS(rhoStar) EQ 1) THEN rmsphiStr.units = '(!4q!X!Ds!N/' + aScale + ')T/e'
values = rmsphi*fScale
rmsphiStr.values = PTR_NEW(values)
rmsphiMin = MIN(values, MAX=rmsphiMAX)
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
IF(N_ELEMENTS(rhoStar) EQ 1) THEN avgWsqStr.units = '(!4q!X!Ds!N/' + aScale + ')!U2!N'
values = avgWsq*eScale
avgWsqStr.values = PTR_NEW(values)
avgWsqMin = MIN(values, MAX=avgWsqMAX)
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
	IF(N_ELEMENTS(rhoStar) EQ 1) THEN pflStr.units = '(!4q!X!Ds!N/' + aScale + ')!u2!Nn!De0!Nc!Ds!N'
	values = pfl[i,*]*eScale
	pflStr.values = PTR_NEW(values)
	pflMin = MIN(values, MAX=pflMAX)
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
	IF(N_ELEMENTS(rhoStar) EQ 1) THEN eflStr.units = '(!4q!X!Ds!N/' + aScale + ')!U2!Nn!De0!NTc!Ds!N'
	values = efl[i,*]*eScale
	eflStr.values = PTR_NEW(values)
	eflMin = MIN(values, MAX=eflMAX)
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
; Read PEFluxes.out file, and add data to output structure
;
nSpecies = IND[0] + 1
IF(N_ELEMENTS(NSpecies_in) EQ 1) THEN nSpecies = nSpecies_In
IF(N_ELEMENTS(nt_in) EQ 1) THEN nm = nSteps_in
nZones = 8
IF(N_ELEMENTS(nZones_in) EQ 1) THEN nZones=nZones_in
FluxStr = GKV_ReadPEFluxes(nSpecies, nZones, nm, 		$
				FluxFile = FluxFile, 		$
				rRange = rRangeIn,		$
				RunID = RunID, FileID = FileID	)	
output = CREATE_STRUCT(output, 'Flux', FluxStr)

;
; Read 'phi' 
;
	phi = {firstTag:0}

	PhixyObj = GKV_ReadTUBExy('phi', nm, Title = '!4u!X', Units='T/e', 	$
					Yang=yang, iSkip=iSkip_in, 		$
					aScale=aScaleIn, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
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
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
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
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
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
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	) 
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
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
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
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	) 
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		) 
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,			$
					rRange=rRange, rTitle=rTitle, 				$
					rMnemonic=rMnemonic, rUnits=rUnits			)
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
					aScale=aScaleIN, rhoStar=rhoStar,			$
					rRange=rRange, rTitle=rTitle, 				$
					rMnemonic=rMnemonic, rUnits=rUnits			)
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
					aScale=aScaleIN, rhoStar=rhoStar,			$
					rRange=rRange, rTitle=rTitle, 				$
					rMnemonic=rMnemonic, rUnits=rUnits			)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		) 
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
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
					aScale=aScaleIN, rhoStar=rhoStar,		$
					rRange=rRange, rTitle=rTitle, 			$
					rMnemonic=rMnemonic, rUnits=rUnits		)
	IF(TypeOf(uparObj)  EQ 11) THEN BEGIN
		uparObj -> Set, CodeName=ObjStr.Codename
		uparObj -> Set, CodePI  =ObjStr.CodePI
		uparObj -> Set, RunID   =ObjStr.RUNID
		uparObj -> Set, FileID  =ObjStr.FileID
		uparObj -> ScaleAxis, 't', const=dt
		upa = CREATE_STRUCT(upa, 'r' , uparObj)
	ENDIF
	
	output = CREATE_STRUCT(output, 'upa', upa)

;
; Read 'mph'
;
;	mph = {firstTag:0}
;	
;	mphxyObj = GKV_ReadTUBExy('mph', nm, Title = 'm!4u!X', units="T/e",	$
;					Yang=yang, iSkip=iSkip_in,		$
;					aScale=aScaleIN, rhoStar=rhoStar,	$
;					rRange=rRange, rTitle=rTitle, 		$
;					rMnemonic=rMnemonic, rUnits=rUnits	) 
;	IF(TypeOf(mphxyObj) EQ 11) THEN BEGIN
;		mphxyObj -> Set, CodeName=ObjStr.Codename
;		mphxyObj -> Set, CodePI  =ObjStr.CodePI
;		mphxyObj -> Set, RunID   =ObjStr.RUNID
;		mphxyObj -> Set, FileID  =ObjStr.FileID
;		upaxyObj -> ScaleAxis, 't', const=dt
;		mph = CREATE_STRUCT(mph, 'xy', mphxyObj)
;	ENDIF
;
;	mphxzObj = GKV_ReadTUBExz('mph', nm, Title = 'm!4u!X', units="T/e",	$
;					Yang=yang, iSkip=iSkip_in,		$
;					aScale=aScaleIN, rhoStar=rhoStar,	$
;					rRange=rRange, rTitle=rTitle, 		$
;					rMnemonic=rMnemonic, rUnits=rUnits	)
;	IF(TypeOf(mphxzObj) EQ 11) THEN BEGIN
;		mphxzObj -> Set, CodeName=ObjStr.Codename
;		mphxzObj -> Set, CodePI  =ObjStr.CodePI
;		mphxzObj -> Set, RunID   =ObjStr.RUNID
;		mphxzObj -> Set, FileID  =ObjStr.FileID
;		mphxzObj -> ScaleAxis, 't', const=dt
;		mph = CREATE_STRUCT(mph, 'xz', mphxzObj)
;	ENDIF
;
;	mphrObj  = GKV_ReadTUBEr( 'mph', nm, Title = 'm!4u!X', units="T/e",	$
;					Yang=yang, iSkip=iSkip_in,		$
;					aScale=aScaleIN, rhoStar=rhoStar,	$
;					rRange=rRange, rTitle=rTitle, 		$
;					rMnemonic=rMnemonic, rUnits=rUnits	)
;	IF(TypeOf(mphrObj)  EQ 11) THEN BEGIN
;		mphrObj -> Set, CodeName=ObjStr.Codename
;		mphrObj -> Set, CodePI  =ObjStr.CodePI
;		mphrObj -> Set, RunID   =ObjStr.RUNID
;		mphrObj -> Set, FileID  =ObjStr.FileID
;		mphrObj -> ScaleAxis, 't', const=dt
;		mph = CREATE_STRUCT(mph, 'r' , mphrObj)
;	ENDIF
;	
;	output = CREATE_STRUCT(output, 'mph', mph)
;
;
; Read 'Tpi'
;
	
	Tpi = {firstTag:0}
	
	TpixyObj = GKV_ReadTUBExy('tpi', nm, Title = 'T!D!9#!Xi!N', units="T",	$
					Yang=yang, iSkip=iSkip_in,		$
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	) 
	IF(TypeOf(TpixyObj) EQ 11) THEN BEGIN
		TpixyObj -> Set, CodeName=ObjStr.Codename
		TpixyObj -> Set, CodePI  =ObjStr.CodePI
		TpixyObj -> Set, RunID   =ObjStr.RUNID
		TpixyObj -> Set, FileID  =ObjStr.FileID
		TpixyObj -> ScaleAxis, 't', const=dt
		Tpi = CREATE_STRUCT(Tpi, 'xy', TpixyObj)
	ENDIF

	TpixzObj = GKV_ReadTUBExz('Tpi', nm, Title = 'T!D!9#!Xi!N', units="T",	$
					Yang=yang, iSkip=iSkip_in,		$
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
	IF(TypeOf(TpixzObj) EQ 11) THEN BEGIN
		TpixzObj -> Set, CodeName=ObjStr.Codename
		TpixzObj -> Set, CodePI  =ObjStr.CodePI
		TpixzObj -> Set, RunID   =ObjStr.RUNID
		TpixzObj -> Set, FileID  =ObjStr.FileID
		TpixzObj -> ScaleAxis, 't', const=dt
		Tpi = CREATE_STRUCT(Tpi, 'xz', TpixzObj)
	ENDIF

	TpirObj  = GKV_ReadTUBEr( 'Tpi', nm, Title = 'T!D!9#!Xi!N', units="T",	$
					Yang=yang, iSkip=iSkip_in,		$
					aScale=aScaleIN, rhoStar=rhoStar,	$
					rRange=rRange, rTitle=rTitle, 		$
					rMnemonic=rMnemonic, rUnits=rUnits	)
	IF(TypeOf(TpirObj)  EQ 11) THEN BEGIN
		TpirObj -> Set, CodeName=ObjStr.Codename
		TpirObj -> Set, CodePI  =ObjStr.CodePI
		TpirObj -> Set, RunID   =ObjStr.RUNID
		TpirObj -> Set, FileID  =ObjStr.FileID
		TpirObj -> ScaleAxis, 't', const=dt
		Tpi = CREATE_STRUCT(Tpi, 'r' , TpirObj)
	ENDIF
	
	output = CREATE_STRUCT(output, 'Tpi', Tpi)

RETURN, output
END	; ****** GEM_Data ****** ;
