FUNCTION GTC_ShearEB, FileName=fileName, Path=path, NT=ntime, DT=delT, Debug=d, Verbose=verbose
;
; Reads the 'sheareb.out' file from GTC.
; Returns a structure with tags corresponding to the
; items stored in this file.
;
; Keywords:
;
;	FileName	An optional 'string' argument containing the FULL name of  
;			the 'sheareb.out' file. If this keyword is not set (or the 
;			requested file cannot be found) then a dialogue box will
;			pop up to allow the user to select the 'history.out' file.
;			When in doubt, ***DON'T*** set this keyword.
;
;	Path		An optional 'string' argument containing the path to the 
;			directory (folder on MACs) where the pop-up dialogue box
;			will initialize the user's search for the the 'sheareb.out'
;			file.  When in doubt, ***DO*** set this keyword.
;
;	NT		An optional argument containing (on input) the number of
;			time slices to be read from 'sheareb.out'.  If set, NT
;			***MUST*** be not be more than the total number of time
;			slices in 'sheareb.out'.  If NT is not set, then GTC_SearEB
;			reads all time slices in 'sheareb.out' using READ_ASCII. It's 
;			 best to set NT if you know how many time slices are present, 
;			as this speeds the reading process.
;
;	DT		An optional argument containing (on input) the interval
;			between time slices (in units of 1/Omega_ci).  Defaults
;			to 20.
;
;	Debug		Set this keyword (i.e, "/Debug") to turn off error trapping.
;			disabling error trapping is most useful during debugging.
;
;	Verbose		Set this keyword (i.e., "/Verbose") so that READ_ASCI prints
;			out a message upon reading ***EACH*** record in 'sheareb.out'
;			This ***REALLY*** slows the read to a crawl, so ***DON'T*** 
;			set 'Verbose' unless you know that 'sheareb.out' has less 
;			than about 1000 records (or so...), or you are really motivated
;			to know how many records are in 'sheareb.out'.
;
; Written by W.M. Nevins
;	7/6/00
;
deBug=0
IF(N_ELEMENTS(d) NE 0) THEN deBug=d
IF(deBug EQ 0) THEN BEGIN						; Set Keyword 'DeBug' to avoid error trap
	Catch, error 							; Set up error trap
	IF error NE 0 then Begin 
		Catch, /Cancel             			; Cancel error trap 
		ok=Error_Message(/Trace)  				; Print out a trace-back in the output log.
		IF(shearID GE 0) THEN FREE_LUN, shearID		; Close any open NCDF files
		IF(N_ELEMENTS(gkvObjs) NE 0) THEN $
   			RETURN, gkvObjs					; Return SOMETHING useful if possible;
   		RETURN, 0                  			;	return 0 if all else fails. 
	ENDIF
ENDIF
;
; find 'sheareb.out' file if needed
;
IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out') ELSE fineIn=fileName
ok = FINDFILE(fileIn, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, 'Bad filename', /Informational
	fileIn=DIALOG_PICKFILE(Path=path, /MustExist, Filter='*.out')
ENDIF	
;
; Open 'sheareb.out' file
;
bufSize = 2L^20
GET_LUN, shearID							; Get a free Logical Unit Number
OPENR, shearID, fileIn, BUFSIZE=bufsize
;
; Read header data from sheareb file
;	
READF, shearID, mpsi
;
; Get number of time steps, interval between time samples
;
nT = 0
IF KEYWORD_SET(ntime) THEN nT=ntime
dT = 20.
IF KEYWORD_SET(delT)  THEN dT=delT

IF(nt EQ 0) THEN BEGIN
	;
	; Close sheareb file
	;
	FREE_LUN, shearID
	shearID=-1
	;
	; Use Read_Ascii to get ***ALL*** data from 'sheareb.out' file
	;
	dataStr = READ_ASCII(fileIN, Data_Start=1, Count=count, Verbose=verbose)
	nT = count/(2L*mpsi)
	sData = REFORM(dataStr.(0), mpsi, 2, nT)
ENDIF ELSE BEGIN
	;
	; Create array to hold data in sheareb file
	;
	sData = FLTARR(mpsi, 2, nT)
	;
	; Read data from sheareb file into 'sData'
	;
	READF, shearID, sData
	;
	; Close sheareb file
	;
	FREE_LUN, shearID
	shearID=-1
ENDELSE
;
; Compute time base
;
t=dT*FINDGEN(nT)
timeUnits = '1/!4X!X!Dci!N'
;
; Set up array of mnemonics and "pretty" names
;
mnemonics = STRARR(2)
mnemonics = [	'Q_E',				'V_ExB'		]

titles = STRARR(2)
titles = [		"Q!DE!N",				"V!DExB!N"		]

units = STRARR(2)
units = [		"n!Di!NT!Di!Nv!Dti!N",	"v!Dti!N"		]
	
nobjects = 2
gkvObjs = OBJARR(nobjects)
str2D = {GKVs2D}
str2D.CodeName = 'GTC'
str2D.CodePI = 'Z. Lin'
Indices = REPLICATE('*', 2)
str2D.Indices = PTR_NEW(Indices)
;
; Add any additional runinfo here ...
;

;
; Set up time 'Grid'
;
tGrid = {Grid}
tGrid.Mnemonic = 't'
tGrid.title = 't'
tGrid.values = PTR_NEW(t)
tGrid.units = timeUnits
tGrid.boundary = 'open'
tGrid.uniform = 1b
tmin = MIN(t, max=tmax)
tGrid.range = [tmin, tmax]
tGrid.irange = [0, nT-1]
;
; Set up radial 'Grid'
;
rGrid  = {Grid}
rGrid.Mnemonic = 'r'
rGrid.title = 'r'
rValues = 0.1 + (0.8/(mpsi-1))*FINDGEN(mpsi)
rGrid.values = PTR_NEW(rValues)
rGrid.units = 'a'
rmin = MIN(rValues, MAX=rmax)
rGrid.range  = [rmin, rmax]
rGrid.irange = [0, (mpsi-1)]
;
; Create GKVs2D objects for sheareb data
;
FOR i=0,1 DO BEGIN
	thisObj = str2D
	thisObj.mnemonic = mnemonics[i]
	thisObj.Title = titles[i]
	values = REFORM(sData[*,i,*])
	thisObj.values = PTR_NEW(values)
	thisObj.units = units[i]
	thisObj.Grid1 = GKVsd_GridCopy(rGrid)
	thisObj.Grid2 = GKVsd_GridCopy(tGrid)
	gkvObjs[i] = OBJ_NEW('GKVs2D', thisObj)
ENDFOR
;
; Create structure 'result' 
;
result = { Q_E:gkvObjs[0], V_ExB:gkvObjs[1] }
;
; Turn off error control
;
IF(debug EQ 0) THEN CATCH, /CANCEL
;
; ... and we're done!
;
RETURN, result
END ; ****** GTC_ShearEB ****** ;




