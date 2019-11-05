FUNCTION GTC_History, FileName=fileName, Path=path, Debug=d
;
; Reads the 'history.out' file from GTC.
; Returns a structure with tags corresponding to the
; items stored in the history file.
;
; Keywords:
;
;	FileName	An optional 'string' argument containing the FULL name of  
;			the 'history.out' file. If this keyword is not set (or the 
;			requested file cannot be found) then a dialogue box will
;			pop up to allow the user to select the 'history.out' file.
;			When in doubt, ***DON'T*** set this keyword.
;
;	Path		An optional 'string' argument containing the path to the 
;			directory (folder on MACs) where the pop-up dialogue box
;			will initialize the user's search for the the 'history.out'
;			file.  When in doubt, ***DO*** set this keyword.
;
;	Debug		Set this keyword (i.e, "/Debug") to turn off error trapping.
;			disabling error trapping is most useful during debugging.
;
; Written by W.M. Nevins
;	7/5/00
;
deBug=0
IF(N_ELEMENTS(d) NE 0) THEN deBug=d
IF(deBug EQ 0) THEN BEGIN						; Set Keyword 'DeBug' to avoid error trap
	Catch, error 							; Set up error trap
	IF error NE 0 then Begin 
		Catch, /Cancel             			; Cancel error trap 
		ok=Error_Message(/Trace)  				; Print out a trace-back in the output log.
		IF(histID GE 0) THEN FREE_LUN, histID		; Close any open NCDF files
		IF(N_ELEMENTS(gkvObj) NE 0) THEN $
   			RETURN, gkvObj					; Return SOMETHING useful if possible;
   		RETURN, 0                  			;	return 0 if all else fails. 
	ENDIF
ENDIF
;
; find 'history.out file if needed
;
IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out') ELSE fineIn=fileName
ok = FINDFILE(fileIn, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, 'Bad filename', /Informational
	fileIn=DIALOG_PICKFILE(Path=path, Filter='*.out')
ENDIF	
;
; Open 'history.out' file
;
bufSize = 2L^20
GET_LUN, histID							; Get a free file unit number
OPENR, histID, fileIn	, BUFSIZE=bufsize
;
; Read header data from history file
;	
READF, histID, nrun, nother, nradial, nmode, ntime, delt
;
; Create array to hold data in history file
;
f=FLTARR(4*nmode+nradial+nother, ntime)
;
; Read data from history file into 'f'
;
READF, histID, f
;
; Close history file
;
FREE_LUN, histID
histID=-1
;
; Compute time base
;
t=delt*FINDGEN(ntime)
timeUnits = '1/!4X!X!Dci!N'
;
; Set up array of mnemonics and "pretty" names
;
mnemonics = STRARR(40)
mnemonics = [	"Exit",	"???",	"orbit_poloidal",	"orbit_flux",	"weight",		$	; [0:4]
			"V_para",	"dE_j",	"dL_j", 			"avg_Q_E",		"potential",	$	; [5:9]
			"marker", 	"dn_i", 	"V_rms",			"phi_rms",		"dKE",		$	; [10:14]
			"Wsq",	"flow",	"df_flow",			"particle_flux","L_flux",		$	; [15:19]
			"chi_i",	"shear1",	"shear2",			"shear3",		"shear4",		$	; [20:24]
			"shear5",	"shear6",	"shear7",			"shear8",		"mode1",		$	; [25:29]
			"mode2",	"mode3",	"mode4",			"mode5",		"mode6",		$	; [30:34]
			"mode7",	"mode8",	"peak chi_i",		"nbegin",		"naverage"		]	; [35:39]

titles = STRARR(40)
titles = [		"Exit",		"???",		"orbit poloidal", 	"orbit flux",	"W!Dj!N",			$	; [0:4]
			"V!D!9#!X!N",	"!4d!XE!Dj!N",	"!4d!XL!Dj!N",			"<Q!DE!N>",	"!4u!X",		$	; [5:9]
			"marker", 		"!4d!Xn!Di!N", 	"<V!S!DExB!R!E2!N>!E1/2!N", $
													"<!4du!X!E2!N>!E1/2!N", $
																"!4d!XE!Di!N",		$	; [10:14]
			"!4R!XW!S!Di!R!E2!N", $
						"N!S!Dp!R!E-1!N!4R!X!Dj!Nv!D!9#!Xj!N", 							$ 	; "1/N!Dp!N!4S!Dj!Nv!9#!Xj!N", $
									"!4R!X!Dj!NW!Dj!Nv!D!9#!Xj!N", $
													"!4C!X!Di!N", $
																"!7P!X!Dr,!9#!X!N",	$	; [15:19]
			"!4v!X!Di!N",	"V!DExB!N",	"V!DExB!N",		"V!DExB!N",	"V!DExB!N",		$	; [20:24]
			"V!DExB!N",	"V!DExB!N",	"V!DExB!N",		"V!DExB!N",	"!4u!X(m,n) mode1",	$	; [25:29]
		"!4u!X(m,n) mode2", "!4u!X(m,n) mode3", "!4u!X(m,n) mode4", "!4u!X(m,n) mode5", "!4u!X(m,n) mode6",		$	; [30:34]
		"!4u!X(m,n) mode7", "!4u!X(m,n) mode8",	"peak chi_i",		"nbegin",		"naverage"			]	; [35:39]

units = [		"",			"",			"",				"",			"",				$	; [0:4]
			"v!Dti!N",		"E!Dj0!N",		"L!Dj0!N",		"n!Di!NT!Di!Nv!Dti!N","T!Di!N/e",		$	; [5:9]
			"",			"n!Di!N",		"v!Dti!N",			"T!Di!N/e",	"n!Di!NT!Di!N",		$	; [10:14]
			"",			"v!Dti!N",		"v!Dti!N",			"n!Di!Nv!Dti!N","n!Di!Nv!S!Dti!R!E2!N", $	; [15:19]
		"a!E2!N!4X!X!Dci!N",	"v!Dti!N",		"v!Dti!N",			"v!Dti!N",		"v!Dti!N",			$	; [20:24]
			"v!Dti!N",		"v!Dti!N",		"v!Dti!N",			"v!Dti!N",		"T!Di!N/e"			]	; [25:29]
	
commandStr = 'result={'
nobjects = 2*nmode + nother
gkvObjs = OBJARR(nobjects)
str1D = {GKVs1D}
str1D.CodeName = 'GTC'
str1D.CodePI = 'Z. Lin'
Indices = REPLICATE('*', 1)
str1D.Indices = PTR_NEW(Indices)
;
; Add any additional runinfo here 
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
tGrid.irange = [0, ntime-1]
iObj = 0
FOR i=4, nother-1 DO BEGIN
	thisObj = str1D
	thisObj.mnemonic = mnemonics[i]
	thisObj.Title = titles[i]
	values = f[i,*]
	thisObj.values = PTR_NEW(values)
	thisObj.units = units[i]
	thisObj.Grid1 = GKVsd_GridCopy(tGrid)
	gkvObjs[iobj] = OBJ_NEW('GKVs1D', thisObj)
	iObjStr = STRING(iObj, FORMAT='(I2)')
	iObjStr = STRTRIM(iObjStr, 2)
	commandStr = commandStr + mnemonics[i] + ':gkvObjs[' + iObjStr + '], '
	iObj = iObj+1
ENDFOR
;
; Form GKVs1D objects for V_ExB vs. k_r
;
kr0 = 2*!PI/0.8
krValues = STRARR(nmode+1)
krValues[0] = 'k!Dr!Na='
FOR i=1, nmode DO krValues[i]=STRTRIM(STRING(i*kr0,FORMAT='(G10.3)'),2)
indices = REPLICATE('*',2)
istart = nother+nradial-2
FOR i=1, nmode DO BEGIN
	thisObj = str1D
	thisObj.mnemonic = mnemonics[20+i]
	thisObj.Title = titles[20+i]
	rvalues = f[istart + 2*i, *]
	ivalues = f[istart + 2*i + 1, *]
	values = COMPLEX(rvalues, ivalues)
	thisObj.values = PTR_NEW(values)
	thisObj.units = units[20+i]
	indices[0] = krValues[0] + krValues[i]
	thisObj.indices = PTR_NEW(indices)
	thisObj.Grid1 = GKVsd_GridCopy(tGrid)
	gkvObjs[iobj] = OBJ_NEW('GKVs1D', thisObj)
	iObjStr = STRING(iObj, FORMAT='(I2)')
	iObjStr = STRTRIM(iObjStr, 2)
	commandStr = commandStr + mnemonics[20+i] + ':gkvObjs[' + iObjStr + '], '
	iObj = iObj+1
ENDFOR
;
; Form GKVs1D objects for (m,n) modes
;
istart = nother+nradial+2*nmode-2
FOR i=1, nmode DO BEGIN
	thisObj = str1D
	thisObj.mnemonic = mnemonics[28+i]
	thisObj.Title = titles[28+i]
	rvalues = f[istart + 2*i, *]
	ivalues = f[istart + 2*i + 1, *]
	values = COMPLEX(rvalues, ivalues)
	thisObj.values = PTR_NEW(values)
	thisObj.units = units[29]
	thisObj.Grid1 = GKVsd_GridCopy(tGrid)
	gkvObjs[iobj] = OBJ_NEW('GKVs1D', thisObj)
	iObjStr = STRING(iObj, FORMAT='(I2)')
	iObjStr = STRTRIM(iObjStr, 2)
	commandStr = commandStr + mnemonics[28+i] + ':gkvObjs[' + iObjStr + '], '
	iObj = iObj+1
ENDFOR
;
; Form GKVs2D object for (radially resolved) heat flux
;
thisObj = {GKVs2D}
rGrid = {Grid}
rGrid.Mnemonic = 'r'
rGrid.Title = 'r'
rGrid.units = 'a'
rvalues = 0.1 + (0.8/nradial)*(0.5+FINDGEN(nradial))
rGrid.values = PTR_NEW(rValues)
rmin = MIN(rvalues, max=rmax)
rGrid.boundary = 'open'
rGrid.uniform = 1B
rGrid.range = [0.1, 0.9]
rGrid.irange = [0, (nradial-1)]
FOR i=0, N_TAGS({GKVsd})-1 DO thisObj.(i) = str1D.(i)
thisObj.mnemonic = mnemonics[nother]
thisObj.Title = titles[nother]
Indices = REPLICATE('*', 2)
thisObj.indices = PTR_NEW(Indices)
values = f[(nother):(nother+nradial), *]
thisObj.values = PTR_NEW(values)
thisObj.units = units[nother]
thisObj.Grid1 = GKVsd_GridCopy(rGrid)
thisObj.Grid2 = GKVsd_GridCopy(tGrid)
gkvObjs[iObj] = OBJ_NEW('GKVs2D', thisObj)
iObjStr = STRING(iObj, FORMAT='(I2)')
iObjStr = STRTRIM(iObjStr, 2)
commandStr = commandStr + mnemonics[nother] + ':gkvObjs[' + iObjStr + ']}
;
; Create structure 'result' 
;
ok = EXECUTE(commandStr)
;
; Turn off error control
;
IF(debug EQ 0) THEN CATCH, /CANCEL
;
; ... and we're done!
;
RETURN, result
END ; ****** GTC_History ****** ;




