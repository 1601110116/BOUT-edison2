PRO GKVs2D::JPEGs, Obj0, obj1, obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9, _EXTRA=extra
;
; Purpose:
;
;		This routine is used to provide a sequence of JPEG
;		images for an animation of 'self' over the independent
;		variable stored in self.Grid2.  It creates a JPEG image of
;		every iSkipth timeslice of 'self', and stores each frame of
;		the result as a separate Jpeg file in a directory
;		selected by the user.
;
; Arguments:	
;
;		Objects (or object arrays) to be plotted (in various colors)
;		over the curve corresponding to "self".
;
; Input KeyWords:
;
;	iSkip			Interval between timeslices to be animated.
;				Defaults to 1 (Optional).
;
;	trange			A two-element array specifying the first and 
;				last time-slices of this animation.  
;				Defaults to self.Grid3.range (Optional).
;
;	'mnemonic'		The mnemonic of self.Grid3 can be used in
;				place of the keyword 'trange' with the same 
;				effect.  Defaults to self.Grid3.range 
;				(Optional).
;
;	Path			Set to path to folder where a new folder
;				containing the sequence of JPEG files is to
;				be stored.  Defaults to the current working
;				directory (Optional).
;
;	DirectoryName		Name of the new folder to be created to
;				contain the sequence of JPEG files.  Enter
;				DirectoryName='none' to store JPEG files 
;				in folder specified by 'Path' keyword 
;				(that is, NOT in a newly created directory)  
;				Defaults to "'Self.mnemonic'_Animation".
;				(Optional)
;
;	FileRoot		Sequence of JPEG files will have the names
;				FileRoot00000.jpg, FILERoot00001.jpg, ...
;				Defaults to self.mnemonic.  (Optional)
;
;	FirstIndex		Index of first JPEG file to be produced.
;				defaults to 00000. (Optional)
;
;	Xsize, Ysize		The size of the frame in 'device' pixels.
;				Defaults to xsize = 500, ysize = 500.
;				(Optional)
;
;	Quality			Sets the quality index in the range of 
;				0 ("terrible") to 100 ("excellent") of the 
;				JPEG images to be produced.  The default is
;				80 to produce "pretty good" quality.  Lower
;				values of "Quality" pjroduce higher compression
;				ratios and smaller files.  Default is Quality = 80.
;				(Optional).
;
;	true			This keyword specifies the index, starting at 1, 
;				of the dimension over which the color is interleaved. 
;				For example, for an image that is pixel interleaved 
;				and has dimensions of (3, m, n), set TRUE to 1. 
;				Specify 2 for row-interleaved images (m, 3, n); and 3 for 
;				band-interleaved images (m, n, 3). If this keyword
;				is not set, then WRITE_JPEG is called with true=1.
;				(Optional)
;
;	ShowLoad		Set this keyword to view images as they are loaded 
;				into the JPEG file.  Default is not to show images.
;				(Optional)
;							
;
; Written by W.M. Nevins
;	1/8/05
;
separator="/"					; Path separator for unix devices
IF (!D.Name EQ "MAC") then separator=":"	; or for MAC's
IF (!D.NAME EQ "WIN") THEN speparator="\"	; or for windows systems
;
; Check for keywords in 'extra'
;
;	iSkip=skip:
;
iskip = 1L
result = GetKeyWord('iskip', extra)
IF(TypeOF(result) NE 7) THEN iskip =LONG(result) > 0L
;
;	trange = [tmin, tmax]
;
trange = self.Grid2.range
result = GetKeyWord('trange', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	IF(  (N_Elements(result) EQ 2) AND ( Query_Real(result) OR Query_Integer(result) )  ) THEN trange = result
ENDIF
;
;	'mnemonic' = trange
;
result = GetKeyWord(self.grid2.mnemonic, extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	IF(  (N_Elements(result) EQ 2) AND ( Query_Real(result) OR Query_Integer(result) )  ) THEN trange = result
ENDIF
tmin = trange[0]
tmax = trange[1]
; 
; Turn trange into a range of indices
;
gridValues = *self.Grid2.values
irange = LONARR(2)
FOR i=0, 1 DO BEGIN
	temp = (gridValues - trange[i])^2
	eps = MIN(temp, indx)
	irange[i] = indx
ENDFOR
imin = irange[0]
imax = irange[1]
;
;	ShowLoad = showLoad
;
showLoad = 0
result = GetKeyWord('ShowLoad', extra)
IF( Query_Integer(result) ) THEN showLoad = result
;
;	Path=path 
;
cd, current = current_working_directory
path = 'current_working_directory'
result = GetKeyWord('path', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN
	;
	; Change to this directory (if it exists)
	;
	IF( STRCMP(result, current_working_directory) ) THEN GOTO, GotPath	; 'result' IS current_working_directory
	CD, result
	CD, CURRENT=path
	IF( STRCMP(path, current_working_directory) ) THEN BEGIN		; 'result' is not a valid directory identifier
		MESSAGE, "specified path is not legal", /INFORMATIONAL
		RETURN	
	ENDIF
	CD, path	; Have a good path
ENDIF
GotPath: 
;
;	DirectoryName = Directory_name
;
Directory_Name = Self.mnemonic + "_Animation"
result = GetKeyWord('DirectoryName', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN	; User provided text
	Directory_Name = STRCOMPRESS(result, /REMOVE_ALL)		; remove all blanks from input text
ENDIF
IF(STRCMP(Directory_Name, 'none', /FOLD_CASE) EQ 0) THEN BEGIN
	FILE_MKDIR, DIrectory_Name
	CD, Directory_Name
ENDIF
;
;	File_Root = FileRoot
;
FileRoot = STRCOMPRESS(self.mnemonic, /REMOVE_ALL)			; Set default to self.mnemonic
result = GetKeyWord('FileRoot', extra)					; Check command line for 'FileRoot' keyword
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN	; User provided text
	FileRoot = STRCOMPRESS(result, /REMOVE_ALL)			; remove all blanks from input text
ENDIF
;
; Add sequence number field to File_Root
;
FileRoot = FileRoot + '00000'
FileRoot = STRCOMPRESS(FileRoot, /REMOVE_ALL)
FileRootLen = STRLEN(FileRoot)
;
; Check for keyword FirstIndex
;
fileIndex = 0
result = GetKeyWord('FirstIndex', extra)
IF Query_Integer(result) THEN fileIndex = result
;
; 	Xsize = xsize, Ysize = ysize
;
xSize = 500
result = GetKeyWord('xsize', extra)
IF(Query_Integer(result)) THEN xSize = result > 100
ySize = 500
result = GetKeyWord('ysize', extra)
IF(Query_Integer(result)) THEN ySize = result > 100
;
;	Quality = quality
;
quality = 80
result = GetKeyWord('Quality', extra)
IF Query_Integer(result) THEN BEGIN
	quality = result < 100
	quality = quality > 0
ENDIF
;
;	True = true
;
true = 1
result = GetKeyWord('true', extra)
IF(Query_Integer(result)) THEN BEGIN
	true = true < 3
	true = true > 1
ENDIF
;
; Parse command line for over-plots (that is, arguments)
;
vrange = self.vrange
numoPlots = 0
IF(N_PARAMS() GT 0) THEN BEGIN
	tmin = (*self.Grid2.Values)[imin]
	tmax = (*self.Grid2.Values)[imax]
	oPlotObjs = obj0
	FOR i=1, N_PARAMS()-1  DO BEGIN
		argString = STRING(i, FORMAT='(I1)')
		commandString = "oPlotObjs = [oPlotObjs, Obj" + argString + "]"
		ok = EXECUTE(commandString)
	ENDFOR
	numoPlots = N_ELEMENTS(oPlotObjs)
	FOR i=0, numoPlots-1 DO BEGIN
		oPlotObjs[i] -> Get, axis=2, range=thisRange
		oPlotObjs[i] -> SignalWindow, axis=2, range=[tmin,tmax]
		oPlotObjs[i] -> Get, vrange=thisVrange
		oPlotObjs[i] -> SignalWindow, axis=2, range=thisRange
		vrange[0] = vrange[0] < thisVrange[0]
		vrange[1] = vrange[1] > thisVrange[1]
	ENDFOR
ENDIF

;
; Save info on device, colors
;
thisDevice = !D.NAME
nColors = !D.TABLE_SIZE
TVLCT, rr, gg, bb, /GET
;
; Render graphic in Z-buffer
;
SET_PLOT, 'Z'
!P.BACKGROUND = 0
!P.COLOR = 1
TVLCT, rr, gg, bb
ERASE, color=0
DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
;
; Visual depth of "Z-buffer" is ALWAYS 8!
;
trueColor = 0
;
; Transfer contents of "extra" to "nextra" if any un-parsed keywords are left in "extra".
;
IF(TypeOF(extra) EQ 8) THEN nextra = extra
;
; begin writting JPEGs
;
FOR i=irange[0], irange[1], iskip DO BEGIN
	FileIndexStr = STRING(fileIndex, FORMAT='(I5)')		; Turn fileIndex into a string
	FileIndexStr = STRTRIM(FileIndexStr, 2)			; Strip out leading (and trailing) blanks
	digits = STRLEN(FileIndexStr)				; Determine number of digits
	STRPUT, FileRoot, FileIndexStr, FileRootLen-digits	; Overwrite tail of 'FileRoot' with new sequence number
	fileIndex = fileIndex + 1				; Increment fileIndex
	FileName = FileRoot + ".jpg"				; Append ".tif" to FileName
	FileName = STRCOMPRESS(FileName, /REMOVE_ALL)		; Strip out any blanks

	IF( KEYWORD_SET(showLoad) ) THEN BEGIN
		SET_PLOT, thisDevice
		!P.BACKGROUND = 0
		!P.COLOR = 1
		ERASE
		thisSlice = self -> Slice(axis=2, Index=i)
		IF(TypeOF(nextra) EQ 8) THEN extra = nextra
		thisSlice -> Draw, vrange=vrange, _EXTRA=extra
		thisSlice -> Trash
		IF(numoPlots GT 0) THEN BEGIN
			FOR j=0, numoPlots-1 DO BEGIN
				thisPlot = oPlotObjs[j] -> Slice(axis=2, index=i)
				thisPlot -> oPlot, color=((j) MOD 10)+2
				thisPlot -> Trash
			ENDFOR
		ENDIF
		SET_PLOT, 'Z'
		!P.BACKGROUND = 0
		!P.COLOR = 1
		TVLCT, rr, gg, bb
		DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
	ENDIF

	thisSlice = self -> Slice(axis=2, Index=i)
	IF(TypeOF(nextra) EQ 8) THEN extra = nextra
	thisSlice -> Draw, vrange=vrange, _EXTRA=extra
	thisSlice -> Trash
	IF(numoPlots GT 0) THEN BEGIN
		FOR j=0, numoPlots-1 DO BEGIN
		thisPlot = oPlotObjs[j] -> Slice(axis=2, index=i)
		thisPlot -> oPlot, color=((j) MOD 10)+2
		thisPlot -> Trash
		ENDFOR
	ENDIF
	
	;
	; Write 'thisImage to the JPEG file.
	;
	IF(trueColor NE 1) THEN BEGIN
		thisImage = TVRD()
		TVLCT, r, g, b, /GET	
		image24 = BYTARR(3, xSize, ySize)
		image24(0,*,*) = r(thisImage)
		image24(1,*,*) = g(thisImage)
		image24(2,*,*) = b(thisImage)
		WRITE_JPEG, FileName, image24, TRUE=1, QUALITY=quality
	ENDIF ELSE BEGIN
		thisImage = TVRD(TRUE=1)
		WRITE_JPEG, FileName, thisimage, QUALITY=quality, TRUE=true
	ENDELSE
	ERASE, COLOR=0
ENDFOR
;
; Get out of Z-Graphics buffer
;
SET_PLOT, thisDevice
!P.BACKGROUND=0
!P.COLOR=1
TVLCT, rr, gg, bb
;
; Change back to 'current_working_directory'
CD, current_working_directory
;
; and we're done ...
;
RETURN
END	; ****** GKVs2D::JPEGs ****** ;
