PRO GKVs3D::MPeg,Quality=quality, _EXTRA=extra
;
;
; Purpose:
;
;		This method makes an MPeg animation from a sequence of 
;		images of slices of 'self' at sequence of values the independent
;		variable stored in self.Grid3.  It creates an image of
;		every iSkipth timeslice of 'self', and stores
;		the result into a file selected by the user.
;
; Arguments:	NONE
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
;	Shade_Surf		Set this KeyWord (i.e., '/Shade_Surf') to
;				make a surface plot instead of an image.
;				Defaults to making 'image' plots (Optional)
;
;	Path			Set to path to folder where the MPeg
;				animation is to stored.  Defaults to the 
;				current working directory (Optional).
;
;	FileName		Name of the file containing the MPeg
;				animation.  Defaults to 'self.mnemonic'.mpg
;				(Optional)
;
;	Xsize, Ysize		The size of the frame in 'device' pixels.
;				Defaults to xsize = 500, ysize = 500.
;				(Optional)
;
;	Quality			Set this keyword to an integer value between 
;				0 (low quality) and 100 (high quality) inclusive 
;				to specify the quality at which the MPEG stream 
;				is to be stored. Higher quality values result in 
;				lower rates of time compression and less motion 
;				prediction which provide higher quality MPEGs but 
;				with substantially larger file size. Lower quality 
;				factors may result in longer MPEG generation times. 
;				The default is 50. (Optional)
;
;	true			This keyword specifies the index, starting at 1, 
;				of the dimension over which the color is interleaved. 
;				For example, for an image that is pixel interleaved 
;				and has dimensions of (3, m, n), set TRUE to 1. 
;				Specify 2 for row-interleaved images (m, 3, n); and 3 for 
;				band-interleaved images (m, n, 3). If this keyword
;				is not set, then WRITE_MPEG is called with true=1.
;				(Optional)
;
;	ShowLoad		Set this keyword to view images as they are loaded 
;				into the MPeg file.  Default is not to show images.
;				(Optional)
;							
;
; Written by W.M. Nevins
;	3/15/01
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
trange = self.Grid3.range
result = GetKeyWord('trange', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	IF(  (N_Elements(result) EQ 2) AND ( Query_Real(result) OR Query_Integer(result) )  ) THEN trange = result
ENDIF
;
;	'mnemonic' = trange
;
result = GetKeyWord(self.grid3.mnemonic, extra)
IF(TypeOF(result) NE 7) THEN BEGIN
	IF(  (N_Elements(result) EQ 2) AND ( Query_Real(result) OR Query_Integer(result) )  ) THEN trange = result
ENDIF
; 
; Turn trange into a range of indices
;
gridValues = *self.Grid3.values
irange = LONARR(2)
FOR i=0, 1 DO BEGIN
	temp = (gridValues - trange[i])^2
	eps = MIN(temp, indx)
	irange[i] = indx
ENDFOR
;
;	Shade_Surf=shadeSurf
;
shadeSurf = 0
result = GetKeyWord('Shade_Surf', extra)
IF( Query_Integer(result) ) THEN shadeSurf = result
;
;	ShowLoad = showLoad
;
showLaod = 0
result = GetKeyWord('ShowLoad', extra)
IF( Query_Integer(result) ) THEN BEGIN
	IF(result GT 0) THEN showLoad = 1
ENDIF
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
;	FileName = FileName
;
FileName = self.mnemonic + ".mpg"
FileName = STRCOMPRESS(FileName, /REMOVE_ALL)
result = GetKeyWord('FileName', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN	; User provided text
	FileName = STRCOMPRESS(result, /REMOVE_ALL)			; remove all blanks from input text
	subStrings = STRSPLIT(FileName, '.', /EXTRACT)			; split user-provided name into substrings on '.'
	nSubStrings = N_ELEMENTS(subStrings)				; Check that user supplied a '.mpg' suffix
	IF(SubStrings[nsubStrings-1] NE 'mpg') THEN FileName=FileName + '.mpg'		; or add one if necessary
ENDIF
;
; Check if file 'FileName' alreaedy exists
;
fileNames = FINDFILE(FileName, COUNT=nFIles)
fileIndex = 0
subStrings = STRSPLIT(FileName, '.', /EXTRACT)
nSubStrings = N_ELEMENTS(subStrings)
filePrefix = STRJOIN(subStrings[0:nSubStrings-2]) + '_'
WHILE nFIles GT 0 DO BEGIN						; fileName already exists, so embed an index number in
	fileIndex = fileIndex + 1					; into FileName
	FileName = filePrefix + STRCOMPRESS(STRING(fileIndex), /REMOVE_ALL) + '.mpg'
	fileNames = FINDFILE(FileName, COUNT=nFIles)
ENDWHILE
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
; Open MPeg file
;
mpegID = MPEG_OPEN([xsize,ysize], FILENAME=FileName, QUALITY=quality)
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
; Save info on device, colors
;
IF KEYWORD_SET(ShowLoad) THEN BEGIN
	WINDOW, /FREE, XSIZE=xSize, YSIZE=ySize, title='Loading MPEG FIle'
	thisWindow = !D.WINDOW
ENDIF
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
; begin writting frames into Z-Buffer (for transfer into MPEG file)
;
frame=-1
FOR i=irange[0], irange[1], iskip DO BEGIN
	frame = frame+1				; Increment frame number

	IF( KEYWORD_SET(showLoad) ) THEN BEGIN
		SET_PLOT, thisDevice
		!P.BACKGROUND = 0
		!P.COLOR = 1
		TVLCT, rr, gg, bb
		ERASE
		IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
			self -> Shade_Surf,	indx1=i, _Extra=nextra
		ENDIF ELSE BEGIN
			self -> Draw, 		indx1=i, _Extra=nextra
		ENDELSE
		SET_PLOT, 'Z'
		!P.BACKGROUND = 0
		!P.COLOR = 1
		TVLCT, rr, gg, bb
		DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
	ENDIF
	
	IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
		self -> Shade_Surf,	indx1=i, _Extra=nextra
	ENDIF ELSE BEGIN
		self -> Draw, 		indx1=i, _Extra=nextra
	ENDELSE
	
	;
	; Get 'Image24' from Z-Buffer
	;
	IF(trueColor NE 1) THEN BEGIN
		thisImage = TVRD()
		TVLCT, r, g, b, /GET	
		image24 = BYTARR(3, xSize, ySize)
		image24(0,*,*) = r(thisImage)
		image24(1,*,*) = g(thisImage)
		image24(2,*,*) = b(thisImage)
		MPEG_PUT, mpegID, IMAGE=image24, FRAME=frame, /ORDER
	ENDIF ELSE BEGIN
		thisImage = TVRD(TRUE=1)
		MPEG_PUT, mpegID, IMAGE=image24, FRAME=frame, /ORDER
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
; Save MPeg animation into FileName
;
MPEG_SAVE, mpegID, FILENAME=FileName
;
; Now, close the MPeg sequence
;
MPEG_CLOSE, mpegID
;
; Close 'showload' window if 'showload; was set
;
IF KEYWORD_SET(ShowLoad) THEN WDELETE, thisWindow
;
; Change back to 'current_working_directory'
;
CD, current_working_directory
;
; and we are done
;
RETURN
END ; ****** GKVs3D::MPeg ****** ;

