PRO GKVs3D::TIFFs, _EXTRA=extra
;
; Purpose:
;
;		This routine is used to provide a sequence of TIFF
;		images for an animation of 'self' over the independent
;		variable stored in self.Grid3.  It creates a TIFF image of
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
;	Path			Set to path to folder where a new folder
;				containing the sequence of TIFF files is to
;				be stored.  Defaults to the current working
;				directory (Optional).
;
;	DirectoryName		Name of the new folder to be created to
;				contain the sequence of TIFF files.  
;				Defaults to "TIFF_Animation".
;				(Optional)
;
;	File_Root		Sequence of TIFF files will have the names
;				File_Root00001, FILE_Root00002, ...
;				Defaults to "TIFF".  (Optional)
;
;	Xsize, Ysize		The size of the frame in 'device' pixels.
;				Defaults to xsize = 500, ysize = 500.
;				(Optional)
;
;	ShowLoad		Set this keyword to view images as they are loaded 
;				into the TIFF file.  Default is not to show images.
;				(Optional)
;							
;
; Written by W.M. Nevins
;	9/20/00
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
Directory_Name = "TIFF_Animation"
result = GetKeyWord('DirectoryName', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN	; User provided text
	Directory_Name = STRCOMPRESS(result, /REMOVE_ALL)			; remove all blanks from input text
ENDIF
FILE_MKDIR, DIrectory_Name
CD, Directory_Name
;
;	File_Root = File_Root
;
File_Root = "TIFF"
result = GetKeyWord('File_Root', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN	; User provided text
	File_Root = STRCOMPRESS(result, /REMOVE_ALL)			; remove all blanks from input text
ENDIF
;
; Set up name of first file
;
FileName = File_Root + '00000'
FileName = STRCOMPRESS(FileName, /REMOVE_ALL)
FileNameLen = STRLEN(FileName)
fileIndex = 0
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
; Check visual depth
;
trueColor = 0
DEVICE, Get_Visual_Depth=thisDepth            
IF(thisDepth gt 8) THEN trueColor = 1 
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
; begin writting TIFFs
;
FOR i=irange[0], irange[1], iskip DO BEGIN
		fileIndex = fileIndex + 1				; Increment fileIndex
		FileIndexStr = STRING(fileIndex, FORMAT='(I5)')		; Turn fileIndex into a string
		FileIndexStr = STRTRIM(FileIndexStr, 2)			; Strip out leading (and trailing) blanks
		digits = STRLEN(FileIndexStr)				; Determine number of digits
		STRPUT, FileName, FileIndexStr, FileNameLen-digits	; Overwrite tail of 'FileName' with new sequence number

	IF( KEYWORD_SET(showLoad) ) THEN BEGIN
		SET_PLOT, thisDevice
		!P.BACKGROUND = 0
		!P.COLOR = 1
		ERASE
		IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
			self -> Shade_Surf,	indx1=i, _Extra=extra
		ENDIF ELSE BEGIN
			self -> Draw, 		indx1=i, _Extra=extra
		ENDELSE
		SET_PLOT, 'Z'
		!P.BACKGROUND = 0
		!P.COLOR = 1
		TVLCT, rr, gg, bb
		DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
	ENDIF
	IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
		self -> Shade_Surf,	indx1=i, _Extra=extra
	ENDIF ELSE BEGIN
		self -> Draw, 		indx1=i, _Extra=extra
	ENDELSE
	thisImage = TVRD()
	TVLCT, r, g, b, /GET
	;
	; Write 'thisImage' to the TIFF file.
	;
	IF(trueColor EQ 1) THEN BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, 		$
				RED=r, GREEN=g, BLUE=b, 					$
				XRESOL=ROUND(!D.X_PX_CM * 2.54), YRESOL=ROUND(!D.X_PX_CM * 2.54)
	ENDIF ELSE BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, 		$
				XRESOL=ROUND(!D.X_PX_CM * 2.54), YRESOL=ROUND(!D.X_PX_CM * 2.54)
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
END	; ****** GKVs3D::TIFFs ****** ;
