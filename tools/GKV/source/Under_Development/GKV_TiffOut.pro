PRO GKVsd::DRAW, _Extra=extra
;
;	Purpose:
;
; 		writes contents of self into currently 
;		open graphics frame
;
;
; 'Virtual compile-time keywords:
;			  title=T_title, Pretty=pretty
;
; Default "Plotting" routine for GKVsd objects (0-D).
; (really just writes values (and error bars if present)
; to center of graphics frame.
;
; First, parse 'extra' for keywords which we wish to intercept 
; (this is the moral equivalent of defining them on the PRO line...
;  but prevents IDL from enforcing keyword abreviation, which would
;  interfer with use of mnemonics as 'run-time' keywords)
;
result = GetKeyWord('title', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN T_title=result	;	title = T_title
result = GetKeyWord('Pretty', extra)
IF(TypeOf(result) NE 7) THEN pretty=result			;	Pretty = pretty
result = GetKeyWord('Thick', extra)				;	Thick = thick
IF(TypeOf(result) NE 7) THEN thick = result
;
; Erase current frame
;
ERASE, 0
;
; Choose an appropriate character size
; (scaled to size of currently open graphics window)
;
newSizes = GKV_CharSize(Old_Sizes = oldSizes)
;
; Reset character size
;
DEVICE, SET_CHARACTER_SIZE = newSizes
;
; Put title at top of Graphics window
;
indices = self -> IndexString([0,1], pretty=pretty)
indexStr = '[' + STRJOIN(indices, ', ') + ']'
IF(KEYWORD_SET(pretty)) THEN BEGIN
	title=self.title + indexStr + " (" + self.units + ")"
ENDIF ELSE BEGIN
	title=self.mnemonic + indexStr + " (" + self.units + ")"
ENDELSE
IF KEYWORD_SET(T_title) THEN title=T_title
charSize = GKV_TitleSize(title, Width = (0.90-0.15), yPosition=ywrite)
xwrite = 0.5*!D.X_SIZE
XYOUTS, xwrite, ywrite, title, ALIGNMENT=0.5, CHARSIZE=charSize, /DEVICE 
;
; Write value of 'self' in center of screen
;
IF(KEYWORD_SET(pretty)) THEN BEGIN
	stringOut = self.title    + indexStr + ' = ' 
ENDIF ELSE BEGIN
	stringOut = self.mnemonic + indexStr + ' = ' 
ENDELSE
value = *self.values
IF( Query_Complex(value) ) THEN BEGIN
	rValue = FLOAT(value)
	iValue = IMAGINARY(value)
	valueString = 	 '( ' + STRCOMPRESS(STRING(rValue, FORMAT='(G10.3)'), /REMOVE_ALL) +	$
			' , ' + STRCOMPRESS(STRING(iValue, FORMAT='(G10.3)'), /REMOVE_ALL) + ' )'
ENDIF ELSE BEGIN
	valueString = STRCOMPRESS(STRING(value, FORMAT='(G10.3)'), /REMOVE_ALL)
ENDELSE
stringOut = stringOut + valueString
IF(PTR_VALID(self.errorBars)) THEN BEGIN
	errorBar = *self.ErrorBars
	errorString = STRCOMPRESS(STRING(errorBar, FORMAT='(G10.3)'), /REMOVE_ALL)
	stringOut = stringOut + ' � ' + errorString
ENDIF
stringOut = stringOut + " (" + self.units + ")"
charSize = GKV_TitleSize(stringOut, Width = 0.8, yPosition=0.5)
XYOUTS, 0.5, 0.5, /NORMAL, COLOR=1, ALIGNMENT=0.5, stringOut, CHARSIZE=charSize, /DEVICE 
;
; Now write code and run info in lower corners
;
maxChars =ROUND(0.4*!D.X_SIZE/!D.X_CH_SIZE)		; Compute maximum allowed characters in a line
xwrite=!D.x_ch_size					; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower left-hand corner of plot window.
codeName = STRMID(self.CodeName, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, CodeName, /Device		; Write CodeName to lower left-hand corne.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
codePI = STRMID(self.CodePI, 0, maxChars)		; Truncate 'CodePI' if necessary
XYOUTS, xwrite, ywrite, CodePI, /Device			; Write CodePI below CodeName.
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window.
runID = STRMID(self.RunID, 0, maxChars)			; Truncate 'runID' if necessary
XYOUTS, xwrite, ywrite, RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
fileID = STRMID(self.FileID, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID
;
; Return character size to input values
;
DEVICE, SET_CHARACTER_SIZE = oldSizes

RETURN
END ; ****** GKVsd::DRAW ****** ;


PRO GKV_TiffOut, argIn, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, 	$
			PACK=pack, fileName=fileNameIn, Path=Path, 		$
			Append=append, xSize=xSizeIn, ySize=ySizeIn
;
; Purpose:
;
;	This routine accepts GKV objects, GKV object arrays,
;	or structures containing GKV objects as input.  It
;	produces a Tiff output file which displays each of
;	these GKV objects.
;
;  Keywords:
;
;	Pack	Setting this keyword (i.e., putting '/pack/ on the command line)
;		will result in 4 plots per page in the output TIFF file.
;		Default is one frame per page. (Optional)
;
;   FileName	Set this keyword to an ascii string containing the name of
;		the desired output file.  (Optional).  Defaults: 
;
;			If the argument is a single GKV object, then the 
;			fileName defaults to 'mnemonic'.tif, where 'mnenmonic'
;			is the mnemonic of the GKV object.
;		
;			If the argument is a structure which includes the tag "Name", 
;			then the outfile name defaults to 'Name'.tif, where 'Name' is
;			the (ascii) contents of the tag "Name".
;
;			If no "fileName keyword is supplied, and the structure argument
;			does not have a tag "Name", then the FileName defaults to
;			"GKV.tif".
;
; 	Path	Set this keyword to the path to the directory where the TIFF file is 
;		to be stored.  Defaults to the current working directory. (Optional)
;
;    Append	Set this keyword (i.e., put "/Append" on the command line) when
;		calling GKV_TIffOut recursively to indicate that frames are to be added
;		to an existing TIFF file.
;
;     xSize	x-dimension (in pixels) of resulting image.  Defaults to 400. (Optional)
;
;     ySize	y-dimension (in pixels) of resulting image.  Defaults to 400. (Optional)
; 
;
;  Written by W.M. Nevins
; 	7/18/02 
;
;
; Set file separator
separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'
;
; Set image size
;
xSize=400
IF(Query_Integer(xSizeIn)) THEN xSize=xSizeIn > 25
ySize=400
IF(Query_Integer(ySizeIn)) THEN ySize=ySizeIn > 25
;
START	:
argInfo=SIZE(argIn)
argDims = argInfo[0]
argType=argInfo[argDims+1]

IF(argType EQ 0) THEN BEGIN
	MESSAGE, "No more arguments. Returning", /INFORMATIONAL
	RETURN
ENDIF
;
; copy argIn to arg to avoid overwritting input
;
arg=argIn
;
; Check if argument is a pointer
;
IF(argType EQ 10) THEN BEGIN	; If arg is a pointer, then
	arg=*arg		; de-reference pointer and 
	GOTO, START		; start again.
ENDIF
;
; If argument is a Structure, replace it with an object array
;
IF(argType EQ 8) THEN BEGIN
	nTags = N_TAGS(arg)
	tagNames = TAG_NAMES(arg)
	objIndex = LONARR(nTags)
	nObjs = 0L
	FOR i=0L, nTags-1 DO BEGIN
		IF(TypeOF(arg.(i)) EQ 11) THEN BEGIN
			IF( OBJ_ISA(arg.(i), "GKVsd") ) THEN BEGIN
				objIndex[nObjs]=i
				nObjs=nObjs+1
			ENDIF
		ENDIF
		IF( STRCMP('Name',TagNames[i], /FOLD_CASE) ) THEN fileRoot=arg.Name
	ENDFOR
	IF(nObjs EQ 0) THEN BEGIN
		MESSAGE, "No GKVsd Objects in this structure", /INFORMATIONAL
		GKV_TiffOut, 	arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, 	$
				fileName=fileNameIn, Path=path, pack=pack, 		$
				Append=append, xSize=xSize, ySize=ySize
		RETURN
	ENDIF 
	temp = OBJARR(nObjs)
	FOR i=0, nObjs-1 DO temp[i] = arg.(objIndex[i])
	arg = temp
	argInfo = SIZE(arg)
	argDims = argInfo[0]
	argType = argInfo[argDims+1]
ENDIF
;
; Check that arg is an object reference
;
IF(TypeOf(arg) NE 11) THEN BEGIN
	MESSAGE, "Argument is not an Object or Structure", /INFORMATIONAL
	GKV_TiffOut, 	arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, 	$
			fileName=fileNameIn, Path=path, pack=pack, 		$
			Append=append, xSize=xSize, ySize=ySize
	RETURN	
ENDIF
;
; Remove any non-GKVsd objects from object array
;
nObjs = N_ELEMENTS(arg)
gkvObjIndex = LONARR(nObjs)
nGkvObjs = 0L
FOR i=0L, nObjs-1 DO BEGIN
	IF( OBJ_ISA(arg[i], 'GKVsd') ) THEN BEGIN
		gkvObjIndex[nGkvObjs] = i
		nGkvObjs = nGkvObjs + 1
	ENDIF
ENDFOR

IF(nGkvObjs EQ 0) THEN BEGIN
	MESSAGE, "No GKVsd Objects in this Object Array", /INFORMATIONAL
	GKV_TiffOut, 	arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, 	$
			fileName=fileNameIn, Path=path, pack=pack, 		$
			Append=append, xSize=xSize, ySize=ySize
	RETURN
ENDIF 
temp = arg[gkvObjIndex[0:nGkvObjs-1]]
arg = temp
nObjs = nGkvObjs
;
; Set 'fileRoot'
;
IF( NOT Query_String(fileRoot) ) THEN arg[0] -> GET, mnemonic=fileRoot
fileRoot = STRCOMPRESS(fileRoot, /REMOVE_ALL)
IF( STRLEN(fileRoot) EQ 0 ) THEN fileRoot = 'GKV'
;
; Check for valid fileName
;
IF(Query_Integer(FileNameIn)) THEN BEGIN
	IF(FileNameIN EQ 1) THEN BEGIN
		fileNameIn = DIALOG_PICKFILE(PATH=path, /WRITE)
		subStrings = STRSPLIT(fileNameIn, separator, /EXTRACT)
		nSubStrings = N_ELEMENTS(SubStrings)
		FileNameIn = subStrings[nSubStrings-1]
		IF(nSubStrings GT 1) THEN path=STRJOIN(subStrings[0:nSubStrings-2], separator)
	ENDIF
ENDIF

IF(Query_String(FileNameIn)) THEN BEGIN
	fileName=FileNameIn
ENDIF ELSE BEGIN
	fileName=fileRoot
ENDELSE
;
; Remove any blanks in fileName
;
fileName = STRCOMPRESS(fileName, /REMOVE_ALL)
;
; add 'tif' extension, if necessary
;
subStrings = STRSPLIT(FileName, '.', /EXTRACT)
nSubStrings = N_ELEMENTS(SubStrings)
suffix = SubStrings[nSubStrings-1]
IF(NOT STRCMP(suffix, 'tif', /FOLD_CASE) ) THEN fileName = fileName + '.tif'
;
; Change directory if a path is supplied
;
CD, CURRENT=current_working_directory
IF(N_ELEMENTS(Path)) THEN CD, path
;
; Check visual depth
;
trueColor = 0
DEVICE, GET_VISUAL_DEPTH=thisDepth            
IF(thisDepth gt 8) THEN trueColor = 1 
;
; Save info on device, colors
;
thisDevice = !D.NAME
nColors = !D.TABLE_SIZE
TVLCT, rr, gg, bb, /GET
;
; Render graphics in Z-buffer
;
SET_PLOT, 'Z'
!P.BACKGROUND = 0
!P.COLOR = 1
TVLCT, rr, gg, bb
ERASE, color=0
DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
;
; Check if this is first frame
;
IF(NOT KEYWORD_SET(Append)) THEN BEGIN
	XYOUTS, 0.5, 0.5, /NORMAL, COLOR=1, ALIGNMENT=0.5, "GKV OUTPUT" 
	thisImage = TVRD()
	TVLCT, r, g, b, /GET
	;
	; Write 'thisImage' to the TIFF file.
	;
	IF(trueColor EQ 1) THEN BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, 				$
				RED=r, GREEN=g, BLUE=b, 					$
				XRESOL=ROUND(!D.X_PX_CM * 2.54), YRESOL=ROUND(!D.X_PX_CM * 2.54 )
	ENDIF ELSE BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, 				$
				XRESOL=ROUND(!D.X_PX_CM * 2.54), YRESOL=ROUND(!D.X_PX_CM * 2.54 )
	ENDELSE
	ERASE, COLOR=0
	Append=1
ENDIF
;
; begin writting into TIFF file
;
FOR i=0, nObjs-1 DO BEGIN
	arg[i] -> Draw
	thisImage = TVRD()
	TVLCT, r, g, b, /GET
	;
	; Write 'thisImage' to the TIFF file.
	;
	IF(trueColor EQ 1) THEN BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, /APPEND,			$
				RED=r, GREEN=g, BLUE=b, 					$
				XRESOL=ROUND(!D.X_PX_CM * 2.54), YRESOL=ROUND(!D.X_PX_CM * 2.54)
	ENDIF ELSE BEGIN
		WRITE_TIFF, fileName, REVERSE(thisImage, 2), 1, /APPEND,	 		$
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
; Make next (recursive) call to GKV_TiffOut
;
GKV_TiffOut, 	arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, 	$
		fileName=fileName, Path=path, pack=pack, 		$
		Append=append, xSize=xSize, ySize=ySize
;
; and we're done ...
;
RETURN
END  ;  ****** GKV_TiffOut ******  ;