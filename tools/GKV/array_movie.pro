Pro Array_Movie, arr, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, arr9, _Extra=extra

;
; Argument is an array of GKV objects to be animated over index of the array
; Written so arguments are as similar to movie.pro as possible
; 
; E Wang 4-20-10
;




irange=(SIZE(arr))[1]
imin = 0
imax=irange-1
result = GetKeyWord('imin', extra)
IF(TypeOf(result) NE 7) THEN imin=result
result = GetKeyWord('imax', extra)
IF(TypeOf(result) NE 7) THEN imax=result
iskip=(imax-imin)/128 + 1
result = GetKeyWord('iskip', extra)
IF(TypeOf(result) NE 7) THEN iskip = result > iskip
timeMnemonic = 't'

nPixels=[400, 400]
result = GetKeyWord('nPixels', extra)
IF(typeOf(result) NE 7) THEN BEGIN
  CASE N_ELEMENTS(result) OF
    1:  nPixels = [FIX(result),    FIX(result)]
    2:  nPixels = [FIX(result[0]), FIX(result[1])] 
    ELSE: MESSAGE, "Couldn't parse 'nPixels = ...'.  Will use defaults", /INFORMATIONAL
  ENDCASE
ENDIF
vrange=[-.5,.5]  ; FIX ME!--------------------------------------------------------
numoPlots = 0

IF (SIZE(arr))[0] EQ 1 THEN BEGIN
IF(N_PARAMS() GT 0) THEN BEGIN
  tmin = imin
  tmax = imax
  oPlotObjs = arr[0]  
  FOR i=1, N_PARAMS()-1  DO BEGIN    
    argString = STRING(i, FORMAT='(I1)')
    commandString = "oPlotObjs = [oPlotObjs, arr" + argString + "]"
    ok = EXECUTE(commandString)
  ENDFOR
  numoPlots = N_ELEMENTS(oPlotObjs)  
  FOR i=0, numoPlots-1 DO BEGIN
    ;oPlotObjs[i] -> Get, axis=2, range=thisRange
    ;oPlotObjs[i] -> SignalWindow, axis=2, range=[tmin,tmax]
    oPlotObjs[i] -> Get, vrange=thisVrange
    ;oPlotObjs[i] -> SignalWindow, axis=2, range=thisRange
    vrange[0] = vrange[0] < thisVrange[0]
    vrange[1] = vrange[1] > thisVrange[1]
  ENDFOR
ENDIF

nframes=(imax - imin)/iskip + 1

XINTERANIMATE, SET=[nPixels[0], nPixels[1], nframes], /SHOWLOAD
iframe=0
FOR i=imin, imax, iskip DO BEGIN  
  arr[i] -> Draw, _EXTRA=extra  
  IF(numoPlots GT 0) THEN BEGIN    
    FOR j=0, N_PARAMS()-2 DO BEGIN     
    ;(oPlotObjs[j])[i] -> oPlot, color=((j) MOD 10)+2
      commandString = "arr" + STRING(j+1, FORMAT='(I1)') +"[i] -> oPlot, color=((j) MOD 10) + 2"
      ok = EXECUTE(commandString)
    ENDFOR
  ENDIF
  XINTERANIMATE, FRAME=iframe, Window=!D.WINDOW
  iframe=iframe+1
ENDFOR
XINTERANIMATE

ENDIF ELSE BEGIN   ;  2D array given


imax = (SIZE(arr))[2] -1
noPlots = (SIZE(arr))[1] 

nframes=(imax - imin)/iskip + 1

XINTERANIMATE, SET=[nPixels[0], nPixels[1], nframes], /SHOWLOAD
iframe=0
FOR i=imin, imax, iskip DO BEGIN  
  arr[0,i] -> Draw, _EXTRA=extra    
    FOR j=1, noPlots-1 DO BEGIN         
    arr[j,i] -> oplot, color=((j)/10 MOD 10) + 1      
    ENDFOR  
  XINTERANIMATE, FRAME=iframe, Window=!D.WINDOW
  iframe=iframe+1
ENDFOR
XINTERANIMATE

ENDELSE

END ; ****** array_movie ****** ;

PRO array_jpegs, arr, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, arr9, _Extra=extra
;
;  JPEGS version of array_movie
;
;  KEYWORDS:
;     iskip, trange, directory_name, fileroot
;
;
;
;

separator="/"         ; Path separator for unix devices
IF (!D.Name EQ "MAC") then separator=":"  ; or for MAC's
IF (!D.NAME EQ "WIN") THEN speparator="\" ; or for windows systems
;
; Check for keywords in 'extra'
;
; iSkip=skip:
;
iskip = 1L
result = GetKeyWord('iskip', extra)
IF(TypeOF(result) NE 7) THEN iskip =LONG(result) > 0L
;
; trange = [tmin, tmax]
;
trange=[0,(SIZE(arr))[2]]
tmin = 0
tmax=trange-1
result = GetKeyWord('trange', extra)
IF(TypeOf(result) NE 7) THEN BEGIN
  IF(  (N_Elements(result) EQ 2) AND ( Query_Real(result) OR Query_Integer(result) )  ) THEN trange = result
ENDIF


irange = trange
;
; Shade_Surf=shadeSurf
;
shadeSurf = 0
result = GetKeyWord('Shade_Surf', extra)
IF( Query_Integer(result) ) THEN shadeSurf = result
;
; ShowLoad = showLoad
;
showLoad = 0
result = GetKeyWord('ShowLoad', extra)
IF( Query_Integer(result) ) THEN showLoad = result
;
; Path=path 
;
cd, current = current_working_directory
path = 'current_working_directory'
result = GetKeyWord('path', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN
  ;
  ; Change to this directory (if it exists)
  ;
  IF( STRCMP(result, current_working_directory) ) THEN GOTO, GotPath  ; 'result' IS current_working_directory
  CD, result
  CD, CURRENT=path
  IF( STRCMP(path, current_working_directory) ) THEN BEGIN    ; 'result' is not a valid directory identifier
    MESSAGE, "specified path is not legal", /INFORMATIONAL
    RETURN  
  ENDIF
  CD, path  ; Have a good path
ENDIF
GotPath: 
;
; DirectoryName = Directory_name
;

arr[0,0] -> get, mnemonic=filemnemonic
Directory_Name = filemnemonic + '_' + STRMID(SYSTIME(0), 4, 3)  + "_" + STRMID(SYSTIME(0), 8, 2)  + "_" + STRMID(SYSTIME(0),20,4)
result = GetKeyWord('DirectoryName', extra)
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN  ; User provided text
  Directory_Name = STRCOMPRESS(result, /REMOVE_ALL)   ; remove all blanks from input text
ENDIF
IF(STRCMP(Directory_Name, 'none', /FOLD_CASE) EQ 0) THEN BEGIN
  FILE_MKDIR, Directory_Name
  CD, Directory_Name
ENDIF
;
; File_Root = FileRoot
;
FileRoot = STRCOMPRESS(filemnemonic, /REMOVE_ALL)      ; Set default to self.mnemonic
result = GetKeyWord('FileRoot', extra)          ; Check command line for 'FileRoot' keyword
IF(  (TypeOf(result) EQ 7) AND (result NE 'undefined')  ) THEN BEGIN  ; User provided text
  FileRoot = STRCOMPRESS(result, /REMOVE_ALL)     ; remove all blanks from input text
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
;   Xsize = xsize, Ysize = ysize
;
xSize = 400
result = GetKeyWord('xsize', extra)
IF(Query_Integer(result)) THEN xSize = result > 100
ySize = 400
result = GetKeyWord('ysize', extra)
IF(Query_Integer(result)) THEN ySize = result > 100
;
; Quality = quality
;
quality = 100
result = GetKeyWord('Quality', extra)
IF Query_Integer(result) THEN BEGIN
  quality = result < 100
  quality = quality > 0
ENDIF
;
; True = true
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
FOR i=irange[0], irange[1]-1, iskip DO BEGIN
  FileIndexStr = STRING(fileIndex, FORMAT='(I5)')   ; Turn fileIndex into a string
  FileIndexStr = STRTRIM(FileIndexStr, 2)     ; Strip out leading (and trailing) blanks
  digits = STRLEN(FileIndexStr)       ; Determine number of digits
  STRPUT, FileRoot, FileIndexStr, FileRootLen-digits  ; Overwrite tail of 'FileRoot' with new sequence number
  fileIndex = fileIndex + 1       ; Increment fileIndex
  FileName = FileRoot + ".jpg"        ; Append ".tif" to FileName
  FileName = STRCOMPRESS(FileName, /REMOVE_ALL)   ; Strip out any blanks
  PRINT, 'Working on image number ', i+1, ' of ', irange[1]

  IF( KEYWORD_SET(showLoad) ) THEN BEGIN
    SET_PLOT, thisDevice
    !P.BACKGROUND = 0
    !P.COLOR = 1
    ERASE
    IF (SIZE(arr))[0] EQ 1 THEN BEGIN
      IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
        arr[0]  -> Shade_Surf, indx1=i, _Extra=nextra
      ENDIF ELSE BEGIN
        arr[0] -> Draw,     indx1=i, _Extra=nextra
      ENDELSE
    ENDIF ELSE BEGIN
      
      IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
        arr[0]  -> Shade_Surf, indx1=i, _Extra=nextra        
      ENDIF ELSE BEGIN
        arr[0,i] -> Draw,     indx1=i, _Extra=nextra        
        FOR j = 1,(SIZE(arr))[2]-2  DO BEGIN
        arr[j,i] -> oplot, indx1=i, color=((j)/10 MOD 10) + 1, _Extra=nextra
        ENDFOR        
      ENDELSE  
    ENDELSE
    SET_PLOT, 'Z'
    !P.BACKGROUND = 0
    !P.COLOR = 1
    TVLCT, rr, gg, bb
    DEVICE, SET_RESOLUTION=[xsize, ysize], SET_COLORS=nColors
  ENDIF
  
  IF (SIZE(arr))[0] EQ 1 THEN BEGIN
      IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN       
        arr[0]  -> Shade_Surf, indx1=i, _Extra=nextra
      ENDIF ELSE BEGIN
        arr[0] -> Draw,     indx1=i, _Extra=nextra
      ENDELSE
  ENDIF ELSE BEGIN            
      
      IF( KEYWORD_SET(ShadeSurf) ) THEN BEGIN
        arr[0]  -> Shade_Surf, indx1=i, _Extra=nextra        
      ENDIF ELSE BEGIN
        arr[0,i] -> Draw,     indx1=i, _Extra=nextra
        FOR j = 0,(SIZE(arr))[1]-2 DO BEGIN
        arr[j+1,i] -> oplot, psym=3, indx1=i, color=((j)/10 MOD 10) + 1, _Extra=nextra      
      ENDFOR
      ENDELSE      
  ENDELSE  
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
END  ; ********  array_jpegs ***** ; 


