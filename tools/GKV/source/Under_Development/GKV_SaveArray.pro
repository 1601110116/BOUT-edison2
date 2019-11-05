;
; *****************************************************************************************************************
; ******************************************     Copyright Notice     *********************************************
; *                                                                                                               *
; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *
; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *
; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *
; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *
; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *
; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * 
; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *
; *                                                                                                               *
; *****************************************************************************************************************
;
; *****************************************************************************************************************
; **********************************************     DISCLAIMER     ***********************************************
; *                                                                                                               *
; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *
; *  Neither the United States Government nor the University of California nor any of their employees, makes      *
; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *
; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *
; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *
; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *
; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *
; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *
; *  States Government or the University of California, and shall not be used for advertising or product          *
; *  endorsement purposes.                                                                                        *
; *                                                                                                               *
; *****************************************************************************************************************
;


FUNCTION GKV_RestoreArray, Path=path, DirectoryName=DirectoryName, FamilyName=familyName, Debug
;
; This function returns an array of GKVsd objects which have previously
; been saved to disk using GKV_SaveArray.  
;
; Input Arguments
;
;	None	
;
; Outupt Arguments	
;
;	None
;
; Input Keywords
;
;	Path		Path to working directory in which the new 
;			directory, "DirectoryName", containing the
;			'.gkv' save files were created.
;			Defaults to current working directory.
;			(Optional)
;			
;	DirectoryName	Name of the directory (within the working
;			directory specified by 'Path') containing the
;			'.gkv' save files.  If not specified the
;			user will pick DirectoryName with "Dialog_Pickfile". 
;			(Optional)
;
;	FamilyName	The Save files should have a name of the form
;			'FamilyName'.'index'.gkv.  If not specified
;			GKV will extract 'FamilyName' from the first
;			file returned by "FINDFILE".
;			(Optional)
;
;	Debug		Set this keyword (i.e., put "/Debug" on the commandline)
;			to print out intermediate information which may be useful
;			in debugging this proceedure
;
;  Written by W.M. Nevins
;	4/12/01
;
separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
CD, current=currentWorkingDirectory
;
; Check for 'path'
;
IF(N_ELEMENTS(path) ) THEN BEGIN
	IF(TypeOf(Path) EQ 7) THEN CD, path
ENDIF
;
; Check for DirectoryName
;
IF( (N_Elements(DirectoryName) EQ 0) OR (TypeOf(DirectoryName) NE 7) ) THEN BEGIN
	DirName=DIALOG_PICKFILE(/DIRECTORY, FILTER='*_arr')
	IF(DirName EQ '') THEN BEGIN
		MESSAGE, "No directory was selected -- Returning", /INFORMATIONAL
		RETURN, 0
	ENDIF
ENDIF ELSE BEGIN
	CD, current=dirName
	len = STRLEN(dirName)
	lastChar = (STRMID(dirName, len-1, 1))[0]
	IF(lastChar NE separator) THEN dirName=dirName + separator
	dirName = dirName + DirectoryName[0]
ENDELSE
cd, dirName	
;
; Check for FamilyName
;
IF( (N_ELEMENTS(FamilyName) EQ 0) OR (TypeOf(FamilyName) NE 7) ) THEN BEGIN
	fileNames = FINDFILE(COUNT=nFiles)
	subStrings = STRSPLIT(fileNames[0], '.', /EXTRACT)
	FamilyName = subStrings[0]
ENDIF
IF(FamilyName EQ "") THEN BEGIN
	MESSAGE, "Couldn't form valid FamilyName -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
fileSearch = familyName +  '.*.gkv'
result = FINDFILE(fileSearch, COUNT=nObjs)
IF(nObjs EQ 0) THEN BEGIN
	MESSAGE, "Couldn't find .gkv files -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Create output array
;
output = OBJARR(nObjs)
;
; Get GKVsd objects to fill array
;
FOR i=0, nObjs-1 DO BEGIN
	index = STRING(i, FORMAT="(i4)")
	index = STRCOMPRESS(index, /REMOVE_ALL)
	fileName = familyName + '.' + index + '.gkv'
	result = FINDFILE(fileName, Count=thisCount)
	IF(thisCount EQ 1) THEN output[i] = GKV_RESTORE(FileName=fileName, Debug=debug)
ENDFOR	
;
; and we're done ...
;
CD, currentWorkingDirectory
RETURN, output
END ; ****** GKV_RestoreArray ****** ;



PRO GKV_SaveArray, ObjArr, Path=path, DirectoryName=DirectoryName, FamilyName=familyName, Debug
;
; This proceedure saves an array of GKVsd objects into to disk.
; It does this by first creating a new directory ("Directoryname_arr")
; and then writing separate .gkv file (see GKVsd::Save) with
; names FamilyName.index.gkv, where 'index' is the zero-based
; array index into ObjArr.
;
; Input Arguments
;
;	ObjArr		An array of GKVsd objects to be saved to disk
;
; Outupt Arguments	
;
;	None
;
; Input Keywords
;
;	Path		Path to working directory in which the new 
;			directory, "DirectoryName_arr", containing the
;			'.gkv' save files will be created.
;			Defaults to current working directory.
;			(Optional)
;			
;	DirectoryName	Name of the new directory to be created within
;			the directory specified by "Path".  Defaults
;			to the mnemonic of ObjArr[0] (or, "GKVsd_arr", if
;			ObjArr[0] has no mnemonic).  
;			(Optional)
; 
;	FamilyName	Family name of the sequence of GKV files.  Defaults
;			to the mnemonic of ObjArr[0] (or, "GKVsd", if
;			ObjArr[0] has no mnemonic).  
;			(Optional)
;
;	Debug		Set this keyword (i.e., put "/Debug" on the commandline)
;			to print out intermediate information which may be useful
;			in debugging this proceedure
;
;  Written by W.M. Nevins
;	4/12/01
;
info=SIZE(ObjArr)
nDims = info[0]

IF( info[nDims+1] NE 11 ) THEN BEGIN	; Argument is not an object reference
	MESSAGE, "Argument is not an object reference -- Returning", /INFORMATIONAL
	RETURN
ENDIF

nObjs = info[nDims+2]
IF( nObjs LT 1 ) THEN BEGIN		; ObjARR has no elements
	MESSAGE, "Object Array has no elements -- Returning", /INFORMATIONAL
ENDIF
;
; Check for "Path"
;
CD, current=currentWorkingDirectory
IF(N_ELEMENTS(Path) NE 0) THEN BEGIN
	IF(TypeOf(path) EQ 7)THEN cd, path
ENDIF
;
; Check for "DirectoryName"
;
ObjArr[0] -> get, mnemonic=dirName					; Use first objects mnemonic as default,
IF(dirName EQ "") THEN dirName = "GKV"					; or 'GKV' if no mnemonic present
dirName = dirName + '_arr'						; Append '_arr' to dirName
IF( N_ELEMENTS(DirectoryName) EQ 1) THEN BEGIN
	IF(TypeOF(DirectoryName) EQ 7) THEN BEGIN
		dirName = STRCOMPRESS(DirectoryName, /REMOVE_ALL)	; Remove all blank space
		subStrings = STRSPLIT(dirName, '_', /EXTRACT)		; Crack 'dirName' on underbars
		nSubStrings = N_ELEMENTS(subStrings)			; and force it to end in '_arr'
		IF(subStrings[nSubStrings-1] NE 'arr') THEN dirName = dirName + '_arr'
	ENDIF
ENDIF
;
; Check to be sure that the name 'DirName' is not already in use
;
fileNames = FINDFILE(COUNT=nFiles)
tryAgain : 
FOR i=0, nFiles-1 DO BEGIN
	thisFile = fileNames[i]
	IF(!D.NAME EQ 'MAC') THEN thisFile = (STRSPLIT(thisFile, ':', /EXTRACT))[0]
	IF(thisFile EQ DirName) THEN BEGIN				; Found a name conflict
		subStrings = STRSPLIT(dirName, '_', /extract)		; Check for embedded index number in 'dirName'
		nStrings = N_ELEMENTS(subStrings)
		ON_IOERROR, noIndex
			j = FIX(subStrings[nStrings-2])
		ON_IOERROR, NULL
		j = j + 1						; 'dirName' has an embedded index number,
		index = STRING(j, FORMAT="(i4)")				;  so we increment by one, and try again
		index = STRCOMPRESS(index, /REMOVE_ALL)
		subStrings[nStrings-2] = index
		dirName = STRJOIN(subStrings, "_")
		GOTO, tryAgain
	noIndex:							; 'dirName' does not have an embedded index number
		index = '1'						;  so we give it one (starting with one)
		index = STRCOMPRESS(index, /REMOVE_ALL)
		newStrings = STRARR(nStrings + 1)
		newStrings[0:nStrings-2] = subStrings[0:nStrings-2]
		newStrings[nStrings-1] = index
		newStrings[nStrings] = subStrings[nStrings-1]
		dirName = STRJOIN(newStrings, "_")
		GOTO, tryAgain
	ENDIF
ENDFOR
;
; Make new directory to hold .gkv files
;
FILE_MKDIR, dirName
cd, dirName
;
; Check for FamilyName
;
ObjArr[0] -> get, mnemonic=famName
IF(famName EQ "") THEN famName = "GKVsd"
IF( N_ELEMENTS(FamilyName) EQ 1) THEN BEGIN
	IF(TypeOF(Familyname)  EQ 7) THEN famName = STRCOMPRESS(FamilyName, /REMOVE_ALL)
ENDIF
;
; Save elements of ObjArr
;
FOR i=0, nObjs-1 DO BEGIN
	index = STRING(i, FORMAT="(i4)")
	index = STRCOMPRESS(index, /REMOVE_ALL)
	fileName = famName + '.' + index
	ObjArr[i] -> Save, FileName=fileName
ENDFOR
;
; Return to 'currentWorkingDirectory'
;
cd, currentWorkingDirectory
;
; and we're done ...
;
RETURN
END ; ****** GKV_SaveArray ****** ;
