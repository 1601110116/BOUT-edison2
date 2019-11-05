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
FUNCTION GKV_StructureHasName, Structure
;
;Purpose:
;
;	This function returns 1 if input structure contains
;	the tag "name", and 0 otherwise
;
; Written by W.M. Nevins
;	4/13/01
;
IF(TypeOf(Structure) NE 8) THEN BEGIN
	MESSAGE, "WARNING: Argument is not a Structure -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
tagNames = TAG_NAMES(Structure)
nTags = N_ELEMENTS(tagNames)
FOR i=0, nTags-1 DO BEGIN
	thisTag = STRUPCASE(tagNames[i])
	IF(thisTag EQ 'NAME') THEN GOTO, GOTIT
ENDFOR
RETURN, 0
GOTIT	:
RETURN, 1
END ; ****** GKV_StructureHasName ****** ;


FUNCTION GKV_RestoreStructure, Path=path, DirectoryName=DirectoryName, Debug=debug
;
; This function returns a structure (possibly containing GKVsd objects or 
; arrays of GKVsd objects) which has previously been saved to disk using 
; GKV_SaveStructure.  
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
;			various save files and subdirectories was created.
;			Defaults to current working directory.
;			(Optional)
;			
;	DirectoryName	Name of the directory (within the working
;			directory specified by 'Path') containing the
;			various save files and subdirectories.  If not specified the
;			user will pick DirectoryName with "Dialog_Pickfile". 
;			(Optional)
;
;	Debug		Set this keyword (i.e., put "/Debug" on the commandline)
;			to print out intermediate information which may be useful
;			in debugging this proceedure
;
;  Written by W.M. Nevins
;	4/13/01
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
	DirName=DIALOG_PICKFILE(/DIRECTORY, FILTER='*_str')
	IF(DirName EQ '') THEN BEGIN
		MESSAGE, "No directory was selected -- Returning", /INFORMATIONAL
		RETURN, 0
	ENDIF 
ENDIF ELSE BEGIN
	CD, current=dirName
	len = STRLEN(dirName)
	STRPUT, lastChar, dirName, len-1
	IF(lastChar NE separator) THEN dirName=dirName + separator
	dirName = dirName + DirectoryName
ENDELSE
cd, dirName	
;
; Check for 'StrInfo.dat'
;
strInfo = FINDFILE('StrInfo.dat', COUNT=nDats)
IF(nDats EQ 0) THEN BEGIN
	MESSAGE, "Couldn't find 'StrInfo.dat' ", /INFORMATIONAL
	output = {name:"GKVsd"}
ENDIF ELSE BEGIN
	RESTORE, strInfo[0], VERBOSE=debug	; This call 'restores' the output structure from previouc
ENDELSE						;  call to GKV_SaveStructure
;
; Check for .gkv files
;
fileNames = FINDFILE('*.gkv', COUNT=nFiles)
IF(nFIles GT 0) THEN BEGIN
	FOR i=0, nFiles-1 DO BEGIN
		thisObj = GKV_RESTORE(FileName=fileNames[i])
		subStrings = STRSPLIT(fileNames[i], '.', /EXTRACT)
		nStrings = N_ELEMENTS(subStrings)
		thisTag = STRJOIN(subStrings[0:nStrings-2], '.')
		output = CREATE_STRUCT(output, thistag, thisObj)
	ENDFOR
ENDIF
;
; Check of '_arr' directories
;
dirNames = FINDFILE("*_arr", COUNT=nDirs)
IF(nDirs GT 0) THEN BEGIN
	FOR i=0, nDirs-1 DO BEGIN
		thisDir = (STRSPLIT(dirNames[i], separator, /extract))[0]
		thisArray = GKV_RestoreArray(DirectoryName=thisDir)
		subStrings = STRSPLIT(dirNames[i], '_', /EXTRACT)
		nStrings = N_ELEMENTS(subStrings)
		thisTag = STRJOIN(subStrings[0:nStrings-2], '_')
		output = CREATE_STRUCT(output, thistag, thisArray)		
	ENDFOR
ENDIF
;
; and we're done ...
;
CD, currentWorkingDirectory
RETURN, output
END ; ****** GKV_RestoreStructure ****** ;



PRO GKV_SaveStructure, Structure, Path=path, DirectoryName=DirectoryName, Debug=debug
;
; This proceedure save a structure (possibly containing GKVsd objects) 
; onto to disk.  It does this by first creating a new directory 
; ("Directoryname_str") and then writing separate SAVE files using
; GKVsd::Save, GKV_SaveArray, and the native IDL SAVE routine as 
; appropriate
;
; Input Arguments
;
;	Structure	A structure which may contain GKVsd objects
;			and/or GKVsd object arrays
;
; Outupt Arguments	
;
;	None
;
; Input Keywords
;
;	Path		Path to working directory in which the new 
;			directory, "DirectoryName_str", containing the
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
;	Debug		Set this keyword (i.e., put "/Debug" on the commandline)
;			to print out intermediate information which may be useful
;			in debugging this proceedure
;
;  Written by W.M. Nevins
;	4/13/01
;

IF( TypeOf(Structure) NE 8 ) THEN BEGIN	; Argument is not a struture
	MESSAGE, "Argument is not a structure -- Returning", /INFORMATIONAL
	RETURN
ENDIF

nTags = N_TAGS(Structure)
IF( nTags LT 1 ) THEN BEGIN		; Structure has no elements
	MESSAGE, "Structure has no tags -- Returning", /INFORMATIONAL
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
DirName = ''
IF(GKV_StructureHasName(Structure)) THEN BEGIN				; If structure has a 'name' tag
	name=structure.name						;  use its value as default directory name
	IF(TypeOf(name) EQ 7) THEN DirName=name
ENDIF
IF(dirName EQ "") THEN dirName = "GKVsd"				;  or 'GKVsd' if no name Tag present
IF(TypeOf(name) NE 7) THEN name = dirName
dirName = dirName + '_str'						; Append '_str' to dirName
IF( N_ELEMENTS(DirectoryName) EQ 1) THEN BEGIN
	IF(TypeOF(DirectoryName) EQ 7) THEN BEGIN
		dirName = STRCOMPRESS(DirectoryName, /REMOVE_ALL)	; Remove all blank space
		subStrings = STRSPLIT(dirName, '_', /EXTRACT)		; Crack 'dirName' on underbars
		nSubStrings = N_ELEMENTS(subStrings)			; and force it to end in '_arr'
		IF(subStrings[nSubStrings-1] NE 'str') THEN BEGIN
			name = dirName
			dirName = dirName + '_str'
		ENDIF ELSE BEGIN
			name = STRJOIN(substrings[0:nSubStrings-2], '_')
		ENDELSE
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
		index = STRING(j, FORMAT="(i4)")			;  so we increment by one, and try again
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
; Make new directory to hold various Save files and directories
;
FILE_MKDIR, dirName
cd, dirName
;
; Create structure to hold tags which are NOT objects or object arrays
;
output = CREATE_STRUCT("Name", name)
;
; Save elements of Structure
;
tagNames = TAG_NAMES(Structure)
FOR i=0, nTags-1 DO BEGIN
	thisName = tagNames[i]
	IF(thisName NE 'NAME') THEN BEGIN
		thisValue = Structure.(i)
		IF(TypeOF(thisValue) EQ 11) THEN BEGIN
			CASE N_ELEMENTS(thisValue) OF
				0	:	
				1	:	thisValue[0] -> Save, fileName=thisName
				else	:	GKV_SaveArray, thisValue, DirectoryName=thisName
			ENDCASE
		ENDIF ELSE BEGIN
			output = CREATE_STRUCT(output, thisName, thisValue)
		ENDELSE
	ENDIF 
ENDFOR
SAVE, output, FILENAME='StrInfo.dat', /COMPRESS, VERBOSE=debug
;
; Return to 'currentWorkingDirectory'
;
cd, currentWorkingDirectory
;
; and we're done ...
;
RETURN
END ; ****** GKV_SaveStructure ****** ;
