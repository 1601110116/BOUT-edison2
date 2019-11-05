FUNCTION GKV_ChangeDir, Path=newPath, New=newDirName, Current=current
;
;  Purpose:
;
;	This function allows the user to select a directory
;	(for example, a new working directory) using the 
;	Dialog_PickFile widget, and then makes this directory
;	the current working directory.
;
; Output:
;
;	The function returns a legel 'path' to the directory
;	selected by the user.
;
;
; Input KeyWords:
;
;	Path	Set this keyword to the path where you wish to begin
;		your search (defaults to the current working directory).
;
;	New	If this keyword is invoked, then this function will create
;		a new directory within the directory selected by the user.
;
;	Current	Set this keyword to prevent GKV_ChangeDir from actually changing
;		the current working directory to that selected by the user.  
;		However, the path to the selected directory is still returned.
;
;  Written by W.M. Nevins
;	4/22/01
;
; Set path
;
CD, current=current_working_directory
path=current_working_directory
IF(TypeOf(newPath) EQ 7) THEN path=newPath
;
; Select directory
;
outPut = DIALOG_PICKFILE(/DIRECTORY, PATH=path)
IF(output EQ "") THEN BEGIN
	MESSAGE, "No Directory selected -- Returning null string"
	RETURN, output
ENDIF
;
; Change to selected directory
;
CD, output
;
; Create new directory if requested by "New" keyword
;
IF(TypeOF(newDirName) EQ 7) THEN BEGIN
	;
	; Check to be sure that the name 'DirName' is not already in use
	;
	fileNames = FINDFILE(COUNT=nFiles)
tryAgain : 
	FOR i=0, nFiles-1 DO BEGIN
		thisFile = fileNames[i]
		IF(!D.NAME EQ 'MAC') THEN thisFile = (STRSPLIT(thisFile, ':', /EXTRACT))[0]
		IF(thisFile EQ newDirName) THEN BEGIN				; Found a name conflict
			subStrings = STRSPLIT(newDirName, '_', /extract)	; Check for embedded index number in 'newDirName'
			nStrings = N_ELEMENTS(subStrings)
			ON_IOERROR, noIndex
				j = FIX(subStrings[nStrings-2])
			ON_IOERROR, NULL
			j = j + 1						; 'newDirName' has an embedded index number,
			index = STRING(j, FORMAT="(i4)")			;  so we increment by one, and try again
			index = STRCOMPRESS(index, /REMOVE_ALL)
			subStrings[nStrings-2] = index
			newDirName = STRJOIN(subStrings, "_")
		GOTO, tryAgain
	noIndex:								; 'dirName' does not have an embedded index number
			MESSAGE, "Name Conflict, adding index number", /INFORMATIONAL
			index = '1'						;  so we give it one (starting with one)
			index = STRCOMPRESS(index, /REMOVE_ALL)
			newStrings = STRARR(nStrings + 1)
			newStrings[0:nStrings-2] = subStrings[0:nStrings-2]
			newStrings[nStrings-1] = index
			newStrings[nStrings] = subStrings[nStrings-1]
			newDirName = STRJOIN(newStrings, "_")
		GOTO, tryAgain
		ENDIF
	ENDFOR
	;
	; Make new directory to hold .gkv files
	;
	FILE_MKDIR, newDirName
	CD, newDirName
	CD, current=output
	PRINT, "Creating Directory ", newDirName
ENDIF
;
; Check for "Current" keyword
;
IF KEYWORD_SET(Current) THEN CD, current_working_directory
;
; and we're done
;
RETURN, output
END ; ****** GKV_ChangeDir ****** ;