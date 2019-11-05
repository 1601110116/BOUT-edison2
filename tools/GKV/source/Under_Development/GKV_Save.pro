Pro GKVsd::Save, GKVsdOjb, FileName=file_name, DeBug=d
;
; Proceedure to SAVE a GKVsd object to disk.
; 
filename=self.mnemonic + '.GKVsd'
IF(KEYWORD_SET(file_name)) THEN BEGIN
	subStrings=STRSPLIT(file_name, '.', /EXTRACT)
	nSubStrings = N_ELEMENTS(subStrings)
	extension = subStrings[nSubStrings-1]
	filename=file_name
	IF(extension NE '.GKVsd') THEN filename=file_name + '.GKVsd'
ENDIF
GKVsdClass = OBJ_CLASS(self)
;

SAVE, GKVsdClass, self, FIleName=filename, /COMPRESS, VERBOSE=d
RETURN
END ; ****** GKVsd::Save ****** ;

FUNCTION GKVsd_RESTORE, FileName=filename, DeBug=D
;
; Proceedure to Restore the GKVsd object previously saved 
; (using GKVsd::Save) in the file 'file_name'
;
IF(NOT KEYWORD_SET(filename)) THEN filename = DIALOG_PICKFILE(Filter='*.GKVsd')
;
RESTORE, filename, /RELAXED_STRUCTURE_ASSIGNMENT, VERBOSE=d
result = OBJ_NEW(GKVsdClass, self)
RETURN, result
END ; ****** GKVsd_RESTORE ****** ;
