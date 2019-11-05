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

FUNCTION Query_Structure, arg
;
; Returns 1 if ARG is a STRUCTURE, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 8 ) then return, 1
return, 0
END


FUNCTION Query_Integer, ARG
;
; Returns 1 if ARG is an integer, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 0 ) then return, 0
if (arg_type lt 4 ) then return, 1
if (arg_type lt 12) then return, 0
if (arg_type lt 15) then return, 1
return, 0
END


Function Query_Real, ARG
;
; Returns 1 if ARG is a real number, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 4 ) then return, 1
if (arg_type eq 5 ) then return, 1
return, 0
END


Function Query_Complex, ARG
;
; Returns 1 if ARG is a Complex number, 0 otherwise
; (see help document for SIZE function for definition of TYPE)
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 6 ) then return, 1
if (arg_type eq 9 ) then return, 1
return, 0
END


Function Query_String, ARG
;
; Returns 1 if ARG is a string, 0 otherwise
;
arg_type=SIZE(ARG, /TYPE)
if (arg_type eq 7 ) then return, 1
return, 0
END


FUNCTION Subscript, string
;
; accepts a string, which is assumed to be Hershey (Vector) font.
; searches for line-shifting characer !N, and replaces it by !L
; 
; First, make sure that all subscripts are obtained with !I (and NOT !L)
;
result=STRSPLIT(string, '!L', /EXTRACT, /REGEX)
result= '!L' + STRJOIN(result, '!I')
;
; Now, replace any !N command with !L
;
result=STRSPLIT(string, '!N', /EXTRACT, /REGEX)
result=STRJOIN(result, '!L')
;
; Finally, append a !L to the front
;
result = '!L' + result  
RETURN, result
END ; ****** Subscript ****** ;


FUNCTION Superscript, string
;
; accepts a string, which is assumed to be Hershey (Vector) font.
; searches for line-shifting characer !N, and replaces it by !U
; 
; 
; First, make sure that all superscripts are obtained with !E (and NOT !U)
;
result=STRSPLIT(string, '!U', /EXTRACT, /REGEX)
result= '!L' + STRJOIN(result, '!E')
;
; Now, replace any !N command with !U
;
result=STRSPLIT(string, '!N', /EXTRACT, /REGEX)
result=STRJOIN(result, '!U')
;
; Finally, append a !U to the front
;
result = '!U' + result  
RETURN, result
END ; ****** Superscript ****** ;


FUNCTION StringClean, string
;
; replaces any !L's in 'string' with !I's; and
; replaces any !U's in 'string' with !E's.
;
result=STRSPLIT(string, '!L', /EXTRACT, /REGEX)
result= STRJOIN(result, '!I')
result=STRSPLIT(result, '!U', /EXTRACT, /REGEX)
result= STRJOIN(result, '!E')
RETURN, result
END ; ****** StringClean ****** ;

FUNCTION GetKeyWord, keyword, structure,Debug=d
;
; Searches 'structure' for the tag 'keyword'.
; Returns value field associated with this tag.
;
; NOTE that 'structure' may be undefined on entry.  In this case, it is
; important that 'structure' only appear as an argument to N_ELEMENTS, TypeOf, or SIZE.
; A null result ('undefined') should be returned.
;
; On return, 'structure' contains all tags and values
; from 'structure' EXCEPT for 'keyword' and its associated value
; if they appeared in 'structure' on entry.
; However, if 'keyword' is the only tag in 'structure' then
; 'structure' is set to -1 on return (IDL 5.3 does not support null structures).
;
; Written by W.M. Nevins
;	3/5/00
;
debug=1
result='undefined'
otherTags=0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
IF(N_ELEMENTS(structure) EQ 0) THEN RETURN, result
arg1Type = TypeOf(keyword)
IF(arg1Type NE 7)  THEN BEGIN
	MESSAGE, 'First argument must be a string', Informational=d
	return, result
ENDIF
arg2Type = typeOF(structure)							; Check for proper argument type
IF(arg2Type NE 8) THEN BEGIN							; Wasn't passed a structure, so
	structure = -1								;	set Nstructure=-1 (as a flag)
	RETURN, result								;	and return 'undefined'
ENDIF
nTags = N_TAGS(structure)
IF(nTags EQ 0) THEN RETURN, result
tagNames = TAG_NAMES(structure)
tagNames = STRTRIM(tagNames, 2)							; Remove both leading and trailing blanks
command_str = 'structure = {'
tagToFind = STRTRIM(keyword)							; Remove both leading and trailing blanks
FOR i=0, ntags-1 DO BEGIN							; Search tags of 'structure'
	IF( STRCMP(tagToFind, tagNames[i], /FOLD_CASE) ) THEN BEGIN
		result = structure.(i)						; IF mulitple occurances of 'keyword' in
	ENDIF ELSE BEGIN							;	'structure', only the LAST occurance is significant
		i_str = STRING(i, FORMAT='(I3)')				; Construct command string for 'Nstructure'
		i_str = STRTRIM(i_str, 2)						; Remove both leading and trailing blanks
		IF(otherTags NE 0) THEN command_str = command_str + ', '
		command_str = command_str + tagNames[i] + ':' + 'structure.(' + i_str + ')'
		othertags = othertags + 1
	ENDELSE
ENDFOR
command_str = command_str + '}'
IF(othertags NE 0) THEN BEGIN
	ok = EXECUTE(command_str)
ENDIF ELSE BEGIN
	structure = -1
ENDELSE
	RETURN, result
END ; ******GetKeyWord ****** ;


FUNCTION GKV_CharSize, Old_Sizes = oldSizes, Spacing = spczng
;
; Purpose:
;
;		This function computes an appropriate size for
;		the Hershey Vector font characters within the
;		currently open graphics window.
;
; Output:
;
;		GKV_CharSize returns a two element integer array
;		containing the recommend charactersitic vertical 
;		height (!D.Y_CH_SIZE) and the line spacing.  
;
; Arguments:
;
;		None.
;
;
; Input Keywords:
;
;		None.
;
;
; Output Keywords:
;
;		Old_Char_Size	On output this keyword is set to a two element
;				integer array containing the initial value of 
;				!D.Y_CH_SIZE and an estimate of the initial value
;				of the linespacing. 
;
; Written by W.M. Nevins
;	7/19/00
;

Old_X_CH_SIZE = !D.X_CH_SIZE
Old_Y_CH_SIZE = !D.Y_CH_SIZE
char_Aspect_Ratio = FLOAT(Old_Y_CH_SIZE)/FLOAT(Old_X_CH_SIZE)
oldSizes = [Old_X_CH_SIZE, Old_Y_CH_SIZE]
X_CH_SIZE = INTARR(2)
;
; Choose an appropriate character size
; scaled to width of the currently open graphics window
;
X_CH_SIZE[0]  = !D.X_VSIZE/60.
Y_CH_SIZE_0 = LONG(char_Aspect_Ratio*X_CH_SIZE[0])
;
; Choose an appropriate character size
; scaled to height of the currently open graphics window
;
Y_CH_SIZE_1 = FLOAT(0.03*!D.Y_VSIZE)
X_CH_SIZE[1] = LONG(Y_CH_SIZE_1/char_Aspect_Ratio)
;
; Choose smallest of these alternatives
;
new_X_CH_SIZE = MIN(X_CH_SIZE)

New_Y_CH_SIZE = LONG(char_Aspect_Ratio*new_X_CH_SIZE)
newSizes = [New_X_CH_SIZE, New_Y_CH_SIZE]
RETURN, newSizes
END	; ****** GKV_CharSize ****** ;


FUNCTION GKV_TitleSize, title, Width=xwidth, yPosition=ypos
;
; Purpose:
;
;		This function computes an appropriate value for
;		the keyword CHARSIZE in XYOUTS for writting titles
;		to GKVsd plots.
;
; Argument:
; (optional)
;
;		The proposed title is accepted as a string argument.
;		When present, the recommend CHARSIZE is chosen such
;		that the title will fit within the available width
;		as well as within the available height.
;		
;  Input Keywords:
;
;			Width		Maximum width available for 'title'
;					in NORMALIZED coordinates.  Defaults to 0.75.
;					(Optional)
;
; Output Keywords:	
;
;			yPosition	Recommended value for y position in
;					DEVICE coordinates for writting the
;					title to a GKVsd plot with XYOUTS
;
; Written by W.M. Nevins
;	7/19/00
;
height = 0.1
width = 0.75
IF(N_ELEMENTS(xwidth) NE 0) THEN width = xwidth
CharSize = FLTARR(2)
;
; Estimate largest Y_CH_SIZE from available height
;
CharSize[0] = (height/3.)*!D.Y_SIZE/!D.Y_CH_SIZE
;
; Estimate largest Y_CH_SIZE from available width
;
CharSize[1] = STR_SIZE(title, width, INITSIZE=CharSize[0])
;
; Choose smaller of these two values
;
result = MIN(CharSize)
;
; Compute recommend yPosition
;
yPos = !D.Y_SIZE - 0.5*(height*FLOAT(!D.Y_SIZE) + result*FLOAT(!D.Y_CH_SIZE))

RETURN, result

END	; ****** GKV_TitleSize ****** ;	


FUNCTION GKV_RemoveParens, stringIn
;
;  This function accepts an ascii string
;  as input, and returns a string in which
;  any leading or trailing parenthesis are 
;  removed.
;  
;  Written by W.M. Nevins
;	9/24/02
;
; First, remove leading/trailing '('
;
subStrings = STRSPLIT(stringIn, '(', /EXTRACT)
temp = STRJOIN(subStrings, '(')
;
; Now, remove trailing/leading ')'
;
subStrings = STRSPLIT(temp, ')', /EXTRACT)
result = STRJOIN(subStrings, ')')
RETURN, result
END  ;  ****** GKV_RemoveParens ******  ;


FUNCTION GKV_InvertString, stringIn
;
; This function assumes that the input string is an
; ascii representation of a fraction, and returns
; an ascii representation of the inverse of that
; fraction.
;
; Written by W.M. Nevins
;	9/24/02
;
; Break string into sub-strings at '/'
;
subStrings1 = STRSPLIT(stringIn, '/', /EXTRACT)
nSubStrings = N_ELEMENTS(subStrings1)
IF(nSubStrings EQ 0) THEN RETURN, ''
IF(nSubStrings GT 2) THEN BEGIN
	MESSAGE, "Too Many /'s", /INFORMATIONAL
	result = '1/(' + GKV_RemoveParens(stringIn) + ')'
	RETURN, result
ENDIF
;
; Remove any leading/trailing parenthesis
;
subStrings = STRARR(nSubStrings)
FOR i=0, nSubStrings-1 DO subStrings[i] = GKV_RemoveParens(subStrings1[i])

CASE nSubStrings OF
	1	:	result = '1/(' + subStrings[0] + ')'
	2	: 	result = '(' + subStrings[1] + ')/(' + subStrings[0] + ')'
ENDCASE
RETURN, result
END  ;  ****** GKV_InvertString  ******  ;
	

PRO GKVdelete, arg, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
;
; purpose of this proceedure is to delete 'arg',
; being sure that pointers are released, objects 
; properly cleaned up, etc.
;
; Written by W.M. Nevins
;	6/7/00
;
nArgs = N_PARAMS()
argType=TypeOf(arg)
;
; Delete first argument
;
CASE argType OF 
	0	:	RETURN
	1	:	BEGIN
				arg=0
			END
	2	:	BEGIN
				arg=0
			END
	3	:	BEGIN
				arg=0
			END
	4	:	BEGIN
				arg=0
			END
	5	:	BEGIN
				arg=0
			END
	6	:	BEGIN
				arg=0
			END
	7	:	BEGIN
				arg=0
			END
	8	:	BEGIN		; arg is a structure
				nTags = N_TAGS(arg)
				FOR i=0, nTags-1 DO GKVdelete, arg.(i)
				arg=0
			END
	9	:	BEGIN
				arg=0
			END
	10	:	BEGIN		; arg is a pointer
				PTR_FREE, arg
			END
	11	:	BEGIN		; arg is an object reference
				IF(OBJ_ISA(arg[0], 'GKVsd')) THEN BEGIN
					n = N_ELEMENTS(arg)
					IF(N EQ 1) THEN BEGIN 
						arg -> Trash
					ENDIF ELSE BEGIN
						FOR i=0,n-1 DO BEGIN
						  IF(OBJ_ISA(arg[i], 'GKVsd')) THEN arg[i] -> trash
						ENDFOR
					ENDELSE
					arg=0
				ENDIF ELSE BEGIN
					OBJ_DESTROY, arg
					arg=0
				ENDELSE
			END
	12	:	BEGIN
				arg=0
			END
	13	:	BEGIN
				arg=0
			END
	14	:	BEGIN
				arg=0
			END
	15	:	BEGIN
				arg=0
			END
ENDCASE
;
; If this gkvDelete ws called with more than one argument,
; use this routine recursively to delete additional arguments
CASE nArgs OF
	1	:	RETURN
	2	:	GKVdelete, arg0
	3	:	GKVdelete, arg0, arg1
	4	:	GKVdelete, arg0, arg1, arg2
	5	:	GKVdelete, arg0, arg1, arg2, arg3
	6	:	GKVdelete, arg0, arg1, arg2, arg3, arg4
	7	:	GKVdelete, arg0, arg1, arg2, arg3, arg4, arg5
	8	:	GKVdelete, arg0, arg1, arg2, arg3, arg4, arg5, arg6
	9	:	GKVdelete, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7
	10	:	GKVdelete, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8
	11	:	GKVdelete, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
	ELSE	: 	MESSAGE, 'Too many arguments... returning', /INFORMATIONAL
ENDCASE

RETURN
END ; ****** GKVdelete ****** ;


FUNCTION GKV_CopyStructure, arg
;
; This function returns a 'deep' copy of its argument, which
; is assumed to be a structure
;
; Written by W.M. Nevins
;  9/30/03
;
IF(TypeOf(arg) NE 8) THEN BEGIN
	MESSAGE, "Argument is not a structure, returning", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Create output structure
;
Output = arg
;
; Get some basic info
;
nTags = N_TAGS(arg)
tagNames = TAG_NAMES(arg)
;
; process tags, putting copies of all GKVsd objects in the OUTPUT structure
;
FOR i=0, nTags-1 DO BEGIN
	value = arg.(i)
	IF( TypeOf(value) EQ 11 ) THEN BEGIN
		n = N_ELEMENTS(value)
		IF(n EQ 1) THEN BEGIN
			value = value -> MakeCopy()
		ENDIF ELSE BEGIN
			array=value
			value = OBJARR(n)
			FOR j=0,n-1 DO value[j] = array[j] -> MakeCopy()
		ENDELSE
	Output.(i) = value
	ENDIF
ENDFOR
RETURN, Output
END  ; ****** FUNCTION GKV_CopyStructure ****** ;
				
