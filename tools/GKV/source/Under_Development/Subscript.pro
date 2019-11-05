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
result= '!L' + STRJOIN(result, '!I')
result=STRSPLIT(string, '!U', /EXTRACT, /REGEX)
result= '!L' + STRJOIN(result, '!E')
RETURN, result
END ; ****** StringClean ****** ;