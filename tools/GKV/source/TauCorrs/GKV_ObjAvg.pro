FUNCTION GKVs1D::ObjAvg, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, $
			 title=utitle, mnemonic=umnemonic, units=uunits, vrange=uvrange, keepSelf=keepSelf
;
; Purpose:
;
;	This function computes the average and standard deviation of
;	the objects passed as argument(s).  
;	It returns an object of the same dimensionality as 'self'
;	whose value is the average of input objects, and whose
;	error bars reflect BOTH the standard deviation, and any
;	existing error bars..  Normally,
;	'self' is NOT included in computing the average and standard
;	deviations.  However, if the keyword "KeepSelf" is set, then
;	'self' be included in computing the average and standard deviation.
;
; Written by W.M. Nevins
;	3/20/01
;
; Check if have arguments, and that the first argument is an object, etc.
;

avg = self -> MakeCopy()
avg -> Restrict
IF( PTR_VALID(avg.ErrorBars) ) THEN BEGIN
	PTR_FREE, avg.ErrorBars
	avg.ErrorBars = PTR_NEW()
ENDIF

nArgs = N_PARAMS()

title = "!12<!X" + avg.title + "!12>!X"
IF(N_ELEMENTS(utitle) EQ 1) 	THEN title=utitle
units = avg.units
IF(N_ELEMENTS(uunits) EQ 1) 	THEN units=uunits
mnemonic = "Avg_" + avg.mnemonic
IF(N_ELEMENTS(umnemonic) EQ 1)	THEN mnemonic=umnemonic
vrange = avg.vrange
IF(N_ELEMENTS(uvrange) EQ 2) 	THEN vrange=uvrange
avg -> Set, title=title, mnemonic=mnemonic, units=units, vrange=vrange


IF (nArgs LT 1) THEN RETURN, avg
;
; Get class of GKVsd object
;
class = OBJ_CLASS(self)
;
; Concatenate arguments into single object array
;
args = STRARR(11)
args = ["arg0", "arg1", "arg2", "arg3", "arg4", "arg5", "arg6", "arg7", "arg8", "arg9", "arg10"]

nElements   = INTARR(nArgs)
totElements = INTARR(nArgs+1)
FOR i=0,nArgs-1 DO BEGIN				; count number of valid GKVsXD objects passed as arguments
	commandString = "nElements[i] = TOTAL( OBJ_ISA(" + args[i] + ", class) )"
	ok = EXECUTE(commandString)
ENDFOR
totElements[1:nArgs] = TOTAL(nElements, /CUMULATIVE)
nObjs = TOTAL(nElements)
obj_Array = OBJARR(nObjs)

FOR i=0,nArgs-1 DO BEGIN				; Load valid GKVsXD objects into 'obj_Array'
	IF(nElements[i] NE 0) THEN BEGIN
		commandString = "addresses = WHERE( OBJ_ISA(" + args[i] + ", class) )"
		ok = EXECUTE(commandString)
		commandString = "obj_Array[totElements[i]:(totElements[i+1]-1)] = " + args[i] + "[addresses]"
		ok = EXECUTE(commandString)
	ENDIF
ENDFOR

nErrors = 0
avgValuesPtr = avg -> GetValues()
avgValues = *avgValuesPtr
PTR_FREE, avgValuesPtr
sqErrors = 0.*avgValues
IF( NOT KEYWORD_SET(keepSelf) ) THEN BEGIN	; If 'keepSelf' is not set, then do NOT
	avgValues = 0.*avgValues		; include 'self' when computing average
	IF( PTR_VALID(avg.ErrorBars) ) THEN BEGIN
		errorPtr = avg -> GetErrors()
		sqErrors = (*errorPtr)^2
		nErrors = 1
		PTR_FREE, errorPtr
	ENDIF
ENDIF
avgValuesSq =avgValues^2

nnobjs=0
FOR i=0,nObjs-1 DO BEGIN				; Compute Sum of objects, sum of objects^2
	obj_array[i] -> get, vrange=vrange
	ok = TOTAL(FINITE(vrange))			; test for bad values
	IF(ok NE 2) THEN GOTO, skip
	nnobjs=nnobjs+1
	objValuesPtr = obj_Array[i] -> GetValues()
	objValues = *objValuesPtr
	PTR_FREE, objValuesPtr
	avgValues = avgValues + objValues
	avgValuesSq = avgValuesSq + objValues^2
	IF( PTR_VALID((obj_Array[i]).ErrorBars) ) THEN BEGIN
		errorPtr = obj_Array[i] -> GetErrors()
		objErrors = *errorPtr	 
		sqErrors = sqErrors + objErrors^2
		nErrors = nErrors + 1
		PTR_FREE, errorPtr
	ENDIF
	skip: 	
ENDFOR

;
; Check for incrementing (i.e., changing) element of indices array
;
oldIndices = *self.indices
nIndices = N_ELEMENTS(oldIndices)
indicesChanged = INTARR(nIndices, nObjs)
FOR i=0, nObjs-1 DO BEGIN
	tempIndices = *((obj_Array[i]).indices)
	indicesChanged[*,i] = (oldIndices NE tempIndices[0:(nindices-1)])
ENDFOR
nChanges = TOTAL(indicesChanged, 2)
j=0
newIndices = STRARR(nIndices)
symbols = STRARR(nIndices)
nSymbols = 0
FOR i=0, nIndices-1 DO BEGIN
	IF(nChanges[i] EQ 0) THEN BEGIN
		newIndices[j] = oldIndices[i]
		j=j+1
	ENDIF ELSE BEGIN
		subStrings = STRSPLIT(oldIndices[i], "=", /EXTRACT)
		symbols[nSymbols] = STRCOMPRESS(subStrings[0], /REMOVE_ALL)
		nSymbols=nSymbols + 1
	ENDELSE
ENDFOR
PTR_FREE, avg.indices
avg -> set, indices = PTR_NEW(newIndices[0:(j-1)])		; Update pointer to 'indices' array
IF( N_ELEMENTS(utitle) EQ 0) THEN BEGIN				; add 'avgOver' as substript on title
	IF(nSymbols GT 0) THEN BEGIN
		avgOver = STRJOIN(symbols[0:nSymbols-1], ",")
		avgOver = STRCOMPRESS(avgOver, /REMOVE_ALL)
		title = title + "!D" + avgOver + "!N!X"
		avg -> set, title=title
	ENDIF
ENDIF


IF KEYWORD_SET(keepSelf) THEN nObjs = nObjs + 1

avgValues = avgValues/nnObjs				; Compute average
avg -> Set, values=PTR_NEW(avgValues)

avgValuesSQ = avgValuesSq/nnObjs			; Compute standard deviation
std = avgValuesSQ - avgValues^2
IF(nErrors GT 0) THEN std = std + sqErrors/nErrors	; add contribution from pre-existing error bars (if any)
std = SQRT(std)
IF PTR_VALID(avg.ErrorBars)  THEN PTR_FREE, avg.ErrorBars
avg -> Set, ErrorBars = PTR_NEW(std)


RETURN, avg

END ;****** GKV_ObjAvg ****** ;





