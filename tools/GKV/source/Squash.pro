FUNCTION GKVs1D::Squash, argObj, Title=titleIn, Mnemonic=mnemonicIn
;
;  Purpose:
;
;	This function returns a GKVs1D object with the data of 'self'
;	as the dependent varialbe, and the data of 'argObj' as the 
;	independent variable.  Plots of this objects can be used to
;	check for relationships beteen data objects.
;
;  Arguments;
;
;	argObj		a GKVsd object of the same dimensionality as 'self' and 
;			having the same independent varialbles.
;
;  Input Keywords
;
;	Title		Title for output object.  Defaults to title of 'self'.
;			(Optional)
;
;	Mnemonic	Mnemonic for output object.  Defaults to mnemonic of 'self'.
;			(Optional)
;
;  Written by W.M. Nevins
;	4/20/01
;
nArgs = N_PARAMS()
IF(nArgs NE 1) THEN BEGIN
	MESSAGE, "Called with to many/few arguments -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
nDims = self -> NumDims()
CASE nDims OF
	1	:	class = "GKVS1D"
	2	:	class = "GKVS2D"
	3	:	class = "GKVS3D"
	4	:	class = "GKVS4D"
ENDCASE
IF(TypeOf(argObj) NE 11) THEN BEGIN
	MESSAGE, "Argument is not an object -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
IF NOT OBJ_ISA(argObj, class) THEN BEGIN
	MESSAGE, "Arg's class differs from self's class -- Returning", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Interpolate argObj onto save grid as self
;
temp = argObj -> INTERPOLATE(self)
;
; Coherce independent variable values into 1-D array
;
inDependentVar = *temp.values
inDepInfo = SIZE(inDependentVar)
inDepDims = inDepInfo[0]
nInDepVals= inDepInfo(inDepDims+2)
inDependentVar = REFORM(inDependentVar, nInDepVals, /OVERWRITE)
;
; Coherce dependent variable values into 1-D array
;
dependentVarPtr = self -> GetValues()
dependentVar = *dependentVarPtr
depInfo = SIZE(dependentVar)
depDims = depInfo[0]
nDepVals= depInfo(depDims+2)
dependentVar = REFORM(dependentVar, nDepVals, /OVERWRITE)
;
; prepare 'result' structure
;
result = {GKVs1D}
nDepTags = N_TAGS({GKVsd})
FOR i=0,nDepTags-1 DO BEGIN
	result.(i) = self.(i)
ENDFOR
IF(N_ELEMENTS(titleIn   ) EQ 1) THEN result.title=titleIn
argOBj -> Get, mnemonic=argMnemonic
result.mnemonic = result.mnemonic + "_vs_" + argMnemonic
IF(N_ELEMENTS(mnemonicIn) EQ 1) THEN result.mnemonic=mnemonicIn
indices = ["*"]
result.Indices = PTR_NEW(indices)
result.values = PTR_NEW(dependentVar)
result.errorBars = PTR_NEW()
;
; Prepare Grid structure
;
result.Grid1.Mnemonic	= temp.mnemonic
result.Grid1.title 	= temp.title
result.Grid1.units	= temp.units
result.Grid1.values	= PTR_NEW(inDependentVar)
result.Grid1.uniform	= 0B
result.Grid1.range	= temp.vrange
result.Grid1.irange	= [0, nInDepVals-1]
;
; Create GKVs1D object for output
;
output = OBJ_NEW("GKVs1D", result)
;
; Clean up
;
temp -> Trash
PTR_FREE, dependentVarPtr
;
; and we're done
;
RETURN, output
END ; ****** Squash ****** ;




