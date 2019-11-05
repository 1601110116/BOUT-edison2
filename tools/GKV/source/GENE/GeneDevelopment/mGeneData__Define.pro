FUNCTION mGeneData::GetData, _Extra=Extra
;
; returns GKV objects using
; data in multiple GeneData objects.
; Keywords are same as for
; GeneData::GetData 
; (and most are simply passed along)
;
; Written by W.M. Nevins
;	2/13/2009
;
Intervals = *self.Intervals
GeneDataObjs = *self.GeneDataObjs
nObjs = N_ELEMENTS(GeneDataObjs)
sampleGKVObj = GeneDataObjs[0] -> GKV_Obj()
tAxis = sampleGKVObj -> NumDims()
sampleGKVObj -> GET, axis=tAxis, GridMnemonic=tMnemonic
;
; parse command line for time interval
;
nExtra = Extra
result = GetKeyWord(tMnemonic, nExtra)
IF(Query_Real(result) + Query_Integer(result)) THEN t=result
;
; If time mnemonic is not on command line, 
; assume user wants entire time interval
; covered by data set
;
IF(N_ELEMENTS(t) EQ 0) THEN BEGIN
	tMin = Intervals[0,0]
	tMax = Intervals[nObjs-1, 1]
	t = [tMin, tMax]
ENDIF
;
; If t[1] is greater than maximum time, set
; to tMax
;
tMax = Intervals[nObjs-1,1]
t[1] = t[1] < tMax
;
; Find time interval containing t[0]
;
lower = t[0] GE Intervals[*,0]
upper = t[0] LE Intervals[*,1]
ok = upper AND lower
lowerIndex = WHERE(ok)
;
IF(N_ELEMENTS(t) EQ 1) THEN BEGIN
	index = lowerIndex[0]
	result = GeneDataObjs[index] -> GetData(_Extra=Extra)
	RETURN, result
ENDIF
;
; Find time interval containing t[1]
;
lower = t[1] GE Intervals[*,0]
upper = t[1] LE Intervals[*,1]
ok = upper AND lower
upperIndex = WHERE(ok)
;
; data overlaps, so endpoints may lie in
; more than one GeneData object. Let's
; use as few GeneDataObjects as possible!
;
lower = MAX(lowerIndex)
upper = MIN(upperIndex)
nnObjs = upper - lower + 1
IF(nnObjs EQ 1) THEN BEGIN
	result = GeneDataObjs[lower] -> GetData(_Extra=Extra)
	RETURN, result
ENDIF
;
; time interval covers multiple GeneDataObjs
; (that is, multiple field files). Form
; array of GKV objects which cover this interval
;
gkvArr = OBJARR(nnObjs)
mExtra = Extra
gkvArr[0] = GeneDataObjs[lower] -> GetData(_Extra=mExtra)
IF(nnObjs GT 2) THEN BEGIN
FOR i = 1, nnObjs-2 DO BEGIN
	nnExtra = nExtra
	gkvArr[i] = GeneDataObjs[lower+i] -> GetData(_Extra=nnExtra) 
ENDFOR
ENDIF
gkvArr[nnObjs-1] = GENEDataObjs[upper] -> GetData(_Extra=Extra)
;
; Now, concatenate this array of GKV objects into a single GKV object.
;
result = gkvarr[0]
FOR i = 1, nnOBjs-1 DO BEGIN
	previousResult = result
	result = previousResult -> Cat(gkvarr[i])
	previousResult -> Trash
ENDFOR
GKVdelete, gkvarr
;
; and we're done!
;
RETURN, result
END  ; ****** mGeneData::GetData ****** ;


FUNCTION mGeneData::MakeCopy
;
; Makes "deep" copy of self
;
; Written by W.M. Nevins
;	2/13/2009
;
class = OBJ_CLASS(self)			; Remember POLYMORPHISM! 
					;	(self may be a subclass of mGeneData)
ok = EXECUTE('copy = { ' + class + ' }'); Create 'copy' of this instance of mGeneData
GeneDataObjs = OBJARR(nObjs)
result = OBJ_NEW(class, GeneDataObjs)
return, result
END ; ***** mGeneData::MakeCopy ****** ;

PRO mGeneData::Trash
;
; Destroys mGeneData objects
;
; Written by W.M. Nevins
;	11/5/2008
OBJ_DESTROY, self
END  ;  ***** GeneData::Trash  *****  ;


PRO mGeneData::INFO
;
; Prints basic information about
; mGeneData class object
;
; Written by W.M. Nevins
;	2/13/2009
;
mGeneDataClass = OBJ_CLASS(self)
Objs = *self.GeneDataObjs
nObjs = N_ELEMENTS(Objs)
number = STRCOMPRESS(STRING(nObjs), /REMOVE_ALL)

PRINT, "Object Class = ", mGeneDataClass
PRINT, "Containing " + number + " GeneData Objects"
Objs[0] -> Info
Print, "Covering the time intervals"
Print, *self.Intervals  
END ; ***** mGeneData::INFO ***** ;


PRO mGeneData::Get, GeneDataObjs=GeneDataObjs, Intervals=Intervals
;
; Get routine for mGeneData object class
; Written by W.M. Nevins
;	2/13/2009
;
GeneDataObjs	= *Self.GeneDataObjs
Intervals	= *Self.Intervals

END ; ***** mGeneData::Get ***** ;

FUNCTION Make_mGeneData, ObjArr
Result = OBJ_NEW("mGeneData", Objarr)
RETURN, Result
END ; ****** mGeneData::Make_mGeneData ****** ;

PRO mGeneData::CleanUp
;
; Clean-up routine for mGeneData object class
; Written by W.M. Nevins
;	2/13/2009
nObjs = N_ELEMENTS(self.GeneDataObjs)
FOR i=0, nObjs-1 DO BEGIN
	thisObj = self.GeneDataObjs[i]
	IF(OBJ_ISA(thisObj, "GeneData")) THEN thisObj -> Trash
ENDFOR
PTR_FREE, self.Intervals, slef.GeneDataObjs
END ; ***** mGeneData::CleanUP ***** ; 

FUNCTION mGeneData::INIT, ObjArray
;
; Initialization routine for mGeneData object class
;
; Written by W.M. Nevins
;	2/13/2009
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF
;
nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN BEGIN
	MESSAGE, "Unable to initialize mGeneData object: no argument supplied", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Find time intervals covered by each GeneData object
; in ObjArray, and store this information in Intervals 
; array
;
nObjs = N_ELEMENTS(ObjArray)
Intervals = FLTARR(nObjs, 2)
FOR i=0, nObjs-1 DO BEGIN
	thisGKV_Obj = ObjArray[i] -> GKV_Obj()
	tAxis = ThisGKV_Obj -> NumDims()
	thisGKV_Obj -> GET, axis=tAxis, gridValues=tValues
	tValues = *tValues
	nt = N_ELEMENTS(tValues)
	Intervals[i,0] = tvalues[0]
	Intervals[i,1] = tValues[nt-1]
ENDFOR
;
; Set Tags of mGeneData object class
; 
self.Intervals = PTR_NEW(Intervals)
self.GeneDataObjs = PTR_NEW(ObjArray)
;
; and, we're done!
;
RETURN, 1	; successful completion of mGeneData object initialization
END ; ***** mGeneData::INIT ***** ;

PRO mGeneData__DEFINE
struct = {	mGeneData,				$	; Define "mGeneData" object class
		Intervals	:	PTR_NEW(),	$	; Pointer to an array containing
								;      time intervals covered by 
								;      corresponding GeneData objects
		GeneDataObjs	:	PTR_NEW()	}	; Pointer to an array containing
								;      GeneData objects 
END  ; ***** mGeneData__Define ***** ;
