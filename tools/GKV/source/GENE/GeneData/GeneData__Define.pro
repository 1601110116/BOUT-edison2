FUNCTION GeneData::MakeCopy
;
; Makes "deep" copy of self
;
; Written by W.M. Nevins
;	11/5/2008
;
class = OBJ_CLASS(self)			; Remember POLYMORPHISM! 
					;	(self may be a subclass of GeneData)
ok = EXECUTE('copy = { ' + class + ' }'); Create 'copy' of this instance of GeneData
copy.GeneDataDir = self.GeneDataDir
IF(PTR_VALID(self.ParamStr)) THEN BEGIN
	paramStr = *self.paramStr
	copy.paramStr = PTR_NEW(paramStr)
ENDIF

IF( OBJ_ISA(self.gkvObj, "GKVsd") ) THEN BEGIN
	copy.gkvObj = self.gkvObj -> MakeCopy(/NoValues)
ENDIF
result = OBJ_NEW(class, copy)
return, result
END ; ***** GeneData::MakeCopy ****** ;

PRO GeneData::Trash
;
; Destroys GeneData objects
;
; Written by W.M. Nevins
;	11/5/2008
IF(OBJ_ISA(self.GkvObj, "GKVsd")) THEN self.GkvObj -> Trash
OBJ_DESTROY, self
END  ;  ***** GeneData::Trash  *****  ;


PRO GeneData::INFO
;
; Prints basic information about
; GeneData class object
;
; Written by W.M. Nevins
;	11/5/2008
GeneDataClass = OBJ_CLASS(self)
paramStrInfo = "absent"
IF( PTR_VALID(self.ParamStr) ) THEN BEGIN
	IF(TypeOF(*self.ParamStr) EQ 8) THEN paramStrInfo = "present"
ENDIF

Print, "Object Class = ", GeneDataClass
PRINT, "Gene Data Directory = ", self.GeneDataDir
Print, "ParamStr ", paramStrInfo
IF( OBJ_ISA(self.GkvObj, "GKVsd") ) THEN BEGIN
	self.GkvObj -> info 
ENDIF ELSE BEGIN
	Print, "  GkvObj absent"
ENDELSE

END ; ***** GeneData::INFO ***** ;



PRO GeneData::Get, 	GeneDataDir=GeneDataDir, ParamStr=Paramstr,	$
			GkvObj=GkvObj, _REF_EXTRA=Extra
;
; Get routine for GeneData object class
; Written by W.M. Nevins
;	11/5/2008
;
GeneDataDir	= Self.GeneDataDir
ParamStr	= *self.ParamStr
GkvObj		= self.GkvObj ; this does not seem to work ...
;
; Forward anything in Extra to the GKVsd::Get routine.
;
self.GkvObj -> Get, _Ref_Extra=Extra

END ; ***** GeneData::Get ***** ;

FUNCTION GeneData::GKV_Obj
RETURN, self.GkvObj
END  ; ****** GeneData::GKV_Obj ****** ;

PRO GeneData::Set, _EXTRA=Extra
;
; Set routine for GeneData object class
; Written by W.M. Nevins
;	11/5/2008
;
result = GetKeyWord("GeneDataDir", Extra)
IF(TypeOf(result) EQ 7) THEN BEGIN
	IF(result NE "undefined") THEN self.GeneDataDir=result
ENDIF

result = GetKeyWord("ParamStr", Extra)
IF(TypeOf(result) EQ 10) THEN self.ParamStr = result
IF(TypeOf(ewaulr) EQ  8) THEN self.ParamStr = PTR_NEW(result)

result = GetKeyWord("GkvObj", Extra)
IF(TypeOf(result) EQ 11) THEN BEGIN
	IF( OBJ_ISA(Self.gkvObj, "GKVsd")) THEN self.gkvObj -> Trash
	Self.gkvObj = result
ENDIF
;
; Forward anything left in Extra to 
; GKVsd::Set to allow user to alter 
; fields in gkvObj
; 
IF(TypeOF(Extra) EQ 8) THEN self.gkvObj -> Set, _Extra=Extra
;
; and we're done
;
END ; ***** GeneData::Set ***** ;


PRO GeneData::CleanUp
;
; Clean-up routine for GeneData object class
; Written by W.M. Nevins
;	11/5/2008
PTR_FREE, self.ParamStr
IF(OBJ_ISA(self.GkvObj, "GKVsd")) THEN self.GkvObj -> Trash
END ; ***** GeneData::CleanUP ***** ; 

FUNCTION GeneData::INIT, GeneDataStr
;
; Initialization routine for GeneData object class
;
; Written by W.M. Nevins
;	11/5/2008
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF
;
nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN BEGIN
	RETURN, 0
ENDIF
;
; Set Tags of GeneData object class
; using data in GeneDataStr 
;
self.GeneDataDir = GeneDataStr.GeneDataDir
self.ParamStr      = GeneDataStr.ParamStr
self.GkvObj        = GeneDataStr.GkvObj
RETURN, 1	; successful completion of GeneData object initialization
END ; ***** GeneData::INIT ***** ;

PRO GeneData__DEFINE
struct = {	GeneData,				$	; Define "GeneData" object class
		GeneDataDir	:	"",		$	; Alpha/Numeric field containing
								;      full address of the directory 
								;      containing the GENE data file(s).
		ParamStr	:	PTR_NEW(),	$	; Pointer to structure containing
								;      info from corresponding 
								;       PARAMETERS.DAT file.
		GkvObj		: OBJ_NEW("GKVsd", {GKVsd}) }	; Corresponding GKV object
END  ; ***** GeneData__Define ***** ;
