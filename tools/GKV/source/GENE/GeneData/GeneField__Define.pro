FUNCTION GeneField::MakeCopy
;
; Makes "deep" copy of self
;
; Written by W.M. Nevins
;	11/5/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
result = self -> GeneData::MakeCopy()
result.fieldFile = self.fieldFile
result.nField	 = self.nField
result.Ifield    = self.iField
RETURN, result
END ; ***** GeneField::MakeCopy ***** ;

PRO GeneField::INFO
;
; Prints basic information about
; GeneField class object
;
; Written by W.M. Nevins
;	11/5/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
PRINT, "Field File = ", self.FieldFile
PRINT, "    nField = ", self.nField
PRINT, "    iField = ", self.iField
self -> GeneData::INFO
END ; ***** GeneField::INFO ***** ;



PRO GeneField::Get, 	FieldFile=FieldFile, 	$
			nField=nField,		$
			iField=iFiled,		$
			_REF_EXTRA=Extra
;
; Get routine for GeneField object class
; Written by W.M. Nevins
;	11/5/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
FieldFile	= Self.FieldFile
nField		= self.nField
iField		= self.iField
self -> GeneData::Get, _REF_EXTRA=Extra
;
END ; ***** GeneField::Get ***** ;

PRO GeneField::Set, _EXTRA=Extra
;
; Set routine for GeneField object class
; Written by W.M. Nevins
;	11/5/2008
;
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
result = GetKeyWord("FieldFile", Extra)
IF(TypeOF(result) EQ 7) THEN BEGIN
	IF(result NE "undefined") THEN self.FieldFile=result
ENDIF

result = GetKeyWord("nField", Extra)
IF(Query_Integer(result)) THEN self.nField = result

result = GetKeyWord("iField", Extra)
IF(Query_Integer(result)) THEN self.iField = result
;
; Forward anything left in Extra to 
; GeneData::Set to allow user to alter 
; remaining fields
;
IF(TypeOf(Extra) EQ 8) THEN self -> GeneData::Set, _Extra=Extra
;
; and we're done
;
END ; ***** GeneField::Set ***** ;


PRO GeneField::CleanUp
;
; Clean-up routine for GeneField object class
; Written by W.M. Nevins
;	11/5/2008
PTR_FREE, self.ParamStr
IF( OBJ_VALID(self.GkvObj) ) THEN BEGIN
	IF( OBJ_ISA(self.GkvObj, "GKVsd") ) THEN self.GkvObj -> Trash
ENDIF
END ; ***** GeneField::CleanUP ***** ; 

FUNCTION GeneField::INIT, GeneFieldStr
;
; Initialization routine for GeneField object class
;
; Written by W.M. Nevins
;	11/5/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
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
; Set Tags of GeneField object class
; using data in GeneFieldStr 
;
self.GeneDataDir	= GeneFieldStr.GeneDataDir
self.ParamStr      	= GeneFieldStr.ParamStr
self.GkvObj        	= GeneFieldStr.GkvObj
self.FieldFile		= "field.dat"
IF(TypeOf(GeneFieldStr.FieldFile) EQ 7) THEN BEGIN
	IF(STRCMP(GeneFieldStr.fieldFile, "") NE 1) THEN	$
	self.FieldFile	= GeneFieldStr.FieldFile
ENDIF
self.nField		= 1b
IF( Query_Integer(GeneFieldStr.nField) ) THEN BEGIN
	IF(self.nField GT 0) THEN	$
	self.nField	= GeneFieldStr.nField
ENDIF		
self.iField		= 1b
IF( Query_Integer(GeneFieldStr.iField) ) THEN BEGIN
	IF(self.iField GT 0) THEN	$
	self.iField	= GeneFieldStr.iField
ENDIF		
RETURN, 1	; successful completion of GeneField object initialization
END ; ***** GeneField::INIT ***** ;

PRO GeneField__DEFINE
;
;  Written by W.M. Nevins
;	11/5/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
struct = {	GeneField,				$	; Define "GENE Field" object class
		INHERITS GeneData,			$	; GeneField is a subclass of GeneData
		FieldFile	:	"field.dat",	$	; Alpha/Numeric field containing
								;      the name of the GENE 
								;      field.dat or mom_xxx.dat file
		nField		:	1b,		$	; Total number of fields in
								;      this FIELD.dat or mom_xxx.dat file
		iField		:	1b		}	; Number of chosen field in 
								;      this FIELD.dat or mom_xxx.dat file	
END  ; ***** GeneField__Define ***** ;
