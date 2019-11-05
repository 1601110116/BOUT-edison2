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

;
; Define Base Class GKVsdNull (Gyrokinetic Visualization signal data for large data sets)
; 
; This class contains information about the code run and the data
; to be analyzed.
;
; Written by W.M. Nevins
;	1/15/00.  
;	Modified 1/29/00 to include error bars
;
; Modified from GKVsd class 3/12/09 E. Wang

PRO GKVsdNull::debug
;
; Purpose:
;
;		This proceedure intentionally creates an error, causing
;		IDL to halt with in this proceedure.  Since 'debug' is
;		a method for all GKVsd objects, you are now able to examine
;		the contents of 'self' using IDL's variable window, 'print',
;		etc.  

y=x		; An intentional error (x is undefined).
END ; ****** GKVsdNull::debug ****** ;

FUNCTION GKVsdNull::MakeCopy, NoValues=novalues, NoErrorBars=noerrorbars, _Extra=extra
;
; Make "deep" copy of self
;
class = OBJ_CLASS(self)				; Remember POLYMORPHISM! 
								;	(self may be a subclass of GKVsd)
ok = EXECUTE('copy = { ' + class + ' }')	; Create 'copy', an instance of 
								;	the 'class' class-structure	
copy.mnemonic	= self.mnemonic		
copy.Title	= self.Title			; Copy fields from "self" into "copy"
Indices		= *self.indices			; Get Indices
copy.indices	= PTR_NEW(Indices)		; 	and make a new pointer for 'copy'
copy.units	= self.units
copy.CodeName	= self.codename			
copy.CodePI	= self.CodePI
copy.RunID		= self.RunID
copy.FileID	= self.FileID
copy.DirID = self.DirID
copy.reali = self.reali
copy.VariableID = self.VariableID
copy.nVars = self.nVars
copy.Grid1 = GKVsd_GridCopy(self.Grid1) 
copy.Grid2 = GKVsd_GridCopy(self.Grid2) 
copy.Grid3 = GKVsd_GridCopy(self.Grid3) 
copy.Grid4 = GKVsd_GridCopy(self.Grid4) 
exestring = "result=OBJ_NEW('" + class + "', copy)"
ok = EXECUTE(exestring)				; Register object
RETURN, result
END ; ***** GKVsdNull::MakeCopy ***** ;


PRO GKVsdNull::Trash
OBJ_DESTROY, self
;heap_gc, /verbose
END ; ****** GKVsdNull::Trash ****** ;



PRO GKVsdNull::Info
;
; Prints information about GKV Data objects
;
GKVsClass = OBJ_CLASS(self)
GKVSCLASS = STRUPCASE(GKVsClass)
Print, GKVsClass
PRINT, " mnemonic", " = ", self.mnemonic
PRINT, "    Title", " = ", self.Title
PRINT, "  Indices", " = [", STRJOIN(*self.Indices, ", "), "]"
PRINT, "    units", " = ", self.units
PRINT, " CodeName", " = ", self.CodeName
PRINT, "   CodePI", " = ", self.CodePI
PRINT, "    RunID", " = ", self.RunID
PRINT, "   FileID", " = ", self.FileID
PRINT, "    DirID", " = ", self.DirID
PRINT, "    VarID", " = ", self.variableID
PRINT, 'Grid1:'
GKVsd_PrintGrid, self.Grid1
IF (self.nVars GT 1) THEN BEGIN
PRINT, 'Grid2:'
GKVsd_PrintGrid, self.Grid2
ENDIF
IF (self.nVars GT 2) THEN BEGIN
PRINT, 'Grid3:'
GKVsd_PrintGrid, self.Grid3
ENDIF
IF (self.nVars GT 3) THEN BEGIN
PRINT, 'Grid4:'
GKVsd_PrintGrid, self.Grid4
ENDIF
RETURN
END ; ****** GKVsdNull::Info ****** ;

PRO GKVsdNull::GET,	mnemonic=mnemonic, Title=Title, Indices=indices, $
			units=units,values=values, vrange=vrange,		$
			ErrorBars=ErrorBars, CodeName=CodeName,		$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, $
			DirID=DirID, reali = reali, VariableID=VariableID, nVars=nVars, $
			axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,        $
      GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,   $
      range=range, irange=irange, _Extra=extra	
;
; Get values of elements of GKVsd Class
;
mnemonic	= self.mnemonic
Title		= self.Title
Indices	= self.Indices
units		= self.units
CodeName	= self.CodeName
CodePI	= self.CodePI
RunID		= self.RunID
FileID	= self.FIleID
reali = self.reali
arg = self.Grid1
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
  1:  IF(axisID EQ 1) THEN  $
      GKVsd_GetGrid, self.Grid1,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  2:  IF(axisID EQ 1) THEN  $
      GKVsd_GetGrid, self.Grid1,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  3:  IF(axisID EQ 1) THEN  $
      GKVsd_GetGrid, self.Grid1,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  7:  IF(axisID EQ self.Grid1.mnemonic) THEN  $
      GKVsd_GetGrid, self.Grid1,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
ELSE   :  ; just continue on else 

ENDCASE
self.Grid1 = arg

arg = self.Grid2
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
  1:  IF(axisID EQ 2) THEN  $
      GKVsd_GetGrid, self.Grid2,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  2:  IF(axisID EQ 2) THEN  $
      GKVsd_GetGrid, self.Grid2,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  3:  IF(axisID EQ 2) THEN  $
      GKVsd_GetGrid, self.Grid2,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  7:  IF(axisID EQ self.Grid2.mnemonic) THEN  $
      GKVsd_GetGrid, self.Grid2,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
ELSE   :  ; just continue on else 

ENDCASE
self.Grid2 = arg

arg = self.Grid3
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
  1:  IF(axisID EQ 3) THEN  $
      GKVsd_GetGrid, self.Grid3,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  2:  IF(axisID EQ 3) THEN  $
      GKVsd_GetGrid, self.Grid3,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
               irange=irange, _Extra=extra
  3:  IF(axisID EQ 3) THEN  $
      GKVsd_GetGrid, self.Grid3,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  7:  IF(axisID EQ self.Grid3.mnemonic) THEN  $
      GKVsd_GetGrid, self.Grid3,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
ELSE   :  ; just continue on else 

ENDCASE
self.Grid3 = arg

arg = self.Grid4
IDinfo = SIZE(axisID)
IDtype = IDinfo[IDinfo[0] + 1]
CASE IDtype OF
  1:  IF(axisID EQ 4) THEN  $
      GKVsd_GetGrid, self.Grid4,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  2:  IF(axisID EQ 4) THEN  $
      GKVsd_GetGrid, self.Grid4,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  3:  IF(axisID EQ 4) THEN  $
      GKVsd_GetGrid, self.Grid4,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
  7:  IF(axisID EQ self.Grid4.mnemonic) THEN  $
      GKVsd_GetGrid, self.Grid4,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,     $
                GridUnits=GridUnits,  GridValues=GridValues,      $
                boundary=boundary, uniform=uniform,range=range,     $
                irange=irange, _Extra=extra
ELSE   :  ; just continue on else 

ENDCASE
self.Grid4 = arg

RETURN
END ; ***** GKVsdNull::GET ***** ;

	
PRO GKVsdNull::SET,	axis=axisID,  GridMnemonic = GridMnemonic, GridTitle=GridTitle,         $
      GridUnits=GridUnits, GridValues=GridValues, boundary=boundary, uniform=uniform,   $
      range=range, irange=irange, reali=reali, $
      mnemonic=mnemonic, Title=Title, Indices=indices,$
			units=units,	$
			ErrorBars=ErrorBars, CodeName=CodeName,		$
			CodePI=CodePI, RunID=RunID, FileID=FIleID, _Extra=extra	
;
; Set values of elements of GKVsd Class
;
type = TypeOf(mnemonic)
IF(type EQ 7) THEN	$
	self.mnemonic = mnemonic

type = TypeOf(Title)
IF(type EQ 7) THEN	$
	self.Title	= Title
	
type= TypeOf(Indices)
IF(type EQ 10) THEN BEGIN
	PTR_FREE, self.Indices
	self.Indices = Indices
ENDIF
IF(type EQ 7) THEN BEGIN
	PTR_FREE, self.Indices
	self.Indices = PTR_NEW(Indices)
ENDIF

type = TypeOf(reali)
IF(type EQ 2) THEN BEGIN
    self.reali = reali
ENDIF
;type = TypeOf(vrange)
;IF(((type GT 0) AND (type LT 7)) OR (type EQ 9)) THEN			$
;	 self.vrange = vrange

type = TypeOf(units)
IF(type EQ 7) THEN	$
	self.units	= units

type = TypeOf(values)
IF(type EQ 10) THEN BEGIN
	PTR_FREE, self.values
	self.values = values
	IF(TypeOf(vrange) EQ 0) THEN BEGIN
		vmin = GKVsd_Min(*values, Max=vmax)
	ENDIF ELSE BEGIN
		vmin=vrange[0]
		vmax=vrange[1]
	ENDELSE
	IF(vmax EQ vmin) THEN vmax = vmin + 1
	vrange = [vmin, vmax]
ENDIF
IF(((type GT 0) AND (type LT 7)) OR (type EQ 9)) THEN BEGIN
	PTR_FREE, self.values
	self.values = PTR_NEW(values)
	IF(TypeOf(vmin) EQ 0) THEN vmin = GKVsd_Min(values, Max=vmax)
	IF(vmax EQ vmin) THEN vmax = vmin + 1
	vrange = [vmin, vmax]
ENDIF


type = TypeOf(CodeName)
IF(type EQ 7) THEN	$
	self.CodeName = CodeName

type = TypeOf(CodePI)
IF(type EQ 7) THEN	$
	self.CodePI = CodePI

type = TypeOf(RunID)
IF(type EQ 7) THEN	$
	self.RunID	= RunID

type = TypeOf(FileID)
IF(type EQ 7) THEN	$
	self.FIleID = FileID


RETURN
END ; ***** GKVsdNull::SET ***** ;


PRO GKVsDNull::CleanUp
PTR_FREE, self.Grid1.values
PTR_FREE, self.Grid2.values
PTR_FREE, self.Grid3.values
PTR_FREE, self.Grid4.values
PTR_FREE, self.indices

RETURN
END ; ***** GKVsdNull::CleanUp ***** ;


FUNCTION GKVsdNull::INIT, signal
;
; Preliminary init for testing object methods ...
;
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF 
;
; Setting self=signal DOESN'T work, so ...
;
; Set tags for GKVsd class
;
	self.mnemonic	= STRCOMPRESS(signal.mnemonic, /REMOVE_ALL)
	self.Title	= STRCOMPRESS(signal.Title,    /REMOVE_ALL)
	self.Indices	= signal.Indices
;        self.vrange	= signal.vrange
	self.units	= STRCOMPRESS(signal.units,    /REMOVE_ALL)
	self.CodeName	= signal.CodeName
	self.CodePI	= signal.CodePI
	self.RunID	= signal.RunID
	self.dirID = signal.DirID
	self.VariableID = signal.VariableID
        self.reali = signal.reali
	self.nVars = signal.nVars
	self.FileID	= signal.FileID
        self.multi      = signal.multi
        self.multiLoc   = signal.multiloc
	NumTags=N_TAGS(self.Grid1)
FOR itag=0, NumTags-1 DO $
  self.Grid1.(itag) = signal.Grid1.(itag)

IF(PTR_VALID(self.grid1.values)) THEN  BEGIN
  ;
  ; Check for uniform grid
  ;
  self.Grid1.uniform = GKVsd_UniformGrid(self.Grid1.values)
  ;
  ; Check if Grid1.range is set ... and set if necessary
  ;
  IF((self.Grid1.range[0] EQ 0.) AND (self.Grid1.range[1] EQ 0.)) THEN $
    self.Grid1.range = [MIN(*self.Grid1.Values, Max=max), max]
  ;
  ; Check if Grid1.irange is set ... and set if necessary
  ;
  IF((self.Grid1.irange[0] EQ 0) AND (self.Grid1.irange[1] EQ 0)) THEN $
    self.Grid1.irange = [0, N_ELEMENTS(*self.Grid1.values)-1]
ENDIF
	
	NumTags=N_TAGS(self.Grid2)
	
FOR itag=0, NumTags-1 DO $
  self.Grid2.(itag) = signal.Grid2.(itag)

IF(PTR_VALID(self.grid2.values)) THEN  BEGIN
  ;
  ; Check for uniform grid
  ;
  self.Grid2.uniform = GKVsd_UniformGrid(self.Grid2.values)
  ;
  ; Check if Grid2.range is set ... and set if necessary
  ;
  IF((self.Grid2.range[0] EQ 0.) AND (self.Grid2.range[1] EQ 0.)) THEN $
    self.Grid2.range = [MIN(*self.Grid2.Values, Max=max), max]
  ;
  ; Check if Grid2.irange is set ... and set if necessary
  ;
  IF((self.Grid2.irange[0] EQ 0) AND (self.Grid2.irange[1] EQ 0)) THEN $
    self.Grid2.irange = [0, N_ELEMENTS(*self.Grid2.values)-1]
ENDIF
	NumTags=N_TAGS(self.Grid3)
FOR itag=0, NumTags-1 DO $
  self.Grid3.(itag) = signal.Grid3.(itag)

IF(PTR_VALID(self.grid3.values)) THEN  BEGIN
  ;
  ; Check for uniform grid
  ;
  self.Grid3.uniform = GKVsd_UniformGrid(self.Grid3.values)
  ;
  ; Check if Grid3.range is set ... and set if necessary
  ;
  IF((self.Grid3.range[0] EQ 0.) AND (self.Grid3.range[1] EQ 0.)) THEN $
    self.Grid3.range = [MIN(*self.Grid3.Values, Max=max), max]
  ;
  ; Check if Grid3.irange is set ... and set if necessary
  ;
  IF((self.Grid3.irange[0] EQ 0) AND (self.Grid3.irange[1] EQ 0)) THEN $
    self.Grid3.irange = [0, N_ELEMENTS(*self.Grid3.values)-1]
ENDIF

NumTags=N_TAGS(self.Grid4)
FOR itag=0, NumTags-1 DO $
  self.Grid4.(itag) = signal.Grid4.(itag)

IF(PTR_VALID(self.grid4.values)) THEN  BEGIN
  ;
  ; Check for uniform grid
  ;
  self.Grid4.uniform = GKVsd_UniformGrid(self.Grid4.values)
  ;
  ; Check if Grid4.range is set ... and set if necessary
  ;
  IF((self.Grid4.range[0] EQ 0.) AND (self.Grid4.range[1] EQ 0.)) THEN $
    self.Grid4.range = [MIN(*self.Grid4.Values, Max=max), max]
  ;
  ; Check if Grid4.irange is set ... and set if necessary
  ;
  IF((self.Grid4.irange[0] EQ 0) AND (self.Grid4.irange[1] EQ 0)) THEN $
    self.Grid4.irange = [0, N_ELEMENTS(*self.Grid4.values)-1]
ENDIF

RETURN, 1
END  ; ***** GKVsdNull:INIT *****;
FUNCTION GKVsdNull::subGrids, _Extra=Extra
;
; Purpose:
;
;	This function is used to select a 
;	subset of the data in self. 
;
; Input Keywords:
;
;	'mnemonic'	Where 'mnemonic' is the mnemonic
;			of one of self's axes.  Set this equal
;			to a 2-element array specifying the 
;			desired range of values from the 
;			selected axis. Default is to include
;			the entire axis. (Optional). 
;
;			If 'mnemonic' is set to a single value,
;			or if the two elements are both closest
;			to the same grid value, then the 
;			corresponding element of the SliceArray
;			is set to 1 (true).
;			  
;
;			Up to nDims mnemonics may be specified
;			on the command line. 
;
;
;
; Output:
;
;	A structure containing the following tags:
;
;
;	IndexArray	An array of length [nDims, 2]
;			containing the indices to the
;			desired minimum and maximum 
;			values of the corresponding axis.
;
;	SliceArray	An array of length [nDims] which
;			is 1 (true) if only a single value
;			was entered for the corresponding
;			axis mnemonic or if both elements of 
;			the desired range were closest to the
;			same grid point; and 0 (false) otherwise.
;
;	GridArray	An array of length [nDims] containing 
;			Grid structures appropriate for the 
;			selected data.
;
; Written by W.M. Nevins
;	11/7/2008
;
nDims      = self.nVars
axisinfo={axis1:'no match', axis2:'no match',axis3:'no match',axis4:'no match'}
IF (N_Elements(Extra) NE 0) THEN BEGIN
axisInfo   = self -> GetAxis(Extra)
ENDIF
IndexArray = LONARR(nDims, 2)
sliceArray = BYTARR(nDims)
;
; Create a structure array to hold 
; self's grid structures
;
myGrids    = REPLICATE(self.grid1, nDims)
IF(nDims GT 1) THEN myGrids[1] = self.grid2
IF(nDims GT 2) THEN myGrids[2] = self.grid3
IF(nDIms GT 3) THEN myGrids[3] = self.grid4
;
; check for new limits to each grid of self
;
FOR i=0, nDims-1 DO BEGIN
	gridLimits = axisInfo.(i)
	thisGrid = myGrids[i]
	gridValues = *thisGrid.values
	IF( TypeOf(gridLimits) NE 7 ) THEN BEGIN
		IF( N_ELEMENTS(gridLimits) EQ 1) THEN BEGIN
			gridLimits = [gridLimits, gridLimits]
			sliceArray[i] = 1
		ENDIF
		resid1 = (gridValues - gridlimits[0])^2
		eps1   = MIN(resid1, iMin)
		resid2 = (gridValues - gridLimits[1])^2
		eps2   = MIN(resid2, iMax)
		IndexArray[i,0] = iMin
		IndexArray[i,1] = iMax                
		myGrids[i] = GKVsd_GridCopy(thisGrid, Irange=[iMin, iMax])                
	ENDIF ELSE BEGIN
		iMin = 0
		iMax = N_ELEMENTS(gridValues) - 1
		myGrids[i] = GKVsd_GridCopy(thisGrid)
	ENDELSE
	IndexArray[i,0] = iMin
	IndexArray[i,1] = iMax
	IF(iMin EQ iMax) THEN sliceArray[i] = 1
ENDFOR

output = {	IndexArray	:	indexArray,	$
		SliceArray	:	sliceArray,	$
		GridArray	:	myGrids		}
RETURN, output	
END ; ******* FUNCTION GKVsdNull::subGrids ******* ;


FUNCTION GKVsdNull::GetAxis, structure, Debug=d
;
; Searches 'structure' for a tag which is the same as the mnemonic of axis1
; Returns a structure with tag 'axis1'.  
; The associated value is whatever value was associated with the mnemonic of axis1.
;
; On return, 'structure' has tag (and value) corresponding to the mnemonic of axis1
; removed.
;
; If the mnemonic of axis1 is the only tag in structure, then structure=-1 on return
; (null structures are not allowed in IDL 5.3)
;
; NOTE:  the argument structure may be undefined on entry.  In this case it is important
; not to do any thing with 'structure (beyond calling N_ELEMENTS, TypeOf, or SIZE), as
; this would result in an error.  Just return with a null result.
; 
; Written by W.M. Nevins
;	2/20/00
;
debug=0
;result={axis1:'no match', axis2:'no match',axis3:'no match',axis4:'no match'}
save1='no match'
save2='no match'
save3='no match'
save4='no match'
otherTags=0
IF(N_ELEMENTS(d) NE 0) THEN debug=d
IF(N_ELEMENTS(structure) EQ 0) THEN RETURN, result
arg1Type = typeOF(structure)							; Check for proper argument type
IF(arg1Type NE 8) THEN BEGIN
	structure = -1
	RETURN, result
ENDIF
nTags = N_TAGS(structure)
IF(nTags EQ 0) THEN RETURN, result
tagNames = TAG_NAMES(structure)
tagNames = STRTRIM(tagNames, 2)							; Remove both leading and trailing blanks
command_str = 'structure = {'
axisMnemonic1 = STRTRIM(self.Grid1.mnemonic)				; Remove both leading and trailing blanks
axisMnemonic2 = STRTRIM(self.Grid2.mnemonic)
axisMnemonic3 = STRTRIM(self.Grid3.mnemonic)
axisMnemonic4 = STRTRIM(self.Grid4.mnemonic)
FOR i=0, ntags-1 DO BEGIN								; Search tags of 'structure'
    IF( STRCMP(axisMnemonic1, tagNames[i], /FOLD_CASE) ) THEN save1=structure.(i)
    IF( STRCMP(axisMnemonic2, tagNames[i], /FOLD_CASE) ) THEN save2=structure.(i)
    IF( STRCMP(axisMnemonic3, tagNames[i], /FOLD_CASE) ) THEN save3=structure.(i)
    IF( STRCMP(axisMnemonic4, tagNames[i], /FOLD_CASE) ) THEN BEGIN save4=structure.(i)
	ENDIF ELSE BEGIN								;	'structure', only the LAST occurance is significant
		i_str = STRING(i, FORMAT='(I3)')				; Construct command string for 'Nstructure'
		i_str = STRTRIM(i_str, 2)						; Remove both leading and trailing blanks
		IF(otherTags NE 0) THEN command_str = command_str + ', '
		command_str = command_str + tagNames[i] + ':' + 'structure.(' + i_str + ')'
		otherTags = otherTags + 1
	ENDELSE
    ENDFOR
result = {axis1:save1,axis2:save2,axis3:save3,axis4:save4}
command_str = command_str + '}'
IF(otherTags NE 0) THEN BEGIN
	ok = EXECUTE(command_str)
ENDIF ELSE BEGIN
	structure = -1
ENDELSE
RETURN, result
END ; ****** GKVsdNull::GetAxis ****** ;




FUNCTION GkvsdNull::GetData, arg, _Extra=Extra
;
; This function acts on GKVsdNull objects, reading 
; data associated with the object from the stored 
; NetCDF file. It returns all, or some 
; logically rectangular subset of the field data 
; as a GKV object. {axis1:result3.axis1, axis2:result3.axis2, axis3:result3.axis3, axis4:'no match'}
; 
; Data subsets are specified by choosing some range of 
; one axis (k_x, k_y, z, or t) using the argument and/or
; keywords described below.
;
;
;	Input Keywords:
;
;	 'mnemonic'	Set the 'mnemonic', of the selected axis 
;			equal to a two-element array, [min, max], to identify the 
;                       desired range of the selected independent variable; 
;			or input only a single value to "slice" field data at this value.;
;
;
; Written by E. Wang
;	3/24/2008
;


;
; Define local variables for procedure
;
CD, CURRENT=current_working_directory

CD, self.DirID
fileIN = self.FileID
vVarDims = self.nVars
varID = self.VariableID
arraystyle = '0'
specIndx='species'
;
; Open the file
;

cdfID = NCDF_OPEN(fileIN, /NoWrite)			; Open netCDF file (sets cdfID to an integer � 0)
cdfInfo = NCDF_INQUIRE(cdfID)				;	and get info regarding contents
separator='/'
IF(!D.NAME EQ 'MAC') THEN separator=':'
IF(!D.NAME EQ 'WIN') THEN separator='\'
subNames = STRSPLIT(fileIN, separator, /EXTRACT)		; Crack 'fileIn' on separators
nSubNames = N_ELEMENTS(subNames)
defaultName = subNames[nSubNames-1]

varInfo = NCDF_VARINQ(cdfID, varID) ; Get info about variables
varName = varInfo.Name
varDims = varInfo.nDims
varDimID= varInfo.Dim
nAtts   = varInfo.nAtts
nDims =cdfInfo.Ndims
dimName = STRARR(nDims)
dimSize = LONARR(nDims)
dimVars = REPLICATE(-1L, nDims)
dimGrid = REPLICATE({GRID}, nDims)


multi = self.multi
multiLoc = self.multiLoc
 complex = self.reali
 nDIms =cdfInfo.Ndims	
 dimName = STRARR(nDims)
 firstDimName = dimName[varDimID[0]]
 firstDimName = STRLOWCASE(firstDimName)
; IF(firstDimName EQ 'ri') THEN BEGIN
;     complex = 1
; ENDIF

;FOR dimID=0,nVarDims-1 DO IF(dimName[Varinfo.Dim[dimID]] EQ specIndx) THEN multiSpecies[varID]=1

FOR dimID = 0, nDims-1 DO BEGIN
	NCDF_DIMINQ, cdfID, dimID, name, dsize		; Get info about dimensions
	dimName[dimID]=name
	dimSize[dimID]=dsize	
	IF(name EQ specIndx) THEN BEGIN 
            nspec=dsize
            ;multi = 1
            ;multiLoc = dimID
        ENDIF       
ENDFOR
 

;
; Create appropriate GKV structure
;

CASE vVarDims OF
			0	: output = {GKVsd}
			1	: output = {GKVs1D}
			2	: output = {GKVs2D}
			3	: output = {GKVs3D}
			4	: output = {GKVs4D}
			ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
		ENDCASE
;
; put in default information
;

output.FileID = self.DirID + self.FileID
output.runID = self.RunID
output.CodeName = self.CodeName
output.CodePI = self.CodePI
output.mnemonic = self.mnemonic
output.Title = self.Title
ind = *self.Indices
output.Indices = PTR_NEW(ind)
output.units = self.units

;
; Parse command line for axis range etc
;

rangeInfo = self -> subGrids(_Extra=Extra)
g1 = REFORM(rangeInfo.IndexArray[0,*]) ; array of beginning and end values supplied by input- g1 is meant to represent grid 1 range
g1[1] = g1[1]+1
c1=g1[1]-g1[0]                         ; NCDF_VARGET requires an inital value and a distance (not final value) c1 is the distance for grid1
dimSize[varDimID[1]]=c1                ; dimSize array is carried over from NetCDF_data()
PTR_FREE, rangeInfo.gridArray[0].values
IF vVarDims GT 1 THEN  BEGIN 
g2 = REFORM(rangeInfo.IndexArray[1,*])
g2[1] = g2[1]+1
c2=g2[1]-g2[0]
dimSize[varDimID[2]]=c2
PTR_FREE, rangeInfo.gridArray[1].values
ENDIF
IF vVarDims GT 2 THEN BEGIN
 g3 = REFORM(rangeInfo.IndexArray[2,*])
g3[1] = g3[1]+1
c3=g3[1]-g3[0]
dimSize[varDimID[3]]=c3
PTR_FREE, rangeInfo.gridArray[2].values
ENDIF
IF vVarDims GT 3 THEN BEGIN 
g4 = REFORM(rangeInfo.IndexArray[3,*])
g4[1] = g4[1]+1
c4=g4[1]-g4[0]
dimSize[varDimID[4]]=c4
PTR_FREE, rangeInfo.gridArray[3].values
ENDIF

IF(multi) THEN BEGIN    
    CASE (vVarDims) OF
        1:   BEGIN 
            CASE (multiLoc) OF
                1  :  offset=[0,g1[0]]
                2  :  offset=[g1[0],0]
            ENDCASE
        END
        2:  BEGIN
            CASE (multiLoc) OF 
                1: offset=[0,g1[0],g2[0]]
                2: offset=[g1[0],0,g2[0]]
                3: offset=[g1[0],g2[0],0]
            ENDCASE
        END
        3:  BEGIN
            CASE(multiLoc) OF 
                1: offset=[0,g1[0],g2[0],g3[0]]
                2: offset=[g1[0],0,g2[0],g3[0]]
                3: offset=[g1[0],g2[0],0,g3[0]]
                4: offset=[g1[0],g2[0],g3[0],0]
            ENDCASE
        END
        4:  BEGIN
            CASE(multiLoc) OF

                1: offset=[0,g1[0],g2[0],g3[0],g4[0]]
                2: offset=[g1[0],0,g2[0],g3[0],g4[0]]
                3: offset=[g1[0],g2[0],0,g3[0],g4[0]]
                4: offset=[g1[0],g2[0],g3[0],0,g4[0]]
                5: offset=[g1[0],g2[0],g3[0],g4[0],0]
            ENDCASE
        END
    ENDCASE
    CASE (vVarDims) OF
        1:   BEGIN 
            CASE (multiLoc) OF
                1  :  count=[1,g1[1]-g1[0]]
                2  :  count=[g1[1]-g1[0],1]
            ENDCASE
        END
        2:  BEGIN
            CASE (multiLoc) OF 
                1: count=[1,g1[1]-g1[0],g2[1]-g2[0]]
                2: count=[g1[1]-g1[0],1,g2[1]-g2[0]]
                3: count=[g1[1]-g1[0],g2[1]-g2[0],1]
            ENDCASE
        END
        3:  BEGIN
            CASE (multiLoc) OF
                1: count=[1,g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0]]
                2: count=[g1[1]-g1[0],1,g2[1]-g2[0],g3[1]-g3[0]]
                3: count=[g1[1]-g1[0],g2[1]-g2[0],1,g3[1]-g3[0]]
                4: count=[g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],1]
            ENDCASE
        END
        4:  BEGIN
            CASE (multiLoc) OF 
                1: count=[1,g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],g4[1]-g4[0]]
                2: count=[g1[1]-g1[0],1,g2[1]-g2[0],g3[1]-g3[0],g4[1]-g4[0]]
                3: count=[g1[1]-g1[0],g2[1]-g2[0],1,g3[1]-g3[0],g4[1]-g4[0]]
                4: count=[g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],1,g4[1]-g4[0]]
                5: count=[g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],g4[1]-g4[0],1]
            ENDCASE
        END
    ENDCASE
ENDIF ELSE BEGIN
IF(complex) THEN BEGIN  ; create offset and count for NCDF_VARGET (command to read in NetCDF data)- first array element is complex

CASE  (vVarDims) OF
    1     : offset=[0,g1[0]]
    2     : offset=[0,g1[0],g2[0]]
    3     : offset=[0,g1[0],g2[0],g3[0]]
    4     : offset=[0,g1[0],g2[0],g3[0],g4[0]]
    ELSE  :
ENDCASE


CASE  (vVarDims) OF
    1     : count=[2,g1[1]-g1[0]]
    2     : count=[2,g1[1]-g1[0],g2[1]-g2[0]]
    3     : count=[2,g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0]]
    4     : count=[2,g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],g4[1]-g4[0]]
    ELSE  :
ENDCASE

ENDIF ELSE BEGIN  ; not complex, create offset and count
CASE  (vVarDims) OF
    1     : offset=[g1[0]]
    2     : offset=[g1[0],g2[0]]
    3     : offset=[g1[0],g2[0],g3[0]]
    4     : offset=[g1[0],g2[0],g3[0],g4[0]]
    ELSE  :
ENDCASE


CASE  (vVarDims) OF
    1     : count=[g1[1]-g1[0]]
    2     : count=[g1[1]-g1[0],g2[1]-g2[0]]
    3     : count=[g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0]]
    4     : count=[g1[1]-g1[0],g2[1]-g2[0],g3[1]-g3[0],g4[1]-g4[0]]
    ELSE  :
ENDCASE
ENDELSE
ENDELSE
slice = rangeInfo.SliceArray

vVardims2 = vVardims
FOR iaxis=1,vVardims DO BEGIN
	IF(slice(iAxis-1)) THEN BEGIN
		vVardims2 = vVardims2 -1 
	ENDIF
ENDFOR

; Get data==================================================

NCDF_VARGET, cdfID, varID, values ;, OFFSET=offset, COUNT=count


; if variable is complex (assuming the first element is 0 or 1 for
; real or complex) values needs to be rewritten
IF(complex) THEN BEGIN
    CASE vVarDims OF 
        1	:	values = REFORM(COMPLEX(values[0,*],       values[1,*]),       dimSize[varDimID[1]])
        2	:	values = REFORM(COMPLEX(values[0,*,*],     values[1,*,*]),     dimSize[varDimID[1]], dimSize[varDimID[2]])
        3	:	values = REFORM(COMPLEX(values[0,*,*,*], values[1,*,*,*]), dimSize[varDimID[1]], dimSize[varDimID[2]], dimSize[varDimID[3]])
        4	:	values = REFORM(COMPLEX(values[0,*,*,*,*], values[1,*,*,*,*]), dimSize[varDimID[1]], dimSize[varDimID[2]], dimSize[varDimID[3]], dimSize[varDimID[4]])
        ELSE	:            ; *** NEED SOME ERROR PROCESSING HERE ***
    ENDCASE
ENDIF 



IF(arrayStyle EQ 'c') THEN BEGIN ; If data was written by a 'c' (or pascal) code,
    CASE varDims OF     ; then the 'values' array must be reformed to 
        1	:         ; to reflect IDL's array storage conventions
        2	: values = REFORM(values, dimSize[varDimID[1]], dimSize[varDimID[0]])
        3	: values = REFORM(values, dimSize[varDimID[2]], dimSize[varDimID[1]], dimSize[varDimID[0]])
        4	: values = REFORM(values, dimSize[varDimID[3]], dimSize[varDimID[2]], dimSize[varDimID[1]], dimSize[varDimID[0]])
        ELSE	:            ; *** NEED SOME ERROR PROCESSING HERE ***
    ENDCASE
ENDIF

IF (multi) THEN BEGIN
CASE (vVarDims) OF
        1:   BEGIN 
            CASE (multiLoc) OF
                1  :  temp = REFORM(values[0,*], dimSize[varDimID[1+complex]])
                2  :  temp = REFORM(values[*,0], dimSize[varDimID[0+complex]])
            ENDCASE
        END
        2	:	BEGIN
            CASE multiLoc OF
                1	:	temp = REFORM(values[0,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
                2	:	temp = REFORM(values[*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]])
                3	:	temp = REFORM(values[*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]])
            ENDCASE
        END
        3	:	BEGIN
            CASE multiLoc OF
                1	:	temp = REFORM(values[0,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
                2	:	temp = REFORM(values[*,0,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
                3	:	temp = REFORM(values[*,*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]])
                4	:	temp = REFORM(values[*,*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
            ENDCASE
        END
        4	:	BEGIN
            CASE multiLoc OF
                1	:	temp = REFORM(values[0,*,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
                2	:	temp = REFORM(values[*,0,*,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
                3	:	temp = REFORM(values[*,*,0,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
                4	:	temp = REFORM(values[*,*,*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[4+complex]])
                5	:	temp = REFORM(values[*,*,*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
            ENDCASE
        ENDCASE
    ENDCASE
    values = temp
    ENDIF

;
; Done editing values- put in structure
;

vmin = GKVsd_Min(values, Max=vmax)
output.values = PTR_NEW(values)
output.vrange = [vmin,vmax]


; create GKV object of original dimensionality (might be sliced based
; on input for a lower final dimensionality)



CASE (vVarDims) OF 
   	0	: output2 = OBJ_NEW('GKVsd' , output)
        1	: output2 = OBJ_NEW('GKVs1D', output)
        2	: output2 = OBJ_NEW('GKVs2D', output)
        3	: output2 = OBJ_NEW('GKVs3D', output)
        4	: output2 = OBJ_NEW('GKVs4D', output) 
ENDCASE


;
; Put Grid(s) into structure
;
output2.Grid1 = GKVsd_GridCopy(self.Grid1, irange=[g1[0],g1[1]-1])
IF (vVarDims GT 1) THEN output2.Grid2 = GKVsd_GridCopy(self.Grid2, irange=[g2[0],g2[1]-1])
IF (vVarDims GT 2) THEN output2.Grid3 = GKVsd_GridCopy(self.Grid3, irange=[g3[0],g3[1]-1])
IF (vVarDims GT 3) THEN output2.Grid4 = GKVsd_GridCopy(self.Grid4, irange=[g4[0],g4[1]-1])

;
; If any Grid is only one element in length, slice the dimension
;

slice = rangeInfo.SliceArray
dAxis = 0

IF (vVarDims GT 3) THEN BEGIN
IF output2.grid4.irange(1)EQ 0 THEN BEGIN
temp = output2 -> Slice(axis=4,index=0)
dAxis=dAxis+1
output2 -> Trash
output2 = temp
ENDIF 
ENDIF

IF (vVarDims GT 2) THEN BEGIN
IF output2.grid3.irange(1)EQ 0 THEN BEGIN
temp = output2 -> Slice(axis=3,index=0)
dAxis=dAxis+1
output2 -> Trash
output2 = temp
ENDIF 
ENDIF

IF (vVarDims GT 1) THEN BEGIN
IF output2.grid2.irange(1)EQ 0 THEN BEGIN
temp = output2 -> Slice(axis=2,index=0)
dAxis=dAxis+1
output2 -> Trash
output2 = temp
ENDIF 
ENDIF

IF output2.grid1.irange(1) EQ 0 THEN BEGIN
temp = output2 -> Slice(axis=1,index=0)
dAxis=dAxis+1
output2 -> Trash
output2 = temp
ENDIF 

;
; Close CDF file, clear memory, return to working directory
;



NCDF_CLOSE, cdfID
cdfID=-1
CD, current_working_directory
RETURN, output2
END



;
; GKVsdNull__Define properties-
;  Null class will contain all the info of the Gkvsd object except the values pointer.  To access the values, nullObj -> GetData() must be run.
;  Differences between Gkvsd objects and Null objects:  
;  Null contains the file directory, Dimension size, 4 Grids (if nDims < 4, extra Grids will be copy of last grid).
;  Eric Wang, 3-17-09  
PRO GKVsdNull__Define
struct = {	GKVsdNull,				$	; Define "Gyrokinetic Visualization signal data" class
		mnemonic	: "",		$	; Alpha/numeric mnemonic (as used in code)
		Title		: "",		$	; "Pretty" name	(use IDL Vector Font)
		Indices		: PTR_NEW(),	$	; Pointer to string array containing values of indices
		units		: "",		$	; Units		(use IDL Vector Font)
   		CodeName	: "",		$	; Name of code which produced data
		CodePI		: "",		$	; Name of PI (or team) responsible for code runs
		RunID		: "",		$	; Run identifier
		DirID   : "",   $ ; Directory of File
                reali   :  1b, $ ; complex or not
		VariableID   :  1b,  $ ; Variable number inside NetCDF files
		nVars    :   1b, $ ; Number of grids to created
                multi    :   1b, $ ; If there are multiple species
                multiLoc :   1b, $ ; Location of species index
		Grid1:{Grid} , $ ; Include one 'Grid' class
		Grid2:{Grid} , $ ; Include 2nd 'Grid' class
		Grid3:{Grid} , $ ; Include 3rd 'Grid' class
		Grid4:{Grid} , $ ; Include 4th 'Grid' class- if data is not 4D, last Grids will be shallow copy of previous grid
		FileID		: ""		}	; Name of file from which data was extracted
								
END ; ***** GKVsd__Define ***** ;

								
		
