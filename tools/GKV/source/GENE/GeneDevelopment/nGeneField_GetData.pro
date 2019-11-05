PRO GKVsd::CopyHeader, gkvStr
;
; Purpose:
;
;	Copies header of self into 
;	argument.
;
; Argument:
;
;	Argument must be a member of the
;	GKVsd class.
;
; Written by W.M. Nevins
;	2/17/2009
;
gkvStr.mnemonic = self.mnemonic
gkvStr.title    = self.title
gkvStr.units    = self.units
gkvStr.CodeName = self.CodeName
gkvStr.CodePI   = self.CodePI
gkvStr.RunID    = self.RunID
gkvStr.FileID   = self.FileID

RETURN
 
END ; ******  GKVsd::CopyHeader ****** ;

PRO GKVs1D::GridSet, arg, grid
;
; This proceedure replaces the grid structure
; specified by 'arg' with the grid structure
; supplied as 'grid' 
;
; Written by W.M. Nevins
;	11/7/2008
CASE arg OF
	1   :	self.grid1 = grid
	2   :	self.grid2 = grid
	3   :	self.grid3 = grid
	4   :	self.grid4 = grid
ENDCASE
RETURN
END ; ****** FUNCTION GKVs4D::GridSet ****** ;


FUNCTION GKVs1D::subGrids, _Extra=Extra
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
nDims      = self -> NumDims()
axisInfo   = self -> GetAxis(Extra)
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
END ; ******* FUNCTION GKVsd::subGrids ******* ;



FUNCTION GeneField::nGetData, arg, _Extra=Extra
;
; This function acts on GeneField objects,
; reads data associated with a particular field 
; (or moment) from the FIELDS.dat our mom_xxx.dat 
; file output by GENE.  It returns all, or some 
; logically rectangular subset of the field data 
; as a GKV object.
; 
; Data subsets are specified by choosing some range of 
; one axis (k_x, k_y, z, or t) using the argument and/or
; keywords described below.
;
;
;	Input Keywords:
;
;	 'mnemonic'	Set the 'mnemonic', (k_x, k_y, z, or t), of the selected axis 
;			equal to a two-element array, [min, max], to identify the 
;                       desired range of the selected independent variable; 
;			or input only a single value to "slice" field data at this value.;
;
;	    Start	Length (in bytes) of header in binary files
;			Defaults to 4. (Optional)
;
;	     EOR	Length (in bytes) of "end of record" mark
;			in binary files. Defaults to 8. (Optional)
;
; Written by W.M. Nevins
;	11/6/2008
; Modified by W.M. Nevins
;	11/8/2008
; to enable reading mom_xxx.dat files
; in addition to field.dat files
;
CD, CURRENT=CurrentWorkingDirectory
path=self.GeneDataDir
CD, path
fileName=self.fieldfile
ok = FINDFILE(fileName, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	MESSAGE, "Couldn't find FIELD file", /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Get a free I/O unit, and open FIELD.dat file
;
GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err
IF(err NE 0) THEN BEGIN
	PRINT, "Error opening FIELD file:
	PRINT, "Error Message = ", !ERR_STRING
	RETURN, 0
ENDIF
;
; Compute size of FIELD data using data from 
; parameter structure.
;
ParamStr = *self.ParamStr
IF(TypeOf(ParamStr) EQ 8) THEN BEGIN
	nkx     = ParamStr.BOX.NX0
	nky	= ParamStr.BOX.NKY0
	nZ	= ParamStr.BOX.NZ0
	nt	= ParamStr.INFO.ITIME/ParamStr.IN_OUT.ISTEP_FIELD
	Endian	= ParamStr.INFO.ENDIANNESS
	Endian	= STRLOWCASE(Endian)
	Endian	= STRCOMPRESS(Endian, /REMOVE_ALL)
	Endian	= STRCOMPRESS(Endian, /REMOVE_ALL)
	Precision=ParamStr.INFO.Precision
	Precision=STRLOWCASE(Precision)
	Precision=STRCOMPRESS(Precision, /REMOVE_ALL)
ENDIF ELSE BEGIN
	MESSAGE, "Invalid Parameter Structure, Returning 0", /INFORMATIONAL
	RETURN, 0
ENDELSE
;
; Some key constants garnared from close
; inspection of a sample field.dat file:
;
start = 4L	; The FIELD.dat file begins with a 4 byte field
EOR = 8L	; Each record is separated by an 8 byte field
;
; These may be operating system dependent, 
; so must allow the user to change 
; "start" and "EOR" from the command line
;
result = GetKeyWord("start", Extra)
IF(Query_Integer(result)) THEN start = result

result = GetKeyWord("EOR", Extra)
IF(Query_Integer(result)) THEN EOR = result
;
; Size of each field block (in bytes)
;
WordSize = 8
dataType = 9
IF( STRCMP(Precision, "single") ) THEN BEGIN
	WordSize=4
	dataType=6
ENDIF
nFieldBytes = LONG64(2L*WordSize*nkx*nky*nz)
;
; Compute number of bytes between time values
;
nFields = LONG64(self.nField)
iSkip = WordSize + EOR + nFields*(nFieldBytes + EOR)
;
; Parse Command line for restrictions to any axis
;
output    = self.gkvObj -> MakeCopy(/NoValues)
rangeInfo = output -> subGrids(_Extra=Extra)
;
; Create structure to hold result
;
outDims = output -> Numdims()
nDims = outDims - TOTAL(rangeInfo.SliceArray)
CASE nDims of
	0:	resultStr = {GKVsd}
	1:	resultStr = {GKVs1D}
	2:	resultStr = {GKVs2D}
	3:	resultStr = {GKVs3D}
	4:	resultStr = {GKVs4D} 
ENDCASE
output -> CopyHeader, resultStr
 
;
; Desired range for each index
;
kx = REFORM(rangeInfo.IndexArray[0,*])
ky = REFORM(rangeInfo.IndexArray[1,*])
z  = REFORM(rangeInfo.IndexArray[2,*])
t  = REFORM(rangeInfo.IndexArray[3,*])

FieldDims =  [	kx[1]-kx[0]+1,	$
		ky[1]-ky[0]+1,	$
		 z[1]- z[0]+1,	$
		 t[1]- t[0]+1	]
;
; form array to hold field values at each time step
;
thisField = MAKE_ARRAY(nkx, nky, nz, TYPE=dataType)
thisField = REFORM(thisfield, nkx, nky, nz, /OVERWRITE)
;
; and another array to hold field values to be returned
;
fieldValues = Make_Array(DIMENSION=FieldDims, type=dataType)
fieldValues = REFORM(fieldValues, fieldDims, /OVERWRITE)
;
; Read Field values
;
iField = self.iField
FieldStart = start + wordSize + EOR + (iField-1)*(nFieldbytes + EOR)
FOR i=t[0],t[1] DO	BEGIN
	thisField=READ_BINARY(ioUnit, DATA_START=fieldStart+i*iSkip, DATA_TYPE=dataType, DATA_DIMS=[nkx,nky,nz], ENDIAN=Endian)
	temp = SHIFT(thisField, nkx/2,0,0)
	FieldValues[*,*,*,i-t[0]] = temp[ kx[0]:kx[1], ky[0]:ky[1], z[0]:z[1] ]
ENDFOR
;
; load field values into resultStr
;
fieldValues = REFORM(FieldValues, /OVERWRITE)
resultStr.values = PTR_NEW(FieldValues)
CASE nDims of
	0:	result = OBJ_NEW("GKVsd",  resultStr)
	1:	result = OBJ_NEW("GKVs1d", resultStr)
	2:	result = OBJ_NEW("GKVs2d", resultStr)
	3:	result = OBJ_NEW("GKVs3d", resultStr)
	4:	result = OBJ_NEW("GKVs4d", resultStr) 
ENDCASE
;
; load new grids into result
;
slice = rangeInfo.SliceArray
newGrids = rangeInfo.GridArray
iAxis=1
output -> Get, indices=IndicesPtr
Indices = *IndicesPtr
Index = WHERE(Indices EQ '*')
FOR i=1,4 DO BEGIN
	IF(slice[i-1]) THEN BEGIN		; This grid has been sliced out		
		grid = newGrids[i-1]	; so just correct corresponding
		gridTitle = grid.title	; element of 'indices' array
		axisValue = *grid.values
		valueStr = STRING(axisValue, FORMAT='(G10.3)')
		Indices[Index[i-1]] = STRCOMPRESS(gridTitle + '=' + valueStr)
	ENDIF ELSE BEGIN
		grid = newGrids[i-1]
		result -> GridSet, iAxis, grid
		iAxis=iAxis+1
	ENDELSE
ENDFOR
result -> Set, indices=PTR_NEW(Indices)
; Slice axes with no range
;
;
; ane we're done!
;
FREE_LUN, ioUnit
CD, CurrentWorkingDirectory
RETURN, result

END ; ****** GeneFields__GetData ****** ;  
