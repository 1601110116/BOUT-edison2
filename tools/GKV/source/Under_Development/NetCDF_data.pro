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

PRO GKVs1D::AxisShift, axis=axisID, debug=d
;
; Proceedure examines grid of associated with axisID.  
; If grid values are discontinuous (at one point) and
; apparantly peroidic (as in Dorland's kx grid),
; Then both this grid and the associated values are shifted
; such that the grid vaues are monotonic.
;
; Written 6/7/00 by
;	W.M. Nevins
;
axis=axisID
IF(TypeOf(axis) EQ 7) THEN axis = Self -> AxisNumber(axisID)
IF(axis EQ 0) THEN RETURN
IF(axis GT 1) THEN BEGIN
	MESSAGE, 'axisID should be 1 to have gotten here!', Informational=d
	RETURN
ENDIF
axisStr = STRING(axis, FORMAT='(i1)')
commandStr = 'grid = self.grid' + axisStr
ok = EXECUTE(commandStr)
axisValues = *Grid.Values
nValues = N_ELEMENTS(axisValues)
delta = SHIFT(axisValues, 1) - axisValues			; Computes difference between nearest neighbors using periodic boundary conditions
deltaNorm = 0.5*TOTAL(ABS(delta))/(nValues-1)		; This should be (about) equal to the typical grid spacing	if we have a single discontinuity
deltaMax = MAX(delta, indexShift)				; 'indexShift' is index to last element of 'axis' before discontinuity.
IF(indexShift EQ nValues) THEN RETURN				; Discontinuity is at boundary (where GKV routines expect it)
axisValues = SHIFT(axisValues, -indexShift)		; Shift axisValues to put discontinuity at boundary.
PTR_FREE, Grid.Values
Grid.values  = PTR_NEW(axisValues)
axisMin=MIN(axisValues, MAX=axisMax)
Grid.range = [axisMin, axisMax]
Grid.irange= [0, nValues-1]
commandStr = 'self.grid' + axisStr + '= Grid'
ok = EXECUTE(commandStr)
;
; Now shift the values array
;
values = *self.Values
values = SHIFT(values, -indexShift)
PTR_FREE, Self.Values
self.Values = PTR_NEW(values)
RETURN 
END ; ****** GKVs1D::AxisShift ****** ;


PRO GKVs2D::AxisShift, axis=axisID, debug=d
;
; Proceedure examines grid of associated with axisID.  
; If grid values are discontinuous (at one point) and
; apparantly peroidic (as in Dorland's kx grid),
; Then both this grid and the associated values are shifted
; such that the grid vaues are monotonic.
;
; Written 6/7/00 by
;	W.M. Nevins
;
axis=axisID
IF(TypeOf(axis) EQ 7) THEN axis = Self -> AxisNumber(axisID)
IF(axis EQ 0) THEN RETURN
IF(axis GT 2) THEN BEGIN
	MESSAGE, 'axisID should be less that 3 to have gotten here!', Informational=d
	RETURN
ENDIF
axisStr = STRING(axis, FORMAT='(i1)')
commandStr = 'grid = self.grid' + axisStr
ok = EXECUTE(commandStr)
axisValues = *Grid.Values
nValues = N_ELEMENTS(axisValues)
delta = SHIFT(axisValues, 1) - axisValues			; Computes difference between nearest neighbors using periodic boundary conditions
deltaNorm = 0.5*TOTAL(ABS(delta))/(nValues-1)		; This should be (about) equal to the typical grid spacing	if we have a single discontinuity
deltaMax = MAX(delta, indexShift)				; 'indexShift' is index to last element of 'axis' before discontinuity.
IF(indexShift EQ nValues) THEN RETURN				; Discontinuity is at boundary (where GKV routines expect it)
axisValues = SHIFT(axisValues, -indexShift)		; Shift axisValues to put discontinuity at boundary.
PTR_FREE, Grid.Values
Grid.values  = PTR_NEW(axisValues)
axisMin=MIN(axisValues, MAX=axisMax)
Grid.range = [axisMin, axisMax]
Grid.irange= [0, nValues-1]
commandStr = 'self.grid' + axisStr + '= Grid'
ok = EXECUTE(commandStr)
;
; Now shift the values array
;
values = *self.Values
CASE axis OF
	1	:	values = SHIFT(values, -indexShift, 0)
	2	:	values = SHIFT(values, 0, -indexShift)
ENDCASE
PTR_FREE, Self.Values
self.Values = PTR_NEW(values)
RETURN 
END ; ****** GKVs2D::AxisShift ****** ;


PRO GKVs3D::AxisShift, axis=axisID, debug=d
;
; Proceedure examines grid of associated with axisID.  
; If grid values are discontinuous (at one point) and
; apparantly peroidic (as in Dorland's kx grid),
; Then both this grid and the associated values are shifted
; such that the grid vaues are monotonic.
;
; Written 6/7/00 by
;	W.M. Nevins
;
axis=axisID
IF(TypeOf(axis) EQ 7) THEN axis = Self -> AxisNumber(axisID)
IF(axis EQ 0) THEN RETURN
IF(axis GT 3) THEN BEGIN
	MESSAGE, 'axisID should be less that 4 to have gotten here!', Informational=d
	RETURN
ENDIF
axisStr = STRING(axis, FORMAT='(i1)')
commandStr = 'grid = self.grid' + axisStr
ok = EXECUTE(commandStr)
axisValues = *Grid.Values
nValues = N_ELEMENTS(axisValues)
delta = SHIFT(axisValues, 1) - axisValues			; Computes difference between nearest neighbors using periodic boundary conditions
deltaNorm = 0.5*TOTAL(ABS(delta))/(nValues-1)		; This should be (about) equal to the typical grid spacing	if we have a single discontinuity
deltaMax = MAX(delta, indexShift)				; 'indexShift' is index to last element of 'axis' before discontinuity.
IF(indexShift EQ nValues) THEN RETURN				; Discontinuity is at boundary (where GKV routines expect it)
axisValues = SHIFT(axisValues, -indexShift)		; Shift axisValues to put discontinuity at boundary.
PTR_FREE, Grid.Values
Grid.values  = PTR_NEW(axisValues)
axisMin=MIN(axisValues, MAX=axisMax)
Grid.range = [axisMin, axisMax]
Grid.irange= [0, nValues-1]
commandStr = 'self.grid' + axisStr + '= Grid'
ok = EXECUTE(commandStr)
;
; Now shift the values array
;
values = *self.Values
CASE axis OF
	1	:	values = SHIFT(values, -indexShift, 0, 0)
	2	:	values = SHIFT(values, 0, -indexShift, 0)
	3	:	values = SHIFT(values, 0, 0, -indexShift)
ENDCASE
PTR_FREE, Self.Values
self.Values = PTR_NEW(values)
RETURN 
END ; ****** GKVs3D::AxisShift ****** ;


PRO GKVs4D::AxisShift, axis=axisID, debug=d
;
; Proceedure examines grid of associated with axisID.  
; If grid values are discontinuous (at one point) and
; apparantly peroidic (as in Dorland's kx grid),
; Then both this grid and the associated values are shifted
; such that the grid vaues are monotonic.
;
; Written 6/7/00 by
;	W.M. Nevins
;
axis=axisID
IF(TypeOf(axis) EQ 7) THEN axis = Self -> AxisNumber(axisID)
IF(axis EQ 0) THEN RETURN
IF(axis GT 4) THEN BEGIN
	MESSAGE, 'axisID should be less that 5 to have gotten here!', Informational=d
	RETURN
ENDIF
axisStr = STRING(axis, FORMAT='(i1)')
commandStr = 'grid = self.grid' + axisStr
ok = EXECUTE(commandStr)
axisValues = *Grid.Values
nValues = N_ELEMENTS(axisValues)
delta = SHIFT(axisValues, 1) - axisValues			; Computes difference between nearest neighbors using periodic boundary conditions
deltaNorm = 0.5*TOTAL(ABS(delta))/(nValues-1)		; This should be (about) equal to the typical grid spacing	if we have a single discontinuity
deltaMax = MAX(delta, indexShift)				; 'indexShift' is index to last element of 'axis' before discontinuity.
IF(indexShift EQ nValues) THEN RETURN				; Discontinuity is at boundary (where GKV routines expect it)
axisValues = SHIFT(axisValues, -indexShift)		; Shift axisValues to put discontinuity at boundary.
PTR_FREE, Grid.Values
Grid.values  = PTR_NEW(axisValues)
axisMin=MIN(axisValues, MAX=axisMax)
Grid.range = [axisMin, axisMax]
Grid.irange= [0, nValues-1]
commandStr = 'self.grid' + axisStr + '= Grid'
ok = EXECUTE(commandStr)
;
; Now shift the values array
;
values = *self.Values
CASE axis OF
	1	:	values = SHIFT(values, -indexShift, 0, 0, 0)
	2	:	values = SHIFT(values, 0, -indexShift, 0, 0)
	3	:	values = SHIFT(values, 0, 0, -indexShift, 0)
	4	:	values = SHIFT(values, 0, 0, 0, -indexShift)
ENDCASE
PTR_FREE, Self.Values
self.Values = PTR_NEW(values)
RETURN 
END ; ****** GKVs4D::AxisShift ****** ;


FUNCTION NetCDF_Data, FileName=fileName, Path=path, ArrayStyle=array_Style, Threshold=thrshd, Species=species, RI=real_or_imaginary, Debug=d
;
; Reads standard netCDF file, and produces GKV objects
;
; KEYWORDS:
;
;	FileName		Name of netCDF file.  Must include either full path, or else the path from the 
;				current working directory.  If FileName is not provided, DIALOG_PICKFILE will 
;				be called to allow the user to select an appropriate netCDF file.
;
;	Path			Sets path to initial directory selected by DIALOG_PICKFILE (which will allow user
;				to change to other directories...).  If no path is provided, DIALOG_PICKFILE will
;				default to current working directory. 
;
;	ArrayStyle		Set to 'c' if netCDF file was written by code (like C or PASCAL) which stores array
;				data in 'column major' format.  Default is 'fortran', appropriate when the netCDF file
;				was written by code (like FORTRAN or IDL) which stores array data in 'row major' format.
;				This issue is discussed in Chapter 5 of 'Building IDL Applications'.
;				Defaults to 'fortran'.
;
;	Threshold		Size (in units of 10^6 elements) of threshold between 'large' and 'small' data sets.
;				User will be prompted before reading a data array of more than 'Threshold*10^6' elements
;				from the netCDF file.  Default is 1 (x10^6).  ***NOT YET IMPLIMENTED***
;
;	Species		A string containing the name (within the netCDF file) of the DIMENSION which indices
;				particle species.  NetCDF_Data will break separate data corresponding to different
;				species into separate GKV objects.  The species index (incremented by 1, so that the
;				first species will be '1', rather than '0') will be appended to the objects mnemonic
;				and added (as a subscript) to the objects title.
;				Defaults to 'species'.
;
;	RI			A string containing the name (within the netCDF file) of the DIMENSION which is used
;				to index the real and imaginary part of a complex variable.  Within IDL's 'row major'
;				format this MUST be the first (that is, most rapidally varying) dimension of the array.
;				***NOTE*** that NCDUMP (which is implimented as a C code) should show 'RI' as the 
;				LAST INDEX of the array within the netCDF file (consistent with C's use of 'column
;				major' format).
;
;	Debug			Set (i.e., put '/Debug' within the argument string) to turn off error trapping.  Turning
;				off error trapping is useful if you are debuggin NetCDF_Data.	
;
;
; Written by W.M. Nevins
;	6/7/00
;
deBug=0
IF(N_ELEMENTS(d) NE 0) THEN deBug=d
IF(deBug EQ 0) THEN BEGIN						; Set Keyword 'DeBug' to avoid error trap
	Catch, error 							; Set up error trap
	IF error NE 0 then Begin 
		Catch, /Cancel             			; Cancel error trap 
		ok=Error_Message(/Trace)  				; Print out a trace-back in the output log.
		IF(cdfID GE 0) THEN NCDF_CLOSE, cdfID		; Close any open NCDF files
		IF(N_ELEMENTS(gkvObjects) NE 0) THEN $
   			RETURN, gkvObjects				; Return SOMETHING useful if possible;
   		RETURN, 0                  			;	return 0 if all else fails. 
	ENDIF
ENDIF
;
; set cdfID to -1 to indicate that no netCDF file is currently open
;
cdfID=-1	

specIndx='species'
IF(KEYWORD_SET(species)) THEN specIndx=STRLOWCASE(species)

ri='ri'
IF(KEYWORD_SET(real_or_imaginary)) THEN ri = real_or_imaginary

arrayStyle='fortran'
IF(KEYWORD_SET(array_Style)) THEN arrayStyle=STRLOWCASE(array_Style)

threshold=1000000L
IF(N_ELEMENTS(thrshd) NE 0) THEN threshold=thrshd*threshold
;
; find netCDF file if needed
;
IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.nc') ELSE fineIn=fileName
ok = FINDFILE(fileIn, Count=nfiles)			; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN				; 	we didn't find a file... 
	fileIn=DIALOG_PICKFILE(Path=path, /MustExist, Filter='*.nc')
ENDIF		
cdfID = NCDF_OPEN(fileIn, /NoWrite)			; Open netCDF file (sets cdfID to an integer „ 0)
cdfInfo = NCDF_INQUIRE(cdfID)				;	and get info regarding contents

nVars = cdfInfo.Nvars						; Number of variables  in this netCDF file
nDIms =cdfInfo.Ndims						; Number of dimensions in this netCDF file
nGatts=cdfInfo.nGatts						; Nubmer of global attributes in this netCDF file
varType = STRARR(nVars)
varSave = BYTARR(nVars)
multiSpecies = BYTARR(nVars)
dimName = STRARR(nDims)
dimSize = LONARR(nDims)
dimVars = REPLICATE(-1L, nDims)
dimGrid = REPLICATE({GRID}, nDims)

;
; Get Global attributes, and set up default GKVsd structure
;
defaultStructure = {GKVsd}
headerElements = N_TAGS({GKVsd})
FOR gattID = 0, nGAtts-1 DO BEGIN
	gattName = NCDF_ATTNAME(cdfID, /GLOBAL, gattID)	; Get name of the gattID^th Global attribute
	lgattName=STRLOWCASE(gattName)				; 	and coherce it to lower case
	NCDF_ATTGET, cdfID, /GLOBAL, gattName, gattValue; Get attribute value
	IF (TypeOf(gattValue) NE 7) THEN gattValue = STRING(gattValue)
	CASE lgattName OF
		'codename'		:	defaultStructure.CodeName	= gattValue
		'title'		:	defaultStructure.CodeName	= gattValue
		'codepi'		:	defaultStructure.CodePI	= gattValue
		'conventions'	:	defaultStructure.CodePI	= gattValue
		'runid'		:	defaultStructure.RunID	= gattValue
		'fileid'		:	defaultStructure.FileID	= gattValue
		ELSE			:
	ENDCASE
ENDFOR

nspec=0
FOR dimID = 0, nDims-1 DO BEGIN
	NCDF_DIMINQ, cdfID, dimID, name, dsize		; Get info about dimensions
	dimName[dimID]=name
	dimSize[dimID]=dsize	
	IF(name EQ specIndx) then nspec=dsize
ENDFOR
FOR dimID = 0, nDims-1 DO BEGIN
;
; Set default {GRID} info
;
	gridStructure = {Grid}								
	gridStructure.Mnemonic	= STRCOMPRESS(dimName[dimID], /REMOVE_ALL)	; Remove any embedded blanks from the grid mnemonic built from 'dimName' 
	gridStructure.Title		= dimName[dimID]
	gridStructure.Values		= PTR_NEW(FINDGEN(dimSize[dimID]))
	gridStructure.Range		= [0., dimSize[dimID]-1.]
	gridStructure.irange		= [0, dimSize[dimID]-1]
;
; Check for corresponding variable	
;
	varID = NCDF_VARID(cdfID, dimName[dimID])					
	dimVars[dimID] = varID

	IF(varID NE -1) THEN BEGIN						
;
; We have found a variable whose name corresponds to that of the dimID^th dimension
;							
		varInfo = NCDF_VARINQ(cdfID, varID)					; Get info about dimension variable
		IF( (N_ELEMENTS(varInfo.ndims) EQ 1) AND (varInfo.Dim[0] EQ dimID) ) THEN BEGIN
;
; ... and the dimensions match.  
;
			varType[varID]='dimension'
			NCDF_VARGET, cdfID, varID, values
			numValues = N_ELEMENTS(values)
			IF(numValues EQ dimSize[dimID]) THEN BEGIN			; Make sure dimensions match
				PTR_FREE, gridStructure.Values				; ... then replace default gridStructure.values with
				gridStructure.Values = PTR_NEW(values)		; 	 the values of this variable.
				vmin=values[0]
				vmax=values[numValues-1]
				gridStructure.Range = [vmin, vmax]
			ENDIF
			nAtts = varInfo.nAtts
			FOR attID = 0, nAtts-1 DO BEGIN					; Look for Mnemonic, title, units, etc.
				attName = NCDF_ATTNAME(cdfID, varID, attID)	;	Get name of attID^th attribute
				lattName = STRLOWCASE(attName)				;	coherce attname to lower case
				lattName = STRCOMPRESS(lattName, /REMOVE_ALL)	; 	and remove any embedded blanks
				NCDF_ATTGET, cdfID, varId, attName, attValue	; Get value of attribute
				sattValue = attValue
				IF(TypeOF(sattValue) NE 7) THEN sattValue = STRING(attValue)
				CASE lattName OF
					'mnemonic'		:  gridStructure.mnemonic	= STRCOMPRESS(sattValue, /REMOVE_ALL)
					'long_name'	:  gridStructure.Title	= sattValue + '!N!X'	; The '!N!X' insures that any font/level changes are local to this text
					'pretty_name'	:  gridStructure.Title	= sattValue + '!N!X'
					'idl_name'		:  gridStructure.Title	= sattValue + '!N!X'
					'units'		:  gridStructure.Units	= sattValue + '!N!X'
					'boundary'		:  gridStructure.Boundary	= sattValue
					ELSE			:
				ENDCASE
			; 
			; End of search for attributes of this dimensional variable 
			;
			ENDFOR
		ENDIF
	;
	; We are finished with this dimensional variable
	;	 
	ENDIF
dimGrid[dimID] = gridStructure
;
; End of search for dimensional variables
;
ENDFOR
;
; Check for 0-D variables, short varibles, and long variables.  Also check for multiSpecies variables (that is, variables dimensioned with specIndx)
;
varNames = STRARR(nVars)
longVars = 0
longVarIndices = 0
FOR varID=0, nVars-1 DO BEGIN
	IF(varType[varID] EQ 'dimension') THEN GOTO, DONE1
	varInfo = NCDF_VARINQ(cdfID, varID)
	varNames[varID] = varInfo.Name
	IF(varinfo.nDims EQ 0) THEN BEGIN
		varType[varID]='0-D'
		GOTO, DONE1
	ENDIF ELSE BEGIN
		varSize=1
		FOR dimID=0, varinfo.nDims-1 DO varSize=varSize*dimSize[varinfo.Dim[dimID]]
		Case (varSize GT threshold) OF
			0	:	BEGIN
					varType[varID] = 'short'
					varSave[varID] = 1
					END
			1	:	BEGIN
					varType[varID] = 'long'
					varSave[varID] = 0
					longVarIndices = [longVarIndices, varID]
					longVars = longVars + 1
					END
		ENDCASE
	ENDELSE
	IF((varinfo.nDims EQ 1) AND (dimName(varInfo.Dim[0]) EQ specIndx)) THEN BEGIN
		varType[varID]='0-D'
	ENDIF ELSE BEGIN
		nVarDims=varInfo.nDims
		FOR dimID=0,nVarDims-1 DO IF(dimName[Varinfo.Dim[dimID]] EQ specIndx) THEN multiSpecies[varID]=1
	ENDELSE
	IF(STRCMP(varInfo.DataType, 'CHAR', /FOLD_CASE)) THEN varType[varID] = 'char'
DONE1	:
ENDFOR
;
; Query user about saving 'long' variables
; 
IF(longVars GT 0) THEN BEGIN
	longVarNames = varNames[longVarIndices[1:longVars]]
	checkedLongVars = XGKV_CHECK_NAMES(Names=longVarNames, Window_Title='GKV::Select Data')
	varSave[longVarIndices[1:longVars]] = checkedLongVars
ENDIF
;
; Find number of GKV objects to be created
;	(one object for each non-dimensional variable)
;
nObjs = 0
FOR varID=0, nVars-1 DO IF((varType[varID] NE 'dimension') AND (varType[varID] NE '0-D')) THEN nObjs = nObjs + 1
nObjs = nObjs + nspec*TOTAL(multiSpecies)
;
; Create Object array for GKV objects to be created
;
gkvObjects = OBJARR(nObjs)

iObjs = 0
FOR varID=0, nVars-1 DO BEGIN
	IF(varSave[varID] EQ 0) THEN GOTO, DONE2
	IF((varType[varID] NE 'dimension') AND (varType[varID] NE '0-D') AND (varType[varID] NE 'char')) THEN BEGIN
		varInfo = NCDF_VARINQ(cdfID, varID)		; Get info about variables
		varName = varInfo.Name
		varDims = varInfo.nDims
		varDimID= varInfo.Dim
		nAtts   = varInfo.nAtts
;
; Check for complex variables (from Bill Dorland, this occurs if first dimension is named 'ri')
;
		complex = 0
		firstDimName = dimName[varDimID[0]]
		firstDimName = STRLOWCASE(firstDimName)
		IF(firstDimName EQ ri) THEN BEGIN
			complex = 1
			varDims = varDims-1						; *** NOTE *** varDims ALREADY takes account of 'complex'
		ENDIF
;
; Check for 'multispecies' variable
;
		vVarDims = varDims
		IF(multiSpecies[varID]) THEN vVarDims = varDims-1	; *** while vVarDims takes account of BOTH 'complex' and 'multiSpecies'
;
; Create GKV structure for this variable
;
		CASE vVarDims OF
			0	: objStructure = {GKVsd}
			1	: objStructure = {GKVs1D}
			2	: objStructure = {GKVs2D}
			3	: objStructure = {GKVs3D}
			4	: objStructure = {GKVs4D}
			ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
		ENDCASE
;
; Copy default information into structure
;
		FOR i=0, headerElements-1 DO objStructure.(i) = defaultStructure.(i)
		objStructure.mnemonic = STRCOMPRESS(varName, /REMOVE_ALL)
		objStructure.Title = varName
		objStructure.Indices = PTR_NEW(REPLICATE('*', vVarDims))		
;
; Check variable attributes
;
		FOR attID=0, nAtts-1 DO BEGIN
			attName = NCDF_ATTNAME(cdfID, varID, attID)	; Get name of attID^th attribute
			lattName = STRLOWCASE(attName)				;	coherce attname to lower case
			lattName=STRCOMPRESS(lattName, /REMOVE_ALL)	;	and remove any embedded blanks
			NCDF_ATTGET, cdfID, varId, attName, attValue	; Get value of attribute
			sattValue = attValue
			IF(TypeOF(sattValue) NE 7) THEN sattValue = STRING(attValue)
			CASE lattName OF
				'mnemonic'		:  gridStructure.mnemonic	= STRCOMPRESS(sattValue, /REMOVE_ALL)
				'long_name'	:  objStructure.Title		= sattValue
				'pretty_name'	:  objStructure.Title 	= sattValue + '!N!X'	; The '!N!X' insures that any font/level changes are local to this text
				'idl_name'		:  objStructure.Title 	= sattValue + '!N!X'
				'units'		:  objStructure.Units 	= sattValue + '!N!X'
				'valid_range'	:  objStructure.vrange	= attValue
				ELSE			:
			ENDCASE
		ENDFOR
		IF( objStructure.Title EQ varName) THEN objStructure.Title = objStructure.Mnemonic	
;
; Get values of this variable
;
		NCDF_VARGET, cdfID, varID, temp
;
; IF values are complex, then need to repack into a complex variable
;
		IF(complex) THEN BEGIN
			CASE varDims OF
				1	:	values = REFORM(COMPLEX(temp[0,*], temp[1,*]), dimSize[varDimID[1]])
				2	:	values = REFORM(COMPLEX(temp[0,*,*], temp[1,*,*]), dimSize[varDimID[1]], dimSize[varDimID[2]])
				3	:	values = REFORM(COMPLEX(temp[0,*,*,*], temp[1,*,*,*]), dimSize[varDimID[1]], dimSize[varDimID[2]], dimSize[varDimID[3]])
				4	:	values = REFORM(COMPLEX(temp[0,*,*,*,*], temp[1,*,*,*,*]), dimSize[varDimID[1]], dimSize[varDimID[2]], dimSize[varDimID[3]], dimSize[varDimID[4]])
				ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
			ENDCASE
		ENDIF ELSE BEGIN
			values = temp
		ENDELSE
		
		IF(arrayStyle EQ 'c') THEN BEGIN			; If data was written by a 'c' (or pascal) code,
			CASE varDims OF						; then the 'values' array must be reformed to 
				1	: 						; to reflect IDL's array storage conventions
				2	: values = REFORM(values, dimSize[varDimID[1]], dimSize[varDimID[0]])
				3	: values = REFORM(values, dimSize[varDimID[2]], dimSize[varDimID[1]], dimSize[varDimID[0]])
				4	: values = REFORM(values, dimSize[varDimID[3]], dimSize[varDimID[2]], dimSize[varDimID[1]], dimSize[varDimID[0]])
				ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
			ENDCASE
		ENDIF
;
; Insert values into objStructure ***NOTE--THERE IS AN ERROR HERE when complex=1!
;
		IF(multiSpecies[varID]) THEN BEGIN
;
; First, find 'nspecLoc', the location of 'specIndx' within this variable (note that nspecLoc is an element of [1, VarDims], NOT [0, VarDims-1], or [1, vVarDims])
;
			FOR dimID=0,varDims-1+complex DO IF(STRCMP(dimName[varDimID[dimID]], specIndx)) THEN nspecLoc = dimID+1
			cNspecLoc = nspecLoc - complex  ; 'cNspecLoc' points to index in (possibly complex) 'values' array corresponding to 'specIndx'
			CASE vVarDims OF
				1	:	BEGIN
					CASE cNspecLoc OF
						1	:	temp = REFORM(values[0,*], dimSize[varDimID[1+complex]])
						2	:	temp = REFORM(values[*,0], dimSize[varDimID[0+complex]])
					ENDCASE
						END
				2	:	BEGIN
					CASE cNspecLoc OF
						1	:	temp = REFORM(values[0,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
						2	:	temp = REFORM(values[*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]])
						3	:	temp = REFORM(values[*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]])
					ENDCASE
						END
				3	:	BEGIN
					CASE cNspecLoc OF
						1	:	temp = REFORM(values[0,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
						2	:	temp = REFORM(values[*,0,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
						3	:	temp = REFORM(values[*,*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]])
						4	:	temp = REFORM(values[*,*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
					ENDCASE
						END
				4	:	BEGIN
					CASE cNspecLoc OF
						1	:	temp = REFORM(values[0,*,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
						2	:	temp = REFORM(values[*,0,*,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
						3	:	temp = REFORM(values[*,*,0,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
						4	:	temp = REFORM(values[*,*,*,0,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[4+complex]])
						5	:	temp = REFORM(values[*,*,*,*,0], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
					ENDCASE
						END
				ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
			ENDCASE		
			objStructure.values = PTR_NEW(temp)
		ENDIF ELSE BEGIN
			objStructure.values = PTR_NEW(values)
		ENDELSE
;
; Add appropriate GRID structures to objStructure
;
		IF(arrayStyle EQ 'c') THEN BEGIN	
			i=VarDims+1						; If data was written by a 'c' (or pascal) code,
			FOR j=1, vVarDims DO BEGIN				;	then the the indices IDL has first correspond to
				i=i-1							;	to those with the highest 'vardim' index.
				IF(multiSpecies[varID]) THEN IF(i EQ nspecLoc) THEN i=i-1
				CASE j OF
					1	: objStructure.Grid1 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					2	: objStructure.Grid2 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					3	: objStructure.Grid3 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					4	: objStructure.Grid4 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
				ENDCASE
			ENDFOR
		ENDIF ELSE	BEGIN
			i=complex
			FOR j=1, vVarDims DO BEGIN				
				i=i+1
				IF(multiSpecies[varID]) THEN IF(i EQ nspecLoc) THEN i=i+1
				CASE j OF
					1	: objStructure.Grid1 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					2	: objStructure.Grid2 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					3	: objStructure.Grid3 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
					4	: objStructure.Grid4 = GKVsd_GridCopy(dimGrid[varDimID[i-1]])
				ENDCASE
			ENDFOR
		ENDELSE
;
; Make appropriate GKVsd object, and store in gkvObjects array
;
		CASE vVarDims OF
			0	: gkvObjects[iObjs] = OBJ_NEW('GKVsd' , objStructure)
			1	: gkvObjects[iObjs] = OBJ_NEW('GKVs1D', objStructure)
			2	: gkvObjects[iObjs] = OBJ_NEW('GKVs2D', objStructure)
			3	: gkvObjects[iObjs] = OBJ_NEW('GKVs3D', objStructure)
			4	: gkvObjects[iObjs] = OBJ_NEW('GKVs4D', objStructure)
		ENDCASE
;
; If this is a multispecies variable, then several objects must be created
;
		IF(multiSpecies[varID]) THEN BEGIN
				specObjs = OBJARR(nspec)
				speciesStr = STRING(1, FORMAT='(I1)')
				mnemonic = objStructure.mnemonic + '_' + speciesStr
				title = objStructure.title + '!Ispec ' + speciesStr  + '!N'
				gkvObjects[iObjs] -> Set, Mnemonic=mnemonic, Title=title 
				specObjs[0] = gkvObjects[iObjs]  
			FOR ispecies=1,nspec-1 DO BEGIN
				speciesStr = STRING(ispecies+1, FORMAT='(I1)')
				mnemonic = objStructure.mnemonic + '_' + speciesStr
				title = objStructure.title + '!Ispec ' + speciesStr + '!N'
				specObjs[ispecies] = gkvObjects[iObjs] -> MakeCopy(/NoValues)
				CASE vVarDims OF
					1	:	BEGIN
						CASE cNspecLoc OF
							1	:	temp = REFORM(values[ispecies,*], dimSize[varDimID[1+complex]])
							2	:	temp = REFORM(values[*,ispecies], dimSize[varDimID[0+complex]])
						ENDCASE
							END
					2	:	BEGIN
						CASE cNspecLoc OF
							1	:	temp = REFORM(values[ispecies,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
							2	:	temp = REFORM(values[*,ispecies,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]])
							3	:	temp = REFORM(values[*,*,ispecies], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]])
						ENDCASE
							END
					3	:	BEGIN
						CASE cNspecLoc OF
							1	:	temp = REFORM(values[ispecies,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
							2	:	temp = REFORM(values[*,ispecies,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
							3	:	temp = REFORM(values[*,*,ispecies,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]])
							4	:	temp = REFORM(values[*,*,*,ispecies], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]])
						ENDCASE
							END
					4	:	BEGIN
						CASE cNspecLoc OF
							1	:	temp = REFORM(values[ispecies,*,*,*,*], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
							2	:	temp = REFORM(values[*,ispecies,*,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
							3	:	temp = REFORM(values[*,*,ispecies,*,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[3+complex]], dimSize[varDimID[4+complex]])
							4	:	temp = REFORM(values[*,*,*,ispecies,*], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[4+complex]])
							5	:	temp = REFORM(values[*,*,*,*,ispecies], dimSize[varDimID[0+complex]], dimSize[varDimID[1+complex]], dimSize[varDimID[2+complex]], dimSize[varDimID[3+complex]])
							ENDCASE
								END
					ELSE	: ; *** NEED SOME ERROR PROCESSING HERE ***
				ENDCASE		
				specObjs[ispecies] -> Set, Mnemonic=mnemonic, Title=title, values=PTR_NEW(temp)						
			ENDFOR
;
; Need to check for non-monotonic 'kx' axis. 
;
			FOR ispecies=0, nspec-1 DO specObjs[ispecies] -> AxisShift, axis='kx'
;
; Load finished GKV objects into gkvObjects array
;
			FOR ispecies=0, nspec-1 DO gkvObjects[iObjs+ispecies] = specObjs[ispecies]
			iObjs = iObjs + nspec
		ENDIF ELSE BEGIN
			gkvObjects[iObjs] -> AxisShift, axis='kx'
			iObjs = iObjs+1
		ENDELSE
	;
	; Finished with this variable
	;
	ENDIF
DONE2:	; ...
;
; finished with loop over variables
;
ENDFOR
;
; We've finished reading the data.  Now close the NETCDF file.
;
NCDF_CLOSE, cdfID
cdfID=-1
;
; Take the OBJARR, gkvObjects, and form an anomymous structure.
;	(use objects's Mnemonic as the tagname
;
IF(iObjs EQ 0) THEN RETURN, 0
FOR i=0,iObjs-1 DO BEGIN
	gkvObjects[i] -> Get, mnemonic=tagName
	tagName = STRCOMPRESS(tagName, /REMOVE_ALL)
	CASE i OF
		0:	result = CREATE_STRUCT(tagName, gkvObjects[0])
		ELSE:	result = CREATE_STRUCT(result, tagName, gkvObjects[i])
	ENDCASE
ENDFOR

RETURN, result

END ; ****** NetCDF_data ****** ;


