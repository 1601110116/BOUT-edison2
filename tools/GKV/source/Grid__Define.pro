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

FUNCTION GKV_GridSame, Grid1, Grid2, All=all, Force=force
;
; Check for same grid structure
; (should be put into "GRID__DEFINE")
;
type1 = TypeOf(Grid1)
type2 = TypeOF(Grid2)
IF(type1 NE 8) THEN RETURN, 0					; argument is not a structure
IF(type2 NE 8) THEN RETURN, 0
IF(Grid1.mnemonic NE Grid2.mnemonic)	THEN RETURN, 0	; Grid mnemonics differ
IF(Grid1.units NE Grid2.units) 		THEN RETURN, 0	; Grid units differ
;IF(Grid1.boundary NE Grid2.boundary)	THEN RETURN, 0	; Grid boundary conditions differ
;IF(Grid1.uniform NE Grid2.uniform)	THEN RETURN, 0	; Grid uniformity differs
min1 = Grid1.irange[0]
min2 = Grid2.irange[0]
max1 = Grid1.irange[1]
max2 = Grid2.irange[1]
values1 = *Grid1.values
values2 = *Grid2.values
IF N_ELEMENTS(all) THEN BEGIN
	min1=0
	min2=0
	max1=N_ELEMENTS(values1) - 1
	max2=N_ELEMENTS(values2) - 1
ENDIF
iforce = 0
IF((max1-min1) NE (max2-min2)) THEN BEGIN
	IF(NOT KEYWORD_SET(Force)) THEN RETURN, 0		; Try to match Grid2 range to Grid1 range
	IF(KEYWORD_SET(All)) THEN RETURN, 0			; Can't force match if must use all grid points
	xmin1 = values1(min1)
	xmax1 = values1(max1)
	err = MIN((values2-xmin1)^2, min2)			; Set min2 to get best match to values1(min1)
	err = MIN((values2-xmax1)^2, max2)			; Set max2 to get best match to values1(max1)
	IF((max1-min1) NE (max2-min2)) THEN RETURN, 0	; Same number of grid points?
	iforce = 1
ENDIF
Err = TOTAL((values1[min1:max1]^2 - values2[min2:max2]^2))/(1+max1-min1)
L = values1[max1] - values1[min1]
dL=L/(max1-min1)
IF(ERR GT 1.e-3*(dL^2)) THEN RETURN, 0
IF(iforce) THEN Grid2.irange = [min2, max2]
RETURN, 1
END ; ****** GKV_GridSame ****** ;


PRO GKVsd_SetGrid, arg, GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
				GridUnits=GridUnits,	GridValues=GridValues, 			$
				boundary=boundary, uniform=uniform,range=range, 		$
				irange=irange, _Extra=extra
;
; Proceedure to set values of elements of a realization of the 'Grid' class
;
argClass = TAG_NAMES(arg, /Structure_Name)
argClass = STRUPCASE(argClass)
IF(argClass NE 'GRID') THEN BEGIN
	MESSAGE, 'GKVsd_SetGrid called with bad argument', /INFORMATIONAL
	RETURN
ENDIF

	IF(N_ELEMENTS(GridMnemonic)) 	THEN	arg.mnemonic	= GridMnemonic
	IF(N_ELEMENTS(GridTitle))		THEN	arg.title		= GridTitle
	IF(N_ELEMENTS(GridUnits))		THEN	arg.units		= GridUnits
; 29-Aug-02 DRM: Bug fix; add "gt 0" to the test of N_ELEMENTS(GridValues):
	IF(N_ELEMENTS(GridValues) GT 0) THEN BEGIN
		vtype = TypeOf(Gridvalues)
		IF(vtype NE 10) THEN Gridvalues = PTR_NEW(Gridvalues)
		PTR_FREE, arg.values
		arg.values	 = Gridvalues
	ENDIF
	IF(N_ELEMENTS(boundary))	THEN	arg.boundary	= boundary
	IF(N_ELEMENTS(uniform))	THEN	arg.uniform	= uniform
	IF(N_ELEMENTS(range) EQ 2) THEN	  arg.range	= range
	IF(N_ELEMENTS(irange) EQ 2) THEN arg.irange	= irange
RETURN
END ; ****** GKVsd_SetGrid ****** ;
	

PRO GKVsd_GetGrid, arg, GridMnemonic = GridMnemonic, GridTitle=GridTitle, 		$
				GridUnits=GridUnits,	GridValues=GridValues, 			$
				boundary=boundary, uniform=uniform,range=range, 		$
				irange=irange, _Extra=extra
;
; Proceedure to get values of elements of a realization of the 'Grid' class
;
argClass = TAG_NAMES(arg, /Structure_Name)
argClass = STRUPCASE(argClass)
IF(argClass NE 'GRID') THEN BEGIN
	MESSAGE, 'GKVsd_GetGrid called with bad argument', /INFORMATIONAL
	RETURN
ENDIF
GridMnemonic	= arg.mnemonic
Gridtitle		= arg.title
Gridunits		= arg.units
Gridvalues		= arg.values
boundary 		= arg.boundary
uniform 		= arg.uniform
range 		= arg.range
irange 		= arg.irange
RETURN
END ; ***** GKVsd_GetGrid ***** ;


PRO GKVsd_PrintGrid, arg
;
; Prints contents of 'Grid' class
;
argClass = TAG_NAMES(arg, /Structure_Name)
argClass = STRUPCASE(argClass)
IF(argClass NE 'GRID') THEN BEGIN
	MESSAGE, 'GKVsd_PrintGrid called with bad argument', /INFORMATIONAL
	RETURN
ENDIF

PRINT, ' Mnemonic = ', arg.mnemonic
PRINT, '    title = ', arg.title
PRINT, '    units = ', arg.units
IF(PTR_VALID(arg.values)) THEN BEGIN
	PRINT, 'valueInfo = ', SIZE(*arg.values)
ENDIF ELSE BEGIN
	PRINT, 'valueInfo = Null Pointer'
ENDELSE
PRINT, ' boundary = ', arg.boundary
PRINT, '  uniform = ', arg.uniform
PRINT, '    range = ', arg.range
PRINT, '   irange = ', arg.irange
RETURN
END ; ****** GKVsd_PrintGrid ****** ;


FUNCTION GKVsd_GridCopy, Grid, Irange=irange
;
; Make 'deep' copy of argument
;
; If the keyword Irange is set, then restrict copy to specified range
;
argClass = TAG_NAMES(Grid, /Structure_Name)
argClass = STRUPCASE(argClass)
IF(argClass NE 'GRID') THEN BEGIN
	MESSAGE, 'GKVsd_GridCopy called with bad argument', /INFORMATIONAL
	RETURN, 0
ENDIF
result = {Grid}
FOR i=0, N_TAGS({Grid})-1 DO result.(i)=Grid.(i)
IF(N_ELEMENTS(irange) EQ 2) THEN BEGIN
	iirange = irange
	imin=iirange[0]
	imax=iirange[1]
	gridValues = (*Grid.values)[imin:imax]
	result.values = PTR_NEW(gridValues)
	result.irange = [0, imax-imin]
	result.range[0] = result.range[0] > gridValues[0]
	result.range[1] = result.range[1] < gridValues[imax-imin]
	RETURN, result
ENDIF
gridValues = *Grid.values
result.values = PTR_NEW(gridValues)
RETURN, result
END ; ****** GKVsd_GridCopy ****** ;


PRO Grid__Define

Grid = {	Grid,					$	; Defining 'Grid' class
		Mnemonic	: '',			$	; ascii Mnemonic 
		title		: '',			$	; "Pretty" Axis title		(use IDL Vector fonts)
		units		: '',			$	; "Pretty" Units		(use IDL Vector fonts)
		values		: PTR_NEW(),		$	; Pointer to an array containing grid values
		boundary	: "",			$	; Boundary conditions in co-ordinate  
								;	(Open, Periodic, Neumann, Dirchelet)  
								;	Separate with "/"for mixed boundary conditions.
		uniform		: 0B,			$	; Set =1 for uniform grid spacing, 0 otherwise			
		range		: FLTARR(2),		$	; Range of grid-values within the plot-window
								;	(in general, this differs from entire data string)
		irange		: LONARR(2)		}	; Range of indices within signal window 
								; 	(in general, this differs from entire data string)

END ; ****** Grid__Define ****** ;
