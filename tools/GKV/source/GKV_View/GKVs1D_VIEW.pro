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

PRO GKVs1D::View, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, _Extra=extra
;
; Purpose:
;
;		This proceedure is part of a 'wrapper' which allows me
;		to use D. Fanning's 'XWINDOW' widget to display GKV
;		objects.  It works for GKVsd Objects of dimensionality
;		1 through 4.
;
; Arguments:
;
;		arg1 through arg10 are optional arguments which will be
;		ploted over 'self'.  Implimented for GKVs1D objects only. 
;		All arguments should be GKVs1D objects. 
;		(Optional)
;
; Keywords:
;
;		Keywords native to this wrapper are:
;
;			Shade_Surf		Set this keyword (i.e., put "/Shade_Surf" on the
;						command line) to obtain a surface plot for objects
;						of 2 (or more) dimensions (default is an 'Image" 
;						plot). (Optional)
;
;			Reverse		Set this keyword (i.e., put "/Reverse" on the
;						command line) to switch the forground and 
;						background colors (see comments at begining
;						of COLOR_SETUP).  (Optional)
;
;		Remaining Keywords are passed first to XWINDOW (see XWINDOW source
;		for a descreption of XWINDOW keywords) and then to the object's
;		drawing routine (see GKVs1D::Draw, GKVs2D::Draw, or GKVs2D::Shade_Surf).
;		In order to preserve the GKVsD objects' drawing routine's runtime 
;		Runtime keywords, it is necessary to spell out all keywords in full
;		(i.e., keyword abreviation will not work!).
;
; Written by W.M. Nevins
;	7/16/00
;
; Save input color information
;
TvLct, rIn, gIn, bIn, /GET
ncolorsIn = !D.TABLE_SIZE
;
; Set bottom 11 colors in color table
;
;result = GetKeyWord('Reverse', extra)
;IF(KEYWORD_SET(result)) THEN BEGIN
;	Color_SetUp, 1, /REVERSE
;ENDIF ELSE BEGIN
;	Color_SetUp, 1
;ENDELSE
bottom = !COLOR_SETUP_NCOLORS  + 1
ncolors = !D.TABLE_SIZE - bottom
;
; Next task is to check for XWINDOW keywords 
; (This is required in order to suppress Keyword abreviation).
;
; Get first few values as needed to set defaults...
;
xsize = 400
result = GetKeyWord('WXSIZE', extra)
IF(result NE 'undefined') THEN xsize = result
ysize = 400
IF(result NE 'undefined') THEN ysize = result
DEVICE, GET_SCREEN_SIZE=screenSize
wxpos = (screenSize(0) - xsize) / 2.
wypos = (screenSize(1) - ysize) / 2.
;
; Reserve some output keywords
;
wid_Out = 0
drawID_Out = 0
top_out = 0
;
; Set up a structure containing the XWINDOW keywords, and set default values
;
;				KEYWORD			VALUE					INDEX
XW_KeyWords = {	$ ;	GROUP_LEADER	:	group,	 	$	;	
				WXSIZE		:	xsize,		$	;	0
				WYSIZE		:	ysize,		$	;	1
				XCOLORS		:	[ncolors,bottom],$	;	2
				WTITLE		:	'GKV VIEW',	$	;	3
				JUST_REGISTER	:	0,			$	;	4
				OUTPUT		:	1,			$	;	5
				NO_CHANGE_CONFIG:	0,			$	;	6
				ERASE			:	0,			$	;	7
				PROTECT		:	0,			$	;	8
				NOMENU		:	0,			$	;	9
				CPMENU		:	1,			$	;	10
				WXPOS			:	wxpos,		$	;	11
				WYPOS			:	wypos,		$	;	12
				BACKGROUND		:	!P.BACKGROUND	}	;	13	*************************************************
;				WID			:	wid_Out,		$	;	14	*	I don't yet have a mechanism for passing	*
;				DRAWID		:	drawID_Out,	$	;	15	*	back output keywords from XWINDOW while 	*
;				TOP			:	top_out		}	;	16	*	still suppressing keyword abreviation.	*
;														*************************************************
; Get values of any XWINDOW keywords 	which were set on command line from the 'extra' structure
; (note that GetKeyWord will remove these keywords from 'extra' in addition to returning the value)
;
; ntags = N_TAGS(XW_KeyWords) 
tagNames = TAG_NAMES(XW_KeyWords)	
FOR i=0,13 DO BEGIN
	result = GetKeyWord(tagNames[i], extra)
	IF(result NE 'undefined') THEN XW_KeyWords.(i) = result
ENDFOR
;
; Now combine XW_KeyWords structure with what is left of extra
;
IF(TypeOF(extra) EQ 8) THEN BEGIN
	extra = CREATE_STRUCT(XW_KeyWords, extra) 	
ENDIF ELSE BEGIN
	extra = XW_KeyWords
ENDELSE		
;
; Get class of GKVsd object
;
class = OBJ_CLASS(self)
;
; Concatenate arguments into 'oPlots' array
;
IF( (class EQ 'GKVS1D') OR (class EQ 'GKVS2D') ) THEN BEGIN
	nArgs = N_PARAMS()					; Find number of arguments
	args = STRARR(11)
	args = ["arg0", "arg1", "arg2", "arg3", "arg4", "arg5", "arg6", "arg7", "arg8", "arg9", "arg10"]
	IF(nArgs GT 0) THEN BEGIN
		nElements   = INTARR(11)
		totElements = INTARR(12)
		FOR i=0,10 DO BEGIN				; count number of valid GKVs1D objects passed as arguments
			nELEMENTS[i]=0
			commandString = "argType = TypeOf(" + args[i] + ")"
			ok = EXECUTE(commandString)
			commandString = "nElements[i] = TOTAL( OBJ_ISA(" + args[i] + ", 'GKVs1D') )"
			IF(argType EQ 11) THEN ok = EXECUTE(commandString)
		ENDFOR
		totElements[1:11] = TOTAL(nElements, /CUMULATIVE)
		nObjs = TOTAL(nElements)
		IF(nObjs GT 0) THEN oPlots = OBJARR(nObjs)
		
		FOR i=0,10 DO BEGIN				; Load valid GKVs1D objects into 'oPlots' array
			IF(nElements[i] NE 0) THEN BEGIN
				commandString = "addresses = WHERE( OBJ_ISA(" + args[i] + ", 'GKVs1D') )"
				ok = EXECUTE(commandString)
				commandString = "oPlots[totElements[i]:(totElements[i+1]-1)] = " + args[i] + "[addresses]"
				ok = EXECUTE(commandString)
			ENDIF
		ENDFOR
		
		IF(nObjs GT 0) THEN $	
			extra = CREATE_STRUCT(extra, 'oPlot', oPlots)	; and, then put 'oPlots' array into 'extra'
	ENDIF
;
; Now call XWINDOW with appropriate arguments
;
ENDIF
IF(class EQ 'GKVS1D') THEN BEGIN
	XWINDOW, 'PlotObj1D', self, _EXTRA=extra
ENDIF ELSE BEGIN
	XWINDOW, 'PlotObj2D', self, _EXTRA=extra
ENDELSE
;
;
RETURN
END	 ; ****** GKVs1D::View ****** ;


PRO PlotObj1D, Obj, _EXTRA=extra
;
;		This proceedure is part of a 'wrapper' which allows me
;		to use D. Fanning's 'XWINDOW' widget to display GKV
;		objects.
;
; Keywords:
;
;		Keywords native to this proceedure are:
;
;			oPlot			For 1-D Set this keyword equal a GKVs1D object,
;						or an array of GKVs1D objects, which you wish to 
;						Plot over original object.
;
;		Remaining Keywords are passed first to XWINDOW (see XWINDOW source
;		for a descreption of XWINDOW keywords) and then to the object's
;		drawing routine (see GKVs1D::Draw, GKVs2D::Draw, or GKVs2D::Shade_Surf).
;		In order to preserve the GKVsD objects' drawing routine's runtime 
;		Runtime keywords, it is necessary to spell out all keywords in full
;		(i.e., keyword abreviation will not work!).
;
; Written by W.M. Nevins
;	7/16/00
;
; Avoid changing contents of 'extra'
;
localExtra = -1
IF(TypeOF(extra) EQ 8) THEN localExtra = extra
;
; First task is to check for the local keyword, oplot.
;
result = GetKeyWord('oPlot', localExtra)
IF(TypeOf(result) EQ 11) THEN oPlot = result
;
; Draw Obj in currently open window
;
IF(typeOf(localExtra) EQ 8) THEN BEGIN
	Obj -> Draw, _EXTRA=localExtra
ENDIF ELSE BEGIN
	Obj -> DRAW
ENDELSE
;
; Check for 'DataPoints' keyword
;
DataPoints = GetKeyWord('DataPoints', localExtra)
IF(TypeOf(DataPoints) EQ 11) THEN BEGIN
	FOR i=0, N_ELEMENTS(DataPoints)-1 DO DataPoints[i] -> oPlot, color=1+(i MOD 9), lineStyle=0, psym=3
ENDIF
;
; Check for 'FieldMap' keyword
;
FieldMap = GetKeyWord('FieldMap', localExtra)
IF(TypeOf(FieldMap) EQ 11) THEN BEGIN
  info = SIZE(FieldMap)
  imax = info[1] - 1
  jmax = info[2] - 1
  FOR i=0, imax DO BEGIN
    FOR j=0,jmax DO BEGIN 
      fieldMap[i,j] -> oPlot, color=1+(i MOD 9), lineStyle=0, psym=3
    ENDFOR
  ENDFOR
ENDIF

;
; Check if any over plots are requested
;
oPlotType = TypeOf(oplot)
IF(oPlottype EQ 11) THEN BEGIN
	;
	; Check for 'NoErrorBars' keyword
	;
	result = GetKeyWord('NoErrorBars', localExtra)
	IF(TypeOF(result) NE 7) THEN noErrorBars = result
	noColor = 0b
	result = GetKeyWord('NoColor', localExtra)
	IF(TypeOf(result) NE 7) THEN noColor = result
	CASE noColor OF
		0	:	FOR i=0, (N_ELEMENTS(oPlot) -1) DO	$
					oPlot[i] -> oPlot, color=2+(i MOD 9), NoErrorBars=noErrorBars
		1	:	FOR i=0, (N_ELEMENTS(oPlot) -1) DO	$
					oPlot[i] -> oPlot, LineStyle=1+(i MOD 5), NoErrorBars=noErrorBars
		else	:
	ENDCASE
	localExtra=extra
	result = GetKeyWord('thick', localExtra)
	IF(typeOf(localExtra) EQ 8) THEN BEGIN				; replot 'obj' ON TOP OF any overplots 
		Obj -> oPlot, thick=2*!P.thick,			$  	; with a heavy line
				NoErrorBars=noErrorBars
	ENDIF ELSE BEGIN
		Obj -> oPlot, thick=2*!P.thick
	ENDELSE
ENDIF
 

RETURN
END	; ****** PlotObj1D ****** ;


PRO PlotObj2D, Obj, _EXTRA=extra
;
;		This proceedure is part of a 'wrapper' which allows me
;		to use D. Fanning's 'XWINDOW' widget to display GKV
;		objects.
;
; Keywords:
;
;		Keywords native to this proceedure are:
;
;			Shade_Surf		Set this keyword (i.e., put "/Shade_Surf" on the
;						command line) to obtain a surface plot for objects
;						of 2 (or more) dimensions (default is an 'Image" 
;						plot).
;
;		Remaining Keywords are passed first to XWINDOW (see XWINDOW source
;		for a descreption of XWINDOW keywords) and then to the object's
;		drawing routine (see GKVs1D::Draw, GKVs2D::Draw, or GKVs2D::Shade_Surf).
;		In order to preserve the GKVsD objects' drawing routine's runtime 
;		keywords, it is necessary to spell out all keywords in full
;		(i.e., keyword abreviation will not work!).
;
; Written by W.M. Nevins
;	7/16/00
;
; Avoid changing contents of 'extra'
;
localExtra = -1
IF(TypeOF(extra) EQ 8) THEN localExtra = extra
;
; Set 'bottom' and 'ncolors'
bottom = !COLOR_SETUP_NCOLORS + 1
ncolors = !D.TABLE_SIZE - bottom
;
; First task is to check for the local keyword, Shade_Surf.
;
shadeSurf = 0
result = GetKeyWord('Shade_Surf', localExtra)
IF Query_Integer(result) THEN shadeSurf = result
;
; Draw Obj in currently open window
; 
IF(shadeSurf NE 0) THEN BEGIN
	IF(TypeOf(localExtra) EQ 8) THEN BEGIN
		Obj -> Shade_Surf, _EXTRA=localExtra
	ENDIF ELSE BEGIN
		Obj -> Shade_Surf
	ENDELSE
ENDIF ELSE BEGIN
	IF(TypeOf(localExtra) EQ 8) THEN BEGIN
		Obj -> Draw, bottom=bottom, ncolors=ncolors, _EXTRA=localExtra
	ENDIF ELSE BEGIN
		Obj -> Draw, bottom=bottom, ncolors=ncolors
	ENDELSE
ENDELSE

RETURN
END	; ****** PlotObj2D ****** ;
