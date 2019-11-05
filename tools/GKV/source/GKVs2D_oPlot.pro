

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
;			oPlot			This keyword should be set to an object array
;						array containing GKV objects
;						to be layen over color contour plot of "self".
;
;			ObjColors		Set this key word to an integer array referencing
;						the colors you wish for the corresponding oPlot 
;						object.
;
;			NoColor			Set this key word (i.e., /NoColor on the command line)
;						if output is to be in black and white. Curves will
;						then be dots, dashes, etc. 
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
;this=that
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
	RETURN
ENDIF 
;
; Not a surface plot, so make color image in guise of a 
; contour plot
;
IF(TypeOf(localExtra) EQ 8) THEN BEGIN
	Obj -> Draw, bottom=bottom, ncolors=ncolors, _EXTRA=localExtra
ENDIF ELSE BEGIN
	Obj -> Draw, bottom=bottom, ncolors=ncolors
ENDELSE
;
; Check if any over plots are requested,
; by examining the local keyword, oplot.
; If so, put contour lines on top of the 
; color image.
;
result = GetKeyWord('oPlot', localExtra)
IF(TypeOf(result) EQ 11) THEN oPlot = result
oPlotType = TypeOf(oplot)
IF(oPlottype EQ 11) THEN BEGIN
	ObjColors = INDGEN(11)
	noColor = 0b
	result = GetKeyWord('NoColor', localExtra)
	IF(TypeOf(result) NE 7) THEN noColor = result
	result = GetKeyWord('ObjColors', localExtra)
	IF Query_Integer(result) THEN ObjColors=result
	CASE noColor OF
		0	:	FOR i=0, (N_ELEMENTS(oPlot) -1) DO	$
					oPlot[i] -> oPlot, color=ObjColors[i Mod 11], _Extra=LocalExtra
		1	:	FOR i=0, (N_ELEMENTS(oPlot) -1) DO	$
					oPlot[i] -> oPlot, LineStyle=1+(i MOD 5), _Extra=LocalExtra
		else	:
	ENDCASE
ENDIF
 
RETURN
END	; ****** PlotObj2D ****** ;


Pro GKVs2D::Draw, _Extra=extra
;
;  'Virtual' keywords:
;
;			  xrange=x_range,  title=T_title, xtitle=x_title, ytitle=y_title, 		$
;			  yrange=y_range, Grid1=Grid_1, Grid2=Grid_2, indx1=indx_1, indx2=indx_2,	$
;			  Pretty=pretty, Polar=polar, Log=log, Vrange=vrange, Bottom
;
; Default Plotting routine for GKV objects (2-D or more).
; Allows any legal (to PLOT) graphics keyword, but gets defaults from
; object definition if not supplied on command line.
;
; In additon, allows user to set range and indices using axis mnemonics
;
; First, parse 'extra' for keywords which we wish to intercept 
; (this is the moral equivalent of defining them on the PRO line...
;  but prevents IDL from enforcing keyword abreviation, which would
;  interfer with use of mnemonics as 'run-time' keywords)
;
result = GetKeyWord('xrange', extra)				; Equivalent to having
IF(TypeOf(result) NE 7) THEN x_range = result		;	xrange = x_range 
result = GetKeyWord('xtitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN x_title=result	;	xtitle = x_title
result = GetKeyWord('title', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN T_title=result	;	title = T_title
result = GetKeyWord('ytitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN y_title=result	;	ytitle = y_title
result = GetKeyWord('yrange', extra)
IF(TypeOf(result) NE 7) THEN y_range=result		;	yrange = y_range
result = GetKeyWord('Grid1', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_1=result	;	Grid1 = Grid_1
IF(TypeOf(result) NE 7) THEN Grid_1=result
result = GetKeyWord('Grid2', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN Grid_2=result	; 	Grid2 = Grid_2
IF(TypeOf(result) NE 7) THEN Grid_2=result
result = GetKeyWord('indx1', extra)
IF(TypeOf(result) NE 7) THEN indx_1=result			; 	indx1 = indx_1
result = GetKeyWord('indx2', extra)
IF(TypeOf(result) NE 7) THEN indx_2=result			; 	indx2	 = indx_2
result = GetKeyWord('Pretty', extra)
IF(TypeOf(result) NE 7) THEN pretty=result			;	Pretty = pretty
result = GetKeyWord('Polar', extra)
IF(TypeOf(result) NE 7) THEN polar=result			;	Polar = polar
result = GetKeyWord('Log', extra)
IF(TypeOF(result) NE 7) THEN Log=result			;	Log = log
result = GetKeyWord('Vrange', extra)
IF(TypeOF(result) NE 7) THEN vrange=result			;	Vrange=vrange
result = GetKeyWord(self.mnemonic, extra)			;	or, set vrange
IF(TypeOF(result) NE 7) THEN vrange=result			;	using 'mnemonic' = [vmin, vmax]  on the command line
bottom = !COLOR_SETUP_NCOLORS + 1
result = GetKeyWord('bottom', extra)
IF(TypeOF(result) NE 7) THEN bottom = result		;	Bottom = bottom
ncolors = !D.TABLE_SIZE - bottom
result = GetKeyWOrd('ncolors', extra)				;	ncolors = ncolors
IF(TypeOF(result) NE 7) THEN ncolors = result		;

Grid1default = 1
Grid2defalut = 2
;
; Find Grid structure for first axis
;
IF(N_ELEMENTS(Grid_1) EQ 0) THEN Grid_1=1			; default value
IF(TypeOF(Grid_1) EQ 7) THEN Grid_1=(self -> AxisNumber(Grid_1)) > 1 
str_1 = STRING(Grid_1, FORMAT='(I1)')
command_str = 'Grid1=self.Grid' + str_1
ok = EXECUTE(command_str)
;
; Find Grid structure for second axis
;
IF(N_ELEMENTS(Grid_2) EQ 0) THEN Grid_2=Grid_1+1	; default value
IF(TypeOF(Grid_2) EQ 7) THEN Grid_2=(self -> AxisNumber(Grid_2)) > (Grid_1+1) 
str_2 = STRING(Grid_2, FORMAT='(I1)')
command_str = 'Grid2=self.Grid' + str_2
ok = EXECUTE(command_str)
ndims = self -> NumDims()			; Get number of dimensions
plotAxis=INTARR(ndims+1)
plotAxis[Grid_1]=1
plotAxis[Grid_2]=2
indx=-1
FOR iaxis=1, ndims DO BEGIN
	IF(plotAxis[iaxis] GT 0) THEN GOTO, DONE01
	plotAxis[iaxis] = indx
	indx=indx-1
DONE01 :  ; moral equivalent of 'continue'? 
ENDFOR
	
;
; Check command line for axis mnemonics
;
axisInfo = self -> GetAxis(extra)
;
; 
FOR iaxis=1, ndims DO BEGIN
	axisValue=axisInfo.(iaxis-1)
	IF(TypeOF(axisValue) NE 7) THEN BEGIN	; Mnemonic for iaxis is used on command line
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -1)) THEN indx_1= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 1) AND (plotAxis[iaxis] EQ -2)) THEN indx_2= self -> AxisIndex(iaxis, AxisValue)
		IF((N_ELEMENTS(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 1)) THEN x_range=AxisValue
		IF((N_Elements(AxisValue) EQ 2) AND (plotAxis[iaxis] EQ 2)) THEN y_range=AxisValue
	ENDIF
ENDFOR
;
; put 'indx_1', 'indx_2' into plot axis if these are set.
;
IF(N_ELEMENTS(indx_1) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -1) THEN plotAxis[iaxis]=-indx_1
	ENDFOR
ENDIF
IF(N_ELEMENTS(indx_2) EQ 1) THEN BEGIN
	FOR iaxis=1, ndims DO BEGIN
		IF(plotAxis[iaxis] EQ -2) THEN plotAxis[iaxis]=-indx_2
	ENDFOR
ENDIF
; 
; Get limits to signal window for each axis
;
imin=Grid1.irange[0]
imax=Grid1.irange[1]
jmin=Grid2.irange[0]
jmax=Grid2.irange[1]
;
; use default plot keyword values from Object definition if they
; are not over ridden on command line.
;
xrange=Grid1.range
IF(N_ELEMENTS(x_range) EQ 2) THEN xrange = x_range

xmin=xrange[0]
xmax=xrange[1]
nx = N_ELEMENTS(*Grid1.values)
temp=(*Grid1.values - xmin)^2			; Finds iimin such that Grid1.values(iimin)
smallest = MIN(temp, iimin)			;	takes value closest to xmin
IF(iimin eq imax) THEN iimin=imin
imin = imin > iimin > 0
xmin=(*Grid1.values)[imin]

temp=(*Grid1.values - xmax)^2			; sets iimax such that Grid1.values(iimax)
smallest = MIN(temp, iimax)			;	takes value closest to xmax
IF(iimax eq imin) THEN iimax=imax
imax = imax < iimax < (nx-1)
xmax=(*Grid1.values)[imax]
xrange = [xmin, xmax]

indices = self -> IndexString(plotAxis, pretty=pretty)
indexStr = '[' + STRJOIN(indices, ', ') + ']'
IF(KEYWORD_SET(pretty)) THEN BEGIN
	title=self.title + indexStr + " (" + self.units + ")"
ENDIF ELSE BEGIN
	title=self.mnemonic + indexStr + " (" + self.units + ")"
ENDELSE
IF KEYWORD_SET(T_title) THEN title=T_title

x_title=Grid1.mnemonic + " (" + Grid1.units + ")"
IF KEYWORD_SET(Pretty)  THEN x_title=Grid1.title + " (" + Grid1.units + ")"
IF KEYWORD_SET(x_title) THEN xtitle=x_title

ytitle=Grid2.mnemonic + " (" + Grid2.units + ")"
IF KEYWORD_SET(Pretty)  THEN ytitle=Grid2.title + " (" + Grid2.units + ")"
IF KEYWORD_SET(y_title) THEN ytitle=y_title

yrange=Grid2.range
IF (N_ELEMENTS(y_range) EQ 2) THEN yrange = y_range

ymin=yrange[0]
ymax=yrange[1]
ny = N_ELEMENTS(*Grid2.values)
temp=(*Grid2.values - ymin)^2			; Finds jjmin such that Grid2.values(iimin)
smallest = MIN(temp, jjmin)			;	takes value closest to ymin
IF(jjmin eq jmax) THEN jjmin=jmin
jmin = jmin > jjmin > 0
ymin=(*Grid2.values)[jmin]

temp=(*Grid2.values - ymax)^2			; sets jjmax such that Grid2.values(iimax)
smallest = MIN(temp, jjmax)			;	takes value closest to ymax
IF(jjmax eq jmin) THEN jjmax=jmax
jmax = jmax < jjmax < (ny-1)
ymax=(*Grid2.values)[jmax]
yrange = [ymin, ymax]

x = (*Grid1.values )[imin:imax]
y = (*Grid2.values )[jmin:jmax]
;
; POLYMORPHISM!  Make this work for subclasses of GKVs2D
;
indx1=0
indx2=0
IF(N_ELEMENTS(indx_1)) THEN indx1=indx_1
IF(N_ELEMENTS(indx_2)) THEN indx2=indx_2
CASE ndims OF
	2:	values = (*self.values)[imin:imax, jmin:jmax]
	3:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1]	
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=3', /Informational
				ENDCASE
				END
				
			2:	values = (*self.values)[indx1, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, "Bad value for Grid_1, ndims=3", /Informational
		ENDCASE
		END 
	4:	BEGIN
		CASE Grid_1 OF
			1:	BEGIN
				CASE Grid_2 OF
					2	:	values = (*self.values)[imin:imax, jmin:jmax, indx1, indx2] 
					3	:	values = (*self.values)[imin:imax, indx1, jmin:jmax, indx2]
					4	:	values = (*self.values)[imin:imax, indx1, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=1, ndims=4', /Informational
				ENDCASE
				END
			2:	BEGIN
				CASE Grid_2 OF
					3	:	values = (*self.values)[indx1, imin:imax, jmin:jmax, indx2] 
					4	:	values = (*self.values)[indx1, imin:imax, indx2, jmin:jmax]
					ELSE	:	MESSAGE, 'Bad value for Grid_2, Grid_1=2, ndims=4', /Informational
				ENDCASE
				END
			3:	 values = (*self.values)[indx1, indx2, imin:imax, jmin:jmax]
			ELSE	:	MESSAGE, 'Bad value for Grid_1, ndims=4', /Informational
		ENDCASE
		END
	ELSE	:	MESSAGE, 'Bad value for ndims', /Informational
ENDCASE
values = FLOAT(REFORM(values, imax-imin+1, jmax-jmin+1, /OVERWRITE))
info=SIZE(values)
index=info[0]+1
vtype=info[index]

IF(KEYWORD_SET(Polar)) THEN BEGIN
	Nrad = N_ELEMENTS(x)
	Rmax = MAX(x)
	dx = Rmax/Nrad
	spacing=[dx,dx]
	bounds=[-Rmax, -Rmax, Rmax, Rmax]
	values = POLAR_SURFACE(values, x, y, /GRID, SPACING=spacing, BOUNDS=bounds)
	info = SIZE(values)
	index=info[0]+1
	vtype=info[index]
	Nx = info[1]
	x = -Rmax + dx*FINDGEN(NX)
	y = x
	xrange = [-Rmax, Rmax]
	yrange = xrange
	xtitle = 'R-R!Io!N (a)'
	ytitle = 'Z (a)'
ENDIF

vmin = self.vrange[0]
vmax = self.vrange[1]
IF(N_ELEMENTS(vrange) EQ 2) THEN BEGIN
	vmin=vrange[0]
	vmax=vrange[1]
ENDIF
vvmin = vmin
vvmax = vmax
IF(KEYWORD_SET(Log)) THEN BEGIN
	values = ALOG10(values)
	vvmax = ALOG10(vmax)
	vvmin = vvmax-12
	IF(vmin GT 0) THEN vvmin = ALOG10(vmin)
	vmin = exp(vvmin*ALOG(10.))
ENDIF

IF(N_ELEMENTS(extra) NE 0) THEN BEGIN
	IF(TypeOf(extra) EQ 8) THEN nextra = extra
ENDIF
;
; Set up 'image' array
;
ncolors = BYTE(ncolors)
bottom  = BYTE(bottom)
IF(ncolors GE (256 - bottom)) THEN ncolors = 256 - bottom
image = BYTSCL(values, Min=vvmin, Max=vvmax, Top=ncolors-1, /NaN) + bottom
;
; Display image
;
TVImage, image, position=[0.15, 0.15, 0.9, 0.8], /Erase, /Minus_One, _Extra=nextra 
;
; Choose an appropriate character size
; (scaled to size of currently open graphics window)
;
newSizes = GKV_CharSize(Old_Sizes = oldSizes)
;
; Reset character size
;
DEVICE, SET_CHARACTER_SIZE = newSizes
;
; Draw axis, etc. over image
;
;
; These 'IF's' are here to insure that the following 'CONTOUR' (which only draws the axises)
;	does not fail when there is no variation in 'values'.
; 
IF(vmin EQ vmax) THEN values[0,0]=values[0,0] +MAX(values) + 1.	
IF( NOT FINITE(vmax)) THEN BEGIN
	values = FLOAT(values) < 100.*vmin
	values[0,0]=values[0,1]+1.
ENDIF
;
; Move call to CONTOUR to end to enable calls following to GKVs2D::oPlot
; (otherwise position, etc. will be left at values set by 'colorbar') 
;  
;CONTOUR, 	values, x, y, position=[0.15, 0.15, 0.9, 0.8], 	$
;		xrange=xrange, yrange=yrange,					$
;		xtitle=xtitle, ytitle=ytitle, XStyle=1, YStyle=1, 	$
;		Ticklen=0.02, /NODATA, /NOErase, _Extra=nextra
;
; Add color bar
;
ColorBar, Ncolors=ncolors, bottom = bottom, range=[vmin,vmax], 	$
		Format='(G10.2)', position=[0.145, 0.85, 0.9, 0.9], Xlog=log
;
; Mark average value on color bar
;
IF (NOT KEYWORD_SET(Log)) THEN BEGIN
	avg=TOTAL(values)/N_elements(values)
	ARROW,  avg, 0.45*(vmax-vmin)+vmin,  avg, vmin, Thick=2., /Data
	XYOUTS, avg, 0.55*(vmax-vmin)+vmin, 'AVG', Alignment=0.5, /Data
ENDIF
;
; Put title at top of Graphics window
;
charSize = GKV_TitleSize(title, Width = (0.9-0.15), yPosition=ywrite)
xwrite = (0.15 + 0.9)/2.*!D.X_SIZE
XYOUTS, xwrite, ywrite, title, ALIGNMENT=0.5, CHARSIZE=charSize, /DEVICE 
;
; Now write code and run info in lower corners
;
result = GetKeyWord("NoFooter", Extra)
IF(Query_Integer(result)) THEN GOTO, NoFooter 

maxChars =ROUND(0.4*!D.X_SIZE/!D.X_CH_SIZE)		; Compute maximum allowed characters in a line
xwrite=!D.x_ch_size					; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower left-hand corner of plot window.
codeName = STRMID(self.CodeName, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, CodeName, /Device		; Write CodeName to lower left-hand corne.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
codePI = STRMID(self.CodePI, 0, maxChars)		; Truncate 'CodePI' if necessary
XYOUTS, xwrite, ywrite, CodePI, /Device			; Write CodePI below CodeName.
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window.
runID = STRMID(self.RunID, 0, maxChars)			; Truncate 'runID' if necessary
XYOUTS, xwrite, ywrite, RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner.
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines.
fileID = STRMID(self.FileID, 0, maxChars)		; Truncate 'CodeName' if necessary
XYOUTS, xwrite, ywrite, FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID

CONTOUR, 	values, x, y, position=[0.15, 0.15, 0.9, 0.8], 		$
		xrange=xrange, yrange=yrange,				$
		xtitle=xtitle, ytitle=ytitle, XStyle=1, YStyle=1, 	$
		Ticklen=0.02, /NODATA, /NOErase; , _Extra=nextra

;
; Return character size to input values
;
NoFooter:  DEVICE, SET_CHARACTER_SIZE = oldSizes
RETURN

END ; ***** GKVs2D::Draw ***** ;


PRO GKVs2D::oPlot, RealOnly=realOnly, _Extra=extra

imin=self.Grid1.irange[0]
imax=self.Grid1.irange[1]
x = (*self.Grid1.Values)[imin:imax]
xRange = [x[0], x[iMax-imin]]
jmin=self.Grid2.irange[0]
jmax=self.Grid2.irange[1]
y = (*self.Grid2.values)[jmin:jmax]
yRange = [y[0], y[jMax-jMin]]
z = (*self.values)[imin:imax, jmin:jmax]
ztype=TypeOf(z)
;
; Check for 'color' in _Extra to pass to ERRPLOT
;
localExtra=-1
IF(TypeOf(extra) NE 0) THEN localExtra=extra
result = GetKeyWord('color', localExtra)
IF Query_Integer(result) THEN BEGIN
	color=result
ENDIF ELSE BEGIN
	color=1
ENDELSE
nLevels = 10
result = GetKeyWord('nLevels', localExtra)
IF Query_Integer(result) THEN nLevels=result
zMin = MIN(z, MAX=zMax)
dz = (zMax - zMin)/(nLevels+1)
Levels = zMIn + dz/2. + dz*FINDGEN(nLevels)
result = GetKeyWord('Levels', localExtra)
IF(Query_Real(result) + Query_Integer(result)) THEN levels=result


;CONTOUR, z, x, y, position=[0.15, 0.15, 0.9, 0.8],	$
;	 xRange=xRange, yRange=yRange, 			$
;	 xStyle=1, yStyle=1, /NoDATA, /NoErase
CONTOUR, z, x, y, LEVELS=Levels, _EXTRA=Extra, /OVERPLOT

RETURN
END ; ***** GKVs2D::oPlot ***** ;
