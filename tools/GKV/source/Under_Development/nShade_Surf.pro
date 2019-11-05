Pro GKVs2D::Shade_Surf, _Extra=extra
;
;  'Virtual' keywords:
;
;			  xrange=x_range,  title=T_title, xtitle=x_title, ytitle=y_title, 		$
;			  yrange=y_range, Grid1=Grid_1, Grid2=Grid_2, indx1=indx_1, indx2=indx_2,	$
;			  zrange=z_range, ztitle=z_title, position=pstn, charsize=chsz,	$
;			  Pretty=pretty, Polar=polar, Log=log, Vrange=vrange
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

result = GetKeyWord('ztitle', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN z_title=result	;	ytitle = y_title
result = GetKeyWord('zrange', extra)
IF(TypeOf(result) NE 7) THEN z_range=result		;	yrange = y_range

result = GetKeyWord('position', extra)
IF(TypeOf(result) EQ 7) THEN $
	IF(result NE 'undefined') THEN pstn=result		;	position = pstn
result = GetKeyWord('charsize', extra)
IF(TypeOf(result) NE 7) THEN chsz=result			;	charsize = chsz


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
IF(TypeOF(result) NE 7) THEN vrange=result			;	using 'mnemonic' = [vmin, vmax]
;
;										; on the command line

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
imax = imax < iimax < nx-1
xmax=(*Grid1.values)[imax]
xrange = [xmin, xmax]

indices = self -> IndexString(plotAxis)
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
	ELSE	:	MESSAGE, 'Bad value for ndmis', /Informational
ENDCASE
values = REFORM(values, imax-imin+1, jmax-jmin+1, /OVERWRITE)
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

ncolors=220
vmin = self.vrange[0]
vmax = self.vrange[1]
vvmin = vmin
vvmax = vmax
IF(N_ELEMENTS(vrange) EQ 2) THEN BEGIN
	vmin=vrange[0]
	vmax=vrange[1]
ENDIF
IF(KEYWORD_SET(Log)) THEN BEGIN
;	values = ALOG10(values)
	vvmax = ALOG10(vmax)
	vvmin = vvmax-12
	IF(vmin GT 0) THEN vvmin = ALOG10(vmin)
	vmin = exp(vvmin*ALOG(10.))
ENDIF

zrange=[vmin, vmax]
IF(N_ELEMENTS(z_range) EQ 2) THEN zrange=z_range

ztitle=self.mnemonic + " (" + self.units + ")"
IF KEYWORD_SET(Pretty)  THEN ztitle=self.title + " (" + self.units + ")"
IF KEYWORD_SET(z_title) THEN ztitle=z_title

position=[0.15, 0.15, 0.9, 0.9] 
IF KEYWORD_SET(pstn) THEN position=pstn
charsize=2.
IF KEYWORD_SET(chsz) THEN charsize=chsz

IF(N_ELEMENTS(extra) NE 0) THEN BEGIN
	IF(TypeOf(extra) EQ 8) THEN nextra = extra
ENDIF

SHADE_SURF, values, x, y, xrange=xrange, yrange=yrange, max_value=vmax, min_value=vmin, zrange=zrange, $
		 title=title,xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, zlog=log, $
		 position=position, charsize=charsize, _Extra=nextra 

xwrite=!D.x_ch_size						; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower left-hand corner of plot window
XYOUTS, xwrite, ywrite, self.CodeName, /Device	; Write CodeName to lower left-hand corner
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines
XYOUTS, xwrite, ywrite, self.CodePI, /Device	; Write CodePI below CodeName
xwrite=!D.x_size -!D.x_ch_size				; Device coordinates for writting to
ywrite=2.0*!D.y_ch_size					; lower right-hand corner of plot window
XYOUTS, xwrite, ywrite, self.RunID,	$
		 Alignment=1., /Device 			; Write RunID to lower right-hand corner
ywrite=ywrite-1.5*!D.y_ch_size				; Move down 1.5 lines
XYOUTS, xwrite, ywrite, self.FileID,	$ 
		Alignment=1. , /Device			; write FileID to below RunID
END ; ****** GKVs2D::Shade_Surf ****** ;
