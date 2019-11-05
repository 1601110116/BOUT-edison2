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


Pro GKVs2D::Movie, Obj0, obj1, obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9, _Extra=extra
;
; Purpose:
;
;		This proceedure creates and displays an animation
;		of the data in 'self' using the IDL tool XINTERANIMATE
;
; Arguments:
;
;			None
;
;
; KeyWords:
;
;	imin		Index within 'self' to the first frame of the animation.
;			Defaults to the first time-slice of 'self'. (Optional)
;
;	imax		Index within 'self' to the final frame of the animation.
;			Defaults to the final time-slice of 'self'. (Optional)
;
;  'Mnemonic'	Alternatively, you can specify the time-interval for the
;			animation in the format "Mnemonic = [start, end]", where
;			"Mnemonic" is the mnemonic for the time-like (third) axis,
;			while 'start' and 'end' are the values of the time-like
;			coordinate at the stare and end of the animation.
;
;	iskip		The animation uses only "iskip"th time-slice.
;			Defaults to the smallest integer such that the total
;			number of frames in the animation is � 128.  
;			NOTE: on many (all?) platforms IDL appears to limit
;			the total number of PIXMAP windows (one PIXMAP window
;			is required for for each frame of the animation) to 128,
;			so the user's choice of iskip is overriden if it would
;			result in more than 128 frame in the animation. (Optional)
;			
;	npixels	The number of pixels to be used in each frame of the animation.
;			If npixels is set to a scalar, a square window of (npixels x npixels)
;			is produced.  If npixels is set to a two-element array, then the
;			a (generally) rectangular windo of (nPixels[0] x nPixels[1]) is produced.
;			Defaults to 400x400.  (Optional)
;
;
; Graphics Keywords:
;
;			Any additional keywords specified on the command line will forwarded to the 
;			plotting routines, allowing the user to customize his animation.
;
; Written by W.M. Nevins
;	10/3/00
;
; Set defaults and get keywords form command line
;
irange=self.Grid2.irange
imin = irange[0]
imax=irange[1]
result = GetKeyWord('imin', extra)
IF(TypeOf(result) NE 7) THEN imin=result
result = GetKeyWord('imax', extra)
IF(TypeOf(result) NE 7) THEN imax=result
iskip=(imax-imin)/128 + 1
result = GetKeyWord('iskip', extra)
IF(TypeOf(result) NE 7) THEN iskip = result > iskip
timeMnemonic = self.Grid2.mnemonic
result = GetKeyWord(timeMnemonic, extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	 CASE N_ELEMENTS(result) OF
	 	1:	BEGIN
	 			imin = self -> AxisIndex(2, result)
	 			iskip= iskip > ((imax-imin)/128 + 1)
	 		END
	 	2:	BEGIN
	 			imin = self -> AxisIndex(2,result[0])
	 			imax = self -> AxisIndex(2,result[1])
	 			iskip= iskip > ((imax-imin)/128 + 1)
	 		END
	 	3:	BEGIN
	 			imin = self -> AxisIndex(2,result[0])
	 			imax = self -> AxisIndex(2,result[1])
	 			dt = (*self.Grid2.values)[imin+1] - (*self.Grid2.values)[imin]
	 			iskip = (result[3]/dt) > ((imax-imin)/128 + 1)
	 		END
		ELSE:	MESSAGE, "Couldn't parse 'Mnemonic = ...' on input line.  Will use defaults", /INFORMATIONAL
	ENDCASE
ENDIF
nPixels=[400, 400]
result = GetKeyWord('nPixels', extra)
IF(typeOf(result) NE 7) THEN BEGIN
	CASE N_ELEMENTS(result) OF
		1:	nPixels = [FIX(result),    FIX(result)]
		2:	nPixels = [FIX(result[0]), FIX(result[1])] 
		ELSE:	MESSAGE, "Couldn't parse 'nPixels = ...'.  Will use defaults", /INFORMATIONAL
	ENDCASE
ENDIF
vrange = self.vrange
numoPlots = 0
IF(N_PARAMS() GT 0) THEN BEGIN
	tmin = (*self.Grid2.Values)[imin]
	tmax = (*self.Grid2.Values)[imax]
	oPlotObjs = obj0
	FOR i=1, N_PARAMS()-1  DO BEGIN
		argString = STRING(i, FORMAT='(I1)')
		commandString = "oPlotObjs = [oPlotObjs, Obj" + argString + "]"
		ok = EXECUTE(commandString)
	ENDFOR
	numoPlots = N_ELEMENTS(oPlotObjs)
	FOR i=0, numoPlots-1 DO BEGIN
		oPlotObjs[i] -> Get, axis=2, range=thisRange
		oPlotObjs[i] -> SignalWindow, axis=2, range=[tmin,tmax]
		oPlotObjs[i] -> Get, vrange=thisVrange
		oPlotObjs[i] -> SignalWindow, axis=2, range=thisRange
		vrange[0] = vrange[0] < thisVrange[0]
		vrange[1] = vrange[1] > thisVrange[1]
	ENDFOR
ENDIF
nframes=(imax - imin)/iskip + 1
;
; Close any existing XINTERANIMATE
;
;XINTERANIMATE, /CLOSE
;
; and start a new onw
;
XINTERANIMATE, SET=[nPixels[0], nPixels[1], nframes], /SHOWLOAD
iframe=0
FOR i=imin, imax, iskip DO BEGIN
  nExtra = Extra
	thisSlice = self -> Slice(axis=2, Index=i)
	thisSlice -> Draw, vrange=vrange, _EXTRA=nExtra
	thisSlice -> Trash
	IF(numoPlots GT 0) THEN BEGIN
		FOR j=0, numoPlots-1 DO BEGIN
		nExtra = Extra
		thisPlot = oPlotObjs[j] -> Slice(axis=2, index=i)
		thisPlot -> oPlot, color=((j) MOD 10)+2, _EXTRA=nExtra
		thisPlot -> Trash
		ENDFOR
	ENDIF
	XINTERANIMATE, FRAME=iframe, Window=!D.WINDOW
	iframe=iframe+1
ENDFOR
XINTERANIMATE
END ; ****** GKVs2D::movie ****** ;



Pro GKVs3D::Movie, _Extra=extra
;
; Purpose:
;
;		This proceedure creates and displays an animation
;		of the data in 'self' using the IDL tool XINTERANIMATE
;
; Arguments:
;
;			None
;
; KeyWords:
;
;	imin		Index within 'self' to the first frame of the animation.
;			Defaults to the first time-slice of 'self'. (Optional)
;
;	imax		Index within 'self' to the final frame of the animation.
;			Defaults to the final time-slice of 'self'. (Optional)
;
;  'Mnemonic'	Alternatively, you can specify the time-interval for the
;			animation in the format "Mnemonic = [start, end]", where
;			"Mnemonic" is the mnemonic for the time-like (third) axis,
;			while 'start' and 'end' are the values of the time-like
;			coordinate at the stare and end of the animation.
;
;	iskip		The animation uses only "iskip"th time-slice.
;			Defaults to the smallest integer such that the total
;			number of frames in the animation is � 128.  
;			NOTE: on many (all?) platforms IDL appears to limit
;			the total number of PIXMAP windows (one PIXMAP window
;			is required for for each frame of the animation) to 128,
;			so the user's choice of iskip is overriden if it would
;			result in more than 128 frame in the animation. (Optional)
;			
;	npixels	The number of pixels to be used in each frame of the animation.
;			If npixels is set to a scalar, a square window of (npixels x npixels)
;			is produced.  If npixels is set to a two-element array, then the
;			a (generally) rectangular windo of (nPixels[0] x nPixels[1]) is produced.
;			Defaults to 400x400.  (Optional)
;
;  Shade_Surf	Set this keyword (i.e., put "/Shade_Surf" on the command line) to produce
;			a surface plot.  Default is to produce an "image" plot (that is, essentially
;			a color contour plot).
;
;
; Graphics Keywords:
;
;			Any additional keywords specified on the command line will forwarded to the 
;			plotting routines, allowing the user to customize his animation.
;
; Written by W.M. Nevins
;	10/3/00
;
; Set defaults and get keywords form command line
;
irange=self.Grid3.irange
imin = irange[0]
imax=irange[1]
result = GetKeyWord('imin', extra)
IF(TypeOf(result) NE 7) THEN imin=result
result = GetKeyWord('imax', extra)
IF(TypeOf(result) NE 7) THEN imax=result
iskip=(imax-imin)/128 + 1

timeMnemonic = self.Grid3.mnemonic
result = GetKeyWord(timeMnemonic, extra)
IF(TypeOf(result) NE 7) THEN BEGIN
	 CASE N_ELEMENTS(result) OF
	 	1:	BEGIN
	 			imin = self -> AxisIndex(3, result)
	 			iskip= (imax-imin)/128 + 1
	 		END
	 	2:	BEGIN
	 			imin = self -> AxisIndex(3,result[0])
	 			imax = self -> AxisIndex(3,result[1])
	 			iskip= (imax-imin)/128 + 1
	 		END
	 	3:	BEGIN
	 			imin = self -> AxisIndex(3,result[0])
	 			imax = self -> AxisIndex(3,result[1])
	 			dt = (*self.Grid3.values)[imin+1] - (*self.Grid3.values)[imin]
	 			iskip = (result[3]/dt) 
	 		END
		ELSE:	MESSAGE, "Couldn't parse 'Mnemonic = ...' on input line.  Will use defaults", /INFORMATIONAL
	ENDCASE
ENDIF
result = GetKeyWord('iskip', extra)
IF(TypeOf(result) NE 7) THEN iskip = result

shadeSurf = 0
result = GetKeyWord("Shade_Surf", extra)
IF(TypeOf(result) NE 7) THEN shadeSurf = KEYWORD_SET(result)
nPixels=[400, 400]
result = GetKeyWord('nPixels', extra)
IF(typeOf(result) NE 7) THEN BEGIN
	CASE N_ELEMENTS(result) OF
		1:	nPixels = [FIX(result),    FIX(result)]
		2:	nPixels = [FIX(result[0]), FIX(result[1])] 
		ELSE:	MESSAGE, "Couldn't parse 'nPixels = ...'.  Will use defaults", /INFORMATIONAL
	ENDCASE
ENDIF
nframes=(imax - imin)/iskip + 1
;
; Close any existing XINTERANIMATE
;
;XINTERANIMATE, /CLOSE
;
; and start a new onw
;
XINTERANIMATE, SET=[nPixels[0], nPixels[1], nframes], /SHOWLOAD
iframe=0
FOR i=imin, imax, iskip DO BEGIN
	IF(KEYWORD_SET(ShadeSurf)) THEN BEGIN
		self -> Shade_Surf, indx1=i, _EXTRA=extra
	ENDIF ELSE BEGIN
		self -> Draw, indx1=i, _EXTRA=extra
	ENDELSE
	XINTERANIMATE, FRAME=iframe, Window=!D.WINDOW
	iframe=iframe+1
ENDFOR
XINTERANIMATE
END ; ****** GKVs3D::movie ****** ;
