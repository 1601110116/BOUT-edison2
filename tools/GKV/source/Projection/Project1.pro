FUNCTION GKVs2D::PROJECT1, _Extra=extra
;
; Purpose:
;
;	Projects 'self' onto first independent variable
;	using specifiec weighting function.
;
; Written by W.M. Nevins
;	4/9/04
;
FORWARD_FUNCTION GKVsd_MIN
outputStr = {GKVs1D}
nTags = N_TAGS({GKVsd})
FOR i=0, nTags-1 DO outputStr.(i) = self.(i)
outputStr.Indices = PTR_NEW(['*'])
outputStr.Grid1 = GKVsd_GridCopy(self.grid1)
iMin = self.Grid1.irange[0]
iMax = self.Grid1.irange[1]
template = self -> slice(axis=1, index=iMin)

TryAgain: iType = TypeOF(ww)
CASE iType OF
	0   :	weight = 1.
	1   :   weight = ww
	2   :   weight = ww
	3   :   weight = ww
	4   :   weight = ww
	5   :   weight = ww
	6   :   weight = ww
	7   :   weight = ww
	8   :   BEGIN
		MESSAGE, "Weight cannot be a structure in this context", /INFORMATIONAL
		RETURN, 0
		END
	9   :   weight = ww
	10  :   BEGIN
		weight = *weight
		GOTO, TryAgain
		END
	11  :   BEGIN
			wObj = ww -> INTERPOLATE(template)
			wptr = wObj -> GetValues()
			weight = *wptr
			wObj -> Trash
			PTR_FREE, wptr
		END
	
	12  :   weight = ww
	13  :   weight = ww
	14  :   weight = ww
	15  :   weight = ww
	ELSE:	BEGIN
		MESSAGE, "Weight has illegal type", /INFORMATIONAL
		RETURN, 0
		END
ENDCASE

values = *self.values
weight = MAKE_ARRAY( (iMax-iMin+1), VALUE=1., /FLOAT)#weight
values = values*weight
values = TOTAL(values, 2)
vMin = GKVsd_Min(values, MAX=vMax)
outputStr.values = PTR_NEW(values)
outputStr.vrange = [vMin, vMax]
output = OBJ_NEW("GKVs1D", outputStr)
RETURN, output
END  ; ****** GKVs2D::PROJECT1.pro ****** ;


FUNCTION GKVs3D::PROJECT1, _Extra=extra
;
; Purpose:
;
;	Projects 'self' onto first independent variable
;	using specifiec weighting function.
;
; Written by W.M. Nevins
;	4/9/04
;
outputStr = {GKVs1D}
nTags = N_TAGS({GKVsd})
FOR i=0, nTags-1 DO outputStr.(i) = self.(i)
outputStr.Indices = PTR_NEW(['*'])
outputStr.Grid1 = GKVsd_GridCopy(self.grid1)
iMin = self.Grid1.irange[0]
iMax = self.Grid1.irange[1]
template = self -> slice(axis=1, index=iMin)
;
; compute weight function
;
weightObj = template -> ProjectWeight3(_Extra=extra)
weightPtr = weightObj -> GetValues()
weight = *WeightPtr
norm = 1./TOTAL(weight)
PTR_FREE, weightPtr
weightObj -> Trash
template  -> trash
;
; Perform projection
;
values = *self.values
type = TypeOf(values)
projection = MAKE_ARRAY(iMax-iMin+1, TYPE=type)
FOR i=imin,imax DO projection[i] = TOTAL(values[i,*,*]*weight)*norm
vMin = GKVsd_Min(projection, MAX=vMax)
outputStr.values = PTR_NEW(projection)
outputStr.vrange = [vMin, vMax]
output = OBJ_NEW("GKVs1D", outputStr)
RETURN, output
END  ; ****** GKVs3D::PROJECT1.pro ****** ;


FUNCTION GKVs2D::PROJECT2, WEIGHT=ww
;
; Purpose:
;
;	Projects 'self' onto second independent variable
;	using specifiec weighting function.
;
; Written by W.M. Nevins
;	4/9/04
;
outputStr = {GKVs1D}
nTags = N_TAGS({GKVsd})
FOR i=0, nTags-1 DO outputStr.(i) = self.(i)
outputStr.Indices = PTR_NEW(['*'])
outputStr.Grid1 = GKVsd_GridCopy(self.grid2)
iMin = self.Grid2.irange[0]
iMax = self.Grid2.irange[1]
template = self -> slice(axis=2, index=iMin)

TryAgain: iType = TypeOF(ww)
CASE iType OF
	0   :	weight = 1.
	1   :   weight = ww
	2   :   weight = ww
	3   :   weight = ww
	4   :   weight = ww
	5   :   weight = ww
	6   :   weight = ww
	7   :   weight = ww
	8   :   BEGIN
		MESSAGE, "Weight cannot be a structure in this context", /INFORMATIONAL
		RETURN, 0
		END
	9   :   weight = ww
	10  :   BEGIN
		weight = *weight
		GOTO, TryAgain
		END
	11  :   BEGIN
			wObj = ww -> INTERPOLATE(template)
			wptr = wObj -> GetValues()
			weight = *wptr
			wObj -> Trash
			PTR_FREE, wptr
		END
	
	12  :   weight = ww
	13  :   weight = ww
	14  :   weight = ww
	15  :   weight = ww
	ELSE:	BEGIN
		MESSAGE, "Weight has illegal type", /INFORMATIONAL
		RETURN, 0
		END
ENDCASE

values = *self.values
weight = weight#MAKE_ARRAY( (iMax-iMin+1), VALUE=1., /FLOAT)
values = values*weight
values = TOTAL(values, 1)
vMin = GKVsd_Min(values, MAX=vMax)
outputStr.values = PTR_NEW(values)
outputStr.vrange = [vMin, vMax]
output = OBJ_NEW("GKVs1D", outputStr)
RETURN, output
END  ; ****** GKVs2D::PROJECT2.pro ****** ;


FUNCTION GKVs3D::PROJECT2, _Extra=extra
;
; Purpose:
;
;	Projects 'self' onto first independent variable
;	using specifiec weighting function.
;
; Written by W.M. Nevins
;	4/9/04
;
outputStr = {GKVs1D}
nTags = N_TAGS({GKVsd})
FOR i=0, nTags-1 DO outputStr.(i) = self.(i)
outputStr.Indices = PTR_NEW(['*'])
outputStr.Grid1 = GKVsd_GridCopy(self.grid2)
iMin = self.Grid2.irange[0]
iMax = self.Grid2.irange[1]
template = self -> slice(axis=2, index=iMin)
;
; compute weight function
;
weightObj = template -> ProjectWeight3(_Extra=extra)
weightPtr = weightObj -> GetValues()
weight = *WeightPtr
norm = 1./TOTAL(weight)
PTR_FREE, weightPtr
weightObj -> Trash
template  -> trash
;
; Perform projection
;
values = *self.values
type = TypeOf(values)
projection = MAKE_ARRAY(iMax-iMin+1, TYPE=type)
FOR i=imin,imax DO projection[i] = TOTAL(values[*,i,*]*weight)*norm
vMin = GKVsd_Min(projection, MAX=vMax)
outputStr.values = PTR_NEW(projection)
outputStr.vrange = [vMin, vMax]
output = OBJ_NEW("GKVs1D", outputStr)
RETURN, output
END  ; ****** GKVs3D::PROJECT2.pro ****** ;


FUNCTION GKVs3D::PROJECT3, _Extra=extra
;
; Purpose:
;
;	Projects 'self' onto first independent variable
;	using specifiec weighting function.
;
; Written by W.M. Nevins
;	4/9/04
;
outputStr = {GKVs1D}
nTags = N_TAGS({GKVsd})
FOR i=0, nTags-1 DO outputStr.(i) = self.(i)
outputStr.Indices = PTR_NEW(['*'])
outputStr.Grid1 = GKVsd_GridCopy(self.grid3)
iMin = self.Grid3.irange[0]
iMax = self.Grid3.irange[1]
template = self -> slice(axis=3, index=iMin)
;
; compute weight function
;
weightObj = template -> ProjectWeight3(_Extra=extra)
weightPtr = weightObj -> GetValues()
weight = *WeightPtr
norm = 1./TOTAL(weight)
PTR_FREE, weightPtr
weightObj -> Trash
template  -> trash
;
; Perform projection
;
values = *self.values
type = TypeOf(values)
projection = MAKE_ARRAY(iMax-iMin+1, TYPE=type)
FOR i=imin,imax DO projection[i] = TOTAL(values[*,*,i]*weight)*norm
vMin = GKVsd_Min(projection, MAX=vMax)
outputStr.values = PTR_NEW(projection)
outputStr.vrange = [vMin, vMax]
output = OBJ_NEW("GKVs1D", outputStr)
RETURN, output
END  ; ****** GKVs3D::PROJECT3.pro ****** ;
