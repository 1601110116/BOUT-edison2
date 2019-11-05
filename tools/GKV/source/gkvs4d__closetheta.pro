PRO GKVs2D::YtoYtilde, theta=theta_in, sHat = shat_in, Params=params
;
; Purpose:
; 
; This proceedure acts on (x,y,t data,
; Removing the  field-line shift to 
; form a locally-orthogonal coordinate
; system at locations off of outboard 
; midplane (theta=0).
;
; Keywords:
;
;  Theta    The value of theta for this data
;
; Params  The "Params" structure returned
;         by GENE_DATA (or multipleGENE_Data)
;         (Optional).
;
; Written by W.M. Nevins
;  5/25/09
;
theta = !PI
IF(N_ELEMENTS(theta_in) EQ 1) THEN theta = theta_in
sHat = 0.786
IF(N_ELEMENTS(sHat_in) EQ 1) THEN sHat = sHat_in
;
; Check for 'params' structure
; if it exists, get defaults values
; from it
;
IF(TypeOF(params) EQ 8) THEN BEGIN
  thisParams = params
  IF(N_ELEMENTS(params) GT 1) THEN thisParams=params[0]
  ky0 = thisParams.box.kymin
  Lx  = thisParams.box.Lx
  kx0 = 2.*!PI/Lx
  nTags = N_TAGS(thisParams)
  paramTags = TAG_NAMES(thisParams)
  FOR i=0, nTags-1 DO $
    IF(STRCMP(paramTags[i], "geometry", /FOLD_CASE)) THEN sHat = thisParams.geometry.sHat
  FOR i=0, nTags-1 DO $
    IF(STRCMP(paramTags[i], "global",   /FOLD_CASE)) THEN sHat = thisParams.Global.sHat
ENDIF
values = *self.values
info = SIZE(values)
nDims = info[0]
;
; This routine won't work on 4D objects ...
;
IF(nDims GT 3) THEN BEGIN
  MESSAGE, "YtoYtidle not implimented for 4D objects", /INFORMATIONAL
  RETURN
ENDIF
result = MAKE_ARRAY(SIZE=info)
xValues = *self.grid1.values
nx = N_ELEMENTS(xValues)
yValues = *self.grid2.values
dy = yValues[1] - yValues[0]
FOR i=0,nx-1 DO BEGIN
  x = xValues[i]
  thisShift = -theta*shat*x/dy
  sgn = FIX(thisShift GE 0) - FIX(thisShift LT 0)
  iShift = FIX(thisShift)
  dShift = sgn*(thisShift - iShift)
  CASE nDims OF
        2  :  result[i,*]   = (1.-dShift)*SHIFT(values[i,*],0,thisShift) $
                             +dShift*SHIFT(values[i,*],0,thisShift+sgn)
        3  :  result[i,*,*] = (1.-dShift)*SHIFT(values[i,*,*],0,thisShift,0) $
                             + dShift*SHIFT(values[i,*,*],0,thisShift+sgn,0)
  ENDCASE
ENDFOR
PTR_FREE, self.values
self.values = PTR_NEW(result)
yTitle = "!Sy!R!U!9A!X!N"
yMnemonic = 'yTilde'
self.grid2.title = yTitle
self.grid2.mnemonic = yMnemonic
RETURN
END ; ****** GKVs3D::YtoYtilde ****** ;

FUNCTION GKVs2D::YtoYtilde, _Extra=Extra
result = self -> MakeCopy()
result -> YtoYtilde, _EXTRA=Extra
RETURN, result
END ; ****** FUNCTION GKVs2D::YtoYtilde ****** ;



PRO GKVs3D::CloseTheta_k, kx0 = kx0In, ky0=ky0In, sHat=sHatIn, params=params
;
; Purpose:
;
; This routine acts on GENE field data in
; (k_x, k_y, theta, t) format and, for each
; (k_x, k_y, t)  "closes" it in theta by
; providing the (redundant) value of the field
; at theta = +pi.
;
; Keywords:
; 
; kx0     The fundamental mode in kx.
;         Defaults to 0.0617324
;         (Optional).
; 
; ky0     The fundamental mode in ky.
;         Defaults to 0.05
;         (Optional)
;
; sHat    The value of the magnetic shear.
;         Defaults to 0.786 (Optional).
;
; Params  The "Params" structure returned
;         by GENE_DATA (or multipleGENE_Data)
;         (Optional).
;
; Written by W.M. Nevins
;  5/24/09
;
; Set default shear
;
kx0 = 0.0617324
ky0 = 0.05
sHat = 0.786
;
; Check for 'params' structure
; if it exists, get defaults values
; from it
;
IF(TypeOF(params) EQ 8) THEN BEGIN
  thisParams = params
  IF(N_ELEMENTS(params) GT 1) THEN thisParams=params[0]
  ky0 = thisParams.box.kymin
  Lx  = thisParams.box.Lx
  kx0 = 2.*!PI/Lx
  nTags = N_TAGS(thisParams)
  paramTags = TAG_NAMES(thisParams)
  FOR i=0, nTags-1 DO $
    IF(STRCMP(paramTags[i], "geometry", /FOLD_CASE)) THEN sHat = thisParams.geometry.sHat
  FOR i=0, nTags-1 DO $
    IF(STRCMP(paramTags[i], "global",   /FOLD_CASE)) THEN sHat = thisParams.Global.sHat
ENDIF
;
; Parse command line
;
IF(N_ELEMENTS(ky0In)  EQ 1) THEN ky0  = ky0In
IF(N_ELEMENTS(kx0In)  EQ 1) THEN kx0  = kx0In 
IF(N_ELEMENTS(sHatIn) EQ 1) THEN sHat = sHatIn
;
; Get theta grid
; and append +pi
;
thetaGrid = self.grid3
thetaMax = thetaGrid.range[1]
IF(ABS(thetaMax-!PI) LT 1.e-6) THEN BEGIN
  MESSAGE, "Theta grid was closed, Returning", /INFORMATIONAL
  RETURN
ENDIF
oldThetaValues = *thetaGrid.values
nThetas = N_ELEMENTS(oldThetaValues)
newThetaValues = FLTARR(nThetas+1)
newThetaValues[0:(nThetas-1)] = oldThetaValues
newThetaValues[nThetas] = !PI
PTR_FREE, thetaGrid.values
thetaGrid.values = PTR_NEW(newThetaValues)
thetaGrid.range = [-!PI,!PI]
thetaGrid.irange = [0, nThetas]
thetaGrid.boundary = "closed"
self.Grid3=thetaGrid
;
; Make an array to hold
; new field values
;
oldValues = *self.values
info = SIZE(oldValues)
nDims = info[0]
;info[nDims+2] = (info[nDims+2]/info[3])(info[3]+1)
info[3] = info[3]+1
newValues = MAKE_ARRAY(SIZE=info)
;
; load old field values
;
CASE nDims OF
      3  :  newValues[*,*,0:(nThetas-1)]   = oldValues
      4  :  newValues[*,*,0:(nThetas-1),*] = oldValues
ENDCASE
;
; Now loop over ky values to compute 
; field values at theta = +pi
Delta_kxky = 2.*!PI*sHat*ky0/kx0
Delta_kxky = LONG(Delta_kxky + 0.4999999999)
nkx = info[1]
nky = info[2]
FOR j=0,nky-1 DO BEGIN
  Delta_i = -j*Delta_kxky
  test = FIX(Delta_i GT 0) - FIX(Delta_i LT 0)
  CASE test OF
        1  :  BEGIN
                iMin  = Delta_i
                iMax  = nkx-1
                iiMin = 0
                iiMax = nkx-1-Delta_i
              END ; newValues[Delta_i:(nkx-1),*,nThetas]   = oldValues[0:(nkx-1-Delta_i),*,0]
        0  :  BEGIN
                iMin  = 0
                iMax  = nkx-1
                iiMin = 0
                iiMax = nkx-1
              END ; newValues[Delta_i:(nkx-1),*,nThetas]   = oldValues[0:(nkx-1),        *,0]
       -1  :  BEGIN
                iMin  = 0
                iMax  = nkx-1+Delta_i
                iiMin = -Delta_i
                iiMax = nkx-1
              END ;  newValues[0:(nkx-1+Delta_i),*,nThetas] = oldValues[-Delta_i:(nkx-1), *,0]
  ENDCASE
   
  CASE nDims OF
      3  :  newValues[iMin:iMax,j,nThetas]   = oldValues[iiMin:iiMax,j,0]
      4  :  newValues[iMin:iMax,j,nThetas,*] = oldValues[iiMin:iiMax,j,0,*]
  ENDCASE
ENDFOR
;
; load new values into 'self'
;
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END ; ****** PRO GKVs3D::CloseTheta_k ****** ;

FUNCTION GKVs3D::CloseTheta_k, _EXTRA=Extra
result = self -> MakeCopy()
result -> CloseTheta_k, _Extra=Extra
RETURN, result
END ; ****** FUNCTION GKVs3D::CloseTheta_k ****** ;


FUNCTION GKVs3D::ExtendTheta, kx0=kx0In, ky0=ky0In, sHat=sHatIn, params=params
;
; Purpose:
;
; This routine acts on GENE field data at
; a fixed k_y in (k_x, theta, t) format
; with theta on the interval [-!Pi,!PI).
; 
; It returns an object containing the field
; on the extended (out to available k_x resolution)
; theta grid.
;
; Keywords:
;
;  kx0    Selects desired value of k_x in the fundamental
;         zone in theta [-!PI,!PI]. Defaults to kx0=0.
;         (Optional).
;         
;  ky0    For 3D objects it is necessary to enter selected
;         value of k_y . If 'self' is a 4D object,
;         then selects this value of k_y. 
;         Defaults to 0.05 (Optional).
;        
;  sHat   The value of the magnetic shear.
;         Defaults to 0.786 (Optional).
;
; Written by W.M. Nevins
;    5/24/09
;
; set defaults
;
kx0 = 0.
ky0 = 0.058
shat = 0.786
;
; Check for 'params' structure
; if it exists, get defaults values
; from it
;
IF(TypeOF(params) EQ 8) THEN BEGIN
  thisParams = params
  IF(N_ELEMENTS(params) GT 1) THEN thisParams=params[0]
  ky0 = thisParams.box.kymin
  sHat = thisParams.geometry.sHat
ENDIF

IF(N_ELEMENTS(kx0In) EQ 1) THEN kx0 = kx0In

ky0 = 0.05
IF(N_ELEMENTS(ky0In) EQ 1) THEN ky0 = ky0In

sHat = 0.786
IF(N_ELEMENTS(sHatIn) EQ 1) THEN sHat = sHatIn
;
; If 'self' is a 4D object, slice at selected value of ky
; 
N_Dims = self -> NumDims()
CASE N_Dims OF
      3  :  thisObj = self -> MakeCopy()
      4  :  BEGIN
      
            thisObj = self -> Slice(k_y=ky0)
            END
ENDCASE
;
; Compute stride in k_x
;
Delta_kx = 2.*!PI*sHat*ky0
;
; Compute theta values
; on the extended theta grid
;
kxRange = thisObj.Grid1.range
iRange = FIX(kxRange/Delta_kx)
nCycles = 1 + iRange[1] - iRange[0]
theta0 = *(thisObj.Grid2.values)
nTheta0 = N_ELEMENTS(theta0)
newThetas = FLTARR(nCycles*nTheta0)
FOR i=iRange[0], iRange[1] DO BEGIN
  jMin = (i-iRange[0])*nTheta0
  jMax = jMin + nTheta0 - 1
  newThetas[jmin:jMax] = theta0 + 2.*!PI*i
ENDFOR
;
; Compute field values on
; the extended theta grid
;
OldValues = *(thisObj.values)
info = SIZE(OldValues)
nt = info[3]
newValues = COMPLEXARR(nCycles*nTheta0, nt)
FOR i=iRange[0], iRange[1] DO BEGIN
  jMin = (i-iRange[0])*nTheta0
  jMax = jMin + nTheta0 - 1
  thisCycle = thisObj -> Slice(k_x=kx0 + i*Delta_kx)
  thisValuePtr = thisCycle -> GetValues()
  newValues[jmin:jMax, *] = *thisValuePtr
  PTR_FREE, thisValuePtr
  thisCycle -> Trash
ENDFOR
;
; Make output object
;
thisObj -> Set, axis=1, GridTitle = "k!Dx0!N", GridMnemonic='k_x0'
output = thisObj -> slice(k_x0=kx0)
output -> Get, values=valuePtr
output -> Get, axis=1, GridValues=GridValuePtr
vMin = GKVsd_Min(newValues, MAX=vMax)
PTR_FREE, valuePtr
output -> Set, values=PTR_NEW(newValues), vrange=[vMin, vMax]
PTR_FREE, GridValuePtr
thetaMin = GKVsd_Min(newThetas, MAX=thetaMax)
iMax = N_ELEMENTS(newThetas) - 1
output -> Set, axis=1, GridValues=PTR_NEW(newThetas)
output -> Set, axis=1, range=[thetaMin, thetaMax], irange=[0, iMax]
thisObj -> Trash
RETURN, output
END ; ***** FUNCTION GKVs3D::ExtendTheta ****** ;


FUNCTION GKVs4D::CloseTheta, _Extra=Extra
result = self -> MakeCopy()
result -> CloseTheta, _Extra=Extra
RETURN, result
END ; ****** FUNCTION GKVs4D::CloseTheta ****** ;


PRO GKVs4D::CloseTheta, _Extra=Extra
;
; Purpose:
;
;  This proceedure acts on 3D plus time GENE data
;  in (x,y,theta,t) format. It "closes" the theta-
;  grid by providing data at theta=pi.  The data at
;  theta=pi is obtained by applying the shear and
;  periodicity rules.
;
;  Keywords:
;
;  Params  Use this keyword to pass the "Params" structure
;          greated by Gene_Data from the corressponding
;          GENE simulation.
;
;  Written by W.M. Nevins
;    5/22/09
;
;ky0  = 0.058
;Lx   = 212.044105
;Shat = 0.786
;Nx   = 256
;Nky  = 16
result = GetKeyWord("params", Extra)
IF(TypeOf(Result) EQ 8) THEN BEGIN
  ky0  = (result.box.kymin)[0]
  Lx   = (result.box.lx)[0]
  Shat = (result.geometry.sHat)[0]
  Nx   = (result.box.nx0)[0]
  Nky  = (result.box.nky0)[0]
ENDIF
result = GetKeyWord("ky0", Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN ky0=result
result = GetKeyWord("Lx", Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN Lx=result
result = GetKeyWord("sHat", Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN sHat=result
result = GetKeyWord("Nx", Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN Nx=result
result = GetKeyWord("Nky", Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN Nky=result

;
; Check if theta grid requires "completion"
;
oldthetaValues = *(self.Grid3.values)
nThetas = N_ELEMENTS(oldThetaValues)
oldThetaMax = oldThetaValues[nThetas-1]
IF( ABS(oldThetaMax - !PI) LT 1.e-6) THEN RETURN
;
; Maximum value of theta was not pi, so grid does need completion
;
Dx = Lx/(Nx)
Ly = 2.*!PI/ky0
Dy = Ly/(2*Nky)
DjDi = 2.*!PI*Shat*Dx/Dy
DjDi = FIX(DjDi + 0.49999999)
oldValues = *self.values
size = SIZE(oldValues)
newSize =size
newSize[3] = Size[3] + 1
newValues = MAKE_ARRAY(SIZE=newSize)
newValues[0:(size[1]-1),0:(size[2]-1),0:(size[3]-1),0:(size[4]-1)] = oldValues
sgn = 1
IF (DjDi LT 0) THEN sgn = -1
FOR i=0L,Size[1]-1 DO BEGIN
  ValuesAtPi = REFORM(oldValues[i,*,0,*])
  ValuesAtPi = SHIFT(ValuesAtPi,DjDi*i,0)
  newValues[i,*, size[3], *] = ValuesAtPi
ENDFOR
PTR_FREE, self.values
self.Values = PTR_NEW(newValues)
;
; add theta= pi point to theta grid
;
newThetas = FLTARR(nThetas + 1)
newThetas[0:nThetas-1] = oldThetavalues
newThetas[nThetas] = !PI
PTR_FREE, self.Grid3.values
self.Grid3.values = PTR_NEW(newThetas)
self.Grid3.range  = [-!PI,!PI]
self.Grid3.irange = [0,nThetas]
RETURN
END ; ****** GKVs4D::CloseTheta ****** ; 


