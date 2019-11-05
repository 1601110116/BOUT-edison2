FUNCTION GKVs1D::drSq, Lx=Lxin
;
; Computes estimate of magnetic diffusion coefficient
; from single field map
;
range = self.grid1.range
Lx=range[1]-range[0]
IF(N_ELEMENTS(Lxin) NE 0) THEN Lx=Lxin
xValues = *self.Grid1.values
;
; Try to remove peroidicity
;
nValues = N_ELEMENTS(xValues)
FOR i=1L, nValues-1 DO BEGIN
  dx = xValues[i] - xValues[i-1]
  WHILE(dx GT  Lx/2.) DO BEGIN
    xValues[i:nValues-1] = xValues[i:nValues-1] - Lx
    dx = xValues[i] - xValues[i-1]
  ENDWHILE
  WHILE(dx LT -Lx/2.) DO BEGIN
    xValues[i:nValues-1] = xValues[i:nValues-1] + Lx
    dx = xValues[i] - xValues[i-1]
  ENDWHILE
ENDFOR
;
; Compute change in radial position, etc.
;
dx = xValues - xValues[0]
dxSq = dx*dx
output = self -> MakeCopy(/NoValues)
output.title = "!4d!Xr!U2!N"
output.mnemonic="drSq"
output.values=PTR_NEW(dxSq)
vMax = MAX(dxSq)
output.vrange=[0,vMax]
PTR_FREE, output.indices
output.indices = PTR_NEW('*')
gridValues = FINDGEN(nValues)
PTR_FREE, output.grid1.values
output.grid1.values = PTR_NEW(gridValues)
output.grid1.range = [0, nValues-1]
output.grid1.title = "n_cycles"
output.grid1.mnemonic = "n_cycles"
output.grid1.units = ""
RETURN, output

END ; ****** drSq ****** ;

FUNCTION Dmag, mapArray
nObjs = N_ELEMENTS(mapArray)
drSqArr = OBJARR(nObjs)
mapArrayInfo = SIZE(mapArray)
imax = mapArrayInfo[1] - 1L
jmax = mapArrayInfo[2] - 1L
ii=0L
FOR i=0L,imax DO BEGIN
  FOR j=0L,jmax DO BEGIN
    drSqArr[ii] = mapArray[i,j] -> drSq() ;maparray[i,j]
    ii = ii+1L
  ENDFOR
ENDFOR
drSqAvg = FLTARR(3000)
thisN   = FLTARR(3000)
FOR ii=0L,nObjs-1 DO BEGIN
  drSqArr[ii] -> Get, values = valuePtr 
  theseValues = *valuePtr
  nValues = N_ELEMENTS(theseValues)
  drSqAvg[0:nValues-1] = drSqAvg[0:nValues-1] + theseValues
  thisN[0:nValues-1] = thisN[0:nValues-1] + 1.
ENDFOR
drSqAvg = drSqAvg/thisN
nCycles = FINDGEN(3000)
nCycles[0] = 1.
dMag = drSqAvg/nCycles
errors = dMag/(thisN-1)
output = drSqArr[0] -> MakeCopy(/NoValues)
output -> set, title = "!13<!4d!Xr!U2!N!13>!X/nCycles"
output -> set, mnemonic = "dMag"
output -> set, values = PTR_NEW(dMag)
output -> set, errorbars =PTR_NEW(errors)
output -> get, axis=1, gridValues = oldptr
PTR_FREE, oldptr
output -> set, axis=1, gridvalues = PTR_NEW(FINDGEN(3000.))
output -> set, axis=1, irange=[0,2999]
output -> set, axis=1,  range=[0,2999.]
gkvdelete, drSqArr
RETURN, output
END ; ****** Dmag ****** ;


PRO FieldMapCleanUP
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
dataObj -> TRASH
RETURN
END ; ****** FieldMapCleanUP ****** ;

FUNCTION oneA_parallel, x, y
;
; Purpose:
;  
;  This function returns the derivatives of
;  A_parallel from single fourier mode in y-direction
;  for purposes of testing mapping routine
; 
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
thisObj = dataObj -> slice(ky=ky0)
thisObj -> get, values = valuePtr
values = *valuePtr
thisObj -> get, axis=1, gridValues=kxptr
kxValues = *kxPtr
thisObj -> Trash
ky=ky0
eye = COMPLEX(0.,1.)
expky = EXP(eye*ky0*y)
expkx = EXP(eye*kxValues*x)
expMatrix = expkx*expky
grad_x = eye*kxValues
grad_y = eye*ky0
appl_x = TOTAL(grad_x*Values*expMatrix)
appl_y = TOTAL(grad_y*Values*expMatrix)
appl_value = TOTAL(values*expMatrix)
appl_x = 2.*FLOAT(appl_x)
appl_y = 2.*FLOAT(appl_y)
appl_value = 2.*FLOAT(appl_value)
result = {  value :  appl_value,  $
               x  :  appl_x,      $
               y  :  appl_y       }
RETURN, result
END ; ****** oneA_parallel ****** ;



FUNCTION numA_parallel, x, y
;
; Purpose:
;
;  This Function returns the derivatives of 
;  A_parallel computed from the fourier
;  representation located in dataObj
;
; Arguments:
;  
;    x  Radial location at which the derivatives of
;       A_parallel are to be evaluated.
;
;    y  Bi-normal location at which the derivatives of 
;       A_parallel are to be evaluated.
;
;  Written by W.M. Nevins
;    5/29/09
;
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
dataObj -> Get, values = ValuePtr
dataObj -> Get, axis = 1, gridValues = kxPtr
dataObj -> Get, axis = 2, gridValues = kyPtr
Values = *valuePtr
kxValues = *kxPtr
kyValues = *kyptr
nx = N_ELEMENTS(kxValues)
ny = N_ELEMENTS(kyValues)
eye = COMPLEX(0.,1.)
expkx = COMPLEXARR(nx)
kx_0 = kxValues[1] - kxValues[0]
expDkx = EXP(eye*kx_0*x)
expkx[0] = EXP(eye*kxValues[0]*x)
FOR i=1,nx-1 DO expkx[i] = expkx[i-1]*expDkx
expky = COMPLEXARR(ny)
ky_0 = kyValues[1] - kyValues[0]
expDky = EXP(eye*ky_0*y)
expky[0] = EXP(eye*kyValues[0]*y)
FOR i=1,ny-1 DO expky[i] = expky[i-1]*expDky
expMatrix = expkx#expky
kValues = expMatrix*Values
grad_x = eye*kxValues#MAKE_ARRAY(ny,VALUE=1.)
grad_y = eye*MAKE_ARRAY(nx, VALUE=1.)#kyValues
appl_value = TOTAL(kValues)
appl_x = TOTAL(grad_x*kValues)
appl_y = TOTAL(grad_y*kValues)
;
; GENE and GS2 do not supply fourier
; coefficeints at negative values of ky.
; We use reality condition to get this 
; information

IF(kyValues[0] EQ 0) THEN BEGIN
  appl_value = appl_value + CONJ(TOTAL(kvalues[*,1:ny-1]))
  appl_x = appl_x + CONJ(TOTAL(grad_x[*,1:ny-1]*kvalues[*,1:ny-1]))
  appl_y = appl_y + CONJ(TOTAL(grad_y[*,1:ny-1]*kvalues[*,1:ny-1]))
ENDIF ELSE BEGIN
  appl_value = 2.*appl_value
  appl_x = 2.*appl_x
  appl_y = 2.*appl_y
ENDELSE

IF(N_ELEMENTS(A0) EQ 1)THEN BEGIN
  appl_value = A0*appl_value
  appl_x = A0*appl_x
  appl_y = a0*appl_y
ENDIF
appl_value = FLOAT(appl_value)
appl_x = FLOAT(appl_x)
appl_y = FLOAT(appl_y)

result = { value :  appl_value,  $
              x  :  appl_x,      $
              y  :  appl_y       }
RETURN, result
END ; ****** FUNCTION numA_parallel ****** ;

FUNCTION FieldMap_SampleAparallel, x, y
;
; Purpose
;
;  This function returns the derivatives of 
;  an analytic version of A_parallel to be
;  used in debugging our field line mapping 
;  proceedure.
;
;  Arguments:
;  
;    x  Radial location at which the derivatives of
;       A_parallel are to be evaluated.
;
;    y  Bi-normal location at which the derivatives of 
;       A_parallel are to be evaluated.
;
;  Written by W.M. Nevins
;    5/29/09
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
;
; We will take A_parallel = A0*COS(kx0*x)*SIN(ky0*y)
;
Aparallel_x = -kx0*A0*SIN(kx0*x)*SIN(ky0*y)
Aparallel_y =  ky0*A0*COS(kx0*x)*COS(ky0*y)
result = {  x  : Aparallel_x,  $
            y  : Aparallel_y   }
RETURN, Result
END ; ****** FUNCTION FieldMap_SampleAparallel ****** ;

Function FieldMap_Resid, z
;
; Purpose:
;
;  This routine computes the residual required for the 
;  implicit field line mapping routine.
;
;  Argument:
;
;    z     A two element array containing the values of 
;          the radial coordinate (z[0]) and the bi-normal
;          coordinate (z[1]) at which the residual will
;          be evaluated.
;
;  Written by W.M. Nevins
;    5/29/09
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
resid = FLTARR(2)
xnew = z[0]
ynew = z[1]
;
; Compute derivatives of A_parallel 
; at desired location. The function
; whose name is contained in the
; variable A_parallel must return
; a structure with the x and y tags
; corresponding to the x and y derivatives
; of A_parallel at this location.
;
CommandLine = "Appl = " + A_parallel + "(xnew, yold)"
ok = EXECUTE(CommandLine)
resid[0] = xnew - (xold + 2.*!PI*qCommon*Appl.y)
resid[1] = ynew - (yold - 2.*!PI*qCommon*Appl.x + 2.*!PI*sHatCommon*xnew)
RETURN, resid
END ; ****** Function FieldMap_Resid ****** ;



FUNCTION FieldMap, params=params, x0=x0_in, y0=y0_in, A_parallel=A_parallel_in, $
              oPlot = oPlot, pSym = pSym_in, nx=nx, ny=ny, A0= A0_in,      $
              ApplObj=ApplObj, Aout=Aout, _Extra=Extra
;
; Purpose:
;
;  This proceedure iterates a field-ine mapping
;  to test implicit area-preserving mapping alogrithm
;
; Keywords:
;
;  Params  The 'params' structure created by Gene_Data or MultipleGene_Data
;
;  x0      Initial x-value for this map
;
;  y0      Initial y-value for this map
;
; Written by W.M. Nevins
;  5/28/09
;
; Check for 'params' structure
; if it exists, get defaults values
; from it
;
COMMON FieldMap_Common1, xold, yold, qCommon, sHatCommon, A_parallel, A0, kx0, ky0, dataObj
ly = 74.1986
ky0 = 0.0846807
shat = 0.8
n_rats = 4
lx = n_rats*ly/(2.*!PI*shat)
kx0 = 2.*!PI/lx
q0 = 1.4
IF(TypeOF(params) EQ 8) THEN BEGIN
  thisParams = params
  IF(N_ELEMENTS(params) GT 1) THEN thisParams=params[0]
  ky0 = thisParams.box.kymin
  ly  = 2.*!PI/ky0
  Lx  = thisParams.box.Lx 
  sHat = thisParams.global.sHat
  q0 = thisParams.global.q0
ENDIF
qCommon = q0
sHatCommon = sHat
;
; Parse command line
;
x0=Lx/SQRT(2.)
IF(N_ELEMENTS(x0_in) EQ 1) THEN x0 = x0_in
y0 = !PI
IF(N_ELEMENTS(y0_in) EQ 1) THEN y0 = y0_in
A_parallel = "FieldMap_SampleAparallel"
 IF(N_ELEMENTS(nx) EQ 1) THEN kx0 = nx*kx0
 IF(N_ELEMENTS(ny) EQ 1) THEN ky0 = ny*ky0
 A0 = 1.
 IF(N_ELEMENTS(A0_in) EQ 1) THEN A0 = A0_in
 IF( TypeOF(applObj) EQ 11) THEN BEGIN
   IF( OBJ_ISA(applObj, "GKVs2D") ) THEN BEGIN
    dataObj = applObj -> MakeCopy()
    A_parallel = "numA_parallel"
   ENDIF
 ENDIF
IF(TypeOf(A_parallel_in) EQ 7) THEN A_parallel = A_parallel_in

;
; Set up arrays to hold results of the mapping
;
maxIterations = 100
result = GetKeyWord("maxIterations", Extra) 
IF(Query_Integer(result)) THEN maxIterations = result
x = FLTARR(maxIterations)
y = FLTARR(maxIterations)
Ax = FLTARR(maxIterations)
Ay = FLTARR(maxIterations)
nIterations = maxIterations
x[0] = x0
y[0] = y0
 ;
 ; Compute guess at result of this mapping
 ;
 const = 2.*!PI*q0
 ;
 ; set error trap
 ;
 CATCH, Error_Status
 IF(Error_Status NE 0) THEN BEGIN
  CATCH, /CANCEL
  MESSAGE, !ERROR_STATE.MSG, /INFORMATIONAL
  MESSAGE, "Stopped mapping after " + STRCOMPRESS(STRING(i),/REMOVE_ALL) + " iterations.", /INFORMATIONAL
  IF(i EQ 0) THEN BEGIN
  ; issue here is that GKVs1D objects
  ; won't "draw" if there is only one
  ; data point, so we add a second one.
    x[1] = x[0]
    y[1] = y[0]
    Ax[1]=Ax[0]
    Ay[1]=Ay[0]
       i = 1
  ENDIF
  x=x[0:i]
  y=y[0:i]
  nIterations = i
  GOTO, plotMap
 ENDIF
FOR i=0,maxIterations-2 DO BEGIN
  xold = x[i]
  yold = y[i]
  commandLine = "Appl = " + A_parallel + "(xold, yold)"
  OK = EXECUTE(commandLine)
  IF(NOT OK) THEN MESSAGE, "Failure in computation of A_parallel"
  Ax[i] = Appl.x
  Ay[i] = Appl.y
  xguess = xold + const*Appl.y
  yguess = yold - const*Appl.x + 2.*!PI*shat*xguess
  zguess = [xguess, yguess]
  z = NEWTON(zguess, "FieldMap_Resid")
  x[i+1] = z[0] MOD lx
  IF(x[i+1] LT 0.) THEN x[i+1] = x[i+1] + lx
  y[i+1] = z[1] MOD ly
  IF(y[i+1] LT 0.) THEN y[i+1] = y[i+1] + ly
ENDFOR
 ;
 ; Plot result of the mapping
 ; 
 PlotMap: psym = 1
;IF(N_ELEMENTS(psym_in) EQ 1) THEN psym=pSym_in
;IF(N_ELEMENTS(oPlot) EQ 0) THEN oPlot = 0
;IF(oPlot EQ 0) THEN BEGIN
;  PLOT, x, y, xrange=[0,lx], yrange=[0,ly], psym=psym, $
;        xTitle = 'x', ytitle='y', xStyle=1, yStyle=1
;ENDIF ELSE BEGIN
;  oPlot, x, y, color=oPlot, psym=psym
;ENDELSE
;
; Turn mapping into GKVs1D object
;
codeName = "GKV"
codePI   = "W.M. Nevins"
FileID   = A_parallel
runID    = ""
indices  = ["*", "y"]
IF(TypeOF(dataObj) EQ 11) THEN BEGIN
  dataObj -> Get, codeName=codeName
  dataObj -> Get, CodePI=codePi
  dataObj -> Get, FileID=FileID
  DataObj -> Get, RunID = RunID
  DataObj -> Get, indices=indicesPtr
  indices = *indicesPtr
  indices[1] = "y"
  DataObj -> Get, axis=2, range=kyRange
  kystring = STRCOMPRESS(STRING(kyrange[0], FORMAT='(F4.2)'), /REMOVE_ALL)
  kystring = "[" + kystring + ", " + STRCOMPRESS(STRING(kyrange[1], FORMAT='(F4.2)'), /REMOVE_ALL)
  FileID = "k!Dy!N=" + kystring + "]"
ENDIF

outstr = {GKVs1D}
outstr.mnemonic = "FieldMap"
outstr.title = "FieldMap"
outstr.CodeName = CodeName
outstr.CodePI = CodePi
outstr.values = PTR_NEW(y)
outstr.vrange = [0,ly]
outstr.indices = PTR_NEW(indices)
outstr.fileID = FileID
outstr.runID  =  RunID
grid1 = {Grid}
grid1.boundary = "periodic (closed)"
grid1.irange = [0,nIterations-1]
grid1.range = [0,lx]
grid1.title='x'
grid1.mnemonic='x'
grid1.values = PTR_NEW(x)
outstr.grid1 = grid1
output = OBJ_NEW("GKVs1D", outstr)
IF(ARG_PRESENT(Aout)) THEN BEGIN
  AxStr = {GKVs1D}
  FOR i=0, 10 DO AxStr.(i) = outStr.(i)
  AxStr.title = "!9d!XA!D!9#!N!X/!9d!Xx"
  AxStr.mnemonic = "Appl_x"
  AxStr.values = PTR_NEW(Ax)
  Vmin = MIN(Ax, MAX=Vmax)
  AxStr.vrange = [Vmin, Vmax]
  AxStr.indices = PTR_NEW(indices)
  Axgrid = {Grid}
  Axgrid.irange = [0,nIterations-1]
  Axgrid.boundary = "open"
  Axgrid.range = FLOAT(Axgrid.irange)
  Axgrid.title = "!12l!X"
  Axgrid.mnemonic = "l"
  Axgrid.units = '2!4p!XqR!D0!N'
  Axgrid.values = PTR_NEW(FINDGEN(nIterations))
  AxStr.Grid1 = Axgrid
  AxObj = OBJ_NEW("GKVs1D", AxStr)
  AyObj = AxObj -> MakeCopy(/NoValues)
  AyObj -> Set, title="!9d!XA!D!9#!N!X/!9d!Xy"
  AyObj -> Set, mnemonic="Appl_y"
  AyObj -> Set, values=PTR_NEW(Ay)
  Vmin = MIN(Ay, MAX=Vmax)
  AyObj -> Set, vrange = [Vmin, Vmax]
  Aout = { x:AxObj, y:AyObj }
ENDIF
RETURN, output
END ; ****** FUNCTION FieldMap ****** ;

FUNCTION MakeFieldMap, params=params, _Extra=Extra
;
ky0=0.0846807
ly = 74.1986
lx = 59.0278
sHat = 0.8
q0 = 1.4
IF(TypeOF(params) EQ 8) THEN BEGIN
  thisParams = params
  IF(N_ELEMENTS(params) GT 1) THEN thisParams=params[0]
  ky0 = thisParams.box.kymin
  ly  = 2.*!PI/ky0
  Lx  = thisParams.box.Lx 
  sHat = thisParams.geometry.sHat
  q0 = thisParams.geometry.q0
ENDIF
;
; Create array of default inital mapping points
; near rational surfaces
;
dxrat = ly/(2.*!PI*sHat)
nRats = FIX(Lx/dxrat + 0.4999999)
nx=nRats
xrat  = dxrat*FINDGEN(nRats)
ny = 10
dy = ly/ny
eps = 0.01/lx
xarr = xrat + eps
yarr = ly/ny*FINDGEN(ny)
;
; Check for arrays of x and y starting values 
; on input line.
;
result = GetKeyWord('xstart', Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN BEGIN
  xarr = result
  nx = N_ELEMENTS(xarr)
ENDIF
result = GetKeyWord('ystart', Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN BEGIN
  yarr = result
  ny = N_ELEMENTS(yarr)
ENDIF
output = OBJARR(nx, ny)
FOR i=0, nx-1 DO BEGIN
  FOR j = 0, ny-1 DO BEGIN
    output[i,j] = FieldMap(x0=xarr[i], y0=yarr[j], params=params, _EXTRA=EXTRA)
  ENDFOR
ENDFOR
output[0,0] -> View, /pretty, psym=1, FieldMap=output, ystyle=1
RETURN, output
END  ; ****** MakeFieldMap ****** ;
