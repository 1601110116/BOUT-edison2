FUNCTION dmag_analysis, _Extra=Extra

;
; Designed to read in directories of poincare surface of section plots
; and produce a 1d object d_mag vs time
; Arguments:  
;     start-     initial time point to start at.  Default=1 
;
;     irange-     Max range of time to be sampled
;
;     iskip-      integer for how many time slices skip between plots.  
;                default = 1
;
;     sourcedir- directory structure should be "prefix#", with # being
;                the time slice of surface of section plot.  In such a
;                case, sourcedir should be set = "prefix"
;     lx-  Length in rho_i of radial domain
;
;


CD, CURRENT=current_working_directory


; check _Extra for 'range' argument
irange = 10
result = GETKEYWORD('irange',Extra)
IF(QUERY_INTEGER(result)) THEN irange = result

; check _Extra for 'sourcedir' argument.  String should read 'G#t'
predir = "G8t"
result=GetKeyWord('sourcedir', Extra)
IF(QUERY_STRING(result) ) THEN predir = result

lx=59.52381
result=GetKeyWord('lx', Extra)
IF(QUERY_REAL(result)) THEN lx = result

iskip = 1
result=GetKeyWord('iskip', Extra)
IF(QUERY_INTEGER(result)) THEN iskip = result

start = 1
result=GetKeyWord('start',Extra)
IF(QUERY_INTEGER(result)) THEN start = result 


range= irange/iskip


outputarr = fltarr(range+1)
outputerror = fltarr(range+1)
skiparr = fltarr(range+1)
surf = "/surface.out"

 
 
count = 0
skipdir = 0
FOR i_time = 1,range+1 DO BEGIN

dir = predir + STRCOMPRESS(STRING(start-1+i_time*iskip), /REMOVE_ALL)

PRINT, "reading directory ", dir
surf = "/surface.out"
surf2 = dir + surf
poin = gyro_poincare(surf2)

IF (TYPEOF(poin) EQ 2) THEN BEGIN
  PRINT, dir, " did not contain surface.out file"
  skiparr[skipdir] = i_time
  skipdir = skipdir + 1
  GOTO, nofile
ENDIF

FOR j  = 0, 99 DO BEGIN
  poin[j,0] -> poincare_norm, lx
ENDFOR

poinarr = REFORM(poin[*,0],10,10)
poinmag = dmag(poinarr)


poinmag -> signalwindow, n_cycles=[1500,3000]
poinmag -> restrict

mystats = poinmag -> stats()

outputarr[i_time-1] = mystats.avg
outputerror[i_time-1] = mystats.avgpm
gkvdelete, poin[*,0]
gkvdelete, poin[*,1]
gkvdelete, poinarr
gkvdelete, mystats
gkvdelete, poinmag
count = count + 1 
CD, current_working_directory
nofile:
ENDFOR


gridStructure = {Grid}
gridStructure.Mnemonic = "t"
gridStructure.Title = "t"
gridStructure.units = "a/c!ds!n"
gridValues=  (INDGEN(range+1)+1) ;INSERT VALUES
gridStructure.Values= PTR_NEW(gridValues)
nPoints = N_ELEMENTS(gridValues)
gridStructure.Range= [0,range]
gridStructure.irange=[0,range]

objStructure = {GKVs1D}
objStructure.Grid1 = gridStructure
objStructure.mnemonic = 'Dmag'
objStructure.title = 'Dmag'
objStructure.indices = PTR_NEW('*')
objStructure.units = ''
objStructure.codename = 'Gyro 8.1'
objStructure.CodePI = 'J. Candy'
;objValues = output[i,1,*]
;objValues = REFORM(objValues)
objStructure.values = PTR_NEW(outputarr)
objStructure.vrange = [MIN(outputarr), MAX(outputarr)]
obj = obj_new("gkvs1d", objstructure)



RETURN, obj
END  ;
