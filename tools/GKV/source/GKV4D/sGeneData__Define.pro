;
; Define sGeneData class for 4D data processing
;
; Written by W.M. Nevins
;	2/10/09
;
FUNCTION make_sGeneData, _Extra=Extra
;
; creates sGeneData object from GeneData object
;
result = OBJ_NEW("sGENEdata", _Extra=Extra)
RETURN, result
END  ;  **** make_sGeneData  **** ;

FUNCTION sGeneData::SubSet, n, _Extra=Extra
;
; Returns nth data subset as a GKVs4D object
;
; Written by W.M. Nevins
;	2/10/09
;
interval = (*self.intervals)[n,*]
interval = REFORM(interval)
thisSubSet = self.data -> GetData(t=interval)
IF(typeOf(Extra) EQ 8) THEN BEGIN
	ThisSubSet -> CloseTheta_k, _Extra=Extra
ENDIF
xydata = thisSubSet -> ktox()
thisSubset -> Trash
xydata -> scaleAxis, 't', /Uniform
;
; Maybe need to re-scale data into different units here?
;
Return, xydata
END  ;  **** sGENEdata::SubSet  ****  ;

PRO sGENEdata::SET, data=GeneData, nSubsets=nSubsets, intervals=intervals

type = TypeOf(GeneData)
IF(type EQ 11) THEN BEGIN
	IF OBJ_ISA(GeneData, "GeneData") THEN BEGIN
		self.data = Genedata
	ENDIF ELSE BEGIN
		MESSAGE, "Data must be GeneData object", /INFORMATIONAL
	ENDELSE
ENDIF ELSE BEGIN
	MESSAGE, "Data must be GeneData object", /INFORMATIONAL
ENDELSE

type = TypeOf(nSubsets)
IF(Query_Integer(nSubsets)) THEN self.nSubsets = nSubsets

type = TypeOf(intervals)
IF(Type NE 0) THEN self.intervals = intervals

END  ;  **** sGeneData::SET ****  ;

PRO sGENEData::GET, data=data, nSubsets=nSubsets, intervals=intervals

data = self.data
nSubsets = self.nSubSets
intervals = self.Intervals

END  ;  **** sGENEdata::GET ****  ;

FUNCTION sGeneData::nSubSets
;
; Returns number of data subsets
;
RETURN, self.nSubSets
END  ;  ******  sGeneData:::nSubSets  ******  ;

PRO sGeneData::Info
;
; Prints information about an sGKVdata object
;
Print, "nSubsets  = ", self.nSubSets
Print, "Intervals = ", self.intervals
self.data -> Info
END  ;  ****** sGeneData::Info  ******  ;

PRO sGeneData::Trash
OBJ_DESTROY, self
END ;  ****** sGeneData::Trash ******  ;

PRO sGeneData::CleanUp
PTR_FREE, self.intervals
Return
END  ;  ****  sGKVdata::CleanUp  ****  ;
 
FUNCTION sGeneData::INIT, GeneDataObj=GeneDataObj, nSubSets=nSubsets, trange=trange_in
;
; Create a new instance of an sGeneData object from 
; information on command line. Generally invoked 
; using:
;
;	FUNCTION make_sGeneData, _Extra=Extra
;
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF 
self.data = GeneDataObj
self.nSubsets = nSubsets
;
; choose intervals to be of equal length in time
; (rather than in number of time-points) as 
; time step can be non-uniform.
;
intervals = FLTARR(nSubSets, 2)
IF(OBJ_ISA(GeneDataObj, "GeneData")) THEN BEGIN
	nullGKVObj = GeneDataObj -> GKV_Obj()
	nullGKVObj -> get, axis=4, irange=irange, gridValues=tPtr
	t=*tPtr
	tMin = t[irange[0]]
	tMax = t[irange[1]]
	trange=[tMin, tMax]
ENDIF

IF(OBJ_ISA(GeneDataObj, "mGeneData")) THEN BEGIN
	GeneDataObj -> Get, Intervals=myIntervals
	tMin = MIN(myIntervals)
	tMax = MAX(myIntervals)
	tRange = [tmin, tmax]
ENDIF

IF(N_ELEMENTS(trange_in) EQ 2) THEN BEGIN
  tRange[0] = tRange[0] > trange_in[0]
  tRange[1] = tRange[1] < tRange_in[1]
ENDIF

dT = (trange[1] - trange[0])/nSubSets
tintervals = trange[0] + dT*FINDGEN(nSubSets+1)
tIntervals[nSubSets] = trange[1]
intervals[0,0] = trange[0]
FOR n=1, nSubSets-1 DO BEGIN
	intervals[n,0] = tIntervals[n]
	intervals[n-1,1] = tIntervals[n]
ENDFOR
intervals[nSubsets-1,1] = tIntervals[nSubSets]

self.Intervals = PTR_NEW(intervals)
 
RETURN, 1 
END  ;  ****  sGeneData::INIT  ****  ;

PRO sGENEdata__Define
struct = {	sGeneData,						$	; define sGKVdata class
		data	:	OBJ_NEW('GeneData', {GeneData}),	$	; use GKVs3D object for starters
		nSubSets:	0L,					$	; number of data subsets
		intervals:	PTR_NEW(0)				}	; pointer to an array of time intervals	
END  ;  **** sGENEdata__Define ****  ;
