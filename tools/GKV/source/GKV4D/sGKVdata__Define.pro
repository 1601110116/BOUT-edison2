;
; Define sGKVdata class for 4D data processing
;
; Written by W.M. Nevins
;	3/25/08
;
FUNCTION makeGKVdata, _Extra=Extra
;
; creates sGKVdata object from GKVs3D object
result = OBJ_NEW("sGKVdata", _Extra=Extra)
RETURN, result
END  ;  **** makeGKVdata  **** ;

FUNCTION sGKVData::SubSet, n
interval = (*self.intervals)[n,*]
thisSubSet = self.data -> MakeCopy()
thisSubSet -> SignalWindow, axis=3, irange=interval
thisSubSet -> Restrict
Return, thisSubSet
END  ;  **** sGKVdata::SubSet  ****  ;

PRO sGKVdata::SET, sGKVdata=sGKVdata, nSubsets=nSubsets, intervals=intervals

type = TypeOf(sGKVdata)
IF(type EQ 11) THEN self.sGKVdata = sGKVdata

type = TypeOf(nSubsets)
IF(Type NE 0) THEN self.nSubsets = nSubsets

type = TypeOf(intervals)
IF(Type NE 0) THEN self.intervals = intervals

END  ;  **** sGKVdata::SET ****  ;

PRO sGKVdata::GET, data=data, nSubsets=nSubsets

data = self.data
nSubsets = self.nSubSets
intervals = self.Intervals

END  ;  **** sGKVdata::GET ****  ;

FUNCTION sGKVdata::nSubSets
;
; Returns number of data subsets
;
RETURN, self.nSubSets
END  ;  ******  sGKVdata:::nSubSets  ******  ;

PRO sGKVdata::Info
;
; Prints information about an sGKVdata object
;
Print, "Data Info = ", self.data
Print, "nSubsets  = ", self.nSubSets
Print, "Intervals = ", self.intervals
END  ;  ****** sGKVdata::Info  ******  ;

PRO sGKVdata::CleanUp
PTR_FREE, self.intervals
Return
END  ;  ****  sGKVdata::CleanUp  ****  ;
 
FUNCTION sGKVdata::INIT, GKVs3Dobj=GKVs3Dobj, nSubSets=nSubsets
;
; Preliminary init for testing object methods ...
;
Catch, error 
IF error NE 0 then Begin 
   Catch, /Cancel             ; Cancel error trap 
   ok=Error_Message(/Trace)   ; print out a trace-back in the output log 
   RETURN, 0                  ; Return 0 for unsuccessful call to INIT function. 
ENDIF 
self.data = GKVs3Dobj
self.nSubsets = nSubsets

intervals = LONARR(nSubSets, 2)
GKVs3Dobj -> get, axis=3, irange=irange
nValues = irange[1] - irange[0] - 1
dN = nValues/nSubSets
FOR n=0, nSubSets-1 DO BEGIN
	intervals[n,0] = irange[0] + n*dN
	intervals[n,1] = intervals[n,0] + dn-1
ENDFOR
self.Intervals = PTR_NEW(intervals)
 
RETURN, 1 
END  ;  ****  sGKVdata::INIT  ****  ;

PRO sGKVdata__Define
struct = {	sGKVdata,			$	; define sGKVdata class
		data	:	OBJ_NEW('GKVs3D', {GKVs3D}),$	; use GKVs3D object for starters
		nSubSets:	0L,		$	; number of data subsets
		intervals:	PTR_NEW(0)	}	; pointer to an array of time intervals	
END  ;  **** sGKVdata__Define ****  ;
