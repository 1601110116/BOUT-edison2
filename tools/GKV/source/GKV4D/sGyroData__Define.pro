;
; Define sGyroData class for 4D data processing
;
; Written by E.Wang
;    6-22-09
;
FUNCTION make_sGyroData, _EXTRA=Extra
;
; create new sGyroData 
;
result = OBJ_NEW("sGyroData", _Extra=Extra)
RETURN, result
END;  ****  make_sGyroData  **** ;


FUNCTION sGyroData::SubSet, n, this=that
;
; returns nth subset of the data as conditioned GKVs4D object
;
interval = (*self.intervals)[n,*]
interval = REFORM(interval)
thisSubSet = self.data -> GetData(t=interval)
output = thisSubSet -> Gyro_Condition(l_perp=1, /closetheta)
thisSubSet -> Trash
RETURN, output
END ;  **** sGyrodata::Subset ****


PRO sGyroData::Set, data=GyroData, nSubsets=nSubsets, intervals=intervals

type = TypeOf(GyroData)
if(type EQ 11) THEN BEGIN
    IF OBJ_ISA(GyroData, "GKVsd_Null") THEN BEGIN
        self.data = GyroData
    ENDIF ELSE BEGIN
        MESSAGE, "Data must be a gkvsd_Null ojbect", /INFORMATIONAL
    ENDELSE
ENDIF ELSE BEGIN
    MESSAGE, "Data must be a gkvsd_Null object", /INFORMATIONAL
ENDELSE

type = TypeOf(nSubsets)
IF(Query_Integer(nSubsets)) THEN self.nSubsets = nSubsets

type = TypeOf(intervals)
IF(Type NE 0) THEN self.intervals = intervals

END  ;  **** sGyroData::Set  **** ; 

PRO sGyroData::Get, data = data, nSubsets=nSubsets, intervals=intervals

data = self.data
nSubsets = self.nSubsets
intervals = self.intervals

END ; **** sGyroData::GET ****

FUNCTION sGyroData::nSubsets
RETURN, self.nSubsets
END 


PRO sGyroData::Info
PRINT, "nSubsets = ", self.nSubSets
PRINT, "Intervals = ", self.intervals
self.data -> Info
END ; **** sgyroData::Info   ****

PRO sGyrodata::Trash
OBJ_DESTROY, self
END ;   **** sGyroData::Trash

PRO sGyroData::Cleanup
PTR_FREE, self.intervals
RETURN
END  ;  ***  sGyroData:Cleanup **** ;


FUNCTION sGyroData::INIT, GyroDataObj=GyroDataObj, nSubSets=nSubsets, trange=trange_in


Catch, error
IF error NE 0 THEN BEGIN
    Catch, /Cancel
    ok=Error_Message(/Trace)
    RETURN, 0
ENDIF
self.data = GyroDataObj
self.nSubsets = nSubsets


intervals = FLTARR(nSubSets, 2 )
IF(OBJ_ISA(GyroDataObj, "GKVsdNull")) THEN BEGIN
    GyroDataObj -> get, axis=4, irange=irange, gridValues=tPtr
    t=*tPtr
    tMin = t[irange[0]]
    tMax = t[irange[1]]
    trange=[tMin,tMax]
ENDIF

IF(N_ELEMENTS(trange_in) EQ 2) THEN trange=trange_in
dT = (trange[1] - trange[0])/nSubSets
tintervals = trange[0] + dT*FINDGEN(nSubSets+1)
intervals[0,0] = trange[0]
FOR n=1, nSubSets-1 DO BEGIN
    intervals[n,0] = tIntervals[n]
    intervals[n-1,1] = tIntervals[n]
ENDFOR
intervals[nSubsets-1,1] = tIntervals[nSubsets]

self.Intervals = PTR_NEW(intervals)
RETURN, 1
END ; **** sGyroData:INIT ****

PRO sGyroData__Define
struct = { sGyroData, $
           data : OBJ_NEW('GKVsdNull', {GKVsdNull}), $
           nSubSets:  0L,                            $
           intervals:  PTR_NEW(0)                   }
END ; ****** PRO sGyroData__Define ****** ;        
