FUNCTION GKVs1D::AvgValuePM, PRINT=print, _Extra=Extra
;
; return structure containing avg and
; error bars over specified time interval
;
; Written by W.M. Nevins
;  7/25/05
;
temp = self -> MakeCopy()
result = GetKeyWord('tRange', Extra)
IF(Query_Real(result) + Query_Integer(result)) THEN BEGIN
	temp -> SignalWindow, axis=1, range=result
ENDIF
temp -> restrict
tempStats = temp -> Stats()
temp -> Trash
result = {	Name	:	"AvgValuePM",		$
		avg	:	tempStats.avg,		$
		avgPM	:	tempStats.avgPM		}
IF(KeyWordSet(print)) THEN PRINT, "avg = ", result.avg, " +/- ", result.avgpm
RETURN, result
END
