FUNCTION GKVs2d::GetGrid, arg, _Extra=extra
;
;
;
CASE N_PARAMS() OF
  0 : axis = self -> AxisIrange(     _Extra=extra)
  1 : axis = self -> AxisIrange(arg, _Extra=extra)
  else  : BEGIN
        MESSAGE, 'DeltaSq called with too many arguments', /INFORMATIONAL
        RETURN, 0
      END
ENDCASE
IF(axis LT 0) THEN BEGIN
  MESSAGE, 'No valid axis identifier', /INFORMATIONAL
  RETURN, 0
ENDIF



gridnum = STRING(axis)
gridnum = STRCOMPRESS(gridnum, /REMOVE_ALL)
command_line = 'thisGrid = self.grid' + gridnum 
ok = EXECUTE(command_line)
result = gkvsd_gridcopy(thisGrid)

RETURN, result
END