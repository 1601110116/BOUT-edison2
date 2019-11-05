FUNCTION GKVs1d::Lupanov, arg
;
; Calculate the Lupanov coefficient of two fieldlines.  arg should be adjacent fieldline
;

my_x1 = *(self.grid1.values)
my_y1 = *(self.values)

arg -> get, values=my_y2_ptr, axis=1, gridvalues=my_x2_ptr

my_x2 = *(my_x2_ptr)
my_y2 = *(my_y2_ptr)



my_x1 = Shift_Poincare_Arr(my_x1)
my_x2 = Shift_Poincare_Arr(my_x2)
my_y1 = Shift_Poincare_Arr(my_y1)
my_y2 = Shift_Poincare_Arr(my_y2)


my_dx = my_x1 - my_x2
my_dy = my_y1 - my_y2

my_dx2 = my_dx * my_dx
my_dy2 = my_dy * my_dy

my_diff = my_dx2 + my_dy2

nPoints=  (SIZE(my_diff))[1]
gridStructure = {Grid}
gridStructure.Mnemonic = "turns"
gridStructure.Title = "turns"
gridStructure.units = ""


gridvals = FINDGEN(nPoints)
gridStructure.Values= PTR_NEW(gridvals)
nPoints = N_ELEMENTS(gridValues)
gridStructure.Range=[1,3000]
gridStructure.irange=[0,2999]

objStructure = {GKVs1D}
objStructure.Grid1 = gridStructure
objStructure.mnemonic = 'Lupanov Exponent'
objStructure.title = 'Lupanov Exponent'
objStructure.indices = PTR_NEW('*')
objStructure.units = ''
objStructure.codename = 'Gyro 8.1'
objStructure.CodePI = 'J. Candy'
objStructure.values = PTR_NEW(my_diff)

result = OBJ_NEW("gkvs1d", objStructure)
RETURN, result
END


FUNCTION Shift_Poincare_Arr, arr
;
; Takes 1-D array in [-.5,.5] periodic range format and shifts data when passing a periodic boundary  
; 
;

nPoints = (SIZE(arr))[1]

FOR i = 1,nPoints-1 DO BEGIN
  dx = arr[i-1] - arr[i]
  IF dx GT .9 THEN BEGIN
    FOR j = i,nPoints-1 DO BEGIN
      arr[j] = arr[j] + 1     
    ENDFOR   
  ENDIF    
  IF dx LT -.9 THEN BEGIN
    FOR j = i,nPoints-1 DO BEGIN
      arr[j] = arr[j] - 1
      ENDFOR  
  ENDIF
ENDFOR

RETURN, arr
END
