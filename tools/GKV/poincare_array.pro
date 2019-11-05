FUNCTION GKVs1D::Poincare_array
;
; Recieves Surface of section plot and returns array whose elements are
; each surface of section plots with number of points increasing with 
; array element.Last element of array is the original object
;



vals = *self.values
vrange=self.vrange
grid = self.grid1 
len = N_ELEMENTS(vals)

output=  objarr(len)

FOR i = 0, len-1 DO BEGIN
temp = self -> makecopy(/novalues, /noerrorbars)
tempvalues = vals[0:i]
tempValuesPtr = PTR_NEW(tempvalues)
temp.values=tempValuesPtr 
temp.vrange = vrange
tempGrid = gkvsd_gridcopy(grid, irange=[0,i])
tempGrid.range = self.grid1.range
PTR_FREE, temp.grid1.values
temp.grid1 = tempGrid

output[i] = temp 


ENDFOR

RETURN, output
END  ; Poincare_array