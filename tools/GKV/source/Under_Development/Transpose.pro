PRO GKVs2D::Transpose
;
; Purpose:
;
;	This proceedure transposes the values (ane errorBars) array
;	and interchanges the grid structures for the two independent 
;	variables.
;
; Input Keywords:
;
;	None
;
; Output Keywords:
;
;	None
;
;  Written by W.M. Nevins
;	6/6/01
;
grid1 = self.Grid2
grid2 = self.Grid1
newValues = TRANSPOSE(*self.values)
IF(PTR_VALID(self.ErrorBars)) THEN BEGIN
	errorBars = TRANSPOSE(*self.ErrorBars)
	PTR_FREE, self.ErrorBars
	self.ErrorBars = PTR_NEW(errorBars)
ENDIF
self.grid1 = grid1
self.grid2 = grid2
PTR_FREE, self.values
self.values = PTR_NEW(newValues)
RETURN
END



PRO GKVs3D::Transpose, _extra=extra
;
; Just a dummy for now to trap calls to transpose on 
; objects of dimensionality higher than 2.
;
; Written by W.M. Nevins
;	6/6/01
;
MESSAGE, "Transpose not implimented at dimensionality higher than 2, Returning", /INFORMATIONAL
RETURN
END  ; ****** GKVs3D::Transpose ****** ;