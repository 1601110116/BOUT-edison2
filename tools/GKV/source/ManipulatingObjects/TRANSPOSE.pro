PRO GKVs2D::TRANSPOSE
;
; This proceedure inverts the order of the
; two independent variables
;
; Written by W.M. Nevins
;  3/29/03
;
grid1=self.grid1
grid2=self.grid2
values=TRANSPOSE(*self.values)
PTR_FREE, self.values
self.values = PTR_NEW(values)
self.grid1=grid2
self.grid2=grid1
RETURN
END  ;  ****** Procedure GKVs2D::TRANSPOSE ****** ;

FUNCTION GKVs2D::TRANSPOSE
;
; This returns a copy of 'self'
; in which the order of the
; two independent variables
; is inverted
;
; Written by W.M. Nevins
;  3/29/03
;
result = self -> MakeCopy()
result -> Transpose
RETURN, result
END  ;  ****** Function GKVs2D::TRANSPOSE ****** ;

PRO GKVs3D::TRANSPOSE
;
;
; This returns a copy of 'self'
; in which the order of the
; two independent variables
; is inverted
;
; Written by W.M. Nevins
;	3/31/05
;
grid1=self.grid1
grid2=self.grid2
tValues = *self.grid3.values
nt = N_ELEMENTS(tValues)
info = SIZE(*self.values)
n1=info[1]
n2=info[2]
info[2]=n1
info[1]=n2
type=info[info[0]+1]
IF(type EQ 5) THEN info[info[0]+1] = 4
IF(type EQ 9) THEN info[info[0]+1] = 6
temp = MAKE_ARRAY(SIZE=info)
FOR i=0L, nT-1 DO temp[*,*,i] = TRANSPOSE( (*self.values)[*,*,i] )
PTR_FREE, self.values
self.values=PTR_NEW(temp)
self.grid1=grid2
self.grid2=grid1
RETURN
END  ;  ****** Procedure GKVs3D::TRANSPOSE ****** ;


FUNCTION GKVs3D::TRANSPOSE
;
; This returns a copy of 'self'
; in which the order of the
; two independent variables
; is inverted
;
; Written by W.M. Nevins
;  3/31/05
;
result = self -> MakeCopy()
result -> Transpose
RETURN, result
END  ;  ****** Function GKVs2D::TRANSPOSE ****** ;
