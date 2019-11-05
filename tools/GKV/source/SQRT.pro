PRO GKVs1D::SQRT, _Extra=EXTRA
;
; This proceedure replaces the self with 
; its square root.
;
;  Arguments
;
;		None
;
;  Keywords
;
;		Any keywords will be passed to 'Set' to allow
;		changes in the title, mnemonic, or units
;
; Written by W.M. Nevins
;	10/6/03
;
FORWARD_FUNCTION GKVsd_MIN
values = *self.values
PTR_FREE, self.values
values = SQRT(values)
self.values = PTR_NEW(values)
IF(PTR_VALID(self.errorBars)) THEN BEGIN
	errorBars = 0.5*(*self.errorBars)/values
	PTR_FREE, self.errorbars
	self.errorBars = PTR_NEW(errorBars)
ENDIF
self.title = '(' + self.title + ')!U1/2!N'
self.mnemonic = 'SQRT_' + self.mnemonic
self.units = '(' + self.units + ')!U1/2!N'
vMin = GKVsd_Min(values, MAX=vMax)
self.vrange = [vMin, vMax]
self -> SET, _Extra=Extra
RETURN
END  ; ****** PRO SQRT ****** ;


FUNCTION GKVs1D::SQRT, _Extra=Extra
;
; This function returns the 
; square root of 'self'
;
;  Arguments
;
;		None
;
;  Keywords
;
;		Any keywords will be passed to 'Set' to allow
;		changes in the title, mnemonic, or units
;
; Written by W.M. Nevins
;	10/6/03
;
result = self -> MakeCopy()
result -> SQRT, _Extra=Extra
RETURN, result
END  ; ****** FUNCTION SQRT ****** ;
