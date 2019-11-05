PRO GKVs1D::Abs, _Extra=EXTRA
;
; This proceedure replaces the self with 
; its absolute value.
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
values = ABS(values)
self.values = PTR_NEW(values)
self.title = '!9!!!X' + self.title + '!9!!!X'
self.mnemonic = 'Abs_' + self.mnemonic
vMin = GKVsd_Min(values, MAX=vMax)
self.vrange = [vMin, vMax]
self -> SET, _Extra=Extra
RETURN
END  ; ****** PRO Abs ****** ;

FUNCTION GKVs1D::Abs, _Extra=Extra
;
; This function returns the absolute
; values of 'self'
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
result -> Abs, _Extra=Extra
RETURN, result
END  ; ****** FUNCTION Abs ****** ;
