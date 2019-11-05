FUNCTION GKVsd::Execute, arg
;
; Purpose:
;
;		Applies any legal IDL function to argument of 'self'
;		and returns the result as a new GKV object.
;
; Arguments:
;
;		Accepts a string argument containing the name of an IDL function.
;
; Output:
;
;		Returns a GKV object of same dimensionality as 'self' whose
;		values are 'arg'(self.values)
;
; Written by W.M. Nevins
;	8/21/00
;
; Check argument
;
IF(TypeOf(arg) NE 7) THEN BEGIN
	MESSAGE, "Requires STRING argument containing the name of an IDL function", /INFORMATIONAL
	RETURN, 0
ENDIF
result = self -> MakeCopy(/NoErrorBars)
SelfValues = *result.values
commandStr = 'ResultValues = ' + arg + '(SelfValues)'
ok = EXECUTE(commandStr)
IF(NOT ok) THEN BEGIN
	MESSAGE, "Error executing " + arg, /INFORMATIONAL
	result -> trash
	RETURN, 0
ENDIF
PTR_FREE, result.values
result.values = PTR_NEW(ResultValues)
result.vrange = [GKVsd_MIN(ResultValues, MAX=vmax), vmax]
result.title = arg + '(' + self.title + ')'
result.mnemonic = arg + '_' + self.mnemonic
result.units = ''
RETURN, result
END ; ****** Execute ****** ;

FUNCTION GKVs2D::OmegaTrap, TransitTime=tt
;
; Purpose:
;
;		Acts on objects containing the potential in a 
;		slice or average perpendicular to B, and computes
;		the trapping frequency.
;
;		Assumes that 'self' is in pg3eq units.
;
; Written by W.M. Nevins
;	8/21/00
;
nDims = self -> NumDims()
IF(nDims LT 2) THEN BEGIN
	MESSAGE, 'OmegaTrap acts on objects of 2-D (or greater)', /INFORMATIONAL
	RETURN, 0
ENDIF
;
; Compute second derivatives of 'self'
;
dselfdx    = self -> dbyd(1)
d2selfdxdy = dselfdx -> dbyd(2)
dselfdx -> trash
d2selfdx2  = self -> d2byd(1)
d2selfdy2  = self -> d2byd(2)
;
; Form (d2'self'/dxdy)^2, store in temp1
;
temp1 = d2selfdxdy -> times(d2selfdxdy)
temp1 -> set, units='(c!Ds!N/L!DT!N)!U2!N'
d2selfdxdy -> trash
;
; Form (d2'self'/dx^2)*(d2'self'/dy^2) and store in temp2
;
temp2 = d2selfdx2 -> times(d2selfdy2)
temp2 -> set, units='(c!Ds!N/L!DT!N)!U2!N'
d2selfdx2 -> trash
d2selfdy2 -> trash
;
; Form the square of the trapping frequency
;
OmegaTrapSq = temp2 -> Minus(temp1)
temp1 -> trash
temp2 -> trash
;
; Take it's square root
;
result = OmegaTrapSq -> Execute('SQRT')
OmegaTrapSq -> Trash
;
; Set title, units
;
title = '!4X!X!Dtrap!N!K'
mnemonic = 'Omega_trap'
units = 'c!Ds!N/L!DT!N'
;
; Now, normalize by multiplying by the transit time
;
transitTime = 1.		; default is no normalization.
IF(N_ELEMENTS(tt) NE 0) THEN BEGIN
	oldResult = result
	result = oldResult -> Times(tt)
	oldResult -> Trash
	title = title + '!4s!X!Dtransit!N'
	units = ''
ENDIF	
;
; Add appropriate title and units
;
result -> Set, title=title, mnemonic=mnemonic, units=units

RETURN, result

END ; ****** OmegaTrap ****** ;