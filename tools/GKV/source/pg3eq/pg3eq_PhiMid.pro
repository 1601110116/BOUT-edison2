FUNCTION pg3eq_PhiMid, pg3eqStR, LsubT=L_t
;
;

LsubT = 0.4025
IF(N_ELEMENTS(L_T) EQ 1) THEN LsubT=L_T

nArgs = N_PARAMS()
IF(NArgs EQ 0) THEN BEGIN
	pg3eqStr = pg3eq_data()
ENDIF

;
temp = pg3eqStr.pot_mid -> EXECUTE('float')
GKVdelete, pg3eqStr

phi = temp -> times(LsubT)
temp -> trash
phi -> ScaleAxis, 't', const=LsubT, units='a/c!Ds!N'
phi -> Set, title='!4u!X', mnemonic='Phi', units='(!4q!X!Ds!N/a)(T/e)'

RETURN, phi
END
 
