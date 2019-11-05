FUNCTION GYRO_Phi_n, 	GYROstr, Path=path, rhoStar=rhoStar, 		$
			rRange=rRangeIn, dn=dn_in
;
; 

CD, CURRENT=currentWorkingDirectory
IF(TypeOF(path) EQ 7) THEN CD, path

nArgs = N_PARAMS()
IF(nArgs EQ 0) THEN GYROStr=NetCDF_Data()

rRange=[0.42,0.58]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar

dn=3
IF(N_ELEMENTS(dn_in) EQ 1) THEN dn=dn_in

thisData=GYROstr.potential -> MakeCopy()
thisData -> signalwindow, axis=1, range=rRange
thisData -> Restrict
phi_n = thisData -> times( rhoStar/SQRT(dn) )
thisData -> trash
temp_r = phi_n -> Avg('r')
temp_r_Sq = temp_r -> AbsSq()
temp_r -> trash
tempSq = phi_n -> AbsSq()
Phi_n -> Trash
tempSq_r = tempSq -> Avg('r')
tempSq -> Trash
Phi_n_r_Sq = temp_r_sq
phi_n_Sq_r = tempSq_r

Phi_n_r_Sq -> set, units='(!4q!X!Ds!N/a)!U2!N(T/e)!U2!N'
Phi_n_Sq_r -> set, units='(!4q!X!Ds!N/a)!U2!N(T/e)!U2!N'
Phi_n_r_Sq -> set, axis=1, gridTitle='N', gridMnemonic='N', gridUnits=''
Phi_n_Sq_r -> set, axis=1, gridTitle='N', gridMnemonic='N', gridUnits=''

IF(nArgs EQ 0) THEN gkvdelete, GYROstr

result = {Name:'GYRO_Phi_n', Phi_n_r_Sq:Phi_n_r_Sq, phi_n_Sq_r:phi_n_Sq_r}
RETURN, result
END
	
