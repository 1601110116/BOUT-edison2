FUNCTION GYRO_Phi, gyroData, path=path, q=qin, aspectRatio=aspectRatioIn, 	$
			rhoStar=rhoStar, r=rin, rRange=rRangeIn, tRange=tRangeIn
;
; 

CD, CURRENT=currentWorkingDirectory
IF(TypeOF(path) EQ 7) THEN CD, path

rRange=[0.4,0.6]
IF(N_ELEMENTS(rRangeIn) EQ 2) THEN rRange=rRangeIn

tRange=[150.,900.]
IF(N_ELEMENTS(tRangeIN) EQ 2) THEN tRange=tRangeIn

phi_k = gyroData -> all_k('n')
phi_k -> signalwindow, axis=1, range=rRange
phi_k -> signalWindow, axis=3, range=tRange
phi_k -> restrict
phi = phi_k -> FFT('n', /INVERSE)
phi_k -> trash
;phistr = phi -> delta(axis='zeta')
;phi -> trash
;phi=phistr.delta

IF(rhoStar LT 1.) THEN rhoStar=1/rhoStar

r=0.5
IF(Query_Real(rin)) THEN r=rin

q=1.4
IF(Query_Real(qin)) THEN q=qin

aspectRatio=2.7775
IF(N_ELEMENTS(aspectRatioIn) EQ 1) THEN aspectRatio=aspectRatioIn

dPerpdRho = (aspectRatio+r)*(rhoStar)/SQRT(1. + ((aspectRatio/r)*q)^2)
print, "Scale Factor = ", dPerpdRho

phi -> get, axis=1, range=rrange
phi -> ScaleAxis, 1, const=rhoStar, offset=-rhoStar*rRange[0], units='!4q!X!Ds!N'
phi -> ScaleAxis, 2, const=dPerpdRho, mnemonic='d_perp', units='!4q!X!Ds!N'
;phi -> window, 'r'
result = phi -> times(rhoStar)
phi -> trash
result -> set, title="!4u!X", mnemonic="Phi", unit='(!4q!X!Ds!N/a)*(T/e)'
RETURN, result
END
	
