FUNCTION AnalysisProtocol,	Data=GEMData, Path=pathIn, File=FileIn, RunID=runIdIn, TRange=trangeIn,	$
				CodePI=CodePIin, FileID=FileIDin, rhoStar=rhoStarIn, aScale=aScaleIn,	$
				CorrFcns=CorrFcns, yCorrFcns=yCorrFcns, rCorrFcns=rCorrFcns, kSpects=kSpects,	$
				debug=debug, Save=save, NoCorrs=noCorrs
;
;  Purpose:
;
;	This proceedure is designed for the analysis of
;	TEM parameter scans from GEM.  It assumes that
;	there is a user present (and paying attention);
;	analyzes data in the selected GEM output files,
;	and returns a structure containing information
;	about this data.
;
; Written by W.M. Nevins
;       9/24/008
; after a simialr GS2 analysis protocal written
;	5/26/02
;
; Get current directory
;
CD, CURRENT=DirectoryIn
IF(TypeOf(pathIn) EQ 7) THEN CD, pathIn
;
; IF No data supplied in call to AnalysisProtocal, then 
; read data from NetCDF file (to be chosen by user)
;
inputData=1b
IF(TypeOF(GEMData) NE 8) THEN BEGIN	; GEMData is not a structure, query user for NetCDF file
	GEMData = GEM_Data(CodePI='CodePIin', FileID=FileIDin,RunID=runIDin,	$
			   rhoStar=rhoStarIn, aScale=aScaleIn)
	InputData=0b
ENDIF

;
; get phi data
;
phi0 = GemData.phi.xy
;
;
; Decompose phi0 into y-average and deviations from same
;
phistr = phi0 -> deltaSq(axis='y')
;
; save phistr
;
IF KEYWORD_SET(save) THEN gkv_saveStructure, phistr
;
; Copy various fluxes into 'Flux' array
;
;Flux = OBJARR(2,3)
;Flux[0,0] = GS2Data.ES_PART_FLUX_1 -> MakeCopy()
;Flux[1,0] = GS2Data.ES_PART_FLUX_2 -> MakeCopy()
;Flux[0,1] = GS2Data.ES_MOM_FLUX_1  -> MakeCopy()
;Flux[1,1] = GS2Data.ES_MOM_FLUX_2  -> MakeCopy()
;Flux[0,2] = GS2Data.ES_HEAT_FLUX_1 -> MakeCopy()
;Flux[1,2] = GS2Data.ES_HEAT_FLUX_2 -> MakeCopy()
;**** NEED TO FIX READ TO HIST.OUT before proceeding *****
;
; reset title and mnemonic of time grid in various fluxes
;
;FOR i=0,1 DO BEGIN
;	FOR j=0,2 DO BEGIN
;		Flux[i,j] -> set, axis=1, Gridtitle='t', GridMnemonic='t'
;	ENDFOR
;ENDFOR
;
; Get trange
;
CASE N_ELEMENTS(trangeIn) OF
	0:	BEGIN
			phistr.avg -> view
			Flux[0,0] -> view, Flux[0,1:2]
			Flux[1,0] -> view, Flux[1,1:2]
			MESSAGE, "Enter trange", /INFORMATIONAL
			trange=FLTARR(2)
			READ, trange, PROMPT="trange = ", FORMAT='(2G)'
		END
	1:	BEGIN
			IF(TypeOf(trangeIn) EQ 7) THEN BEGIN
				IF(STRUPCASE(trangeIn) EQ 'ALL') THEN FLUX[0,0] -> GET, axis='t', GridRange=trange
			ENDIF
		END
	2:	BEGIN
			IF(Query_Real(trangeIn)) THEN BEGIN
				trange=trangeIn
			ENDIF ELSE BEGIN
				MESSAGE, "Bad trange on command line", /INFORAMTIONAL
				phistr.avg -> view
				Flux[0,0] -> view, Flux[0,1:2]
				Flux[1,0] -> view, Flux[1,1:2]
				MESSASGE, "Enter new trange", /INFORMATIONAL
				trange=FLTARR(2)
				READ, trange, PROMPT="trange = ", FORMAT='(2G)'
			ENDELSE
		END
ENDCASE
;
; Restrict data and interpolate onto a uniform time grid
;
FOR i=0,1 DO BEGIN
	phistr.(i) -> signalwindow, t=trange
	phistr.(i) -> restrict
	FOR j=0,2 DO BEGIN
		Flux[i,j] -> signalwindow, t=trange
		Flux[i,j] -> restrict
		Flux[i,j] -> ScaleAxis, 't', /Uniform
	ENDFOR
ENDFOR
;
; Compute average flux, and estimated error in these averages
;
AvgFlux = FLTARR(2,3,2)
FOR i=0,1 DO BEGIN
	FOR j=0,2 DO BEGIN
		theseStats = Flux[i,j] -> stats()
		AvgFlux[i,j,0] = theseStats.avg
		AvgFlux[i,j,1] = theseStats.avgPM
	ENDFOR
ENDFOR
result = CREATE_STRUCT('Name', 'GS2Analysis', 'tRange', tRange, 'AvgFlux', avgFlux)
;
; Perform correlation analysis on phi0
;
IF KEYWORD_SET(noCorrs) THEN BEGIN
	IF KEYWORD_SET(save) THEN gkv_saveStructure, result
	IF(NOT InputData) THEN gkvDelete, GS2Data
	CD, directoryIn	
	RETURN, result
ENDIF
IF KEYWORD_SET(debug) THEN WINDOW, 0
phiCorrs = phistr.delta -> taucorrs('x', CorrFcns=CorrFcns, yCorrFcns=yCorrFcns, 	$
					 rCorrFcns=rCorrFcns, kSpects=kSpects, debug=debug)
result = CREATE_STRUCT(result, "PhiCorrs", phiCorrs)
;
; Compute radial averages of taucorr, rcorr, ycorr, and yhalf
;
tauCorrAvg = phiCorrs.tauCorr -> Avg('x', /IgnoreErrors)
rCorrAvg   = phiCorrs.rCorr   -> Avg('x', /IgnoreErrors)
yCorrAvg   = phiCorrs.yCorr   -> Avg('x', /IgnoreErrors)
yHalfAvg   = phiCorrs.yHalf   -> Avg('x', /IgnoreErrors)
result = CREATE_STRUCT(result, 'tauCorr', tauCorrAvg, 'rCorr', rCorrAvg, 'yCorr', yCorrAvg, 'yHalf', yHalfAvg)
;
; Now, filter avg of phi, using computed tauCorrAvg and rCorrAvg
;
tauCorr = tauCorrAvg -> GetValues()
  rCorr =   rCorrAvg -> GetValues()
temp1 = phistr.avg -> Filter('t', dT=tauCorr)
AvgPhi=      temp1 -> Filter('x', dL=  rCorr)
temp1 -> Trash
ExBShear = avgPhi -> d2byd('x', units='c!Ds!N/a', title='ExBShear', mnemonic='ExBShear')
result = CREATE_STRUCT(result, 'ExBShear', ExBShear)
shearSpect = ExBShear -> xspect()
shearSpect -> set, units='(c!Ds!N/a)*(!4q!X!Ds!N)'
shearSpect_w = ShearSpect -> INT('omega', units='(c!Ds!N/a)!U2!N*(!4q!X!Ds!N)')
shearspect_wk= shearSpect_w -> INT('k_x', units='(c!Ds!N/a)!U2!N')
result = CREATE_STRUCT(result, 'ShearSpect', shearSpect, 'ShearSpect_w', ShearSpect_w, 'Shear_Spect_wk', ShearSpect_wk -> GetValues())

IF KEYWORD_SET(save) THEN gkv_saveStructure, result
;
; Clean up ...
;
shearspect_wk -> trash
IF(NOT InputData) THEN gkvDelete, GS2Data

CD, directoryIn
RETURN, result
END  ;  ****** AnalysisProtocal ******  ;
