FUNCTION GKV4D_SampleFunction, thisSet
; sample command file to act on GYRO data in (r, y, t) format. 
;
x_0=85.
phiStr = thisSet -> DeltaSq(axis=2)
ref_0 = phiStr.delta -> slice(x=x_0)
theta_0=-!PI/2.
theta_1= 0.
theta_2= !PI/2.
dPhi_0 = ref_0 -> slice(theta=theta_0)
dPhi_1 = ref_0 -> slice(theta=theta_1)
dPhi_2 = ref_0 -> slice(theta=theta_2)

ref_0 -> Trash

Corrs_0 = phiStr.delta -> GKVs4D::Xcorr(ref=dPhi_0, /Norm)
Corrs_1 = phiStr.delta -> GKVs4D::Xcorr(ref=dPhi_1, /Norm)
Corrs_2 = phiStr.delta -> GKVs4D::Xcorr(ref=dPhi_2, /Norm)

dPhi_0 -> Trash
dPhi_1 -> Trash
dPhi_2 -> Trash
phistr.delta -> Trash


CAT = {	avg	: 	phistr.avg,	$
	Iavg	:	phistr.Iavg,	$
	I	: 	phistr.I,	$
	I_x	:	phistr.I_x	}
	
AVG = {	Corr_0		:	Corrs_0,	$
	Corr_1		:	Corrs_1,	$
	Corr_2		:	Corrs_2} ;	$
;	Corr_3		:	Corrs_3,	$
;	Corr_4		:	Corrs_4	}
;
; and, we're done ...
;
result={	CAT	: CAT,		$
		AVG	: AVG		}
RETURN, Result
END  ;  ****** GKV4D_SampleFunction  ******  ;

FUNCTION GKV4D_Execute, thisSet, userFunction
IF(TypeOf(userFunction) EQ 0) THEN userFunction="GKV4D_SampleFunction"
OK = EXECUTE("Result = " + userFunction + "(thisSet)")
RETURN, Result
END  ;  ****** GKV4D_Execute  ******  ;

FUNCTION GKV4D, sGKVdata=sGKVdata, UserFunction=commandFile, _Extra=Extra
;
; 
; Purpose:
;
; This function processes subsets of 
; a data set as directed by commandFile
;
; Keywords:
;
; sGKVdata	Set this keyword to an sGKVdata object.
;		sGKVdata objects will deal with selecting
;		data subsets and providing GKVs3D objects
;		for processing. REQUIRED
;
; userFunction	This keyword identifies a user-supplied function
;		containing the GKV commands to be
;		excecuted in processing the data in "thisSet".
;		This file must return a strucdture which
;		contains two substructures,
;		CAT (containing objects to be concatenated)
;		and AVG (containiing objects to be averaged).
;		Defaults to the sample function in this .pro file.
;
; any additional keywords will be put into the structure "Extra" and
; passed to SubSet for use in conditioning the data.
;		
;
; Written by W.M. Nevins
; 	3/25/08
;
; for the moment, dataInfo will be the full name of a file containing 
; the data to be processed.
;
;
nSets = sGKVdata -> nSubsets()
firstSet = sGKVdata -> SubSet(0, _Extra=Extra)
result = GKV4D_Execute(firstSet, CommandFile)
firstSet -> Trash
IF(nSets EQ 1) THEN RETURN, result

CAT = result.CAT
AVG = result.AVG
;
; Establish an error trap
; always return SOMETHING!
;
;CATCH, ERROR_STATUS
;IF(ERROR_STATUS NE 0) THEN BEGIN
;  PRINT, "Error index: ", ERROR_STATUS
;  PRINT, "Error message: ", !ERROR_STATE.MSG
;  result.cat = cat
;  result.avg=avg
;  CATCH, /CANCEL
;  RETURN, result
;ENDIF

nCATs = N_TAGS(CAT)
nAVGs = N_TAGS(AVG)
titles = STRARR(nAVGs)
mnemonics = STRARR(nAVGs)
FOR iTag=0, nAVGs-1 DO BEGIN
	AVG.(iTag)[0] -> Get, title=thisTitle
	titles[iTag] = thistitle
	AVG.(iTag)[0] -> Get, mnemonic=thisMnemonic
	mnemonics[iTag]=thisMnemonic
ENDFOR

FOR iSet = 1, nSets-1 DO BEGIN
	PRINT, iset
	thisSet = sGKVdata -> SubSet(iSet)
	thisResult = GKV4D_Execute(thisSet, CommandFile)
	thisSet -> Trash
	thisCAT = thisResult.CAT
	thisAVG = ThisResult.AVG
	FOR iTag = 0, nCATs-1 DO BEGIN
		temp = CAT.(iTag)
		nObjs = N_ELEMENTS(temp)
		FOR iObj=0, nObjs-1 DO $
		  CAT.(iTag)[iObj] = temp[iObj] -> Cat(thisCAT.(iTag)[iObj])
		gkvdelete, temp, thisCat.(iTag)
	ENDFOR
	FOR iTag = 0, nAVGs-1 DO BEGIN
		temp = AVG.(iTag)
		nObjs = N_ELEMENTS(temp)
		FOR iObj=0, nObjs-1 DO $
		  AVG.(itag)[iObj] = temp[iObj] -> PLUS(thisAVG.(iTag)[iObj])
		gkvdelete, temp,thisAvg.(iTag) 
	ENDFOR
	GKVdelete, thisResult	
ENDFOR

FOR iTag = 0,nAVGs-1 DO BEGIN
  nObjs = N_ELEMENTS(AVG.(iTag))
  FOR iObj=0, nObjs-1 DO BEGIN
	  temp = AVG.(iTag)[iObj] -> Over(nSets)
    temp -> Set, title=titles[iTag], units=""
    temp -> Set, mnemonic=mnemonics[iTag]
    AVG.(iTag)[iObj] -> Trash
    AVG.(iTag)[iObj] = temp
  ENDFOR
ENDFOR


result.cat=cat
result.avg=avg
CATCH, /Cancel
Return, result
END ; ** GKV4D **  ;
