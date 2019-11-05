FUNCTION GKV_AFitBasis, x, M
result = FLTARR(2, N_ELEMENTS(x))
result[0,*] = x^2
result[1,*] = x - 0.5*x^2 + (1.-x)*ALOG(1.-x)
result = REFORM(result)
RETURN, result
END ; ****** GKV_AFitBasis ****** ;

FUNCTION GKV_AFitBasisPrime, x, M
result = FLTARR(2, N_ELEMENTS(x))
result[0,*] =2.*x
result[1,*] = -(x + ALOG(1.-x))
result = REFORM(result)
RETURN, result
END ; ****** GKV_AFitBasisPrime ****** ;

FUNCTION GKVs2D::omegaTrapSq, 	DeBug=d, fraction=f, pFraction=p, transitTime=tt, 	$
						Shades=shades, OmegaTrapSq=omegaTrapSq, 			$
						levelFactor=lf, Areas=aa, areaFit=afit, 			$
						OmegaSq=omSq, OmegaFit=omSqFit
;
; Purpose:
;
;		This proceedure computes the (non-local) ExB trapping frequency
;		distribution in the regions about local maxima in the potential
;
;	Object:
;
;			The object which this proceedure acts on is assumed to contain
;			the electrostatic potential vs. two orthogonal coordinates
;			perpendicular to B.
;
;
; Written by W.M. Nevins
;	8/22/00
;

;
; Get potential
;
phi = *( self -> getValues() )
phiSq = phi^2
info = SIZE(phiSq)
nx = info[1]
ny = info[2]
n_total = nx*ny
omegaTrapSq = FLTARR(nx,ny)
;
; Get Grid spacing (assume a uniform grid within current signal windows),
; and 
;
irange = self.Grid1.irange
x  = (*self.Grid1.values)[irange[0]:irange[1]]
dx = x[1] - x[0]
jrange = self.Grid2.irange
y  = (*self.Grid2.values)[jrange[0]:jrange[1]]
dy = y[1] - y[0]
;
; Find all maxima of phisq
;
map = LONARR(nx, ny)
;
; Make linear arrays of displacements
;
x_offsets=[-1,-1,-1, 0, 0, 1, 1, 1]
y_offsets=[-1, 0, 1,-1, 1,-1, 0, 1]
;
; find all maximia using only array operations for greater speed
;
FOR ii=0,7 DO map = map + ( phisq GT SHIFT(phisq ,x_offsets(ii), y_offsets(ii)) )
map = LONG(map gt 7) 
;
; MAP elements corresponding to maxima of phisq are now marked with 1. 
; All other elements of MAP are zero.
;
; We now use HISTOGRAM to obtain the REVERSE_INDICES array
;
n_maxima = TOTAL(map)
Result   = HISTOGRAM(map, binsize=1, min=1, max=1, REVERSE_INDICES=R)
n_max    = LONG(R(1) - R(0))
IF( n_max ne n_maxima) THEN BEGIN
	MESSAGE, "REVERSE_INDICES array does not have expected structure", /INFORMATIONAL
	RETURN, 0
ENDIF

reverseIndices = R[2:n_max+1]
;
; The ith element of 'reverseIndices' is the (one dimensional) index within MAP (and phiSQ)
; of the ith maximum of phiSq.  
;  
maxOrder = REVERSE(SORT(phiSq[reverseIndices]))
;
; ... and the ith element of maxOrder is the (one dimensional) index into reverseIndices
; of the (i+1)th largest maximum of phiSq.  It follows that the magnitude of the (i+1)th 
; largest maximum of phisq is phiSq[reverseIndices[maxOrder[i]]].  That is, ...
;
PhiSqMax = phiSq[ reverseIndices[maxOrder] ]
;
; While the 'x' and 'y' indices of the (i+1)th maximum are given by
;
iyIndex = reverseIndices[maxOrder]/nx
ixIndex = reverseIndices[maxOrder] MOD nx
;
; Now, number each maximum in MAP by it's order (1 for largest, 2 for next largest, etc.)
;
FOR max_number=0L, n_max-1 DO BEGIN
	index = reverseIndices[maxOrder[max_number]]
	map[index] = max_number + 1
ENDFOR
;
; MAP elements corresponding to maxima of ARRAY are now marked with 
; a unique Max_Number -- an integer between 1 and n_max.
; All other elements are still zero
;
;
; Prepare for loop in which maxima numbers are propagated 'down' in phiSq until
; 
;
fraction = 0.995
IF KEYWORD_SET(f) THEN fraction=f
old_marked=0L
n_marked = LONG(TOTAL(map GT 0))

iteration_number = 0L
;
IF KEYWORD_SET(d) THEN PRINT, "OmegaTrapSq:: finished marking Maxima. N_max = ", n_max

start_of_loop: iteration_number = iteration_number+1L
;
; Find all neighbors of marked cells,
; mark with 'MaxOrder' of local maxima 
; (which is stored in non-zero elements of MAP)
;
FOR ii=0,7 DO BEGIN
	map = map + (map eq 0)*( phiSq LT SHIFT(phiSq, x_offsets[ii], y_offsets[ii]) )*SHIFT(map, x_offsets[ii], y_offsets[ii])
ENDFOR

IF(iteration_number gt 1000) THEN BEGIN
	PRINT, "OmegaTrapSq:: Too many iterations, iteration_number = ", iteration_number
	RETURN, 0
ENDIF
old_marked = n_marked
n_marked = LONG(TOTAL(map gt 0))

IF KEYWORD_SET(d) THEN BEGIN
	PRINT, "OmegaTrapSq:: Find_Peak iteration number ", iteration_number
	PRINT, "OmegaTrapSq::        Total cells marked =", n_marked
	PRINT, "OmegaTrapSq:: Total # of cells in array =", n_total
ENDIF
IF (n_marked GT fraction*n_total) THEN GOTO, finished_loop		; Stop iteration when "fraction" of cells are marked (defaults to 99.5%)
IF (n_marked GT old_marked      ) THEN GOTO, start_of_loop		; Stop iteration when algorithm fails to find any more cells to mark
;
; (essentially) all elements of the Map array are now filled with the order number of the nearest (in some sense) maximum.
;
finished_loop: 
;
; Now, for each maxima we compute the area (about this maxima) enclosed by the curves of constant phiSq
;
pFraction = 0.8					; pFraction is the fraction (in phiSq) of the peak in PhiSq over which
IF KEYWORD_SET(p) THEN pFraction = p	; we expect to find closed phiSq contours... (this may be the Achilles heel of this algorithm...)
transitTime = 1.
IF KEYWORD_SET(tt) THEN transitTime = tt
phiSqMin = FLTARR(n_max)
levelFactor = 5.
IF KEYWORD_SET(lf) THEN levelFactor = lf
IF ARG_PRESENT(aa) THEN BEGIN
	aa = OBJARR(5)
	areaStr 			= {GKVs1D}
	areaStr.mnemonic 	= 'Area'
	areaStr.Title		= 'Area'
	indices = ['*']
	areaStr.Indices 	= PTR_NEW(indices)
	areaStr.Units 		= '(' + self.Grid1.units + ')*(' + self.Grid2.units + ')'
	areaStr.CodeName 	= self.CodeName
	areaStr.CodePI		= self.CodePI
	areaStr.RunID		= self.RunID
	areaStr.FileID		= self.FileID
	areaStr.Grid1.Mnemonic	= 'DphiSq'
	areaStr.Grid1.Title		= '!4Du!X!U2!N'
	areaStr.Grid1.units		= '(' + self.units + ')!U2!N'
	areaStr.Grid1.boundary	= 'open'
	areaStr.Grid1.uniform		= 0B
ENDIF
IF ARG_PRESENT(afit) THEN BEGIN
	afit = OBJARR(5)
	aFitStr 			= areaStr
	aFitStr.mnemonic 	= 'Area_fit'
	aFITStr.Title		= 'Area Fit'
	indices = ['*']
	aFitStr.Indices 	= PTR_NEW(indices)
ENDIF
IF ARG_PRESENT(omSq) THEN BEGIN
	omSq = OBJARR(5)
	omSqStr 			= {GKVs1D}
	omSqStr.mnemonic 	= 'OmegaSQ_ExB'
	omSqStr.Title		= '!4X!X!S!U2!R!DE!9X!XB!N'
	indices = ['*']
	omSqStr.Indices 	= PTR_NEW(indices)
	units = '(c!Ds!N/L!DT!N)!U2!N'
	IF KEYWORD_SET(tt) THEN units = ''
	omSqStr.Units 		= units
	omSqStr.CodeName 	= self.CodeName
	omSqStr.CodePI		= self.CodePI
	omSqStr.RunID		= self.RunID
	omSqStr.FileID		= self.FileID
	omSqStr.Grid1.Mnemonic	= 'DphiSq'
	omSqStr.Grid1.Title		= '!4Du!X!U2!N'
	omSqStr.Grid1.units		= '(' + self.units + ')!U2!N'
	omSqStr.Grid1.boundary	= 'open'
	omSqStr.Grid1.uniform		= 0B
ENDIF
IF ARG_PRESENT(omSqFit) THEN BEGIN
	omSqfit = OBJARR(5)
	omSqFitStr 		= omSqStr
	omSqFitStr.mnemonic 	= 'OmegaSQ_ExB_fit'
	omSqFITStr.Title	= '!4X!X!S!U2!R!DE!9X!XB!N Fit'
	indices = ['*']
	omSqFitStr.Indices 	= PTR_NEW(indices)
ENDIF

FOR max_number=1L, n_max DO BEGIN
	maxArea = TOTAL( map EQ max_Number )		; area associate with this maximum
;
; Use HISTOGRAM to obtain a reverse indices array pointing to the elements in phi, phiSq, map, etc.
; associated with this maxima
;
	Result   = HISTOGRAM(map, binsize=1, min=Max_Number, max=Max_Number, REVERSE_INDICES=Rmax)
	max_Area    = LONG(Rmax(1) - Rmax(0))
	IF( max_Area ne maxArea) THEN BEGIN
		MESSAGE, "Rmax array does not have expected structure", /INFORMATIONAL
		RETURN, 0
	ENDIF
	RIMax = Rmax[2:max_area+1]
	phiSqMin[max_number-1] = MIN(phiSq[RIMax])
	nlevels = LONG(SQRT(maxArea)/levelFactor) > 3		;	estimate of the number of meaningful contour levels...
	dPhiSq = pFraction*(PhiSqMax[max_number-1] - phiSqMin[max_number-1])/(nlevels-1)
	phiSqLevels = PhiSqMax[max_number-1] - dPhiSQ*FINDGEN(nlevels)
	area = FLTARR(nlevels)
	FOR ilevel=1, nlevels-1 DO Area[ilevel] = TOTAL( (phiSq GE phiSqLevels[ilevel])*(map EQ max_number) )*dx*dy
	centeredPhiSqLevels = phiSqLevels + dPhiSq/2.
	centeredPhiSqLevels[0] = PhiSqMax[max_number-1]
	temp = !PI*dPhiSq/(area - SHIFT(area, 1))
	omegaSq = temp^2/centeredPhiSqLevels
	centeredPhiSqLevels = [centeredPhiSqLevels, 0.]
	omegaSq = [omegaSq, 0.]
;
; Now, use the differential expression for omegaTrap to compute the trapping frequency squared right at the maximum.
;	First compute indices of neighboring points assuming periodic boundary conditions
;
	ix      = ixIndex[max_number-1]
	ixPlus  = (ix+1) MOD nx
	ixMinus = ix - 1
	IF(ixMinus LT 0) THEN ixMinus = ixMinus + nx
	iy      = iyIndex[max_number-1]
	iyPlus  = (iy+1) MOD ny
	iyMinus = iy - 1
	IF(iyMinus LT 0) THEN iyMinus = iyMinus + ny
	
	phixx =	(    phi[ixPlus, iy] - 2.*phi[ix, iy] + phi[ixMinus, iy])/(dx^2)
	phiyy =	(    phi[ix, iyPlus]	- 2.*phi[ix, iy] + phi[ix, iyMinus])/(dy^2)
	phixy =	(   (phi[ixPlus, iyPlus ] - phi[ixMinus, iyPlus ])/(2.*dx)			$
			  - (phi[ixPlus, iyMinus] - phi[ixMinus, iyMinus])/(2.*dx)   )/(2.*dy)
	omegaSq[0] = ABS(phixx*phiyy - phixy^2)
;
; Now, normalize omegaSq to transit time
;
	omegaSq = omegaSq*transitTime^2
	tempOmegaSq = INTERPOL(OmegaSq, centeredPhiSqLevels, phiSq[[RIMax]])
	omegaTrapSq[RIMax] = tempOmegaSQ
	IF(max_number GT 5) THEN GOTO, doneIt
;
; for five largest maximia, return GKVs1D objects w/ areas and trapping frequencies vs. dPhiSq
;
	IF ARG_PRESENT(aa) THEN BEGIN
		thisAreaStr = areaStr
		thisAreaStr.values = PTR_NEW(area)
		gridValues = PhiSqMax[max_number-1] - phiSqLevels
		thisAreaStr.Grid1.values = PTR_NEW(gridValues)
		thisAreaStr.Grid1.irange = [0, nlevels-1]
		thisAreaStr.Grid1.range  = [0., gridValues[nlevels-1]]
		aa[max_number-1] = OBJ_NEW('GKVs1D', thisAreaStr)
	ENDIF
	IF ARG_PRESENT(afit) THEN BEGIN
		thisAFitStr = aFitStr
		;
		; form 'xValues' -- which is dPhiSQ normalized so that well-depth is one
		;
		xValues = (PhiSqMax[max_number-1] - phiSqLevels)/(PhiSqMax[max_number-1] - PhiSqMin[max_number-1])
		;
		; Subtract off linear piece, using known bounce frequency at top of well as boundary condition
		;
		areaPrime_0 = !PI*(PhiSqMax[max_number-1] - PhiSqMin[max_number-1])/SQRT((phixx*phiyy - phixy^2)*PhiSqMax[max_number-1])
		yvalues = area - areaPrime_0*xValues
		yvalues[0] = 0.
		;
		; Initial guess of coeffiients--assume that only logarithmic part contributes, and match value at max area
		;
		aCoef = [0., yvalues[nlevels-1]/(pfraction - 0.5*pfraction^2 + (1-pfraction)*ALOG(1.-pfraction))]
		coef = SVDFIT(xValues, yValues, A=aCoef, FUNCTION_NAME="GKV_AFitBasis")
		IF KEYWORD_SET(d) THEN print, coef
		;
		; Form (relatively fine) grid to display fit to area as a function of dPhiSq
		dx = 0.01
		xx  = dx*FINDGEN(101)
		;
		; Get basis functions evaluated on xx grid, and compute fitted y-values on xx grid
		;
		aBasis = GKV_AFitBasis(xx)
		yy  = xx*areaPrime_0 + coef[0]*aBasis[0,*] + coef[1]*aBasis[1,*]
		;
		; Form xxValues -- the unnormalized well depth
		;
		xxValues = xx*(PhiSqMax[max_number-1] - PhiSqMin[max_number-1])
		;
		; Now, load results into thisAFitStr ... 
		;
		thisAFitStr.values = PTR_NEW(yy)
		thisAFitStr.Grid1.values = PTR_NEW(xxValues)
		thisAFitStr.Grid1.irange = [0, 100]
		thisAFitStr.Grid1.range  = [0., xValues[100]]
		afit[max_number-1] = OBJ_NEW('GKVs1D', thisAFitStr)
	ENDIF
	IF ARG_PRESENT(omSq) THEN BEGIN
		thisOmSqStr = omSqStr
		thisOmSqStr.values = PTR_NEW(omegaSq)
		gridValues = PhiSqMax[max_number-1] - centeredPhiSqLevels
		thisOmSqStr.Grid1.values = PTR_NEW(gridValues)
		thisOmSqStr.Grid1.irange = [0, nlevels]
		thisOmSqStr.Grid1.range  = [0., PhiSqMax[max_number-1]]
		omSq[max_number-1] = OBJ_NEW('GKVs1D', thisOmSqStr)
	ENDIF
	IF (ARG_PRESENT(omSqFit) AND ARG_PRESENT(aFit)) THEN BEGIN
		thisOmSqFitStr = omSqFitStr
		;
		; Get basis function for aPrime, and form derivatiive of area wrt phiSq (aPrime)
		;
		aBasisPrime = GKV_AFitBasisPrime(xx)
		aPrime = ( areaPrime_0 + coef[0]*aBasisPrime[0,*] + coef[1]*aBasisPrime[1,*] )/ (PhiSqMax[max_number-1] - PhiSqMin[max_number-1])
		;
		; Form phiSqValues (as opposed to well depth), and compute square of 
		; trapping frequency based on fit to area as fcn of dPhiSq
		;
		phiSqValues = PhiSqMax[max_number-1] - xxValues
		omSqFitValues = (!PI/aPrime)^2/phiSqValues
		;
		; Normalize to transit time
		;
		omSqFitValues = omSqFitValues*transitTime^2
		thisOmSqFitStr.values = PTR_NEW(omSqFitValues)
		thisOmSqFitStr.Grid1.values = PTR_NEW(xxValues)
		thisOmSqFitStr.Grid1.irange = [0, 100]
		thisOmSqFitStr.Grid1.range  = [0., PhiSqMax[max_number-1]]
		omSqFit[max_number-1] = OBJ_NEW('GKVs1D', thisOmSqFitStr)
	ENDIF
doneIt:
ENDFOR
;
omegaSqMax = MAX(omegaTrapSq, /NaN)
IF KEYWORD_SET(d) THEN BEGIN
	bottom = !COLOR_SETUP_NCOLORS
	tableSize = !D.TABLE_SIZE
	dOmegaSQ = omegaSqMax/(tableSize-bottom-1)
	shades = BYTE(bottom + BYTE(omegaTrapSq/dOmegaSQ))
	Shade_Surf, phiSq, x, y, SHADES=shades
ENDIF
;
; Now, lets make a GKVs2D object out of OmegaTrapSq
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
result -> set, mnemonic='OmegaSq_ExB', title='!4X!X!S!U2!R!DE!9X!XB!N'
units = '(c!Ds!N/L!DT!N)!U2!N'
IF KEYWORD_SET(tt) THEN units = ''
result -> set, units=units, vrange=[0,omegaSqMax]
result.values = PTR_NEW(omegaTrapSq)

RETURN, result
END ; ****** GKVs2D::omegaTrapSq ****** ;