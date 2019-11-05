FUNCTION GKV_AFitBasis, x, M
result = FLTARR(2, N_ELEMENTS(x))
result[0,*] = x^2
result[1,*] = x - 0.5*x^2 + (1.-x)*ALOG(1.-x)
;  result[1,*] = x^3
result = REFORM(result)
RETURN, result
END ; ****** GKV_AFitBasis ****** ;

FUNCTION GKV_AFitBasisPrime, x, M
result = FLTARR(2, N_ELEMENTS(x))
result[0,*] =2.*x
result[1,*] = -(x + ALOG(1.-x))
;  result[1,*] = 3.*x^2
result = REFORM(result)
RETURN, result
END ; ****** GKV_AFitBasisPrime ****** ;

FUNCTION GKVs2D::omegaTrapSq, 	DeBug=d, fraction=f, pFraction=p, transitTime=tt, 	$
						Shades=shades, OmegaTrapSq=omegaTrapSq, 			$
						levelFactor=lf, Areas=aa, areaFit=afit, 			$
						OmegaSq=omSq, OmegaFit=omSqFit, 				$
						levelSeparatrix=levelSeparatrix, CenterDiff=centerDiff
; 
; Purpose:
;
;		GIven the electrostatic potential (as 'self') in the plane perpendicular to B,
;		this proceedure computes, and returns (as a GKVs2D object) the square of the  
;		(non-local) ExB trapping frequency for (nearly)all points in this plane.    
;		Optional keywords also allow access to information regarding the profile of 
;		Omega_ExB about the five largest maxima in 'self'.
;	
;		The object which this proceedure acts on is assumed to contain
;		the electrostatic potential vs. two orthogonal coordinates
;		in a plane perpendicular to B.
;
;
;	Arguments:
;
;			None
;
;
;	Output:
;
;			This function returns a GKVs2D object containing an estimate of the 
;			ExB trapping frequency in the same plane as the eletrostatic potential
;			supplied in 'self'.
;
;
;	Input Keywords:
;
;		TransitTime	Normally, the (parallel to B) transit time of a typical ion.
;				If 'TransitTime' is set, then the values in the object returned
;				by this function is normalized to the transit time [that is,
;				the object returned contains (Omega_ExB*transitTime)^2, rather
;				than Omega_ExB^2]. (Optional)
;
;		Fraction	The algorithm for dividing the plane perpendicular to B into regions
;				about each potential maximum requires that more than 'Fraction' of the
;				grid points in this plane be allocate to some potential maximum for
;				convergence.  Defaults to 99.5% (Optional).
;
;		pFraction	The algorithm for computing the (non-local) ExB trapping frequency
;				numerically computes the area about each potential maximum as a function
;				of the potential difference from this maximum.  The ExB trapping frequency
;				is then simply related to the derivative of this area wrt thi potential 
;				difference.  The numerical evaluation of A(Dphi) can be inaccurate near 
;				the edges of the potential well due to errors in allocating grid cells among
;				regions about potential maximia.  To avoid this, we take form the ExB trapping
;				frequency only over the first 'pFraction' of the potential well (and then 
;				interpolate from there to Omega_ExB=0 at phi=0).  Defaults to 80%.
;				(Optional)
;
;		levelFactor	Used in determining the number of 'bins' for computing A(Dphi).  Generally,
;				the annular width of the region about a potential maximum corresponding to
;				one 'bin' in Dphi will be about 'levelFactor' grid cells across.  Increasing 
;				levelFactor reduces the number of bins (reducing the resolution of the resulting
;				estimate of Omega_ExB), while decreasing levelFactor increases the 'shot' noise
;				in our estimate of A(Dphi) which results from the absence of any interpolation
;				of boundary points (grid cells are taken to be entirely inside or outside of 
;				the desired area).  Defaults to 5 (from experience...).  
;				(Optional)
;
;		CenterDiff	Set this keyword (i.e., but '/CenterDiff' on the command line) to use the differential
;				approximation to Omega_ExB to evaluate Omega_ExB^2 at the extrema of Phi.
;				Default is to use only dA/dPhi method, as the differential approximation can
;				produce discontinuous results when applied to data with a large range of length
;				scales.  (Optional)
;
;		DeBug		Set this keyword (i.e., put '/Debug' on the command line) to enable extra
;				output which may prove useful in the event that this routine malfunctions.
;				(Optional)
;
;
;	Output Keywords:
;
;		Shades		Set this keyword to any variable name and, on return, this variable will 
;				contain a byte array of the same size as the 'signalwindow' of 'self' containing
;				an 'image' of the OmegaTrapSq computed.  These 'shades' can be use, for example,
;				to shade a surface plot of phi (or phi squared, etc.).
;				(Optional)
;
;		OmegaTrapSq	Set this keyword to any variable name and, on return, this variable will contain 
;				a floating point array of the same size as the 'signalwindow' of 'self' containing
;				the values of OmegaTrapSq compute by this function.
;
;		Areas		Set this keyword to any variable name and, on return, this variable will contain 
;				an Object Array.  The elements of this array are five GKVs1D objects containing  
;				the Area vs. dPhi for the five largest potential maxima in 'self'.  
;				(Optional)
;
;		AreaFit		Set this keyword to any variable name and, on return, this variable will contain 
;				an Object Array.  The elements of this array are five GKVs1D objects containing
;				a (smooth) analytic fit to A(Dphi).
;				(Optional)
;
;		OmegaSq		Set this keyword to any variable name and, on return, this variable will contain 
;				an Object Array.  The elements of this array are five GKVs1D objects containing
;				estimates of the square of the ExB trapping frequency vs. Dphi based on numerical
;				differencing of the 'raw' A(Dphi) for the five largest potential extrema in 'self'.
;				N.B., numerical differencing of the raw A(Dphi) is the algorithm actually employed
;				in computing the values of OmegaTrapSq returned by this Function.
;				(Optional)
;
;		OmegaFit	Set this keyword to any variable name and, on return, this variable will contain 
;				an Object Array.  The elements of this array are five GKVs1D objects containing
;				estimates of the square of the  ExB trapping frequency based on the analytic
;				derivative of the (smooth) analytic fit to A(Dphi).  This information is provided
;				to allow evaluation of the quality of the estimates of the ExB trapping frequency 
;				based on numerical differentiation.
;				(Optional)
;
;
; Written by W.M. Nevins
;	8/22/00
;

;
; Get potential
;
phi = *( self -> getValues() )
absPhi = ABS(phi)
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
PhiSqMax =  phiSq[ reverseIndices[maxOrder] ]
phiMax   = absPhi[ reverseIndices[maxOrder] ]
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
; Form array of phiSq levels
;
LargestMax = MAX(PhiSqMax)
PhiSqLevels = LargestMax*(1. - 1./100.*FINDGEN(101))
areaAtLevels = LONARR(101)
FOR iphiSq=0,100 DO areaAtLevels[iphiSq] = TOTAL(phiSq GT phiSqLevels(iphiSq)) + TOTAL(PhiSqMax LT phiSqLevels(iphiSq))

fraction = 0.995
IF KEYWORD_SET(f) THEN fraction=f
old_marked=0L
n_marked = LONG(TOTAL(map GT 0))
;
; Loop over phiSq levels
;
FOR iphiSq=1,100 DO BEGIN

	iteration_number = 0L
;
	IF KEYWORD_SET(d) THEN PRINT, "OmegaTrapSq:: finished marking Maxima. N_max = ", n_max, ", iPhiSq = ", iPhiSq

start_of_loop: iteration_number = iteration_number+1L
;
; Find all neighbors of marked cells, 
; whose value of phiSq is greater than PhiSqLevels[iphiSq].
; mark with 'MaxOrder' of local maxima 
; (which is stored in non-zero elements of MAP)
;
	FOR ii=0,7 DO BEGIN
		map = map	+ 	 (map eq 0)							$
					*( phiSq LT SHIFT(phiSq, x_offsets[ii], y_offsets[ii]) )	$
					*SHIFT(map, x_offsets[ii], y_offsets[ii])			$
					*(phiSQ GE phiSqLevels[iPhiSq])
	ENDFOR

	IF(iteration_number gt 100) THEN BEGIN
		PRINT, "OmegaTrapSq:: Too many iterations, iteration_number = ", iteration_number
		RETURN, 0
	ENDIF
	old_marked = n_marked
	n_marked = LONG(TOTAL(map gt 0))

	IF KEYWORD_SET(d) THEN BEGIN
		PRINT, "OmegaTrapSq:: Find_Peak iteration number ", iteration_number
		PRINT, "OmegaTrapSq::        Total cells marked =", n_marked
		PRINT, "OmegaTrapSq:: Total # of cells at this PhiSqLevel =", areaAtLevels[iPhiSq]
	ENDIF
	IF (n_marked GT fraction*areaAtLevels[iPhiSq]) THEN GOTO, finished_loop		; Stop iteration when "fraction" of cells are marked (defaults to 99.5%)
	IF (n_marked GT old_marked                   ) THEN GOTO, start_of_loop		; Stop iteration when algorithm fails to find any more cells to mark
;
; (essentially) all elements of the Map array are now filled with the order number of the nearest (in some sense) maximum.
;
finished_loop: 
;
; Now Deal with 'collisions' between maxima
;
	IF KEYWORD_SET(levelSeparatrix) THEN BEGIN
		iteration_number = 0L
Start_Loop_2 : iteration_number = iteration_number+1L
		;
		; Find any neighbors in current iPhi band who have already been marked.
		; Change mark to the index of largest maxima (i.e., the LOWER index).
		; This insures proper location of separatrix between extrema.
		;
		oMap = map
		FOR ii=0,7 DO BEGIN
			map = map	+	 (SHIFT(map, x_offsets[ii], y_offsets[ii]) - map)			$
						*(SHIFT(map, x_offsets[ii], y_offsets[ii]) GT 0)			$
						*(SHIFT(map, x_offsets[ii], y_offsets[ii]) LT map)			$
						*(phiSQ GE phiSqLevels[iPhiSq])							$
						*(phiSQ LE phiSqLevels[iPhiSq-1])						
		ENDFOR
		n_changes = TOTAL(oMap - map)
		IF(iteration_number gt 100)	THEN GOTO, finished_loop_2		; Stop after 100 iterations, or
		IF (n_changes GT 0) 			THEN GOTO, Start_loop_2		; Stop iteration when "fraction" of cells are marked (defaults to 99.5%)
finished_loop_2	:
	ENDIF
	;
	; End loop over iPhi levels
	;
ENDFOR
;
; Now, for each maxima we compute the area (about this maxima) enclosed by the curves of constant Phi
;
pFraction = 0.8					; pFraction is the fraction (in phi) of the peak in PhiSq over which
IF KEYWORD_SET(p) THEN pFraction = p	; we expect to find closed phi contours... (this may be the Achilles heel of this algorithm...)
transitTime = 1.
IF KEYWORD_SET(tt) THEN transitTime = tt
phiMin = FLTARR(n_max)
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
	areaStr.Grid1.Mnemonic	= 'Dphi'
	areaStr.Grid1.Title		= '!4Du!X'
	areaStr.Grid1.units		= self.units 
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
	omSqStr.Grid1.Mnemonic	= 'Dphi'
	omSqStr.Grid1.Title		= '!4Du!X'
	omSqStr.Grid1.units		= self.units
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
	IF(maxArea EQ 0) THEN GOTO, DoneIt
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
	phiMin[max_number-1] = MIN(absPhi[RIMax])
	nlevels = LONG(SQRT(maxArea)/levelFactor)	;	estimate of the number of meaningful contour levels...
	CASE nlevels OF
		0:	BEGIN
				omegaSQ = 0.	; This value will be overwritten by differential approximation.
				centeredPhiLevels = PhiMax[max_number-1]
			END
		1:	BEGIN
				omegaSQ = (2.*!PI*(PhiMax[max_number-1] - PhiMin[max_number-1])/maxArea)^2	; This value will be overwritten by differential approximation.
				centeredPhiLevels = PhiMax[max_number-1]
			END
		ELSE:BEGIN
				dPhi = pFraction*(PhiMax[max_number-1] - phiMin[max_number-1])/(nlevels-1)
				phiLevels = PhiMax[max_number-1] - dPhi*FINDGEN(nlevels)
				area = FLTARR(nlevels)
				FOR ilevel=1, nlevels-1 DO 	$
					Area[ilevel] = TOTAL( (absPhi GE phiLevels[ilevel])*(map EQ max_number) )*dx*dy
				centeredPhiLevels = phiLevels + dPhi/2.
				centeredPhiLevels[0] = PhiMax[max_number-1]
				temp = 2.*!PI*dPhi/(area - SHIFT(area, 1))
				omegaSq = temp^2
				omegaSq[0] = omegaSq[1]	; Just extrapolate in for now (differential approximation doesn't seem to be working...)
			END
	ENDCASE
	centeredPhiLevels = [centeredPhiLevels, phiMin[max_number-1]]
	omegaSq = [omegaSq, 0.]
	
	IF KEYWORD_SET(centerDiff) THEN BEGIN
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
	ENDIF
;
; Now, normalize omegaSq to transit time
;
	omegaSq = omegaSq*transitTime^2
	tempOmegaSq = INTERPOL(OmegaSq, centeredPhiLevels, absPhi[[RIMax]])
	omegaTrapSq[RIMax] = tempOmegaSQ
	IF(max_number GT 5) THEN GOTO, doneIt
;
; for five largest maximia, return GKVs1D objects w/ areas and trapping frequencies vs. dPhiSq
;
	IF ARG_PRESENT(aa) THEN BEGIN
		thisAreaStr = areaStr
		thisAreaStr.values = PTR_NEW(area)
		gridValues = PhiMax[max_number-1] - phiLevels
		thisAreaStr.Grid1.values = PTR_NEW(gridValues)
		thisAreaStr.Grid1.irange = [0, nlevels-1]
		thisAreaStr.Grid1.range  = [0., gridValues[nlevels-1]]
		aa[max_number-1] = OBJ_NEW('GKVs1D', thisAreaStr)
	ENDIF
	IF ARG_PRESENT(afit) THEN BEGIN
		thisAFitStr = aFitStr
		;
		; form 'xValues' -- which is dPhi normalized such that the well-depth is one
		; (that is, xvalues=0 at the potential peak, and xvalues=1 at the minimum potential between peaks)
		;
		xValues = (PhiMax[max_number-1] - phiLevels)/(PhiMax[max_number-1] - PhiMin[max_number-1])
		;
		; Subtract off linear piece, using the ***known*** (from differential approximation above)
		; bounce frequency at top of well as boundary condition
		;
		areaPrime_0 = 2.*!PI*(PhiMax[max_number-1] - PhiMin[max_number-1])/SQRT((phixx*phiyy - phixy^2))
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
		;
		xx  = 0.01*FINDGEN(101)
		;
		; Get basis functions evaluated on xx grid, and compute fitted y-values on xx grid
		;
		aBasis = GKV_AFitBasis(xx)
		yy  = xx*areaPrime_0 + coef[0]*aBasis[0,*] + coef[1]*aBasis[1,*]
		;
		; Form xxValues -- the unnormalized well depth
		;
		xxValues = xx*(PhiMax[max_number-1] - PhiMin[max_number-1])
		;
		; Now, load results into thisAFitStr ... 
		;
		thisAFitStr.values = PTR_NEW(yy)
		thisAFitStr.Grid1.values = PTR_NEW(xxValues)
		thisAFitStr.Grid1.irange = [0, 100]
		thisAFitStr.Grid1.range  = [0., xxValues[100]]
		afit[max_number-1] = OBJ_NEW('GKVs1D', thisAFitStr)
	ENDIF
	IF ARG_PRESENT(omSq) THEN BEGIN
		thisOmSqStr = omSqStr
		thisOmSqStr.values = PTR_NEW(omegaSq)
		gridValues = PhiMax[max_number-1] - centeredPhiLevels
		thisOmSqStr.Grid1.values = PTR_NEW(gridValues)
		thisOmSqStr.Grid1.irange = [0, nlevels]
		thisOmSqStr.Grid1.range  = [0., PhiMax[max_number-1]]
		omSq[max_number-1] = OBJ_NEW('GKVs1D', thisOmSqStr)
	ENDIF
	IF (ARG_PRESENT(omSqFit) AND ARG_PRESENT(aFit)) THEN BEGIN
		thisOmSqFitStr = omSqFitStr
		;
		; Get basis function for aPrime, and form derivatiive of area wrt phiSq (aPrime)
		;
		aBasisPrime = GKV_AFitBasisPrime(xx)
		aPrime = ( areaPrime_0 + coef[0]*aBasisPrime[0,*] + coef[1]*aBasisPrime[1,*] )/ (PhiMax[max_number-1] - PhiMin[max_number-1])
		;
		; Form square of trapping frequency based on fit to area as fcn of dPhi
		;
		omSqFitValues = (2.*!PI/aPrime)^2
		;
		; Normalize to transit time
		;
		omSqFitValues = omSqFitValues*transitTime^2
		thisOmSqFitStr.values = PTR_NEW(omSqFitValues)
		thisOmSqFitStr.Grid1.values = PTR_NEW(xxValues)
		thisOmSqFitStr.Grid1.irange = [0, 100]
		thisOmSqFitStr.Grid1.range  = [0., xxValues[100]]
		omSqFit[max_number-1] = OBJ_NEW('GKVs1D', thisOmSqFitStr)
	ENDIF
doneIt:
ENDFOR
;
	omegaSqMax = MAX(omegaTrapSq, /NaN)
IF( ARG_PRESENT(shades) OR KEYWORD_SET(d) ) THEN BEGIN
	bottom = !COLOR_SETUP_NCOLORS
	tableSize = !D.TABLE_SIZE
	dOmegaSQ = omegaSqMax/(tableSize-bottom-1)
	shades = LONG(bottom + omegaTrapSq/dOmegaSQ)
	IF KEYWORD_SET(d) THEN Shade_Surf, phiSq, x, y, SHADES=shades
ENDIF
;
; Now, lets make a GKVs2D object out of OmegaTrapSq
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
FOR i=1,2 DO result -> GridRestrict, i
result -> set, mnemonic='OmegaSq_ExB', title='!4X!X!S!U2!R!DE!9X!XB!N'
units = '(c!Ds!N/L!DT!N)!U2!N'
IF KEYWORD_SET(tt) THEN units = ''
result -> set, units=units, vrange=[0,omegaSqMax]
result.values = PTR_NEW(omegaTrapSq)

RETURN, result
END ; ****** GKVs2D::omegaTrapSq ****** ;

FUNCTION GKVs3D::omegaTrapSq, 	DeBug=d, fraction=f, pFraction=p, transitTime=tt, 	$
						Shades=shades, OmegaTrapSq=omegaTrapSq, 			$
						levelFactor=lf, levelSeparatrix=levelSeparatrix
;
; Purpose:
;
;		GIven the electrostatic potential (as 'self') in the plane perpendicular to B,
;		as a function of time, this proceedure computes, and returns (as a GKVs2D object)   
;		the square of the (non-local) ExB trapping frequency for (nearly) all points in     
;		this plane.  Optional keywords also allow access to information regarding the profile of 
;		Omega_ExB about the five largest maxima in 'self'.
;	
;		The object which this proceedure acts on is assumed to contain
;		the electrostatic potential vs. two orthogonal coordinates and time
;		in a plane perpendicular to B.
;
;
;	Arguments:
;
;			None
;
;
;	Output:
;
;			This function returns a GKVs3D object containing an estimate of the 
;			ExB trapping frequency vs. time in the same plane as the eletrostatic
;			potential supplied in 'self'.
;
;
;	Input Keywords:
;
;		TransitTime	Normally, the (parallel to B) transit time of a typical ion.
;					If 'TransitTime' is set, then the values in the object returned
;					by this function is normalized to the transit time [that is,
;					the object returned contains (Omega_ExB*transitTime)^2, rather
;					than Omega_ExB^2]. (Optional)
;
;		Fraction		The algorithm for dividing the plane perpendicular to B into regions
;					about each potential maximum requires that more than 'Fraction' of the
;					grid points in this plane be allocate to some potential maximum for
;					convergence.  Defaults to 99.5% (Optional).
;
;		pFraction		The algorithm for computing the (non-local) ExB trapping frequency
;					numerically computes the area about each potential maximum as a function
;					of the potential difference from this maximum.  The ExB trapping frequency
;					is then simply related to the derivative of this area wrt thi potential 
;					difference.  The numerical evaluation of A(Dphi) can be inaccurate near 
;					the edges of the potential well due to errors in allocating grid cells among
;					regions about potential maximia.  To avoid this, we take form the ExB trapping
;					frequency only over the first 'pFraction' of the potential well (and then 
;					interpolate from there to Omega_ExB=0 at phi=0).  Defaults to 80%.
;					(Optional)
;
;		levelFactor	Used in determining the number of 'bins' for computing A(Dphi).  Generally,
;					the annular width of the region about a potential maximum corresponding to
;					one 'bin' in Dphi will be about 'levelFactor' grid cells across.  Increasing 
;					levelFactor reduces the number of bins (reducing the resolution of the resulting
;					estimate of Omega_ExB), while decreasing levelFactor increases the 'shot' noise
;					in our estimate of A(Dphi) which results from the absence of any interpolation
;					of boundary points (grid cells are taken to be entirely inside or outside of 
;					the desired area).  Defaults to 5 (from experience...).  
;					(Optional)
;
;		levelSeparatrix	Set this keyword (il.e., put '/levelSepartrix on the command line) to perform
;					extra loop which insures that the separatrix between neighboring extrema is 
;					(nearly) a level surface of phi.  Default is not to insure that the separatricies
;					are level surfaces of phi.
;					(Optional)
;
;
;		DeBug			Set this keyword (i.e., put '/Debug' on the command line) to enable extra
;					output which may prove useful in the event that this routine malfunctions.
;					(Optional)
;
;
;	Output Keywords:
;
;		Shades		Set this keyword to any variable name and, on return, this variable will 
;					contain a byte array of the same size as the 'signalwindow' of 'self' containing
;					an 'image' of the OmegaTrapSq computed.  These 'shades' can be use, for example,
;					to shade a surface plot of phi (or phi squared, etc.).
;					(Optional)
;
;		OmegaTrapSq	Set this keyword to any variable name and, on return, this variable will contain 
;					a floating point array of the same size as the 'signalwindow' of 'self' containing
;					the values of OmegaTrapSq compute by this function.
;
;
; Written by W.M. Nevins
;	8/25/00
;
pFraction = 0.8					; pFraction is the fraction (in phi) of the peak in PhiSq over which
IF KEYWORD_SET(p) THEN pFraction = p	; we expect to find closed phi contours... (this may be the Achilles heel of this algorithm...)
transitTime = 1.
IF KEYWORD_SET(tt) THEN transitTime = tt
levelFactor = 5.
IF KEYWORD_SET(lf) THEN levelFactor = lf

;
; Get potential
;
phiValues = *( self -> getValues() )
info = SIZE(phiValues)
nx = info[1]
ny = info[2]
nt = info[3]
n_total = nx*ny*nt
omegaTrapSq = FLTARR(nx,ny, nt)
maxOmegaSqValues = FLTARR(nt)
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
krange = self.Grid3.irange
t  = (*self.Grid3.values)[krange[0]:krange[1]]
dt = t[1] - t[0]

;
; Make linear arrays of displacements
;
x_offsets=[-1,-1,-1, 0, 0, 1, 1, 1]
y_offsets=[-1, 0, 1,-1, 1,-1, 0, 1]
;
; Begin loop over time-level
;
FOR it = 0, nt-1 DO BEGIN
	;
	; First get phi values at this time level
	;
	phi = phiValues(*,*,it)
	absPhi = ABS(phi)
	phiSq = phi^2
	localOmegaTrapSQ = REPLICATE(0., nx, ny)
	;
	; Make 'map' array for this time level
	;
	map = REPLICATE(0L, nx, ny)
	;
	; Find all maximia using only array operations for greater speed
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
	; of the ith maximum of phiSq (at this time-level).  
	;  
	maxOrder = REVERSE(SORT(phiSq[reverseIndices]))
	;
	; ... and the ith element of maxOrder is the (one dimensional) index into reverseIndices
	; of the (i+1)th largest maximum of phiSq.  It follows that the magnitude of the (i+1)th 
	; largest maximum of phisq is phiSq[reverseIndices[maxOrder[i]]].  That is, ...
	;
	PhiSqMax =  phiSq[ reverseIndices[maxOrder] ]
	phiMax   = absPhi[ reverseIndices[maxOrder] ]
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
	; Form array of phiSq levels
	;
	LargestMax = MAX(PhiSqMax)
	PhiSqLevels = LargestMax*(1. - 1./100.*FINDGEN(101))
	areaAtLevels = LONARR(101)
	FOR iphiSq=0,100 DO areaAtLevels[iphiSq] = TOTAL(phiSq GT phiSqLevels(iphiSq)) + TOTAL(PhiSqMax LT phiSqLevels(iphiSq))

	fraction = 0.995
	IF KEYWORD_SET(f) THEN fraction=f
	old_marked=0L
	n_marked = LONG(TOTAL(map GT 0))
	;
	; Loop over phiSq levels
	;
	FOR iphiSq=1,100 DO BEGIN

		iteration_number = 0L
	;

start_of_loop: iteration_number = iteration_number+1L
		;
		; Find all neighbors of marked cells at lower values of phiSq, 
		; whose value of phiSq is still greater than PhiSqLevels[iphiSq].
		; mark with 'MaxOrder' of local maxima 
		; (which is stored in non-zero elements of MAP)
		;
		FOR ii=0,7 DO BEGIN
			map = map	+ 	 (map eq 0)										$
						*( phiSq LT SHIFT(phiSq, x_offsets[ii], y_offsets[ii]) )		$
						*(phiSQ GT phiSqLevels[iPhiSq])							$
						*SHIFT(map, x_offsets[ii], y_offsets[ii])					
		ENDFOR

		old_marked = n_marked
		n_marked = LONG(TOTAL(map gt 0))
		IF(iteration_number gt 100) 				THEN GOTO, finished_loop		; Stop after 100 iterations, or
		IF (n_marked GT fraction*areaAtLevels[iPhiSq]) 	THEN GOTO, finished_loop		; Stop iteration when "fraction" of cells are marked (defaults to 99.5%), or
		IF (n_marked GT old_marked) 				THEN GOTO, start_of_loop		; Stop iteration when algorithm fails to find any more cells to mark
		;
		; (essentially) all elements of the Map array with phiSq > phiSqlevels(iPhi)
		; are now filled with the order number of the nearest maximum.
		;
finished_loop: 
	;
	; Now Deal with 'collisions' between maxima
	;
		IF KEYWORD_SET(levelSeparatrix) THEN BEGIN
			iteration_number = 0L
Start_Loop_2 : 	iteration_number = iteration_number+1L
			;
			; Find any neighbors in current iPhi band who have already been marked.
			; Change mark to the index of largest maxima (i.e., the LOWER index).
			; This insures proper location of separatrix between extrema.
			;
			oMap = map
			FOR ii=0,7 DO BEGIN
				map = map	+	 (SHIFT(map, x_offsets[ii], y_offsets[ii]) - map)			$
							*(SHIFT(map, x_offsets[ii], y_offsets[ii]) GT 0)			$
							*(SHIFT(map, x_offsets[ii], y_offsets[ii]) LT map)			$
							*(phiSQ GT phiSqLevels[iPhiSq])							$
							*(phiSQ LT phiSqLevels[iPhiSq-1])						
			ENDFOR
			n_changes = TOTAL(oMap - map)
			IF(iteration_number gt 100)	THEN GOTO, finished_loop_2		; Stop after 100 iterations, or
			IF (n_changes GT 0) 			THEN GOTO, Start_loop_2		; Stop iteration when "fraction" of cells are marked (defaults to 99.5%)
finished_loop_2	:
		ENDIF
		;
		; End loop over iPhi levels
		;
	ENDFOR
	;
	; Now, for each maxima we compute the area (about this maxima) enclosed by the curves of constant Phi
	;
	phiMin = FLTARR(n_max)
	FOR max_number=1L, n_max DO BEGIN
		maxArea = TOTAL( map EQ max_Number )		; area associate with this maximum
		IF(maxArea EQ 0) THEN GOTO, DoneIt
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
		phiMin[max_number-1] = MIN(absPhi[RIMax])
		nlevels = LONG(SQRT(maxArea)/levelFactor)	;	estimate of the number of meaningful contour levels...
		CASE nlevels OF
			0:	BEGIN
					omegaSQ = 0.	; This value will be overwritten by differential approximation.
					centeredPhiLevels = PhiMax[max_number-1]
				END
			1:	BEGIN
					omegaSQ = 0.	; This value will be overwritten by differential approximation.
					centeredPhiLevels = PhiMax[max_number-1]
				END
			ELSE:BEGIN
					dPhi = pFraction*(PhiMax[max_number-1] - phiMin[max_number-1])/(nlevels-1)
					phiLevels = PhiMax[max_number-1] - dPhi*FINDGEN(nlevels)
					area = FLTARR(nlevels)
					FOR ilevel=1, nlevels-1 DO 	$
						Area[ilevel] = TOTAL( (absPhi GE phiLevels[ilevel])*(map EQ max_number) )*dx*dy
					centeredPhiLevels = phiLevels + dPhi/2.
					centeredPhiLevels[0] = PhiMax[max_number-1]
					temp = 2.*!PI*dPhi/(area - SHIFT(area, 1))
					omegaSq = temp^2
					omegaSq[0] = omegaSq[1]
				END
		ENDCASE

		centeredPhiLevels = [centeredPhiLevels, phiMin[max_number-1]]
		omegaSq = [omegaSq, 0.]
		;
		; Now, use the differential expression for omegaTrap to compute the trapping frequency squared right at the maximum.
		;	First compute indices of neighboring points assuming periodic boundary conditions
		;
	;	ix      = ixIndex[max_number-1]
	;	ixPlus  = (ix+1) MOD nx
	;	ixMinus = ix - 1
	;	IF(ixMinus LT 0) THEN ixMinus = ixMinus + nx
	;	iy      = iyIndex[max_number-1]
	;	iyPlus  = (iy+1) MOD ny
	;	iyMinus = iy - 1
	;	IF(iyMinus LT 0) THEN iyMinus = iyMinus + ny
	;
	;	phixx =	(    phi[ixPlus, iy] - 2.*phi[ix, iy] + phi[ixMinus, iy])/(dx^2)
	;	phiyy =	(    phi[ix, iyPlus]	- 2.*phi[ix, iy] + phi[ix, iyMinus])/(dy^2)
	;	phixy =	(   (phi[ixPlus, iyPlus ] - phi[ixMinus, iyPlus ])/(2.*dx)			$
	;			  - (phi[ixPlus, iyMinus] - phi[ixMinus, iyMinus])/(2.*dx)   )/(2.*dy)
	;	omegaSq[0] = ABS(phixx*phiyy - phixy^2)
		;
		; Now, normalize omegaSq to transit time
		;
		omegaSq = omegaSq*transitTime^2
		tempOmegaSq = INTERPOL(OmegaSq, centeredPhiLevels, absPhi[[RIMax]])
		localOmegaTrapSQ[RIMax] = tempOmegaSQ
		;
		; End of loop over maxima at this time level
		;
		doneIt:
	ENDFOR
	maxOmegaSqValues[it] = MAX(localomegaTrapSq, /NaN)	
	omegaTrapSq[*,*,it] = localOmegaTrapSQ	
	;
	; End of loop over time levels
	;
ENDFOR
;
; Now, lets make a GKVs3D object out of OmegaTrapSq
;
result = self -> MakeCopy(/NoValues, /NoErrorBars)
FOR i=1,3 DO result -> GridRestrict, i
result -> set, mnemonic='OmegaSq_ExB', title='!4X!X!S!U2!R!DE!9X!XB!N'
units = '(c!Ds!N/L!DT!N)!U2!N'
IF KEYWORD_SET(tt) THEN units = ''
result -> set, units=units, vrange=[0,MAX(maxOmegaSqValues)]
result.values = PTR_NEW(omegaTrapSq)


RETURN, result
END ; ****** GKVs3D::omegaTrapSq ****** ;


;
; *****************************************************************************************************************
; ******************************************     Copyright Notice     *********************************************
; *                                                                                                               *
; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *
; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *
; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *
; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *
; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *
; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * 
; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *
; *                                                                                                               *
; *****************************************************************************************************************
;
; *****************************************************************************************************************
; **********************************************     DISCLAIMER     ***********************************************
; *                                                                                                               *
; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *
; *  Neither the United States Government nor the University of California nor any of their employees, makes      *
; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *
; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *
; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *
; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *
; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *
; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *
; *  States Government or the University of California, and shall not be used for advertising or product          *
; *  endorsement purposes.                                                                                        *
; *                                                                                                               *
; *****************************************************************************************************************
;