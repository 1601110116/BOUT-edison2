
FUNCTION GKVs1D::PDF, 	nStd = nn, xmin=xl, xmax = xu, nbins = nb, 	$
			Weight=Weight, Reverse_Indices=R
;
; Purpose:
;
;	Creates a probability distribution function object (GKVpdf) from
;	from the data of 'self'.  The resulting GKVpdf object describes the
;	pdf of the bare values in 'self' (there is no effort to group events,
;	etc.  If this is desired, see GKVs2D::GKVpdf)
;
; Input Keywords:
;
;	nStd	Determine range of output independent variable by specifying the 
;		number of standart deviations (of 'self.values') desired.
;		Defaults to 8. (Optional)
;
;	xmin,	Determine range of output independent varialbe by specifying 
;	xmax	'xmin' and 'xmax'.  Default is to center range at mean of
;		'self.values', with a range of ±4*sigma. (Optional).
;
;	nBins	Specify number of bins to use when calling "HIISTOGRAM" to
;		form the PDF.  Default is to use SQRT(nPoints), 
;		where 'nPoiints' is the number of points is the number of 
;		elements in 'self.values'. (Optional)
;
; Reverse_Indices Set this keyword (to an as yet undefined symbol) to have
;		the Reverse_Indices array produced by the HISTOGRAM function
;		(used to compute the PDF) returned. 
;
;
; Written by W.M. Nevins
;	6/16/01
; Modified by W.M. Nevins
; to include possibilty of weighted averages
;	4/16/2008
;
; Modified by W.M. Nevins
; to return the REVERSE_INDICES
; array.
;	1/29/2009
;
valuePtr = self -> GetValues(/Open)
values = *valuePtr
PTR_FREE, valuePtr
IF(Query_COMPLEX(values)) THEN values=FLOAT(values)
nPoints = N_ELEMENTS(values)
moments = MOMENT(values)
std = SQRT(moments[1])
nBins = FIX(SQRT(nPoints))
IF QUERY_INTEGER(nb) THEN nBins = nb
PRINT, 'number of bins for PDF = ', nBins
nStd = 4.
IF(N_ELEMENTS(nn) EQ 1) THEN nStd = nn/2. > 1.
xmin = moments[0] - nStd*std
IF (N_ELEMENTS(xl) EQ 1) THEN xmin = xl
xmax = moments[0] + nStd*std
IF(N_ELEMENTS(xu) EQ 1) THEN xmax = xu
binSize = (xmax - xmin)/(nBins -1.)
binCenters = xmin + binSize/2 + binSize*FINDGEN(nBINS)
;
; Create structure for creation of output pdf object
;
out = {GKVs1D}
FOR i=0, N_TAGS({GKVsd}) DO out.(i) = self.(i)
out.title = 'PDF(' + self.title + ')'
out.mnemonic = self.mnemonic + '_pdf'
out.units='1/(' + self.units + ')'
indices = *(self.indices)
out.indices = PTR_NEW(indices)
grid1 = {Grid}
grid1.values = PTR_NEW(binCenters)
grid1.range = [xmin, xmax]
grid1.irange = [0,nBins-1]
grid1.units = self.units
grid1.title = self.title
grid1.mnemonic = self.mnemonic
out.grid1 = grid1
;
; Check for "Weight" keyword.
;
IF(TypeOf(Weight) NE 11)THEN BEGIN
	pdf = HISTOGRAM(values, NBINS=nbins, MAX=xmax, min=XMIN, REVERSE_INDICES=R)
	pdf = FLOAT(pdf)
	errorBars = (pdf/SQRT(pdf-1))/(nPoints*binSize)
	pdf = pdf/(nPoints*binSize)
	out.values = PTR_NEW(pdf)
	vmax = MAX(pdf)
	out.vrange = [0., vmax]
	out.ErrorBars = PTR_NEW(errorBars)
	output = OBJ_NEW('GKVs1D', out)
ENDIF ELSE BEGIN
	wPtr = WEIGHT -> GetValues()
	wValues = *wPtr
	wInfo = SIZE(wValues)
	wElements = wInfo[wInfo[0]+2]
	wValues = REFORM(wValues, wElements, /OVERWRITE)
	pdf = HISTOGRAM(values, NBINS=nbins, MAX=xmax, min=XMIN, REVERSE_INDICES=R)
	pdf = FLOAT(pdf)
;	errorBars = (pdf/SQRT(pdf-1))/(nPoints*binSize)
;	pdf = pdf/(nPoints*binSize)
	w = pdf
	FOR i=0,nBins-1 DO BEGIN
		rMin = r[i]
		rMax = r[i+1]-1
		IF(rMax GT rMin) THEN BEGIN
			w[i] = TOTAL( wValues[ R[rMin:rMax] ] )
		ENDIF ELSE BEGIN
			w[i] = 0.
		ENDELSE
	ENDFOR
	w = w/(binSize*wElements)
;	out.values = PTR_NEW(pdf)
;	vmax = MAX(pdf)
;	out.vrange = [0., vmax]
;	out.ErrorBars = PTR_NEW(errorBars)
;	pdfObj = OBJ_NEW('GKVs1D', out)
	wOut = {GKVs1D}
	FOR i=0, N_TAGS({GKVs1D})-1 DO wOut.(i) = out.(i)
	wMax = MAX(w, MIN=wMin)
	wOut.values = PTR_NEW(w)
	wOut.vrange = [wmin, wmax]
	WEIGHT -> GET, title=title, mnemonic=mnemonic, units=units
	wOut.title=title
	wOut.mnemonic=mnemonic
	wOut.units = units
	wErrors = w/SQRT(pdf-1)
	wOut.ErrorBars = PTR_NEW(wErrors)
	wObj = OBJ_NEW('GKVs1D', wOut)
	output = wObj
;	output = {	NAME	:	"pdf",		$
;			pdf	:	pdfObj,		$
;			wPDF	:	wObj		}
ENDELSE
RETURN, output
END ; ****** GKVs1D::GKVpdf ****** ;
