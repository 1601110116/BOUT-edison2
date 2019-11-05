FUNCTION GKVs1D::MakeEvenOdd, iaxis
;
;
; Purpose:
;
; decomposes "self" into even and odd functions
; in the independent variable specifiec by 'iaxis'
; about the centroid of this coordinate
;
; Results are returned as a structure containing
; even and odd parts of 'self'
;
; Written by W.M. Nevins
;     2/28/2006 
;
; Compute even/odd parts
;
values = *Self.values
evenValues = 0.5*(values + REVERSE(values,iaxis))
 oddValues = 0.5*(values - REVERSE(values,iaxis))
nDims = self -> NumDims()
evenSq = FLOAT(evenValues*CONJ(evenValues))
 oddSq = FLOAT( oddValues*CONJ( oddValues))
CASE nDims OF
	1 :	BEGIN
			EvenIntensity = evenSq
			 OddIntensity =  oddSq
		END
	2 :	BEGIN
			EvenIntensity = TOTAL(evenSq, iaxis)
			 OddIntensity = TOTAL( oddSq, iaxis)
		END
	3 :	BEGIN
	 CASE iaxis OF 
	   1 : BEGIN
      temp = TOTAL(evenSq, 2)
      EvenIntensity = TOTAL(temp, 2)
      temp = TOTAL( oddSq, 2)
       OddIntensity = TOTAL(temp, 2)
         END 
	   2 : BEGIN
      temp = TOTAL(evenSq, 1)
      EvenIntensity = TOTAL(temp, 2)
      temp = TOTAL( oddSq, 1)
       OddIntensity = TOTAL(temp, 2)
         END 
	   3 : BEGIN
			temp = TOTAL(evenSq, 1)
			EvenIntensity = TOTAL(temp, 1)
			temp = TOTAL( oddSq, 1)
			 OddIntensity = TOTAL(temp, 1)
			   END 
			    ENDCASE 
		END
	4: BEGIN
	 CASE iaxis of
	   1: BEGIN
     temp = TOTAL(evenSq, 2)
     temp2 = TOTAL(temp, 2)
     EvenIntensity = TOTAL(temp2,2)
     temp = TOTAL( oddSq, 2)
     temp2 = TOTAL(temp,2)
     OddIntensity = TOTAL(temp2, 2)
     END
	   2: BEGIN
     temp = TOTAL(evenSq, 1)
     temp2 = TOTAL(temp, 2)
     EvenIntensity = TOTAL(temp2,2)
     temp = TOTAL( oddSq, 1)
     temp2 = TOTAL(temp,2)
     OddIntensity = TOTAL(temp2, 2)
     END
	   3: BEGIN
     temp = TOTAL(evenSq, 1)
     temp2 = TOTAL(temp, 1)
     EvenIntensity = TOTAL(temp2,2)
     temp = TOTAL( oddSq, 1)
     temp2 = TOTAL(temp,1)
     OddIntensity = TOTAL(temp2, 2)
     END
	   4: BEGIN
	   temp = TOTAL(evenSq, 1)
	   temp2 = TOTAL(temp, 1)
	   EvenIntensity = TOTAL(temp2,1)
	   temp = TOTAL( oddSq, 1)
	   temp2 = TOTAL(temp,1)
     OddIntensity = TOTAL(temp2, 1)
     END
     ENDCASE
	  END
	  
ENDCASE
;
; create output objects
;
even = self -> MakeCopy(/NoValues)
even.values = PTR_NEW(evenValues)
even.title = self.title + "!Deven!N"
eMin = GKVsd_MIN(evenValues, MAX=eMAX)
even.vrange=[eMin, eMax]

odd = even -> MakeCopy(/NoValues)
odd.title = self.title + "!Dodd!N"
odd.values=PTR_NEW(oddValues)
oMin = GKVsd_MIN(oddValues, MAX=oMAX)
odd.vrange=[oMin, oMax]

IevenStr = {GKVs1D}
FOR i=0,10 DO IevenStr.(i) = self.(i)
IevenStr.mnemonic = "I_even"
IevenStr.title = "!12I!X!Deven!N{" + self.title + "}"
IevenStr.Indices = PTR_NEW("*")
IevenStr.Values = PTR_NEW(EvenIntensity)
Vmin = MIN(EvenIntensity, Max=Vmax)
IevenStr.vrange=[Vmin, Vmax]

CASE nDims OF
	1 :	IevenStr.Grid1 = GKVsd_GridCopy(self.Grid1)
	2 : 	CASE iaxis OF
			1 :	IevenStr.Grid1 = GKVsd_GridCopy(self.Grid1)
			2 :	IevenStr.Grid1 = GKVsd_GridCopy(self.Grid2)
		ENDCASE
	3 :	CASE iaxis OF 
	    1 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid1)
      2 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid2)	    
	    3 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid3)    
	    ENDCASE
	4: CASE iaxis OF 
      1 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid1)
      2 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid2)     
      3 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid3)    
      4 : IevenStr.Grid1 = GKVsd_GridCopy(self.Grid4)  
      ENDCASE
ENDCASE
Ieven = OBJ_NEW("GKVs1D", IevenStr)

Iodd = Ieven -> MakeCopy(/NoValues)
Iodd -> Set, Title="!12I!X!Dodd!N{" + self.title + "}"
Iodd -> Set, Mnemonic="I_odd"
Iodd -> Set, Values=PTR_NEW(OddIntensity)
Vmin = MIN(OddIntensity, Max=Vmax)
Iodd -> Set, vRange=[Vmin, Vmax]

Itotal = Ieven -> Plus(Iodd)
Itotal -> Set, Title="!12I!X!Dtotal!N{" + self.title + "}"
Itotal -> Set, Mnemonic="I_total"


result = {	Name:	"EvenOdd",	$
		Even:	even,		$
		 Odd:	 odd,		$
		Ieven:	Ieven,		$
		Iodd:	Iodd,		$
		Itotal:	Itotal		}		
return, result
END ; ****** GKVs1D::MakeEvenOdd ****** ; 

FUNCTION GKVs1D::EvenOdd, arg, _Extra=extra
;
; Purpose:
;
;	Returns even and odd parts of 'self'
;	with respect to specified axis.
;
;
; Arguments:
;
;			The (optional) argument is any legal axis identifier.
;			That is, either an integer between 1 and nDims, or
;			a STRING containing an axis mnemonic.
;
; Keywords:
;
;	     Axis	If no argument is provided, then this keyword may be 
;			used to identify the axis. Set axis 
;			equal to any legal axis identifier (see above).
;
;	 mnemonic	Set the mnemonic of the selected axis equal to a two element
;			array specifying the desired range in the selected independent 
;			variable.  This both
;			specifies the independent variable in which "self" is to
;			be decomposed into even and odd parts, and sets the centroid
;			about which this decomposition will be performed to the
;			center of the specified range. 
;
;	   irange	Set 'irange' to a two-element (integer) array to reset the signal
;			window of the selected independent variable.  The value of irange
;			is interpreted as an index into the grid.values array.
;
;	    range	Set 'range' to a two-element (floating point) array to set the
;			range in the independent variable.
;
;	   icenter	Set 'icenter' to integer to specify the index of the 
;			desired centroid in the selected independent variable.
;			The value of icenter is interpreted as an index into 
;			the grid.values array.
;
;	    center	Set 'center' to the desired floating point value of the
;			desired centroid in the selected independent variable.
;
;	    iOffset	Set to an integer to correct offset in in the selected 
;			independent variable 
;
; Written by W.M. Nevins
;	2/22/06
;
;

input = self -> MakeCopy()
;
; Find axis identifier
;
input -> GET

CASE N_PARAMS() OF
	0	:	iaxis = input -> AxisIrange(     axisValue=axisValue, _Extra=extra)
	1	:	iaxis = input -> AxisIrange(arg, _Extra=extra)
	else	:	BEGIN
				MESSAGE, 'EvenOdd called with too many arguments', /INFORMATIONAL
				RETURN, 0
			END
ENDCASE
IF(iaxis LT 0) THEN BEGIN
	MESSAGE, 'No valid axis identifier', /INFORMATIONAL
	RETURN, 0
ENDIF
input -> Restrict
;
; Get grid structure
;
axisString = STRING(iaxis, FORMAT='(i1)')
commandString = 'Grid = input.Grid' + axisString
ok = EXECUTE(commandString)
gridValues = *(Grid.values)
IF(Grid.uniform EQ 0) THEN BEGIN
	MESSAGE, "Selected Grid is non-uniform, results may be inconsistent", /INFORMATIONAL
	input -> ScaleAxis, iaxis, /Uniform
	commandString = 'Grid = input.Grid' + axisString
	ok = EXECUTE(commandString)
	gridValues = *(Grid.values)
ENDIF
dx = gridValues[1] - gridValues[0] 
;
; Get irange, imin, imax
;
irange= grid.irange
imin = irange[0]
imax = irange[1]
nPoints = imax - imin + 1
;
; Check for "iOffset" keyword
;
iOffset = 0
result = GetKeyWord("iOffSet", Extra)
IF(Query_Integer(result)) THEN iOffset = result
;
; Check for "center" keyWord
;
iCenter = nPoints/2 + (nPoints MOD 2) - 1
center  = GridValues[iCenter] + dx/2.*( 1-(Npoints MOD 2) )
result = GetKeyWord("Center", Extra)
IF(Query_Real(axisValue) + Query_Integer(axisValue)) THEN result=axisValue
IF(Query_Real(result) + Query_Integer(result)) THEN BEGIN
	center=result
	epsilon = MIN((gridValues-center)^2, iCenter)
	iCenter = iCenter + iOffset
	iHalf = MIN([iMax - iCenter, iCenter - iMin])
	iMin = iCenter - iHalf
	iMax = iCenter + iHalf
	iRange = [iMin, iMax]
	input -> signalwindow, axis=iaxis, irange=irange
	input -> restrict
ENDIF
;
; check for 'iCenter' keyword
;
result = GetKeyWord('iCenter', Extra)
IF(Query_Integer(result)) THEN BEGIN
	iCenter = result + iOffset
	iHalf = MIN([iMax - iCenter, iCenter - iMin])
	iMin = iCenter - iHalf
	iMax = iCenter + iHalf
	iRange = [iMin, iMax]
	input -> signalwindow, axis=iaxis, irange=irange
	input -> restrict
ENDIF
;
; input is now centered about reference point, 
;
result = input -> MakeEvenOdd(iaxis)
input -> Trash
RETURN, result
END ; ****** GKVs1D::EvenOdd ****** ;
