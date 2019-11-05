FUNCTION Tube_Norm_Field, 	fieldStr, rhoStar_T, ScaleTime, 	$
				Lxin, Lyin, Lzin, rhoIn
;
; This function accepts the field structure produced by 
; Tube_data as an argument, and re-norms both the dependent
; variable and the time-base.
;
; Argument:
;
; 1st argument:	A "field" data structure as produced
;		by Tube_Data (e.g., "phi", "rho", "dne",
;		"dni", "apa", or "upa").
;
; 2nd argument:	The value of rhoStar to be used in rescaling
;		the dependent variable.
;
; 3rd argument:	The constant to be used in rescaling 
;		the time base.
;
; 4th argument:	If non-zero, the the dimension of 
;		the flux tube in the x-direction
;		in units of rho. 
;
; 5th argument:	If non-zero, the the dimension of 
;		the flux tube in the y-direction
;		in units of rho 
;
; 6th argument:	If non-zero, the the dimension of 
;		the flux tube in the z-direction
;		in units of L_T. NOT CURRENTLY IN USE.
;
; 7th argument:	If specified, the symbol (rho_e or rho_s) 
;		to be used in specifying units
;
;
; Written by W.M. Nevins
;	3/14/06
; Modified to include mass ratio for
; examining electron plasmas
;	12/19/06
;
nArgs = N_PARAMS()
IF(nArgs GT 7) THEN BEGIN
	MESSAGE, "Called with wrong number of Arguments", /INFORMATIONAL
	RETURN, 0
ENDIF
IF(nArgs LT 2) THEN BEGIN
	MESSAGE, "Called with wrong number of Arguments", /INFORMATIONAL
	RETURN, 0
ENDIF

IF(TypeOf(fieldStr) NE 8) THEN BEGIN
	MESSAGE, '1st argument is not a structure', /INFORMATIONAL
	return, 0
ENDIF

rhoStar = RhoStar_T
IF(rhoStar GT 1) THEN rhoStar = 1./rhoStar

IF(nArgs EQ 2) THEN ScaleTime=rhoStar

Lx=0
IF(nArgs GT 3) THEN Lx = Lxin

Ly=0
IF(nArgs GT 4) THEN Ly = Lyin

Lz=0
IF(nArgs GT 5) THEN Lz = Lzin

rhoSym = "!4q!X"
IF(nArgs GT 6) THEN rhoSym = rhoIn

nTags = N_TAGS(fieldStr)
tagNames = TAG_NAMES(fieldStr)
IF(nTags LE 1) THEN RETURN, 0
FOR i=1, nTags-1 DO BEGIN
	fieldStr.(i) -> scaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"
	fieldStr.(i) -> Get, values=valuePtr, units=units
	newValues = (*valuePtr)/rhoStar
	PTR_FREE, valuePtr
	newUnits = "(" + rhoSym + "/L!DT!N)(" + units + ")"
	fieldStr.(i) -> Set, values=PTR_NEW(newValues), units=newUnits
	IF( STRCMP(TagNames[i], "r", /FOLD_CASE) ) THEN BEGIN
		fieldStr.(i) -> Get, axis=1, irange=irange
		IF(Lx EQ 0) THEN Lx = FLOAT(irange[1])
		dx = Lx/irange[1]
		fieldStr.(i) -> ScaleAxis, 1, const=dx, units=rhoSym
	ENDIF
	IF(STRCMP(TagNames[i], "xz", /FOLD_CASE) ) THEN BEGIN
		fieldStr.(i) -> Get, axis=1, irange=irange
		IF(Lx EQ 0) THEN Lx = FLOAT(irange[1])
		dx = Lx/irange[1]
	
		fieldStr.(i) -> ScaleAxis, 1, const=dx, units=rhoSym
		fieldStr.(i) -> Get, axis=2, irange=irange
;		IF(Lz EQ 0) THEN Lz = FLOAT(irange[1])
;		dz = Lz/irange[1]
;		fieldStr.(i) -> ScaleAxis, 2, const=dz, units="L!DT!N"
	ENDIF
	IF( STRCMP(TagNames[i], "xy", /FOLD_CASE) ) THEN BEGIN
		fieldStr.(i) -> Get, axis=1, irange=irange
		IF(Lx EQ 0) THEN Lx = FLOAT(irange[1])
		dx = Lx/irange[1]
		fieldStr.(i) -> ScaleAxis, 1, const=dx, units=rhoSym
		fieldStr.(i) -> Get, axis=2, irange=irange
		IF(Ly EQ 0) THEN Ly = FLOAT(irange[1])
		dy = Ly/irange[1]
		fieldStr.(i) -> ScaleAxis, 2, const=dy, units=rhoSym
	ENDIF
ENDFOR

RETURN, 1
END ; ****** Tube_Norm_Field ****** ;

FUNCTION Tube_Norm, TubeStr, _Extra=Extra
;
; This function accepts the data structure produced
; by TubeData as input, and returns the objects from
; the substructure .hist with both tha amplitude and 
; time base re-scaled to gyrokinetic units
;
; Argument:
;
;		The data structure returned by Tube_Data.
;		If none is provided, then Tube_Norm
;		will invoke Tube_Data to allow the user
;		to select the relevant data.
;
; Keywords:
;
;	rhoStar_T	The value of rho/L_T to be used
;			in rescaling both timebase and 
;			data.
;
;	MassRatio	Ratio of the mass of kinetic species
;			to the proton mass.  Defaults to 1.
;			(that is, a proton plasma). Optional
;
;	Lx		Flux tube dimension in x-direction
;			as specified in the GEM.IN file.
;			Defaults such that grid-spacing in 
;			the x-direction is rho_e.
;	Ly		Flux tube dimension in y-direction
;			as specified in the GEM.IN file.
;			Defaults such that grid-spacing in 
;			the y-direction is rho_e.
;	Lz		Flux tube dimension in z-direction
;			as specified in the GEM.IN file.
;			Defaults such that grid-spacing in 
;			the z-direction is one unit.
;
; Written By W.M. Nevins
;	1/5/06
; Modified to include mass ratio for
; examining electron plasmas
;	12/19/06
;
;first get rhoStar_T
;
result = GetKeyWord("rhoStar_T", Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN BEGIN
	rhoStar_T=FLOAT(result)
	rhoStar_T = rhoStar_T < 1./RhoStar_T
ENDIF ELSE BEGIN
	MESSAGE, "RhoStar_T not provided, returning", /INFORMATIONAL
	RETURN, 0
ENDELSE
rhoStarSq = rhoStar_T*rhoStar_T
;
; Next get mass ratio
;
MassRatio=1.
Result = GetKeyWord("MassRatio", Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN BEGIN
	MassRatio = FLOAT(result)
	MassRatio = MassRatio > 1./MassRatio
ENDIF
rhoSym = "!4q!X!De!N"
IF(MassRatio EQ 1.) THEN rhoSym = "!4q!X!Ds!N"
;
; Next get flux-tube dimensions
;
Lx=0
Result = GetKeyWord("Lx", Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN BEGIN
	Lx = FLOAT(result)
ENDIF
Ly=0
Result = GetKeyWord("Lx", Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN BEGIN
	Ly = FLOAT(result)
ENDIF
Lz=0
Result = GetKeyWord("Lx", Extra)
IF(Query_Integer(result) + Query_Real(result)) THEN BEGIN
	Lz = FLOAT(result)
ENDIF
;
; Check number of arguments
;
nArgs = N_PARAMS(0)
IF(nArgs GT 1) THEN BEGIN
	MESSAGE, "Called with more than one argument, returning", /INFORMATIONAL
	RETURN, 0
ENDIF

IF(nArgs EQ 0) THEN BEGIN
	IF(TypeOf(Extra) EQ 8) THEN TubeStr = Tube_Data(_Extra=Extra)
	IF(TypeOf(Extra) NE 8) THEN TubeStr = Tube_Data()
ENDIF

ScaleTime = rhoStar_T*SQRT(MassRatio)
Lx = Lx*SQRT(MassRatio)
Ly = Ly*SQRT(MassRatio)
Lz = Lz*RhoStar_T 

histIn = tubeStr.hist

TE = histIn.TE -> Over(rhoStarSq)
TE -> Set, units="(" + rhoSym + "/L!DT!N)!U2!NT"
TE -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"

nSpecies = N_ELEMENTS(histIN.KE)
KE = OBJARR(nSpecies)
FOR ispecies = 0, nSpecies-1 DO BEGIN
	KE[iSpecies] = histIn.KE[iSpecies] -> over(RhoStarSq)
	KE[iSpecies] -> Set, units="(" + rhoSym + "/L!DT!N)!U2!NT"
	KE[iSpecies] -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"
ENDFOR

FE = histIn.FE -> Over(rhoStarSq/MassRatio)
FE -> Set, units="(" + rhoSym + "/L!DT!N)!U2!NT"
FE -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"

rmsPhi = histIn.rmsPhi -> Over(rhoStar_T/SQRT(MassRatio))
rmsPhi -> Set, units="(" + rhosym + "/L!DT!N)(T/e)"
rmsPhi -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"

histIn.avgWsq -> get, title=title, mnemonic=mnemonic
Temp = histIn.avgWsq -> AbsSq()
avgWsq = Temp -> Over(rhoStarSq/Massratio)
avgWsq -> Set, units="(" + rhoSym + "/L!DT!N)!U2!N"
avgWsq -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"
avgWsq -> Set, title=title, mnemonic=mnemonic
Temp -> Trash

pFlux = OBJARR(nSpecies)
FOR ispecies = 0, nSpecies-1 DO BEGIN
	pFlux[iSpecies] = histIn.pFlux[iSpecies] -> over(RhoStarSq/SQRT(MassRatio))
	pFlux[iSpecies] -> Set, units="(" + rhoSym + "/L!DT!N)!4q!X!Ds!Nc!Ds!N"
	pFlux[iSpecies] -> ScaleAxis, 't', const=Scaletime, units="L!DT!N/v!Dt!N"
ENDFOR

eFlux = OBJARR(nSpecies)
FOR ispecies = 0, nSpecies-1 DO BEGIN
	eFlux[iSpecies] = histIn.eFlux[iSpecies] -> over(RhoStarSq/SQRT(MassRatio))
	eFlux[iSpecies] -> Set, units="(" + rhoSym + "/L!DT!N)!4q!X!Ds!Nc!Ds!N"
	eFlux[iSpecies] -> ScaleAxis, 't', const=Scaletime, units="L!DT!N/v!Dt!N"
ENDFOR

nModes = N_ELEMENTS(histIn.mode)
mode = OBJARR(nModes)
FOR iMode = 0, nModes-1 DO BEGIN
	mode[iMode] = histIn.mode[iMode] -> over(RhoStar_T/SQRT(MassRatio))
	mode[iMode] -> Set, units="(" + rhoSym + "/L!DT!N)(T/e)"
	mode[iMode] -> ScaleAxis, 't', const=ScaleTime, units="L!DT!N/v!Dt!N"
ENDFOR


output = {	Name:	"TubeNorm",	$
		TE:	TE,		$
		KE:	KE,		$
		FE:	FE,		$
		rmsPhi:	rmsPhi,		$
		avgWsq:	avgWsq,		$
		pFlux:	pFlux,		$
		eFlux:	eFlux,		$
		mode:	mode		}

ok = Tube_Norm_Field(tubeStr.phi, rhoStar_T/SQRT(MassRatio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "Phi", tubeStr.phi)

ok = Tube_Norm_Field(tubeStr.rho, rhoStar_T/SQRT(Massratio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "rho", tubeStr.rho)

ok = Tube_Norm_Field(tubeStr.dne, rhoStar_T/SQRT(Massratio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "dne", tubeStr.dne)

ok = Tube_Norm_Field(tubeStr.dni, rhoStar_T/SQRT(Massratio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "dni", tubeStr.dni)

ok = Tube_Norm_Field(tubeStr.apa, rhoStar_T/SQRT(Massratio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "apa", tubeStr.apa)

ok = Tube_Norm_Field(tubeStr.upa, rhoStar_T/SQRT(Massratio), ScaleTime, Lx, Ly, Lz, rhoSym)
IF(ok) THEN output = CREATE_STRUCT(output, "upa", tubeStr.upa)


RETURN, output


END ; ****** Tube_Norm ****** ;
