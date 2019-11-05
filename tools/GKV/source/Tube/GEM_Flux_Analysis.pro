FUNCTION GKVs2D::SquashArray, arg
;
; Purpose:
;
;	Acts on fluctuation intensity on same
;	grid as objects to be "squashed". 
;	Fluctuation intensity is "self",
;	objects to be "squashed" are in
;	'arg' which should be an object
;	array.
;
;  Argument;
;
;	An 1D array of objects on same grid as "self"
;
;  Output:
;
;	An object array of same dimension as 'arg'
;	containing 'arg' squashed with self.
;
; Written by W.M. Nevins
;	8/12/08
;
nObjs = N_ELEMENTS(arg)
output = OBJARR(nObjs)
FOR i=0,nObjs-1 DO output[i] = arg[i] -> Squash(self)
RETURN, output 
END ; ****** SquashArray ****** ;


FUNCTION GKVs2D::GEM_Flux_Analysis, arg, _Extra=Extra
;
;  Purpose:
;
;	Acts on fluctuation intensity vs. (x,t) 
;	as output from GKVs3D::DeltaSq. Argument
;	should be the "flux" structure produced
;	by GEM_Data (stored in the tag "FLUX").
;	This routine averages the fluctuation
;	intensity radially over the same radial
;	bins used in computing GEM's fluxes. 
;
; Argument:
;
;	"FLUX" structure found under the tag "FLUX"
;	in strucutre produced by GEM_Data. 
;
;  Keywords:
;	
;	ForceUnits	Set this keyword to force output
;			intensity to have the same
;			time units as 'arg'	
;
;  Output:
;
;	Intensity averaged over same radial bins
;
;	For each flux in Flux Structure, a "SQUASH"
;	plot showing flux vs intensity.
;
;  Written by W.M. Nevins
;	9/12/08
;
; Check for valid argument
;
 IF(TypeOF(arg) NE 8) THEN BEGIN
	MESSAGE, 'Argument not a structure, RETURNING', /INFORMATIONAL
	return, 0
ENDIF
FluxTags = TAG_NAMES(arg)
nTags = N_TAGS(arg)
IF(nTags NE 3) THEN BEGIN
	MESSAGE, 'Argument does not have 3 tags, RETURNING', /INFORMATIONAL
	RETURN, 0
ENDIF
nSpecies = N_ELEMENTS(arg.energy)
particle = arg.particle
energy = arg.energy
firstEnergy = energy[0]
;
; Create averaged intensity object. that is,
; a 2D object containing 'self' averaged over same
; radial bins as flux arrays contained in 'arg' but
; sampled onto the same time grid as 'self'
;
avgIstr = {GKVs2D}
FOR i=0,10 DO avgIstr.(i) = self.(i)
AvgIstr.Grid1 = GKVsd_GridCopy(firstEnergy.grid1)
AvgIstr.Grid2 = GKVsd_GridCopy(self.grid2)
result = GetKeyWord('ForceUnits', Extra)
IF(result NE 'undefined') THEN BEGIN
	IF(result EQ 1) THEN AvgIstr.Grid2.units = firstEnergy.grid2.units
ENDIF 
xValues = *(avgIstr.grid1.values)
nx = N_ELEMENTS(xValues)
xMin = xValues[0]
xMax = xValues[nx-1]
tValues = *(AvgIstr.grid2.values)
nt = N_ELEMENTS(tValues)
Dx = xMax - xMin
binSeps = xMin + Dx/(nx)*FINDGEN(nx+1)
values = FLTARR(nx, nt)
FOR i=0, nx-1 DO BEGIN
	avgI_Obj = self -> Avg(axis=1, range=binSeps[i:i+1])
	values[i,*] = *(avgI_Obj.values)
	avgI_Obj -> TRASH
ENDFOR
avgIstr.values = PTR_NEW(values)
avgIstr.errorBars = PTR_NEW()
indices = *(self.indices)
avgIstr.indices = PTR_NEW(indices)
avgI = OBJ_NEW('GKVs2D', avgIstr)
output = {	Name	:	"GEM_Flux_Analysis",	$
		avgI	:	avgI			}	
;
; now "squash" fluxes with intensity
;
xText = STRARR(nx)
FOR ix=0,nx-1 DO BEGIN
	xRange = STRING(binSeps[ix:ix+1], FORMAT='("x=[", F6.2,", ", F6.2, "]")')
	xtext[ix] = STRCOMPRESS(xRange)
ENDFOR
particleOut = OBJARR(nSpecies, nx)
energyOut   = OBJARR(nSpecies, nx)
FOR iSpecies=0, nSpecies-1 DO BEGIN
	FOR i=0,nx-1 DO BEGIN
		ITemp = avgI -> slice(axis=1, value=xValues[i])
		pTemp = particle[iSpecies] -> slice(axis=1, value=xValues[i])
		pTemp1= pTemp -> Interpolate(ITemp)
		eTemp =   energy[iSpecies] -> slice(axis=1, value=xValues[i])
		eTemp1= eTemp -> Interpolate(ITemp)
		particleOut[iSpecies, i] = pTemp1 -> Squash(ITemp)
		  energyout[iSpecies, i] = eTemp1 -> Squash(ITemp)
		ParticleOut[iSpecies, i] -> set, indices=PTR_NEW(['*', xtext[i]])
		  energyOut[iSpecies, i] -> set, indices=PTR_NEW(['*', xtext[i]])
		gkvDelete, Itemp, pTemp, pTemp1, eTemp, eTemp1
	ENDFOR
ENDFOR
output = CREATE_STRUCT(output, 'Particle', particleOut, 'Energy', energyOut)

RETURN, output
END ; ****** GKVs2D::GEM_Flux_Analysis ****** ;  

 


