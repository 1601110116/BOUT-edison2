FUNCTION UmanskyData, inputStructure
;
; This function accepts an IDL structure containing BOUT data
; and returns a structure containing GKV objects.
;
rGrid={grid}
rGrid.mnemonic="r"
rGrid.title="r"
rGrid.units="cm"
nValues = inputStructure.ngx
rmin=15.
rmax=45.
rValues = rmin + (rmax-rmin)*FINDGEN(nValues)/(nValues-1.)
rGrid.values=PTR_NEW(rValues)
rgrid.uniform = 1b
rGrid.boundary = "open"
rGrid.range=[rmin,rmax]
rGrid.iRange=[0,nValues-1]

yGrid = {grid}
yGrid.mnemonic='y'
yGrid.title='y'
yGrid.units='cm'
ymin=0.
ymax=1800.
nValues=inputStructure.ngy
yValues = ymin + (ymax-ymin)*FINDGEN(nValues)/(nValues-1.)
yGrid.values=PTR_NEW(yValues)
yGrid.uniform=1b
yGrid.boundary="periodic (open)"
yGrid.range=[ymin,ymax]
yGrid.iRange = [0,nValues-1]

thetaGrid={grid}
thetaGrid.mnemonic="theta"
thetaGrid.title="!4h!3"
thetaGrid.units=""
thetaMin=0.
thetaMax=!PI/4.
nValues=inputStructure.ngz
thetaValues=thetaMin + (thetaMax-thetaMin)*FINDGEN(nValues)/(nValues-1.)
thetaGrid.values=PTR_NEW(thetaValues)
thetaGrid.boundary="periodic (closed)"
thetaGrid.uniform=1b
thetaGrid.range=[thetaMin, thetaMax]
thetaGrid.irange=[0,nValues-1]

tGrid={grid}
tGrid.mnemonic='t'
tGrid.title='t'
tGrid.units='s'
omega_ci=inputStructure.wci
tValues=inputStructure.t_array/omega_ci
tMin=MIN(tValues)
tmax=MAX(tValues)
nValues=N_ELEMENTS(tValues)
tGrid.values=PTR_NEW(tValues)
tGrid.boundary="open"
tGrid.uniform=1b
tGrid.range=[tmin,tmax]
tGrid.irange=[0,nValues-1]

nStr={GKVs4D}
nStr.mnemonic="n_i"
nStr.title="n!Di!N"
nStr.indices=PTR_NEW(["*","*","*","*"])
nStr.units="cm-3"
nStr.values=PTR_NEW(inputStructure.ni_xyzt*inputStructure.ni_x)
nStr.codeName="BOUT"
nStr.codePI="M. Umansky"
nStr.grid1=rGrid
nStr.grid2=yGrid
nStr.grid3=thetaGrid
nStr.grid4=tGrid

nObj=OBJ_NEW("GKVs4D", nStr)

phiObj = nObj -> MakeCopy(/NoValues)
phiObj -> set, title="!4u!3", mnemonic="phi", units='V'
phiValues = inputStructure.phi_xyzt*inputStructure.Te_x
phiObj -> set, values=PTR_NEW(phiValues)

output = { n_i : nObj, phi: phiObj }
return, output
END