FORWARD_FUNCTION Omega_GAM
FORWARD_FUNCTION GKV_Create_Object
FUNCTION GKV_TempestData, FileID=FileID, RunID=RunID, PATH=path, R_0=R_0in,	$
			tau_e = tau_ein
;
; Reads IDL "dat" files supplied by XU and returns GKV objects 
; containing resulting data.
;
CD, Current=WorkingDirectory
IF(N_ELEMENTS(path) NE 0) THEN CD, path
; find .dat file:
thisFile = DIALOG_PICKFILE(GET_PATH=filePath)
RESTORE, thisFile
FilePieces = STR_SEP(thisFile, "/")
nPieces = N_ELEMENTS(FilePieces)
shortName = FilePieces[nPieces-1]
q=epsilon*Btor/Bpol
phiInfo = SIZE(phi_xyt)
nx=phiInfo(1)
ny=phiInfo(2)
nt=phiInfo(3)
phi_xyt = REFORM(phi_xyt, 1, nx,ny,nt, /OVERWRITE)
xValues = REFORM(  xr, /OVERWRITE)
IF(N_ELEMENTS(xValues) GT nx) THEN xValues=xValues[0:nx-1] 
yValues = REFORM(xthe, /OVERWRITE)
IF(N_ELEMENTS(yValues) GT ny) THEN yValues=yValues[0:ny-1] 
tValues = REFORM(tArr, /OVERWRITE)
IF(N_ELEMENTS(tValues) GT nt) THEN tValues=tValues[0:nt-1]
phiObj=GKV_Create_Object(phi_xyt, tValues=tValues, title="!4u!X", Mnemonic="phi")  
phiObj -> Set, units="V", CodeName="Tempest", CodePI="X.Q. Xu"
PhiObj -> ScaleAxis, theta=xValues, title="x", mnemonic="x", units="m"
PhiObj -> ScaleAxis,  zeta=yValues, title="!4h!X", mnemonic="theta", units=""
RunID_String = "q="+STRING(q, FORMAT='(F5.3)')+", epsilon=" +String(epsilon, FORMAT='(F5.3)')
phiObj -> Set, FileID=shortName, RunID=RunID_String
IF(N_ELEMENTS(FileID) NE 0) THEN PhiObj -> Set, FileID=FileID
IF(N_ELEMENTS( RunID) NE 0) THEN PhiObj -> Set,  RunID=FileID
phiObj -> Set, Axis='t', gridUnits='s'
phiObj -> get, axis=1, range=range
Lx = (range[1] - range[0])*FLOAT(nx+1)/FLOAT(nx)
kx = 2.*!PI/Lx
rho_i = 1.02*SQRT(2.*Mass*Ti)/(1.e4*Btor)
kx = kx*rho_i
R_0 = 1.71
IF(N_ELEMENTS(R_0in) NE 0) THEN R_0=R_0in
v_th = 9.79e3*SQRT(2.*Ti/Mass)
Omega_transit = v_th/(R_0*q)
tau_e=1.
IF(N_ELEMENTS(tau_ein) NE 0) THEN tau_e=tau_ein
const1 = 1.+(23.+16.*tau_e+4.*tau_e^2)/(q^2*(7.+4.*tau_e)^2)
Omega_GAMhat = SQRT(7.+4.*tau_e)/2.*q*const1
Omega_GAM = Omega_GAMhat*Omega_transit

ComplexOmegaHat=Omega_GAM(q=q, tau_e=1, krho=kx)
ComplexOmega = ComplexOmegaHat*Omega_transit
eye = COMPLEX(0.,1.)

const2 = 1./( 1. + 2.*(23./4.+4.*tau_e+tau_e^2)/(q^2*(7./2.+2.*tau_e)^2) )
const3 = Omega_GAMhat^4 + (1.+2.*tau_e)*Omega_GAMhat^2
const4 = Omega_GAMhat^6/64. + (1.+0.375*tau_e)*(0.125*Omega_GAMhat^4+0.75*Omega_GAMhat^2)
gamma_GAMhat = EXP(-omega_GAMhat^2)*const3 + 0.25*(kx*q)^2*exp(-0.25*Omega_GAMhat^2)*const4
gamma_GAMhat = -0.5*SQRT(!PI)*q^2*const2*gamma_GAMhat
gamma_GAM = gamma_GAMhat*Omega_transit

S=0.5*(3.27+SQRT(epsilon)+0.722*epsilon+(0.692-0.722)*epsilon/q^2)
resid =  1./(1.+S*q^2/SQRT(epsilon))

temp1   = phiObj -> slice(theta=0)
phidata = temp1 -> slice(x=3.*Lx/4.)
norm = phiData -> slice(t=0)
nPhiData = phiData -> over(norm)
nPhiData -> scaleAxis, 't', /uniform
nPhiData -> Get, ErrorBars=ErrorBars
PTR_FREE, ErrorBars
ErrorBars = PTR_NEW()
nPhiData -> Set, ErrorBars=ErrorBars, title="!4u!X/!4u!X!D0!N", units=""
avgPhi = nPhiData -> avg('t')
avgPhi -> Get, ErrorBars=ErrorBars
PTR_FREE, ErrorBars
ErrorBars = PTR_NEW()
avgPhi -> Set, ErrorBars=ErrorBars
dPhidata = nphiData -> DeTrend()
norm -> trash
temp1 -> trash
temp1 = nphidata -> times(0.)
temp1 -> set, units=""
avgPhi -> get, title=title, mnemonic=mnemonic
avgPhi -> set, units=""
avgPhi = temp1 -> Plus(avgPhi)
avgPhi -> set, title=title, mnemonic=mnemonic, units=""
temp2 = temp1 -> plus('t')
temp3 = temp2 -> times(-eye*ComplexOmega)
temp4 = temp3 -> Execute("exp")
temp5 = temp4 -> Execute("FLOAT")
phidata   -> trash
residObj = temp1 -> PLUS(resid)
residObj -> set, title="!4u!D!9$!X!N/!4u!X!D0!N", mnemonic="Phi", units="",	$
			runID="Xiao, Catto Theory", FileID=FileID

temp6 = temp5 -> Times(1.-resid)
temp7 = temp6 -> Plus(resid)
FileIDfield = "k!Dx!N!4q!X!Di!N=" + STRING(kx, FORMAT='(F5.3)')
temp7 -> set, title="!4u!X/!4u!X!D0!N", mnemonic="Phi", units="",	$
			runID="Sugama, Watanabe Theory", FileID=FileIDfield
nPhiData -> Set, FileID=FileIDfield
dPhiData -> Set, FileID=FileIDfield
temp1 -> trash
temp2 -> trash
temp3 -> trash
temp4 -> trash
temp5 -> trash
temp6 -> trash


Output = {	Name	: "GKV_TempestData",	$
		File	:	thisFile,	$
		epsilon	:	epsilon,	$
		Bpol	:	Bpol,		$
		Btor	:	Btor,		$
		q	:	q,		$
		Mass	:	Mass,		$
		Ti	:	Ti,		$
		kx	:	kx,		$
		Lx	:	Lx,		$
		Omega_GAM:	Omega_GAM,	$
		Gamma_GAM:	Gamma_GAM,	$
		COmega	:	ComplexOmega,	$
		Omegatransit:	Omega_Transit,	$
		Phi	:	PhiObj,		$
		phidata	:	nphidata,	$
		theory	:	temp7,		$
		resid	:	residObj,	$
		avgPhi	:	avgPhi,		$
		dPhiData:	dPhiData	}
		
CD, WorkingDirectory
RETURN, output
END  ;  ****** GKV_TempestData ******  ; 
