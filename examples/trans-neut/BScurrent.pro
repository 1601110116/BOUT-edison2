 pro bscurrent,Ni_in,Ti_in,Te_in,Jpar_BS0,time=time,path=path,grid=grid,q95_input=q95_input,Aratio=Aratio

;;Bootstrap current is calculated based on the references 
;;O. Sauter, C. Angioni, and Y. R. Lin-Liu, Phys. Plasmas 6, 2834 (1999);
;; and ibid. 9, 5140 (2002)::some corrections to PoP 6, 2834 (1999).
;;

if ~keyword_set(path) then path='data'
if ~keyword_set(grid) then grid=file_import("data/circle.grd.hl2a.nc")
if ~keyword_set(q95_input) then q95_input = 5.
if ~keyword_set(Aratio) then Aratio = 0.35

if ~keyword_set(time) then BEGIN

print,'Error!, time has to be given as same as the time when P0 was taken!'

time= 0
print, 'time set to be ZERO!'
ENDIF

g=grid

print,'collecting normalized parameters'
Ni_x=collect(path=path,var="Ni_x")
Te_x=collect(path=path,var="Te_x")
Ti_x=collect(path=path,var="Ti_x")
Lbar=collect(path=path,var="Lbar")
tbar=collect(path=path,var="tbar")

PI = 3.14159265
MU0 = 4.0e-7*PI     
KB = 1.38065e-23      ;     // Boltamann constant
eV_K = 11605.0        ;         // 1eV = 11605K
Zi = 1.               ;    charge state
AA = 2.
density_unit= 1.0e19       ; 1/m^3

Ni=reform(Ni_in[*,*,*,time])
Ti=reform(Ti_in[*,*,*,time])
Te=reform(Te_in[*,*,*,time])

;;;Normalization
  g2=g
  g2.psixy /= Lbar*Lbar*g.bmag; normalized to g2
  g2.Bxy /= g.bmag
  g2.Btxy /= g.bmag
  g2.Rxy /= lbar

      q95=q95_input
      pei= Ni*(Te+Ti)
      Pe = Ni*Te
      Pi = Ni*Ti
      ZZ = Zi

      lambda_ei = 24.-alog(sqrt(Ni_x)/Te_x);                          // in unit cm
      lambda_ii = 23.-alog(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/Ti_x^1.5);    // unit cm

      nueix     = 2.91e-6*Ni_x*lambda_ei/Te_x^1.5;             // unit 1/s
      nuiix     = 4.78e-8*Ni_x*lambda_ii/Ti_x^1.5/sqrt(AA);         // unit 1/s


      V_th_e= 4.19e5*sqrt(Te*Te_x)*tbar/Lbar;
      V_th_i= 9.79e3*sqrt(Ti*Te_x/AA)*tbar/Lbar;

      nu_estar = 100.*nueix * q95*tbar / (V_th_e) / Aratio^(1.5)
      nu_istar = 100.*nuiix * q95*tbar / (V_th_i) / Aratio^(1.5)

;;; 2D multiply 3D vars

      ft = BS_ft(100,g2);

      s=size(Ni,/dim)   
      ;;;s[0] X, S[1] Y, S[2] Z
      nx = s[0]
      ny = s[1]
      ;nz = s[2]

      ss=size(Ni)
      IF ss[0] EQ 2 THEN BEGIN
         nz = 1
      ENDIF

      IF ss[0] EQ 3 THEN BEGIN
         nz = ss[3]
      ENDIF


      f31 = FLTARR(nx,ny,nz)
      f32ee = FLTARR(nx,ny,nz)
      f32ei = FLTARR(nx,ny,nz)
      f34 = FLTARR(nx,ny,nz)
      BSal0 = FLTARR(nx,ny,nz)
      BSal = FLTARR(nx,ny,nz)
      Jpar_BS0 = FLTARR(nx,ny,nz)

      FOR k=0,nz-1 DO BEGIN


      f31[*,*,k] = ft[*,*] / (1.+(1.-0.1*ft[*,*])*sqrt(nu_estar[*,*,k]) + 0.5*(1.-ft[*,*])*nu_estar[*,*,k]/Zi);
      f32ee[*,*,k] = ft[*,*] / (1.+0.26*(1.-ft[*,*])*sqrt(nu_estar[*,*,k]) + 0.18*(1.-0.37*ft[*,*])*nu_estar[*,*,k]/sqrt(Zi));
      f32ei[*,*,k] = ft[*,*] / (1.+(1.+0.6*ft[*,*])*sqrt(nu_estar[*,*,k]) + 0.85*(1.-0.37*ft[*,*])*nu_estar[*,*,k]*(1.+Zi));
      f34[*,*,k] = ft[*,*] / (1.+(1.-0.1*ft[*,*])*sqrt(nu_estar[*,*,k]) + 0.5*(1.-0.5*ft[*,*])*nu_estar[*,*,k]/Zi);
      ENDFOR

      L31 = F31(f31,Zi) ;
      L32 = F32ee(f32ee,Zi)+F32ei(f32ei,Zi) ;
      L34 = F31(f34,Zi) ;

      FOR k=0,nz-1 DO BEGIN

         BSal0[*,*,k] = - (1.17*(1.-ft[*,*]))/(1.-0.22*ft[*,*]-0.19*ft[*,*]^2);
         BSal[*,*,k] = (BSal0[*,*,k]+0.25*(1-ft[*,*]^2)*sqrt(nu_istar[*,*,k]))/(1.+0.5*sqrt(nu_istar[*,*,k])) + 0.31*nu_istar[*,*,k]^2*ft[*,*]^6;
         BSal[*,*,k] *= 1./(1.+0.15*nu_istar[*,*,k]^2*ft[*,*]^6);
      ENDFOR

      dpeidx = dydx(pei,g2.psixy)
      dTedx = dydx(Te,g2.psixy)
      dTidx = dydx(Ti,g2.psixy)
      ;Jpar_BS0 = L31* DDX(pei)/Pe  + L32*DDX(Te)/Te + L34*DDX(Ti)/(Zi*Te)*BSal;
      Jpar_BS0 = L31* dpeidx/Pe  + L32*dTedx/Te + L34*dTidx/(Zi*Te)*BSal;

      FOR k=0,nz-1 DO BEGIN
         Jpar_BS0[*,*,k] *=  -g2.Rxy[*,*]*g2.Btxy[*,*]*Pe[*,*,k]*(MU0*KB*Ni_x*density_unit*Te_x*eV_K)/(g2.Bxy[*,*]*g2.Bxy[*,*])/(g.bmag*g.bmag) ; //NB:   J_hat = MU0*Lbar * J / mesh->Bxy;

      ENDFOR

END


FUNCTION DYDX, var, xval
  
  s = SIZE(var,/dim)

  ; 3D [x,y,z]
  nx = s[0]
  ny = s[1]
  
  ss=size(var)
  IF ss[0] EQ 2 THEN BEGIN
     nz = 1
  ENDIF

  IF ss[0] EQ 3 THEN BEGIN
     nz = ss[3]
  ENDIF



  result = FLTARR(nx,ny,nz)

    FOR j=0,ny-1 DO BEGIN
      FOR k=0,nz-1 DO BEGIN

       result[*,j,k] = deriv(xval[*,j],var[*,j,k])

    ENDFOR
ENDFOR
  
  RETURN, result
END



FUNCTION BS_ft, index,g


  result1=0. ;
  
  ;dxlam = 1./g.bmag/index     ;  maybe g.bmag should be normalized
  dxlam = 1./index     ;  maybe g.bmag should be normalized

  xlam = 0.;

  FOR i=0, index DO BEGIN
    
      result1 += xlam*dxlam/sqrt(1.-xlam*g.Bxy) ;
      xlam += dxlam;
    
 ENDFOR

  result = 1.- 0.75*g.Bxy*g.Bxy * result1 ;

  RETURN, result

END

FUNCTION F31, input, Zi


  result = ( 1 + 1.4/(Zi+1.) ) * input;
  result -= 1.9/(Zi+1.) * input*input;
  result += 0.3/(Zi+1.) * input*input*input;
  result += 0.2/(Zi+1.) * input*input*input*input;

  RETURN, result

END

FUNCTION F32ee,input, Zi


  result = (0.05+0.62*Zi)/(Zi*(1+0.44*Zi))*(input-input*input*input*input);
  result +=1./(1.+0.22*Zi)*( input*input-input*input*input*input-1.2*(input*input*input-input*input*input*input) );
  result += 1.2/(1.+0.5*Zi)*input*input*input*input;

  RETURN, result
END

FUNCTION F32ei, input, Zi

  
  result = -(0.56+1.93*Zi)/(Zi*(1+0.44*Zi))*(input-input*input*input*input);
  result += 4.95/(1.+2.48*Zi)*( input*input-input*input*input*input-0.55*(input*input*input-input*input*input*input) );
  result -= 1.2/(1.+0.5*Zi)*input*input*input*input;

  RETURN, result

END
