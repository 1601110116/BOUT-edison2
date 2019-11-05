;===========================================================================;
pro generate_grid, z, n=n, zmax=zmax, debug=debug
;
COMMON grdcmn, grid

if not keyword_set(N) then n=64
if not keyword_set(Zmax) then Zmax=1e0

dz=Zmax/(N-1)
z=dz*FINDGEN(n)


;;-array of mode vectors on this grid

;;-first method
; N21 = N/2 + 1 ; Midpoint+1 is the most negative frequency subscript:
; k = INDGEN(N) ; The array of subscripts:
; k[N21] = N21 -N + FINDGEN(N21-2) ; Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):
; k = k/(!PI*dz) ; Compute T0 frequency:
; k1=SHIFT(k, -N21) ; Shift so that the most negative frequency is first:

;;;-second method
;kMin=-!PI/dz
;kMax=!PI/dz
;k2=kmin+(kMax-kMin)*FINDGEN(n)/(n-1)

;;grid={n:n, z:z, dz:dz, k:k1}
grid={n:n, z:z, dz:dz, zmax:zmax}

;
if keyword_set(DEBUG) then STOP
end



function Tprofile, z, option=option
;-generate T(z) profile
COMMON grdcmn, grid

if not keyword_set(OPTION) then OPTION=1

z=grid.z
midz=MEDIAN(grid.z)

case option of

    1: begin
        Tz=exp(-1e2*(grid.z-midz)^2)
    end

    10: begin
        Tz=exp(-1e2*((grid.z-midz)/midz)^2)
    end

    11: begin
        Tz=exp(-1e3*(grid.z-0.9*midz)^2)+exp(-1e3*(grid.z-1.3*midz)^2)
    end

    12: begin
        Tz=exp(-1e3*(grid.z-0.9*midz)^2)+exp(-1e3*(grid.z-1.3*midz)^2)+exp(-1e3*(grid.z-1.0*midz)^2)
    end

    13: begin
        Tz=1*exp(-1e3*(grid.z-1.3*midz)^2)+3*exp(-1e3*(grid.z-1.0*midz)^2)
    end

    2: begin
        Tz=tanh(1e1*(grid.z-0.5))+1e0
    end

    21: begin
        Tz=tanh(1e2*(grid.z-0.5))+1e0
    end

    3: begin
        Tz= (-tanh(1e1*(grid.z-0.5))+1)/2
    end

    31: begin
        Tz= (-tanh(1e2*(grid.z-0.5))+1)/2
    end

    32: begin
        Tz= (-tanh(1e2*(grid.z-0.65))+1)/2
    end

endcase



;
return, Tz
end


function signum, k, approx=approx, imin=imin, imax=imax
;
; The sign function - exact and approximate
;
;Test:
;IDL> kmin=1e-4 & kmax=1e9 & k=10^(alog10(kmin)+alog10(kmax/kmin)*findgen(101)/100)
;IDL> plot, k, signum(k,/app),/xl, chars=2 & oplot, k,k*0+1, lin=2
;========================================;

if keyword_set(APPROX) then begin

    if not keyword_set(IMIN) then IMIN=0
    if not keyword_set(IMAX) then IMAX=7

    res=0d0
    alp=5d0
    bet=1.04d0


    for n=iMin, iMax do begin
    ;;for n=3,3 do begin
    ;;for n=2,2 do begin

        alpn=alp^n

        res_now=bet*alpn*k/(1d0*k^2+alpn^2)
        res=res+res_now

    endfor

endif else begin

    eps=1d-20
    res=k/(ABS(k)+eps)

endelse

;
;
return, res
end




function qFourier, Tz, local=local, approx=approx, imin=imin, imax=imax, debug=debug
;
; Calculate q(z) profile by Fourier transform
;---------------------------------------------
COMMON grdcmn, grid


f=FFT(Tz,/double)

n=N_ELEMENTS(Tz)
if (n mod 2) eq 1 then nmid=n/2 else nmid=n/2-1

u=SHIFT(dindgen(n)-nmid,-nmid)

if keyword_set(LOCAL) then begin
 
    ;;-this is the just derivative dT/dz calculated spectrally
    qz=-(2*!DPI/n)*FFT((complex(0,1)*u*f),/inv,/double)/grid.dz 

endif else begin

    ;;qz=-(2*!DPI/n)*FFT((complex(0,1)*Signum(u,approx=approx)*f),/inv,/double)
    qz=-(2*!DPI)*FFT((complex(0,1)*Signum(u,approx=approx, imin=imin, imax=imax)*f),/inv,/double)/grid.zmax
    ;;qz=-FFT((complex(0,1)*Signum(u,approx=approx)*f),/inv,/double)

endelse

;
if keyword_set(DEBUG) then STOP
return, qz
end




function qLorentz, Tz, local=local, imin=imin, imax=imax, debug=debug
;
; Calculate q(z) by finite-difference using representation by Lorentzians
;---------------------------------------------
COMMON grdcmn, grid

if not keyword_set(IMIN) then IMIN=0
if not keyword_set(IMAX) then IMAX=7


nz=N_ELEMENTS(Tz)
qz=fltarr(nz)  ;;-total qz


if keyword_set(LOCAL) then begin

    qz=-DERIV(grid.z, Tz)

endif else begin

    alp=5d0
    bet=1.04d0


    for n=iMin, iMax do begin

        alpn=alp^n

        z = grid.z
        ;;Pcoef=(1/40.)*REPLICATE(-1d0,nz) ;;-d2/dz2 term - fudge-factor 1/40 needed???
        Pcoef=((grid.zmax/(2*!DPI))^2)*REPLICATE(-1d0,nz) ;;-d2/dz2 term - fudge-factor 1/40 needed???
        ;;Pcoef=REPLICATE(-1d0,nz)
        Qcoef=REPLICATE(alpn^2,nz)
        Rcoef=-bet*alpn*DERIV(z,Tz)

        ;;STOP

        ;;-solve Helmholtz equation
        Helm, nz, z, Pcoef, Qcoef, Rcoef, res_now
        qz=qz+res_now
        
    endfor

endelse


;
if keyword_set(DEBUG) then STOP
return, qz
end




pro glf, nsize=nsize, option=option, approx=approx, zmax=zmax, local=local, norm=norm
;
;
;
COMMON grdcmn, grid

Generate_Grid, n=nsize, zmax=zmax

;;-generate some T(z) profile
tz=TProfile(option=option)

;;-calculate q(z) by exact spectral method
qzs=QFourier(tz, local=local)

;;-calculate q(z) by spectral method with Lorentzian approximation
qzl=QFourier(tz, local=local, /approx)

;;-calculate q(z) by finite-difference method with Lorentzian approximation
qzl2=QLorentz(tz, local=local)


tek_color

if keyword_set(NORM) then yrange=[-1,1]

!p.multi=[0,1,2,0,0]
var=tz
if keyword_set(NORM) then var=var/MAX(ABS(var))
plot, grid.z, var, yr=yrange, xtit='z', tit="T(z)", chars=2

var=qzs
if keyword_set(NORM) then var=var/MAX(ABS(var))
plot, grid.z, var,/nodata, xtit='z', tit="q(z)", chars=2
oplot, grid.z, var, col=2

var=qzl
if keyword_set(NORM) then var=var/MAX(ABS(var))
oplot, grid.z, var, col=3

var=qzl2
if keyword_set(NORM) then var=var/MAX(ABS(var)) else var=var
oplot, grid.z, var, col=4, lin=2



XYOUTS, 0.25, 0.40,/norm, "q(z) spectral exact", chars=1, col=2
XYOUTS, 0.25, 0.38,/norm, "q(z) spectral w/ Lorentzians", chars=1, col=3
XYOUTS, 0.25, 0.36,/norm, "q(z) finite-difference w/ Lorentzians", chars=1, col=4
!p.multi=0


;
;
;
end


pro testq, nsize=nsize, opt=opt, zmax=zmax
;
;
;


Generate_Grid, n=nsize, zmax=zmax
tz=TProfile(z, opt=opt)

;;-exact spectral result
qf0=QFOURIER(tz)

imin=0
imax=7

;;-approximate spectral result
qf=QFOURIER(tz,/app, imin=imin, imax=imax)

;;-finite-difference result
ql=QLORENTZ(tz, imin=imin, imax=imax)





plot,  z, qf0, xtit='z', chars=2
oplot, z, qf, col=2
oplot, z, ql, col=3,lin=2


XYOUTS, /norm, 0.55, 0.25, 'spectral - exact', chars=2
XYOUTS, /norm, 0.55, 0.23, 'spectral - Lorentz', chars=2, col=2
XYOUTS, /norm, 0.55, 0.21, 'finite-difference - Lorentz', chars=2, col=3

print, "error=", MAX(ABS(ql-qf))


;
;
;
end
