pro helm, nz, z, Pcoef, Qcoef, Rcoef, f, BC=BC, debug=debug
;
; Solve Helmholtz equation: [P(z)*d2/dz2 + Q(z)] f(z) = R(z)
;
; Inputs: array size nz
;         real arrays z, Pcoef=P(z), Qcoef=Q(z), Rcoef=R(z)
;
; Output: array f
;===============================================================;


;;-1 for Neumann, 2 for Dirichlet
if not keyword_set(BC) then bc={left:2, right:2}

;;-form a real (complex?) matrix
m=dblarr(nz,nz)
rhs=dblarr(nz)



for i=1,nz-2 do begin

    ;;-2nd derivative term
    dz=0.5*(z[i+1]-z[i-1])
    m[i,i]   = -2d0*Pcoef[i]/dz^2
    m[i+1,i] =  1d0*Pcoef[i]/dz^2
    m[i-1,i] =  1d0*Pcoef[i]/dz^2


    ;;-free term
    m[i,i]=m[i,i]+Qcoef[i]


    ;;-right-hand side
    rhs[i]=Rcoef[i]

endfor


if (bc.left eq 1) then begin
    ;;-use zero-gradient B.C. on left end
    i=0
    m[i,i]=1d0
    m[i+1,i]=-1d0
    rhs[i]=0d0
endif else begin
    ;;-use zero-value B.C. on left end
    i=0
    m[i,i]=1d0
    rhs[i]=0d0
endelse


if (bc.right eq 1) then begin
    ;;-use zero-gradient B.C. on right end
    i=nz-1
    m[i,i]=1d0
    m[i-1,i]=-1d0
    rhs[i]=0d0
endif else begin
    ;;-use zero-value B.C. on right end
    i=nz-1
    m[i,i]=1d0
    rhs[i]=0d0
endelse

;;STOP

;;-solve REAL linear system
LUDC, m, index  
f = LUSOL(m, index, rhs, /DOUBLE)  


;
;
;
if keyword_set(DEBUG) then STOP
end


pro helmtest, nz=nz, opt=opt, debug=debug
;
;
;

if not keyword_set(NZ) then nz=10

if not keyword_set(OPT) then OPT=1


CASE OPT of 
    1: begin
        z=1d0*findgen(nz)/(nz-1)
        Pcoef=1d0 + dblarr(nz)
        Qcoef=2d0 + dblarr(nz)
        g=exp(-10*(z-MEDIAN(z))^2)*cos(z*1e1)
    end

    11: begin
        z=2d0*findgen(nz)/(nz-1)
        Pcoef=1d0 + dblarr(nz)
        Qcoef=2d0 + dblarr(nz)
        g=exp(-10*(z-MEDIAN(z))^2)*cos(z*1e1)
    end

    2: begin
        z=1d0*findgen(nz)/(nz-1)
        Pcoef=1d0 + dblarr(nz)
        Qcoef=0d0 + dblarr(nz)
        g=exp(-20*(z-MEDIAN(z))^2)*cos(z*1e2)
    end


    3: begin
        z=1d0*findgen(nz)/(nz-1)
        Pcoef=1d0 + dblarr(nz)
        Qcoef=0d0 + dblarr(nz)

        midZ=MEDIAN(z)
        envel=exp(-20*(z-midZ)^2)
        envel(where(ABS(z-midZ) lt 0.25))=1.0
        g=envel*SIN(z*1e2)
    end

    ELSE: begin
        print, "Allowed options 1-3, using option 1..."

        z=1d0*findgen(nz)/(nz-1)
        Pcoef=1d0 + dblarr(nz)
        Qcoef=2d0 + dblarr(nz)
        g=exp(-10*(z-MEDIAN(z))^2)*cos(z*1e1)        
    end

ENDCASE

Helm, nz, z, Pcoef, Qcoef, g, f, debug=debug

!p.multi=[0,1,2,0,0]
tek_color
plot, z,f, tit='Solution', xtit='z', chars=2, yticks=4
plot, z,g, tit='RHS', xtit='z', chars=2
oplot, z, Pcoef*(DERIV(z,DERIV(z,f)))+Qcoef*f, lin=2, col=2
!p.multi=0

;
;
;
end
