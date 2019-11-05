FUNCTION D_GAM_SK, arg, Complex=Complex
;
; The zeros of this function yield 
; the GAM dispersion relation. After
;
; H. Sugama and T.H. Watanabe, J. Plasma Physics (2006)
;
; Written by W.M. Nevins
;	3/23/2006
COMMON DD_GAM, q, Tau, sigma
eye = COMPLEX(0.d, 1.d)
one = COMPLEX(1.d, 0.d)

CASE N_ELEMENTS(complex) OF
	1:	omega = arg
	else:	omega  = COMPLEX(arg[0], arg[1])
ENDCASE

omega2 = omega*omega
omega3 = omega2*omega
omega4 = omega3*omega
Z = Z_plasma(omega)
term_1 = 2.*omega3 + 3.*omega + (2.*omega4 + 2.*omega2 + 1.)*Z
term_2 = (omega/2.)*(2.*omega + (2.*omega2 + 1)*Z)^2
denom  = (tau + 1. + omega*Z)
D_GAM = omega + (q*q/2.)*(term_1 - term_2/denom)
IF(N_ELEMENTS(complex) EQ 1)  THEN RETURN, D_GAM
result = [DOUBLE(D_GAM), IMAGINARY(D_GAM)]
RETURN, result
END ; ****** D_GAM_SK ****** ;


FUNCTION GKVs1D_Omega_GAM, q=q_in, tau_e=tau_ein, krho=krho_in
;
; Returns a a GKVs1D object containing the GAM frequency vs. q
;
obj = {GKVs1D}
obj.title = "!4x!X!DGAM!N"
obj.mnemonic = "Omega_GAM"
obj.indices = PTR_NEW(["*"])
obj.units = "v!Dti!N/qR"
obj.CodeName = "GKV"
obj.CodePI = "W.M. Nevins"
obj.runID = "GAM dispersion relation:
obj.FileID= "after Sugama and Watanabe"
grid = {grid}
grid.mnemonic = 'q'
grid.title = 'q'
grid.units = ''

q = 1. + 0.04*FINDGEN(101)
IF(N_ELEMENTS(q_in) GT 1) THEN q=q_in
grid.values = PTR_NEW(q)
qMin = MIN(q, MAX=qMax)
iMax = N_ELEMENTS(q) - 1
grid.range = [qMin, qMax]
grid.irange= [0, iMax]
obj.Grid1 = grid
FORWARD_FUNCTION Omega_GAM
omega = Omega_GAM(q=q, tau_e=tau_ein, krho=krho_in)
obj.values = PTR_NEW(omega)
Omega = OBJ_NEW("GKVs1D", obj)
temp = omega -> EXECUTE("IMAGINARY")
gamma = temp -> times(-1.)
temp -> Trash
gamma -> set, title="-!4c!X!DGAM!N", mnemonic="Gamma_GAM", units="v!Dti!N/qR"
nOmega_GAM = Omega -> over('q')
Omega -> GET, title=title
nOmega_GAM -> Set, title=title, units="v!Dti!N/R"
result = { 	Name: 	"Omega_GAM",	$
		Omega: 	Omega,		$
		Gamma:	gamma,		$
		nOmega:	nOmega_GAM}
RETURN, result
END ; ****** GKVs1D_Omega_GAM ****** ;

FUNCTION Omega_GAM, q=q_in, tau_e=tau_ein, krho=krho_in
;
; Computes GAM dispersion relation after 
;
; Sugama and Watanabe, J. Plasma Phys. <72>, 825 (2006).
; (but see Erratum by same authors in J. Plasma Phys. (2007).
;
; Note that the theory of S&W includes only the first transit
; harmonic. Must use theory of Gao (not included here) to 
; get effects of all transit harmonics. Generally, you need
; to keep about q/2 transit harmonics.
;
; Inputs:
;
;     q	     safety factor (optional). 
;            Defaults to 1.
;
;     tau_e  electron to ion temperature ratio, T_e/T_i
;            (optional). Defaults to 1.
;
;    kRho    Product of radial wave number and ion
;            gyroradius [defined as SQRT(2*T_i*m_i/eB)].
;            (optional).  Defaults to 0.
;
; Modified 1/17/08 by W.M. Nevins
; to include comments on units, etc.
;
; Modified 1/19/08 by W.M. Nevins
; to include corrections to Gamma
; from Erratum by Sugama & Watanabe in J. Plasma Phys. (2007).
;  
;
q=1.
IF(Query_real(q_in)+Query_Integer(q_in)) THEN q=q_in
tau_e=1.
IF(Query_real(tau_ein)+Query_Integer(tua_ein)) THEN tau_e=tau_ein
;
; Compute real part of GAM frequency after S&W, Eq. (2.9)
;
factor_1 = SQRT(7.+4.*tau_e)/2.
factor_2 = 1. + (23.+16.*tau_e+4.*tau_e^2)/(q*(7.+4.*tau_e))^2
omega = factor_1*q*factor_2
;
; comparing with S&W, Eq. (2.9), units must be (v_ti/qR)
; where v_ti=SQRT(2*T_I/m_i) [see text below Eq. (2.7)]
;
; Now compute imaginary part of GAM frequency
; (that is, the GAM damping rate) after S&W, Eq. (2.10)
;
factor_3 = 1. + 2.*(23./4.+4.*tau_e+tau_e^2)/(q*(7./2.+2.*tau_e))^2
omegaSq = omega*omega
factor_4 = omegaSq*omegaSq + (1.+2.*tau_e)*omegaSq
gamma = -SQRT(!PI)/2.*q*q/factor_3*EXP(-omegaSq)*factor_4
IF(N_ELEMENTS(krho_in) EQ 0) THEN RETURN, COMPLEX(omega, gamma)
;
; Check for finite krho
;
krho=0.
IF(Query_real(krho_in)+Query_Integer(krho_in)) THEN krho=krho_in
;
; Expression for Gamma is corrected as per Erratum 
; [Sugama and Watanabe, J. Plasma Physics (2007)]
;
; factor_5=omegaSQ^3/64.+(1.+3./8*tau_e)*((omegaSQ^2)/8.+3.*omegaSq/4.)
; Corrected in Erratum to:
factor_5 = omegaSQ^3/128.+(1.+tau_e)*(omegaSQ^2)/16.		$
           + (3./8.+7./16.*tau_e+5./32.*tau_e^2)*omegaSq

; factor_6=EXP(-omegaSq)*factor_4 + 0.25*(krho*q)^2*EXP(-omegaSq/4.)*factor_5
; Corrected in Erratum to:
factor_6=EXP(-omegaSq)*factor_4 + 0.25*(krho*q)^2*EXP(-omegaSq/4.)*factor_5

gamma = -SQRT(!PI)/2.*q*q/factor_3*factor_6
;
; again comparing with S&W, Eq. (2.10), gamma is also 
; in units of (v_ti/qR) AND 
;      kRho = k_r*v_ti/Omega_ci
;           = k_r*SQRT(2*T_i*m_i)/eB
; 
RETURN, COMPLEX(omega, gamma)
END ; ****** Omega_GAM ****** ;


FUNCTION myERF, arg
;
; Computes Error function from Dawson's integral
;
eye = COMPLEX(0.,1.)
z = eye*arg
result = -2.*eye/SQRT(!PI)*exp(z^2)*Dawson(z)
RETURN, result
END ; ****** myERF ****** ;

FUNCTION Dawson, arg
;
; Computes Dawson's integral
; following "Numerical Recipies in C"
; (see Sec. 6.1)
;
one = COMPLEX(1.0d, 0.0d)
;IF(ABS(arg) LT 0.2) THEN BEGIN
	a1 = 2.0/3.0
	a2 = 0.4
	a3 = 2.0/7.0
	argSq = arg*arg
	result1 = arg*(one - a1*argSq*(one - a2*argSq*(one - a3*argSq) ) )
;	RETURN, result
;ENDIF
h = 0.4d
nMax = 10
c = DCOMPLEXARR(nMax+1)
FOR i=1,nMax DO c[i] = EXP(-((2.0*i-1.0)*h)^2)
xx = DOUBLE(arg)
n0 = 2.*FIX(0.5*xx/h+0.5)
xp = arg - n0*h
e1=exp(2.*xp*h)
e2=e1*e1
d1 = n0 + 1
d2=d1-2.
sum=COMPLEX(0.d,0.d)
FOR i=1, nMax DO BEGIN
	sum = sum + c[i]*( e1/d1 + one/(d2*e1) )
	d1 = d1 + 2.	
	d2 = d2 - 2.
	e1 = e1*e2
ENDFOR
result2 = EXP(-xp*xp)/SQRT(!PI)*sum
result = (ABS(arg) LT 0.2)*result1 + (ABS(arg) GE 0.2)*result2
RETURN, result
END ; ****** Dawson ****** ;


FUNCTION myERFC, arg, itMax=itMax_in
;
; The error function compliment supplied with IDL 6.0
; does not seem to work with complex argument.  This is
; my poor effort to make something that does better.
;
; This function returns an asymptotic approximation
; to EXP(arg^2)*ERFC(arg)
;
; Written by W.M. Nevins
;	3/20/06
;
; Asymptotic expansion after 
; Abramowitz and Stegun Eq. 7.1.23
;
RETURN, EXP(arg^2)*ERFC(arg)
itMax=1
IF(Query_Integer(itMax_in)) THEN itMax = itMax_in
ONE = DCOMPLEX(1.,0.)
zero= DCOMPLEX(0.,0.)
erfc = one
term = one
over2zSQ = one/(2.*arg^2)
FOR i=1, itMax DO BEGIN
	term = -term*(2*i-1)*over2zSq
	erfc = erfc + term
ENDFOR
erfc = (ONE/arg)*erfc
;RETURN, erfc
END ; ****** myERFC ****** ;

FUNCTION Z_Plasma, arg, sigma=sigma_in
;
; Computes the plasma dispersion function
; after Ichimaru, "Basic Principles of Plasma Physics"
; see Sec. 4.1.  See also Chap. 7 of Abramowitz and Stegun
;
;
one = COMPLEX(1.,0.)
eye = COMPLEX(0.,1.)
zeta = arg/SQRT(2.)
sigma = 1
IF(Query_Real(sigma_in) + Query_Integer(sigma_in)) THEN sigma=sigma_in
;
; For small argument use IDL's ERFC
;
zetaSq = zeta*zeta 
Z0 = one - eye*SQRT(!PI)*zeta*EXP(-zetaSq)*(ERFC(-eye*zeta) -2.*sigma)
;
; or express in terms of Dawson's integral
;
Z1 = one - arg*Dawson(zeta) + eye*SQRT(!PI)*zeta*EXP(-zetaSQ)
;
; For large arguments use asymptotic expansion
; (after Ichimaru, "Basic Principles of Plasma Physics", Chp. 4
;
argSq = arg*arg
sum = COMPLEX(0.d, 0.d)*arg
term = 1./argSq
nMax = 3
FOR i=1,nMax DO BEGIN
	sum = sum + (2.*i-1.)*term
	term = term/argSq
ENDFOR
Z2 = eye*SQRT(!PI)*zeta*EXP(-zetaSq) - sum
;
; Choose which version of plasma Z-function to return
;
test = ABS(zeta)
z = (test LT 5)*z0 + (test GE 5)*z2
RETURN, z
END ; ****** Z_Plasma ****** ;   

FUNCTION K_Diamond, p, sigma=sigma_in
;
; Computes the function K(p) from 
; Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
sigma = 0
eye = COMPLEX(0.,1.)
IF(Query_Integer(sigma_in) + Query_Real(sigma_in)) THEN sigma=sigma_in
;
; For small arguments just use IDL's ERFC
;
K1 = 1./SQRT(!PI)*EXP(p*p)*(ERFC(p) - 2*sigma)
;
; For large arguments use asymptotic expansion to ERFC
; from Abromowitz and Stegun, 7.1.23
;
e1  = 1./(2.*p*p)
sum = COMPLEX(0.d,0.d)*p + 1.
term= COMPLEX(1.d, 0.d)*e1
FOR i=1,4 DO BEGIN
	sum = sum + (2.*i-1)*term
	term= -e1*term 
ENDFOR
K2 = (1./!PI)*sum/p ; - 2.*sigma/SQRT(!PI)*EXP(p*p)
test = ABS(p)
k = (test LT 5)*K1 + (test GE 5)*k2
return, K
END ; ****** K_Diamond ****** ;

FUNCTION altK_diamond, p
;
;
; Computes the function K(p) from 
;
; Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
; from Dawson's integral evaluated 
; following "Numerical Recipies in C"
; (see Sec. 6.1)
;
one = COMPLEX(1.0d, 0.0d)
eye = COMPLEX(0.0d, 1.0d)
arg = eye*p
; 
; Compute using small argument expansion
;
a1 = 2.0/3.0
a2 = 0.4
a3 = 2.0/7.0
argSq = arg*arg
result1 = arg*(one - a1*argSq*(one - a2*argSq*(one - a3*argSq) ) )
;
; Compute using method from "Numerical Recipes ...:
;
h = 0.4d
nMax = 10
c = DCOMPLEXARR(nMax+1)
FOR i=1,nMax DO c[i] = EXP(-((2.0*i-1.0)*h)^2)
xx = DOUBLE(arg)
n0 = 2.*FIX(0.5*xx/h+0.5)
xp = arg - n0*h
e1=exp(2.*xp*h)
e2=e1*e1
d1 = n0 + 1
d2=d1-2.
sum=COMPLEX(0.d,0.d)
FOR i=1, nMax DO BEGIN
	sum = sum + c[i]*( e1/d1 + one/(d2*e1) )
	d1 = d1 + 2.	
	d2 = d2 - 2.
	e1 = e1*e2
ENDFOR
result2 = EXP(-xp*xp)/SQRT(!PI)*sum
result = (ABS(arg) LT 0.2)*result1 + (ABS(arg) GE 0.2)*result2
result = result/SQRT(!PI)
RETURN, result
END ; ****** altK_Diamond ****** ;

FUNCTION D_DIAMOND, p, tau=tau_in, sigma=sigma_in
;
; Computes the function "D" from 
;
; Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
; Arguments:
;
;	p	The Laplace transform varialble
;
; Keywords:
;
;	tau	The temperature ratio, T_i/T_e.
;		Defaults to 1.0
;
;	sigma	Set to 1 to include the residual from the pole in ERFC
;		(this is about insuring that ERFC got analytically continued
;		in the proper manner). Default=0
;
; Written by W.M. Nevins
;	3/17/06

tau=1.
IF(Query_Integer(tau_in) + Query_Real(tau_in)) THEN tau=tau_in

Sigma=0
IF(Query_Integer(sigma_in) + Query_Real(sigma_in)) THEN sigma=sigma_in

one = DCOMPLEX(1.,0.)
eye = DCOMPLEX(0.,1.)

K = k_diamond(p, sigma=sigma)

psQ = p*p
firstTerm   = (one - 2.*pSq + 2.*pSq*pSq)*K + 3.*p - 2.*p*pSq
numerator   = (2.*p + (one - 2.*pSq)*K)^2
denominator = one - p*K + tau

result = firstTerm + (p/2.)*numerator/denominator
RETURN, result
END ; ****** D_Diamond ****** ;

FUNCTION D_GAM, arg
;
; The dispersion relation of the Geodesic-Acoustic
; mode are the zeros of this function. See:
;
; Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
; Argument:
;
;	p	The Laplace transform variable
;
; Written by W.M. Nevins
;	3/17/06
COMMON DD_GAM, q, Tau, sigma
p = DCOMPLEX(arg[0], arg[1])
;p = arg
D=D_Diamond(p, tau=tau, sigma=sigma)
one = DCOMPLEX(1.,0.)
D_GAM = p - 0.5*one*q^2*D
result = [DOUBLE(D_gam), IMAGINARY(D_GAM)]
;result = D_GAM
RETURN, result
END ; ****** D_GAM ****** ;


FUNCTION D_GAM_complex, arg, sigma=sigma_in, q=q_in, tau=tau_in
;
; The dispersion relation of the Geodesic-Acoustic
; mode are the zeros of this function. See:
;
; Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
; Argument:
;
;	p	The Laplace transform variable
;
; Written by W.M. Nevins
;	3/17/06
COMMON DD_GAM, q, Tau, sigma
q=1.
IF(Query_Integer(  q_in) + Query_Real(  q_in)) THEN q=q_in
tau=1.
IF(Query_integer(tau_in) + Query_Real(tau_in)) THEN tau=tau_in
Sigma=0.
IF(Query_Integer(sigma_in) + Query_Real(sigma_in)) THEN sigma=sigma_in
p = arg
D=D_Diamond(p, tau=tau, sigma=sigma)
one = DCOMPLEX(1.,0.)
D_GAM = p - 0.5*one*q^2*D
result = D_GAM
RETURN, result
END ; ****** D_GAM ****** ;

FUNCTION omega_GAM_SK, q=q_in, tau=tau_in, sigma=sigma_in, omega_1=omega_in
;
; Computes the GAM dispersion relation following:
;
; 	Sugama and Watanabe, J. Phys. Plasmas, (2006)
;
; Arguments:  None
;
; KeyWords:
;
;	q	Safety factor. 
;		Defaults to 1.0
;
;	tau	Temperature ratio, T_i/T_e.
;		Defaults to 1.0.
;
;	omega_1	Initial guess for root.  
;		Defaults to q*SQRT(7./4. + 1./tau) - 0.01i	
;
; Written by W.M. Nevins
;	3/23/2006
COMMON DD_GAM, q, Tau, sigma
q=1.
IF(Query_Integer(  q_in) + Query_Real(  q_in)) THEN q=q_in
tau=1.
IF(Query_integer(tau_in) + Query_Real(tau_in)) THEN tau=tau_in
Sigma=0.
IF(Query_Integer(sigma_in) + Query_Real(sigma_in)) THEN sigma=sigma_in
;
; Estimate value of root
;   
omega_1 = Omega_GAM(q=q, tau_e=1./tau)
arg_1 = [DOUBLE(omega_1), IMAGINARY(omega_1)]
IF(Query_Integer(omega_in) + Query_Real(omega_in)) THEN arg_1=omega_in
IF(Query_Complex(omega_in)) THEN BEGIN
	arg_1[0] = DOUBLE(omega_in)
	arg_1[1] = IMAGINARY(omega_in)
ENDIF
;arg_1 = DCOMPLEX(-0.1, omega_1)
PRINT, "omega_1 = ", arg_1
PRINT, "D_GAM(omega_1) = ", D_GAM_SK(arg_1)
result = NEWTON(arg_1, "D_GAM_SK", tolf=1.e-6, /DOUBLE)
PRINT, "NEWTON = ", result
;IF(Query_Real(result)) THEN result = COMPLEX(result, 0.)
PRINT, "D_GAM(omega_final) = ", D_GAM_SK(result)
;return, result
RETURN, COMPLEX(result[0], result[1])
END ; ****** omega_GAM ****** ;

FUNCTION p_GAM, q=q_in, tau=tau_in, sigma=sigma_in, p_1=p_in
;
; Computes the GAM dispersion relation following:
;
; 	Lebedev et al, Phys. Plasmas <3> (8), 3023 (Aug. 1996)
;
; Arguments:  None
;
; KeyWords:
;
;	q	Safety factor. 
;		Defaults to 1.0
;
;	tau	Temperature ratio, T_i/T_e.
;		Defaults to 1.0.	
;
; Written by W.M. Nevins
;	3/17/06
COMMON DD_GAM, q, Tau, sigma
q=1.
IF(Query_Integer(  q_in) + Query_Real(  q_in)) THEN q=q_in
tau=1.
IF(Query_integer(tau_in) + Query_Real(tau_in)) THEN tau=tau_in
Sigma=0.
IF(Query_Integer(sigma_in) + Query_Real(sigma_in)) THEN sigma=sigma_in
;
; Estimate value of root
;   
omega_1 = q*SQRT(7./4. + 1./tau)
;omega_1=2.38
arg_1 = [-0.1, omega_1]
IF(Query_Integer(p_in) + Query_Real(p_in)) THEN arg_1=p_in
IF(Query_Complex(p_in)) THEN BEGIN
	arg_1[0] = DOUBLE(p_in)
	arg_1[1] = IMAGINARY(p_in)
ENDIF
;arg_1 = DCOMPLEX(-0.1, omega_1)
PRINT, "p_1 = ", arg_1
PRINT, "D_GAM(p_1) = ", D_GAM(arg_1)
result = NEWTON(arg_1, "D_GAM", tolf=1.e-6, /DOUBLE)
PRINT, "NEWTON = ", result
IF(Query_Real(result)) THEN result = COMPLEX(result, 0.)
PRINT, "D_GAM(p_final) = ", D_GAM(result)
;return, result
RETURN, COMPLEX(result[0], result[1])
END ; ****** p_GAM ****** ;
