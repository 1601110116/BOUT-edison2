; Import the experimentally measured radial electric field, Epsi, into 
; a grid file, where Epsi is the flux function given as
; Epsi = -dPhi/dpsi = Er/RBp. Measured Epsi are often fitted with a B-spline,
; so this code makes a call external to a fortran B-spline routine. 
; 
; Inputs:
; (1) tt - a vector containing the "knots" of the B-spline fit
; (2) coef - a vector containing the B-spline fitting coefficients
; (3) file_path - a string describing the path of the grid file to be modified

pro Epsi_2_grid, tt, coef, file_path

; Set up a psi array where we evaluate the fit. Need to be careful about going 
; beyond psi=1.0.  Some fits are not done past psi=1.0.  You cannot tell that 
; from the fit coefficients.

; import the grid file you wish to import Epsi to
g=file_import(file_path)

tmp = size(g.psixy, /dimensions)
pol_size = tmp[1] ; # of poloidal grid points
pol_ind = pol_size / 2 ; index of midpoint in poloidal grid
delvar, tmp

; flux coordinate, psi (in physical units: Wb, or T * m^2)
psi = g.psixy[*, pol_ind]

; normalized psi
psi_n = (psi - g.PSI_AXIS) / (g.PSI_BNDRY - g.PSI_AXIS)

; To evaluate the spline, we will need to make a call external. For this to
; work, we expect to have some routines in a directory that we can find.
; Need to set environmental variable to point to this.

; do NOT appear to need(???)
;setenv,'IDLUTILS_DIR=./idlutils'

; Determine the electric field by evaluating the B-Spline
Epsi = slatec_bvalu(psi_n, tt, coef)

; However, [Epsi] = kiloRadians/s, so multiply by 1000 to get to Radians/s
Epsi = 1000d * Epsi

; BOUT++'s elm-pb module will only accept a 2-D array for Epsi. However, as
; Epsi is a flux function, it is uniform in the poloidal direction.
Epsi_2d = Epsi # replicate(1, pol_size) 

; Now, import the Epsi to the grid 
handle = file_open(file_path, /write)
s = file_write(handle, 'Epsi', Epsi_2d) ; s = 0 => no error, s = 1 => error
file_close, handle

; Calculate the diamagnetic contribution to Epsi
; Epsi_diamag = (1 / 2 * ni0 * Z * e) * dP/dpsi, where it is assumed Pi = Pe = 0.5*P
ni0 = g.ni0[*, pol_ind] ; in units of 1e20 m^-3
e = 1.602e1 ; electron charge normalized to 1e-20 (because ni0 normalized to 1e20)
Z = 1d ; ion charge
Epsi_diamag = (0.5 / (ni0 * Z * e)) * deriv(psi, g.pressure[*, pol_ind]) ; note: derivative wrt
									; psi, not psi_n 
; Plot total Epsi and Diamagnetic contribution to Epsi
y_min = min([min(Epsi), min(Epsi_diamag)])
y_max = max([max(Epsi), max(Epsi_diamag)])
character_size = 1.5

plot, psi_n, Epsi, $
	title = 'Total Electric Field, E!D!7w!3!N, and the Diamagnetic Contribution, E!D!7w,!3Dia!3!N', $
	xtitle = 'Normalized Flux, !7w!3!Dn!N', $
	subtitle = 'Solid = Total, Dashed = Diamagnetic Contribution', $
	ytitle = 'Electric Field (1/s)', $
	yrange = [y_min, y_max], $
	charsize = character_size
oplot, psi_n, Epsi_diamag, linestyle = 2

end
