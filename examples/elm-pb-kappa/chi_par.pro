;Calculate chi_par and kappa_par

g = file_import('data/cbm18_8_y064_x516_090309.nc') ;
P = g.pressure;
n = 1e19;
Ti = P / (2 * n * 1.6e-19)
Te = Ti;
logLambda = 24 - alog(sqrt(n) / 1000);
AA = 2;
Mi = AA * 1.6726e-27;
Me = 9.109e-31;
taui = fltarr(g.nx, g.ny);
taui = 2.09e7 * Ti ^ 1.5 / n / logLambda * sqrt(AA);

taue = fltarr(g.nx, g.ny);
taue = 3.44e5 * Te ^ 1.5 / n / logLambda;

kappa_par_i = fltarr(g.nx, g.ny)      ;
kappa_par_i = 3.9  * Ti * 1.6e-19 * taui / Mi;

kappa_par_e = fltarr(g.nx, g.ny);
kappa_par_e = 3.2  * Te * 1.6e-19 * taue / Me;

chi_par = fltarr(g.nx, g.ny);
chi_par = (kappa_par_i + kappa_par_e) / 2;
Lbar = collect(path='data', var='Lbar');
Tbar = collect(path='data', var='Tbar');
p0 = plot(g.rxy(*,32), chi_par(*,32), xtitle="R(m)", yrange=[0,9e6], ytitle="$\chi_{||}(m^2/s)$")

vth_i = 9.79e3 * sqrt(Ti / AA)
vth_e = 4.19e5 * sqrt(Te)

q = fltarr(g.nx, g.ny)
for i = 0, g.ny-1 do q(*,i) = -g.shiftangle / (2 * !pi)

kappa_par_i_fs = vth_i * q * Lbar
kappa_par_e_fs = vth_e * q * Lbar

kappa_par_i = 1 / (1 / kappa_par_i + 1 / kappa_par_i_fs)
kappa_par_e = 1 / (1 / kappa_par_e + 1 / kappa_par_e_fs)

chi_par = (kappa_par_i + kappa_par_e) / 2
p1 = plot(g.rxy(*,32), chi_par(*,32), /overplot)

