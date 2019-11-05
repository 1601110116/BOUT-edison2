from equilibrium_parsers import PfileParser
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
parser = PfileParser('p1120907032.01010')
psinorm = parser.get_psinorm()
ne, ne_unit = parser.get_ne()
ni, ni_unit = parser.get_ni()
te, te_unit = parser.get_te()
ti, ti_unit = parser.get_ti()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('Profiles from pfile')
line_ne = ax1.plot(psinorm, np.asarray(ne), label=r'$n_{e}$')
line_ni = ax1.plot(psinorm, np.asarray(ni), label=r'$n_{i}$')
ax1.set_xlabel(r'$\psi_{nor}$')
ax1.set_ylabel(r'$n\ \left(\mathrm{m^{-3}}\right)$')
ax1.legend()

line_te = ax2.plot(psinorm, np.asarray(te), label=r'$T_{e}$')
line_ti = ax2.plot(psinorm, np.asarray(ti), label=r'$T_{i}$')
ax2.set_xlabel(r'$\psi_{nor}$')
ax2.set_ylabel(r'$T\ \left(\mathrm{keV}\right)$')
ax2.legend()
plt.savefig('exp_n_T_profiles.png')
plt.show()

from equilibrium_parsers import CXRSParser

parser = CXRSParser('cxrs1120907032.01010_v20140623.txt')
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
fig.suptitle(r'$\mathrm{CXRS}\ B^{5+}\ \mathrm{profiles}$')

psinorm, tB, tB_err, tB_unit = parser.get_tB()
axs[0][0].errorbar(psinorm, np.asarray(tB), yerr=np.asarray(tB_err), errorevery=5,
                   capsize=2.5)
axs[0][0].set_xlabel(r'$\psi_{nor}$')
axs[0][0].set_ylabel(r'$T_{B^{5+}}\ \left(\mathrm{keV}\right)$')

psinorm, nB, nB_err, nB_unit = parser.get_nB()
axs[0][1].errorbar(psinorm, np.asarray(nB), yerr=np.asarray(nB_err), errorevery=5,
                   capsize=2.5)
axs[0][1].set_xlabel(r'$\psi_{nor}$')
axs[0][1].set_ylabel(r'$n_{B^{5+}}\ \left(\mathrm{keV}\right)$')

psinorm, vpolB, vpolB_err, vpolB_unit = parser.get_vpolB()
line_vpol = axs[1][0].errorbar(psinorm, np.asarray(vpolB), yerr=np.asarray(vpolB_err), color='blue',
                               label='$v_{pol,B^{5+}}$', errorevery=5, capsize=2.5, ecolor='blue')
psinorm, vtorB, vtorB_err, vtorB_unit = parser.get_vtorB()
line_vtor = axs[1][0].errorbar(psinorm, np.asarray(vtorB), yerr=np.asarray(vtorB_err), color='orange',
                               label='$v_{tor,B^{5+}}$', errorevery=5, capsize=2.5, ecolor='orange')
axs[1][0].set_xlabel(r'$\psi_{nor}$')
axs[1][0].set_ylabel(r'$v_{B^{5+}}\ \left(\mathrm{km/s}\right)$')
axs[1][0].legend()

psinorm, Er, Er_err, Er_unit = parser.get_Er()
axs[1][1].errorbar(psinorm, np.asarray(Er), yerr=np.asarray(Er_err), errorevery=5,
                   capsize=2.5)
axs[1][1].set_xlabel(r'$\psi_{nor}$')
axs[1][1].set_ylabel(r'$E_{r}\ \left(\mathrm{kV/m}\right)$')
plt.savefig('CXRS_profiles.png')
plt.show()
