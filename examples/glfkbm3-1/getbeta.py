
import matplotlib.pyplot as plt
import pickle
from numpy import *
from boutdata import *
from boututils import *

g = file_import('data/cyclone_516x64.nc')
n0 = g['NiHat'] * 1.e19
t0 = g['TeHat'] * 2000 * 1.6e-19
b0 = g['Bxy']
pe = n0 * t0
mu0 = 4.e-7 * pi

beta = pe / (b0 * b0 / (2 * mu0))
print beta[240, 32]
plt.plot(beta[:, 32])
plt.xlabel(r'x', fontsize = 18)
plt.ylabel(r'$\beta$', fontsize = 18)
plt.title(r'$\beta$ when $n_0=10^{19}$', fontsize = 18)
plt.show()

