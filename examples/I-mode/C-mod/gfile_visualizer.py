from boututils.file_import import file_import
import numpy as np
import matplotlib.pyplot as plt
g = file_import('bout.grd_C-Mod.nc')
psi = (g['psixy']-g['psi_axis'])/(g['psi_bndry']-g['psi_axis'])
x = psi[:, 38]
y = np.arange(64)
X, Y = np.meshgrid(y, x)
Jpar0 = g['Jpar0']
plt.contourf(X, Y, Jpar0)
plt.colorbar()
plt.show()