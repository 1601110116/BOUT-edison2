"""
all fields in CXRS file are set to 0 in the pf region in the output
src_name: the .nc grid file to be modified
dst_name: the name of the output file
CXRS_name is the p-file used
"""

from shutil import copyfile
from equilibrium_parsers import CXRSParser
from boututils.file_import import file_import
import scipy.interpolate as spi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from boututils.datafile import DataFile

src_name: str = "bout.grd_C-Mod_expNi.nc"
dst_name: str = "bout.grd_C-Mod_expNiEr.nc"
CXRS_name: str = "cxrs1120907032.01010_v20140623.txt"
figsize = (6, 6)
fontsize = 24
CXRS_Er_kVm = True
pf_0_Er = False

# read Er profile from CXRS file
# 'C' means CXRS
Cp = CXRSParser(CXRS_name)
Er_psinormC, ErC, Er_errC, Er_unitC = Cp.get_Er()

copyfile(src_name, dst_name)
src_file = file_import(src_name)
psinormxy = (src_file["psixy"] - src_file["psi_axis"]) \
    / (src_file["psi_bndry"] - src_file["psi_axis"])
# psi of the grid in the second branch, i.e. no pf region
psinormx = psinormxy[:, src_file["npol"][1]]

# interpolate the profiles from CXRS file to the grids in g-file
tmp = spi.splrep(Er_psinormC, ErC)
E_rx = spi.splev(psinormx, tmp)

# reconcile the unit with the .cxx file
if CXRS_Er_kVm:
    E_rx = E_rx * 1e3

# use the interpolated profile for all y
ny = src_file["ny"]
E_r = np.tile(E_rx[:, np.newaxis], (1, ny))

# set all values in the pf region to be 0
ixseps1 = src_file["ixseps1"]
jyseps1_1 = src_file["jyseps1_1"]
jyseps2_2 = src_file["jyseps2_2"]
if pf_0_Er:
    E_r[: ixseps1, :jyseps1_1 + 1] = 0.0
    E_r[: ixseps1, jyseps2_2 + 1:] = 0.0
else:
    E_r[: ixseps1, :jyseps1_1 + 1] = E_rx[ixseps1]
    E_r[: ixseps1, jyseps2_2 + 1:] = E_rx[ixseps1]

# contour the result
Rxy = src_file["Rxy"]
Zxy = src_file["Zxy"]
cmap = plt.get_cmap("jet")
fig = plt.figure(figsize=figsize)
ct = plt.contourf(Rxy, Zxy, E_r, levels=60, cmap=cmap, antialiased=True)
xticks = mpl.ticker.MaxNLocator(nbins=5).tick_values(Rxy.min(), Rxy.max())
yticks = mpl.ticker.MaxNLocator(nbins=10).tick_values(Zxy.min(), Zxy.max())
cbar = plt.colorbar(ct, fraction=0.08, aspect=40)
cbar.ax.tick_params(labelsize=fontsize)
plt.title(r'$E_{r}\left(V/m\right)$', fontsize=fontsize)

plt.savefig("ErinGrid.png")
plt.show()

# save the profiles in dst_file
with DataFile(dst_name, True) as dst_file:
    dst_file.write("E_r", E_r)
