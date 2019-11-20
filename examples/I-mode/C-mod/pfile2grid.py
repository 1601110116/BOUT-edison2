"""
src_name: the .nc grid file to be modified
dst_name: the name of the output file
pfile_name is the p-file used
density: same of the 'density' in BOUT.inp
"""
from shutil import copyfile
from equilibrium_parsers import PfileParser
from boututils.file_import import file_import
import scipy.interpolate as spi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from boututils.datafile import DataFile

src_name: str = "bout.grd_C-Mod.nc"
dst_name: str = "bout.grd_C-Mod_expNi.nc"
pfile_name: str = "p1120907032.01010"
density = 1.0e19
figsize = (10, 12)
fontsize = 24
pfile_keV = True  # if the unit of Te is keV in p-file

# read profiles from p-file
pp = PfileParser(pfile_name)
psinorm_p = pp.get_psinorm()
Ni_p, Ni_unit_p = pp.get_ni()
Ne_p, Ne_unit_p = pp.get_ne()
Ti_p, Ti_unit_p = pp.get_ti()
Te_p, Te_p_unit = pp.get_te()

copyfile(src_name, dst_name)
src_file = file_import(src_name)
psinormxy = (src_file["psixy"] - src_file["psi_axis"]) \
    / (src_file["psi_bndry"] - src_file["psi_axis"])
# psi of the grid in the second branch, i.e. no pf region
psinormx = psinormxy[:, src_file["npol"][1]]

# interpolate the profiles from p-file to the grids in g-file
tmp = spi.splrep(psinorm_p, Ni_p)
Niexpx = spi.splev(psinormx, tmp)
tmp = spi.splrep(psinorm_p, Ne_p)
Neexpx = spi.splev(psinormx, tmp)
tmp = spi.splrep(psinorm_p, Ti_p)
Tiexpx = spi.splev(psinormx, tmp)
tmp = spi.splrep(psinorm_p, Te_p)
Teexpx = spi.splev(psinormx, tmp)

# reconcile the unit with .cxx file
Niexpx = Niexpx / density
Neexpx = Neexpx / density
if pfile_keV:
    Tiexpx = Tiexpx * 1e3
    Teexpx = Teexpx * 1e3

# use the interpolated profile for all y
ny = src_file["ny"]
Niexp = np.tile(Niexpx[:, np.newaxis], (1, ny))
Neexp = np.tile(Neexpx[:, np.newaxis], (1, ny))
Tiexp = np.tile(Tiexpx[:, np.newaxis], (1, ny))
Teexp = np.tile(Teexpx[:, np.newaxis], (1, ny))

# set all values in private flux region to be the same as LCFS
ixseps1 = src_file["ixseps1"]
jyseps1_1 = src_file["jyseps1_1"]
jyseps2_2 = src_file["jyseps2_2"]
Niexp[: ixseps1, :jyseps1_1+1] = Niexpx[ixseps1]
Neexp[: ixseps1, :jyseps1_1+1] = Neexpx[ixseps1]
Tiexp[: ixseps1, :jyseps1_1+1] = Tiexpx[ixseps1]
Teexp[: ixseps1, :jyseps1_1+1] = Teexpx[ixseps1]
Niexp[: ixseps1, jyseps2_2+1:] = Niexpx[ixseps1]
Neexp[: ixseps1, jyseps2_2+1:] = Neexpx[ixseps1]
Tiexp[: ixseps1, jyseps2_2+1:] = Tiexpx[ixseps1]
Teexp[: ixseps1, jyseps2_2+1:] = Teexpx[ixseps1]

# contour the results
Rxy = src_file["Rxy"]
Zxy = src_file["Zxy"]
cmap = plt.get_cmap("jet")
fig, axs = plt.subplots(2, 2, figsize=figsize, sharey='row', sharex='col')
contours = np.ndarray(shape=(2, 2), dtype=mpl.contour.QuadContourSet)
contours[0][0] = axs[0][0].contourf(Rxy, Zxy, Niexp, levels=60, cmap=cmap, antialiased=True)
contours[0][1] = axs[0][1].contourf(Rxy, Zxy, Neexp, levels=60, cmap=cmap, antialiased=True)
contours[1][0] = axs[1][0].contourf(Rxy, Zxy, Tiexp, levels=60, cmap=cmap, antialiased=True)
contours[1][1] = axs[1][1].contourf(Rxy, Zxy, Teexp, levels=60, cmap=cmap, antialiased=True)
axs[0][0].set_title(r'$n_{i}\left(%s%s\right)$' % (density, Ni_unit_p), fontsize=fontsize)
axs[0][1].set_title(r'$n_{e}\left(%s%s\right)$' % (density, Ne_unit_p), fontsize=fontsize)
axs[1][0].set_title(r'$T_{i}\left(eV\right)$', fontsize=fontsize)
axs[1][1].set_title(r'$T_{e}\left(eV\right)$', fontsize=fontsize)
xticks = mpl.ticker.MaxNLocator(nbins=5).tick_values(Rxy.min(), Rxy.max())
yticks = mpl.ticker.MaxNLocator(nbins=10).tick_values(Zxy.min(), Zxy.max())
cbars = np.ndarray(shape=(2, 2), dtype=mpl.colorbar.ColorbarBase)
for i in range(2):
    for j in range(2):
        axs[i][j].tick_params(labelsize=fontsize)
        axs[i][j].set_xticks(xticks)
        axs[i][j].set_yticks(yticks)
        cbars[i][j] = fig.colorbar(contours[i][j], ax=axs[i][j], fraction=0.08, aspect=40)
        cbars[i][j].ax.tick_params(labelsize=fontsize)
plt.savefig("nTinGrid.png")
plt.show()

# save the profiles in dst_file
with DataFile(dst_name, True) as dst_file:
    dst_file.write("Niexp", Niexp)
    dst_file.write("Neexp", Neexp)
    dst_file.write("Tiexp", Tiexp)
    dst_file.write("Teexp", Teexp)

