"""
grd_file_name is the .nc grid file to be modified
p_file_name is the p-file used to modify grd_file
"""

src_name: str = "bout.grd_C-Mod.nc"
dst_name: str = "bout.grd_C-Mod_expNi.nc"
pfile_name: str = "p1120907032.01010"

from shutil import copyfile
copyfile(src_name, dst_name)
from boututils.datafile import DataFile
from equilibrium_parsers import PfileParser
pp = PfileParser(pfile_name)
with open(pfile_name, 'r'):
    psi_p = pp.get_psinorm()
    Ni_p, Ni_p_unit = pp.get_ni()
    Te_p, Te_p_unit = pp.get_te()

# with DataFile(src_name) as src_file:
#     pass
#
f = DataFile("bout.grd_C-Mod.nc")
# var = f.impl.handle.variables["Ti0"][:]
# # b = var[:]
# Bt = f.read("Btxy")
# Bp = f.read("Bpxy")
# R = f.read("Rxy")
# htht = f.read("hthe")
# q = Bt * htht / (Bp)
# pass
