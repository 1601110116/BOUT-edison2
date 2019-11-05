from boutdata.collect import collect
from boututils.datafile import DataFile

Ni = collect("Ni")
shape = Ni.shape

from boututils.plotdata import plotdata
plotdata(Ni[100, :, 0, :])

# from boututils.showdata import showdata
# showdata(Ni[-1, :, 0, :])
pass
import _boutcore_build