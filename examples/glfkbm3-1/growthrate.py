
import matplotlib.pyplot as plt
import pickle
from numpy import *
from boutdata import *

p = collect('P', path='data', info = True)
dcp = mean(p, axis = 3)
rmsp = sqrt(mean(p ** 2, axis = 3) - dcp ** 2)
pickle.dump(rmsp, open('rmsp.py.dat', 'w'))
plt.figure()
plt.contourf(rmsp[-1, :, :], 128)
plt.colorbar()
plt.figure()
plt.plot(gradient(log(rmsp[:, rmsp[-1, :500, 32].argmax(), 32])))
plt.show()

