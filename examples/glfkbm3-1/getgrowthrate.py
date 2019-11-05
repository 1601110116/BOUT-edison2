
import matplotlib.pyplot as plt
import pickle
from numpy import *
from boutdata import *

rmsp = pickle.load(open('rmsp.py.dat', 'r'))
plt.figure()
plt.contourf(rmsp[-1, :, :], 128)
plt.colorbar()
plt.figure()
plt.plot(gradient(log(rmsp[:, rmsp[-1, :500, 32].argmax(), 32])))
plt.show()

