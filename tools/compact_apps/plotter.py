import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as text
from matplotlib.lines import Line2D
import matplotlib.font_manager
import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size=14)

fig = plt.figure(figsize=(6,4))
plot = fig.add_subplot(111)

fig.subplots_adjust(left=0.16,right=0.94)
fig.subplots_adjust(bottom=0.12,top=0.92)

data = np.loadtxt("output/data.512")

index = data[:,0]
value = data[:,1]

#------------------------------------------------------------
plot.set_xlabel('Index')
plot.set_ylabel('Value')
plt.title('TDMA')
plot.grid(True)

plot.plot(cores,gflop,'-o',color='darkgreen', lw=2.5)

a=np.array([1,2,3,4])
tick_locs = a
tick_lbls = ['1','2','3','4']
plt.xticks(tick_locs, tick_lbls)

plot.set_xlim([0.5,4.5])
plot.set_ylim([0,80])
   
plt.savefig('TDMA.pdf')
plt.savefig('TDMA.eps')
#------------------------------------------------------------





