"""
plots implied timescales and percent of kenetic variance captured from tica
to help choose a tica lag time
"""

import pyemma
print(pyemma.__version__)
import pyemma.coordinates as coor
import numpy as np
import mdtraj as md
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
## set path
start_time = time.time()
indir = '/scratch/jap12009/gamma/mix'
topfile = '/scratch/jap12009/gamma/gamma1.pdb'
save_file = 'pgTICs-ca-con-8.pkl'
## create list of trajectories and Colvar files
traj_list = [indir + "/mix_p80win{:2.1f}.xtc".format(i) for i in np.arange(0.7,15.5,0.1).tolist()]
## define topology
f = coor.featurizer(topfile)
f.add_distances_ca()
## load trajectories and colvar files
inp = coor.source(traj_list, f)
## tica
lags = [10, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
cumvar3 = []
timescales1 = []
timescales2 = []
timescales3 = []

for i in lags:
     tica = coor.tica(inp, lag=i, dim=6, commute_map=True, kinetic_map=False, skip=20000, stride=10 )
     cumvar3.append(tica.cumvar[2])
     timescales1.append(tica.timescales[0])
     timescales2.append(tica.timescales[1])
     timescales3.append(tica.timescales[2])
      

lags_ns = [int(i/1000) for i in lags]

cv = [ i * 100 for i in cumvar3]
ts1 = [ i / 100 for i in timescales1]
ts2 = [ i / 100 for i in timescales2]
ts3 = [ i / 100 for i in timescales3]


plt.plot(lags_ns, cv)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('Cumlative Variance (%)', fontsize=20)
plt.tight_layout()
plt.savefig('mix_lag_vs_cumvar2.eps')
plt.clf()

plt.plot(lags_ns, ts1)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 0 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('mix_lag_vs_time1-2.eps')
plt.clf()

plt.plot(lags_ns, ts2)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 1 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('mix_lag_vs_time2-2.eps')
plt.clf()

plt.plot(lags_ns, ts3)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 2 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('mix_lag_vs_time3-2.eps')
plt.clf()

print('script took', (time.time() - start_time)/60, ' minutes')
