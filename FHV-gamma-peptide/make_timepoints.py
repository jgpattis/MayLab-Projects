import pyemma
print(pyemma.__version__)
import pyemma.coordinates as coor
import numpy as np
import mdtraj as md
import pickle
import matplotlib.pyplot as plt
import time
## set path
start_time = time.time()
indir = '/scratch/jap12009/gamma/pg'
topfile = '/scratch/jap12009/gamma/gamma1.pdb'
save_file = 'pgTICs-ca-con-8.pkl'
## create list of trajectories and Colvar files
traj_list = [indir + "/p60win{:2.1f}.xtc".format(i) for i in np.arange(0.8,15.5,0.1)]
## define topology
f = coor.featurizer(topfile)
f.add_distances_ca()
## load trajectories and colvar files
inp = coor.load(traj_list, f)
## tica

inp1 = [i[20000::10] for i in inp]

cumvar3 = []
timescales1 = []
timescales2 = []
timescales3 = []

print('length of inp1 is ', len(inp1))
print('length of inp1 is ', inp1[0].shape)

#for i in range(5000,39999,5000):
for i in [ 2000, 4000, 6000, 8000, 10000, 12000]
     k = int(i/10)
     inp2 = [j[:k,:] for j in inp1]
     print('shape of inp2 0 is ', inp2[0].shape)
     tica = coor.tica(inp2, lag=200, dim=6, kinetic_map=False)
     cumvar3.append(tica.cumvar[2])
     timescales1.append(tica.timescales[0])
     timescales2.append(tica.timescales[1])
     timescales3.append(tica.timescales[2])
      
tica = coor.tica(inp1, lag=200, dim=6, kinetic_map=False)
cumvar3.append(tica.cumvar[2])
timescales1.append(tica.timescales[0])
timescales2.append(tica.timescales[1])
timescales3.append(tica.timescales[2])


cv = [ i * 100 for i in cumvar3]
ts1 = [ i / 100 for i in timescales1]
ts2 = [ i / 100 for i in timescales2]
ts3 = [ i / 100 for i in timescales3]

lags_ns = range(5, 41, 5)

plt.plot(lags_ns, cv)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('Cumlative Variance (%)', fontsize=20)
plt.tight_layout()
plt.savefig('pg_lag_vs_cumvar.eps')
plt.clf()

plt.plot(lags_ns, ts1)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 1 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('pg_lag_vs_time1.eps')
plt.clf()

plt.plot(lags_ns, ts2)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 2 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('pg_lag_vs_time2.eps')
plt.clf()

plt.plot(lags_ns, ts3)
plt.xlabel('lag time (ns)', fontsize=20)
plt.ylabel('timescale 3 (ns)', fontsize=20)
plt.tight_layout()
plt.savefig('pg_lag_vs_time3.eps')
plt.clf()

print('script took', (time.time() - start_time)/60, ' minutes')
