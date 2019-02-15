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
indir = '/scratch/jap12009/gamma/mix'
topfile = '/scratch/jap12009/gamma/gamma1.pdb'
save_file = 'pgTICs-ca-con-8.pkl'
## create list of trajectories and Colvar files
traj_list = [indir + "/mix_p80win{:2.1f}.xtc".format(i) for i in np.arange(0.7,15.5,0.1).tolist()]
## define topology
f = coor.featurizer(topfile)
f.add_distances_ca()
## load trajectories and colvar files
inp = coor.load(traj_list, f)
## tica

inp1 = [i[20000::10] for i in inp]


print('length of inp1 is ', len(inp1))
print('length of inp1 is ', inp1[0].shape)

for i in range(5000,59999,5000):
     k = int(i/10)
     inp2 = [j[:k,:] for j in inp1]
     print('shape of inp2 0 is ', inp2[0].shape)
     tica = coor.tica(inp2, lag=200, dim=15, kinetic_map=False)
     pickle.dump(tica, open('mix_tica_' + str(i) + '.pickle', 'wb'))
      
tica = coor.tica(inp1, lag=200, dim=15, kinetic_map=False)

pickle.dump(tica, open('mix_tica_full.pickle', 'wb'))


print('script took', (time.time() - start_time)/60, ' minutes')
