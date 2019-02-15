"""
run tica on umbrella sampling data
"""

import pyemma
print(pyemma.__version__)
import pyemma.coordinates as coor
import numpy as np
import mdtraj as md
import pickle
## set path
indir = '/scratch/shn13007/mixfolder'
topfile = '/scratch/shn13007/mixfolder/protein.pdb'
save_file = 'mixTICs-ca-con-8.pkl'
## create list of trajectories and Colvar files
traj_list = [indir + "/mix_p80win{:2.1f}.xtc".format(i) for i in np.arange(0.7,15.5,0.1)]
## define topology
f = coor.featurizer(topfile)
f.add_distances_ca()
## load trajectories and colvar files
inp = coor.source(traj_list, f)
## tica
tica = coor.tica(inp, lag=8000, dim=6, kinetic_map=True, skip=20000, stride=10 )
print('cumlative variance of ca con test nw with tau 4 ns is ', tica.cumvar)
print('Timescals are ', tica.timescales)
## make IC 1 start as a negative number(on the left)
for i in range(2):
        if tica.eigenvectors[0, i] > 0:
                tica.eigenvectors[:, i] *= -1
## get IC's froms tica object
Y = tica.get_output(stride=10,skip=20000)
## Save to a file to use later
pickle.dump(Y, open(save_file,'wb'))
