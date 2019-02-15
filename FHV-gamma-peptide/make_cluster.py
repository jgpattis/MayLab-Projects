"""
Cluster structures to pull out representative structures from
dirrerent parts of the landscape
"""
import numpy as np
import pickle
import pyemma.plots as mplt
import matplotlib
import matplotlib.pyplot as plt
import pyemma.coordinates as coor
import mdtraj as md
## set path
indir = '/scratch/shn13007/mixfolder'
tica_file = '/jason_tram/mixTICs-ca-con-8.pkl'
## outfile names
out = 'thirty_100.pkl'
ics = 3  # look at this many ICs
num_clu = 30
tica = pickle.load(open(indir + tica_file, 'rb'))
for i in tica:
        i[:,0] *= -1

tica1 = [ i[:,:ics] for i in tica]

cl = coor.cluster_mini_batch_kmeans(tica1, k=num_clu, max_iter=100)

pickle.dump(cl, open(out, 'wb'))

fig1, ax1 = mplt.plot_free_energy(np.vstack(tica1)[:,0], np.vstack(tica1)[:,1]);
plt.scatter(cl.clustercenters[:,0], cl.clustercenters[:,1], c='k')
plt.xlabel('tIC 1', fontsize=18)
plt.ylabel('tIC 2', fontsize=18)
plt.savefig('cluster_centers1v2.pdf')
plt.clf()

fig1, ax1 = mplt.plot_free_energy(np.vstack(tica1)[:,0], np.vstack(tica1)[:,2]);
plt.scatter(cl.clustercenters[:,0], cl.clustercenters[:,2], c='k') 
plt.xlabel('tIC 1', fontsize=18)
plt.ylabel('tIC 3', fontsize=18)
plt.savefig('cluster_centers1v3.pdf')
plt.clf()
a = np.arange(0.8,15.5,0.1)
b = list(range(num_clu))
ind = cl.sample_indexes_by_cluster(b, 1)
ind_traj = [ i[0][0] for i in ind]
short_list = [ a[i] for i in ind_traj]
topfile = '/scratch/shn13007/pgfolder/protein.pdb'
## create list of trajectories and Colvar files
traj_list = [indir + "/p60win{:2.1f}.xtc".format(i) for i in short_list]
f = coor.featurizer(topfile)
inp = coor.load(traj_list, f)
inp2 = [ i[20000::10] for i in inp]
frame_list = [ i[0][1] for i in ind]
#inp3 = inp2.tolist()
#frame_list2 = frame_list.tolist()
#coor.save_traj(inp2, , 'pg_cluster_centers.xtc', top=topfile)
for i, (j, k) in enumerate(zip(inp2, frame_list)):
	coor.save_traj(j, k, 'pg_cluster_centers' + str(i) + '_' + str(j) + '_' + str(k) + '.pdb', top=topfile)

	 
