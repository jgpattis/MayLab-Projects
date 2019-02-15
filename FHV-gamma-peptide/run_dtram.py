import numpy as np
import pickle
import pyemma.coordinates as coor
import pyemma.thermo as thermo
## set path
indir = '/scratch/shn13007/mixfolder'
## outfile names
clust_out = 'mix_dtram_clust_50_kj_200.pkl'
out = 'mix_dtram_50_kj_200.pkl'
## load colvars 
center = np.arange(0.7,15.5,0.1)
center2 = center.tolist()
colvar_list = [indir + "/comboCOLVAR{:2.1f}".format(i) for i in center]
col = [np.loadtxt(f, skiprows=1) for f in colvar_list]
length = len(center)
print(length)
force = 119.503
force_list = [ 500 ] * length
## Start doing stuff
max_centers = 1500
dmin = 0.015
kt = 0.596
#cv = list(col[20000:,1])
cv = [ i[20000::10,1] for i in col]
cv2 = [ i.copy(order='C') for i in cv]

us_cluster = coor.cluster_regspace( cv2, max_centers= max_centers, dmin= dmin)
w = thermo.estimate_umbrella_sampling( cv2, us_cluster.dtrajs, center2, force_list, kT=2.496, maxiter=50000, lag=200, dt_traj='10 ps', save_convergence_info=200, estimator='dtram')

pickle.dump(us_cluster, open(clust_out, 'wb'))
pickle.dump(w, open(out,'wb'))

