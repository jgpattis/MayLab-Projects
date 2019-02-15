"""
Calculates 3 dimentional divergence allong the collective variable and independent component 2 and 3
Plots observed and consensus distribution of largest dirvergence window
"""

import numpy as np
import pickle
import pyemma.coordinates as coor
import pyemma.thermo as thermo
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class import_colvar:
     def __init__(self, system):
          if system == 'pg':
               self.indir = '/scratch/jap12009/gamma/pg' 
               self.center = np.arange(0.8,15.5,0.1).tolist()
               self.colvar_string = "/new60nsCOLVAR{:2.1f}"
          elif system == 'pc':
               self.indir = '/scratch/jap12009/gamma/pc'
               self.center = [ 0.7, 0.8, 0.9, 1.0, 1.3, 1.5, 1.6,  1.7, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 3.0, 3.3, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5,  10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.7, 13.0, 13.3, 13.5, 13.6, 13.7, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.7, 14.9, 15.0, 15.1, 15.3]
               self.colvar_string = "/new80nsCOLVAR{:2.1f}"
          elif system == 'mix':
               self.indir = '/scratch/jap12009/gamma/mix'
               self.center = np.arange(0.7,15.5,0.1).tolist()
               self.colvar_string = "/comboCOLVAR{:2.1f}"
          else:
               print('You made a typo on the system name')
          self.colvar_list = [self.indir + self.colvar_string.format(i) for i in self.center]
          self.col = [np.loadtxt(f, skiprows=1) for f in self.colvar_list]
          self.outname = system

def trim(end, data):
     if end == 'full':
          temp = [ i[20000::10,1] for i in data]
          temp2 = []
          for i in temp:
               temp2.append(np.expand_dims(i, axis=1))
          return [ i.copy(order='C') for i in temp2]
     else:
          temp = [ i[20000:end:10,1] for i in data]
          temp2 = []
          for i in temp:
               temp2.append(np.expand_dims(i, axis=1))
          return [ i.copy(order='C') for i in temp2]

class cluster:
     def __init__(self, ic, data=None, pickle_file=None):
          if data is None:
               self.cluster = pickle.load(open(pickle_file, 'rb'))
          else:
               self.cluster = coor.cluster_regspace(data, max_centers=1000, dmin=0.025)
               pickle.dump(self.cluster, open(ic.outname + '_cl_full.pickle', 'wb'))

     def cluster_trimmed(self, data1):
          self.cluster.dtrajs = [ self.cluster.assign(i) for i in data1 ]

def run_tram(trimmed_data, clust, ic):
     force = 119.503
     force_list = [500] * len(trimmed_data)
     KT = 2.496
     lag=200
     return thermo.estimate_umbrella_sampling(trimmed_data, clust.cluster.dtrajs, ic.center, force_list, kT=KT, maxiter=100000, maxerr=1.0E-4,  dt_traj='10 ps', lag=lag, save_convergence_info=200, estimator='dtram')

def find_diff(full, short):
     tram_f = pickle.load(open(full, 'rb'))
     tram_4 = pickle.load(open(short, 'rb'))
     both = np.intersect1d(tram_4.active_set, tram_f.active_set)
     tram_fb = tram_f.f_full_state[both]
     tram_4b = tram_4.f_full_state[both]
     min_pos = np.argmin(tram_fb)
     tram_f_al = tram_fb - (tram_fb[min_pos])
     tram_4_al = tram_4b - (tram_4b[min_pos])
     dif_sum_norm = (np.sum(np.absolute(np.subtract(tram_f_al, tram_4_al)))) / len(tram_f_al)
     return dif_sum_norm, tram_f_al, tram_4_al, both

## plot time vs dif_sum_norm

## plot aligned pmf

start_time= time.time()

indir = '/scratch/jap12009/gamma/mix'
topfile = '/scratch/jap12009/gamma/gamma1.pdb'
save_file = 'pgTICs-ca-con-8.pkl'
## create list of trajectories and Colvar files
#traj_list = [indir + "/p60win{:2.1f}.xtc".format(i) for i in np.arange(0.8,15.5,0.1)]
## define topology
#f = coor.featurizer(topfile)
#f.add_distances_ca()
## load trajectories and colvar files
#inp = coor.source(traj_list, f)

#tica = coor.tica(inp, lag=500, dim=3, commute_map=True, kinetic_map=False, skip=20000, stride=10 )

col = import_colvar('mix')

col_skip = [ i[20000::10, 1] for i in col.col]

#print('shape of col_skip is ', len(col_skip[0]))
print('min and max of col_skip ', np.min(col_skip), np.max(col_skip))

clust_col_skip_obj = coor.cluster_regspace(col_skip, max_centers=1000, dmin=0.025)

clust_col_skip_dtraj = clust_col_skip_obj.dtrajs
#clust_col_skip_dtraj =  pickle.load(open('clust_col_skip_dtraj_cl_full.pickle', 'rb'))
pickle.dump(clust_col_skip_obj, open('clust_col_skip_obj_cl_full.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(clust_col_skip_dtraj, open('clust_col_skip_dtraj_cl_full_try2.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

print('length of clust_col_skip_dtraj is ', len(clust_col_skip_dtraj))
print('length of clust_col_skip_dtraj[0] is ', len(clust_col_skip_dtraj[0]))

#Y = tica.get_output()
tica =  pickle.load(open('mix_tica_full.pickle', 'rb'))
Y = tica.get_output()
print('shape of tica is ', len(Y[0]))
#cl_f = pickle.load(open('pg_cl_full_ax1.pickle', 'rb'))

Y2 = [ i[:,1:3] for i in Y]

cluster_tic = coor.cluster_kmeans(Y2, k=10, max_iter=100)

clust_out = cluster_tic.dtrajs
#clust_out =  pickle.load(open('tic_cl_full_dtraj.pickle', 'rb'))
pickle.dump(clust_out, open('tic_cl_full_dtraj3_try2.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

pickle.dump(cluster_tic.clustercenters, open('tic_cl_full_centers3_try2.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

D = [0] * len(clust_col_skip_dtraj)

for i in range(len(clust_col_skip_dtraj)):
     D[i] = [0] * len(clust_col_skip_dtraj[0])
     for j in range(len(clust_col_skip_dtraj[0])):
          temp = clust_col_skip_dtraj[i][j] * 10 + clust_out[i][j]
          D[i][j] = np.int(temp)

print('shape of D is ', len(D))
print('shape of D0 is ', len(D[0]))

temp2 = []
for i in col_skip:
     temp2.append(np.expand_dims(i, axis=1))
new_traj = [ i.copy(order='C') for i in temp2]

D2 = [ np.asarray(i, dtype=int, order='C') for i in D ]
force = 119.503
force_list = [500] * len(new_traj)
KT = 2.496
lag=200
center = np.arange(0.7,15.5,0.1).tolist()
print('center is ', center)
tram = thermo.estimate_umbrella_sampling(new_traj, D2, center, force_list, kT=KT, maxiter=500000, maxerr=5.0E-4,  dt_traj='10 ps', lag=lag, save_convergence_info=200, estimator='dtram')

print('tram f shape is ', len(tram.f))

pickle.dump(tram.f, open('tram_f_mix_3dim_5e4_3_try2.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(tram, open('tram_mix_3dim_e4_3_2try2.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

def manual_kl_divergence(P, Q, scalar=True):
    if len(P.shape) == 1:
        P = np.array([P])
        Q = np.array([Q])
    vec = []
    for row in range(P.shape[0]):
        temp = 0
        for i, entry in enumerate(P[row]):
            if entry * Q[row][i] != 0:  # i.e. one or both is not zero
                temp += entry * np.log(entry / Q[row][i])
        vec.append(temp)
    result = np.array(vec)
    # sometimes there are very tiny negative values close to zero
    result = [np.max([i, 0]) for i in result]
    if scalar:
        return np.sum(result)
    else:
        return result

def js_divergence(P, Q, scalar=True):
    M = np.mean([P, Q], axis=0)
    return (0.5 * manual_kl_divergence(P, M, scalar=scalar) +
            0.5 * manual_kl_divergence(Q, M, scalar=scalar))

#center = np.arange(0.8,15.5,0.1).tolist()

#final = pickle.load(open('pg_tram_full_ax1_5.pickle', 'rb'))

def calc_diva(tram):
     div_c = [0] * len(tram.models)
     div_o = [0] * len(tram.models)
     div = [0] * len(tram.models)
     jsd = [0] * len(tram.models)
     for i in range(len(tram.models)):
          div_c[i] = tram.models[i].pi
          div_o[i] = tram.count_matrices[i].sum(axis=1)/tram.count_matrices[i].sum()
          div[i] = manual_kl_divergence(div_c[i], div_o[i])
          jsd[i] = js_divergence(div_c[i], div_o[i])
     return div, jsd

#center = np.arange(0.8,15.5,0.1).tolist()

#final = pickle.load(open('pg_tram_full_ax1_5.pickle', 'rb'))

kla, jsa = calc_diva(tram)

pickle.dump(kla, open('kl_mix_multi_5e4_active3.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

pickle.dump(jsa, open('js_mix_multi_5e4_active3.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

#print('length of kl is ', len(kl))

p = [i/16 for i in center]

plt.plot(p, kla)
plt.xlabel('% helicity', fontsize=16)
plt.ylabel('KL divergence', fontsize=16)
plt.tight_layout()
plt.savefig('mix_final_multi_kl_divergence3_try2.eps')
plt.clf()

plt.plot(p, jsa)
plt.xlabel('% helicity', fontsize=16)
plt.ylabel('JS divergence', fontsize=16)
plt.tight_layout()
plt.savefig('mix_final_multi_JS_divergence3_try2.eps')
plt.clf()

#dif, final, four, both = find_diff('pg_tram_full_ax1.pickle', 'pg_tram_40_ax1.pickle')

#print('diff is ', dif)

#centers = cl_f.clustercenters[both]

#plt.plot(centers, final, label='final')
#plt.plot(centers, four, label='40 ns')
#plt.legend()
#plt.xlabel('S alpha', fontsize=16)
#plt.ylabel('pmf (kj/mol)', fontsize=16)
#plt.tight_layout()
#plt.savefig('pmf_final_vs_40.eps')


print('script took', (time.time() - start_time)/60, ' minutes')

