"""
compares the difference in the pmf between earlier timepoints and the final pmf
for both wham and tram
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pickle
import pyemma.coordinates as coor
import pyemma.thermo as thermo
import time
import matplotlib.pyplot as plt
from functools import reduce
from scipy.optimize import minimize_scalar
from sklearn.metrics import mean_squared_error as sklearn_mse

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
     #tram_f = pickle.load(open(full, 'rb'))
     #tram_4 = pickle.load(open(short, 'rb'))
     #both = np.intersect1d(tram_4.active_set, tram_f.active_set)
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

cl_f = pickle.load(open('mix_cl_full_ax1.pickle', 'rb'))
#dif, final, four, both = find_diff('pg_tram_full_ax1.pickle', 'pg_tram_40_ax1.pickle')

tram_10 = pickle.load(open('mix_tram_10ns_ax1_e5.pickle', 'rb'))
tram_12 = pickle.load(open('mix_tram_12-5ns_ax1_e5.pickle', 'rb'))
tram_15 = pickle.load(open('mix_tram_15ns_ax1_e5.pickle', 'rb'))
tram_17 = pickle.load(open('mix_tram_17-5ns_ax1_e5.pickle', 'rb'))
tram_20 = pickle.load(open('mix_tram_20ns_ax1_e5.pickle', 'rb'))
tram_22 = pickle.load(open('mix_tram_22-5ns_ax1_e5.pickle', 'rb'))
tram_25 = pickle.load(open('mix_tram_25ns_ax1_e5.pickle', 'rb'))
tram_27 = pickle.load(open('mix_tram_27-5ns_ax1_e5.pickle', 'rb'))
tram_30 = pickle.load(open('mix_tram_30ns_ax1_e5.pickle', 'rb'))
tram_32 = pickle.load(open('mix_tram_32-5ns_ax1_e5.pickle', 'rb'))
tram_35 = pickle.load(open('mix_tram_35ns_ax1_e5.pickle', 'rb'))
tram_37 = pickle.load(open('mix_tram_37-5ns_ax1_e5.pickle', 'rb'))
tram_40 = pickle.load(open('mix_tram_40ns_ax1_e5.pickle', 'rb'))
tram_42 = pickle.load(open('mix_tram_42-5ns_ax1_e5.pickle', 'rb'))
tram_45 = pickle.load(open('mix_tram_45ns_ax1_e5.pickle', 'rb'))
tram_47 = pickle.load(open('mix_tram_47-5ns_ax1_e5.pickle', 'rb'))
tram_50 = pickle.load(open('mix_tram_50ns_ax1_e5.pickle', 'rb'))
tram_52 = pickle.load(open('mix_tram_52-5ns_ax1_e5.pickle', 'rb'))
tram_55 = pickle.load(open('mix_tram_55ns_ax1_e5.pickle', 'rb'))
tram_57 = pickle.load(open('mix_tram_57-5ns_ax1_e5.pickle', 'rb'))
tram_full = pickle.load(open('mix_tram_full_ax1_5.pickle', 'rb'))

wham_10 = pickle.load(open('mix_wham_10ns_ax1_e4.pickle', 'rb'))
wham_12 = pickle.load(open('mix_wham_12-5ns_ax1_e4.pickle', 'rb'))
wham_15 = pickle.load(open('mix_wham_15ns_ax1_e4.pickle', 'rb'))
wham_17 = pickle.load(open('mix_wham_17-5ns_ax1_e4.pickle', 'rb'))
wham_20 = pickle.load(open('mix_wham_20ns_ax1_e4.pickle', 'rb'))
wham_22 = pickle.load(open('mix_wham_22-5ns_ax1_e4.pickle', 'rb'))
wham_25 = pickle.load(open('mix_wham_25ns_ax1_e4.pickle', 'rb'))
wham_27 = pickle.load(open('mix_wham_27-5ns_ax1_e4.pickle', 'rb'))
wham_30 = pickle.load(open('mix_wham_30ns_ax1_e4.pickle', 'rb'))
wham_32 = pickle.load(open('mix_wham_32-5ns_ax1_e4.pickle', 'rb'))
wham_35 = pickle.load(open('mix_wham_35ns_ax1_e4.pickle', 'rb'))
wham_37 = pickle.load(open('mix_wham_37-5ns_ax1_e4.pickle', 'rb'))
wham_40 = pickle.load(open('mix_wham_40ns_ax1_e4.pickle', 'rb'))
wham_42 = pickle.load(open('mix_wham_42-5ns_ax1_e4.pickle', 'rb'))
wham_45 = pickle.load(open('mix_wham_45ns_ax1_e4.pickle', 'rb'))
wham_47 = pickle.load(open('mix_wham_47-5ns_ax1_e4.pickle', 'rb'))
wham_50 = pickle.load(open('mix_wham_50ns_ax1_e4.pickle', 'rb'))
wham_52 = pickle.load(open('mix_wham_52-5ns_ax1_e4.pickle', 'rb'))
wham_55 = pickle.load(open('mix_wham_55ns_ax1_e4.pickle', 'rb'))
wham_57 = pickle.load(open('mix_wham_57-5ns_ax1_e4.pickle', 'rb'))
wham_full = pickle.load(open('mix_wham_full_ax1_e4.pickle', 'rb'))

both = reduce(np.intersect1d, (tram_10.active_set, tram_12.active_set, tram_15.active_set, tram_17.active_set, tram_20.active_set, tram_22.active_set, tram_25.active_set, tram_27.active_set, tram_30.active_set, tram_32.active_set, tram_35.active_set, tram_37.active_set, tram_40.active_set, tram_42.active_set, tram_45.active_set, tram_47.active_set, tram_50.active_set, tram_52.active_set, tram_55.active_set, tram_57.active_set, tram_full.active_set, wham_10.active_set, wham_12.active_set, wham_15.active_set, wham_17.active_set, wham_20.active_set, wham_22.active_set, wham_25.active_set, wham_27.active_set, wham_30.active_set, wham_32.active_set, wham_35.active_set, wham_37.active_set, wham_40.active_set, wham_42.active_set, wham_45.active_set, wham_47.active_set, wham_50.active_set, wham_52.active_set, wham_55.active_set, wham_57.active_set, wham_full.active_set))


tram_fb = tram_full.f_full_state[both]
tram_10b = tram_10.f_full_state[both]
tram_12b = tram_12.f_full_state[both]
tram_15b = tram_15.f_full_state[both]
tram_17b = tram_17.f_full_state[both]
tram_20b = tram_20.f_full_state[both]
tram_22b = tram_22.f_full_state[both]
tram_25b = tram_25.f_full_state[both]
tram_27b = tram_27.f_full_state[both]
tram_30b = tram_30.f_full_state[both]
tram_32b = tram_32.f_full_state[both]
tram_35b = tram_35.f_full_state[both]
tram_37b = tram_37.f_full_state[both]
tram_40b = tram_40.f_full_state[both]
tram_42b = tram_42.f_full_state[both]
tram_45b = tram_45.f_full_state[both]
tram_47b = tram_47.f_full_state[both]
tram_50b = tram_50.f_full_state[both]
tram_52b = tram_52.f_full_state[both]
tram_55b = tram_55.f_full_state[both]
tram_57b = tram_57.f_full_state[both]

min_pos = np.argmin(tram_fb)
tram_f_al = tram_fb - (tram_fb[min_pos])
tram_10_al = tram_10b - (tram_10b[min_pos])
tram_12_al = tram_12b - (tram_12b[min_pos])
tram_15_al = tram_15b - (tram_15b[min_pos])
tram_17_al = tram_17b - (tram_17b[min_pos])
tram_20_al = tram_20b - (tram_20b[min_pos])
tram_22_al = tram_22b - (tram_22b[min_pos])
tram_25_al = tram_25b - (tram_25b[min_pos])
tram_27_al = tram_27b - (tram_27b[min_pos])
tram_30_al = tram_30b - (tram_30b[min_pos])
tram_32_al = tram_32b - (tram_32b[min_pos])
tram_35_al = tram_35b - (tram_35b[min_pos])
tram_37_al = tram_37b - (tram_37b[min_pos])
tram_40_al = tram_40b - (tram_40b[min_pos])
tram_42_al = tram_42b - (tram_42b[min_pos])
tram_45_al = tram_45b - (tram_45b[min_pos])
tram_47_al = tram_47b - (tram_47b[min_pos])
tram_50_al = tram_50b - (tram_50b[min_pos])
tram_52_al = tram_52b - (tram_52b[min_pos])
tram_55_al = tram_55b - (tram_55b[min_pos])
tram_57_al = tram_57b - (tram_57b[min_pos])

wham_fb = wham_full.f_full_state[both]
wham_10b = wham_10.f_full_state[both]
wham_12b = wham_12.f_full_state[both]
wham_15b = wham_15.f_full_state[both]
wham_17b = wham_17.f_full_state[both]
wham_20b = wham_20.f_full_state[both]
wham_22b = wham_22.f_full_state[both]
wham_25b = wham_25.f_full_state[both]
wham_27b = wham_27.f_full_state[both]
wham_30b = wham_30.f_full_state[both]
wham_32b = wham_32.f_full_state[both]
wham_35b = wham_35.f_full_state[both]
wham_37b = wham_37.f_full_state[both]
wham_40b = wham_40.f_full_state[both]
wham_42b = wham_42.f_full_state[both]
wham_45b = wham_45.f_full_state[both]
wham_47b = wham_47.f_full_state[both]
wham_50b = wham_50.f_full_state[both]
wham_52b = wham_52.f_full_state[both]
wham_55b = wham_55.f_full_state[both]
wham_57b = wham_57.f_full_state[both]

wham_f_al = wham_fb - (wham_fb[min_pos])
wham_10_al = wham_10b - (wham_10b[min_pos])
wham_12_al = wham_12b - (wham_12b[min_pos])
wham_15_al = wham_15b - (wham_15b[min_pos])
wham_17_al = wham_17b - (wham_17b[min_pos])
wham_20_al = wham_20b - (wham_20b[min_pos])
wham_22_al = wham_22b - (wham_22b[min_pos])
wham_25_al = wham_25b - (wham_25b[min_pos])
wham_27_al = wham_27b - (wham_27b[min_pos])
wham_30_al = wham_30b - (wham_30b[min_pos])
wham_32_al = wham_32b - (wham_32b[min_pos])
wham_35_al = wham_35b - (wham_35b[min_pos])
wham_37_al = wham_37b - (wham_37b[min_pos])
wham_40_al = wham_40b - (wham_40b[min_pos])
wham_42_al = wham_42b - (wham_42b[min_pos])
wham_45_al = wham_45b - (wham_45b[min_pos])
wham_47_al = wham_47b - (wham_47b[min_pos])
wham_50_al = wham_50b - (wham_50b[min_pos])
wham_52_al = wham_52b - (wham_52b[min_pos])
wham_55_al = wham_55b - (wham_55b[min_pos])
wham_57_al = wham_57b - (wham_57b[min_pos])

num_list = [10, 12, 15, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50, 52, 55, 57]
#tram_list = [tram_10_al, tram_12_al, tram_15_al, tram_17_al, tram_20_al, tram_22_al, tram_25_al, tram_27_al, tram_30_al, tram_32_al, tram_35_al, tram_37_al, tram_40_al, tram_42_al, tram_]

#tram_list_al = [ 'tram_' + str(i) + '_al' for i in num_list]
#wham_list_al = [ 'wham_' + str(i) + '_al' for i in num_list]

tram_f_k = tram_f_al * 0.6
tram_10_k = tram_10_al * 0.6
tram_12_k = tram_12_al * 0.6
tram_15_k = tram_15_al * 0.6
tram_17_k = tram_17_al * 0.6
tram_20_k = tram_20_al * 0.6
tram_22_k = tram_22_al * 0.6
tram_25_k = tram_25_al * 0.6
tram_27_k = tram_27_al * 0.6
tram_30_k = tram_30_al * 0.6
tram_32_k = tram_32_al * 0.6
tram_35_k = tram_35_al * 0.6
tram_37_k = tram_37_al * 0.6
tram_40_k = tram_40_al * 0.6
tram_42_k = tram_42_al * 0.6
tram_45_k = tram_45_al * 0.6
tram_47_k = tram_47_al * 0.6
tram_50_k = tram_50_al * 0.6
tram_52_k = tram_52_al * 0.6
tram_55_k = tram_55_al * 0.6
tram_57_k = tram_57_al * 0.6

wham_f_k = wham_f_al * 0.6
wham_10_k = wham_10_al * 0.6
wham_12_k = wham_12_al * 0.6
wham_15_k = wham_15_al * 0.6
wham_17_k = wham_17_al * 0.6
wham_20_k = wham_20_al * 0.6
wham_22_k = wham_22_al * 0.6
wham_25_k = wham_25_al * 0.6
wham_27_k = wham_27_al * 0.6
wham_30_k = wham_30_al * 0.6
wham_32_k = wham_32_al * 0.6
wham_35_k = wham_35_al * 0.6
wham_37_k = wham_37_al * 0.6
wham_40_k = wham_40_al * 0.6
wham_42_k = wham_42_al * 0.6
wham_45_k = wham_45_al * 0.6
wham_47_k = wham_47_al * 0.6
wham_50_k = wham_50_al * 0.6
wham_52_k = wham_52_al * 0.6
wham_55_k = wham_55_al * 0.6
wham_57_k = wham_57_al * 0.6

#tram_list_k = [ 'tram_' + str(i) + '_k' for i in num_list]
#wham_list_k = [ 'wham_' + str(i) + '_k' for i in num_list]

tram_list = [tram_10_k, tram_12_k, tram_15_k, tram_17_k, tram_20_k, tram_22_k, tram_25_k, tram_27_k, tram_30_k, tram_32_k, tram_35_k, tram_37_k, tram_40_k, tram_42_k, tram_45_k, tram_47_k, tram_50_k, tram_52_k, tram_55_k, tram_57_k]

wham_list = [wham_10_k, wham_12_k, wham_15_k, wham_17_k, wham_20_k, wham_22_k, wham_25_k, wham_27_k, wham_30_k, wham_32_k, wham_35_k, wham_37_k, wham_40_k, wham_42_k, wham_45_k, wham_47_k, wham_50_k, wham_52_k, wham_55_k, wham_57_k]

r_list_tram = [0] * len(tram_list)
r_list_wham = [0] * len(wham_list)
diff_list_tram = [0] * len(tram_list)
diff_list_wham = [0] * len(tram_list)

for k,i in enumerate(tram_list):
    def difference(x):
        new = i + x
        return np.asscalar(sklearn_mse(tram_f_k, new))

    res = minimize_scalar(difference)

    r_list_tram[k] = i + res.x

for k,i in enumerate(wham_list):
    def difference(x):
        new = i + x
        return np.asscalar(sklearn_mse(wham_f_k, new))

    res = minimize_scalar(difference)

    r_list_wham[k] = i + res.x

for k,i in enumerate(r_list_tram):
    diff_list_tram[k] = np.sum(np.absolute(np.subtract(tram_f_k, i)))

for k,i in enumerate(r_list_wham):
    diff_list_wham[k] = np.sum(np.absolute(np.subtract(wham_f_k, i)))


nums = [10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5]

centers = cl_f.clustercenters[both]

p = centers/16

fig1, ax1 = plt.subplots()
#ax1.plot(nums, [t10, t12, t15, t17, t20, t22, t25, t27, t30, t32, t35, t37, t40, t42, t45, t47, t50, t52, t55, t57], label='dTRAM')
#ax1.plot(nums, [w10, w12, w15, w17, w20, w22, w25, w27, w30, w32, w35, w37, w40, w42, w45, w47, w50, w52, w55, w57], label='WHAM')
ax1.plot(nums, diff_list_tram, label='dTRAM')
ax1.plot(nums, diff_list_wham, label='WHAM')
ax1.set_xlabel('Trajectory length (ns)', fontsize=16)
ax1.set_ylabel('Difference in PMF (kcal/mol)', fontsize=16)
ax1.legend()
fig1.tight_layout()
fig1.savefig('difference_tram_wham_mix_rmse.eps')

#fig2, ax2 = plt.subplots()
#ax2.plot( p, wham_10_k, label='10 ns')
#ax2.plot( p, wham_15_k, label='15 ns')
#ax2.plot( p, wham_20_k, label='20 ns')
#ax2.plot( p, wham_25_k, label='25 ns')
#ax2.plot( p, wham_30_k, label='30 ns')
#ax2.plot( p, wham_35_k, label='35 ns')
#ax2.plot( p, wham_f_k, label='final 40 ns')
#ax2.set_xlabel('Percent Helicity', fontsize=16)
#ax2.set_ylabel('PMF (kcal/mol)', fontsize=16)
#ax2.legend()
#ax2.set_title('WHAM', fontsize=20)
#fig2.tight_layout()
#fig2.savefig('wham_pmf_comparison.eps')

#fig3, ax3 = plt.subplots()
#ax3.plot( p, tram_10_k, label='10 ns')
#ax3.plot( p, tram_15_k, label='15 ns')
#ax3.plot( p, tram_20_k, label='20 ns')
#ax3.plot( p, tram_25_k, label='25 ns')
#ax3.plot( p, tram_30_k, label='30 ns')
#ax3.plot( p, tram_35_k, label='35 ns')
#ax3.plot( p, tram_f_k, label='final 40 ns')
#ax3.set_xlabel('Percent Helicity', fontsize=16)
#ax3.set_ylabel('PMF (kcal/mol)', fontsize=16)
#ax3.legend()
#ax3.set_title('TRAM', fontsize=20)
#fig3.tight_layout()
#fig3.savefig('tram_pmf_comparison.eps')

#plt.plot(centers, final, label='final')
#plt.plot(centers, four, label='40 ns')
#plt.legend()
#plt.xlabel('S alpha', fontsize=16)
#plt.ylabel('pmf (kj/mol)', fontsize=16)
#plt.tight_layout()
#plt.savefig('pmf_final_vs_40.eps')


print('script took', (time.time() - start_time)/60, ' minutes')

