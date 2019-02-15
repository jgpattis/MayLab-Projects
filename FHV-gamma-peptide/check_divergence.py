"""
Calculates 1 dimentional divergence allong the collective variable
Plots observed and consensus distribution of largest dirvergence window
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
import msmbuilder.utils
## set path
start_time = time.time()
indir = '/scratch/jap12009/gamma/mix'
topfile = '/scratch/jap12009/gamma/gamma1.pdb'

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

def calc_div(tram):
     div_c = [0] * len(tram.models)
     div_o = [0] * len(tram.models)
     div = [0] * len(tram.models)
     jsd = [0] * len(tram.models)
     for i in range(len(tram.models)):
          div_c[i] = tram.models[i].pi
          div_o[i] = tram.count_matrices_full[i].sum(axis=1)/tram.count_matrices_full[i].sum()
          div[i] = manual_kl_divergence(div_c[i], div_o[i])
          jsd[i] = js_divergence(div_c[i], div_o[i])
     return div, jsd

center = np.arange(0.7,15.5,0.1).tolist()

final = pickle.load(open('mix_tram_full_ax1_5.pickle', 'rb'))

kl, js = calc_div(final)

p = [i/16 for i in center]

plt.plot(p, kl)
plt.xlabel('% helicity', fontsize=16)
plt.ylabel('KL divergence', fontsize=16)
plt.tight_layout()
plt.savefig('mix_final_kl_divergence.eps')
plt.clf()

plt.plot(p, js)
plt.xlabel('% helicity', fontsize=16)
plt.ylabel('JS divergence', fontsize=16)
plt.tight_layout()
plt.savefig('mix_final_JS_divergence.eps')
plt.clf()

cl = pickle.load(open('mix_cl_full_ax1.pickle', 'rb'))
cl_percent = [i/16 for i in cl.clustercenters]
max_js = np.argmax(js)

print('max window of 1d js divergence is ', max_js)
print('max 1d js divergence is ', js[max_js])
print('shape of f full state is ', final.f_full_state.shape)

div_c = final.models[max_js].pi

div_o = final.count_matrices[max_js].sum(axis=1)/final.count_matrices[max_js].sum()


pickle.dump(final.active_set, open('active_set_1d_mix.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

pickle.dump(div_c, open('div_c_' + str(max_js) + '_1d_mix.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

pickle.dump(div_o, open('div_o_' + str(max_js) + '_1d_mix.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

#div_o[div_o < 0.001] = np.nan
#div_c[div_c < 0.001] = np.nan
cl_per_arr = np.array(cl_percent).reshape(451,)
sort_per = np.argsort(cl_per_arr)
print('shape of array is ', cl_per_arr.shape)
print('shape of sort is ',sort_per.shape)
print('shape of div_o is ',div_o.shape)
x = cl_per_arr[sort_per]

plt.plot(cl_per_arr[sort_per], div_o[sort_per], label='Observed distribution')
plt.plot(cl_per_arr[sort_per], div_c[sort_per], label='Consensus distribution')
plt.xlim(0.93,0.97)
plt.xlabel('% helicity', fontsize=16)
plt.ylabel('Probability', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig('obs_vs_con_1d_window_' + str(max_js) + '_mix.eps')
plt.clf()
