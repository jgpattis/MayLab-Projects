from msmbuilder.io import load_trajs, save_trajs, save_generic, load_generic
from pyemma.msm import bayesian_markov_model
import time
import matplotlib
matplotlib.use('Agg')
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
sns.set_style('ticks')
colors = sns.color_palette()
import pickle
start_time = time.time()
## Load
meta, ktrajs = load_trajs('../../ktrajs_cen_30_100_5')


dtrajs = list(ktrajs.values())
type(dtrajs[0])
## Fit
msm = bayesian_markov_model(dtrajs, lag=16, nsamples=100000)

print('done with bmm')

## Load
kmeans = load_generic('../../kcenters_30_100_5.pickl')
msm2 = load_generic('../msm_kcen_30_100_5_16.pickl')
meta, ttrajs = load_trajs('../../../ttrajs_a0_30')
txx = np.concatenate(list(ttrajs.values()))
a1 = ttrajs[14]

print('done with load')

print('active_count_fraction is ', msm.active_count_fraction)
print('active_state_fraction is ', msm.active_state_fraction)
print('shape of transition matrix is ', msm.P.shape)
print('mean of eigenvalues is ', msm.sample_mean('eigenvalues', 7))
print('std of eigenvalues is ', msm.sample_std('eigenvalues', 7))
print('mean of pi is ', msm.sample_mean('pi'))
print('std of pi is ', msm.sample_std('pi'))
print('shape of right eigenvectors is ', msm.sample_std('eigenvectors_right', 2).shape)
print('shape of sample f of eigenvalues is ', len(msm.sample_f('eigenvalues', 1)))


## Plot microstates
def plot_microstates(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    pop = msm.sample_std('pi')
    eig = msm.sample_std('eigenvectors_right', 2)[:,1]
    print('pop is ', pop)
    print('eig is ', eig)
    scale = 100 / np.max(pop)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g')
    im = ax.scatter(kmeans.cluster_centers_[msm2.state_labels_, 0],
               kmeans.cluster_centers_[msm2.state_labels_, 1],
               s=scale * pop + add_a_bit,
               c=eig,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='First Eigenvector STD', fontsize=16)

def plot_microstates13(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    pop = msm.sample_std('pi')
    eig = msm.sample_std('eigenvectors_right', 2)[:,1]
    scale = 100 / np.max(pop)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    im = ax.scatter(kmeans.cluster_centers_[msm2.state_labels_, 0],
               kmeans.cluster_centers_[msm2.state_labels_, 2],
               s=scale * pop + add_a_bit,
               c=eig,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='First Eigenvector STD', fontsize=16)

def plot_microstates_pop(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    pop_std = msm.sample_std('pi')
    pop = msm2.populations_
    scale = 100 / np.max(pop)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g')
    im = ax.scatter(kmeans.cluster_centers_[msm2.state_labels_, 0],
               kmeans.cluster_centers_[msm2.state_labels_, 1],
               s=scale * pop + add_a_bit,
               c=pop_std,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='population STD', fontsize=16)

def plot_microstates13_pop(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    pop_std = msm.sample_std('pi')
    pop = msm2.populations_
    scale = 100 / np.max(pop)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    im = ax.scatter(kmeans.cluster_centers_[msm2.state_labels_, 0],
               kmeans.cluster_centers_[msm2.state_labels_, 2],
               s=scale * pop + add_a_bit,
               c=pop_std,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='population STD', fontsize=16)

## Plot
fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates_pop(ax)
fig.tight_layout()
fig.savefig('msm-microstates_30_100_16_baysian_pop_pyemma.eps')
fig.clf()

fig1, ax1 = plt.subplots(figsize=(7, 5))
plot_microstates13_pop(ax1)
fig1.tight_layout()
fig1.savefig('msm-microstates_30_100_16_baysian_pop13_pyemma.eps')
fig1.clf()

eig = msm.sample_f('eigenvalues', 4)
eig1 = [i[1] for i in eig]
eig2 = [i[2] for i in eig]
eig3 = [i[3] for i in eig]

n1, bins1, patches1 = plt.hist(eig1, bins=50, histtype='step', label='eig 2', color='b')
n2, bins2, patches2 = plt.hist(eig2, bins=50, histtype='step', label='eig 3', color='r')
n3, bins3, patches3 = plt.hist(eig3, bins=50, histtype='step', label='eig 4', color='g')
line1 = plt.scatter(msm2.eigenvalues_[1], 0, marker='x', s=200, c='b', label='eig 2 max like')
line2 = plt.scatter(msm2.eigenvalues_[2], 0, marker='x', s=200, c='r', label='eig 3 max like')
line3 = plt.scatter(msm2.eigenvalues_[3], 0, marker='x', s=200, c='g', label='eig 4 max like')
plt.xlabel('Eigenvalue', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig('eigenval_hist_pyemma.eps')
plt.clf()

pop = msm.sample_f('pi')

avg_pi = msm.sample_mean('pi')

plt.hist(avg_pi, bins=50, histtype='step', label='distribution of pi', color='b')
plt.xlabel('Probability', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig('average_population_distribution_pyemma.eps')
plt.clf()

sorted_pi = np.argsort(avg_pi)

for i in sorted_pi[:4]:
     k = [j[i] for j in pop]
     plt.hist(k, bins=50, histtype='step', label='state ' + str(i), color='b')
     line1 = plt.scatter(msm2.populations_[i], 0, marker='x', s=200, c='r', label='Max like')
     plt.xlabel('Probability', fontsize=16)
     plt.ylabel('Frequency', fontsize=16)
     plt.legend()
     plt.tight_layout()
     plt.savefig('population_hist_pyemma_lowest_pi_state_' + str(i) + '.eps')
     plt.clf()

for i in sorted_pi[-4:]:
     k = [j[i] for j in pop]
     plt.hist(k, bins=50, histtype='step', label='state ' + str(i), color='b')
     line1 = plt.scatter(msm2.populations_[i], 0, marker='x', s=200, c='r', label='Max like')
     plt.xlabel('Probability', fontsize=16)
     plt.ylabel('Frequency', fontsize=16)
     plt.legend()
     plt.tight_layout()
     plt.savefig('population_hist_pyemma_highest_pi_state_' + str(i) + '.eps')
     plt.clf()

std_pi = msm.sample_std('pi')

sorted_std = np.argsort(std_pi)

for i  in sorted_std[:4]:
     k = [j[i] for j in pop]
     plt.hist(k, bins=50, histtype='step', label='state ' + str(i), color='b')
     line1 = plt.scatter(msm2.populations_[i], 0, marker='x', s=200, c='r', label='Max like')
     plt.xlabel('Probability', fontsize=16)
     plt.ylabel('Frequency', fontsize=16)
     plt.legend()
     plt.tight_layout()
     plt.savefig('population_hist_pyemma_lowest_pi_std_state_' + str(i) + '.eps')
     plt.clf()

for i in sorted_std[-4:]:
     k = [j[i] for j in pop]
     plt.hist(k, bins=50, histtype='step', label='state ' + str(i), color='b')
     line1 = plt.scatter(msm2.populations_[i], 0, marker='x', s=200, c='r', label='Max like')
     plt.xlabel('Probability', fontsize=16)
     plt.ylabel('Frequency', fontsize=16)
     plt.legend()
     plt.tight_layout()
     plt.savefig('population_hist_pyemma_highest_pi_std_state_' + str(i) + '.eps')
     plt.clf()

print('sorted_std is ', sorted_std)
np.savetxt('sorted_std.txt', sorted_std, fmt='%6d')

path_list = np.loadtxt('states_in_top5.txt', dtype='int')
limit_to_path = sorted_std[path_list]
print('sorted_std limit to path is ', limit_to_path)
np.savetxt('on_path_sorted_std.txt', limit_to_path, fmt='%6d')

print('script took', (time.time() - start_time)/60, ' minutes')

msm.save('pyemma_30_100_16_baysian.h5')
