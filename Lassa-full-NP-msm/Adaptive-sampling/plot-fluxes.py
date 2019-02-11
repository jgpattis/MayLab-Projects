import matplotlib
matplotlib.use('Agg')
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from msmbuilder.io import load_trajs, load_generic
from msmbuilder.io.sampling import sample_states
from sklearn.neighbors import KDTree
import msmexplorer as msme
from msmbuilder.tpt import net_fluxes, fluxes
from msmbuilder.tpt import paths

sns.set_style('ticks')
colors = sns.color_palette()

## Load
kmeans = load_generic('../kcenters_30_100_5.pickl')
msm = load_generic('msm_kcen_30_100_5_16.pickl')
meta, ttrajs = load_trajs('../../ttrajs_a0_30')
txx = np.concatenate(list(ttrajs.values()))
a1 = ttrajs[14]

## Plot microstates
def plot_microstates(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g')
    ax.scatter(middles2[two, 0],
               middles2[two, 1],
               #s=scale * msm.populations_ + add_a_bit,
               #c=count_sum,
               #cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    #plt.colorbar(label='First Dynamical Eigenvector', ax=ax, fontsize=16)


def plot_microstates13(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    ax.scatter(middles2[two, 0],
               middles2[two, 2],
               #s=scale * msm.populations_ + add_a_bit,
               #c=count_sum,
               #cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)

ttrajs5 = {}
for k, v in ttrajs.items():
     ttrajs5[k] = v[:,:5]

clusters = [0] * 100
selected = [0] * 100
selected_txx = [0] * 100
middles0 = [0] * 100

for i in range(100):
     selected[i] = [0] * len(ttrajs5)
     for j in range(len(ttrajs5)):
          test = np.where(kmeans.labels_[j] == i)
          selected[i][j] = ttrajs5[j][test]
     selected_txx[i] = np.concatenate(selected[i])
     middles0[i] = selected_txx[i].mean(axis=0)

middles = np.asarray(middles0)
middles2 = middles[msm.state_labels_,:]

start_tic = ttrajs5[14][0]

tree = KDTree(middles2)

dist, start_ind = tree.query(start_tic, k=1)

print('distance is ', dist)
print('start index is ', start_ind[0][0])


pop_ind = msm.populations_.argmax()

print('end index is ', pop_ind)

two = [start_ind[0][0], pop_ind]

## Plot
fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates(ax)
fig.tight_layout()
fig.savefig('msm-microstates_two.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates13(ax)
fig.tight_layout()
fig.savefig('msm-microstates_two_1v3.eps')
fig.clf()

pos = dict(zip(range(len(middles2)), middles2[:,:2]))

pos13 = dict(zip(range(len(middles2)), middles2[:,0:3:2]))

ax = msme.plot_tpaths(msm, start_ind[0], [pop_ind], pos=pos)
plt.tight_layout()
plt.savefig('tpath.eps')
plt.clf()

ax = msme.plot_tpaths(msm, start_ind[0], [pop_ind], pos=pos13)
plt.tight_layout()
plt.savefig('tpath_13.eps')
plt.clf()


net_flux_matrix = net_fluxes(start_ind[0], [pop_ind], msm)

flux_matrix = fluxes(start_ind[0], [pop_ind], msm)

paths_list, fluxes1 = paths(start_ind[0], [pop_ind], net_flux_matrix)

paths_list_d, fluxes1_d = paths(start_ind[0], [pop_ind], net_flux_matrix, remove_path='bottleneck')


print('num states in path 1', len(paths_list[0]))

test_3 = np.unique(np.concatenate((paths_list[0],paths_list[1], paths_list[2]), 0))
test_5 = np.unique(np.concatenate((paths_list[0],paths_list[1], paths_list[2], paths_list[3], paths_list[4]), 0))

test_3d = np.unique(np.concatenate((paths_list_d[0],paths_list_d[1], paths_list_d[2]), 0))
test_5d = np.unique(np.concatenate((paths_list_d[0],paths_list_d[1], paths_list_d[2], paths_list_d[3], paths_list_d[4]), 0))

tot_flux = fluxes1.sum()

print('sum of fluxes is ', tot_flux)

print('flux of path 1 ', fluxes1[0] / tot_flux)

print('flux of path 2 ', fluxes1[1] / tot_flux)

print('cum flux of path 1 + 2 ', (fluxes1[:2].sum()) / tot_flux)

print('flux of path 3 ', fluxes1[2] / tot_flux)
 
print('cum flux of path 1 - 3', (fluxes1[:3].sum()) / tot_flux)

tot_fluxd = fluxes1_d.sum()

print('sum of fluxesd is ', tot_fluxd)

print('flux of pathd 1 ', fluxes1_d[0] / tot_fluxd)

print('flux of path 2 ', fluxes1_d[1] / tot_fluxd)

print('cum flux of path 1 + 2 ', (fluxes1_d[:2].sum()) / tot_fluxd)

print('flux of path 3 ', fluxes1_d[2] / tot_fluxd)
 
print('cum flux of path 1 - 3', (fluxes1_d[:3].sum()) / tot_fluxd)

plt.imshow(flux_matrix, origin='lower', interpolation='none', cmap='viridis_r')
plt.colorbar()
plt.xlabel('State', fontsize=16)
plt.ylabel('State', fontsize=16)
plt.savefig('flux_matrix_r.eps')
plt.tight_layout()
plt.clf()

plt.imshow(net_flux_matrix, origin='lower', interpolation='none', cmap='viridis_r')
plt.colorbar()
plt.xlabel('State', fontsize=16)
plt.ylabel('State', fontsize=16)
plt.tight_layout()
plt.savefig('net_flux_matrix_r.eps')
plt.clf()

states_in_top5 = np.unique(np.concatenate((paths_list[0],paths_list[1],paths_list[2],paths_list[3],paths_list[4])))
np.savetxt('states_in_top5.txt', states_in_top5, fmt='%6d')

def plot_microstate_paths(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g', label='start')
    ax.scatter(middles2[paths_list[5], 0], middles2[paths_list[5], 1], label='path6 ' + str(np.round(fluxes1[5] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[4], 0], middles2[paths_list[4], 1], label='path5 ' + str(np.round(fluxes1[4] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[3], 0], middles2[paths_list[3], 1], label='path4 ' + str(np.round(fluxes1[3] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[2], 0], middles2[paths_list[2], 1], label='path3 ' + str(np.round(fluxes1[2] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[1], 0], middles2[paths_list[1], 1], label='path2 ' + str(np.round(fluxes1[1] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[0], 0], middles2[paths_list[0], 1], label='path1 ' + str(np.round(fluxes1[0] / tot_flux, decimals=3)))
    ax.legend()
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    #plt.colorbar(label='First Dynamical Eigenvector', ax=ax, fontsize=16)


def plot_microstate_paths13(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    ax.scatter(middles2[paths_list[5], 0], middles2[paths_list[5], 2], label='path6 ' + str(np.round(fluxes1[5] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[4], 0], middles2[paths_list[4], 2], label='path5 ' + str(np.round(fluxes1[4] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[3], 0], middles2[paths_list[3], 2], label='path4 ' + str(np.round(fluxes1[3] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[2], 0], middles2[paths_list[2], 2], label='path3  ' + str(np.round(fluxes1[2] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[1], 0], middles2[paths_list[1], 2], label='path2 ' + str(np.round(fluxes1[1] / tot_flux, decimals=3)))
    ax.scatter(middles2[paths_list[0], 0], middles2[paths_list[0], 2], label='path1 ' + str(np.round(fluxes1[0] / tot_flux, decimals=3)))
    ax.legend()
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstate_paths(ax)
fig.tight_layout()
fig.savefig('msm-microstate_paths.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstate_paths13(ax)
fig.tight_layout()
fig.savefig('msm-microstate_paths_1v3.eps')
fig.clf()

def plot_microstate_pathsd(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g', label='start')
    ax.scatter(middles2[paths_list_d[5], 0], middles2[paths_list_d[5], 1], label='path5 ' + str(np.round(fluxes1_d[5] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[4], 0], middles2[paths_list_d[4], 1], label='path4 ' + str(np.round(fluxes1_d[4] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[3], 0], middles2[paths_list_d[3], 1], label='path3 ' + str(np.round(fluxes1_d[3] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[2], 0], middles2[paths_list_d[2], 1], label='path2 ' + str(np.round(fluxes1_d[2] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[1], 0], middles2[paths_list_d[1], 1], label='path1 ' + str(np.round(fluxes1_d[1] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[0], 0], middles2[paths_list_d[0], 1], label='path0 ' + str(np.round(fluxes1_d[0] / tot_fluxd, decimals=3)))
    ax.legend()
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    #plt.colorbar(label='First Dynamical Eigenvector', ax=ax, fontsize=16)


def plot_microstate_paths13d(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    ax.scatter(middles2[paths_list_d[5], 0], middles2[paths_list_d[5], 2], label='path5 ' + str(np.round(fluxes1_d[5] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[4], 0], middles2[paths_list_d[4], 2], label='path4 ' + str(np.round(fluxes1_d[4] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[3], 0], middles2[paths_list_d[3], 2], label='path3 ' + str(np.round(fluxes1_d[3] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[2], 0], middles2[paths_list_d[2], 2], label='path2 ' + str(np.round(fluxes1_d[2] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[1], 0], middles2[paths_list_d[1], 2], label='path1 ' + str(np.round(fluxes1_d[1] / tot_fluxd, decimals=3)))
    ax.scatter(middles2[paths_list_d[0], 0], middles2[paths_list_d[0], 2], label='path0 ' + str(np.round(fluxes1_d[0] / tot_fluxd, decimals=3)))
    ax.legend()
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstate_pathsd(ax)
fig.tight_layout()
fig.savefig('msm-microstate_paths_bottleneck.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstate_paths13d(ax)
fig.tight_layout()
fig.savefig('msm-microstate_paths_bottleneck_1v3.eps')
fig.clf()

