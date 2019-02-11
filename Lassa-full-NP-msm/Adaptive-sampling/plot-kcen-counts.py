import matplotlib
matplotlib.use('Agg')
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from msmbuilder.io import load_trajs, load_generic
from msmbuilder.io.sampling import sample_states

sns.set_style('ticks')
colors = sns.color_palette()

## Load
kmeans = load_generic('../kcenters_30_100_5.pickl')
msm = load_generic('msm_kcen_30_100_5_16.pickl')
meta, ttrajs = load_trajs('../../ttrajs_a0_30')
txx = np.concatenate(list(ttrajs.values()))
a1 = ttrajs[1]

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
    im = ax.scatter(middles[msm.state_labels_, 0],
               middles[msm.state_labels_, 1],
               s=scale * msm.populations_ + add_a_bit,
               c=count_sum,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='Count sum', fontsize=16)

def plot_microstates13(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    im = ax.scatter(middles[msm.state_labels_, 0],
               middles[msm.state_labels_, 2],
               s=scale * msm.populations_ + add_a_bit,
               c=count_sum,
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='Count sum', fontsize=16)

def plot_microstates_new(ax):
    ax.hexbin(txx[:, 0], txx[:, 1],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,1], marker="x", s=200, c='g')
    im = ax.scatter(middles[score[:10], 0],
               middles[score[:10], 1],
               c=count_sum[score[:10]],
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 2", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='Count sum', fontsize=16)

def plot_microstates_new13(ax):
    ax.hexbin(txx[:, 0], txx[:, 2],
              cmap='Greys',
              mincnt=1,
              bins='log',
              )

    scale = 100 / np.max(msm.populations_)
    add_a_bit = 25
    ax.scatter(a1[0,0],a1[0,2], marker="x", s=200, c='g')
    im = ax.scatter(middles[score[:10], 0],
               middles[score[:10], 2],
               c=count_sum[score[:10]],
               cmap='RdBu'
               )
    ax.set_xlabel("tIC 1", fontsize=16)
    ax.set_ylabel("tIC 3", fontsize=16)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(label='Count sum', fontsize=16)

count = msm.countsmat_

plt.imshow(count, origin='lower', interpolation='none', cmap='viridis_r')
plt.colorbar()
plt.xlabel('State', fontsize=16)
plt.ylabel('State', fontsize=16)
plt.savefig('counts_cen_30_100_5_16_r.eps')
plt.clf()

plt.imshow(count, origin='lower', interpolation='none', cmap='viridis')
plt.colorbar()
plt.xlabel('State', fontsize=16)
plt.ylabel('State', fontsize=16)
plt.savefig('counts_cen_30_100_5_16.eps')
plt.clf()

print('msm.mapping ', msm.mapping_)

count_sum = np.sum(count, axis=0)

score = np.argsort(count_sum)

print('score is ', score)

print('min is ', count_sum[score[0]])

s2 = score.tolist()

ttrajs5 = {}
for k, v in ttrajs.items():
     ttrajs5[k] = v[:,:5]

clusters = [0] * 100
selected = [0] * 100
selected_txx = [0] * 100
middles0 = [0] * 100
cluster_counts = [0] * 100

for i in range(100):
     selected[i] = [0] * len(ttrajs5)
     for j in range(len(ttrajs5)):
          test = np.where(kmeans.labels_[j] == i)
          selected[i][j] = ttrajs5[j][test]
     selected_txx[i] = np.concatenate(selected[i])
     middles0[i] = selected_txx[i].mean(axis=0)
     cluster_counts[i] = len(selected_txx[i])

middles = np.asarray(middles0)

out = sample_states(ttrajs5, middles[score[:10]], k=5)

print('out is ', out)

out2 = sample_states(ttrajs5, middles[score[:10]], k=1)

np.savetxt('out_round9_cen_70_100_16.txt', out2, fmt='%6d')

state = range(100)

plt.plot(state, count_sum)
plt.xlabel('State', fontsize=16)
plt.ylabel('Count Sum', fontsize=16)
plt.savefig('counts_cen_sum_30_100_5_16.eps')
plt.clf()

count_sum_sort = [0] * 100
count_sum_sort = count_sum[score]

plt.plot(state, count_sum_sort)
plt.xlabel('State', fontsize=16)
plt.ylabel('Count Sum', fontsize=16)
plt.savefig('counts_cen_sum_30_100_5_16_sorted.eps')
plt.clf()

cluster_counts2 = np.array(cluster_counts)
cluster_count_sum = [0] * 100
cluster_count_sum = cluster_counts2[score]

plt.plot(state, cluster_count_sum)
plt.xlabel('State', fontsize=16)
plt.ylabel('Count Sum', fontsize=16)
plt.savefig('counts_cen_sum_from_cluster_30_100_5_16_sorted.eps')
plt.clf()

## Plot
fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates(ax)
fig.tight_layout()
fig.savefig('msm-microstates_cen_30_100_5_16_sum_counts_middle.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates13(ax)
fig.tight_layout()
fig.savefig('msm-microstates_cen_30_100_5_16_sum_counts_1v3_middle.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates_new(ax)
fig.tight_layout()
fig.savefig('msm-microstates_cen_30_100_5_16_new_middle.eps')
fig.clf()

fig, ax = plt.subplots(figsize=(7, 5))
plot_microstates_new13(ax)
fig.tight_layout()
fig.savefig('msm-microstates_cen_30_100_5_16_new_1v3_middle.eps')
fig.clf()

dia = count.diagonal()

plt.plot(state, dia)
plt.xlabel('State', fontsize=16)
plt.ylabel('Count Sum', fontsize=16)
plt.savefig('counts_diagnal_cen_30_100_5_16.eps')
plt.clf()

count -= count.diagonal() * np.eye(*count.shape)

count_sum2 = np.sum(count, axis=0)

plt.plot(state, count_sum2)
plt.xlabel('State', fontsize=16)
plt.ylabel('Count Sum not diagnal', fontsize=16)
plt.savefig('counts_sum_not_diagnal_cen_30_100_5_16.eps')
plt.clf()

