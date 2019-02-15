## FAST script

from msmbuilder.io import gather_metadata, save_meta, NumberedRunsParser, load_meta, itertrajs, save_generic, backup, save_trajs
from msmbuilder.cluster import KCenters

import pickle
import mdtraj as md
import numpy as np

##parameters

num_clusters = 50
alpha = 0.5
rmsd_target = '/scratch/jap12009/msm/fast/try1/monomer2_4us_ZN.pdb'
spawn = 10

f = open('round.txt')
lines = f.read()
f.close
round_num = int(lines)

## Construct and save the dataframe
parser = NumberedRunsParser(
    traj_fmt="trj-{run}.xtc",
    top_fn="/scratch/jap12009/msm/fast/try1/frame0nw_startingAPO.pdb",
    step_ps=240,
)
meta = gather_metadata("/scratch/jap12009/msm/fast/try1/trj/trj-*.xtc", parser)
save_meta(meta)

## Set up parameters for clustering
kcen = KCenters(
    n_clusters=num_clusters,
    metric='rmsd',
)

## Try to limit RAM usage
def guestimate_stride():
    total_data = meta['nframes'].sum()
    want = kcen.n_clusters * 20
    stride = max(1, total_data // want)
    print("Since we have", total_data, "frames, we're going to stride by",
          stride, "during fitting, because this is probably adequate for",
          kcen.n_clusters, "clusters")
    return stride


## Fit
kcen.fit([traj for _, traj in itertrajs(meta, stride=guestimate_stride())])
print(kcen.summarize())

## Save
save_generic(kcen, 'clusterer' + str(round_num) +'.pickl')


## Save centroids
def frame(traj_i, frame_i):
    # Note: kmedoids does 0-based, contiguous integers so we use .iloc
    row = meta.iloc[traj_i]
    return md.load_frame(row['traj_fn'], frame_i, top=row['top_fn'])


centroids = md.join((frame(ti, fi) for ti, fi in kcen.cluster_ids_),
                    check_topology=False)

centroids_fn = 'centroids_' + str(round_num) + '.xtc'
backup(centroids_fn)
centroids.save("centroids_" + str(round_num) + ".xtc")

## Count based function ## Find counts then normalize
label_list = np.concatenate(kcen.labels_)

cluster_counts = [0] * num_clusters
for i in label_list:
    cluster_counts[i] = cluster_counts[i] + 1

max_count = max(cluster_counts)
min_count = min(cluster_counts)
norm_counts = [0] * num_clusters
range_count = max_count - min_count
j = 0
for i in range(num_clusters):
    norm_counts[j] = (float(max_count) - cluster_counts[i])/ range_count
    j = j + 1

## Structure based function ## Determine function then normalize
target = md.load(rmsd_target)
rms = md.rmsd(centroids, target)

norm_rms = [0] * num_clusters
max_rms = max(rms)
min_rms = min(rms)
range_rms = max_rms - min_rms
k = 0
for i in range(num_clusters):
    norm_rms[k] = float((max_rms - rms[i])/ range_rms)
    k = k + 1

## Combine count and structure based scores then sort by score
score = np.zeros((num_clusters,2))

for i in range(num_clusters):
    score[i] = i
    score[i][1] = norm_rms[i] + alpha * norm_counts[i]

sorted_score = sorted(score, key=lambda x: x[1], reverse=True)

## Find IDs of cluster centers that need to be pulled out
spawn1 = np.zeros((spawn,2))
for i in range(spawn):
    j = int(sorted_score[i][0])
    spawn1[i] = kcen.cluster_ids_[j]

spawn2 = np.copy(spawn1)
for i in range(spawn):
    spawn2[i][1] = (spawn2[i][1] + 1)

full = np.zeros((num_clusters,7))
for i in range(num_clusters):
    j = int(sorted_score[i][0])
    full[i][0:2] = kcen.cluster_ids_[j]
    full[i][2] = sorted_score[i][1]
    full[i][3] = rms[j]
    full[i][4] = norm_rms[j]
    full[i][5] = cluster_counts[j]
    full[i][6] = norm_counts[j]

full1 = np.copy(full)
for i in range(num_clusters):
    full1[i][1] = full1[i][1] + 1

np.savetxt('output-' + str(round_num) + '.txt', spawn2, fmt='%6i %6i', delimiter=' ')
np.savetxt('full_score-' + str(round_num) + '.txt', full1, fmt='%6i %6i %6f %6f %6f %6f %6f', delimiter=' ')
