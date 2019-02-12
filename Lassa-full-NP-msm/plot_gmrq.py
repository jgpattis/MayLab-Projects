import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

lags = [ 5, 15, 30, 45, 60, 75]
components = [ 3, 4, 5, 6, 7, 15]
clusters = [ 25, 50, 75, 100, 125, 150, 200, 300]

r = pd.read_pickle('gmrq_try10.pickl')


all_list = [0] * len(components)
all_list2 = [0] * len(components)
median_list = [0] * len(components)
std_list = [0] * len(components)
med_test = [0] * len(components)
for i in range(len(components)):
    all_list[i] = r[r['components'] == components[i]]
    all_list2[i] = [0] * len(clusters)
    median_list[i] = [0] * len(clusters)
    std_list[i] = [0] * len(clusters)
    med_test[i] = [0] * len(clusters)
    for j in range(len(clusters)):
        all_list2[i][j] = all_list[i][all_list[i]['n_clusters'] == clusters[j]]
        median_list[i][j] = all_list2[i][j].groupby('tica_lag').aggregate(np.nanmedian).drop('fold', axis=1)
        std_list[i][j] = all_list2[i][j].groupby('tica_lag').aggregate(np.nanstd).drop('fold', axis=1)
        med_test[i][j] = [0] * len(lags)
        for k in range(len(lags)):
            med_test[i][j][k] = median_list[i][j]['test_score'][lags[k]]

med_array = np.array(med_test)
best = np.where(med_array == med_array.max())

fig, axs = plt.subplots(len(components),len(clusters), figsize=(12,9))

lags2 = [i * 0.24 for i in lags]
for j in range(len(components)):
    for k in range(len(clusters)):
        axs[j][k].errorbar(lags2, median_list[j][k]['train_score'], yerr= std_list[j][k]['train_score'], c='b')
        axs[j][k].errorbar(lags2, median_list[j][k]['test_score'], yerr= std_list[j][k]['test_score'], c='r')
        #axs[j][k].set_ylim(3,4.1)       #plot first then find a good range

axs[best[0][0]][best[1][0]].scatter(lags2[best[2][0]], m1[best], c='g', s=100, marker='*')       
 
for n in range(len(components)):
    axs[n][0].set_ylabel('GMRQ-4', fontsize=12)
            
for l in range(len(components)):
    axs[l][len(clusters)-1].yaxis.set_label_position('right')
    axs[l][len(clusters)-1].set_ylabel(str(components[l]) + ' Comp', fontsize=12)

for m in range(len(clusters)):
    axs[0][m].set_title(str(clusters[m]) + ' Cl', fontsize=12)

for p in range(len(clusters)):
    axs[len(components)-1][p].set_xlabel('tICA lag (ns)', fontsize=12)
    
fig.tight_layout()

fig.savefig('gmrq_plot_try10.eps')
