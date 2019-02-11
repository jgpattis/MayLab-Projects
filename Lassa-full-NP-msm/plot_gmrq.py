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
        for k in range(len(lags)):
            med_test[i][j] = [0] * len(lags)
            temp = median_list[i][j]['test_score']
            med_test[i][j][k] = temp[lags[k]]

