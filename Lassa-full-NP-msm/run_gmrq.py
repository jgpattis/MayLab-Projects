from msmbuilder.io import load_trajs
from msmbuilder.decomposition import tICA
from sklearn.pipeline import Pipeline
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import MiniBatchKMeans
from sklearn.cross_validation import ShuffleSplit
from msmbuilder.msm import MarkovStateModel
import pandas as pd

meta, ftrajs = load_trajs('ftrajs')

ftrajs_list = list(ftrajs.values())

components = [ 3, 4, 5, 6, 7]
lags = [ 5, 15, 30, 45, 60, 75]
cl = [ 25, 50, 75, 100, 125, 150, 200, 300]

model = Pipeline([
    ('tica', tICA(kinetic_mapping=True)), 
    ('cluster', MiniBatchKMeans()), 
    ('msm', MarkovStateModel(lag_time=16, n_timescales=3, verbose=False))
    #### use if percent retained is lower than 80 percent
    ##('msm', MarkovStateModel(lag_time=16, n_timescales=3, ergodic_cuttoff=off, prior_conts=0.00001, verbose=False))
])

shuffle_split = ShuffleSplit(len(ftrajs_list), n_iter=15, test_size=0.5)
print(len(ttrajs))
results = []

for lag in lags:
    for n_clusters in cl:
        for comp in components:
            model.set_params(cluster__n_clusters=n_clusters, tica__lag_time=lag, tica__n_components=comp)
                for fold, (train, validate) in enumerate(shuffle_split):
                    try:
                        train_ttrajs = [ ftrajs_list[i] for i in train ]
                        val_ttrajs = [ ftrajs_list[i] for i in validate ]

                        model.fit(train_ttrajs)

                        train_score = model.score(train_ttrajs)
                        test_score = model.score(val_ttrajs)

                        msm_timescales = model.named_steps['msm'].timescales_
                        pr = model.named_steps['msm'].percent_retained_
                        print('I did something ', fold)
                        print('clusters lag components train test', n_clusters, lag, comp, train_score, test_score)
                        results.append({
                            'train_score': train_score,
                            'test_score': test_score,
                            'n_clusters': n_clusters,
                            'components': comp,
                            'msm_timescale0': msm_timescales[0],
                            'msm_timescale1': msm_timescales[1],
                            'percent_retained' : pr,
                            'tica_lag': lag,
                            'fold': fold
                         })

                     except:
                         print('something went wrong! -Jason')
                         pass

results1 = pd.DataFrame(results)

## Save
print(results1.head())
results1.to_pickle('gmrq_try10.pickl')
