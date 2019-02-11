These scripts assume that trajectories have already:
    -been featurized
    -run through tica
    -clustered using k-centers clustering
using the msmbuilder turorial protocal

plot-kcen-counts.py   will determine areas of tica space that have not been well sampled
    spawning new simulations here will help discover new states

plot-fluxes.py   will calculate the top maximum flux pathways

microstate_bayes_pyemma_and_plot.py   will use a baysian msm to determine the uncertainty in stationary population along the top 5 maximum flux pathways

automate_seed.sh  will spawn new simulations from the points selected in plot-kcen-counts.py (low counts) or microstate_bayes_pyemma_and_plot.py (high uncerainty)
