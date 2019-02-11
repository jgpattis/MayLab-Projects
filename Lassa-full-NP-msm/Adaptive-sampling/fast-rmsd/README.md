Directory setup

/topdir/         - top directory which contains all other subdirectories

/topdir/water/   - contains solvated trajectories to pull out and spawn new simulations

/topdir/trj/     - contains desolvated trajectories to help RMSD and cluster analysis run faster

/topdir/grompp/  - contains topology files so gromacs can compile pdb files into run input files (tpr)

/topdir/round-0/ - initial 6 simulations from starting structure

/topdir/round-1/ - created by fast_master.sh new round of simulations are run within spawn-$x of this directory


Scripts

fast_master.sh   - submit to debug queue of cluster
                   submits fast_rmsd.py
                   pulls out new starting structures and grompps them into tprs
                   submits new jobs

fast_rmsd.py     - performs k-centers clustering
                   performs a rmsd calculation to find structures close to monomer2_4us_ZN.pdb
                   scores clusters based on cluster counts and rmsd
                   passes score back to fast_master.sh

submit.sh        - submitted automaticly by fast_master.py
                   runs simulation of new spawns

pulloutwater.sh  - concatonates trajectories and strips out waters
                   places trajectory in /topdir/water/ and /topdir/trj/
                   keeps an index of where trajectories came from


Counters

round.txt        - keeps track of what round it is

index.txt        - how many trajectories are in /topdir/trj/

index_log.txt    - log with path to solvated trajectory and tpr for all trajectories

count.txt        - restart count for individual simulations
