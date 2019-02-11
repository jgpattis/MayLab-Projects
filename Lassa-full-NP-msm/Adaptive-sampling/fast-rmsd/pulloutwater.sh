#!/bin/bash
#SBATCH -o %j.out
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -p general

module load intelics/2013.1.039-full gromacs/5.0.1-ics

topdir=`cat path.txt`
parent=/scratch/jap12009/msm/fast/try1
outname=monomer
gromppfolder=/scratch/jap12009/msm/fast/grompp
round=`cat $parent/round.txt`
input=`cat tpr.txt`
frame=`cat frame.txt`
trjfolder=/scratch/jap12009/msm/fast/try1/trj
water=/scratch/jap12009/msm/fast/try1/water

name1=trj_water.xtc #full
name2=trj_nw.xtc # no water
name3=trj_skip240.xtc #skip 10

srun gmx_mpi trjcat -f $topdir/part*/*.xtc -o $topdir/$name1

echo 24 | srun --mpi=pmi2 gmx_mpi trjconv -s $topdir/part1/run.tpr -f $topdir/$name1 -n $gromppfolder/index.ndx -o $topdir/$name2 -pbc whole

echo 24 | srun --mpi=pmi2 gmx_mpi trjconv -s $topdir/part1/run.tpr -f $topdir/$name2 -n $gromppfolder/index.ndx -o $topdir/$name3 -pbc nojump -skip 240

index=`cat $parent/index.txt`

cp $topdir/trj_water.xtc $water/trj-$index-water.xtc
cp $topdir/trj_skip240.xtc $trjfolder/trj-$index.xtc
echo "$index   $topdir   $input   $frame" >> $parent/index_log.txt
let index=index+1
echo $index > $parent/index.txt

exit
