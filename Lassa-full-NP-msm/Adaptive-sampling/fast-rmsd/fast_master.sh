#!/bin/bash
#SBATCH -J fast_master
#SBATCH -o %j.out
#SBATCH -n 1
#SBATCH -t 00:29:00
#SBATCH -p debug
#SBATCH --exclude=cn[01-136]

## Inputs
spawns=11 # spawns variable sould be one higher then spawn wanted
topdir=/scratch/jap12009/msm/fast/try1
trjdir=/scratch/jap12009/msm/fast/try1/water
grompp_folder=/scratch/jap12009/msm/fast/grompp
tpr=/scratch/jap12009/msm/fast/try1/round-0/spawn-1/start_fast_1.tpr

if [ ! -f $topdir/round.txt ]; then 
	echo 1 > $topdir/round.txt
fi

round=`cat $topdir/round.txt | awk '{print $1}'`

echo "round is $round"

module load tcl/8.6.6.8606 sqlite/3.18.0 python/3.6.1

if [ -f $topdir/output-$round.txt ]; then
	mv $topdir/output-$round.txt $topdir/backup-output-$round.txt
fi

if [ -f $topdir/full_score-$round.txt ]; then
        mv $topdir/full_score-$round.txt $topdir/backup-full_score-$round.txt
fi

python $topdir/fast_rmsd.py

while [ ! -f $topdir/output-$round.txt ]; do
	sleep 30
done

module load intelics/2013.1.039-full gromacs/5.0.1-ics

if [ ! -d $topdir/round-$round ]; then
        mkdir $topdir/round-$round
fi

counter3=1
while [ $counter3 -lt $spawns ]; do
	mkdir $topdir/round-$round/spawn-$counter3
        mkdir $topdir/round-$round/spawn-$counter3/v1
	mkdir $topdir/round-$round/spawn-$counter3/v2
	let counter3=counter3+1
done

counter=1
while [ $counter -lt $spawns ]; do
	trj=$( head -n $counter $topdir/output-$round.txt | tail -n 1 | awk '{print $1}' )
	frame=$( head -n $counter $topdir/output-$round.txt | tail -n 1 | awk '{print $2}' )
	(( time_ps= frame * 240 ))
	echo 0 | srun -n 1 gmx_mpi trjconv -s $tpr -f $trjdir/trj-$trj-water.xtc -dump $time_ps -o $topdir/round-$round/spawn-$counter/spawn-$round-$counter.pdb
	echo "$topdir/round-$round/spawn-$counter/v1/" > $topdir/round-$round/spawn-$counter/v1/path.txt
	echo "$topdir/round-$round/spawn-$counter/spawn-$round-$counter.pdb" > $topdir/round-$round/spawn-$counter/v1/tpr.txt
	if [ -f $topdir/round-$round/spawn-$counter/v1/frame.txt ]; then
		mv $topdir/round-$round/spawn-$counter/v1/frame.txt $topdir/round-$round/spawn-$counter/v1/backup_frame.txt
	fi
        echo "$trj   $frame   $time_ps" > $topdir/round-$round/spawn-$counter/v1/frame.txt
	cp $topdir/submit.sh $topdir/round-$round/spawn-$counter/v1/
	cd $topdir/round-$round/spawn-$counter/v1/
        sbatch -J ro-$round-sp-$counter-v1 $topdir/round-$round/spawn-$counter/v1/submit.sh
	echo "$topdir/round-$round/spawn-$counter/v2/" > $topdir/round-$round/spawn-$counter/v2/path.txt
	echo "$topdir/round-$round/spawn-$counter/spawn-$round-$counter.pdb" > $topdir/round-$round/spawn-$counter/v2/tpr.txt
	if [ -f $topdir/round-$round/spawn-$counter/v2/frame.txt ]; then
		mv $topdir/round-$round/spawn-$counter/v2/frame.txt $topdir/round-$round/spawn-$counter/v2/backup_frame.txt
	fi
	echo "$trj   $frame   $time_ps" > $topdir/round-$round/spawn-$counter/v2/frame.txt
        cp $topdir/submit.sh $topdir/round-$round/spawn-$counter/v2/
	cd $topdir/round-$round/spawn-$counter/v2/
        sbatch -J ro-$round-sp-$counter-v2 $topdir/round-$round/spawn-$counter/v2/submit.sh
	let counter=counter+1
done

let round=round+1
echo $round > $topdir/round.txt

exit
