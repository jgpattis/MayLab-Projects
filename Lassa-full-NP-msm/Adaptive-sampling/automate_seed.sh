#!/bin/bash
#SBATCH -J automate_seed
#SBATCH -o %j.out
#SBATCH -n 1
#SBATCH -t 00:29:00
#SBATCH -p debug
#SBATCH --exclude=cn[01-136]

## Inputs
spawns=11 # spawns variable sould be one higher then spawn wanted
topdir=/scratch/jap12009/msm/fast/msm_counts
curr=/scratch/jap12009/msm/try1/pyContacts/analysis_1_18/pyContacts/tica_a0/cluster/msm/auto     # change manually
grompp_folder=/scratch/jap12009/monomer3/starting_grompp/test
out_file=out_round27_cen_70_100_16.txt      # change manually
round=28           # change manually

echo "round is $round"

module purge
module load intelics/2013.1.039-full gromacs/5.0.1-ics

if [ ! -d $topdir/round-$round ]; then
        mkdir $topdir/round-$round
fi

number_list="one two three four five six seven eight nine ten"

counter=1
while [ $counter -lt $spawns ]; do
        trj=$( head -n $counter $curr/$out_file | tail -n 1 | awk '{print $1}' )
        frame=$( head -n $counter $curr/$out_file | tail -n 1 | awk '{print $2}' )
        counterword=$( echo $number_list | awk '{print $var}' var="$counter" )
        if [ ! -d $topdir/round-$round/$counterword-$trj-$frame ]; then
                mkdir $topdir/round-$round/$counterword-$trj-$frame
        fi
        if [ ! -d $topdir/round-$round/$counterword-$trj-$frame-v1 ]; then
                mkdir $topdir/round-$round/$counterword-$trj-$frame-v1
        fi
        if [ "$trj" -lt 2 ]; then
                let counter=counter+1
                continue
        fi
        keytpr=$( head -n $trj $grompp_folder/key.txt | tail -n 1 | awk '{print $3}' )
        keytrj=$( head -n $trj $grompp_folder/key.txt | tail -n 1 | awk '{print $2}' )
        watername=$( head -n $trj $grompp_folder/key.txt | tail -n 1 | awk '{print $4}' )
        (( time_ps= (frame + 1) * 240 ))
        srun echo 0 | gmx_mpi trjconv -s $keytrj/$keytpr -f $keytrj/$watername -b $time_ps -dump $time_ps -o $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb
        cp $keytrj/index.ndx $topdir/round-$round/$counterword-$trj-$frame
        cp $keytrj/index.ndx $topdir/round-$round/$counterword-$trj-$frame-v1
        if [ ! -f $topdir/round-$round/$counterword-$trj-$frame/topol.top ]; then
                cp $grompp_folder/topol.top $grompp_folder/topol_Protein_chain_B.itp $grompp_folder/topol_Other_chain_B2.itp $topdir/round-$round/$counterword-$trj-$frame
                na=$( grep NA $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb | wc -l )
                cl=$( grep CL $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb | wc -l )
                SOL1=$( grep SOL $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb | wc -l )
                SOL=$(( $SOL1 / 3 ))
                echo "SOL         $SOL" >> $topdir/round-$round/$counterword-$trj-$frame/topol.top
                echo "NA               $na" >> $topdir/round-$round/$counterword-$trj-$frame/topol.top
                echo "CL               $cl" >> $topdir/round-$round/$counterword-$trj-$frame/topol.top
        fi
        srun gmx_mpi grompp -f $topdir/test_diffv.mdp -c $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb -n $topdir/round-$round/$counterword-$trj-$frame/index.ndx -maxwarn 1 -p $topdir/round-$round/$counterword-$trj-$frame/topol.top -o $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.tpr
        srun gmx_mpi grompp -f $topdir/test_diffv.mdp -c $topdir/round-$round/$counterword-$trj-$frame/$counterword-$trj-$frame.pdb -n $topdir/round-$round/$counterword-$trj-$frame/index.ndx -maxwarn 1 -p $topdir/round-$round/$counterword-$trj-$frame/topol.top -o $topdir/round-$round/$counterword-$trj-$frame-v1/$counterword-$trj-$frame-v1.tpr
        if [ "$counter" -lt 3 ]; then           ## used for submitting to different partitions
             cp $topdir/submit.sh $topdir/round-$round/$counterword-$trj-$frame
             cd $topdir/round-$round/$counterword-$trj-$frame
             sbatch -J Round-$round submit.sh
             cp $topdir/submit.sh $topdir/round-$round/$counterword-$trj-$frame-v1
             cd $topdir/round-$round/$counterword-$trj-$frame-v1
             sbatch -J Round-$round submit.sh
        else
             cp $topdir/submit.sh $topdir/round-$round/$counterword-$trj-$frame
             cd $topdir/round-$round/$counterword-$trj-$frame
             sbatch -J Round-$round submit.sh
             cp $topdir/submit.sh $topdir/round-$round/$counterword-$trj-$frame-v1
             cd $topdir/round-$round/$counterword-$trj-$frame-v1
             sbatch -J Round-$round submit.sh
        fi
        cd $curr
        let counter=counter+1
done

module purge
