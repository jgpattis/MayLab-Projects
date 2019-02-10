#!/bin/bash
#SBATCH -o %j.out
#SBATCH -N 2
#SBATCH -n 48
#SBATCH --ntasks-per-node=24
#SBATCH -t 11:15:00
#SBATCH -p general
#SBATCH --exclude=cn[01-136]

module purge
module load intelics/2013.1.039-full gromacs/5.0.1-ics

topdir=`cat path.txt`
script=submit.sh
input=`cat tpr.txt`
parent=/scratch/jap12009/msm/fast/try1
maxiter=6
outname=monomer
gromppfolder=/scratch/jap12009/msm/fast/grompp 
round=`cat $parent/round.txt`

echo "round is $round"
echo "input is $input"
 
if [ ! -f $topdir/count.txt ]; then
        echo 1 > $topdir/count.txt
fi

count=`cat $topdir/count.txt | awk '{print $1}'`

     mkdir -p $topdir/part$count
     currout=$topdir/part$count
     let prev=$count-1
     prevout=$topdir/part$prev
     let next=$count+1
    
     if [ $count -eq "1" ]; then
           echo "started part $count"
           cd $currout
           srun -n 1 --mpi=pmi2 gmx_mpi grompp -f $gromppfolder/test_diffv.mdp -c $input -p $gromppfolder/topol.top -maxwarn 3 -n $gromppfolder/index.ndx -o $currout/run.tpr
           srun --mpi=pmi2 gmx_mpi mdrun -s $currout/run.tpr -deffnm $currout/$outname$count -cpt 100 -cpo $currout/part$count.cpt -c $currout/$outname$count.pdb -maxh 11
     else
           echo "restart part $count"
           cd $currout
	   srun --mpi=pmi2 gmx_mpi mdrun -s $topdir/part1/run.tpr -deffnm $currout/$outname$count -c $currout/$outname$count.pdb -cpi $prevout/part$prev -cpt 100 -cpo $currout/part$count.cpt -maxh 11
     fi
     echo making files$count
 let count=$count+1
 let oldcount=$count-1
 echo $count > $topdir/count.txt
 echo starting new$count
 if [ -f $currout/$outname$oldcount.pdb ]; then
    echo "Reached 20ns Job Complete"
    cp $parent/pulloutwater.sh $topdir
    cd $topdir
    sbatch pulloutwater.sh
 elif [ $count -gt $maxiter ]; then 
    echo "Reached $maxiter iterations.  Job Complete"
 else
   mkdir -p $topdir/part$count
   cd $topdir 
   sbatch -J ro-$round-cou-$count $script
 fi
