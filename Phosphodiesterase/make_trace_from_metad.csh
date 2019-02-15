#!/bin/bash
#SBATCH -J convert
#SBATCH -o %j.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:29:00
#SBATCH -p debug

### This script makes a trace of the center of mass of a small molecule
###  as it is pushed out of the binding pocket from metadynamics
### The TCL script created can be loaded into VMD

#source /apps2/Modules/3.2.6/init/tcsh

module load intelics/2016.1-full plumed/2.2.1 gromacs/5.0.1-plumed

topdir=/scratch/pde6_pull/new/noion/trial1_half

name1=PDE6_water.xtc #full

name2=PDE6_nw.xtc # no water whole

name3=PDE6_nw2.xtc # pbc nojump

name4=PDE6_fit2t5tpr.xtc #skip dt 100

gmx_mpi trjcat -f $topdir/part*/*.xtc -o $topdir/$name1

echo 24 | gmx_mpi trjconv -s $topdir/part1/run.tpr -f $topdir/$name1 -n $topdir/index.ndx -o $topdir/$name2 -pbc whole -e 40248 

echo 24 | gmx_mpi trjconv -s $topdir/part1/run.tpr -f $topdir/$name2 -n $topdir/index.ndx -o $topdir/$name3 -pbc nojump

echo 4 24 | gmx_mpi trjconv -s $topdir/../trial5_half/part1/run.tpr -f $topdir/$name3 -n $topdir/index.ndx -o $topdir/$name4 -fit rot+trans -dt 100

echo 13 | gmx_mpi traj -s $topdir/../trial5_half/part1/run.tpr -n $topdir/index.ndx -f $topdir/$name4 -ox drugcom4.xvg -com -tu ns -nojump

rm noheader4.xvg
rm noheaderadjusted.xvg
rm both.xvg
rm init.tcl
rm init2.tcl
rm pde6_1pathway_sm.tcl

grep -v '#' drugcom4.xvg | grep -v '@' > noheader4.xvg

x1=`head -n 1 noheader4.xvg | awk '{print$2}'`

y1=`head -n 1 noheader4.xvg | awk '{print$3}'`

z1=`head -n 1 noheader4.xvg | awk '{print$4}'`

x5=`head -n 1 ../trial5_half/noheader4.xvg | awk '{print$2}'`

y5=`head -n 1 ../trial5_half/noheader4.xvg | awk '{print$3}'`

z5=`head -n 1 ../trial5_half/noheader4.xvg | awk '{print$4}'`

x=`echo '$x5 - $x1' | bc -l`

y=`echo '$y5 - $y1' | bc -l`

z=`echo '$z5 - $z1' | bc -l`

echo "-1      0      0      0" > noheaderadjusted.xvg

cat smooth_pde6_t1_300.dat >> noheaderadjusted.xvg

paste smooth_pde6_t1_300.dat noheaderadjusted.xvg > both.xvg

awk < both.xvg -v a="$x" -v b="$y" -v c="$z" '{print "graphics top line {" 10*$6-a " " 10*$7-b " " 10*$8-c "} {" 10*$2-a " " 10*$3-b " " 10*$4-c "} width 2 style solid\n"}' >  init.tcl

len=`wc -l init.tcl | awk '{print $1}'`

echo 'set center {0 0 0}' > init2.tcl

minus=$((len-1))

minus2=$((len-2))

tail -n $minus init.tcl >> init2.tcl

head -n $minus2 init2.tcl > pde6_1pathway_sm.tcl
