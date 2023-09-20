#!/bin/sh
#PBS -l nodes=36:ppn=4
###PBS -l nodes=2:ppn=4,walltime=00:30:00

#PBS -m ea -M Travis_Whyte@baylor.edu

module add mvapich2-intel/2.2-2

cd $PBS_O_WORKDIR

num=`cat $PBS_NODEFILE | wc -l`
echo "Requested processors: $num"
echo "Node(s):"
uniq $PBS_NODEFILE
echo

echo "Cleaning scratch files"
##rm /data/barals/qqcd/scratch-workinglagfib10/es.dat
##rm /data/barals/qqcd/scratch-workinglagfib10/hfes.dat
rm /data/lashombp/qqcd/scratch2/DUBLIN.OUT
##rm /data/whytet/qqcd/scratch/CFGSPROPS.LOG
rm /data/lashombp/qqcd/scratch2/EIGEN_VALS.LOG
##rm /data/whytet/qqcd/scratch-2/trueresidual.dat
##rm /data/whytet/qqcd/scratch-2/residual.dat
##rm /data/whytet/qqcd/scratch-2/eigresidual.dat
##rm /data/whytet/qqcd/scratch-2/wresidual.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/es.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/hfes.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/ns.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/hfespoly.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/pp.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/pp7.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/pp4.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/hfesps.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch2/esps.dat

echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout2 >quick2.out
OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
