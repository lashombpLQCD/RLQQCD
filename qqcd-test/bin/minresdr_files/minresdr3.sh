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
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/DUBLIN.OUT
##rm /data/whytet/qqcd/scratch/CFGSPROPS.LOG
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/EIGEN_VALS.LOG
##rm /data/whytet/qqcd/scratch-3/trueresidual.dat
##rm /data/whytet/qqcd/scratch-3/residual.dat
##rm /data/whytet/qqcd/scratch-3/eigresidual.dat
##rm /data/whytet/qqcd/scratch-3/wresidual.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/es.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/hfes.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/ns.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/hfespoly.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/pp.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/pp7.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/pp4.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/hfesps.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch3/esps.dat



echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout3 >quick3.out
OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
