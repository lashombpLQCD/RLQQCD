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
##rm /mata/barals/qqcd/scratch-workinglagfib10/es.dat
##rm /mata/barals/qqcd/scratch-workinglagfib10/hfes.dat
rm /mata/lashombp/qqcd/scratch5/DUBLIN.OUT
##rm /mata/whytet/qqcd/scratch/CFGSPROPS.LOG
rm /mata/lashombp/qqcd/scratch5/EIGEN_VALS.LOG
##rm /mata/whytet/qqcd/scratch-5/trueresidual.dat
##rm /mata/whytet/qqcd/scratch-5/residual.dat
##rm /mata/whytet/qqcd/scratch-5/eigresidual.dat
##rm /mata/whytet/qqcd/scratch-5/wresidual.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/es.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/hfes.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/ns.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/hfespoly.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/pp.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/pp7.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/pp4.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/hfesps.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch5/esps.dat


echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout5 >quick5.out
OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
