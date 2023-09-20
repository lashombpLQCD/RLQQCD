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
rm /data/lashombp/qqcd/submethodsminresdr/DUBLIN.OUT
##rm /data/lashombp/qqcd/scratch/CFGSPROPS.LOG
rm /data/lashombp/qqcd/submethodsminresdr/EIGEN_VALS.LOG
##rm /data/lashombp/qqcd/submethodsminresdr/trueresidual.dat
##rm /data/lashombp/qqcd/submethodsminresdr/residual.dat
##rm /data/lashombp/qqcd/submethodsminresdr/eigresidual.dat
##rm /data/lashombp/qqcd/submethodsminresdr/wresidual.dat
##rm /data/lashombp/qqcd/submethodsminresdr/orthog.dat
##rm /data/lashombp/qqcd/submethodsminresdr/vn.dat
##rm /data/lashombp/qqcd/submethodsminresdr/htest.dat
##rm /data/lashombp/qqcd/submethodsminresdr/loweigresidual.dat
##rm /data/lashombp/qqcd/submethodsminresdr/tctest.dat
##rm /data/lashombp/qqcd/submethodsminresdr/ritzvectsorthog.dat
##rm /data/lashombp/qqcd/submethodsminresdr/ritztest.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/es.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/hfes.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/ns.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/hfespoly.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/pp.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/pp7.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/pp4.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/hfesps.dat
rm /data/lashombp/qqcd/submethodsminresdr/scratch1/esps.dat



echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout1 >quick1.out
OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
