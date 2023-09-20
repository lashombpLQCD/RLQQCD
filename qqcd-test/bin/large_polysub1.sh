#!/bin/sh
#PBS -l nodes=27:ppn=4
###PBS -l nodes=2:ppn=4,walltime=00:30:00
###PBS -l walltime=00:10:00
###PBS -m ea -M Paul_Lashomb@baylor.edu
###PBS -l nodes=27:ppn=4

module purge 
module add intel/19.2 
#module add mvapich2-intel/2.2-2

export FI_PROVIDER=verbs 

cd $PBS_O_WORKDIR

num=`cat $PBS_NODEFILE | wc -l`
echo "Requested processors: $num"
echo "Node(s):"
uniq $PBS_NODEFILE
echo

echo "Cleaning scratch files"
##rm /data/barals/qqcd/scratch-workinglagfib10/es.dat
##rm /data/barals/qqcd/scratch-workinglagfib10/hfes.dat



rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/DUBLIN.OUT
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/CFGSPROPS.LOG
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/EIGEN_VALS.LOG

rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/EIGEN_VALS.LOG
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/CFGSPROPS.LOG
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/DUBLIN.OUT
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/es.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/esps.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/hfes.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/hfespoly.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/hfesps.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/linearresidual.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/ns.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/pp4.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/pp7.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/pp.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/ps.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/newpp.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/newhfespoly.dat
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/newpp_noise-1.LOG 
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/oldpp_noise-1.LOG 
rm /data/lashombp/qqcd/preconditioning_12121216/scratch1/prepp_noise-1.LOG 



##rm /data/lashombp/qqcd/submethodsppminresdr/trueresidual.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/residual.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/eigresidual.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/wresidual.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/orthog.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/vn.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/htest.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/loweigresidual.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/tctest.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/ritzvectsorthog.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/ritztest.dat



##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/es.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/hfes.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/ns.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/hfespoly.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/pp.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/pp7.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/pp4.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/hfesps.dat
##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/esps.dat

##rm /data/lashombp/qqcd/polysub_12121216_10noises/scratch1/linearresidual.dat


echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout1 >quick1.out

### added -verbose before 

OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
