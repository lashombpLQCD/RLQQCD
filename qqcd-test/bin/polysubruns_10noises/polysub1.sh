#!/bin/bash
#SBATCH -p flat-quadrant   		# Queue (partition) name 
#SBATCH -n 108 					# Total number of tasks
#SBATCH -t 48:00:00	 			# Time (HH:MM:SS)
#SBATCH -N 27					# Total number of nodes 
#SBATCH --ntasks-per-node=4     # Number of tasks per node 
#SBATCH -o polysub1_run.o%j
#SBATCH -e polysub1_run.e%j
#SBATCH --mail-user=paul_lashomb@baylor.edu
#SBATCH --mail-type=all

###export IBRUN_TASKS_PER_NODE = 4



###module add mvapich2-intel/2.2-2 # not necessary as modules already added

###cd $PBS_O_WORKDIR # old PBS 

###cd $SLURM_SUBMIT_DIR # new Slurm equiv

###num=`cat $PBS_NODEFILE | wc -l` # old PBS


###num = `cat $SLURM_JOB_NODELIST | wc -l` # new Slurm equiv
###echo $SLURM_JOB_NODELIST

###num = `$SLURM_JOB_NODELIST`


###echo "Requested processors: $num"
###echo "Node(s):"

###uniq $PBS_NODEFILE # old PBS
###uniq $SLURM_JOB_NODELIST # new Slurm equiv
echo

echo "Cleaning scratch files"
##rm /data/barals/qqcd/scratch-workinglagfib10/es.dat
##rm /data/barals/qqcd/scratch-workinglagfib10/hfes.dat

###----------Kodiak----------

##rm /data/lashombp/qqcd/submethodsppminresdr/DUBLIN.OUT
##rm /data/lashombp/qqcd/scratch/CFGSPROPS.LOG
##rm /data/lashombp/qqcd/submethodsppminresdr/EIGEN_VALS.LOG



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


##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/es.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/hfes.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/ns.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/hfespoly.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/pp.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/pp7.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/pp4.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/hfesps.dat
##rm /data/lashombp/qqcd/submethodsppminresdr/scratch1/esps.dat


###-----------Stampede2-------------

dats="/work/00526/tg457765/stampede2/qqcd/polysub/scratch1/*.dat"
logs="/work/00526/tg457765/stampede2/qqcd/polysub/scratch1/*.LOG"
outs="/work/00526/tg457765/stampede2/qqcd/polysub/scratch1/*.OUT"
pths="/work/00526/tg457765/stampede2/qqcd/polysub/scratch1/"



### Removes .dat files
for filename in $dats; do 
	[ -e "$filename" ] || continue 
	rm $filename   
done 

### Removes .LOG files
for filename in $logs; do 
	[ -e "$filename" ] || continue 
		if [ "$filename" == "$pths""hopping-matrix.LOG" ]; then 
		echo "skipped hopping-matrix.LOG" 
		else 
		rm $filename
		fi  
done 

### Removes .OUT files
for filename in $outs; do 
	[ -e "$filename" ] || continue 
	rm $filename 
done 






###file1="DUBLIN.OUT"
###file2="CFGSPROPS.LOG"
###file3="EIGEN_VALS.LOG"



###i="4"
###n="15"

###file4="linearresidual.dat"
####file5="es.dat"
###file6="hfes.dat"
###file7="ns.dat"
###file8="hfespoly.dat"
###file9="pp.dat"
###file10="pp7.dat"
###file11="newpp.dat"
###file12="pp4.dat"
###file13="hfesps.dat"
###file14="esps.dat"

###while [ $i -lt $n ]
###do
###if [ -f "$polypath1








echo

echo "Job starting at `date`"
echo

START=`date '+%s'`

###export MV2_ENABLE_AFFINITY=0


###mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout1 >quick1.out # old PBS


### -n specifies the number of process to run
### arbitrary allows you to specify the list of hosts
### -w requests a specific list of hosts
### $SLURM_JOB_NODELIST is the file that contains the allocated hostnames 


###export OMP_NUM_THREADS=1
###OMP_NUM_THREADS=1



ibrun ./qqcd-testout >quick1.out # Stampede2 uses ibrun for MPI, not mpiexec  

### ibrun ./qqcd-test1 >quick1.out
### ibrun -m arbitrary -w $SLURM_JOB_NODELIST ./qqcd-test1 >quick1.out


END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
