# RLQQCD

The Randy Lewis qqcd program in FORTRAN. This version uses the new method of forming the subtraction polynomial, as outlined in out [High-degree subtraction polynomials](https://arxiv.org/pdf/2306.06188.pdf) paper.


### Basic File Structure
- `qqcd/andyworkinglagfib/qqcd-test/`: Base directory of program 
  - `cfgsprops/`
    - `cfgsprops.f90`: Performs all calculations from generating/loading gauge configurations to calling subroutines in `quark/` to perform subtraction calculations. 
    - `quark/`: Contains all of the subtraction source files 
  - `glue/`: Where gauge configurations are generated, saved, and read in
  - `user/`: Settings for a particular run (e.g. configuration directories, lattice size, etc.). The FORTRAN program starts in `cfgspropsmain.f90`. 
  - `bin/`: Directory for compiled program and batch submission files for running on HPC cluster
    
### Basic Operation 
- Set desired user settings in the `qqcd/andyworking/qqcd-test/user` directory. 
- Set directories and compile the program by modifying and executing the *make* file in `qqcd/andyworklagfib/qqcd-test/sbin#`. (sbin# set in `qqcd/andyworking/qqcd-test/user/cfgspropsmain.f90`) 
- Each job is submitted using a .sh batch submission file inside `qqcd/andyworkinglagfib/qqcd-test/bin` containing the compiled program, `qqcd-testout#`.


## Common Operations
### Printing the Hopping Matrix 
The quark matrix comes in the form $`M = I - kB`$, where $`B`$ is the hopping matrix. The hopping matrix from FORTRAN can be exported to a .LOG file, which can then can then be converted to a .mat file for Matlab using `hoppingconvert.m` in our FORTRAN utilities repository: [https://github.com/lashombpLQCD/LQCDUtilities](https://github.com/lashombpLQCD/LQCDUtilities). 

To print the .LOG file, use the following steps: 
1. In `qqcd/andyworkinglagfib/qqcd-test/user/latsize.f90`, set the dimensions of the lattice used (nx, ny, nz, nt) and how many processors are used for each dimension (npx, npy, npz, npt).
2. In `qqcd/andyworkinglagfib/qqcd-test/user/cfgspropsmain.f90`, set the variable `rwdir` to the desired directory for printing the .LOG file.
3. In `qqcd/andyworkinglagfib/qqcd-test/user/sbin#`, edit the *make* file to match the directories used, and then compile the program by running the `make` command.
4. Edit your job submission file in `qqcd/andyworklagfig/qqcd-test/bin` if needed, and submit the job to the cluster.
5. Repeat each step for the remaining configurations. 


Once the .mat file has been made from the .LOG file, import the .mat file into Matlab by simply loading the .mat file and creating the quark matrix via, 
```
load h01_24242424.mat

n = size(B,1);
I = eye(n,n);
k = 0.1570;    % Use whatever value of kappa used in FORTRAN  
A = I - kB;

clear B; 
```

### Running the Subtraction Routines 
The main subtraction routines are within the file `qqcd/andyworkinglagfib/qqcd-test/cfgsprops/quark/disconloops.f90`. The Monte Carlo trace calculation is performed within the subroutine `subroutine testFUNC`. 
