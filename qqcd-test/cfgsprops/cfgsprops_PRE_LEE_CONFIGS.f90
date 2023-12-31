! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! cfgsprops.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This is the main module for the generation of gauge field configurations
! and the construction of quark propagators.
! It also computes the average plaquette, tadpole factors, Wilson loops,
! and the pion, rho and nucleon correlators.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module cfgsprops

!    use MPI
    use kinds
    use latdims
    use basics
    use lagfib
    use lattice
    use gaugetools
    use debug
    use heatbath
    use diracops
    use quark
    use disconloops
    use gmresrhs
    use shift
    use vevbleft
    use printops !BS
    implicit none
    private


   EXTERNAL IARGC,GETARG
! Use the following line if the MPI module is not available.
   include 'mpif.h'

! Define access to subroutines.
    public  :: generator
    private :: sweep, somewloops, perm, check4ones

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine generator(rwdir,newcalc,ifilesave,cfgfile,iseed,ncold,nsweep, &
                      n0sweep,n1sweep,nhit,gaction,alat,beta,tad,tadupdate, &
                      tadpole,tadlandau,landimp,alpha,omega,thetamax,nwilo, &
                      wlat,Rmax,npaths,nfuzz,epsfuzz,src,bc,resmax,itermin, &
                      nsmear,asmear,nkappa,kappa,omegaM3R,ntmqcd, &
                      mtmqcd,cSWtmqcd,inverter,nGMRES,opSST,tSST,pSST,nSST, &
                      mSST,cSWSST,invSST,numnoises,kappaloop,cSWloop,docorr, &
                      nkappacorr,ntmqcdcorr,ntmqcdloop,mtmqcdloop,nSSTcorr,&
                      nucleonparticle,gitest,hoption,wilsonmm)

! Main subroutine for gauge field generation and propagator construction.

    integer(kind=KI), intent(in)    :: newcalc, iseed, ncold, nsweep, nhit, &
                                       tadlandau, nwilo, npaths, nfuzz, &
                                       gitest, nsmear, nkappa, ntmqcd, &
                                       itermin, docorr, opSST, & !moved inverter to intent(inout) TW 11/17/15                                      
                                       tSST, invSST, numnoises, nkappacorr, &
                                       ntmqcdcorr, ntmqcdloop, &
                                       nucleonparticle,hoption,wilsonmm
    integer(kind=KI), intent(inout) :: n0sweep, n1sweep, tadupdate, ifilesave,inverter            
    character(len=*), intent(in),    dimension(:)    :: rwdir
    integer(kind=KI), intent(inout)                  :: tadpole
    character(len=*), intent(inout)                  :: cfgfile
    integer(kind=KI), intent(inout), dimension(:)    :: gaction, nGMRES
    real(kind=KR),    intent(inout)                  :: beta, epsfuzz
    real(kind=KR),    intent(inout), dimension(:)    :: alat, tad, Rmax
    integer(kind=KI), intent(in),    dimension(:)    :: landimp, bc, &
                                                        pSST, nSST, nSSTcorr
    real(kind=KR),    intent(in)                     :: alpha, omega, thetamax,&
                                                        resmax, asmear, &
                                                        cSWloop
    real(kind=KR),    intent(in),    dimension(2)    :: kappaloop
    real(kind=KR),    intent(inout)                  :: cSWtmqcd, cSWSST
    integer(kind=KI), intent(in),    dimension(:,:)  :: wlat
    integer(kind=KI), intent(in),    dimension(0:,:) :: src
    real(kind=KR),    intent(inout), dimension(:)    :: kappa
    real(kind=KR),    intent(in),    dimension(:,:)  :: mtmqcd, mSST,mtmqcdloop
    real(kind=KR2),   intent(in)                     :: omegaM3R

!*LABELLING CONVENTION FOR LATTICE SITES:
! Because next-nearest neighbour interactions may be present in some
! spacetime directions, the lattice will be blocked so there are two
! "partial" global lattices per direction.  (For coding simplicity,
! this is done even when next-nearest neighbour interactions are not used.)
! These are referred to as the two "blocks" of the lattice in each
! direction, so the lattice has 16 blocks in total.
!
! 2-dimensional example of lattice site labelling scheme
! (isg,even/odd,block):
!         (7,o,3) (7,o,4) (7,e,3) (7,e,4) (8,o,3) (8,o,4) (8,e,3) (8,e,4)
!         (7,o,1) (7,o,2) (7,e,1) (7,e,2) (8,o,1) (8,o,2) (8,e,1) (8,e,2)
!         (5,e,3) (5,e,4) (5,o,3) (5,o,4) (6,e,3) (6,e,4) (6,o,3) (6,o,4)
! y       (5,e,1) (5,e,2) (5,o,1) (5,o,2) (6,e,1) (6,e,2) (6,o,1) (6,o,2)
! /\      (3,o,3) (3,o,4) (3,e,3) (3,e,4) (4,o,3) (4,o,4) (4,e,3) (4,e,4)
! |       (3,o,1) (3,o,2) (3,e,1) (3,e,2) (4,o,1) (4,o,2) (4,e,1) (4,e,2)
! |___\ x (1,e,3) (1,e,4) (1,o,3) (1,o,4) (2,e,3) (2,e,4) (2,o,3) (2,o,4)
!     /   (1,e,1) (1,e,2) (1,o,1) (1,o,2) (2,e,1) (2,e,2) (2,o,1) (2,o,2)

!*SOME BASIC PARAMETERS THAT ARE DEFINED IN module basics:
! nps = total number of processes.
! nvhalf is the number of even (or odd) lattices sites per process per block.
! nbmax is the number of spacetime sites in the largest boundary
!       shared by two processes.  (The boundary is 3-D.)
!*LINK VARIABLES:
!    real(kind=KR), dimension(18,ntotal,4,2,16)     :: u
    real(kind=KR),allocatable, dimension(:,:,:,:,:)     :: u
! Link variables are stored as u(a,b,c,d,e) where
! a=1,18 runs over Re and Im parts of the 3x3 colour matrix as follows:
!             /  u(1,i)+i*u(2,i)    u(3,i)+i*u(4,i)    u(5,i)+i*u(6,i)  \
!        u = |   u(7,i)+i*u(8,i)    u(9,i)+i*u(10,i)  u(11,i)+i*u(12,i)  |
!             \ u(13,i)+i*u(14,i)  u(15,i)+i*u(16,i)  u(17,i)+i*u(18,i) /
!        where i=sqrt(-1).
! b=1,ntotal runs over even or odd sites on one block of a sublattice, plus
!            room for the largest boundary shared by two processes.
! c=1,4 runs over the x,y,z,t spacetime directions.
! d=1,2 runs over the even(1) and odd(2) lattice sites.
! e=1,16 runs over the blocks of a sublattice.


!*FERMION SPINORS:
!    real(kind=KR), dimension(6,ntotal,4,2,8,nshifts) :: xe, xo
!    real(kind=KR), dimension(6,ntotal,4,2,8)     :: be, bo
    real(kind=KR),allocatable, dimension(:,:,:,:,:,:) :: xe, xo
    real(kind=KR),allocatable, dimension(:,:,:,:,:)     :: be, bo
! Source vectors are stored in two vectors, be(a,b,c,d,e) and bo(a,b,c,d,e),
! and resultant propagators are stored as xe(a,b,c,d,e,f) and xo(a,b,c,d,e,f),
! according to their checkerboard "colour" on the global lattice. 
! The arguments match those of the link variables (see above), with the
! first argument ordered as follows,
!            / v(1)+i*v(2) \
!        v = | v(3)+i*v(4) |, where i=sqrt(-1),
!            \ v(5)+i*v(6) /
! and the third argument, c=1,4 runs over Dirac index instead of x,y,z,t.
! The final argument of xe,xo is to accommodate multiple hopping parameters.

!*SOME BASIC VARIABLE DEFINITIONS:
    integer(kind=KI)                         :: numprocs, myid
    integer(kind=KI)                         :: is, ishift, ntm
    integer(kind=KI)                         :: isignal, ngmresrhs
    integer(kind=KI), dimension(4,2)         :: nn
    logical,          dimension(4)           :: ldiv
    integer(kind=KI), dimension(4)           :: nms
    integer(kind=KI), dimension(nvhalf,4,2)  :: lv, lvbc
    integer(kind=KI), dimension(nbmax,4,2,2) :: ib
    logical,          dimension(16,4)        :: lbd
    integer(kind=KI), dimension(16,4)        :: iblv
    integer(kind=KI), dimension(2,8)         :: vecbl
    integer(kind=KI), dimension(2,16)        :: vecblinv
! numprocs = total number of processes assigned by MPI (must equal nps).
! myid = number (i.e. single-integer address) of current process
!        (value = 0, 1, 2, ...).
! nn(j,1) = single-integer address, on the grid of processes, of the
!           neighbour to myid in the +j direction.
! nn(j,2) = single-integer address, on the grid of processes, of the
!           neighbour to myid in the -j direction.
! ldiv(mu) is true if there is more than one process in the mu direction,
!          otherwise it is false.
! nms(mu) is number of even (or odd) boundary sites per block between
!         processes in the mu'th direction IF there is more than one
!         process in this direction.
! lv(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
! lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                     If there is only one process in the mu direction,
!                     then lvbc still points to the buffer whenever
!                     non-periodic boundary conditions are needed.
! ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                    block of the sublattice.
! ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                    block of the sublattice.
! lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!             direction, otherwise it is false.
! iblv(ibl,mu) is the number of the neighbouring block (to the block
!              numbered ibl) in the mu direction.
!              Both ibl and iblv() run from 1 through 16.
! vecbl(1,:) are the blocks of the lattice where all sites are globally "even".
! vecbl(2,:) are the blocks of the lattice where all sites are globally "odd".
! vecblinv() is an inverse for vecbl().

!*LOCAL VARIABLE DEFINITIONS:
    integer(kind=KI) :: ierr, icfgsave, icount, ikappa, isweep, myidin, &
                        firstcfg, MRT, MRT2, icsrc, idsrc, mrvectype, &
                        cgvectype, itmflag, checksum1, checksum2, idebug, &
                        SSTvectype, iSST, iSSTflag, iloopflag, nLS, iLS, &
                        icfgin, ibit, nsmrsnk, irhs, rhs
    real(kind=KR)                       :: rtest, itest, plaq, cSWin, kapin, &
                                           muin
    real(kind=KR), dimension(4)         :: tadmu, tadrun
    real(kind=KR), dimension(4,4,4)     :: coact
    real(kind=KR), dimension(2)         :: kaptm
    real(kind=KR), dimension(nshifts)   :: mutemp
    real(kind=KR)                       :: mudelta, muWilson
    character(len=128)                  :: lagfibfile,lagfibstore,propfile
    real(kind=KR2)                      :: timeinit, timefin, sweeptime
    character(len=1), dimension(2)      :: LSsymbol
    character(len=5)                    :: trailer
!   integer(kind=KI), parameter         :: kmax=9
    real(kind=KR), dimension(nt,kmax,4) :: pion, rho, decay
    real(kind=KR), dimension(2*kmax,4)  :: pionSST
    real(kind=KR), dimension(nshifts)   :: sigma

! The last index is ii=1,2 for be and bo stored in btemp....
!    real(kind=KR), dimension(6,ntotal,2,2,8,2)       :: btemp
    real(kind=KR),allocatable, dimension(:,:,:,:,:,:):: btemp
!    real(kind=KR), dimension(6,ntotal,2,2,8,1,2)     :: xtemp
    real(kind=KR), allocatable,dimension(:,:,:,:,:,:,:)     :: xtemp
   
    integer(kind=KI)                    :: sweepnum, sweeptotal
    integer(kind=KI)                    :: iparticle, ir, irloop, irloop1, irloop2
    real(kind=KR)                       :: time1, time2,time11,time22
    integer(kind=KI),dimension(4)       :: idg5
    integer(kind=KI),dimension(4,4)     :: signg5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(kind=KI)                     :: nrhs,exists
    integer(kind=KI)                     :: ii,jj,kk,ll,mm,iaa,ibb,icc,idd,iee,coll,evodd,roww,idag, &
                                                           prc

    real(kind=KR2),     dimension(nxyzt,3,4)           :: realz2noise,imagz2noise
    real(kind=KR2)                                     :: rz2, iz2,indxx
    real(kind=KR2),  dimension(6,ntotal,4,2,8)         :: evenVin, evenVout
    real(kind=KR2),  dimension(6,ntotal,4,2,8)         :: oddVin, oddVout
    



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Added for testing the status of allocatable arrays 
    integer(kind=KI)                                   :: statu,statxe,statxeo,&                                                          statbe,statbo,statxt,&                                                          statbt 




!call npsvalue
!print *, "nps=",nps
   !BS set nsweep equal to one for newcalc =2 


    idg5(:)=(/3,4,1,2/)
    signg5(1,:)=(/1,1,-1,-1/)
    signg5(2,:)=(/1,1,-1,-1/)
    signg5(3,:)=(/-1,-1,1,1/)
    signg5(4,:)=(/-1,-1,1,1/)    
!BSprint *,"hello b"!BS


!print*, "At Line 249 before MPI_INIT"  







! Initialize MPI
     call MPI_INIT(ierr)
!BSprint *,"hello c"!BS



!print*, "At Line 263 after MPI_INIT" 




! Ensure that nx/4/npx, ny/4/npy, nz/4/npz, nt/4/npt are integers.
    rtest = real(nx,KR)/4.0/real(npx,KR)
!    print *,"rtest=",rtest
    itest = real(nx/4/npx,KR)
!    print *,"itest=",itest
    if (abs(rtest-itest)>0.000001_KR) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nx/4/npx is not an integer: nx,npx = ", nx, npx
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'nx/4/npx is not an integer: nx,npx'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif
    rtest = real(ny,KR)/4.0/real(npy,KR)
    itest = real(ny/4/npy,KR)
    if (abs(rtest-itest)>0.000001_KR) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "ny/4/npy is not an integer: ny,npy = ", ny, npy
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'ny/4/npy is not an integer: ny,npy'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif
    rtest = real(nz,KR)/4.0/real(npz,KR)
    itest = real(nz/4/npz,KR)
    if (abs(rtest-itest)>0.000001_KR) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nz/4/npz is not an integer: nz,npz = ", nz, npz
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'nz/4/npz is not an integer: nz,npz'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif
    rtest = real(nt,KR)/4.0/real(npt,KR)
    itest = real(nt/4/npt,KR)
    if (abs(rtest-itest)>0.000001_KR) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nt/4/npt is not an integer: nt,npt = ", nt, npt
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'nt/4/npt is not an integer: nt,npt'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif

! Verify that nk1 has been set to a sufficiently large value.
    if (nk1<nkappa) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nk1 is less than nkappa: nk1,nkappa = ", nk1, nkappa
     close(unit=8,status="keep")
    endif

! Don't allow the qmrgam5 inverter to be used for twisted mass QCD.
    rtest = abs(mtmqcd(1,2)) + abs(mtmqcd(2,2)) + abs(mtmqcd(3,2)) &
          + abs(mtmqcd(4,2)) + abs(mtmqcd(5,2)) + abs(mtmqcd(6,2)) &
          + abs(mtmqcd(7,2)) + abs(mtmqcd(8,2)) + abs(mtmqcd(9,2))
    if (inverter==3 .and. rtest>0.0001_KR) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "wrong inverter for tmQCD: rtest = ", rtest
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'wrong inverter for tmQCD: rtest'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif

! Verify that the array set for GMRES(n),GMRES-DR(n,k) is sufficiently large.
    if (ntmqcd/=0 .and. inverter>3 .and. nGMRES(1)>nmaxGMRES) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nGMRES is larger than nmaxGMRES."
      write(unit=8,fmt=*) "nGMRES = ", nGMRES
      write(unit=8,fmt=*) "nmaxGMRES = ", nmaxGMRES
     close(unit=8,status="keep")
     print *, "Calling MPI_ABORT near 'nGMRES is larger than nmaxGMRES'"
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     stop
    endif

! Verify that k<n and k<nGMRES(2) in GMRES-DR(n,k).
    if (ntmqcd/=0 .and. inverter==5 .and. &
       (nGMRES(2)>nGMRES(1).or.nGMRES(2)>kmaxGMRES)) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nGMRES(2) is too large."
      write(unit=8,fmt=*) "nGMRES = ", nGMRES
      write(unit=8,fmt=*) "kmaxGMRES = ", kmaxGMRES
     close(unit=8,status="keep")
    endif

! Verify that k<n and k<nGMRES(2) in GMRES-DRshift(n,k)
    if (ntmqcd/=0 .and. inverter==6 .and. &
       (nGMRES(2)>nGMRES(1).or.nGMRES(2)>kmaxGMRES)) then
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "nGMRES(2) is too large."
      write(unit=8,fmt=*) "nGMRES = ", nGMRES
      write(unit=8,fmt=*) "kmaxGMRES = ", kmaxGMRES
     close(unit=8,status="keep")
    endif

! Initiate MPI, confirm number of MPI processes and determine myid.
    myid = 0
!   if (nps==1) then
!    myid = 0
!   else
! MPI INIT is moved to top so that MPI_ABORT is initalized.
!    call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
     if (numprocs/=nps) then
      open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
           ,form="formatted")
       write(unit=8,fmt=*) "numprocs,nps = ",numprocs,nps
      close(unit=8,status="keep")
      call MPI_FINALIZE(ierr)
      stop
     endif
!   endif


    ! Added status variables -PL 7-15-22 
    allocate(u(18,ntotal,4,2,16),stat=statu)
    allocate(xe(6,ntotal,4,2,8,nshifts),stat=statxe)
    allocate(xo(6,ntotal,4,2,8,nshifts),stat=statxo)
    allocate(be(6,ntotal,4,2,8),stat=statbe)
    allocate(bo(6,ntotal,4,2,8),stat=statbo)
    allocate(btemp(6,ntotal,2,2,8,2),stat=statbt)
    allocate(xtemp(6,ntotal,2,2,8,1,2),stat=statxt)



! Record the user-defined parameters.
if(myid==0) then
        inquire(file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", exist=exists)
  if (.not. exists) then

open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",status="new",     &
        action="write",form="formatted")
   close(unit=8,status="keep")
 endif
endif










    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write",&
          form="formatted",status="old")
      write(unit=8,fmt="(a1)")           " "
      write(unit=8,fmt="(a12)")          "USER INPUTS:"
      write(unit=8,fmt="(a1)")           " "
      write(unit=8,fmt="(a18,4i10)")     "nx,ny,nz,nt =     ", nx,ny,nz,nt
      write(unit=8,fmt="(a18,4i10)")     "npx,npy,npz,npt = ", npx,npy,npz,npt
      write(unit=8,fmt="(a18,i10)")      "nk1 =             ", nk1
      write(unit=8,fmt="(a18,3i10)")     "KI, KR, KR2 =     ", KI, KR, KR2
      write(unit=8,fmt="(a12,i10)")      "newcalc =   ", newcalc
      write(unit=8,fmt="(a12,i10)")      "ifilesave = ", ifilesave
      write(unit=8,fmt="(a12,a64)")      "cfgfile =   ", cfgfile
      write(unit=8,fmt="(a12,i10)")      "iseed =     ", iseed
      write(unit=8,fmt="(a12,i10)")      "ncold =     ", ncold
      write(unit=8,fmt="(a12,i10)")      "nsweep =    ", nsweep
      write(unit=8,fmt="(a12,i10)")      "n0sweep =   ", n0sweep
      write(unit=8,fmt="(a12,i10)")      "n1sweep =   ", n1sweep
      write(unit=8,fmt="(a12,i10)")      "nhit =      ", nhit
      write(unit=8,fmt="(a12,4i10)")     "gaction =   ", gaction
      write(unit=8,fmt="(a12,4es17.10)") "alat =      ", alat
      write(unit=8,fmt="(a12,es17.10)")  "beta =      ", beta
      write(unit=8,fmt="(a12,4es17.10)") "tad =       ", tad
      write(unit=8,fmt="(a12,i10)")      "tadupdate = ", tadupdate
      write(unit=8,fmt="(a12,i10)")      "tadpole =   ", tadpole
      write(unit=8,fmt="(a12,i10)")      "tadlandau = ", tadlandau
      write(unit=8,fmt="(a12,4i10)")     "landimp =   ", landimp
      write(unit=8,fmt="(a12,es17.10)")  "alpha =     ", alpha
      write(unit=8,fmt="(a12,es17.10)")  "omega =     ", omega
      write(unit=8,fmt="(a12,es17.10)")  "thetamax =  ", thetamax
      write(unit=8,fmt="(a12,i10)")      "nwilo =     ", nwilo
      if (nwilo/=0) then
       do icount = 1,abs(nwilo)
        write(unit=8,fmt="(a14,i6,4i10)")"i,wlat(i,:) = ", icount,wlat(icount,:)
        write(unit=8,fmt="(a14,i6,es17.10)")"i,Rmax(i) =   ",icount,Rmax(icount)
       enddo ! icount
      endif
      write(unit=8,fmt="(a12,i10)")      "npaths =    ", npaths
      write(unit=8,fmt="(a12,i10)")      "nfuzz =     ", nfuzz
      write(unit=8,fmt="(a12,es17.10)")  "epsfuzz =   ", epsfuzz
      write(unit=8,fmt="(a12,4i10)")     "src(0,:) =  ", src(0,:)
      write(unit=8,fmt="(a12,4i10)")     "bc =        ", bc
      write(unit=8,fmt="(a12,es17.10)")  "resmax =    ", resmax
      write(unit=8,fmt="(a12,i10)")      "itermin =   ", itermin
      write(unit=8,fmt="(a12,i10)")      "nsmear =    ", nsmear
      write(unit=8,fmt="(a12,es17.10)")  "asmear =    ", asmear
      write(unit=8,fmt="(a12,i10)")      "nkappa =    ", nkappa
      write(unit=8,fmt="(a12,i10)")      "wilsonmm =  ", wilsonmm
      write(unit=8,fmt="(a12,i10)")      "hopton =    ", hoption
      if (nkappa>0) then
       do icount = 1,nkappa
        write(unit=8,fmt="(a14,i6,es17.10)") &
                                         "i,kappa(i) =  ", icount,kappa(icount)
       enddo ! icount
      endif
      write(unit=8,fmt="(a12,es17.10)")  "omegaM3R =  ", omegaM3R
      write(unit=8,fmt="(a12,i10)")      "ntmqcd =    ", ntmqcd
      if (abs(ntmqcd)>0) then
       do icount = 1,abs(ntmqcd)
        write(unit=8,fmt="(a16,i6,2es17.10)") &
                                         "i,mtmqcd(i,:) = ", icount, &
                                                             mtmqcd(icount,:)
       enddo ! icount
       do icount = 1,abs(ntmqcd)
        write(unit=8,fmt="(a13,i6,4i10)") "i,src(i,:) = ", icount,src(icount,:)
       enddo ! icount
      endif

      write(unit=8, fmt="(a16,i10)") "ntmqcdloop =   ", ntmqcdloop
      if(abs(ntmqcdloop) > 0) then
       do icount = 1, abs(ntmqcdloop)
        write(unit=8, fmt="(a16,i6,2es17.10)") &
                     "i,mtmqcdloop(i,:) =", icount,mtmqcdloop(icount,:)
       enddo !icount
      endif
    

      write(unit=8,fmt="(a12,es17.10)")  "cSWtmqcd =  ", cSWtmqcd
      write(unit=8,fmt="(a12,i10)")      "inverter =  ", inverter
      write(unit=8,fmt="(a12,2i10)")     "nGMRES =    ", nGMRES
      write(unit=8,fmt="(a12,i10)")      "opSST =     ", opSST
      write(unit=8,fmt="(a12,i10)")      "tSST =      ", tSST
      write(unit=8,fmt="(a12,3i10)")     "pSST =      ", pSST
      write(unit=8,fmt="(a12,2i10)")     "nSST =      ", nSST
      if (abs(nSST(2))>0) then
       do icount = 1,abs(nSST(2))
        write(unit=8,fmt="(a16,i6,2es17.10)") &
                                       "i,mSST(i,:) =   ", icount,mSST(icount,:)
       enddo ! icount
      endif
      write(unit=8,fmt="(a12,es17.10)")  "cSWSST =    ", cSWSST
      write(unit=8,fmt="(a12,i10)")      "invSST =    ", invSST
      write(unit=8,fmt="(a12,i10)")      "numnoises = ", numnoises
      write(unit=8,fmt="(a12,2es17.10)") "kappaloop = ", kappaloop(1), kappaloop(2)
      write(unit=8,fmt="(a12,es17.10)")  "cSWloop =   ", cSWloop
      write(unit=8,fmt="(a12,i10)")      "docorr =    ", docorr
      write(unit=8,fmt="(a12,i10)")      "nkappacorr =", nkappacorr
      write(unit=8,fmt="(a12,i10)")      "ntmqcdcorr =", ntmqcdcorr
      write(unit=8,fmt="(a12,i10)")      "ntmqcdloop =", ntmqcdloop
      write(unit=8,fmt="(a12,i10)")      "nucleonparticle=", nucleonparticle
      write(unit=8,fmt="(a12,i10)")      "gitest =    ", gitest
      write(unit=8,fmt="(a12,i10)")      "nirloop=    ", nirloop
      write(unit=8,fmt="(a12)")          "            "
      write(unit=8,fmt="(a12)")          "OUTPUTS:    "
      write(unit=8,fmt="(a12)")          "            "
     close(unit=8,status="keep")
    endif

    ! This file will record the eigenvalues from gmres-dr or lan-dr
    if (myid==0) then
     open(unit=9,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",action="write",&
          form="formatted",status="new")
     close(unit=9,status="keep")
    endif
!      print *,"start"
    !Timing the whole run
     call cpu_time(time1)
!      print *,"end"
! Need to initialize the "sweepcount" to be passed into 
! twistdiscon. 
     sweepnum = 0

     if (newcalc==3) then
       sweeptotal = (nsweep - n0sweep)/n1sweep + 1 
     else
       sweeptotal = nsweep       
     endif ! newcalc

! Fix a few parameters.
    select case(newcalc)
     case(0)
      tadupdate = 0
      ifilesave = 4
     case(1)
      tadupdate = 0
      if (opSST==0) then
       ifilesave = 4
      else
       ifilesave = max(ifilesave,3)
      endif
     case(2)
      tadupdate = 0
      if (nkappa==0.and.ntmqcd==0.and.opSST==0) then
       ifilesave = 4
      elseif ((nkappa/=0.or.ntmqcd/=0).and.opSST==0) then
       ifilesave = max(ifilesave,2)
       if (ifilesave==3) ifilesave=4
      else
       ifilesave = max(ifilesave,1)
      endif
     case(3)
      if (tadupdate>0) tadpole=1
     case default
      open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
           ,form="formatted")
       write(unit=8,fmt=*) "improper value, newcalc = ", newcalc
      close(unit=8,status="keep")
      print *, "Calling MPI_ABORT near 'improper value, newcalc'"
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      stop
    end select
!-Set MRT to be MPI_REAL or MPI_DOUBLE_PRECISION.  It must match KR.
    if (KR==2.or.KR==8) then
     MRT = MPI_DOUBLE_PRECISION
    else
     MRT = MPI_REAL
    endif
!-Set MRT2 to be MPI_REAL or MPI_DOUBLE_PRECISION.  It must match KR2.
    if (KR2==2.or.KR2==8) then
     MRT2 = MPI_DOUBLE_PRECISION
    else
     MRT2 = MPI_REAL
    endif
! For GMRES, cgvectype is the "n" in GMRES(n).
! For BiCGstab, cgvectype=1 => first call to fermprop.
! For CGNE and QMR, cgvectype is irrelevant.
    if (inverter>3) then
     cgvectype = nGMRES(1)
    else
     cgvectype = 1
    endif
    if (invSST>3) then
     SSTvectype = nGMRES(1)
    else
     SSTvectype = 1
    endif
! For minimal residual, mrvectype=1 => b_e=0,  b_o/=0
! For minimal residual, mrvectype=2 => b_e/=0, b_o=0
! For minimal residual, mrvectype=3 => b_e/=0, b_o/=0 (not a point source)
    checksum1 = src(0,1) + src(0,2) + src(0,3) + src(0,4)
    checksum2 = (src(0,1) + src(0,2) + src(0,3) + src(0,4))/2
    if (checksum1==2*checksum2) then
     mrvectype = 2
    else
     mrvectype = 1
    endif
! Decide whether the clover coefficient is intended to be exactly zero or not.
    if (abs(cSWtmqcd)<0.0001) then
     itmflag = -1
    else
     itmflag = -2
    endif
    if (abs(cSWSST)<0.0001) then
     iSSTflag = -1
    else
     iSSTflag = -2
    endif
    if (abs(cSWloop)<0.0001) then
     iloopflag = -1
    else
     iloopflag = -2
    endif
! Local source, smeared source or both?  Define filename trailers.
    if (nsmear==0) then
     nLS = 1
     LSsymbol(1) = "L"
    elseif (nsmear>0) then
     nLS = 1
     LSsymbol(1) = "S"
    else
     nLS = 2
     LSsymbol(1) = "L"
     LSsymbol(2) = "S"
    endif




!if (myid==0) then 
!    print*, "At Line 636" 
!endif






! Setup the links.
!    print *,"strat"
    call initlat(myid,bc,nn,ldiv,nms,lv,lvbc,ib,lbd,iblv,vecbl,vecblinv)
!    print *,"end"
! Initiate the random number generator, and file which stores latest seeds.
    if (newcalc==3.or.numnoises>0.or.gitest==1) then
     lagfibfile = "cfgfiles/LAGFIBxxx.SAVE"
     write(unit=lagfibfile(16:18),fmt="(i3.3)") myid
     call initlagfib(myid,iseed)
if (.false.) then ! change true to false for generating lagfib files when newcalc = 3 TW                                     
     if (ncold==0.or.numnoises>0.or.gitest==1) then
      open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),action="read", &
           form="unformatted",status="old",position="rewind")
       read(unit=8) iwlagfib,jrlagfib,krlagfib
      close(unit=8,status="keep")
     endif
endif
    endif




! Set beta, anisotropies and tadpole factors.
! Helpful statements -WW
    if (newcalc==3.and.ncold==1) then
     u = 0.0_KR
     u(1,:,:,:,:) = 1.0_KR
     u(9,:,:,:,:) = 1.0_KR
     u(17,:,:,:,:) = 1.0_KR
     icfgsave = 0
    else




!if (myid==0) then 
!    print*, "At Line 680" 
!endif







     write(unit=cfgfile(16:18),fmt="(i3.3)") myid
     print *,"myid,cfgfile=",myid,cfgfile
     
     call GFread(u,trim(rwdir(myid+1))//trim(cfgfile),gaction,beta, &
                 alat,tad,icfgsave,myidin)
     
     if(myid.eq.0) print *,"icfgsave=",icfgsave
     call printlog("Just did GFread",myid,rwdir)
     if (myidin/=myid) then
      open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
           ,form="formatted")
       write(unit=8,fmt=*) "subroutine cfgsprops: myid =",myid," AND ",myidin
      close(unit=8,status="keep")
      call MPI_FINALIZE(ierr)
      stop
     endif
     if (newcalc<3) icfgsave=icfgsave-1 !...so the first cfg is not skipped.
    endif

! Determine the action's coefficients.
    call fixaction(gaction,beta,alat,tad,coact)

! Perform the desired number of gauge field updates.
    if (newcalc<3) then
     n0sweep = 1
     n1sweep = 1
    endif
    tadrun = 0
    icount = -n0sweep
    if (myid==0) timeinit = MPI_WTIME()
    firstcfg = 1

    if(myid==0) then
     open(unit=12,file=trim(rwdir(myid+1))//"DUBLIN.OUT",action="write",&
          form="formatted",status="new")
     close(unit=12,status="keep")
    endif ! myid

    if(myid==0) call printlog("Entering sweep loop",myid,rwdir)
    do isweep = 1,nsweep

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
            action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt=*)  "isweep = ",isweep
       close(unit=8,status="keep")
     endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     icount = icount + 1
     if (newcalc==3) call sweep(gaction,coact,nhit,u,myid,nn,ldiv,nms,lv,ib, &
                                lbd,iblv,MRT)
!-for newcalc=3, save selected configurations (and random number information).
!-for newcalc<3, read selected configurations.
     if (icount==0) then
     sweepnum = sweepnum + 1

      call printlog("Updated sweepnum",myid,rwdir)
      if(myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write",&
            form="formatted",status="old",position="append")
        write(unit=8, fmt="(a15,i5)") "sweepnum = ",sweepnum
       close(unit=8,status="keep")
      endif ! myid

! This is just purely infomational. Used as a check of the scalar form factor.
      if(myid==0) then
       open(unit=12,file=trim(rwdir(myid+1))//"DUBLIN.OUT",action="write",&
            form="formatted",status="old",position="append")
        write(unit=12, fmt="(a15,i5)") "sweepnum = ",sweepnum
       close(unit=12,status="keep")
      endif ! myid


      icount = -n1sweep
!-delete the configuration and propagator files from the previously saved sweep,
! if desired.
      
      if (ifilesave<4) then
       if (firstcfg==0) then
        if (ifilesave==0) then
         open(unit=9,file=trim(rwdir(myid+1))//trim(cfgfile), &
              action="write",status="old",form="unformatted")
         close(unit=9,status="delete")
        endif
        if (ifilesave<3.and.nkappa>0) then
         do iLS = 1,nLS
          do idsrc = 1,4
           do icsrc = 1,3
            trailer = ".xxk"//LSsymbol(iLS)
            write(unit=trailer(2:2),fmt="(i1.1)") icsrc
            write(unit=trailer(3:3),fmt="(i1.1)") idsrc
            do ikappa = 1,nkappa
             write(unit=trailer(4:4),fmt="(i1.1)") ikappa
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer
             open(unit=9,file=propfile, &
                  action="write",status="old",form="unformatted")
             close(unit=9,status="delete")
            enddo ! ikappa
           enddo ! icsrc
          enddo ! idsrc
         enddo ! iLS
        endif
        if (ifilesave<3.and.ntmqcd/=0) then
         do iLS = 1,nLS
          do idsrc = 1,4
           do icsrc = 1,3
            trailer = ".xxk"//LSsymbol(iLS)
            write(unit=trailer(2:2),fmt="(i1.1)") icsrc
            write(unit=trailer(3:3),fmt="(i1.1)") idsrc
            do ikappa = 1,abs(ntmqcd)
             write(unit=trailer(4:4),fmt="(i1.1)") ikappa
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmU"
             open(unit=9,file=propfile, &
                  action="write",status="old",form="unformatted")
             close(unit=9,status="delete")
             if (ntmqcd<0) then
              propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmD"
              open(unit=9,file=propfile, &
                   action="write",status="old",form="unformatted")
              close(unit=9,status="delete")
             endif
            enddo ! ikappa
           enddo ! icsrc
          enddo ! idsrc
         enddo ! iLS
        endif ! (ifilesave<3.and.ntmqcd/=0)
        if ((ifilesave<2.or.ifilesave==3).and.opSST/=0) then
         do iLS = 1,nLS
          do idsrc = 1,4
           do icsrc = 1,3
            trailer = ".xxk"//LSsymbol(iLS)
            write(unit=trailer(2:2),fmt="(i1.1)") icsrc
            write(unit=trailer(3:3),fmt="(i1.1)") idsrc
            if (opSST/=0) then
             do ikappa = 1,abs(nSST(1))
              do iSST = 1,abs(nSST(2))
               ibit = (ikappa-1)*abs(nSST(2)) + iSST
               write(unit=trailer(4:4),fmt="(i1.1)") ibit
               propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmUtmU"
               open(unit=9,file=propfile, &
                    action="write",status="old",form="unformatted")
               close(unit=9,status="delete")
               if (nSST(2)<0) then
                propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmUtmD"
                open(unit=9,file=propfile, &
                     action="write",status="old",form="unformatted")
                close(unit=9,status="delete")
               endif
               if (nSST(1)<0) then
                propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmDtmU"
                open(unit=9,file=propfile, &
                     action="write",status="old",form="unformatted")
                close(unit=9,status="delete")
                if (nSST(2)<0) then
                 propfile =trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmDtmD"
                 open(unit=9,file=propfile, &
                      action="write",status="old",form="unformatted")
                 close(unit=9,status="delete")
                endif
               endif
              enddo ! iSST
             enddo ! ikappa
            endif
           enddo ! icsrc
          enddo ! idsrc
         enddo ! iLS
        endif ! ((ifilesave<2.or.ifilesave==3).and.opSST/=0)
       else
        firstcfg = 0
       endif ! (firstcfg==0)
      endif ! (ifilesave<4)
!-read/write the current configuration files.
      icfgsave = icfgsave + 1

! If we are reading in configs then sweepnum will be the number of the
! config read in.
      sweepnum = icfgsave





!if (myid==0) then 
!    print*, "At Line 875" 
!endif










      write(unit=cfgfile(11:15),fmt="(i5.5)") icfgsave
      write(unit=cfgfile(16:18),fmt="(i3.3)") myid
      if (newcalc==3) then
       call GFwrite(u,trim(rwdir(myid+1))//trim(cfgfile),gaction,beta, &
                    alat,tad,icfgsave,myid)
       open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
            action="write",form="unformatted",status="replace")
        write(unit=8) iwlagfib,jrlagfib,krlagfib
       close(unit=8,status="keep")
     If(myid==0) then
      call printlog("Just did a write",myid,rwdir)
     endif ! myid
      else
       call GFread(u,trim(rwdir(myid+1))//trim(cfgfile),gaction,beta, &
                   alat,tad,icfgsave,myidin)
     If(myid==0) then
      call printlog("Just did a read",myid,rwdir)
     endif ! myid
      endif ! (newcalc==3)
!-perform a local gauge transformation if requested.
      if (gitest==1) then
       call gaugerot(u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
            action="write",form="unformatted",status="replace")
        write(unit=8) iwlagfib,jrlagfib,krlagfib
       close(unit=8,status="keep")
      endif ! (gitest==1)
!-tadpoles and average plaquette calculations.
      if (tadpole==1) then
       call aveplaq(plaq,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       if (tadlandau==1) then
        call tadLcalc(tadmu,u,alpha,omega,thetamax,landimp,alat,tad,myid, &
                      nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       else
        call tadPcalc(tadmu,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       endif
       if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
             action="write",form="formatted",status="old",position="append")
         write(unit=8,fmt="(i10,a13, es17.10)") icfgsave, " plaq:       ", plaq
         write(unit=8,fmt="(i10,a13,4es17.10)") icfgsave, " ux,uy,uz,ut:", tadmu
        close(unit=8,status="keep")
       endif
      endif ! (tadpole==1)
!-tadpole updates.


      if (tadupdate>0) then
       tadrun = tadrun + tadmu
       if (modulo(isweep,tadupdate)==0) then
        tad = 0.5_KR*( tad + tadrun/real(tadupdate,KR) )
        tadrun = 0
        call fixaction(gaction,beta,alat,tad,coact)
        if (myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
              action="write",form="formatted",status="old",position="append")
          write(unit=8,fmt="(a17,4es17.10)") "ADJUSTMENT: tad =",tad
         close(unit=8,status="keep")
        endif
       endif
      endif


!-Wilson loop calculations.
      if (nwilo/=0) &
       call somewloops(nwilo,wlat,Rmax,npaths,nfuzz,epsfuzz,u,myid,nn,ldiv, &
                       nms,lv,ib,lbd,iblv,icfgsave,rwdir,tad,cfgfile,MRT)


!-Wilson fermion propagator construction with or without multi-massing.

      if(nirloop==1) ir=1
      if(nirloop==2) ir=-1



      if (newcalc>1 .and. nkappa>0) then
!
! This is the paramter which allows one to avoid the M3R for the
! Wilson case inside fermprop. It directs the program instead to
! inverter=6 inverters (which do not multi-mass).
!
       if(wilsonmm==1) then

        do iLS = 1,nLS
         pion = 0.0_KR
         rho = 0.0_KR
         decay = 0.0_KR
         kaptm(2) = 0.0_KR
         trailer = ".xxk"//LSsymbol(iLS)
         do idsrc = 1,4
          write(unit=trailer(3:3),fmt="(i1.1)") idsrc
          do icsrc = 1,3
           write(unit=trailer(2:2),fmt="(i1.1)") icsrc
           call pointsource(be,bo,src(0,:),icsrc,idsrc,myid,iblv,vecblinv)
           if (LSsymbol(iLS)=="S") then
            call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear, &
                       bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            be(:,:,:,:,:) = xe(:,:,:,:,:,1)
            bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
           endif

           call fermpropR(rwdir,be,bo,nkappa,kappa,cSWtmqcd,inverter,coact,bc, &
                        resmax,itermin,omegaM3R,u,mrvectype,nGMRES(2),xe,xo, &
                        myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv, &
                        0,hoption,MRT,MRT2,ir)

          
           do ikappa = 1,nkappa
            write(unit=trailer(4:4),fmt="(i1.1)") ikappa
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer
            call FPwrite(xe(:,:,:,:,:,ikappa),xo(:,:,:,:,:,ikappa),vecblinv, &
                         propfile,kappa(ikappa),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                         icfgsave,myid)
            call pionrho(pion(:,ikappa,1),rho(:,ikappa,1),decay(:,ikappa,1), &
                      xe(:,:,:,:,:,ikappa),xo(:,:,:,:,:,ikappa),idsrc,ntmqcd,myid,MRT)
           enddo ! ikappa
          enddo ! icsrc
         enddo ! idsrc
         if(myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
              action="write",form="formatted",status="old",position="append")
         do ikappa = 1,nkappa
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,pion =  ", &
                                               ikappa," ",pion(:,ikappa,1)
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,rho =   ", &
                                               ikappa," ",rho(:,ikappa,1)
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,decay = ", &
                                               ikappa," ",decay(:,ikappa,1)
         enddo ! ikappa
         close(unit=8,status="keep")
         endif
        enddo ! iLS
   
       else !wilsonmm

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! output matrix values
if (.false.) then !changed to true TW 10/16/17
! Code to get out Heo, Hoe, and 1/k^2-HeoHoe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
    ! idag = 0 is for Hsingle/Hdouble to not be daggered
    idag = 0

    ! identify the correct number of 
    if (myid==0) write(*,*) "nvhalf = ", nvhalf
    if (myid==0) write(*,*) "numprocs= ", numprocs

    do prc=0,numprocs-1 
      do evodd=1,2
        do jj=1,nvhalf!ntotal
          do mm=1,8
            do ll=1,2
              do ii=1,5,2
                do kk=1,4
                  ! initialize vectors
                  evenVin=0.0_KR2
                  evenVout=0.0_KR2

                  oddVin=0.0_KR2
                  oddVout=0.0_KR2

                  ! since x is split in even/odd initialization is akward
                  ! multipy unity vector by M
                  if (evodd==1) then
                    if (myid == prc) then
                      evenVin(ii,jj,kk,ll,mm) = 1.0_KR2
                    endif
                    ! gonna do H_oe; must mulitply by xe
                    call Hsingle(oddVout,u,evenVin,idag,coact,bc,2_KI,vecbl,vecblinv, &
                             myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

                    ! gonna do HeoHoe multiplication for simplified matrix
                    !call Hdouble(evenVout,u,evenVin,idag,coact,kappa(1),bc,vecbl,vecblinv,myid,nn,ldiv, &
                    !        nms,lvbc,ib,lbd,iblv,MRT)

                    evenVin = 0.0_KR2
                    if (myid == prc) then
                      evenVin(ii,jj,kk,ll,mm) = 1.0_KR2
                    endif
                  else
                    if (myid==prc) then
                      oddVin(ii,jj,kk,ll,mm)  = 1.0_KR2
                    endif

                    ! gonna do H_eo; must mulitply by xo
                    call Hsingle(evenVout,u,oddVin,idag,coact,bc,1_KI,vecbl,vecblinv, &
                             myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

                    oddVin = 0.0_KR2
                    if (myid==prc) then
                      oddVin(ii,jj,kk,ll,mm)  = 1.0_KR2
                    endif
                  endif

                  ! convert from RL notation to (site,color,dirac) notation
                  ! This version of changevector is a direct copy of the one
                  ! vevbleft.f90. Done because of the difficulty of compilation.
                  call changevector(realz2noise,imagz2noise,evenVin,oddVin(:,:nvhalf,:,:,:),numprocs,MRT,myid)

                  if (myid==0) then
                    coll = 1 ! used to identify the correct colum of M being computed
                    ! look for where the 1 is to identify the column being computed
                    do iaa=1,nxyzt ! position
                      do ibb=1,3   ! color
                        do icc=1,4 ! dirac
                          rz2=realz2noise(iaa,ibb,icc)
                          iz2=imagz2noise(iaa,ibb,icc)
                          if (ABS(rz2) > .1) then
                            roww = coll
                          endif
                          coll=coll+1
                        enddo
                      enddo
                    enddo
                    coll = roww
                  endif

                  if (nps .ne. 1) then
                    call MPI_BCAST(coll,1,MRT2,0,MPI_COMM_WORLD,ierr)
                  endif

                  if (coll < nxyzt*3*4+1) then
                    ! initialize output for vector conversion
                    ! convert to readable notation
                    call changevector(realz2noise,imagz2noise,evenVout,oddVout(:,:nvhalf,:,:,:),numprocs,MRT,myid)

                    ! index that runs through the converted vectors
                    roww=1
                    ! fill M matrix
                    do iaa=1,nxyzt
                      do ibb=1,3
                        do icc=1,4
                          rz2=realz2noise(iaa,ibb,icc)
                          iz2=imagz2noise(iaa,ibb,icc)

                          ! Output Hoe/Heo in [row col real imag] format
                          if (myid ==0) then
                            if (ABS(rz2) > .000001 .or. ABS(iz2) > .000001) then
                              write(*,fmt='(I6,A,I6,A,F20.16,A,F20.16)') roww,' ',coll,' ',rz2,' ',iz2
                            endif
                          endif

                          roww=roww+1
                        enddo
                      enddo
                    enddo
                  endif

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

  if (myid==0) write(*,*) "poop "
  call MPI_FINALIZE(ierr)
  stop
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if (.false.) then        
        do iLS = 1,nLS
         pion = 0.0_KR
         rho = 0.0_KR
         decay = 0.0_KR
         kaptm(2) = 0.0_KR
         trailer = ".xxk"//LSsymbol(iLS)
         do ikappa=1,nkappa
          xe=0.0_KR2
          xo=0.0_KR2
          kaptm(1) = kappa(ikappa) 
          do idsrc = 1,4
           write(unit=trailer(3:3),fmt="(i1.1)") idsrc
           do icsrc = 1,3
            write(unit=trailer(2:2),fmt="(i1.1)") icsrc
            call setrhs(idsrc,icsrc) !solve using multiple right-hand sides set up
            call pointsource(be,bo,src(0,:),icsrc,idsrc,myid,iblv,vecblinv)
            if (LSsymbol(iLS)=="S") then
             call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear, &
                        bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             be(:,:,:,:,:) = xe(:,:,:,:,:,1)
             bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
            endif
!
! We can use either fermpropR or fermprop_wilson here. The fermpropR runs the
! Wilson fermions as a twsited mass pair with mu=0, whereas the fermprop_wilson
! is just a simplified version of the original fermprop.
!
            call fermpropR(rwdir,be,bo,itmflag,kaptm,cSWtmqcd,inverter,coact,bc, &
                        resmax,itermin,omegaM3R,u,cgvectype,nGMRES(2),xe,xo, &
                        myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv, &
                        1,hoption,MRT,MRT2,ir)

            !call fermprop_wilson(rwdir,be,bo,kaptm(1),inverter,coact,bc, &
            !            resmax,itermin,u,cgvectype,nGMRES(2),xe(:,:,:,:,:,1),xo(:,:,:,:,:,1), &
            !            myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv, &
            !            hoption,MRT,MRT2)
           

           write(unit=trailer(4:4),fmt="(i1.1)") ikappa
           propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer
           call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv, &
                        propfile,kappa(ikappa),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                        icfgsave,myid)
           call pionrho(pion(:,ikappa,1),rho(:,ikappa,1),decay(:,ikappa,1), &
                      xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),idsrc,ntmqcd,myid,MRT)
          
          enddo ! icsrc
         enddo ! idsrc
        enddo !ikappa
        if(myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
             action="write",form="formatted",status="old",position="append")
         do ikappa = 1,nkappa
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,pion =  ", &
                                               ikappa," ",pion(:,ikappa,1)
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,rho =   ", &
                                               ikappa," ",rho(:,ikappa,1)
          write(unit=8,fmt="(a15,i1,a1,500es17.10)") "ikappa,decay = ", &
                                               ikappa," ",decay(:,ikappa,1)
         enddo ! ikappa
        close(unit=8,status="keep")
        endif
       enddo ! iLS
endif
       endif !if(wilsonmm==1)
      endif ! (newcalc>1 .and. nkappa>0)

   
!-twisted mass Wilson/clover fermion propagator construction.
      if (newcalc>1 .and. ntmqcd/=0 .and. inverter/=6) then

       do iLS = 1,nLS
        pion = 0.0_KR
        rho = 0.0_KR
        decay = 0.0_KR
        trailer = ".xxk"//LSsymbol(iLS)
        do idsrc = 1,4
         write(unit=trailer(3:3),fmt="(i1.1)") idsrc
         do icsrc = 1,3
          write(unit=trailer(2:2),fmt="(i1.1)") icsrc
          do ikappa = 1,abs(ntmqcd)
!           print *,"ntmqcd=",ntmqcd           
           write(unit=trailer(4:4),fmt="(i1.1)") ikappa
           kaptm(1) = mtmqcd(ikappa,1)
           kaptm(2) = mtmqcd(ikappa,2)
           propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmU"
           call pointsource(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)

           if (LSsymbol(iLS)=="S") then
            call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear, &
                       bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            be(:,:,:,:,:) = xe(:,:,:,:,:,1)
            bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
           endif

           call fermprop(rwdir,be,bo,itmflag,kaptm,cSWtmqcd,inverter,coact,bc,&
                      resmax,itermin,omegaM3R,u,cgvectype,nGMRES(2),xe,xo, &
                      myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv, &
                      abs(ntmqcd),hoption,MRT,MRT2)

           call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv, &
                        propfile,kaptm(1),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                        icfgsave,myid)
           call pionrho(pion(:,ikappa,1),rho(:,ikappa,1),decay(:,ikappa,1), &
                        xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),idsrc,ntmqcd,myid,MRT)

           if (ntmqcd<0) then
            kaptm(2) = -1.0_KR*mtmqcd(ikappa,2)
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmD"
            call pointsource(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)
            if (LSsymbol(iLS)=="S") then
             call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear,&
                        bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             be(:,:,:,:,:) = xe(:,:,:,:,:,1)
             bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
            endif
 
            call setrhs(idsrc,icsrc)

            call fermprop(rwdir,be,bo,itmflag,kaptm,cSWtmqcd,inverter,coact, &
                          bc,resmax,itermin,omegaM3R,u,cgvectype,nGMRES(2), &
                          xe,xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl, &
                          vecblinv,ntmqcd,hoption,MRT,MRT2)
            call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv, &
                         propfile,kaptm(1),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                         icfgsave,myid)
            call pionrho(pion(:,ikappa,2),rho(:,ikappa,2),decay(:,ikappa,2), &
                         xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),idsrc,ntmqcd,myid,MRT)
           endif
          enddo ! ikappa
         enddo ! icsrc
        enddo ! idsrc

        if (myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
              action="write",form="formatted",status="old",position="append")
          do ikappa = 1,abs(ntmqcd)
           write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,pion [tmU] =  ", &
                                                ikappa," ",pion(:,ikappa,1)
           write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,rho [tmU] =   ", &
                                                ikappa," ",rho(:,ikappa,1)
           write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,decay [tmU] = ", &
                                                ikappa," ",decay(:,ikappa,1)
           if (ntmqcd<0) then
            write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,pion [tmD] =  ",&
                                                 ikappa," ",pion(:,ikappa,2)
            write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,rho [tmD] =   ",&
                                                 ikappa," ",rho(:,ikappa,2)
            write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,decay [tmD] = ",&
                                                 ikappa," ",decay(:,ikappa,2)
           endif
          enddo ! ikappa
         close(unit=8,status="keep")
        endif
       enddo ! iLS
     ! enddo ! irhs
      endif ! (newcalc>1 .and. ntmqcd/=0 .and. inverter/=6)

!-twisted mass , single-mass Wilson fermion propagator construction.
! inversion with multiple RHS

    if (newcalc>1 .and. ntmqcd/=0 .and. inverter==6) then

    !call printlog("WARNING :: NOT DOING PROPAGATOR FOR THIS CONFIGURATION!!!",myid,rwdir)
    !if (myid==0) then
    !   open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
    !        action="write",form="formatted",status="old",position="append")
    !   write(unit=8,fmt="(a21,1i3)") "Configuration Number=", icfgsave
    !   close(unit=8,status="keep")
    !endif ! myid

! If logic == .false. , then you are skipping the nucleon propagator part.

    if (.true.) then ! ABDOU ~ If you don't need to do props, turn to false.
! right here Carl

!  call cpu_time(time1)

    do ikappa = 1,abs(ntmqcd)
      do iLS = 1,nLS
        pion = 0.0_KR
        rho = 0.0_KR
        decay = 0.0_KR
        trailer = ".xxk"//LSsymbol(iLS)
         do idsrc = 1,4
           write(unit=trailer(3:3),fmt="(i1.1)") idsrc
           do icsrc = 1,3
             write(unit=trailer(2:2),fmt="(i1.1)") icsrc
             write(unit=trailer(4:4),fmt="(i1.1)") ikappa
             kaptm(1) = mtmqcd(ikappa,1)
             kaptm(2) = mtmqcd(ikappa,2)

             call pointsource(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)

! zsmear is a zero momentum source. It projects out the scalar
! signal for the Dublin confrence. 7-18-05 

!            call zsmear(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)
!            if (myid==9) then
!               call check4ones(be,bo,icsrc,idsrc)
!            endif !  myid
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  stop
          

             call gamma5vector(be,abs(ntmqcd),myid)
             call gamma5vector(bo,abs(ntmqcd),myid)

             if (LSsymbol(iLS)=="S") then
               call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear, &
                          bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
               be(:,:,:,:,:) = xe(:,:,:,:,:,1)
               bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
             endif

             
             call setrhs(idsrc,icsrc)
               if (myid==0) then
                  print *,"idsrc+= ,icsrc=+",idsrc,icsrc
               endif

             call fermprop(rwdir,be,bo,itmflag,kaptm,cSWtmqcd,inverter,coact,bc,&
                           resmax,itermin,omegaM3R,u,cgvectype,nGMRES(2),xe, &
                           xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv, &
                           abs(ntmqcd),hoption,MRT,MRT2)


             write(unit=trailer(4:4),fmt="(i1.1)") ikappa
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmU"
             call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv, &
                          propfile,kaptm(1),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                          icfgsave,myid)

             call pionrho(pion(:,ikappa,1),rho(:,ikappa,1),decay(:,ikappa,1), &
                          xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),idsrc,ntmqcd,myid,MRT)
              
          enddo ! icsrc
        enddo ! idsrc
         
! Do tmD sperately for single-mass, mulitple RHS source vectors

        do idsrc = 1,4
         write(unit=trailer(3:3),fmt="(i1.1)") idsrc
         do icsrc = 1,3
          write(unit=trailer(2:2),fmt="(i1.1)") icsrc
           if (ntmqcd<0) then
            kaptm(2) = -1.0_KR*mtmqcd(ikappa,2)

             call pointsource(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)

! zsmear is a zero momentum source. It projects out the scalar
! signal for the Dublin confrence. 7-18-05

!            call zsmear(be,bo,src(ikappa,:),icsrc,idsrc,myid,iblv,vecblinv)


             call gamma5vector(be,ntmqcd,myid)
             call gamma5vector(bo,ntmqcd,myid)

            if (LSsymbol(iLS)=="S") then
             call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear,&
                        bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             be(:,:,:,:,:) = xe(:,:,:,:,:,1)
             bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
            endif ! (LSymbol=="S")
 
            call setrhs(idsrc,icsrc)
               if (myid==0) then
                  print *,"idsrc= ,icsrc=",idsrc,icsrc
               endif
            call fermprop(rwdir,be,bo,itmflag,kaptm,cSWtmqcd,inverter,coact, &
                          bc,resmax,itermin,omegaM3R,u,cgvectype,nGMRES(2), &
                          xe,xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl,vecblinv,&
                          ntmqcd,hoption,MRT,MRT2)


               write(unit=trailer(4:4),fmt="(i1.1)") ikappa
               propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer
               propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmD"
               call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv, &
                            propfile,kaptm(1),kaptm(2),cSWtmqcd,icsrc,idsrc, &
                            icfgsave,myid)
               call pionrho(pion(:,ikappa,2),rho(:,ikappa,2),decay(:,ikappa,2), &
                            xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),idsrc,ntmqcd,myid,MRT)
           endif ! (ntmqcd<0)
         enddo ! icsrc
        enddo ! idsrc
!     enddo ! ikappa
      
        if (myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
           action="write",form="formatted",status="old",position="append")
            write(unit=8,fmt="(a29,i1,a1,500es17.10)") "ikappa,pion [tmU] =  ", &
                                                        ikappa," ",pion(:,ikappa,1)
            write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,rho [tmU] =   ", &
                                                        ikappa," ",rho(:,ikappa,1)
            write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,decay [tmU] = ", &
                                                        ikappa," ",decay(:,ikappa,1)
            if (ntmqcd<0) then
              write(unit=8,fmt="(a29,i1,a1,500es17.10)") "ikappa,pion [tmD] =  ",&
                                                          ikappa," ",pion(:,ikappa,2)
              write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,rho [tmD] =   ",&
                                                          ikappa," ",rho(:,ikappa,2)
              write(unit=8,fmt="(a21,i1,a1,500es17.10)") "ikappa,decay [tmD] = ",&
                                                          ikappa," ",decay(:,ikappa,2)
            endif ! (ntmqcd<0)
         close(unit=8,status="keep")
        endif ! myid
       enddo ! iLS

       enddo ! ikappa

   endif ! .true. !!!! ****DEAN CTB*****

!    call cpu_time(time2)
!    if (myid==0) print *, " time to do nucleon propagators=", time2 - time1  

      endif ! (newcalc>1 .and. ntmqcd/=0 .and. inverter==6)

!-connected current insertions: sequential source technique.
      if (newcalc>0 .and. opSST/=0) then
       do iLS = 1,nLS
        pionSST = 0.0_KR
        trailer = ".xxk"//LSsymbol(iLS)
        do idsrc = 1,4
         write(unit=trailer(3:3),fmt="(i1.1)") idsrc
         do icsrc = 1,3
          write(unit=trailer(2:2),fmt="(i1.1)") icsrc
          do ikappa = 1,abs(nSST(1))
           write(unit=trailer(4:4),fmt="(i1.1)") ikappa
           propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmU"
           call FPread(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                       kapin,muin,cSWin,icsrc,idsrc,icfgin,myidin)
           call SSTsource(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),src(0,:),opSST,tSST, &
                          pSST,be,bo,bc,myid,iblv,vecblinv)
           if (opSST<0) then
            call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear, &
                       bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            be(:,:,:,:,:) = xe(:,:,:,:,:,1)
            bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
           endif
           do iSST = 1,abs(nSST(2))
            ibit = (ikappa-1)*abs(nSST(2)) + iSST
            write(unit=trailer(4:4),fmt="(i1.1)") ibit
            kaptm(1) = mSST(ikappa,1)
            kaptm(2) = mSST(ikappa,2)
            call fermprop(rwdir,be,bo,iSSTflag,kaptm,cSWSST,invSST,coact, &
                          bc,resmax,itermin,omegaM3R,u,SSTvectype,nGMRES(2), &
                          xe,xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl, &
                          vecblinv,0,hoption,MRT,MRT2)
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmUtmU"
            call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                         kaptm(1),kaptm(2),cSWSST,icsrc,idsrc,icfgsave,myid)
            call SSTpion(pionSST(ibit,1),xe(:,:,:,:,:,1),xo(:,:,:,:,:,1), &
                         src(0,:),icsrc,idsrc,myid,iblv,vecblinv,MRT)
            if (nSST(2)<0) then
             kaptm(2) = -1.0_KR*mSST(ikappa,2)
             call fermprop(rwdir,be,bo,iSSTflag,kaptm,cSWSST,invSST,coact,bc, &
                           resmax,itermin,omegaM3R,u,SSTvectype,nGMRES(2),xe, &
                           xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl, &
                           vecblinv,0,hoption,MRT,MRT2)
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmUtmD"
             call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                          kaptm(1),kaptm(2),cSWSST,icsrc,idsrc,icfgsave,myid)
             call SSTpion(pionSST(ibit,2),xe(:,:,:,:,:,1),xo(:,:,:,:,:,1), &
                          src(0,:),icsrc,idsrc,myid,iblv,vecblinv,MRT)
            endif
           enddo ! iSST
           if (nSST(1)<0) then
            write(unit=trailer(4:4),fmt="(i1.1)") ikappa
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmD"
            call FPread(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                        kapin,muin,cSWin,icsrc,idsrc,icfgin,myidin)
            call SSTsource(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),src(0,:),opSST, &
                           tSST,pSST,be,bo,bc,myid,iblv,vecblinv)
            if (opSST<0) then
             call smear(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),u,be,bo,asmear,nsmear,&
                        bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             be(:,:,:,:,:) = xe(:,:,:,:,:,1)
             bo(:,:,:,:,:) = xo(:,:,:,:,:,1)
            endif
            do iSST = 1,abs(nSST(2))
             ibit = (ikappa-1)*abs(nSST(2)) + iSST
             write(unit=trailer(4:4),fmt="(i1.1)") ibit
             kaptm(1) = mSST(ikappa,1)
             kaptm(2) = -1.0_KR*mSST(ikappa,2)
             call fermprop(rwdir,be,bo,iSSTflag,kaptm,cSWSST,invSST,coact,bc, &
                           resmax,itermin,omegaM3R,u,SSTvectype,nGMRES(2),xe, &
                           xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd,iblv,vecbl, &
                           vecblinv,0,hoption,MRT,MRT2)
             propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmDtmD"
             call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                          kaptm(1),kaptm(2),cSWSST,icsrc,idsrc,icfgsave,myid)
             call SSTpion(pionSST(ibit,3),xe(:,:,:,:,:,1),xo(:,:,:,:,:,1), &
                          src(0,:),icsrc,idsrc,myid,iblv,vecblinv,MRT)
             if (nSST(2)<0) then
              kaptm(2) = mSST(ikappa,2)
              call fermprop(rwdir,be,bo,iSSTflag,kaptm,cSWSST,invSST,coact, &
                            bc,resmax,itermin,omegaM3R,u,SSTvectype, &
                            nGMRES(2),xe,xo,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
                            iblv,vecbl,vecblinv,0,hoption,MRT,MRT2)
              propfile = trim(rwdir(myid+1))//trim(cfgfile)//trailer//"tmDtmU"
              call FPwrite(xe(:,:,:,:,:,1),xo(:,:,:,:,:,1),vecblinv,propfile, &
                           kaptm(1),kaptm(2),cSWSST,icsrc,idsrc,icfgsave,myid)
              call SSTpion(pionSST(ibit,4),xe(:,:,:,:,:,1),xo(:,:,:,:,:,1), &
                           src(0,:),icsrc,idsrc,myid,iblv,vecblinv,MRT)
             endif
            enddo ! iSST
           endif
          enddo ! ikappa
         enddo ! icsrc
        enddo ! idsrc
        if (myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
              action="write",form="formatted",status="old",position="append")
          do ikappa = 1,abs(nSST(1))
           do iSST = 1,abs(nSST(2))
            ibit = (ikappa-1)*abs(nSST(2)) + iSST
            write(unit=8,fmt="(a37,i2,i2,a2,a1,es17.10)") &
                  "kappa1,kappa2,smear,pion [tmUtmU] =  ",ikappa,iSST, &
                  LSsymbol(iLS)," ",pionSST(ibit,1)
            if (nSST(2)<0) then
             write(unit=8,fmt="(a37,i2,i2,a2,a1,es17.10)") &
                   "kappa1,kappa2,smear,pion [tmUtmD] =  ",ikappa,iSST, &
                   LSsymbol(iLS)," ",pionSST(ibit,2)
            endif
            if (nSST(1)<0) then
             write(unit=8,fmt="(a37,i2,i2,a2,a1,es17.10)") &
                   "kappa1,kappa2,smear,pion [tmDtmD] =  ",ikappa,iSST, &
                   LSsymbol(iLS)," ",pionSST(ibit,3)
             if (nSST(2)<0) then
              write(unit=8,fmt="(a37,i2,i2,a2,a1,es17.10)") &
                    "kappa1,kappa2,smear,pion [tmDtmU] =  ",ikappa,iSST, &
                    LSsymbol(iLS)," ",pionSST(ibit,4)
             endif
            endif
           enddo ! iSST
          enddo ! ikappa
         close(unit=8,status="keep")
        endif
       enddo ! iLS
      endif ! (newcalc>0 .and. opSST/=0)

!-disconnected loop construction: stochastic estimator technique (STE).

! NOTE: ~ here we have the option to do the (STE) for regular Wilson
!         case when inverter < 6 and for the twisted Wilson case
!         when inverter==6.

! WARNING, WARNING, WARNING, WARNING
! I need to change the inverter logic such that I can execute
! either discon OR twistdiscon independant of inverter.

      if (inverter/=6) then
        if (numnoises>0) then
          call discon(numnoises,be,bo,rwdir,iloopflag,kappaloop(1),cSWloop,coact, &
                      bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv, &
                      ib,lbd,iblv,vecbl,vecblinv,hoption,MRT,MRT2)
          open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
               action="write",form="unformatted",status="replace")
            write(unit=8) iwlagfib,jrlagfib,krlagfib
          close(unit=8,status="keep")
        endif ! (numnoises>0)
!-correlation function construction.
       if (docorr>0) then
         if (modulo(docorr,3)==1) then
           nsmrsnk = 0
         elseif (modulo(docorr,3)==2) then
           nsmrsnk = nsmear
         else
           nsmrsnk = -nsmear
         endif
         if (docorr<4 .or. docorr>6) &
           call mesoncorr(nkappacorr,ntmqcdcorr,nSSTcorr,src(0,:),nsmear,nsmrsnk, &
                          asmear,bc,rwdir,cfgfile,u,vecbl,vecblinv,myid,nn,ldiv, &
                          nms,lvbc,ib,lbd,iblv,MRT,be,bo,xe(:,:,:,:,:,1), &
                          xo(:,:,:,:,:,1),icfgsave)
         if (docorr>3) &
           call corr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                     MRT)
       endif ! (docorr>0)
!-calculate the total time for this number of sweeps.

           if (myid==0) then
             timefin = MPI_WTIME()
             sweeptime = timefin - timeinit
             open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
             action="write",form="formatted",status="old",position="append")
               write(unit=8,fmt="(a8,i10,a20,es17.10,a8)") &
                     "TIMING: ",isweep," sweeps performed in",sweeptime," seconds"
               close(unit=8,status="keep")
           endif ! myid
!-end the set of calculations for this configuration.

         endif ! (inverter/=6)

! Use the stochastic estimator technique to do the loop construcion for tmU and
! tmD.

! Recall mtmqcd(ikappa,1) is the kappa value for the parity maximal twist of a given mass.
! Recall mtmqcd(ikappa,2) is the mu value for the parity maximal twist of a given mass.

      if (inverter==6) then  
        if (ntmqcdloop /= 0) then
          mudelta = 0.0_KR

! I am going to modify the gauge fields here to do a test of the inversion. -WW

!         u=0.0_KR
!         u(1,:,3,:,:)=1.0_KR 
!         u(9,:,3,:,:)=1.0_KR
!         u(17,:,3,:,:)=1.0_KR
!         u(1,:,4,:,:)=1.0_KR
!         u(9,:,4,:,:)=1.0_KR
!         u(17,:,4,:,:)=1.0_KR
!         u(:,:,1,:,:) = 0.0_KR
!         u(:,:,2,:,:) = 0.0_KR


      do ikappa = 1,abs(ntmqcdloop)

        if(nirloop.eq.1) then
          irloop1=1
          irloop2=1
        else if(nirloop.eq.2) then
          irloop1=2
          irloop2=2
        else if(nirloop.eq.3) then
          irloop1=1
          irloop2=2
        endif

        do irloop=irloop1,irloop2
         if(irloop.eq.1) ir=1
         if(irloop.eq.2) ir=-1

              !if(nirloop.eq.3.and.ir.eq.-1) then
              !  open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
              !       action="read",form="unformatted",status="old",position="rewind")
              !       read(unit=8) iwlagfib,jrlagfib,krlagfib
              !  close(unit=8,status="keep")
              !endif

              if(mod(ikappa,2).eq.0) then
                open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
                     action="read",form="unformatted",status="old",position="rewind")
                     read(unit=8) iwlagfib,jrlagfib,krlagfib
                close(unit=8,status="keep")
              endif


             if (numnoises>0) then
                ntm = abs(ntmqcdloop)
                if(myid==0) call printlog("Doing twisted case", myid,rwdir)
                mudelta = atan(2.0_KR*mtmqcdloop(ikappa,1)*mtmqcdloop(ikappa,2)) 
                call twistdiscon(numnoises,be,bo,rwdir,iloopflag,mtmqcdloop(ikappa,1),mtmqcdloop(ikappa,2),&
                                 cSWloop,coact,bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv, &
                                 ib,lbd,iblv,vecbl,vecblinv,hoption,MRT,MRT2,numprocs,mudelta,&
                                 ikappa,ntm,nGMRES(2),sweepnum,sweeptotal,icfgsave,ir)






!BS                
! Note here is where I can comment out the lagfib updates. -WW
       !      if(.true.) then
              if(.not.(nirloop.eq.3.and.ir.eq.1)) then
                open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
                     action="write",form="unformatted",status="replace")
                     write(unit=8) iwlagfib,jrlagfib,krlagfib
                close(unit=8,status="keep")
              endif


!             if(.true.) then
 !             if(.not.(nirloop.eq.3.and.ir.eq.1)) then
  !              open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
   !                  action="write",form="unformatted",status="replace")
    !                 write(unit=8) iwlagfib,jrlagfib,krlagfib
     !           close(unit=8,status="keep")
      !        endif














!BS


              if(mod(ikappa,2).eq.0) then
                open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
                     action="write",form="unformatted",status="replace")
                     write(unit=8) iwlagfib,jrlagfib,krlagfib
                close(unit=8,status="keep")
              endif

             if (ntmqcdloop<0) then
                ntm = ntmqcdloop
! NOTE~ mudelta changes sign when we do the "down" part of the maximal twist.
                call twistdiscon(numnoises,be,bo,rwdir,iloopflag,mtmqcdloop(ikappa,1),-mtmqcdloop(ikappa,2),&
                                 cSWloop,coact,bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv, &
                                 ib,lbd,iblv,vecbl,vecblinv,hoption,MRT,MRT2,numprocs,-mudelta,&
                                 ikappa,ntm,nGMRES(2),sweepnum,sweeptotal,icfgsave,ir)
                open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
                action="write",form="unformatted",status="replace")
                write(unit=8) iwlagfib,jrlagfib,krlagfib
                close(unit=8,status="keep")
             endif ! (ntmqcdloop<0)
           endif ! (numnoises>0) 
          enddo ! irloop
         enddo ! ikappa
       else
! Use the twisted routine in the Wilson case.

! ABDOU ~ take out the ikappa loop if you are multi-massing. This
!         will affect the kappaloop variable...probably need to change this in here
!         and in the file ~/qqcd/user/cfgspropsmain.f90

         if (numnoises>0) then
           do ikappa = 1,1 !2 !****CTB**** 
            !call printlog("WARNING :: Hey Abdou, we are only doing one loop",myid,rwdir) 
            !do ikappa = 1,1
               mudelta  = 0.0_KR
               muWilson = 0.0_KR
               ntm = abs(ntmqcdloop)
               call printlog("Doing Wilson case", myid,rwdir)

           !    ! twisted mass Wilson version by Dean
           !    call twistdiscon(numnoises,be,bo,rwdir,iloopflag,kappaloop(ikappa),muWilson,&
           !                     cSWloop,coact,bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv, &
           !                     ib,lbd,iblv,vecbl,vecblinv,hoption,MRT,MRT2,numprocs,mudelta,&
           !                     ikappa,ntm,nGMRES(2),sweepnum,sweeptotal,icfgsave,1)

               ! wilson version by Victor
             !  call eigdiscon(u,kappa,numnoises,coact,bc,vecbl,vecblinv,myid,ntm,nn,ldiv, &
             !                 nms,lvbc,lv,ib,lbd,iblv,rwdir,MRT,MRT2,numprocs)


! BS Just to make sure testFUNC is not called while generating configurations
  if (newcalc /= 3) then
!print *,myid,": calling testFUNC"
               call testFUNC(u,kappa,numnoises,coact,bc,vecbl,vecblinv,myid,ntm,nn,ldiv, &
                              nms,lvbc,lv,ib,lbd,iblv,rwdir,MRT,MRT2,numprocs) !added inverter TW 11/17/15

end if
               ! wilson perturb version by R. Lewis
 !              call discon(numnoises,be,bo,rwdir,iloopflag,kappa(1),cSWloop,coact,bc, &
 !                          resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
 !                          iblv,vecbl,vecblinv,hoption,MRT,MRT2)
            enddo ! ikappa
         endif ! check numnoises
       endif ! check Wilson or TWmass





      !call printlog("WARNING :: NOT WRITING LOG FILES!!!",myid,rwdir)
      if (.true.) then



!-correlation function construction for twisted mass doublet value; tmU and tmD.
          if (docorr>0) then
            if (modulo(docorr,3)==1) then
              nsmrsnk = 0
            elseif (modulo(docorr,3)==2) then
              nsmrsnk = nsmear
            else
              nsmrsnk = -nsmear
            endif
            if (docorr<4 .or. docorr>6) &

! Determine each correlator for each solution coming from the previously solved
! system of equations. 
         
              call mesoncorr(nkappacorr,ntmqcdcorr,nSSTcorr,src(0,:),nsmear,nsmrsnk, &
                             asmear,bc,rwdir,cfgfile,u,vecbl,vecblinv,myid,nn,ldiv, &
                             nms,lvbc,ib,lbd,iblv,MRT,be,bo,xe(:,:,:,:,:,1), &
                             xo(:,:,:,:,:,1),icfgsave)

! The following comment is for subroutine nucleoncorr:
! This structure corresponds to the particle nucleon you wish to create. If you would like
! the proton (nucleonparticle=1), the neutron (nucleonparticle=2), and for both (nucleonparticle=3)

            if (docorr>3) then
             if(nucleonparticle==1) then ! proton
               call nucleoncorr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                                1,MRT,icfgsave)
             elseif(nucleonparticle==2) then ! neutron
               call nucleoncorr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                                0,MRT,icfgsave)
             elseif(nucleonparticle==3) then ! both proton and neutron
                call nucleoncorr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                                 1,MRT,icfgsave)
                call nucleoncorr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                                 0,MRT,icfgsave)
             else ! (Wilson case)
              call corr(nkappacorr,ntmqcdcorr,nsmear,src(0,:),bc,myid,rwdir,cfgfile,&
                        MRT)
             endif ! nucleonparticle
            endif ! docorr>3
          endif ! (docorr>0)
!-calculate the total time for this number of sweeps.
            if (myid==0) then
              timefin = MPI_WTIME()
              sweeptime = timefin - timeinit
              open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
              action="write",form="formatted",status="old",position="append")
                write(unit=8,fmt="(a8,i10,a9,i10,a20,es17.10,a8)") &
                      "TIMING: ",isweep,"for mass ", ikappa, &
                               " sweeps performed in",sweeptime," seconds"
              close(unit=8,status="keep")
            endif ! myid
!-end the set of calculations for this configuration.
 



       endif ! .false !!!! *****DEAN CTB*****





      endif ! (inverter==6)

     endif ! (icount==0)
    enddo ! isweep

! Compute disconnected loops before and after a local gauge transformation.
    if (gitest==1 .and. numnoises>0 .and. inverter/=6) then
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
        action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt="(a65)") &
             "BLGT: disconnected loops before a local gauge transformation"
      close(unit=8,status="keep")
     endif ! myid
     idebug = 0
     call z2source(be,bo,nvhalf)
write(*,*) 'XB'
     call discon(idebug,be,bo,rwdir,iloopflag,kappaloop(1),cSWloop,coact,bc, &
                 resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
                 iblv,vecbl,vecblinv,hoption,MRT,MRT2)
     open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
          action="write",form="unformatted",status="replace")
      write(unit=8) iwlagfib,jrlagfib,krlagfib
     close(unit=8,status="keep")
!-perform the gauge transformation...
     call qqcdrot(u,be,bo,myid,nn,ldiv,nms,lv,ib,lbd,iblv,vecbl,MRT)
!-after the gauge transformation.
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
        action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt="(a59)") &
             "ALGT: disconnected loops after a local gauge transformation"
      close(unit=8,status="keep")
     endif ! myid
     idebug = 0
write(*,*) 'XC'
     call discon(idebug,be,bo,rwdir,iloopflag,kappaloop(1),cSWloop,coact,bc, &
                 resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
                 iblv,vecbl,vecblinv,hoption,MRT,MRT2)
     open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
          action="write",form="unformatted",status="replace")
      write(unit=8) iwlagfib,jrlagfib,krlagfib
     close(unit=8,status="keep")

! Compute disconnected loops before and after a local gauge transformation
!         for twisted disconnected loops
   else if (inverter.eq.0) then! (gitest==1.and.numnoises>0.and.inverter/=6)
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
        action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt="(a60)") &
             "BLGT: twisted disconnected loops before a local gauge transformation"
      close(unit=8,status="keep")
     endif
     idebug = 0
     call twistedvalues
     do ishift = 1,nshifts
       call z2source(be,bo,nvhalf)
       call twistdiscon(idebug,be,bo,rwdir,iloopflag,mtmqcd(ishift,1),mtmqcd(ishift,2),&
                   cSWloop,coact,bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
                   iblv,vecbl,vecblinv,hoption,MRT,MRT2,numprocs,mudelta,&
                   ishift,ntmqcd,nGMRES(2),sweepnum,sweeptotal,icfgsave,1)
       open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
            action="write",form="unformatted",status="replace")
      write(unit=8) iwlagfib,jrlagfib,krlagfib
     close(unit=8,status="keep")
!-perform the gauge transformation...
     call qqcdrot(u,be,bo,myid,nn,ldiv,nms,lv,ib,lbd,iblv,vecbl,MRT)
!-after the gauge transformation.
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
        action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt="(a59)") &
             "ALGT: disconnected loops after a local gauge transformation"
      close(unit=8,status="keep")
     endif
     idebug = 0
     call twistdiscon(idebug,be,bo,rwdir,iloopflag,kappa(ishift),mutemp(ishift),&
                 cSWloop,coact,bc,resmax,itermin,cgvectype,u,myid,nn,ldiv,nms,lvbc,lv,ib,lbd, &
                 iblv,vecbl,vecblinv,hoption,MRT,MRT2,numprocs,mudelta,&
                 ishift,ntmqcd,nGMRES(2),sweepnum,sweeptotal,icfgsave,1)
     open(unit=8,file=trim(rwdir(myid+1))//trim(lagfibfile),      &
          action="write",form="unformatted",status="replace")
      write(unit=8) iwlagfib,jrlagfib,krlagfib
     close(unit=8,status="keep")
   enddo ! ishift
   endif ! inverter/=6

!if (myid ==0) then
!print *,"calling printH"
!endif
!     call printH(rwdir,u,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)



    ! Added status arguments -PL 7-15-22
    deallocate(u,stat=statu)
    deallocate(xe,stat=statxe)
    deallocate(xo,stat=statxo)
    deallocate(be,stat=statbe)
    deallocate(bo,stat=statbo)
    deallocate(btemp,stat=statbt)
    deallocate(xtemp,stat=statxt)

    call cpu_time(time2)


     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG", &
        action="write",form="formatted",status="old",position="append")
       write(unit=8,fmt="(a50,d20.10)") &
             "Time for the entire run(sec.)=",time2-time1
      close(unit=8,status="keep")
     endif


! Conclude MPI.
   if (nps/=1) call MPI_FINALIZE(ierr)
!BS    call MPI_FINALIZE(ierr)

 end subroutine generator

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine sweep(gaction,coact,nhit,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)

    integer(kind=KI), intent(in),    dimension(:)         :: gaction
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in)                          :: nhit, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                      :: ieo, ibl, mu
    real(kind=KR),   dimension(18,nvhalf) :: c

! Perform the pseudo-heatbath updates.  (c contains the staples.)
     do ieo = 1,2
      do ibl = 1,16
       do mu = 1,4
        call staple(mu,ieo,ibl,coact,gaction,c,u,myid,nn,ldiv,nms,lv,ib,lbd, &
                    iblv,MRT)
        call phbup(mu,ieo,ibl,c,u,nhit)
       enddo ! mu
      enddo ! ibl
     enddo ! ieo

 end subroutine sweep

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine somewloops(nwilo,wlat,Rmax,npaths,nfuzz,epsfuzz,u,myid,nn,ldiv, &
                       nms,lv,ib,lbd,iblv,icfgsave,rwdir,tad,cfgfile,MRT)
! Compute, and save to a file, the Wilson loops requested.

    integer(kind=KI), intent(in) :: nwilo, npaths, nfuzz, myid, icfgsave, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: wlat, nn, iblv
    real(kind=KR),    intent(in),    dimension(:)         :: Rmax, tad
    real(kind=KR),    intent(in)                          :: epsfuzz
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    character(len=*), intent(in),    dimension(:)         :: rwdir
    character(len=*), intent(in)                          :: cfgfile

    integer(kind=KI) :: Rx, Ry, Rz, Rsum, it, nperms, iperms, fperms, &
                        npermsbit, xxdir, yydir, zzdir, ttdir, iwilo, &
                        myidin, icfgsavein
    real(kind=KR)    :: tadtsq, tadxsq, tadysq, tadzsq, tadx2Rx, Rmaxrnd, beta
    real(kind=KR),    dimension(4)                        :: alat, tadin
    integer(kind=KI), dimension(4)                        :: gaction
    integer(kind=KI), parameter                      :: iuser=10000, juser=100
    integer(kind=KI), dimension(juser)                    :: pathi, pathf
    real(kind=KR),    dimension(max(nx,ny,nz,nt))         :: wloop, wloopbit
! USER: For large Wilson loops, increase the value of iuser or juser.
    integer(kind=KI), dimension(iuser,juser)              :: permvecs

! Consider each orientation of the lattice.
    do iwilo = 1,abs(nwilo)
     xxdir = wlat(iwilo,1)
     yydir = wlat(iwilo,2)
     zzdir = wlat(iwilo,3)
     ttdir = wlat(iwilo,4)
     Rmaxrnd = Rmax(iwilo) + 0.0001_KR

! Fuzz the spatial links, if nfuzz>0.
     if (nfuzz>0) call fuzz(ttdir,u,nfuzz,epsfuzz,myid,nn,ldiv,nms,lv,ib, &
                            lbd,iblv,MRT)

! Construct the Wilson loops that exist in a single "spatial" direction.
     if (nwilo<0) then
      tadtsq = 1.0_KR
      tadxsq = 1.0_KR
     else
      tadtsq = tad(ttdir)**2
      tadxsq = tad(xxdir)**2
     endif
     Ry = 0
     Rz = 0
     do Rx = 1,int(Rmaxrnd)
      tadx2Rx = tadxsq**Rx
      pathi = xxdir
      pathf = xxdir
      call wilsonloop(ttdir,Rx,pathi,Rx,pathf,wloop,u,myid,nn,ldiv,nms,lv, &
                      ib,lbd,iblv,MRT)
      do it = 1,max(nx,ny,nz,nt)
       wloop(it) = wloop(it)/tadx2Rx/tadtsq**it
      enddo ! it
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
            action="write",form="formatted",status="old",position="append")
        write(unit=8,fmt="(i10,a12,4i1,3i3,500es17.10)") icfgsave, &
             " wloop (1D):", xxdir, yydir, zzdir, ttdir, Rx, Ry, Rz, wloop
       close(unit=8,status="keep")
      endif
     enddo ! Rx
 
! Construct the Wilson loops that require two "spatial" directions.
     if (nwilo<0) then
      tadtsq = 1.0_KR
      tadxsq = 1.0_KR
      tadysq = 1.0_KR
     else
      tadtsq = tad(ttdir)**2
      tadxsq = tad(xxdir)**2
      tadysq = tad(yydir)**2
     endif
     Rz = 0
     do Rx = 1,int(Rmaxrnd)
      nperms = 1
      do Ry = 1,Rx
       Rsum = Rx + Ry
       nperms = nperms*Rsum/Ry
       if (nperms<=npaths.and.sqrt(real(Rx**2+Ry**2,KR))<Rmaxrnd) then
        call perm(xxdir,Rx,yydir,Ry,zzdir,Rz,nperms,permvecs)
        wloop = 0.0_KR
        do iperms = 1,nperms
         pathi(:) = permvecs(iperms,:)
         do fperms = 1,nperms
          pathf(:) = permvecs(fperms,:)
          call wilsonloop(ttdir,Rsum,pathi,Rsum,pathf,wloopbit,u,myid,nn, &
                          ldiv,nms,lv,ib,lbd,iblv,MRT)
          wloop = wloop + wloopbit
         enddo ! fperms
        enddo ! iperms
        do it = 1,max(nx,ny,nz,nt)
         wloop(it) = wloop(it)/nperms**2/tadxsq**Rx/tadysq**Ry/tadtsq**it
        enddo ! it
        if (myid==0) then
         open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
              action="write",form="formatted",status="old",position="append")
           write(unit=8,fmt="(i10,a12,4i1,3i3,500es17.10)") icfgsave, &
                " wloop (2D):", xxdir, yydir, zzdir, ttdir, Rx, Ry, Rz, wloop
         close(unit=8,status="keep")
        endif
       endif
      enddo ! Ry
     enddo ! Rx
 
! Construct the Wilson loops that require three "spatial" directions.
     if (nwilo<0) then
      tadtsq = 1.0_KR
      tadxsq = 1.0_KR
      tadysq = 1.0_KR
      tadzsq = 1.0_KR
     else
      tadtsq = tad(ttdir)**2
      tadxsq = tad(xxdir)**2
      tadysq = tad(yydir)**2
      tadzsq = tad(zzdir)**2
     endif
     do Rx = 1,int(Rmaxrnd)
      npermsbit = 1
      do Ry = 1,Rx
       npermsbit = npermsbit*(Rx+Ry)/Ry
       nperms = npermsbit
       do Rz = 1,Ry
        nperms = nperms*(Rx+Ry+Rz)/Rz
        if (nperms<=npaths.and.sqrt(real(Rx**2+Ry**2+Rz**2,KR))<Rmaxrnd) then
         Rsum = Rx + Ry + Rz
         call perm(xxdir,Rx,yydir,Ry,zzdir,Rz,nperms,permvecs)
         wloop = 0.0_KR
         do iperms = 1,nperms
          pathi(:) = permvecs(iperms,:)
          do fperms = 1,nperms
           pathf(:) = permvecs(fperms,:)
           call wilsonloop(ttdir,Rsum,pathi,Rsum,pathf,wloopbit,u,myid,nn, &
                           ldiv,nms,lv,ib,lbd,iblv,MRT)
           wloop = wloop + wloopbit
          enddo ! fperms
         enddo ! iperms
         do it = 1,max(nx,ny,nz,nt)
          wloop(it) = wloop(it)/nperms**2/tadxsq**Rx/tadysq**Ry/tadzsq**Rz &
                      /tadtsq**it
         enddo ! it
         if (myid==0) then
          open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
               action="write",form="formatted",status="old",position="append")
           write(unit=8,fmt="(i10,a12,4i1,3i3,500es17.10)") icfgsave, &
                " wloop (3D):", xxdir, yydir, zzdir, ttdir, Rx, Ry, Rz, wloop
          close(unit=8,status="keep")
         endif
        endif
       enddo ! Rz
      enddo ! Ry
     enddo ! Rx

! If the links have been fuzzed, then read in the original links.
     if (nfuzz>0) call GFread(u,trim(rwdir(myid+1))//trim(cfgfile), &
                              gaction,beta,alat,tadin,icfgsavein,myidin)

    enddo ! iwilo

 end subroutine somewloops

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine perm(mu,Nmu,nu,Nnu,la,Nla,ionumvecs,permvecs)
! Construct a set of "i" vectors, called permvecs(i,j), with j=Nmu+Nnu+Nla.
! Each vector has Nmu entries equal to mu, Nnu entries equal to nu 
! and Nla entries equal to la.  Each vector is a unique permutation
! of the entries, so there are (Nmu+Nnu+Nla)!/Nmu!/Nnu!/Nla! vectors.
! NOTE: On output, ionumvecs = (Nmu+Nnu+Nla)!/Nmu!/Nnu!/Nla!.

    integer(kind=KI), intent(in)                   :: mu, Nmu, nu, Nnu, la, Nla
    integer(kind=KI), intent(out)                  :: ionumvecs
    integer(kind=KI), intent(out),  dimension(:,:) :: permvecs

    integer(kind=KI) :: numvecs, lenvecs, ieo, jeo, ieosave, icount, &
                        i, j, k, jloop, inu, ila
    integer(kind=KI), dimension(2,npathsmax,twoRmax) :: temp

! Begin with Nmu links of type mu.
    numvecs = 1
    lenvecs = Nmu
    ieo = 1
    jeo = 2
    temp(ieo,1,:) = mu

! Insert the next link of type nu in all possible locations.
    if (Nnu>0) then
     do inu = 1,Nnu
      icount = 0
      do i = 1,numvecs
!-find the rightmost nu-type link.
       j = lenvecs
       lastnu: do
        if (j==0) exit lastnu
        if (temp(ieo,i,j)==nu) exit lastnu
        j = j - 1
       enddo lastnu
!-create new permutations.
       do jloop = j,lenvecs
        icount = icount + 1
        if (jloop>0) then
         do k = 1,jloop
          temp(jeo,icount,k) = temp(ieo,i,k)
         enddo ! k
        endif
        temp(jeo,icount,jloop+1) = nu
        if (jloop<lenvecs) then
         do k = jloop+1,lenvecs
          temp(jeo,icount,k+1) = temp(ieo,i,k)
         enddo ! k
        endif
       enddo ! jloop
!-end do loops.
      enddo ! i
      ieosave = ieo
      ieo = jeo
      jeo = ieosave
      lenvecs = lenvecs + 1
      numvecs = icount
     enddo ! inu
    endif

! Insert the next link of type la in all possible locations.
    if (Nla>0) then
     do ila = 1,Nla
      icount = 0
      do i = 1,numvecs
!-find the rightmost la-type link.
       j = lenvecs
       lastla: do
        if (j==0) exit lastla
        if (temp(ieo,i,j)==la) exit lastla
        j = j - 1
       enddo lastla
!-create new permutations.
       do jloop = j,lenvecs
        icount = icount + 1
        if (jloop>0) then
         do k = 1,jloop
          temp(jeo,icount,k) = temp(ieo,i,k)
         enddo ! k
        endif
        temp(jeo,icount,jloop+1) = la
        if (jloop<lenvecs) then
         do k = jloop+1,lenvecs
          temp(jeo,icount,k+1) = temp(ieo,i,k)
         enddo ! k
        endif
       enddo ! jloop
!-end do loops.
      enddo ! i
      ieosave = ieo
      ieo = jeo
      jeo = ieosave
      lenvecs = lenvecs + 1
      numvecs = icount
     enddo ! ila
    endif

! Write the final permutations to the output matrix.
    permvecs = 0
    do j = 1,Nmu+Nnu+Nla
     do i = 1,numvecs
      permvecs(i,j) = temp(ieo,i,j)
     enddo ! i
    enddo ! j
    ionumvecs = numvecs

 end subroutine perm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine gamma5(a,b,kappa)

 real(kind=KR), intent(inout), dimension(:,:,:,:,:)        :: a,b
 real(kind=KR), intent(in),    dimension(:)                :: kappa
 real(kind=KR),                dimension(6,ntotal,2,2,8,2) :: temp
 integer(kind=KI)                                          :: fac

           fac = 1.0_KR/(2.0_KR*kappa(1)*kappa(2))

           temp(:,:,1,:,:,1) = a(:,:,1,:,:)
           temp(:,:,2,:,:,1) = a(:,:,2,:,:)
           temp(:,:,1,:,:,2) = b(:,:,1,:,:)
           temp(:,:,2,:,:,2) = b(:,:,2,:,:)

           a(:,:,1,:,:) = -fac*a(:,:,3,:,:)
           a(:,:,2,:,:) = -fac*a(:,:,4,:,:)
           a(:,:,3,:,:) =  fac*temp(:,:,1,:,:,1)
           a(:,:,4,:,:) =  fac*temp(:,:,2,:,:,1)

           b(:,:,1,:,:) = -fac*b(:,:,3,:,:)
           b(:,:,2,:,:) = -fac*b(:,:,4,:,:)
           b(:,:,3,:,:) =  fac*temp(:,:,1,:,:,2)
           b(:,:,4,:,:) =  fac*temp(:,:,2,:,:,2)

 
 end subroutine gamma5
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine check4ones(be,bo,icsrc,idsrc)
   
    real(kind=KR), intent(in), dimension(6,ntotal,4,2,8)     :: be, bo
    integer(kind=KI), intent(in) :: icsrc, idsrc
    integer(kind=KI) :: i,j,k,ic

      if (icsrc==1) then
        ic =1
      elseif(icsrc==2) then
        ic=3
      elseif(icsrc==3)  then
        ic=5
      endif

      print *, "Checking for ones"
      do i =1,nvhalf
          do j=1,2
            do k=1,8
              if(abs(be(ic,i,idsrc,j,k)) == 1.0_KR) then
                print "(a32,5i3)", "be(icsrc,i,idsrc,j,k)=", icsrc,i,idsrc,j,k
              endif
!              if(abs(bo(ic,i,idsrc,j,k)) /= 1.0_KR) then
!                print *, "bo"
!                print *, "icsrc,idsrc,j,k=", icsrc,idsrc,j,k
!             endif
            enddo ! k
           enddo ! j
        enddo ! i
   end subroutine check4ones
 end module cfgsprops

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
