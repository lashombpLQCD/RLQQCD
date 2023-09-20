!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! inverters.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module performs the fermion matrix inversion.
! Useful references: 
!
! "Accelerating Wilson fermion matrix inversions by means of the
! stabilized biconjugate gradient algorithm"
! Frommer, Hannemann, Nockel, Lippert & Schilling, Int.J.Mod.Phys.C5,1073(1994)
!
! "Progress on lattice QCD algorithms"
! Ph. de Forcrand, Nucl.Phys.Proc.Suppl.47:228(1996)
!
! "Krylov space solvers for shifted linear systems"
! B. Jegerlehner, hep-lat/9612014.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

 module inverters

!   use MPI
    use kinds
    use latdims
    use basics
    use diracops
    use pseudolapack
    use gmresrhs
    implicit none
    private

! Use the following line if the MPI module is not available.
    include 'mpif.h'

! Define access to subroutines.
    public  :: bicgstab, m3r, cgne, qmr, qmrgam5, gmres, ppgmresproject, &
               ppgmresdr, gmresdr5EIG, gmresdrEIG, gmresdr, mmgmresdr,&
               mmgmresproject, gmresdrshift,ppmmgmresdr,ppmmgmresproject, &
               gmresproject, LUinversion, psc,landr,cg_proj,landr_proj,cg, &
               vecdot, vecdot1, minresdr5EIG, minresdr
    private :: ludcmp, lubksb, mprove, vecdotmyid, changevectorX, &
               qrfactorizationLAPACK, leastsquaresLAPACK, eigmodesLAPACK 

 contains

!!!!!!!



 
 subroutine leastsquaresLAPACK(a,m,n,b,x,myid)
 real(kind=KR2),   intent(in), dimension(:,:,:) :: a
 real(kind=KR2),   intent(in), dimension(:,:)   :: b
 real(kind=KR2),   intent(out), dimension(:,:)  :: x
 integer(kind=KI), intent(in)                   :: m,n,myid

 complex*16, allocatable, dimension(:,:) :: a_lap
 complex*16, allocatable, dimension(:) :: work
 complex*16, allocatable, dimension(:,:) :: b_lap

! complex*16, dimension(m,n) :: a_lap
! complex*16, dimension(10*m*n) :: work
! complex*16, dimension(m,1) :: b_lap

 CHARACTER*1 trans
 integer(kind=KI) info, lda, ldb, lwork, nrhs

 integer(kind=KI) i,j,k



 allocate(a_lap(m,n))
 allocate(work(10*m*n))
 allocate(b_lap(m,1))


 trans = 'N'
 nrhs = 1
 lda = m
 ldb = m
 lwork = 10*m*n

 do i=1,m
   do j=1,n
     a_lap(i,j) = cmplx(a(1,i,j),a(2,i,j),KR2)
   enddo
 enddo

 do i=1,m
   b_lap(i,1) = cmplx(b(1,i),b(2,i),KR2)
 enddo 

 call zgels(trans,m,n,nrhs,a_lap,lda,b_lap,ldb,work,lwork,info)

 do i=1,n
   x(1,i) = real(b_lap(i,1),KR2)
   x(2,i) = aimag(b_lap(i,1))
 enddo

 deallocate(a_lap)
 deallocate(work)
 deallocate(b_lap)
 end subroutine leastsquaresLAPACK


! mode = 1 calculate rite evectors, otherwise not
subroutine eigmodesLAPACK(mat,m,mode,eval,evec)
 real(kind=KR2),   intent(in), dimension(:,:,:)  :: mat
 integer(kind=KI), intent(in)                    :: m
 integer(kind=KI), intent(in)                    :: mode
 real(kind=KR2),   intent(out), dimension(:,:)   :: eval
 real(kind=KR2),   intent(out), dimension(:,:,:) :: evec

 ! local LAPACK variables
 complex*16, allocatable, dimension(:,:) :: a
 complex*16, allocatable, dimension(:,:) :: vl
 complex*16, allocatable, dimension(:,:) :: vr
 complex*16, allocatable, dimension(:)   :: w
 complex*16, allocatable, dimension(:)   :: work
!complex(kind=KR2), allocatable, dimension(:,:) :: a!BS changed comples to KR2
!complex(kind=KR2), allocatable, dimension(:,:) :: vl
!complex(kind=KR2), allocatable, dimension(:,:) :: vr
!complex(kind=KR2), allocatable, dimension(:)   :: w
! complex(kind=KR2), allocatable, dimension(:)   :: work
 ! real*16,    allocatable, dimension(:)   :: rwork
 real(kind=KR2),    allocatable, dimension(:)   :: rwork

 CHARACTER*1 jobvl
 CHARACTER*1 jobvr

 integer(kind=KI) info, lda, ldvl, ldvr, lwork, n

 real(kind=KR2), parameter :: DOUBLEPRECISIONLIMIT = 0 !1e-14_KR2

 ! loop variables
 integer(kind=KI) i,j,k

 ! set LAPACK integer values
 n     = m
 lda   = m
 ldvl  = m
 ldvr  = m
 lwork = 10*n ! notes give as >= max(1,2*n)
 jobvl = 'N'
 if (mode == 1) then
   jobvr = 'V'! else
   jobvr = 'N'
 endif

 ! allocate appropriate space
 allocate(a(lda,n))
 allocate(vl(ldvl,n))
 allocate(vr(ldvr,n))
 allocate(w(n))
 allocate(work(lwork))
 allocate(rwork(2*n))

 ! copy into LAPACK type arrays
 do i=1,lda
   do j=1,n
       a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
          if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(0.0,mat(2,i,j),KR2)
     else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(mat(1,i,j),0.0,KR2)
     else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
     else if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(0.0,0.0,KR2)
     else
       write(*,*) 'AHHHHHHHHHHHHH!!!!!!!!!!!',DOUBLEPRECISIONLIMIT 
     endif
   enddo
 enddo



 ! allocate appropriate space
 !allocate(a(lda,n))
 !allocate(vl(ldvl,n))
 !allocate(vr(ldvr,n))
 !allocate(w(n))
 !allocate(work(lwork))
 !allocate(rwork(2*n))

 ! copy into LAPACK type arrays
 !do i=1,lda
 !  do j=1,n
      ! a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
  !        if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >=
  !        DOUBLEPRECISIONLIMIT) then
  !     a(i,j) = cmplx(0.0,mat(2,i,j),KR2)
  !   else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <=
  !   DOUBLEPRECISIONLIMIT) then
  !     a(i,j) = cmplx(mat(1,i,j),0.0,KR2)
  !   else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >=
  !   DOUBLEPRECISIONLIMIT) then
  !     a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
  !   else if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <=
  !   DOUBLEPRECISIONLIMIT) then
  !     a(i,j) = cmplx(0.0,0.0,KR2)
  !   else
  !     write(*,*) 'AHHHHHHHHHHHHH!!!!!!!!!!!',DOUBLEPRECISIONLIMIT
  !   endif
  ! enddo
! enddo


    print * ,'about to call zgeev','  mode=',mode
 ! call LAPACK routine
 call zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
    
    print * ,'just  called zgeev','mode=',mode
 !   print * , 'w=',w,'vr =',vr,'info',info
 ! error checking
 if (info > 0) then
   write(*,*) 'ERROR: QR algorithm failed - zgeev'
   stop
 endif

 ! move evalues into final storage position
 do i=1,n
   eval(1,i) = real(w(i),KR2)
   eval(2,i) = aimag(w(i))
 enddo

 ! move evectors into final storage position
 do i=1,ldvr
   do j=1,n
     evec(1,i,j) = real(vr(i,j),KR2)
     evec(2,i,j) = aimag(vr(i,j))
   enddo
 enddo

 deallocate(a)
 deallocate(vl)
 deallocate(vr)
 deallocate(w)
 deallocate(work)
 deallocate(rwork)
 end subroutine eigmodesLAPACK



 subroutine qrfactorizationLAPACK(a,m,n,qfact,rfact,myid)
 integer(kind=KI), intent(in)                    :: m, n,myid 
 real(kind=KR2),   intent(in),  dimension(:,:,:) :: a
 real(kind=KR2),   intent(out), dimension(:,:,:) :: qfact
 real(kind=KR2),   intent(out), dimension(:,:,:) :: rfact

 complex*16, allocatable, dimension(:,:) :: a_lap
 complex*16, allocatable, dimension(:)   :: tau
 complex*16, allocatable, dimension(:)   :: work
 integer(kind=KI) :: lda, info, lwork

 integer(kind=KI) :: i,j,k
 
 allocate(a_lap(m,n))
 allocate(tau(min(m,n)))
 allocate(work(10*n))

 lda = m
 lwork = 10*n
 k = min(m,n)

 do i=1,m
   do j=1,n
     a_lap(i,j) = cmplx(a(1,i,j),a(2,i,j),KR2)
   enddo
 enddo

 call zgeqrf(m,n,a_lap,lda,tau,work,lwork,info)

 rfact = 0.0_KR2
 do i=1,m
   do j=i,n
     rfact(1,i,j) = real(a_lap(i,j),KR2)
     rfact(2,i,j) = aimag(a_lap(i,j))
   enddo
 enddo

 call zungqr(m,n,k,a_lap,lda,tau,work,lwork,info)

 qfact = 0.0_KR2
 do i=1,m
   do j=1,n
     qfact(1,i,j) = real(a_lap(i,j),KR2)
     qfact(2,i,j) = aimag(a_lap(i,j))
   enddo
 enddo

 deallocate(a_lap)
 deallocate(tau)
 deallocate(work)
 end subroutine qrfactorizationLAPACK
!!!!!!!

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

 subroutine bicgstab(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv,iflag, &
                     kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                     lbd,iblv,MRT,MRT2)
! Stabilized biconjugate gradient matrix inverter.
! Solves M*x=phi for the vector x.
! The notation of Frommer et al, Int.J.Mod.Phys.C5,1073(1994) is followed.
! (except that their vector "s" is stored in the vector "r" to conserve memory)
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by bicgstab.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(in)             :: itermin, iflag, myid, MRT, MRT2
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR),  dimension(6,ntotal,4,2,8) :: r, p
    real(kind=KR),  dimension(6,nvhalf,4,2,8) :: rhat, v, t
    real(kind=KR2), dimension(2)              :: rho, oldrho, alpha, beta, &
                                                 omega, betaomega, &
                                                 alphatmp, top, bottom, rsq,&
                                                 r0srt,exitnorm
    real(kind=KR2)                            :: betaomegabit, betabit, &
                                                 alphabit, omegabit
    integer(kind=KI)                          :: i, icri, idag

! Initializations.
    itercount = 0
    idag = 0
    call Hdbletm(t,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - t(:,i,:,:,:)
    enddo ! i
    do i = 1,nvhalf
     rhat(:,i,:,:,:) = r(:,i,:,:,:)
    enddo ! i
    call vecdot(r,r,r0srt,MRT2)
    r0srt(1) =sqrt(r0srt(1))


    rho(1) = 1.0_KR2
    rho(2) = 0.0_KR2
    alpha = rho
    omega = rho
    v = 0.0_KR
    p = 0.0_KR
    
! Main loop.
    maindo: do 
     oldrho = rho
     call vecdot(rhat,r,rho,MRT2)
     betaomegabit = 1.0_KR2/(oldrho(1)**2+oldrho(2)**2)
     betaomega(1) = betaomegabit*(rho(1)*oldrho(1)*alpha(1) &
                                 +rho(1)*oldrho(2)*alpha(2) &
                                 -rho(2)*oldrho(1)*alpha(2) &
                                 +rho(2)*oldrho(2)*alpha(1))
     betaomega(2) = betaomegabit*(rho(2)*oldrho(1)*alpha(1) &
                                 -rho(1)*oldrho(2)*alpha(1) &
                                 +rho(1)*oldrho(1)*alpha(2) &
                                 +rho(2)*oldrho(2)*alpha(2))
     betabit = 1.0_KR2/(omega(1)**2+omega(2)**2)
     beta(1) = betabit*(betaomega(1)*omega(1)+betaomega(2)*omega(2))
     beta(2) = betabit*(betaomega(2)*omega(1)-betaomega(1)*omega(2))
     do i = 1,nvhalf
      t(:,i,:,:,:) = p(:,i,:,:,:)
     enddo ! i
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       p(icri  ,i,:,:,:) = r(icri  ,i,:,:,:)                                  &
                 + beta(1)*t(icri  ,i,:,:,:) - betaomega(1)*v(icri  ,i,:,:,:) &
                 - beta(2)*t(icri+1,i,:,:,:) + betaomega(2)*v(icri+1,i,:,:,:)
       p(icri+1,i,:,:,:) = r(icri+1,i,:,:,:)                                  &
                 + beta(1)*t(icri+1,i,:,:,:) - betaomega(1)*v(icri+1,i,:,:,:) &
                 + beta(2)*t(icri  ,i,:,:,:) - betaomega(2)*v(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     call Hdbletm(v,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                  myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(rhat,v,alphatmp,MRT2)
     alphabit = 1.0_KR2/(alphatmp(1)**2+alphatmp(2)**2)
     alpha(1) = alphabit*(rho(1)*alphatmp(1)+rho(2)*alphatmp(2))
     alpha(2) = alphabit*(rho(2)*alphatmp(1)-rho(1)*alphatmp(2))
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       r(icri  ,i,:,:,:) = r(icri  ,i,:,:,:)                              &
                - alpha(1)*v(icri  ,i,:,:,:) + alpha(2)*v(icri+1,i,:,:,:)
       r(icri+1,i,:,:,:) = r(icri+1,i,:,:,:)                              &
                - alpha(1)*v(icri+1,i,:,:,:) - alpha(2)*v(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     call Hdbletm(t,u,GeeGooinv,r,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                  myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(t,r,top,MRT2)
     call vecdot(t,t,bottom,MRT2)
     omegabit = 1.0_KR2/(bottom(1)**2+bottom(2)**2)
     omega(1) = omegabit*(top(1)*bottom(1)+top(2)*bottom(2))
     omega(2) = omegabit*(top(2)*bottom(1)-top(1)*bottom(2))
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       x(icri  ,i,:,:,:) = x(icri  ,i,:,:,:)                              &
                + omega(1)*r(icri  ,i,:,:,:) + alpha(1)*p(icri  ,i,:,:,:) &
                - omega(2)*r(icri+1,i,:,:,:) - alpha(2)*p(icri+1,i,:,:,:)
       x(icri+1,i,:,:,:) = x(icri+1,i,:,:,:)                              &
                + omega(1)*r(icri+1,i,:,:,:) + alpha(1)*p(icri+1,i,:,:,:) &
                + omega(2)*r(icri  ,i,:,:,:) + alpha(2)*p(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       r(icri  ,i,:,:,:) = r(icri  ,i,:,:,:) &
                - omega(1)*t(icri  ,i,:,:,:) + omega(2)*t(icri+1,i,:,:,:)
       r(icri+1,i,:,:,:) = r(icri+1,i,:,:,:) &
                - omega(1)*t(icri+1,i,:,:,:) - omega(2)*t(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     itercount = itercount + 1
     call vecdot(r,r,rsq,MRT2)
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       write(unit=8,fmt="(a12,i9,es17.10)") "bicgstab",itercount,rsq(1)/r0srt(1)
      close(unit=8,status="keep")
     endif
     exitnorm(1) = sqrt(rsq(1))
     exitnorm(1) = exitnorm(1)/r0srt(1)
     if (exitnorm(1)<resmax.and.itercount>=itermin) exit maindo
    enddo maindo

 end subroutine bicgstab

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

 subroutine m3r(rwdir,phi,x,resmax,itermin,itercount,omega,u,nkappa,kappa, &
                coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv, &
                MRT,MRT2)
! Multi-mass minimal residual matrix inverter.
! Solves M_i*x=phi for the vector x, where M_i is the Dirac matrix with kappa_i.
! Ref's: Glassner,Gusken,Lippert,Ritzenhofer,Schilling,Frommer,hep-lat/9605008.
!        Ying,Dong,Liu,Nucl.Phys.B(Proc.Suppl.)53(1997)993.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   omega is an overrelaxation parameter.  Ying, Dong and Liu use omega=1.1,
!         but I have found this not to converge for some configurations when
!         the quark mass is light.  omega=1.0 means no overrelaxation.

!   u() contains the gauge fields for this sublattice.
!   nkappa is the number of hopping parameters.
!   kappa() contains the hopping parameters.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by this subroutine.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:)   :: phi
    real(kind=KR),    intent(out), dimension(:,:,:,:,:,:) :: x
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(in)           :: itermin, nkappa, myid, MRT, MRT2
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR2),   intent(in)                          :: omega
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:)   :: u
    real(kind=KR),    intent(in),  dimension(:)           :: kappa
    real(kind=KR),    intent(in),  dimension(:,:,:)       :: coact
    integer(kind=KI), intent(in),  dimension(:)           :: bc, nms
    integer(kind=KI), intent(in),  dimension(:,:)         :: vecbl, vecblinv
    integer(kind=KI), intent(in),  dimension(:,:)         :: nn, iblv
    logical,          intent(in),  dimension(:)           :: ldiv
    integer(kind=KI), intent(in),  dimension(:,:,:)       :: lvbc
    integer(kind=KI), intent(in),  dimension(:,:,:,:)     :: ib
    logical,          intent(in),  dimension(:,:)         :: lbd

    real(kind=KR),  dimension(6,ntotal,4,2,8) :: r
    real(kind=KR),  dimension(6,nvhalf,4,2,8) :: p
    real(kind=KR2), dimension(2)              :: alphatop, alphabot, alpha, rsq,tv
!   integer(kind=KI), parameter               :: kmax=9
    real(kind=KR2), dimension(2,kmax)         :: ffac
    real(kind=KR2), dimension(kmax)           :: sigma
    real(kind=KR2)                            :: denRe, denIm, oldfRe, oldfIm,resnum,normnum
    integer(kind=KI)                          :: i, icri, ikappa, idag

! Initializations.
    itercount = 0
    idag = 0
    x = 0.0_KR
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:)
    enddo ! i

    call vecdot(r,r,tv,MRT)
    normnum=sqrt(tv(1))

    ffac(1,:) = 1.0_KR2
    ffac(2,:) = 0.0_KR2
    do ikappa = 2,nkappa
     sigma(ikappa) = 1.0_KR2/kappa(ikappa)**2 - 1.0_KR2/kappa(1)**2
    enddo ! ikappa

! Main loop.
    maindo: do 
     call Hdouble(p,u,r,idag,coact,kappa(1),bc,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(p,r,alphatop,MRT2)
     call vecdot(p,p,alphabot,MRT2)
     alpha = omega/alphabot(1)*alphatop
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       x(icri  ,i,:,:,:,1) = x(icri  ,i,:,:,:,1)        &
                           + alpha(1)*r(icri  ,i,:,:,:) &
                           - alpha(2)*r(icri+1,i,:,:,:)
       x(icri+1,i,:,:,:,1) = x(icri+1,i,:,:,:,1)        &
                           + alpha(1)*r(icri+1,i,:,:,:) &
                           + alpha(2)*r(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     do ikappa = 2,nkappa
      denRe = 1.0_KR2 + sigma(ikappa)*alpha(1)
      denIm = sigma(ikappa)*alpha(2)
      oldfRe = ffac(1,ikappa)
      oldfIm = ffac(2,ikappa)
      ffac(1,ikappa) = (oldfRe*denRe+oldfIm*denIm)/(denRe**2+denIm**2)
      ffac(2,ikappa) = (oldfIm*denRe-oldfRe*denIm)/(denRe**2+denIm**2)
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        x(icri  ,i,:,:,:,ikappa) = x(icri  ,i,:,:,:,ikappa)                  &
                                 + ffac(1,ikappa)*alpha(1)*r(icri  ,i,:,:,:) &
                                 - ffac(1,ikappa)*alpha(2)*r(icri+1,i,:,:,:) &
                                 - ffac(2,ikappa)*alpha(1)*r(icri+1,i,:,:,:) &
                                 - ffac(2,ikappa)*alpha(2)*r(icri  ,i,:,:,:)
        x(icri+1,i,:,:,:,ikappa) = x(icri+1,i,:,:,:,ikappa)                  &
                                 + ffac(1,ikappa)*alpha(1)*r(icri+1,i,:,:,:) &
                                 + ffac(1,ikappa)*alpha(2)*r(icri  ,i,:,:,:) &
                                 + ffac(2,ikappa)*alpha(1)*r(icri  ,i,:,:,:) &
                                 - ffac(2,ikappa)*alpha(2)*r(icri+1,i,:,:,:)
       enddo ! icri
      enddo ! i
     enddo ! ikappa
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       r(icri  ,i,:,:,:) = r(icri  ,i,:,:,:) - alpha(1)*p(icri  ,i,:,:,:) &
                                             + alpha(2)*p(icri+1,i,:,:,:)
       r(icri+1,i,:,:,:) = r(icri+1,i,:,:,:) - alpha(1)*p(icri+1,i,:,:,:) &
                                             - alpha(2)*p(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     itercount = itercount + 1
     call vecdot(r,r,rsq,MRT2)
     resnum=sqrt(rsq(1))/normnum
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       !write(unit=8,fmt="(a12,i9,es17.10)") "m3r",itercount,rsq(1) !Randy
       write(unit=8,fmt="(a12,i9,es17.10)") "m3r",itercount,resnum  !Abdou
      close(unit=8,status="keep")
     endif
     !if (real(rsq(1),KR)<resmax.and.itercount>=itermin) exit maindo (Randy)
     if (resnum<resmax.and.itercount>=itermin) exit maindo !Abdou
    enddo maindo

 end subroutine m3r

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

 subroutine cgne(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv,iflag,kappa, &
                 coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv, &
                 MRT,MRT2)
! Conjugate gradient matrix inverter acting on the "normal equations".
! Solves M*x=phi for the vector x.  ["Normal equations" means the code
! actually works with M^dagger*M*x=M^dagger*phi].
! References: Numerical Recipes, pages 77 and 78;
!             de Forcrand, Nucl.Phys.B(Proc.Suppl.)47,228(1996).
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by cgne.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(in)             :: itermin, iflag, myid, MRT, MRT2
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR),  dimension(6,ntotal,4,2,8) :: r, p
    real(kind=KR),  dimension(6,nvhalf,4,2,8) :: Mvec
    real(kind=KR2), dimension(2)              :: Mdagrsq, Mpsq, rsq, r0srt,exitnorm
    real(kind=KR2)                            :: Mdagrsqold, alpha, beta
    integer(kind=KI)                          :: idag, i,iterprint

! Initialize matrices that will only get defined via a subroutine call.
    p = 0.0_KR
    Mvec = 0.0_KR

! Initialization steps for the cgne.
    itercount = 0
    idag = 0
    call Hdbletm(p,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - p(:,i,:,:,:)
    enddo ! i

    call vecdot(r,r,r0srt,MRT2)
    r0srt(1) =sqrt(r0srt(1)) 

    idag = 1
    call Hdbletm(Mvec,u,GeeGooinv,r,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     p(:,i,:,:,:) = Mvec(:,i,:,:,:)
    enddo ! i
    call vecdot(Mvec,Mvec,Mdagrsq,MRT2)

! Main loop.
    maindo: do 
     idag = 0
     call Hdbletm(Mvec,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl, &
                  vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(Mvec,Mvec,Mpsq,MRT2)
     alpha = Mdagrsq(1)/Mpsq(1)
     do i = 1,nvhalf
      r(:,i,:,:,:) = r(:,i,:,:,:) - alpha*Mvec(:,i,:,:,:)
      x(:,i,:,:,:) = x(:,i,:,:,:) + alpha*p(:,i,:,:,:)
     enddo ! i
     itercount = itercount + 1
     call vecdot(r,r,rsq,MRT2)
     iterprint = 5
     if (myid==0) then
      if(mod(itercount,iterprint) ==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       write(unit=8,fmt="(a12,i9,es17.10)") "cgne",itercount,sqrt(rsq(1))/r0srt(1)
      close(unit=8,status="keep")
      endif ! mod
     endif
     exitnorm(1) = rsq(1)
     exitnorm(1) = sqrt(exitnorm(1))/r0srt(1)
     if (exitnorm(1) <resmax.and.itercount>=itermin) exit maindo
     Mdagrsqold = Mdagrsq(1)
     idag = 1
     call Hdbletm(Mvec,u,GeeGooinv,r,idag,coact,kappa,iflag,bc,vecbl, &
                  vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(Mvec,Mvec,Mdagrsq,MRT2)
     beta = Mdagrsq(1)/Mdagrsqold
     do i = 1,nvhalf
      p(:,i,:,:,:) = Mvec(:,i,:,:,:) + beta*p(:,i,:,:,:)
     enddo ! i
    enddo maindo

 end subroutine cgne

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine qmr(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv,iflag,kappa, &
                coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv, &
                MRT,MRT2)
! Quasi minimal residual matrix inverter (no look-ahead).
! Solves M*x=phi for the vector x.
! References: www.netlib.org/templates (I follow this netlib reference)
!             Frommer et al, hep-lat/9504020.
!             Ying, Dong and Yiu, Nucl. Phys. B (Proc.Suppl.) 53, 993 (1997).
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by qmr.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(in)             :: itermin, iflag, myid, MRT, MRT2
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR),  dimension(6,ntotal,4,2,8) :: p, q, pold, qold
    real(kind=KR),  dimension(6,nvhalf,4,2,8) :: vtilde, ptilde, v, r, &
                                                 d, s, w, wtilde, qtilde
    real(kind=KR2), dimension(2)              :: vecsq, delta, beta, eta, &
                                                 etaold, eps
    real(kind=KR2)                            :: theta, gam, gammaold, &
                                                 thetaold, deltabit
    real(kind=KR)                             :: rho, rhoold, xi
    integer(kind=KI)                          :: i, icri, idag

! Initializations.
    itercount = 0
    idag = 0
    call Hdbletm(v,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - v(:,i,:,:,:)
    enddo ! i
    vtilde = r
    call vecdot(vtilde,vtilde,vecsq,MRT2)
    rho = sqrt(vecsq(1))
    wtilde = r
    call vecdot(wtilde,wtilde,vecsq,MRT2)
    xi = sqrt(vecsq(1))
    gam = 1.0_KR2
    eta(1) = -1.0_KR2
    eta(2) = 0.0_KR2

! Main loop.
    maindo: do 
     do i = 1,nvhalf
      v(:,i,:,:,:) = vtilde(:,i,:,:,:)/rho
      w(:,i,:,:,:) = wtilde(:,i,:,:,:)/xi
     enddo ! i
     call vecdot(w,v,delta,MRT2)
     if (itercount==0) then
      do i = 1,nvhalf
       p(:,i,:,:,:) = v(:,i,:,:,:)
       q(:,i,:,:,:) = w(:,i,:,:,:)
      enddo ! i
     else
      vecsq(1) = (delta(1)*eps(1)+delta(2)*eps(2))/(eps(1)**2+eps(2)**2)
      vecsq(2) = (delta(2)*eps(1)-delta(1)*eps(2))/(eps(1)**2+eps(2)**2)
      pold = p
      qold = q
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        p(icri  ,i,:,:,:) = v(icri  ,i,:,:,:)-vecsq(1)*xi*pold(icri  ,i,:,:,:) &
                                             +vecsq(2)*xi*pold(icri+1,i,:,:,:)
        p(icri+1,i,:,:,:) = v(icri+1,i,:,:,:)-vecsq(2)*xi*pold(icri  ,i,:,:,:) &
                                             -vecsq(1)*xi*pold(icri+1,i,:,:,:)
        q(icri  ,i,:,:,:) = w(icri  ,i,:,:,:)-vecsq(1)*rho*qold(icri  ,i,:,:,:)&
                                             -vecsq(2)*rho*qold(icri+1,i,:,:,:)
        q(icri+1,i,:,:,:) = w(icri+1,i,:,:,:)+vecsq(2)*rho*qold(icri  ,i,:,:,:)&
                                             -vecsq(1)*rho*qold(icri+1,i,:,:,:)
       enddo ! icri
      enddo ! i
     endif
     idag = 0
     call Hdbletm(ptilde,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl, &
                  vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vecdot(q,ptilde,eps,MRT2)
     deltabit = 1.0_KR2/(delta(1)**2+delta(2)**2)
     beta(1) = deltabit*(eps(1)*delta(1)+eps(2)*delta(2))
     beta(2) = deltabit*(eps(2)*delta(1)-eps(1)*delta(2))
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       vtilde(icri  ,i,:,:,:) = ptilde(icri  ,i,:,:,:) &
                  - beta(1)*v(icri  ,i,:,:,:) + beta(2)*v(icri+1,i,:,:,:)
       vtilde(icri+1,i,:,:,:) = ptilde(icri+1,i,:,:,:) &
                  - beta(2)*v(icri  ,i,:,:,:) - beta(1)*v(icri+1,i,:,:,:)
      enddo ! icri
     enddo ! i
     call vecdot(vtilde,vtilde,vecsq,MRT2)
     rhoold = rho
     rho = sqrt(vecsq(1))
     idag = 1
     call Hdbletm(qtilde,u,GeeGooinv,q,idag,coact,kappa,iflag,bc,vecbl, &
                  vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       wtilde(icri  ,i,:,:,:) = qtilde(icri  ,i,:,:,:) &
                  - beta(1)*w(icri  ,i,:,:,:) - beta(2)*w(icri+1,i,:,:,:)
       wtilde(icri+1,i,:,:,:) = qtilde(icri+1,i,:,:,:) &
                  - beta(1)*w(icri+1,i,:,:,:) + beta(2)*w(icri  ,i,:,:,:)
      enddo ! icri
     enddo ! i
     call vecdot(wtilde,wtilde,vecsq,MRT2)
     xi = sqrt(vecsq(1))
     thetaold = theta
     theta = rho/(gam*sqrt(beta(1)**2+beta(2)**2))
     gammaold = gam
     gam = 1.0_KR2/sqrt(1.0_KR2+theta**2)
     etaold = eta
     eta(1) = -(gam/gammaold)**2*rhoold/(beta(1)**2+beta(2)**2) &
              *(etaold(1)*beta(1)+etaold(2)*beta(2))
     eta(2) = -(gam/gammaold)**2*rhoold/(beta(1)**2+beta(2)**2) &
              *(etaold(2)*beta(1)-etaold(1)*beta(2))
     if (itercount==0) then
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        d(icri  ,i,:,:,:) = eta(1)*p(icri  ,i,:,:,:) - eta(2)*p(icri+1,i,:,:,:)
        d(icri+1,i,:,:,:) = eta(1)*p(icri+1,i,:,:,:) + eta(2)*p(icri  ,i,:,:,:)
        s(icri  ,i,:,:,:) = eta(1)*ptilde(icri  ,i,:,:,:) &
                          - eta(2)*ptilde(icri+1,i,:,:,:)
        s(icri+1,i,:,:,:) = eta(1)*ptilde(icri+1,i,:,:,:) &
                          + eta(2)*ptilde(icri  ,i,:,:,:)
       enddo ! icri
      enddo ! i
     else
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        d(icri  ,i,:,:,:) = eta(1)*p(icri  ,i,:,:,:) - eta(2)*p(icri+1,i,:,:,:)&
                          + (thetaold*gam)**2*d(icri  ,i,:,:,:)
        d(icri+1,i,:,:,:) = eta(1)*p(icri+1,i,:,:,:) + eta(2)*p(icri  ,i,:,:,:)&
                          + (thetaold*gam)**2*d(icri+1,i,:,:,:)
        s(icri  ,i,:,:,:) = eta(1)*ptilde(icri  ,i,:,:,:) &
                          - eta(2)*ptilde(icri+1,i,:,:,:) &
                          + (thetaold*gam)**2*s(icri  ,i,:,:,:)
        s(icri+1,i,:,:,:) = eta(1)*ptilde(icri+1,i,:,:,:) &
                          + eta(2)*ptilde(icri  ,i,:,:,:) &
                          + (thetaold*gam)**2*s(icri+1,i,:,:,:)
       enddo ! icri
      enddo ! i
     endif
     do i = 1,nvhalf
      x(:,i,:,:,:) = x(:,i,:,:,:) + d(:,i,:,:,:)
      r(:,i,:,:,:) = r(:,i,:,:,:) - s(:,i,:,:,:)
     enddo ! i
     itercount = itercount + 1
     call vecdot(r,r,vecsq,MRT2)
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       write(unit=8,fmt="(a12,i9,es17.10)") "qmr",itercount,vecsq(1)
      close(unit=8,status="keep")
     endif
     if (real(vecsq(1),KR)<resmax.and.itercount>=itermin) exit maindo
    enddo maindo

 end subroutine qmr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine qmrgam5(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv,iflag, &
                    kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                    lbd,iblv,MRT,MRT2)
! Quasi minimal residual matrix inverter with gamma5 symmetry (no look-ahead).
! Solves M*x=phi for the vector x, BUT ONLY FOR AN MATRIX THAT SATISFIES
! M = gam5*M^dagger*gam5.
! In particular, THIS SUBROUTINE IS NOT VALID FOR TWISTED MASS QCD.
! References: www.netlib.org/templates (I follow this netlib reference)
!             Frommer et al, hep-lat/9504020.
!             Ying, Dong and Yiu, Nucl. Phys. B (Proc.Suppl.) 53, 993 (1997).
!             I follow the netlib reference, but with the gamma5 symmetry
!             employed as in Frommer et al.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by qmr.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(in)             :: itermin, iflag, myid, MRT, MRT2
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR),  dimension(6,ntotal,4,2,8) :: p
    real(kind=KR),  dimension(6,nvhalf,4,2,8) :: vtilde, ptilde, v, r, &
                                                 d, s, extra
    real(kind=KR2), dimension(2)              :: vecsq, delta, beta, eta, &
                                                 etaold, eps
    real(kind=KR2)                            :: theta, gam, gammaold, &
                                                 thetaold, deltabit
    real(kind=KR)                             :: rho, rhoold
    integer(kind=KI)                          :: i, icri, idag

! Initializations.
    itercount = 0
    idag = 0
    call Hdbletm(v,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - v(:,i,:,:,:)
    enddo ! i
    vtilde = r
    call vecdot(vtilde,vtilde,vecsq,MRT2)
    rho = sqrt(vecsq(1))
    gam = 1.0_KR2
    eta(1) = -1.0_KR2
    eta(2) = 0.0_KR2

! Main loop.
    maindo: do 
     do i = 1,nvhalf
      v(:,i,:,:,:) = vtilde(:,i,:,:,:)/rho
     enddo ! i
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       extra(icri  ,i,1,:,:) =  v(icri+1,i,3,:,:)
       extra(icri+1,i,1,:,:) = -v(icri  ,i,3,:,:)
       extra(icri  ,i,2,:,:) =  v(icri+1,i,4,:,:)
       extra(icri+1,i,2,:,:) = -v(icri  ,i,4,:,:)
       extra(icri  ,i,3,:,:) = -v(icri+1,i,1,:,:)
       extra(icri+1,i,3,:,:) =  v(icri  ,i,1,:,:)
       extra(icri  ,i,4,:,:) = -v(icri+1,i,2,:,:)
       extra(icri+1,i,4,:,:) =  v(icri  ,i,2,:,:)
      enddo ! icri
     enddo ! i
     call vecdot(extra,v,delta,MRT2)
     if (itercount==0) then
      do i = 1,nvhalf
       p(:,i,:,:,:) = v(:,i,:,:,:)
      enddo ! i
     else
      vecsq(1) = rho*(delta(1)*eps(1)+delta(2)*eps(2))/(eps(1)**2+eps(2)**2)
      vecsq(2) = rho*(delta(2)*eps(1)-delta(1)*eps(2))/(eps(1)**2+eps(2)**2)
      do i = 1,nvhalf
       extra(:,i,:,:,:) = p(:,i,:,:,:)
      enddo ! i
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        p(icri  ,i,:,:,:) = v(icri  ,i,:,:,:)-vecsq(1)*extra(icri  ,i,:,:,:) &
                                             +vecsq(2)*extra(icri+1,i,:,:,:)
        p(icri+1,i,:,:,:) = v(icri+1,i,:,:,:)-vecsq(2)*extra(icri  ,i,:,:,:) &
                                             -vecsq(1)*extra(icri+1,i,:,:,:)
       enddo ! icri
      enddo ! i
     endif
     call Hdbletm(ptilde,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl, &
                  vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       extra(icri  ,i,1,:,:) =  p(icri+1,i,3,:,:)
       extra(icri+1,i,1,:,:) = -p(icri  ,i,3,:,:)
       extra(icri  ,i,2,:,:) =  p(icri+1,i,4,:,:)
       extra(icri+1,i,2,:,:) = -p(icri  ,i,4,:,:)
       extra(icri  ,i,3,:,:) = -p(icri+1,i,1,:,:)
       extra(icri+1,i,3,:,:) =  p(icri  ,i,1,:,:)
       extra(icri  ,i,4,:,:) = -p(icri+1,i,2,:,:)
       extra(icri+1,i,4,:,:) =  p(icri  ,i,2,:,:)
      enddo ! icri
     enddo ! i
     call vecdot(extra,ptilde,eps,MRT2)
     deltabit = 1.0_KR2/(delta(1)**2+delta(2)**2)
     beta(1) = deltabit*(eps(1)*delta(1)+eps(2)*delta(2))
     beta(2) = deltabit*(eps(2)*delta(1)-eps(1)*delta(2))
     do i = 1,nvhalf
      do icri = 1,5,2 ! 6=nri*nc
       vtilde(icri  ,i,:,:,:) = ptilde(icri  ,i,:,:,:) &
                  - beta(1)*v(icri  ,i,:,:,:) + beta(2)*v(icri+1,i,:,:,:)
       vtilde(icri+1,i,:,:,:) = ptilde(icri+1,i,:,:,:) &
                  - beta(2)*v(icri  ,i,:,:,:) - beta(1)*v(icri+1,i,:,:,:)
      enddo ! icri
     enddo ! i
     call vecdot(vtilde,vtilde,vecsq,MRT2)
     rhoold = rho
     rho = sqrt(vecsq(1))
     thetaold = theta
     theta = rho/(gam*sqrt(beta(1)**2+beta(2)**2))
     gammaold = gam
     gam = 1.0_KR2/sqrt(1.0_KR2+theta**2)
     etaold = eta
     eta(1) = -(gam/gammaold)**2*rhoold/(beta(1)**2+beta(2)**2) &
              *(etaold(1)*beta(1)+etaold(2)*beta(2))
     eta(2) = -(gam/gammaold)**2*rhoold/(beta(1)**2+beta(2)**2) &
              *(etaold(2)*beta(1)-etaold(1)*beta(2))
     if (itercount==0) then
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        d(icri  ,i,:,:,:) = eta(1)*p(icri  ,i,:,:,:) - eta(2)*p(icri+1,i,:,:,:)
        d(icri+1,i,:,:,:) = eta(1)*p(icri+1,i,:,:,:) + eta(2)*p(icri  ,i,:,:,:)
        s(icri  ,i,:,:,:) = eta(1)*ptilde(icri  ,i,:,:,:) &
                          - eta(2)*ptilde(icri+1,i,:,:,:)
        s(icri+1,i,:,:,:) = eta(1)*ptilde(icri+1,i,:,:,:) &
                          + eta(2)*ptilde(icri  ,i,:,:,:)
       enddo ! icri
      enddo ! i
     else
      do i = 1,nvhalf
       do icri = 1,5,2 ! 6=nri*nc
        d(icri  ,i,:,:,:) = eta(1)*p(icri  ,i,:,:,:) - eta(2)*p(icri+1,i,:,:,:)&
                          + (thetaold*gam)**2*d(icri  ,i,:,:,:)
        d(icri+1,i,:,:,:) = eta(1)*p(icri+1,i,:,:,:) + eta(2)*p(icri  ,i,:,:,:)&
                          + (thetaold*gam)**2*d(icri+1,i,:,:,:)
        s(icri  ,i,:,:,:) = eta(1)*ptilde(icri  ,i,:,:,:) &
                          - eta(2)*ptilde(icri+1,i,:,:,:) &
                          + (thetaold*gam)**2*s(icri  ,i,:,:,:)
        s(icri+1,i,:,:,:) = eta(1)*ptilde(icri+1,i,:,:,:) &
                          + eta(2)*ptilde(icri  ,i,:,:,:) &
                          + (thetaold*gam)**2*s(icri+1,i,:,:,:)
       enddo ! icri
      enddo ! i
     endif
     do i = 1,nvhalf
      x(:,i,:,:,:) = x(:,i,:,:,:) + d(:,i,:,:,:)
      r(:,i,:,:,:) = r(:,i,:,:,:) - s(:,i,:,:,:)
     enddo ! i
     itercount = itercount + 1
     call vecdot(r,r,vecsq,MRT2)
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       write(unit=8,fmt="(a12,i9,es17.10)") "qmrgam5",itercount,vecsq(1)
      close(unit=8,status="keep")
     endif
     if (real(vecsq(1),KR)<resmax.and.itercount>=itermin) exit maindo
    enddo maindo

 end subroutine qmrgam5

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine gmres(rwdir,phi,x,nGMRES,resmax,itermin,itercount,LUcount,u, &
                  GeeGooinv,iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES(n) matrix inverter.
! Solves M*x=phi for the vector x.
! The notation of Saad, "Iterative Methods for Sparse Linear Systems",
! algorthms 6.9 and 6.11.  (This published book can also be found in the web.)
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   nGMRES is the "n" in GMRES(n) -- the number of iterations between restarts.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmres.
!   LUcount is the number of extra LU inversions used.

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in)     :: nGMRES, itermin, iflag, myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount, LUcount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

!    real(kind=KR),  dimension(6,ntotal,4,2,8,nmaxGMRES) :: v
    real(kind=KR),allocatable,dimension(:,:,:,:,:,:) :: v
    real(kind=KR),  dimension(6,nvhalf,4,2,8)           :: r, w
    real(kind=KR),  dimension(2,nmaxGMRES+1,nmaxGMRES)  :: h
    real(kind=KR),  dimension(2,nmaxGMRES,nmaxGMRES)    :: amini, alud
    real(kind=KR),  dimension(2,nmaxGMRES)              :: bmini, y
    integer(kind=KI)                          :: idag, i, j, k, l, icri, mGMRES
    real(kind=KR2), dimension(2)                        :: vecsq
    real(kind=KR2)                                      :: beta
    real(kind=KR)                                       :: rowsign
    integer(kind=KI), dimension(nmaxGMRES)              :: indx

! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.

  allocate(v(6,ntotal,4,2,8,nmaxGMRES))


! Initializations.
    itercount = 0
    LUcount = 0
    idag = 0

! This outer loop includes the restarting of GMRES (Saad's algorithm 6.11).
    maindo: do 

! Initialize GMRES at each restart (Saad's algorithm 6.9).
     call Hdbletm(w,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                  myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     do i = 1,nvhalf
      r(:,i,:,:,:) = phi(:,i,:,:,:) - w(:,i,:,:,:)
     enddo ! i
     call vecdot(r,r,vecsq,MRT2)
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       write(unit=8,fmt="(a12,i9,es17.10)") "gmres",itercount,vecsq(1)
      close(unit=8,status="keep")
     endif
     if (real(vecsq(1),KR)<resmax.and.itercount>=itermin) exit maindo
     beta = sqrt(vecsq(1))
     do i = 1,nvhalf
      v(:,i,:,:,:,1) = r(:,i,:,:,:)/beta
     enddo ! i

! Generate the Arnoldi basis and the matrix Hbar (Saad's algorithm 6.9).
     h = 0.0_KR
     j = 0
     jdo: do
      j = j + 1
      call Hdbletm(w,u,GeeGooinv,v(:,:,:,:,:,j),idag,coact,kappa,iflag,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      do i = 1,j
       call vecdot(v(:,:,:,:,:,i),w,h(:,i,j),MRT2)
       do k = 1,nvhalf
        do icri = 1,5,2 ! 6=nri*nc
         w(icri  ,k,:,:,:) = w(icri  ,k,:,:,:) - h(1,i,j)*v(icri  ,k,:,:,:,i) &
                                               + h(2,i,j)*v(icri+1,k,:,:,:,i)
         w(icri+1,k,:,:,:) = w(icri+1,k,:,:,:) - h(1,i,j)*v(icri+1,k,:,:,:,i) &
                                               - h(2,i,j)*v(icri  ,k,:,:,:,i)
        enddo ! icri
       enddo ! k
      enddo ! i
      call vecdot(w,w,vecsq,MRT2)
      if (vecsq(1)==0.0_KR2) then
       exit jdo
      else
       h(1,j+1,j) = sqrt(vecsq(1))
      endif
      do k = 1,nvhalf
       v(:,k,:,:,:,j+1) = w(:,k,:,:,:)/h(1,j+1,j)
      enddo ! k
      if (j==nGMRES) exit jdo
     enddo jdo
     mGMRES = j

! Compute the vector y that minimizes the length of vector beta*e1 - Hbar*y.
!1:define a corresponding linear algebra problem: matrix=amini and vector=bmini.
     amini = 0.0_KR
     do l = 1,mGMRES
      do k = 1,mGMRES
       do j = 1,mGMRES+1
        amini(1,k,l) = amini(1,k,l) + h(1,j,k)*h(1,j,l) + h(2,j,k)*h(2,j,l)
        amini(2,k,l) = amini(2,k,l) + h(1,j,k)*h(2,j,l) - h(2,j,k)*h(1,j,l)
       enddo ! j
      enddo ! k
      bmini(1,l) = beta*h(1,1,l)
      bmini(2,l) = -beta*h(2,1,l)
     enddo ! l
     y = bmini
!2:copy amini into alud, since ludcmp will destroy its input matrix.
     alud = amini
!3:decompose the matrix alud into LU form.
     call ludcmp(alud,mGMRES,indx,rowsign)
!4:solve the linear algebra problem.
     call lubksb(alud,mGMRES,indx,y)
     call mprove(amini,alud,mGMRES,indx,bmini,y,LUcount)

! Update the approximate solution to the original linear algebra problem.
     do j = 1,mGMRES
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        x(icri  ,k,:,:,:) = x(icri  ,k,:,:,:) + y(1,j)*v(icri  ,k,:,:,:,j) &
                                              - y(2,j)*v(icri+1,k,:,:,:,j)
        x(icri+1,k,:,:,:) = x(icri+1,k,:,:,:) + y(1,j)*v(icri+1,k,:,:,:,j) &
                                              + y(2,j)*v(icri  ,k,:,:,:,j)
       enddo ! k
      enddo ! icri
     enddo ! j

! Record the conclusion of this iteration and begin a new one.
     itercount = itercount + 1
    enddo maindo

    deallocate(v)
 end subroutine gmres

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine gmresdr5EIG(rwdir,be,bo,xe,xo,GMRES,rtol,itermin,itercount,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2, idag
    real(kind=KR),    intent(in)                          :: rtol
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j,  jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, icri, m, k,kk, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ind
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv,rninit,rn,vn, vnf,xrn,resii,htest1,htest2
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag,dabbs
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w
    real(kind=KR2),   dimension(2,201,201)                  :: orthoge,orthog
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv,d
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: ff, punty
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: dd,th, tmpVec
    real(kind=KR2),   dimension(2,kmaxGMRES)                :: rho
    real(kind=KR2),   dimension(kmaxGMRES)                  :: rna, sita
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,gon,rr,gondag
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: h, h2, h3, hprint, hh,g,gg,greal, &
                                                               tmpmat,hnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             ::re,ro,xte,xto,xre,xro
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: worke,worko
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpE, tmpE2
    real(kind=KR2),   dimension(6)              :: vte, vto
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpO, tmpO2
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpOO, tmpEE,xtmpEE,xtmpOO,xeTMP,xoTMP !BS
    real(kind=KR2), dimension(2) :: tmp1,tmp2,tmp3,tmp4,xtmp1,xtmp2,dote,doto!BS
    real(kind=KR2), dimension(2,nmaxGMRES)      :: vnj
   
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htempe,htmpo
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempe,wve,fe,wxe
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempo,wvo,fo,wxo!BS
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow,gblclr 
    integer(kind=KI) :: didmaindo, exists
     if (myid==0) then
 print *,'just entered gmresdr5EIG'
 endif

if (myid == 0) then
   print *, 'the value of kappa in gmresdr5EIG is = ',kappa(1)
endif

        
  ! initialize variables
  m = GMRES(1)  ! size of Krylov Subspace
  k = GMRES(2)  ! number of eigenvector/values to deflate
  !m = 10 ! change back to above after testing
  !k = 5
  
  !idag = 0      ! flag to do M*x NOT Mdag*x

!if (myid ==0) then
!print *, 'rtol =',rtol
!endif




  ! ignore initial guess and use 0
  xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2
  xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

  icycle = 1    ! initialize cycle counter
  j = 1         ! intiialize iteration counter
  h = 0.0_KR2   ! initialize h matrix

  ! calculate inital residual by ignoring intial guess and use x=0
  ! r = b - Ax = b - A*0 = b
  re = be
  ro = bo
  call vecdot(re,re,tmp1,MRT2)
  call vecdot(ro,ro,tmp2,MRT2)
  rninit = tmp1 + tmp2
  rninit(1) = sqrt(rninit(1))
  rninit(2) = 0.0_KR2
  rn = rninit
  vn = rn

xrn = rn

  ! first vector
  ve(:6,:ntotal,:4,:2,:8,1) = (1.0_KR2 / vn(1)) * re(:6,:ntotal,:4,:2,:8)
  vo(:6,:ntotal,:4,:2,:8,1) = (1.0_KR2 / vn(1)) * ro(:6,:ntotal,:4,:2,:8)

  ! right hand side to future least square problem
  c = 0.0_KR2
  c(1,1) = vn(1)

!BS      if (myid==0) then
 !print *,'starting while loop'
 !endif

  ! begin cycles
 ! do while((rn(1)/rninit(1)) > rtol .AND. icycle <= kcyclim)
  do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 2000))!becomes zero in 29 cycles
 ! do while(((rn(1)/rninit(1)) > 1e-9) .AND. (icycle <= 15))
  !do while((resi > 1e-9) .AND. (icycle <= 15))
 ! do while(icycle <=150)! to make sure if it works
    ! begin gmres(m) iterations
    do while(j <= m)
      ! First build H_oe * v_e.
     
!!!!  Added line by BS. ordering gamama5 accordingly !!!!!
      do ieo=1,2
        do ibleo=1,8
          call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO, 5,ieo,ibleo)
        enddo
      enddo
!!!!  Added until this line !!!!     
     
!!!! commented out by BS  !!!!
!  gblclr = 2   
 !call Hsingle(wvo,u,ve(:,:,:,:,:,j),idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
      !             nms,lvbc,ib,lbd,iblv,MRT)
      ! Next build H_eo * v_o.
     ! gblclr = 1
     ! call Hsingle(wve,u,vo(:,:,:,:,:,j),idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
      !             nms,lvbc,ib,lbd,iblv,MRT)
! tmpE = ve(:,:,:,:,:,j) - kappa(1) * wve
 !     tmpO = vo(:,:,:,:,:,j) - kappa(1) * wvo
! ****BS comment out until this line *******

     gblclr = 2

   call Hsingle(wvo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
      ! Next build H_eo * v_o.
      gblclr = 1
      call Hsingle(wve,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)


      wve = tmpEE - kappa(1) * wve
      wvo = tmpOO - kappa(1) * wvo
! BS ********start comment out
      ! multiply by gama5 here
      ! BS 10/4/2015 displaced gamma 5 multiplication  
     ! do ieo=1,2
      !  do ibleo=1,8
       !   call gammamult(tmpE, tmpO, wve, wvo, 5,ieo,ibleo)
       ! enddo
     ! enddo
  ! BS ****end comment out

      call vecdot(wve,wve,tmp1,MRT2)
      call vecdot(wvo,wvo,tmp2,MRT2)
      vnf = tmp1 + tmp2
      vnf(1) = sqrt(vnf(1))
      vnf(2) = 0.0_KR2

      do i=1,j
        call vecdot(ve(:,:,:,:,:,i),wve,tmp1,MRT2)
        call vecdot(vo(:,:,:,:,:,i),wvo,tmp2,MRT2)
        h(:,i,j) = tmp1(:) + tmp2(:) !commented out TW 2/26/18
        !h(1,i,j) = tmp1(1) + tmp2(1)
        !h(2,i,j) = tmp1(2) + tmp2(2)
                 
        do icri = 1,5,2 ! 6=nri*nc
          do jj = 1,nvhalf
            wve(icri  ,jj,:,:,:) = wve(icri  ,jj,:,:,:) &
                                 - h(1,i,j)*ve(icri  ,jj,:,:,:,i) &
                                 + h(2,i,j)*ve(icri+1,jj,:,:,:,i)
            wvo(icri  ,jj,:,:,:) = wvo(icri  ,jj,:,:,:) &
                                 - h(1,i,j)*vo(icri  ,jj,:,:,:,i) &
                                 + h(2,i,j)*vo(icri+1,jj,:,:,:,i)
            wve(icri+1,jj,:,:,:) = wve(icri+1,jj,:,:,:) &
                                 - h(2,i,j)*ve(icri  ,jj,:,:,:,i) &
                                 - h(1,i,j)*ve(icri+1,jj,:,:,:,i)
            wvo(icri+1,jj,:,:,:) = wvo(icri+1,jj,:,:,:) &
                                 - h(2,i,j)*vo(icri  ,jj,:,:,:,i) &
                                 - h(1,i,j)*vo(icri+1,jj,:,:,:,i)
          enddo
        enddo
      enddo

      call vecdot(wve,wve,tmp1,MRT2)
      call vecdot(wvo,wvo,tmp2,MRT2)
      vn = tmp1 + tmp2
      vn(1) = sqrt(vn(1))
      vn(2) = 0.0_KR2

      vnj(1,j) = vn(1)

      ! --- reorthogonalization section ---
      if (vn(1) < (1.1_KR * vnf(1))) then
        do i=1,j
          call vecdot(ve(:,:,:,:,:,i),wve,tmp2,MRT2)
          call vecdot(vo(:,:,:,:,:,i),wvo,tmp3,MRT2)
          tmp1 = tmp2 + tmp3 !commented out for testing TW 2/26/18
          ! tmp1(1) = tmp2(1) + tmp3(1) !added TW same date
           !tmp1(2) = tmp2(2) + tmp3(2) !added TW same date
          do icri = 1,5,2 ! 6=nri*nc
            do jj = 1,nvhalf
              wve(icri  ,jj,:,:,:) = wve(icri  ,jj,:,:,:) &
                                   - tmp1(1)*ve(icri  ,jj,:,:,:,i) &
                                   + tmp1(2)*ve(icri+1,jj,:,:,:,i)
              wvo(icri  ,jj,:,:,:) = wvo(icri  ,jj,:,:,:) &
                                   - tmp1(1)*vo(icri  ,jj,:,:,:,i) &
                                   + tmp1(2)*vo(icri+1,jj,:,:,:,i)
              wve(icri+1,jj,:,:,:) = wve(icri+1,jj,:,:,:) &
                                   - tmp1(2)*ve(icri  ,jj,:,:,:,i) &
                                   - tmp1(1)*ve(icri+1,jj,:,:,:,i)
              wvo(icri+1,jj,:,:,:) = wvo(icri+1,jj,:,:,:) &
                                   - tmp1(2)*vo(icri  ,jj,:,:,:,i) &
                                   - tmp1(1)*vo(icri+1,jj,:,:,:,i)
            enddo
          enddo
          h(:,i,j) = h(:,i,j) + tmp1(:) !commented out for testing TW 2/26/18
          ! h(1,i,j) = h(1,i,j) + tmp1(1) !added TW same date
          ! h(2,i,j) = h(2,i,j) + tmp1(2) !added TW same date
        enddo
        call vecdot(wve,wve,tmp1,MRT2)
        call vecdot(wvo,wvo,tmp2,MRT2)
        vn = tmp1 + tmp2
        vn(1) = sqrt(vn(1))
        vn(2) = 0.0_KR2
        vnj(2,j) = vn(1) !for output of vnh.dat
      endif
      ! --- --- --- --- --- --- --- --- ---

      h(:,j+1,j) = vn !changed from vn to vn(1) TW 2/26/18
      ve(:6,:ntotal,:4,:2,:8,j+1) = (1.0_KR2 / h(1,j+1,j)) * wve(:6,:ntotal,:4,:2,:8)
      vo(:6,:ntotal,:4,:2,:8,j+1) = (1.0_KR2 / h(1,j+1,j)) * wvo(:6,:ntotal,:4,:2,:8)

      j = j + 1
    enddo

 if (.false.) then
  if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !print out orthog. of basis vectors for the krylov subspace TW 2/18/18                     
   do i=1,m+1
    do j=1,m+1

      dote=0.0_KR2
      doto=0.0_KR2

     call vecdot(ve(:,:,:,:,:,i),ve(:,:,:,:,:,j),dote,MRT2)
     call vecdot(vo(:,:,:,:,:,i),vo(:,:,:,:,:,j),doto,MRT2)


     orthog(1,i,j)=dote(1) + doto(1)
     orthog(2,i,j)=dote(2) + doto(2)

   enddo!j
  enddo!i


        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"orthog.dat", exist=exists)
            if (.not. exists) then

               open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="new",&
               action="write",form="formatted")
               close(unit=48,status="keep")
            endif
       endif


     do i=1,m+1
       do j=1,m+1
            if (myid==0) then
                 open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=48,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12,es22.12)")"icycle=",icycle,"V'V(",i,",",j,")=",orthog(1,i,j),orthog(2,i,j)
                 close(unit=48,status="keep")
            endif
      enddo
    enddo



  endif !icycle
endif!true or false

 if (.false.) then
  if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then 
        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"vn.dat", exist=exists)
            if (.not. exists) then

               open(unit=58,file=trim(rwdir(myid+1))//"vn.dat",status="new",&
               action="write",form="formatted")
               close(unit=58,status="keep")
            endif
       endif
     do i = 1,j
       if (myid==0) then
              open(unit=58,file=trim(rwdir(myid+1))//"vn.dat",status="old",action="write",&
              form="formatted",position="append")
              write(unit=58,fmt="(a9,i7,a10,i3,a10,es22.12,a10,es22.12)")"icycle=",icycle,"j=",i,"vn=",vnj(1,i),"vn after reorthog =",vnj(2,i)
              close(unit=58,status="keep")
            endif
      enddo
  endif !icycle
endif !true/false

if (.true.) then
  if (icycle==1 .OR. icycle==2 .OR. icycle==10 .OR. icycle==50 .OR. icycle==100 .OR. icycle==200 .OR. icycle==250 .OR. icycle==1000 .OR. icycle==2000) then
        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"htest.dat", exist=exists)
            if (.not. exists) then

               open(unit=98,file=trim(rwdir(myid+1))//"htest.dat",status="new",&
               action="write",form="formatted")
               close(unit=98,status="keep")
            endif
       endif
     do i = 1,m
      do kk = 1,m
       if (myid==0) then
              open(unit=98,file=trim(rwdir(myid+1))//"htest.dat",status="old",action="write",&
              form="formatted",position="append")
              write(unit=98,fmt="(es22.12,es22.12)")h(1,kk,i),h(2,kk,i)
              close(unit=98,status="keep")
            endif
      enddo !kk
     enddo !i
  endif !icycle
endif !true/false

    call leastsquaresLAPACK(h,m+1,m,c,d,myid)
    call matvecmult(h,m+1,m,d,m,srv)
    srv(:,:m+1) = c(:,:m+1) - srv(:,:m+1)


    ! Setup and sovle linear equations problem
    ! x(:) = x(:) + v(:,1:m)*d(1:m)
    do i=1,m
      do icri = 1,5,2 ! 6=nri*nc
        do jj = 1,nvhalf
          xe(icri  ,jj,:,:,:) = xe(icri  ,jj,:,:,:) &
                              + d(1,i)*ve(icri  ,jj,:,:,:,i) &
                              - d(2,i)*ve(icri+1,jj,:,:,:,i)
          xo(icri  ,jj,:,:,:) = xo(icri  ,jj,:,:,:) &
                              + d(1,i)*vo(icri  ,jj,:,:,:,i) &
                              - d(2,i)*vo(icri+1,jj,:,:,:,i)
          xe(icri+1,jj,:,:,:) = xe(icri+1,jj,:,:,:) &
                              + d(2,i)*ve(icri  ,jj,:,:,:,i) &
                              + d(1,i)*ve(icri+1,jj,:,:,:,i)
          xo(icri+1,jj,:,:,:) = xo(icri+1,jj,:,:,:) &
                              + d(2,i)*vo(icri  ,jj,:,:,:,i) &
                              + d(1,i)*vo(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo


!BS  4/9/016 to calculate residual using directmethod


      do ieo=1,2
        do ibleo=1,8
          call gammamult( xe(:,:,:,:,:), xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
        enddo
      enddo




      gblclr = 2

      call Hsingle(wxo,u,xtmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
      ! Next build H_eo * v_o.
      gblclr = 1
      call Hsingle(wxe,u,xtmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)


      wxe = xtmpEE - kappa(1) * wxe
      wxo = xtmpOO - kappa(1) * wxo


      xre = wxe - be
      xro = wxo - bo



      call vecdot(xre,xre,xtmp1,MRT2)
      call vecdot(xro,xro,xtmp2,MRT2)

      xrn = xtmp1+xtmp2
      xrn(1)=sqrt(xrn(1))
      xrn(2)=0.0_KR2



    ! r = v(:,1:m+1)*srv(1:m+1)
    re = 0.0_KR2
    ro = 0.0_KR2
    do icri = 1,5,2 ! 6=nri*nc
      do jj = 1,nvhalf
        do i=1,m+1
          re(icri  ,jj,:,:,:) = re(icri,jj,:,:,:)   &
                              + srv(1,i)*ve(icri  ,jj,:,:,:,i) &
                              - srv(2,i)*ve(icri+1,jj,:,:,:,i)
          ro(icri  ,jj,:,:,:) = ro(icri,jj,:,:,:)   &
                              + srv(1,i)*vo(icri  ,jj,:,:,:,i) &
                              - srv(2,i)*vo(icri+1,jj,:,:,:,i)
          re(icri+1,jj,:,:,:) = re(icri+1,jj,:,:,:)  &
                              + srv(2,i)*ve(icri  ,jj,:,:,:,i) &
                              + srv(1,i)*ve(icri+1,jj,:,:,:,i)
          ro(icri+1,jj,:,:,:) = ro(icri+1,jj,:,:,:)  &
                              + srv(2,i)*vo(icri  ,jj,:,:,:,i) &
                              + srv(1,i)*vo(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo
   

    call vecdot(re,re,tmp1,MRT2)
    call vecdot(ro,ro,tmp2,MRT2)
    rn = tmp1 + tmp2
    rn(1) = sqrt(rn(1))
    rn(2) = 0.0_KR2

    !do i = 1,m
     ! do kk = 1,m
      !  h(1,kk,i) = (1.0_KR2 / 2.0_KR2) * (h(1,kk,i) + h(1,i,kk))
       ! h(2,kk,i) = (1.0_KR2 / 2.0_KR2) * (h(2,kk,i) - h(2,i,kk))
      !enddo
    !enddo


    ! prepare for next cycle
    hh(:2,1:m,1:m) = h(:2,1:m,1:m)
    do i=1,m
      do jj=1,m
        hcht(1,i,jj) = h(1,jj,i)
        hcht(2,i,jj) = -1.0_KR2 * h(2,jj,i)
      enddo
    enddo

    ff=0.0_KR2
    ff(1,m,1)=1.0_KR2

    call linearsolver(m,1,hcht,ipiv,ff)

    do i=1,m
      hh(1,i,m) = hh(1,i,m)-h(2,m+1,m)**2*ff(1,i,1)-2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(2,i,1)+h(1,m+1,m)**2*ff(1,i,1)
      hh(2,i,m) = hh(2,i,m)-h(2,m+1,m)**2*ff(2,i,1)+2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(1,i,1)+h(1,m+1,m)**2*ff(2,i,1)
    enddo

! BS 1/21/2016 for testing purpose remove after done
!m=3
     
!hh(1,1) = 1
!hh(2,1) = 0
!hh(3,1) = 0
!hh(1,2) = 0
!hh(2,2) = 2
!hh(3,2) = 0
!hh(1,3) = 0
!hh(2,3) = 0
!hh(3,3) = 4



!allocate(hh(3,3))



    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg) !BS 1/19/2016
   ! call eigmodesLAPACK(hh,m,1,dd,gg)

    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)

    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo

    ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      call matvecmult(h,m,m,g(:,:,i),m,tmpVec) !I think this is redundant -TW
      do jj=1,m
        call cmpxmult(rho(:,i),g(:,jj,i),tmp1)
        tmpVec(:,jj) = tmpVec(:,jj) - tmp1
      enddo
      call vecdagvec(tmpVec,m,tmpVec,m,tmp1)
      tmp1(1) = sqrt(tmp1(1))
      tmp1(2) = 0.0_KR2

      rna(i) = sqrt((tmp1(1)*tmp1(1)) + (h(1,m+1,m)*h(1,m+1,m) + h(2,m+1,m)*h(2,m+1,m)) &
                * (g(1,m,i)*g(1,m,i) + g(2,m,i)*g(2,m,i)))
    enddo
!BS changed 


    call sort(rna,k,ind)

    do i = 1,k
       sita(i) = rna(ind(i))
    enddo

    do i =1,k
       rna(i) = sita(i)
    enddo

if (.true.) then

    do i=1,k  !BS 5/4/2016

        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"residual.dat", exist=exists)
            if (.not. exists) then
               open(unit=73,file=trim(rwdir(myid+1))//"residual.dat",status="new",&
               action="write",form="formatted")
               close(unit=73,status="keep")
            endif
       endif



       if (myid==0) then
         open(unit=73,file=trim(rwdir(myid+1))//"residual.dat",status="old",action="write",&
            form="formatted",position="append")
            write(unit=73,fmt="(i7,a6,i7,a6,es19.12)") icycle,"   ",i," ",rna(i)
            close(unit=73,status="keep")
      endif

   enddo !i

endif !true or false


    do i=1,k
      gg(:,:,i) = g(:,:,ind(i))
    enddo

    do i=1,k
      gg(:,m+1,i) = 0.0_KR2
    enddo
!Chris

    beta(1) = h(1,m+1,m)
    beta(2) = h(2,m+1,m)

    do j = 1,m
        punty(1,j,1) = -beta(1)*ff(1,j,1)+beta(2)*ff(2,j,1)
        punty(2,j,1) = -beta(2)*ff(1,j,1)-beta(1)*ff(2,j,1)
    enddo

    do j = 1,m
        gg(:,j,k+1) = punty(:,j,1)
    enddo
    gg(1,m+1,k+1) = 1.0_KR2

!end Chris

   ! gg(:,:,k+1) = srv

    call qrfactorizationLAPACK(gg,m+1,k+1,gon,rr,myid)

    do i=1,m+1
      do jj=1,m+1
        gondag(1,i,jj) = gon(1,jj,i)
        gondag(2,i,jj) =-1.0_KR2 * gon(2,jj,i)
      enddo
    enddo

    ! hcnew = gon'*h*gon(1:m,1:k) 
    call matmult(gondag,k+1,m+1,h,m+1,m,tmpmat)
    call matmult(tmpmat,k+1,m,gon,m,k,hcnew)

    h(:,k+1,:m) = 0.0_KR2

    i = 1
    do while ((i <= k) .AND. (rna(i) < 1E-12)) !added i<=k to match matlab TW 8218     
      hcnew(:,i+1:k+1,i) = 0.0_KR2! BS travis suggests to be consistenst with matlab
      i = i + 1
    enddo

    ! form right eigenvectors; evector is in shift module
    do j=1,k
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vte(:) = 0.0_KR2
              vto(:) = 0.0_KR2
              do kk=1,m
                do icri = 1,5,2
                  vte(icri)   = vte(icri) &
                        + ve(icri  ,i,id,ieo,ibleo,kk)*gg(1,kk,j) &
                        - ve(icri+1,i,id,ieo,ibleo,kk)*gg(2,kk,j) 
                  vto(icri)   = vto(icri) &
                        + vo(icri  ,i,id,ieo,ibleo,kk)*gg(1,kk,j) &
                        - vo(icri+1,i,id,ieo,ibleo,kk)*gg(2,kk,j) 
                  vte(icri+1) = vte(icri+1) &
                        + ve(icri  ,i,id,ieo,ibleo,kk)*gg(2,kk,j) &
                        + ve(icri+1,i,id,ieo,ibleo,kk)*gg(1,kk,j) 
                  vto(icri+1) = vto(icri+1) &
                        + vo(icri  ,i,id,ieo,ibleo,kk)*gg(2,kk,j) &
                        + vo(icri+1,i,id,ieo,ibleo,kk)*gg(1,kk,j) 
                enddo
              enddo
              evectore(:,i,id,ieo,ibleo,j) = vte(:)
              evectoro(:,i,id,ieo,ibleo,j) = vto(:)
            enddo
          enddo
        enddo
      enddo
      evalue(:,j) = th(:,ind(j))
            !if (myid==0) then
             !  write(*,*) 'evalue5EIG:',j,'j:',evalue
            !endif

    enddo
       !     do while(j <= m)
 !BS          if (myid==0) then
    !BS           write(*,*) 'evalue5EIG:',j,'j:',evalue
       !BS    endif
       ! enddo
!!***** added TW 1/23/18 to get residuals of lowest lying eigenpairs*******
    if (.true.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,k
           do ieo=1,2           
            do ibleo=1,8
             call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
            enddo
           enddo

         gblclr = 2
         call Hsingle(tmpOO,u,tmpE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
         ! Next build H_eo * v_o.
         gblclr = 1
         call Hsingle(tmpEE,u,tmpO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
         xeTMP = tmpE - kappa(1) * tmpEE
         xoTMP = tmpO - kappa(1) * tmpOO

         tmpE = xeTMP - (evalue(1,i)*evectore(:,:,:,:,:,i))
         tmpO = xoTMP - (evalue(1,i)*evectoro(:,:,:,:,:,i))

         call vecdot(tmpE,tmpE,tmp1,MRT2)
         call vecdot(tmpO,tmpO,tmp2,MRT2)
         tmp1 = tmp1 + tmp2
         call cmpxsqrt(tmp1, tmp2)

         if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"loweigresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=88,file=trim(rwdir(myid+1))//"loweigresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=88,status="keep")
           endif
         endif

         if (myid==0) then
           open(unit=88,file=trim(rwdir(myid+1))//"loweigresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=88,fmt="(i4,a6,i7,a6,es19.12,a6,es19.12,a6,es19.12)")icycle," ", i," ",evalue(1,i)," ",evalue(2,i)," ",tmp2(1)
           close(unit=88,status="keep")
         endif


      enddo !i
     endif !icycle
    endif !.true./.false

!******end add TW 1/23/18********************

    do i=1,k
      do jj=1,k+1
        h(:,jj,i) = hcnew(:,jj,i)
      enddo
    enddo

    call matvecmult(gondag,k+1,m+1,srv,m+1,c)

    do i=k+2,m+1
      c(:,i) = 0.0_KR2
    enddo

! The section is here Dr. Wilcox -- Travis 8/11/18
    do jj=1,k+1
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vte(:) = 0.0_KR2
              vto(:) = 0.0_KR2
              do kk=1,m+1
                do icri = 1,5,2
                  vte(icri)   = vte(icri) &
                      + ve(icri  ,i,id,ieo,ibleo,kk)*gon(1,kk,jj) &
                      - ve(icri+1,i,id,ieo,ibleo,kk)*gon(2,kk,jj) 
                  vto(icri)   = vto(icri) &
                      + vo(icri  ,i,id,ieo,ibleo,kk)*gon(1,kk,jj) &
                      - vo(icri+1,i,id,ieo,ibleo,kk)*gon(2,kk,jj) 
                  vte(icri+1) = vte(icri+1) &
                      + ve(icri  ,i,id,ieo,ibleo,kk)*gon(2,kk,jj) &
                      + ve(icri+1,i,id,ieo,ibleo,kk)*gon(1,kk,jj) 
                  vto(icri+1) = vto(icri+1) &
                      + vo(icri  ,i,id,ieo,ibleo,kk)*gon(2,kk,jj) &
                      + vo(icri+1,i,id,ieo,ibleo,kk)*gon(1,kk,jj) 
                enddo
              enddo
              worke(:,i,id,ieo,ibleo,jj) = vte(:)
              worko(:,i,id,ieo,ibleo,jj) = vto(:)
            enddo
          enddo
        enddo
      enddo
    enddo

    do i=1,k+1
      ve(:,:,:,:,:,i) = worke(:,:,:,:,:,i)
      vo(:,:,:,:,:,i) = worko(:,:,:,:,:,i)
    enddo

    do i=1,k
      call vecdot(ve(:,:,:,:,:,i),ve(:,:,:,:,:,k+1),tmp2,MRT2)
      call vecdot(vo(:,:,:,:,:,i),vo(:,:,:,:,:,k+1),tmp3,MRT2)
      tmp1 = tmp2 + tmp3

      do icri = 1,5,2 ! 6=nri*nc
        !do jj = 1,ntotal
        do jj = 1,nvhalf
          ve(icri  ,jj,:,:,:,k+1) = ve(icri  ,jj,:,:,:,k+1) &
                                  - tmp1(1)*ve(icri  ,jj,:,:,:,i) &
                                  + tmp1(2)*ve(icri+1,jj,:,:,:,i)
          vo(icri  ,jj,:,:,:,k+1) = vo(icri  ,jj,:,:,:,k+1) &
                                  - tmp1(1)*vo(icri  ,jj,:,:,:,i) &
                                  + tmp1(2)*vo(icri+1,jj,:,:,:,i)
          ve(icri+1,jj,:,:,:,k+1) = ve(icri+1,jj,:,:,:,k+1) &
                                  - tmp1(2)*ve(icri  ,jj,:,:,:,i) &
                                  - tmp1(1)*ve(icri+1,jj,:,:,:,i)
          vo(icri+1,jj,:,:,:,k+1) = vo(icri+1,jj,:,:,:,k+1) &
                                  - tmp1(2)*vo(icri  ,jj,:,:,:,i) &
                                  - tmp1(1)*vo(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo


    call vecdot(ve(:,:,:,:,:,k+1),ve(:,:,:,:,:,k+1),tmp2,MRT2)
    call vecdot(vo(:,:,:,:,:,k+1),vo(:,:,:,:,:,k+1),tmp3,MRT2)
    tmp1 = tmp2 + tmp3
    tmp1(1) = sqrt(tmp1(1))

    do jj = 1,nvhalf
      ve(:,jj,:,:,:,k+1) = (1.0_KR2 / tmp1(1))*ve(:,jj,:,:,:,k+1)
      vo(:,jj,:,:,:,k+1) = (1.0_KR2 / tmp1(1))*vo(:,jj,:,:,:,k+1)
    enddo

    !if (myid==0) then
    !  write(*,*) 'cycle:',icycle,'resnorm:',rn(1)/rninit(1)
    !endif
!     if (myid==0) then
 !print *,'about to end gmresdr5EIG'
 !endif
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmrs5EIG",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmrs5EIG",icycle,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)

       !BS   write(unit=8,fmt="(a12,i9,es17.10)") "igmresdr5EIG",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif

    j = k+1
    icycle = icycle + 1
  enddo
      
          if (myid==0) then !BS 
            do i =1,k 
                 write(*,fmt="(a12,i6,f19.11,a5,f19.11)") "evalue5EIG:",i,evalue(1,i),"",evalue(2,i)
 !BS               print * , 'evalue5EIG:',evalue
            enddo 
         endif
            

          if (myid==0.AND.(icycle==150)) then !BS
             print * ,'rn/rninit = ',rn(1)/rninit(1)
          endif
          
  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmrs5EIG",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmrs5EIG",icycle-1,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)
 
 !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr5EG",icycle-1,rn(1)/rninit(1)
    close(unit=8,status="keep")
  endif
end subroutine gmresdr5EIG

 subroutine gmresdrEIG(rwdir,b,x,GMRES,rtol,itermin,itercount,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: b
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2, idag
    real(kind=KR),    intent(in)                          :: rtol
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j,  jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, icri, m, k,kk, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ind
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv,rninit,rn,vn, vnf,xrn,rina
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag,dabbs
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w
    real(kind=KR),       dimension(2,201,201)                     :: jamuna!BS
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv,d
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: ff
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: dd,th,tmpVec,resi
    real(kind=KR2),   dimension(2,kmaxGMRES)                :: rho
    real(kind=KR2),   dimension(kmaxGMRES)                  :: rna,sita
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,gon,rr,gondag
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: h, h2, h3, hprint, hh,g,gg,greal, &
                                                               tmpmat,hnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: r, xt
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: work
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpE, tmpE2
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpO, tmpO2
    real(kind=KR), dimension(6,ntotal,4,2,8,nmaxGMRES+1)  ::muna,tina,matmul,matmully
    real(kind=KR), dimension(6,ntotal,4,2,8)  ::getemp,wv,f,xtmp,mina 
   ! real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(6)              :: vt
    real(kind=KR2), dimension(2) :: tmp1,tmp2,tmp3,tmp4,broda
    real(kind=KR2)   :: temp11,temp22,temp33,normmatmully
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htemp
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr
    integer(kind=KI) :: didmaindo,pickle
    integer(kind=KI) :: exists,vena  ! initialize variables
  m = GMRES(1)  ! size of Krylov Subspace
  k = GMRES(2)  ! number of eigenvector/values to deflate

 ! print *,'k=',k
  !idag = 0      ! flag to do M*x NOT Mdag*x

  ! ignore initial guess and use 0
  x(:6,:ntotal,:4,:2,:8) = 0.0_KR2

  icycle = 1    ! initialize cycle counter
  j = 1         ! intiialize iteration counter
  h = 0.0_KR2   ! initialize h matrix

  ! calculate inital residual by ignoring intial guess and use x=0
  ! r = Ax - b = A*0 - b = b
  r = b
  call vecdot(r,r,rninit,MRT2)
  rninit(1) = sqrt(rninit(1))
  rninit(2) = 0.0_KR2
  rn = rninit
  vn = rn

  ! first vector
  v(:6,:ntotal,:4,:2,:8,1) = (1.0_KR2 / vn(1)) * r(:6,:ntotal,:4,:2,:8)

  ! right hand side to future least square problem
  c = 0.0_KR2
  c(1,1) = vn(1)

  if (myid==0) then
     print *,"gmresdrprinttest1"
  endif
  ! begin cycles
  do while(((rn(1)/rninit(1)) > 1e-150_KR2) .AND. (icycle <= 150))!becomes zero in 29 cycles
 ! do while(icycle <= 39)!BS 1/22/2016 to make sure enough evalue pass tolerance
    ! begin gmres(m) iterations
    do while(j <= m)
      call Hdbletm(wv,u,GeeGooinv,v(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      call vecdot(wv,wv,vnf,MRT2)
      vnf(1) = sqrt(vnf(1))
      vnf(2) = 0.0_KR2

      do i=1,j
        call vecdot(v(:,:,:,:,:,i),wv,h(:,i,j),MRT2)
        do icri = 1,5,2 ! 6=nri*nc
          do jj = 1,nvhalf
            wv(icri  ,jj,:,:,:) = wv(icri  ,jj,:,:,:) &
                                  - h(1,i,j)*v(icri  ,jj,:,:,:,i) &
                                  + h(2,i,j)*v(icri+1,jj,:,:,:,i)
            wv(icri+1,jj,:,:,:) = wv(icri+1,jj,:,:,:) &
                                  - h(2,i,j)*v(icri  ,jj,:,:,:,i) &
                                  - h(1,i,j)*v(icri+1,jj,:,:,:,i)
          enddo
        enddo
      enddo

      call vecdot(wv,wv,vn,MRT2)
      vn(1) = sqrt(vn(1))
      vn(2) = 0.0_KR2

      ! --- reorthogonalization section ---
      if (vn(1) < (1.1_KR * vnf(1))) then
        do i=1,j
          call vecdot(v(:,:,:,:,:,i),wv,tmp1,MRT2)
          do icri = 1,5,2 ! 6=nri*nc
            do jj = 1,nvhalf
              wv(icri  ,jj,:,:,:) = wv(icri  ,jj,:,:,:) &
                                  - tmp1(1)*v(icri  ,jj,:,:,:,i) &
                                  + tmp1(2)*v(icri+1,jj,:,:,:,i)
              wv(icri+1,jj,:,:,:) = wv(icri+1,jj,:,:,:) &
                                  - tmp1(2)*v(icri  ,jj,:,:,:,i) &
                                  - tmp1(1)*v(icri+1,jj,:,:,:,i)
            enddo
          enddo
          h(:,i,j) = h(:,i,j) + tmp1(:)
        enddo
        call vecdot(wv,wv,vn,MRT2)
        vn(1) = sqrt(vn(1))
        vn(2) = 0.0_KR2
      endif
      ! --- --- --- --- --- --- --- --- ---

      h(:,j+1,j) = vn
      v(:6,:ntotal,:4,:2,:8,j+1) = (1.0_KR2 / h(1,j+1,j)) * wv(:6,:ntotal,:4,:2,:8)

      j = j + 1
    enddo


!------------------************************----------------------------

!BS 5/12/2016 print jamuna to check orthogonality condition for all
!icycles,creates orthonormal.dat !look at scratch file,remember to set nps=1
!though

!   call vecdothigherrank (muna,muna,jamuna,m+1,MRT2)! vecdothigherrank was
!   supposed to give jamuna replacing the algorithms below
 
if (.false.) then
  do i=1,m+1
   do j=1,m+1

      broda=0.0_KR2

     call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,j),broda,MRT2)

     jamuna(1,i,j)=broda(1)
     jamuna(2,i,j)=broda(2)

   enddo!j
  enddo!i


        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"orthonormal.dat", exist=exists)
            if (.not. exists) then

               open(unit=47,file=trim(rwdir(myid+1))//"orthonormal.dat",status="new",&
               action="write",form="formatted")
               close(unit=47,status="keep")
            endif
       endif


     do i=1,m+1
       do j=1,m+1
            if (myid==0) then
                 open(unit=47,file=trim(rwdir(myid+1))//"orthonormal.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=47,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12,es22.12)") "icycle=",icycle,"V'V(",i,",",j,")=",jamuna(1,i,j),jamuna(2,i,j)
                 close(unit=47,status="keep")
            endif
      enddo
    enddo




endif!true or false



!BS 05/12/016 print of jamuna completed succesfully....change (.false.) to (.true.) if want to create
!orthonormal.dat

!------------------------*******************************---------------------------















!-------------------********************************************--------------------------

!BS 05/12/2016 creates orthogonali.dat to see ||AV(1:m)-V(1:(m+1)*h(m+1,m)||

   if (.true.) then
 !    if (icycle>=28) then

          muna(:6,:ntotal,:4,:2,:8,:(m+1)) =  v(:6,:ntotal,:4,:2,:8,:(m+1))
          do i=1,m
             call Hdbletm(wv,u,GeeGooinv,v(:,:,:,:,:,i),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                          myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             tina(:,:,:,:,:,i) = wv(:,:,:,:,:)
          enddo

          call matrixmultiplylike(muna,h,m,matmul,MRT2)

if (myid ==0) then
     print *,"tina(6,128,4,2,8,5) =", tina(6,128,4,2,8,5)
     print *,"tina(6,128,4,2,7,5) =", tina(6,128,4,2,7,5)
     print *,"tina(6,128,4,2,6,5) =", tina(6,128,4,2,6,5)
     print *,"tina(6,128,4,2,5,5) =", tina(6,128,4,2,5,5)
     print *,"tina(6,128,4,2,4,5) =", tina(6,128,4,2,4,5)
     print *,"tina(6,128,4,2,3,5) =", tina(6,128,4,2,3,5)
     print *,"tina(6,18,4,2,3,5) =",  tina(6,18,4,2,3,5)
     print *,"matmul(6,128,4,2,8,5) =", matmul(6,128,4,2,8,5)
     print *,"matmul(6,128,4,2,7,5) =", matmul(6,128,4,2,7,5)
     print *,"matmul(6,128,4,2,6,5) =", matmul(6,128,4,2,6,5)
     print *,"matmul(6,128,4,2,5,5) =", matmul(6,128,4,2,5,5)
     print *,"matmul(6,128,4,2,4,5) =", matmul(6,128,4,2,4,5)
     print *,"matmul(6,128,4,2,3,5) =", matmul(6,128,4,2,3,5)
     print *,"matmul(6,18,4,2,3,5) =",  matmul(6,18,4,2,3,5)
endif


         normmatmully=0.0_KR2

          do ibleo = 1,8
           do ieo = 1,2
            do id = 1,4
             do isite = 1,nvhalf
              do icri =1,6
               do j = 1,m

           
               matmully(icri,isite,id,ieo,ibleo,j) = matmul(icri,isite,id,ieo,ibleo,j)  - tina(icri,isite,id,ieo,ibleo,j)
               normmatmully = normmatmully + matmully(icri,isite,id,ieo,ibleo,j)**2
               enddo!j
              enddo!icri
             enddo!isite
            enddo!id
           enddo!ieo
          enddo!ibleo

        if (myid==0) then
           print *,"normmatmully = ",normmatmully
        endif   
if (.false.) then
        if (myid==0) then
            inquire(file=trim(rwdir(myid+1))//"orthogonali.dat", exist=exists)
            if (.not. exists) then

               open(unit=48,file=trim(rwdir(myid+1))//"orthogonali.dat",status="new",&
               action="write",form="formatted")
               close(unit=48,status="keep")
            endif
        endif



        if (myid==0) then
            open(unit=48,file=trim(rwdir(myid+1))//"orthogonali.dat",status="old",action="write",&
                 form="formatted",position="append")
                do ibleo = 1,8
                 do ieo = 1,2
                  do id = 1,4
                   do isite = 1,nvhalf
                    do icri =1,5,2
                     do j = 1,m

                        write(unit=48,fmt="(i6,i6,i6,i6,i6,i6,i6,es22.12,es22.12)") icycle,&
                              icri,isite,id,ieo,ibleo,j,matmully(icri,isite,id,ieo,ibleo,m),&
                               matmully(icri+1,isite,id,ieo,ibleo,m)
                     enddo!j
                    enddo!icri
                   enddo!isite
                  enddo!id
                 enddo!ieo
                enddo!ibleo
            close(unit=48,status="keep")
        endif! myid

endif
!     endif!if cycle==28
   endif! true or false




!BS 05/12/016 end creating orthogonali.dat

!----------------------------*************************-------------------------







    call leastsquaresLAPACK(h,m+1,m,c,d,myid)
    call matvecmult(h,m+1,m,d,m,srv)
    srv(:,:m+1) = c(:,:m+1) - srv(:,:m+1)

    ! x(:) = x(:) + v(:,1:m)*d(1:m)
    do i=1,m
      do icri = 1,5,2 ! 6=nri*nc
        do jj = 1,nvhalf
          x(icri  ,jj,:,:,:) = x(icri  ,jj,:,:,:) &
                                + d(1,i)*v(icri  ,jj,:,:,:,i) &
                                - d(2,i)*v(icri+1,jj,:,:,:,i)
          x(icri+1,jj,:,:,:) = x(icri+1,jj,:,:,:) &
                                + d(2,i)*v(icri  ,jj,:,:,:,i) &
                                + d(1,i)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo





!BS added this part for  linear equation residual 




      call Hdbletm(wv,u,GeeGooinv,x(:,:,:,:,:),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)


      xtmp = wv - b

     
      call vecdot(xtmp,xtmp,xrn,MRT2)
            xrn(1)=sqrt(xrn(1))
            xrn(2)=0.0_KR2

      


!BS residual for linear equation ends here











    ! r = v(:,1:m+1)*srv(1:m+1)
    r = 0.0_KR2
    do icri = 1,5,2 ! 6=nri*nc
      do jj = 1,nvhalf
        do i=1,m+1
          r(icri  ,jj,:,:,:) = r(icri,jj,:,:,:)   &
                              + srv(1,i)*v(icri  ,jj,:,:,:,i) &
                              - srv(2,i)*v(icri+1,jj,:,:,:,i)
          r(icri+1,jj,:,:,:) = r(icri+1,jj,:,:,:)  &
                              + srv(2,i)*v(icri  ,jj,:,:,:,i) &
                              + srv(1,i)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo

    call vecdot(r,r,rn,MRT2)
    rn(1) = sqrt(rn(1))
    rn(2) = 0.0_KR2

    ! prepare for next cycle
    hh(:2,1:m,1:m) = h(:2,1:m,1:m)
    do i=1,m
      do jj=1,m
        hcht(1,i,jj) = h(1,jj,i)
        hcht(2,i,jj) = -1.0_KR2 * h(2,jj,i)
      enddo
    enddo

    ff=0.0_KR2
    ff(1,m,1)=1.0_KR2

    call linearsolver(m,1,hcht,ipiv,ff)

    do i=1,m
      hh(1,i,m) = hh(1,i,m)-h(2,m+1,m)**2*ff(1,i,1)-2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(2,i,1)+h(1,m+1,m)**2*ff(1,i,1)
      hh(2,i,m) = hh(2,i,m)-h(2,m+1,m)**2*ff(2,i,1)+2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(1,i,1)+h(1,m+1,m)**2*ff(2,i,1)
    enddo

    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg)!BS 1/19/2016
    !call eigmodesLAPACK(hh,m,1,dd,gg)

    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)


    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo


   ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      do jj=1,m
        call cmpxmult(rho(:,i),g(:,jj,i),tmp1)
        tmpVec(:,jj) = tmpVec(:,jj) - tmp1
      enddo
      call vecdagvec(tmpVec,m,tmpVec,m,tmp1)
      tmp1(1) = sqrt(tmp1(1))
      tmp1(2) = 0.0_KR2

      rna(i) = sqrt((tmp1(1)*tmp1(1)) + (((h(1,m+1,m)*h(1,m+1,m)) + (h(2,m+1,m)*h(2,m+1,m))) &
                * ((g(1,m,i)*g(1,m,i)) + (g(2,m,i)*g(2,m,i)))))!BS changed + to *



    enddo!i



    call sort(rna,k,ind)

    do i=1,k
       sita(i)=rna(ind(i))!BS 5/4/2016
    enddo



    do i=1,k
       rna(i)=sita(i)!BS 5/4/2016
    enddo





    do i=1,k  !BS 5/4/2016

        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"residual.dat", exist=exists)
            if (.not. exists) then

               open(unit=33,file=trim(rwdir(myid+1))//"residual.dat",status="new",&
               action="write",form="formatted")
               close(unit=33,status="keep")
            endif
       endif




             if (myid==0) then
            open(unit=33,file=trim(rwdir(myid+1))//"residual.dat",status="old",action="write",&
            form="formatted",position="append")
                 write(unit=33,fmt="(i7,a6,i7,a6,es19.12)") icycle,"   ",i," ",rna(i)
           close(unit=33,status="keep")
          endif

   enddo !i

















    do i=1,k
      gg(:,:,i) = g(:,:,ind(i))
    enddo

    do i=1,k
      gg(:,m+1,i) = 0.0_KR2
    enddo

    gg(:,:,k+1) = srv

    call qrfactorizationLAPACK(gg,m+1,k+1,gon,rr,myid)

    do i=1,m+1
      do jj=1,m+1
        gondag(1,i,jj) = gon(1,jj,i)
        gondag(2,i,jj) =-1.0_KR2 * gon(2,jj,i)
      enddo
    enddo

    ! hcnew = gon'*h*gon(1:m,1:k) 
    call matmult(gondag,k+1,m+1,h,m+1,m,tmpmat)
    call matmult(tmpmat,k+1,m,gon,m,k,hcnew)

    h(:,k+1,:m) = 0.0_KR2

    i = 1
    do while (rna(i) < 1E-12)
      hnew(:,i+1:k+1,i) = 0.0_KR2
      i = i + 1
    enddo

    ! form right eigenvectors; evector is in shift module
    do j=1,k
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vt(:) = 0.0_KR2
              do kk=1,m
                do icri = 1,5,2
                  vt(icri)   = vt(icri) &
                        + v(icri  ,i,id,ieo,ibleo,kk)*gg(1,kk,j) &
                        - v(icri+1,i,id,ieo,ibleo,kk)*gg(2,kk,j) 
                  vt(icri+1) = vt(icri+1) &
                        + v(icri  ,i,id,ieo,ibleo,kk)*gg(2,kk,j) &
                        + v(icri+1,i,id,ieo,ibleo,kk)*gg(1,kk,j) 
                enddo
              enddo
              evector(:,i,id,ieo,ibleo,j) = vt(:)
            enddo
          enddo
        enddo
      enddo
      evalue(:,j) = th(:,ind(j))
! if (myid==0) then
 !              write(*,*) 'evalueEIG:',j,'j:',evalue
  !          endif

    enddo
       !BS do while(j <= m)
          !BS if (myid==0) then
            !BS   write(*,*) 'evalueEIG:',j,'j:',evalue
         !BS  endif
       !BS enddo

    do i=1,k
      do jj=1,k+1
        h(:,jj,i) = hcnew(:,jj,i)
      enddo
    enddo

    call matvecmult(gondag,k+1,m+1,srv,m+1,c)

    do i=k+2,m+1
      c(:,i) = 0.0_KR2
    enddo

    do jj=1,k+1
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vt(:) = 0.0_KR2
              do kk=1,m+1
                do icri = 1,5,2
                  vt(icri)   = vt(icri) &
                      + v(icri  ,i,id,ieo,ibleo,kk)*gon(1,kk,jj) &
                      - v(icri+1,i,id,ieo,ibleo,kk)*gon(2,kk,jj) 
                  vt(icri+1) = vt(icri+1) &
                      + v(icri  ,i,id,ieo,ibleo,kk)*gon(2,kk,jj) &
                      + v(icri+1,i,id,ieo,ibleo,kk)*gon(1,kk,jj) 
                enddo
              enddo
              work(:,i,id,ieo,ibleo,jj) = vt(:)
            enddo
          enddo
        enddo
      enddo
    enddo

    do i=1,k+1
      v(:,:,:,:,:,i) = work(:,:,:,:,:,i)
    enddo

    do i=1,k
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,k+1),tmp1,MRT2)

      do icri = 1,5,2 ! 6=nri*nc
        !do jj = 1,ntotal
        do jj = 1,nvhalf
          v(icri  ,jj,:,:,:,k+1) = v(icri  ,jj,:,:,:,k+1) &
                             - tmp1(1)*v(icri  ,jj,:,:,:,i) &
                             + tmp1(2)*v(icri+1,jj,:,:,:,i)
          v(icri+1,jj,:,:,:,k+1) = v(icri+1,jj,:,:,:,k+1) &
                             - tmp1(2)*v(icri  ,jj,:,:,:,i) &
                             - tmp1(1)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo

    call vecdot(v(:,:,:,:,:,k+1),v(:,:,:,:,:,k+1),tmp1,MRT2)
    tmp1(1) = sqrt(tmp1(1))

    !do jj = 1,ntotal
    do jj = 1,nvhalf
      v(:,jj,:,:,:,k+1) = (1.0_KR2 / tmp1(1))*v(:,jj,:,:,:,k+1)
    enddo

!if (myid==0) then
!  write(*,*) 'cycle:',icycle,'resnorm:',rn(1)/rninit(1)
!endif

      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", form="formatted",status="old",position="append")
!      BS  write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
       !BS write(unit=8,fmt="(a9,i4,a3,a5,es10.6,es10.6,es17.10)") "gmresdr",icycle,"rn","rnit",rn(1),rninit(1),rn(1)/rninit(1)
                write(unit=8,fmt="(a9,i4,a4,es11.4,a6,es11.4,es17.10)") "gmrsEIG",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
                    write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)")"gmrsEIG",icycle,"xrn=",xrn(1),"xrnit=",rninit(1),xrn(1)/rninit(1)

     !  write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif



do vena=1,k


!BS added this part for true residual



      call Hdbletm(wv,u,GeeGooinv,evector(:,:,:,:,:,vena),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              do icri =1,5,2 


    mina(icri,i,id,ieo,ibleo) = evalue(1,vena)*evector(icri,i,id,ieo,ibleo,vena)&
                                   - evalue(2,vena)*evector(icri+1,i,id,ieo,ibleo,vena)
    mina(icri+1,i,id,ieo,ibleo) = evalue(1,vena)*evector(icri+1,i,id,ieo,ibleo,vena)&
                                   + evalue(2,vena)*evector(icri,i,id,ieo,ibleo,vena)

              enddo
            enddo
         enddo
     enddo
enddo



      xtmp = wv - mina


      call vecdot(xtmp,xtmp,xrn,MRT2)
            rina(1)=sqrt(xrn(1))
            rina(2)=0.0_KR2

!BS ************print of trueresidual.dat starts*************
if (.false.) then 


        if (myid==0) then
            inquire(file=trim(rwdir(myid+1))//"trueresidual.dat", exist=exists)
            if (.not. exists) then

               print *, "File does not exist. Creating it."
               open(unit=46,file=trim(rwdir(myid+1))//"trueresidual.dat",status="new",&
               action="write",form="formatted")
               close(unit=46,status="keep")
            endif
       endif


             if (myid==0) then
            open(unit=46,file=trim(rwdir(myid+1))//"trueresidual.dat",status="old",action="write",&
            form="formatted",position="append")
                 write(unit=46,fmt="(i5,a5,i5,a5,es20.12)") icycle,"   ",vena," ",rina(1)
           close(unit=46,status="keep")
          endif

endif!true or false
!BS **********print of trueresidual.dat ends***********
enddo






    j = k+1
    icycle = icycle + 1
  enddo
        
        !  if (myid==0) then !BS
              
              !BS   write(*,fmt="(a12,f19.11)") "evalueEIG:",evalue
          !    print *,'evalueEIG:',evalue
         ! endif
     if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
   !BS  write(unit=8,fmt="(a12,i9,es17.10)")
   !"--gmresdr",icycle-1,rn(1)/rninit(1)
                write(unit=8,fmt="(a9,i4,a4,es11.4,a6,es11.4,es17.10)") "gmrsEIG",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)

    close(unit=8,status="keep")
  endif
  


  if (myid==0) then
    do i =1,k 
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
   !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr",icycle-1,rn(1)/rninit(1)
                write(unit=8,fmt="(i9,F17.12,F17.12)") i,evalue(1,i),evalue(2,i)

    close(unit=8,status="keep")
    enddo  
  endif
end subroutine gmresdrEIG






















!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine gmresdr(rwdir,phi,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=phi for the vector x.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresdr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo, mvp
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w, work
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,ztmp, q
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht,matss
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3, hprint, t, hhh
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws, ev
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h , xt ,h1
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(2,kmaxGMRES,kmaxGMRES)      :: greal
   
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp, bopart
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
    real(kind=KR2), dimension(6,ntotal,4,2,8,1)  :: try1, try2, try3, try4,emp 
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
    integer(kind=KI) :: didmaindo,roww
    integer(kind=KI), dimension(nmaxGMRES)              :: sortev

    real(kind=KR), dimension(nxyzt,3,4)   :: realz2noise, imagz2noise


! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.



! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)
    didmaindo = 0
    icycle = 1
    idag = 0
    ss = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR
    mvp = 0
 ! htemp = 0.0_KR
 ! site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf
 !      do idirac=1,4
 !        do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !              !print *, "site,icolorr =", site,icolorr
 !              !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps
!
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
!
!                call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!
!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
! To print single rhs source vector use ..
!
!               !irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!               !       print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! idirac
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock

!*****MORGAN'S STEP 1: Start.
! Compute r=phi-M*x and v=r/|r| and beta=|r|.
  call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    mvp = mvp + 1
!     try1 = 0.0_KR2
!     try2 = 0.0_KR2
!     try3 = 0.0_KR2
!     try4 = 0.0_KR2
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
     v(:,i,:,:,:,1) = r(:,i,:,:,:)
     xt(:,i,:,:,:) = x(:,i,:,:,:)
    enddo ! i
!    call Hdbletm(try2(:,:,:,:,:,1),u,GeeGooinv,try1(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(try3(:,:,:,:,:,1),u,GeeGooinv,try2(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(try4(:,:,:,:,:,1),u,GeeGooinv,try3(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(v(:,:,:,:,:,1),u,GeeGooinv,try4(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)


!    do i = 1,nvhalf
!     v(:,i,:,:,:,1) = h1(:,i,:,:,:)
!    enddo!!
    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))
    normnum = beta(1)

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.
    itercount = 0
    j = 0
    maindo: do




     if ( icycle > kcyclim) exit maindo

     j = j + 1
     jp1 = j + 1
     itercount = itercount + 1

!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.

!     try1 = 0.0_KR2
!     try2 = 0.0_KR2
!     try3 = 0.0_KR2
!     try4 = 0.0_KR
     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                  iblv,MRT)

     mvp = mvp + 1
!     call Hdbletm(try2(:,:,:,:,:,1),u,GeeGooinv,try1(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(try3(:,:,:,:,:,1),u,GeeGooinv,try2(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(try4(:,:,:,:,:,1),u,GeeGooinv,try3(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,try4(:,:,:,:,:,1),idag, &
!                  coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc, &
!                  ib,lbd,iblv,MRT)

!     do i=1,nvhalf
!      v(:,i,:,:,:,jp1) = h1(:,i,:,:,:)
!     enddo!i
     do i = 1,j
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)
      hc(:,i,j) = beta(:)

!     print *, "i,j, hc(:,i,j)=", i,j, hc(:,i,j)

      hc2(:,i,j) = hc(:,i,j)
      hc3(:,i,j) = hc(:,i,j)
      hcht(1,j,i) = hc(1,i,j)
      hcht(2,j,i) = -hc(2,i,j)
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                              - beta(1)*v(icri  ,k,:,:,:,i) &
                              + beta(2)*v(icri+1,k,:,:,:,i)
        v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                              - beta(2)*v(icri  ,k,:,:,:,i) &
                              - beta(1)*v(icri+1,k,:,:,:,i)
       enddo ! k
      enddo ! icri
     enddo ! i



     call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
     hc(1,jp1,j) = sqrt(beta(1))
     hc(2,jp1,j) = 0.0_KR2
     hc2(:,jp1,j) = hc(:,jp1,j)
     hc3(:,jp1,j) = hc(:,jp1,j)
     hcht(1,j,jp1) = hc(1,jp1,j)
     hcht(2,j,jp1) = -hc(2,jp1,j)
     const = 1.0_KR2/sqrt(beta(1))
     v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
     c(:,jp1) = 0.0_KR2
     c2(:,jp1) = c(:,jp1)


! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.
     if (icycle/=1) then
      do jj = 1,kDR
       do i = jj+1,kDR+1
        tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
               + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
        tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
               + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
        tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
               - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
        tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
               - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
        hc(:,jj,j) = tv1(:)
        hc(:,i,j) = tv2(:)
       enddo ! i
      enddo ! jj
      if (j>kDR+1) then
       do i = kDR+1,j-1
        tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
               + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
        tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
               + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
        tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
               - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
        tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
               - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
        hc(:,i,j) = tv1(:)
        hc(:,i+1,j) = tv2(:)
       enddo ! i
      endif
     elseif (j/=1) then
      do i = 1,j-1
       tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
              + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
       tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
              + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
       tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
              - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
       tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
              - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
       hc(:,i,j) = tv1(:)
       hc(:,i+1,j) = tv2(:)
      enddo ! i
     endif
     amags = hc(1,j,j)**2 + hc(2,j,j)**2
     tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
     tv(2) = 0.0_KR2
     gc(1,j) = sqrt(amags)/tv(1)
     gc(2,j) = 0.0_KR2
     con2 = gc(1,j)/amags
     gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
     gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
     hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
     hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
     hc(:,j+1,j) = 0.0_KR2
     tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
     tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
     tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
     tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
     c(:,j) = tv1(:)
     c(:,j+1) = tv2(:)
     do i = 1,j
      ss(:,i) = c(:,i)
     enddo ! i
     con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
     const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
     ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
     ss(1,j) = const
     if (j/=1) then
      do i = 1,j-1
       ir = j - i + 1
       irm1 = ir - 1
       do jj = 1,irm1
        const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
        ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
        ss(1,jj) = const
       enddo ! jj
       con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
       const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
       ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
       ss(1,irm1) = const
      enddo ! i
     endif

! Form the approximate new solution x = xt + V*ss.
     xb = 0.0_KR2
     do jj = 1,j
      do icri = 1,5,2 ! 6=nri*nc
       do i = 1,nvhalf
        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) + ss(1,jj)*v(icri  ,i,:,:,:,jj)&
                                                - ss(2,jj)*v(icri+1,i,:,:,:,jj)
        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) + ss(2,jj)*v(icri  ,i,:,:,:,jj)&
                                                + ss(1,jj)*v(icri+1,i,:,:,:,jj)
       enddo ! i
      enddo ! icri
     enddo ! jj
     do i = 1,nvhalf
      x(:,i,:,:,:) = xt(:,i,:,:,:) + xb(:,i,:,:,:)
     enddo ! i

! Define a small residual vector, srv = c-Hbar*ss, which corresponds to the
! kDR+1 column of the new V that will be formed.
     do i = 1,nDR+1
      srv(:,i) = c2(:,i)
     enddo ! i
     do jj = 1,nDR
      do i = 1,nDR+1
       srv(1,i) = srv(1,i) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
       srv(2,i) = srv(2,i) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
      enddo ! i
     enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
     if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=phi-M*x.

      call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      mvp = mvp + 1
      do i = 1,nvhalf
       r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
      enddo ! i

      beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
!        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr-norm",mvp,beta(1)/normnum
       close(unit=8,status="keep")
      endif

      if ((beta(1)/normnum)<resmax .and. itercount>=itermin) then
         if(didmaindo==0) then
         print *, "Everyone needs more resmax!"
         endif ! didmaindo
         exit maindo
      endif ! beta
      didmaindo=1
      do i = 1,nvhalf
       xt(:,i,:,:,:) = x(:,i,:,:,:)
      enddo ! i

!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)
      do i = 1,nDR
       em(:,i,1) = 0.0_KR
      enddo ! i
      em(1,nDR,1) = 1.0_KR2
      nrhs = 1





      call linearsolver(nDR,nrhs,hcht,ipiv,em)


      do i = 1,nDR
       hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                    - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                    - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
       hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                    + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                    - em(2,i,1)*hc2(2,nDR+1,nDR)**2
      enddo ! i

      ilo = 1
      ihi = nDR
      call hessenberg(hc2,nDR,ilo,ihi,tau,work)




      do jj = 1,nDR
       do i = jj,nDR
        z(:,i,jj) = hc2(:,i,jj)
       enddo ! i
      enddo ! jj
      call qgenerator(nDR,ilo,ihi,z,tau,work)
      ischur = 1
      call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

         if (myid==0) then
            print *, "------------"
            print *, "evaluesgmresdr:", w
            print *, "------------"
         endif 

!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan -- 
! so reorder the harmonic Ritz values and rearrange the Schur form.
      mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
      do i = 2,nDR
       mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
       is = 0
       ritzloop: do
        is = is + 1
        if (is>i-1) exit ritzloop
        if (mag(i)<mag(is)) then
         tval = mag(i)
         do ivb = i-1,is,-1
          mag(ivb+1) = mag(ivb)
         enddo ! ivb
         mag(is) = tval
         exit ritzloop
        endif
       enddo ritzloop
      enddo ! i
      do i = 1,nDR
       myselect(i) = .false.
       if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
      enddo ! i

      call orgschur(myselect,nDR,hc2,z,w,idis)

!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.
      do i = 1,kDR
       z(:,nDR+1,i) = 0.0_KR2
      enddo ! i
      do i = 1,nDR+1
       z(:,i,kDR+1) = srv(:,i)
      enddo ! i
      do jj = 1,kDR
       tv = 0.0_KR2
       do i = 1,nDR+1
        tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
        tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
       enddo ! i
       do i = 1,nDR+1
        z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
        z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
       enddo ! i
      enddo ! jj
      jj = nDR + 1
      rv = twonorm(z(:,:,kDR+1),jj)
      rv = 1.0_KR2/rv
      do i = 1,nDR+1
       z(:,i,kDR+1) = rv*z(:,i,kDR+1)
      enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.
      do jj = 1,kDR
       do ii = 1,nDR+1
        ws(:,ii,jj) = 0.0_KR2
       enddo ! ii
       do ii = 1,nDR
        do i = 1,nDR+1
         ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                 - z(2,ii,jj)*hc3(2,i,ii)
         ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                 + z(2,ii,jj)*hc3(1,i,ii)
        enddo ! i
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hcnew(:,ii,jj) = 0.0_KR2
        do i = 1,nDR+1
         hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                         + z(2,i,ii)*ws(2,i,jj)
         hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                         - z(2,i,ii)*ws(1,i,jj)
        enddo ! i
       enddo ! ii
      enddo ! jj

      do jj = 1,nDR
       do ii = 1,nDR
        hcht(:,ii,jj) = 0.0_KR2
        hc2(:,ii,jj) = 0.0_KR2
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hc(:,ii,jj) = hcnew(:,ii,jj)
        hc2(:,ii,jj) = hcnew(:,ii,jj)
        hc3(:,ii,jj) = hcnew(:,ii,jj)
       enddo ! ii
       do ii = 1,kDR+1
        hcht(1,jj,ii) = hcnew(1,ii,jj)
        hcht(2,jj,ii) = -hcnew(2,ii,jj)
       enddo ! ii
      enddo ! jj
      do ii = 1,kDR+1
       c(:,ii) = 0.0_KR2
       do i = 1,nDR+1
        c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i) + z(2,i,ii)*srv(2,i)
        c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i) - z(2,i,ii)*srv(1,i)
       enddo ! i
       c2(:,ii) = c(:,ii)
      enddo ! ii

      vtemp = 0.0_KR

      do ibleo = 1,8
       do ieo = 1,2
        do id = 1,4
         do i = 1,nvhalf
          do jj = 1,kDR+1
           vt(:,jj) = 0.0_KR2
           do k = 1,nDR+1
            do icri = 1,5,2 ! 6=nri*nc
             vt(icri  ,jj) = vt(icri  ,jj) &
                           + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
             vt(icri+1,jj) = vt(icri+1,jj) &
                           + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
            enddo ! icri
           enddo ! k
          enddo ! jj
          do jj = 1,kDR+1
           v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
           vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
          enddo ! jj
         enddo ! i
        enddo ! id
       enddo ! ieo
      enddo ! ibleo
         !if(myid==0) then
         ! print *, "vtemp in this shizat-only 1:kDR", vtemp
         !endif ! myid

!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.
      do jj = 1,kDR
       call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
       do icri = 1,5,2 ! 6=nri*nc
        do i = 1,nvhalf
         v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                 - beta(1)*v(icri  ,i,:,:,:,jj) &
                                 + beta(2)*v(icri+1,i,:,:,:,jj)
         v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                 - beta(2)*v(icri  ,i,:,:,:,jj) &
                                 - beta(1)*v(icri+1,i,:,:,:,jj)
        enddo ! i
       enddo ! icri
      enddo ! jj
      call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)
      const = 1.0_KR2/sqrt(beta(1))
      do i = 1,nvhalf
       v(:,i,:,:,:,kDR+1)     = const*v(:,i,:,:,:,kDR+1)
      enddo ! i

! Need to have the vtemp vector for the gmresproj routine....

      do jj = 1,nvhalf
       vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1)
      enddo ! jj

! Rotations for newly formed hc() matrix.
      do jj = 1,kDR
       do i = jj+1,kDR+1
        amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
        con2 = 1.0_KR2/amags
        tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
        tv(2) = 0.0_KR2
        gca(1,i,jj) = sqrt(amags)/tv(1)
        gca(2,i,jj) = 0.0_KR2
        gsa(1,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
        gsa(2,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
        do j = jj,kDR
         tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                         + gsa(2,i,jj)*hc(2,i,j)
         tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                         - gsa(2,i,jj)*hc(1,i,j)
         tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                        + gsa(2,i,jj)*hc(2,jj,j)
         tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                        - gsa(2,i,jj)*hc(1,jj,j)
         hc(:,jj,j) = tv1(:)
         hc(:,i,j) = tv2(:)
        enddo ! j
        tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
        tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
        tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
        tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
        c(:,jj) = tv1(:)
        c(:,i) = tv2(:)
       enddo ! i
      enddo ! jj
      j = kDR
      icycle = icycle + 1
     endif
    enddo maindo

     !     if (myid==0) then
     !        print *, "evaluesgmresdr:-Dean wuz here"
     !        do ii =1,kDR
     !          print *, w(:,ii)
     !        enddo ! ii
     !     endif 
   
   call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
    enddo ! i
        beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))

!       print *, "final residue=",beta(1)/normnum 

end subroutine gmresdr
!-----------------------------------------------------------------------------
subroutine mmgmresdr(rwdir,phi,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=phi for the vector x.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresdr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w, work
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,ztmp, q
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht,matss
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3, hprint, t, hhh
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws, ev
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h , xt ,h1,&
                                                               inter,inter1,tt
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(2,kmaxGMRES,kmaxGMRES)      :: greal
   
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp, bopart
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
    real(kind=KR2), dimension(6,ntotal,4,2,8,1)  :: try1, try2, try3, try4,emp 
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
    integer(kind=KI) :: didmaindo,roww
    integer(kind=KI), dimension(nmaxGMRES)              :: sortev

    real(kind=KR), dimension(nxyzt,3,4)   :: realz2noise, imagz2noise


! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.



! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)
    didmaindo = 0
    icycle = 1
    idag = 0
    ss = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR

 ! htemp = 0.0_KR
 ! site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf
 !      do idirac=1,4
 !        do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !              !print *, "site,icolorr =", site,icolorr
 !              !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps
!
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
!
!                call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!
!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
! To print single rhs source vector use ..
!
!               !irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!               !       print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! idirac
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock

!*****MORGAN'S STEP 1: Start.

    idag = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute r=phi-M*x and v=r/|r| and beta=|r|.
 !!!!!test the M^*M case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    idag = 1
!  do i=1,nvhalf
!   tt(:,i,:,:,:)=phi(:,i,:,:,:)
!  enddo!i
!  call Hdbletm(inter1,u,GeeGooinv,tt,idag,coact,kappa,iflag,bc, &
!                vecbl,vecblinv, &
!                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   idag = 0 
  call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)!check M^dagger*M

!   idag = 1
!  call Hdbletm(h,u,GeeGooinv,inter,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
!                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)!check M^dagger*M
!   idag = 0
!     try1 = 0.0_KR2
!     try2 = 0.0_KR2
!     try3 = 0.0_KR2
!     try4 = 0.0_KR2
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
     v(:,i,:,:,:,1) = r(:,i,:,:,:)
     xt(:,i,:,:,:) = x(:,i,:,:,:)
    enddo ! i
!    call Hdbletm(try2(:,:,:,:,:,1),u,GeeGooinv,try1(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(try3(:,:,:,:,:,1),u,GeeGooinv,try2(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(try4(:,:,:,:,:,1),u,GeeGooinv,try3(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)
!    call Hdbletm(v(:,:,:,:,:,1),u,GeeGooinv,try4(:,:,:,:,:,1),idag,coact, &
!                 kappa,iflag,bc,vecbl,vecblinv, myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                 iblv,MRT)


!    do i = 1,nvhalf
!     v(:,i,:,:,:,1) = h1(:,i,:,:,:)
!    enddo!!
    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))
    normnum = beta(1)

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.
    itercount = 0
    j = 0
    maindo: do




     if ( icycle > kcyclim) exit maindo

     j = j + 1
     jp1 = j + 1
     itercount = itercount + 1

!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.

!     try1 = 0.0_KR2
!     try2 = 0.0_KR2
!     try3 = 0.0_KR2
!     try4 = 0.0_KR2
     try1=0
     idag = 0
     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                  iblv,MRT)
!     idag = 1
!     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,try1(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)

!     call Hdbletm(try2(:,:,:,:,:,1),u,GeeGooinv,try1(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(try3(:,:,:,:,:,1),u,GeeGooinv,try2(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(try4(:,:,:,:,:,1),u,GeeGooinv,try3(:,:,:,:,:,1),idag,coact, &
!                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                  iblv,MRT)
!     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,try4(:,:,:,:,:,1),idag, &
!                  coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc, &
!                  ib,lbd,iblv,MRT)

!     do i=1,nvhalf
!      v(:,i,:,:,:,jp1) = h1(:,i,:,:,:)
!     enddo!i
     do i = 1,j
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)
      hc(:,i,j) = beta(:)

!     print *, "i,j, hc(:,i,j)=", i,j, hc(:,i,j)

      hc2(:,i,j) = hc(:,i,j)
      hc3(:,i,j) = hc(:,i,j)
      hcht(1,j,i) = hc(1,i,j)
      hcht(2,j,i) = -hc(2,i,j)
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                              - beta(1)*v(icri  ,k,:,:,:,i) &
                              + beta(2)*v(icri+1,k,:,:,:,i)
        v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                              - beta(2)*v(icri  ,k,:,:,:,i) &
                              - beta(1)*v(icri+1,k,:,:,:,i)
       enddo ! k
      enddo ! icri
     enddo ! i



     call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
     hc(1,jp1,j) = sqrt(beta(1))
     hc(2,jp1,j) = 0.0_KR2
     hc2(:,jp1,j) = hc(:,jp1,j)
     hc3(:,jp1,j) = hc(:,jp1,j)
     hcht(1,j,jp1) = hc(1,jp1,j)
     hcht(2,j,jp1) = -hc(2,jp1,j)
     const = 1.0_KR2/sqrt(beta(1))
     v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
     c(:,jp1) = 0.0_KR2
     c2(:,jp1) = c(:,jp1)


! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.
     if (icycle/=1) then
      do jj = 1,kDR
       do i = jj+1,kDR+1
        tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
               + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
        tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
               + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
        tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
               - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
        tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
               - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
        hc(:,jj,j) = tv1(:)
        hc(:,i,j) = tv2(:)
       enddo ! i
      enddo ! jj
      if (j>kDR+1) then
       do i = kDR+1,j-1
        tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
               + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
        tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
               + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
        tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
               - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
        tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
               - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
        hc(:,i,j) = tv1(:)
        hc(:,i+1,j) = tv2(:)
       enddo ! i
      endif
     elseif (j/=1) then
      do i = 1,j-1
       tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
              + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
       tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
              + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
       tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
              - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
       tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
              - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
       hc(:,i,j) = tv1(:)
       hc(:,i+1,j) = tv2(:)
      enddo ! i
     endif
     amags = hc(1,j,j)**2 + hc(2,j,j)**2
     tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
     tv(2) = 0.0_KR2
     gc(1,j) = sqrt(amags)/tv(1)
     gc(2,j) = 0.0_KR2
     con2 = gc(1,j)/amags
     gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
     gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
     hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
     hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
     hc(:,j+1,j) = 0.0_KR2
     tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
     tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
     tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
     tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
     c(:,j) = tv1(:)
     c(:,j+1) = tv2(:)
     do i = 1,j
      ss(:,i) = c(:,i)
     enddo ! i
     con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
     const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
     ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
     ss(1,j) = const
     if (j/=1) then
      do i = 1,j-1
       ir = j - i + 1
       irm1 = ir - 1
       do jj = 1,irm1
        const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
        ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
        ss(1,jj) = const
       enddo ! jj
       con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
       const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
       ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
       ss(1,irm1) = const
      enddo ! i
     endif

! Form the approximate new solution x = xt + V*ss.
     xb = 0.0_KR2
     do jj = 1,j
      do icri = 1,5,2 ! 6=nri*nc
       do i = 1,nvhalf
        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) + ss(1,jj)*v(icri  ,i,:,:,:,jj)&
                                                - ss(2,jj)*v(icri+1,i,:,:,:,jj)
        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) + ss(2,jj)*v(icri  ,i,:,:,:,jj)&
                                                + ss(1,jj)*v(icri+1,i,:,:,:,jj)
       enddo ! i
      enddo ! icri
     enddo ! jj
     do i = 1,nvhalf
      x(:,i,:,:,:) = xt(:,i,:,:,:) + xb(:,i,:,:,:)
     enddo ! i

! Define a small residual vector, srv = c-Hbar*ss, which corresponds to the
! kDR+1 column of the new V that will be formed.
     do i = 1,nDR+1
      srv(:,i) = c2(:,i)
     enddo ! i
     do jj = 1,nDR
      do i = 1,nDR+1
       srv(1,i) = srv(1,i) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
       srv(2,i) = srv(2,i) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
      enddo ! i
     enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
     if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=phi-M*x.

      idag=0
       inter = 0
      call Hdbletm(h,u,GeeGooinv,xb,idag,coact,kappa,iflag, &
                   bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!      idag=1
!      call Hdbletm(h,u,GeeGooinv,inter,idag,coact,kappa,iflag, &
!                   bc,vecbl,vecblinv, &
!                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      do i = 1,nvhalf
       r(:,i,:,:,:) = r(:,i,:,:,:) - h(:,i,:,:,:)
      enddo ! i

     idag=0

      beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
!        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr-norm",itercount,beta(1)/normnum
       close(unit=8,status="keep")
      endif

      if ((beta(1)/normnum)<resmax .and. itercount>=itermin) then
         if(didmaindo==0) then
         print *, "Everyone needs more resmax!"
         endif ! didmaindo
         exit maindo
      endif ! beta
      didmaindo=1
      do i = 1,nvhalf
       xt(:,i,:,:,:) = x(:,i,:,:,:)
      enddo ! i

!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)
      do i = 1,nDR
       em(:,i,1) = 0.0_KR
      enddo ! i
      em(1,nDR,1) = 1.0_KR2
      nrhs = 1





      call linearsolver(nDR,nrhs,hcht,ipiv,em)


      do i = 1,nDR
       hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                    - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                    - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
       hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                    + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                    - em(2,i,1)*hc2(2,nDR+1,nDR)**2
      enddo ! i

      ilo = 1
      ihi = nDR
      call hessenberg(hc2,nDR,ilo,ihi,tau,work)




      do jj = 1,nDR
       do i = jj,nDR
        z(:,i,jj) = hc2(:,i,jj)
       enddo ! i
      enddo ! jj
      call qgenerator(nDR,ilo,ihi,z,tau,work)
      ischur = 1
      call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

         if (myid==0) then
            print *, "------------"
            print *, "evaluesgmresdr:", w
            print *, "------------"
         endif 

!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan -- 
! so reorder the harmonic Ritz values and rearrange the Schur form.
      mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
      do i = 2,nDR
       mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
       is = 0
       ritzloop: do
        is = is + 1
        if (is>i-1) exit ritzloop
        if (mag(i)<mag(is)) then
         tval = mag(i)
         do ivb = i-1,is,-1
          mag(ivb+1) = mag(ivb)
         enddo ! ivb
         mag(is) = tval
         exit ritzloop
        endif
       enddo ritzloop
      enddo ! i
      do i = 1,nDR
       myselect(i) = .false.
       if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
      enddo ! i

      call orgschur(myselect,nDR,hc2,z,w,idis)

!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.
      do i = 1,kDR
       z(:,nDR+1,i) = 0.0_KR2
      enddo ! i
      do i = 1,nDR+1
       z(:,i,kDR+1) = srv(:,i)
      enddo ! i
      do jj = 1,kDR
       tv = 0.0_KR2
       do i = 1,nDR+1
        tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
        tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
       enddo ! i
       do i = 1,nDR+1
        z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
        z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
       enddo ! i
      enddo ! jj
      jj = nDR + 1
      rv = twonorm(z(:,:,kDR+1),jj)
      rv = 1.0_KR2/rv
      do i = 1,nDR+1
       z(:,i,kDR+1) = rv*z(:,i,kDR+1)
      enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.
      do jj = 1,kDR
       do ii = 1,nDR+1
        ws(:,ii,jj) = 0.0_KR2
       enddo ! ii
       do ii = 1,nDR
        do i = 1,nDR+1
         ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                 - z(2,ii,jj)*hc3(2,i,ii)
         ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                 + z(2,ii,jj)*hc3(1,i,ii)
        enddo ! i
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hcnew(:,ii,jj) = 0.0_KR2
        do i = 1,nDR+1
         hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                         + z(2,i,ii)*ws(2,i,jj)
         hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                         - z(2,i,ii)*ws(1,i,jj)
        enddo ! i
       enddo ! ii
      enddo ! jj

      do jj = 1,nDR
       do ii = 1,nDR
        hcht(:,ii,jj) = 0.0_KR2
        hc2(:,ii,jj) = 0.0_KR2
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hc(:,ii,jj) = hcnew(:,ii,jj)
        hc2(:,ii,jj) = hcnew(:,ii,jj)
        hc3(:,ii,jj) = hcnew(:,ii,jj)
       enddo ! ii
       do ii = 1,kDR+1
        hcht(1,jj,ii) = hcnew(1,ii,jj)
        hcht(2,jj,ii) = -hcnew(2,ii,jj)
       enddo ! ii
      enddo ! jj
      do ii = 1,kDR+1
       c(:,ii) = 0.0_KR2
       do i = 1,nDR+1
        c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i) + z(2,i,ii)*srv(2,i)
        c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i) - z(2,i,ii)*srv(1,i)
       enddo ! i
       c2(:,ii) = c(:,ii)
      enddo ! ii

      vtemp = 0.0_KR

      do ibleo = 1,8
       do ieo = 1,2
        do id = 1,4
         do i = 1,nvhalf
          do jj = 1,kDR+1
           vt(:,jj) = 0.0_KR2
           do k = 1,nDR+1
            do icri = 1,5,2 ! 6=nri*nc
             vt(icri  ,jj) = vt(icri  ,jj) &
                           + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
             vt(icri+1,jj) = vt(icri+1,jj) &
                           + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
            enddo ! icri
           enddo ! k
          enddo ! jj
          do jj = 1,kDR+1
           v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
           vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
          enddo ! jj
         enddo ! i
        enddo ! id
       enddo ! ieo
      enddo ! ibleo
         !if(myid==0) then
         ! print *, "vtemp in this shizat-only 1:kDR", vtemp
         !endif ! myid

!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.
      do jj = 1,kDR
       call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
       do icri = 1,5,2 ! 6=nri*nc
        do i = 1,nvhalf
         v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                 - beta(1)*v(icri  ,i,:,:,:,jj) &
                                 + beta(2)*v(icri+1,i,:,:,:,jj)
         v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                 - beta(2)*v(icri  ,i,:,:,:,jj) &
                                 - beta(1)*v(icri+1,i,:,:,:,jj)
        enddo ! i
       enddo ! icri
      enddo ! jj
      call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)
      const = 1.0_KR2/sqrt(beta(1))
      do i = 1,nvhalf
       v(:,i,:,:,:,kDR+1)     = const*v(:,i,:,:,:,kDR+1)
      enddo ! i

! Need to have the vtemp vector for the gmresproj routine....

      do jj = 1,nvhalf
       vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1)
      enddo ! jj

! Rotations for newly formed hc() matrix.
      do jj = 1,kDR
       do i = jj+1,kDR+1
        amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
        con2 = 1.0_KR2/amags
        tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
        tv(2) = 0.0_KR2
        gca(1,i,jj) = sqrt(amags)/tv(1)
        gca(2,i,jj) = 0.0_KR2
        gsa(1,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
        gsa(2,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
        do j = jj,kDR
         tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                         + gsa(2,i,jj)*hc(2,i,j)
         tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                         - gsa(2,i,jj)*hc(1,i,j)
         tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                        + gsa(2,i,jj)*hc(2,jj,j)
         tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                        - gsa(2,i,jj)*hc(1,jj,j)
         hc(:,jj,j) = tv1(:)
         hc(:,i,j) = tv2(:)
        enddo ! j
        tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
        tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
        tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
        tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
        c(:,jj) = tv1(:)
        c(:,i) = tv2(:)
       enddo ! i
      enddo ! jj
      j = kDR
      icycle = icycle + 1
     endif
    enddo maindo

     !     if (myid==0) then
     !        print *, "evaluesgmresdr:-Dean wuz here"
     !        do ii =1,kDR
     !          print *, w(:,ii)
     !        enddo ! ii
     !     endif 
   
!   call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
!                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!    do i = 1,nvhalf
!     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
!    enddo ! i
!        beta =0.0_KR
!      call vecdot(r,r,beta,MRT2)
!      beta(1) = sqrt(beta(1))

!       print *, "final residue=",beta(1)/normnum 

end subroutine mmgmresdr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine psc(rwdir,b,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=b for the vector x.
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globlly-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresdrshift.
 
    use shift 

    character(len=*), intent(in),    dimension(:)             :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: b 
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)     :: x
    integer(kind=KI), intent(in),    dimension(:)             :: GMRES, bc, nms
    integer(kind=KI), intent(in)                              :: itermin, iflag, &
                                                                 myid, MRT, MRT2
                                                             
    real(kind=KR),    intent(in)                              :: resmax
    integer(kind=KI), intent(out)                             :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)             :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)         :: coact
    integer(kind=KI), intent(in),    dimension(:,:)           :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)           :: nn, iblv
    logical,          intent(in),    dimension(:)             :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)         :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)       :: ib
    logical,          intent(in),    dimension(:,:)           :: lbd
  
    integer(kind=KI), dimension(2)                            :: GMRESproj
    real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)       :: tempxvshift
    real(kind=KR2)                                            :: rn, rnc, rnt, rrr
    real(kind=KR2),   dimension(2)                            :: gam
    real(kind=KR2),   dimension(2,nshifts)                    :: gamtemp
    integer(kind=KI)                                          :: mvp, nDR, kDR, is, &
                                                                 icri, jj, idag, i, k, &
                                                                 ishift, ierr, j,l,m,n
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)       :: xshift

! NOTE~ the dimension of gdr must be big enough to hold all the norms
!       for the specified size of matrix

    real(kind=KR),    dimension(gdrsize,nshifts)              :: gdr 
!   real(kind=KR),    dimension(15000,nshifts)              :: gdr 
    real(kind=KR)                                             :: projresmax
    real(kind=KR),    dimension(2)                            :: beta
    real(kind=KR2),   dimension(nshifts)                      :: sigma

! This is still PSC

! Shift sigmamu to base mu above  (mtmqcd(1,2))
 !  print *,'calling psc from inverters.f90'

! Inititilizations  
 
! NOTE ~ BIG NOTE - might need to pass in idag into psc if that part of gamma5mult is needed

    idag = 0
    mvp = 0
    gdr = 0.0_KR
      v = 0.0_KR
    

    nDR = GMRES(1)
    kDR = GMRES(2)

    if (isignal == 1) then

      xshift(:,:,:,:,:,1) = x(:,:,:,:,:)
      vtemp = 0.0_KR2
      v = 0.0_KR


      call gmresdr(rwdir,b,xshift(:,:,:,:,:,1),GMRES,resmax,itermin, &
                     itercount,u,GeeGooinv,iflag,kappa,coact,bc,vecbl, &
                     vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)


      beg(:,1:nvhalf,:,:,:) = vtemp(:,1:nvhalf,:,:,:,kDR+1)
     
    else ! (isignal == 1)

    itercount = 0
    projresmax = 1.0e-8_KR

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta, MRT2)


! For comparison with gmresdr(m,k) we need to compare with  
!     gmres(m-k) - proj(k)

    GMRESproj(1) = GMRES(1) - GMRES(2)
    GMRESproj(2) = GMRES(2)
!     GMRESproj(1) = 150
!     GMRESproj(2) = 100

    call gmresproject(rwdir,b,xshift,GMRESproj,projresmax,itermin, &
                      itercount,u,GeeGooinv,iflag,kappa,coact,bc,vecbl, &
                      vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv, &
                      MRT,MRT2,isignal,mvp,gdr)


    endif ! (isignal /= 1)

    x(:,:,:,:,:) = xshift(:,:,:,:,:,1)

end subroutine psc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine ppgmresdr(rwdir,phi,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) with polynomial preconditioning
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j, k, jp1, jj, p, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo, mvp
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv
!!!!!!!!!!!!!!!!!!!!!!!!!!solving the coeffcients for the polynomial!!!!!!!!!!!!
!    integer(kind=KI), dimension(2)                          :: ipiv2!order of po                                                                     lynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w, work
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,ztmp, q
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht,matss
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3, hprint, t, hhh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    real(kind=KR2),   dimension(2,6,6)                      :: lsmat
!    real(kind=KR2),   dimension(2,6,1)                      :: cls
!    real(kind=KR2),   dimension(2,6)                        :: co!coefficients
!!!!!!!!!!!parameters in determing the preconditioning polynomial!!!!!!!!!!!!!!!

    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws, ev
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h , xt
!!!!!!!!!!!Some intermediate parameters used for the preconditioning!!!!!!!!!!!!
!    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: w1, z1 ,z2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
!    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: v, vprime
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(2,kmaxGMRES,kmaxGMRES)      :: greal
   
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp, bopart
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
!    real(kind=KR2), dimension(6,ntotal,4,2,8,6) :: try!6 is the degree of P(A)*A
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
    integer(kind=KI) :: didmaindo,roww
    integer(kind=KI), dimension(nmaxGMRES)              :: sortev

    real(kind=KR), dimension(nxyzt,3,4)   :: realz2noise, imagz2noise


! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.



! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)
    p = 6!the order of the polynomial
    didmaindo = 0
    icycle = 1
    idag = 0
    ss = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR
    w1 = 0.0_KR2
    z1 = 0.0_KR2
    z2 = 0.0_KR2
    try = 0.0_KR2
    mvp = 0 
! htemp = 0.0_KR
 ! site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf
 !      do idirac=1,4
 !        do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !              !print *, "site,icolorr =", site,icolorr
 !              !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps
!
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
!
!                call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!
!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
! To print single rhs source vector use ..
!
!               !irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!               !       print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! idirac
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock

!*****MORGAN'S STEP 1: Start.
! Compute r=phi-M*x and v=r/|r| and beta=|r|.
    call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    mvp = mvp + 1
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
     try(:,i,:,:,:,1) = r(:,i,:,:,:)
     xt(:,i,:,:,:) = x(:,i,:,:,:)
    enddo ! i
    !Determine the polynomial.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nvhalf
     vprime(:,k,:,:,:,1) = phi(:,k,:,:,:)
    enddo !k
    do i = 1,p
     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     mvp = mvp + 1
    enddo !i
    
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
!      print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
     enddo!j
    enddo!i
        



   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),phi(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)    
!    co = 0.0_KR2    
!    co(1,1) = 4
  ! if(myid==0) then
  !  do i=1,p
  !   print *, "i,result(:,i)=", i, co(:,i)
 !   enddo!i  
 !  endif!myid
!!!!
!!!!!Times the polynomial to the residue
        
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                         -co(2,1)*try(icri+1,k,:,:,:,1)
      y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                           +co(2,1)*try(icri,k,:,:,:,1)
     enddo!k
    enddo!icri

!    print *,"original component1:real=",w1(1,1,1,1,1)
!    print *,"original component2:imaginary=",w1(2,1,1,1,1)
 !   print *,"conditioned component1:real=",y(1,1,1,1,1)
!    print *,"conditioned component2:imaginary=",y(2,1,1,1,1)


   do i=1,p-1 
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
    mvp = mvp + 1
     do icri=1,5,2
      do k=1,nvhalf
       y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                             +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
       y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
   enddo!i
!     do k=1,nvhalf
!      r(:,k,:,:,:) = z2(:,k,:,:,:)
!     enddo!k
   
    
    do k=1,nvhalf
     v(:,k,:,:,:,1) = y(:,k,:,:,:,1) ! Use y=P(A)*r to generate V_(m+1)&Hbar_m
    enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))
    normnum = beta(1)

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.
    itercount = 0
    j = 0
    maindo: do




     if ( icycle > kcyclim) exit maindo

     j = j + 1
     jp1 = j + 1
     itercount = itercount + 1

!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
    try = 0.0_KR2
    y = 0.0_KR2
    do i=1,nvhalf
     try(:,i,:,:,:,1) = v(:,i,:,:,:,j)
    enddo!i
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                         -co(2,1)*try(icri+1,k,:,:,:,1)
      y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                           +co(2,1)*try(icri,k,:,:,:,1)
     enddo!k
    enddo!icri
    do i=1,p-1 
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     mvp = mvp + 1
     do icri=1,5,2
      do k=1,nvhalf
       y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                             +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
       y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i

    call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,y(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mvp = mvp + 1
      do i = 1,j
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)
      hc(:,i,j) = beta(:)

!     print *, "i,j, hc(:,i,j)=", i,j, hc(:,i,j)

      hc2(:,i,j) = hc(:,i,j)
      hc3(:,i,j) = hc(:,i,j)
      hcht(1,j,i) = hc(1,i,j)
      hcht(2,j,i) = -hc(2,i,j)
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                              - beta(1)*v(icri  ,k,:,:,:,i) &
                              + beta(2)*v(icri+1,k,:,:,:,i)
        v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                              - beta(2)*v(icri  ,k,:,:,:,i) &
                              - beta(1)*v(icri+1,k,:,:,:,i)
       enddo ! k
      enddo ! icri
     enddo ! i



     call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
     hc(1,jp1,j) = sqrt(beta(1))
     hc(2,jp1,j) = 0.0_KR2
     hc2(:,jp1,j) = hc(:,jp1,j)
     hc3(:,jp1,j) = hc(:,jp1,j)
     hcht(1,j,jp1) = hc(1,jp1,j)
     hcht(2,j,jp1) = -hc(2,jp1,j)
     const = 1.0_KR2/sqrt(beta(1))
     v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
     c(:,jp1) = 0.0_KR2
     c2(:,jp1) = c(:,jp1)


! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.
     if (icycle/=1) then
      do jj = 1,kDR
       do i = jj+1,kDR+1
        tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
               + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
        tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
               + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
        tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
               - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
        tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
               - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
        hc(:,jj,j) = tv1(:)
        hc(:,i,j) = tv2(:)
       enddo ! i
      enddo ! jj
      if (j>kDR+1) then
       do i = kDR+1,j-1
        tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
               + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
        tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
               + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
        tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
               - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
        tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
               - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
        hc(:,i,j) = tv1(:)
        hc(:,i+1,j) = tv2(:)
       enddo ! i
      endif
     elseif (j/=1) then
      do i = 1,j-1
       tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
              + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
       tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
              + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
       tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
              - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
       tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
              - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
       hc(:,i,j) = tv1(:)
       hc(:,i+1,j) = tv2(:)
      enddo ! i
     endif
     amags = hc(1,j,j)**2 + hc(2,j,j)**2
     tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
     tv(2) = 0.0_KR2
     gc(1,j) = sqrt(amags)/tv(1)
     gc(2,j) = 0.0_KR2
     con2 = gc(1,j)/amags
     gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
     gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
     hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
     hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
     hc(:,j+1,j) = 0.0_KR2
     tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
     tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
     tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
     tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
     c(:,j) = tv1(:)
     c(:,j+1) = tv2(:)
     do i = 1,j
      ss(:,i) = c(:,i)
     enddo ! i
     con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
     const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
     ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
     ss(1,j) = const
     if (j/=1) then
      do i = 1,j-1
       ir = j - i + 1
       irm1 = ir - 1
       do jj = 1,irm1
        const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
        ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
        ss(1,jj) = const
       enddo ! jj
       con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
       const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
       ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
       ss(1,irm1) = const
      enddo ! i
     endif

! Form the approximate new solution x = xt + V*ss.
     xb = 0.0_KR2
     do jj = 1,j
      do icri = 1,5,2 ! 6=nri*nc
       do i = 1,nvhalf
        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) + ss(1,jj)*v(icri  ,i,:,:,:,jj)&
                                                - ss(2,jj)*v(icri+1,i,:,:,:,jj)
        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) + ss(2,jj)*v(icri  ,i,:,:,:,jj)&
                                                + ss(1,jj)*v(icri+1,i,:,:,:,jj)
       enddo ! i
      enddo ! icri
     enddo ! jj
     do i = 1,nvhalf
      x(:,i,:,:,:) = xt(:,i,:,:,:) + xb(:,i,:,:,:)
     enddo ! i

! Define a small residual vector, srv = c-Hbar*ss, which corresponds to the
! kDR+1 column of the new V that will be formed.
     do i = 1,nDR+1
      srv(:,i) = c2(:,i)
     enddo ! i
     do jj = 1,nDR
      do i = 1,nDR+1
       srv(1,i) = srv(1,i) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
       srv(2,i) = srv(2,i) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
      enddo ! i
     enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
     if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=phi-M*x.

      call Hdbletm(h,u,GeeGooinv,xb,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      mvp = mvp + 1
      do i = 1,nvhalf
       r(:,i,:,:,:) = r(:,i,:,:,:) - h(:,i,:,:,:)
      enddo ! i
!!!!!!!!!!!!!!!!!!!!!! r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k=1,nvhalf
!     w1(:,k,:,:,:) = r(:,k,:,:,:)
!    enddo!k
    
!    do icri=1,5,2
!     do k=1,nvhalf
!      y(icri,k,:,:,:) = co(1,1)*w1(icri,k,:,:,:)-co(2,1)*w1(icri+1,k,:,:,:)
!      y(icri+1,k,:,:,:) = co(1,1)*w1(icri+1,k,:,:,:)+co(2,1)*w1(icri,k,:,:,:)
!     enddo!k
!    enddo!icri

!    do k=1,nvhalf
!     z1(:,k,:,:,:) = w1(:,k,:,:,:)
!    enddo!k

!    do i=2,p
!     call Hdbletm(z1(:,:,:,:,:),u,GeeGooinv,z1(:,:,:,:,:),idag,coact,kappa, &
!                  iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
!                  lbd,iblv,MRT )  !z1=M*z
!     do icri=1,5,2
!      do k=1,nvhalf
!       y(icri,k,:,:,:) = y(icri,k,:,:,:)+co(1,i)*z1(icri,k,:,:,:) &
!                         -co(2,i)*z1(icri+1,k,:,:,:)
!       y(icri+1,k,:,:,:) = y(icri+1,k,:,:,:)+co(1,i)*z1(icri+1,k,:,:,:) &
!                           +co(2,i)*z1(icri,k,:,:,:)   !y=P(A)*r
!      enddo!k
!     enddo!icri
!    enddo!i

!    do k=1,nvhalf
!     r(:,k,:,:,:)=y(:,k,:,:,:)
!    enddo!k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
!        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr-norm",mvp,beta(1)/normnum
       close(unit=8,status="keep")
      endif

      if ((beta(1)/normnum)<resmax .and. itercount>=itermin) then
         if(didmaindo==0) then
         print *, "Everyone needs more resmax!"
         endif ! didmaindo
         exit maindo
      endif ! beta
      didmaindo=1
      do i = 1,nvhalf
       xt(:,i,:,:,:) = x(:,i,:,:,:)
      enddo ! i

!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)
      do i = 1,nDR
       em(:,i,1) = 0.0_KR
      enddo ! i
      em(1,nDR,1) = 1.0_KR2
      nrhs = 1





      call linearsolver(nDR,nrhs,hcht,ipiv,em)


      do i = 1,nDR
       hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                    - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                    - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
       hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                    + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                    - em(2,i,1)*hc2(2,nDR+1,nDR)**2
      enddo ! i

      ilo = 1
      ihi = nDR
      call hessenberg(hc2,nDR,ilo,ihi,tau,work)




      do jj = 1,nDR
       do i = jj,nDR
        z(:,i,jj) = hc2(:,i,jj)
       enddo ! i
      enddo ! jj
      call qgenerator(nDR,ilo,ihi,z,tau,work)
      ischur = 1
      call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

         if (myid==0) then
            print *, "------------"
            print *, "evaluesgmresdr:", w
            print *, "------------"
         endif 

!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan -- 
! so reorder the harmonic Ritz values and rearrange the Schur form.
      mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
      do i = 2,nDR
       mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
       is = 0
       ritzloop: do
        is = is + 1
        if (is>i-1) exit ritzloop
        if (mag(i)<mag(is)) then
         tval = mag(i)
         do ivb = i-1,is,-1
          mag(ivb+1) = mag(ivb)
         enddo ! ivb
         mag(is) = tval
         exit ritzloop
        endif
       enddo ritzloop
      enddo ! i
      do i = 1,nDR
       myselect(i) = .false.
       if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
      enddo ! i

      call orgschur(myselect,nDR,hc2,z,w,idis)

!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.
      do i = 1,kDR
       z(:,nDR+1,i) = 0.0_KR2
      enddo ! i
      do i = 1,nDR+1
       z(:,i,kDR+1) = srv(:,i)
      enddo ! i
      do jj = 1,kDR
       tv = 0.0_KR2
       do i = 1,nDR+1
        tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
        tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
       enddo ! i
       do i = 1,nDR+1
        z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
        z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
       enddo ! i
      enddo ! jj
      jj = nDR + 1
      rv = twonorm(z(:,:,kDR+1),jj)
      rv = 1.0_KR2/rv
      do i = 1,nDR+1
       z(:,i,kDR+1) = rv*z(:,i,kDR+1)
      enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.
      do jj = 1,kDR
       do ii = 1,nDR+1
        ws(:,ii,jj) = 0.0_KR2
       enddo ! ii
       do ii = 1,nDR
        do i = 1,nDR+1
         ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                 - z(2,ii,jj)*hc3(2,i,ii)
         ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                 + z(2,ii,jj)*hc3(1,i,ii)
        enddo ! i
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hcnew(:,ii,jj) = 0.0_KR2
        do i = 1,nDR+1
         hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                         + z(2,i,ii)*ws(2,i,jj)
         hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                         - z(2,i,ii)*ws(1,i,jj)
        enddo ! i
       enddo ! ii
      enddo ! jj

      do jj = 1,nDR
       do ii = 1,nDR
        hcht(:,ii,jj) = 0.0_KR2
        hc2(:,ii,jj) = 0.0_KR2
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hc(:,ii,jj) = hcnew(:,ii,jj)
        hc2(:,ii,jj) = hcnew(:,ii,jj)
        hc3(:,ii,jj) = hcnew(:,ii,jj)
       enddo ! ii
       do ii = 1,kDR+1
        hcht(1,jj,ii) = hcnew(1,ii,jj)
        hcht(2,jj,ii) = -hcnew(2,ii,jj)
       enddo ! ii
      enddo ! jj
      do ii = 1,kDR+1
       c(:,ii) = 0.0_KR2
       do i = 1,nDR+1
        c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i) + z(2,i,ii)*srv(2,i)
        c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i) - z(2,i,ii)*srv(1,i)
       enddo ! i
       c2(:,ii) = c(:,ii)
      enddo ! ii

      vtemp = 0.0_KR

      do ibleo = 1,8
       do ieo = 1,2
        do id = 1,4
         do i = 1,nvhalf
          do jj = 1,kDR+1
           vt(:,jj) = 0.0_KR2
           do k = 1,nDR+1
            do icri = 1,5,2 ! 6=nri*nc
             vt(icri  ,jj) = vt(icri  ,jj) &
                           + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
             vt(icri+1,jj) = vt(icri+1,jj) &
                           + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
            enddo ! icri
           enddo ! k
          enddo ! jj
          do jj = 1,kDR+1
           v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
           vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
          enddo ! jj
         enddo ! i
        enddo ! id
       enddo ! ieo
      enddo ! ibleo
         !if(myid==0) then
         ! print *, "vtemp in this shizat-only 1:kDR", vtemp
         !endif ! myid

!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.
      do jj = 1,kDR
       call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
       do icri = 1,5,2 ! 6=nri*nc
        do i = 1,nvhalf
         v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                 - beta(1)*v(icri  ,i,:,:,:,jj) &
                                 + beta(2)*v(icri+1,i,:,:,:,jj)
         v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                 - beta(2)*v(icri  ,i,:,:,:,jj) &
                                 - beta(1)*v(icri+1,i,:,:,:,jj)
        enddo ! i
       enddo ! icri
      enddo ! jj
      call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)
      const = 1.0_KR2/sqrt(beta(1))
      do i = 1,nvhalf
       v(:,i,:,:,:,kDR+1)     = const*v(:,i,:,:,:,kDR+1)
      enddo ! i

! Need to have the vtemp vector for the gmresproj routine....

      do jj = 1,nvhalf
       vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1)
      enddo ! jj

! Rotations for newly formed hc() matrix.
      do jj = 1,kDR
       do i = jj+1,kDR+1
        amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
        con2 = 1.0_KR2/amags
        tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
        tv(2) = 0.0_KR2
        gca(1,i,jj) = sqrt(amags)/tv(1)
        gca(2,i,jj) = 0.0_KR2
        gsa(1,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
        gsa(2,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
        do j = jj,kDR
         tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                         + gsa(2,i,jj)*hc(2,i,j)
         tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                         - gsa(2,i,jj)*hc(1,i,j)
         tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                        + gsa(2,i,jj)*hc(2,jj,j)
         tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                        - gsa(2,i,jj)*hc(1,jj,j)
         hc(:,jj,j) = tv1(:)
         hc(:,i,j) = tv2(:)
        enddo ! j
        tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
        tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
        tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
        tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
        c(:,jj) = tv1(:)
        c(:,i) = tv2(:)
       enddo ! i
      enddo ! jj
      j = kDR
      icycle = icycle + 1
     endif
    enddo maindo

     !     if (myid==0) then
     !        print *, "evaluesgmresdr:-Dean wuz here"
     !        do ii =1,kDR
     !          print *, w(:,ii)
     !        enddo ! ii
     !     endif 

 end subroutine ppgmresdr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine ppmmgmresdr(rwdir,phi,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) with polynomial preconditioning
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: icycle, i, j, k, jp1, jj, p, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv
!!!!!!!!!!!!!!!!!!!!!!!!!!solving the coeffcients for the polynomial!!!!!!!!!!!!
!    integer(kind=KI), dimension(2)                          :: ipiv2!order of po                                                                     lynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w, work
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,ztmp, q
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht,matss
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3, hprint, t, hhh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    real(kind=KR2),   dimension(2,6,6)                      :: lsmat
!    real(kind=KR2),   dimension(2,6,1)                      :: cls
!    real(kind=KR2),   dimension(2,6)                        :: co!coefficients
!!!!!!!!!!!parameters in determing the preconditioning polynomial!!!!!!!!!!!!!!!

    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws, ev
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h , xt ,&
                                                               inter,inter1,tt
!!!!!!!!!!!Some intermediate parameters used for the preconditioning!!!!!!!!!!!!
!    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: w1, z1 ,z2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
!    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: v, vprime
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(2,kmaxGMRES,kmaxGMRES)      :: greal
   
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp, bopart
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
!    real(kind=KR2), dimension(6,ntotal,4,2,8,6) :: try!6 is the degree of P(A)*A
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
    integer(kind=KI) :: didmaindo,roww
    integer(kind=KI), dimension(nmaxGMRES)              :: sortev

    real(kind=KR), dimension(nxyzt,3,4)   :: realz2noise, imagz2noise


! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.



! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)
    p = 6!the order of the polynomial
    didmaindo = 0
    icycle = 1
    idag = 0
    ss = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR
    w1 = 0.0_KR2
    z1 = 0.0_KR2
    z2 = 0.0_KR2
    try = 0.0_KR2
 ! htemp = 0.0_KR
 ! site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf
 !      do idirac=1,4
 !        do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !              !print *, "site,icolorr =", site,icolorr
 !              !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps
!
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
!
!                call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!
!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
! To print single rhs source vector use ..
!
!               !irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!               !       print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! idirac
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock

!*****MORGAN'S STEP 1: Start.
! Compute r=phi-M*x and v=r/|r| and beta=|r|.
    idag = 1
  do i=1,nvhalf
   tt(:,i,:,:,:)=phi(:,i,:,:,:)
  enddo!i
  call Hdbletm(inter1,u,GeeGooinv,tt,idag,coact,kappa,iflag,bc, &
                vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)!inter1 will be the new
                                                       !  phi
    idag = 0 
   call Hdbletm(inter,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)!check M^dagger*M

   idag = 1
  call Hdbletm(h,u,GeeGooinv,inter,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)!check M^dagger*M
   idag = 0
!     try1 = 0.0_KR2
!     try2 = 0.0_KR2
!     try3 = 0.0_KR2
!     try4 = 0.0_KR2
    do i = 1,nvhalf
     r(:,i,:,:,:) = inter1(:,i,:,:,:) - h(:,i,:,:,:)
     v(:,i,:,:,:,1) = r(:,i,:,:,:)
     xt(:,i,:,:,:) = x(:,i,:,:,:)
    enddo ! i
    !Determine the polynomial.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nvhalf
     vprime(:,k,:,:,:,1) = inter1(:,k,:,:,:)
    enddo !k
    do i = 1,p
     idag=0
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     idag=1
     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag, &
                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    
    enddo !i
    
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
!      print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
     enddo!j
    enddo!i
        



   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),inter1(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)    
!    co = 0.0_KR2    
!    co(1,1) = 4
   if(myid==0) then
    do i=1,p
     print *, "i,result(:,i)=", i, co(:,i)
    enddo!i  
   endif!myid
!!!!
!!!!!Times the polynomial to the residue
        
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                         -co(2,1)*try(icri+1,k,:,:,:,1)
      y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                           +co(2,1)*try(icri,k,:,:,:,1)
     enddo!k
    enddo!icri

!    print *,"original component1:real=",w1(1,1,1,1,1)
!    print *,"original component2:imaginary=",w1(2,1,1,1,1)
 !   print *,"conditioned component1:real=",y(1,1,1,1,1)
!    print *,"conditioned component2:imaginary=",y(2,1,1,1,1)


   do i=1,p-1 
     idag=0
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z

     idag=1
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     do icri=1,5,2
      do k=1,nvhalf
       y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                             +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
       y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
   enddo!i
!     do k=1,nvhalf
!      r(:,k,:,:,:) = z2(:,k,:,:,:)
!     enddo!k
   
    
    do k=1,nvhalf
     v(:,k,:,:,:,1) = y(:,k,:,:,:,1) ! Use y=P(A)*r to generate V_(m+1)&Hbar_m
    enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))
    normnum = beta(1)

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.
    itercount = 0
    j = 0
    maindo: do




     if ( icycle > kcyclim) exit maindo

     j = j + 1
     jp1 = j + 1
     itercount = itercount + 1

!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
    try = 0.0_KR2
    y = 0.0_KR2
    do i=1,nvhalf
     try(:,i,:,:,:,1) = v(:,i,:,:,:,j)
    enddo!i
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                         -co(2,1)*try(icri+1,k,:,:,:,1)
      y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                           +co(2,1)*try(icri,k,:,:,:,1)
     enddo!k
    enddo!icri
    do i=1,p-1 
     idag=0
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     idag=1
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z

     do icri=1,5,2
      do k=1,nvhalf
       y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                             +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
       y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i

    idag=0
    call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,y(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)
    idag=1
    call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1,j
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)
      hc(:,i,j) = beta(:)

!     print *, "i,j, hc(:,i,j)=", i,j, hc(:,i,j)

      hc2(:,i,j) = hc(:,i,j)
      hc3(:,i,j) = hc(:,i,j)
      hcht(1,j,i) = hc(1,i,j)
      hcht(2,j,i) = -hc(2,i,j)
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                              - beta(1)*v(icri  ,k,:,:,:,i) &
                              + beta(2)*v(icri+1,k,:,:,:,i)
        v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                              - beta(2)*v(icri  ,k,:,:,:,i) &
                              - beta(1)*v(icri+1,k,:,:,:,i)
       enddo ! k
      enddo ! icri
     enddo ! i



     call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
     hc(1,jp1,j) = sqrt(beta(1))
     hc(2,jp1,j) = 0.0_KR2
     hc2(:,jp1,j) = hc(:,jp1,j)
     hc3(:,jp1,j) = hc(:,jp1,j)
     hcht(1,j,jp1) = hc(1,jp1,j)
     hcht(2,j,jp1) = -hc(2,jp1,j)
     const = 1.0_KR2/sqrt(beta(1))
     v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
     c(:,jp1) = 0.0_KR2
     c2(:,jp1) = c(:,jp1)


! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.
     if (icycle/=1) then
      do jj = 1,kDR
       do i = jj+1,kDR+1
        tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
               + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
        tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
               + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
        tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
               - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
        tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
               - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
        hc(:,jj,j) = tv1(:)
        hc(:,i,j) = tv2(:)
       enddo ! i
      enddo ! jj
      if (j>kDR+1) then
       do i = kDR+1,j-1
        tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
               + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
        tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
               + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
        tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
               - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
        tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
               - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
        hc(:,i,j) = tv1(:)
        hc(:,i+1,j) = tv2(:)
       enddo ! i
      endif
     elseif (j/=1) then
      do i = 1,j-1
       tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
              + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
       tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
              + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
       tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
              - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
       tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
              - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
       hc(:,i,j) = tv1(:)
       hc(:,i+1,j) = tv2(:)
      enddo ! i
     endif
     amags = hc(1,j,j)**2 + hc(2,j,j)**2
     tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
     tv(2) = 0.0_KR2
     gc(1,j) = sqrt(amags)/tv(1)
     gc(2,j) = 0.0_KR2
     con2 = gc(1,j)/amags
     gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
     gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
     hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
     hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
     hc(:,j+1,j) = 0.0_KR2
     tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
     tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
     tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
     tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
     c(:,j) = tv1(:)
     c(:,j+1) = tv2(:)
     do i = 1,j
      ss(:,i) = c(:,i)
     enddo ! i
     con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
     const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
     ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
     ss(1,j) = const
     if (j/=1) then
      do i = 1,j-1
       ir = j - i + 1
       irm1 = ir - 1
       do jj = 1,irm1
        const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
        ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
        ss(1,jj) = const
       enddo ! jj
       con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
       const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
       ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
       ss(1,irm1) = const
      enddo ! i
     endif

! Form the approximate new solution x = xt + V*ss.
     xb = 0.0_KR2
     do jj = 1,j
      do icri = 1,5,2 ! 6=nri*nc
       do i = 1,nvhalf
        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) + ss(1,jj)*v(icri  ,i,:,:,:,jj)&
                                                - ss(2,jj)*v(icri+1,i,:,:,:,jj)
        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) + ss(2,jj)*v(icri  ,i,:,:,:,jj)&
                                                + ss(1,jj)*v(icri+1,i,:,:,:,jj)
       enddo ! i
      enddo ! icri
     enddo ! jj
     do i = 1,nvhalf
      x(:,i,:,:,:) = xt(:,i,:,:,:) + xb(:,i,:,:,:)
     enddo ! i

! Define a small residual vector, srv = c-Hbar*ss, which corresponds to the
! kDR+1 column of the new V that will be formed.
     do i = 1,nDR+1
      srv(:,i) = c2(:,i)
     enddo ! i
     do jj = 1,nDR
      do i = 1,nDR+1
       srv(1,i) = srv(1,i) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
       srv(2,i) = srv(2,i) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
      enddo ! i
     enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
     if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=phi-M*x.

      idag=0
      call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,xb,idag,&
                   coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      call Hdbletm(h,u,GeeGooinv,test(:,:,:,:,:,1),idag,&
                   coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      do i = 1,nvhalf
       r(:,i,:,:,:) = r(:,i,:,:,:) - h(:,i,:,:,:)
      enddo ! i
!!!!!!!!!!!!!!!!!!!!!! r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k=1,nvhalf
!     w1(:,k,:,:,:) = r(:,k,:,:,:)
!    enddo!k
    
!    do icri=1,5,2
!     do k=1,nvhalf
!      y(icri,k,:,:,:) = co(1,1)*w1(icri,k,:,:,:)-co(2,1)*w1(icri+1,k,:,:,:)
!      y(icri+1,k,:,:,:) = co(1,1)*w1(icri+1,k,:,:,:)+co(2,1)*w1(icri,k,:,:,:)
!     enddo!k
!    enddo!icri

!    do k=1,nvhalf
!     z1(:,k,:,:,:) = w1(:,k,:,:,:)
!    enddo!k

!    do i=2,p
!     call Hdbletm(z1(:,:,:,:,:),u,GeeGooinv,z1(:,:,:,:,:),idag,coact,kappa, &
!                  iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
!                  lbd,iblv,MRT )  !z1=M*z
!     do icri=1,5,2
!      do k=1,nvhalf
!       y(icri,k,:,:,:) = y(icri,k,:,:,:)+co(1,i)*z1(icri,k,:,:,:) &
!                         -co(2,i)*z1(icri+1,k,:,:,:)
!       y(icri+1,k,:,:,:) = y(icri+1,k,:,:,:)+co(1,i)*z1(icri+1,k,:,:,:) &
!                           +co(2,i)*z1(icri,k,:,:,:)   !y=P(A)*r
!      enddo!k
!     enddo!icri
!    enddo!i

!    do k=1,nvhalf
!     r(:,k,:,:,:)=y(:,k,:,:,:)
!    enddo!k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
!        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr-norm",itercount,beta(1)/normnum
       close(unit=8,status="keep")
      endif

      if ((beta(1)/normnum)<resmax .and. itercount>=itermin) then
         if(didmaindo==0) then
         print *, "Everyone needs more resmax!"
         endif ! didmaindo
         exit maindo
      endif ! beta
      didmaindo=1
      do i = 1,nvhalf
       xt(:,i,:,:,:) = x(:,i,:,:,:)
      enddo ! i

!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)
      do i = 1,nDR
       em(:,i,1) = 0.0_KR
      enddo ! i
      em(1,nDR,1) = 1.0_KR2
      nrhs = 1





      call linearsolver(nDR,nrhs,hcht,ipiv,em)


      do i = 1,nDR
       hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                    - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                    - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
       hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                    + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                    - em(2,i,1)*hc2(2,nDR+1,nDR)**2
      enddo ! i

      ilo = 1
      ihi = nDR
      call hessenberg(hc2,nDR,ilo,ihi,tau,work)




      do jj = 1,nDR
       do i = jj,nDR
        z(:,i,jj) = hc2(:,i,jj)
       enddo ! i
      enddo ! jj
      call qgenerator(nDR,ilo,ihi,z,tau,work)
      ischur = 1
      call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

         if (myid==0) then
            print *, "------------"
            print *, "evaluesgmresdr:", w
            print *, "------------"
         endif 

!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan -- 
! so reorder the harmonic Ritz values and rearrange the Schur form.
      mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
      do i = 2,nDR
       mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
       is = 0
       ritzloop: do
        is = is + 1
        if (is>i-1) exit ritzloop
        if (mag(i)<mag(is)) then
         tval = mag(i)
         do ivb = i-1,is,-1
          mag(ivb+1) = mag(ivb)
         enddo ! ivb
         mag(is) = tval
         exit ritzloop
        endif
       enddo ritzloop
      enddo ! i
      do i = 1,nDR
       myselect(i) = .false.
       if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
      enddo ! i

      call orgschur(myselect,nDR,hc2,z,w,idis)

!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.
      do i = 1,kDR
       z(:,nDR+1,i) = 0.0_KR2
      enddo ! i
      do i = 1,nDR+1
       z(:,i,kDR+1) = srv(:,i)
      enddo ! i
      do jj = 1,kDR
       tv = 0.0_KR2
       do i = 1,nDR+1
        tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
        tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
       enddo ! i
       do i = 1,nDR+1
        z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
        z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
       enddo ! i
      enddo ! jj
      jj = nDR + 1
      rv = twonorm(z(:,:,kDR+1),jj)
      rv = 1.0_KR2/rv
      do i = 1,nDR+1
       z(:,i,kDR+1) = rv*z(:,i,kDR+1)
      enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.
      do jj = 1,kDR
       do ii = 1,nDR+1
        ws(:,ii,jj) = 0.0_KR2
       enddo ! ii
       do ii = 1,nDR
        do i = 1,nDR+1
         ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                 - z(2,ii,jj)*hc3(2,i,ii)
         ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                 + z(2,ii,jj)*hc3(1,i,ii)
        enddo ! i
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hcnew(:,ii,jj) = 0.0_KR2
        do i = 1,nDR+1
         hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                         + z(2,i,ii)*ws(2,i,jj)
         hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                         - z(2,i,ii)*ws(1,i,jj)
        enddo ! i
       enddo ! ii
      enddo ! jj

      do jj = 1,nDR
       do ii = 1,nDR
        hcht(:,ii,jj) = 0.0_KR2
        hc2(:,ii,jj) = 0.0_KR2
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hc(:,ii,jj) = hcnew(:,ii,jj)
        hc2(:,ii,jj) = hcnew(:,ii,jj)
        hc3(:,ii,jj) = hcnew(:,ii,jj)
       enddo ! ii
       do ii = 1,kDR+1
        hcht(1,jj,ii) = hcnew(1,ii,jj)
        hcht(2,jj,ii) = -hcnew(2,ii,jj)
       enddo ! ii
      enddo ! jj
      do ii = 1,kDR+1
       c(:,ii) = 0.0_KR2
       do i = 1,nDR+1
        c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i) + z(2,i,ii)*srv(2,i)
        c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i) - z(2,i,ii)*srv(1,i)
       enddo ! i
       c2(:,ii) = c(:,ii)
      enddo ! ii

      vtemp = 0.0_KR

      do ibleo = 1,8
       do ieo = 1,2
        do id = 1,4
         do i = 1,nvhalf
          do jj = 1,kDR+1
           vt(:,jj) = 0.0_KR2
           do k = 1,nDR+1
            do icri = 1,5,2 ! 6=nri*nc
             vt(icri  ,jj) = vt(icri  ,jj) &
                           + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
             vt(icri+1,jj) = vt(icri+1,jj) &
                           + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
            enddo ! icri
           enddo ! k
          enddo ! jj
          do jj = 1,kDR+1
           v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
           vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
          enddo ! jj
         enddo ! i
        enddo ! id
       enddo ! ieo
      enddo ! ibleo
         !if(myid==0) then
         ! print *, "vtemp in this shizat-only 1:kDR", vtemp
         !endif ! myid

!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.
      do jj = 1,kDR
       call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
       do icri = 1,5,2 ! 6=nri*nc
        do i = 1,nvhalf
         v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                 - beta(1)*v(icri  ,i,:,:,:,jj) &
                                 + beta(2)*v(icri+1,i,:,:,:,jj)
         v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                 - beta(2)*v(icri  ,i,:,:,:,jj) &
                                 - beta(1)*v(icri+1,i,:,:,:,jj)
        enddo ! i
       enddo ! icri
      enddo ! jj
      call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)
      const = 1.0_KR2/sqrt(beta(1))
      do i = 1,nvhalf
       v(:,i,:,:,:,kDR+1)     = const*v(:,i,:,:,:,kDR+1)
      enddo ! i

! Need to have the vtemp vector for the gmresproj routine....

      do jj = 1,nvhalf
       vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1)
      enddo ! jj

! Rotations for newly formed hc() matrix.
      do jj = 1,kDR
       do i = jj+1,kDR+1
        amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
        con2 = 1.0_KR2/amags
        tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
        tv(2) = 0.0_KR2
        gca(1,i,jj) = sqrt(amags)/tv(1)
        gca(2,i,jj) = 0.0_KR2
        gsa(1,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
        gsa(2,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
        do j = jj,kDR
         tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                         + gsa(2,i,jj)*hc(2,i,j)
         tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                         - gsa(2,i,jj)*hc(1,i,j)
         tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                        + gsa(2,i,jj)*hc(2,jj,j)
         tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                        - gsa(2,i,jj)*hc(1,jj,j)
         hc(:,jj,j) = tv1(:)
         hc(:,i,j) = tv2(:)
        enddo ! j
        tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
        tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
        tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
        tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
        c(:,jj) = tv1(:)
        c(:,i) = tv2(:)
       enddo ! i
      enddo ! jj
      j = kDR
      icycle = icycle + 1
     endif
    enddo maindo

     !     if (myid==0) then
     !        print *, "evaluesgmresdr:-Dean wuz here"
     !        do ii =1,kDR
     !          print *, w(:,ii)
     !        enddo ! ii
     !     endif 

 end subroutine ppmmgmresdr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 subroutine gmresdrshift(rwdir,b,xshift,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,mvp,gdr)
! GMRES-DR(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=b for the vector x.
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresdrshift.

! This is DRSHIFT
 
    use shift 

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount, mvp
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
    ! real(kind=KR2),   intent(inout), dimension(:,:,:,:,:,:) :: vtemp
   
 
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp

    integer(kind=KI) :: icycle, i, j, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa,tempnrhs

    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv, alpha
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si  
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    ! real(kind=KR2),   dimension(1500,nshifts)               :: gdr
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    ! real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
! NOTE ~ if you wish to use this routine efficently in memory space
!        comment the ::xt,xbt line and uncomment the same line in
!        qqcd/cfgsprops/quark/shift.f90
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
!   real(kInd=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha
!   real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: beg
    integer(kind=KI)                                        :: isignal, isite, icolorir, idirac,&
                                                               iblock, site, printnum, icolorr, irow, &
                                                               ishift, ilat


! This is still DRSHIFT

!   if(myid==0) then
!     print *, "Made it to DRSHIFT"
!   endif ! myid


! Shift sigmamu to base mu above  (mtmqcd(1,2))

!   sigma = 0.0_KR

!   sigma(1) =  0.00_KR/kappa(1)**2
!   sigma(2) = -0.002_KR/kappa(1)**2
!   sigma(3) = -0.006_KR/kappa(1)**2

!   sigma(1) =  0.00_KR/kappa(1)**2
!   sigma(2) = -0.099_KR/kappa(1)**2
!   sigma(3) = -0.199_KR/kappa(1)**2

!   sigma(1) =  0.00_KR/kappa(1)**2
!   sigma(2) = -0.30_KR/kappa(1)**2
!   sigma(3) = -0.50_KR/kappa(1)**2

!   call twistedsigma

!   do ishift = 1,nshifts
!     sigma(ishift) =  sigmamu(ishift)/kappa(1)**2
!   enddo ! ishift

! Leading dimensions used for LAPACK routines
    ldh = nmaxGMRES + 1
    ldz = nmaxGMRES + 1
    ldhcht = nmaxGMRES
    lwork = nmaxGMRES
    lzwork = nmaxGMRES+1
 
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)
    icycle = 1
    idag = 0
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR
! Print statment init - temporary 
    beta1 = 0.0_KR
    beta2 = 0.0_KR
    beta3 = 0.0_KR


! Initialize cmult

    do is = 1,nshifts
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0
    enddo ! is 

! Verify that nk1 has been set to a sufficiently large value.

    ! if (nk1<nkappa) then
    !   open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
    !        ,form="formatted")
    !   write(unit=8,fmt=*) "nk1 is less than nkappa: nk1,nkappa = ", nk1,nkappa 
    !   close(unit=8,status="keep")
    !endif

!*****MORGAN'S STEP 1: Start.
! Compute r=b-M*x and v=r/|r| and beta=|r|.
! HEY! Does this routine assume kappa(2) is mu value?
    ! if (myid==0) then
    !   print *, "before first, h ="
    !   call checkNonZero(h(:,:,:,:,:),ntotal)
    ! endif

   htemp = 0.0_KR
   site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf 
 !      do idirac=1,4
 !        do icolorir=1,5,2
  !              site = ieo + 16*(isite - 1) + 2*(iblock - 1) 
 !               icolorr = icolorir/2 +1
                !print *, "site,icolorr =", site,icolorr
                !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps

!                getemp = 0.0_KR
!                getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
!                getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
        
!                call gamma5mult(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)

! To print single rhs source vector use ..

 !               irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
 !                      print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)

 !            enddo ! icolorir
 !         enddo ! idirac 
 !      enddo ! isite
 !   enddo ! ieo
 ! enddo ! iblock

    call gamma5mult(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    ! if (myid==0) then
    !   rint *, "after first, h ="
    !   call checkNonZero(h(:,:,:,:,:),ntotal)
    ! endif

    do i = 1,nvhalf
       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)
       v(:,i,:,:,:,1) = r(:,i,:,:,:)
       do is=1,nshifts
          xt(:,i,:,:,:,is) = xshift(:,i,:,:,:,1)
       enddo ! is
    enddo ! i

! Copy the initial resiudal into the initial residual for each shift

    do is=1,nshifts
     rshift(:,:,:,:,:,is) = r(:,:,:,:,:)  
    enddo ! is

    ! if (myid==0) then
    !    print *, "r at top ="
    !    call checkLargeValue(r,nvhalf)
    ! endif

    beta = 0.0_KR2
    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))

!   if (myid==0) then
!      print *, "foobeta", beta(1)
!   endif

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
    ! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.

    itercount = 0
    j = 0
    mvp = 0

    maindo: do

     if ( icycle > kcyclim) exit maindo

       j = j + 1
       jp1 = j + 1
       itercount = itercount + 1

       ! if (myid==0) then
       !   print *, "maindo: j=",j
       ! endif
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.

       ! if (myid==0) then
       !   print *, "STEP 2A, rank=",myid
       ! endif

! HEY! Might have to shift on kappa(2)
! (or the vectors v(j) and v(jp1))
       call gamma5mult(v(:,:,:,:,:,jp1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT)
       mvp = mvp + 1
 
       do i = 1,j

          ! if (myid==0) then
          !   print *,"before 2.beta", beta(1), beta(2)
          !endif

          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif

          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
          hcht(1,j,i) = hc(1,i,j)
          hcht(2,j,i) = -hc(2,i,j)

          do icri = 1,5,2 ! 6=nri*nc 
             do k = 1,nvhalf
                v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,k,:,:,:,i) &
                                        + beta(2)*v(icri+1,k,:,:,:,i)
                v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,k,:,:,:,i) &
                                        - beta(1)*v(icri+1,k,:,:,:,i)
             enddo ! k
          enddo ! icri
       enddo ! i

       ! ... shift the diagonal entries of hc for first shift localKappa(2)
       ! (sigma(2))

       ! if (myid==0) then
       !   print *, "around 3.beta sigma(1)=", sigma(1)
       ! endif

       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)
       hcht(1,j,j) = hc(1,j,j)
       hcht(2,j,j) = -hc(2,j,j)

       ! ... first shift on diag matrix complete
 
       ! if (myid==0) then
       !   print *,"before 3.beta", beta(1), beta(2)
       ! endif
     
       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)

       ! if (myid==0) then
       !   print *,"after 3.beta", sqrt(beta(1)), beta(2)
       ! endif

       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)
       hcht(1,j,jp1) = hc(1,jp1,j)
       hcht(2,j,jp1) = -hc(2,jp1,j)
       const = 1.0_KR2/sqrt(beta(1))
       v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)

       ! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.

      !if (myid == 0) then
      !  print *, "c2 =", c2
      !  print *, "hc2 =", hc2(:,1:j+1,1:j)  
      !endif

       if (icycle/=1) then

          do jj = 1,kDR
             do i = jj+1,kDR+1
                tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
                       + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
                tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
                       + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
                tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
                       - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
                tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
                       - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
                hc(:,jj,j) = tv1(:)
                hc(:,i,j) = tv2(:)
             enddo ! i
          enddo ! jj
          if (j>kDR+1) then
             do i = kDR+1,j-1
                tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                       + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
                tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                       + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
                tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                       - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
                tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                       - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
                hc(:,i,j) = tv1(:)
                hc(:,i+1,j) = tv2(:)
             enddo ! i
          endif ! (j>kDR+1)

       elseif (j/=1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (icycle/=1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

       con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
       const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
       ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
       ss(1,j) = const

       if (j/=1) then
          do i = 1,j-1
             ir = j - i + 1
             irm1 = ir - 1
             do jj = 1,irm1
                const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
                ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
                ss(1,jj) = const
             enddo ! jj
             con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
             const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
             ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
             ss(1,irm1) = const
          enddo ! i
       endif ! (j/=1)

       ! ... Define new variable d to assist in shifting masses
       ! d is the "short" solution vector to small problem

       d(1,:,1) = ss(1,:) 
       d(2,:,1) = ss(2,:)

       ! ... form srv for first solution vector

       st = 0.0_KR

      !if (myid == 0 ) then
      !  print *, "first d=", d
      !endif

       do ii = 1,j+1
          do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1) 
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
          enddo ! jj
          srv(1,ii,1) = c2(1,ii) - st(1,ii)
          srv(2,ii,1) = c2(2,ii) - st(2,ii)
       enddo ! ii

       beta(1) = 0.0_KR

       do jj = 1,j+1
          beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
       enddo ! jj

       ! ... form gdr

        gdr(mvp,1) = sqrt(beta(1))
        !if (myid==0) then
        !  print *, "j =", j
        !  print *, "gdr for first shift=", gdr(mvp,1)
        !endif

       ! ... To keep subspaces parallel for different shifts vector
       !     needs to be rotated. Use a QR factorization to do the
       !     ortognal rotation.
       !     BIG Shifting loop

       do is = 2,nshifts

          do ii=1,j+1
             do jj=1,j
                hcs(1,ii,jj) = hc2(1,ii,jj)
                hcs(2,ii,jj) = hc2(2,ii,jj)
                hcs2(1,ii,jj) = hc2(1,ii,jj)
                hcs2(2,ii,jj) = hc2(2,ii,jj)
             enddo ! jj
          enddo ! ii

          do jj=1,j
             ! HEY! replaced :'s with 1's (first index), e.g., hcs(:,jj,jj)
             hcs(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is)
             hcs2(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is) 
          enddo ! jj

          ! Copy the 12complex arrays into true complex arrays for use with
          ! lapack routines

          call real2complex_mat(hcs, nmaxGMRES+1, nmaxGMRES+1, zhcs)

          !if (myid == 0) then
          !   print *, "before zgeqrf"
          !endif

          call zgeqrf(j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)

          !if (myid == 0) then
          !   print *, "after zgeqrf"
          !endif

          ! ... store R (upper triangular) in rr

          do ii=1,jp1 
             do jj=1,j
                zrr(ii,jj) = zhcs(ii,jj)
             enddo ! jj
          enddo ! ii

          call zungqr(j+1,j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)

          ! Copy the complex zhcs array back to hcs 

          call complex2real_mat(zhcs, nmaxGMRES+1, nmaxGMRES+1, hcs)

          ! ... hcs after this call is the qq part of the qr factorization

          ! Now zero out crot (keeps shifts parrallel) and srvrot

           do ii=1,jp1
              zcrot(ii)=0.0_KR
              zsrvrot(ii)=0.0_KR
           enddo ! ii

          do ii=1,jp1
             do jj=1,jp1

                ztemp1 =DCMPLX(cmult(1,is), cmult(2,is))
                ztemp2 = DCMPLX(c2(1,jj), c2(2,jj))
                ztemp3 = DCMPLX(srv(1,jj,1), srv(2,jj,1))

                ! if (myid == 0) then
                !   print *, "ii, jj, ztemp1, ztemp2  = ", ii, jj, ztemp1, ztemp2
                ! endif

                zcrot(ii) = zcrot(ii) + ztemp1 * CONJG(zhcs(jj,ii)) * ztemp2
                zsrvrot(ii) = zsrvrot(ii) + CONJG(zhcs(jj,ii)) * ztemp3

                ! if ((myid/=myid) .and. (ABS(zcrot(ii)) < 1.0d-004)) then
                !   print *, "HEY!"
                !   print *, "ii, jj, ztemp1, ztemp2  = ", ii, jj, ztemp1, ztemp2
                !   print *, "zhcs = "
                !   call printArray(zhcs)
                !   call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                ! endif

             enddo ! jj
          enddo ! ii

          ! ... construct alpha

          ! if ((myid == 0) .and. (is==3)) then
          !   print *, "jp1, zcrot(jp1), zsrvrot(jp1) = ", jp1, zcrot(jp1), zsrvrot(jp1)
          ! endif
             
          zalpha = zcrot(jp1)/zsrvrot(jp1)

          alph(1, is) = REAL(zalpha)
          alph(2, is) = AIMAG(zalpha)

          ! if (myid==0)  then
          !   print *,"is, zalpha, alph = ", is, zalpha, alph(:,is)
          ! endif

          do jj=1,j
             ztemp1 = zalpha * zsrvrot(jj)
             cmas(1,jj) = REAL(zcrot(jj)) - REAL(ztemp1)    ! alpha*srvrot(1,jj)
             cmas(2,jj) = AIMAG(zcrot(jj)) - AIMAG(ztemp1)   ! alpha*srvrot(2,jj)
          enddo ! jj

          ! if ((myid==0) .and. (is==3))  then
          !   print *,"cmas = ", cmas
          ! endif

          ! ... solve linear eqns problem d(1:j,is)=rr(1:j,1:j)\cmas

          do ii=1,j
             d(1,ii,is)=cmas(1,ii)
             d(2,ii,is)=cmas(2,ii)
          enddo ! ii

          call real2complex_mat(d,nmaxGMRES,nshifts,zd)

          ! do i=1,2
          !   d(1,j,is) = d(1,j,is)/rr(i,j,j)
          !   d(2,j,is) = d(2,j,is)/rr(i,j,j)       
          ! enddo ! i

          zd(j,is) = zd(j,is)/zrr(j,j)

          if (j /= 1) then 
             do i=1,j-1
                ir = j-i +1
                irm1 = ir -1
                call zaxpy(irm1,-zd(ir,is),zrr(1,ir),1,zd(1,is),1)
                zd(irm1,is) = zd(irm1,is)/zrr(irm1,irm1)
             enddo ! i
          endif ! j/=1

          call complex2real_mat(zd, nmaxGMRES, nshifts, d)

          do ii=1,jp1
             st(1,ii) = 0.0_KR 
             st(2,ii) = 0.0_KR
             srvis(1,ii) = 0.0_KR
             srvis(2,ii) = 0.0_KR
          enddo ! ii

          do ii=1,jp1
             do jj=1,j
                st(1,ii) = st(1,ii) + hcs2(1,ii,jj)*d(1,jj,is) &
                                    - hcs2(2,ii,jj)*d(2,jj,is)
                st(2,ii) = st(2,ii) + hcs2(1,ii,jj)*d(2,jj,is) &
                                    + hcs2(2,ii,jj)*d(1,jj,is)
             enddo ! jj
             srvis(1,ii) = cmult(1,is)*c2(1,ii)-cmult(2,is)*c2(2,ii) &
                         - st(1,ii)
             srvis(2,ii) = cmult(1,is)*c2(2,ii)+cmult(2,is)*c2(1,ii) &
                         - st(2,ii)
          enddo ! ii

          ! ... form the norm of srvis and put in gdr

          beta(1) =0.0_KR

          do jj=1,j+1
             beta(1) = beta(1) + srvis(1,jj)*srvis(1,jj) &
                               + srvis(2,jj)*srvis(2,jj)
          enddo ! jj

          gdr(mvp,is) = sqrt(beta(1))

          !if ((myid == 0) .and. (is==2))  then
          !   print *, "is, j =", is, j 
          !   print *, "gdr in shift=", gdr(:,is) 
          !   print *, "itercount = ", itercount
          !   print *, "zcrot = ", zcrot(jp1)
          !   print *, "zsrvrot = ", zsrvrot(jp1)
          !endif

       enddo ! BIG is loop
 
       ! Form the approximate new solution x = xt + xb.
       ! First zero xb then xb = V*d(1,2,:,is), and x =xt +xb
       ! ... Need to keep track of each solution for each shift.
       !     Loop over shifts while creating the soln vector.

       if (j>=nmaxGMRES) then 
          do is = 1,nshifts
             cmult(1,is) = alph(1,is) 
             cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i
           
             xb = 0.0_KR2
             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj
             do i = 1,nvhalf
                xshift(:,i,:,:,:,is) = xt(:,i,:,:,:,is) + xb(:,i,:,:,:)
                xbt(:,i,:,:,:,is) = xb(:,i,:,:,:)
             enddo ! i
          enddo ! is
       endif ! (j>=nmaxGMRES)
  
       ! Define a small residual vector, srv = c-Hbar*ss, which corresponds
       ! to the kDR+1 column of the new V that will be formed.

        do i = 1,nDR+1
          srv(:,i,1) = c2(:,i)
       enddo ! i

       do jj = 1,nDR
          do i = 1,nDR+1
             srv(1,i,1) = srv(1,i,1) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
             srv(2,i,1) = srv(2,i,1) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
          enddo ! i
       enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
        if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=b-M*x.

          ! if (myid==0) then
          !  print *, "STEP 2B, rank=", myid
          !endif

          ! if (myid == 0) then
          !   print *,"ntotal=", ntotal
          !   call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          ! endif

          do is=1,nshifts
           
             !if ((is==2) .and. (myid==0)) then
             !   print *, "is, before rshift =", is
             !   ! call checkLargeValue(h(:,:,:,:,:),ntotal)
             !   call checkLargeValue(rshift(:,:,:,:,:,is),nvhalf)
             !endif

             call gamma5mult(h,u,GeeGooinv,xbt(:,:,:,:,:,is),idag,coact,kappa, &
                             iflag,bc,vecbl, vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

            do i = 1,nvhalf

               rshift(:,i,:,:,:,is) = rshift(:,i,:,:,:,is) - h(:,i,:,:,:) &
                                                +sigma(is)*xbt(:,i,:,:,:,is) 
            enddo ! i

             !if ((is==2) .and. (myid==0)) then
             !   print *, "is, after rshift =", is
             !   call checkLargeValue(rshift(:,:,:,:,:,is),nvhalf)
             !endif

            beta = 0.0_KR

            call vecdot(rshift(:,:,:,:,:,is),rshift(:,:,:,:,:,is),beta,MRT2)
          ! ~NOTE might have to take sqrt of beta 1 
            beta(1) = sqrt(beta(1))
            betashift(1,is) = beta(1)
          enddo ! is

!         beta = beta1

          ! if (myid==0) then
          !   print *, "icycle,j,norm = ", icycle,j,sqrt(beta(1))
          ! endif ! print

          if (myid==0) then
             open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                  form="formatted",status="old",position="append")
! I commented out the betashift of higher dimension because I wanted to do 
! only one mass.
             write(unit=8,fmt=*) "gmresdrshift",itercount,betashift(1,1) !,betashift(1,2),betashift(1,3),&
                                                                         !  betashift(1,4)
             close(unit=8,status="keep")
          endif
        ! if (myid ==0) then
        !    open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
        !         form="formatted",status="old",position="append")
        !    write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
        !          betashift(1,1) !, " ", betashift(1,2), " ", betashift(1,3), " ", betashift(1,4)
        !    close(unit=8,status="keep")
        ! endif ! myid
           

          if (betashift(1,1)<resmax .and. itercount>=itermin) exit maindo

          do is=1,nshifts
           do i = 1,nvhalf
              xt(:,i,:,:,:,is) = xshift(:,i,:,:,:,is)
           enddo ! i
          enddo ! is
!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)

          ! if (myid==0) then
          !   print *, "STEP 2C, rank=", myid
          ! endif

          do i = 1,nDR
             em(:,i,1) = 0.0_KR
          enddo ! i
          em(1,nDR,1) = 1.0_KR2
          !nrhs = 1
           tempnrhs = 1

          call linearsolver(nDR,tempnrhs,hcht,ipiv,em)

          do i = 1,nDR
             hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                          - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                          - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
             hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                          + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                          - em(2,i,1)*hc2(2,nDR+1,nDR)**2
          enddo ! i

          ilo = 1
          ihi = nDR
          call hessenberg(hc2,nDR,ilo,ihi,tau,work)
          do jj = 1,nDR
             do i = jj,nDR
                z(:,i,jj) = hc2(:,i,jj)
            enddo ! i
          enddo ! jj
          call qgenerator(nDR,ilo,ihi,z,tau,work)
          ischur = 1
          call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

          ! Print evalues


!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan --
! so reorder the harmonic Ritz values and rearrange the Schur form.

          ! if (myid==0) then
          !   print *, "STEP 3, rank=", myid
          ! endif

          mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
          do i = 2,nDR
             mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
             is = 0
             ritzloop: do
                is = is + 1
                if (is>i-1) exit ritzloop
                if (mag(i)<mag(is)) then
                   tval = mag(i)
                   do ivb = i-1,is,-1
                      mag(ivb+1) = mag(ivb)
                   enddo ! ivb
                   mag(is) = tval
                   exit ritzloop
                endif ! (mag(i)<mag(is))
             enddo ritzloop
          enddo ! i

          do i = 1,nDR
             myselect(i) = .false.
             if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
          enddo ! i

          call orgschur(myselect,nDR,hc2,z,w,idis)
 
!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.

          ! print *, "STEP 4, rank=", myid

          do i = 1,kDR
             z(:,nDR+1,i) = 0.0_KR2
          enddo ! i
          do i = 1,nDR+1
             z(:,i,kDR+1) = srv(:,i,1)
          enddo ! i
          do jj = 1,kDR
             tv = 0.0_KR2
             do i = 1,nDR+1
                tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
                tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
             enddo ! i
             do i = 1,nDR+1
                z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
                z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
             enddo ! i
          enddo ! jj
          jj = nDR + 1
          rv = twonorm(z(:,:,kDR+1),jj)
          rv = 1.0_KR2/rv
          do i = 1,nDR+1
             z(:,i,kDR+1) = rv*z(:,i,kDR+1)
          enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.

          ! print *, "STEP 5, rank=", myid

          do jj = 1,kDR
             do ii = 1,nDR+1
                ws(:,ii,jj) = 0.0_KR2
             enddo ! ii
             do ii = 1,nDR
                do i = 1,nDR+1
                   ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                           - z(2,ii,jj)*hc3(2,i,ii)
                   ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                           + z(2,ii,jj)*hc3(1,i,ii)
                enddo ! i
             enddo ! ii
          enddo ! jj
          do jj = 1,kDR
             do ii = 1,kDR+1
                hcnew(:,ii,jj) = 0.0_KR2
                do i = 1,nDR+1
                   hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                                   + z(2,i,ii)*ws(2,i,jj)
                   hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                                   - z(2,i,ii)*ws(1,i,jj)
                enddo ! i
             enddo ! ii
          enddo ! jj
          do jj = 1,nDR
             do ii = 1,nDR
                hcht(:,ii,jj) = 0.0_KR2
                hc2(:,ii,jj) = 0.0_KR2
             enddo ! ii
          enddo ! jj
          
    !     if (myid == 0) then
    !      print *, "hcnew =", hcnew(:,1:3,1:2)
    !      print *, "z =", z(:,1:4,1:3)
    !     endif

          do jj = 1,kDR
             do ii = 1,kDR+1
                hc(:,ii,jj) = hcnew(:,ii,jj)
                hc2(:,ii,jj) = hcnew(:,ii,jj)
                hc3(:,ii,jj) = hcnew(:,ii,jj)
             enddo ! ii
             do ii = 1,kDR+1
                hcht(1,jj,ii) = hcnew(1,ii,jj)
                hcht(2,jj,ii) = -hcnew(2,ii,jj)
             enddo ! ii
          enddo ! jj
          do ii = 1,kDR+1
             c(:,ii) = 0.0_KR2
             do i = 1,nDR+1
                c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i,1) + z(2,i,ii)*srv(2,i,1)
                c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i,1) - z(2,i,ii)*srv(1,i,1)
             enddo ! i
             c2(:,ii) = c(:,ii)
          enddo ! ii
       
         !if (myid == 0) then
! PRINTC2
         !   print *, "c2 = ", c2
        !endif

! Zero out vtemp that is to be passed to gmresproject.
          
          vtemp = 0.0_KR

          do ibleo = 1,8
             do ieo = 1,2
                do id = 1,4
                   do i = 1,nvhalf
                      do jj = 1,kDR+1
                         vt(:,jj) = 0.0_KR2
                         do k = 1,nDR+1
                            do icri = 1,5,2 ! 6=nri*nc
                               vt(icri  ,jj) = vt(icri  ,jj) &
                                             + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                                             - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
                               vt(icri+1,jj) = vt(icri+1,jj) &
                                             + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                                             + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
                            enddo ! icri
                         enddo ! k
                      enddo ! jj
                      do jj = 1,kDR+1
                         v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
                         vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
                      enddo ! jj
                   enddo ! i
                enddo ! id
             enddo ! ieo
          enddo ! ibleo

! Above a  copy of the first k-Harmonic Ritz vectors in vt
! are copied into vtemp for the projections section for
! subsequent right hand sides.

!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.

          ! print *, "STEP 6, rank=", myid

          do jj = 1,kDR
             call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
             do icri = 1,5,2 ! 6=nri*nc
                do i = 1,nvhalf
                   v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                           - beta(1)*v(icri  ,i,:,:,:,jj) &
                                           + beta(2)*v(icri+1,i,:,:,:,jj)
                   v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                           - beta(2)*v(icri  ,i,:,:,:,jj) &
                                           - beta(1)*v(icri+1,i,:,:,:,jj)
                enddo ! i
             enddo ! icri
          enddo ! jj

          call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)

          const = 1.0_KR2/sqrt(beta(1))
          do i = 1,nvhalf
             v(:,i,:,:,:,kDR+1) = const*v(:,i,:,:,:,kDR+1)
          enddo ! i
! Need a copy of v for the projection section. This helps
! multi-shifting and  multiple right hand sides occur together in
! gmresproject.
          do jj = 1,nvhalf
           vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1) 
          enddo ! jj

!do jj = 1,kDR+1
! do iblock =1,8
!   do ieo = 1,2
!     do isite=1,nvhalf 
!       do idirac=1,4
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1

! To print single rhs source vector use ..

!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, vtemp(icolorir,isite,idirac,ieo,iblock,jj), vtemp(icolorir+1,isite,idirac,ieo,iblock,jj)

!             enddo ! icolorir
!          enddo ! idirac 
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock
! enddo ! jj

! Will pass vt into subroutine gmresproject. 
! Rotations for newly formed hc() matrix.

          do jj = 1,kDR
             do i = jj+1,kDR+1
                amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
                con2 = 1.0_KR2/amags
                tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
                tv(2) = 0.0_KR2
                gca(1,i,jj) = sqrt(amags)/tv(1)
                gca(2,i,jj) = 0.0_KR2
                gsa(1,i,jj) = gca(1,i,jj)*con2 &
                            *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
                gsa(2,i,jj) = gca(1,i,jj)*con2 &
                            *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
                do j = jj,kDR
                   tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                                   + gsa(2,i,jj)*hc(2,i,j)
                   tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                                   - gsa(2,i,jj)*hc(1,i,j)
                   tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                                  + gsa(2,i,jj)*hc(2,jj,j)
                   tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                                  - gsa(2,i,jj)*hc(1,jj,j)
                   hc(:,jj,j) = tv1(:)
                   hc(:,i,j) = tv2(:)
                enddo ! j
                tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
                tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
                tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
                tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
                c(:,jj) = tv1(:)
                c(:,i) = tv2(:)
             enddo ! i
          enddo ! jj

          j = kDR
          icycle = icycle + 1
! NOTE ~ end endif
       endif ! (j>=nDR)
     
    enddo maindo

          if (myid==0) then
             print *, "evaluesgmresdrshift:"
             do ii =1,nDR
               print *, w(:,ii)
             enddo ! ii
          endif 
 

 end subroutine gmresdrshift

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine gmresproject(rwdir,b,xshift,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)
! GMRES-PROJECT(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
!
! Gmresproject projects over the approximate eigenvectors found in the 
! deflation section of gmresdr (gmresdrshift). These rojected evectors
! are used in the basis to solve the following right-hand side with the
! usual extraction methods.  
!
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresproject.

! This is PROJ
 
    use shift

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    ! real(kind=KR2),   intent(in),    dimension(:,:,:)     :: hcnew
    ! real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: vtemp
!   real(kind=KR2),   intent(inout), dimension(:,:,:,:,:) :: beg
    integer(kind=KI), intent(in)                          :: isignal
    integer(kind=KI), intent(inout)                       :: mvp
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
 
    integer(kind=KI) :: icycle, i, j, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa, ifreq,temprhs !,nrhs
 
    ! real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)     :: xshift
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               con2, rv, amags
    real(kind=KR2)                                          :: rn, rinit, rnnn 
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info,rnt
    real(kind=KR2),   dimension(2)                          :: beta, alpha, &
                                                               tv, tv1 ,tv2
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(kcyclim,nshifts)            :: rnale
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhc2
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    ! real(kind=KR2),   dimension(2,nmaxGMRES)                :: gc, gs
    real(kind=KR2),   dimension(2)                          :: gam

    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha

    integer(kind=KI)                                        ::  isite, icolorir, idirac, iblock, &
                                                                site, icolorr, irow, ishift

! This is still PROJ 

 
! Shift sigmamu to base mu above (mtmqcd(1,2))
 
    sigma = 0.0_KR

! DEAN ~ HEY! I need to take out the sigma(1) part because I am not 
!        shifting at all and the residuals should be r = b-Ax
 
! init all x to 0

      xshift = 0.0_KR
 
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)


    ldh = nmaxGMRES+ 1
    ldz = nmaxGMRES+ 1
    lzwork = nmaxGMRES+1

!   if (myid==0) then
!     print *, "nDR = ", nDR
!     print *, "kDR = ", kDR
!   endif ! myid

    icycle = 1
    ifreq = 1
    idag = 0
    alph = 0.0_KR
    zcrot = 0.0_KR
    zsrvrot = 0.0_KR
    cmas = 0.0_KR
    srv = 0.0_KR
    srvis = 0.0_KR
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcs = 0.0_KR
    hcht = 0.0_KR
    zrr = 0.0_KR
    c2 = 0.0_KR
    c = 0.0_KR
    v = 0.0_KR
    xb = 0.0_KR
    d = 0.0_KR
    mvp = 0


! Initialize cmult
     
!   do is = 1,nshifts
       is =1
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0.0_KR
!   enddo ! is
 
    ! do is=1,nshifts
    !   xshift(:,:,:,:,:,is) = x(:,:,:,:,:,:)
    ! enddo ! is
 
! Compute r=b-M*x and v=r/|r| and beta=|r|.

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta,MRT2)
   !if (myid == 0) then
   ! print *, "norm o b in project =", sqrt(beta(1))
   !endif

! do iblock =1,8
!   do ieo = 1,2
!     do idirac=1,4
!       do isite=1,nvhalf 
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1
!
! To print single rhs source vector use ..
!
!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, b(icolorir,isite,idirac,ieo,iblock), b(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! isite
!       enddo ! idirac 
!    enddo ! ieo
!  enddo ! iblock

    call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

    mvp = mvp + 1
! NOTE ~ I don't think that I need this ....

    do ii = 1,nvhalf
        rshift(:,ii,:,:,:,1) = b(:,ii,:,:,:) - h(:,ii,:,:,:) +sigma(1)*xshift(:,ii,:,:,:,1)
    enddo ! ii

! For error correction put error in one direction and use gmresproject to 
! correct and solve the later right hand sides.

    if (isignal == 1) then
!     do is=1,nshifts
       rshift(:,:,:,:,:,1) = beg(:,:,:,:,:)
!     enddo ! is
    else
      beg(:,:,:,:,:) = rshift(:,:,:,:,:,1)
     !beg(:,:,:,:,:) = b(:,:,:,:,:)
    endif
!   call checkNonZero(beg,nvhalf)
! Create rinit for logic passing used in gmresprojet    

   beta = 0.0_KR

   call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

   rinit = sqrt(beta(1))
   rn = rinit

   call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 
   mvp = mvp + 1
! should b be vtemp because I am solving the error solution???
    
    do i = 1,nvhalf
!       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)

! NOTE~ need to add sigma(1)*xshift(...,1) to end of this even though sigma1=0

        r(:,i,:,:,:) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(1)*xshift(:,i,:,:,:,1)
!      v(:,i,:,:,:,:,1) = r(:,i,:,:,:,:)

!      do is=1,nshifts
!         xt(:,i,:,:,:,is) = x(:,i,:,:,:,:)
!      enddo ! is
    enddo ! i
 
! Copy the initial resiudal into the initial residual for each shift

! r and rshift(...,1) are the source vector at this point

!   do is=2,nshifts
       is =1
       rshift(:,:,:,:,:,is) = r(:,:,:,:,:)
!      rshift(:,:,:,:,:,is) = rshift(:,:,:,:,:,1)
!   enddo ! is
 
! Need to zero out the first shift after creating the first res.
! so that the solution from the projection section is not initiated
! with something other than zeros.

! ~NOTE this should not be done here, but not effetcing hopefully

    xshift(:,:,:,:,:,1) = 0.0_KR
   !if (myid == 0) then
   ! print *, "xshift2  in proj="
!    call checkNonZero(xshift(:,:,:,:,:,2),ntotal)
   !endif

!....Need logic here to determine the cycle in which it leaves...

    k = kDR

    cycledo: do

     if (rn/rinit <= resmax .or. icycle > kcyclim ) exit cycledo

     !if (myid==0) then
    !!   print *, "At start of loop PROJ rn/rinit = ", rn, rinit, rn/rinit
    ! endif

! DEAN~ By uncommenting next line - take out projection 
!     if (icycle -1 ==-1 )then

! The next if statment allows the projection step to occur.
     if (icycle-1 == ((icycle-1)/ifreq)*ifreq) then
 
      do i=1,k+1
        call vecdot(vtemp(:,:,:,:,:,i),rshift(:,:,:,:,:,1),beta,MRT2)
         c(1,i) = beta(1)
         c(2,i) = beta(2)
         c2(1,i) = c(1,i)
         c2(2,i) = c(2,i)
      enddo ! i

      do ii=1,k+1
        do jj=1,k 
          hc2(1,ii,jj) = hcnew(1,ii,jj)
          hc2(2,ii,jj) = hcnew(2,ii,jj)
        enddo ! jj
      enddo ! ii
 
      do jj =1,k
       hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(1)
      enddo ! jj

      do jj = 1,k
         do i = jj+1,k+1
            amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
            con2 = 1.0_KR2/amags
            tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
            tv(2) = 0.0_KR2
            gca(1,i,jj) = sqrt(amags)/tv(1)
            gca(2,i,jj) = 0.0_KR2
            gsa(1,i,jj) = gca(1,i,jj)*con2 &
                          *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
            gsa(2,i,jj) = gca(1,i,jj)*con2 &
                        *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
            do j = jj,k
               tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
                                             + gsa(2,i,jj)*hc2(2,i,j)
               tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
                                             - gsa(2,i,jj)*hc2(1,i,j)
               tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
                                            + gsa(2,i,jj)*hc2(2,jj,j)
               tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
                                            - gsa(2,i,jj)*hc2(1,jj,j)
               hc2(:,jj,j) = tv1(:)
               hc2(:,i,j) = tv2(:)
            enddo ! j
            tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
            tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
            tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
            tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
            c(:,jj) = tv1(:)
            c(:,i) = tv2(:)
         enddo ! i
      enddo ! jj

    ! Solve linear equation
! ~ NOTE check dimension of ss...does it need to be zeroed out?

      do i = 1,k
         ss(:,i) = c(:,i)
      enddo ! i

      con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
      const = con2*(ss(1,k)*hc2(1,k,k)+ss(2,k)*hc2(2,k,k))
      ss(2,k) = con2*(ss(2,k)*hc2(1,k,k)-ss(1,k)*hc2(2,k,k))
      ss(1,k) = const

      if (k/=1) then
         do i = 1,k-1
            ir = k - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc2(1,jj,ir) + ss(2,ir)*hc2(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc2(2,jj,ir) - ss(2,ir)*hc2(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc2(1,irm1,irm1)+ss(2,irm1)*hc2(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc2(1,irm1,irm1)-ss(1,irm1)*hc2(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (k/=1)

    ! ... Define new variable d to assist in shifting masses
    ! d is the "short" solution vector to small problem

      do jj=1,k
        d(1,jj,1) = ss(1,jj)
        d(2,jj,1) = ss(2,jj)
      enddo ! jj

! Put this in to take out of algorithum...temporarily. ONCE
!   IT WORKS NEED TO TAKE OUT THE is LOOP!!!

!     is=1
!     if(is==0) then
!
!     do is=2,nshifts
!
!        do ii=1,k+1
!           do jj=1,k
!              hc2(1,ii,jj) = hcnew(1,ii,jj)
!              hc2(2,ii,jj) = hcnew(2,ii,jj)
!           enddo ! jj
!        enddo ! ii
!
!        do ii=1,k+1
!          st(:,ii) = 0.0_KR
!        enddo ! ii
!
! NOTE ~ if first shift is not zero we need hc2 to be shifted...
!
!        do ii = 1,k+1
!           do jj = 1,k
!             st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
!             st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
!           enddo ! jj
!           srv(1,ii,1) = c2(1,ii) - st(1,ii)
!           srv(2,ii,1) = c2(2,ii) - st(2,ii)
!        enddo ! ii
!
!        do jj=1,k
!          hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(is)
!        enddo ! jj
!
! NOTE ~ d here is a "work" vector give differenet temp name..
!        NOT the short solution vector!
!
!        do ii=1,k+1
!           d(1,ii,is) = cmult(1,is)*st(1,ii) - cmult(2,is)*st(2,ii)
!           d(2,ii,is) = cmult(1,is)*st(2,ii) + cmult(2,is)*st(1,ii)
!        enddo ! ii
!
!       !if (myid == 0) then
!       ! print *, "is, ourrhshereis =", is, d(:,1:2,is)
!       ! print *, "hc2 before rotation =", hc2(:,1:3,1:2)
!       ! print *, " k before =", k
!       !endif ! myid
!
!        do jj = 1,k
!           do i = jj+1,k
!              amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
!              con2 = 1.0_KR2/amags
!              tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
!              tv(2) = 0.0_KR2
!              gca(1,i,jj) = sqrt(amags)/tv(1)
!              gca(2,i,jj) = 0.0_KR2
!              gsa(1,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
!              gsa(2,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
!              do j = jj,k
!                 tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
!                                                  + gsa(2,i,jj)*hc2(2,i,j)
!                 tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
!                                                  - gsa(2,i,jj)*hc2(1,i,j)
!                 tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
!                                                 + gsa(2,i,jj)*hc2(2,jj,j)
!                 tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
!                                                 - gsa(2,i,jj)*hc2(1,jj,j)
!                 hc2(:,jj,j) = tv1(:)
!                 hc2(:,i,j) = tv2(:)
!              enddo ! j
!              tv1(1) = gca(1,i,jj)*d(1,jj,is) + gsa(1,i,jj)*d(1,i,is) + gsa(2,i,jj)*d(2,i,is)
!              tv1(2) = gca(1,i,jj)*d(2,jj,is) + gsa(1,i,jj)*d(2,i,is) - gsa(2,i,jj)*d(1,i,is)
!              tv2(1) = gca(1,i,jj)*d(1,i,is) - gsa(1,i,jj)*d(1,jj,is) + gsa(2,i,jj)*d(2,jj,is)
!              tv2(2) = gca(1,i,jj)*d(2,i,is) - gsa(1,i,jj)*d(2,jj,is) - gsa(2,i,jj)*d(1,jj,is)
!!             d(:,jj,is) = tv1(:)
!              d(:,i,is) = tv2(:)
!           enddo ! i
!        enddo ! jj
!
! NOTE~ why are we zeroing out lower triangluar part of hc2
!       when these values should be allready zero.
!
!        do jj=1,k-1
!           do ii=jj+1,k
!              hc2(1,ii,jj) = 0.0_KR2
!              hc2(2,ii,jj) = 0.0_KR2
!           enddo ! ii
!        enddo ! jj
!
!        do ii=1,k
!           hc2(1,k+1,ii) = 0.0_KR2
!           hc2(2,k+1,ii) = 0.0_KR2
!        enddo ! ii
!
!      ! if (myid == 0) then
!      !  print *, "is, d in between =", is, d(:,1:2,is)
!      !  print *, "k in da' middle =", k
!      ! endif ! myid 
! NOTE~ back solve these linear equations
!
!     con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
!     const = con2*(d(1,k,is)*hc2(1,k,k)+d(2,k,is)*hc2(2,k,k))
!     d(2,k,is) = con2*(d(2,k,is)*hc2(1,k,k)-d(1,k,is)*hc2(2,k,k))
!     d(1,k,is) = const
!
!     if (k/=1) then
!        do i = 1,k-1
!           ir = k - i + 1
!           irm1 = ir - 1
!           do jj = 1,irm1
!              const = d(1,jj,is) - d(1,ir,is)*hc2(1,jj,ir) + d(2,ir,is)*hc2(2,jj,ir)
!              d(2,jj,is) = d(2,jj,is) - d(1,ir,is)*hc2(2,jj,ir) - d(2,ir,is)*hc2(1,jj,ir)
!              d(1,jj,is) = const
!           enddo ! jj
!           con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
!           const = con2*(d(1,irm1,is)*hc2(1,irm1,irm1)+d(2,irm1,is)*hc2(2,irm1,irm1))
!           d(2,irm1,is) = con2*(d(2,irm1,is)*hc2(1,irm1,irm1)-d(1,irm1,is)*hc2(2,irm1,irm1))
!           d(1,irm1,is) = const
!        enddo ! i
!     endif ! (k/=1)
!
!    !if (myid == 0) then
!    ! print *, "is, oursolnvector =" , is, d(:,1:2,is)
!    ! print *, "k after =", k
!    !endif ! myid
!
! NOTE ~ took out least squares solution upon DR Morgan request..
!
!        call real2complex_mat(hc2, nmaxGMRES+1, nmaxGMRES, zhc2)
!
!        call zgels('N', k+1, k, 1, zhc2, nmaxGMRES+1, zd(1,is), nmaxGMRES, zwork, lzwork, info)
!
!        call complex2real_mat(zd,nmaxGMRES,nshifts,d)
!        call complex2real_mat(zhc2, nmaxGMRES+1,nmaxGMRES, hc2)
!        
!     enddo ! is
!     
!     endif ! (is==0)
     
! Form the approximate new solution x = xt + xb.
! First zero xb then xb = V*d(1,2,:,is), and x =xt +xb
! ... Need to keep track of each solution for each shift.
!     Loop over shifts while creating the soln vector.
 
! DEAN ~ ONCE THIS ALGORITHUM WORKS, NEED TO TAKE OUT THE
!        is  LOOP AND HARD CODE IN is==1

 !          do is = 1,nshifts
               is = 1
               do i=1,k
                  sr(i) = d(1,i,is)
                  si(i) = d(2,i,is)
               enddo ! i
               do jj=1,nvhalf
                  xb(:,jj,:,:,:) = 0.0_KR2
               enddo ! jj
               do jj = 1,k
                  do icri = 1,5,2 ! 6=nri*nc
                     do i = 1,nvhalf
                        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                           + sr(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           - si(jj)*vtemp(icri+1,i,:,:,:,jj)
                        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                           + si(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           + sr(jj)*vtemp(icri+1,i,:,:,:,jj)
                     enddo ! i
                  enddo ! icri
               enddo ! jj
               do i = 1,nvhalf
                  xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
               enddo ! i

! This call to Hdbletm is purely for informational purpose.

               call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag, &
                               bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
               
               do ii = 1,nvhalf
                  rshift(:,ii,:,:,:,is) = beg(:,ii,:,:,:) - h(:,ii,:,:,:) &
                                                +sigma(is)*xshift(:,ii,:,:,:,is)
               enddo ! ii
!           enddo ! is

!DD note ~ equivqlent to exiting the 505 loop.

   endif ! (icycle-1 == ((icycle-1)/ifreq)*ifreq) 

  !if (myid == 0) then
  ! print *, "cycles of gmres"
  !endif

! Normalize rshift(:,:,:,:,:,1) and put into the first coln. of V, the 
! orthonormal matrix whose colns span the Krylov subspace.

! NOTE ~ may not need all the res norms for shifts >=2

   betashift = 0.0_KR

!    do is =1,nshifts
     is=1

       call vecdot(rshift(:,:,:,:,:,1),rshift(:,:,:,:,:,1), beta,MRT2)
       betashift(1,is) = beta(1)
       betashift(2,is) = beta(2)


!    enddo ! is

   const = 1.0_KR/sqrt(betashift(1,1))


! With the normalized residual form the first coln of V.
   do jj=1,nvhalf
     v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
   enddo ! jj  

   ztemp1 = DCMPLX(betashift(1,1),betashift(2,1))
   ztemp1 = sqrt(ztemp1)
   c(1,1) = REAL(ztemp1)
   c(2,1) = AIMAG(ztemp1)
   c2(1,1) = c(1,1)
   c2(2,1) = c(2,1)
   
! Perform GMRES between projections...

     do j = 1,nDR
       jp1 = j + 1
       itercount = itercount + 1
 
! LEFT OFF HERE!
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
 
       call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT)

       mvp = mvp + 1

       do i = 1,j
 
          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif
 
          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
 
          do icri = 1,5,2 ! 6=nri*nc
             do jj = 1,nvhalf
                v(icri  ,jj,:,:,:,jp1) = v(icri  ,jj,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,jj,:,:,:,i) &
                                        + beta(2)*v(icri+1,jj,:,:,:,i)
                v(icri+1,jj,:,:,:,jp1) = v(icri+1,jj,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,jj,:,:,:,i) &
                                        - beta(1)*v(icri+1,jj,:,:,:,i)
             enddo ! jj 
          enddo ! icri
       enddo ! i
 
 
       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)

       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
 
       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)

       const = 1.0_KR2/sqrt(beta(1))
       do jj=1,nvhalf
          v(:,jj,:,:,:,jp1) = const*v(:,jj,:,:,:,jp1)
       enddo

       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)
 
! DD note ~ I need to find a way of doing Givens rotations
! in the pseudocomplex routines.....

       if (j /= 1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (j /= 1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

      con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
      const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
      ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
      ss(1,j) = const

      if (j/=1) then
         do i = 1,j-1
            ir = j - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (j/=1)

      ! ... Define new variable d to assist in shifting masses
      ! d is the "short" solution vector to small problem

      do i=1,j
         d(1,i,1) = ss(1,i)
         d(2,i,1) = ss(2,i)
      enddo ! i

      do i=1,jp1
        st(:,i) = 0.0_KR2
      enddo ! i
  
      do ii = 1,jp1
         do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
         enddo ! jj
         srv(1,ii,1) = c2(1,ii) - st(1,ii)
         srv(2,ii,1) = c2(2,ii) - st(2,ii)
      enddo ! ii

      beta(1) = 0.0_KR

      do jj = 1,j+1
         beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
      enddo ! jj

      ! ... form gdr

      gdr(mvp,1) = sqrt(beta(1))

       ! ... To keep subspaces parallel for different shifts vector
       !     needs to be rotated. Use a QR factorization to do the
       !     ortognal rotation.
       !     BIG Shifting loop


!      if (is==0) then
!
!      do is = 2,nshifts
!
!         do ii=1,jp1
!            do jj=1,j
!               hcs(1,ii,jj) = hc2(1,ii,jj)
!               hcs(2,ii,jj) = hc2(2,ii,jj)
!               hcs2(1,ii,jj) = hc2(1,ii,jj)
!               hcs2(2,ii,jj) = hc2(2,ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         do jj=1,j
!            hcs(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is)
!            hcs2(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is) 
!         enddo ! jj
!
!         ! Copy the 12complex arrays into true complex arrays for use with
!         ! lapack routines
!
!         call real2complex_mat(hcs, nmaxGMRES+1, nmaxGMRES+1, zhcs)
!
!         call zgeqrf(j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! ... store R (upper triangular) in rr
!
!         do ii=1,jp1 
!            do jj=1,j
!               zrr(ii,jj) = zhcs(ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         call zungqr(j+1,j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! Copy the complex zhcs array back to hcs 
!
!         call complex2real_mat(zhcs, nmaxGMRES+1, nmaxGMRES+1, hcs)
!
!         ! ... hcs after this call is the qq part of the qr factorization
!
!         ! Now zero out crot (keeps shifts parrallel) and srvrot
!
!          do ii=1,jp1
!             zcrot(ii)=0.0_KR
!             zsrvrot(ii)=0.0_KR
!          enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,jp1
!
!               ztemp1 = DCMPLX(cmult(1,is), cmult(2,is))
!               ztemp2 = DCMPLX(c2(1,jj), c2(2,jj))
!               ztemp3 = DCMPLX(srv(1,jj,1), srv(2,jj,1))
!
!               zcrot(ii) = zcrot(ii) + ztemp1 * CONJG(zhcs(jj,ii)) * ztemp2
!               zsrvrot(ii) = zsrvrot(ii) + CONJG(zhcs(jj,ii)) * ztemp3
!
!            enddo ! jj
!         enddo ! ii
!
!         ! ... construct alpha
!
!           if ((myid == 0) .and. (is==3)) then
!             print *, "jp1, zcrot(jp1), zsrvrot(jp1) = ", jp1, zcrot(jp1), zsrvrot(jp1)
!           endif
!            
!         zalpha = zcrot(jp1)/zsrvrot(jp1)
!
!         alph(1, is) = REAL(zalpha)
!         alph(2, is) = AIMAG(zalpha)
!
!           if (myid==0)  then
!             print *,"is, zalpha, alph = ", is, zalpha, alph(:,is)
!           endif
!
!         do jj=1,j
!            ztemp1 = zalpha * zsrvrot(jj)
!            cmas(1,jj) = REAL(zcrot(jj)) - REAL(ztemp1)    ! alpha*srvrot(1,jj)
!            cmas(2,jj) = AIMAG(zcrot(jj)) - AIMAG(ztemp1)   ! alpha*srvrot(2,jj)
!         enddo ! jj
!
!           if ((myid==0) .and. (is==3))  then
!             print *,"cmas = ", cmas
!           endif
!
!         ! ... solve linear eqns problem d(1:j,is)=rr(1:j,1:j)\cmas
!
!         do ii=1,j
!            d(1,ii,is)=cmas(1,ii)
!            d(2,ii,is)=cmas(2,ii)
!         enddo ! ii
!
!         call real2complex_mat(d,nmaxGMRES,nshifts,zd)
!
!         zd(j,is) = zd(j,is)/zrr(j,j)
!
!         if (j /= 1) then 
!            do i=1,j-1
!               ir = j-i +1
!               irm1 = ir -1
!               call zaxpy(irm1,-zd(ir,is),zrr(1,ir),1,zd(1,is),1)
!               zd(irm1,is) = zd(irm1,is)/zrr(irm1,irm1)
!            enddo ! i
!         endif ! j/=1
!
!         call complex2real_mat(zd, nmaxGMRES, nshifts, d)
!
!         do ii=1,jp1
!            st(1,ii) = 0.0_KR 
!            st(2,ii) = 0.0_KR
!            srvis(1,ii) = 0.0_KR
!            srvis(2,ii) = 0.0_KR
!         enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,j
!               st(1,ii) = st(1,ii) + hcs2(1,ii,jj)*d(1,jj,is) &
!                                   - hcs2(2,ii,jj)*d(2,jj,is)
!               st(2,ii) = st(2,ii) + hcs2(1,ii,jj)*d(2,jj,is) &
!                                   + hcs2(2,ii,jj)*d(1,jj,is)
!            enddo ! jj
!            srvis(1,ii) = cmult(1,is)*c2(1,ii)-cmult(2,is)*c2(2,ii) &
!                        - st(1,ii)
!            srvis(2,ii) = cmult(1,is)*c2(2,ii)+cmult(2,is)*c2(1,ii) &
!                        - st(2,ii)
!         enddo ! ii
!
!         ! ... form the norm of srvis and put in gdr
!
!         beta(1) =0.0_KR
!
!         do jj=1,j+1
!            beta(1) = beta(1) + srvis(1,jj)*srvis(1,jj) &
!                              + srvis(2,jj)*srvis(2,jj)
!         enddo ! jj
!
!         gdr(mvp,is) = sqrt(beta(1))
!
!         if(myid==0.and.is==2) then 
!           print *, "mvp, gdr(mvp,2)=", mvp, gdr(mvp,2)
!         endif ! myid
!
!      enddo ! BIG is loop
!
!      endif ! (is==0)

       if (j>=nDR) then 

        !if(myid==0) then
        !  print *,"j=nDR", j,nDR
        !endif ! myid

! DEAN~ NEED TO HARD CODE IN is==1

!         do is = 1,nshifts
          is=1
             ! cmult(1,is) = alph(1,is) 
             ! cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i

             do jj=1,nvhalf
               xb(:,jj,:,:,:) = 0.0_KR2
             enddo ! jj

             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj

             do i = 1,nvhalf
                ! HEY!
                xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
             enddo ! i

            !if (myid == 0) then
            ! print *, "idag in proj=", idag
            !endif

             call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

              mvp = mvp + 1
             do i=1,nvhalf
                rshift(:,i,:,:,:,is) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(is)*xshift(:,i,:,:,:,is)
             enddo ! i 
           

!       enddo ! is
       endif ! (j>=nDR)
  
       if (j >= nDR) then
          betashift = 0.0_KR2

! DEAN HARD CODE IN is=1

!         do is=1,nshifts
          is=1
             call vecdot(rshift(:,:,:,:,:,is), rshift(:,:,:,:,:,is), beta, MRT2)
          
! ~ NOTE ... this should be a sqrt of beta.....

             betashift(1,is) = sqrt(beta(1))
             betashift(2,is) = 0.0_KR2
             rnale(icycle,is) = sqrt(beta(1))
!            cmult(1,is) = alph(1,is)
!            cmult(2,is) = alph(2,is)
!         enddo ! is

! use betashift (1,1) to check convergence of residual and ultimately exit.

         rn = betashift(1,1)

          !BS  if(myid==0) then
           !  print *, "Print before write to LOG",rn
           !endif ! myid

          ! if (j >= nDR) then 
             if (myid==0) then
               open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                     form="formatted",status="old",position="append")
             !BS432    write(unit=8,fmt=*) "gmresdrproject-gdr",mvp,rn/rinit
!                  write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,gdr(mvp,1)/rinit!,gdr(mvp,2), gdr(mvp,3), gdr(mvp,4)
                 !print *, betashift(1,1), betashift(1,2), betashift(1,3)
                close(unit=8,status="keep")
             endif
          ! endif ! (j >= nDR)

      !   if (myid ==0) then
      !      open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
      !           form="formatted",status="old",position="append")
!     !      write(unit=8,fmt=*) itercount,betashift(1,1),betashift(1,2),betashift(1,3)
      !      write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
      !            gdr(mvp,1)/rinit!, " ", gdr(mvp,2), " ", gdr(mvp,3), " ", gdr(mvp,4)
      !      close(unit=8,status="keep")
      !   endif ! myid

          rshift = 0.0_KR2

!         if (myid == 0) then
!            print *,"Checking rshift, v = ", v
!         endif

          do ii=1,nDR+1
             do icri = 1,5,2 ! 6=nri*nc
                do jj=1,nvhalf
                   rshift(icri  ,jj,:,:,:,1) = rshift(icri  ,jj,:,:,:,1) + v(icri  ,jj,:,:,:,ii)*srv(1,ii,1) &
                                                                         - v(icri+1,jj,:,:,:,ii)*srv(2,ii,1)
                   rshift(icri+1,jj,:,:,:,1) = rshift(icri+1,jj,:,:,:,1) + v(icri  ,jj,:,:,:,ii)*srv(2,ii,1) &
                                                                         + v(icri+1,jj,:,:,:,ii)*srv(1,ii,1)
                enddo ! jj
             enddo ! icri
          enddo ! ii

!          rn=gdr(mvp,1)

        !if (myid==0) then
        !   print *, "At bottom of loop PROJ rn/rinit = ", rn, rinit, rn/rinit, "j,mvp=",j,mvp
        !endif
        !if (myid==0) then
        !   print *, "end of gmres "
        !endif

! Don't need this vecdot

         call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

! NOTE~ since I took sqrt of betashift above don't need to here,

         const = 1.0_KR/betashift(1,1)

         do jj=1,nvhalf
            v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
         enddo ! jj

         icycle = icycle + 1

      endif ! (j>=nDR)
 
    enddo ! j BIG J LOOP

 enddo cycledo

 !  if (isignal == 1) then 
 !    if (myid == 0) then
 !        print *, "vk+1 gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  else ! isignal
 !    if (myid == 0) then
 !        print *, "RHS gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  endif ! isiganl

  ! if(myid==0) then
  !   print *, "Leaveing PROJ"
  ! endif ! myid


 end subroutine gmresproject
!-----------------------------------------------------------------------------
subroutine mmgmresproject(rwdir,b,xshift,GMRES,resmax,itermin, &
                       itercount,u,GeeGooinv, &
                       iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                       lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)
! GMRES-PROJECT(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
!
! Gmresproject projects over the approximate eigenvectors found in the 
! deflation section of gmresdr (gmresdrshift). These rojected evectors
! are used in the basis to solve the following right-hand side with the
! usual extraction methods.  
!
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresproject.

! This is PROJ
 
    use shift

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    ! real(kind=KR2),   intent(in),    dimension(:,:,:)     :: hcnew
    ! real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: vtemp
!   real(kind=KR2),   intent(inout), dimension(:,:,:,:,:) :: beg
    integer(kind=KI), intent(in)                          :: isignal
    integer(kind=KI), intent(inout)                       :: mvp
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
 
    integer(kind=KI) :: icycle, i, j, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa, ifreq,temprhs !,nrhs
 
    ! real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)     :: xshift
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               con2, rv, amags
    real(kind=KR2)                                          :: rn, rinit, rnnn 
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info,rnt
    real(kind=KR2),   dimension(2)                          :: beta, alpha, &
                                                               tv, tv1 ,tv2
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(kcyclim,nshifts)            :: rnale
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhc2
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h ,tes
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,1)          ::tes2    
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    ! real(kind=KR2),   dimension(2,nmaxGMRES)                :: gc, gs
    real(kind=KR2),   dimension(2)                          :: gam

    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha

    integer(kind=KI)                                        ::  isite, icolorir, idirac, iblock, &
                                                                site, icolorr, irow, ishift

! This is still PROJ 

 
! Shift sigmamu to base mu above (mtmqcd(1,2))
 
    sigma = 0.0_KR

! DEAN ~ HEY! I need to take out the sigma(1) part because I am not 
!        shifting at all and the residuals should be r = b-Ax
 
! init all x to 0

      xshift = 0.0_KR
 
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)


    ldh = nmaxGMRES+ 1
    ldz = nmaxGMRES+ 1
    lzwork = nmaxGMRES+1

!   if (myid==0) then
!     print *, "nDR = ", nDR
!     print *, "kDR = ", kDR
!   endif ! myid

    icycle = 1
    ifreq = 1
    idag = 0
    alph = 0.0_KR
    zcrot = 0.0_KR
    zsrvrot = 0.0_KR
    cmas = 0.0_KR
    srv = 0.0_KR
    srvis = 0.0_KR
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcs = 0.0_KR
    hcht = 0.0_KR
    zrr = 0.0_KR
    c2 = 0.0_KR
    c = 0.0_KR
    v = 0.0_KR
    xb = 0.0_KR
    d = 0.0_KR



! Initialize cmult
     
!   do is = 1,nshifts
       is =1
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0.0_KR
!   enddo ! is
 
    ! do is=1,nshifts
    !   xshift(:,:,:,:,:,is) = x(:,:,:,:,:,:)
    ! enddo ! is
 
! Compute r=b-M*x and v=r/|r| and beta=|r|.

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta,MRT2)
   !if (myid == 0) then
   ! print *, "norm o b in project =", sqrt(beta(1))
   !endif

! do iblock =1,8
!   do ieo = 1,2
!     do idirac=1,4
!       do isite=1,nvhalf 
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1
!
! To print single rhs source vector use ..
!
!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, b(icolorir,isite,idirac,ieo,iblock), b(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! isite
!       enddo ! idirac 
!    enddo ! ieo
!  enddo ! iblock

    call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

! NOTE ~ I don't think that I need this ....

    do ii = 1,nvhalf
        rshift(:,ii,:,:,:,1) = b(:,ii,:,:,:) - h(:,ii,:,:,:) +sigma(1)*xshift(:,ii,:,:,:,1)
    enddo ! ii

! For error correction put error in one direction and use gmresproject to 
! correct and solve the later right hand sides.

    if (isignal == 1) then
!     do is=1,nshifts
       rshift(:,:,:,:,:,1) = beg(:,:,:,:,:)
!     enddo ! is
    else
      beg(:,:,:,:,:) = rshift(:,:,:,:,:,1)
     !beg(:,:,:,:,:) = b(:,:,:,:,:)
    endif
!   call checkNonZero(beg,nvhalf)
! Create rinit for logic passing used in gmresprojet    

   beta = 0.0_KR

   call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

   rinit = sqrt(beta(1))
   rn = rinit

   call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 
! should b be vtemp because I am solving the error solution???
    
    do i = 1,nvhalf
!       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)

! NOTE~ need to add sigma(1)*xshift(...,1) to end of this even though sigma1=0

        r(:,i,:,:,:) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(1)*xshift(:,i,:,:,:,1)
!      v(:,i,:,:,:,:,1) = r(:,i,:,:,:,:)

!      do is=1,nshifts
!         xt(:,i,:,:,:,is) = x(:,i,:,:,:,:)
!      enddo ! is
    enddo ! i
 
! Copy the initial resiudal into the initial residual for each shift

! r and rshift(...,1) are the source vector at this point

!   do is=2,nshifts
       is =1
       rshift(:,:,:,:,:,is) = r(:,:,:,:,:)
!      rshift(:,:,:,:,:,is) = rshift(:,:,:,:,:,1)
!   enddo ! is
 
! Need to zero out the first shift after creating the first res.
! so that the solution from the projection section is not initiated
! with something other than zeros.

! ~NOTE this should not be done here, but not effetcing hopefully

    xshift(:,:,:,:,:,1) = 0.0_KR
   !if (myid == 0) then
   ! print *, "xshift2  in proj="
!    call checkNonZero(xshift(:,:,:,:,:,2),ntotal)
   !endif

!....Need logic here to determine the cycle in which it leaves...

    k = kDR

    cycledo: do

     if (rn/rinit <= resmax .or. icycle > kcyclim ) exit cycledo

     !if (myid==0) then
    !!   print *, "At start of loop PROJ rn/rinit = ", rn, rinit, rn/rinit
    ! endif

! DEAN~ By uncommenting next line - take out projection 
!     if (icycle -1 ==-1 )then

! The next if statment allows the projection step to occur.
     if (icycle-1 == ((icycle-1)/ifreq)*ifreq) then
!---------------------MM------------------------------------------------------ 
   idag = 1
   call Hdbletm(tes2(:,:,:,:,:,1),u,GeeGooinv,rshift(:,:,:,:,:,1), &
                idag,coact,kappa,iflag,bc,vecbl, &
                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   do i= 1,nvhalf
    rshift(:,i,:,:,:,1)=tes2(:,i,:,:,:,1)
   enddo!i
!------------------------------------------------------------------------------

      do i=1,k+1
        call vecdot(vtemp(:,:,:,:,:,i),rshift(:,:,:,:,:,1),beta,MRT2)
         c(1,i) = beta(1)
         c(2,i) = beta(2)
         c2(1,i) = c(1,i)
         c2(2,i) = c(2,i)
      enddo ! i

      do ii=1,k+1
        do jj=1,k 
          hc2(1,ii,jj) = hcnew(1,ii,jj)
          hc2(2,ii,jj) = hcnew(2,ii,jj)
        enddo ! jj
      enddo ! ii
 
      do jj =1,k
       hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(1)
      enddo ! jj

      do jj = 1,k
         do i = jj+1,k+1
            amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
            con2 = 1.0_KR2/amags
            tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
            tv(2) = 0.0_KR2
            gca(1,i,jj) = sqrt(amags)/tv(1)
            gca(2,i,jj) = 0.0_KR2
            gsa(1,i,jj) = gca(1,i,jj)*con2 &
                          *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
            gsa(2,i,jj) = gca(1,i,jj)*con2 &
                        *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
            do j = jj,k
               tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
                                             + gsa(2,i,jj)*hc2(2,i,j)
               tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
                                             - gsa(2,i,jj)*hc2(1,i,j)
               tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
                                            + gsa(2,i,jj)*hc2(2,jj,j)
               tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
                                            - gsa(2,i,jj)*hc2(1,jj,j)
               hc2(:,jj,j) = tv1(:)
               hc2(:,i,j) = tv2(:)
            enddo ! j
            tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
            tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
            tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
            tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
            c(:,jj) = tv1(:)
            c(:,i) = tv2(:)
         enddo ! i
      enddo ! jj

    ! Solve linear equation
! ~ NOTE check dimension of ss...does it need to be zeroed out?

      do i = 1,k
         ss(:,i) = c(:,i)
      enddo ! i

      con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
      const = con2*(ss(1,k)*hc2(1,k,k)+ss(2,k)*hc2(2,k,k))
      ss(2,k) = con2*(ss(2,k)*hc2(1,k,k)-ss(1,k)*hc2(2,k,k))
      ss(1,k) = const

      if (k/=1) then
         do i = 1,k-1
            ir = k - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc2(1,jj,ir) + ss(2,ir)*hc2(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc2(2,jj,ir) - ss(2,ir)*hc2(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc2(1,irm1,irm1)+ss(2,irm1)*hc2(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc2(1,irm1,irm1)-ss(1,irm1)*hc2(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (k/=1)

    ! ... Define new variable d to assist in shifting masses
    ! d is the "short" solution vector to small problem

      do jj=1,k
        d(1,jj,1) = ss(1,jj)
        d(2,jj,1) = ss(2,jj)
      enddo ! jj

! Put this in to take out of algorithum...temporarily. ONCE
!   IT WORKS NEED TO TAKE OUT THE is LOOP!!! !          do is = 1,nshifts
               is = 1
               do i=1,k
                  sr(i) = d(1,i,is)
                  si(i) = d(2,i,is)
               enddo ! i
               do jj=1,nvhalf
                  xb(:,jj,:,:,:) = 0.0_KR2
               enddo ! jj
               do jj = 1,k
                  do icri = 1,5,2 ! 6=nri*nc
                     do i = 1,nvhalf
                        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                           + sr(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           - si(jj)*vtemp(icri+1,i,:,:,:,jj)
                        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                           + si(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           + sr(jj)*vtemp(icri+1,i,:,:,:,jj)
                     enddo ! i
                  enddo ! icri
               enddo ! jj
               do i = 1,nvhalf
                  xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
               enddo ! i

! This call to Hdbletm is purely for informational purpose.

               call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag, &
                               bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
               do ii = 1,nvhalf
                  rshift(:,ii,:,:,:,is) = beg(:,ii,:,:,:) - h(:,ii,:,:,:) &
                                                +sigma(is)*xshift(:,ii,:,:,:,is)
               enddo ! ii
!           enddo ! is

!DD note ~ equivqlent to exiting the 505 loop.

   endif ! (icycle-1 == ((icycle-1)/ifreq)*ifreq) 

  !if (myid == 0) then
  ! print *, "cycles of gmres"
  !endif

! Normalize rshift(:,:,:,:,:,1) and put into the first coln. of V, the 
! orthonormal matrix whose colns span the Krylov subspace.

! NOTE ~ may not need all the res norms for shifts >=2

   betashift = 0.0_KR

!    do is =1,nshifts
     is=1

       call vecdot(rshift(:,:,:,:,:,1),rshift(:,:,:,:,:,1), beta,MRT2)
       betashift(1,is) = beta(1)
       betashift(2,is) = beta(2)


!    enddo ! is

   const = 1.0_KR/sqrt(betashift(1,1))


! With the normalized residual form the first coln of V.
   do jj=1,nvhalf
     v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
   enddo ! jj  

   ztemp1 = DCMPLX(betashift(1,1),betashift(2,1))
   ztemp1 = sqrt(ztemp1)
   c(1,1) = REAL(ztemp1)
   c(2,1) = AIMAG(ztemp1)
   c2(1,1) = c(1,1)
   c2(2,1) = c(2,1)
   
! Perform GMRES between projections...

     do j = 1,nDR
       jp1 = j + 1
       itercount = itercount + 1
 
! LEFT OFF HERE!
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
 
       idag = 0
       call Hdbletm(tes2(:,:,:,:,:,1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT)

       idag = 1
       call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,tes2(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT)

       mvp = mvp + 2

       do i = 1,j
 
          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif
 
          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
 
          do icri = 1,5,2 ! 6=nri*nc
             do jj = 1,nvhalf
                v(icri  ,jj,:,:,:,jp1) = v(icri  ,jj,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,jj,:,:,:,i) &
                                        + beta(2)*v(icri+1,jj,:,:,:,i)
                v(icri+1,jj,:,:,:,jp1) = v(icri+1,jj,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,jj,:,:,:,i) &
                                        - beta(1)*v(icri+1,jj,:,:,:,i)
             enddo ! jj 
          enddo ! icri
       enddo ! i
 
 
       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)

       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
 
       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)

       const = 1.0_KR2/sqrt(beta(1))
       do jj=1,nvhalf
          v(:,jj,:,:,:,jp1) = const*v(:,jj,:,:,:,jp1)
       enddo

       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)
 
! DD note ~ I need to find a way of doing Givens rotations
! in the pseudocomplex routines.....

       if (j /= 1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (j /= 1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

      con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
      const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
      ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
      ss(1,j) = const

      if (j/=1) then
         do i = 1,j-1
            ir = j - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (j/=1)

      ! ... Define new variable d to assist in shifting masses
      ! d is the "short" solution vector to small problem

      do i=1,j
         d(1,i,1) = ss(1,i)
         d(2,i,1) = ss(2,i)
      enddo ! i

      do i=1,jp1
        st(:,i) = 0.0_KR2
      enddo ! i
  
      do ii = 1,jp1
         do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
         enddo ! jj
         srv(1,ii,1) = c2(1,ii) - st(1,ii)
         srv(2,ii,1) = c2(2,ii) - st(2,ii)
      enddo ! ii

      beta(1) = 0.0_KR

      do jj = 1,j+1
         beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
      enddo ! jj

      ! ... form gdr

      gdr(mvp,1) = sqrt(beta(1))
!       

       if (j>=nDR) then 

        !if(myid==0) then
        !  print *,"j=nDR", j,nDR
        !endif ! myid

! DEAN~ NEED TO HARD CODE IN is==1

!         do is = 1,nshifts
          is=1
             ! cmult(1,is) = alph(1,is) 
             ! cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i

             do jj=1,nvhalf
               xb(:,jj,:,:,:) = 0.0_KR2
             enddo ! jj

             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj

             do i = 1,nvhalf
                ! HEY!
                xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
             enddo ! i

            !if (myid == 0) then
            ! print *, "idag in proj=", idag
            !endif


             idag=0
             call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

             do i=1,nvhalf
                rshift(:,i,:,:,:,is) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(is)*xshift(:,i,:,:,:,is)
             enddo ! i 

              idag=1
              call Hdbletm(tes2(:,:,:,:,:,1),u,GeeGooinv,&
                    rshift(:,:,:,:,:,is),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
              idag=0
          do i=1,nvhalf
           rshift(:,i,:,:,:,1)=tes2(:,i,:,:,:,1)
          enddo!i

!       enddo ! is
       endif ! (j>=nDR)
  
       if (j >= nDR) then
          betashift = 0.0_KR2

! DEAN HARD CODE IN is=1

!         do is=1,nshifts
          is=1
             call vecdot(rshift(:,:,:,:,:,is), rshift(:,:,:,:,:,is), beta, MRT2)
          
! ~ NOTE ... this should be a sqrt of beta.....

             betashift(1,is) = sqrt(beta(1))
             betashift(2,is) = 0.0_KR2
             rnale(icycle,is) = sqrt(beta(1))
!            cmult(1,is) = alph(1,is)
!            cmult(2,is) = alph(2,is)
!         enddo ! is

! use betashift (1,1) to check convergence of residual and ultimately exit.

         rn = betashift(1,1)

           !BSif(myid==0) then
            ! print *, "Print before write to LOG",rn
           !endif ! myid

          ! if (j >= nDR) then 
             if (myid==0) then
               open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                     form="formatted",status="old",position="append")
               !BS432  write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,rn/rinit
               !BS432   write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,gdr(mvp,1)/rinit!,gdr(mvp,2), gdr(mvp,3), gdr(mvp,4)
                 !print *, betashift(1,1), betashift(1,2), betashift(1,3)
                close(unit=8,status="keep")
             endif
          ! endif ! (j >= nDR)

      !   if (myid ==0) then
      !      open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
      !           form="formatted",status="old",position="append")
!     !      write(unit=8,fmt=*) itercount,betashift(1,1),betashift(1,2),betashift(1,3)
      !      write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
      !            gdr(mvp,1)/rinit!, " ", gdr(mvp,2), " ", gdr(mvp,3), " ", gdr(mvp,4)
      !      close(unit=8,status="keep")
      !   endif ! myid

          rshift = 0.0_KR2

!         if (myid == 0) then
!            print *,"Checking rshift, v = ", v
!         endif

          do ii=1,nDR+1
             do icri = 1,5,2 ! 6=nri*nc
                do jj=1,nvhalf
                   rshift(icri  ,jj,:,:,:,1) = rshift(icri  ,jj,:,:,:,1) + v(icri  ,jj,:,:,:,ii)*srv(1,ii,1) &
                                                                         - v(icri+1,jj,:,:,:,ii)*srv(2,ii,1)
                   rshift(icri+1,jj,:,:,:,1) = rshift(icri+1,jj,:,:,:,1) + v(icri  ,jj,:,:,:,ii)*srv(2,ii,1) &
                                                                         + v(icri+1,jj,:,:,:,ii)*srv(1,ii,1)
                enddo ! jj
             enddo ! icri
          enddo ! ii

!          rn=gdr(mvp,1)

        !if (myid==0) then
        !   print *, "At bottom of loop PROJ rn/rinit = ", rn, rinit, rn/rinit, "j,mvp=",j,mvp
        !endif
        !if (myid==0) then
        !   print *, "end of gmres "
        !endif

! Don't need this vecdot

         call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

! NOTE~ since I took sqrt of betashift above don't need to here,

         const = 1.0_KR/betashift(1,1)

         do jj=1,nvhalf
            v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
         enddo ! jj

         icycle = icycle + 1

      endif ! (j>=nDR)
 
    enddo ! j BIG J LOOP

 enddo cycledo

 !  if (isignal == 1) then 
 !    if (myid == 0) then
 !        print *, "vk+1 gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  else ! isignal
 !    if (myid == 0) then
 !        print *, "RHS gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  endif ! isiganl

  ! if(myid==0) then
  !   print *, "Leaveing PROJ"
  ! endif ! myid


 end subroutine mmgmresproject

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine ppmmgmresproject(rwdir,b,xshift,GMRES,resmax, &
                            itermin,itercount,u,GeeGooinv, &
                       iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                           lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)
! GMRES-PROJECT(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
!
! Gmresproject projects over the approximate eigenvectors found in the 
! deflation section of gmresdr (gmresdrshift). These rojected evectors
! are used in the basis to solve the following right-hand side with the
! usual extraction methods.  
!
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresproject.

! This is PROJ
 
    use shift

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    ! real(kind=KR2),   intent(in),    dimension(:,:,:)     :: hcnew
    ! real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: vtemp
!   real(kind=KR2),   intent(inout), dimension(:,:,:,:,:) :: beg
    integer(kind=KI), intent(in)                          :: isignal
    integer(kind=KI), intent(inout)                       :: mvp
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
 
    integer(kind=KI) :: icycle, i, j, k,kk, jp1, jj, p, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa, ifreq,temprhs !,nrhs
 
    ! real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)     :: xshift
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               con2, rv, amags
    real(kind=KR2)                                          :: rn, rinit, rnnn 
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info,rnt
    real(kind=KR2),   dimension(2)                          :: beta, alpha, &
                                                               tv, tv1 ,tv2
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(kcyclim,nshifts)            :: rnale
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhc2
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h, tes
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,1)       :: tes2
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    ! real(kind=KR2),   dimension(2,nmaxGMRES)                :: gc, gs
    real(kind=KR2),   dimension(2)                          :: gam

    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha

    integer(kind=KI)                                        ::  isite, icolorir, idirac, iblock, &
                                                                site, icolorr, irow, ishift

! This is still PROJ 

 
! Shift sigmamu to base mu above (mtmqcd(1,2))
    p = 6!order of the polynomial
    sigma = 0.0_KR
    y = 0.0_KR2!put y=0
    try = 0.0_KR2!put try=0
!    vprime = 0.0_KR2
! DEAN ~ HEY! I need to take out the sigma(1) part because I am not 
!        shifting at all and the residuals should be r = b-Ax
 
! init all x to 0

      xshift = 0.0_KR
 
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)


    ldh = nmaxGMRES+ 1
    ldz = nmaxGMRES+ 1
    lzwork = nmaxGMRES+1

!   if (myid==0) then
!     print *, "nDR = ", nDR
!     print *, "kDR = ", kDR
!   endif ! myid

    icycle = 1
    ifreq = 1
    idag = 0
    alph = 0.0_KR
    zcrot = 0.0_KR
    zsrvrot = 0.0_KR
    cmas = 0.0_KR
    srv = 0.0_KR
    srvis = 0.0_KR
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcs = 0.0_KR
    hcht = 0.0_KR
    zrr = 0.0_KR
    c2 = 0.0_KR
    c = 0.0_KR
    v = 0.0_KR
    xb = 0.0_KR
    d = 0.0_KR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!poly parameters!!!!!!!!!!!!!!!!!!!!!
!    lsmat = 0.0_KR2
!    co = 0.0_KR2
!    cls = 0.0_KR2
!    ipiv2 = 0.0_KI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize cmult
     
!   do is = 1,nshifts
       is =1
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0.0_KR
!   enddo ! is
 
    ! do is=1,nshifts
    !   xshift(:,:,:,:,:,is) = x(:,:,:,:,:,:)
    ! enddo ! is
 
! Compute r=b-M*x and v=r/|r| and beta=|r|.

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta,MRT2)
   !if (myid == 0) then
   ! print *, "norm o b in project =", sqrt(beta(1))
   !endif

! do iblock =1,8
!   do ieo = 1,2
!     do idirac=1,4
!       do isite=1,nvhalf 
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1
!
! To print single rhs source vector use ..
!
!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, b(icolorir,isite,idirac,ieo,iblock), b(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! isite
!       enddo ! idirac 
!    enddo ! ieo
!  enddo ! iblock

    call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

! NOTE ~ I don't think that I need this ....

    do ii = 1,nvhalf
        rshift(:,ii,:,:,:,1) = b(:,ii,:,:,:) - h(:,ii,:,:,:) +sigma(1)*xshift(:,ii,:,:,:,1)
    enddo ! ii

! For error correction put error in one direction and use gmresproject to 
! correct and solve the later right hand sides.

    if (isignal == 1) then
!     do is=1,nshifts
       rshift(:,:,:,:,:,1) = beg(:,:,:,:,:)
!     enddo ! is
    else
      beg(:,:,:,:,:) = rshift(:,:,:,:,:,1)
     !beg(:,:,:,:,:) = b(:,:,:,:,:)
    endif
!   call checkNonZero(beg,nvhalf)
! Create rinit for logic passing used in gmresprojet    

   beta = 0.0_KR

   call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

   rinit = sqrt(beta(1))
   rn = rinit

   call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 
! should b be vtemp because I am solving the error solution???
    
    do i = 1,nvhalf
!       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)

! NOTE~ need to add sigma(1)*xshift(...,1) to end of this even though sigma1=0

        r(:,i,:,:,:) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(1)*xshift(:,i,:,:,:,1)
!      v(:,i,:,:,:,:,1) = r(:,i,:,:,:,:)

!      do is=1,nshifts
!         xt(:,i,:,:,:,is) = x(:,i,:,:,:,:)
!      enddo ! is
    enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!generate the poly!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do kk = 1,nvhalf
!     vprime(:,kk,:,:,:,1) = b(:,kk,:,:,:)
!    enddo !kk
!    do i = 1,p
!     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
!                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
!                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!    enddo !i
    
!    do i=2,p+1
!     do j=2,p+1
!      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
!      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
!     enddo!j
!    enddo!i
        
!   do i=2,p+1
!     call vecdot(vprime(:,:,:,:,:,i),b(:,:,:,:,:),beta,MRT2)
!     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
!   enddo!i
    
!    call linearsolver(p,1,lsmat,ipiv2,cls)
!    co(:,:) = cls(:,:,1)    
   if(myid==0) then
    do i=1,p
     print *, "i,result from the project(:,i)=", i, co(:,i)
    enddo!i  
   endif!myid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copy the initial resiudal into the initial residual for each shift

! r and rshift(...,1) are the source vector at this point

!   do is=2,nshifts
       is =1
       rshift(:,:,:,:,:,is) = r(:,:,:,:,:)


!      rshift(:,:,:,:,:,is) = rshift(:,:,:,:,:,1)
!   enddo ! is
 
! Need to zero out the first shift after creating the first res.
! so that the solution from the projection section is not initiated
! with something other than zeros.

! ~NOTE this should not be done here, but not effetcing hopefully

    xshift(:,:,:,:,:,1) = 0.0_KR
   !if (myid == 0) then
   ! print *, "xshift2  in proj="
!    call checkNonZero(xshift(:,:,:,:,:,2),ntotal)
   !endif

!....Need logic here to determine the cycle in which it leaves...

    k = kDR

    cycledo: do

     if (rn/rinit <= resmax .or. icycle > kcyclim ) exit cycledo

     !if (myid==0) then
    !!   print *, "At start of loop PROJ rn/rinit = ", rn, rinit, rn/rinit
    ! endif

! DEAN~ By uncommenting next line - take out projection 
!     if (icycle -1 ==-1 )then

! The next if statment allows the projection step to occur.
     if (icycle-1 == ((icycle-1)/ifreq)*ifreq) then

!---------------------MM------------------------------------------------------ 
   idag = 1
   call Hdbletm(tes2(:,:,:,:,:,1),u,GeeGooinv,rshift(:,:,:,:,:,1), &
                idag,coact,kappa,iflag,bc,vecbl, &
                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   do i= 1,nvhalf
    rshift(:,i,:,:,:,1)=tes2(:,i,:,:,:,1)
   enddo!i
!------------------------------------------------------------------------------


 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii = 1,nvhalf
          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
      enddo!ii

      do icri=1,5,2
       do kk=1,nvhalf
        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                           -co(2,1)*try(icri+1,kk,:,:,:,1)
        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                             +co(2,1)*try(icri,kk,:,:,:,1)
       enddo!kk
      enddo!icri

    do i=1,p-1 
     idag=0
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     idag=1
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
 
       do icri=1,5,2
        do kk=1,nvhalf
         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
        enddo!k
       enddo!icri
    enddo!i
      do kk=1,nvhalf
       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,k+1
        call vecdot(vtemp(:,:,:,:,:,i),rshift(:,:,:,:,:,1),beta,MRT2)
         c(1,i) = beta(1)
         c(2,i) = beta(2)
         c2(1,i) = c(1,i)
         c2(2,i) = c(2,i)
      enddo ! i

      do ii=1,k+1
        do jj=1,k 
          hc2(1,ii,jj) = hcnew(1,ii,jj)
          hc2(2,ii,jj) = hcnew(2,ii,jj)
        enddo ! jj
      enddo ! ii
 
      do jj =1,k
       hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(1)
      enddo ! jj

      do jj = 1,k
         do i = jj+1,k+1
            amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
            con2 = 1.0_KR2/amags
            tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
            tv(2) = 0.0_KR2
            gca(1,i,jj) = sqrt(amags)/tv(1)
            gca(2,i,jj) = 0.0_KR2
            gsa(1,i,jj) = gca(1,i,jj)*con2 &
                          *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
            gsa(2,i,jj) = gca(1,i,jj)*con2 &
                        *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
            do j = jj,k
               tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
                                             + gsa(2,i,jj)*hc2(2,i,j)
               tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
                                             - gsa(2,i,jj)*hc2(1,i,j)
               tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
                                            + gsa(2,i,jj)*hc2(2,jj,j)
               tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
                                            - gsa(2,i,jj)*hc2(1,jj,j)
               hc2(:,jj,j) = tv1(:)
               hc2(:,i,j) = tv2(:)
            enddo ! j
            tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
            tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
            tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
            tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
            c(:,jj) = tv1(:)
            c(:,i) = tv2(:)
         enddo ! i
      enddo ! jj

    ! Solve linear equation
! ~ NOTE check dimension of ss...does it need to be zeroed out?

      do i = 1,k
         ss(:,i) = c(:,i)
      enddo ! i

      con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
      const = con2*(ss(1,k)*hc2(1,k,k)+ss(2,k)*hc2(2,k,k))
      ss(2,k) = con2*(ss(2,k)*hc2(1,k,k)-ss(1,k)*hc2(2,k,k))
      ss(1,k) = const

      if (k/=1) then
         do i = 1,k-1
            ir = k - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc2(1,jj,ir) + ss(2,ir)*hc2(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc2(2,jj,ir) - ss(2,ir)*hc2(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc2(1,irm1,irm1)+ss(2,irm1)*hc2(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc2(1,irm1,irm1)-ss(1,irm1)*hc2(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (k/=1)

    ! ... Define new variable d to assist in shifting masses
    ! d is the "short" solution vector to small problem

      do jj=1,k
        d(1,jj,1) = ss(1,jj)
        d(2,jj,1) = ss(2,jj)
      enddo ! jj

! Put this in to take out of algorithum...temporarily. ONCE
!   IT WORKS NEED TO TAKE OUT THE is LOOP!!!

!     is=1
!     if(is==0) then
!
!     do is=2,nshifts
!
!        do ii=1,k+1
!           do jj=1,k
!              hc2(1,ii,jj) = hcnew(1,ii,jj)
!              hc2(2,ii,jj) = hcnew(2,ii,jj)
!           enddo ! jj
!        enddo ! ii
!
!        do ii=1,k+1
!          st(:,ii) = 0.0_KR
!        enddo ! ii
!
! NOTE ~ if first shift is not zero we need hc2 to be shifted...
!
!        do ii = 1,k+1
!           do jj = 1,k
!             st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
!             st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
!           enddo ! jj
!           srv(1,ii,1) = c2(1,ii) - st(1,ii)
!           srv(2,ii,1) = c2(2,ii) - st(2,ii)
!        enddo ! ii
!
!        do jj=1,k
!          hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(is)
!        enddo ! jj
!
! NOTE ~ d here is a "work" vector give differenet temp name..
!        NOT the short solution vector!
!
!        do ii=1,k+1
!           d(1,ii,is) = cmult(1,is)*st(1,ii) - cmult(2,is)*st(2,ii)
!           d(2,ii,is) = cmult(1,is)*st(2,ii) + cmult(2,is)*st(1,ii)
!        enddo ! ii
!
!       !if (myid == 0) then
!       ! print *, "is, ourrhshereis =", is, d(:,1:2,is)
!       ! print *, "hc2 before rotation =", hc2(:,1:3,1:2)
!       ! print *, " k before =", k
!       !endif ! myid
!
!        do jj = 1,k
!           do i = jj+1,k
!              amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
!              con2 = 1.0_KR2/amags
!              tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
!              tv(2) = 0.0_KR2
!              gca(1,i,jj) = sqrt(amags)/tv(1)
!              gca(2,i,jj) = 0.0_KR2
!              gsa(1,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
!              gsa(2,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
!              do j = jj,k
!                 tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
!                                                  + gsa(2,i,jj)*hc2(2,i,j)
!                 tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
!                                                  - gsa(2,i,jj)*hc2(1,i,j)
!                 tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
!                                                 + gsa(2,i,jj)*hc2(2,jj,j)
!                 tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
!                                                 - gsa(2,i,jj)*hc2(1,jj,j)
!                 hc2(:,jj,j) = tv1(:)
!                 hc2(:,i,j) = tv2(:)
!              enddo ! j
!              tv1(1) = gca(1,i,jj)*d(1,jj,is) + gsa(1,i,jj)*d(1,i,is) + gsa(2,i,jj)*d(2,i,is)
!              tv1(2) = gca(1,i,jj)*d(2,jj,is) + gsa(1,i,jj)*d(2,i,is) - gsa(2,i,jj)*d(1,i,is)
!              tv2(1) = gca(1,i,jj)*d(1,i,is) - gsa(1,i,jj)*d(1,jj,is) + gsa(2,i,jj)*d(2,jj,is)
!              tv2(2) = gca(1,i,jj)*d(2,i,is) - gsa(1,i,jj)*d(2,jj,is) - gsa(2,i,jj)*d(1,jj,is)
!!             d(:,jj,is) = tv1(:)
!              d(:,i,is) = tv2(:)
!           enddo ! i
!        enddo ! jj
!
! NOTE~ why are we zeroing out lower triangluar part of hc2
!       when these values should be allready zero.
!
!        do jj=1,k-1
!           do ii=jj+1,k
!              hc2(1,ii,jj) = 0.0_KR2
!              hc2(2,ii,jj) = 0.0_KR2
!           enddo ! ii
!        enddo ! jj
!
!        do ii=1,k
!           hc2(1,k+1,ii) = 0.0_KR2
!           hc2(2,k+1,ii) = 0.0_KR2
!        enddo ! ii
!
!      ! if (myid == 0) then
!      !  print *, "is, d in between =", is, d(:,1:2,is)
!      !  print *, "k in da' middle =", k
!      ! endif ! myid 
! NOTE~ back solve these linear equations
!
!     con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
!     const = con2*(d(1,k,is)*hc2(1,k,k)+d(2,k,is)*hc2(2,k,k))
!     d(2,k,is) = con2*(d(2,k,is)*hc2(1,k,k)-d(1,k,is)*hc2(2,k,k))
!     d(1,k,is) = const
!
!     if (k/=1) then
!        do i = 1,k-1
!           ir = k - i + 1
!           irm1 = ir - 1
!           do jj = 1,irm1
!              const = d(1,jj,is) - d(1,ir,is)*hc2(1,jj,ir) + d(2,ir,is)*hc2(2,jj,ir)
!              d(2,jj,is) = d(2,jj,is) - d(1,ir,is)*hc2(2,jj,ir) - d(2,ir,is)*hc2(1,jj,ir)
!              d(1,jj,is) = const
!           enddo ! jj
!           con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
!           const = con2*(d(1,irm1,is)*hc2(1,irm1,irm1)+d(2,irm1,is)*hc2(2,irm1,irm1))
!           d(2,irm1,is) = con2*(d(2,irm1,is)*hc2(1,irm1,irm1)-d(1,irm1,is)*hc2(2,irm1,irm1))
!           d(1,irm1,is) = const
!        enddo ! i
!     endif ! (k/=1)
!
!    !if (myid == 0) then
!    ! print *, "is, oursolnvector =" , is, d(:,1:2,is)
!    ! print *, "k after =", k
!    !endif ! myid
!
! NOTE ~ took out least squares solution upon DR Morgan request..
!
!        call real2complex_mat(hc2, nmaxGMRES+1, nmaxGMRES, zhc2)
!
!        call zgels('N', k+1, k, 1, zhc2, nmaxGMRES+1, zd(1,is), nmaxGMRES, zwork, lzwork, info)
!
!        call complex2real_mat(zd,nmaxGMRES,nshifts,d)
!        call complex2real_mat(zhc2, nmaxGMRES+1,nmaxGMRES, hc2)
!        
!     enddo ! is
!     
!     endif ! (is==0)
     
! Form the approximate new solution x = xt + xb.
! First zero xb then xb = V*d(1,2,:,is), and x =xt +xb
! ... Need to keep track of each solution for each shift.
!     Loop over shifts while creating the soln vector.
 
! DEAN ~ ONCE THIS ALGORITHUM WORKS, NEED TO TAKE OUT THE
!        is  LOOP AND HARD CODE IN is==1

 !          do is = 1,nshifts
               is = 1
               do i=1,k
                  sr(i) = d(1,i,is)
                  si(i) = d(2,i,is)
               enddo ! i
               do jj=1,nvhalf
                  xb(:,jj,:,:,:) = 0.0_KR2
               enddo ! jj
               do jj = 1,k
                  do icri = 1,5,2 ! 6=nri*nc
                     do i = 1,nvhalf
                        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                           + sr(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           - si(jj)*vtemp(icri+1,i,:,:,:,jj)
                        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                           + si(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           + sr(jj)*vtemp(icri+1,i,:,:,:,jj)
                     enddo ! i
                  enddo ! icri
               enddo ! jj
               do i = 1,nvhalf
                  xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
               enddo ! i

! This call to Hdbletm is purely for informational purpose.

               call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag, &
                               bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
               do ii = 1,nvhalf
                  rshift(:,ii,:,:,:,is) = beg(:,ii,:,:,:) - h(:,ii,:,:,:) &
                                                +sigma(is)*xshift(:,ii,:,:,:,is)
               enddo ! ii
!           enddo ! is

!DD note ~ equivqlent to exiting the 505 loop.

   endif ! (icycle-1 == ((icycle-1)/ifreq)*ifreq) 

  !if (myid == 0) then
  ! print *, "cycles of gmres"
  !endif

! Normalize rshift(:,:,:,:,:,1) and put into the first coln. of V, the 
! orthonormal matrix whose colns span the Krylov subspace.

! NOTE ~ may not need all the res norms for shifts >=2

 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii = 1,nvhalf
          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
      enddo!ii

      do icri=1,5,2
       do kk=1,nvhalf
        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                           -co(2,1)*try(icri+1,kk,:,:,:,1)
        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                             +co(2,1)*try(icri,kk,:,:,:,1)
       enddo!kk
      enddo!icri

   do i=1,p-1 
     idag=0
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     idag=1
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z

       do icri=1,5,2
        do kk=1,nvhalf
         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
        enddo!k
       enddo!icri
   enddo!i
      do kk=1,nvhalf
       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   betashift = 0.0_KR

!    do is =1,nshifts
     is=1

       call vecdot(rshift(:,:,:,:,:,1),rshift(:,:,:,:,:,1), beta,MRT2)
       betashift(1,is) = beta(1)
       betashift(2,is) = beta(2)


!    enddo ! is

   const = 1.0_KR/sqrt(betashift(1,1))


! With the normalized residual form the first coln of V.
   do jj=1,nvhalf
     v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
   enddo ! jj  

   ztemp1 = DCMPLX(betashift(1,1),betashift(2,1))
   ztemp1 = sqrt(ztemp1)
   c(1,1) = REAL(ztemp1)
   c(2,1) = AIMAG(ztemp1)
   c2(1,1) = c(1,1)
   c2(2,1) = c(2,1)
   
! Perform GMRES between projections...

     do j = 1,nDR
       jp1 = j + 1
       itercount = itercount + 1
 
! LEFT OFF HERE!
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!v_j=P(A)*v_j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    try = 0.0_KR2
    y = 0.0_KR2
    do i=1,nvhalf
     try(:,i,:,:,:,1) = v(:,i,:,:,:,j)
    enddo!i
    do icri=1,5,2
     do kk=1,nvhalf
      y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                         -co(2,1)*try(icri+1,kk,:,:,:,1)
      y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                           +co(2,1)*try(icri,kk,:,:,:,1)
     enddo!kk
    enddo!icri
    do i=1,p-1 
     call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     idag=1
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z

     do icri=1,5,2
      do kk=1,nvhalf
       y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                             +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
       y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i

    idag=0
    call Hdbletm(test(:,:,:,:,:,1),u,GeeGooinv,y(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)
    idag=1
    call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,test(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)

     mvp = mvp+2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
!


       do i = 1,j
 
          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif
 
          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
 
          do icri = 1,5,2 ! 6=nri*nc
             do jj = 1,nvhalf
                v(icri  ,jj,:,:,:,jp1) = v(icri  ,jj,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,jj,:,:,:,i) &
                                        + beta(2)*v(icri+1,jj,:,:,:,i)
                v(icri+1,jj,:,:,:,jp1) = v(icri+1,jj,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,jj,:,:,:,i) &
                                        - beta(1)*v(icri+1,jj,:,:,:,i)
             enddo ! jj 
          enddo ! icri
       enddo ! i
 
 
       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)

       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
 
       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)

       const = 1.0_KR2/sqrt(beta(1))
       do jj=1,nvhalf
          v(:,jj,:,:,:,jp1) = const*v(:,jj,:,:,:,jp1)
       enddo

       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)
 
! DD note ~ I need to find a way of doing Givens rotations
! in the pseudocomplex routines.....

       if (j /= 1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (j /= 1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

      con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
      const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
      ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
      ss(1,j) = const

      if (j/=1) then
         do i = 1,j-1
            ir = j - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (j/=1)

      ! ... Define new variable d to assist in shifting masses
      ! d is the "short" solution vector to small problem

      do i=1,j
         d(1,i,1) = ss(1,i)
         d(2,i,1) = ss(2,i)
      enddo ! i

      do i=1,jp1
        st(:,i) = 0.0_KR2
      enddo ! i
  
      do ii = 1,jp1
         do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
         enddo ! jj
         srv(1,ii,1) = c2(1,ii) - st(1,ii)
         srv(2,ii,1) = c2(2,ii) - st(2,ii)
      enddo ! ii

      beta(1) = 0.0_KR

      do jj = 1,j+1
         beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
      enddo ! jj

      ! ... form gdr

      gdr(mvp,1) = sqrt(beta(1))

       ! ... To keep subspaces parallel for different shifts vector
       !     needs to be rotated. Use a QR factorization to do the
       !     ortognal rotation.
       !     BIG Shifting loop


!      if (is==0) then
!
!      do is = 2,nshifts
!
!         do ii=1,jp1
!            do jj=1,j
!               hcs(1,ii,jj) = hc2(1,ii,jj)
!               hcs(2,ii,jj) = hc2(2,ii,jj)
!               hcs2(1,ii,jj) = hc2(1,ii,jj)
!               hcs2(2,ii,jj) = hc2(2,ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         do jj=1,j
!            hcs(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is)
!            hcs2(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is) 
!         enddo ! jj
!
!         ! Copy the 12complex arrays into true complex arrays for use with
!         ! lapack routines
!
!         call real2complex_mat(hcs, nmaxGMRES+1, nmaxGMRES+1, zhcs)
!
!         call zgeqrf(j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! ... store R (upper triangular) in rr
!
!         do ii=1,jp1 
!            do jj=1,j
!               zrr(ii,jj) = zhcs(ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         call zungqr(j+1,j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! Copy the complex zhcs array back to hcs 
!
!         call complex2real_mat(zhcs, nmaxGMRES+1, nmaxGMRES+1, hcs)
!
!         ! ... hcs after this call is the qq part of the qr factorization
!
!         ! Now zero out crot (keeps shifts parrallel) and srvrot
!
!          do ii=1,jp1
!             zcrot(ii)=0.0_KR
!             zsrvrot(ii)=0.0_KR
!          enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,jp1
!
!               ztemp1 = DCMPLX(cmult(1,is), cmult(2,is))
!               ztemp2 = DCMPLX(c2(1,jj), c2(2,jj))
!               ztemp3 = DCMPLX(srv(1,jj,1), srv(2,jj,1))
!
!               zcrot(ii) = zcrot(ii) + ztemp1 * CONJG(zhcs(jj,ii)) * ztemp2
!               zsrvrot(ii) = zsrvrot(ii) + CONJG(zhcs(jj,ii)) * ztemp3
!
!            enddo ! jj
!         enddo ! ii
!
!         ! ... construct alpha
!
!           if ((myid == 0) .and. (is==3)) then
!             print *, "jp1, zcrot(jp1), zsrvrot(jp1) = ", jp1, zcrot(jp1), zsrvrot(jp1)
!           endif
!            
!         zalpha = zcrot(jp1)/zsrvrot(jp1)
!
!         alph(1, is) = REAL(zalpha)
!         alph(2, is) = AIMAG(zalpha)
!
!           if (myid==0)  then
!             print *,"is, zalpha, alph = ", is, zalpha, alph(:,is)
!           endif
!
!         do jj=1,j
!            ztemp1 = zalpha * zsrvrot(jj)
!            cmas(1,jj) = REAL(zcrot(jj)) - REAL(ztemp1)    ! alpha*srvrot(1,jj)
!            cmas(2,jj) = AIMAG(zcrot(jj)) - AIMAG(ztemp1)   ! alpha*srvrot(2,jj)
!         enddo ! jj
!
!           if ((myid==0) .and. (is==3))  then
!             print *,"cmas = ", cmas
!           endif
!
!         ! ... solve linear eqns problem d(1:j,is)=rr(1:j,1:j)\cmas
!
!         do ii=1,j
!            d(1,ii,is)=cmas(1,ii)
!            d(2,ii,is)=cmas(2,ii)
!         enddo ! ii
!
!         call real2complex_mat(d,nmaxGMRES,nshifts,zd)
!
!         zd(j,is) = zd(j,is)/zrr(j,j)
!
!         if (j /= 1) then 
!            do i=1,j-1
!               ir = j-i +1
!               irm1 = ir -1
!               call zaxpy(irm1,-zd(ir,is),zrr(1,ir),1,zd(1,is),1)
!               zd(irm1,is) = zd(irm1,is)/zrr(irm1,irm1)
!            enddo ! i
!         endif ! j/=1
!
!         call complex2real_mat(zd, nmaxGMRES, nshifts, d)
!
!         do ii=1,jp1
!            st(1,ii) = 0.0_KR 
!            st(2,ii) = 0.0_KR
!            srvis(1,ii) = 0.0_KR
!            srvis(2,ii) = 0.0_KR
!         enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,j
!               st(1,ii) = st(1,ii) + hcs2(1,ii,jj)*d(1,jj,is) &
!                                   - hcs2(2,ii,jj)*d(2,jj,is)
!               st(2,ii) = st(2,ii) + hcs2(1,ii,jj)*d(2,jj,is) &
!                                   + hcs2(2,ii,jj)*d(1,jj,is)
!            enddo ! jj
!            srvis(1,ii) = cmult(1,is)*c2(1,ii)-cmult(2,is)*c2(2,ii) &
!                        - st(1,ii)
!            srvis(2,ii) = cmult(1,is)*c2(2,ii)+cmult(2,is)*c2(1,ii) &
!                        - st(2,ii)
!         enddo ! ii
!
!         ! ... form the norm of srvis and put in gdr
!
!         beta(1) =0.0_KR
!
!         do jj=1,j+1
!            beta(1) = beta(1) + srvis(1,jj)*srvis(1,jj) &
!                              + srvis(2,jj)*srvis(2,jj)
!         enddo ! jj
!
!         gdr(mvp,is) = sqrt(beta(1))
!
!         if(myid==0.and.is==2) then 
!           print *, "mvp, gdr(mvp,2)=", mvp, gdr(mvp,2)
!         endif ! myid
!
!      enddo ! BIG is loop
!
!      endif ! (is==0)

       if (j>=nDR) then 

        !if(myid==0) then
        !  print *,"j=nDR", j,nDR
        !endif ! myid

! DEAN~ NEED TO HARD CODE IN is==1

!         do is = 1,nshifts
          is=1
             ! cmult(1,is) = alph(1,is) 
             ! cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i

             do jj=1,nvhalf
               xb(:,jj,:,:,:) = 0.0_KR2
             enddo ! jj

             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj

             do i = 1,nvhalf
                ! HEY!
                xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
             enddo ! i

            !if (myid == 0) then
            ! print *, "idag in proj=", idag
            !endif

             call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,&
                         iflag,bc,vecbl,vecblinv,myid,nn,ldiv, &
                         nms,lvbc,ib,lbd,iblv,MRT)

             do i=1,nvhalf
                rshift(:,i,:,:,:,is) = beg(:,i,:,:,:) - h(:,i,:,:,:) &
                                 + sigma(is)*xshift(:,i,:,:,:,is)
             enddo ! i 
           
              call Hdbletm(tes2(:,:,:,:,:,1),u,GeeGooinv,&
                    rshift(:,:,:,:,:,is),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          do i=1,nvhalf
           rshift(:,i,:,:,:,1)=tes2(:,i,:,:,:,1)
          enddo!i

!       enddo ! is
       endif ! (j>=nDR)
  
       if (j >= nDR) then
          betashift = 0.0_KR2

! DEAN HARD CODE IN is=1
 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do ii = 1,nvhalf
!          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
!      enddo!ii

!      do icri=1,5,2
!       do kk=1,nvhalf
!        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
!                           -co(2,1)*try(icri+1,kk,:,:,:,1)
!        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
!                             +co(2,1)*try(icri,kk,:,:,:,1)
!       enddo!kk
!      enddo!icri

!      do i=1,p-1 
!       call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag, &
!                    coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
!                    lvbc,ib,lbd,iblv,MRT )  !z1=M*z
!       do icri=1,5,2
!        do kk=1,nvhalf
!         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
!                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
!                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
!         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
!                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
!                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
!        enddo!k
!       enddo!icri
!      enddo!i
!      do kk=1,nvhalf
!       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
!      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         do is=1,nshifts
          is=1
             call vecdot(rshift(:,:,:,:,:,is), rshift(:,:,:,:,:,is), beta, MRT2)
          
! ~ NOTE ... this should be a sqrt of beta.....

             betashift(1,is) = sqrt(beta(1))
             betashift(2,is) = 0.0_KR2
             rnale(icycle,is) = sqrt(beta(1))
!            cmult(1,is) = alph(1,is)
!            cmult(2,is) = alph(2,is)
!         enddo ! is

! use betashift (1,1) to check convergence of residual and ultimately exit.

         rn = betashift(1,1)

           !if(myid==0) then
           !  print *, "Print before write to LOG"
           !endif ! myid

          ! if (j >= nDR) then 
             if (myid==0) then
                open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                     form="formatted",status="old",position="append")
                !BS432  write(unit=8,fmt=*) "gmresdrproject-rn",itercount,rn/rinit
                !BS432  write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,gdr(mvp,1)/rinit!,gdr(mvp,2), gdr(mvp,3), gdr(mvp,4)
                 !print *, betashift(1,1), betashift(1,2), betashift(1,3)
                close(unit=8,status="keep")
             endif
          ! endif ! (j >= nDR)

      !   if (myid ==0) then
      !      open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
      !           form="formatted",status="old",position="append")
!     !      write(unit=8,fmt=*) itercount,betashift(1,1),betashift(1,2),betashift(1,3)
      !      write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
      !            gdr(mvp,1)/rinit!, " ", gdr(mvp,2), " ", gdr(mvp,3), " ", gdr(mvp,4)
      !      close(unit=8,status="keep")
      !   endif ! myid

          rshift = 0.0_KR2

!         if (myid == 0) then
!            print *,"Checking rshift, v = ", v
!         endif

          do ii=1,nDR+1
             do icri = 1,5,2 ! 6=nri*nc
                do jj=1,nvhalf
                 rshift(icri  ,jj,:,:,:,1) = rshift(icri  ,jj,:,:,:,1) & 
                                            + v(icri,jj,:,:,:,ii)*srv(1,ii,1) &
                                            - v(icri+1,jj,:,:,:,ii)*srv(2,ii,1)
                 rshift(icri+1,jj,:,:,:,1) = rshift(icri+1,jj,:,:,:,1) & 
                                             +v(icri,jj,:,:,:,ii)*srv(2,ii,1) &
                                             + v(icri+1,jj,:,:,:,ii)*srv(1,ii,1)
                enddo ! jj
             enddo ! icri
          enddo ! ii

!          rn=gdr(mvp,1)

        !if (myid==0) then
        !   print *, "At bottom of loop PROJ rn/rinit = ", rn, rinit, rn/rinit, "j,mvp=",j,mvp
        !endif
        !if (myid==0) then
        !   print *, "end of gmres "
        !endif

! Don't need this vecdot

         call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

! NOTE~ since I took sqrt of betashift above don't need to here,

         const = 1.0_KR/betashift(1,1)

         do jj=1,nvhalf
            v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
         enddo ! jj

         icycle = icycle + 1

      endif ! (j>=nDR)
 
    enddo ! j BIG J LOOP

 enddo cycledo

 !  if (isignal == 1) then 
 !    if (myid == 0) then
 !        print *, "vk+1 gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  else ! isignal
 !    if (myid == 0) then
 !        print *, "RHS gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  endif ! isiganl

  ! if(myid==0) then
  !   print *, "Leaveing PROJ"
  ! endif ! myid


 end subroutine ppmmgmresproject



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine ppgmresproject(rwdir,b,xshift,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)
! GMRES-PROJECT(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
!
! Gmresproject projects over the approximate eigenvectors found in the 
! deflation section of gmresdr (gmresdrshift). These rojected evectors
! are used in the basis to solve the following right-hand side with the
! usual extraction methods.  
!
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresproject.

! This is PROJ
 
    use shift

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    ! real(kind=KR2),   intent(in),    dimension(:,:,:)     :: hcnew
    ! real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: vtemp
!   real(kind=KR2),   intent(inout), dimension(:,:,:,:,:) :: beg
    integer(kind=KI), intent(in)                          :: isignal
    integer(kind=KI), intent(inout)                       :: mvp
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
 
    integer(kind=KI) :: icycle, i, j, k,kk, jp1, jj, p, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa, ifreq,temprhs !,nrhs
 
    ! real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)     :: xshift
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               con2, rv, amags
    real(kind=KR2)                                          :: rn, rinit, rnnn 
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info,rnt
    real(kind=KR2),   dimension(2)                          :: beta, alpha, &
                                                               tv, tv1 ,tv2
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(kcyclim,nshifts)            :: rnale
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhc2
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    ! real(kind=KR2),   dimension(2,nmaxGMRES)                :: gc, gs
    real(kind=KR2),   dimension(2)                          :: gam

    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha

    integer(kind=KI)                                        ::  isite, icolorir, idirac, iblock, &
                                                                site, icolorr, irow, ishift

! This is still PROJ 

 
! Shift sigmamu to base mu above (mtmqcd(1,2))
    p = 6!order of the polynomial
    sigma = 0.0_KR
    y = 0.0_KR2!put y=0
    try = 0.0_KR2!put try=0
!    vprime = 0.0_KR2
! DEAN ~ HEY! I need to take out the sigma(1) part because I am not 
!        shifting at all and the residuals should be r = b-Ax
 
! init all x to 0

      xshift = 0.0_KR
 
      mvp = 0
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)


    ldh = nmaxGMRES+ 1
    ldz = nmaxGMRES+ 1
    lzwork = nmaxGMRES+1

!   if (myid==0) then
!     print *, "nDR = ", nDR
!     print *, "kDR = ", kDR
!   endif ! myid

    icycle = 1
    ifreq = 1
    idag = 0
    alph = 0.0_KR
    zcrot = 0.0_KR
    zsrvrot = 0.0_KR
    cmas = 0.0_KR
    srv = 0.0_KR
    srvis = 0.0_KR
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcs = 0.0_KR
    hcht = 0.0_KR
    zrr = 0.0_KR
    c2 = 0.0_KR
    c = 0.0_KR
    v = 0.0_KR
    xb = 0.0_KR
    d = 0.0_KR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!poly parameters!!!!!!!!!!!!!!!!!!!!!
!    lsmat = 0.0_KR2
!    co = 0.0_KR2
!    cls = 0.0_KR2
!    ipiv2 = 0.0_KI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize cmult
     
!   do is = 1,nshifts
       is =1
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0.0_KR
!   enddo ! is
 
    ! do is=1,nshifts
    !   xshift(:,:,:,:,:,is) = x(:,:,:,:,:,:)
    ! enddo ! is
 
! Compute r=b-M*x and v=r/|r| and beta=|r|.

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta,MRT2)
   !if (myid == 0) then
   ! print *, "norm o b in project =", sqrt(beta(1))
   !endif

! do iblock =1,8
!   do ieo = 1,2
!     do idirac=1,4
!       do isite=1,nvhalf 
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1
!
! To print single rhs source vector use ..
!
!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, b(icolorir,isite,idirac,ieo,iblock), b(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! isite
!       enddo ! idirac 
!    enddo ! ieo
!  enddo ! iblock

    call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

    mvp = mvp + 1
! NOTE ~ I don't think that I need this ....

    do ii = 1,nvhalf
        rshift(:,ii,:,:,:,1) = b(:,ii,:,:,:) - h(:,ii,:,:,:) +sigma(1)*xshift(:,ii,:,:,:,1)
    enddo ! ii

! For error correction put error in one direction and use gmresproject to 
! correct and solve the later right hand sides.

    if (isignal == 1) then
!     do is=1,nshifts
       rshift(:,:,:,:,:,1) = beg(:,:,:,:,:)
!     enddo ! is
    else
      beg(:,:,:,:,:) = rshift(:,:,:,:,:,1)
     !beg(:,:,:,:,:) = b(:,:,:,:,:)
    endif
!   call checkNonZero(beg,nvhalf)
! Create rinit for logic passing used in gmresprojet    

   beta = 0.0_KR

   call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

   rinit = sqrt(beta(1))
   rn = rinit

   call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 
    mvp = mvp + 1
! should b be vtemp because I am solving the error solution???
    
    do i = 1,nvhalf
!       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)

! NOTE~ need to add sigma(1)*xshift(...,1) to end of this even though sigma1=0

        r(:,i,:,:,:) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(1)*xshift(:,i,:,:,:,1)
!      v(:,i,:,:,:,:,1) = r(:,i,:,:,:,:)

!      do is=1,nshifts
!         xt(:,i,:,:,:,is) = x(:,i,:,:,:,:)
!      enddo ! is
    enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!generate the poly!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do kk = 1,nvhalf
!     vprime(:,kk,:,:,:,1) = b(:,kk,:,:,:)
!    enddo !kk
!    do i = 1,p
!     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
!                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
!                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!    enddo !i
    
!    do i=2,p+1
!     do j=2,p+1
!      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
!      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
!     enddo!j
!    enddo!i
        
!   do i=2,p+1
!     call vecdot(vprime(:,:,:,:,:,i),b(:,:,:,:,:),beta,MRT2)
!     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
!   enddo!i
    
!    call linearsolver(p,1,lsmat,ipiv2,cls)
!    co(:,:) = cls(:,:,1)    
   if(myid==0) then
    do i=1,p
     print *, "i,result from the project(:,i)=", i, co(:,i)
    enddo!i  
   endif!myid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copy the initial resiudal into the initial residual for each shift

! r and rshift(...,1) are the source vector at this point

!   do is=2,nshifts
       is =1
       rshift(:,:,:,:,:,is) = r(:,:,:,:,:)


!      rshift(:,:,:,:,:,is) = rshift(:,:,:,:,:,1)
!   enddo ! is
 
! Need to zero out the first shift after creating the first res.
! so that the solution from the projection section is not initiated
! with something other than zeros.

! ~NOTE this should not be done here, but not effetcing hopefully

    xshift(:,:,:,:,:,1) = 0.0_KR
   !if (myid == 0) then
   ! print *, "xshift2  in proj="
!    call checkNonZero(xshift(:,:,:,:,:,2),ntotal)
   !endif

!....Need logic here to determine the cycle in which it leaves...

    k = kDR

    cycledo: do

     if (rn/rinit <= resmax .or. icycle > kcyclim ) exit cycledo

     !if (myid==0) then
    !!   print *, "At start of loop PROJ rn/rinit = ", rn, rinit, rn/rinit
    ! endif

! DEAN~ By uncommenting next line - take out projection 
!     if (icycle -1 ==-1 )then

! The next if statment allows the projection step to occur.
     if (icycle-1 == ((icycle-1)/ifreq)*ifreq) then
 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii = 1,nvhalf
          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
      enddo!ii

      do icri=1,5,2
       do kk=1,nvhalf
        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                           -co(2,1)*try(icri+1,kk,:,:,:,1)
        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                             +co(2,1)*try(icri,kk,:,:,:,1)
       enddo!kk
      enddo!icri

      do i=1,p-1 
       call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag, &
                    coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT )  !z1=M*z
       mvp = mvp + 1
       do icri=1,5,2
        do kk=1,nvhalf
         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
        enddo!k
       enddo!icri
      enddo!i
      do kk=1,nvhalf
       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,k+1
        call vecdot(vtemp(:,:,:,:,:,i),rshift(:,:,:,:,:,1),beta,MRT2)
         c(1,i) = beta(1)
         c(2,i) = beta(2)
         c2(1,i) = c(1,i)
         c2(2,i) = c(2,i)
      enddo ! i

      do ii=1,k+1
        do jj=1,k 
          hc2(1,ii,jj) = hcnew(1,ii,jj)
          hc2(2,ii,jj) = hcnew(2,ii,jj)
        enddo ! jj
      enddo ! ii
 
      do jj =1,k
       hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(1)
      enddo ! jj

      do jj = 1,k
         do i = jj+1,k+1
            amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
            con2 = 1.0_KR2/amags
            tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
            tv(2) = 0.0_KR2
            gca(1,i,jj) = sqrt(amags)/tv(1)
            gca(2,i,jj) = 0.0_KR2
            gsa(1,i,jj) = gca(1,i,jj)*con2 &
                          *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
            gsa(2,i,jj) = gca(1,i,jj)*con2 &
                        *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
            do j = jj,k
               tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
                                             + gsa(2,i,jj)*hc2(2,i,j)
               tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
                                             - gsa(2,i,jj)*hc2(1,i,j)
               tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
                                            + gsa(2,i,jj)*hc2(2,jj,j)
               tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
                                            - gsa(2,i,jj)*hc2(1,jj,j)
               hc2(:,jj,j) = tv1(:)
               hc2(:,i,j) = tv2(:)
            enddo ! j
            tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
            tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
            tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
            tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
            c(:,jj) = tv1(:)
            c(:,i) = tv2(:)
         enddo ! i
      enddo ! jj

    ! Solve linear equation
! ~ NOTE check dimension of ss...does it need to be zeroed out?

      do i = 1,k
         ss(:,i) = c(:,i)
      enddo ! i

      con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
      const = con2*(ss(1,k)*hc2(1,k,k)+ss(2,k)*hc2(2,k,k))
      ss(2,k) = con2*(ss(2,k)*hc2(1,k,k)-ss(1,k)*hc2(2,k,k))
      ss(1,k) = const

      if (k/=1) then
         do i = 1,k-1
            ir = k - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc2(1,jj,ir) + ss(2,ir)*hc2(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc2(2,jj,ir) - ss(2,ir)*hc2(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc2(1,irm1,irm1)+ss(2,irm1)*hc2(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc2(1,irm1,irm1)-ss(1,irm1)*hc2(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (k/=1)

    ! ... Define new variable d to assist in shifting masses
    ! d is the "short" solution vector to small problem

      do jj=1,k
        d(1,jj,1) = ss(1,jj)
        d(2,jj,1) = ss(2,jj)
      enddo ! jj

! Put this in to take out of algorithum...temporarily. ONCE
!   IT WORKS NEED TO TAKE OUT THE is LOOP!!!

!     is=1
!     if(is==0) then
!
!     do is=2,nshifts
!
!        do ii=1,k+1
!           do jj=1,k
!              hc2(1,ii,jj) = hcnew(1,ii,jj)
!              hc2(2,ii,jj) = hcnew(2,ii,jj)
!           enddo ! jj
!        enddo ! ii
!
!        do ii=1,k+1
!          st(:,ii) = 0.0_KR
!        enddo ! ii
!
! NOTE ~ if first shift is not zero we need hc2 to be shifted...
!
!        do ii = 1,k+1
!           do jj = 1,k
!             st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
!             st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
!           enddo ! jj
!           srv(1,ii,1) = c2(1,ii) - st(1,ii)
!           srv(2,ii,1) = c2(2,ii) - st(2,ii)
!        enddo ! ii
!
!        do jj=1,k
!          hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(is)
!        enddo ! jj
!
! NOTE ~ d here is a "work" vector give differenet temp name..
!        NOT the short solution vector!
!
!        do ii=1,k+1
!           d(1,ii,is) = cmult(1,is)*st(1,ii) - cmult(2,is)*st(2,ii)
!           d(2,ii,is) = cmult(1,is)*st(2,ii) + cmult(2,is)*st(1,ii)
!        enddo ! ii
!
!       !if (myid == 0) then
!       ! print *, "is, ourrhshereis =", is, d(:,1:2,is)
!       ! print *, "hc2 before rotation =", hc2(:,1:3,1:2)
!       ! print *, " k before =", k
!       !endif ! myid
!
!        do jj = 1,k
!           do i = jj+1,k
!              amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
!              con2 = 1.0_KR2/amags
!              tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
!              tv(2) = 0.0_KR2
!              gca(1,i,jj) = sqrt(amags)/tv(1)
!              gca(2,i,jj) = 0.0_KR2
!              gsa(1,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
!              gsa(2,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
!              do j = jj,k
!                 tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
!                                                  + gsa(2,i,jj)*hc2(2,i,j)
!                 tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
!                                                  - gsa(2,i,jj)*hc2(1,i,j)
!                 tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
!                                                 + gsa(2,i,jj)*hc2(2,jj,j)
!                 tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
!                                                 - gsa(2,i,jj)*hc2(1,jj,j)
!                 hc2(:,jj,j) = tv1(:)
!                 hc2(:,i,j) = tv2(:)
!              enddo ! j
!              tv1(1) = gca(1,i,jj)*d(1,jj,is) + gsa(1,i,jj)*d(1,i,is) + gsa(2,i,jj)*d(2,i,is)
!              tv1(2) = gca(1,i,jj)*d(2,jj,is) + gsa(1,i,jj)*d(2,i,is) - gsa(2,i,jj)*d(1,i,is)
!              tv2(1) = gca(1,i,jj)*d(1,i,is) - gsa(1,i,jj)*d(1,jj,is) + gsa(2,i,jj)*d(2,jj,is)
!              tv2(2) = gca(1,i,jj)*d(2,i,is) - gsa(1,i,jj)*d(2,jj,is) - gsa(2,i,jj)*d(1,jj,is)
!!             d(:,jj,is) = tv1(:)
!              d(:,i,is) = tv2(:)
!           enddo ! i
!        enddo ! jj
!
! NOTE~ why are we zeroing out lower triangluar part of hc2
!       when these values should be allready zero.
!
!        do jj=1,k-1
!           do ii=jj+1,k
!              hc2(1,ii,jj) = 0.0_KR2
!              hc2(2,ii,jj) = 0.0_KR2
!           enddo ! ii
!        enddo ! jj
!
!        do ii=1,k
!           hc2(1,k+1,ii) = 0.0_KR2
!           hc2(2,k+1,ii) = 0.0_KR2
!        enddo ! ii
!
!      ! if (myid == 0) then
!      !  print *, "is, d in between =", is, d(:,1:2,is)
!      !  print *, "k in da' middle =", k
!      ! endif ! myid 
! NOTE~ back solve these linear equations
!
!     con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
!     const = con2*(d(1,k,is)*hc2(1,k,k)+d(2,k,is)*hc2(2,k,k))
!     d(2,k,is) = con2*(d(2,k,is)*hc2(1,k,k)-d(1,k,is)*hc2(2,k,k))
!     d(1,k,is) = const
!
!     if (k/=1) then
!        do i = 1,k-1
!           ir = k - i + 1
!           irm1 = ir - 1
!           do jj = 1,irm1
!              const = d(1,jj,is) - d(1,ir,is)*hc2(1,jj,ir) + d(2,ir,is)*hc2(2,jj,ir)
!              d(2,jj,is) = d(2,jj,is) - d(1,ir,is)*hc2(2,jj,ir) - d(2,ir,is)*hc2(1,jj,ir)
!              d(1,jj,is) = const
!           enddo ! jj
!           con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
!           const = con2*(d(1,irm1,is)*hc2(1,irm1,irm1)+d(2,irm1,is)*hc2(2,irm1,irm1))
!           d(2,irm1,is) = con2*(d(2,irm1,is)*hc2(1,irm1,irm1)-d(1,irm1,is)*hc2(2,irm1,irm1))
!           d(1,irm1,is) = const
!        enddo ! i
!     endif ! (k/=1)
!
!    !if (myid == 0) then
!    ! print *, "is, oursolnvector =" , is, d(:,1:2,is)
!    ! print *, "k after =", k
!    !endif ! myid
!
! NOTE ~ took out least squares solution upon DR Morgan request..
!
!        call real2complex_mat(hc2, nmaxGMRES+1, nmaxGMRES, zhc2)
!
!        call zgels('N', k+1, k, 1, zhc2, nmaxGMRES+1, zd(1,is), nmaxGMRES, zwork, lzwork, info)
!
!        call complex2real_mat(zd,nmaxGMRES,nshifts,d)
!        call complex2real_mat(zhc2, nmaxGMRES+1,nmaxGMRES, hc2)
!        
!     enddo ! is
!     
!     endif ! (is==0)
     
! Form the approximate new solution x = xt + xb.
! First zero xb then xb = V*d(1,2,:,is), and x =xt +xb
! ... Need to keep track of each solution for each shift.
!     Loop over shifts while creating the soln vector.
 
! DEAN ~ ONCE THIS ALGORITHUM WORKS, NEED TO TAKE OUT THE
!        is  LOOP AND HARD CODE IN is==1

 !          do is = 1,nshifts
               is = 1
               do i=1,k
                  sr(i) = d(1,i,is)
                  si(i) = d(2,i,is)
               enddo ! i
               do jj=1,nvhalf
                  xb(:,jj,:,:,:) = 0.0_KR2
               enddo ! jj
               do jj = 1,k
                  do icri = 1,5,2 ! 6=nri*nc
                     do i = 1,nvhalf
                        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                           + sr(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           - si(jj)*vtemp(icri+1,i,:,:,:,jj)
                        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                           + si(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           + sr(jj)*vtemp(icri+1,i,:,:,:,jj)
                     enddo ! i
                  enddo ! icri
               enddo ! jj
               do i = 1,nvhalf
                  xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
               enddo ! i

! This call to Hdbletm is purely for informational purpose.

               call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag, &
                               bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
               do ii = 1,nvhalf
                  rshift(:,ii,:,:,:,is) = beg(:,ii,:,:,:) - h(:,ii,:,:,:) &
                                                +sigma(is)*xshift(:,ii,:,:,:,is)
               enddo ! ii
!           enddo ! is

!DD note ~ equivqlent to exiting the 505 loop.

   endif ! (icycle-1 == ((icycle-1)/ifreq)*ifreq) 

  !if (myid == 0) then
  ! print *, "cycles of gmres"
  !endif

! Normalize rshift(:,:,:,:,:,1) and put into the first coln. of V, the 
! orthonormal matrix whose colns span the Krylov subspace.

! NOTE ~ may not need all the res norms for shifts >=2

 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii = 1,nvhalf
          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
      enddo!ii

      do icri=1,5,2
       do kk=1,nvhalf
        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                           -co(2,1)*try(icri+1,kk,:,:,:,1)
        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                             +co(2,1)*try(icri,kk,:,:,:,1)
       enddo!kk
      enddo!icri

      do i=1,p-1 
       call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag, &
                    coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT )  !z1=M*z
       mvp = mvp + 1
       do icri=1,5,2
        do kk=1,nvhalf
         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
        enddo!k
       enddo!icri
      enddo!i
      do kk=1,nvhalf
       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   betashift = 0.0_KR

!    do is =1,nshifts
     is=1

       call vecdot(rshift(:,:,:,:,:,1),rshift(:,:,:,:,:,1), beta,MRT2)
       betashift(1,is) = beta(1)
       betashift(2,is) = beta(2)


!    enddo ! is

   const = 1.0_KR/sqrt(betashift(1,1))


! With the normalized residual form the first coln of V.
   do jj=1,nvhalf
     v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
   enddo ! jj  

   ztemp1 = DCMPLX(betashift(1,1),betashift(2,1))
   ztemp1 = sqrt(ztemp1)
   c(1,1) = REAL(ztemp1)
   c(2,1) = AIMAG(ztemp1)
   c2(1,1) = c(1,1)
   c2(2,1) = c(2,1)
   
! Perform GMRES between projections...

     do j = 1,nDR
       jp1 = j + 1
       itercount = itercount + 1
 
! LEFT OFF HERE!
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!v_j=P(A)*v_j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    try = 0.0_KR2
    y = 0.0_KR2
    do i=1,nvhalf
     try(:,i,:,:,:,1) = v(:,i,:,:,:,j)
    enddo!i
    do icri=1,5,2
     do kk=1,nvhalf
      y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
                         -co(2,1)*try(icri+1,kk,:,:,:,1)
      y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
                           +co(2,1)*try(icri,kk,:,:,:,1)
     enddo!kk
    enddo!icri
    do i=1,p-1 
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
    mvp = mvp + 1
     do icri=1,5,2
      do kk=1,nvhalf
       y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
                             +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
       y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i

    call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,y(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)

     mvp = mvp+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
!


       do i = 1,j
 
          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif
 
          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
 
          do icri = 1,5,2 ! 6=nri*nc
             do jj = 1,nvhalf
                v(icri  ,jj,:,:,:,jp1) = v(icri  ,jj,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,jj,:,:,:,i) &
                                        + beta(2)*v(icri+1,jj,:,:,:,i)
                v(icri+1,jj,:,:,:,jp1) = v(icri+1,jj,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,jj,:,:,:,i) &
                                        - beta(1)*v(icri+1,jj,:,:,:,i)
             enddo ! jj 
          enddo ! icri
       enddo ! i
 
 
       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)

       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
 
       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)

       const = 1.0_KR2/sqrt(beta(1))
       do jj=1,nvhalf
          v(:,jj,:,:,:,jp1) = const*v(:,jj,:,:,:,jp1)
       enddo

       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)
 
! DD note ~ I need to find a way of doing Givens rotations
! in the pseudocomplex routines.....

       if (j /= 1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (j /= 1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

      con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
      const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
      ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
      ss(1,j) = const

      if (j/=1) then
         do i = 1,j-1
            ir = j - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (j/=1)

      ! ... Define new variable d to assist in shifting masses
      ! d is the "short" solution vector to small problem

      do i=1,j
         d(1,i,1) = ss(1,i)
         d(2,i,1) = ss(2,i)
      enddo ! i

      do i=1,jp1
        st(:,i) = 0.0_KR2
      enddo ! i
  
      do ii = 1,jp1
         do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
         enddo ! jj
         srv(1,ii,1) = c2(1,ii) - st(1,ii)
         srv(2,ii,1) = c2(2,ii) - st(2,ii)
      enddo ! ii

      beta(1) = 0.0_KR

      do jj = 1,j+1
         beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
      enddo ! jj

      ! ... form gdr

      gdr(mvp,1) = sqrt(beta(1))

       ! ... To keep subspaces parallel for different shifts vector
       !     needs to be rotated. Use a QR factorization to do the
       !     ortognal rotation.
       !     BIG Shifting loop


!      if (is==0) then
!
!      do is = 2,nshifts
!
!         do ii=1,jp1
!            do jj=1,j
!               hcs(1,ii,jj) = hc2(1,ii,jj)
!               hcs(2,ii,jj) = hc2(2,ii,jj)
!               hcs2(1,ii,jj) = hc2(1,ii,jj)
!               hcs2(2,ii,jj) = hc2(2,ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         do jj=1,j
!            hcs(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is)
!            hcs2(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is) 
!         enddo ! jj
!
!         ! Copy the 12complex arrays into true complex arrays for use with
!         ! lapack routines
!
!         call real2complex_mat(hcs, nmaxGMRES+1, nmaxGMRES+1, zhcs)
!
!         call zgeqrf(j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! ... store R (upper triangular) in rr
!
!         do ii=1,jp1 
!            do jj=1,j
!               zrr(ii,jj) = zhcs(ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         call zungqr(j+1,j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! Copy the complex zhcs array back to hcs 
!
!         call complex2real_mat(zhcs, nmaxGMRES+1, nmaxGMRES+1, hcs)
!
!         ! ... hcs after this call is the qq part of the qr factorization
!
!         ! Now zero out crot (keeps shifts parrallel) and srvrot
!
!          do ii=1,jp1
!             zcrot(ii)=0.0_KR
!             zsrvrot(ii)=0.0_KR
!          enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,jp1
!
!               ztemp1 = DCMPLX(cmult(1,is), cmult(2,is))
!               ztemp2 = DCMPLX(c2(1,jj), c2(2,jj))
!               ztemp3 = DCMPLX(srv(1,jj,1), srv(2,jj,1))
!
!               zcrot(ii) = zcrot(ii) + ztemp1 * CONJG(zhcs(jj,ii)) * ztemp2
!               zsrvrot(ii) = zsrvrot(ii) + CONJG(zhcs(jj,ii)) * ztemp3
!
!            enddo ! jj
!         enddo ! ii
!
!         ! ... construct alpha
!
!           if ((myid == 0) .and. (is==3)) then
!             print *, "jp1, zcrot(jp1), zsrvrot(jp1) = ", jp1, zcrot(jp1), zsrvrot(jp1)
!           endif
!            
!         zalpha = zcrot(jp1)/zsrvrot(jp1)
!
!         alph(1, is) = REAL(zalpha)
!         alph(2, is) = AIMAG(zalpha)
!
!           if (myid==0)  then
!             print *,"is, zalpha, alph = ", is, zalpha, alph(:,is)
!           endif
!
!         do jj=1,j
!            ztemp1 = zalpha * zsrvrot(jj)
!            cmas(1,jj) = REAL(zcrot(jj)) - REAL(ztemp1)    ! alpha*srvrot(1,jj)
!            cmas(2,jj) = AIMAG(zcrot(jj)) - AIMAG(ztemp1)   ! alpha*srvrot(2,jj)
!         enddo ! jj
!
!           if ((myid==0) .and. (is==3))  then
!             print *,"cmas = ", cmas
!           endif
!
!         ! ... solve linear eqns problem d(1:j,is)=rr(1:j,1:j)\cmas
!
!         do ii=1,j
!            d(1,ii,is)=cmas(1,ii)
!            d(2,ii,is)=cmas(2,ii)
!         enddo ! ii
!
!         call real2complex_mat(d,nmaxGMRES,nshifts,zd)
!
!         zd(j,is) = zd(j,is)/zrr(j,j)
!
!         if (j /= 1) then 
!            do i=1,j-1
!               ir = j-i +1
!               irm1 = ir -1
!               call zaxpy(irm1,-zd(ir,is),zrr(1,ir),1,zd(1,is),1)
!               zd(irm1,is) = zd(irm1,is)/zrr(irm1,irm1)
!            enddo ! i
!         endif ! j/=1
!
!         call complex2real_mat(zd, nmaxGMRES, nshifts, d)
!
!         do ii=1,jp1
!            st(1,ii) = 0.0_KR 
!            st(2,ii) = 0.0_KR
!            srvis(1,ii) = 0.0_KR
!            srvis(2,ii) = 0.0_KR
!         enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,j
!               st(1,ii) = st(1,ii) + hcs2(1,ii,jj)*d(1,jj,is) &
!                                   - hcs2(2,ii,jj)*d(2,jj,is)
!               st(2,ii) = st(2,ii) + hcs2(1,ii,jj)*d(2,jj,is) &
!                                   + hcs2(2,ii,jj)*d(1,jj,is)
!            enddo ! jj
!            srvis(1,ii) = cmult(1,is)*c2(1,ii)-cmult(2,is)*c2(2,ii) &
!                        - st(1,ii)
!            srvis(2,ii) = cmult(1,is)*c2(2,ii)+cmult(2,is)*c2(1,ii) &
!                        - st(2,ii)
!         enddo ! ii
!
!         ! ... form the norm of srvis and put in gdr
!
!         beta(1) =0.0_KR
!
!         do jj=1,j+1
!            beta(1) = beta(1) + srvis(1,jj)*srvis(1,jj) &
!                              + srvis(2,jj)*srvis(2,jj)
!         enddo ! jj
!
!         gdr(mvp,is) = sqrt(beta(1))
!
!         if(myid==0.and.is==2) then 
!           print *, "mvp, gdr(mvp,2)=", mvp, gdr(mvp,2)
!         endif ! myid
!
!      enddo ! BIG is loop
!
!      endif ! (is==0)

       if (j>=nDR) then 

        !if(myid==0) then
        !  print *,"j=nDR", j,nDR
        !endif ! myid

! DEAN~ NEED TO HARD CODE IN is==1

!         do is = 1,nshifts
          is=1
             ! cmult(1,is) = alph(1,is) 
             ! cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i

             do jj=1,nvhalf
               xb(:,jj,:,:,:) = 0.0_KR2
             enddo ! jj

             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj

             do i = 1,nvhalf
                ! HEY!
                xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
             enddo ! i

            !if (myid == 0) then
            ! print *, "idag in proj=", idag
            !endif

             call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,&
                         iflag,bc,vecbl,vecblinv,myid,nn,ldiv, &
                         nms,lvbc,ib,lbd,iblv,MRT)
             mvp = mvp + 1
             do i=1,nvhalf
                rshift(:,i,:,:,:,is) = beg(:,i,:,:,:) - h(:,i,:,:,:) &
                                 + sigma(is)*xshift(:,i,:,:,:,is)
             enddo ! i 
           

!       enddo ! is
       endif ! (j>=nDR)
  
       if (j >= nDR) then
          betashift = 0.0_KR2

! DEAN HARD CODE IN is=1
 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do ii = 1,nvhalf
!          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
!      enddo!ii

!      do icri=1,5,2
!       do kk=1,nvhalf
!        y(icri,kk,:,:,:,1) = co(1,1)*try(icri,kk,:,:,:,1) &
!                           -co(2,1)*try(icri+1,kk,:,:,:,1)
!        y(icri+1,kk,:,:,:,1) = co(1,1)*try(icri+1,kk,:,:,:,1) &
!                             +co(2,1)*try(icri,kk,:,:,:,1)
!       enddo!kk
!      enddo!icri

!      do i=1,p-1 
!       call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag, &
!                    coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
!                    lvbc,ib,lbd,iblv,MRT )  !z1=M*z
!       do icri=1,5,2
!        do kk=1,nvhalf
!         y(icri  ,kk,:,:,:,1) = y(icri ,kk,:,:,:,1) &
!                               +co(1,i+1)*try(icri,kk,:,:,:,i+1) &
!                               -co(2,i+1)*try(icri+1,kk,:,:,:,i+1)
!         y(icri+1,kk,:,:,:,1) = y(icri+1,kk,:,:,:,1) &
!                               +co(1,i+1)*try(icri+1,kk,:,:,:,i+1) &
!                               +co(2,i+1)*try(icri,kk,:,:,:,i+1)   !y=P(A)*r
!        enddo!k
!       enddo!icri
!      enddo!i
!      do kk=1,nvhalf
!       rshift(:,kk,:,:,:,1) = y(:,kk,:,:,:,1) ! Let rshift=P(A)*rshift
!      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         do is=1,nshifts
          is=1
             call vecdot(rshift(:,:,:,:,:,is), rshift(:,:,:,:,:,is), beta, MRT2)
          
! ~ NOTE ... this should be a sqrt of beta.....

             betashift(1,is) = sqrt(beta(1))
             betashift(2,is) = 0.0_KR2
             rnale(icycle,is) = sqrt(beta(1))
!            cmult(1,is) = alph(1,is)
!            cmult(2,is) = alph(2,is)
!         enddo ! is

! use betashift (1,1) to check convergence of residual and ultimately exit.

         rn = betashift(1,1)

           !if(myid==0) then
           !  print *, "Print before write to LOG"
           !endif ! myid

          ! if (j >= nDR) then 
             if (myid==0) then
                open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                     form="formatted",status="old",position="append")
               !BS432   write(unit=8,fmt=*) "gmresdrproject-rn",mvp,rn/rinit
!                  write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,gdr(mvp,1)/rinit!,gdr(mvp,2), gdr(mvp,3), gdr(mvp,4)
                 !print *, betashift(1,1), betashift(1,2), betashift(1,3)
                close(unit=8,status="keep")
             endif
          ! endif ! (j >= nDR)

      !   if (myid ==0) then
      !      open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
      !           form="formatted",status="old",position="append")
!     !      write(unit=8,fmt=*) itercount,betashift(1,1),betashift(1,2),betashift(1,3)
      !      write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
      !            gdr(mvp,1)/rinit!, " ", gdr(mvp,2), " ", gdr(mvp,3), " ", gdr(mvp,4)
      !      close(unit=8,status="keep")
      !   endif ! myid

          rshift = 0.0_KR2

!         if (myid == 0) then
!            print *,"Checking rshift, v = ", v
!         endif

          do ii=1,nDR+1
             do icri = 1,5,2 ! 6=nri*nc
                do jj=1,nvhalf
                 rshift(icri  ,jj,:,:,:,1) = rshift(icri  ,jj,:,:,:,1) & 
                                            + v(icri,jj,:,:,:,ii)*srv(1,ii,1) &
                                            - v(icri+1,jj,:,:,:,ii)*srv(2,ii,1)
                 rshift(icri+1,jj,:,:,:,1) = rshift(icri+1,jj,:,:,:,1) & 
                                             +v(icri,jj,:,:,:,ii)*srv(2,ii,1) &
                                             + v(icri+1,jj,:,:,:,ii)*srv(1,ii,1)
                enddo ! jj
             enddo ! icri
          enddo ! ii

!          rn=gdr(mvp,1)

        !if (myid==0) then
        !   print *, "At bottom of loop PROJ rn/rinit = ", rn, rinit, rn/rinit, "j,mvp=",j,mvp
        !endif
        !if (myid==0) then
        !   print *, "end of gmres "
        !endif

! Don't need this vecdot

         call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

! NOTE~ since I took sqrt of betashift above don't need to here,

         const = 1.0_KR/betashift(1,1)

         do jj=1,nvhalf
            v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
         enddo ! jj

         icycle = icycle + 1

      endif ! (j>=nDR)
 
    enddo ! j BIG J LOOP

 enddo cycledo

 !  if (isignal == 1) then 
 !    if (myid == 0) then
 !        print *, "vk+1 gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  else ! isignal
 !    if (myid == 0) then
 !        print *, "RHS gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  endif ! isiganl

  ! if(myid==0) then
  !   print *, "Leaveing PROJ"
  ! endif ! myid


 end subroutine ppgmresproject

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


!--------------------------------------------------------------------------
  subroutine landr(rwdir,phi,x,LDR,resmax,itermin,itercount,u,GeeGooinv, &
                   iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                   lvbc,ib,lbd,iblv,hoption,MRT,MRT2)


! Direct implementation of the Lanczos algorithm
! LAN-DR(n,k) is similar to GMRES-DR(n,k)matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=phi for the vector x but for M hermetian. In this case, the
! Arnoldi algorithm is replaced by the Lanczos algorithm and the matrix
! H is real, symmetric. See section 3.2 in Morgan's paper.
!  
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   LDR(1)=n is the maximum dimension of the subspace.
!   LDR(2)=k is the number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by landr.

! Notes:
    !With Lanczos algorithm we might be able to save on storage. For the moment,
    !I'll implement the algorithm directly and keep a space for all the vectors
    !defined as "vc" below. We can work on the storage later. 

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:)         :: LDR
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             hoption,myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    ! Local variables
    integer(kind=KI)::ierr,nDR,kDR,icycle,icri,idag,i,j,k,ii,jj,kk,isite,jmin
    real(kind=KR2),dimension(6,ntotal,4,2,8,nmaxGMRES+1)::vc
    real(kind=KR2),dimension(6,ntotal,4,2,8)::r,h,w,ht,dphi
    real(kind=KR2),dimension(2)::tv
    real(kind=KR2),dimension(nmaxGMRES)::ym
    real(kind=KR2),dimension(nmaxGMRES+1,nmaxGMRES)::hbar
    real(kind=KR2),dimension(nmaxGMRES+1,nmaxGMRES)::hbarnew
    real(kind=KR2),dimension(nmaxGMRES+1,nmaxGMRES)::hbarnew_temp
    real(kind=KR2),dimension(nmaxGMRES,kmaxGMRES)::pk
    real(kind=KR2),dimension(nmaxGMRES+1,kmaxGMRES+1)::pkp1
    real(kind=KR2)::normnum,resnum,normnum_hermetian,resnum_hermetian,&
                    resnumh_approx,const,beta,cr,sumvalue    

    !Linear System of Equations
    character(len=1)::uplo,id
    real(kind=KR2),allocatable,dimension(:,:)    ::zl
    real(kind=KR2),allocatable,dimension(:,:)    ::bl
    integer(kind=KI),allocatable,dimension(:)    ::ipivl
    real(kind=KR2),allocatable,dimension(:)      ::workl
    integer(kind=KI)::info,ilwork,nrhs,jstart,kj 

    !Eigenvalue part using lapack subroutine
    character::jobz
    real(kind=KR2),allocatable,dimension(:)  ::dh,tmpeigen
    real(kind=KR2),allocatable,dimension(:,:)::zh,zhtmp
    real(kind=KR2),dimension(kcyclim,nmaxGMRES)::eigresnorm1
    character(len=5)::eigtrailer


    !Check that we are dealing with a hermetian system
    if((hoption.ne.1).and.(hoption.ne.2)) then 
      if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a30,i5)")"hoption=",hoption
        write(unit=8,fmt="(a70)")"Error: calling lan-dr for non-hermetian system"
        close(unit=8,status="keep")
      endif
      !call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      stop
     endif




    !Adjust the right-hand side to the appropriate choice
    dphi=phi


    !Multiply the right-hand side by gamma5 if solving gamma5*M*x=gamma5*phi
    if(hoption==1) then 
      call gam5x(dphi)
    endif

    !Multiply the right-hand side by M^dagger if solving M^dagger*M*x=M*^dagger*phi
    if(hoption==2) then
      idag=1
      call Hdbletm(h,u,GeeGooinv,dphi,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      dphi= h
    endif
   
    !Assume initial guess to be zero
    x=0.0_KR2
  
    eigtrailer="eigxx"
    
    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a30)")"started landr"
     write(unit=8,fmt='(a10,d20.10,a10,d20.10)')"kappa=",kappa(1),"mu=",kappa(2) 
     close(unit=8,status="keep")
    endif


    if(myid==0) then
     open(unit=12,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",status="old",&
          action="write",form="formatted",position="append")
     write(unit=12,fmt='(a40)')"-------LAN-DR---------"
     write(unit=12,fmt='(a10,d20.10,a10,d20.10)')"kappa=",kappa(1),"mu=",kappa(2) 
     write(unit=12,fmt='(a5,a5,a20,a5,a20)')"eig","#","eig value","cycle","eig-res-norm"
     close(unit=12,status="keep")
    endif


    nDR=LDR(1)
    kDR=LDR(2)


    !Set some parameters used in solving the small linear system
    !and the eigenvalue problems by the lapack libraray routines.
    uplo='U' 
    nrhs=1
    jobz='V' 
    id='d'

    !initializations
    itercount=0
    vc=0.0_KR2
    vtemp=0.0_KR2
    r=0.0_KR2
    h=0.0_KR2
    w=0.0_KR2
    ym=0.0_KR2
    cr=0.0_KR2
    hbar=0.0_KR2
    hbarnew=0.0_KR2
    tcnew=0.0_KR2



    !Assume zero initial guess
    x=0.0_KR
    
    !STEP 1: Start
 
    !Compute r0=dphi and vc_1=r0/||r0||_2 and the actual residual norm 

    do i = 1,nvhalf
     r(:,i,:,:,:) = dphi(:,i,:,:,:) 
    enddo ! i

    call vecdot(r,r,tv,MRT2)
    
    normnum_hermetian=sqrt(tv(1))
    beta=normnum_hermetian
    const=1.0_KR2/normnum_hermetian
    cr=beta
    do i=1,nvhalf
     vc(:,i,:,:,:,1)=const*r(:,i,:,:,:)
    enddo !i





    call vecdot(phi,phi,tv,MRT2)
    normnum=sqrt(tv(1))

    resnum=1.0_KR
    resnum_hermetian=1.0_KR
    resnumh_approx=1.0_KR

   if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a70,i10,3d20.10)")"landr-mvp,rresnum,rresnum_hermetian,rresnum_herme_approx",&
                                           itercount,resnum,resnum_hermetian,resnumh_approx
     close(unit=8,status="keep")
   endif

    icycle=1

    maindo: do

    if(icycle.eq.1) then
     jstart=1
    else
     jstart=kDR+1
    endif

    jdo: do j=jstart,nDR


     !STEP 1: First cycle is a purely Lanczos (see Saad's first edition, algorithm 6.15 )
     if( ((icycle.eq.1).and.(j.eq.1)).or.((icycle.ne.1).and.(j.eq.kDR+1)) ) then 
      jmin=1
     else
      jmin=j-1
     endif


     if(hoption==2) then
    
      idag=0 
      call Hdbletm(ht,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,&
                   bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      itercount=itercount+1

      idag=1
      call Hdbletm(w,u,GeeGooinv,ht,idag,coact,kappa,iflag,&
                   bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      itercount=itercount+1

     else
      idag=0
      call Hdbletm(w,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,&
                    bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      call gam5x(w)
      itercount=itercount+1
      
     endif
     
     do i=jmin,j


      if((icycle.eq.1).or.(j.ne.kDR+1)) then

       if(i.eq.j-1) then
          hbar(i,j)=hbar(j,i)
       else
         call vecdot(vc(:,:,:,:,:,i),w,tv,MRT2)
         hbar(i,j)=tv(1)
       endif      
       do k=1,nvhalf
        do icri=1,5,2
         w(icri,k,:,:,:)=w(icri,k,:,:,:)-hbar(i,j)*vc(icri,k,:,:,:,i) !+&
                         !tv(2)*vc(icri+1,k,:,:,:,i)
         w(icri+1,k,:,:,:)=w(icri+1,k,:,:,:)-hbar(i,j)*vc(icri+1,k,:,:,:,i) !-&
                          !tv(2)*vc(icri,k,:,:,:,i)
         
        enddo!icri
       enddo !k

      else

       if(i.eq.j) then
         call vecdot(vc(:,:,:,:,:,i),w,tv,MRT2)
         hbar(i,j)=tv(1)
       else
         hbar(i,j)=hbar(j,i)
       endif      
       do k=1,nvhalf
        do icri=1,5,2
         w(icri,k,:,:,:)=w(icri,k,:,:,:)-hbar(i,j)*vc(icri,k,:,:,:,i) !+&
                         !tv(2)*vc(icri+1,k,:,:,:,i)
         w(icri+1,k,:,:,:)=w(icri+1,k,:,:,:)-hbar(i,j)*vc(icri+1,k,:,:,:,i)!-&
                         !tv(2)*vc(icri,k,:,:,:,i)
         
        enddo!icri
       enddo !k

       endif
     enddo !i




  
     call vecdot(w,w,tv,MRT2)



     !Re-orthogonalization
     if(icycle.eq.1) then
       kj=j
     else
       kj=kDR
     endif

    !switch off partial re-orthogonalization for a moment.
    !Also re-orthogonalization should happen after the first cycle.

    if(.true.) then     
     do i=1,kj
      call vecdot(vc(:,:,:,:,:,i),w,tv,MRT2)
      do k=1,nvhalf
       do icri=1,5,2
        w(icri,k,:,:,:)=w(icri,k,:,:,:)-tv(1)*vc(icri,k,:,:,:,i)+&
                        tv(2)*vc(icri+1,k,:,:,:,i)
        w(icri+1,k,:,:,:)=w(icri+1,k,:,:,:)-tv(1)*vc(icri+1,k,:,:,:,i)-&
                        tv(2)*vc(icri,k,:,:,:,i)

       enddo !icri
      enddo !k
     enddo !i     
     call vecdot(w,w,tv,MRT2)
    endif !if(.true.)


     hbar(j+1,j)=sqrt(tv(1))
     if(abs(hbar(j+1,j)).LT.0.0000000000000001) then
      nDR=j
      exit jdo
     endif
     
     const=1.0_KR2/hbar(j+1,j)

     do i=1,nvhalf
      vc(:,i,:,:,:,j+1)=const*w(:,i,:,:,:)
     enddo !i

     enddo jdo


     !Solving the linear system
     ilwork=nDR*nDR
     allocate(zl(nDR,nDR))
     allocate(bl(nDR,nrhs))
     allocate(ipivl(nDR))
     allocate(workl(ilwork))

     zl=0.0_KR2
     bl=0.0_KR2
     if(icycle.eq.1) then
      bl(1,1)=cr
     else
      bl(kDR+1,1)=cr
     endif
     
     do i=1,nDR
      do ii=1,nDR
       zl(i,ii)=hbar(i,ii)
      enddo !ii
     enddo !i

     
     !call DSYSV(uplo,nDR,nrhs,zl,nDR,ipivl,bl,nDR,workl,ilwork,info)
     
     
     ym(1:nDR)=bl(1:nDR,1)
     
    !update the solution with real ym
    do j=1,nDR
     do i=1,nvhalf
      x(:,i,:,:,:)=x(:,i,:,:,:)+ym(j)*vc(:,i,:,:,:,j)
     enddo !i
    enddo !j

    !Compute the residual norm of the hermetian propblem and the actual residual norm for M*x=phi:
    !To get the residual norm of the original problem note that if solving M^dagger*M*x=M^dagger*phi 
    !, we need to calculate the residual vector using the exact formula: r=phi-M*x since M^dagger*r 
    !has a different norm than r. However, for the system gamma5*M*x=gamma5*phi, the approximate 
    !formula can be used since r and gamma5*r has the same residual norm

    
    !r=-hbar(m+1,m)*e_nDR^T*y_nDR*vc(nDR+1) and check for convergence 
    const=-hbar(nDR+1,nDR)*ym(nDR)
    resnumh_approx=abs(const)/normnum_hermetian
 

 
    idag=0   
    call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   

    do i = 1,nvhalf
     ht(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
    enddo


    call vecdot(ht,ht,tv,MRT2)

    resnum = sqrt(tv(1))/normnum

    !Calculate the exact residual norm for the hermetian system not from the approximate formula
    !which is subject to roundoff errors.
    idag=1   
    call Hdbletm(ht,u,GeeGooinv,h,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   

    do i = 1,nvhalf
     ht(:,i,:,:,:) = dphi(:,i,:,:,:) - ht(:,i,:,:,:)
    enddo


    call vecdot(ht,ht,tv,MRT2)

    resnum_hermetian = sqrt(tv(1))/normnum_hermetian


   if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a70,i10,3d20.10)")"landr-mvp,rresnum,rresnum_hermetian,rresnum_herme_approx",&
                                           itercount,resnum,resnum_hermetian,resnumh_approx
     close(unit=8,status="keep")
   endif



    if(resnumh_approx.LT.resmax) then

      deallocate(zl)
      deallocate(bl)
      deallocate(ipivl)
      deallocate(workl)

      if(myid==0) then
      open(unit=12,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",status="old",&
         action="write",form="formatted",position="append")
      write(unit=12,fmt='(a40)')"-------Eigenvalues---------"
      do i=1,nDR
       write(unit=eigtrailer(4:5),fmt='(i2.2)')i
       write(unit=12,fmt='(a5,i5,d20.10,i5,d20.10)')eigtrailer,i,heigen(i),icycle-1,eigresnorm1(icycle-1,i)
      enddo
      close(unit=12,status="keep")
      endif

      exit maindo

    endif

    !Continue if not converged
    cr=abs(const)

   !Build the eigenvector part:

   !Solving the eigenvalue part 
   allocate(zh(nDR,nDR))
   allocate(dh(nDR))
   allocate(tmpeigen(nDR))
   allocate(zhtmp(nDR,nDR))

   do i=1,nDR
    do ii=1,nDR
     zh(i,ii)=hbar(i,ii)
    enddo !ii
   enddo !i


   !call DSYEV(jobz,uplo,nDR,zh,nDR,dh,workl,ilwork,info)

   !Sort eigenvalues in ascending order according to their absolute values
   do i=1,nDR
    tmpeigen(i)=abs(dh(i))
   enddo

   ipivl=0
   call indexx(nDR,tmpeigen,ipivl)

   do i=1,nDR
    tmpeigen(i)=dh(i)
   enddo
  
   zhtmp=zh
   do i=1,nDR
    dh(i)=tmpeigen(ipivl(i))
    heigen(i)=dh(i)
    zh(:,i)=zhtmp(:,ipivl(i))
   enddo
            

   !Calculate the residual norm for the eigenvalue problem using the approximate formula:
   !residual norm for eigenvector k = abs(hbar(mDR+1,mDR)*zh(mDR,k))
   !Note: this formula is affected by roundoff erros. The more accurate way is to calculate
   !the norm of A*v-lamda*v

   do i=1,nDR
    eigresnorm1(icycle,i)=abs(hbar(nDR+1,nDR)*zh(nDR,i))
   enddo
   



   if(myid==0.and..false.) then
    open(unit=12,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",status="old",&
         action="write",form="formatted",position="append")
    write(unit=12,fmt='(a40)')"-------Eigenvalues---------"
    do i=1,nDR
      write(unit=eigtrailer(4:5),fmt='(i2.2)')i
      write(unit=12,fmt='(a5,i5,d20.10,i5,d20.10)')eigtrailer,i,heigen(i),icycle,eigresnorm1(icycle,i)
    enddo
    close(unit=12,status="keep")
   endif


   !Construct the hbarnew_k

  !define pk
  pk=0.0_KR2
  do i=1,kDR
   do ii=1,nDR
    pk(ii,i)=zh(ii,i)
   enddo !ii
  enddo !i

  !define pkp1
  pkp1=0.0_KR2
  do i=1, kDR
   do ii=1,nDR
    pkp1(ii,i)=zh(ii,i)
   enddo !ii
  enddo !i
  pkp1(nDR+1,kDR+1)=1.0_KR2
  

   hbarnew=0.0_KR2
   hbarnew_temp=0.0_KR2


   do i=1,nDR+1
    do ii=1,kDR
     hbarnew_temp(i,ii)=0.0_KR2
     do kk=1,nDR
      hbarnew_temp(i,ii)=hbarnew_temp(i,ii)+&
                         hbar(i,kk)*pk(kk,ii)
     enddo !kk
    enddo !ii
   enddo !i

   do i=1,kDR+1
    do ii=1,kDR
     hbarnew(i,ii)=0.0_KR2
     do kk=1,nDR+1
      hbarnew(i,ii)=hbarnew(i,ii)+&
                    pkp1(kk,i)*hbarnew_temp(kk,ii)
     enddo !kk
    enddo !ii
   enddo !i
   

   !Construct the new V from the old V and the eigenvectors

   vtemp=0.0_KR2
   do i=1,kDR+1
    vtemp(:,:,:,:,:,i)=0.0_KR2
    do isite=1,nvhalf
     do ii=1,nDR+1
      vtemp(:,isite,:,:,:,i)=vtemp(:,isite,:,:,:,i)+&
                             vc(:,isite,:,:,:,ii)*pkp1(ii,i)
     enddo !ii
    enddo !isite
   enddo !i

     
   !orthogonalize the kDR+1 vector againist the prvious kDR vectors


   if(.true.) then
   do i=1,kDR
    call vecdot(vtemp(:,:,:,:,:,i),vtemp(:,:,:,:,:,kDR+1),tv,MRT2)
    do isite=1,nvhalf
     do icri=1,5,2
      vtemp(icri,isite,:,:,:,kDR+1)=vtemp(icri,isite,:,:,:,kDR+1)&
          -tv(1)*vtemp(icri,isite,:,:,:,i)+tv(2)*vtemp(icri+1,isite,:,:,:,i) 

      vtemp(icri+1,isite,:,:,:,kDR+1)=vtemp(icri+1,isite,:,:,:,kDR+1)&
          -tv(1)*vtemp(icri+1,isite,:,:,:,i)-tv(2)*vtemp(icri,isite,:,:,:,i)
     enddo !icri
    enddo !isite
   enddo !i
   endif





    hbar(1:kDR+1,1:kDR)=hbarnew(1:kDR+1,1:kDR)
    vc(:,1:nvhalf,:,:,:,1:kDR+1)=vtemp(:,1:nvhalf,:,:,:,1:kDR+1)
    tcnew(1:kDR+1,1:kDR)=hbarnew(1:kDR+1,1:kDR)

    icycle = icycle + 1


    deallocate(zl)
    deallocate(bl)
    deallocate(ipivl)
    deallocate(workl)
    deallocate(zh)
    deallocate(dh)
    deallocate(tmpeigen)
    deallocate(zhtmp)

   enddo maindo
end subroutine landr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine cg(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv, &
                  iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                  lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
! LAN-DR(n,k) is similar to GMRES-DR(n,k)matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=phi for the vector x but for M hermetian.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   LDR(1)=n is the maximum dimension of the subspace.
!   LDR(2)=k is the number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by landr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             hoption,myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    ! Local variables
    integer(kind=KI)::ierr,icycle,j,k,icri,i,mvp,idag,ifreq,isite
    real(kind=KR2),dimension(6,ntotal,4,2,8)::h1,h,r,p,dphi
    real(kind=KR2),dimension(2)::tv1,tv2
    real(kind=KR2)::const,normnum,resnum,normnum_hermetian,resnum_hermetian,&
                    alpha,beta
   
    !Check that we are dealing with a hermetian system
    if((hoption.ne.1).and.(hoption.ne.2)) then 
      if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a30,i5)")"hoption=",hoption
        write(unit=8,fmt="(a70)")"Error: calling cg for non-hermetian system"
        close(unit=8,status="keep")
      endif
      !call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      stop
    endif



    !Adjust the right-hand side to the appropriate choice
    dphi=phi


    !Multiply the right-hand side by gamma5 if solving gamma5*M*x=gamma5*phi
    if(hoption==1) then 
      call gam5x(dphi)
    endif

    !Multiply the right-hand side by M^dagger if solving M^dagger*M*x=M*^dagger*phi
    if(hoption==2) then
      idag=1
      call Hdbletm(h,u,GeeGooinv,dphi,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      dphi=h
    endif
   
    !Assume initial guess to be zero
    x=0.0_KR2

    !calculate the actual residual norm every 30 iterations (this requires one or two extra 
    !matrix-products).  
    ifreq=30
  
    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a50)") "-------------------------------------------"
     write(unit=8,fmt="(a30,2d20.10)")"started cg",kappa
     write(unit=8,fmt="(a50)")"--------------------------------------------"
     close(unit=8,status="keep")
    endif

    icycle=1
    h=0.0_KR2
    r=0.0_KR2
    p=0.0_KR2
    alpha=0.0_KR2
    beta=0.0_KR2
    tv1=0.0_KR2
    tv2=0.0_KR2
    mvp=0
    itercount=0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

   call vecdot(phi,phi,tv1,MRT2)
   normnum=sqrt(tv1(1))
   
   call vecdot(dphi,dphi,tv1,MRT2)
   normnum_hermetian=sqrt(tv1(1)) 

   resnum=1.0_KR
   resnum_hermetian=1.0_KR
  


    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a60,i10,2d20.10)")"cg-mvp,resnum,resnum_hermetian",&
                                      mvp,resnum,resnum_hermetian
     close(unit=8,status="keep")
    endif

    !Compute r0=dphi
    do i = 1,nvhalf
     r(:,i,:,:,:) = dphi(:,i,:,:,:) 
     p(:,i,:,:,:) = r(:,i,:,:,:)
    enddo ! i

  
    j=0
    maindo: do

      if(hoption==2) then
        idag=0
        call Hdbletm(h1,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

        idag=1
        call Hdbletm(h,u,GeeGooinv,h1,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

        mvp=mvp+2
      else
        idag=0
        call Hdbletm(h,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        call gam5x(h)
        mvp=mvp+1
      endif

      call vecdot(p,h,tv1,MRT2)
      call vecdot(r,r,tv2,MRT2)
      
      alpha=tv2(1)/tv1(1)

      do i=1,nvhalf
        x(:,i,:,:,:)=x(:,i,:,:,:)+alpha*p(:,i,:,:,:)
      enddo !i


      do i=1,nvhalf
        r(:,i,:,:,:)  =r(:,i,:,:,:)-alpha*h(:,i,:,:,:)
      enddo !i

      tv1=0.0_KR2
      call vecdot(r,r,tv1,MRT2)
      beta=tv1(1)/tv2(1)

      resnum_hermetian=sqrt(tv1(1))/normnum_hermetian

      !Calculate the residual norm  of the original problem every 30 iterations
      if(mvp.eq.((mvp/ifreq)*ifreq)) then
    
          idag=0
          call Hdbletm(h1,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          
          do isite=1,nvhalf
           h1(:,isite,:,:,:)=phi(:,isite,:,:,:)-h1(:,isite,:,:,:)
          enddo

          call vecdot(h1,h1,tv1,MRT2)
          resnum=sqrt(tv1(1))/normnum
      
          if (myid==0) then
            open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
            write(unit=8,fmt="(a60,i10,2d20.10)")"cg-mvp,rresnum,rresnum_hermetian",&
                                             mvp,resnum,resnum_hermetian
            close(unit=8,status="keep")
          endif
          if(resnum_hermetian.LT.resmax) exit maindo

      endif       

      do i=1,nvhalf
       p(:,i,:,:,:)=r(:,i,:,:,:)+beta*p(:,i,:,:,:)
      enddo !i    
      j=j+1
 
    enddo maindo
   
end subroutine cg
! - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 subroutine cg_proj(rwdir,phi,x,LDR,resmax,itermin,itercount,u,GeeGooinv, &
                  iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                  lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
! LAN-DR(n,k) is similar to GMRES-DR(n,k)matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=phi for the vector x but for M hermetian.
! INPUT:
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   LDR(1)=n is the maximum dimension of the subspace.
!   LDR(2)=k is the number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by landr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:)         :: LDR
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             hoption,myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    ! Local variables
    integer(kind=KI)::ierr,nDR,kDR,icycle,j,k,icri,i,mvp,idag,ifreq,isite
    real(kind=KR2),dimension(6,ntotal,4,2,8)::h1,h,r,p,dphi
    real(kind=KR2),dimension(2)::tv1,tv2
    real(kind=KR2)::const,normnum,resnum,normnum_hermetian,resnum_hermetian,alpha,beta
    real(kind=KR2),dimension(2,kmaxGMRES)::ym,ym_approx
    real(kind=KR2),dimension(2,kmaxGMRES,1)::proj_rhs
    real(kind=KR2),dimension(2,kmaxGMRES,kmaxGMRES)::proj_lhs
    integer(kind=KI),dimension(kmaxGMRES)::ipiv_proj
    real(kind=KR2),dimension(2,kmaxGMRES,kmaxGMRES)::dotprods
   
    !Check that we are dealing with a hermetian system
    if((hoption.ne.1).and.(hoption.ne.2)) then 
      if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a30,i5)")"hoption=",hoption
        write(unit=8,fmt="(a70)")"Error: calling cg_proj for non-hermetian system"
        close(unit=8,status="keep")
      endif
      !call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      stop
    endif



    !Adjust the right-hand side to the appropriate choice
    dphi=phi


    !Multiply the right-hand side by gamma5 if solving gamma5*M*x=gamma5*phi
    if(hoption==1) then 
      call gam5x(dphi)
    endif

    !Multiply the right-hand side by M^dagger if solving M^dagger*M*x=M*^dagger*phi
    if(hoption==2) then
      idag=1
      call Hdbletm(h,u,GeeGooinv,dphi,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      dphi=h
    endif
   
    !Assume initial guess to be zero
    x=0.0_KR2

    !calculate the actual residual norm every 30 iterations (this requires one or two extra 
    !matrix-products).  
    ifreq=30
  
    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a50)") "-------------------------------------------"
     write(unit=8,fmt="(a30,2d20.10)")"started cg_proj",kappa
     write(unit=8,fmt="(a50)")"--------------------------------------------"
     close(unit=8,status="keep")
    endif


    nDR = LDR(1)
    kDR = LDR(2)
    icycle=1
    h=0.0_KR2
    r=0.0_KR2
    p=0.0_KR2
    alpha=0.0_KR2
    beta=0.0_KR2
    tv1=0.0_KR2
    tv2=0.0_KR2
    mvp=0
    itercount=0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 

   call vecdot(phi,phi,tv1,MRT2)
   normnum=sqrt(tv1(1))
   
   call vecdot(dphi,dphi,tv1,MRT2)
   normnum_hermetian=sqrt(tv1(1)) 

   resnum=1.0_KR
   resnum_hermetian=1.0_KR

   if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
     write(unit=8,fmt="(a60,i10,2d20.10)")"cg-proj,rresnum,rresnum_hermetian",&
                                      mvp,resnum,resnum_hermetian
     close(unit=8,status="keep")
   endif

    !Compute r0=dphi
    do i = 1,nvhalf
     r(:,i,:,:,:) = dphi(:,i,:,:,:) 
    enddo ! i

    !MinRes Projection
    if(.true.) then

     !Build the r.h.s. v^dagger*r0
     do j=1,kDR
      call vecdot(vtemp(:,:,:,:,:,j),r,proj_rhs(:,j,1),MRT2)
      const=1.0_KR2/heigen(j)
      ym(:,j)=const*proj_rhs(:,j,1)
     enddo


     !projection
     do j=1,kDR
      do isite=1,nvhalf
       do icri=1,5,2
        x(icri,isite,:,:,:)=x(icri,isite,:,:,:)+ym(1,j)*vtemp(icri,isite,:,:,:,j)&
                           -ym(2,j)*vtemp(icri+1,isite,:,:,:,j)
        x(icri+1,isite,:,:,:)=x(icri+1,isite,:,:,:)+ym(1,j)*vtemp(icri+1,isite,:,:,:,j)&
                             +ym(2,j)*vtemp(icri,isite,:,:,:,j)
       enddo !icri
      enddo !isite
     enddo !j

    endif !if(.true.)
  

    !STEP 1: 
    !Compute r0 after the projection

    if(hoption==2) then

      idag=0
      call Hdbletm(h1,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      idag=1
      call Hdbletm(h,u,GeeGooinv,h1,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT) 

    else

      idag=0
      call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      
      call gam5x(h)      

    endif


    do i = 1,nvhalf
     r(:,i,:,:,:) = dphi(:,i,:,:,:) - h(:,i,:,:,:)
     p(:,i,:,:,:) = r(:,i,:,:,:)
    enddo ! i

  
    j=0
    maindo: do

      if(hoption==2) then
        idag=0
        call Hdbletm(h1,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

        idag=1
        call Hdbletm(h,u,GeeGooinv,h1,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

        mvp=mvp+2
      else
        idag=0
        call Hdbletm(h,u,GeeGooinv,p,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        call gam5x(h)
        mvp=mvp+1
      endif

      call vecdot(p,h,tv1,MRT2)
      call vecdot(r,r,tv2,MRT2)
      
      alpha=tv2(1)/tv1(1)

      do i=1,nvhalf
        x(:,i,:,:,:)=x(:,i,:,:,:)+alpha*p(:,i,:,:,:)
      enddo !i


      do i=1,nvhalf
        r(:,i,:,:,:)  =r(:,i,:,:,:)-alpha*h(:,i,:,:,:)
      enddo !i

      tv1=0.0_KR2
      call vecdot(r,r,tv1,MRT2)
      beta=tv1(1)/tv2(1)

      resnum_hermetian=sqrt(tv1(1))/normnum_hermetian

      !Calculate the residual norm  of the original problem every 30 iterations
      if(mvp.eq.((mvp/ifreq)*ifreq)) then
    
          idag=0
          call Hdbletm(h1,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          
          do isite=1,nvhalf
           h1(:,isite,:,:,:)=phi(:,isite,:,:,:)-h1(:,isite,:,:,:)
          enddo

          call vecdot(h1,h1,tv1,MRT2)
          resnum=sqrt(tv1(1))/normnum
      
          if (myid==0) then
            open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
            write(unit=8,fmt="(a60,i10,2d20.10)")"cg_proj,rresnum,rresnum_hermetian",&
                                             mvp,resnum,resnum_hermetian
            close(unit=8,status="keep")
          endif
          if(resnum.LT.resmax) exit maindo

      endif       

      do i=1,nvhalf
       p(:,i,:,:,:)=r(:,i,:,:,:)+beta*p(:,i,:,:,:)
      enddo !i    
      j=j+1
 
    enddo maindo
   
 end subroutine cg_proj
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine landr_proj(rwdir,phi,x,nGMRES,resmax,itermin,itercount,u,GeeGooinv, &
                       iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                       lvbc,ib,lbd,iblv,hoption,MRT,MRT2)

    use shift
    
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in), dimension(:)            :: nGMRES
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             hoption,myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    !Local variables
    real(kind=KR)::resmax_proj

    resmax_proj=1.0e-08
   
    if(isignal == 1 ) then
      !call cg(rwdir,phi,x,resmax,itermin,itercount,u,GeeGooinv, &
      !            iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
      !            lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
      call landr(rwdir,phi,x,nGMRES,resmax,itermin,itercount,u,GeeGooinv, &
                  iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                  lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
    
    else

       call cg_proj(rwdir,phi,x,nGMRES,resmax_proj,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
    endif

 end subroutine landr_proj

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine minresdr5EIG(rwdir,be,bo,xe,xo,ldr,resmax,itermin,itercount,u,GeeGooinv, &
                      iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                      lvbc,ib,lbd,iblv,MRT,MRT2)

! ldr(1)=mdr, ldr(2)=kdr
! MinRes-DR(mdr,kdr) that solves the non-hermetian problem M*x=phi by converting
! it into a hermetain problem. For Wilson fermions this can be done through
! multiplication with gamma5 and solving gamma5*M*x=gamma5*phi. For all cases
! we can use M^dagge*M*x=M^dagger*phi.
!
! INPUT:
!   rwdir is the read and write directory
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   mdr is the maximum dimension of the subspace.
!   kdr is the number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by landr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI),intent(in), dimension(:)             :: ldr
    integer(kind=KI), intent(in)                          :: itermin,iflag,idag,&
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    !Local variables
    real(kind=KR2), dimension(6,ntotal,4,2,8) :: dphie,dphio,h,r0,r0e,r0o,we,wo,&
                                                tmpOO,tmpEE,xtmpEE,xtmpOO,wxe,wxo,&
                                                xre,xro,tmpE,tmpO,xeTMP,xoTMP
    real(kind=KR2), dimension(6)              :: vte,vto
    !real(kind=KR), allocatable,dimension(:,:,:,:,:,:) :: ve,vo
    real(kind=KR2), dimension(6,ntotal,4,2,8,nmaxGMRES+1) ::vtempe,vtempo
    real(kind=KR2),allocatable,dimension(:,:)::tc,tcs,tcs2,tbar
    real(kind=KR2),allocatable,dimension(:)::cvec,cvec2,cvecs,srv
    real(kind=KR2),allocatable,dimension(:)::ym
    real(kind=KR2), dimension(2)::tv,tv1,tv2,xrn,xtmp1,xtmp2,dote,doto,dot,tmp1,tmp2          
    real(kind=KR2) ::normnum,normnum_hermetian,const,cr,f,g,cs,sn,f1,rn,rn2,tmp,&
                    srvnorm,residual_norm
    real(kind=KR2), dimension(201,201) :: orthog ! orthog should have dims (mdr+1,mdr+1) TW                               
    real(kind=KR2), dimension(160,160) :: evectorthog !evectorthog should have dims (kdr,kdr) TW   
    real(kind=KR2), dimension(2,160)   :: tmpevect1,tmpevect2 !these should have dims (2,kdr) TW
    real(kind=KR2), dimension(201)     :: tmpvec                                                             
    real(kind=KR2), dimension(160)     :: tmpevect ! should have dim(kdr) TW 
    integer(kind=KI)::isite,icri,j,i,jj,ii,mdr,kdr,icycle,jstart,jmin,id,ibleo,ieo,kk,gblclr
    real(kind=KR2),allocatable,dimension(:,:)::lslhs,lsrhs
    real(kind=KR2),allocatable,dimension(:)::work
    integer(kind=KI)::nrhs,lwork,info,lworke,exists
    character(len=1)::trans,jobvl,jobvr
    integer(kind=KI),allocatable,dimension(:) ::ipivl
    real(kind=KR2),allocatable,dimension(:,:)  ::vec1
    real(kind=KR2),allocatable,dimension(:,:)  ::tch
    real(kind=KR2),allocatable,dimension(:) ::wr,wi,wv,tmpwr,tmpwi
    real(kind=KR2),allocatable,dimension(:,:)::vl,vr,vro
    real(kind=KR2),allocatable,dimension(:)  ::worke
    real(kind=KR2),allocatable,dimension(:,:)::pk


    !Parameters
    mdr=ldr(1)
    kdr=ldr(2)
    nrhs=1
    lwork=10*(mdr+1)
    trans='n'
    jobvl='n' !don't compute left-eigenvectors
    jobvr='v' !compute right-eigenvectors
    lworke=10*mdr

    if (myid == 0) then
      print *, 'value of mdr = ', mdr
      print *, 'value of kdr = ', kdr
    endif

    !Allocations
    !allocate(ve(6,ntotal,4,2,8,mdr+1)) !removed since ve/vo defined in shift-TW    
    !allocate(vo(6,ntotal,4,2,8,mdr+1))
    allocate(tc(mdr+1,mdr))
    allocate(tcs(mdr+1,mdr))
    allocate(tcs2(mdr+1,mdr))
    allocate(tbar(mdr+1,mdr))
    allocate(cvec(mdr+1))
    allocate(cvec2(mdr+1))
    allocate(cvecs(mdr+1))
    allocate(srv(mdr+1))
    allocate(ym(mdr))
    allocate(lslhs(mdr+1,mdr))
    allocate(lsrhs(mdr+1,nrhs))
    allocate(work(lwork))
    allocate(ipivl(mdr))
    allocate(vec1(mdr,nrhs))
    allocate(tch(mdr,mdr))
    allocate(wr(mdr))
    allocate(wi(mdr))
    allocate(wv(mdr))
    allocate(tmpwr(mdr))
    allocate(tmpwi(mdr))
    allocate(vl(mdr,mdr))
    allocate(vr(mdr,mdr))
    allocate(vro(mdr,mdr))
    allocate(worke(lworke))
    allocate(pk(mdr+1,kdr+1))


    !Intitalizations
    dphie=0.0_KR2
    dphio=0.0_KR2
    !h=0.0_KR
    r0e=0.0_KR2
    r0o=0.0_KR2
    !w=0.0_KR
    wo=0.0_KR2
    we=0.0_KR2
    ve=0.0_KR2
    vo=0.0_KR2
    tc=0.0_KR2
    tcs=0.0_KR2
    tcs2=0.0_KR2
    cvec=0.0_KR2
    cvec2=0.0_KR2
    srv=0.0_KR2
    ym=0.0_KR2
    lslhs=0.0_KR2
    lsrhs=0.0_KR2
    vec1=0.0_KR2
    tch=0.0_KR2
    wr=0.0_KR2
    wi=0.0_KR2
    wv=0.0_KR2
    vl=0.0_KR2
    vr=0.0_KR2
    vro=0.0_KR2
    itercount=0
    pk=0.0_KR2
    !vtempe=0.0_KR2
    !vtempo=0.0_KR2
    tcnew=0.0_KR2



    !Copy the right-hand sidei
    !dphi = phi !ORIGINAL
    dphie=be
    dphio=bo

    r0e = be ! added for testing TW 8618
    r0o = bo
    !Set the initial guess to zero
    !x=0.0_KR

    !Compute the initial residual norm assuming zero initial guess

    !call vecdot(phi,phi,tv,MRT2)  !ORIGINAL
    !normnum=sqrt(tv(1))           !ORIGINAL
    call vecdot(be,be,tv1,MRT2)
    call vecdot(bo,bo,tv2,MRT2)
    tv = tv1 + tv2
    normnum=sqrt(tv(1))

    !if (myid == 0) then
     ! print *, 'at start of minresdr5EIG, value of normnum = ', normnum
    !endif

    !Compute the right hand side of the hermetian problem
    ! --------------commented out TW 8318 -----------------!
    !if(hoption==1) then
     !call gam5x(dphi)
    !elseif(hoption==2) then
     !idag=1
     !call Hdbletm(h,u,GeeGooinv,dphi,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
     !            myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !dphi= h
    !endif
    ! -------------end comment out TW 8318 ---------------!

    !Compute the residual of the hermetian problem and vc_1
    !do isite=1,nvhalf !ORIGINAL
     !r0(:,isite,:,:,:)=dphi(:,isite,:,:,:)  !ORIGINAL
    !enddo !isite !ORIGINAL

    !call vecdot(r0,r0,tv,MRT2) !ORIGINAL
    !normnum_hermetian=sqrt(tv(1)) !ORIGINAL
    !do isite=1,nvhalf    !ORIGINAL
    ! vc(:,isite,:,:,:,1)=const*r0(:,isite,:,:,:) !ORIGINAL
    !enddo !isite !ORIGINAL

    call vecdot(r0e,r0e,tv1,MRT2)
    call vecdot(r0o,r0o,tv2,MRT2)
    tv = tv1 + tv2
    normnum_hermetian=sqrt(tv(1))
    xrn(1) = sqrt(tv(1))  ! this value will be output to CFGSPROPS.LOG
                          ! repeatedly. We are not calculating the residual
                          ! of the linear eqns to cut down on mvp's
                          ! for data collection
    rn2 = sqrt(tv(1))     ! same as above

    if (myid == 0) then
       print *, 'value of normnum_hermitian = ',normnum_hermetian
    endif

    const=1.0_KR2/normnum_hermetian
 
    !Trying this out to make it like the gmresdr5EIG algorithm TW
    ! ve(:6,:ntotal,:4,:2,:8,1)=const*r0e(:6,:ntotal,:4,:2,:8)
    ! vo(:6,:ntotal,:4,:2,:8,1)=const*r0o(:6,:ntotal,:4,:2,:8)
    do isite=1,nvhalf    !ORIGINAL
    ! vc(:,isite,:,:,:,1)=const*r0(:,isite,:,:,:) !ORIGINAL
      ve(:,isite,:,:,:,1)=const*r0e(:,isite,:,:,:)
      vo(:,isite,:,:,:,1)=const*r0o(:,isite,:,:,:)
    enddo !isite !ORIGINAL

    !if (myid == 0) then !if added by TW
     !do isite =1,nvhalf
      !print *, 've = ',ve(:,isite,:,:,:,1)
      !print *, 'vo = ',vo(:,isite,:,:,:,1)
     !enddo
    !endif

    !if (myid == 0) then 
     ! print *, 'at beginning of minresdr5EIG, value of const = ', const
    !endif


    icycle=1


    maindo: do

    if(icycle==1) then
     jstart=1
    else
     jstart=kdr+1
    endif

    jdo: do j=jstart,mdr
     ! insert gamma5 mult, followed by Heo/Hoe mult here - TW 8/3/18
     ! ------------- COMMENTED OUT TW 8318 -----------------------------------------!                                                
     !if(hoption==1) then !REPLACE THIS WITH OUR WAY OF MULTIPLYING BY GAMMA5 TW 8318        
      !idag=0
      !call Hdbletm(w,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
       !            myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      !call gam5x(w) !moved gam5x(w) to above, and operate on vc
      !itercount=itercount+1
     !elseif(hoption==2) then
      !idag=0
      !call Hdbletm(h,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
       !            myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      !idag=1
      !call Hdbletm(w,u,GeeGooinv,h,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
        !           myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      !itercount=itercount+2
     !endif !hoption==1

       do ieo=1,2
        do ibleo=1,8
          call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO,5,ieo,ibleo)
        enddo
      enddo

      gblclr = 2
      call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                   nms,lvbc,ib,lbd,iblv,MRT)

      gblclr = 1
      call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                   nms,lvbc,ib,lbd,iblv,MRT)

      we = tmpEE - kappa(1) * we
      wo = tmpOO - kappa(1) * wo

      itercount=itercount+1
     !endif !hoption==1
     ! ------------ END COMMENT OUT TW 8318 ---------------------------------------!                                       
     !if (myid == 0) then
      !print *, 'we(1,1,1,1,1) =', we(1,1,1,1,1)
      !print *, 'wo(1,1,1,1,1) =', wo(1,1,1,1,1)
     !endif

     if( ((icycle==1).and.(j==1)) .or. ((icycle>1).and.(j==kdr+1)) ) then
      jmin=1
     else
      jmin=j-1
     endif

     if(j.gt.1) then
      do i=jmin,j-1
       tc(i,j)=tc(j,i)
       !if (myid == 0) then 
        ! print *, 'tc(',i,j,') = ',tc(i,j)
       !endif
        do isite=1,nvhalf
         we(:,isite,:,:,:)=we(:,isite,:,:,:)-tc(i,j)*ve(:,isite,:,:,:,i)
         wo(:,isite,:,:,:)=wo(:,isite,:,:,:)-tc(i,j)*vo(:,isite,:,:,:,i)
         !w(:,isite,:,:,:)=w(:,isite,:,:,:)-tc(i,j)*vc(:,isite,:,:,:,i) !ORIGINAL              
        enddo !isite
      enddo
     endif !j.gt.1

     !call vecdot(vc(:,:,:,:,:,j),w,tv,MRT2) !ORIGINAL
     !tc(j,j)=tv(1) !ORIGINAL
     call vecdot(ve(:,:,:,:,:,j),we,tv1,MRT2)
     call vecdot(vo(:,:,:,:,:,j),wo,tv2,MRT2)
      tv = tv1 + tv2
      !if (myid == 0) then
       !print *, 'tv(1) into tc(j,j) = ', tv(1)
       !print *, 'tv(2) into tc(j,j) = ', tv(2)
      !endif
      tc(j,j)=tv(1) !should this be sqrt? - TW 8/9/18
      do isite=1,nvhalf
       we(:,isite,:,:,:)=we(:,isite,:,:,:)-tc(j,j)*ve(:,isite,:,:,:,j)
       wo(:,isite,:,:,:)=wo(:,isite,:,:,:)-tc(j,j)*vo(:,isite,:,:,:,j)
       !w(:,isite,:,:,:)=w(:,isite,:,:,:)-tc(j,j)*vc(:,isite,:,:,:,j) !ORIGINAL
      enddo !isite

     !Re-Orthogonalization part againist eigenvectors
     if(icycle>1) then
      do i=1,kdr
       !call vecdot(vc(:,:,:,:,:,i),w,tv,MRT2) !ORIGINAL
       call vecdot(ve(:,:,:,:,:,i),we,tv1,MRT2)
       call vecdot(vo(:,:,:,:,:,i),wo,tv2,MRT2)
       tv = tv1 + tv2
       !if (myid == 0) then
        ! print *, ' tv1(1) in reorthog = ', tv1(1)
        ! print *, ' tv1(2) in reorthog = ', tv1(2)
        ! print *, ' tv2(1) in reorthog = ', tv2(1)
        ! print *, ' tv2(2) in reorthog = ', tv2(2)
       !endif

       do isite=1,nvhalf 
        do icri=1,5,2
         !w(icri,isite,:,:,:)=w(icri,isite,:,:,:)-tv(1)*vc(icri,isite,:,:,:,i)+& !ORIGINAL
         !                       tv(2)*vc(icri+1,isite,:,:,:,i)
         !w(icri+1,isite,:,:,:)=w(icri+1,isite,:,:,:)-tv(1)*vc(icri+1,isite,:,:,:,i)-& !ORIGINAL             
         !                      tv(2)*vc(icri,isite,:,:,:,i)

         ! changing this to be like gmresdr5EIG. see reorthog section in
         ! gmresdr5EIG --- TW 8/11/18
         !we(icri,isite,:,:,:) =we(icri,isite,:,:,:)-tv1(1)*ve(icri,isite,:,:,:,i)+&
         !                      tv1(2)*ve(icri+1,isite,:,:,:,i)
         !we(icri+1,isite,:,:,:)=we(icri+1,isite,:,:,:)-tv1(1)*ve(icri+1,isite,:,:,:,i)-&
         !                      tv1(2)*ve(icri,isite,:,:,:,i)
         !wo(icri,isite,:,:,:)=wo(icri,isite,:,:,:)-tv2(1)*vo(icri,isite,:,:,:,i)+&
         !                      tv2(2)*vo(icri+1,isite,:,:,:,i)
         !wo(icri+1,isite,:,:,:)=wo(icri+1,isite,:,:,:)-tv2(1)*vo(icri+1,isite,:,:,:,i)-&
         !                      tv2(2)*vo(icri,isite,:,:,:,i)
         ! --- end comment out TW 8/11/18

         we(icri,isite,:,:,:) =we(icri,isite,:,:,:) &
                               -tv(1)*ve(icri,isite,:,:,:,i) &
                               +tv(2)*ve(icri+1,isite,:,:,:,i)
         we(icri+1,isite,:,:,:)=we(icri+1,isite,:,:,:) &
                                -tv(1)*ve(icri+1,isite,:,:,:,i) &
                                -tv(2)*ve(icri,isite,:,:,:,i)
         wo(icri,isite,:,:,:)=wo(icri,isite,:,:,:) &
                              -tv(1)*vo(icri,isite,:,:,:,i) &
                              +tv(2)*vo(icri+1,isite,:,:,:,i)
         wo(icri+1,isite,:,:,:)=wo(icri+1,isite,:,:,:) &
                                -tv(1)*vo(icri+1,isite,:,:,:,i) &
                                -tv(2)*vo(icri,isite,:,:,:,i)
        enddo
       enddo
      enddo
     endif

     call vecdot(we,we,tv1,MRT2)
     call vecdot(wo,wo,tv2,MRT2)
     tv = tv1 + tv2
     !call vecdot(w,w,tv,MRT2) !ORIGINAL    
     if(tv(1).lt.1.e-12) then
      mdr=j
      exit jdo
     endif
     tc(j+1,j)=sqrt(tv(1))
     !if (myid == 0) then
      !  print *,'tc(',j+1,j,') = ', tc(j+1,j)
     !endif
     const=1.0_KR2/tc(j+1,j)
     !if(myid == 0) then
      ! print *,'const before end jdo =',const
     !endif

     !do isite=1,nvhalf !ORIGINAL
      !vc(:,isite,:,:,:,j+1)=const*w(:,isite,:,:,:)!ORIGINAL
     !enddo !isite !ORIGINAL
    
     ! --------- Inserted to be like gmresdr5EIG - TW
     !ve(:6,:ntotal,:4,:2,:8,j+1) = const*we(:6,:ntotal,:4,:2,:8)
     !vo(:6,:ntotal,:4,:2,:8,j+1) = const*wo(:6,:ntotal,:4,:2,:8)

     do isite=1,nvhalf !ORIGINAL
      !vc(:,isite,:,:,:,j+1)=const*w(:,isite,:,:,:)!ORIGINAL
      ve(:,isite,:,:,:,j+1)=const*we(:,isite,:,:,:)
      vo(:,isite,:,:,:,j+1)=const*wo(:,isite,:,:,:)
     enddo !isite !ORIGINAL

    enddo jdo

 if (.false.) then
  if (icycle==1 .OR. icycle==2 .OR. icycle==10 .OR. icycle==100 .OR. icycle==200 .OR. icycle==500 .OR.  icycle==1000 .OR. icycle==2000) then !print out orthog. of basis vectors for the krylov subspace TW 2/18/18
   do i=1,mdr+1
    do kk=1,mdr+1

      dote=0.0_KR2
      doto=0.0_KR2

     call vecdot(ve(:,:,:,:,:,i),ve(:,:,:,:,:,kk),dote,MRT2)
     call vecdot(vo(:,:,:,:,:,i),vo(:,:,:,:,:,kk),doto,MRT2)


     orthog(i,kk)=dote(1) + doto(1)
     !orthog(i,kk)=dote(2) + doto(2)

   enddo!j
  enddo!i


        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"orthog.dat", exist=exists)
            if (.not. exists) then

               open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="new",&
               action="write",form="formatted")
               close(unit=48,status="keep")
            endif
       endif


     do i=1,mdr+1
       do kk=1,mdr+1
            if (myid==0) then
                 open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=48,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12)")"icycle=",icycle,"V'V(",i,",",kk,")=",orthog(i,kk)
                 close(unit=48,status="keep")
            endif
      enddo
    enddo



  endif !icycle
endif!true or false


if (.false.) then
  if (icycle==1 .OR. icycle==2 .OR. icycle==10 .OR. icycle==50 .OR. icycle==100 .OR. icycle==200 .OR. icycle==250 .OR. icycle==1000 .OR. icycle==2000) then
        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"tctest.dat", exist=exists)
            if (.not. exists) then

               open(unit=87,file=trim(rwdir(myid+1))//"tctest.dat",status="new",&
               action="write",form="formatted")
               close(unit=87,status="keep")
            endif
       endif
     do i = 1,mdr
      do kk = 1,mdr
       if (myid==0) then
              open(unit=87,file=trim(rwdir(myid+1))//"tctest.dat",status="old",action="write",&
              form="formatted",position="append")
              write(unit=87,fmt="(es22.12)")tc(kk,i)
              close(unit=87,status="keep")
            endif
      enddo !kk
     enddo !i
  endif !icycle
endif !true/false

    !Solve the least square problem using givens rotations:
    !------------------------------------------------------
    tcs=tc
    !cvec=0.0_KR2

    !compute cvec=V_{m+1}^dagger*r0
   if (icycle == 1) then !inserted for testing -TW 8/23/18
    do i=1,kdr+1 !Dr, Wilcox,should this be i=1,kdr+1? see the linear eqns part of the minresdr alg in the paper -Travis (was mdr+1, changed for test)                                  
     !call vecdot(vc(:,:,:,:,:,i),r0,tv,MRT2) !ORIGINAL
     !cvec(i)=tv(1) !ORIGINAL
     call vecdot(ve(:,:,:,:,:,i),r0e,tv1,MRT2) !might want to call vecdagvec here      
     call vecdot(vo(:,:,:,:,:,i),r0o,tv2,MRT2)
     tv = tv1 + tv2
     cvec(i)=tv(1)
    enddo !i

    do i = kdr+2,mdr+1
     cvec(i) = 0.0_KR2
    enddo
   endif !icycle, !end of test -TW 8/23/18

    lslhs=tc
    lsrhs(:,1)=cvec(:)
    cvecs=cvec
 
    !if (myid == 0) then
     !do i = 1,mdr+1
      ! print *, 'elements of cvec =', cvec(i)
     !enddo
    !endif

    !-----Dr. Wilcox, this next portion was commented out in the original form
    !     of the subroutine -Travis

    ! WARNINIG: The following Givens rotations part need to be fixed to work for
    ! all cycles
    ! Its current form is valid only for the first cycle. For next cycles, note
    ! that tc has
    ! a full k+1 x k+1 leading matrix.

    ! Start of Givens rotations
    !do i=1,mdr
    ! tcs2=tcs
    ! cvec2=cvec
    ! f=tcs(i,i)
    ! g=tcs(i+1,i)
    ! call rotate2(f,g,cs,sn,f1)
    ! tcs(i,i)=f1
    ! tcs(i+1,i)=0.0_KR2
    ! do jj=1,mdr
    !  if(jj.ne.i) then
    !    tcs(i,jj)=cs*tcs2(i,jj)+sn*tcs2(i+1,jj)
    !    tcs(i+1,jj)=-sn*tcs2(i,jj)+cs*tcs2(i+1,jj)
    !  endif
    ! enddo
    ! cvec(i)=cs*cvec2(i)+sn*cvec2(i+1)
    ! cvec(i+1)= -sn*cvec2(i)+cs*cvec2(i+1)
    !enddo !i


    !ym(mdr)=cvec(mdr)/tcs(mdr,mdr)
    !do i=mdr-1,1,-1
    ! tmp=0.0_KR
    ! do jj=i+1,mdr
    !  tmp=tmp+tcs(i,jj)*ym(jj)
    ! enddo
    ! ym(i)=(cvec(i)-tmp)/tcs(i,i)
    !enddo

    ! End of Givens rotations
    ! --- end the commented portion out - Travis

   !Solve the least-squares problem by using lapack subroutine dgels
   call dgels(trans,mdr+1,mdr,nrhs,lslhs,mdr+1,lsrhs,mdr+1,work,lwork,info) 
   ym(1:mdr)=lsrhs(1:mdr,1)

   !if (myid == 0) then
    ! do i = 1,mdr 
     !  print *, 'elements of ym = ', ym(i)
    ! enddo
   !endif
   !Compute the components of the short residual vector
   do i=1,mdr+1
    tmp=0.0_KR2
    do j=1,mdr
     tmp=tmp+tc(i,j)*ym(j)
    enddo !j
     srv(i)= cvecs(i)- tmp
   enddo !i

   srvnorm=0.0_KR2
   do i=1,mdr+1
    srvnorm=srvnorm+srv(i)*srv(i)
   enddo
   !if (myid == 0) then
    !print *,'srvnorm before sqrt = ', srvnorm
   !endif

   srvnorm=sqrt(srvnorm)

   !if (myid == 0) then
    !  print *, 'srvnorm = ', srvnorm
   !endif




   !Update the solution and compute the residual norm
   do i=1,mdr
    do isite=1,nvhalf
     !x(:,isite,:,:,:)=x(:,isite,:,:,:)+ym(i)*vc(:,isite,:,:,:,i) !ORIGINAL
     xe(:,isite,:,:,:)=xe(:,isite,:,:,:)+ym(i)*ve(:,isite,:,:,:,i)
     xo(:,isite,:,:,:)=xo(:,isite,:,:,:)+ym(i)*vo(:,isite,:,:,:,i)
    enddo !isite
   enddo !i

   !Short cut
   rn=abs(lsrhs(mdr+1,1))
   !if (myid == 0) then
    !  print *, 'rn = ', rn
   !endif
   !cr=rn

   !--------Dr. Wilcox, I commented out this original resiudal calc -Travis
   !direct, WILL NEED TO CHANGE THIS RESIDUAL CALC - TW 8/3/18
   !idag=0
   !call gam5x(x) !added for the computation of r = M*g5*x-b TW 8318
   !call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
    !            myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

   !do isite=1,nvhalf
    !w(:,isite,:,:,:)=phi(:,isite,:,:,:)-h(:,isite,:,:,:)
   !enddo

   !Residual norm of the original problem
   !call vecdot(w,w,tv,MRT2)
   !residual_norm=sqrt(tv(1))

   ! WILL NEED TO CHANGE THIS RESIDUAL CALC AS WELL - TW 8/3/18
   !if(hoption==1) then
    !call gam5x(w)
    !call gam5x(w)
    !do isite=1,nvhalf
    ! r0(:,isite,:,:,:)=w(:,isite,:,:,:)
    !enddo
    !call vecdot(r0,r0,tv,MRT2)
    !rn2=sqrt(tv(1)) !LEFT UNCOMMENTED FOR TESTING - TW 8/3/18
   !elseif(hoption==2) then
    !idag=1
    !call Hdbletm(r0,u,GeeGooinv,w,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
     !            myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    !call vecdot(r0,r0,tv,MRT2)
    !rn2=sqrt(tv(1))
   !endif

    !start calculation of nnew residual vector!!! -TW 8/20/18

 if (.false.) then !residual calculation, does not determine stopping criteria -Travis          
   do i = 1,mdr+1
   tmp = 0.0_KR2
    do kk = 1,mdr
     tmp = tmp + tc(i,kk)*ym(kk)
    enddo
   tmpvec(i) = tmp
   enddo
    
   do i = 1,mdr+1
     do isite = 1,nvhalf
       r0e(:,isite,:,:,:) = r0e(:,isite,:,:,:) - tmpvec(i)*ve(:,isite,:,:,:,i)
       r0o(:,isite,:,:,:) = r0o(:,isite,:,:,:) - tmpvec(i)*vo(:,isite,:,:,:,i)
     enddo
   enddo
 

    call vecdot(r0e,r0e,tv1,MRT2)
    call vecdot(r0o,r0o,tv2,MRT2)
    tv = tv1 + tv2
    tv(2) = 0.0_KR2
    rn2 = sqrt(tv(1))
 endif !true/false

   ! start residual calculation for linear eqns---- Travis
   if (.false.) then
      do ieo=1,2
        do ibleo=1,8
          call gammamult( xe(:,:,:,:,:),xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
        enddo
      enddo
 
      gblclr = 2
      call Hsingle(wxo,u,xtmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
      ! Next build H_eo * v_o.
      gblclr = 1
      call Hsingle(wxe,u,xtmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

      wxe = xtmpEE - kappa(1) * wxe
      wxo = xtmpOO - kappa(1) * wxo


      xre = wxe - be
      xro = wxo - bo



      call vecdot(xre,xre,xtmp1,MRT2)
      call vecdot(xro,xro,xtmp2,MRT2)

      xrn = xtmp1+xtmp2
      if (myid == 0) then
         print *, 'xrn(1) = ', xrn(1)
         print *, 'xrn(2) = ', xrn(2)
      endif
      xrn(1)=sqrt(xrn(1))
      xrn(2)=0.0_KR2
    endif !true/false 


      residual_norm = rn2
      !if(myid == 0) then
        !print *, 'xrn(1) = ', xrn(1)
       ! print *, 'resdiual_norm = ', residual_norm
      !endif

   if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
    write(unit=8,fmt="(a40,i10,es17.10,es17.10,es17.10,es17.10)")"minresdr5EIG,itercount,rn,xrn,srvnorm,residual",&
         itercount,rn/normnum_hermetian,xrn(1)/normnum_hermetian,srvnorm/normnum_hermetian,residual_norm/normnum
    close(unit=8,status="keep")
   endif

   if( (srvnorm/normnum_hermetian.LT. resmax) .OR. (icycle > 2000) ) then

     if(myid==0) then
      open(unit=10,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",action="write",&
      form="formatted",status="old",position="append")
      write(unit=10,fmt='(a30,2d17.10)')"#kappa,mu=",kappa(1),kappa(2)
      write(unit=10,fmt='(4a20)')"#eigen num","wr","wi","eig-res-norm"
     endif
     !printout the harmonic ritz vectors and their residual
     do i=1,kdr
      !w=0.0_KR !ORIGINAL
      we=0.0_KR2
      wo=0.0_KR2
      do j=1,kdr+1
       do isite=1,nvhalf
        !w(:,isite,:,:,:)=w(:,isite,:,:,:)+vc(:,isite,:,:,:,j)*tc(j,i) !ORIGINAL
        we(:,isite,:,:,:)=we(:,isite,:,:,:)+ve(:,isite,:,:,:,j)*tc(j,i)
        wo(:,isite,:,:,:)=wo(:,isite,:,:,:)+vo(:,isite,:,:,:,j)*tc(j,i)
       enddo
      enddo
      do isite=1,nvhalf
       !w(:,isite,:,:,:)=w(:,isite,:,:,:)-wr(i)*vc(:,isite,:,:,:,i) !ORIGINAL
       we(:,isite,:,:,:)=we(:,isite,:,:,:)-wr(i)*ve(:,isite,:,:,:,i)
       wo(:,isite,:,:,:)=wo(:,isite,:,:,:)-wr(i)*vo(:,isite,:,:,:,i)
      enddo
      !call vecdot(w,w,tv,MRT2) !ORIGINAL
      call vecdot(we,we,tv1,MRT2)
      call vecdot(wo,wo,tv2,MRT2)
      tv = tv1 + tv2
      tv(1) = sqrt(tv(1))
      tv(2) = 0.0_KR2
      if(myid==0) then
       write(unit=10,fmt='(i20,3d20.10)')i,wr(i),wi(i),sqrt(tv(1))
      endif
     enddo !i

     if(myid==0) close(unit=10,status="keep")

     !Dellocations
     !deallocate(ve)
     !deallocate(vo)
     deallocate(tc) 
     deallocate(tcs) 
     deallocate(tcs2)
     deallocate(tbar)
     deallocate(cvec)
     deallocate(cvec2)
     deallocate(cvecs)
     deallocate(srv)
     deallocate(ym)
     deallocate(lslhs)
     deallocate(lsrhs)
     deallocate(work)
     deallocate(ipivl)
     deallocate(vec1)
     deallocate(tch) 
     deallocate(wr)
     deallocate(wi)
     deallocate(wv)
     deallocate(tmpwr)
     deallocate(tmpwi)
     deallocate(vl)
     deallocate(vr)
     deallocate(vro)
     deallocate(worke)
     deallocate(pk)
     exit maindo

    endif







   !Compute the k smallest harmonic Ritz values and vectors:

   !Compute H_m^{-T}*e_m
   do i=1,mdr
    do j=1,mdr
     tch(i,j)=tc(j,i) !ORIGINAL, commented out TW 8/20/18
     ! tch(i,j) = tc(i,j) !commented out for testing! -TW 5:22pm 8/23/18
    enddo
   enddo

   !begin comment out for appropriate matrix for ritz pair calculation!!! -TW
   !8/20/18

   vec1=0.0_KR2
   vec1(mdr,1)=1.0_KR2

   call DGESV(mdr,nrhs,tch,mdr,ipivl,vec1,mdr,info) !lapack equivalent of linearsolver                                 

   !build the matrix for which to compute the eigenvectors
   tch=0.0_KR2

   do i=1,mdr
    tch(i,mdr)=tc(i,mdr)+tc(mdr+1,mdr)**2*vec1(i,1)
   enddo !i

   do i=1,mdr
    do j=1,mdr-1
     tch(i,j)=tc(i,j)
    enddo
   enddo
 
   !end comment out TW 8/20/18

   !below print out for testing -TW
   !if (myid == 0) then 
    ! print *, 'begin printing of tch'
     !do i = 1,mdr
      ! do kk = 1,mdr
       !   print *, 'tch(',i,kk,') = ',tch(i,kk)
      ! enddo
     !enddo
   !endif

   !we might need to call a hermitian solver here TW 8318
   call DGEEV(jobvl,jobvr,mdr,tch,mdr,wr,wi,vl,mdr,vr,mdr,worke,lworke,info) !lapack equivalent of eigencalc                                                


   !Sort eigenvalues in ascending order according to their absolute values
   do i=1,mdr
    tmpwr(i)=sqrt(wr(i)**2+wi(i)**2)
   enddo

   ipivl=0
   call indexx(mdr,tmpwr,ipivl)

   do i=1,mdr 
    tmpwr(i)=wr(i)
    tmpwi(i)=wi(i)
   enddo

   vro = 0.0_KR2
   vro=vr
   do i=1,mdr 
    wr(i)=tmpwr(ipivl(i))                                   
    wi(i)=tmpwi(ipivl(i))
    vr(:,i)=vro(:,ipivl(i))
   enddo

if (.false.) then
  if (icycle==1 .OR. icycle==2 .OR. icycle==10 .OR. icycle==50 .OR. icycle==100 .OR. icycle==200 .OR. icycle==250 .OR. icycle==1000 .OR. icycle==2000)  then
        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"ritztest.dat", exist=exists)
            if (.not. exists) then

               open(unit=83,file=trim(rwdir(myid+1))//"ritztest.dat",status="new",&
               action="write",form="formatted")
               close(unit=83,status="keep")
            endif
       endif
     do i = 1,mdr
      do kk = 1,mdr
       if (myid==0) then
              open(unit=83,file=trim(rwdir(myid+1))//"ritztest.dat",status="old",action="write",&
              form="formatted",position="append")
              write(unit=83,fmt="(es22.12)")vr(kk,i)
              close(unit=83,status="keep")
            endif
      enddo !kk
     enddo !i
  endif !icycle
endif !true/false
   !--- added for testing TW
   !tmpevect = 0.0_KR2
   !tmpevect(:) = vr(:,1)
   !tmpevect1(1,:) = tmpevect(:)
   !tmpevect1(2,:) = 0.0_KR2
   !call vecdotvec(tmpevect1,mdr,tmpevect1,mdr,dot)
   !if (myid == 0) then
    ! print *,'vr(:,1)*vr(:,1) =', dot(1)
    ! print *,'vr(:,1)*vr(:,1) =', dot(2)
   !endif
   ! ----end test TW


   do i=1,kdr
    heigen(i)=wr(i)
   enddo
   
   !if (myid == 0) then
    !do i = 1,kdr 
     !print *, 'real eval = ',wr(i)
     !print *, 'imag eval = ',wi(i)
    !enddo
   !endif

   !Building orthonormal basis for the kdr smallest eigenvectors: - this is the
   !orginal!! -TW 8/21/18
   ! This is standard gram-schmidt procedure, I don't understand why the ritz
   ! vectors are not orthogonal - Travis
   do j=1,kdr
    wv = 0.0_KR2
    wv(:)=vr(:,j)
    do i=1,j-1 
     tmp=0.0_KR2  
     do ii=1,mdr
      tmp=tmp+vr(ii,i)*wv(ii)
     enddo
     wv(:)=wv(:)-tmp*vr(:,i)
    enddo !i
    tmp=0.0_KR2
    do ii=1,mdr
     tmp=tmp+wv(ii)*wv(ii)
    enddo
    const=1.0_KR2/sqrt(tmp)
    vr(:,j)=const*wv(:)
   enddo !j
   ! end the original portion - TW 8/21/18
   ! end comment out for testing 8/21/18 -TW


!check orthonormality of kdr evects here TW 8/10/18

 if (.false.) then !I cant tell if the evects are not orthonormal or if I have
 ! just not implemented this test correctly --Travis
  if (icycle==1 .OR. icycle==2 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !print out orthog. of ritz vectors TW 8/10/18
   do kk = 1,kdr
   !tmpevect(:) = vr(:,kk)
    do jj = 1,kdr
    tmp = 0.0_KR2
     do i = 1,mdr
      tmp = tmp + vr(i,jj)*vr(i,kk)
     enddo
     evectorthog(jj,kk)=tmp
    enddo
   enddo



        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"ritzvectsorthog.dat", exist=exists)
            if (.not. exists) then

               open(unit=76,file=trim(rwdir(myid+1))//"ritzvectsorthog.dat",status="new",&
               action="write",form="formatted")
               close(unit=76,status="keep")
            endif
       endif


     do i=1,kdr
       do kk=1,kdr
            if (myid==0) then
                 open(unit=76,file=trim(rwdir(myid+1))//"ritzvectsorthog.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=76,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12)")"icycle=",icycle,"ev'ev(",i,",",kk,")=",evectorthog(i,kk)
                 close(unit=76,status="keep")
            endif
      enddo
    enddo



  endif !icycle
endif!true or false

   !build even and odd right evects, added by TW 
     do i=1,kdr
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do isite = 1,nvhalf
              vte(:) = 0.0_KR2
              vto(:) = 0.0_KR2
              do kk=1,mdr
                do icri = 1,5,2
                  vte(icri)   = vte(icri) &
                        + ve(icri  ,isite,id,ieo,ibleo,kk)*vr(kk,i) !&
                        !- ve(icri+1,isite,id,ieo,ibleo,kk)*vr(kk,i) !need to remove this line                      
                  vto(icri)   = vto(icri) &
                        + vo(icri  ,isite,id,ieo,ibleo,kk)*vr(kk,i) !&
                        !- vo(icri+1,isite,id,ieo,ibleo,kk)*vr(kk,i) !need to remove this line                  
                  vte(icri+1) = vte(icri+1) &
                        + ve(icri+1,isite,id,ieo,ibleo,kk)*vr(kk,i) !&
                        !+ ve(icri  ,isite,id,ieo,ibleo,kk)*vr(kk,i) & !need to
                        !remove this line
                  vto(icri+1) = vto(icri+1) &
                        + vo(icri+1,isite,id,ieo,ibleo,kk)*vr(kk,i) !&
                        !+ vo(icri  ,isite,id,ieo,ibleo,kk)*vr(kk,i) & !need to
                        !remove this line
                enddo
              enddo
              evectore(:,isite,id,ieo,ibleo,i) = vte(:)
              evectoro(:,isite,id,ieo,ibleo,i) = vto(:)
            enddo
          enddo
        enddo
      enddo
      evalue(1,i) = wr(i)
      evalue(2,i) = wi(i)
    enddo
    !end evector get - TW

    if (.false.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==500 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,kdr
           do ieo=1,2
            do ibleo=1,8
             call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
            enddo
           enddo

         gblclr = 2
         call Hsingle(tmpOO,u,tmpE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
         ! Next build H_eo * v_o.
         gblclr = 1
         call Hsingle(tmpEE,u,tmpO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
         xeTMP = tmpE - kappa(1) * tmpEE
         xoTMP = tmpO - kappa(1) * tmpOO

         tmpE = xeTMP - (evalue(1,i)*evectore(:,:,:,:,:,i))
         tmpO = xoTMP - (evalue(1,i)*evectoro(:,:,:,:,:,i))

         call vecdot(tmpE,tmpE,tmp1,MRT2)
         call vecdot(tmpO,tmpO,tmp2,MRT2)
         tmp1 = tmp1 + tmp2
         call cmpxsqrt(tmp1, tmp2)

         if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"loweigresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=88,file=trim(rwdir(myid+1))//"loweigresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=88,status="keep")
           endif
         endif

         if (myid==0) then
           open(unit=88,file=trim(rwdir(myid+1))//"loweigresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=88,fmt="(i4,a6,i7,a6,es19.12,a6,es19.12,a6,es19.12)")icycle," ", i," ",evalue(1,i)," ",evalue(2,i)," ",tmp2(1)
           close(unit=88,status="keep")
         endif


      enddo !i
     endif !icycle
    endif !.true./.false

   !Build the matrix pk(mdr+1,kdr+1)
   pk=0.0_KR2
   do i=1,kdr
    do j=1,mdr
     pk(j,i)=vr(j,i)
    enddo
    pk(mdr+1,i)=0.0_KR2
   enddo !i

  ! do i = 1,mdr
   ! pk(i,kdr+1) = 0.0_KR2
  ! enddo
  ! pk(mdr+1,kdr+1) = 1.0_KR2 !added TW 8/21/18


   !below is the old way of forming the k+1 column of P
   !it can become unstable after linnear eqns converge
   !if this doesnt work, we need to use the approach in the 
   ! minres-dr/lan-dr paper for k+1 column (same as in gmresdr5EIG)
   ! - Travis

   !begin comment out of ORIGINAL to implement the algorithm in the lan-dr
   ! minres-dr paper - Travis 8/20/18
   cvec2=0.0_KR2
   do i=1,mdr+1
    tmp=0.0_KR2
    do jj=1,mdr
     tmp=tmp+tc(i,jj)*ym(jj)
    enddo
   cvec2(i)=cvecs(i)-tmp
   enddo

   ! end comment out Travis 8/20/18

   !begin addition for forming k+1 column of p -TW 8/20/18
   !cvec2 = 0.0_KR2
   !do i = 1,mdr
    !do kk = 1,mdr
     ! tbar(i,kk) = tc(kk,i)
    !enddo
   !enddo

   !vec1 = 0.0_KR2
   !vec1(mdr,1) = 1.0_KR2
   !call DGESV(mdr,nrhs,tbar,mdr,ipivl,vec1,mdr,info) !lapack equivalent of linearsolver

   !do i = 1,mdr
    !cvec2(i) = -tc(mdr+1,mdr)*vec1(i,1)
   !enddo
   
   ! end addition for forming k+1 column of p - TW 8/20/18

   ! comment out for testing - TW 8/23/18

   ! added for testing -TW 8/23/18 4:28pm
   !do i = 1,mdr+1
    !cvec2(i) = pk(i,kdr+1)
   !enddo

   do i=1,kdr
    tmp=0.0_KR2
    do j=1,mdr+1
     tmp=tmp+pk(j,i)*cvec2(j)
    enddo
    do j=1,mdr+1
     cvec2(j)=cvec2(j)-tmp*pk(j,i)
    enddo
   enddo

   tmp=0.0_KR2
   do j=1,mdr+1
    tmp=tmp+cvec2(j)*cvec2(j)
   enddo
   const=1.0_KR2/sqrt(tmp)
   cvec2(:)=const*cvec2(:)
   do j=1,mdr+1
    pk(j,kdr+1)=cvec2(j)
   enddo

   !form new rhs for least squares -TW 8/23/18
   do i = 1,kdr+1
    tmp = 0.0_KR2
    do kk = 1,mdr+1
     tmp = tmp + pk(kk,i)*srv(kk)
    enddo
    cvec(i) = tmp
   enddo
  
   do i = kdr+2,mdr+1
    cvec(i) = 0.0_KR2
   enddo
   !end of new rhs for least squares -TW 8/23/18

   !Build V_(k+1)^(new) and tcbar_k^new
   tcs2=0.0_KR2

   do i=1,mdr+1
    do j=1,kdr
     tmp=0.0_KR2
     do ii=1,mdr
      tmp=tmp+tc(i,ii)*pk(ii,j)
     enddo
     tcs2(i,j)=tmp
    enddo
   enddo

   do i=1,kdr+1
    do j=1,kdr
     tcnew(i,j)=0.0_KR2
     do ii=1,mdr+1
      tcnew(i,j) = tcnew(i,j) + pk(ii,i)*tcs2(ii,j)
     enddo
    enddo
   enddo

   do i=1,kdr+1
    do j=1,kdr
     tc(i,j)=tcnew(i,j)
    enddo
   enddo

   !This is the analogous portion to the part in gmresdr5EIG Dr. Wilcox
   !I added a lot of portions for testing different things, but the uncommented
   !portion is the structure of the original with even/odd - Travis
   !vtemp=0.0_KR !ORIGINAL
   vtempe=0.0_KR2
   vtempo=0.0_KR2
   do j=1,kdr+1
    !vtemp(:,:,:,:,:,j)=0.0_KR !ORIGINAL
    vtempe(:,:,:,:,:,j)= 0.0_KR2
    vtempo(:,:,:,:,:,j)= 0.0_KR2
      do i=1,mdr+1
     !vtemp(:,:,:,:,:,j) = vtemp(:,:,:,:,:,j) + vc(:,:,:,:,:,i)*pk(i,j) !ORIGINAL

     ! my first additions below -- TW 
     vtempe(:,:,:,:,:,j) = vtempe(:,:,:,:,:,j) + ve(:,:,:,:,:,i)*pk(i,j)
     vtempo(:,:,:,:,:,j) = vtempo(:,:,:,:,:,j) + vo(:,:,:,:,:,i)*pk(i,j)
     ! ---end comment out of my first additions - TW


    
        enddo !i   !added tw 8/11/18, testing
   enddo !j

   !vc(:,:,:,:,:,1:kdr+1)=vtemp(:,:,:,:,:,1:kdr+1)  !ORIGINAL
   !ve(:,:,:,:,:,1:kdr+1)=vtempe(:,:,:,:,:,1:kdr+1) !my first additions -TW
   !vo(:,:,:,:,:,1:kdr+1)=vtempo(:,:,:,:,:,1:kdr+1) !my first additions -TW

   do i = 1,kdr+1 !added for testing TW
    ve(:,:,:,:,:,i)=vtempe(:,:,:,:,:,i) 
    vo(:,:,:,:,:,i)=vtempo(:,:,:,:,:,i)
   enddo ! end test do TW

   !re-orthogonalize the v_k+1 vector againist the first k vectors
   do i=1,kdr
    !call vecdot(vc(:,:,:,:,:,i),vc(:,:,:,:,:,kdr+1),tv,MRT2)  !ORIGINAL
    call vecdot(ve(:,:,:,:,:,i),ve(:,:,:,:,:,kdr+1),tv1,MRT2)
    call vecdot(vo(:,:,:,:,:,i),vo(:,:,:,:,:,kdr+1),tv2,MRT2)
    tv = tv1 + tv2
    do isite=1,nvhalf
     do icri=1,5,2
     ! vc(icri,isite,:,:,:,kdr+1) = vc(icri,isite,:,:,:,kdr+1) - & !ORIGINAL                    
     !   tv(1)*vc(icri,isite,:,:,:,i)+tv(2)*vc(icri+1,isite,:,:,:,i) !ORIGINAL            

     ! vc(icri+1,isite,:,:,:,kdr+1) = vc(icri+1,isite,:,:,:,kdr+1) - & !ORIGINAL                 
     !   tv(1)*vc(icri+1,isite,:,:,:,i)+tv(2)*vc(icri,isite,:,:,:,i) !ORIGINAL


      ve(icri,isite,:,:,:,kdr+1) = ve(icri,isite,:,:,:,kdr+1)  &
                                   -tv(1)*ve(icri,isite,:,:,:,i) &
                                   +tv(2)*ve(icri+1,isite,:,:,:,i)

      ve(icri+1,isite,:,:,:,kdr+1) = ve(icri+1,isite,:,:,:,kdr+1)  &
                                     -tv(1)*ve(icri+1,isite,:,:,:,i) &
                                     -tv(2)*ve(icri,isite,:,:,:,i) !changed from + to - TW                  

      vo(icri,isite,:,:,:,kdr+1) = vo(icri,isite,:,:,:,kdr+1)  &
                                  -tv(1)*vo(icri,isite,:,:,:,i) &
                                  +tv(2)*vo(icri+1,isite,:,:,:,i)

      vo(icri+1,isite,:,:,:,kdr+1) = vo(icri+1,isite,:,:,:,kdr+1)  &
                                     -tv(1)*vo(icri+1,isite,:,:,:,i) &
                                     -tv(2)*vo(icri,isite,:,:,:,i) !changed from + to - TW        
     enddo
    enddo
   enddo !i

   !call vecdot(vc(:,:,:,:,:,kdr+1),vc(:,:,:,:,:,kdr+1),tv,MRT2) !ORIGINAL
   !const=1.0_KR/sqrt(tv(1)) !ORIGINAL
   !vc(:,:,:,:,:,kdr+1)=const*vc(:,:,:,:,:,kdr+1) !ORIGINAL
   call vecdot(ve(:,:,:,:,:,kdr+1),ve(:,:,:,:,:,kdr+1),tv1,MRT2)
   call vecdot(vo(:,:,:,:,:,kdr+1),vo(:,:,:,:,:,kdr+1),tv2,MRT2)
   tv = tv1 + tv2
   tv(1) = sqrt(tv(1))
   tv(2) = 0.0_KR2
   const=1.0_KR/tv(1)
   do isite = 1,nvhalf !added by TW 8/11/18, see end of gmresdr5EIG
   !ve(:,:,:,:,:,kdr+1)=const*ve(:,:,:,:,:,kdr+1) !my first change, not in do statement           
   !vo(:,:,:,:,:,kdr+1)=const*vo(:,:,:,:,:,kdr+1) !my first change, not in do statement          
    ve(:,isite,:,:,:,kdr+1)=const*ve(:,isite,:,:,:,kdr+1)
    vo(:,isite,:,:,:,kdr+1)=const*vo(:,isite,:,:,:,kdr+1)
   enddo !isite 

   icycle=icycle+1
  enddo maindo

end subroutine minresdr5EIG
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine minresdr(rwdir,phi,x,ldr,resmax,itermin,itercount,u,GeeGooinv, &
                      iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                      lvbc,ib,lbd,iblv,hoption,MRT,MRT2)

! ldr(1)=mdr, ldr(2)=kdr
! MinRes-DR(mdr,kdr) that solves the non-hermetian problem M*x=phi by converting
! it into a hermetain problem. For Wilson fermions this can be done through
! multiplication with gamma5 and solving gamma5*M*x=gamma5*phi. For all cases
! we can use M^dagge*M*x=M^dagger*phi.
!
! INPUT:
!   rwdir is the read and write directory
!   phi() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   mdr is the maximum dimension of the subspace.
!   kdr is the number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by landr.

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI),intent(in), dimension(:)             :: ldr
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             hoption,myid, MRT,MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    !Local variables
    real(kind=KR), dimension(6,ntotal,4,2,8) :: dphi,h,r0,w
    real(kind=KR), dimension(201,201) :: orthog
    real(kind=KR), allocatable,dimension(:,:,:,:,:,:)::vc
    real(kind=KR),allocatable,dimension(:,:)::tc,tcs,tcs2
    real(kind=KR),allocatable,dimension(:)::cvec,cvec2,cvecs,srv
    real(kind=KR),allocatable,dimension(:)::ym
    real(kind=KR), dimension(2)::tv,dot
    real(kind=KR) ::normnum,normnum_hermetian,const,cr,f,g,cs,sn,f1,rn,rn2,tmp,&
                    srvnorm,residual_norm
    integer(kind=KI)::isite,icri,j,i,jj,ii,kk,idag,mdr,kdr,icycle,jstart,jmin
    real(kind=KR),allocatable,dimension(:,:)::lslhs,lsrhs
    real(kind=KR),allocatable,dimension(:)::work
    integer(kind=KI)::nrhs,lwork,info,lworke,exists
    character(len=1)::trans,jobvl,jobvr
    integer(kind=KI),allocatable,dimension(:) ::ipivl
    real(kind=KR),allocatable,dimension(:,:)  ::vec1
    real(kind=KR),allocatable,dimension(:,:)  ::tch
    real(kind=KR),allocatable,dimension(:) ::wr,wi,wv,tmpwr,tmpwi
    real(kind=KR),allocatable,dimension(:,:)::vl,vr,vro
    real(kind=KR),allocatable,dimension(:)  ::worke
    real(kind=KR),allocatable,dimension(:,:)::pk


    !Parameters
    mdr=ldr(1)
    kdr=ldr(2)
    nrhs=1
    lwork=10*(mdr+1)
    trans='n'
    jobvl='n' !don't compute left-eigenvectors
    jobvr='v' !compute right-eigenvectors
    lworke=10*mdr



    !Allocations
    allocate(vc(6,ntotal,4,2,8,mdr+1))
    allocate(tc(mdr+1,mdr))
    allocate(tcs(mdr+1,mdr))
    allocate(tcs2(mdr+1,mdr))
    allocate(cvec(mdr+1))
    allocate(cvec2(mdr+1))
    allocate(cvecs(mdr+1))
    allocate(srv(mdr+1))
    allocate(ym(mdr))
    allocate(lslhs(mdr+1,mdr))
    allocate(lsrhs(mdr+1,nrhs))
    allocate(work(lwork))
    allocate(ipivl(mdr))
    allocate(vec1(mdr,nrhs))
    allocate(tch(mdr,mdr))
    allocate(wr(mdr))
    allocate(wi(mdr))
    allocate(wv(mdr))
    allocate(tmpwr(mdr))
    allocate(tmpwi(mdr))
    allocate(vl(mdr,mdr))
    allocate(vr(mdr,mdr))
    allocate(vro(mdr,mdr))
    allocate(worke(lworke))
    allocate(pk(mdr+1,kdr+1))


    !Intitalizations
    dphi=0.0_KR
    h=0.0_KR
    r0=0.0_KR
    w=0.0_KR
    vc=0.0_KR
    tc=0.0_KR
    tcs=0.0_KR
    tcs2=0.0_KR
    cvec=0.0_KR
    cvec2=0.0_KR
    srv=0.0_KR
    ym=0.0_KR
    lslhs=0.0_KR
    lsrhs=0.0_KR
    vec1=0.0_KR
    tch=0.0_KR
    wr=0.0_KR
    wi=0.0_KR
    wv=0.0_KR
    vl=0.0_KR
    vr=0.0_KR
    vro=0.0_KR
    itercount=0
    pk=0.0_KR
    vtemp=0.0_KR
    tcnew=0.0_KR



    !Copy the right-hand side
    dphi=phi

    !Set the initial guess to zero
    x=0.0_KR

    !Compute the initial residual norm assuming zero initial guess
    call vecdot(phi,phi,tv,MRT2)
    normnum=sqrt(tv(1))

    !Compute the right hand side of the hermetian problem
    if(hoption==1) then
     call gam5x(dphi)
    elseif(hoption==2) then
     idag=1
     call Hdbletm(h,u,GeeGooinv,dphi,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     dphi= h
    endif

    !Compute the residual of the hermetian problem and vc_1
    do isite=1,nvhalf
     r0(:,isite,:,:,:)=dphi(:,isite,:,:,:)
    enddo !isite

    call vecdot(r0,r0,tv,MRT2)
    normnum_hermetian=sqrt(tv(1))
    const=1.0_KR/normnum_hermetian
    do isite=1,nvhalf
     vc(:,isite,:,:,:,1)=const*r0(:,isite,:,:,:)
    enddo !isite


    icycle=1


    maindo: do

    if(icycle==1) then
     jstart=1
    else
     jstart=kdr+1
    endif

    jdo: do j=jstart,mdr

     if(hoption==1) then
      idag=0
      call Hdbletm(w,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      call gam5x(w)
      itercount=itercount+1
     elseif(hoption==2) then
      idag=0
      call Hdbletm(h,u,GeeGooinv,vc(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      idag=1
      call Hdbletm(w,u,GeeGooinv,h,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      itercount=itercount+2
     endif !hoption==1

     if( ((icycle==1).and.(j==1)) .or. ((icycle>1).and.(j==kdr+1)) ) then
      jmin=1
     else
      jmin=j-1
     endif

     if(j.gt.1) then
      do i=jmin,j-1
       tc(i,j)=tc(j,i)
       do isite=1,nvhalf
        w(:,isite,:,:,:)=w(:,isite,:,:,:)-tc(i,j)*vc(:,isite,:,:,:,i)
       enddo !isite
      enddo
     endif !j.gt.1

     call vecdot(vc(:,:,:,:,:,j),w,tv,MRT2)
     tc(j,j)=tv(1)

     do isite=1,nvhalf
      w(:,isite,:,:,:)=w(:,isite,:,:,:)-tc(j,j)*vc(:,isite,:,:,:,j)
     enddo !isite

     !Re-Orthogonalization part againist eigenvectors
     if(icycle>1) then
      do i=1,kdr
       call vecdot(vc(:,:,:,:,:,i),w,tv,MRT2)

       do isite=1,nvhalf
        do icri=1,5,2
         w(icri,isite,:,:,:)=w(icri,isite,:,:,:)-tv(1)*vc(icri,isite,:,:,:,i)+&
                               tv(2)*vc(icri+1,isite,:,:,:,i)
         w(icri+1,isite,:,:,:)=w(icri+1,isite,:,:,:)-tv(1)*vc(icri+1,isite,:,:,:,i)-&
                               tv(2)*vc(icri,isite,:,:,:,i)
        enddo
       enddo

      enddo
     endif

     call vecdot(w,w,tv,MRT2)
     if(tv(1).lt.1.e-12) then
      mdr=j
      exit jdo
     endif
     tc(j+1,j)=sqrt(tv(1))
     const=1.0_KR/tc(j+1,j)
     do isite=1,nvhalf
      vc(:,isite,:,:,:,j+1)=const*w(:,isite,:,:,:)
     enddo !isite

    enddo jdo

    !-----inserting output for othog of basis vectors -TW 8/11/18
 if (.true.) then
  if (icycle==1 .OR. icycle==2) then !print out orthog. of basis vectors for the krylov subspace TW 2/18/18
   do i=1,mdr+1
    do kk=1,mdr+1

      !dote=0.0_KR2
      !doto=0.0_KR2
      dot = 0.0_KR

     !call vecdot(ve(:,:,:,:,:,i),ve(:,:,:,:,:,kk),dote,MRT2)
     !call vecdot(vo(:,:,:,:,:,i),vo(:,:,:,:,:,kk),doto,MRT2)
     call vecdot(vc(:,:,:,:,:,i),vc(:,:,:,:,:,kk),dot,MRT2)

     !orthog(i,kk)=dote(1) + doto(1)
     !orthog(i,kk)=dote(2) + doto(2)
     orthog(i,kk) = dot(1)

   enddo!j
  enddo!i


        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"orthog.dat", exist=exists)
            if (.not. exists) then

               open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="new",&
               action="write",form="formatted")
               close(unit=48,status="keep")
            endif
       endif


     do i=1,mdr+1
       do kk=1,mdr+1
            if (myid==0) then
                 open(unit=48,file=trim(rwdir(myid+1))//"orthog.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=48,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12)")"icycle=",icycle,"V'V(",i,",",kk,")=",orthog(i,kk)
                 close(unit=48,status="keep")
            endif
      enddo
    enddo



  endif !icycle
endif!true or false


    !Solve the least square problem using givens rotations:
    !------------------------------------------------------
    tcs=tc
    cvec=0.0_KR2

    !compute cvec=V_{m+1}^dagger*r0
    do i=1,mdr+1
     call vecdot(vc(:,:,:,:,:,i),r0,tv,MRT2)
     cvec(i)=tv(1)
    enddo !i

    lslhs=tc
    lsrhs(:,1)=cvec(:)
    cvecs=cvec

    ! WARNINIG: The following Givens rotations part need to be fixed to work for
    ! all cycles
    ! Its current form is valid only for the first cycle. For next cycles, note
    ! that tc has
    ! a full k+1 x k+1 leading matrix.

    ! Start of Givens rotations
    !do i=1,mdr
    ! tcs2=tcs
    ! cvec2=cvec
    ! f=tcs(i,i)
    ! g=tcs(i+1,i)
    ! call rotate2(f,g,cs,sn,f1)
    ! tcs(i,i)=f1
    ! tcs(i+1,i)=0.0_KR2
    ! do jj=1,mdr
    !  if(jj.ne.i) then
    !    tcs(i,jj)=cs*tcs2(i,jj)+sn*tcs2(i+1,jj)
    !    tcs(i+1,jj)=-sn*tcs2(i,jj)+cs*tcs2(i+1,jj)
    !  endif
    ! enddo
    ! cvec(i)=cs*cvec2(i)+sn*cvec2(i+1)
    ! cvec(i+1)= -sn*cvec2(i)+cs*cvec2(i+1)
    !enddo !i


    !ym(mdr)=cvec(mdr)/tcs(mdr,mdr)
    !do i=mdr-1,1,-1
    ! tmp=0.0_KR
    ! do jj=i+1,mdr
    !  tmp=tmp+tcs(i,jj)*ym(jj)
    ! enddo
    ! ym(i)=(cvec(i)-tmp)/tcs(i,i)
    !enddo

    ! End of Givens rotations


   !Solve the least-squares problem by using lapack subroutine dgels

   call dgels(trans,mdr+1,mdr,nrhs,lslhs,mdr+1,lsrhs,mdr+1,work,lwork,info)
   ym(1:mdr)=lsrhs(1:mdr,1)

   !Compute the components of the short residual vector
   do i=1,mdr+1
    tmp=0.0_KR
    do j=1,mdr
     tmp=tmp+tc(i,j)*ym(j)
    enddo !j
     srv(i)= cvecs(i)- tmp
   enddo !i

   srvnorm=0.0_KR
   do i=1,mdr+1
    srvnorm=srvnorm+srv(i)*srv(i)
   enddo

   srvnorm=sqrt(srvnorm)






   !Update the solution and compute the residual norm
   do i=1,mdr
    do isite=1,nvhalf
     x(:,isite,:,:,:)=x(:,isite,:,:,:)+ym(i)*vc(:,isite,:,:,:,i)
    enddo !isite
   enddo !i

   !Short cut
   rn=abs(lsrhs(mdr+1,1))
   !cr=rn

   !direct
   idag=0
   call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

   do isite=1,nvhalf
    w(:,isite,:,:,:)=phi(:,isite,:,:,:)-h(:,isite,:,:,:)
   enddo

   !Residual norm of the original problem
   call vecdot(w,w,tv,MRT2)
   residual_norm=sqrt(tv(1))

   if(hoption==1) then
    call gam5x(w)
    do isite=1,nvhalf
     r0(:,isite,:,:,:)=w(:,isite,:,:,:)
    enddo
    call vecdot(r0,r0,tv,MRT2)
    rn2=sqrt(tv(1))
   elseif(hoption==2) then
    idag=1
    call Hdbletm(r0,u,GeeGooinv,w,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    call vecdot(r0,r0,tv,MRT2)
    rn2=sqrt(tv(1))
   endif


   if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
    write(unit=8,fmt="(a40,i10,4d20.10)")"minres-dr,itercount,rn,rn2,srvnorm,residual",&
         itercount,rn/normnum_hermetian,rn2/normnum_hermetian,srvnorm/normnum_hermetian,residual_norm/normnum
    close(unit=8,status="keep")
   endif

   if( (residual_norm/normnum.LT. resmax) .OR. (itercount > 5000) ) then

     if(myid==0) then
      open(unit=10,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",action="write",&
      form="formatted",status="old",position="append")
      write(unit=10,fmt='(a30,2d17.10)')"#kappa,mu=",kappa(1),kappa(2)
      write(unit=10,fmt='(4a20)')"#eigen num","wr","wi","eig-res-norm"
     endif
     !printout the harmonic ritz vectors and their residual
     do i=1,kdr
      w=0.0_KR
      do j=1,kdr+1
       do isite=1,nvhalf
        w(:,isite,:,:,:)=w(:,isite,:,:,:)+vc(:,isite,:,:,:,j)*tc(j,i)
       enddo
      enddo
      do isite=1,nvhalf
       w(:,isite,:,:,:)=w(:,isite,:,:,:)-wr(i)*vc(:,isite,:,:,:,i)
      enddo
      call vecdot(w,w,tv,MRT2)
      if(myid==0) then
       write(unit=10,fmt='(i20,3d20.10)')i,wr(i),wi(i),sqrt(tv(1))
      endif
     enddo !i

     if(myid==0) close(unit=10,status="keep")

     !Dellocations
     deallocate(vc)
     deallocate(tc)
     deallocate(tcs)
     deallocate(tcs2)
     deallocate(cvec)
     deallocate(cvec2)
     deallocate(cvecs)
     deallocate(srv)
     deallocate(ym)
     deallocate(lslhs)
     deallocate(lsrhs)
     deallocate(work)
     deallocate(ipivl)
     deallocate(vec1)
     deallocate(tch)
     deallocate(wr)
     deallocate(wi)
     deallocate(wv)
     deallocate(tmpwr)
     deallocate(tmpwi)
     deallocate(vl)
     deallocate(vr)
     deallocate(vro)
     deallocate(worke)
     deallocate(pk)
     exit maindo

    endif







   !Compute the k smallest harmonic Ritz values and vectors:

   !Compute H_m^{-T}*e_m
   do i=1,mdr
    do j=1,mdr
     tch(i,j)=tc(j,i) !tc is real
    enddo
   enddo

   vec1=0.0_KR
   vec1(mdr,1)=1.0_KR

   call DGESV(mdr,nrhs,tch,mdr,ipivl,vec1,mdr,info)

   !build the matrix for which to compute the eigenvectors
   tch=0.0_KR

   do i=1,mdr
    tch(i,mdr)=tc(i,mdr)+tc(mdr+1,mdr)**2*vec1(i,1)
   enddo !i

   do i=1,mdr
    do j=1,mdr-1
     tch(i,j)=tc(i,j)
    enddo
   enddo




   call DGEEV(jobvl,jobvr,mdr,tch,mdr,wr,wi,vl,mdr,vr,mdr,worke,lworke,info)


   !Sort eigenvalues in ascending order according to their absolute values
   do i=1,mdr
    tmpwr(i)=sqrt(wr(i)**2+wi(i)**2)
   enddo

   ipivl=0
   call indexx(mdr,tmpwr,ipivl)

   do i=1,mdr
    tmpwr(i)=wr(i)
    tmpwi(i)=wi(i)
   enddo

   vro=vr
   do i=1,mdr
    wr(i)=tmpwr(ipivl(i))
    wi(i)=tmpwi(ipivl(i))
    vr(:,i)=vro(:,ipivl(i))
   enddo

   do i=1,kdr
    heigen(i)=wr(i)
   enddo

   !Building orthonormal basis for the kdr smallest eigenvectors:
   do j=1,kdr
    wv(:)=vr(:,j)
    do i=1,j-1
     tmp=0.0_KR
     do ii=1,mdr
      tmp=tmp+vr(ii,i)*wv(ii)
     enddo
     wv(:)=wv(:)-tmp*vr(:,i)
    enddo !i
    tmp=0.0_KR
    do ii=1,mdr
     tmp=tmp+wv(ii)*wv(ii)
    enddo
    const=1.0_KR/sqrt(tmp)
    vr(:,j)=const*wv(:)
   enddo !j




   !Build the matrix pk(mdr+1,kdr+1)
   pk=0.0_KR
   do i=1,kdr
    do j=1,mdr
     pk(j,i)=vr(j,i)
    enddo
    pk(mdr+1,i)=0.0_KR
   enddo !i

   cvec2=0.0_KR
   do i=1,mdr+1
    tmp=0.0_KR
    do jj=1,mdr
     tmp=tmp+tc(i,jj)*ym(jj)
    enddo
   cvec2(i)=cvecs(i)-tmp
   enddo

   do i=1,kdr
    tmp=0.0_KR
    do j=1,mdr+1
     tmp=tmp+pk(j,i)*cvec2(j)
    enddo
    do j=1,mdr+1
     cvec2(j)=cvec2(j)-tmp*pk(j,i)
    enddo
   enddo

   tmp=0.0_KR
   do j=1,mdr+1
    tmp=tmp+cvec2(j)*cvec2(j)
   enddo
   const=1.0_KR/sqrt(tmp)
   cvec2(:)=const*cvec2(:)

   do j=1,mdr+1
    pk(j,kdr+1)=cvec2(j)
   enddo

   !Build V_(k+1)^(new) and tcbar_k^new
   tcs2=0.0_KR

   do i=1,mdr+1
    do j=1,kdr
     tmp=0.0_KR
     do ii=1,mdr
      tmp=tmp+tc(i,ii)*pk(ii,j)
     enddo
     tcs2(i,j)=tmp
    enddo
   enddo

   do i=1,kdr+1
    do j=1,kdr
     tcnew(i,j)=0.0_KR
     do ii=1,mdr+1
      tcnew(i,j) = tcnew(i,j) + pk(ii,i)*tcs2(ii,j)
     enddo
    enddo
   enddo

   do i=1,kdr+1
    do j=1,kdr
     tc(i,j)=tcnew(i,j)
    enddo
   enddo


   vtemp=0.0_KR
   do j=1,kdr+1
    vtemp(:,:,:,:,:,j)=0.0_KR
    do i=1,mdr+1
     vtemp(:,:,:,:,:,j) = vtemp(:,:,:,:,:,j) + vc(:,:,:,:,:,i)*pk(i,j)
    enddo
   enddo

   vc(:,:,:,:,:,1:kdr+1)=vtemp(:,:,:,:,:,1:kdr+1)


   !re-orthogonalize the v_k+1 vector againist the first k vectors
   do i=1,kdr
    call vecdot(vc(:,:,:,:,:,i),vc(:,:,:,:,:,kdr+1),tv,MRT2)
    do isite=1,nvhalf
     do icri=1,5,2
      vc(icri,isite,:,:,:,kdr+1) = vc(icri,isite,:,:,:,kdr+1) - &
        tv(1)*vc(icri,isite,:,:,:,i)+tv(2)*vc(icri+1,isite,:,:,:,i)

      vc(icri+1,isite,:,:,:,kdr+1) = vc(icri+1,isite,:,:,:,kdr+1) - &
        tv(1)*vc(icri+1,isite,:,:,:,i)-tv(2)*vc(icri,isite,:,:,:,i)
     enddo
    enddo
   enddo !i

   call vecdot(vc(:,:,:,:,:,kdr+1),vc(:,:,:,:,:,kdr+1),tv,MRT2)
   const=1.0_KR/sqrt(tv(1))
   vc(:,:,:,:,:,kdr+1)=const*vc(:,:,:,:,:,kdr+1)


   icycle=icycle+1
  enddo maindo

end subroutine minresdr
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -










 
 subroutine LUinversion(a,n,ainv,counter)
! Construct the inverse of the square complex matrix "a".
! INPUT:
!   "a" is the square complex matrix to be inverted.
!   expected size: a(2,n,n)
!   "n" is the number of rows (and columns) in the matrix "a".
! OUTPUT:
!   "ainv" is the desired inverse of the matrix "a".
!   "counter" is the number of extra LU inversions required.

    real(kind=KR),    intent(in),  dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                    :: n
    real(kind=KR),    intent(out), dimension(:,:,:) :: ainv
    integer(kind=KI), intent(inout)                 :: counter

    integer(kind=KI)                         :: i, j, iri
    real(kind=KR)                            :: rowsign
    integer(kind=KI), parameter              :: nmax=6
    integer(kind=KI), dimension(nmax)        :: indx
    real(kind=KR),    dimension(2,nmax,nmax) :: alud, unity

! Set unity to the identity.
    unity = 0.0_KR
    do i = 1,n
     unity(1,i,i) = 1.0_KR
     ainv(1,i,i) = 1.0_KR
    enddo ! i

! Copy "a" into alud, since ludcmp will destroy its input matrix.
    do j = 1,n
     do i = 1,n
      do iri = 1,2
       alud(iri,i,j) = a(iri,i,j)
      enddo ! iri
     enddo ! i
    enddo ! j

! Decompose the matrix alud into LU form.
    call ludcmp(alud,n,indx,rowsign)

! Find the inverse, column by column.
    do j = 1,n
     call lubksb(alud,n,indx,ainv(:,:,j))
     call mprove(a,alud,n,indx,unity(:,:,j),ainv(:,:,j),counter)
    enddo ! j

 end subroutine LUinversion

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine ludcmp(a,n,indx,d)
! Adapted from Numerical Recipes, pages 38 and 39.
! Given a matrix a(ReIm,1:n,1:n) with physical dimension np by np, this routine
! replaces it by the LU decomposition of a rowwise permutation of itself.
! "a" and "n" are input.  "a" is output, arranged as in equation (2.3.14);
! indx(1:n) is an output vector that records the row permutation effected
! by the partial pivoting; "d" is output as +/-1 depending on whether the
! number of row interchanges was even or odd, respectively.  This routine
! is used in combination with lubksb to solve linear equations or invert a
! matrix.

    real(kind=KR),    intent(inout), dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n
    integer(kind=KI), intent(out),   dimension(:)     :: indx
    real(kind=KR),    intent(out)                     :: d

    integer(kind=KI)                  :: i, j, k, imax
    real(kind=KR)                     :: aatry, aamax, ajjsq
    real(kind=KR),    dimension(2)    :: csum, dum
    integer(kind=KI), parameter       :: nmax=max(6,nmaxGMRES)
    real(kind=KR),    dimension(nmax) :: vv

! No row interchanges yet.
    d = 1.0_KR

! Loop over rows to get the implicit scaling information.
    do i = 1,n
     aamax = 0.0_KR
     do j = 1,n
      aatry = sqrt( a(1,i,j)**2 + a(2,i,j)**2 )
      if (aatry>aamax) aamax=aatry
     enddo ! j
     if (aamax==0.0_KR) then
      open(unit=8,file="INVERTERS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine ludcmp: aamax =", aamax
      close(unit=8,status="keep")
     endif
     vv(i) = 1.0_KR/aamax
    enddo ! i

! Loop over columns of Crout's method.
    do j = 1,n
!-this is equation (2.3.12) except for i=j.
     do i = 1,j-1
      csum(1) = a(1,i,j)
      csum(2) = a(2,i,j)
      do k = 1,i-1
       csum(1) = csum(1) - a(1,i,k)*a(1,k,j) + a(2,i,k)*a(2,k,j)
       csum(2) = csum(2) - a(1,i,k)*a(2,k,j) - a(2,i,k)*a(1,k,j)
      enddo ! k
      a(1,i,j) = csum(1)
      a(2,i,j) = csum(2)
     enddo ! i
!-initialize for the search for largest pivot element.
     aamax = 0.0_KR
!-this is i=j of equation (2.3.12) and i=j+1,...,N of equation (2.3.13).
     do i = j,n
      csum(1) = a(1,i,j)
      csum(2) = a(2,i,j)
      do k = 1,j-1
       csum(1) = csum(1) - a(1,i,k)*a(1,k,j) + a(2,i,k)*a(2,k,j)
       csum(2) = csum(2) - a(1,i,k)*a(2,k,j) - a(2,i,k)*a(1,k,j)
      enddo ! k
      a(1,i,j) = csum(1)
      a(2,i,j) = csum(2)
!-figure of merit for the pivot.
      dum(1) = vv(i)*sqrt(csum(1)**2+csum(2)**2)
!-is it better than the best so far?
      if (dum(1)>=aamax) then
       imax = i
       aamax = dum(1)
      endif
     enddo ! i
!-do we need to interchange rows?
     if (j/=imax) then
!-if yes, then do so...
      do k = 1,n
       dum(1) = a(1,imax,k)
       dum(2) = a(2,imax,k)
       a(1,imax,k) = a(1,j,k)
       a(2,imax,k) = a(2,j,k)
       a(1,j,k) = dum(1)
       a(2,j,k) = dum(2)
      enddo ! k
!-...and change the parity of d.
      d = -d
!-also interchange the scale factor.
      vv(imax) = vv(j)
     endif
     indx(j) = imax
!-if the pivot element is zero the matrix is singular (at least to the
! precision of the algorithm).  For some applications on singular matrices,
! it is desirable to substitute TINY for zero.
     ajjsq = a(1,j,j)**2 + a(2,j,j)**2
     if (ajjsq==0.0_KR) then
      open(unit=8,file="INVERTERS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine ludcmp: ajjsq =", ajjsq
      close(unit=8,status="keep")
     endif
!-now, finally, divide by the pivot element.
     if (j/=n) then
      dum(1) = a(1,j,j)/ajjsq
      dum(2) = -a(2,j,j)/ajjsq
      do i = j+1,n
       aatry = a(1,i,j)*dum(1) - a(2,i,j)*dum(2)
       a(2,i,j) = a(1,i,j)*dum(2) + a(2,i,j)*dum(1)
       a(1,i,j) = aatry
      enddo ! i
     endif
    enddo ! j

 end subroutine ludcmp

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine lubksb(a,n,indx,b)
! Adapted from Numerical Recipes, page 39.
! Solves the set of n linear equations A.X=B.  Here "a" is input, not as the
! matrix A but rather as its LU decomposition, determined by the routine
! ludcmp.  indx is input as the permutation vector returned by ludcmp.
! b(ReIm,1:n) is input as the right-hand side vector B, and returns with the
! solution vector X.  "a", n, np, and indx are not modified by this routine
! and can be left in place for successive calls with different right-hand
! sides b.  This routine takes into account the possibility that b will begin
! with many zero elements, so it is efficient for use in matrix inversion.

    real(kind=KR),    intent(in),    dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n
    integer(kind=KI), intent(in),    dimension(:)     :: indx
    real(kind=KR),    intent(inout), dimension(:,:)   :: b

    integer(kind=KI)            :: i, j, ii, ll
    real(kind=KR), dimension(2) :: csum
    real(kind=KR)               :: csummag, aiisqinv

! When ii is set to a positive value, it will become the index of the first
! nonvanishing element of b.  We now do the forward substitution, equation
! (2.3.6).  The only new wrinkle is to unscramble the permutation as we go.
    ii = 0
    do i = 1,n
     ll = indx(i)
     csum(1) = b(1,ll)
     csum(2) = b(2,ll)
     csummag = sqrt(csum(1)**2+csum(2)**2)
     b(:,ll) = b(:,i)
     if (ii/=0) then
      do j = ii,i-1
       csum(1) = csum(1) - a(1,i,j)*b(1,j) + a(2,i,j)*b(2,j)
       csum(2) = csum(2) - a(1,i,j)*b(2,j) - a(2,i,j)*b(1,j)
      enddo ! j
     elseif (csummag/=0.0_KR) then
!-a nonzero element was encountered, so from now on we will have to do
! the sums in the loop above.
      ii = i
     endif
     b(1,i) = csum(1)
     b(2,i) = csum(2)
    enddo ! i

! Now we do the backsubstitution, equation (2.3.7).
    do i = n,1,-1
     csum(1) = b(1,i)
     csum(2) = b(2,i)
     do j = i+1,n
      csum(1) = csum(1) - a(1,i,j)*b(1,j) + a(2,i,j)*b(2,j)
      csum(2) = csum(2) - a(1,i,j)*b(2,j) - a(2,i,j)*b(1,j)
     enddo ! j
     aiisqinv = 1.0_KR/(a(1,i,i)**2+a(2,i,i)**2)
     b(1,i) = ( csum(1)*a(1,i,i) + csum(2)*a(2,i,i) )*aiisqinv
     b(2,i) = ( csum(2)*a(1,i,i) - csum(1)*a(2,i,i) )*aiisqinv
    enddo ! i

 end subroutine lubksb

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mprove(a,alud,n,indx,b,x,counter)
! Adapted from Numerical Recipes, page 48.
! Improves a solution vector x(ReIm,1:n) of the linear set of equations A.X=B.
! The matrix a(ReIm,1:n,1:n) and the vectors b(ReIm,1:n) and x(ReIm,1:n) are
! input, as is the dimension n.  Also input is alud, the LU decomposition
! of "a" as returned by ludcmp, and the vector indx also returned by that
! routine.  On output, only x(ReIm,1:n) is modified, to an improved set of
! values.

    real(kind=KR),    intent(in),    dimension(:,:,:) :: a, alud
    integer(kind=KI), intent(in)                      :: n
    integer(kind=KI), intent(in),    dimension(:)     :: indx
    real(kind=KR),    intent(in),    dimension(:,:)   :: b
    real(kind=KR),    intent(inout), dimension(:,:)   :: x
    integer(kind=KI), intent(inout)                   :: counter

    integer(kind=KI)                    :: i, j
    integer(kind=KI), parameter         :: nmax=max(6,nmaxGMRES)
    real(kind=KR),    dimension(2,nmax) :: r
    real(kind=KR)                       :: rmag
    real(kind=KR),    parameter         :: rmax=1.0e-18_KR
    real(kind=KR2),   dimension(2)      :: sdp

! Iterate until desired accuracy is obtained.
    loopy: do

!-calculate the right-hand side, accumulating the residual in double precision.
     rmag = 0.0_KR
     do i = 1,n
      sdp(1) = -b(1,i)
      sdp(2) = -b(2,i)
      do j = 1,n
       sdp(1) = sdp(1) + real(a(1,i,j),KR2)*real(x(1,j),KR2) &
                       - real(a(2,i,j),KR2)*real(x(2,j),KR2)
       sdp(2) = sdp(2) + real(a(1,i,j),KR2)*real(x(2,j),KR2) &
                       + real(a(2,i,j),KR2)*real(x(1,j),KR2)
      enddo ! j
      r(:,i) = sdp(:)
      rmag = rmag + r(1,i)**2 + r(2,i)**2
     enddo ! i

!-if the residual is small enough then exit the subroutine, otherwise
! solve for the error term and subtract it from the old solution.
     if (rmag<rmax) then
      exit loopy
     else
      counter = counter + 1
      call lubksb(alud,n,indx,r)
      do i = 1,n
       x(1,i) = x(1,i) - r(1,i)
       x(2,i) = x(2,i) - r(2,i)
      enddo ! i
     endif

    enddo loopy

 end subroutine mprove

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 





!BS  Subroutine for matrix multiplication is formed 



  subroutine vecdothigherrank(a,b,adb,m,MRT2)
! Calculate the dot product of two vectors, a and b, where the dot product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.

!    real(kind=KR),    intent(in),  dimension(:,:,:,:,:,:) :: a, b
 !   integer(kind=KI), intent(in)                          :: MRT2

    real(kind=KR),    intent(in),  dimension(:,:,:,:,:,:) :: a, b
    real(kind=KR2),   intent(out), dimension(2,m+1,m+1)         :: adb
    integer(kind=KI), intent(in)                            :: m
    integer(kind=KI), intent(in)                          :: MRT2
    integer(kind=KI)                         ::  i,j 
    real(kind=KR2),     allocatable,  dimension(:)         :: broda
    real(kind=KR2),     allocatable,  dimension(:,:,:,:,:)         :: aunt
    real(kind=KR2),     allocatable,  dimension(:,:,:,:,:)         :: munt
     




  do i=1,m
   do j=1,m
      broda=0.0_KR2
      aunt=0.0_KR2
      munt=0.0_KR2
       aunt(:,:,:,:,:)=a(:,:,:,:,:,i)
       munt(:,:,:,:,:)=b(:,:,:,:,:,j)
     print *,"gmresdrprinttest-insidehigherrank-before"
     call vecdot(a(:,:,:,:,:,i),b(:,:,:,:,:,j),broda,MRT2)
     if ((i==5).and.(j==4)) then

        print *,"gmresdrprinttest-insidehigherrank-later"
         print *,"broda1and2=",broda(1),broda(2)

     endif
     adb(1,i,j)=broda(1)
     adb(2,i,j)=broda(2)
   print *,"adb(1,5,4)= ",adb(1,5,4),adb(2,5,4)

   enddo
  enddo


     print *,"gmresdrprinttest-insidehigherrank2"






 end subroutine vecdothigherrank





!BS end subroutine vecdothigherrank


!BS  subroutine matrixmultiplylike
  subroutine matrixmultiplylike(a,h,m,matmul,MRT2)
! Calculate the dot product of two vectors, a and b, where the dot product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.

    
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:,:) :: a
    real(kind=KR),    intent(in),  dimension(:,:,:)       :: h
    real(kind=KR),    intent(out),  dimension(:,:,:,:,:,:) :: matmul
    integer(kind=KI), intent(in)                          :: MRT2,m

    integer(kind=KI)             :: isite, id, ieo, ibleo, ierr,j,Aham,icri
    real(kind=KR2), dimension(2) :: adbbit
     real(kind=KR2)               :: temp1,temp2
    do ibleo = 1,8
     do ieo = 1,2
      do id = 1,4
       do isite = 1,nvhalf
        do icri =1,5,2
         do j = 1,m
          do Aham = 1,m+1

            matmul(icri,isite,id,ieo,ibleo,j)   =temp1             &
                                                + a(icri,isite,id,ieo,ibleo,Aham)*h(1,Aham,j)   &
                                                - a(icri+1,isite,id,ieo,ibleo,Aham)*h(2,Aham,j) 
            matmul(icri+1,isite,id,ieo,ibleo,j) = temp2           &
                                                + a(icri,isite,id,ieo,ibleo,Aham)*h(2,Aham,j)   &
                                                + a(icri+1,isite,id,ieo,ibleo,Aham)*h(1,Aham,j)


            temp1=matmul(icri,isite,id,ieo,ibleo,j)
            temp2=matmul(icri+1,isite,id,ieo,ibleo,j)
 

 
          enddo !Aham
         enddo !j
        enddo !icri
       enddo  !isite 
      enddo !id
     enddo  !ieo
    enddo   !ibleo





! Sum the contributions from all processes.
    if (nps==1) then
    print*,"matrixmultiplylike is working" 
    else
print *,"Error in subroutine matrrixmultiplylike , nps not equal to one"

!     call MPI_REDUCE(adbbit(1),adb(1),2,MRT2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 !    call MPI_BCAST(adb(1),2,MRT2,0,MPI_COMM_WORLD,ierr)
    endif

 end subroutine matrixmultiplylike

!BS  end subroutine matrixmultiplylike


















  subroutine vecdot(a,b,adb,MRT2)
! Calculate the dot product of two vectors, a and b, where the dot product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.

    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: a, b
    real(kind=KR2),   intent(out), dimension(:)         :: adb
    integer(kind=KI), intent(in)                        :: MRT2

    integer(kind=KI)             :: isite, id, ieo, ibleo, ierr
    real(kind=KR2), dimension(2) :: adbbit

    adbbit = 0.0_KR2
    do ibleo = 1,8
     do ieo = 1,2
      do id = 1,4
       do isite = 1,nvhalf
        adbbit(1) = adbbit(1) & 
                  + real(a(1,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(2,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(4,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo) &
                       + a(6,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo),KR2)
        adbbit(2) = adbbit(2) &
                  + real(a(1,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       - a(2,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       - a(4,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo) &
                       - a(6,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo),KR2)
       enddo ! isite
      enddo ! id
     enddo ! ieo
    enddo ! ibleo

! Sum the contributions from all processes.
    if (nps==1) then
     adb = adbbit
    else
     call MPI_REDUCE(adbbit(1),adb(1),2,MRT2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(adb(1),2,MRT2,0,MPI_COMM_WORLD,ierr)
    endif

 end subroutine vecdot


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine vecdot1(a,b,adb,MRT2)
! Calculate the dot product of two vectors, a and b, where the dot product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.

    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: a, b
    real(kind=KR2),   intent(out), dimension(:)         :: adb
    integer(kind=KI), intent(in)                        :: MRT2

    integer(kind=KI)             :: isite, id, ieo, ibleo, ierr
    real(kind=KR2), dimension(2) :: adbbit

    adbbit = 0.0_KR2
    do ibleo = 1,8
     do ieo = 1,2
      do id = 1,4
       do isite = 1,ntotal
        adbbit(1) = adbbit(1) & 
                  + real(a(1,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(2,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(4,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo) &
                       + a(6,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo),KR2)
        adbbit(2) = adbbit(2) &
                  + real(a(1,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       - a(2,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       - a(4,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo) &
                       - a(6,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo),KR2)
       enddo ! isite
      enddo ! id
     enddo ! ieo
    enddo ! ibleo

! Sum the contributions from all processes.
    if (nps==1) then
     adb = adbbit
    else
     call MPI_REDUCE(adbbit(1),adb(1),2,MRT2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(adb(1),2,MRT2,0,MPI_COMM_WORLD,ierr)
    endif

 end subroutine vecdot1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine vecdotmyid(a,b,adb,MRT2,myid)
! Calculate the dot product of two vectors, a and b, where the dot product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.

    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: a, b
    real(kind=KR2),   intent(out), dimension(:)         :: adb
    integer(kind=KI), intent(in)                        :: MRT2,myid

    integer(kind=KI)             :: isite, id, ieo, ibleo, ierr
    real(kind=KR2), dimension(2) :: adbbit

    adbbit = 0.0_KR2
    do ibleo = 1,8
     do ieo = 1,2
      do id = 1,4
       do isite = 1,nvhalf
        adbbit(1) = adbbit(1) & 
                  + real(a(1,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(2,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(4,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo) &
                       + a(6,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo),KR2)
        adbbit(2) = adbbit(2) &
                  + real(a(1,isite,id,ieo,ibleo)*b(2,isite,id,ieo,ibleo) &
                       - a(2,isite,id,ieo,ibleo)*b(1,isite,id,ieo,ibleo) &
                       + a(3,isite,id,ieo,ibleo)*b(4,isite,id,ieo,ibleo) &
                       - a(4,isite,id,ieo,ibleo)*b(3,isite,id,ieo,ibleo) &
                       + a(5,isite,id,ieo,ibleo)*b(6,isite,id,ieo,ibleo) &
                       - a(6,isite,id,ieo,ibleo)*b(5,isite,id,ieo,ibleo),KR2)
       enddo ! isite
      enddo ! id
     enddo ! ieo
    enddo ! ibleo

! Sum the contributions from all processes.
    if (nps==1) then
     adb = adbbit
    else
     print *,"before MPI REDUCE, Process :",myid
     call MPI_REDUCE(adbbit(1),adb(1),2,MRT2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     print *,"after MPI REDUCE, process :",myid
     call MPI_BCAST(adb(1),2,MRT2,0,MPI_COMM_WORLD,ierr)
     print *,"after MPI BCAST, process :",myid
    endif

 end subroutine vecdotmyid

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   

subroutine real2complex_vec(r, m, z)

    real(kind=KR),    intent(in),  dimension(:,:)       :: r
    complex(kind=KCC), intent(out), dimension(:)         :: z
    integer(kind=KI), intent(in)                        :: m

    integer(kind=KI)                                    :: i

    do i=1,m
       z(i) = DCMPLX(r(1,i), r(2,i))
    enddo ! i

end subroutine real2complex_vec

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine real2complex_mat(r, m, n, z)

    real(kind=KR),    intent(in),  dimension(:,:,:)     :: r
    complex(kind=KCC), intent(out), dimension(:,:)       :: z
    integer(kind=KI), intent(in)                        :: m, n

    integer(kind=KI)                                    :: i, j

    do i=1,m
       do j=1,n
          z(i,j) = DCMPLX(r(1,i,j), r(2,i,j))
       enddo ! j
    enddo ! i

end subroutine real2complex_mat

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine complex2real_vec(z, m, r)

    real(kind=KR),    intent(out),  dimension(:,:)      :: r
    complex(kind=KCC), intent(in), dimension(:)          :: z
    integer(kind=KI), intent(in)                        :: m

    integer(kind=KI)                                    :: i

    do i=1,m
       r(1,i) = REAL(z(i))
       r(2,i) = AIMAG(z(i))
    enddo ! i

end subroutine complex2real_vec 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine complex2real_mat(z, m, n, r)

    real(kind=KR),    intent(out),  dimension(:,:,:)    :: r
    complex(kind=KCC), intent(in), dimension(:,:)        :: z
    integer(kind=KI), intent(in)                        :: m, n

    integer(kind=KI)                                    :: i, j

    do i=1,m
       do j=1,n
          r(1,i,j) = REAL(z(i,j))
          r(2,i,j) = AIMAG(z(i,j))
       enddo ! j
    enddo ! i

end subroutine complex2real_mat

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine mult_complex(ar, ai, br, bi, cr, ci)

    real(kind=KR),     intent(in)     :: ar, ai, br, bi
    real(kind=KR),     intent(out)    :: cr, ci

    complex(kind=KCC)                  :: za, zb, zc

    za = DCMPLX(ar, ai)
    zb = DCMPLX(br, bi)

    zc = za * zb

    cr = REAL(zc)
    ci = AIMAG(zc)

end subroutine mult_complex

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine printArray(a)

  complex(kind=KCC), intent(in), dimension(:,:)    :: a
  integer :: i, j

  print "(31e8.3)", (a(i,:),i=1,31)

end subroutine printArray

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine checkNonZero2(x,n)

   real(kind=KR2), intent(in), dimension(:,:,:,:,:) :: x
   integer, intent(in) :: n
 ! real(kind=KR2), intent(in) :: lg

   integer :: icolorir,isite,idirac,ioe,iblock

   do icolorir=1,6
      do isite=1,n
         do idirac=1,4
            do ioe=1,2
               do iblock=1,8
                  if ((abs(x(icolorir,isite,idirac,ioe,iblock)) /= 0.0_KR) .or.(abs(x(icolorir+1,isite,idirac,ioe,iblock)) /= 0.0_KR)) then
                     print "(a,5i4,a3,e15.7)", "icolorr,isite,idirac,ioe,iblock= ", &
                                           icolorir,isite,idirac,ioe,iblock, &
                                           " x= ",x(icolorir,isite,idirac,ioe,&
iblock)
                     print "(a,5i4,a3,e15.7)", "icolorr+1,isite,idirac,ioe,iblock= ", &
                                           icolorir,isite,idirac,ioe,iblock, &
                                           " x= ",x(icolorir+1,isite,idirac,ioe,iblock)
                   endif
               enddo ! iblock
            enddo ! ioe
         enddo ! idirac
      enddo ! isite
   enddo ! icolorir

end subroutine checkNonZero2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine checkNonZero(x,n,ib,ieo,iss,id,ic,site,icolorr)

   real(kind=KR2), intent(in), dimension(:,:,:,:,:) :: x
   integer, intent(in) :: n, ib,iss,ieo,id,ic,site,icolorr

   integer :: icolorir,isite,idirac,ie,iblock,icolorc,irow,icol,is

   icol = icolorr + 3*(iss -1) + 24*(id -1) + 96*(ieo -1) + 192*(ib - 1) 
    do iblock=1,8
      do ie =1,2
        do isite=1,n
          do idirac=1,4
             do icolorir=1,5,2
                ! is = isite + n*(iblock -1)
                  icolorc = icolorir/2 + 1
                  irow =  icolorc + 3*(isite - 1) + 24*(idirac -1) + 96*(ie -1) + 192*(iblock -1)
                 !icol =  icolorc + 3*(idirac -1) + 12*(isite - 1) + 96*(ie -1) + 192*(iblock -1)
                  if ((abs(x(icolorir,isite,idirac,ie,iblock)) /= 0.0).or. (abs(x(icolorir+1,isite,idirac,ie,iblock))/= 0.0)) then
                        print *, irow, icol, x(icolorir,isite,idirac,ie,iblock), x(icolorir+1,isite,idirac,ie,iblock)
                 endif
               enddo ! icolorir
           enddo ! idirac
        enddo ! isite
      enddo ! ie
   enddo ! iblock

end subroutine checkNonZero

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

subroutine printmatrix(x,n,isite,iblock,icolorir,idirac,site)

   real(kind=KR2), intent(in), dimension(:,:,:,:,:) :: x
   integer(kind=KI), intent(inout)                  :: isite,iblock,icolorir,idirac
   integer, intent(in)                              :: n,site

   integer                                          :: i,j,k,l,m,ioe
   real(kind=KR)                                    :: rhoe,ihoe
   
  
   ioe = 1

   print "(a,i4,a,5i4)", "printmatrix: site=",site,"icolorir,isite,idirac,ioe,iblock", icolorir,isite,idirac,ioe,iblock
  
!   ib = (site - isite)/n - 1
!   do m=1,ib
!     do i=1,5,2
!        do k=1,4
!           if (abs(x(i,j,k,2,iblock)) /= 0.0) then
!           site = (m-1)*n + j
           
          
            rhoe = x(icolorir  ,isite,idirac,ioe,iblock)
            ihoe = x(icolorir+1,isite,idirac,ioe,iblock)
       
!            if ((abs(rhoe) /= 0.0_KR) .or. (abs(ihoe) /= 0.0_KR)) then
               if (icolorir==1) then
                  print *, "site, color, dirac =", site, 1, idirac
                  print *, "REAL(hoe) , IMAG(hoe) =", rhoe, ihoe
               elseif(icolorir==3) then
                  print *, "site, color, dirac =", site, 2, idirac
                  print *, "REAL(hoe) , IMAG(hoe) =", rhoe, ihoe
               elseif(icolorir==5) then
                  print *, "site, color, dirac =", site, 3, idirac
                  print *, "REAL(hoe) , IMAG(hoe) =", rhoe, ihoe
               endif !(i==1)
!            endif ! abs()
!        enddo ! k
!     enddo ! k
!  enddo ! m 

end subroutine printmatrix

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

!subroutine matrix(x,isig,u,GeeGooinv,idag,coact,kappa,iflag,bc,vecbl, &
!                 vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

! This subroutine is strickly for debugging. We use it to pass the matrix
! and the RHS to a matlab program that the above alogrotihums
! were developed from.

! If isig =1 -> print matrix
!    isig =2 -> print rhs

!   real(kind=KR2), intent(in), dimension(:,:,:,:,:)      :: x
!   integer(kind=KI), intent(in)                          :: isig
!   real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
!   real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
!   real(kind=KR),    intent(in),    dimension(:)         :: kappa
!   integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
!   integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
!   integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
!   logical,          intent(in),    dimension(:)         :: ldiv
!   integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
 !  integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 !  integer(kind=KI), intent(in)                          :: idag, iflag 
!   integer(kind=KI), intent(in)                          :: myid, MRT
 !  logical,          intent(in),    dimension(:,:)       :: lbd
!
!  
!   real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp
!   real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
!   integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
!
!
! if (isig=1) then
!
!   htemp = 0.0_KR
!   site = 0.0_KR
!   do iblock =1,8
 !    do ieo = 1,2
 !      do isite=1,nvhalf
 !        do idirac=1,4
 !          do icolorir=1,5,2
!                site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
 !
 !               call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
 !                               vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 !
 !               call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
!            enddo ! icolorir
!          enddo ! idirac
 !       enddo ! isite
 !     enddo ! ieo
 !   enddo ! iblock
!
!  else ! isig

!    htemp = 0.0_KR
!    site = 0.0_KR
!    do iblock =1,8
!      do ieo = 1,2
!        do isite=1,nvhalf
 !         do idirac=1,4
 !           do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !               irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
 !                      print *, irow, x(icolorir,isite,idirac,ieo,iblock), x(icolorir+1,isite,idirac,ieo,iblock)
 !
 !           enddo ! icolorir
 !         enddo ! idirac
  !      enddo ! isite
  !    enddo ! ieo
  !  enddo ! iblock
 ! endif ! isig

! end subroutine matrix

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

  subroutine changevectorX(realz2noise,imagz2noise,bepart,bopart,numprocs,MRT,myid)

  real(kind=KR),    intent(out),    dimension(nxyzt,3,4)           :: realz2noise,imagz2noise
  real(kind=KR),    intent(in),     dimension(6,ntotal,4,2,8)        :: bepart
  real(kind=KR),    intent(in),     dimension(6,nvhalf,4,2,8)        :: bopart
  integer(kind=KI), intent(in)                                       :: myid,numprocs,MRT

  real(kind=KR),                    dimension(nxyzt,2,3,4)       :: evenz2,oddz2
  integer(kind=KI),                 dimension(nxyzt)                 :: hasdata

  real(kind=KR),                    dimension(6,ntotal,4,2,8)       :: betemp
  real(kind=KR),                    dimension(6,nvhalf,4,2,8)       :: botemp
  real(kind=KR),                    dimension(6,ntotal,4,2,16)      :: rtempb,itempb

  integer(kind=KI)                                                   :: iblock,ieo,i,&
                                                                        isite,icolor,idirac,site,&
                                                                        thesite, ic
  integer(kind=KI)                                                   :: kc,proc,count,ierr
  integer(kind=KI)                                                   :: counte, counto
  integer(kind=KI),                 dimension(2)                     :: ix
  integer(kind=KI)                                                   :: iy,iz,it
  integer(kind=KI) :: j,k,l,m,n
  integer(kind=KI) :: ieo1,ieo2,itbit,izbit,iybit,ixbit,itbit2,&
                      izbit2,iybit2,ixbit2,ixbit3,iblbit
  integer(kind=KI), dimension(4)  :: ip,np

  evenz2 = 0.0_KR
  oddz2 = 0.0_KR
  hasdata = 0
  realz2noise = 0.0_KR
  imagz2noise = 0.0_KR

  counte = (6*ntotal*4*2*8)
  counto = (6*nvhalf*4*2*8)

!!!! WARNING!!!!
! ABDOU ~ not sure in here, but it is VERY IMPORTANT that the shift index is handled with 
!         care. You must make sure that the looping order is correct such that we don't
!         "jumble" the data in array space.

  if (myid==0) then
      do proc=0,numprocs-1
         if (proc==0) then
             betemp = bepart
             botemp = bopart
         else
             call MPI_RECV(betemp, counte, MRT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
             call MPI_RECV(botemp, counto, MRT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         endif

! Working code
         do isite=1,nvhalf
            do ieo = 1,2
               do ic=1,3
                  do idirac = 1,4
                     rtempb(ic,isite,idirac,ieo,2)   = botemp(ic+(ic-1),isite,idirac,ieo,1)
                     rtempb(ic,isite,idirac,ieo,3)   = botemp(ic+(ic-1),isite,idirac,ieo,2)
                     rtempb(ic,isite,idirac,ieo,5)   = botemp(ic+(ic-1),isite,idirac,ieo,3)
                     rtempb(ic,isite,idirac,ieo,8)   = botemp(ic+(ic-1),isite,idirac,ieo,4)
                     rtempb(ic,isite,idirac,ieo,9)   = botemp(ic+(ic-1),isite,idirac,ieo,5)
                     rtempb(ic,isite,idirac,ieo,12)   = botemp(ic+(ic-1),isite,idirac,ieo,6)
                     rtempb(ic,isite,idirac,ieo,14)   = botemp(ic+(ic-1),isite,idirac,ieo,7)
                     rtempb(ic,isite,idirac,ieo,15)   = botemp(ic+(ic-1),isite,idirac,ieo,8)

                     itempb(ic,isite,idirac,ieo,2)   = botemp(2*ic,isite,idirac,ieo,1)
                     itempb(ic,isite,idirac,ieo,3)   = botemp(2*ic,isite,idirac,ieo,2)
                     itempb(ic,isite,idirac,ieo,5)   = botemp(2*ic,isite,idirac,ieo,3)
                     itempb(ic,isite,idirac,ieo,8)   = botemp(2*ic,isite,idirac,ieo,4)
                     itempb(ic,isite,idirac,ieo,9)   = botemp(2*ic,isite,idirac,ieo,5)
                     itempb(ic,isite,idirac,ieo,12)   = botemp(2*ic,isite,idirac,ieo,6)
                     itempb(ic,isite,idirac,ieo,14)   = botemp(2*ic,isite,idirac,ieo,7)
                     itempb(ic,isite,idirac,ieo,15)   = botemp(2*ic,isite,idirac,ieo,8)

                     rtempb(ic,isite,idirac,ieo,1)   = betemp(ic+(ic-1),isite,idirac,ieo,1)
                     rtempb(ic,isite,idirac,ieo,4)   = betemp(ic+(ic-1),isite,idirac,ieo,2)
                     rtempb(ic,isite,idirac,ieo,6)   = betemp(ic+(ic-1),isite,idirac,ieo,3)
                     rtempb(ic,isite,idirac,ieo,7)   = betemp(ic+(ic-1),isite,idirac,ieo,4)
                     rtempb(ic,isite,idirac,ieo,10)   = betemp(ic+(ic-1),isite,idirac,ieo,5)
                     rtempb(ic,isite,idirac,ieo,11)   = betemp(ic+(ic-1),isite,idirac,ieo,6)
                     rtempb(ic,isite,idirac,ieo,13)   = betemp(ic+(ic-1),isite,idirac,ieo,7)
                     rtempb(ic,isite,idirac,ieo,16)   = betemp(ic+(ic-1),isite,idirac,ieo,8)

                     itempb(ic,isite,idirac,ieo,1)   = betemp(2*ic,isite,idirac,ieo,1)
                     itempb(ic,isite,idirac,ieo,4)   = betemp(2*ic,isite,idirac,ieo,2)
                     itempb(ic,isite,idirac,ieo,6)   = betemp(2*ic,isite,idirac,ieo,3)
                     itempb(ic,isite,idirac,ieo,7)   = betemp(2*ic,isite,idirac,ieo,4)
                     itempb(ic,isite,idirac,ieo,10)   = betemp(2*ic,isite,idirac,ieo,5)
                     itempb(ic,isite,idirac,ieo,11)   = betemp(2*ic,isite,idirac,ieo,6)
                     itempb(ic,isite,idirac,ieo,13)   = betemp(2*ic,isite,idirac,ieo,7)
                     itempb(ic,isite,idirac,ieo,16)   = betemp(2*ic,isite,idirac,ieo,8)


                  enddo ! idirac
               enddo ! ic
            enddo ! ieo
         enddo ! isite

! The mapping of sites is detemined by the picture for the blocking in qqcd/cfqsprops/cfgsprops.f90.
! Notice that the "block" index is not linear, meaning that the site (first index) jumps from one to 
! two when the "block" changes to the right of the intial entry in the diagram. Our job is to
! map this non-linear fashion into the site which fills nxyzt.

  np(1) = npx
  np(2) = npy
  np(3) = npz
  np(4) = npt
  call atoc(proc,np,ip)

! Begin main loop, constructing ix,iy,iz,it from isite,ieo,ibl.
  isite = 0
  ieo1 = 2
  ieo2 = 1
  do itbit = 2,nt/npt,2
     itbit2 = itbit + ip(4)*nt/npt
     ieo1 = 3 - ieo1
     ieo2 = 3 - ieo2
     do izbit = 2,nz/npz,2
        izbit2 = izbit + ip(3)*nz/npz
        ieo1 = 3 - ieo1
        ieo2 = 3 - ieo2
        do iybit = 2,ny/npy,2
           iybit2 = iybit + ip(2)*ny/npy
           ieo1 = 3 - ieo1
           ieo2 = 3 - ieo2
           do ixbit = 4,nx/npx,4
              ixbit2 = ixbit + ip(1)*nx/npx
              isite = isite + 1
              do ieo = 1,2
                 do iblock = 1,16
                    if (iblock>8) then
                        it = itbit2
                    else
                        it = itbit2 - 1
                    endif ! (iblock>8)
                    iblbit = 1 + modulo(iblock-1,8)
                    if (iblbit>4) then
                        iz = izbit2
                    else
                        iz = izbit2 - 1
                    endif ! (iblbit>4)
                    iblbit = 1 + modulo(iblbit-1,4)
                    if (iblbit>2) then
                        iy = iybit2
                    else
                        iy = iybit2 - 1
                    endif ! (iblbit>2)
                    if (modulo(iblock,2)==1) then
                        ixbit3 = ixbit2 - 1
                    else
                        ixbit3 = ixbit2
                    endif ! (modulo==1)
                    ix(ieo1) = ixbit3 - 2
                    ix(ieo2) = ixbit3

                    thesite = ix(ieo) + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nz
                    do icolor = 1,3
                       do idirac = 1,4
                          realz2noise(thesite,icolor,idirac) = rtempb(icolor,isite,idirac,ieo,iblock)
                          imagz2noise(thesite,icolor,idirac) = itempb(icolor,isite,idirac,ieo,iblock)
                       enddo ! idirac
                    enddo ! icolor
                 enddo ! iblock
              enddo ! ieo
           enddo ! ixbit
        enddo ! iybit
     enddo ! izbit
  enddo ! itbit

      end do ! proc
  else
      call MPI_SEND(bepart, counte, MRT, 0, 0, MPI_COMM_WORLD, ierr)
      call MPI_SEND(bopart, counto, MRT, 0, 0, MPI_COMM_WORLD, ierr)
  endif

  end subroutine changevectorX

  subroutine atoc(ia,nc,ic)
    integer(kind=KI), intent(in)                   :: ia
    integer(kind=KI), intent(in),  dimension(:)    :: nc
    integer(kind=KI), intent(out), dimension(:)    :: ic

    ic(1) = modulo(ia,nc(1))
    ic(2) = modulo(ia/nc(1),nc(2))
    ic(3) = modulo(ia/(nc(1)*nc(2)),nc(3))
    ic(4) = modulo(ia/(nc(1)*nc(2)*nc(3)),nc(4))

  end subroutine atoc



 end module inverters
