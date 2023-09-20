!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! preconditioning.f90, Paul Lashomb, paul_lashomb@baylor.edu 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module contains tools for performing polynomial preconditioning of 
! the Wilson-Dirac matrix, as well as constructing polynomials for use in 
! the polynomial noise subtraction methods.    
!
!
!
! Comments: 
! 
! This module offers functionality for the "new" polynomial method, as 
! opposed to the "power" method of forming the polynomial. The required 
! subroutines for the "new" polynomial method are, 
!
! "new" poly: 
!     modlejacomp   -- performs Leja ordering of harmonic ritz values 
!     gmresEIGritz  -- obtains harmonic ritz values from gmres
!     gmres5EIGritz -- gamma5 hermitian form of gmresEIGritz  
!     newpoly       -- applies polynomial p(M). 
!     newphi        -- applies phi(M) = Mp(M). Used in preconditioning only 
!    
!
! In addition, there are some subroutines for the original "power" and 
! "Newton" method that are now deprecated and have been replaced with 
! the "new" polynomial method. Note that this is equivalent to the polynomial 
! methods hardcoded into ppgmresdr (inverters.f90), and both  xipolygenerator 
! and ppaverage2(disconloops.f90). The "old" polynomial method subroutines are, 
!
! "old" poly: 
!     powerbasis    -- constructs the coefficients for the "power" polynomial    
!     applyppgamma5 -- applies the polynomial from the "power" polynomial 
!     newtonbasis   -- plannned, but was obsolete after "new" polynomial created  
!                      and remains unfinished 
!
! Additional comments: 
!
! There are three subroutines -- vecdot, leastsquaresLAPACK, and 
! qrleastsquaresLAPACK -- which are duplicated from inverters.f90 but stored 
! as its own private subroutines to avoid circular references to inverters.f90. 
! It was decided that duplicating the very short subroutines was more readable 
! than adding gmresEIGritz, gmres5EIGritz, and powerbasis all in inverters.f90, 
! as each referenced these subroutines.    
!  
!
!
! Useful references: 
!
! "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
! Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
! "Polynomial Preconditioned GMRES and GMRES-DR", 
! Quan Liu, Ronald B Morgan, and Walter Wilcox, SIAM Journal on Sci. Comp. (2015) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     


 module preconditioning

!   use MPI
    use kinds
    use latdims
    use basics
    use diracops
    use pseudolapack
    use gmresrhs
!    use quark 
!    use printops 

    implicit none
    private

! Use the following line if the MPI module is not available.
    include 'mpif.h'

! Define access to subroutines.
!    public  :: ppminresdr5EIG,&
!               gmresEIGritz, gmres5EIGritz, newppgmresdr5EIG, newppminresdr5EIG 

    public  :: modlejacomp, newpoly, newdoublepoly, newphi, newdoublephi, & 
               powerbasis, applyppgamma5, gmresEIGritz, doublegmresEIGritz, &
               doublegmres5EIGritz, gmres5EIGritz, genpoly, gendoublepoly  
    private :: vecdot, leastsquaresLAPACK, qrfactorizationLAPACK  
 
 contains

!!!!!!!





! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine modlejacomp(th,thq,pp)

!===========================================================================
!   A subroutine to sort a set of complex elements using Modified Leja 
!   ordering 
!===========================================================================
!   INPUT:
!
!       th      = A set of either Ritz values or Harmonic Ritz values 
!       MRT     = MPIREALTYPE
!---------------------------------------------------------------------------
!   OUTPUT:
!
!       thq     = The set th after modifed Leja ordering
!---------------------------------------------------------------------------
!   References:
!
!       Z. Bai, D. Hu, and L. Reichel, A Newton basis GMRES implementation, 
!       IMA J. Numer. Anal., 14 (1994), pp. 563–581.
!
!---------------------------------------------------------------------------

!use shift

!============================Declarations================================

! IN/OUT

real(kind=KR2),        intent(in),      dimension(:,:)              :: th
real(kind=KR2),        intent(out),     dimension(:,:)              :: thq
integer(kind=KI),      intent(in)                                   :: pp

! LOCAL

real(kind=KR2),                         dimension(2,pp)             :: thr
real(kind=KR2),                         dimension(pp)               :: pr
real(kind=KR2)                                                      :: tmp1,tmp2,&
                                                                       prj
integer(kind=KI)                                                    :: i,j,ii

!===========================Initializations==============================

    !pp = size(th,2); 

    thr(:,:) = 0.0_KR2 
    thr(:,1) = th(:,1) 

    pr(:) = 0.0_KR2

!===========================Main Algorithm===============================


    do i = 2,pp

        tmp1 = sqrt(th(1,i)**2 + th(2,i)**2)
        tmp2 = sqrt(thr(1,1)**2 + thr(2,1)**2)

        if( tmp1 > tmp2 ) then
            thr(:,1) = th(:,i)
        endif

    enddo


    do j = 2,pp

        do i = 1,pp

            pr(i) = 0.0_KR2
            tmp1 = 0.0_KR2
            
            do ii = 1,j-1
                tmp1 = sqrt( ( th(1,i) - thr(1,ii) )**2 & 
                           + ( th(2,i) - thr(2,ii) )**2 )
                pr(i) = pr(i) + log(tmp1)
            enddo
        enddo

        thr(:,j) = th(:,1)

        prj = pr(1)


        do i = 2,pp

            if( pr(i) > prj ) then
                thr(:,j) = th(:,i)
                prj = pr(i)
            endif

        enddo

    enddo

    !thq = transpose(thr)
    thq = thr 
!===========================End of Algorithm=============================



endsubroutine modlejacomp
 !!!!!!





! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine newpoly(phie,phio,vINe,vINo,deg,hRV,MRT,ieo,ibleo,gblclr,&
                    kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                    bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2,tEE,tOO)

!===========================================================================
!   A subroutine to apply the polynomial p(M) to an input vector v and 
!   return the result (i.e. p(M)v) according to Algorithm 3 from the 
!   referenced paper. Note that this polynomial is what is used as the 
!   polynomial preconditioner, as well as the polynomial from the noise 
!   subtraction methods, i.e. M p(M)y = b.   
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg     = Degree of polynomial is deg-1
!       hRV    = Harmonic Ritz values for p(M)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Currently intended for application of new poly method
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================

integer(kind=KI),      intent(in)                                   :: deg
real(kind=KR2),        intent(in),     dimension(:,:,:,:,:)         :: vINe,vINo
real(kind=KR2),        intent(in),     dimension(2,deg)             :: hRV
real(kind=KR2),                        dimension(2,deg)             :: harv 
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                        dimension(6,ntotal,4,2,8)    :: tmpEE,tmpOO,&
                                                                       we,wo
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polye
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polyo
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prode 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prodo 
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phie                    
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phio 
!real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tmpe
!real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tmpo 
!real(kind=KR2)                                                      :: a,b
real(kind=KR2)                                                      :: realtol 
real(kind=KR2),                         dimension(2,deg)            :: smod ! square modulus  


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffe 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffo 
real(kind=KR2),                         dimension(2)                :: tmpe,tmpo,tmp10                         

real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: tEE, tOO 


character(len=128) :: newppfile 
integer(kind=KI)                                                    :: irhs 
integer(kind=KI) :: exists



!===========================Initializations==============================


polye(:6,:ntotal,:4,:2,:8) = 0.0_KR2
polyo(:6,:ntotal,:4,:2,:8) = 0.0_KR2
phie(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
phio(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
prode = vINe
prodo = vINo
smod(:2,:deg) = 0.0_KR2 
harv(:2,:deg) = hRV(:2,:deg) 


! Reciprocal of harmonic ritz values: 
! smod = 1 / (a + ib) = (a-ib) /(a**2 + b**2)  

do i = 1,deg 
    smod(1,i) =  1/(harv(1,i)**2 + harv(2,i)**2) * harv(1,i) 
    smod(2,i) = -1/(harv(1,i)**2 + harv(2,i)**2) * harv(2,i)  
enddo


!realtol = 1.0E-8 ! A tolerance to check if Im[rv(i)] == 0 


! If using Mg5 system, harmonic ritz values are real 
if (herm==1) then 
    do i=1,deg 
        harv(2,i) = 0.0_KR2 
    enddo
endif 


i = 1


!===========================Main Algorithm===============================

!write(*,*) "In main algorithm" 

do while( i <= deg-1 ) 
    !if( harv(2,i) .lt. realtol ) then 


    if (.true.) then !Comment out

     
     if (herm==1) then 
  
        do icri = 1,5,2
            do k = 1,nvhalf 
                polye(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) + (1/harv(1,i))*prode(icri  ,k,:,:,:) 
                polyo(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) + (1/harv(1,i))*prodo(icri  ,k,:,:,:) 
                polye(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:) + (1/harv(1,i))*prode(icri+1,k,:,:,:) 
                polyo(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) + (1/harv(1,i))*prodo(icri+1,k,:,:,:) 
            enddo
        enddo
     else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                polye(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) &
                                      + smod(1,i)*prode(icri  ,k,:,:,:) &
                                      - smod(2,i)*prode(icri+1,k,:,:,:)  
                polyo(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) & 
                                      + smod(1,i)*prodo(icri  ,k,:,:,:) & 
                                      - smod(2,i)*prodo(icri+1,k,:,:,:)  
                polye(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:)   & 
                                      + smod(2,i)*prode(icri  ,k,:,:,:) &    
                                      + smod(1,i)*prode(icri+1,k,:,:,:) 
                polyo(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) & 
                                      + smod(2,i)*prodo(icri  ,k,:,:,:) & 
                                      + smod(1,i)*prodo(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 
     endif     
    

    endif ! Comment out






    if (.true.) then ! Comment out



        ! Matrix multiplication for w = A*prod ( = A * tmp )

        tmpEE = prode 
        tmpOO = prodo 

        if (herm==1) then ! For M gamma5    
            do ieo=1,2
                do ibleo=1,8 
                    call gammamult(prode,prodo,tmpEE,tmpOO,5,ieo,ibleo)
                enddo
            enddo
        endif 


        gblclr = 2

        call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        gblclr = 1

        call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        we = tmpEE - kappa(1) * we
        wo = tmpOO - kappa(1) * wo





    if (.true.) then 

    if (herm==1) then         
        do icri = 1,5,2
            do k = 1,nvhalf 
                prode(icri  ,k,:,:,:) = prode(icri  ,k,:,:,:) - (1/harv(1,i))*we(icri  ,k,:,:,:) 
                prodo(icri  ,k,:,:,:) = prodo(icri  ,k,:,:,:) - (1/harv(1,i))*wo(icri  ,k,:,:,:) 
                prode(icri+1,k,:,:,:) = prode(icri+1,k,:,:,:) - (1/harv(1,i))*we(icri+1,k,:,:,:) 
                prodo(icri+1,k,:,:,:) = prodo(icri+1,k,:,:,:) - (1/harv(1,i))*wo(icri+1,k,:,:,:) 
            enddo
        enddo
    else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                prode(icri  ,k,:,:,:) = prode(icri  ,k,:,:,:) &
                                      - smod(1,i)*we(icri  ,k,:,:,:) &
                                      + smod(2,i)*we(icri+1,k,:,:,:)  
                prodo(icri  ,k,:,:,:) = prodo(icri  ,k,:,:,:) & 
                                      - smod(1,i)*wo(icri  ,k,:,:,:) & 
                                      + smod(2,i)*wo(icri+1,k,:,:,:)  
                prode(icri+1,k,:,:,:) = prode(icri+1,k,:,:,:)   & 
                                      - smod(2,i)*we(icri  ,k,:,:,:) &    
                                      - smod(1,i)*we(icri+1,k,:,:,:) 
                prodo(icri+1,k,:,:,:) = prodo(icri+1,k,:,:,:) & 
                                      - smod(2,i)*wo(icri  ,k,:,:,:) & 
                                      - smod(1,i)*wo(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 
    endif 




    endif 






         
        !if (.false.) then 

        !do icri = 1,5,2
        !    do k = 1,nvhalf 
        !        prode(icri  ,k,:,:,:) = - smod(1,i)*we(icri  ,k,:,:,:) &
        !                              + smod(2,i)*we(icri+1,k,:,:,:)  
        !        prodo(icri  ,k,:,:,:) = - smod(1,i)*wo(icri  ,k,:,:,:) & 
        !                              + smod(2,i)*wo(icri+1,k,:,:,:)  
        !        prode(icri+1,k,:,:,:) = - smod(2,i)*we(icri  ,k,:,:,:) &    
        !                              - smod(1,i)*we(icri+1,k,:,:,:) 
        !        prodo(icri+1,k,:,:,:) = - smod(2,i)*wo(icri  ,k,:,:,:) & 
        !                              - smod(1,i)*wo(icri+1,k,:,:,:)  
        !    enddo !k 
        !enddo !icri 


        ! End of matrix multiplication 

        !prode = prode - (1/harv(1,i))*we
        !prodo = prodo - (1/harv(1,i))*wo
       
        !endif 






    endif ! Comment out




        phie = polye
        phio = polyo 

        i = i + 1
!    else
!
!        !write(*,*) "In else" 
!        a = harv(1,i)  
!        b = harv(2,i)  
!        
!        ! Matrix multiplication for w = A*prod         
!
!
!        if (herm==1) then ! For M gamma5
!            we = prode 
!            wo = prodo
!            do ieo=1,2
!                do ibleo=1,8 
!                    call gammamult(prode,prodo,we,wo,5,ieo,ibleo)
!                enddo
!            enddo
!        else              ! For M 
!            we = prode 
!            wo = prodo 
!        endif 
!
!        gblclr = 2
!
!        call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
!                     nms,lvbc,ib,lbd,iblv,MRT)
!
!        gblclr = 1
!
!        call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
!                     nms,lvbc,ib,lbd,iblv,MRT)
!
!        we = tmpEE - kappa(1) * we
!        wo = tmpOO - kappa(1) * wo
!        
!        ! End of matrix multiplication 
!
!
!        tmpe = 2*a*prode - we
!        tmpo = 2*a*prodo - wo
!
!        polye = polye + (1/(a**2 + b**2)) * tmpe 
!        polyo = polyo + (1/(a**2 + b**2)) * tmpo 
!
!        if( i .le. deg - 2) then 
!            ! Matrix multiplication for w = A*tmp         
!
!
!
!            if (herm==1) then ! For M gamma5
!                we = tmpe 
!                wo = tmpo
!                do ieo=1,2
!                    do ibleo=1,8 
!                        call gammamult(tmpe,tmpo,we,wo,5,ieo,ibleo)
!                    enddo
!                enddo
!            else              ! For M 
!                we = tmpe 
!                wo = tmpo 
!            endif 
!
!
!            gblclr = 2
!
!            call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
!                         nms,lvbc,ib,lbd,iblv,MRT)
!
!            gblclr = 1
!
!            call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
!                         nms,lvbc,ib,lbd,iblv,MRT)
!
!            we = tmpEE - kappa(1) * we
!            wo = tmpOO - kappa(1) * wo
!        
!            ! End of matrix multiplication 
!
!            prode = prode - (1/(a**2 + b**2))*we
!            prodo = prodo - (1/(a**2 + b**2))*wo  
!        endif 
!        i = i + 2
    !endif 
enddo


!if( harv(2,deg) .lt. realtol ) then 


if (.true.) then ! Comment out


if (herm==1) then 
    do icri = 1,5,2
        do k = 1,nvhalf 
            phie(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) + (1/harv(1,deg))*prode(icri  ,k,:,:,:) 
            phio(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) + (1/harv(1,deg))*prodo(icri  ,k,:,:,:) 
            phie(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:) + (1/harv(1,deg))*prode(icri+1,k,:,:,:) 
            phio(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) + (1/harv(1,deg))*prodo(icri+1,k,:,:,:) 
        enddo
    enddo
else
    do icri = 1,5,2
        do k = 1,nvhalf
            phie(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) &
                               + smod(1,deg)*prode(icri  ,k,:,:,:) &
                               - smod(2,deg)*prode(icri+1,k,:,:,:)  
            phio(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) & 
                               + smod(1,deg)*prodo(icri  ,k,:,:,:) & 
                               - smod(2,deg)*prodo(icri+1,k,:,:,:)  
            phie(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:)   & 
                               + smod(2,deg)*prode(icri  ,k,:,:,:) &    
                               + smod(1,deg)*prode(icri+1,k,:,:,:) 
            phio(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) & 
                               + smod(2,deg)*prodo(icri  ,k,:,:,:) & 
                               + smod(1,deg)*prodo(icri+1,k,:,:,:)  
        enddo !k 
    enddo !icri 
endif 




endif ! Comment out 

        !polye = polye + (1/harv(1,deg))*prode
        !polyo = polyo + (1/harv(1,deg))*prodo

!endif 

!ye = polye 
!yo = polyo 


!write(*,*) "Finished!" 

!===========================End of Algorithm=============================



tEE = 0.0_KR2
tOO = 0.0_KR2 

!phie = we
!phio = wo

tEE = phie
tOO = phio


endsubroutine newpoly
 !!!!!!




! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine newdoublepoly(phie,phio,vINe,vINo,deg1,hRV1,deg2,hRV2,MRT,ieo,&
                          ibleo,gblclr,kappa,u,coact,vecbl,vecblinv,idag,myid,&
                          nn,iblv,bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2)

!===========================================================================
!   A subroutine to apply the double polynomial to an input vector v and 
!   return the result (i.e. p_1(M) * p_2 (phi_1(M)) * v) according to 
!   a modified form of Algorithm 1 from the referenced paper. Note that this  
!   polynomial is what is used as the polynomial preconditioner for double 
!   polynomial preconditioning, as well as the polynomial for subtraction 
!   using double polynomials. These allow for high-degree polynomials 
!   without requiring as large a subspace and without as much orthogonalization 
!   costs. 
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg1  = Degree of polynomial p_1(M) is deg-1
!       deg2  = Degree of polynomial p_2(M) is deg-1 
!       hRV1  = Harmonic Ritz values for p_1(M) 
!       hRV2  = Harmonic Ritz values for p_2(M)
!       MRT   = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!       
!       Currently intended for application of new double polynomial method
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================

integer(kind=KI),      intent(in)                                   :: deg1,deg2
real(kind=KR2),        intent(in),     dimension(:,:,:,:,:)         :: vINe,vINo
real(kind=KR2),        intent(in),     dimension(2,deg1)            :: hRV1
real(kind=KR2),        intent(in),     dimension(2,deg2)            :: hRV2
real(kind=KR2),                        dimension(2,deg1)            :: harv1
real(kind=KR2),                        dimension(2,deg2)            :: harv2
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                        dimension(6,ntotal,4,2,8)    :: tmpEE,tmpOO,&
                                                                       we,wo
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polye
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polyo
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prode 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prodo 
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phie,phio                    
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: phi1e,phi1o 
real(kind=KR2)                                                      :: realtol 
real(kind=KR2),                         dimension(2,deg1)           :: smod1 ! square modulus  
real(kind=KR2),                         dimension(2,deg2)           :: smod2 ! square modulus  

real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffe 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffo 
real(kind=KR2),                         dimension(2)                :: tmpe,tmpo,tmp10                         

!real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: tEE,tOO 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tEE,tOO 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tempEE,tempOO 


character(len=128) :: newppfile 
integer(kind=KI)                                                    :: irhs 
integer(kind=KI) :: exists



!===========================Initializations==============================


polye(:6,:ntotal,:4,:2,:8) = 0.0_KR2
polyo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

phie(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
phio(:6,:ntotal,:4,:2,:8)  = 0.0_KR2

phi1e(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
phi1o(:6,:ntotal,:4,:2,:8)  = 0.0_KR2

prode = vINe
prodo = vINo

smod1(:2,:deg1) = 0.0_KR2 
smod2(:2,:deg2) = 0.0_KR2 

harv1(:2,:deg1) = hRV1(:2,:deg1) 
harv2(:2,:deg2) = hRV2(:2,:deg2) 



! Reciprocal of harmonic ritz values: 
! smod = 1 / (a + ib) = (a-ib) /(a**2 + b**2)  

do i = 1,deg1 
    smod1(1,i) =  1/(harv1(1,i)**2 + harv1(2,i)**2) * harv1(1,i) 
    smod1(2,i) = -1/(harv1(1,i)**2 + harv1(2,i)**2) * harv1(2,i)  
enddo

do i = 1,deg2 
    smod2(1,i) =  1/(harv2(1,i)**2 + harv2(2,i)**2) * harv2(1,i) 
    smod2(2,i) = -1/(harv2(1,i)**2 + harv2(2,i)**2) * harv2(2,i)  
enddo

!realtol = 1.0E-8 ! A tolerance to check if Im[rv(i)] == 0 


! If using Mg5 system, harmonic ritz values are real 
if (herm==1) then 
    do i=1,deg1 
        harv1(2,i) = 0.0_KR2 
    enddo
endif 

if (herm==1) then 
    do i=1,deg2 
        harv2(2,i) = 0.0_KR2 
    enddo
endif 

i = 1



!===========================Main Algorithm===============================

!write(*,*) "In main algorithm" 

do while( i <= deg2-1 ) 
    !if( harv(2,i) .lt. realtol ) then 


    if (.true.) then !Comment out

     
     if (herm==1) then 
  
        do icri = 1,5,2
            do k = 1,nvhalf 
                polye(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) + (1/harv2(1,i))*prode(icri  ,k,:,:,:) 
                polyo(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) + (1/harv2(1,i))*prodo(icri  ,k,:,:,:) 
                polye(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:) + (1/harv2(1,i))*prode(icri+1,k,:,:,:) 
                polyo(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) + (1/harv2(1,i))*prodo(icri+1,k,:,:,:) 
            enddo
        enddo
     else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                polye(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) &
                                      + smod2(1,i)*prode(icri  ,k,:,:,:) &
                                      - smod2(2,i)*prode(icri+1,k,:,:,:)  
                polyo(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) & 
                                      + smod2(1,i)*prodo(icri  ,k,:,:,:) & 
                                      - smod2(2,i)*prodo(icri+1,k,:,:,:)  
                polye(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:)   & 
                                      + smod2(2,i)*prode(icri  ,k,:,:,:) &    
                                      + smod2(1,i)*prode(icri+1,k,:,:,:) 
                polyo(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) & 
                                      + smod2(2,i)*prodo(icri  ,k,:,:,:) & 
                                      + smod2(1,i)*prodo(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 
     endif     
    

    endif ! Comment out




    ! Replaces the matrix application M from newpoly with application of phi_1(M) 


    !tmpEE = prode 
    !tmpOO = prodo 

    !if (herm==1) then ! For M gamma5    
    !    do ieo=1,2
    !        do ibleo=1,8 
    !            call gammamult(prode,prodo,tmpEE,tmpOO,5,ieo,ibleo)
    !        enddo
    !    enddo
    !endif 


    call newphi(phi1e,phi1o,prode,prodo,deg1,harv1,MRT,ieo,ibleo,gblclr,&
                kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2) 



    if (.true.) then 

    if (herm==1) then         
        do icri = 1,5,2
            do k = 1,nvhalf 
                prode(icri  ,k,:,:,:) = prode(icri  ,k,:,:,:) - (1/harv2(1,i))*phi1e(icri  ,k,:,:,:) 
                prodo(icri  ,k,:,:,:) = prodo(icri  ,k,:,:,:) - (1/harv2(1,i))*phi1o(icri  ,k,:,:,:) 
                prode(icri+1,k,:,:,:) = prode(icri+1,k,:,:,:) - (1/harv2(1,i))*phi1e(icri+1,k,:,:,:) 
                prodo(icri+1,k,:,:,:) = prodo(icri+1,k,:,:,:) - (1/harv2(1,i))*phi1o(icri+1,k,:,:,:) 
            enddo
        enddo
    else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                prode(icri  ,k,:,:,:) = prode(icri  ,k,:,:,:) &
                                      - smod2(1,i)*phi1e(icri  ,k,:,:,:) &
                                      + smod2(2,i)*phi1e(icri+1,k,:,:,:)  
                prodo(icri  ,k,:,:,:) = prodo(icri  ,k,:,:,:) & 
                                      - smod2(1,i)*phi1o(icri  ,k,:,:,:) & 
                                      + smod2(2,i)*phi1o(icri+1,k,:,:,:)  
                prode(icri+1,k,:,:,:) = prode(icri+1,k,:,:,:)   & 
                                      - smod2(2,i)*phi1e(icri  ,k,:,:,:) &    
                                      - smod2(1,i)*phi1e(icri+1,k,:,:,:) 
                prodo(icri+1,k,:,:,:) = prodo(icri+1,k,:,:,:) & 
                                      - smod2(2,i)*phi1o(icri  ,k,:,:,:) & 
                                      - smod2(1,i)*phi1o(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 

    endif ! Hermitian test 



    endif ! Comment out




    phie = polye
    phio = polyo 



    i = i + 1

enddo



if (.true.) then ! Comment out


if (herm==1) then 
    do icri = 1,5,2
        do k = 1,nvhalf 
            phie(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) + (1/harv2(1,deg2))*prode(icri  ,k,:,:,:) 
            phio(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) + (1/harv2(1,deg2))*prodo(icri  ,k,:,:,:) 
            phie(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:) + (1/harv2(1,deg2))*prode(icri+1,k,:,:,:) 
            phio(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) + (1/harv2(1,deg2))*prodo(icri+1,k,:,:,:) 
        enddo
    enddo
else
    do icri = 1,5,2
        do k = 1,nvhalf
            phie(icri  ,k,:,:,:) = polye(icri  ,k,:,:,:) &
                               + smod2(1,deg2)*prode(icri  ,k,:,:,:) &
                               - smod2(2,deg2)*prode(icri+1,k,:,:,:)  
            phio(icri  ,k,:,:,:) = polyo(icri  ,k,:,:,:) & 
                               + smod2(1,deg2)*prodo(icri  ,k,:,:,:) & 
                               - smod2(2,deg2)*prodo(icri+1,k,:,:,:)  
            phie(icri+1,k,:,:,:) = polye(icri+1,k,:,:,:)   & 
                               + smod2(2,deg2)*prode(icri  ,k,:,:,:) &    
                               + smod2(1,deg2)*prode(icri+1,k,:,:,:) 
            phio(icri+1,k,:,:,:) = polyo(icri+1,k,:,:,:) & 
                               + smod2(2,deg2)*prodo(icri  ,k,:,:,:) & 
                               + smod2(1,deg2)*prodo(icri+1,k,:,:,:)  
        enddo !k 
    enddo !icri 
endif 




endif ! Comment out 


!===========================End of Algorithm=============================



tEE = 0.0_KR2
tOO = 0.0_KR2 


tEE = phie
tOO = phio

tempEE = 0.0_KR2 
tempOO = 0.0_KR2 


! Applies the final p_1(M) to the vector p_2(phi_1(M))v  

call newpoly( tEE,tOO,phie,phio,deg1,harv1,MRT,ieo,ibleo,gblclr,kappa,u,&
              coact,vecbl,vecblinv,idag,myid,nn,iblv,bc,nms,ib,lbd,ldiv,&
              lvbc,herm,MRT2,tempEE,tempOO )

phie = tEE 
phio = tOO 







endsubroutine newdoublepoly
 !!!!!!



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine newphi(phie,phio,vINe,vINo,deg,hRV,MRT,ieo,ibleo,gblclr,&
                   kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                   bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2)

!===========================================================================
!   A subroutine to apply the polynomial phi(M) to an input vector v and 
!   return the result (i.e. phi(M)v) according to Algorithm 1 from the 
!   referenced paper. Note that this polynomial is the same as the 
!   polynomial preconditioned system, i.e. M p(M)y = phi(M)x = b.   
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg     = Degree of polynomial is deg-1
!       hRV    = Harmonic Ritz values for p(M)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Currently intended for polynomial preconditioning involving new 
!       polynomial method. 
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================

integer(kind=KI),      intent(in)                                   :: deg
real(kind=KR2),        intent(in),     dimension(:,:,:,:,:)         :: vINe,vINo
real(kind=KR2),        intent(in),     dimension(2,deg)             :: hRV
real(kind=KR2),                        dimension(2,deg)             :: harv 
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                        dimension(6,ntotal,4,2,8)    :: tmpEE,tmpOO,&
                                                                       piE,piO
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polye
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: polyo
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prode 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: prodo 
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phie                    
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phio 
!real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tmpe
!real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tmpo 
!real(kind=KR2)                                                      :: a,b
real(kind=KR2)                                                      :: realtol 
real(kind=KR2),                         dimension(2,deg)            :: smod ! square modulus  


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffe 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffo 
real(kind=KR2),                         dimension(2)                :: tmpe,tmpo,tmp10                         

!real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: tEE, tOO 


character(len=128) :: newppfile 
integer(kind=KI)                                                    :: irhs 
integer(kind=KI) :: exists



!===========================Initializations==============================


polye(:6,:ntotal,:4,:2,:8) = 0.0_KR2
polyo(:6,:ntotal,:4,:2,:8) = 0.0_KR2
phie(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
phio(:6,:ntotal,:4,:2,:8)  = 0.0_KR2
piE = vINe
piO = vINo
smod(:2,:deg) = 0.0_KR2 
harv(:2,:deg) = hRV(:2,:deg) 


! Reciprocal of harmonic ritz values: 
! smod = 1 / (a + ib) = (a-ib) /(a**2 + b**2)  

do i = 1,deg 
    smod(1,i) =  1/(harv(1,i)**2 + harv(2,i)**2) * harv(1,i) 
    smod(2,i) = -1/(harv(1,i)**2 + harv(2,i)**2) * harv(2,i)  
enddo



! If using Mg5 system, harmonic ritz values are real 
if (herm==1) then 
    do i=1,deg 
        harv(2,i) = 0.0_KR2 
    enddo
endif 



i = 1



!===========================Main Algorithm===============================

!write(*,*) "In main algorithm" 

do while( i <= deg ) 
    !if( harv(2,i) .lt. realtol ) then 



! Line 4

    tmpEE = piE 
    tmpOO = piO 

    if (herm==1) then ! For M gamma5    
        do ieo=1,2
            do ibleo=1,8 
                call gammamult(piE,piO,tmpEE,tmpOO,5,ieo,ibleo)
            enddo
        enddo
    endif 


    gblclr = 2

    call Hsingle(prodo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                 nms,lvbc,ib,lbd,iblv,MRT)

    gblclr = 1

    call Hsingle(prode,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                 nms,lvbc,ib,lbd,iblv,MRT)

    prode = tmpEE - kappa(1) * prode
    prodo = tmpOO - kappa(1) * prodo




! Line 5

    if (.true.) then !Comment out

     
     if (herm==1) then 
  
        do icri = 1,5,2
            do k = 1,nvhalf 
                piE(icri  ,k,:,:,:) = piE(icri  ,k,:,:,:) - (1/harv(1,i))*prode(icri  ,k,:,:,:) 
                piO(icri  ,k,:,:,:) = piO(icri  ,k,:,:,:) - (1/harv(1,i))*prodo(icri  ,k,:,:,:) 
                piE(icri+1,k,:,:,:) = piE(icri+1,k,:,:,:) - (1/harv(1,i))*prode(icri+1,k,:,:,:) 
                piO(icri+1,k,:,:,:) = piO(icri+1,k,:,:,:) - (1/harv(1,i))*prodo(icri+1,k,:,:,:) 
            enddo
        enddo
     else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                piE(icri  ,k,:,:,:) = piE(icri  ,k,:,:,:) &
                                      - smod(1,i)*prode(icri  ,k,:,:,:) &
                                      + smod(2,i)*prode(icri+1,k,:,:,:)  
                piO(icri  ,k,:,:,:) = piO(icri  ,k,:,:,:) & 
                                      - smod(1,i)*prodo(icri  ,k,:,:,:) & 
                                      + smod(2,i)*prodo(icri+1,k,:,:,:)  
                piE(icri+1,k,:,:,:) = piE(icri+1,k,:,:,:)   & 
                                      - smod(2,i)*prode(icri  ,k,:,:,:) &    
                                      - smod(1,i)*prode(icri+1,k,:,:,:) 
                piO(icri+1,k,:,:,:) = piO(icri+1,k,:,:,:) & 
                                      - smod(2,i)*prodo(icri  ,k,:,:,:) & 
                                      - smod(1,i)*prodo(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 
     endif     
    

    endif ! Comment out


! Line 6 

    i = i + 1


enddo



! Line 15 

if (.true.) then ! Comment out

do icri = 1,5,2
    do k = 1,nvhalf 
        phie(icri  ,k,:,:,:) = vINe(icri  ,k,:,:,:) - piE(icri  ,k,:,:,:) 
        phio(icri  ,k,:,:,:) = vINo(icri  ,k,:,:,:) - piO(icri  ,k,:,:,:) 
        phie(icri+1,k,:,:,:) = vINe(icri+1,k,:,:,:) - piE(icri+1,k,:,:,:) 
        phio(icri+1,k,:,:,:) = vINo(icri+1,k,:,:,:) - piO(icri+1,k,:,:,:) 
    enddo !k 
enddo !icri 

endif ! Comment out 


!===========================End of Algorithm=============================



!tEE = 0.0_KR2
!tOO = 0.0_KR2 

!tEE = phie
!tOO = phio


 endsubroutine newphi
 !!!!!!






! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine newdoublephi(phi2e,phi2o,vINe,vINo,deg1,hRV1,deg2,hRV2,MRT,ieo,&
                         ibleo,gblclr,kappa,u,coact,vecbl,vecblinv,idag,myid,& 
                         nn,iblv,bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2)

!===========================================================================
!   A subroutine to apply the polynomial phi_2(phi_1(M)) to an input vector v 
!   and return the result (i.e. phi_2(phi_1(M))v) according to Algorithm 1 
!   from the referenced paper. Note that this polynomial is the same as the 
!   double polynomial preconditioned system, i.e., 
!   phi_2(phi_1(M))y = M p_1(M) p_2(phi_1(M))y = b.   
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg     = Degree of polynomial is deg-1
!       hRV    = Harmonic Ritz values for p(M)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Currently intended for polynomial preconditioning involving new 
!       polynomial method. 
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================

integer(kind=KI),      intent(in)                                   :: deg1,deg2
real(kind=KR2),        intent(in),     dimension(:,:,:,:,:)         :: vINe,vINo
real(kind=KR2),        intent(in),     dimension(2,deg1)            :: hRV1
real(kind=KR2),        intent(in),     dimension(2,deg2)            :: hRV2
real(kind=KR2),                        dimension(2,deg1)            :: harv1 
real(kind=KR2),                        dimension(2,deg2)            :: harv2 
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                        dimension(6,ntotal,4,2,8)    :: tmpEE,tmpOO,&
                                                                       piE,piO
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: phi1e,phi1o                    
real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: phi2e,phi2o 
real(kind=KR2)                                                      :: realtol 
real(kind=KR2),                         dimension(2,deg1)           :: smod1 ! square modulus  
real(kind=KR2),                         dimension(2,deg2)           :: smod2 ! square modulus  


real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffe 
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: diffo 
real(kind=KR2),                         dimension(2)                :: tmpe,tmpo,tmp10                         

!real(kind=KR2),        intent(out),     dimension(6,ntotal,4,2,8)   :: tEE, tOO 


character(len=128) :: newppfile 
integer(kind=KI)                                                    :: irhs 
integer(kind=KI) :: exists



!===========================Initializations==============================

phi1e(:6,:ntotal,:4,:2,:8) = 0.0_KR2
phi1o(:6,:ntotal,:4,:2,:8) = 0.0_KR2

phi2e(:6,:ntotal,:4,:2,:8) = 0.0_KR2
phi2o(:6,:ntotal,:4,:2,:8) = 0.0_KR2

piE = vINe
piO = vINo

smod2(:2,:deg2) = 0.0_KR2 

harv1(:2,:deg1) = hRV1(:2,:deg1) 
harv2(:2,:deg2) = hRV2(:2,:deg2) 


! Reciprocal of harmonic ritz values: 
! smod = 1 / (a + ib) = (a-ib) /(a**2 + b**2)  

do i = 1,deg2 
    smod2(1,i) =  1/(harv2(1,i)**2 + harv2(2,i)**2) * harv2(1,i) 
    smod2(2,i) = -1/(harv2(1,i)**2 + harv2(2,i)**2) * harv2(2,i)  
enddo



! If using Mg5 system, harmonic ritz values are real 

if (herm==1) then 
    do i=1,deg1 
        harv1(2,i) = 0.0_KR2 
    enddo
endif 

if (herm==1) then 
    do i=1,deg2 
        harv2(2,i) = 0.0_KR2 
    enddo
endif 


i = 1



!===========================Main Algorithm===============================

!write(*,*) "In main algorithm" 

do while( i <= deg2 ) 
    !if( harv(2,i) .lt. realtol ) then 



! Line 4

    tmpEE = piE 
    tmpOO = piO 

    ! Not sure if I need this for hermitian 
    !if (herm==1) then ! For M gamma5    
    !    do ieo=1,2
    !        do ibleo=1,8 
    !            call gammamult(piE,piO,tmpEE,tmpOO,5,ieo,ibleo)
    !        enddo
    !    enddo
    !endif 




! Replaces matrix application M with application of phi_1(M) 

    call newphi(phi1e,phi1o,tmpEE,tmpOO,deg1,harv1,MRT,ieo,ibleo,gblclr,&
                kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2) 
    


! Line 5

    if (.true.) then !Comment out

     
     if (herm==1) then 
  
        do icri = 1,5,2
            do k = 1,nvhalf 
                piE(icri  ,k,:,:,:) = piE(icri  ,k,:,:,:) - (1/harv2(1,i))*phi1e(icri  ,k,:,:,:) 
                piO(icri  ,k,:,:,:) = piO(icri  ,k,:,:,:) - (1/harv2(1,i))*phi1o(icri  ,k,:,:,:) 
                piE(icri+1,k,:,:,:) = piE(icri+1,k,:,:,:) - (1/harv2(1,i))*phi1e(icri+1,k,:,:,:) 
                piO(icri+1,k,:,:,:) = piO(icri+1,k,:,:,:) - (1/harv2(1,i))*phi1o(icri+1,k,:,:,:) 
            enddo
        enddo
     else 
        do icri = 1,5,2
            do k = 1,nvhalf 
                piE(icri  ,k,:,:,:) = piE(icri  ,k,:,:,:) &
                                      - smod2(1,i)*phi1e(icri  ,k,:,:,:) &
                                      + smod2(2,i)*phi1e(icri+1,k,:,:,:)  
                piO(icri  ,k,:,:,:) = piO(icri  ,k,:,:,:) & 
                                      - smod2(1,i)*phi1o(icri  ,k,:,:,:) & 
                                      + smod2(2,i)*phi1o(icri+1,k,:,:,:)  
                piE(icri+1,k,:,:,:) = piE(icri+1,k,:,:,:)   & 
                                      - smod2(2,i)*phi1e(icri  ,k,:,:,:) &    
                                      - smod2(1,i)*phi1e(icri+1,k,:,:,:) 
                piO(icri+1,k,:,:,:) = piO(icri+1,k,:,:,:) & 
                                      - smod2(2,i)*phi1o(icri  ,k,:,:,:) & 
                                      - smod2(1,i)*phi1o(icri+1,k,:,:,:)  
            enddo !k 
        enddo !icri 
     endif     
    

    endif ! Comment out


! Line 6 

    i = i + 1


enddo



! Line 15 

if (.true.) then ! Comment out

do icri = 1,5,2
    do k = 1,nvhalf 
        phi2e(icri  ,k,:,:,:) = vINe(icri  ,k,:,:,:) - piE(icri  ,k,:,:,:) 
        phi2o(icri  ,k,:,:,:) = vINo(icri  ,k,:,:,:) - piO(icri  ,k,:,:,:) 
        phi2e(icri+1,k,:,:,:) = vINe(icri+1,k,:,:,:) - piE(icri+1,k,:,:,:) 
        phi2o(icri+1,k,:,:,:) = vINo(icri+1,k,:,:,:) - piO(icri+1,k,:,:,:) 
    enddo !k 
enddo !icri 

endif ! Comment out 


!===========================End of Algorithm=============================



!tEE = 0.0_KR2
!tOO = 0.0_KR2 

!tEE = phie
!tOO = phio


 endsubroutine newdoublephi
 !!!!!!




! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine genpoly(rwdir,lorv,deg,be,bo,MRT,ieo,ibleo,gblclr,kappa,u,coact,vecbl,vecblinv,&
                    idag,myid,nn,iblv,bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2)

!===========================================================================
!   A subroutine to generate the polynomial phi(M) and p(M). 
!   Runs a single cycle of GMRES, Modified Leja orders them, and returns them  
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg     = Degree of polynomial is deg-1
!       hRV    = Harmonic Ritz values for p(M)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Currently intended for application of new poly method
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================
character(len=*),      intent(in),      dimension(:)                :: rwdir
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: be, bo
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------


    ! Params for gmresEIGritz call 
    
integer(kind=KI),   dimension(2)                :: paramsGMRES 
real(kind=KR),      dimension(18,nvhalf,8,2,16) :: GeeGooinv  !Not sure what this is
integer(kind=KI)                                :: iflag  !not sure what this is



!------------------------------------------------------
!------------------------------------------------------

! SET TO DESIRED POLYNOMIAL DEGREE 
integer(kind=KI),      intent(in)                                   :: deg
real(kind=KR2),     dimension(2,1000)                               :: rv ! Harmonic ritz values 
real(kind=KR2),     intent(out),         dimension(2,1000)          :: lorv !Leja-ordered harmonic ritz values

!-------------------------------------------------------
!------------------------------------------------------


!    integer(kind=KI)                                :: i ! index 
real(kind=KR2),     dimension(6,ntotal,4,2,8)   :: xejunk, xojunk 




!===========================Initializations==============================

    paramsGMRES(1) = deg ! size of Krylov subspace
    paramsGMRES(2) = 1   ! number of eigenvector/values to deflate 

    xejunk = 0.0_KR2 
    xojunk = 0.0_KR2

!===========================Main Algorithm===============================


    ! changed be,bo to onesE,onesO
    call gmresEIGritz(rwdir,be,bo,xejunk,xojunk,paramsGMRES,u,GeeGooinv, &
                      iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                      lvbc,ib,lbd,iblv,MRT,MRT2,rv)


    if (.true.) then 
        if (myid==0) then 
            do i = 1,deg
                print*, 'rv:        ', rv(1,i), rv(2,i)  
            enddo 
        endif
    endif 
    ! Call modlejacomp to perform modified Leja ordering on Harmonic ritz values
    
    call modlejacomp(rv,lorv,deg)    

    if (.true.) then 
        if (myid==0) then 
            do i = 1,deg
                print*, 'lorv:        ', lorv(1,i), lorv(2,i) 
            enddo
        endif 
    endif 

!===========================End of Algorithm=============================




endsubroutine genpoly 
 !!!!!!




! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine gendoublepoly(rwdir,lorv1,lorv2,deg1,deg2,be,bo,MRT,ieo,ibleo,gblclr,&
                          kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                          bc,nms,ib,lbd,ldiv,lvbc,herm,MRT2)

!===========================================================================
!   A subroutine to generate the double polynomial phi_2(phi_1(M)) and 
!   p_1(M) p_2(phi_1(A)). 
!   Runs a single cycle of GMRES on M to obtain harmonic ritz values, 
!   Modified Leja orders them to form phi_1(M), runs another cycle of GMRES 
!   on phi_1(M) to obtain harmonic ritz values, Modified Leja orders them, 
!   and returns them to form phi_2(phi_1(M)) and p_1(M) p_2(phi_1(A)). 
!
!   References: 
!
!   "Toward Efficient and Stable Polynomial Preconditioning for GMRES.",
!   Jennifer A. Loe and Ronald B. Morgan,arXiv: Numerical Analysis (2020)
!
!===========================================================================
!   INPUT:
!
!       be    = Even vector z to apply p(M) to
!       bo    = Odd vector z to apply p(M) to
!       deg     = Degree of polynomial is deg-1
!       hRV    = Harmonic Ritz values for p(M)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       polye   = vIne with GMRES polynomial p(M)
!                 applied to it. (i.e. polye = p(M) * vIne)
!       polyo   = vIno with GMRES polynomial p(M)
!                 applied to it. (i.e. polyo = p(M) * vIno)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Currently intended for application of new poly method
!
!--------------------------------------------------------------------------
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================
character(len=*),      intent(in),      dimension(:)                :: rwdir
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: be, bo
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                        dimension(6,ntotal,4,2,8)    :: tmpEE,tmpOO,&
                                                                       we,wo
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,k,&
                                                                       icri
integer(kind=KI),      intent(in)                                   :: herm

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------

    ! Params for gmresEIGritz call 
    
integer(kind=KI),   dimension(2)                :: paramsGMRES1, paramsGMRES2 
real(kind=KR),      dimension(18,nvhalf,8,2,16) :: GeeGooinv  !Not sure what this is
integer(kind=KI)                                :: iflag  !not sure what this is



!------------------------------------------------------
!------------------------------------------------------

! SET TO DESIRED POLYNOMIAL DEGREE 
integer(kind=KI),      intent(in)                                   :: deg1,deg2 
real(kind=KR2),                         dimension(2,1000)           :: rv1,rv2     !Harmonic ritz values 
real(kind=KR2),                         dimension(2,1000)           :: lorv1,lorv2 !Modified Leja-ordered harmonic ritz values

!-------------------------------------------------------
!------------------------------------------------------


!    integer(kind=KI)                                :: i ! index 
real(kind=KR2),     dimension(6,ntotal,4,2,8)   :: xejunk, xojunk 



!===========================Initializations==============================

    paramsGMRES1(1) = deg1 ! dimension of Krylov subspace  
    paramsGMRES1(2) = 1    ! number of eigenvector/values to deflate (doesn't matter as deflation                           
                           ! portion is skipped and harmonic ritz values obtained before it  

    paramsGMRES2(1) = deg2 ! dimension of second Krylov subspace 
    paramsGMRES2(2) = 1    ! same as paramsGMRES1(2) 

	
    xejunk = 0.0_KR2 
    xojunk = 0.0_KR2

!===========================Main Algorithm===============================


    !%%%%% FORMING FIRST POLYNOMIAL %%%%%    
    ! Finding harmonic ritz values for A to use in forming phi_1(A) and p_1(A)  ( phi_1(A) y = b, where x = p_1(A) y)  


    xejunk(:6,:ntotal,:4,:2,:8) = 0.0_KR2
    xojunk(:6,:ntotal,:4,:2,:8) = 0.0_KR2      


    call gmresEIGritz(rwdir,be,bo,xejunk,xojunk,paramsGMRES1,u,GeeGooinv, &
                      iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                      lvbc,ib,lbd,iblv,MRT,MRT2,rv1)

    if (myid==0) then 
        do i = 1,deg1 
            print*, 'rv1:        ', rv1(1,i), rv1(2,i)  
        enddo 
    endif

    ! Call modlejacomp to perform modified Leja ordering on Harmonic ritz values
    
    call modlejacomp(rv1,lorv1,deg1)    

    if (myid==0) then 
        do i = 1,deg1 
            print*, 'lorv1:        ', lorv1(1,i), lorv1(2,i)  
        enddo 
    endif



    !%%%%% FORMING SECOND POLYNOMIAL %%%%%
    ! Finding harmonic ritz values for phi_1(A) to use in forming phi_2(phi_1(A)) and p_2(phi_1(A))  ( phi_2(phi_1(A)) z = b, where x = p_1(A) p_2(phi_1(A)) z) 

    xejunk(:6,:ntotal,:4,:2,:8) = 0.0_KR2
    xojunk(:6,:ntotal,:4,:2,:8) = 0.0_KR2      

    call doublegmresEIGritz(rwdir,be,bo,xejunk,xojunk,paramsGMRES2,lorv1,u,GeeGooinv, &
                      iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                      lvbc,ib,lbd,iblv,MRT,MRT2,rv2)

    if (myid==0) then 
        do i = 1,deg2 
            print*, 'rv2:        ', rv2(1,i), rv2(2,i)  
        enddo 
    endif

    ! Call modlejacomp to perform modified Leja ordering on Harmonic ritz values
    
    call modlejacomp(rv2,lorv2,deg2)    

    if (myid==0) then 
        do i = 1,deg2  
            print*, 'lorv2:        ', lorv2(1,i), lorv2(2,i)  
        enddo 
    endif


!===========================End of Algorithm=============================


endsubroutine gendoublepoly 
 !!!!!!





! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine powerbasis(phie,phio,p,co,MRT,MRT2,ieo,ibleo,gblclr,kappa,u,coact,&
                       vecbl,vecblinv,idag,myid,nn,iblv,bc,nms,ib,lbd,ldiv,lvbc)

!==========================================================================
!   A MinRes-DR(mdr,kdr) polynomial preconditioner constructed through
!   solving the Normal Equations, V'Vg=V'b. In paper by Dr. Morgan on polynomial
!   preconditioned GMRES-DR, V=AY. Note that this method is identical to the 
!   original method used in the "power" method of noise subtraction that was 
!   limited to degree 7 polynomials, as well as what was hard-coded for 
!   preconditioning in ppgmresdr. The "power" method has since been replaced 
!   with the "new" polynomial method, outlined in the subroutine newpoly.   
!
!   References: 
!   
!   "Polynomial Preconditioned GMRES and GMRES-DR", 
!   Quan Liu, Ronald B Morgan, and Walter Wilcox, SIAM Journal on Sci. Comp. (2015) 
!
!==========================================================================
!  INPUT:
!       phie    = the "even" part of the right hand side of M(gamma5)x=phi
!       phio    = the "odd" part of the right hand side of M(gamma5)x=phi 
!       p       = polynomial degree.
!       MRT     = MPIREALTYPE
!       MRT2    = MPIDBLETYPE
!--------------------------------------------------------------------------
!  OUTPUT:
!       co      = coefficients for polynomial 
!--------------------------------------------------------------------------

!use shift


!===========================Initializations==============================

real(kind=KR),         intent(in),      dimension(6,ntotal,4,2,8)   :: phie,phio
integer(kind=KI),      intent(in)                                   :: p
real(kind=KR2),                         dimension(p)                :: beta,&
                                                                       betaE,betaO
real(kind=KR2),                         dimension(2,p,p)            :: lsmat,&
                                                                       lsmatE,&
                                                                       lsmatO
real(kind=KR2),                         dimension(2,p,1)            :: cls,clsE,&
                                                                       clsO
real(kind=KR2),        intent(out),     dimension(2,p)              :: co
real(kind=KR2),                         dimension(2,p)              :: coE,coO
integer(kind=KI),      intent(in)                                   :: MRT,MRT2
real(kind=KR2),                         dimension(6,ntotal,4,2,8)   :: tmpEE,tmpOO,&
                                                                       we,wo
real(kind=KR2),                         dimension(6,ntotal,4,2,8,p+1) :: vprimee,&
                                                                         vprimeo
integer(kind=KI)                                                    :: ieo,ibleo,&
                                                                       gblclr,i,j,k

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u 
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact          
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv 
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
integer(kind=KI),                       dimension(2)                :: ipiv2
!------------------------------------------------------------------------


we = 0.0_KR2
wo = 0.0_KR2
ipiv2 = 0.0_KI

!=======================Determining the polynomial=======================


!lsmat cls have 2 for first argument because they are for the even (1) and odd
!(2) parts

!Try using the "full" AY matrix with the "full" b vector to find the
!coefficients. Also, try using the "AY" constructed from just the even (odd)
!parts and use only the even (odd) RHS. Check and see if the coefficients are
!all the same. If they're all the same, then just use the even part rather than
!using both the even and odd parts.


!-------------------------Constructing V for V'Vg=V'b----------------------

    do k = 1,nvhalf
        vprimee(:,k,:,:,:,1) = phie(:,k,:,:,:)
        vprimeo(:,k,:,:,:,1) = phio(:,k,:,:,:)
    enddo !k


    do i = 1,p
        do ieo=1,2
            do ibleo=1,8
                call gammamult( vprimee(:,:,:,:,:,i),vprimeo(:,:,:,:,:,i),&
                                tmpEE,tmpOO,5,ieo,ibleo)
            enddo !ibleo
        enddo !ieo
            

        gblclr = 2

        call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        gblclr = 1

        call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        we = tmpEE - kappa(1) * we
        wo = tmpOO - kappa(1) * wo


        vprimee(:,:,:,:,:,i+1) = we(:,:,:,:,:)
        vprimeo(:,:,:,:,:,i+1) = wo(:,:,:,:,:)

    enddo !i


!----------------Constructing Matrix for Normal Equations------------------

    do i=2,p+1
        do j=2,p+1

            call vecdot(vprimee(:,:,:,:,:,i),vprimee(:,:,:,:,:,j),betaE,MRT2)
            call vecdot(vprimeo(:,:,:,:,:,i),vprimeo(:,:,:,:,:,j),betaO,MRT2)

    !START OF "FULL" MATRIX V'V (V = AY)
            beta(:) = betaE(:) + betaO(:)
            lsmat(:,i-1,j-1) = beta(:)
    !END

    !START OF "EVEN" ("ODD") MATRIX Ve'Ve (Vo'Vo)
            !lsmatE(:,i-1,j-1) = betaE(:)  !lsmat(2,p,p) ,cls(2,p,1)
            !lsmatO(:,i-1,j-1) = betaO(:)  !print *, "i,j, lsmat(:,i,j)=",
            !i-1,j-1, lsmat(:,i-1,j-1)
    !END


        enddo!j
    enddo!i


!--------------------Constructing RHS for Normal Equations------------------


    do i=2,p+1

        call vecdot(vprimee(:,:,:,:,:,i),phie(:,:,:,:,:),betaE,MRT2)
        call vecdot(vprimeo(:,:,:,:,:,i),phio(:,:,:,:,:),betaO,MRT2)


    !START OF "FULL" MATRIX RHS
       beta(:) = betaE(:) + betaO(:)
       cls(:,i-1,1) = beta(:)
    !END

    !START OF "EVEN" ("ODD") MATRIX RHS
       !clsE(:,i-1,1) = betaE(:)  !print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
       !clsO(:,i-1,1) = betaO(:)
    !END

    enddo!i


!---------------------Solving (lsmat)co =(cls) --------------------------

    !I think the even/odd pieces are all correct except for the linear solver.
    !It needs to be over the entire matrix to get all of the coefficients. The
    !matrix lsmat, though, needs to be changed slightly to make it right.

    !START COEFFICIENTS FROM "FULL" MATRIX
        call linearsolver(p,1,lsmat,ipiv2,cls)
        co(:,:) = cls(:,:,1)
    !END

    !Print out the coefficients temporarily 

   ! do i=1,p
   !     print *, 'p, co(:,p)=',i, co(:,i)  
   ! enddo !i

    !START COEFFICIENTS FROM "EVEN" ("ODD") MATRIX
        !call linearsolver(p,1,lsmatE,ipiv2,clsE)
        !call linearsolver(p,1,lsmatO,ipiv2,clsO)
        !coE(:,:) = clsE(:,:,1)
        !coO(:,:) = clsO(:,:,1)
    !END


!===========================End of Algorithm=============================


 endsubroutine powerbasis
 !!!!!!


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine applyppgamma5(trye,tryo,p,co,ye,yo,MRT,ieo,ibleo,gblclr,&
                          kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                          bc,nms,ib,lbd,ldiv,lvbc)

!===========================================================================
!   A subroutine to implement the polynomial preconditioner p(Ag5) by applying
!   it to a vector z and returning it (i.e. y = p(Ag5)z). This was an earlier 
!   implementation of the polynomial for subtraction to be used with the 
!   "power" basis and the "Newton" basis. See the following references for 
!   more details in section 2.1 for the "power" method, and section 3.3 for 
!   "Newton polynomial" method. While the "power" method was implemented, 
!   the "Newton" basis polynomial was not completed, as the "new" polynomial 
!   method, outlined in subroutine newpoly, replaced the method.  
!
!   References: 
!   
!   "Polynomial Preconditioned GMRES and GMRES-DR", 
!   Quan Liu, Ronald B Morgan, and Walter Wilcox, SIAM Journal on Sci. Comp. (2015) 
!
!===========================================================================
!   INPUT:
!
!       trye    = Even vector z to apply p(A) to
!       tryo    = Odd vector z to apply p(A) to
!       p       = Degree of polynomial
!       co      = Coefficients of polynomial p(A)
!       MRT     = MPIREALTYPE

!---------------------------------------------------------------------------
!   OUTPUT:
!
!       ye      = trye with polynomial preconditioner p(A)
!                 applied to it. (i.e. ye = p(A) * trye)
!       yo      = tryo with polynomial preconditioner p(A)
!                 applied to it. (i.e. yo = p(A) * tryo)
!---------------------------------------------------------------------------
!   FUTURE VERSIONS:
!
!       Only able to handle power basis, currently. Will be modified to use
!       switch between power basis and Newton basis in the future.
!
!---------------------------------------------------------------------------
!   MISC COMMENTS: 
!
!       Be careful in the application of co(:,i) to each order of M on b. 
!       The color index (i.e. the first index of ye and yo) goes over both 
!       real and complex elements, so it must be done with care when 
!       multiplying by co(:,i) which has both a real co(1,i) and imaginary 
!       co(2,i) part. This is where the "extra" negative sign comes from in 
!       front of trye and tryo when finding ye and yo. The first component 
!       looks like (for either the even or odd part), 
!
!              / try(1)+i*try(2) \
!        try = | try(3)+i*try(4) |, where i=sqrt(-1)
!              \ try(5)+i*try(6) /
!
!       so multiplying a complex number (a + ib) by v to get y gives,
!       (the . just indicates scalar multiplication)  
!
!                       / (a.try(1) - b.try(2)) + i (a.try(2) + b.try(1)) \
!        y = (a+ib).v = | (a.try(3) - b.try(4)) + i (a.try(4) + b.try(3)) | 
!                       \ (a.try(5) - b.try(6)) + i (a.try(6) + b.try(5)) /
!         
!                       / y(1)+i*y(2) \
!                     = | y(3)+i*y(4) |
!                       \ y(5)+i*y(6) /
!
!
!
!---------------------------------------------------------------------------


!use shift

!===========================Initializations==============================

real(kind=KR),         intent(out),    dimension(6,ntotal,4,2,8)   :: ye,yo
real(kind=KR),                         dimension(6,ntotal,4,2,8)   :: re,ro
integer(kind=KI),      intent(in)                                  :: p
real(kind=KR),                         dimension(6,ntotal,4,2,8,p) :: trye,tryo
real(kind=KR2),        intent(in),     dimension(2,p)              :: co
integer(kind=KI),      intent(in)                                  :: MRT
real(kind=KR2),                        dimension(6,ntotal,4,2,8)   :: tmpEE,tmpOO,&
                                                                      we,wo
integer(kind=KI)                                                   :: ieo,ibleo,&
                                                                      gblclr,i,k,&
                                                                      icri

!---------------------------NEEDED FOR HSINGLE()-------------------------
real(kind=KR),         intent(in),      dimension(:)                :: kappa
real(kind=KR),         intent(in),      dimension(:,:,:,:,:)        :: u
real(kind=KR),         intent(in),      dimension(:,:,:)            :: coact
integer(kind=KI),      intent(in),      dimension(:,:)              :: vecbl,&
                                                                       vecblinv
integer(kind=KI),      intent(in)                                   :: idag,myid
integer(kind=KI),      intent(in),      dimension(:,:)              :: nn,iblv
integer(kind=KI),      intent(in),      dimension(:)                :: bc,nms
integer(kind=KI),      intent(in),      dimension(:,:,:,:)          :: ib
logical,               intent(in),      dimension(:,:)              :: lbd
logical,               intent(in),      dimension(:)                :: ldiv
integer(kind=KI),      intent(in),      dimension(:,:,:)            :: lvbc
!------------------------------------------------------------------------



!===========================Main Algorithm===============================


    do icri=1,5,2
        do k=1,nvhalf

            ye(icri,k,:,:,:) = co(1,1)*trye(icri,k,:,:,:,1) &
                                -co(2,1)*trye(icri+1,k,:,:,:,1)
            ye(icri+1,k,:,:,:) = co(1,1)*trye(icri+1,k,:,:,:,1) &
                                +co(2,1)*trye(icri,k,:,:,:,1)

            yo(icri,k,:,:,:) = co(1,1)*tryo(icri,k,:,:,:,1) &
                                -co(2,1)*tryo(icri+1,k,:,:,:,1)
            yo(icri+1,k,:,:,:) = co(1,1)*tryo(icri+1,k,:,:,:,1) &
                                +co(2,1)*tryo(icri,k,:,:,:,1)

        enddo!k
    enddo!icri



    do i=1,p-1

        do ieo=1,2
            do ibleo=1,8
                call gammamult(trye(:,:,:,:,:,i),tryo(:,:,:,:,:,i),&
                               tmpEE,tmpOO,5,ieo,ibleo)
            enddo !ibleo
        enddo !ieo

        gblclr = 2

        call Hsingle(wo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        gblclr = 1

        call Hsingle(we,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)

        we = tmpEE - kappa(1) * we
        wo = tmpOO - kappa(1) * wo


        trye(:,:,:,:,:,i+1) = we(:,:,:,:,:)
        tryo(:,:,:,:,:,i+1) = wo(:,:,:,:,:)


        do icri=1,5,2
            do k=1,nvhalf

                ye(icri  ,k,:,:,:) = ye(icri ,k,:,:,:) &
                                    +co(1,i+1)*trye(icri,k,:,:,:,i+1) &
                                    -co(2,i+1)*trye(icri+1,k,:,:,:,i+1)
                ye(icri+1,k,:,:,:) = ye(icri+1,k,:,:,:) &
                                    +co(1,i+1)*trye(icri+1,k,:,:,:,i+1) &
                                    +co(2,i+1)*trye(icri,k,:,:,:,i+1)
                                    !y=P(A)*r


                yo(icri  ,k,:,:,:) = yo(icri ,k,:,:,:) &
                                    +co(1,i+1)*tryo(icri,k,:,:,:,i+1) &
                                    -co(2,i+1)*tryo(icri+1,k,:,:,:,i+1)
                yo(icri+1,k,:,:,:) = yo(icri+1,k,:,:,:) &
                                    +co(1,i+1)*tryo(icri+1,k,:,:,:,i+1) &
                                    +co(2,i+1)*tryo(icri,k,:,:,:,i+1)

            enddo!k
        enddo!icri
    enddo !i


!===========================End of Algorithm=============================


!!!ORIGINAL CODE FOR PRECONDITIONED GMRESDR

!    do icri=1,5,2
!        do k=1,nvhalf
!            y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
!                                -co(2,1)*try(icri+1,k,:,:,:,1)
!            y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
!                                +co(2,1)*try(icri,k,:,:,:,1)
!     enddo!k
!    enddo!icri

!!    print *,"original component1:real=",w1(1,1,1,1,1)
!!    print *,"original component2:imaginary=",w1(2,1,1,1,1)
!!    print *,"conditioned component1:real=",y(1,1,1,1,1)
!!    print *,"conditioned component2:imaginary=",y(2,1,1,1,1)


!    do i=1,p-1
!        call
!        Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
!                kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
!                lbd,iblv,MRT )  !z1=M*z

!        mvp = mvp + 1
!        do icri=1,5,2
!            do k=1,nvhalf
!                y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
!                                    +co(1,i+1)*try(icri,k,:,:,:,i+1) &
!                                    -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
!                y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
!                                    +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
!                                    +co(2,i+1)*try(icri,k,:,:,:,i+1)
!                                    !y=P(A)*r
!            enddo!k
!        enddo!icri
!    enddo!i

!!     do k=1,nvhalf
!!      r(:,k,:,:,:) = z2(:,k,:,:,:)
!!     enddo!k

!!!!!END OF ORIGINAL CODE FROM PRECONDITIONED GMRESDR


 endsubroutine applyppgamma5
 !!!!!!



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine gmresEIGritz(rwdir,be,bo,xe,xo,GMRES,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,hrv)

!===========================================================================
!   An implementation of gmres for obtaining ritz values to be used in 
!   polynomial preconditioning. Code copied and pasted from gmresdr5EIG
!   in inverters.f90 with gamma5 multiplication removed and hardcoded for 
!   a single cycle. This is for the non-hermitian system -PL
!
!    
!   Some things could be removed near the end since there is no restart and 
!   deflation, so it's wasted calculation. Single cycle implemented by 
!   simply removing while loop over main algorithm (i.e. part of algorithm 
!   constituting a single cycle) Solution vector x can be entirely removed 
!   as it does not update or return a solution vector, and only obtains 
!   harmonic ritz values. 
!===========================================================================

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: iflag, &
                                                             myid, MRT, MRT2, idag
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
    logical,          dimension(GMRES(1))                  :: myselect
    integer(kind=KI), dimension(GMRES(1))                  :: ipiv,ind
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv,rninit,rn,vn, vnf,xrn,resii,htest1,htest2
    real(kind=KR2),   dimension(GMRES(1))                  :: mag,dabbs
    real(kind=KR2),   dimension(2,GMRES(1))                :: ss, gs, gc, &
                                                               tau, w
    real(kind=KR2),   dimension(2,201,201)                  :: orthoge,orthog
    real(kind=KR2),   dimension(2,GMRES(1)+1)              :: c, c2, srv,d
    real(kind=KR2),   dimension(2,GMRES(1),1)              :: ff, punty
    real(kind=KR2),   dimension(2,GMRES(1))                :: dd, tmpVec
    real(kind=KR2),   dimension(2,GMRES(2))                :: rho
    real(kind=KR2),   dimension(GMRES(2))                  :: rna, sita
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(1)+1)  :: z,gon,rr,gondag
    real(kind=KR2),   dimension(2,GMRES(1),GMRES(1)+1)    :: hcht
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(1))    :: h, h2, h3, hprint, hh,g,gg,greal, &
                                                               tmpmat,hnew
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(2))    :: ws
    real(kind=KR2),   dimension(2,GMRES(2)+1,GMRES(2))    :: gca, gsa
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             ::re,ro,xte,xto,xre,xro
    real(kind=KR2),   dimension(6,ntotal,4,2,8,GMRES(1)+1) :: worke,worko
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpE, tmpE2
    real(kind=KR2),   dimension(6)              :: vte, vto
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpO, tmpO2
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpOO, tmpEE,xtmpEE,xtmpOO,xeTMP,xoTMP !BS
    real(kind=KR2), dimension(2) :: tmp1,tmp2,tmp3,tmp4,xtmp1,xtmp2,dote,doto!BS
    real(kind=KR2), dimension(2,GMRES(1))      :: vnj
   
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htempe,htmpo
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempe,wve,fe,wxe
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempo,wvo,fo,wxo!BS
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow,gblclr 
    integer(kind=KI) :: didmaindo, exists

    ! Getting the Harmonic ritz values 

    real(kind=KR2),  dimension(2,GMRES(1))                      :: th 
    real(kind=KR2),   intent(out), dimension(2,GMRES(1))        :: hrv 


if (myid==0) then
 print *,'just entered gmresEIGritz'
endif
   
!if (myid == 0) then
!   print *, 'the value of kappa in gmresEIGritz is = ',kappa(1)
!endif



        
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
  !xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2
  !xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

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



 if (myid==0) then
 print *,'starting while loop'
 endif




  ! begin cycles
 ! do while((rn(1)/rninit(1)) > rtol .AND. icycle <= kcyclim)
  !do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 2000))!becomes zero in 29 cycles
 ! do while(((rn(1)/rninit(1)) > 1e-9) .AND. (icycle <= 15))

!BEGIN MAIN ALGORITHM ! Commented out while loop so only a single cycle is performed -PL 1/12/21
!  do while((xrn(1) > 1e-9_KR2) .AND. (icycle <= 500)) !changed to this for residual norm print out TW 11/14/18                    
 ! do while(icycle <=150)! to make sure if it works
    ! begin gmres(m) iterations



!    if (myid == 0) then
!    print *, 'after main do while'
!    endif




    do while(j <= m)
      ! First build H_oe * v_e.
     
!!!!  Added line by BS. ordering gamama5 accordingly !!!!!
      
      !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 
      !do ieo=1,2
      !  do ibleo=1,8
      !    call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO, 5,ieo,ibleo)
      !  enddo
      !enddo
      
    
    tmpEE = ve(:,:,:,:,:,j) 
    tmpOO = vo(:,:,:,:,:,j) 

    
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




!if (myid==0) then
!print*, 'basis built'
!endif




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

if (.false.) then
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








    ! IGNORE SOLUTION UPDATE 
    if (.false.) then 


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

    !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 

    !  do ieo=1,2
    !    do ibleo=1,8
    !      call gammamult( xe(:,:,:,:,:), xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
    !    enddo
    !  enddo


    xtmpEE = xe(:,:,:,:,:) 
    xtmpOO = xo(:,:,:,:,:) 




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


      endif ! END SOLUTION UPDATE 


!if (myid == 0) then
!print *, 'before writing of residual norm'
!endif



if (.false.) then !to get residual norm of linear equations -TW 11/14/18



!if (myid==0) then
!print *, 'inside true/false'
!endif



        if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"linearresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=34,status="keep")
           endif
         endif



!if(myid==0) then
!print*, 'before write of residuals'
!endif



         if (myid==0) then
           open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=34,fmt="(i4,a6,es19.12)")icycle," ",xrn(1)
           close(unit=34,status="keep")
         endif
endif



if (.false.) then ! Ignore residual calculation 

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

    
    if (myid==0) then 
        print*, 'rn from single cycle gmres = ', rn(1) 
    endif 

endif



 


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

    if (myid == 0) then 
        print *, "Before linearsolver" 
    endif 

    call linearsolver(m,1,hcht,ipiv,ff)

    if (myid == 0) then 
        print *, "After linearsolver" 
    endif 


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


!
!allocate(hh(3,3))


    if (myid == 0) then 
        print *, "Before eigencalc" 
    endif 

    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg) !BS 1/19/2016
   ! call eigmodesLAPACK(hh,m,1,dd,gg)


    if (myid == 0) then 
        print *, "After eigencalc" 
    endif 


    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)

    !Saving the Harmonic Ritz values to output -PL 2-11-21

    hrv = 0.0_KR2

    do i=1,m
        hrv(:,i) = dd(:,ind(i))
    enddo


    ! Writing eigenvalues to file 

! if (myid==0) then
!   open(unit=12,file=trim(rwdir(myid+1))//"gmresEIGritz_EIGVALS.LOG",status="new",&
!        action="write",form="formatted",position="append")
!   write(unit=12,fmt='(a60)')"-------GMRES-DR (MATLAB ver) (full computed)---------"
!   write(unit=12,fmt='(a22,i4)') "               numEV: ",numEV
!   do eig=1,numEV
!     write(unit=12,fmt='(i5,a5,f25.15,a5,f25.15)')  eig,"     ",eigValExpand(1,eig), & 
!                                                        "     ",eigValExpand(2,eig)
!   enddo
!   close(unit=12,status="keep")
! endif




    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo





    !IGNORING RESTART AND DEFLATION 

    if (.false.) then 





    ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      call matvecmult(h,m,m,g(:,:,i),m,tmpVec) 
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

if (.false.) then

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
    if (.false.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,k
           
        !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21    
        !   do ieo=1,2           
        !    do ibleo=1,8
        !     call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
        !    enddo
        !   enddo


        tmpE = evectore(:,:,:,:,:,i) 
        tmpO = evectoro(:,:,:,:,:,i) 


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
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)

       !BS   write(unit=8,fmt="(a12,i9,es17.10)") "igmresdr5EIG",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif

    j = k+1
    icycle = icycle + 1
  !enddo !!END OF MAIN ALGORITHM !Removed while loop so only a single cycle is performed -PL 1/12/21

  
      
          if (myid==0) then !BS 
            do i =1,k 
                 write(*,fmt="(a12,i6,f19.11,a5,f19.11)") "evalueEIGritz:",i,evalue(1,i),"",evalue(2,i)
 !BS               print * , 'evalue5EIG:',evalue
            enddo 
         endif
            

          if (myid==0.AND.(icycle==150)) then !BS
             print * ,'rn/rninit = ',rn(1)/rninit(1)
          endif
          
  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle-1,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)
 
 !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr5EG",icycle-1,rn(1)/rninit(1)
    close(unit=8,status="keep")
  endif






endif ! ENCLOSING RESTART AND DEFLATION PORTION 





if(myid==0) then
print*, 'finished gmresEIGritz'
endif





end subroutine gmresEIGritz

!!!!!!



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine doublegmresEIGritz(rwdir,be,bo,xe,xo,GMRES,harv,u,GeeGooinv,& 
                               iflag,idag,kappa,coact,bc,vecbl,vecblinv,&
							   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2,hrv)

!===========================================================================
!   An implementation of gmres for obtaining ritz values to be used in 
!   polynomial preconditioning. Code copied and pasted from gmresdr5EIG
!   in inverters.f90 with gamma5 multiplication removed and hardcoded for 
!   a single cycle. This is for the non-hermitian system -PL
!
!    
!   Some things could be removed near the end since there is no restart and 
!   deflation, so it's wasted calculation. Single cycle implemented by 
!   simply removing while loop over main algorithm (i.e. part of algorithm 
!   constituting a single cycle) Solution vector x can be entirely removed 
!   as it does not update or return a solution vector, and only obtains 
!   harmonic ritz values. 
!===========================================================================

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: iflag, &
                                                             myid, MRT, MRT2, idag
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
    logical,          dimension(GMRES(1))                  :: myselect
    integer(kind=KI), dimension(GMRES(1))                  :: ipiv,ind
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv,rninit,rn,vn, vnf,xrn,resii,htest1,htest2
    real(kind=KR2),   dimension(GMRES(1))                  :: mag,dabbs
    real(kind=KR2),   dimension(2,GMRES(1))                :: ss, gs, gc, &
                                                               tau, w
    real(kind=KR2),   dimension(2,201,201)                  :: orthoge,orthog
    real(kind=KR2),   dimension(2,GMRES(1)+1)              :: c, c2, srv,d
    real(kind=KR2),   dimension(2,GMRES(1),1)              :: ff, punty
    real(kind=KR2),   dimension(2,GMRES(1))                :: dd, tmpVec
    real(kind=KR2),   dimension(2,GMRES(2))                :: rho
    real(kind=KR2),   dimension(GMRES(2))                  :: rna, sita
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(1)+1)  :: z,gon,rr,gondag
    real(kind=KR2),   dimension(2,GMRES(1),GMRES(1)+1)    :: hcht
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(1))    :: h, h2, h3, hprint, hh,g,gg,greal, &
                                                               tmpmat,hnew
    real(kind=KR2),   dimension(2,GMRES(1)+1,GMRES(2))    :: ws
    real(kind=KR2),   dimension(2,GMRES(2)+1,GMRES(2))    :: gca, gsa
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             ::re,ro,xte,xto,xre,xro
    real(kind=KR2),   dimension(6,ntotal,4,2,8,GMRES(1)+1) :: worke,worko
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpE, tmpE2
    real(kind=KR2),   dimension(6)              :: vte, vto
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpO, tmpO2
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpOO, tmpEE,xtmpEE,xtmpOO,xeTMP,xoTMP !BS
    real(kind=KR2), dimension(2) :: tmp1,tmp2,tmp3,tmp4,xtmp1,xtmp2,dote,doto!BS
    real(kind=KR2), dimension(2,GMRES(1))      :: vnj
   
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htempe,htmpo
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempe,wve,fe,wxe
    real(kind=KR), dimension(6,ntotal,4,2,8)  :: getempo,wvo,fo,wxo!BS
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow,gblclr 
    integer(kind=KI) :: didmaindo, exists

    ! Getting the Harmonic ritz values for phi_2(A) 

    real(kind=KR2),  dimension(2,GMRES(1))                      :: th 
    real(kind=KR2),   intent(out), dimension(2,GMRES(1))        :: hrv 
	
	! Harmonic ritz values for polynomial preconditioned operator phi_1(A) 
    real(kind=KR2),  intent(in),   dimension(2,GMRES(1))             :: harv 


if (myid==0) then
 print *,'just entered gmresEIGritz'
endif
   
!if (myid == 0) then
!   print *, 'the value of kappa in gmresEIGritz is = ',kappa(1)
!endif



        
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
  !xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2
  !xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

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



 if (myid==0) then
 print *,'starting while loop'
 endif




  ! begin cycles
 ! do while((rn(1)/rninit(1)) > rtol .AND. icycle <= kcyclim)
  !do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 2000))!becomes zero in 29 cycles
 ! do while(((rn(1)/rninit(1)) > 1e-9) .AND. (icycle <= 15))

!BEGIN MAIN ALGORITHM ! Commented out while loop so only a single cycle is performed -PL 1/12/21
!  do while((xrn(1) > 1e-9_KR2) .AND. (icycle <= 500)) !changed to this for residual norm print out TW 11/14/18                    
 ! do while(icycle <=150)! to make sure if it works
    ! begin gmres(m) iterations



!    if (myid == 0) then
!    print *, 'after main do while'
!    endif




    do while(j <= m)
      ! First build H_oe * v_e.
     
!!!!  Added line by BS. ordering gamama5 accordingly !!!!!
      
      !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 
      !do ieo=1,2
      !  do ibleo=1,8
      !    call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO, 5,ieo,ibleo)
      !  enddo
      !enddo
      
    
    tmpEE = ve(:,:,:,:,:,j) 
    tmpOO = vo(:,:,:,:,:,j) 

    
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



    ! 


    ! Replacing matrix application A with preconditioned system phi_1(A) 

    if (.false.) then 
    gblclr = 2

    call Hsingle(wvo,u,tmpEE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                      nms,lvbc,ib,lbd,iblv,MRT)
    ! Next build H_eo * v_o.
    gblclr = 1
    call Hsingle(wve,u,tmpOO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv,&
                     nms,lvbc,ib,lbd,iblv,MRT)


    wve = tmpEE - kappa(1) * wve
    wvo = tmpOO - kappa(1) * wvo
	  
    endif ! Falsed out single matrix application 

    


    call newphi(wve,wvo,tmpEE,tmpOO,m,harv,MRT,ieo,ibleo,gblclr,&
                kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                bc,nms,ib,lbd,ldiv,lvbc,0,MRT2) 
	  
	  
	  
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




!if (myid==0) then
!print*, 'basis built'
!endif




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

if (.false.) then
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








    ! IGNORE SOLUTION UPDATE 
    if (.false.) then 


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

    !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 

    !  do ieo=1,2
    !    do ibleo=1,8
    !      call gammamult( xe(:,:,:,:,:), xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
    !    enddo
    !  enddo


    xtmpEE = xe(:,:,:,:,:) 
    xtmpOO = xo(:,:,:,:,:) 




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


      endif ! END SOLUTION UPDATE 


!if (myid == 0) then
!print *, 'before writing of residual norm'
!endif



if (.false.) then !to get residual norm of linear equations -TW 11/14/18



!if (myid==0) then
!print *, 'inside true/false'
!endif



        if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"linearresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=34,status="keep")
           endif
         endif



!if(myid==0) then
!print*, 'before write of residuals'
!endif



         if (myid==0) then
           open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=34,fmt="(i4,a6,es19.12)")icycle," ",xrn(1)
           close(unit=34,status="keep")
         endif
endif



if (.false.) then ! Ignore residual calculation 

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

    
    if (myid==0) then 
        print*, 'rn from single cycle gmres = ', rn(1) 
    endif 

endif



 


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

    if (myid == 0) then 
        print *, "Before linearsolver" 
    endif 

    call linearsolver(m,1,hcht,ipiv,ff)

    if (myid == 0) then 
        print *, "After linearsolver" 
    endif 


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


!
!allocate(hh(3,3))


    if (myid == 0) then 
        print *, "Before eigencalc" 
    endif 

    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg) !BS 1/19/2016
   ! call eigmodesLAPACK(hh,m,1,dd,gg)


    if (myid == 0) then 
        print *, "After eigencalc" 
    endif 


    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)

    !Saving the Harmonic Ritz values to output -PL 2-11-21

    hrv = 0.0_KR2

    do i=1,m
        hrv(:,i) = dd(:,ind(i))
    enddo


    ! Writing eigenvalues to file 

! if (myid==0) then
!   open(unit=12,file=trim(rwdir(myid+1))//"gmresEIGritz_EIGVALS.LOG",status="new",&
!        action="write",form="formatted",position="append")
!   write(unit=12,fmt='(a60)')"-------GMRES-DR (MATLAB ver) (full computed)---------"
!   write(unit=12,fmt='(a22,i4)') "               numEV: ",numEV
!   do eig=1,numEV
!     write(unit=12,fmt='(i5,a5,f25.15,a5,f25.15)')  eig,"     ",eigValExpand(1,eig), & 
!                                                        "     ",eigValExpand(2,eig)
!   enddo
!   close(unit=12,status="keep")
! endif




    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo





    !IGNORING RESTART AND DEFLATION 

    if (.false.) then 





    ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      call matvecmult(h,m,m,g(:,:,i),m,tmpVec) 
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

if (.false.) then

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
    if (.false.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,k
           
        !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21    
        !   do ieo=1,2           
        !    do ibleo=1,8
        !     call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
        !    enddo
        !   enddo


        tmpE = evectore(:,:,:,:,:,i) 
        tmpO = evectoro(:,:,:,:,:,i) 


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
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)

       !BS   write(unit=8,fmt="(a12,i9,es17.10)") "igmresdr5EIG",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif

    j = k+1
    icycle = icycle + 1
  !enddo !!END OF MAIN ALGORITHM !Removed while loop so only a single cycle is performed -PL 1/12/21

  
      
          if (myid==0) then !BS 
            do i =1,k 
                 write(*,fmt="(a12,i6,f19.11,a5,f19.11)") "evalueEIGritz:",i,evalue(1,i),"",evalue(2,i)
 !BS               print * , 'evalue5EIG:',evalue
            enddo 
         endif
            

          if (myid==0.AND.(icycle==150)) then !BS
             print * ,'rn/rninit = ',rn(1)/rninit(1)
          endif
          
  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmresEIGritz",icycle-1,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)
 
 !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr5EG",icycle-1,rn(1)/rninit(1)
    close(unit=8,status="keep")
  endif






endif ! ENCLOSING RESTART AND DEFLATION PORTION 





if(myid==0) then
print*, 'finished gmresEIGritz'
endif





end subroutine doublegmresEIGritz

!!!!!!




!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine gmres5EIGritz(rwdir,be,bo,xe,xo,GMRES,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,hrv)

!===========================================================================
!   An implementation of gmres for obtaining harmonic ritz values to be used in 
!   polynomial preconditioning. Code copied and pasted from gmresdr5EIG
!   in inverters.f90 with gamma5 multiplication left in and hardcoded for 
!   a single cycle. This is for the hermitian system  -PL
!
!    
!   Some things could be removed near the end since there is no restart and 
!   deflation, so it's wasted calculation. Single cycle implemented by 
!   simply removing while loop over main algorithm (i.e. part of algorithm 
!   constituting a single cycle) Solution vector x can be entirely removed 
!   as it does not update or return a solution vector, and only obtains 
!   harmonic ritz values. 
!===========================================================================

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: iflag, &
                                                             myid, MRT, MRT2, idag
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
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: dd, tmpVec
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

    ! Getting the Harmonic ritz values 

    real(kind=KR2),  dimension(2,nmaxGMRES)                 :: th 
    real(kind=KR2),   intent(out), dimension(2,GMRES(1))        :: hrv 


if (myid==0) then
 print *,'just entered gmres5EIGritz'
endif
   
!if (myid == 0) then
!   print *, 'the value of kappa in gmresEIGritz is = ',kappa(1)
!endif



        
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
  !xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2
  !xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

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



 if (myid==0) then
 print *,'starting gmres5EIGritz while loop'
 endif




  ! begin cycles
 ! do while((rn(1)/rninit(1)) > rtol .AND. icycle <= kcyclim)
  !do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 2000))!becomes zero in 29 cycles
 ! do while(((rn(1)/rninit(1)) > 1e-9) .AND. (icycle <= 15))

!BEGIN MAIN ALGORITHM ! Commented out while loop so only a single cycle is performed -PL 1/12/21
!  do while((xrn(1) > 1e-9_KR2) .AND. (icycle <= 500)) !changed to this for residual norm print out TW 11/14/18                    
 ! do while(icycle <=150)! to make sure if it works
    ! begin gmres(m) iterations



!    if (myid == 0) then
!    print *, 'after main do while'
!    endif




    do while(j <= m)
      ! First build H_oe * v_e.
     
!!!!  Added line by BS. ordering gamama5 accordingly !!!!!
      
      !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 
      do ieo=1,2
        do ibleo=1,8
          call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO, 5,ieo,ibleo)
        enddo
      enddo
      
    
    !tmpEE = ve(:,:,:,:,:,j) 
    !tmpOO = vo(:,:,:,:,:,j) 

    
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




!if (myid==0) then
!print*, 'basis built'
!endif




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

if (.false.) then
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







    ! IGNORE SOLUTION UPDATE
    if (.false.) then 


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

    !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 

      do ieo=1,2
        do ibleo=1,8
          call gammamult( xe(:,:,:,:,:), xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
        enddo
      enddo


    !xtmpEE = xe(:,:,:,:,:) 
    !xtmpOO = xo(:,:,:,:,:) 




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

      
      endif ! END SOLUTION UPDATE 


!if (myid == 0) then
!print *, 'before writing of residual norm'
!endif



if (.false.) then !to get residual norm of linear equations -TW 11/14/18



!if (myid==0) then
!print *, 'inside true/false'
!endif
 
        if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"linearresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=34,status="keep")
           endif
         endif



!if(myid==0) then
!print*, 'before write of residuals'
!endif



         if (myid==0) then
           open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=34,fmt="(i4,a6,es19.12)")icycle," ",xrn(1)
           close(unit=34,status="keep")
         endif
endif



if (.false.) then ! IGNORE RESIDUAL CALCULATION 

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


endif ! END RESIDUAL CALCULATION 

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


!
!allocate(hh(3,3))



    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg) !BS 1/19/2016
   ! call eigmodesLAPACK(hh,m,1,dd,gg)

    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)

    !Saving the Harmonic Ritz values to output -PL 2-11-21

    hrv = 0.0_KR2

    do i=1,m
        hrv(:,i) = dd(:,ind(i))
    enddo

    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo



    ! IGNORING RESTART AND DEFLATION 
    if (.false.) then 







    ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      call matvecmult(h,m,m,g(:,:,i),m,tmpVec) 
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

if (.false.) then

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
    if (.false.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,k
           
        !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21    
           do ieo=1,2           
            do ibleo=1,8
             call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
            enddo
           enddo


        !tmpE = evectore(:,:,:,:,:,i) 
        !tmpO = evectoro(:,:,:,:,:,i) 


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
        write(unit=8,fmt="(a13,i4,a6,es11.4,a6,es11.4,es17.10)") "gmres5EIGritz",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
        write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "gmres5EIGritz",icycle,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)

       !BS   write(unit=8,fmt="(a12,i9,es17.10)") "igmresdr5EIG",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif

    j = k+1
    icycle = icycle + 1
  !enddo !!END OF MAIN ALGORITHM !Removed while loop so only a single cycle is performed -PL 1/12/21

  
      
          if (myid==0) then !BS 
            do i =1,k 
                 write(*,fmt="(a15,i6,f19.11,a5,f19.11)") "evalue5EIGritz:",i,evalue(1,i),"",evalue(2,i)
 !BS               print * , 'evalue5EIG:',evalue
            enddo 
         endif
            

          if (myid==0.AND.(icycle==150)) then !BS
             print * ,'rn/rninit = ',rn(1)/rninit(1)
          endif
          
  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
          write(unit=8,fmt="(a13,i4,a6,es11.4,a6,es11.4,es17.10)") "gmres5EIGritz",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
          write(unit=8,fmt="(a13,i4,a6,es11.4,a6,es11.4,es17.10)") "gmres5EIGritz",icycle-1,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)
 
 !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr5EG",icycle-1,rn(1)/rninit(1)
    close(unit=8,status="keep")
  endif





endif ! END RESTART AND DEFLATION PORTION 




if(myid==0) then
print*, 'finished gmres5EIGritz'
endif





end subroutine gmres5EIGritz

!!!!!!





!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine doublegmres5EIGritz(rwdir,be,bo,xe,xo,GMRES,harv,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,hrv)

!===========================================================================
!   An implementation of gmres for obtaining harmonic ritz values to be used in 
!   polynomial preconditioning. Code copied and pasted from gmresdr5EIG
!   in inverters.f90 with gamma5 multiplication left in and hardcoded for 
!   a single cycle. This is for the hermitian system  -PL
!
!    
!   Some things could be removed near the end since there is no restart and 
!   deflation, so it's wasted calculation. Single cycle implemented by 
!   simply removing while loop over main algorithm (i.e. part of algorithm 
!   constituting a single cycle) Solution vector x can be entirely removed 
!   as it does not update or return a solution vector, and only obtains 
!   harmonic ritz values. 
!===========================================================================

    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: be,bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe,xo
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: iflag, &
                                                             myid, MRT, MRT2, idag
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
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: dd, tmpVec
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

    ! Getting the Harmonic ritz values 

    real(kind=KR2),  dimension(2,nmaxGMRES)                 :: th 
    real(kind=KR2),  intent(out), dimension(2,GMRES(1))     :: hrv 

    ! Harmonic ritz values for polynomial preconditioned operator phi_1(Agamma5) 
    real(kind=KR2),  intent(in),   dimension(2,GMRES(1))    :: harv 


if (myid==0) then
 print *,'just entered doublegmres5EIGritz'
endif
   
!if (myid == 0) then
!   print *, 'the value of kappa in gmresEIGritz is = ',kappa(1)
!endif



        
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
  !xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2
  !xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2

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



 if (myid==0) then
 print *,'starting doublegmres5EIGritz while loop'
 endif




  ! begin cycles
 ! do while((rn(1)/rninit(1)) > rtol .AND. icycle <= kcyclim)
  !do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 2000))!becomes zero in 29 cycles
 ! do while(((rn(1)/rninit(1)) > 1e-9) .AND. (icycle <= 15))

!BEGIN MAIN ALGORITHM ! Commented out while loop so only a single cycle is performed -PL 1/12/21
!  do while((xrn(1) > 1e-9_KR2) .AND. (icycle <= 500)) !changed to this for residual norm print out TW 11/14/18                    
 ! do while(icycle <=150)! to make sure if it works
    ! begin gmres(m) iterations



!    if (myid == 0) then
!    print *, 'after main do while'
!    endif




    do while(j <= m)
      ! First build H_oe * v_e.
     
!!!!  Added line by BS. ordering gamama5 accordingly !!!!!
      



    ! Replacing single matrix application Ag5 with preconditioned system phi_1(Ag5) 

    if (.false.) then 


      !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 
      do ieo=1,2
        do ibleo=1,8
          call gammamult( ve(:,:,:,:,:,j), vo(:,:,:,:,:,j),tmpEE,tmpOO, 5,ieo,ibleo)
        enddo
      enddo
      
    
    !tmpEE = ve(:,:,:,:,:,j) 
    !tmpOO = vo(:,:,:,:,:,j) 

    
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


    endif ! Falsed out single matrix application 


    tmpEE = ve(:,:,:,:,:,j) 
    tmpOO = vo(:,:,:,:,:,j) 


    ! Applying phi_1(Ag5) instead of Ag5 to v 
    call newphi(wve,wvo,tmpEE,tmpOO,m,harv,MRT,ieo,ibleo,gblclr,&
                kappa,u,coact,vecbl,vecblinv,idag,myid,nn,iblv,&
                bc,nms,ib,lbd,ldiv,lvbc,1,MRT2) 





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




!if (myid==0) then
!print*, 'basis built'
!endif




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

if (.false.) then
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







    ! IGNORE SOLUTION UPDATE
    if (.false.) then 


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

    !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21 

      do ieo=1,2
        do ibleo=1,8
          call gammamult( xe(:,:,:,:,:), xo(:,:,:,:,:),xtmpEE,xtmpOO,5,ieo,ibleo)
        enddo
      enddo


    !xtmpEE = xe(:,:,:,:,:) 
    !xtmpOO = xo(:,:,:,:,:) 




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

      
      endif ! END SOLUTION UPDATE 


!if (myid == 0) then
!print *, 'before writing of residual norm'
!endif



if (.false.) then !to get residual norm of linear equations -TW 11/14/18



!if (myid==0) then
!print *, 'inside true/false'
!endif
 
        if (myid==0) then
          inquire(file=trim(rwdir(myid+1))//"linearresidual.dat", exist=exists)
           if (.not. exists) then
            open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="new",&
            action="write",form="formatted")
            close(unit=34,status="keep")
           endif
         endif



!if(myid==0) then
!print*, 'before write of residuals'
!endif



         if (myid==0) then
           open(unit=34,file=trim(rwdir(myid+1))//"linearresidual.dat",status="old",action="write",&
           form="formatted",position="append")
           write(unit=34,fmt="(i4,a6,es19.12)")icycle," ",xrn(1)
           close(unit=34,status="keep")
         endif
endif



if (.false.) then ! IGNORE RESIDUAL CALCULATION 

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


endif ! END RESIDUAL CALCULATION 

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


!
!allocate(hh(3,3))



    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg) !BS 1/19/2016
   ! call eigmodesLAPACK(hh,m,1,dd,gg)

    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)

    !Saving the Harmonic Ritz values to output -PL 2-11-21

    hrv = 0.0_KR2

    do i=1,m
        hrv(:,i) = dd(:,ind(i))
    enddo

    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo




    ! IGNORING RESTART AND DEFLATION 
    if (.false.) then 







    ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      call matvecmult(h,m,m,g(:,:,i),m,tmpVec) 
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

if (.false.) then

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
    if (.false.) then
       if (icycle==1 .OR. icycle==10 .OR. icycle==100 .OR. icycle==1000 .OR. icycle==2000) then !change to .false. if you dont need/want it
        do i = 1,k
           
        !! Commented out from gmresdr5EIG to remove Hermitian forcing -PL 1/12/21    
           do ieo=1,2           
            do ibleo=1,8
             call gammamult(evectore(:,:,:,:,:,i),evectoro(:,:,:,:,:,i), tmpE, tmpO,5,ieo,ibleo)
            enddo
           enddo


        !tmpE = evectore(:,:,:,:,:,i) 
        !tmpO = evectoro(:,:,:,:,:,i) 


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
        write(unit=8,fmt="(a19,i4,a6,es11.4,a6,es11.4,es17.10)") "doublegmres5EIGritz",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
        write(unit=8,fmt="(a19,i4,a6,es11.4,a6,es11.4,es17.10)") "doublegmres5EIGritz",icycle,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)

       !BS   write(unit=8,fmt="(a12,i9,es17.10)") "igmresdr5EIG",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif

    j = k+1
    icycle = icycle + 1
  !enddo !!END OF MAIN ALGORITHM !Removed while loop so only a single cycle is performed -PL 1/12/21

  
      
          if (myid==0) then !BS 
            do i =1,k 
                 write(*,fmt="(a21,i6,f19.11,a5,f19.11)") "evaluedouble5EIGritz:",i,evalue(1,i),"",evalue(2,i)
 !BS               print * , 'evalue5EIG:',evalue
            enddo 
         endif
            

          if (myid==0.AND.(icycle==150)) then !BS
             print * ,'rn/rninit = ',rn(1)/rninit(1)
          endif
          
  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
          write(unit=8,fmt="(a19,i4,a6,es11.4,a6,es11.4,es17.10)") "doublegmres5EIGritz",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
          write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)") "doublegmres5EIGritz",icycle-1,"xrn=",xrn(1),"rnit=",rninit(1),xrn(1)/rninit(1)
 
 !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr5EG",icycle-1,rn(1)/rninit(1)
    close(unit=8,status="keep")
  endif





endif ! END RESTART AND DEFLATION PORTION 




if(myid==0) then
print*, 'finished gmres5EIGritz'
endif





end subroutine doublegmres5EIGritz

!!!!!!





!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




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


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



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




!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




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




 end module preconditioning

