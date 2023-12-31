! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! pseudolapack.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! A collection of linear algebra routines.
!
! The code is adapted from f77 routines in LAPACK.
! The official reference for LAPACK is the LAPACK Users' Guide Third Edition,
! updated 22 Aug 1999.  It is available at
! http://www.netlib.org/lapack/lug/lapack_lug.html .
! The authors are Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra,
! J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney and D. Sorensen.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module pseudolapack

    use kinds
    use latdims
    implicit none
    private

    real(kind=KR2), private, parameter :: dlamchs=2.3e-308_KR2, &
                                          dlamche=1.2e-16_KR2, &
                                          dlamchp=2.3e-16_KR2, &
                                          dlamch1=2.0e-146_KR2, &
                                          dlamch2=dlamch1

! Define access to subroutines and functions.
     public  :: linearsolver, eigencalc, svdcalc, orgschur, hessenberg, &
                qgenerator, evalues, twonorm,indexx
     private :: LUfactor, schurswap, balance, Tvectors, evectors, &
                bidiagonalize, bidiagqgen, qgenkernel, svdkernel, QRdouble, &
                makereflector, multreflector, mvutil1, mvutil2, mvutil3, &
                mvutil4, mvutil5, mvutil6, mmutil1, mmutil2, rowswap, &
                svd2by2, sing2by2, rotate2, rotatecomp, compdiv, pythag2, &
                pythag3, safesqrt

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine linearsolver(n,nrhs,a,ipiv,b)
! Compute the solution to a complex system of linear equations a*x = b,
! where "a" is an n-by-n matrix and x and b are n-by-nrhs matrices.
! The LU decomposition with partial pivoting and row interchanges is used to
! factor "a" as a=P*L*U, where P is a permutation matrix, L is unit lower
! triangular, and U is upper triangular.  The factored form of "a" is then
! used to solve the system of equations a*x=b.
! This subroutine corresponds to SUBROUTINE ZGESV of LAPACK.
! INPUT:
!   n is the number of linear equations, i.e., the order of the matrix "a".
!   nrhs is the number of right hand sides, i.e., the number of columns of the
!        matrix b.
!   a(1,:,:),a(2,:,:) are input as the Re,Im parts of the n by n matrix.
!   b(1,:,:),b(2,:,:) are input as the Re,Im parts of the n by nrhs matrix of b.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the factors L and U from the factorization a=P*L*U;
!                     the unit diagonal elements of L are not stored.
!   b(1,:,:),b(2,:,:) are the Re,Im parts of the n by nrhs solution matrix of x.
!   ipiv contains the pivot indices that define the permutation matrix P;
!        row i of the matrix was interchanged with row ipiv(i).

    integer(kind=KI), intent(in)                      :: n, nrhs
    integer(kind=KI), intent(out),   dimension(:)     :: ipiv
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a, b

! Compute the LU factorization of "a".
    call LUfactor(n,n,a,ipiv)
! Solve the system a*x=b, overwriting b with x.
    call mvutil6(n,nrhs,a,ipiv,b)

 end subroutine linearsolver

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine eigencalc(a,n,ivec,w,vr)
! Compute, for an n by n complex nonsymmetric matrix "a", the eigenvalues and,
! optionally, the right eigenvectors.  The computed eigenvectors are normalized
! to have Euclidean norm equal to 1 and largest component real.
! This subroutine corresponds to SUBROUTINE ZGEEV of LAPACK.
! INPUT:
!   a(1,:,:),a(2,:,:) are input as the Re,Im parts of the n by n matrix.
!                     expected size: a(2,n,n)
!   n is the order of the matrix.  It must be positive.
!   ivec=1 to compute eigenvectors, otherwise they are not computed.
! OUTPUT:
!   a() gets destroyed on output.
!   w(1,:),w(2,:) contain the Re,Im parts of the eigenvalues.
!                 expected size: w(2,n)
!   vr(1,:,:),vr(2,:,:) are the Re,Im parts of the right eigenvectors if ivec=1.
!                       In particular, vr(:,:,j) is the eigenvector
!                       corresponding to eigenvalue w(:,j).
!                       expected size: vr(2,n,n)

    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n, ivec
    real(kind=KR2),   intent(out),   dimension(:,:)   :: w
    real(kind=KR2),   intent(out),   dimension(:,:,:) :: vr

    integer(kind=KI)                         :: i, j, k, ilo, ihi, myid
    real(kind=KR2)                           :: scl, dmax, vrRe, time1, time2
    real(kind=KR2), dimension(2)             :: tmp
    real(kind=KR2), dimension(2,2*nmaxGMRES) :: workspace
    real(kind=KR2), dimension(nmaxGMRES)     :: ascale
    real(kind=KR2), dimension(2,nmaxGMRES)   :: tau

!if(myid==0) then !TW 2/19/17
 !  print *, 'Entering eigencalc' !TW
!endif !TW

! Balance the matrix.
    call balance(a,n,ilo,ihi,ascale)!BS 1/13/2016 Dr Morgn''s suggestion
! Reduce the matrix to upper Hessenberg form.
!ilo = 1
!ihi = 200
    call hessenberg(a,n,ilo,ihi,tau,workspace)

! Compute eigenvalues and (if requested) eigenvectors.
    if (ivec==0) then
     call evalues(ivec,n,ilo,ihi,a,w,vr,workspace)
    elseif (ivec==1) then
     do j = 1,n
      do i = j,n
       vr(:,i,j) = a(:,i,j)
      enddo ! i
     enddo ! j
!-generate the required complex unitary matrix Q (here stored as vr).
     call qgenerator(n,ilo,ihi,vr,tau,workspace)
    ! if (myid==0) then !TW added time and print statements lines 115-122
     !call cpu_time(time1)!3/2/17
    ! endif
     call evalues(ivec,n,ilo,ihi,a,w,vr,workspace)
    ! if(myid==0) then
     !call cpu_time(time2)
      ! print *, 'Time for subroutine evalues is', time2-time1,' seconds'
    ! endif
!-compute left and/or right eigenvectors.
     call Tvectors(n,a,vr,workspace,tau(1,:))
!-undo balancing of right eigenvectors.
     call evectors(n,ilo,ihi,ascale,vr)
!-normalize right eigenvectors and make largest component real.
     do i = 1,n
      scl = 1.0_KR2/twonorm(vr(:,:,i),n)
      vr(:,:,i) = scl*vr(:,:,i)
      do k = 1,n
       if (abs(vr(1,k,i))<dlamch2) vr(1,k,i)=0.0_KR2
       if (abs(vr(2,k,i))<dlamch2) vr(2,k,i)=0.0_KR2
       ascale(k) = vr(1,k,i)**2 + vr(2,k,i)**2
      enddo ! k
      k = 1
      if (n>1) then
       dmax = ascale(1)
       do j = 2,n
        if (ascale(j)>dmax) then
         k = j
         dmax = ascale(j)
        endif
       enddo ! j
      endif
      tmp(1) = vr(1,k,i)/sqrt(ascale(k))
      tmp(2) = -vr(2,k,i)/sqrt(ascale(k))
      do j = 1,n
       vrRe = tmp(1)*vr(1,j,i) - tmp(2)*vr(2,j,i)
       vr(2,j,i) = tmp(1)*vr(2,j,i) + tmp(2)*vr(1,j,i)
       vr(1,j,i) = vrRe
      enddo ! j
      vr(2,k,i) = 0.0_KR2
     enddo ! i
    endif

 end subroutine eigencalc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine svdcalc(m,n,a,s,u)
! Compute the singular value decomposition (SVD) of a complex m by n matrix "a"
! and the left singular vectors. The SVD is written: a = u * s * v^dagger
! where s is an m by n matrix which is zero except for its min(m,n) diagonal
! elements, u is an m by m unitary matrix, and v is an n by n unitary matrix.
! The diagonal elements of s are the singular values of "a"; they are real and
! non-negative, and are returned in descending order.  The first min(m,n)
! columns of u are the left singular vectors of "a".
! This subroutine corresponds to SUBROUTINE ZGESVD of LAPACK.
! INPUT:
!   m is the number of rows in the input matrix "a".
!   n is the number of columns in the input matrix "a".
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the input matrix.
!                     expected size: a(2,m,n)
! OUTPUT:
!   a() is destroyed on output.
!   s() contains the singular values of "a", sorted so that s(i) >= s(i+1).
!   u(1,:,:),u(2,:,:) are the Re,Im parts of the m by m unitary matrix u.

    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    real(kind=KR2),   intent(out),   dimension(:)     :: s
    real(kind=KR2),   intent(out),   dimension(:,:,:) :: u

    character(len=1)                         :: uplo
    integer(kind=KI)                         :: i, j, irwork
    real(kind=KR2)                           :: eps, smlnum, bignum
    real(kind=KR2), dimension(nmaxGMRES)     :: rwork
    real(kind=KR2), dimension(2,2*nmaxGMRES) :: work

! Set the constants to control overflow and underflow.
    eps = dlamchp
    smlnum = sqrt(dlamchs)/eps
    bignum = 1.0_KR2/smlnum

! "a" has at least as many rows as columns.
    if (m>=n) then
!-bidiagonalize "a".
     call bidiagonalize(m,n,a,s,rwork,work(:,1:n),work(:,n+1:2*n))
!-copy result to u and generate left bidiagonalizing vectors in u.
     do j = 1,n
      do i = j,m
       u(:,i,j) = a(:,i,j)
      enddo ! i
     enddo ! j
     call bidiagqgen(m,n,u,work(:,1:n),work(:,n+1:2*n))
!-perform bidiagonal QR iteration, computing left singular vectors in u and
! computing right singular vectors in vt.
     uplo = "U"
     irwork = 1 + n
     call svdkernel(uplo,n,m,s,rwork(1:irwork-1),u,rwork(irwork:))

    else
! "a" has more columns than rows.
!-bidiagonalize "a".
     call bidiagonalize(m,n,a,s,rwork,work(:,1:m),work(:,m+1:2*m))
!-copy result to u and generate left bidiagonalizing vectors in u.
     do j = 1,m
      do i = j,m
       u(:,i,j) = a(:,i,j)
      enddo ! i
     enddo ! j
     call bidiagqgen(m,n,u,work(:,1:m),work(:,m+1:2*m))
!-perform bidiagonal QR iteration, computing left singular vectors in u and
! computing right singular vectors in vt.
     uplo = "L"
     irwork = 1 + m
     call svdkernel(uplo,m,m,s,rwork(1:irwork-1),u,rwork(irwork:))
    endif

 end subroutine svdcalc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine orgschur(eselect,n,t,q,w,m)
! Reorder the Schur factorization of a complex matrix A = Q*T*Q^H, so that a
! selected cluster of eigenvalues appears in the leading positions on the
! diagonal of the upper triangular matrix T, and the leading columns of Q form
! an orthonormal basis of the corresponding right invariant subspace.
! This subroutine corresponds to SUBROUTINE ZTRSEN of LAPACK.
! INPUT:
!   eselect specifies the eigenvalues in the selected cluster.  To select the
!           j'th eigenvalue, eselect(j) must be set to .true..
!   n is the order of the matrix T.
!   t is the upper triangular matrix T.
!   q is the matrix Q of Schur vectors.
! OUTPUT:
!   t is overwritten by the reordered matrix with selected eigenvalues as the
!     leading diagonal elements.
!   q is its input value postmultiplied by the unitary transformation matrix
!     which reorders t; the leading m columns of Q form an orthonormal basis
!     for the specified invariant subspace.
!   w contains the reordered eigenvalues of T, in the same order as they appear
!     on the diagonal of T.
!   m is the dimension of the specified invariant subspace (0 <= m <= n).

    logical,          intent(in),    dimension(:)     :: eselect
    integer(kind=KI), intent(in)                      :: n
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: t, q
    real(kind=KR2),   intent(out),   dimension(:,:)   :: w
    integer(kind=KI), intent(out)                     :: m

    integer(kind=KI) :: k, ks

! Set m to the number of selected eigenvalues.
    m = 0
    do k = 1,n
     if (eselect(k)) m=m+1
    enddo ! k

    if (m/=n .and. m/=0) then
! Collect the selected eigenvalues at the top left corner of T.
     ks = 0
     do k = 1,n
      if (eselect(k)) then
       ks = ks + 1
! Swap the k'th eigenvalue to position ks.
       if (k/=ks) call schurswap(n,t,q,k,ks)
      endif
     enddo ! k
    endif

! Copy reordered eigenvalues to w.
    do k = 1,n
     w(:,k) = t(:,k,k)
    enddo ! k

 end subroutine orgschur

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine LUfactor(m,n,a,ipiv)
! Compute an LU factorization of a general m by n matrix "a" using partial
! pivoting with row interchanges.  The factorization has the form a=P*L*U
! where P is a permutation matrix, L is lower triangular with unit diagonal
! elements (lower trapezoidal if m > n), and U is upper triangular (upper
! trapezoidal if m < n).
! This subroutine corresponds to SUBROUTINE ZGETRF and ZGETF2 of LAPACK.
! INPUT:
!   m is the number of rows of the matrix "a".
!   n is the number of columns of the matrix "a".
!   a(1,:,:),a(2,:,:) are input as the Re,Im parts of the initial m by n matrix.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the factors L and U from the factorization a=P*L*U;
!                     the unit diagonal elements of L are not stored.
!   ipiv contains the pivot indices; for 1 <= i <= min(m,n), row i of the
!        matrix was interchanged with row ipiv(i).

    integer(kind=KI), intent(in)                      :: m, n
    integer(kind=KI), intent(out),   dimension(:)     :: ipiv
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a

    integer(kind=KI)             :: i, j, jp, row, col
    real(kind=KR)                :: smax, stry
    real(kind=KR2), dimension(2) :: minusone, atemp

! A useful definition.
    minusone(1) = -1.0_KR2
    minusone(2) = 0.0_KR2

    do j = 1,min(m,n)
! Find pivot and test for singularity.
     jp = j
     if (m>j) then
      smax = abs(a(1,j,j)) + abs(a(2,j,j))
      do i = j+1,m
       stry = abs(a(1,i,j)) + abs(a(2,i,j))
       if (stry>smax) then
        jp = i
        smax = stry
       endif
      enddo ! i
     endif
     ipiv(j) = jp
     if (a(1,jp,j)/=0.0_KR2 .or. a(2,jp,j)/=0.0_KR2) then
! Apply the interchange to columns 1:n.
      if (jp/=j) then
       do i = 1,n
        atemp(:) = a(:,j,i)
        a(:,j,i) = a(:,jp,i)
        a(:,jp,i) = atemp(:)
       enddo ! i
      endif
! Compute elements j+1:m of j'th column.
      if (j<m) then
       atemp(1) = a(1,j,j)/(a(1,j,j)**2+a(2,j,j)**2)
       atemp(2) = -a(2,j,j)/(a(1,j,j)**2+a(2,j,j)**2)
       do i = j+1,m
        stry = atemp(1)*a(1,i,j) - atemp(2)*a(2,i,j)
        a(2,i,j) = atemp(1)*a(2,i,j) + atemp(2)*a(1,i,j)
        a(1,i,j) = stry
       enddo ! i
      endif
     else
      open(unit=8,file="PSEUDOLAPACK.ERROR",action="write",status="replace" &
          ,form="formatted")
       write(unit=8,fmt=*) "Singularity in LUfactor."
      close(unit=8,status="keep")
      stop
     endif
! Update trailing submatrix.
     row = m - j
     col = n - j
     if (j<min(m,n)) &
      call mvutil5(row,col,minusone,a(:,j+1:m,j),a(:,j,j+1:n),a(:,j+1:m,j+1:n))
    enddo ! j

 end subroutine LUfactor

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine schurswap(n,t,q,ifst,ilst)
! Reorder the Schur factorization of a complex matrix A = Q*T*Q^H, so that the
! diagonal element of T with row index ifst is moved to row ilst.  The Schur
! form T is reordered by a unitary similarity transformation Z^H*T*Z, and
! optionally the matrix Q of Schur vectors is updated by postmultplying it
! with Z.
! This subroutine corresponds to SUBROUTINE ZTREXC of LAPACK.
! INPUT:
!   n is the order of the matrix T.
!   t is the upper triangular matrix T.
!   q is the matrix Q od Schur vectors.
!   ifst,ilst: The element with row index IFST is moved to row ILST by a
!              sequence of transpositions between adjacent elements.
! OUTPUT:
!   t is the reordered upper triangular matrix.
!   q is its input value postmultiplied by the unitary transformation matrix Z
!     which reorders T.

    integer(kind=KI), intent(in)                      :: n, ifst, ilst
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: t, q

    integer(kind=KI)             :: i, k, m1, m2, m3
    real(kind=KR2)               :: cs
    real(kind=KR2), dimension(2) :: sn, t11, t22, temp

    if (n>1 .and. ifst/=ilst) then
     if (ifst<ilst) then
! Move the ifst'th diagonal element forward down the diagonal.
      m1 = 0
      m2 = -1
      m3 = 1
     else
! Move the ifst'th diagonal element backward up the diagonal.
      m1 = -1
      m2 = 0
      m3 = -1
     endif
     do k = ifst+m1,ilst+m2,m3
! Interchange the k'th and (k+1)'th diagonal elements.
      t11(:) = t(:,k,k)
      t22(:) = t(:,k+1,k+1)
! Determine the transformation to perform the interchange.
      temp = t22 - t11
      call rotatecomp(t(:,k,k+1),temp,cs,sn,temp)
! Apply transformation to the matrix T.
      if (k+2<=n) then
       do i = k+2,n
        temp(1) = cs*t(1,k,i) + sn(1)*t(1,k+1,i) - sn(2)*t(2,k+1,i)
        temp(2) = cs*t(2,k,i) + sn(1)*t(2,k+1,i) + sn(2)*t(1,k+1,i)
        t(1,k+1,i) = cs*t(1,k+1,i) - sn(1)*t(1,k,i) - sn(2)*t(2,k,i)
        t(2,k+1,i) = cs*t(2,k+1,i) - sn(1)*t(2,k,i) + sn(2)*t(1,k,i)
        t(:,k,i) = temp(:)
       enddo ! i
      endif
      do i = 1,k-1
       temp(1) = cs*t(1,i,k) + sn(1)*t(1,i,k+1) + sn(2)*t(2,i,k+1)
       temp(2) = cs*t(2,i,k) + sn(1)*t(2,i,k+1) - sn(2)*t(1,i,k+1)
       t(1,i,k+1) = cs*t(1,i,k+1) - sn(1)*t(1,i,k) + sn(2)*t(2,i,k)
       t(2,i,k+1) = cs*t(2,i,k+1) - sn(1)*t(2,i,k) - sn(2)*t(1,i,k)
       t(:,i,k) = temp(:)
      enddo ! i
      t(:,k,k) = t22(:)
      t(:,k+1,k+1) = t11(:)
! Accumulate transformation in the matrix Q.
      do i = 1,n
       temp(1) = cs*q(1,i,k) + sn(1)*q(1,i,k+1) + sn(2)*q(2,i,k+1)
       temp(2) = cs*q(2,i,k) + sn(1)*q(2,i,k+1) - sn(2)*q(1,i,k+1)
       q(1,i,k+1) = cs*q(1,i,k+1) - sn(1)*q(1,i,k) + sn(2)*q(2,i,k)
       q(2,i,k+1) = cs*q(2,i,k+1) - sn(1)*q(2,i,k) - sn(2)*q(1,i,k)
       q(:,i,k) = temp(:)
      enddo ! i
     enddo ! k
    endif

 end subroutine schurswap

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine balance(a,n,ilo,ihi,scales)
! Balance a complex matrix.
! Balancing may reduce the 1-norm of the matrix and improve the accuracy of
! the computed eigenvalues and/or eigenvectors.
!
! STEP ONE is to move off-diagonal zero entries to the lower left corner of
! the matrix, using permutations among row-and-column pairs.  The complete set
! of permutations is denoted by the similarity transformation P.
! STEP TWO is to apply a diagonal similarity transformation (D) to the "centre"
! rows and "centre" columns to make the rows and columns as close in norm as
! possible.
!
!    STEP 1          / T1  X   Y \   STEP 2   / T1     X*D          Y    \
! A  ----->  P*A*P = |  0  B   Z |   ----->   |  0  inv(D)*B*D  inv(D)*Z |
!                    \  0  0  T2 /            \  0      0          T2    /
!
! This subroutine corresponds to SUBROUTINE ZGEBAL of LAPACK.
! INPUT:
!   a(1,:,:),a(2,:,:) are input as the Re,Im parts of the n by n matrix.
!                     expected size: a(2,n,n)
!   n is the order of the matrix.  n>1 is assumed.
! OUTPUT:
!   a(1,:,:), a(2,:,:) are output as the Re,Im parts of the balanced matrix.
!                      expected size: a(2,n,n)
!   ilo and ihi are set to integers such that on exit
!               a(i,j) = 0 if i > j and j = 1,...,ilo-1 or I = ihi+1,...,n.
!   scales is details of the permutations and scaling factors applied to "a".
!          If P(j) is the index of the row and column interchanged with row and
!          column j and D(j) is the scaling factor applied to row and column j,
!          then scales(j) = P(j), for j = 1,...,ilo-1
!                         = D(j), for j = ilo,...,ihi
!                         = P(j), for j = ihi+1,...,n.
!          The order in which the interchanges are made is n to ihi+1,
!          then 1 to ilo-1. 
!          expected size: scales(n)

    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n
    integer(kind=KI), intent(out)                     :: ilo, ihi
    real(kind=KR2),   intent(out),   dimension(:)     :: scales

    logical                      :: conv
    integer(kind=KI)             :: i, j, ii, iflag, ica, ira
    real(kind=KR2), dimension(2) :: bit
    real(kind=KR2), parameter    :: sclfac=8.0_KR2, factor=0.95_KR2
    real(kind=KR2)               :: sfmin1, sfmin2, sfmax1, sfmax2, c, r, ca, &
                                    ra, catry, ratry, g, f, s

! Initialization.
    sfmin1 = dlamchs/dlamchp
    ilo = 1
    ihi = n

!*STEP ONE: Permute to put zeros in the lower left corner, if possible.

! First, re-order rows such that rows with more zeros are below rows with fewer.
    xrowswap: do
     j = ihi + 1
     jrow: do
      j = j - 1
      if (j<1) exit xrowswap
      iflag = 0
      i = 0
      icol: do
       i = i + 1
       if (i>ihi) exit icol
       if (i/=j .and. (a(1,j,i)/=0.0_KR2 .or. a(2,j,i)/=0.0_KR2)) then
        iflag = 1
        exit icol
       endif
      enddo icol
      if (iflag==0) then
       scales(ihi) = j
       if (j/=ihi) then
        do ii = 1,ihi
         bit(:) = a(:,ii,j)
         a(:,ii,j) = a(:,ii,ihi)
         a(:,ii,ihi) = bit(:)
        enddo ! ii
        do ii = ilo,n
         bit(:) = a(:,j,ii)
         a(:,j,ii) = a(:,ihi,ii)
         a(:,ihi,ii) = bit(:)
        enddo ! ii
       endif
       if (ihi==1) then
        exit xrowswap
       else
        exit jrow
       endif
      endif
     enddo jrow
     ihi = ihi - 1
    enddo xrowswap

! Now, re-order centre columns so that columns with more zeros are to the left.
    if (ihi>1) then
     colswap: do
      j = ilo - 1
      jcol: do
       j = j + 1
       if (j>ihi) exit colswap
       iflag = 0
       i = ilo - 1
       irow: do
        i = i + 1
        if (i>ihi) exit irow
        if (i/=j .and. (a(1,i,j)/=0.0_KR2 .or. a(2,i,j)/=0.0_KR2)) then
         iflag = 1
         exit irow
        endif
       enddo irow
       if (iflag==0) then
        scales(ilo) = j
        if (j/=ilo) then
         do ii = 1,ihi
          bit(:) = a(:,ii,j)
          a(:,ii,j) = a(:,ii,ilo)
          a(:,ii,ilo) = bit(:)
         enddo ! ii
         do ii = ilo,n
          bit(:) = a(:,j,ii)
          a(:,j,ii) = a(:,ilo,ii)
          a(:,ilo,ii) = bit(:)
         enddo ! ii
        endif
        exit jcol
       endif
      enddo jcol
      ilo = ilo + 1
     enddo colswap
     do i = ilo,ihi
      scales(i) = 1.0_KR2
     enddo ! i
    endif

!*STEP TWO: Balance the submatrix in rows ilo to ihi.
    if (ihi>1) then
     sfmax1 = 1.0_KR2/sfmin1
     sfmin2 = sfmin1*sclfac
     sfmax2 = 1.0_KR2/sfmin2
!-iterative loop for norm reduction.
     convergence: do
      conv = .true.
      do i = ilo,ihi
       c = 0.0_KR2
       r = 0.0_KR2
       do j = ilo,ihi
        if (j/=i) then
         c = c + abs(a(1,j,i)) + abs(a(2,j,i))
         r = r + abs(a(1,i,j)) + abs(a(2,i,j))
        endif
       enddo ! j
       ica = 1
       ca = abs(a(1,ica,i)) + abs(a(2,ica,i))
       if (ihi>1) then
        do j = 2,ihi
         catry = abs(a(1,j,i)) + abs(a(2,j,i))
         if (catry>ca) then
          ica = j
          ca = catry
         endif
        enddo ! j
       endif
       ira = ilo
       ra = abs(a(1,i,ira)) + abs(a(2,i,ira))
       if (ilo<n) then
        do j = ilo+1,n
         ratry = abs(a(1,i,j)) + abs(a(2,i,j))
         if (ratry>ra) then
          ira = j
          ra = ratry
         endif
        enddo ! j
       endif
!-guard against c=0 and/or r=0 due to underflow.
       if (c/=0.0_KR2 .and. r/=0.0_KR2) then
        g = r/sclfac
        f = 1.0_KR2
        s = c + r
        loloop: do
         if (c>=g .or. max(f,c,ca)>=sfmax2 .or. min(r,g,ra)<=sfmin2) exit loloop
         f = f*sclfac
         c = c*sclfac
         ca = ca*sclfac
         r = r/sclfac
         g = g/sclfac
         ra = ra/sclfac
        enddo loloop
        g = c/sclfac
        hiloop: do
         if (g<r .or. max(r,ra)>=sfmax2 .or. min(f,c,g,ca)<=sfmin2) exit hiloop
         f = f/sclfac
         c = c/sclfac
         g = g/sclfac
         ca = ca/sclfac
         r = r*sclfac
         ra = ra*sclfac
        enddo hiloop
!-now balance.
        iflag = 0
        if (c+r<factor*s) then
         if (f<1.0_KR2 .and. scales(i)<1.0_KR2) then
          if (f*scales(i)>sfmin1) then
           iflag = 1
          endif
         elseif (f>1.0_KR2 .and. scales(i)>1.0_KR2) then
          if (f*scales(i)<sfmax1/f) then
           iflag = 1
          endif
         endif
        endif
        if (iflag==1) then
         g = 1.0_KR2/f
         scales(i) = scales(i)*f
         conv = .false.
         do j = ilo,n
          a(:,i,j) = g*a(:,i,j)
         enddo ! j
         do j = 1,ihi
          a(:,j,i) = f*a(:,j,i)
         enddo ! j
        endif
       endif
      enddo ! i
      if (conv) exit convergence
     enddo convergence
    endif

 end subroutine balance

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine hessenberg(a,n,ilo,ihi,tau,workspace)
! Reduce a complex general matrix "a" to upper Hessenberg form h by
! unitary similarity transformation:  Q^dagger * a * Q = h .
! where Q = h_(ilo) * h_(ilo+1) * h_(ilo+2) * ... * h_(ihi-1)
!       The "elementary reflector" is h_i = 1 - tau_i * v_i * v_i^dagger
!       v_i is a complex vector with Re part v(1,:) and Im part v(2,:).
!           v(:,j) = 0 for j<=i and for j>ihi,
!           v(1,i+1) = 1 and v(2,i+1) = 0,
!           v(:,i+2), v(:,i+3), ..., v(:,ihi) are stored in
!           a(:,i+2,i), a(:,i+3,i), ..., a(:,ihi,i) for each vector v_i.
! This subroutine corresponds to SUBROUTINES ZGEHRD and ZGEHD2 of LAPACK.
! INPUT:
!   a(1,:,:),a(2,:,:) are input as the Re,Im parts of the n by n matrix.
!            expected size: a(2,n,n)
!   n is the order of the matrix.  n>1 is assumed.
!   ilo and ihi: It is assumed that "a" is already upper triangular in rows
!                and columns 1, 2, 3, ..., ilo-1 and ihi+1, ihi+2, ..., n, as
!                output by subroutine balance, otherwise choose ilo=1 and ihi=n.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are output as the Re,Im parts of the hessenberg matrix.
!            expected size: a(2,n,n)
!            Outputted entries in "a" below the first subdiagonal would have
!            been wasted, so some are used to store the "v_i" vectors in H_i,
!            as defined above.
!   tau(:,:) contains the scalar factors of the elementary reflectors in Q.
!   workspace(:,:) is used for scratch space.

    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n, ilo, ihi
    real(kind=KR2),   intent(out),   dimension(:,:)   :: tau, workspace

    integer(kind=KI)             :: i, row, col, dagger, myn
    real(kind=KR2), dimension(2) :: one, zero, tautemp, alpha

    tau = 0.0_KR2
    one(1) = 1.0_KR2
    one(2) = 0.0_KR2
    zero = 0.0_KR2
    do i = ilo,ihi-1
! Compute elementary reflector h_i to annihilate a(:,i+2,i)...a(:,ihi,i).
     alpha(:) = a(:,i+1,i)
     myn = ihi - i
     call makereflector(myn,alpha,a(:,min(i+2,n):n,i),tau(:,i))
     a(1,i+1,i) = 1.0_KR2
     a(2,i+1,i) = 0.0_KR2
! Apply h_i to a(:,ii,jj) for 1<=ii<=ihi and i+1<=jj<=ihi, from the right.
     if (tau(1,i)/=0.0_KR2 .or. tau(2,i)/=0.0_KR2) then
      dagger = 0
      row = ihi
      col = ihi - i
      call mvutil1(dagger,row,col,one,a(:,1:ihi,i+1:ihi),a(:,i+1:n,i),zero, &
                   workspace)
      tautemp(1) = -tau(1,i)
      tautemp(2) = -tau(2,i)
      call mvutil2(row,col,tautemp,workspace,a(:,i+1:n,i),a(:,1:ihi,i+1:ihi))
     endif
! Apply h_i^dagger to a(:,ii,jj) for i+1<=ii<=ihi and i+1<=jj<=n, from the left.
     if (tau(1,i)/=0.0_KR2 .or. tau(2,i)/=0.0_KR2) then
      dagger = 1
      row = ihi - i
      col = n - i
      call mvutil1(dagger,row,col,one,a(:,i+1:ihi,i+1:n),a(:,i+1:n,i),zero, &
                   workspace)
      tautemp(1) = -tau(1,i)
      tautemp(2) = tau(2,i)
      call mvutil2(row,col,tautemp,a(:,i+1:n,i),workspace,a(:,i+1:ihi,i+1:n))
     endif
     a(:,i+1,i) = alpha(:)
    enddo ! i

 end subroutine hessenberg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine qgenerator(n,ilo,ihi,a,tau,workspace)
! Generate a complex unitary matrix Q which is defined as the product of
! ihi-ilo elementary reflectors of order n, received from subroutine hessenberg.
! Q = h_(ilo) h_(ilo+1) . . . h_(ihi-1).
! This subroutine corresponds to SUBROUTINES ZUNGHR,ZUNGQR,ZUNG2R of LAPACK.
! INPUT:
!   n is the order of the matrix Q.
!   ilo,ihi are as defined in subroutine hessenberg.
!   a() contains the vectors which define the elementary reflectors.
!   tau(:,i) is the scalar factor, tau_i, for elementary reflector h_i.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the n by n unitary matrix Q.
!   workspace() is used for scratch space.

    integer(kind=KI), intent(in)                      :: n, ilo, ihi
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    real(kind=KR2),   intent(in),    dimension(:,:)   :: tau
    real(kind=KR2),   intent(out),   dimension(:,:)   :: workspace

    integer(kind=KI)             :: dagger, nh, i, j, row, col
    real(kind=KR2)               :: aRe
    real(kind=KR2), dimension(2) :: one, zero, tautemp

    dagger = 1
    nh = ihi - ilo
    one(1) = 1.0_KR2
    one(2) = 0.0_KR2
    zero(1) = 0.0_KR2
    zero(2) = 0.0_KR2
    if (n/=0) then
! Shift the vectors which define the elementary reflectors one column to the
! right, and set the first ilo and the last n-ihi rows and columns to those
! of the unit matrix.
     do j = ihi,ilo+1,-1
      do i = 1,j-1
       a(:,i,j) = 0.0_KR2
      enddo ! i
      do i = j+1,ihi
       a(:,i,j) = a(:,i,j-1)
      enddo ! i
      do i = ihi+1,n
       a(:,i,j) = 0.0_KR2
      enddo ! i
     enddo ! j
     do j = 1,ilo
      do i = 1,n
       a(:,i,j) = 0.0_KR2
      enddo ! i
      a(1,j,j) = 1.0_KR2
      a(2,j,j) = 0.0_KR2
     enddo ! j
     do j = ihi+1,n
      do i = 1,n
       a(:,i,j) = 0.0_KR2
      enddo ! i
      a(1,j,j) = 1.0_KR2
      a(2,j,j) = 0.0_KR2
     enddo ! j
! Generate Q(ilo+1:ihi,ilo+1:ihi).
     if (nh>0) then
! Apply h_i to a(:,ilo+1:ihi,ilo+1:ihi) from the left.
      do i = ihi,ilo+1,-1
       if (i<ihi) then
        a(1,i,i) = 1.0_KR2
        a(2,i,i) = 0.0_KR2
        if (tau(1,i-1)/=0.0_KR2 .or. tau(2,i-1)/=0.0_KR2) then
         row = ihi - i + 1
         col = ihi - i
         call mvutil1(dagger,row,col,one,a(:,i:ihi,i+1:ihi),a(:,i:ihi,i),zero, &
                      workspace)
         tautemp(1) = -tau(1,i-1)
         tautemp(2) = -tau(2,i-1)
         call mvutil2(row,col,tautemp,a(:,i:ihi,i),workspace,a(:,i:ihi,i+1:ihi))
        endif
        do j = i+1,ihi
         aRe = -tau(1,i-1)*a(1,j,i) + tau(2,i-1)*a(2,j,i)
         a(2,j,i) = -tau(1,i-1)*a(2,j,i) - tau(2,i-1)*a(1,j,i)
         a(1,j,i) = aRe
        enddo ! j
       endif
       a(1,i,i) = 1.0_KR2 - tau(1,i-1)
       a(2,i,i) = 0.0_KR2 - tau(2,i-1)
! Set a(:,ilo+1:i-1,i) to zero.
       do j = ilo+1,i-1
        a(:,j,i) = 0.0_KR2
       enddo ! j
      enddo ! i
     endif
    endif

 end subroutine qgenerator

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine evalues(ischur,n,ilo,ihi,h,w,z,workspace)
! Compute the eigenvalues of a complex upper Hessenberg matrix h, and,
! optionally, the matrices T and Z from the Schur decomposition h = Z T Z^H,
! where T is an upper triangular matrix (the Schur form), and Z is the unitary
! matrix of Schur vectors.
! Optionally Z may be postmultiplied into an input unitary matrix Q, so that
! this routine can give the Schur factorization of a matrix A which has been
! reduced to the Hessenberg form h by the unitary matrix Q:
! A = Q h Q^H = (QZ) T (QZ)^H.
! This subroutine corresponds to SUBROUTINE ZHSEQR of LAPACK.
! INPUT:
!   ischur=0 => compute eigenvalues only (no Schur vectors are computed).
!   ischur=1 => compute eigenvalues and the Schur form T.
!   n is the order of the matrix h.
!   ilo,ihi are as defined in subroutine hessenberg.
!   h is the upper Hessenberg matrix h.
!   z must contain an n by n matrix Q which is assumed to be equal to the
!     unit matrix except for the submatrix a(:,ilo:ihi,ilo,ihi).
!     (if ischur=0 then z is irrelevant.)
! OUTPUT:
!   h contains the upper triangular matrix T from the Schur decomposition
!     (the Schur form) if ischur=1.
!   w contains the computed eigenvalues.  If ischur=1, the eigenvalues are
!     stored in the same order as on the diagonal of the Schur form returned
!     in h, with w(:,i) = h(:,i,i).
!   z contains Q*Z.
!   workspace() is used for scratch space.

    integer(kind=KI), intent(in)                      :: ischur, n, ilo, ihi
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: h, z
    real(kind=KR2),   intent(out),   dimension(:,:)   :: w, workspace

    character(len=1)            :: side
    integer(kind=KI)            :: i, j, k, l, nh, ns, i1, i2, iflag, maxb, &
                                   itn, its, ii, ivec, myilo, nv, dagger, &
                                   row, col, nr, myid
    integer(kind=KI), parameter :: nsmax=15
    real(kind=KR2)              :: rtemp, hRe, unfl, ovfl, ulp, smlnum, &
                                      tst1, mysum, rtempii
    real(kind=KR2), dimension(2)                     :: one, temp, tau, tautemp
    real(kind=KR2), dimension(2,nmaxGMRES)           :: v, vv
    real(kind=KR2), dimension(2,nmaxGMRES,nmaxGMRES) :: s

!if(myid==0) then !TW 2/19/17
 !  print *, 'On the', itn, 'th itertion'!TW 2/17/17
!endif

! Definitions.
    iflag = 0
    one(1) = 1.0_KR2
    one(2) = 0.0_KR2

! Store the eigenvalues isolated by subroutine balance.
    do i = 1,ilo-1
     w(:,i) = h(:,i,i)
    enddo ! i
    do i = ihi+1,n
     w(:,i) = h(:,i,i)
    enddo ! i

    if (n/=0) then
     if (ilo==ihi) then
      w(:,ilo) = h(:,ilo,ilo)
     else
! Set rows and columns (ilo to ihi) to zero below the first subdiagonal.
      do j = ilo,ihi-2
       do i = j+2,n
        h(:,i,j) = 0.0_KR2
       enddo ! i
      enddo ! j
      nh = ihi - ilo + 1
! i1 and i2 are the indices of the first row and last column of h to which
! transformations must be applied. If eigenvalues only are being computed,
! i1 and i2 are re-set inside the main loop.
      if (ischur==1) then
       i1 = 1
       i2 = n
      else
       i1 = ilo
       i2 = ihi
      endif
! Ensure that the subdiagonal elements are real.
      do i = ilo+1,ihi
       temp(:) = h(:,i,i-1)
       if (temp(2)/=0.0_KR2) then
        rtemp = pythag2(temp(1),temp(2))
        h(1,i,i-1) = rtemp
        h(2,i,i-1) = 0.0_KR2
        temp = temp/rtemp
        if (i2>i) then
         do j = i+1,i2
          hRe = temp(1)*h(1,i,j) + temp(2)*h(2,i,j)
          h(2,i,j) = temp(1)*h(2,i,j) - temp(2)*h(1,i,j)
          h(1,i,j) = hRe
         enddo ! j
        endif
        do j = i1,i-1
         hRe = temp(1)*h(1,j,i) - temp(2)*h(2,j,i)
         h(2,j,i) = temp(1)*h(2,j,i) + temp(2)*h(1,j,i)
         h(1,j,i) = hRe
        enddo ! j
        if (i<ihi) then
         hRe = temp(1)*h(1,i+1,i) - temp(2)*h(2,i+1,i)
         h(2,i+1,i) = temp(1)*h(2,i+1,i) + temp(2)*h(1,i+1,i)
         h(1,i+1,i) = hRe
        endif
        if (ischur==1) then
         do j = ilo,ihi
          hRe = temp(1)*z(1,j,i) - temp(2)*z(2,j,i)
          z(2,j,i) = temp(1)*z(2,j,i) + temp(2)*z(1,j,i)
          z(1,j,i) = hRe
         enddo ! j
        endif
       endif
      enddo ! i
! Determine the order of the multi-shift QR algorithm to be used.
      
      ns = 6
      maxb = 250 ! Note: ZHSEQR.F used maxb=50 but I use a larger number so
!                        I can avoid ever using the multi-shift algorithm,
!                        since it relies on the choice "itn=30*nh" which seems
!                        not large enough for my requests to converge.
!                        If the multi-shift algorithm becomes desirable,
!                        I will increase itn and then convergence might be
!                        reached.


! ==============================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!- TESTING THIS VALUE. CHANGED FROM 900 to 1000 -PL 

     !maxb = 900
     maxb = 1000 

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! =============================================

     if (nh<maxb) then

!-use the standard double-shift algorithm.
       call QRdouble(ischur,n,ilo,ihi,h,w,ilo,ihi,z)
       !call QRdoubleLAPACK(ischur,n,ilo,ihi,h,w,ilo,ihi,z)
      else
!-use a multi-shift algorithm.
       unfl = dlamchs
       ovfl = 1.0_KR2/unfl
       ulp = dlamchp
       smlnum = unfl*(nh/ulp)
       itn = 20*nh
!-The main multi-shift loop begins here. i is the loop index and decreases
! from ihi to ilo in steps of at most maxb.  Each iteration of the loop works
! with the active submatrix in rows and columns l to i.  Eigenvalues i+1 to
! ihi have already converged.  Either l = ilo, or h(l,l-1) is negligible so
! that the matrix splits.
       i = ihi
       mainloop: do
        if (i<ilo) exit mainloop
!-perform multiple-shift QR iterations on rows and columns ilo to i until a
! submatrix of order at most maxb splits off at the bottom because a
! subdiagonal element has become negligible.
        l = ilo
        its = -1
        itsloop: do
         its = its + 1
         if (its>itn) exit itsloop
!-look for a single small subdiagonal element.
         k = i + 1
         kloop: do
          k = k - 1
          if (k<=l) exit kloop
          tst1 = abs(h(1,k-1,k-1)) + abs(h(2,k-1,k-1)) &
               + abs(h(1,k,k)) + abs(h(2,k,k))
          if (tst1==0.0_KR2 .and. l/=i+1) then
           do j = l,i
            mysum = 0.0_KR2
            do ii = l,min(i,j+1)
             mysum = mysum + sqrt(h(1,ii,j)**2+h(2,ii,j)**2)
            enddo ! ii
            tst1 = max(tst1,mysum)
           enddo ! j
          endif
          if (abs(h(1,k,k-1))<=max(ulp*tst1,smlnum)) exit kloop
         enddo kloop
         l = k
         if (l>ilo) then
          h(:,l,l-1) = 0.0_KR2
         endif
! Exit from loop if a submatrix of order <= maxb has split off.
         if (l>=i-maxb+1) then
          iflag = 1
          exit itsloop
         endif
! Now the active submatrix is in rows and columns l to i. If eigenvalues only
! are being computed, only the active submatrix need be transformed.
         if (ischur==0) then
          i1 = l
          i2 = i
         endif
         if (its==20 .or. its==30) then
! Exceptional shifts.
          do ii = i-ns+1,i
           w(1,ii) = 1.5_KR2*(abs(h(1,ii,ii-1))+abs(h(1,ii,ii)))
           w(2,ii) = 0.0_KR2
          enddo ! ii
         else
! Use eigenvalues of trailing submatrix of order ns as shifts.
          do j = 1,ns
           do ii = 1,ns
            s(:,ii,j) = h(:,i-ns+ii,i-ns+j)
           enddo ! ii
          enddo ! j
          ivec = 0
          myilo = 1
          call QRdouble(ivec,ns,myilo,ns,s,w(:,i-ns+1:n),myilo,ns,z)
          !call QRdoubleLAPACK(ivec,ns,myilo,ns,s,w(:,i-ns+1:n),myilo,ns,z)
         endif
! Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns)), where G is the
! Hessenberg submatrix h(:,l:i,l:i) and w is the vector of shifts (stored in w).
! The result is stored in the local array v.
         v(:,1) = one
         do ii = 2,ns+1
          v(:,ii) = 0.0_KR2
         enddo ! ii
         nv = 1
         do j = i-ns+1,i
          do ii = 1,nv+1
           vv(:,ii) = v(:,ii)
          enddo ! ii
          dagger = 0
          row = nv + 1
          col = nv
          call mvutil1(dagger,row,col,one,h(:,l:l+row-1,l:l+col-1),vv,-w(:,j),v)
          nv = nv + 1
! Scale v(:,1:nv) so that max(abs(v(:,i))) = 1.
! If v is zero, reset it to the unit vector.
          rtemp = abs(v(1,1)) + abs(v(2,1))
          do ii = 2,nv
           rtempii = abs(v(1,ii)) + abs(v(2,ii))
           if (rtempii>rtemp) then
            rtemp = rtempii
           endif
          enddo ! ii
          if (rtemp==0.0_KR2) then
           v(:,1) = one(:)
           do ii = 2,nv
            v(:,ii) = 0.0_KR2
           enddo ! ii
          else
           rtemp = max(rtemp,smlnum)
           do ii = 1,nv
            v(:,ii) = v(:,ii)/rtemp
           enddo ! ii
          endif
         enddo ! j
! Multiple-shift QR step.
         do k = l,i-1
! The first iteration of this loop determines a reflection G from the vector
! v and applies it from left and right to h, thus creating a nonzero bulge
! below the subdiagonal.
! Each subsequent iteration determines a reflection G to restore the
! Hessenberg form in the (k-1)th column, and thus chases the bulge one step
! toward the bottom of the active submatrix. nr is the order of G.
          nr = min(ns+1,i-k+1)
          if (k>l) then
           do ii = 1,nr
            v(:,ii) = h(:,k-2+ii,k-1)
           enddo ! ii
          endif
          call makereflector(nr,v(:,1),v(:,2:nr),tau)
          if (k>l) then
           h(:,k,k-1) = v(:,1)
           do ii = k+1,i
            h(:,ii,k-1) = 0.0_KR2
           enddo ! ii
          endif
          v(:,1) = one(:)
! Apply G^dagger from the left to transform the rows of the matrix in columns k
! to i2.
          side = "L"
          row = nr
          col = i2 - k + 1
          tautemp(1) = tau(1)
          tautemp(2) = -tau(2)
          call multreflector(side,row,col,v,tautemp,h(:,k:nr+k-1,k:i2), &
                             workspace)
! Apply G from the right to transform the columns of the matrix in rows i1 to
! min(k+nr,i).
          side = "R"
          row = min(k+nr,i) - i1 + 1
          col = nr
          call multreflector(side,row,col,v,tau,h(:,i1:row+i1-1,k:nr+k-1), &
                             workspace)
          if (ischur==1) then
! Accumulate transformations in the matrix z.
          side = "R"
          call multreflector(side,nh,nr,v,tau,z(:,ilo:nh+ilo-1,k:nr+k-1), &
                             workspace)
          endif
         enddo ! k
! Ensure that h(:,i,i-1) is real.
         temp(:) = h(:,i,i-1)
         if (temp(2)/=0.0_KR2) then
          rtemp = pythag2(temp(1),temp(2))
          h(1,i,i-1) = rtemp
          h(2,i,i-1) = 0.0_KR2
          temp = temp/rtemp
          if (i2>i) then
           do ii = i+1,i2
            hRe = temp(1)*h(1,i,ii) + temp(2)*h(2,i,ii)
            h(2,i,ii) = temp(1)*h(2,i,ii) - temp(2)*h(1,i,ii)
            h(1,i,ii) = hRe
           enddo ! ii
          endif
          do ii = i1,i-1
           hRe = temp(1)*h(1,ii,i) - temp(2)*h(2,ii,i)
           h(2,ii,i) = temp(1)*h(2,ii,i) + temp(2)*h(1,ii,i)
           h(1,ii,i) = hRe
          enddo ! ii
          if (ischur==1) then
           do ii = ilo,ilo+nh-1
            hRe = temp(1)*z(1,ii,i) - temp(2)*z(2,ii,i)
            z(2,ii,i) = temp(1)*z(2,ii,i) + temp(2)*z(1,ii,i)
            z(1,ii,i) = hRe
           enddo ! ii
          endif
         endif
        enddo itsloop
! Failure to converge in remaining number of iterations.
        if (iflag==1) then
         iflag = 0
        else
         open(unit=8,file="PSEUDOLAPACK.ERROR",action="write",status="replace" &
             ,form="formatted")
          write(unit=8,fmt=*) "Failed to converge in evalues. Is itn too small?"
         close(unit=8,status="keep")
         stop
        endif
! A submatrix of order <= maxb in rows and columns l to i has split off.
! Use the double-shift QR algorithm to handle it.
        call QRdouble(ischur,n,l,i,h,w,ilo,ihi,z)
        !call QRdoubleLAPACK(ischur,n,l,i,h,w,ilo,ihi,z)
! Decrement number of remaining iterations, and return to start of the main
! loop with a new value of i.
        itn = itn - its
       ! if (myid==0) then
        !   print *, 'On the', itn, 'th itertion'!TW 2/17/17
       ! endif
        i = l - 1
       enddo mainloop
      endif
     endif
    endif

 end subroutine evalues

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine Tvectors(n,t,vr,work,rwork)
! Compute all of the right eigenvectors of a complex upper triangular matrix t
! and return the product Q*x, where x is the matrix of right eigenvectors and
! Q is an inputted unitary matrix.  If t was obtained from the Schur
! factorization of an original matrix A = Q*t*Q^dagger, then Q*x is the matrix
! of right eigenvectors of A.
! This subroutine corresponds to SUBROUTINE ZTREVC of LAPACK.
! INPUT:
!   n is the order of the matrix t.
!   t is an upper triangular matrix.  t is modified, but restored on exit.
!   vr contains an n-by-n matrix Q (usually the unitary matrix Q of Schur
!      vectors returned by subroutine evalues).
!   work(1,:),work(2,:) are the Re,Im parts of a scratch vector.
!   rwork(:) is another scratch vector.
! OUTPUT:
!   vr contains the matrix Q*x.

    integer(kind=KI), intent(in)                      :: n
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: t, vr
    real(kind=KR2),   intent(out),   dimension(:,:)   :: work
    real(kind=KR2),   intent(out),   dimension(:)     :: rwork

    integer(kind=KI)             :: i, j, k, is, kk, ii, myn, dagger, row, col
    real(kind=KR2)               :: unfl, ovfl, ulp, smlnum, smin, bit, remax, &
                                    smax
    real(kind=KR2), dimension(2) :: alpha, myscale

    if (n/=0) then
! Set the constants to control overflow and underflow.
     unfl = dlamchs
     ovfl = 1.0_KR2/unfl
     ulp = dlamchp
     smlnum = unfl*(n/ulp)
! Store the diagonal elements of t in working array work.
     do i = 1,n
      work(:,i+n) = t(:,i,i)
     enddo ! i
! Compute 1-norm of each column of strictly upper triangular part of t to
! control overflow in triangular solver.
     rwork = 0.0_KR2
     do j = 2,n
      do i = 1,j-1
       rwork(j) = rwork(j) + abs(t(1,i,j)) + abs(t(2,i,j))
      enddo ! i
     enddo ! j
! Compute right eigenvectors.
     is = n
     do kk = n,1,-1
      smin = max(ulp*(abs(t(1,kk,kk))+abs(t(2,kk,kk))),smlnum)
      work(1,1) = 1.0_KR2
      work(2,1) = 0.0_KR2
!-form right-hand side.
      do k = 1,kk-1
       work(:,k) = -t(:,k,kk)
      enddo ! k
!-solve the triangular system: (t(1:kk-1,1:kk-1) - t(kk,kk))*x = myscale*work.
      do k = 1,kk-1
       t(:,k,k) = t(:,k,k) - t(:,kk,kk)
       if (abs(t(1,k,k))+abs(t(2,k,k))<smin) then
        t(1,k,k) = smin
        t(2,k,k) = 0.0_KR2
       endif
      enddo ! k
      if (kk>1) then
       myn = kk - 1
       call mvutil3(myn,t,work,myscale(1),rwork)
       myscale(2) = 0.0_KR2
       work(:,kk) = myscale(:)
      endif
!-copy the vector Q*x to vr and normalize.
      if (kk>1) then
       dagger = 0
       row = n
       col = kk - 1
       alpha(1) = 1.0_KR2
       alpha(2) = 0.0_KR2
       call mvutil1(dagger,row,col,alpha,vr,work,myscale,vr(:,:,kk))
      endif
      ii = 1
      if (n>1) then 
       smax = abs(vr(1,1,kk)) + abs(vr(2,1,kk))
       do i = 2,n
        bit = abs(vr(1,i,kk)) + abs(vr(2,i,kk))
        if (bit>smax) then
         ii = i
         smax = bit
        endif
       enddo ! i
      endif
      remax = 1.0_KR2/(abs(vr(1,ii,kk))+abs(vr(2,ii,kk)))
      vr(:,:,kk) = remax*vr(:,:,kk)
!-set back the original diagonal elements of t.
      do k = 1,kk-1
       t(:,k,k) = work(:,k+n)
      enddo ! k
      is = is - 1
     enddo ! kk
    endif

 end subroutine Tvectors

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine evectors(n,ilo,ihi,myscale,v)
! Form the right eigenvectors of a general complex matrix by backward
! transformation on the computed eigenvectors of the balanced matrix outputted
! by subroutine balance.
! This subroutine corresponds to SUBROUTINE ZGEBAK of LAPACK.
! INPUT:
!   n is the order of the square matrix v.
!   ilo,ihi are the integers determined by subroutine balance.
!   myscale contains details of the permutation and scaling factors, as 
!           returned by subroutine balance.
!   v is the matrix of right eigenvectors to be transformed, as returned by
!     subroutine Tvectors.
! OUTPUT:
!   v is overwritten by the transformed eigenvectors.

    integer(kind=KI), intent(in)                      :: n, ilo, ihi
    real(kind=KR2),   intent(in),    dimension(:)     :: myscale
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: v

    integer(kind=KI)             :: i, j, k, ii
    real(kind=KR2), dimension(2) :: vtemp

    if (n/=0) then
! Backward balance.
   !  print *,'ilo=',ilo,'ili=',ihi !BS 1/19/2016
     if (ilo/=ihi) then
      do i = ilo,ihi
       v(:,i,:) = myscale(i)*v(:,i,:)
      enddo ! i
     endif
! Backward permutation.
     do ii = 1,n
      i = ii
      if (i<ilo .or. i>ihi) then
       if (i<ilo) i=ilo-ii
       k = myscale(i)
       if (k/=i) then
        do j = 1,n
         vtemp(:) = v(:,i,j)
         v(:,i,j) = v(:,k,j)
         v(:,k,j) = vtemp(:)
        enddo ! j
       endif
      endif
     enddo ! ii
    endif

 end subroutine evectors

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine bidiagonalize(m,n,a,d,e,tauq,taup)
! Reduce a general complex m by n matrix "a" to upper (if m>=n) or lower (if
! m<n) bidiagonal form "b" by a unitary transformation: q^H * a * p = b.
! The matrices q and p are represented as products of elementary reflectors:
! If m >= n, then q = H(1) H(2) . . . H(n)  and  p = G(1) G(2) . . . G(n-1)
!  Each H(i) and G(i) has the form: H(i) = I - tauq*v*v^dagger,
!                                   G(i) = I - taup*u*u^dagger,
!  where tauq and taup are complex scalars, and v and u are complex vectors.
! If m < n, then q = H(1) H(2) . . . H(m-1)  and  p = G(1) G(2) . . . G(m)
!  Each H(i) and G(i) has the form: H(i) = I - tauq*v*v^dagger,
!                                   G(i) = I - taup*u*u^dagger,
!  where tauq and taup are complex scalars, and v and u are complex vectors.
! This subroutine corresponds to SUBROUTINES ZGEBRD,ZGEBD2 of LAPACK.
! INPUT:
!   m is the number of rows in the matrix "a".
!   n is the number of columns in the matrix "a".
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the m by n matrix to be reduced.
! OUTPUT:
!   a() If m >= n, the diagonal and the first superdiagonal are overwritten
!        with the upper bidiagonal matrix b; the elements below the diagonal,
!        along with the array tauq, represent the unitary matrix q as a product
!        of elementary reflectors, and the elements above the first
!        superdiagonal, along with the array taup, represent the unitary matrix
!        p as a product of elementary reflectors.
!       If m < n, the diagonal and the first subdiagonal are overwritten with
!        the lower bidiagonal matrix b; the elements below the first
!        subdiagonal, along with the array tauq, represent the unitary matrix q
!        as a product of elementary reflectors, and the elements above the
!        diagonal, along with the array taup, represent the unitary matrix p as
!        a product of elementary reflectors.
!   d contains the diagonal elements of the bidiagonal matrix b: d(i) = a(i,i).
!   e contains the off-diagonal elements of the bidiagonal matrix b:
!        If m >= n, e(i) = a(i,i+1) for i = 1,2,...,n-1.
!        If m < n, e(i) = a(i+1,i) for i = 1,2,...,m-1.
!   tauq contains the scalar factors of the elementary reflectors which
!        represent the unitary matrix Q.
!        expected size: tauq(2,min(m,n))
!   taup contains the scalar factors of the elementary reflectors which
!        represent the unitary matrix Q.
!        expected size: taup(2,min(m,n))

    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(out),   dimension(:)     :: d, e
    real(kind=KR2),   intent(out),   dimension(:,:)   :: tauq, taup
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a

    integer(kind=KI)                       :: i, j, row, col, dagger
    real(kind=KR2), dimension(2)           :: alpha, one, zero
    real(kind=KR2), dimension(2,nmaxGMRES) :: workspace

! Useful definitions.
    one(1) = 1.0_KR2
    one(2) = 0.0_KR2
    zero = 0.0_KR2

! Reduce to upper bidiagonal form.
    if (m>=n) then
!-generate elementary reflector h(:,i) to annihilate a(:,i+1:m,i).
     do i = 1,n
      row = m - i + 1
      alpha(:) = a(:,i,i)
      call makereflector(row,alpha,a(:,min(i+1,m):m,i),tauq(:,i))
      d(i) = alpha(1)
      a(:,i,i) = one(:)
!-apply h(:,i)^dagger to a(:,i:m,i+1:n) from the left.
      if (tauq(1,i)/=0.0_KR2 .or. tauq(2,i)/=0.0_KR2) then
       dagger = 1
       row = m - i + 1
       col = n - i
       call mvutil1(dagger,row,col,one,a(:,i:m,i+1:n),a(:,i:m,i),zero,workspace)
       alpha(1) = -tauq(1,i)
       alpha(2) = tauq(2,i)
       call mvutil2(row,col,alpha,a(:,i:m,i),workspace,a(:,i:m,i+1:n))
      endif
      a(1,i,i) = d(i)
      a(2,i,i) = 0.0_KR2
!-generate elementary reflector g(i) to annihilate a(:,i,i+2:n).
      if (i<n) then
       do j = i+1,n
        a(2,i,j) = -a(2,i,j)
       enddo ! j
       row = n - i
       alpha(:) = a(:,i,i+1)
       call makereflector(row,alpha,a(:,i,min(i+2,n):n),taup(:,i))
       e(i) = alpha(1)
       a(:,i,i+1) = one(:)
!-apply g(i) to a(:,i+1:m,i+1:n) from the right.
       if (taup(1,i)/=0.0_KR2 .and. taup(2,i)/=0.0_KR2) then
        dagger = 0
        row = m - i
        col = n - i
        call mvutil1(dagger,row,col,one,a(:,i+1:m,i+1:n),a(:,i,i+1:n),zero, &
                     workspace)
        alpha(1) = -taup(1,i)
        alpha(2) = -taup(2,i)
        call mvutil2(row,col,alpha,workspace,a(:,i,i+1:n),a(:,i+1:m,i+1:n))
       endif
       do j = i+1,n
        a(2,i,j) = -a(2,i,j)
       enddo ! j
       a(1,i,i+1) = e(i)
       a(2,i,i+1) = 0.0_KR2
      else
       taup(:,i) = 0.0_KR2
      endif
     enddo ! i

! Reduce to lower bidiagonal form.
    else
!-generate elementary reflector g(i) to annihilate a(:,i,i+1:n).
     do i = 1,m
      do j = i,n
       a(2,i,j) = -a(2,i,j)
      enddo ! j
      row = n - i + 1
      alpha(:) = a(:,i,i)
      call makereflector(row,alpha,a(:,i,min(i+1,n):n),taup(:,i))
      d(i) = alpha(1)
      a(:,i,i) = one(:)
!-apply g(i) to a(:,i+1:m,i:n) from the right.
      if (taup(1,i)/=0.0_KR2 .or. taup(2,i)/=0.0_KR2) then
       dagger = 0
       row = m - i
       col = n - i + 1
       call mvutil1(dagger,row,col,one,a(:,min(i+1,m):m,i:n),a(:,i,i:n),zero, &
                    workspace)
       alpha(:) = -taup(:,i)
       call mvutil2(row,col,alpha,workspace,a(:,i,i:n),a(:,min(i+1,m):m,i:n))
      endif
      do j = i,n
       a(2,i,j) = -a(2,i,j)
      enddo ! j
      a(1,i,i) = d(i)
      a(2,i,i) = 0.0_KR2
!-generate elementary reflector h(i) to annihilate a(:,i+2:m,i).
      if (i<m) then
       row = m - i
       alpha(:) = a(:,i+1,i)
       call makereflector(row,alpha,a(:,min(i+2,m):m,i),tauq(:,i))
       e(i) = alpha(1)
       a(:,i+1,i) = one(:)
!-apply H(i)^dagger to a(:,i+1:m,i+1:n) from the left.
       if (tauq(1,i)/=0.0_KR2 .or. tauq(2,i)/=0.0_KR2) then
        dagger = 1
        row = m - i
        col = n - i
        call mvutil1(dagger,row,col,one,a(:,i+1:m,i+1:n),a(:,i+1:m,i),zero, &
                     workspace)
        alpha(1) = -tauq(1,i)
        alpha(2) = tauq(2,i)
        call mvutil2(row,col,alpha,a(:,i+1:m,i),workspace,a(:,i+1:m,i+1:n))
       endif
       a(1,i+1,i) = e(i)
       a(2,i+1,i) = 0.0_KR2
      else
       tauq(:,i) = 0.0_KR2
      endif
     enddo ! i
    endif

 end subroutine bidiagonalize

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine bidiagqgen(n,k,a,tau,workspace)
! Generate the complex unitary matrix Q determined by subroutine bidiagonalize
! when reducing a complex matrix "a" to bidiagonal form: a = Q * b * p^dagger.
! Q is defined as a product of elementary reflectors h_i.
! "a" is assumed to have been an n by k matrix, and Q is of order n:
!  if n >= k, Q = h_1 h_2... h_k and subroutine bidiagqgen returns the first n
!             columns of Q.
!  if m < k, Q = h_1 h_2... h_(m-1) and subroutine bidiagqgen returns Q as an
!  n by n matrix.
! This subroutine corresponds to SUBROUTINES ZUNGBR of LAPACK.
! INPUT:
!   n is the number of rows = number of columns of Q to be computed.
!   k is the number of columns in the original m by k matrix.
!   a() contains the elementary reflectors: h_i is in its i'th column.
!   tau(:,i) is the scalar factor, tau_i, for elementary reflector h_i.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the n by n matrix Q.
!   workspace() is used for scratch space.

    integer(kind=KI), intent(in)                      :: n, k
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    real(kind=KR2),   intent(in),    dimension(:,:)   :: tau
    real(kind=KR2),   intent(out),   dimension(:,:)   :: workspace

    integer(kind=KI) :: i, j, row, col, refl

    if (n>=k) then
     row = n
     col = n
     refl = k
     call qgenkernel(row,col,refl,a,tau,workspace)
    else
! Shift the vectors which define the elementary reflectors one column to the
! right, and set the first row and column of Q to those of the unit matrix.
     do j = n,2,-1
      a(:,1,j) = 0.0_KR2
      do i = j+1,n
       a(:,i,j) = a(:,i,j-1)
      enddo ! i
     enddo ! j
     a(1,1,1) = 1.0_KR2
     a(2,1,1) = 0.0_KR2
     do i = 2,n
      a(:,i,1) = 0.0_KR2
     enddo ! i
     if (n>1) then
! Form Q(2:m,2:m).
      row = n - 1
      col = n - 1
      refl = n - 1
      call qgenkernel(row,col,refl,a(:,2:n,2:n),tau,workspace)
     endif
    endif

 end subroutine bidiagqgen

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine qgenkernel(m,n,k,a,tau,workspace)
! Generate an m by n complex matrix Q with orthonormal columns, which is
! defined as the first n columns of a product of k elementary reflectors of
! order m: Q = h_1 h_2... h_k.
! This subroutine corresponds to SUBROUTINES ZUNGQR,ZUNG2R of LAPACK.
! INPUT:
!   m is the number of rows in the matrix Q.
!   n is the number of columns in the matrix Q.
!   k is the number of elementary reflectors whose product defines the matrix Q.
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the vectors which define the
!                     elementary reflectors: the i'th column is h_i.
!   tau(:,i) is the scalar factor, tau_i, for elementary reflector h_i.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the m by n matrix Q.
!   workspace() is used for scratch space.

    integer(kind=KI), intent(in)                      :: m, n, k
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    real(kind=KR2),   intent(in),    dimension(:,:)   :: tau
    real(kind=KR2),   intent(out),   dimension(:,:)   :: workspace

    integer(kind=KI)             :: dagger, i, j, row, col
    real(kind=KR2)               :: aRe
    real(kind=KR2), dimension(2) :: one, zero, tautemp

    dagger = 1
    one(1) = 1.0_KR2
    one(2) = 0.0_KR2
    zero = 0.0_KR2
    if (n>0) then
! Initialize columns k+1:n to columns of the unit matrix.
     do j = k+1,n
      do i = 1,m
       a(:,i,j) = 0.0_KR2
      enddo ! i
       a(:,j,j) = one(:)
     enddo ! j
! Apply h_i to a(:,i:m,i:n) from the left.
      do i = k,1,-1
       if (i<n) then
        a(:,i,i) = one(:)
        if (tau(1,i)/=0.0_KR2 .or. tau(2,i)/=0.0_KR2) then
         row = m - i + 1
         col = n - i
         call mvutil1(dagger,row,col,one,a(:,i:m,i+1:n),a(:,i:m,i),zero, &
                      workspace)
         tautemp(1) = -tau(1,i)
         tautemp(2) = -tau(2,i)
         call mvutil2(row,col,tautemp,a(:,i:m,i),workspace,a(:,i:m,i+1:n))
        endif
       endif
       if (i<m) then
        do j = i+1,m
         aRe = -tau(1,i)*a(1,j,i) + tau(2,i)*a(2,j,i)
         a(2,j,i) = -tau(1,i)*a(2,j,i) - tau(2,i)*a(1,j,i)
         a(1,j,i) = aRe
        enddo ! j
       endif
       a(:,i,i) = one(:) - tau(:,i)
! Set a(:,1:i-1,i) to zero.
       do j = 1,i-1
        a(:,j,i) = 0.0_KR2
       enddo ! j
      enddo ! i
    endif

 end subroutine qgenkernel

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine svdkernel(uplo,n,nru,d,e,u,rwork)
! Compute the singular value decomposition (SVD) of a real n by n (upper or
! lower) bidiagonal matrix b:  b = Q * s * P^T, where s is a diagonal matrix
! with non-negative diagonal elements (the singular values of b), and Q and P
! are orthogonal matrices.  The routine computes s and u*Q for a given complex
! input matrix u.
! This subroutine corresponds to SUBROUTINE ZBDSQR of LAPACK.
! INPUT:
!   uplo="U" if b is upper diagonal; uplo="L" if b is lower diagonal.
!   n is the order of the matrix b.
!   nru is the number of rows being requested from the matrix u.
!       NOTE: This code assumes nru>0.
!   d contains the diagonal elements of the bidiagonal matrix b.
!   e contains the off-diagonal elements of the bidiagonal matrix b.
!   u(1,:,:),u(2,:,:) are the Re,Im parts of an nru by n matrix.
! OUTPUT:
!   d is overwritten by the singular values of b in decreasing order.
!   e is destroyed.
!   u is overwritten by u*Q.
!   rwork is a scratch vector.

    character(len=*), intent(in)                      :: uplo
    integer(kind=KI), intent(in)                      :: n, nru
    real(kind=KR2),   intent(inout), dimension(:)     :: d, e
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: u
    real(kind=KR2),   intent(out),   dimension(:)     :: rwork

    character(len=1)             :: side, fwdbwd
    integer(kind=KI)             :: i, j, m, nm1, nm12, nm13, idir, maxit, &
                                    iter, oldll, oldm, isub, iflag, mainflag, &
                                    ll, lll, col
    real(kind=KR2)               :: eps, unfl, cs, sn, r, tol, tolmul, smax, &
                                    sminl, smin, mu, sminoa, thresh, abss, &
                                    abse, sigmx, sigmn, cosl, sinl, cosr, &
                                    sinr, sminlo, shift, sll, oldcs, oldsn, &
                                    f, g, h
    real(kind=KR2), dimension(2) :: utemp
    integer(kind=KI), parameter  :: maxitr=6
! maxitr controls the maximum number of passes of the algorithm through its
! inner loop. The algorithms stops (and so fails to converge) if the number
! of passes through the inner loop exceeds maxitr*n^2.

    if (n>1) then
     nm1 = n - 1
     nm12 = 2*nm1
     nm13 = nm12 + nm1
     idir = 0
     eps = dlamche
     unfl = dlamchs
! If matrix lower bidiagonal, rotate to be upper bidiagonal by applying Givens
! rotations on the left.
     if (uplo=="L") then
      do i = 1,n-1
       call rotate2(d(i),e(i),cs,sn,r)
       d(i) = r
       e(i) = sn*d(i+1)
       d(i+1) = cs*d(i+1)
       rwork(i) = cs
       rwork(nm1+i) = sn
      enddo ! i
! Update singular vectors.
      side = "R"
      fwdbwd = "F"
      call mmutil1(side,fwdbwd,nru,n,rwork(1:n-1),rwork(n:2*n-2),u)
     endif
! Compute singular values to relative accuracy tol.  (By setting tol to be
! negative, algorithm will compute singular values to absolute accuracy
! abs(tol)*norm(input matrix)).
     tolmul = max(10.0_KR2,min(100.0_KR2,1.0_KR2/sqrt(sqrt(sqrt(eps)))))
     tol = tolmul*eps
! Compute approximate maximum and minimum singular values.
     smax = 0.0_KR2
     do i = 1,n
      smax = max(smax,abs(d(i)))
     enddo ! i
     do i = 1,n-1
      smax = max(smax,abs(e(i)))
     enddo ! i
     sminl = 0.0_KR2
     if (tol>=0.0_KR2) then
! Relative accuracy desired.
      sminoa = abs(d(1))
      if (sminoa/=0.0_KR2) then
       mu = sminoa
       i = 1
       iloop: do
        i = i + 1
        if (i>n) exit iloop
        mu = abs(d(i))*(mu/(mu+abs(e(i-1))))
        sminoa = min(sminoa,mu)
        if (sminoa==0.0_KR2) exit iloop
       enddo iloop
      endif
      sminoa = sminoa/sqrt(real(n,KR2))
      thresh = max(tol*sminoa,maxitr*n*n*unfl)
     else
! Absolute accuracy desired.
      thresh = max(abs(tol)*smax,maxitr*n*n*unfl)
     endif
! Prepare for main iteration loop for the singular values (maxit is the maximum
! number of passes through the inner loop permitted before nonconvergence is
! signalled.)
     maxit = maxitr*n*n
     iter = 0
     oldll = -1
     oldm = -1
! m points to the last element of the unconverged part of the matrix.
     m = n

! Begin main iteration loop.
     mainloop: do
      mainflag = 0
! Check for convergence.
      if (m<=1) exit mainloop
! Check to see if the maximum iteration count has been exceeded.
      if (iter>maxit) then
       open(unit=8,file="PSEUDOLAPACK.ERROR",action="write",status="replace" &
           ,form="formatted")
        write(unit=8,fmt=*) "Failed to converge in svdkernel."
       close(unit=8,status="keep")
       stop
      endif
! Find diagonal block of matrix to work on.
      if (tol<0.0_KR2 .and. abs(d(m))<=thresh) d(m)=0.0_KR2
      smax = abs(d(m))
      smin = smax
      iflag = 0
      lll = 0
      l3loop1: do
       lll = lll + 1
       if (lll>m-1) exit l3loop1
       ll = m - lll
       abss = abs(d(ll))
       abse = abs(e(ll))
       if (tol<0.0_KR2 .and. abss<=thresh) d(ll)=0.0_KR2
       if (abse<=thresh) then
        iflag = 1
        exit l3loop1
       endif
       smin = min(smin,abss)
       smax = max(smax,abss,abse)
      enddo l3loop1
      if (iflag==0) then
       ll = 0
      else
       e(ll) = 0.0_KR2
! Matrix splits since e(ll) = 0.
       if (ll==m-1) then
! Convergence of bottom singular value, return to top of loop (via mainflag=1).
        m = m - 1
        mainflag = 1
       endif
      endif

      if (mainflag==0) then
       ll = ll + 1
! e(ll) through e(m-1) are nonzero, e(ll-1) is zero.
       if (ll==m-1) then
! 2 by 2 block, handle separately.
        call svd2by2(d(m-1),e(m-1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl)
        d(m-1) = sigmx
        e(m-1) = 0.0_KR2
        d(m) = sigmn
! Compute singular vectors.
        do i = 1,nru
         utemp(:) = cosl*u(:,i,m-1) + sinl*u(:,i,m)
         u(:,i,m) = cosl*u(:,i,m) - sinl*u(:,i,m-1)
         u(:,i,m-1) = utemp(:)
        enddo ! i
        m = m - 2
        mainflag = 1
       endif
      endif

! If working on new submatrix, choose shift direction (from larger end diagonal
! element towards smaller).
      if (mainflag==0) then
       if (ll>oldm .or. m<oldll) then
        if (abs(d(ll))>=abs(d(m))) then
! Chase bulge from top (big end) to bottom (small end).
         idir = 1
        else
! Chase bulge from bottom (big end) to top (small end).
         idir = 2
        endif
       endif
      endif

! Apply convergence tests.
      if (mainflag==0) then
       if (idir==1) then
! Run convergence test in forward direction.
! First apply standard test to bottom of matrix.
        if (abs(e(m-1))<=abs(tol)*abs(d(m)) .or. &
           (tol<0.0_KR2 .and. abs(e(m-1))<=thresh)) then
         e(m-1) = 0.0_KR2
         mainflag = 1
        endif
        if (mainflag==0 .and. tol>=0.0_KR2) then
! If relative accuracy desired, apply convergence criterion forward.
         mu = abs(d(ll))
         sminl = mu
         lll = ll - 1
         l3loop2: do
          lll = lll + 1
          if (lll>m-1) exit l3loop2
          if (abs(e(lll))<=tol*mu) then
           e(lll) = 0.0_KR2
           mainflag = 1
           exit l3loop2
          endif
          sminlo = sminl
          mu = abs(d(lll+1))*(mu/(mu+abs(e(lll))))
          sminl = min(sminl,mu)
         enddo l3loop2
        endif
       else
! Run convergence test in backward direction.
! First apply standard test to top of matrix.
        if (abs(e(ll))<=abs(tol)*abs(d(ll)) .or. &
           (tol<0.0_KR2 .and. abs(e(ll))<=thresh)) then
         e(ll) = 0.0_KR2
         mainflag = 1
        endif
        if (mainflag==0 .and. tol>=0.0_KR2) then
! If relative accuracy desired, apply convergence criterion backward.
         mu = abs(d(m))
         sminl = mu
         lll = m
         l3loop3: do
          lll = lll - 1
          if (lll<ll) exit l3loop3
          if (abs(e(lll))<=tol*mu) then
           e(lll) = 0.0_KR2
           mainflag = 1
           exit l3loop3
          endif
          sminlo = sminl
          mu = abs(d(lll))*(mu/(mu+abs(e(lll))))
          sminl = min(sminl,mu)
         enddo l3loop3
        endif
       endif
      endif

      if (mainflag==0) then
       oldll = ll
       oldm = m
! Compute shift.  First, test if shifting would ruin relative accuracy,
! and if so set the shift to zero.
       if (tol>=0.0_KR2 .and. n*tol*(sminl/smax)<=max(eps,0.01_KR2*tol)) then
! Use a zero shift to avoid loss of relative accuracy.
        shift = 0.0_KR2
       else
! Compute the shift from 2-by-2 block at end of matrix.
        if (idir==1) then
         sll = abs(d(ll))
         call sing2by2(d(m-1),e(m-1),d(m),shift,r)
        else
         sll = abs(d(m))
         call sing2by2(d(ll),e(ll),d(ll+1),shift,r)
        endif
! Test if shift negligible, and if so set to zero.
        if (sll>0.0_KR2) then
         if ((shift/sll)**2<eps) shift=0.0_KR2
        endif
       endif
! Increment iteration count.
       iter = iter + m - ll
! If shift = 0, do simplified QR iteration,
       if (shift==0.0_KR2) then
        if (idir==1) then
! Chase bulge from top to bottom.
! Save cosines and sines for later singular vector updates.
         cs = 1.0_KR2
         oldcs = 1.0_KR2
         do i = ll,m-1
          call rotate2(d(i)*cs,e(i),cs,sn,r)
          if (i>ll) e(i-1)=oldsn*r
          call rotate2(oldcs*r,d(i+1)*sn,oldcs,oldsn,d(i))
          rwork(i-ll+1) = cs
          rwork(i-ll+1+nm1) = sn
          rwork(i-ll+1+nm12) = oldcs
          rwork(i-ll+1+nm13) = oldsn
         enddo ! i
         h = d(m)*cs
         d(m) = h*oldcs
         e(m-1) = h*oldsn
! Update singular vectors.
         side = "R"
         fwdbwd = "F"
         col = m - ll + 1
         call mmutil1(side,fwdbwd,nru,col,rwork(nm12+1:nm12+col), &
                      rwork(nm13+1:nm13+col),u(:,:,ll:n))
! Test convergence.
         if (abs(e(m-1))<=thresh) e(m-1)=0.0_KR2
        else
! Chase bulge from bottom to top.
! Save cosines and sines for later singular vector updates.
         cs = 1.0_KR2
         oldcs = 1.0_KR2
         do i = m,ll+1,-1
          call rotate2(d(i)*cs,e(i-1),cs,sn,r)
          if (i<m) e(i)=oldsn*r
          call rotate2(oldcs*r,d(i-1)*sn,oldcs,oldsn,d(i))
          rwork(i-ll) = cs
          rwork(i-ll+nm1) = -sn
          rwork(i-ll+nm12) = oldcs
          rwork(i-ll+nm13) = -oldsn
         enddo ! i
         h = d(ll)*cs
         d(ll) = h*oldcs
         e(ll) = h*oldsn
! Update singular vectors.
         side = "R"
         fwdbwd = "B"
         col = m - ll + 1
         call mmutil1(side,fwdbwd,nru,col,rwork(1:col),rwork(n:n-1+col), &
                      u(:,:,ll:n))
! Test convergence.
         if (abs(e(ll))<=thresh) e(ll)=0.0_KR2
        endif
       else
! Use nonzero shift.
        if (idir==1) then
! Chase bulge from top to bottom.
! Save cosines and sines for later singular vector updates.
         f = (abs(d(ll))-shift)*(sign(1.0_KR2,d(ll))+shift/d(ll))
         g = e(ll)
         do i = ll,m-1
          call rotate2(f,g,cosr,sinr,r)
          if (i>ll) e(i-1)=r
          f = cosr*d(i) + sinr*e(i)
          e(i) = cosr*e(i) - sinr*d(i)
          g = sinr*d(i+1)
          d(i+1) = cosr*d(i+1)
          call rotate2(f,g,cosl,sinl,r)
          d(i) = r
          f = cosl*e(i) + sinl*d(i+1)
          d(i+1) = cosl*d(i+1) - sinl*e(i)
          if (i<m-1) then
           g = sinl*e(i+1)
           e(i+1) = cosl*e(i+1)
          endif
          rwork(i-ll+1) = cosr
          rwork(i-ll+1+nm1) = sinr
          rwork(i-ll+1+nm12) = cosl
          rwork(i-ll+1+nm13) = sinl
         enddo ! i
         e(m-1) = f
! Update singular vectors.
         side = "R"
         fwdbwd = "F"
         col = m - ll + 1
         call mmutil1(side,fwdbwd,nru,col,rwork(nm12+1:nm12+col), &
                      rwork(nm13+1:nm13+col),u(:,:,ll:n))
! Test convergence.
         if (abs(e(m-1))<=thresh) e(m-1)=0.0_KR2
        else
! Chase bulge from bottom to top.
! Save cosines and sines for later singular vector updates.
         f = (abs(d(m))-shift)*(sign(1.0_KR2,d(m))+shift/d(m))
         g = e(m-1)
         do i = m,ll+1,-1
          call rotate2(f,g,cosr,sinr,r)
          if (i<m) e(i)=r
          f = cosr*d(i) + sinr*e(i-1)
          e(i-1) = cosr*e(i-1) - sinr*d(i)
          g = sinr*d(i-1)
          d(i-1) = cosr*d(i-1)
          call rotate2(f,g,cosl,sinl,r)
          d(i) = r
          f = cosl*e(i-1) + sinl*d(i-1)
          d(i-1) = cosl*d(i-1) - sinl*e(i-1)
          if (i>ll+1) then
           g = sinl*e(i-2)
           e(i-2) = cosl*e(i-2)
          endif
          rwork(i-ll) = cosr
          rwork(i-ll+nm1) = -sinr
          rwork(i-ll+nm12) = cosl
          rwork(i-ll+nm13) = -sinl
         enddo ! i
         e(ll) = f
! Test convergence.
         if (abs(e(ll))<=thresh) e(ll)=0.0_KR2
! Update singular vectors if desired.
         side = "R"
         fwdbwd = "B"
         col = m - ll + 1
         call mmutil1(side,fwdbwd,nru,col,rwork(1:col),rwork(n:n-1+col), &
                      u(:,:,ll:n))
        endif
       endif
      endif

! QR iteration finished, go back and check convergence.
     enddo mainloop
    endif

! All singular values have converged, so make them positive.
    do i = 1,n
     if (d(i)<0.0_KR2) d(i)=-d(i)
    enddo ! i

! Sort the singular values into decreasing order (insertion sort on singular
! values, but only one transposition per singular vector).
    do i = 1,n-1
! Scan for smallest d(i).
     isub = 1
     smin = d(1)
     do j = 2,n+1-i
      if (d(j)<=smin) then
       isub = j
       smin = d(j)
      endif
     enddo ! j
! Swap singular values and vectors.
     if (isub/=n+1-i) then
      d(isub) = d(n+1-i)
      d(n+1-i) = smin
      do j = 1,nru
       utemp(:) = u(:,j,isub)
       u(:,j,isub) = u(:,j,n+1-i)
       u(:,j,n+1-i) = utemp(:)
      enddo ! j
     endif
    enddo ! i

 end subroutine svdkernel

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine QRdoubleLAPACK(ischur,n,ilo,ihi,h,w,iloz,ihiz,z)
! Apply the standard double-shift QR algorithm to the Hessenberg submatrix
! in rows and columns ilo to ihi.
! This subroutine corresponds to SUBROUTINE ZLAHQR of LAPACK.
! INPUT:
!   ischur=1 means eigenvectors are requested, ischur=0 means only eigenvalues.
!   n is the order of the matrix h.
!   ilo,ihi: It is assumed that h is already upper triangular in rows and
!            columns ihi+1:n, and that h(ilo,ilo-1) = 0 (unless ilo = 1).
!            Subroutine QRdouble works primarily with the Hessenberg submatrix
!            in rows and columns ilo to ihi, but applies transformations to all
!            of h if ischur=1.
!   h is the upper Hessenberg matrix.
!   iloz,ihiz specify the rows of Z to which transformations must be applied if
!             ischur=1.
!   z is only used if ischur=1.  In that case, Z must contain the current matrix
!     Z of transformations accumulated by subroutine evalues.
! OUTPUT:
!   h is only meaningful on exit if ischur=1.  In that case, h is upper
!     triangular in rows and columns ilo:ihi, with any 2-by-2 diagonal blocks
!     in standard form.
!   w contains the computed eigenvalues between ilo and ihi.
!     If ischur=1, the eigenvalues are stored in the same order as on the
!     diagonal of the Schur form returned in h, with w(:,i) = h(:,i,i).
!   z is updated if ischur=1.  Transformations are applied only to the submatrix
!     z(iloz:ihiz,ilo:ihi).

 integer(kind=KI), intent(in)                      :: ischur, n, ilo, ihi, &
                                                      iloz, ihiz
 real(kind=KR2),   intent(inout), dimension(:,:,:) :: h, z
 real(kind=KR2),   intent(out),   dimension(:,:)   :: w

 ! counter variables
 integer(kind=KI)   :: i,j

 ! Lapack variales
 Logical    :: wantt
 Logical    :: wantz
 integer    :: ldh
 integer    :: ldz
 integer    :: info
 complex*16, allocatable, dimension(:,:) :: h_LAPACK
 complex*16, allocatable, dimension(:)   :: w_LAPACK
 complex*16, allocatable, dimension(:,:) :: z_LAPACK

 ! define lapack variables
 ldh = n
 ldz = n

 ! allocate space
 allocate(h_LAPACK(ldh,n))
 allocate(w_LAPACK(n))
 allocate(z_LAPACK(ldz,n))

 if (ischur == 1) then
   wantz = .true.
   wantt = .true.
 else
   wantz = .false.
   wantt = .false.
 endif

 ! copy RL to Lapack notation
 do i=1,ldh
   do j=1,n
     h_LAPACK(i,j) = cmplx(h(1,i,j),h(2,i,j),KR2)
   enddo
 enddo

 do i=1,ldz
   do j=1,n
     z_LAPACK(i,j) = cmplx(z(1,i,j),z(2,i,j),KR2)
   enddo
 enddo

 call zlahqr(wantt, wantz, n, ilo, ihi, h_LAPACK, ldh, w_LAPACK, iloz, ihiz, z_LAPACK, ldz, info)

 ! copy Lapack to RL notation
 do i=1,ldh
   do j=1,n
     h(1,i,j) = real(h_LAPACK(i,j),KR2)
     h(2,i,j) = aimag(h_LAPACK(i,j))
   enddo
 enddo

 do i=1,ldz
   do j=1,n
     z(1,i,j) = real(z_LAPACK(i,j),KR2)
     z(2,i,j) = aimag(z_LAPACK(i,j))
   enddo
 enddo

 do i=1,n
   w(1,i) = real(w_LAPACK(i),KR2)
   w(2,i) = aimag(w_LAPACK(i))
 enddo

 ! deallocate space
 deallocate(h_LAPACK)
 deallocate(w_LAPACK)
 deallocate(z_LAPACK)
 end subroutine QRdoubleLAPACK

 subroutine QRdouble(ischur,n,ilo,ihi,h,w,iloz,ihiz,z)
! Apply the standard double-shift QR algorithm to the Hessenberg submatrix
! in rows and columns ilo to ihi.
! This subroutine corresponds to SUBROUTINE ZLAHQR of LAPACK.
! INPUT:
!   ischur=1 means eigenvectors are requested, ischur=0 means only eigenvalues.
!   n is the order of the matrix h.
!   ilo,ihi: It is assumed that h is already upper triangular in rows and
!            columns ihi+1:n, and that h(ilo,ilo-1) = 0 (unless ilo = 1).
!            Subroutine QRdouble works primarily with the Hessenberg submatrix
!            in rows and columns ilo to ihi, but applies transformations to all
!            of h if ischur=1.
!   h is the upper Hessenberg matrix.
!   iloz,ihiz specify the rows of Z to which transformations must be applied if
!             ischur=1.
!   z is only used if ischur=1.  In that case, Z must contain the current matrix
!     Z of transformations accumulated by subroutine evalues.
! OUTPUT:
!   h is only meaningful on exit if ischur=1.  In that case, h is upper
!     triangular in rows and columns ilo:ihi, with any 2-by-2 diagonal blocks
!     in standard form.
!   w contains the computed eigenvalues between ilo and ihi.
!     If ischur=1, the eigenvalues are stored in the same order as on the
!     diagonal of the Schur form returned in h, with w(:,i) = h(:,i,i).
!   z is updated if ischur=1.  Transformations are applied only to the submatrix
!     z(iloz:ihiz,ilo:ihi).

    integer(kind=KI), intent(in)                      :: ischur, n, ilo, ihi, &
                                                         iloz, ihiz
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: h, z
    real(kind=KR2),   intent(out),   dimension(:,:)   :: w

    integer(kind=KI)               :: i, j, k, l, m, ii, i1, i2, itn, its, &
                                      iflag, jflag, myn, myid
    real(kind=KR2)                 :: ulp, smlnum, tst1, myrsum, s, t2, hRe, &
                                      rtemp
    real(kind=KR2), dimension(2)   :: t, u, x, y, ytemp, tbit, h10, h11, h21, &
                                      h22, h11s, t1, temp, mysum
    real(kind=KR2), dimension(2,2) :: v

  !  if (myid==0) then
   !    print *, 'Entering QRdouble'
   ! endif

    if (n/=0) then
     if (ilo==ihi) then
      w(:,ilo) = h(:,ilo,ilo)
     else
! Set machine-dependent constants for the stopping criterion.
! If norm(h) <= sqrt(ovfl), overflow should not occur.
      ulp = dlamchp
      smlnum = dlamchs/ulp
! i1 and i2 are the indices of the first row and last column of h to which
! transformations must be applied.  If eigenvalues only are being computed,
! i1 and i2 are set inside the main loop.
      if (ischur==1) then
       i1 = 1
       i2 = n
      endif
! itn is the total number of QR iterations allowed.
      itn = 30*(ihi-ilo+1)!TW changed from 30 to 20 2/27/2017 
! The main loop begins here.  i is the loop index and decreases from ihi to ilo
! in steps of 1.  Each iteration of the loop works with the active submatrix in
! rows and columns l to i.  Eigenvalues i+1 to ihi have already converged.
! Either l = ilo, or h(l,l-1) is negligible so that the matrix splits.
      i = ihi
      maindo: do
       if (i<ilo) exit maindo
! Perform QR iterations on rows and columns ilo to i until a submatrix of order
! 1 splits off at the bottom because a subdiagonal element has become
! negligible.
       l = ilo
       its = -1
       itsdo: do
        its = its + 1
        if (its>itn) exit itsdo
! Look for a single small subdiagonal element.
        k = i + 1
        kdo: do
         k = k - 1
         if (k<l+1) exit kdo
         tst1 = abs(h(1,k-1,k-1)) + abs(h(2,k-1,k-1)) &
              + abs(h(1,k,k)) + abs(h(2,k,k))
         if (tst1==0.0_KR2 .and. i+1/=l) then
          do j = l,i
           myrsum = 0.0_KR2
           do ii = l,min(i,j+1)
            myrsum = myrsum + sqrt(h(1,ii,j)**2+h(2,ii,j)**2)
           enddo ! ii
           tst1 = max(tst1,myrsum)
          enddo ! j
         endif
         if (abs(h(1,k,k-1))<=max(ulp*tst1,smlnum)) exit kdo
        enddo kdo
        l = k
        if (l>ilo) then
! h(l,l-1) is negligible.
         h(:,l,l-1) = 0.0_KR2
        endif
! Exit from loop if a submatrix of order 1 has split off.
        if (l<i) then
         iflag = 0
        else
         iflag = 1
         exit itsdo
        endif
! Now the active submatrix is in rows and columns l to i. If eigenvalues only
! are being computed, only the active submatrix need be transformed.
        if (ischur==0) then
         i1 = l
         i2 = i
        endif
        if (its==10 .or. its==20) then
! Exceptional shift.
         s = 0.75_KR2*abs(h(1,i,i-1))
         t(1) = s + h(1,i,i)
         t(2) = h(2,i,i)
        else
! Wilkinson's shift.
         t(:) = h(:,i,i)
         u(:) = h(:,i-1,i)*h(1,i,i-1)
         if (u(1)/=0.0_KR2 .or. u(2)/=0.0_KR2) then
          x(:) = 0.5_KR2*(h(:,i-1,i-1)-t(:))
          ytemp(1) = x(1)**2 - x(2)**2 + u(1)
          ytemp(2) = 2.0_KR2*x(1)*x(2) + u(2)
          call safesqrt(ytemp,y)
          if (x(1)*y(1)+x(2)*y(2)<0.0_KR2) y=-y
          ytemp = x + y
          call compdiv(u,ytemp,tbit)
          t = t - tbit
         endif
        endif
! Look for two consecutive small subdiagonal elements.
        jflag = 0
        m = i
        mdo: do
         m = m - 1
         if (m<l+1) exit mdo
! Determine the effect of starting the single-shift QR iteration at row m, and
! see if this would make h(:,m,m-1) negligible.
         h11(:) = h(:,m,m)
         h22(:) = h(:,m+1,m+1)
         h11s = h11 - t
         h21(:) = h(:,m+1,m)
         s = abs(h11s(1)) + abs(h11s(2)) + sqrt(h21(1)**2+h21(2)**2)
         h11s = h11s/s
         h21 = h21/s
         v(:,1) = h11s(:)
         v(:,2) = h21(:)
         h10(:) = h(:,m,m-1)
         tst1 = (abs(h11s(1))+abs(h11s(2)))*(abs(h11(1))+abs(h11(2)) &
                +abs(h22(1))+abs(h22(2)))
         myrsum = sqrt( (h10(1)*h21(1)-h10(2)*h21(2))**2 &
                     + (h10(1)*h21(2)+h10(2)*h21(1))**2 )
         if (myrsum<=ulp*tst1) then
          jflag = 1
          exit mdo
         endif
        enddo mdo
        if (jflag==0) then
         h11(:) = h(:,l,l)
         h22(:) = h(:,l+1,l+1)
         h11s = h11 - t
         h21(:) = h(:,l+1,l)
         s = abs(h11s(1)) + abs(h11s(2)) + sqrt(h21(1)**2+h21(2)**2)
         h11s = h11s/s
         h21 = h21/s
         v(:,1) = h11s(:)
         v(:,2) = h21(:)
        endif
! Single-shift QR step.
        do k = m,i-1
! The first iteration of this loop determines a reflection g from the vector
! v and applies it from left and right to h, thus creating a nonzero bulge
! below the subdiagonal.
! Each subsequent iteration determines a reflection g to restore the Hessenberg
! form in the (k-1)th column, and thus chases the bulge one step toward the
! bottom of the active submatrix.
! v(:,2) is always real before the call to subroutine makereflector, and hence
! after the call t2 ( = t1*v(:,2) ) is also real.
         if (k>m) then
          v(:,1) = h(:,k,k-1)
          v(:,2) = h(:,k+1,k-1)
         endif
         myn = 2
         call makereflector(myn,v(:,1),v(:,2:2),t1)
         if (k>m) then
          h(:,k,k-1) = v(:,1)
          h(:,k+1,k-1) = 0.0_KR2
         endif
         t2 = t1(1)*v(1,2) - t1(2)*v(2,2)
! Apply g from the left to transform the rows of the matrix in columns k to i2.
         do j = k,i2
          mysum(1) = t1(1)*h(1,k,j) + t1(2)*h(2,k,j) + t2*h(1,k+1,j)
          mysum(2) = t1(1)*h(2,k,j) - t1(2)*h(1,k,j) + t2*h(2,k+1,j)
          h(:,k,j) = h(:,k,j) - mysum(:)
          h(1,k+1,j) = h(1,k+1,j) - mysum(1)*v(1,2) + mysum(2)*v(2,2)
          h(2,k+1,j) = h(2,k+1,j) - mysum(1)*v(2,2) - mysum(2)*v(1,2)
         enddo ! j
! Apply g from the right to transform the columns of the matrix in rows i1 to
! min(k+2,i).
         do j = i1,min(k+2,i)
          mysum(1) = t1(1)*h(1,j,k) - t1(2)*h(2,j,k) + t2*h(1,j,k+1)
          mysum(2) = t1(1)*h(2,j,k) + t1(2)*h(1,j,k) + t2*h(2,j,k+1)
          h(:,j,k) = h(:,j,k) - mysum(:)
          h(1,j,k+1) = h(1,j,k+1) - mysum(1)*v(1,2) - mysum(2)*v(2,2)
          h(2,j,k+1) = h(2,j,k+1) - mysum(2)*v(1,2) + mysum(1)*v(2,2)
         enddo ! j
! Accumulate transformations in the matrix z.
         if (ischur==1) then
          do j = iloz, ihiz
           mysum(1) = t1(1)*z(1,j,k) - t1(2)*z(2,j,k) + t2*z(1,j,k+1)
           mysum(2) = t1(1)*z(2,j,k) + t1(2)*z(1,j,k) + t2*z(2,j,k+1)
           z(:,j,k) = z(:,j,k) - mysum(:)
           z(1,j,k+1) = z(1,j,k+1) - mysum(1)*v(1,2) - mysum(2)*v(2,2)
           z(2,j,k+1) = z(2,j,k+1) - mysum(2)*v(1,2) + mysum(1)*v(2,2)
          enddo ! j
         endif
! If the QR step was started at row m > l because two consecutive small
! subdiagonals were found, then extra scaling must be performed to ensure that
! h(:,m,m-1) remains real.
         if (k==m .and. m>l) then
          temp(1) = 1.0_KR2 - t1(1)
          temp(2) = 0.0_KR2 - t1(2)
          myrsum = sqrt(temp(1)**2+temp(2)**2)
          temp = temp/myrsum
          hRe = h(1,m+1,m)*temp(1) + h(2,m+1,m)*temp(2)
          h(2,m+1,m) = h(2,m+1,m)*temp(1) - h(1,m+1,m)*temp(2)
          h(1,m+1,m) = hRe
          if (m+2<=i) then
           hRe = h(1,m+2,m+1)*temp(1) - h(2,m+2,m+1)*temp(2)
           h(2,m+2,m+1) = h(2,m+2,m+1)*temp(1) + h(1,m+2,m+1)*temp(2)
           h(1,m+2,m+1) = hRe
          endif
          do j = m,i
           if (j/=m+1) then
            if (i2>j) then
             do ii = j+1,i2
              hRe = temp(1)*h(1,j,ii) - temp(2)*h(2,j,ii)
              h(2,j,ii) = temp(1)*h(2,j,ii) + temp(2)*h(1,j,ii)
              h(1,j,ii) = hRe
             enddo ! ii
            endif
            do ii = i1,j-1
             hRe = temp(1)*h(1,ii,j) + temp(2)*h(2,ii,j)
             h(2,ii,j) = temp(1)*h(2,ii,j) - temp(2)*h(1,ii,j)
             h(1,ii,j) = hRe
            enddo ! ii
            if (ischur==1) then
             do ii = iloz,ihiz
              hRe = temp(1)*z(1,ii,j) + temp(2)*z(2,ii,j)
              z(2,ii,j) = temp(1)*z(2,ii,j) - temp(2)*z(1,ii,j)
              z(1,ii,j) = hRe
             enddo ! ii
            endif
           endif
          enddo ! j
         endif
        enddo ! k
! Ensure that h(:,i,i-1) is real.
        temp(:) = h(:,i,i-1)
        if (temp(2)/=0.0_KR2) then
         rtemp = sqrt(temp(1)**2+temp(2)**2)
         h(1,i,i-1) = rtemp
         h(2,i,i-1) = 0.0_KR2
         temp = temp/rtemp
         if (i2>i) then
          do ii = i+1,i2
           hRe = temp(1)*h(1,i,ii) + temp(2)*h(2,i,ii)
           h(2,i,ii) = temp(1)*h(2,i,ii) - temp(2)*h(1,i,ii)
           h(1,i,ii) = hRe
          enddo ! ii
         endif
         do ii = i1,i-1
          hRe = temp(1)*h(1,ii,i) - temp(2)*h(2,ii,i)
          h(2,ii,i) = temp(1)*h(2,ii,i) + temp(2)*h(1,ii,i)
          h(1,ii,i) = hRe
         enddo ! ii
         if (ischur==1) then
          do ii = iloz,ihiz
           hRe = temp(1)*z(1,ii,i) - temp(2)*z(2,ii,i)
           z(2,ii,i) = temp(1)*z(2,ii,i) + temp(2)*z(1,ii,i)
           z(1,ii,i) = hRe
          enddo ! ii
         endif
        endif
       enddo itsdo
! Failure to converge in remaining number of iterations.
       if (iflag==0) then
        open(unit=8,file="PSEUDOLAPACK.ERROR",action="write",status="replace" &
            ,form="formatted")
         write(unit=8,fmt=*) "Failed to converge in QRdouble."
        close(unit=8,status="keep")
        stop
       endif
! h(:,i,i-1) is negligible: one eigenvalue has converged.
       w(:,i) = h(:,i,i)
! Decrement number of remaining iterations, and return to start of the main
! loop with new value of i.
       itn = itn - its
     !  if (myid==0) then !TW 2/27/2017
      !    print *, 'On the', itn, 'th iteration of QRdouble'!TW 2/27/2017
       !endif !TW 2/27/2017
       i = l - 1
      enddo maindo
     endif
    endif

 end subroutine QRdouble

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine makereflector(n,alpha,x,tau)
! Generate a complex elementary reflector H of order n such that
! H^dagger * / alpha \ = / beta \, and H^dagger * H = 1.
!            \   x   /   \  0   /
! where alpha is a complex scalar,
!       beta is a real scalar,
!       x is an (n-1)-component complex vector.
! H is written in the form H = 1 - tau * / 1 \ * ( 1 v^dagger ),
!                                        \ v /
! where tau is a complex scalar,
!       v is an (n-1)-component complex vector.
! Note that H is not Hermitian.
! If x=0 and alpha is real then the code chooses tau=0, otherwise
! 1 <= Re(tau) <= 2 and abs(tau-1) <= 1.
! This subroutine corresponds to SUBROUTINE ZLARFG of LAPACK.
! INPUT:
!   n,
!   alpha(:) with 2 entries expected,
!   x(:,:) with (2,n-1) entries expected.
! OUTPUT:
!   alpha() is overwritten with the value of beta,
!   x() is overwritten with the vector v,
!   tau.

    integer(kind=KI), intent(in)                    :: n
    real(kind=KR2),   intent(inout), dimension(:)   :: alpha
    real(kind=KR2),   intent(inout), dimension(:,:) :: x
    real(kind=KR2),   intent(out),   dimension(:)   :: tau

    real(kind=KR2)               :: safmin, xnorm, beta, rsafmn, xRe
    real(kind=KR2), dimension(2) :: a, b
    integer(kind=KI)             :: knt, i

    safmin = dlamchs/dlamche
    if (n<1) then
     tau = 0.0_KR2
    else
     xnorm = twonorm(x,n-1)
     if (xnorm==0.0_KR2 .and. alpha(2)==0.0_KR2) then
      tau = 0.0_KR2
     else
      beta = -sign(pythag3(alpha(1),alpha(2),xnorm),alpha(1)) 
      rsafmn = 1.0_KR2/safmin
      if (abs(beta)<safmin) then
! xnorm and beta may be inaccurate; scale x and recompute them.
       knt = 0
       bigloop: do
        knt = knt + 1
        do i = 1,n-1
         x(:,i) = rsafmn*x(:,i)
        enddo ! i
        beta = rsafmn*beta
        alpha = rsafmn*alpha
        if (abs(beta)>=safmin) exit bigloop
       enddo bigloop
! The new beta is no larger than 1 and no smaller than safmin.
       xnorm = twonorm(x,n-1)
       beta = -sign(pythag3(alpha(1),alpha(2),xnorm),alpha(1)) 
       tau(1) = (beta-alpha(1))/beta
       tau(2) = -alpha(2)/beta
       a(1) = 1.0_KR2
       a(2) = 0.0_KR2
       b(1) = alpha(1) - beta
       b(2) = alpha(2)
       call compdiv(a,b,alpha)
       do i = 1,n-1
        xRe = alpha(1)*x(1,i) - alpha(2)*x(2,i)
        x(2,i) = alpha(1)*x(2,i) + alpha(2)*x(1,i)
        x(1,i) = xRe
       enddo ! i
! If alpha is subnormal, it may lose relative accuracy.
       alpha(1) = beta
       alpha(2) = 0.0_KR2
       do i = 1,knt
        alpha(:) = safmin*alpha(:)
       enddo ! i
      else
       tau(1) = (beta-alpha(1))/beta
       tau(2) = -alpha(2)/beta
       a(1) = 1.0_KR2
       a(2) = 0.0_KR2
       b(1) = alpha(1) - beta
       b(2) = alpha(2)
       call compdiv(a,b,alpha)
       do i = 1,n-1
        xRe = alpha(1)*x(1,i) - alpha(2)*x(2,i)
        x(2,i) = alpha(1)*x(2,i) + alpha(2)*x(1,i)
        x(1,i) = xRe
       enddo ! i
       alpha(1) = beta
       alpha(2) = 0.0_KR2
      endif
     endif
    endif

 end subroutine makereflector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine multreflector(side,m,n,v,tau,c,workspace)
! Apply a complex elementary reflector h to a complex m by n matrix c,
! from either the left or the right.  h is represented in the form
! h = 1 - tau * v * v'
! where tau is a complex scalar and v is a complex vector.
! If tau = 0, then h is taken to be the unit matrix.
! This subroutine corresponds to SUBROUTINE ZLARFX of LAPACK.
! INPUT:
!   side="L" will construct h*c; side="R" will construct c*h.
!   m is the number of rows in the matrix c.
!   n is the number of columns in the matrix c.
!   v(1,:),v(2,:) are the Re,Im parts of vector v in the representation of h.
!       expected size: v(:,1:m) if side="L"; v(:,1:n) if side="R".
!   tau(1),tau(2) are the Re,Im parts of tau in the representation of h.
!   c(1,:,:),c(2,:,:) are the Re,Im parts of the input matrix c.
!   workspace(1,:),workspace(2,:) are the Re,Im parts of a scratch vector.
!       expected size: (:,1:m) if side="L"; (:,1:n) if side="R".
! OUTPUT:
!   c() is overwritten by the solution matrix h*c or c*h.

    character(len=*), intent(in)                      :: side
    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(in),    dimension(:,:)   :: v
    real(kind=KR2),   intent(in),    dimension(:)     :: tau
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: c
    real(kind=KR2),   intent(inout), dimension(:,:)   :: workspace

    integer(kind=KI)             :: dagger
    real(kind=KR2), dimension(2) :: one, zero

    one(1) = 1.0_KR2
    one(2) = 0.0_KR2
    zero(:) = 0.0_KR2
    if (tau(1)/=0.0_KR2 .or. tau(2)/=0.0_KR2) then
     if (side=="L") then
      dagger = 1
      call mvutil1(dagger,m,n,one,c,v,zero,workspace)
      call mvutil2(m,n,-tau,v,workspace,c)
     else
      dagger = 0
      call mvutil1(dagger,m,n,one,c,v,zero,workspace)
      call mvutil2(m,n,-tau,workspace,v,c)
     endif
    endif

 end subroutine multreflector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil1(dagger,m,n,alpha,a,x,beta,y)
! Perform one of the following matrix-vector operations:
! y = beta*y + alpha*a*x        [dagger=0],
! y = beta*y + alpha*a^dagger*x [dagger=1],
! where alpha and beta are complex scalars, x and y are complex vectors
! and "a" is a complex m by n matrix.
! This subroutine corresponds to SUBROUTINE ZGEMV of LAPACK.
! INPUT:
!   dagger is a character that specifies "a" or a^dagger as defined above.
!   m is the number of rows in the matrix (cannot be negative).
!   n is the number of columns in the matrix (cannot be negative).
!   alpha defines the complex scalar defined above.
!   "a" is the matrix defined above.
!   x is the vector defined above.
!   beta defines the complex scalar defined above.
!   y is the original vector defined above.  It will get destroyed.
! OUTPUT:
!   y is over-written with the solution.

    integer(kind=KI), intent(in)                      :: dagger, m, n
    real(kind=KR2),   intent(in),    dimension(:)     :: alpha, beta
    real(kind=KR2),   intent(in),    dimension(:,:,:) :: a
    real(kind=KR2),   intent(in),    dimension(:,:)   :: x
    real(kind=KR2),   intent(inout), dimension(:,:)   :: y

    integer(kind=KI)             :: lenx, leny, i, j
    real(kind=KR2), dimension(2) :: temp
    real(kind=KR2)               :: yRe

    if (m>0 .and. n>0 .and. (alpha(1)/=0.0_KR2 .or. alpha(2)/=0.0_KR2 .or. &
                             beta(1)/=1.0_KR2 .or. beta(2)/=0.0_KR2)) then
     if (dagger==0) then
      lenx = n
      leny = m
     else
      lenx = m
      leny = n
     endif
! First construct y = beta*y.
     if (beta(1)/=1.0_KR2 .or. beta(2)/=0.0_KR2) then
      if (beta(1)==0.0_KR2 .and. beta(2)==0.0_KR2) then
       y = 0.0_KR2
      else
       do i = 1,leny
        yRe = beta(1)*y(1,i) - beta(2)*y(2,i)
        y(2,i) = beta(1)*y(2,i) + beta(2)*y(1,i)
        y(1,i) = yRe
       enddo ! i
      endif
     endif
! Now construct y = y + alpha*a*x, if requested.
     if (alpha(1)/=0.0_KR2 .or. alpha(2)/=0.0_KR2) then
      if (dagger==0) then
       do j = 1,n
        if (x(1,j)/=0.0_KR2 .or. x(2,j)/=0.0_KR2) then
         temp(1) = alpha(1)*x(1,j) - alpha(2)*x(2,j)
         temp(2) = alpha(1)*x(2,j) + alpha(2)*x(1,j)
         do i = 1,m
          y(1,i) = y(1,i) + temp(1)*a(1,i,j) - temp(2)*a(2,i,j)
          y(2,i) = y(2,i) + temp(1)*a(2,i,j) + temp(2)*a(1,i,j)
         enddo ! i
        endif
       enddo ! j
! Now construct y = y + alpha*a^dagger*x, if requested.
      else
       do j = 1,n
        temp = 0.0_KR2
        do i = 1,m
         temp(1) = temp(1) + a(1,i,j)*x(1,i) + a(2,i,j)*x(2,i)
         temp(2) = temp(2) + a(1,i,j)*x(2,i) - a(2,i,j)*x(1,i)
        enddo ! i
        y(1,j) = y(1,j) + alpha(1)*temp(1) - alpha(2)*temp(2)
        y(2,j) = y(2,j) + alpha(1)*temp(2) + alpha(2)*temp(1)
       enddo ! j
      endif
     endif
    endif

 end subroutine mvutil1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil2(m,n,alpha,x,y,a)
! Perform the following matrix-vector operation:
! a = a + alpha*x*y^dagger
! where alpha is a complex scalar, x and y are complex vectors
! and "a" is a complex m by n matrix.
! This subroutine corresponds to SUBROUTINE ZGERC of LAPACK.

    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(in),    dimension(:)     :: alpha
    real(kind=KR2),   intent(in),    dimension(:,:)   :: x, y
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a

    integer(kind=KI)             :: i, j
    real(kind=KR2), dimension(2) :: temp

    if (m/=0 .and. n/=0 .and. (alpha(1)/=0.0_KR2 .or. alpha(2)/=0.0_KR2) ) then
     do j = 1,n
      if (y(1,j)/=0.0_KR2 .or. y(2,j)/=0.0_KR2) then
       temp(1) = alpha(1)*y(1,j) + alpha(2)*y(2,j)
       temp(2) = alpha(2)*y(1,j) - alpha(1)*y(2,j)
       do i = 1,m
        a(1,i,j) = a(1,i,j) + x(1,i)*temp(1) - x(2,i)*temp(2)
        a(2,i,j) = a(2,i,j) + x(1,i)*temp(2) + x(2,i)*temp(1)
       enddo ! i
      endif
     enddo ! j
    endif

 end subroutine mvutil2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil3(n,a,x,myscale,cnorm)
! Solve the triangular system a*x = s*b with scaling to prevent overflow.
! Here "a" is an upper triangular matrix, x and b are n-component vectors, and
! s is a scaling factor, usually less than or equal to 1, chosen so that the
! components of x will be less than the overflow threshold.
! This subroutine corresponds to SUBROUTINE ZLATRS of LAPACK.
! INPUT:
!   n is the order of the matrix "a".
!   a is the upper triangular matrix.  The leading n by n upper triangular part
!     of the array contains the upper triangular matrix, and the strictly lower
!     triangular part is not referenced.
!   x is the right hand side b of the triangular system.
!   myscale is the scaling factor s for the triangular system a*x = s*b.
!           If myscale = 0, the matrix "a" is singular or badly scaled, and
!           the vector x is an exact or approximate solution to a*x = 0.
!   cnorm(j) contains the norm of the off-diagonal part of the j-th column of
!            "a".  cnorm(j) must be greater than or equal to the infinity-norm.
! OUTPUT:
!   x is overwritten by the solution vector x.
!   cnorm() is updated.

    integer(kind=KI), intent(in)                      :: n
    real(kind=KR2),   intent(in),    dimension(:,:,:) :: a
    real(kind=KR2),   intent(inout), dimension(:,:)   :: x
    real(kind=KR2),   intent(out)                     :: myscale
    real(kind=KR2),   intent(inout), dimension(:)     :: cnorm

    integer(kind=KI)             :: imax, i, j, ii
    real(kind=KR2), dimension(2) :: tjjs, xold
    real(kind=KR2)               :: smlnum, bignum, dmax, trydmax, tscal, xj, &
                                    xmax, tryxmax, xbnd, grow, tjj, myrec, tmax

    if (n/=0) then
! Determine machine dependent parameters to control overflow.
     smlnum = dlamchs/dlamchp
     bignum = 1.0_KR2/smlnum
     myscale = 1.0_KR2
! Scale the column norms by tscal if the maximum element in cnorm is greater
! than bignum/2.
     imax = 1
     if (n>1) then
      dmax = abs(cnorm(1))
      do i = 2,n
       trydmax = abs(cnorm(i))
       if (trydmax>dmax) then
        imax = i
        dmax = trydmax
       endif
      enddo ! i
     endif
     tmax = cnorm(imax)
     if (tmax<=0.5_KR2*bignum) then
      tscal = 1.0_KR2
     else
      tscal = 0.5_KR2/(smlnum*tmax)
      cnorm = tscal*cnorm
     endif
! Compute a bound on the computed solution vector to see if subroutine mvutil4
! can be used.
     xmax = 0.0_KR2
     do j = 1,n
      xmax = max(xmax,abs(0.5_KR2*x(1,j))+abs(0.5_KR2*x(2,j)))
     enddo ! j
     xbnd = xmax
! If "a" is non-unit triangular, then compute grow = 1/g(j) and xbnd = 1/m(j).
! Initially, g(0) = max{x(i), i=1,...,n}.
     if (tscal==1.0_KR2) then
      grow = 0.5_KR2/max(xbnd,smlnum)
      xbnd = grow
      j = n + 1
      jdo: do
       j = j - 1
       if (j<1 .or. grow<=smlnum) exit jdo
       tjjs(:) = a(:,j,j)
       tjj = abs(tjjs(1)) + abs(tjjs(2))
       if (tjj>=smlnum) then
        xbnd = min(xbnd,min(1.0_KR2,tjj)*grow)
       else
        xbnd = 0.0_KR2
       endif
       if (tjj+cnorm(j)>=smlnum) then
        grow = grow*(tjj/(tjj+cnorm(j)))
       else
        grow = 0.0_KR2
       endif
      enddo jdo
      grow = xbnd
     else
      grow = 0.0_KR2
     endif
! Use subroutine mvutil4 if the reciprocal of the bound on elements of x is
! not too small.
     if (grow*tscal>smlnum) then
      call mvutil4(n,a,x)
     else
! An alternative to subroutine mvutil4, with scaling of intermediate results.
      if (xmax>0.5_KR2*bignum) then
! Scale x so that its components are less than or equal to bignum in absolute
! value.
       myscale = (0.5_KR2*bignum)/xmax
       do i = 1,n
        x(:,i) = myscale*x(:,i)
       enddo ! i
       xmax = bignum
      else
       xmax = 2.0_KR2*xmax
      endif
! Solve a*x = b.
      do j = n,1,-1
!-compute x(:,j) = b(:,j)/a(:,j,j), scaling x if necessary.
       xj = abs(x(1,j)) + abs(x(2,j))
       tjjs(:) = a(:,j,j)*tscal
       tjj = abs(tjjs(1)) + abs(tjjs(2))
       if (tjj>smlnum) then
!-abs(a(j,j)) > smlnum:
        if (tjj<1.0_KR2) then
         if (xj>tjj*bignum) then
!-scale x by 1/b(j).
          myrec = 1.0_KR2/xj
          x = myrec*x
          myscale = myscale*myrec
          xmax = xmax*myrec
         endif
        endif
        xold(:) = x(:,j)
        call compdiv(xold(:),tjjs,x(:,j))
        xj = abs(x(1,j)) + abs(x(2,j))
       elseif (tjj>0.0_KR2) then
!-0<abs(a(j,j)) <= smlnum:
        if (xj>tjj*bignum) then
!-scale x by (1/abs(x(j)))*abs(a(j,j))*bignum to avoid overflow when dividing
! by a(j,j).
         myrec = (tjj*bignum)/xj
         if (cnorm(j)>1.0_KR2) then
!-scale by 1/cnorm(j) to avoid overflow when multiplying x(j) times column j.
          myrec = myrec/cnorm(j)
         endif
         x = myrec*x
         myscale = myscale*myrec
         xmax = xmax*myrec
        endif
        xold(:) = x(:,j)
        call compdiv(xold(:),tjjs,x(:,j))
        xj = abs(x(1,j)) + abs(x(2,j))
       else
!-a(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and myscale = 0, and compute a
! solution to a*x = 0.
        x = 0.0_KR2
        x(1,j) = 1.0_KR2
        xj = 1.0_KR2
        myscale = 0.0_KR2
        xmax = 0.0_KR2
       endif
!-scale x if necessary to avoid overflow when adding a multiple of column j.
       if (xj>1.0_KR2) then
        myrec = 1.0_KR2/xj
        if (cnorm(j)>(bignum-xmax)*myrec) then
!-scale x by 1/(2*abs(x(j))).
         myrec = myrec*0.5_KR2
         x = myrec*x
         myscale = myscale*myrec
        endif
       elseif (xj*cnorm(j)>bignum-xmax) then
!-scale x by 1/2.
        x = 0.5_KR2*x
        myscale = myscale*0.5_KR2
       endif
       if (j>1) then
! Compute the update x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
        do i = 1,j-1
         x(:,i) = x(:,i) - tscal*(x(1,j)*a(1,i,j)-x(2,j)*a(2,i,j))
        enddo ! i
        i = 1
        if (j>2) then
         xmax = abs(x(1,1)) + abs(x(2,1))
         do ii = 2,j-1
          tryxmax = abs(x(1,i)) + abs(x(2,i))
          if (tryxmax>xmax) then
           i = ii
           xmax = tryxmax
          endif
         enddo ! ii
        endif
       endif
      enddo ! j
      myscale = myscale/tscal
     endif
! Scale the column norms by 1/tscal for return.
     if (tscal/=1.0_KR2) then
      cnorm = cnorm/tscal
     endif
    endif

 end subroutine mvutil3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil4(n,a,x)
! Solve the system of equations a*x = b, where b and x are n-component vectors
! and "a" is an n by n upper triangular matrix.
! No test for singularity or near-singularity is included in this routine.
! Such tests must be performed before calling this routine.
! This subroutine corresponds to SUBROUTINE ZTRSV of LAPACK.
! INPUT:
!   n is the order of the matrix "a".
!   "a" is the n by n upper triangular matrix.  The strictly lower triangular
!       part of "a" is not referenced.
!   x is the n-component righthand side vector b.
! OUTPUT:
!   x is overwritten with the solution vector x.

    integer(kind=KI), intent(in)                      :: n
    real(kind=KR2),   intent(in),    dimension(:,:,:) :: a
    real(kind=KR2),   intent(inout), dimension(:,:)   :: x

    integer(kind=KI)             :: i, j
    real(kind=KR2), dimension(2) :: temp

    if (n>0) then
     do j = n,1,-1
      if (x(1,j)/=0.0_KR2 .or. x(2,j)/=0.0_KR2) then
       call compdiv(x(:,j),a(:,j,j),temp)
       x(:,j) = temp(:)
       do i = j-1,1,-1
        x(1,i) = x(1,i) - temp(1)*a(1,i,j) + temp(2)*a(2,i,j)
        x(2,i) = x(2,i) - temp(1)*a(2,i,j) - temp(2)*a(1,i,j)
       enddo ! i
      endif
     enddo ! j
    endif

 end subroutine mvutil4

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil5(m,n,alpha,x,y,a)
! Perform the following matrix-vector operation:
! a = a + alpha*x*y^T
! where alpha is a complex scalar, x and y are complex vectors
! and "a" is a complex m by n matrix.
! This subroutine corresponds to SUBROUTINE ZGERU of LAPACK.

    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(in),    dimension(:)     :: alpha
    real(kind=KR2),   intent(in),    dimension(:,:)   :: x, y
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a

    integer(kind=KI)             :: i, j
    real(kind=KR2), dimension(2) :: temp

    if (m/=0 .and. n/=0 .and. (alpha(1)/=0.0_KR2 .or. alpha(2)/=0.0_KR2) ) then
     do j = 1,n
      if (y(1,j)/=0.0_KR2 .or. y(2,j)/=0.0_KR2) then
       temp(1) = alpha(1)*y(1,j) - alpha(2)*y(2,j)
       temp(2) = alpha(2)*y(1,j) + alpha(1)*y(2,j)
       do i = 1,m
        a(1,i,j) = a(1,i,j) + x(1,i)*temp(1) - x(2,i)*temp(2)
        a(2,i,j) = a(2,i,j) + x(1,i)*temp(2) + x(2,i)*temp(1)
       enddo ! i
      endif
     enddo ! j
    endif

 end subroutine mvutil5

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvutil6(n,nrhs,a,ipiv,b)
! Solve the system of equations a*x = b using the LU factorization computed by
! subroutine LUfactor.
! This subroutine corresponds to SUBROUTINE ZGETRS of LAPACK.
! INPUT:
!   n is the order of the matrix "a".
!   nrhs is the number of right hand sides, i.e., the number of columns of the
!        matrix b.
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the factors L and U from the
!                     factorization a=P*L*U as computed by subroutine LUfactor.
!   ipiv contains the pivot indices from subroutine LUfactor; for 1<=i<=N,
!        row i of the matrix was interchanged with row ipiv(i).
!   b(1,:,:),b(2,:,:) are the Re,Im parts of the right-hand side matrix "b".
! OUTPUT:
!   b(1,:,:),b(2,:,:) are overwritten with the Re,Im parts of the solution "x".

    integer(kind=KI), intent(in)                      :: n, nrhs
    real(kind=KR2),   intent(in),    dimension(:,:,:) :: a
    integer(kind=KI), intent(in),    dimension(:)     :: ipiv
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: b

    integer(kind=KI)             :: nmin
    real(kind=KR2), dimension(2) :: alpha
    character(len=2)             :: tridiag

! A useful definition.
    alpha(1) = 1.0_KR2
    alpha(2) = 0.0_KR2

! Apply row interchanges to the right hand sides.
    nmin = 1
    call rowswap(nrhs,b,nmin,n,ipiv)

! Solve L*x=b, overwriting b with x.
    tridiag = "LU"
    call mmutil2(tridiag,n,nrhs,alpha,a,b)

! Solve U*x=b, overwriting b with x.
    tridiag = "UN"
    call mmutil2(tridiag,n,nrhs,alpha,a,b)

 end subroutine mvutil6

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mmutil1(side,fwdbwd,m,n,c,s,a)
! Perform the transformation a=p*a or a=a*p^T, where "a" is an m by n complex
! matrix and p is an orthogonal matrix consisting of a sequence of plane
! rotations.
! This subroutine corresponds to SUBROUTINE ZLASR of LAPACK.
! INPUT:
!   side: If side="L" then compute a=p*a.
!         If side="R" then compute a=a*p^T.
!   fwdbwd: If fwdbwd="F" then p = p(z-1)*...*p(2)*p(1).
!           If fwdbwd="B" then p = p(1)*p(2)*...*p(z-1).
!   m is the number of rows in the matrix "a".
!   n is the number of columns in the matrix "a".
!   c(),s() contain the cosine and sine that define the matrix p(k).
!           The 2 by 2 plane rotation part of the matrix p(k), called r(k),
!           is assumed to be of the form r(k) = / c(k)  s(k) \.
!                                               \-s(k)  c(k) /
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the m by n matrix.
! OUTPUT:
!   a() is overwritten by p*a if side="R" or by a*p^T if side="L".

    character(len=*), intent(in)                      :: side, fwdbwd
    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(in),    dimension(:)     :: c, s
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a

    integer(kind=KI)             :: i, j
    real(kind=KR2)               :: ctemp, stemp
    real(kind=KR2), dimension(2) :: temp

    if (side=="L") then
     if (fwdbwd=="F") then
      do j = 1,m-1
       ctemp = c(j)
       stemp = s(j)
       if (ctemp/=1.0_KR2 .or. stemp/=0.0_KR2) then
        do i = 1,n
         temp(:) = a(:,j+1,i)
         a(:,j+1,i) = ctemp*temp(:) - stemp*a(:,j,i)
         a(:,j,i) = stemp*temp(:) + ctemp*a(:,j,i)
        enddo ! i
       endif
      enddo ! j
     elseif (fwdbwd=="B") then
      do j = m-1,1,-1
       ctemp = c(j)
       stemp = s(j)
       if (ctemp/=1.0_KR2 .or. stemp/=0.0_KR2) then
        do i = 1,n
         temp(:) = a(:,j+1,i)
         a(:,j+1,i) = ctemp*temp(:) - stemp*a(:,j,i)
         a(:,j,i) = stemp*temp(:) + ctemp*a(:,j,i)
        enddo ! i
       endif
      enddo ! j
     endif
    elseif (side=="R") then
     if (fwdbwd=="F") then
      do j = 1,n-1
       ctemp = c(j)
       stemp = s(j)
       if (ctemp/=1.0_KR2 .or. stemp/=0.0_KR2) then
        do i = 1,m
         temp(:) = a(:,i,j+1)
         a(:,i,j+1) = ctemp*temp(:) - stemp*a(:,i,j)
         a(:,i,j) = stemp*temp(:) + ctemp*a(:,i,j)
        enddo ! i
       endif
      enddo ! j
     elseif (fwdbwd=="B") then
      do j = n-1,1,-1
       ctemp = c(j)
       stemp = s(j)
       if (ctemp/=1.0_KR2 .or. stemp/=0.0_KR2) then
        do i = 1,m
         temp(:) = a(:,i,j+1)
         a(:,i,j+1) = ctemp*temp(:) - stemp*a(:,i,j)
         a(:,i,j) = stemp*temp(:) + ctemp*a(:,i,j)
        enddo ! i
       endif
      enddo ! j
     endif
    endif

 end subroutine mmutil1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mmutil2(tridiag,m,n,alpha,a,b)
! Solve the matrix equation a*x = alpha*b, where alpha is a scalar, x and b are
! m by n matrices, and "a" is an upper or lower triangular matrix.
! The matrix x is overwritten on x.
! This subroutine corresponds to SUBROUTINE ZTRSM of LAPACK.
! INPUT:
!   tridiag = "LU" means that "a" is lower triangular and unit triangular.
!   tridiag = "UN" means that "a" is upper triangular and non-unit triangular.
!   m is the number of rows in the matrix "b".
!   n is the number of columns in the matrix "b".
!   alpha is the complex scalar, alpha.  When alpha is zero then "a" is not
!         referenced and b need not be set before entry.
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the matrix "a".
!                     The third argument runs from 1 to m.
!                     When tridiag="LU", the leading m by m lower triangular
!                     part of the array "a" must contain the lower triangular
!                     matrix and neither the diagonal nor the strictly upper
!                     triangular part of "a" is referenced; the diagonal is
!                     assumed to be unity.
!                     When tridiag="UN", the leading m by m upper triangular
!                     part of the array "a" must contain the upper triangular
!                     matrix and the strictly lower triangular part of "a" is
!                     not referenced.
!   b(1,:,:),b(2,:,:) are the Re,Im parts of the matrix b.
! OUTPUT:
!   b(1,:,:),b(2,:,:) are overwritten by the solution matrix x.

    character(len=*), intent(in)                      :: tridiag
    integer(kind=KI), intent(in)                      :: m, n
    real(kind=KR2),   intent(in),    dimension(:)     :: alpha
    real(kind=KR2),   intent(in),    dimension(:,:,:) :: a
    real(kind=KR2),   intent(inout), dimension(:,:,:) :: b

    integer(kind=KI) :: i, j, k
    real(kind=KR2)   :: temp, denom

    if (alpha(1)==0.0_KR2 .and. alpha(2)==0.0_KR2) then
     do j = 1,n
      do i = 1,m
       b(:,i,j) = 0.0_KR2
      enddo ! i
     enddo ! j
    else
     if (tridiag=="UN") then
      do j = 1,n
       if (alpha(1)/=1.0_KR2 .or. alpha(2)/=0.0_KR2) then
        do i = 1,m
         temp = alpha(1)*b(1,i,j) - alpha(2)*b(2,i,j)
         b(2,i,j) = alpha(1)*b(2,i,j) + alpha(2)*b(1,i,j)
         b(1,i,j) = temp
        enddo ! i
       endif
       do k = m,1,-1
        if (b(1,k,j)/=0.0_KR2 .or. b(2,k,j)/=0.0_KR2) then
         denom = 1.0_KR2/(a(1,k,k)**2+a(2,k,k)**2)
         temp = denom*(b(1,k,j)*a(1,k,k)+b(2,k,j)*a(2,k,k))
         b(2,k,j) = denom*(b(2,k,j)*a(1,k,k)-b(1,k,j)*a(2,k,k))
         b(1,k,j) = temp
         do i = 1,k-1
          temp = b(1,i,j) - b(1,k,j)*a(1,i,k) + b(2,k,j)*a(2,i,k)
          b(2,i,j) = b(2,i,j) - b(1,k,j)*a(2,i,k) - b(2,k,j)*a(1,i,k)
          b(1,i,j) = temp
         enddo ! i
        endif
       enddo ! k
      enddo ! j
     else
      do j = 1,n
       if (alpha(1)/=1.0_KR2 .or. alpha(2)/=0.0_KR2) then
        do i = 1,m
         temp = alpha(1)*b(1,i,j) - alpha(2)*b(2,i,j)
         b(2,i,j) = alpha(1)*b(2,i,j) + alpha(2)*b(1,i,j)
         b(1,i,j) = temp
        enddo ! i
       endif
       do k = 1,m
        if (b(1,k,j)/=0.0_KR2 .or. b(2,k,j)/=0.0_KR2) then
         do i = k+1,m
          temp = b(1,i,j) - b(1,k,j)*a(1,i,k) + b(2,k,j)*a(2,i,k)
          b(2,i,j) = b(2,i,j) - b(1,k,j)*a(2,i,k) - b(2,k,j)*a(1,i,k)
          b(1,i,j) = temp
         enddo ! i
        endif
       enddo ! k
      enddo ! j
     endif
    endif

 end subroutine mmutil2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine rowswap(n,a,k1,k2,ipiv)
! Perform a series of row interchanges on the matrix "a".
! One row interchange is initiated for each of rows k1 through k2 of "a".
! This subroutine corresponds to SUBROUTINE ZLASWP of LAPACK.
! INPUT:
!   n is the number of columns in the matrix "a".
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the matrix of column dimension "a".
!   k1 is the first element of ipiv for which a row interchange will be done.
!   k2 is the last element of ipiv for which a row interchange will be done.
!   ipiv is the vector of pivot indices.  Only the elements in positions
!        k1 through k2 of ipiv are accessed.  ipiv(l) implies that rows k and l
!        are to be interchanged.
! OUTPUT:
!   a(1,:,:),a(2,:,:) are the Re,Im parts of the permuted matrix.

    real(kind=KR2),   intent(inout), dimension(:,:,:) :: a
    integer(kind=KI), intent(in)                      :: n, k1, k2
    integer(kind=KI), intent(in),    dimension(:)     :: ipiv

    integer(kind=KI)             :: n32, i, j, k, ix, ip
    real(kind=KR2), dimension(2) :: temp

    n32 = (n/32)*32

    if (n32/=0) then
     do j = 1,n32,32
      ix = k1
      do i = k1,k2
       ip = ipiv(ix)
       if (ip/=i) then
        do k = j,j+31
         temp(:) = a(:,i,k)
         a(:,i,k) = a(:,ip,k)
         a(:,ip,k) = temp(:)
        enddo ! k
       endif
       ix = ix + 1
      enddo ! i
     enddo ! j
    endif

    if (n32/=n) then
     n32 = n32 + 1
     ix = k1
     do i = k1,k2
      ip = ipiv(ix)
      if (ip/=i) then
       do k = n32,n
        temp(:) = a(:,i,k)
        a(:,i,k) = a(:,ip,k)
        a(:,ip,k) = temp(:)
       enddo ! k
      endif
      ix = ix + 1
     enddo ! i
    endif

 end subroutine rowswap

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine svd2by2(f,g,h,ssmin,ssmax,snr,csr,snl,csl)
! Compute the singular value decomposition of a 2 by 2 triangular matrix,
!     / f  g \
!     \ 0  h /.
!  to arrive at
!     / csl  snl \ / f  g \ / csr -snr \  =  / ssmax   0   \
!     \-snl  csl / \ 0  h / \ snr  csr /     \  0    ssmin /.
! This subroutine corresponds to SUBROUTINE DLASV2 of LAPACK.
! INPUT:
!   f is the (1,1) element of the 2 by 2 matrix.
!   g is the (1,2) element of the 2 by 2 matrix.
!   h is the (2,2) element of the 2 by 2 matrix.
! OUTPUT:
!   abs(ssmin) is the smaller singular value.
!   abs(ssmax) is the larger singular value.
!   (csl,snl) is a unit left singular vector for the singular value abs(ssmax).
!   (csr,snr) is a unit right singular vector for the singular value abs(ssmax).

    real(kind=KR2), intent(in)  :: f, g, h
    real(kind=KR2), intent(out) :: ssmin, ssmax, snr, csr, snl, csl

    logical          :: gasmal, swap
    integer(kind=KI) :: pmax
    real(kind=KR2)   :: a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m, mm, r, &
                        s, slt, srt, t, temp, tsign, tt

    ft = f
    fa = abs(ft)
    ht = h
    ha = abs(h)

! pmax points to the maximum absolute element of matrix.
! pmax = 1 if f has largest absolute value.
! pmax = 2 if g has largest absolute value.
! pmax = 3 if h has largest absolute value.
    pmax = 1
    swap = (ha>fa)
    if (swap) then
     pmax = 3
     temp = ft
     ft = ht
     ht = temp
     temp = fa
     fa = ha
     ha = temp
    endif

! Now fa>=ha.
    gt = g
    ga = abs(gt)
    if (ga==0.0_KR2) then
! Diagonal matrix.
     ssmin = ha
     ssmax = fa
     clt = 1.0_KR2
     crt = 1.0_KR2
     slt = 0.0_KR2
     srt = 0.0_KR2
    else
     gasmal = .true.
     if (ga>fa) then
      pmax = 2
! Case of very large ga.
      if (fa/ga<dlamche) then
       gasmal = .false.
       ssmax = ga
       if (ha>1.0_KR2) then
        ssmin = fa/(ga/ha)
       else
        ssmin = (fa/ga)*ha
       endif
       clt = 1.0_KR2
       slt = ht/gt
       srt = 1.0_KR2
       crt = ft/gt
      endif
     endif
     if (gasmal) then
! Normal case.
      d = fa - ha
      if (d==fa) then
! Copes with infinite f or h.
       l = 1.0_KR2
      else
       l = d/fa
      endif
! Note that 0 <= l <= 1.
      m = gt/ft
! Note that abs(M) <= 1/macheps.
      t = 2.0_KR2 - l
! Note that t >= 1.
      mm = m**2
      tt = t**2
      s = sqrt(tt+mm)
! Note that 1 <= S <= 1 + 1/macheps.
      if (l==0.0_KR2) then
       r = abs(m)
      else
       r = sqrt(l**2+mm)
      endif
! Note that 0 <= R <= 1 + 1/macheps.
      a = 0.5_KR2*(s+r)
! Note that 1 <= a <= 1 + abs(m).
      ssmin = ha/a
      ssmax = fa*a
      if (mm==0.0_KR2) then
! Note that m is very tiny.
       if (l==0.0_KR2) then
        t = sign(2.0_KR2,ft)*sign(1.0_KR2,gt)
       else
        t = gt/sign(d,ft) + m/t
       endif
      else
       t = (m/(s+t)+m/(r+l))*(1.0_KR2+a)
      endif
      l = sqrt(t**2+4.0_KR2)
      crt = 2.0_KR2/l
      srt = t/l
      clt = (crt+srt*m)/a
      slt = (ht/ft)*srt/a
     endif
    endif

    if (swap) then
     csl = srt
     snl = crt
     csr = slt
     snr = clt
    else
     csl = clt
     snl = slt
     csr = crt
     snr = srt
    endif

! Correct signs of ssmax and ssmin.
    if (pmax==1) tsign=sign(1.0_KR2,csr)*sign(1.0_KR2,csl)*sign(1.0_KR2,f)
    if (pmax==2) tsign=sign(1.0_KR2,snr)*sign(1.0_KR2,csl)*sign(1.0_KR2,g)
    if (pmax==3) tsign=sign(1.0_KR2,snr)*sign(1.0_KR2,snl)*sign(1.0_KR2,h)
    ssmax = sign(ssmax,tsign)
    ssmin = sign(ssmin,tsign*sign(1.0_KR2,f)*sign(1.0_KR2,h))

 end subroutine svd2by2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine sing2by2(f,g,h,ssmin,ssmax)
! Compute the singular values of the 2 by 2 triangular matrix,  / f  g \
!                                                               \ 0  h /.
! This subroutine corresponds to SUBROUTINE DLAS2 of LAPACK.
! INPUT:
!   f is the (1,1) element of the 2 by 2 matrix.
!   g is the (1,2) element of the 2 by 2 matrix.
!   h is the (2,2) element of the 2 by 2 matrix.
! OUTPUT:
!   ssmin is the smaller singular value.
!   ssmax is the larger singular value.

    real(kind=KR2), intent(in)  :: f, g, h
    real(kind=KR2), intent(out) :: ssmin, ssmax

    real(kind=KR2) :: as, at, au, c, fa, fhmn, fhmx, ga, ha

    fa = abs(f)
    ga = abs(g)
    ha = abs(h)
    fhmn = min(fa,ha)
    fhmx = max(fa,ha)
    if (fhmn==0.0_KR2) then
     ssmin = 0.0_KR2
     if (fhmx==0.0_KR2) then
      ssmax = ga
     else
      ssmax = max(fhmx,ga)*sqrt(1.0_KR2+(min(fhmx,ga)/max(fhmx,ga))**2)
     endif
    else
     if (ga<fhmx) then
      as = 1.0_KR2 + fhmn/fhmx
      at = (fhmx-fhmn)/fhmx
      au = (ga/fhmx)**2
      c = 2.0_KR2/(sqrt(as**2+au)+sqrt(at**2+au))
      ssmin = fhmn*c
      ssmax = fhmx/c
     else
      au = fhmx/ga
      if (au==0.0_KR2) then
! Avoid possible harmful underflow if exponent range asymmetric (true ssmin
! may not underflow even if au underflows).
       ssmin = (fhmn*fhmx)/ga
       ssmax = ga
      else
       as = 1.0_KR2 + fhmn/fhmx
       at = (fhmx-fhmn)/fhmx
       c = 1.0_KR2/(sqrt(1.0_KR2+(as*au)**2)+sqrt(1.0_KR2+(at*au)**2))
       ssmin = (fhmn*c)*au
       ssmin = 2*ssmin
       ssmax = ga/(2.0_KR2*c)
      endif
     endif
    endif

 end subroutine sing2by2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine rotate2(f,g,cs,sn,r)
! Generate a plane rotation so that
!     / cs  sn \  / f \  =  / r \   where cs^2 + sn^2 = 1.
!     \-sn  cs /  \ g /     \ 0 /
! If f exceeds g in magnitude, cs will be positive.
! This subroutine corresponds to SUBROUTINE DLARTG of LAPACK.
! INPUT:
!   f is the first component of the vector to be rotated.
!   g is the second component of the vector to be rotated.
! OUTPUT:
!   cs is the cosine of the rotation.
!   sn is the sine of the rotation.
!   r is the nonzero component of the rotated vector.

    real(kind=KR2), intent(in)  :: f, g
    real(kind=KR2), intent(out) :: cs, sn, r

    integer(kind=KI) :: i, icount
    real(kind=KR2)   :: safmin, eps, safmn2, safmx2, f1, g1, myscale

    safmin = dlamchs
    eps = dlamche
    safmn2 = dlamch1
    safmx2 = 1.0_KR2/safmn2
    if (g==0.0_KR2) then
     cs = 1.0_KR2
     sn = 0.0_KR2
     r = f
    elseif (f==0.0_KR2) then
     cs = 0.0_KR2
     sn = 1.0_KR2
     r = g
    else
     f1 = f
     g1 = g
     myscale = max(abs(f1),abs(g1))
     if (myscale>=safmx2) then
      icount = 0
      firstloop: do
       icount = icount + 1
       f1 = f1*safmn2
       g1 = g1*safmn2
       myscale = max(abs(f1),abs(g1))
       if (myscale<safmx2) exit firstloop
      enddo firstloop
      r = sqrt(f1**2+g1**2)
      cs = f1/r
      sn = g1/r
      do i = 1,icount
       r = r*safmx2
      enddo ! i
     elseif (myscale<=safmn2) then
      icount = 0
      secondloop: do
       icount = icount + 1
       f1 = f1*safmx2
       g1 = g1*safmx2
       myscale = max(abs(f1),abs(g1))
       if (myscale>safmn2) exit secondloop
      enddo secondloop
      r = sqrt(f1**2+g1**2)
      cs = f1/r
      sn = g1/r
      do i = 1,icount
       r = r*safmn2
      enddo ! i
     else
      r = sqrt(f1**2+g1**2)
      cs = f1/r
      sn = g1/r
     endif
     if (abs(f)>abs(g) .and. cs<0.0_KR2) then
      cs = -cs
      sn = -sn
      r = -r
     endif
    endif

 end subroutine rotate2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine rotatecomp(f,g,cs,sn,r)
! Generate a plane rotation so that
!     / cs  sn \  / f \  =  / r \   where cs^2 + |sn|^2 = 1.
!     \-sn  cs /  \ g /     \ 0 /
! If g=0, then cs=1 and sn=0.
! If f=0, then cs=0 and sn is chosen so that r is real.
! This subroutine corresponds to SUBROUTINE ZLARTG of LAPACK.
! INPUT:
!   f is the first component of the vector to be rotated.
!   g is the second component of the vector to be rotated.
! OUTPUT:
!   cs is the cosine of the rotation.
!   sn is the sine of the rotation.
!   r is the nonzero component of the rotated vector.

    real(kind=KR2), intent(in),  dimension(:) :: f, g
    real(kind=KR2), intent(out)               :: cs
    real(kind=KR2), intent(out), dimension(:) :: sn, r

    integer(kind=KI)             :: i, iflag, icount
    real(kind=KR2)               :: safmin, eps, safmn2, safmx2, myscale, &
                                    temp, f2, g2, d, f2s, g2s, dr, di
    real(kind=KR2), dimension(2) :: fs, gs, ff

    iflag = 0
    safmin = dlamchs
    eps = dlamche
    safmn2 = dlamch1
    safmx2 = 1.0_KR2/safmn2
    myscale = max(abs(f(1)),abs(f(2)),abs(g(1)),abs(g(2)))
    fs = f
    gs = g
    icount = 0
    if (myscale>=safmx2) then
     firstloop: do
      icount = icount + 1
      fs = fs*safmn2
      gs = gs*safmn2
      myscale = myscale*safmn2
      if (myscale<safmx2) exit firstloop
     enddo firstloop
    elseif (myscale<=safmn2) then
     if (g(1)==0.0_KR2 .and. g(2)==0.0_KR2) then
      cs = 1.0_KR2
      sn = 0.0_KR2
      r = f
      iflag = 1
     endif
     if (iflag==0) then
      secondloop: do
       icount = icount - 1
       fs = fs*safmx2
       gs = gs*safmx2
       myscale = myscale*safmx2
       if (myscale>safmn2) exit secondloop
      enddo secondloop
     endif
    endif
    if (iflag==0) then
     f2 = fs(1)**2 + fs(2)**2
     g2 = gs(1)**2 + gs(2)**2
     if (f2<=max(g2,1.0_KR2)*safmin) then
! This is a rare case: f is very small.
      if (f(1)==0.0_KR2 .and. f(2)==0.0_KR2) then
       cs = 0.0_KR2
       r = pythag2(g(1),g(2))
! Do complex/real division explicitly with two real divisions.
       d = pythag2(gs(1),gs(2))
       sn(1) = gs(1)/d
       sn(2) = -gs(2)/d
       iflag = 1
      endif
      if (iflag==0) then
       f2s = pythag2(fs(1),fs(2))
! g2 and g2s are accurate: g2 is at least safmin, and g2s is at least safmn2.
       g2s = sqrt(g2)
       cs = f2s/g2s
! Make sure abs(ff) = 1.
! Do complex/real division explicitly with 2 real divisions.
       if (max(abs(f(1)),abs(f(2)))>1.0_KR2) then
        d = pythag2(f(1),f(2))
        ff = f/d
       else
        dr = safmx2*f(1)
        di = safmx2*f(2)
        d = pythag2(dr,di)
        ff(1) = dr/d
        ff(2) = di/d
       endif
       sn(1) = (ff(1)*gs(1)+ff(2)*gs(2))/g2s
       sn(2) = (ff(2)*gs(1)-ff(1)*gs(2))/g2s
       r(1) = cs*f(1) + sn(1)*g(1) - sn(2)*g(2)
       r(2) = cs*f(2) + sn(1)*g(2) + sn(2)*g(1)
      endif
     else
! This is the most common case.
! Neither f2 nor f2/g2 are less than safmin.
! f2s cannot overflow and it is accurate.
      f2s = sqrt(1.0_KR2+g2/f2)
! Do the f2s(real)*fs(complex) multiply with two real multiplies.
      r = f2s*fs
      cs = 1.0_KR2/f2s
      d = f2 + g2
! Do complex/real division explicitly with two real divisions.
      sn = r/d
      temp = sn(1)*gs(1) + sn(2)*gs(2)
      sn(2) = sn(2)*gs(1) - sn(1)*gs(2)
      sn(1) = temp
      if (icount/=0) then
       if (icount>0) then
        do i = 1,icount
         r = r*safmx2
        enddo ! i
       else
        do i = 1,-icount
         r = r*safmn2
        enddo ! i
       endif
      endif
     endif
    endif

 end subroutine rotatecomp

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 function twonorm(x,n) result(xnorm)
! Compute the 2-norm, also called Euclidean norm, of a complex vector.
! This function corresponds to FUNCTION DZNRM2 of LAPACK.
! INPUT:
!   x(1,:),x(2,:) are the Re,Im parts of the complex vector.
!   n is the length of the complex vector.
! OUTPUT:
!   xnorm is the 2-norm of the complex vector.

    real(kind=KR2),   intent(in), dimension(:,:) :: x
    integer(kind=KI), intent(in)                 :: n
    real(kind=KR2)                               :: xnorm

    real(kind=KR2)   :: bigscale, ssq, temp
    integer(kind=KI) :: ix

    if (n<1) then
     xnorm = 0.0_KR2
    else
     bigscale = 0.0_KR2
     ssq = 1.0_KR2
     do ix = 1,n
      if (x(1,ix)/=0.0_KR2) then
       temp = abs(x(1,ix))
       if (bigscale<temp) then
        ssq = 1.0_KR2 + ssq*(bigscale/temp)**2
        bigscale = temp
       else
        ssq = ssq + (temp/bigscale)**2
       endif
      endif
      if (x(2,ix)/=0.0_KR2) then
       temp = abs(x(2,ix))
       if (bigscale<temp) then
        ssq = 1.0_KR2 + ssq*(bigscale/temp)**2
        bigscale = temp
       else
        ssq = ssq + (temp/bigscale)**2
       endif
      endif
     enddo ! ix
     xnorm = bigscale*sqrt(ssq)
    endif

 end function twonorm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine compdiv(a,b,p)
! Compute p(1)+i*p(2) = (a(1)+i*a(2))/(b(1)+i*b(2)), where i=sqrt(-1)
! while avoiding unnecessary floating overflow.
! This subroutine corresponds to SUBROUTINE ZLADIV (and thus DLADIV) of LAPACK.

    real(kind=KR2), intent(in),  dimension(:) :: a, b
    real(kind=KR2), intent(out), dimension(:) :: p

    real(kind=KR2) :: c, d

    if (abs(b(2))<abs(b(1))) then
     c = b(2)/b(1)
     d = b(1) + b(2)*c
     p(1) = (a(1)+a(2)*c)/d
     p(2) = (a(2)-a(1)*c)/d
    else
     c = b(1)/b(2)
     d = b(2) + b(1)*c
     p(1) = (a(2)+a(1)*c)/d
     p(2) = (-a(1)+a(2)*c)/d
    endif

 end subroutine compdiv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 function pythag2(x,y) result(r)
! Compute r=sqrt(x^2+y^2) while avoiding unnecessary floating overflow.
! This function corresponds to FUNCTION DLAPY2 of LAPACK.

    real(kind=KR2), intent(in) :: x, y
    real(kind=KR2)             :: r

    real(kind=KR2) :: xabs, yabs, w, z

    xabs = abs(x)
    yabs = abs(y)
    w = max(xabs,yabs)
    z = min(xabs,yabs)
    if (z==0.0_KR2) then
     r = w
    else
     r = w*sqrt( 1.0_KR2 + (z/w)**2 )
    endif

 end function pythag2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 function pythag3(x,y,z) result(r)
! Compute r=sqrt(x^2+y^2+z^2) while avoiding unnecessary floating overflow.
! This function corresponds to FUNCTION DLAPY3 of LAPACK.

    real(kind=KR2), intent(in) :: x, y, z
    real(kind=KR2)             :: r

    real(kind=KR2) :: xabs, yabs, zabs, w

    xabs = abs(x)
    yabs = abs(y)
    zabs = abs(z)
    w = max(xabs,yabs,zabs)
    if (w==0.0_KR2) then
     r = 0.0_KR2
    else
     r = w*sqrt( (xabs/w)**2 + (yabs/w)**2 + (zabs/w)**2 )
    endif

 end function pythag3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine safesqrt(x,y)
! Take the square root of a complex number without arithmetic exception.
! This subroutine is not from LAPACK; they use complex fortran variables and
! thus avoid the issue.
! INPUT:
!   x(1),x(2) are the Re,Im parts of the original number.
! OUTPUT:
!   y(1),y(2) are the Re,Im parts of the square root of the original number.
! Which root is chosen: y(1) is not negative.
!                       If y(1)=0, then y(2) is not negative.

    real(kind=KR2), intent(in),  dimension(:) :: x
    real(kind=KR2), intent(out), dimension(:) :: y

    y(1) = sqrt(0.5_KR2*(sqrt(x(1)**2+x(2)**2)+x(1)))
    if (y(1)/=0.0_KR2) then
     y(2) = 0.5_KR2*x(2)/y(1)
    else
     y(2) = sqrt(0.5_KR2*(sqrt(x(1)**2+x(2)**2)-x(1)))
     if (y(2)/=0.0_KR2) then
      y(1) = 0.5_KR2*x(2)/y(2)
     else
      y(:) = 0.0_KR2
     endif
    endif

 end subroutine safesqrt

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine indexx(n,arr,indx)

!Indexes an array arr, i.e. outputs the array indx of length n such that 
!arr(indx(j)) is in ascending order for j=1,2,3,..,n. The input quantity
!arr is not changed.

real(kind=KR2),dimension(:),intent(in)::arr
integer(kind=KI),dimension(:),intent(out)::indx
integer(kind=KI),parameter::nn=15,nstack=50

real(kind=KR2)::a
integer(kind=KI)::n,k,i,j,indext,jstack,l,r,itemp
integer(kind=KI),dimension(nstack)::istack

do i=1,n
 indx(i)=i
enddo !i

jstack=0
l=1
r=n

do
  if(r-l < nn) then
    do j=l+1,r
     indext=indx(j)
     a=arr(indext)
     do i=j-1,l,-1
        if(arr(indx(i)) <= a) exit
        indx(i+1)=indx(i)
     enddo !i
     indx(i+1)=indext
    enddo !j

    if(jstack==0) return
    r=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+r)/2
    
    itemp=indx(k)
    indx(k)=indx(l+1)
    indx(l+1)=itemp

    if(arr(indx(l)).GT.arr(indx(r))) then
      itemp=indx(l)
      indx(l)=indx(r)
      indx(r)=itemp
    endif

    if(arr(indx(l+1)).GT.arr(indx(r))) then
      itemp=indx(l+1)
      indx(l+1)=indx(r)
      indx(r)=itemp
    endif

    if(arr(indx(l)).GT.arr(indx(l+1))) then
      itemp=indx(l)
      indx(l)=indx(l+1)
      indx(l+1)=itemp
    endif

    i=l+1
    j=r
    indext=indx(l+1)
    a=arr(indext)

    do 

      do
         i=i+1
         if(arr(indx(i)) .GE. a) exit
      enddo

      do
         j=j-1
         if(arr(indx(j)).LE. a) exit
      enddo

      if(j.LT.i) exit
      
      itemp=indx(i)
      indx(i)=indx(j)
      indx(j)=itemp

    enddo

    indx(l+1)=indx(j)
    indx(j)=indext
    jstack=jstack+2

    if(jstack > nstack) then
      print *, "Error in indexx: nstack too small"
      stop
    endif

    if( (r-i+1) .GE. (j-l) ) then
      istack(jstack) = r
      istack(jstack-1)=i
      r=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    endif
   endif
  enddo

end subroutine indexx     





! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end module pseudolapack

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
