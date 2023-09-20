! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! kinds.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! For f90, 32-bit is KR =4 and KI =4.
! For f90, 64-bit is KR2=8 and KI2=8.
! For F,   32-bit is KR =1 and KI =3.
! For F,   64-bit is KR2=2 and KI2=4.
! NOTE: module lagfib requires 32-bit integers: KI=4 for f90; KI=3 for F.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 module kinds
    implicit none
    integer, public, parameter :: KI=4, KR=8, KR2=8, KCC=8! !changed KR and KR2 from 8 to 4 for testing purposes TW 3/26/17                      
 end module kinds
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
