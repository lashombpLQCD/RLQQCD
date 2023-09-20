        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZGERC__genmod
          INTERFACE 
            SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
              COMPLEX(KIND=8) :: A(LDA,*)
            END SUBROUTINE ZGERC
          END INTERFACE 
        END MODULE ZGERC__genmod
