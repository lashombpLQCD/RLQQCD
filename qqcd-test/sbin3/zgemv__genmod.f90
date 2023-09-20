        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZGEMV__genmod
          INTERFACE 
            SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE ZGEMV
          END INTERFACE 
        END MODULE ZGEMV__genmod
