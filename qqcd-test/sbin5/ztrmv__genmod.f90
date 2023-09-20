        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZTRMV__genmod
          INTERFACE 
            SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE ZTRMV
          END INTERFACE 
        END MODULE ZTRMV__genmod
