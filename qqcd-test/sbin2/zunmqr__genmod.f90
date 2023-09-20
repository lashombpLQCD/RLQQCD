        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZUNMQR__genmod
          INTERFACE 
            SUBROUTINE ZUNMQR(SIDE,TRANS,M,N,K,A,LDA,TAU,C,LDC,WORK,    &
     &LWORK,INFO)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: TAU(*)
              COMPLEX(KIND=8) :: C(LDC,*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZUNMQR
          END INTERFACE 
        END MODULE ZUNMQR__genmod
