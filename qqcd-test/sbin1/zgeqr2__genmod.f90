        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZGEQR2__genmod
          INTERFACE 
            SUBROUTINE ZGEQR2(M,N,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: TAU(*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZGEQR2
          END INTERFACE 
        END MODULE ZGEQR2__genmod
