        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZGEQRF__genmod
          INTERFACE 
            SUBROUTINE ZGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: TAU(*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZGEQRF
          END INTERFACE 
        END MODULE ZGEQRF__genmod
