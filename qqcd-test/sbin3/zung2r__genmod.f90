        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZUNG2R__genmod
          INTERFACE 
            SUBROUTINE ZUNG2R(M,N,K,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: TAU(*)
              COMPLEX(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZUNG2R
          END INTERFACE 
        END MODULE ZUNG2R__genmod
