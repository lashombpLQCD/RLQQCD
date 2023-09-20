        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZLARFT__genmod
          INTERFACE 
            SUBROUTINE ZLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: DIRECT
              CHARACTER(LEN=1) :: STOREV
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: V(LDV,*)
              COMPLEX(KIND=8) :: TAU(*)
              COMPLEX(KIND=8) :: T(LDT,*)
            END SUBROUTINE ZLARFT
          END INTERFACE 
        END MODULE ZLARFT__genmod
