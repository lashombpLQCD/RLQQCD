        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZLARFB__genmod
          INTERFACE 
            SUBROUTINE ZLARFB(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,LDT&
     &,C,LDC,WORK,LDWORK)
              INTEGER(KIND=4) :: LDWORK
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: SIDE
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIRECT
              CHARACTER(LEN=1) :: STOREV
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: V(LDV,*)
              COMPLEX(KIND=8) :: T(LDT,*)
              COMPLEX(KIND=8) :: C(LDC,*)
              COMPLEX(KIND=8) :: WORK(LDWORK,*)
            END SUBROUTINE ZLARFB
          END INTERFACE 
        END MODULE ZLARFB__genmod
