        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZLASCL__genmod
          INTERFACE 
            SUBROUTINE ZLASCL(TYPE,KL,KU,CFROM,CTO,M,N,A,LDA,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TYPE
              INTEGER(KIND=4) :: KL
              INTEGER(KIND=4) :: KU
              REAL(KIND=8) :: CFROM
              REAL(KIND=8) :: CTO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ZLASCL
          END INTERFACE 
        END MODULE ZLASCL__genmod
