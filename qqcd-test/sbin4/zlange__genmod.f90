        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:44 2016
        MODULE ZLANGE__genmod
          INTERFACE 
            FUNCTION ZLANGE(NORM,M,N,A,LDA,WORK)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: NORM
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: WORK(*)
              REAL(KIND=8) :: ZLANGE
            END FUNCTION ZLANGE
          END INTERFACE 
        END MODULE ZLANGE__genmod
