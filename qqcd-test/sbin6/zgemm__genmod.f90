        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 15 14:42:43 2016
        MODULE ZGEMM__genmod
          INTERFACE 
            SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,&
     &C,LDC)
              INTEGER(KIND=4) :: LDC
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANSA
              CHARACTER(LEN=1) :: TRANSB
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=8) :: ALPHA
              COMPLEX(KIND=8) :: A(LDA,*)
              COMPLEX(KIND=8) :: B(LDB,*)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C(LDC,*)
            END SUBROUTINE ZGEMM
          END INTERFACE 
        END MODULE ZGEMM__genmod
