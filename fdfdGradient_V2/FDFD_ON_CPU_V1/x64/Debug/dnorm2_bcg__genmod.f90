        !COMPILER-GENERATED INTERFACE MODULE: Wed Jan 25 14:21:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DNORM2_BCG__genmod
          INTERFACE 
            FUNCTION DNORM2_BCG(N,ZX)
              INTEGER(KIND=4), INTENT(IN) :: N
              COMPLEX(KIND=8), INTENT(IN) :: ZX(N)
              REAL(KIND=8) :: DNORM2_BCG
            END FUNCTION DNORM2_BCG
          END INTERFACE 
        END MODULE DNORM2_BCG__genmod
