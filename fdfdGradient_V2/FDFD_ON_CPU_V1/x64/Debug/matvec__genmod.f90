        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 14 07:20:40 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATVEC__genmod
          INTERFACE 
            SUBROUTINE MATVEC(N,X,Y)
              USE GLOBAL
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: X(NX,NY,NZ,3)
              COMPLEX(KIND=8) :: Y(NX,NY,NZ,3)
            END SUBROUTINE MATVEC
          END INTERFACE 
        END MODULE MATVEC__genmod
