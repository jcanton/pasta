MODULE  sparse_matrix_profiles 


   IMPLICIT NONE


   TYPE  CSR_MUMPS_Matrix
   
       INTEGER,      DIMENSION(:), POINTER :: i
       INTEGER,      DIMENSION(:), POINTER :: i_mumps    
       INTEGER,      DIMENSION(:), POINTER :: j
       REAL(KIND=8), DIMENSION(:), POINTER :: e ! matrix elements

   END TYPE  CSR_MUMPS_Matrix


   TYPE  CSR_MUMPS_Complex_Matrix
   
       INTEGER,         DIMENSION(:), POINTER :: i
       INTEGER,         DIMENSION(:), POINTER :: i_mumps    
       INTEGER,         DIMENSION(:), POINTER :: j
       COMPLEX(KIND=8), DIMENSION(:), POINTER :: e ! matrix elements

   END TYPE  CSR_MUMPS_Complex_Matrix

 
   TYPE  MSR_Matrix

       INTEGER,      DIMENSION(:), POINTER :: i ! Different from CSR !   
       INTEGER,      DIMENSION(:), POINTER :: j
       REAL(KIND=8), DIMENSION(:), POINTER :: e ! matrix elements

   END TYPE  MSR_Matrix


END MODULE  sparse_matrix_profiles
