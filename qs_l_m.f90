MODULE qs_L_M


! AXISYMMETRIC

   
   USE sparse_matrix_profiles
   
   USE prep_mesh_p1p2_sp
   
  
   CONTAINS


SUBROUTINE Dirichlet_L_M (js_D, AA,  symmetric)
!==============================================

!  This entry point modifies the elements of
!  the sparse matrix AA (in CSR format) to 
!  impose Dirichlet boundary conditions
   
   IMPLICIT NONE

   INTEGER, DIMENSION(:),  INTENT(IN)    :: js_D
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
  
   REAL (KIND=8), PARAMETER :: irons = 1.0d30 
   INTEGER :: n, i, p
   
   IF (PRESENT(symmetric)) THEN
   
      DO n = 1, SIZE(js_D);  i = js_D(n)
    
         DO p = AA%i(i), AA%i(i+1) - 1
        
            IF (AA%j(p) == i) THEN 
       
               AA%e(p) = irons
            
               EXIT
       
            ENDIF
       
         ENDDO
   
      ENDDO
  
   ELSE
    
      ! The possible symmetric character of matrix AA%e is destroyed  
    
      DO n = 1, SIZE(js_D);  i = js_D(n)
    
         DO p = AA%i(i), AA%i(i+1) - 1
        
            IF (AA%j(p) == i) THEN
               AA%e(p) = 1
            ELSE
               AA%e(p) = 0
            ENDIF
     
         ENDDO
   
      ENDDO
   
   ENDIF
  
END SUBROUTINE Dirichlet_L_M


!------------------------------------------------------------------------------


SUBROUTINE qs_0y0_L_M (alpha, AA,  symmetric)
!============================================

!  alpha < w, y _ >   ===>   AA%e    

   USE Gauss_points
   USE Gauss_points_L

   IMPLICIT NONE
   
   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
  
   INTEGER :: m, l, ni, nj, i, j, p
   REAL(KIND=8) :: al, x

   AA%e = 0

   DO m = 1, me

      DO l = 1, l_G_L

         al = alpha * jac_py_L(l,m)

         DO ni = 1, n_w_L;  i = jj(ni, m)

            DO nj = 1, n_w_L;  j = jj(nj, m)
               
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
             
               x = ww_L(ni,l) * al * ww_L(nj,l)
                 
               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_0y0_L_M


!------------------------------------------------------------------------------


SUBROUTINE qs_1y1_L_M (alpha, AA,  symmetric)
!===========================================

!  alpha << (Dw), y (D_) >>   ===>   AA%e    

   USE Gauss_points
   USE Gauss_points_L

   IMPLICIT NONE
  
   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
 
   INTEGER :: m, l, ni, nj, i, j, p
   REAL(KIND=8) :: al, x
   
   AA%e = 0
   
   DO m = 1, me

      DO l = 1, l_G_L

         al = alpha * jac_py_L(l,m)

         DO ni = 1, n_w_L;  i = jj(ni, m)

            DO nj = 1, n_w_L;  j = jj(nj, m)
              
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
    
               x = al * SUM(dw_L(:,ni,l,m) * dw_L(:,nj,l,m))
              
               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_1y1_L_M


END MODULE qs_L_M
