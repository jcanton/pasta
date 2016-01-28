MODULE qs_sp_M


! AXISYMMETRIC


!  ww(ni,l)      : Parabolic shape function
!  Dw_re(k,n,l)  : Derivative w.r.t x_k of ww on the reference simplex
!  d...          : Derivative in the physical space
!  pp_w(l)       : Weight of Gauss points of the Parabolic approximation
!  dwl(k,n)      : dw(n)/dxk * Jac_det(m)   ---> SUM(MNR(k,:,m)*Dw_re(:,n,l))
!
!  Don't forget to multiply by Jacobian determinant for mass terms
!  Don't forget to divide by Jacobian determinant for stiffness terms
!
!  For stiffness matrix only we have:
!  M^TM_j(k,h,m) = SUM(MNR(k,:,m)*MNR(h,:,m))/JAC(m)
!
!   ===> dd_ij_ml = SUM_SUM(Dw_re(:,ni,l) * MTM_j(:,:,m) * Dw_re(:,nj,l))

  USE sparse_matrix_profiles

  USE prep_mesh_p1p2_sp  


  CONTAINS
 
 
SUBROUTINE Dirichlet_M (js_D, AA,  symmetric)
!============================================

!  Mdification of elements of matrix AA to
!  impose Dirichlet boundary conditions at nodes js_D

   USE Gauss_points

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
  
  
END SUBROUTINE Dirichlet_M


!------------------------------------------------------------------------------


SUBROUTINE qs_1_y_sp_M (alpha, AA,  symmetric)
!=============================================

! AXISYMMETRIC :  supplementary term to the negative Laplacian
!                 acting on a vector field, for the Poisson equation 
!                 of the nonaxial components of a vector unknown:   
! 
!                 + alpha 1/y^2 _
!                 
!                 which, after multiplication by y, becomes
!
!                 + alpha 1/y _
!                 
!                 and in weak form:
!
!                 + alpha < w, (1/y) _ >   ===>   AA
!               
! The value of  alpha  is  1  in axisymmetric problems
! and can be equal to  m^2  for the m-th Fourier mode in 
! nonaxisymmetric problems in axisymmetric domains 
!           
!  ===>   AA%e   

             
   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
 
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, alpha_m
   INTEGER      :: l, m, ni, nj, i, j, p


   AA%e = 0

   
   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
      ENDDO

   ENDDO


   DO m = 1, me

      alpha_m = alpha * JAC(m)

      DO l = 1, l_G

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)
                
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
               
               x  =  alpha_m * w_w_p(ni,nj,l) / yy_G(l,m) ! 1/y^2

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_1_y_sp_M


!------------------------------------------------------------------------------


SUBROUTINE qs_0y0_sp_M (alpha, AA,  symmetric)
!=============================================

!  alpha < w, y _ >    ===>   AA
!                          
!  ===>   AA%e         

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
  
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, alpha_m
   INTEGER      :: l, m, ni, nj, i, j, p


   AA%e = 0


   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
      ENDDO

   ENDDO


   DO m = 1, me

      alpha_m = alpha * JAC(m)

      DO l = 1, l_G

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)
               
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
             
               x  =  alpha_m * w_w_p(ni,nj,l) * yy_G(l,m)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_0y0_sp_M


!------------------------------------------------------------------------------


SUBROUTINE qs_00_sp_M (alpha, AA,  symmetric)
!============================================

!  alpha < w, _ >   ===>   AA    [no multiplication by y = R]
!                          
!  ===>   AA%e         

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   LOGICAL, OPTIONAL,      INTENT(IN)    :: symmetric
  
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, alpha_m
   INTEGER      :: l, m, ni, nj, i, j, p


   AA%e = 0


   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
      ENDDO

   ENDDO


   DO m = 1, me

      alpha_m = alpha * JAC(m)

      DO l = 1, l_G

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)
               
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
             
               x  =  alpha_m * w_w_p(ni,nj,l) ! * yy_G(l,m)  
                                              ! modifica rispetto
                                              ! a  qs_0y0_sp_M

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_00_sp_M


!------------------------------------------------------------------------------



SUBROUTINE qs_1y1_sp_M (alpha, AA,  PLUS,  symmetric)
!====================================================

!  +  alpha  [ << (Dw), y D_ >>  
!              +  IF (PLUS)  PLUS * < w, (1/y) _ > ]  ===>   AA
!                          
!  ===>   AA%e         

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),           INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   REAL(KIND=8), OPTIONAL, INTENT(IN)    :: PLUS
   LOGICAL,      OPTIONAL, INTENT(IN)    :: symmetric
  
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: dwl
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p
   
   REAL(KIND=8) :: x, alpha_p
   INTEGER      :: k, l, m, ni, nj, i, j, p
 
  
!   AA%e = 0


   IF (PRESENT(PLUS)) THEN  

      DO ni = 1, n_w

         DO nj = 1, n_w
      
            w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
         ENDDO

      ENDDO

   ENDIF
   


   DO m = 1, me

      DO l = 1, l_G

         DO k = 1, k_d
            DO ni = 1, n_w
               dwl(k,ni) = SUM(MNR(k,:,m) * Dw_re(:,ni,l))
            ENDDO
         ENDDO


         alpha_p = alpha * pp_w(l) * yy_G(l,m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)
 
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
             
               x  =  alpha_p * SUM(dwl(:,ni) * dwl(:,nj))/JAC(m)

               IF (PRESENT(PLUS)) &
                  x  =  x  +  alpha * PLUS * JAC(m) * w_w_p(ni,nj,l) / yy_G(l,m) ! 1/y^2

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_1y1_sp_M



!------------------------------------------------------------------------------


SUBROUTINE qs_y_diff_mass_sp_M (alpha, beta, AA,  PLUS,  symmetric)
!==================================================================

!  +  alpha  [ << (Dw), y D_ >>  
!              +  IF (PLUS)  PLUS * < w, (1/y) _ > ] 
!  +  beta  < w, y _ >    ===>   AA
!                          
!  ===>   AA%e         

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8),           INTENT(IN)    :: alpha, beta
   TYPE(CSR_MUMPS_Matrix), INTENT(INOUT) :: AA  
   REAL(KIND=8), OPTIONAL, INTENT(IN)    :: PLUS
   LOGICAL,      OPTIONAL, INTENT(IN)    :: symmetric
 
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: dwl
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p

   REAL(KIND=8) :: x, alpha_p,  beta_m,  dd_ij_ml
   INTEGER      :: k, l, m, ni, nj, i, j, p
  

   AA%e = 0


   DO ni = 1, n_w

      DO nj = 1, n_w
      
         w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
      ENDDO

   ENDDO
 

   DO m = 1, me

      beta_m = beta * JAC(m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO ni = 1, n_w
               dwl(k,ni) = SUM(MNR(k,:,m) * Dw_re(:,ni,l))
            ENDDO
         ENDDO


         alpha_p = alpha * pp_w(l) * yy_G(l,m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)
 
               IF (PRESENT(symmetric)  .AND.  i > j) CYCLE
             
               dd_ij_ml = SUM(dwl(:,ni) * dwl(:,nj))/JAC(m)

               x  =  alpha_p * dd_ij_ml  +  beta_m * w_w_p(ni,nj,l) * yy_G(l,m)
               
               IF (PRESENT(PLUS))  &
                  x  =  x  +  alpha * PLUS * JAC(m) * w_w_p(ni,nj,l) / yy_G(l,m) ! 1/y^2
               
               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_y_diff_mass_sp_M


!------------------------------------------------------------------------------


SUBROUTINE qs_y_skew_adv_sp_M (gg, AA)
!=====================================
!
!  < w, y (g.D)_ >  +  < w, [y(D.g) + g(2)] _ > / 2    ===>   AA
!
! AXISYMMETRIC :  Modification to include the additional term 
!                 of the expanded divergence in cylindrical 
!                 coordinates for axisymmetric problems
!                 
!                 + 1/2 g(2)/y _
!                 
!                 which, after multiplication by y, becomes
!
!                 + 1/2 g(2) _
!
!                 and in weak form:
!
!                 + < w, g(2) _ > / 2 
!                          
!                          
!  ===>   AA%e        

! WARNING:  This ENTRY POINT is meaningful and yields a valid 
!           matrix only provided that the vector field gg has the
!           second component equal to zero for y = 0 (i.e. on x axis) 
!           whenever the computational domain has a part of its 
!           boundary reaching the axis x. 

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)           :: gl
   REAL(KIND=8), DIMENSION(n_w)           :: gl_p
   INTEGER,      DIMENSION(n_w)           :: jjm

   INTEGER      :: k, l, m, n, ni, nj, i, j, p
   REAL(KIND=8) :: dgl, x


   AA%e = 0


   DO ni = 1, n_w
   
      DO nj = 1, n_w
         
         w_w_p(ni,nj,:) = ww(ni,:) * ww(nj,:) * pp_w / 2
               
      ENDDO
   
   ENDDO


   DO m = 1, me

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO ni = 1, n_w
            gl_p(ni) = SUM(gl * dwl(:,ni)) * pp_w(l) * yy_G(l,m) 
         ENDDO

         dgl = SUM(ggm * dwl) * yy_G(l,m) 

         DO ni = 1, n_w;  i = jjm(ni)

            DO nj = 1, n_w;  j = jjm(nj)

               x  =  ww(ni,l) * gl_p(nj)  & 
             
                  +  w_w_p(ni,nj,l) * (dgl + gl(2)*JAC(m))

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_y_skew_adv_sp_M



!------------------------------------------------------------------------------


SUBROUTINE qs_0y01_sp_M (gg, AA,  PLUS)
!======================================
!
!  < w, y (g.D)_ >  +  IF (PLUS)  PLUS * < w, g(2)_ >   ===>   AA
!
! AXISYMMETRIC :  Modification to include the additional term 
!                 of the advection equation for the angular
!                 component of velocity                      
!                          
!  ===>   AA%e        

! WARNING:  When  PLUS is true, this  ENTRY POINT
!           is meaningful and provides a valid matrix AA
!           only provided that the vector field gg has the
!           second component equal to zero for y = 0 (i.e. on x axis) 
!           whenever part of the boundaro of the computational 
!           domain lies on the axis x. 

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: AA  
   REAL(KIND=8), OPTIONAL,       INTENT(IN)    :: PLUS

   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)           :: gl
   REAL(KIND=8), DIMENSION(n_w)           :: gl_p
   INTEGER,      DIMENSION(n_w)           :: jjm

   INTEGER      :: k, l, m, ni, nj, i, j, p, n
   REAL(KIND=8) :: x
   
   
   AA%e = 0


   IF (PRESENT(PLUS)) THEN 

      DO ni = 1, n_w

         DO nj = 1, n_w
      
            w_w_p(ni,nj, :) = ww(ni,:) * ww(nj,:) * pp_w 
      
         ENDDO

      ENDDO

   ENDIF
   

   DO m = 1, me

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO n = 1, n_w
            gl_p(n) = SUM(gl * dwl(:,n)) * pp_w(l) * yy_G(l,m)
         ENDDO

         DO ni = 1, n_w;  i = jjm(ni)

            DO nj = 1, n_w;  j = jjm(nj)

               x  =  ww(ni,l) * gl_p(nj)  
          
               IF (PRESENT(PLUS))  &
                  x  =  x  +  PLUS * w_w_p(ni,nj,l) * gl(2) * JAC(m)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_0y01_sp_M



!------------------------------------------------------------------------------


SUBROUTINE qs_skew_adv_sp_M_crt (gg, AA)
!=======================================

!  < w, y (g.D)_ >  +  < w, y (D.g) _ > / 2    ===>   AA
!                          
!  ===>   AA%e        

   USE Gauss_points

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: AA  

   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: w_w_p
   REAL(KIND=8), DIMENSION(k_d, n_w)      :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)           :: gl
   REAL(KIND=8), DIMENSION(n_w)           :: gl_p
   INTEGER,      DIMENSION(n_w)           :: jjm

   INTEGER      :: k, l, m, ni, nj, i, j, p, n
   REAL(KIND=8) :: dgl, x

   
   AA%e = 0


   DO ni = 1, n_w
   
      DO nj = 1, n_w
         
         w_w_p(ni,nj,:) = ww(ni,:) * ww(nj,:) * pp_w / 2
               
      ENDDO
   
   ENDDO


   DO m = 1, me

      jjm = jj(:,m)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         dgl = SUM(ggm * dwl)

         DO n = 1, n_w
            gl_p(n) = SUM(gl * dwl(:,n)) * pp_w(l) * yy_G(l,m)
         ENDDO

         DO ni = 1, n_w;  i = jjm(ni)

            DO nj = 1, n_w;  j = jjm(nj)

               x  =  ww(ni,l) * gl_p(nj)  +  dgl * w_w_p(ni,nj,l) * yy_G(l,m)

               DO p = AA%i(i),  AA%i(i+1) - 1
                  IF (AA%j(p) == j) THEN;  AA%e(p) = AA%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qs_skew_adv_sp_M_crt



END MODULE qs_sp_M
