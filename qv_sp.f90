MODULE  qv_sp


! AXISYMMETRIC


CONTAINS


!------------------------------------------------------------------------------


SUBROUTINE qv_0y0_sp (m0, jj, gg, alpha,  v0)
!====================================

!  alpha << w, y g >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:)              :: v0

   REAL(KIND=8) :: gkl
   INTEGER      :: mm, m, l, k


   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, 3

            gkl = SUM(gg(k, jj(:,m)) * ww(:,l)) * JAC(m) * pp_w(l) * yy_G(l,m)

            v0(k, jj(:,m)) = v0(k, jj(:,m)) + alpha * ww(:,l) * gkl

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_0y0_sp

!------------------------------------------------------------------------------

SUBROUTINE qv_0y0_dR_sp (m0, jj, gg, alpha,  v0)
!====================================

!  alpha << w, g >>   ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:)              :: v0

   REAL(KIND=8) :: gkl
   INTEGER      :: mm, m, l, k


   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, 3

            gkl = SUM(gg(k, jj(:,m)) * ww(:,l)) * JAC(m) * pp_w(l)

            v0(k, jj(:,m)) = v0(k, jj(:,m)) + alpha * ww(:,l) * gkl

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_0y0_dR_sp

!------------------------------------------------------------------------------

SUBROUTINE qv_01_sp (m0, jj, ff,  v0)
!====================================

!  << w, Df >>    ===>    v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(k_d)      :: dfl, yl_p
   REAL(KIND=8), DIMENSION(n_w)      :: ffm
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   INTEGER :: mm, m, l, k, n

!   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
     
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            dfl(k) = SUM(ffm * dwl(k,:))
         ENDDO

         yl_p  =  dfl * pp_w(l)
         
         DO k = 1, k_d
            v0(k, jjm)  =  v0(k, jjm)  +  ww(:,l) * yl_p(k)
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_01_sp

!------------------------------------------------------------------------------

SUBROUTINE qv_y_10_hybrid_sp (m0, jj, jj_L, ff_L,  v0)
!===================================================

!  << D.w, f_L >>    ===>    v0   VECCHIO

! CANTON 

!  << y w, Df_L >>    integrando per parti

!  << y D.w  +  2 w_R , f_L >>    ===>    v0

! CANTON

   USE Gauss_points
  
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj, jj_L
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff_L
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(n_w_L)    :: ffm_L   
   INTEGER,      DIMENSION(n_w_L)    :: jjm_L
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L
  
   REAL(KIND=8) :: fl  
   
   INTEGER :: mm, m, l, k, n


   SELECT CASE (k_d)

      CASE (2)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(5,:) + ww(6,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:) + ww(4,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(4,:) + ww(5,:))

      CASE (3)
               
         w_L(1,:) = ww(1,:) + 0.5*(ww(n_w-2,:) + ww(n_w-1,:) + ww(n_w,:))
         w_L(2,:) = ww(2,:) + 0.5*(ww(6,:)     + ww(n_w-3,:) + ww(n_w,:))
         w_L(3,:) = ww(3,:) + 0.5*(ww(5,:)     + ww(n_w-3,:) + ww(n_w-1,:))
         w_L(4,:) = ww(4,:) + 0.5*(ww(5,:)     + ww(6,:)     + ww(n_w-2,:))
 
   END SELECT

!   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)
      
      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      ffm_L = ff_L(jjm_L)
  
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
         
         fl = SUM(ffm_L * w_L(:,l)) * pp_w(l) 
         
         v0(:, jjm) = v0(:, jjm)  +  dwl * fl * yy_G(l,m)  ! new * yy_G(l,m)
       
         v0(2, jjm) = v0(2, jjm)  +  ww(:,l) * fl          ! new  tutta la linea
       
      ENDDO

   ENDDO

END SUBROUTINE qv_y_10_hybrid_sp

!------------------------------------------------------------------------------


SUBROUTINE qv_10_sp (m0, jj, ff,  v0)
!====================================

!  << D.w, f >>    ===>    v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(n_w)      :: ffm   
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: fl  
   
   INTEGER :: mm, m, l, k, n

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
  
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
         
         fl = SUM(ffm * ww(:,l)) * pp_w(l)
       
         v0(:, jjm) = v0(:, jjm)  +  dwl * fl
     
      ENDDO

   ENDDO

END SUBROUTINE qv_10_sp


!------------------------------------------------------------------------------


SUBROUTINE qv_mass_grad_sp (m0, jj, gg, ff,  v0)
!===============================================

!     << w, g >>
!  +  << w, Df >>    ===>    v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: gl, dfl, yl_p
   REAL(KIND=8), DIMENSION(n_w)      :: ffm
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   INTEGER :: mm, m, l, k, n

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            dfl(k) = SUM(ffm * dwl(k,:))
         ENDDO

         DO k = 1, k_d
            yl_p(k)  =  (gl(k) * JAC(m)  +  dfl(k)) * pp_w(l)
         ENDDO

         DO k = 1, k_d
            v0(k, jjm)  =  v0(k, jjm)  +  ww(:,l) * yl_p(k)
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_mass_grad_sp


!------------------------------------------------------------------------------


SUBROUTINE qv_mass_grad_Madv_sp (m0, jj, gg, ff, vv,  v0)
!========================================================

!     << w, g >>
!  +  << w, Df >>  
!  -  << w, (v.D)v >>    ===>   v0
!
!  MINUS the QUADRATIC term
!  

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, vv
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, vvm, dwl
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dvl
   REAL(KIND=8), DIMENSION(k_d)      ::  gl, vl, dfl, yl_p
   REAL(KIND=8), DIMENSION(n_w)      :: ffm
   INTEGER,      DIMENSION(n_w)      :: jjm
 
   INTEGER :: mm, m, l, k, k1, n

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
      ggm = gg(:, jjm)
      vvm = vv(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            vl(k) = SUM(vvm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            dfl(k) = SUM(ffm * dwl(k,:))
            DO k1 = 1, k_d
               dvl(k,k1) = SUM(vvm(k,:) * dwl(k1,:))
            ENDDO
         ENDDO

         DO k = 1, k_d
            yl_p(k)  =  ( gl(k) * JAC(m)  +  dfl(k)   &
                           -  SUM(vl * dvl(k,:)) ) * pp_w(l)
         ENDDO

         DO k = 1, k_d
            v0(k, jjm)  =  v0(k, jjm)  +  ww(:,l) * yl_p(k)
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_mass_grad_Madv_sp


!------------------------------------------------------------------------------


SUBROUTINE qv_mass_Gibp_sp (m0, jj, gg, ff,  v0)
!===============================================

!                  Gradient integrated by parts

!     << w, g >>
!  +  << D.w, f >>    ===>    v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: gl, yl_p
   REAL(KIND=8), DIMENSION(n_w)      :: ffm   
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: fl  
   
   INTEGER :: mm, m, l, k, n

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
         
         fl = SUM(ffm * ww(:,l)) * pp_w(l)
         
         DO k = 1, k_d
            yl_p(k)  =  gl(k) * JAC(m) * pp_w(l)
         ENDDO


         DO k = 1, k_d
            v0(k, jjm)  =  v0(k, jjm)  +  ww(:,l) * yl_p(k)  &
                                       + dwl(k,:) * fl
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_mass_Gibp_sp


!------------------------------------------------------------------------------


SUBROUTINE qv_0y01_sp (m0, jj, gg,  v0)
!======================================
!
!  << w, y (g.D)g >>    ===>   v0

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(3, n_w)   :: ggm
   REAL(KIND=8), DIMENSION(3, k_d)   :: dgl
   REAL(KIND=8), DIMENSION(3)        :: gl  
   INTEGER,      DIMENSION(n_w)      :: jjm
 
   REAL(KIND=8) :: yl_p
   INTEGER :: mm, m, l, k3, k, n

!   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k3 = 1, 3
            gl(k3) = SUM(ggm(k3,:) * ww(:,l)) 
         ENDDO
         
         DO k = 1, k_d          
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k3 = 1, 3
            DO k = 1, k_d
               dgl(k3, k) = SUM(ggm(k3,:) * dwl(k,:))
            ENDDO
         ENDDO
       
         DO k3 = 1, 3
         
            yl_p = SUM(gl(1:2) * dgl(k3,:)) * pp_w(l) * yy_G(l,m)
         
            v0(k3, jjm)  =  v0(k3, jjm)  +  ww(:,l) * yl_p  
      
            SELECT CASE (k3) 
      
               CASE(2);  v0(2, jjm) = v0(2, jjm)  -  ww(:,l) * gl(3)**2 * JAC(m) * pp_w(l)
            
               CASE(3);  v0(3, jjm) = v0(3, jjm)  +  ww(:,l) * gl(3) * gl(2) * JAC(m) * pp_w(l)
      
            END SELECT
      
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_0y01_sp


!------------------------------------------------------------------------------


SUBROUTINE qv_mass_Gibp_Madv_sp (m0, jj, gg, ff, vv,  v0)
!========================================================

!                  Gradient integrated by parts

!     << w, g >>
!  +  << D.w, f >>  
!  -  << w, (v.D)v >>    ===>   v0
!
!  MINUS the QUADRATIC term
!  

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg, vv
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, vvm, dwl
   REAL(KIND=8), DIMENSION(k_d, k_d) :: dvl
   REAL(KIND=8), DIMENSION(k_d)      :: gl, vl, yl_p
   REAL(KIND=8), DIMENSION(n_w)      :: ffm
   INTEGER,      DIMENSION(n_w)      :: jjm
 
   REAL(KIND=8) :: fl
 
   INTEGER :: mm, m, l, k, k1, n

   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ffm = ff(jjm)
      ggm = gg(:, jjm)
      vvm = vv(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            gl(k) = SUM(ggm(k,:) * ww(:,l))
            vl(k) = SUM(vvm(k,:) * ww(:,l))
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 1, k_d
            DO k1 = 1, k_d
               dvl(k,k1) = SUM(vvm(k,:) * dwl(k1,:))
            ENDDO
         ENDDO
         
         fl = SUM(ffm * ww(:,l)) * pp_w(l)
        
         DO k = 1, k_d
            yl_p(k)  =  ( gl(k) * JAC(m)    &
                           -  SUM(vl * dvl(k,:)) ) * pp_w(l)
         ENDDO

         DO k = 1, k_d
            v0(k, jjm)  =  v0(k, jjm)  +  ww(:,l) * yl_p(k)  &
                                       + dwl(k,:) * fl
         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qv_mass_Gibp_Madv_sp


!------------------------------------------------------------------------------



END MODULE qv_sp
