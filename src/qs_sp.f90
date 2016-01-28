MODULE  qs_sp


! AXISYMMETRIC 


   USE prep_mesh_p1p2_sp


CONTAINS

SUBROUTINE Dirichlet (js_D, us_D,  uu,  symmetric)
!=================================================

   IMPLICIT NONE
   INTEGER,      DIMENSION(:), INTENT(IN)  :: js_D
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: us_D
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: uu
   LOGICAL,      OPTIONAL,     INTENT(IN)  :: symmetric
   
   REAL (KIND=8), PARAMETER :: irons = 1.0d30  

   IF (PRESENT(symmetric)) THEN
   
      IF (symmetric)  uu(js_D) = irons * us_D

   ELSE

      uu(js_D) = us_D
   
   ENDIF

END SUBROUTINE Dirichlet


!------------------------------------------------------------------------------


SUBROUTINE qs_0y0_sp (ff,  u0)
!=============================

!  < w, y f >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u0

   REAL(KIND=8) :: fl
   INTEGER      :: m, l

   u0 = 0

   DO m = 1, me

      DO l = 1, l_G

         fl = SUM(ff(jj(:,m)) * ww(:,l)) * JAC(m) * pp_w(l) * yy_G(l,m) 

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fl

      ENDDO

   ENDDO

END SUBROUTINE qs_0y0_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_000_sp (ff, hh,  u0)
!=================================

!  < w, f h >   ===>   u0


   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ff, hh
   REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u0

   REAL(KIND=8) :: fhl
   INTEGER      :: m, l

   u0 = 0

   DO m = 1, me

      DO l = 1, l_G

         fhl = SUM(ff(jj(:,m)) * ww(:,l))  & 
              * SUM(hh(jj(:,m)) * ww(:,l))  & 
              * JAC(m) * pp_w(l) 

         u0(jj(:,m)) = u0(jj(:,m)) + ww(:,l) * fhl

      ENDDO

   ENDDO

END SUBROUTINE qs_000_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_0y1_sp (gg,  u0)
!=============================

!  < w, y (D.g) + g(2) >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl, ggm
   INTEGER,      DIMENSION(n_w)      :: jjm
  
   INTEGER      :: m, l, n, k
   REAL(KIND=8) :: dgl_p

   u0 = 0

   DO m = 1, me

      jjm = jj(:,m)
      
      ggm = gg(:, jjm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
         
         dgl_p = (SUM(ggm * dwl) * yy_G(l,m)  &
                + SUM(ggm(2,:) * ww(:,l)) * JAC(m)) * pp_w(l)
         
         u0(jjm) = u0(jjm) + ww(:,l) * dgl_p
      
      ENDDO

   ENDDO


END SUBROUTINE qs_0y1_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_0y01_sp (gg, ff,  u0)
!==================================

!  +  << w, y g.Df >>    ===>    u0

   USE Gauss_points

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: gl, dfl
   REAL(KIND=8), DIMENSION(n_w)      :: ffm
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: yl_p
   INTEGER :: m, l, k, n

   u0 = 0

   DO m = 1, me

      jjm = jj(:,m)

      ggm = gg(:,jjm)
      ffm = ff(jjm)
     
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
                              ! errore:  Jacobiano da eliminare
         ! yl_p = SUM(gl * dfl) * JAC(m) * pp_w(l) * yy_G(l,m)
         
         yl_p = SUM(gl * dfl) * pp_w(l) * yy_G(l,m)
         
         u0(jjm) = u0(jjm) + ww(:,l) * yl_p

      ENDDO

   ENDDO

END SUBROUTINE qs_0y01_sp


!------------------------------------------------------------------------------


SUBROUTINE qs_y_mass_Madv_Mg2u_sp (ff, gg, uu,  u0)
!==================================================

!  < w,  y f  -  y (g.D)u  -  g(2)u >   ===>   u0

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, uu
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: ggm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: gl,  dul
   REAL(KIND=8), DIMENSION(n_w)      :: ffm, uum
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: fl,  g_dul,  g2_ul
   INTEGER      :: m, n, l, k

   u0 = 0

   DO m = 1, me

      jjm = jj(:,m)
     
      ffm = ff(jjm) 
      uum = uu(jjm)
      ggm = gg(:, jjm)
 
      DO l = 1, l_G

         fl = SUM(ffm * ww(:,l)) * JAC(m) * yy_G(l,m) 

         DO k = 1, k_d
            
            gl(k) = SUM(ggm(k,:) * ww(:,l))
             
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO

         ENDDO

         DO k = 1, k_d
            dul(k) = SUM(uum * dwl(k,:))
         ENDDO
                                            
         g_dul = SUM(gl * dul) * yy_G(l,m)  
        
         g2_ul = gl(2) * SUM(uum * ww(:,l)) * JAC(m) 

         u0(jjm) = u0(jjm) + ww(:,l) * (fl - g_dul - g2_ul) * pp_w(l)

      ENDDO

   ENDDO

END SUBROUTINE qs_y_mass_Madv_Mg2u_sp


!------------------------------------------------------------------------------

!------------------------------------------------------------------------------


SUBROUTINE qs_0y1_sp_c (gg,  u0)
!===============================

!  < w, y k.D x g >   ===>   u0       ( 2d only )

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl, ggm
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   INTEGER :: m, l, k, n
   REAL(KIND=8) :: c_gl_p


   IF (k_d /= 2) THEN
   
       WRITE (*,*) 'Program qs_0y1_sp_c is valid only in two dimensions' 
   
       STOP 
       
   ENDIF


   u0 = 0

   DO m = 1, me
    
      jjm = jj(:,m)
      
      ggm = gg(:, jjm)
    
      DO l = 1, l_G
    
         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
      
         c_gl_p = SUM(ggm(2,:) * dwl(1,:)   &
                    - ggm(1,:) * dwl(2,:)) * pp_w(l) * yy_G(l,m)

         u0(jjm) = u0(jjm) + ww(:,l) * c_gl_p

      ENDDO
 
   ENDDO


END SUBROUTINE qs_0y1_sp_c


!------------------------------------------------------------------------------


SUBROUTINE qs_01_sp_c (gg,  u0)
!==============================

!  < w, k.D x g >   ===>   u0    [no multiplication by y = R] (2d only)

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl, ggm
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   INTEGER :: m, l, k, n
   REAL(KIND=8) :: c_gl_p


   IF (k_d /= 2) THEN
   
       WRITE (*,*) 'Program qs_01_sp_c_cartesian is valid only in two dimensions' 
   
       STOP 
       
   ENDIF


   u0 = 0

   DO m = 1, me
    
      jjm = jj(:,m)
      
      ggm = gg(:, jjm)
    
      DO l = 1, l_G
    
         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
      
         c_gl_p = SUM(ggm(2,:) * dwl(1,:)  &
                    - ggm(1,:) * dwl(2,:)) * pp_w(l) ! * yy_G(l,m) 
                                                     ! modifica rispetto
                                                     ! a  qs_0y1_sp_c

         u0(jjm) = u0(jjm) + ww(:,l) * c_gl_p

      ENDDO
 
   ENDDO


END SUBROUTINE qs_01_sp_c


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


SUBROUTINE qs_0y0_sp_s (ms_N, fs,  u0)
!=====================================

!  < ws, y fs >_s   ===>   u0   incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:), INTENT(IN)    :: ms_N
   REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: fs
   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: u0

   INTEGER :: mm, ms, ls
   REAL(KIND=8) :: fls

   DO mm = 1, SIZE(ms_N);  ms = ms_N(mm)
     
      DO ls = 1, l_Gs

         fls = SUM(fs(iis(:,ms)) * wws(:,ls)) * JACs(ms) * pp_ws(ls) * yy_Gs(ls,ms)

         u0(jjs(:,ms)) = u0(jjs(:,ms)) + wws(:,ls) * fls

      ENDDO
   
   ENDDO


END SUBROUTINE qs_0y0_sp_s


!------------------------------------------------------------------------------


SUBROUTINE qs_0y1_sp_s (ms_N, gs,  u0)
!=====================================

!  < ws, y n.g_s >_s   ===>   u0   incremental accumulation of boundary terms
  
   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms_N
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:)                :: u0

   REAL(KIND=8), DIMENSION(k_d) :: gls
   REAL(KIND=8) :: x
   INTEGER :: mm, ms, ls, k
   
   DO mm = 1, SIZE(ms_N);  ms = ms_N(mm)

      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k, iis(:,ms)) * wws(:,ls))
         ENDDO

         x = SUM(gls * normals(:,ls,ms)) * JACs(ms) * pp_ws(ls) * yy_Gs(ls,ms)
         
         u0(jjs(:,ms)) = u0(jjs(:,ms)) + wws(:,ls) * x
         !u0(jjs(:,ms)) = wws(:,ls) * x
                
      ENDDO

   ENDDO


END SUBROUTINE qs_0y1_sp_s

!------------------------------------------------------------------------------


SUBROUTINE qs_0y1_sp_sl (ms_N, gs,  u0)
!=====================================

!  < ws, y n.g_s >_sl   ===>   u0   incremental accumulation of boundary terms
 
   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms_N
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gs
   REAL(KIND=8), DIMENSION(:)                :: u0

   REAL(KIND=8), DIMENSION(k_d) :: gls
   REAL(KIND=8) :: x, slval
   INTEGER :: mm, ms, ls, k

   slval = 0d0

   DO mm = 1, SIZE(ms_N);  ms = ms_N(mm)

      DO ls = 1, l_Gs

         DO k = 1, k_d
            gls(k) = SUM(gs(k, iis(:,ms)) * wws(:,ls))
         ENDDO

         x = SUM(gls * normals(:,ls,ms)) * JACs(ms) * pp_ws(ls) * yy_Gs(ls,ms)

         slval = slval + x

         u0(jjs(:,ms)) = u0(jjs(:,ms)) + wws(:,ls) * slval

      ENDDO

   ENDDO


END SUBROUTINE qs_0y1_sp_sl

!------------------------------------------------------------------------------


SUBROUTINE qs_y_mass_ROT21_sp (ff, vv, uu,  u0)
!==============================================

!  < w,  y f  +  y [v x Dxu]_3 >        ===>   u0
!
!  < w,  y f  +  y v.Du  +  v(2) u >    ===>   u0

   USE Gauss_points

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, uu
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: vv
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: u0
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: vvm, dwl
   REAL(KIND=8), DIMENSION(k_d)      :: vl,  dul
   REAL(KIND=8), DIMENSION(n_w)      :: ffm, uum
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: fl,  v_dul,  v2_ul
   INTEGER      :: m, n, l, k

   u0 = 0

   DO m = 1, me

      jjm = jj(:,m)
     
      ffm = ff(jjm) 
      uum = uu(jjm)
      vvm = vv(:, jjm)
 
      DO l = 1, l_G

         fl = SUM(ffm * ww(:,l)) * JAC(m) * yy_G(l,m) 

         DO k = 1, k_d
            
            vl(k) = SUM(vvm(k,:) * ww(:,l))
             
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO

         ENDDO

         DO k = 1, k_d
            dul(k) = SUM(uum * dwl(k,:))
         ENDDO
                                            
         v_dul = SUM(vl * dul) * yy_G(l,m)  
        
         v2_ul = vl(2) * SUM(uum * ww(:,l)) * JAC(m) 

         u0(jjm) = u0(jjm) + ww(:,l) * (fl + v_dul + v2_ul) * pp_w(l)

      ENDDO

   ENDDO

END SUBROUTINE qs_y_mass_ROT21_sp


END MODULE qs_sp
