MODULE qc_sp_M
!  COUPLED VELOCITY-PRESSURE SYSTEM OF EQUATIONS

!  THE USE OF THIS PROGRAM ASSUMES THAT THE LARGE SPARSE
!  MATRIX OF THE COUPLED SYSTEM IS STORED IN CSR FORMAT


!  ww(ni,l)      : Parabolic shape function
!  Dw_re(k,n,l)  : Derivative w.r.t x_k of ww on the reference simplex
!  d...          : Derivative in the physical space
!  pp_w(l)       : Weight of Gauss points of the Parabolic approximation
!  dwl(k,n)      : dw(n)/dxk * Jac_det(m)   ---> SUM(MNR(k,:,m)*Dw_re(:,n,l))
!
!  Don't forget to multiply by the Jacobian determinant for mass terms
!  Don't forget to divide by the Jacobian determinant for stiffness terms
!
!  For stiffness matrix only we have:
!  M^TM_j(k,h,m) = SUM(MNR(k,:,m)*MNR(h,:,m))/JAC(m)
!
!   ===> dd_ij_ml = SUM_SUM(Dw_re(:,ni,l) * MTM_j(:,:,m) * Dw_re(:,nj,l))


  USE sparse_matrix_profiles

  USE dynamic_structures


CONTAINS


SUBROUTINE ComputeJacobianMatrix(np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, CC,  Re, UU)
!===============================================================================

   IMPLICIT NONE
   
   !-----------------------------------------------------------------------!
   INTEGER,                            INTENT(IN) :: np 
   INTEGER,            DIMENSION(:),   INTENT(IN) :: mm
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj_L
   INTEGER,   POINTER, DIMENSION(:)               :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:),   INTENT(IN) :: js_D 
   LOGICAL,                            INTENT(IN) :: DESINGULARIZE
   
   TYPE(CSR_MUMPS_Matrix),          INTENT(INOUT) :: CC

   REAL(KIND=8),                           INTENT(IN) :: Re
   REAL(KIND=8), OPTIONAL, DIMENSION(:,:), INTENT(IN) :: UU
   !-----------------------------------------------------------------------!

   REAL(KIND=8), PARAMETER :: zero = 0,  one = 1
   INTEGER :: Nx, p
   

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (weak)  
   !             CC  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   Nx = SIZE(CC%i) - 1

   CC%e = 0

   ! CALL qc_1y1_sp_M (mm, jj,     1d0/Re,  CC) ! + stiffness (ROT-DIV)

   CALL qc_1y1_sp_gg_M (mm, jj,  1d0/Re,  CC) ! + stiffness (GRAD:GRAD)

   CALL qc_1y0_sp_M (mm, jj, jj_L, -one,  CC) ! + pressure gradient (ibp)  

   CALL qc_0y1_sp_M (mm, jj, jj_L, -one,  CC) ! - velocity divergence
      

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  [(U.G)_ + (_.G)U)]  + 1/Re *  K_  +  G_ (ibp)  
   !-------------ADD THE ITERATION DEPENDENDT PART------------------- 


   IF (PRESENT(UU)) THEN 
      CALL qc_oseen2y_sp_M (mm, jj, UU,  CC) ! linearized terms
   END IF
   
   CALL Dirichlet_c_M (np, js_Axis, js_D,  CC)

   IF (DESINGULARIZE) THEN
      ! reduction of the row of the last equation to 
      ! the diagonal element alone, set equal to 1   
      DO p = CC%i(Nx), CC%i(Nx + 1) - 1
         CC%e(p) = 0
         IF (CC%j(p) == Nx) CC%e(p) = 1
      ENDDO
   ENDIF
   
END SUBROUTINE ComputeJacobianMatrix
 

SUBROUTINE ComputeJacobianMatrix_OLD(np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, CC,  Re_U)
!==========================================================================================

   IMPLICIT NONE
   
   !-----------------------------------------------------------------------!
   INTEGER,                            INTENT(IN) :: np 
   INTEGER,            DIMENSION(:),   INTENT(IN) :: mm
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj_L
   INTEGER,   POINTER, DIMENSION(:)               :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:),   INTENT(IN) :: js_D 
   LOGICAL,                            INTENT(IN) :: DESINGULARIZE
   
   TYPE(CSR_MUMPS_Matrix),          INTENT(INOUT) :: CC
   
   REAL(KIND=8), OPTIONAL, DIMENSION(:,:), INTENT(IN) :: Re_U
   !-----------------------------------------------------------------------!

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: CCe_save
   REAL(KIND=8), PARAMETER :: zero = 0,  one = 1
   INTEGER :: Nx, p
   LOGICAL, SAVE :: initialized = .FALSE.
   

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (weak)  
   !             CC  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   Nx = SIZE(CC%i) - 1

   IF (.NOT. initialized) THEN

      CC%e = zero  !  remind the incremental accumulation in CC
      !==========================================================

      initialized = .TRUE.

      CALL qc_1y1_sp_M (mm, jj,        one,  CC) ! + stiffness (ROT-DIV)

      CALL qc_1y0_sp_M (mm, jj, jj_L, -one,  CC) ! + pressure gradient (ibp)  

      CALL qc_0y1_sp_M (mm, jj, jj_L, -one,  CC) ! - velocity divergence
      
      ALLOCATE (CCe_save(SIZE(CC%e)))

      CCe_save = CC%e  ! STORE THE CONSTANT CONTRIBUTION (Stokes operator)
      
   END IF      

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (ibp)  
   !-------------ADD THE ITERATION DEPENDENDT PART------------------- 

   CC%e = CCe_save

   IF (PRESENT(Re_U)) THEN 
      CALL qc_oseen2y_sp_M (mm, jj, Re_U,  CC) ! linearized terms
   END IF
   
   CALL Dirichlet_c_M (np, js_Axis, js_D,  CC)

   IF (DESINGULARIZE) THEN    
      ! reduction of the row of the last equation to 
      ! the diagonal element alone, set equal to 1   
      DO p = CC%i(Nx), CC%i(Nx + 1) - 1
         CC%e(p) = 0
         IF (CC%j(p) == Nx) CC%e(p) = 1
      ENDDO
   ENDIF
   
END SUBROUTINE ComputeJacobianMatrix_OLD

!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

SUBROUTINE extract_Dirichlet_c (np, js_Axis, js_D, xx,  old_us_D)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   INTEGER,             DIMENSION(:), INTENT(IN)  :: js_Axis
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D
   REAL(KIND=8),        DIMENSION(:), INTENT(IN)  :: xx
   TYPE(dyn_real_line), DIMENSION(:)              :: old_us_D

   INTEGER :: k

   DO k = 1, 3

      old_us_D(k)%DRL = xx(js_D(k)%DIL + (k-1) * np)
   
   ENDDO
   
END SUBROUTINE extract_Dirichlet_c

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_DIFF (np, js_Axis, js_D, us_D, old_us_D,  xx)
!===========================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   INTEGER,             DIMENSION(:), INTENT(IN)  :: js_Axis
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: us_D
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: old_us_D
   REAL(KIND=8),        DIMENSION(:)              :: xx

   INTEGER :: k

   DO k = 1, 3

      IF (k == 2  .OR.  k == 3) xx(js_Axis + (k-1)*np) = 0 ! Theoretically not needed

      xx(js_D(k)%DIL + (k-1) * np) = us_D(k)%DRL - old_us_D(k)%DRL
   
   ENDDO
   
END SUBROUTINE Dirichlet_c_DIFF

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c (np, js_Axis, js_D, us_D,  xx)
!====================================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   INTEGER,             DIMENSION(:), INTENT(IN)  :: js_Axis
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D  
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: us_D  
   REAL(KIND=8),        DIMENSION(:), INTENT(OUT) :: xx

   INTEGER :: k

   DO k = 1, 3
  
      IF (k == 2  .OR.  k == 3) xx(js_Axis + (k-1)*np) = 0
  
      xx(js_D(k)%DIL + (k-1)*np) = us_D(k)%DRL
  
   ENDDO 

END SUBROUTINE Dirichlet_c

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_M_MASS (np, js_Axis, js_D,  CC)
!================================================

!  Modification of elements of selected rows of
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes js_Axis, jsx_D, jsy_D and jst_D

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                          INTENT(IN)    :: np 
   INTEGER,            DIMENSION(:), INTENT(IN)    :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN)    :: js_D   
   TYPE(CSR_MUMPS_Matrix),           INTENT(INOUT) :: CC  

   INTEGER :: k, n, i_, p

   
   DO k = 2, 3

      DO n = 1, SIZE(js_Axis)
    
         i_ = js_Axis(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 0d0
                
         ENDDO  
         
      ENDDO
   
   ENDDO
   
   
   DO k = 1, 3

      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 0d0
         
         ENDDO
         
      ENDDO
   
   ENDDO
   
   
END SUBROUTINE Dirichlet_c_M_MASS

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_M (np, js_Axis, js_D,  CC)
!================================================

!  Modification of elements of selected rows of
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes js_Axis, jsx_D, jsy_D and jst_D

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                          INTENT(IN)    :: np 
   INTEGER,            DIMENSION(:), INTENT(IN)    :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN)    :: js_D   
   TYPE(CSR_MUMPS_Matrix),           INTENT(INOUT) :: CC  

   INTEGER :: k, n, i_, p

   
   DO k = 2, 3

      DO n = 1, SIZE(js_Axis)
    
         i_ = js_Axis(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 1
                
         ENDDO  
         
      ENDDO
   
   ENDDO
   
   
   DO k = 1, 3

      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 1
         
         ENDDO
         
      ENDDO
   
   ENDDO
   
   
END SUBROUTINE Dirichlet_c_M

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_rc_M (np, js_Axis, js_D, diagE,  CC)
!========================================================

!  Modification of elements of selected rows AND columns
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes  js_D(1)%DIL  and  js_D(2)%DIL  

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                          INTENT(IN)    :: np
   INTEGER,            DIMENSION(:), INTENT(IN)    :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN)    :: js_D
   REAL(KIND=8),                     INTENT(IN)    :: diagE
   TYPE(CSR_MUMPS_Matrix),           INTENT(INOUT) :: CC 

   INTEGER :: k, n, i_, p

   DO k = 2, 3

      DO n = 1, SIZE(js_Axis)
    
         i_ = js_Axis(n) + (k-1)*np
    
         ! column
         WHERE ( CC%j == i_ )

            CC%e = 0

         ENDWHERE
    
         ! row
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = 1
                
         ENDDO  
         
      ENDDO
   
   ENDDO
   
  
   DO k = 1, 3
  
      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n)  +  (k-1) * np

         ! column
         WHERE ( CC%j == i_ )

            CC%e = 0

         ENDWHERE
    
         ! row
         DO p = CC%i(i_), CC%i(i_+1) - 1
         
            CC%e(p) = 0
         
            IF (CC%j(p) == i_) CC%e(p) = diagE
         
         ENDDO
    
      ENDDO
  
   ENDDO 
   
  
END SUBROUTINE Dirichlet_rc_M

!------------------------------------------------------------------------------

SUBROUTINE qc_0y0_zero_sp_M (m0, jj, alpha,  CC) 
!===============================================


!  alpha << w, _ >>   ===>   CC  
!  
!  mass matrix but ONLY on Diagonal blocks 
!
!  The last (third) block which is zero (no cumulation)
!
!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  
 
   REAL(KIND=8), DIMENSION(n_w, n_w, l_G) :: a_w_w_p

   REAL(KIND=8) :: x
   INTEGER      :: np, mm, l, m, ni, nj, i, j, p, i_, j_


   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_0y0_zero_sp_M  is implemented only in 2D'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF


   np = MAXVAL(jj)  

   
   DO ni = 1, n_w

      DO nj = 1, n_w
      
         a_w_w_p(ni, nj, :) = alpha * ww(ni,:) * ww(nj,:) * pp_w
      
      ENDDO

   ENDDO


   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G
      
         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               x = a_w_w_p(ni, nj, l) * JAC(m) * yy_G(l,m)

               ! diagonal block of the first block row
            
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! diagonal block of the second block row
               
               i_ = i + np;   j_ = j + np 
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! diagonal block of the third block row
               
               i_ = i + 2*np;   j_ = j + 2*np 
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_0y0_zero_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_1y1_sp_M (m0, jj, alpha,  CC) 
!==========================================

!  +  alpha [<< (Dxw), y Dx_ >>  +  < (D.w), y D._ >]   ===>   CC  
!                          
!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY
!
!  +  alpha  times
!
!  dwx/dx . dvx/dx  +  dwx/dy . dvx/dy   |   dwx/dx . dvy/dy  -  dwx/dy . dvy/dx
!
! -dwy/dx . dvx/dy  +  dwy/dy . dvx/dx   |   dwy/dx . dvy/dx  +  dwy/dy . dvy/dy 
!
!
!  dwx/dx . d_x/dx  +  dwx/dy . d_x/dy   |   dwx/dx . d_y/dy  -  dwx/dy . d_y/dx
!
! -dwy/dx . d_x/dy  +  dwy/dy . d_x/dx   |   dwy/dx . d_y/dx  +  dwy/dy . d_y/dy 
!

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8) :: alpha_pyJ,  x
   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_

   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_1y1_sp_M  is implemented only for axisymmetric'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF

   np = MAXVAL(jj)  
   
!  WRITE (*,*) 'np = ', np, 'Echo from qc_1y1_sp_M'
      

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         alpha_pyJ = alpha * pp_w(l) * yy_G(l,m) * JAC(m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               ! FIRST BLOCK-ROW
               !================
                   
               ! diagonal block (1,1)
               
               ! y dwx/dx . dvx/dx  +  y dwx/dy . dvx/dy
               ! y dwx/dx . d_x/dx  +  y dwx/dy . d_x/dy
                          
               x = alpha_pyJ * SUM(dwl(:,ni) * dwl(:,nj)) / JAC(m)**2

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off-diagonal block (1,2)  
                     
               j_ = j + np
              
               ! y dwx/dx (dvy/dy + vy/R)  -  y dwx/dy . dvy/dx
               ! y dwx/dx (d_y/dy + _y/R)  -  y dwx/dy . d_y/dx
             
               x = alpha_pyJ * ( (dwl(1,ni)/JAC(m)) * (dwl(2,nj)/JAC(m) + ww(nj,l)/yy_G(l,m))  & 
                               - dwl(2,ni) * dwl(1,nj) / JAC(m)**2 )
              
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off-diagonal block (1,3)  ZERO
             
              
               ! SECOND BLOCK-ROW
               !=================
              
               i_ = i + np;
               
               ! off-diagonal block (2,1) 
               
               ! y (dwy/dy + wy/R) dvx/dx  -  y dwy/dx . dvx/dy  
               ! y (dwy/dy + wy/R) d_x/dx  -  y dwy/dx . d_x/dy  
              
               x = alpha_pyJ * ( (dwl(2,ni)/JAC(m) + ww(ni,l)/yy_G(l,m)) * dwl(1,nj)/JAC(m)  &
                                - dwl(1,ni) * dwl(2,nj) / JAC(m)**2 )
                
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO                                     

               ! diagonal block (2,2)
              
               j_ = j + np 
              
               ! y dwy/dx . dvy/dx  +  y (dwy/dy + wy/R) (dvy/dy + vy/R) 
               ! y dwy/dx . d_y/dx  +  y (dwy/dy + wy/R) (d_y/dy + _y/R)
             
               x = alpha_pyJ * ( dwl(1,ni) * dwl(1,nj) / JAC(m)**2  &
                                + (dwl(2,ni)/JAC(m) + ww(ni,l)/yy_G(l,m))  &
                                * (dwl(2,nj)/JAC(m) + ww(nj,l)/yy_G(l,m)) )
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! off-diagonal block (2,3) ZERO 
              
                  
               ! THIRD BLOCK-ROW
               !================
              
               i_ = i + 2*np;
               
               ! off-diagonal block (3,1) ZERO 
               
               ! off-diagonal block (3,2) ZERO 
            
               ! diagonal block (3,3) 
               
               j_ = j + 2*np 
              
               ! x = idem, vedi sopra  
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_1y1_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_1y1_sp_gg_M (m0, jj, alpha,  CC) 
!==========================================

!  +  alpha [<< (Dw), y D_ >> ]   ===>   CC  
!                          
!  ===>   CC   cumulative       

!  TWO-DIMENSIONAL VERSION ONLY
!

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER, DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8),            INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),  INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8) :: alpha_pyJ,  x
   INTEGER      :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_

   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_1y1_sp_gg_M  is implemented only for axisymmetric'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF

   np = MAXVAL(jj)  
   
!  WRITE (*,*) 'np = ', np, 'Echo from qc_1y1_sp_M'
      

   DO mm = 1, SIZE(m0);  m = m0(mm)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         alpha_pyJ = alpha * pp_w(l) * yy_G(l,m) * JAC(m)

         DO ni = 1, n_w;  i = jj(ni, m)

            DO nj = 1, n_w;  j = jj(nj, m)

               ! FIRST BLOCK-ROW
               !================
                   
               ! diagonal block (1,1)
               
               ! y dwx/dx . dvx/dx  +  y dwx/dy . dvx/dy
               ! y dwx/dx . d_x/dx  +  y dwx/dy . d_x/dy
                          
               !x = alpha_pyJ * SUM(dwl(:,ni) * dwl(:,nj)) / JAC(m)**2
               x = alpha_pyJ *  ( dwl(1,ni) * dwl(1,nj) / JAC(m)**2  &
                                + dwl(2,ni) * dwl(2,nj) / JAC(m)**2  )

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off-diagonal block (1,2)    ZERO
                     
               
               ! off-diagonal block (1,3)  ZERO
             
              
               ! SECOND BLOCK-ROW
               !=================
              
               i_ = i + np;
               
               ! off-diagonal block (2,1) ZERO  

               ! diagonal block (2,2)
              
               j_ = j + np 
              
               ! y dwy/dx . dvy/dx  +  y (dwy/dy . dvy/dy + wy/R . vy/R) 
               ! y dwy/dx . d_y/dx  +  y (dwy/dy . d_y/dy + wy/R . _y/R)
             
               x = alpha_pyJ *  (  dwl(1,ni) * dwl(1,nj) / JAC(m)**2  &
                                +  dwl(2,ni) * dwl(2,nj) / JAC(m)**2  &
                                +  ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2)
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! off-diagonal block (2,3) ZERO 
              
                  
               ! THIRD BLOCK-ROW
               !================
              
               i_ = i + 2*np;
               
               ! off-diagonal block (3,1) ZERO 
               
               ! off-diagonal block (3,2) ZERO 
            
               ! diagonal block (3,3) 
               
               j_ = j + 2*np 
              
               ! x = idem, vedi sopra:
               !
               ! y dwt/dx . dvt/dx  +  y (dwt/dy . dvt/dy + wt/R . vt/R) 
               ! y dwt/dx . d_t/dx  +  y (dwt/dy . d_t/dy + wt/R . _t/R)
               !
               ! x = alpha_pyJ *    (dwl(1,ni) * dwl(1,nj) / JAC(m)**2  &
               !                  +  dwl(2,ni) * dwl(2,nj) / JAC(m)**2  &
               !                  +  ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2)

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_1y1_sp_gg_M

!------------------------------------------------------------------------------

SUBROUTINE qc_1y1_sp_gg (m0, jj, gg, alpha,  v0)   !  BEWARE: NOT  _M 
!======================================   

!  << (Dw), y Dg >> 

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: gg
   REAL(KIND=8),                 INTENT(IN)  :: alpha
   REAL(KIND=8), DIMENSION(:,:)              :: v0

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   REAL(KIND=8), DIMENSION(3, n_w)   :: ggm
   REAL(KIND=8), DIMENSION(3, k_d)   :: dgl
   REAL(KIND=8), DIMENSION(3)        :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dwdgl_k
   INTEGER,      DIMENSION(n_w)      :: jjm
   
   REAL(KIND=8) :: z1, z2
 
   INTEGER :: mm, m, n, k, k1, l
      
!   v0 = 0

   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)

      ggm = gg(:,jjm)
  
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO k = 2, 3  ! g_x == g_z == g(k=1) non serve
            gl(k) = SUM(ggm(k,:) * ww(:,l))
         ENDDO

         DO k = 1, 3
            DO k1 = 1, k_d ! k_d == 2
               dgl(k,k1) = SUM(ggm(k,:) * dwl(k1,:)) 
            ENDDO
         ENDDO

         z1 = yy_G(l,m) * pp_w(l) / JAC(m)

         z2 = JAC(m) * pp_w(l) / yy_G(l,m)

         DO k = 1, 3
          
            DO n = 1, n_w
               dwdgl_k(n) = SUM(dwl(:,n) * dgl(k,:))
            ENDDO

            v0(k, jjm)  =  v0(k, jjm)  +  alpha * dwdgl_k * z1

            SELECT CASE (k)

               CASE(2); v0(2, jjm)  =  v0(2, jjm)  +  alpha * ww(:,l) * gl(2) * z2

               CASE(3); v0(3, jjm)  =  v0(3, jjm)  +  alpha * ww(:,l) * gl(3) * z2

            END SELECT

         ENDDO
         
      ENDDO

   ENDDO

END SUBROUTINE qc_1y1_sp_gg

!------------------------------------------------------------------------------

SUBROUTINE qc_advecy_sp_M (m0, jj, gg,  CC) 
!===========================================

!  ACCUMULATES CONTRIBUTIONS ONLY TO DIAGONAL BLOCKS

!  +  << w, y (g.D)_ >>   ===>   CC    /ALL VECTORS/
!    
!  ACCUMULATES CONTRIBUTIONS TO ALL FOUR VELOCITY BLOCKS
!                      
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl   
   REAL(KIND=8), DIMENSION(3, n_w)   :: ggm
   REAL(KIND=8), DIMENSION(3)        :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dgl_py
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k3, k, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) ::xd ! xd --> diagonal contribution

      WRITE(*,*) '***************************************'
      WRITE(*,*) '*** Warning:                        ***'
      WRITE(*,*) '*** implemented by Jacopo Canton    ***'
      WRITE(*,*) '*** use with care and without trust ***'
      WRITE(*,*) '***************************************'

   np = SIZE(gg, 2)

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
    
         ! nell'operatore d'advezione intervengono solo le prime            
         ! due componenti g(1) = gx = gz  e  g(2) = gy = gR 
     
         DO n = 1, n_w
            dgl_py(n) = SUM(gl(1:2) * dwl(:,n)) * pp_w(l) * yy_G(l,m)
         ENDDO

         DO ni = 1, n_w;  i = jjm(ni)
                           
            DO nj = 1, n_w;  j = jjm(nj)
               
               xd = ww(ni,l) * dgl_py(nj) 
               
               ! FIRST BLOCK ROW   
               !================    
               
               ! block (1,1)
                
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xd;  EXIT;  ENDIF
               ENDDO

               ! block (1,2)  ZERO

               ! block (1,3)  ZERO
             
                  
               ! SECOND BLOCK ROW
               !=================
               i_ = i + np 
         
               ! block (2,1)  ZERO

               ! block (2,2)

               j_ = j + np
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd;  EXIT;  ENDIF
               ENDDO

               ! block (2,3)  ZERO


               ! THIRD BLOCK ROW
               !================
               i_ = i + 2*np 
            
               ! block (3,1)  ZERO
            
               ! block (3,2)  ZERO
            
               ! block (3,3)
              
               j_ = j + 2*np
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd;  EXIT;  ENDIF
               ENDDO

            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_advecy_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_oseen2y_sp_M (m0, jj, gg,  CC) 
!===========================================

!  ACCUMULATES CONTRIBUTIONS ONLY TO DIAGONAL BLOCKS

!  +  << w, y (g.D)_ >>  +  << w, y (_.D)g >>   ===>   CC    /ALL VECTORS/
!    
!  ACCUMULATES CONTRIBUTIONS TO ALL FOUR VELOCITY BLOCKS
!                      
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl   
   REAL(KIND=8), DIMENSION(3, n_w)   :: ggm
   REAL(KIND=8), DIMENSION(k_d, 3)   :: dglo_py
   REAL(KIND=8), DIMENSION(3)        :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dgl_py
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k3, k, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) :: wij, xd, xo  ! xd --> diagonal contribution
                                ! xo --> off diagonal contribution   

   np = SIZE(gg, 2)

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
    
         ! nell'operatore d'advezione intervengono solo le prime            
         ! due componenti g(1) = gx = gz  e  g(2) = gy = gR 
     
         DO n = 1, n_w
            dgl_py(n) = SUM(gl(1:2) * dwl(:,n)) * pp_w(l) * yy_G(l,m)
         ENDDO

         DO k = 1, k_d
            DO k3 = 1, 3
               dglo_py(k, k3) = SUM(dwl(k,:) * ggm(k3,:)) * pp_w(l) * yy_G(l,m)        
            ENDDO
         ENDDO
         
         DO ni = 1, n_w;  i = jjm(ni)
                           
            DO nj = 1, n_w;  j = jjm(nj)
               
               xd = ww(ni,l) * dgl_py(nj) 
               
               wij = ww(ni,l) * ww(nj,l)
               
               ! FIRST BLOCK ROW   
               !================    
               
               ! block (1,1)   xo =  _ dgx/dx
                
               xo = wij * dglo_py(1,1) 
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xd + xo;  EXIT;  ENDIF
               ENDDO

               ! block (1,2)   xo =  _ dgx/dy              

               j_ = j + np

               xo = wij * dglo_py(2,1)      
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO
                              
               ! block (1,3)  ZERO
             
                  
               ! SECOND BLOCK ROW
               !=================
               i_ = i + np 
         
               ! block (2,1)   xo =  _ dgy/dx 
       
               xo = wij * dglo_py(1,2)     
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO

               ! block (2,2)   xo =  _ dgy/dy 

               j_ = j + np
            
               xo = wij * dglo_py(2,2) 
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd + xo;  EXIT;  ENDIF
               ENDDO

               ! block (2,3)   xo = - 2 _ gt/R  
                                    ! MINUS
               j_ = j + 2*np
            
               xo = - 2 * wij * gl(3) * pp_w(l) * JAC(m)
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO


               ! THIRD BLOCK ROW
               !================
               i_ = i + 2*np 
            
               ! block (3,1)   xo = _ dgt/dx
            
               xo = wij * dglo_py(1,3) 
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO
               
               ! block (3,2)   xo = _ (dgt/dy + gt/R)  
            
               j_ = j + np
            
               xo = wij * (dglo_py(2,3)  +  gl(3) * pp_w(l) * JAC(m)) 
            
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO

            
               ! block (3,3)   xo = _ gy/R 
              
               j_ = j + 2*np
            
               xo = wij * gl(2) * pp_w(l) * JAC(m)
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd + xo;  EXIT;  ENDIF
               ENDDO

            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_oseen2y_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_1y0_sp_M (m0, jj, jj_L, alpha,  CC)
!================================================

!  +  alpha ( < D.w, y _L >  +  < w(2), _L > )   ===>   CC
!                          
!  ===>   CC   cumulative      

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, jj_L
   REAL(KIND=8),                 INTENT(IN)    :: alpha 
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER      :: np, mm, m, l, n, k, ni, nj, i, j, p, i_, j_ 
   REAL(KIND=8) :: x
 
 
   np = MAXVAL(jj)  


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


   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO

         DO ni = 1, n_w;   i = jjm(ni);   i_ = i + np
              
            DO nj = 1, n_w_L;   j = jjm_L(nj);   j_ = j + 3*np

               ! first rectangular off-diagonal block of the last block-column  
                 
               x = alpha * dwl(1,ni) * w_L(nj,l) * pp_w(l) * yy_G(l,m)
                 
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
               ! second rectangular off-diagonal block of the last block-column  
   
               !  +  alpha ( < D.w, y _L >  +  < w(2), _L > )   ===>   CC  
            
               ! x1 = alpha * dwl(2,ni) * w_L(nj,l) * pp_w(l) * yy_G(l,m)
               ! x2 = alpha *  ww(ni,l) * w_L(nj,l) * pp_w(l) * JAC(m)
              
               x = alpha * (dwl(2,ni) * yy_G(l,m)  +  ww(ni,l) * JAC(m))  &                        
                         * w_L(nj,l) * pp_w(l)
               
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
             
              ! The third rectangular block of the fourth block-column is zero  
   
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_1y0_sp_M

!------------------------------------------------------------------------------

SUBROUTINE qc_0y1_sp_M (m0, jj, jj_L, alpha,  CC)
!================================================

!  +  alpha  ( < w_L, y D._ >  +  < w_L, __(2) > )  ===>   CC   
!                          
!  ===>   CC   cumulative      

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)    :: m0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jj, jj_L
   REAL(KIND=8),                 INTENT(IN)    :: alpha
   TYPE(CSR_MUMPS_Matrix),       INTENT(INOUT) :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER      :: np, mm, m, l, n, k, ni, nj, i, j, p, i_, j_ 
   REAL(KIND=8) :: x  

   np = MAXVAL(jj)  
  

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


   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)
      jjm_L = jj_L(:,m)

      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l))
            ENDDO
         ENDDO
 
         DO ni = 1, n_w_L;   i = jjm_L(ni);   i_ = i + 3*np

            DO nj = 1, n_w;   j = jjm(nj);   j_ = j + np
               
               ! first rectangular off-diagonal block of the bottom block-row  
              
               x = alpha * w_L(ni,l) * pp_w(l) * yy_G(l,m) * dwl(1,nj) 
                 
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
  
               ! second rectangular off-diagonal block of the bottom block-row
                 
               ! x1 = alpha * w_L(ni,l) * pp_w(l) * yy_G(l,m) * dwl(2,nj)             
               ! x2 = alpha * w_L(ni,l) * pp_w(l) *  ww(nj,l) * JAC(m)
                
               x = alpha * w_L(ni,l) * pp_w(l)  &
                         * (dwl(2,nj) * yy_G(l,m)  +  ww(nj,l) * JAC(m))
                 
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
               ! The third rectangular block of the fourth block-row is zero 
              
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_0y1_sp_M

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE qc_ny0_sp_s (ms0, jjs, iis, fs,  v0)
!==============================================

!  << n.ws, y fs >>_s   ===>   v0   incremental accumulation of boundary terms

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)    :: jjs, iis
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: fs
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: v0

   INTEGER :: mm, ms, ls, k
   REAL(KIND=8) :: fls, x

   DO mm = 1, SIZE(ms0);  ms = ms0(mm)
     
      DO ls = 1, l_Gs

         fls = SUM(fs(iis(:,ms)) * wws(:,ls)) 

         x = fls * JACs(ms) * pp_ws(ls) * yy_Gs(ls,ms) 

         DO k = 1, k_d     

            v0(k, jjs(:,ms)) = v0(k, jjs(:,ms)) + (normals(k, ls,ms)*x) * wws(:,ls) 

         ENDDO

      ENDDO
   
   ENDDO

       
END SUBROUTINE qc_ny0_sp_s

!------------------------------------------------------------------------------

SUBROUTINE qc_ty0_sp_s (ms0, jjs, iis, gzs,  v0)
!===============================================

!  << nxws, y gzs >>_s   ===>   v0   incremental accumulation of boundary terms
  
   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN)  :: ms0
   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jjs, iis
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: gzs
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: v0

   REAL(KIND=8) :: gzls, x
   INTEGER :: mm, ms, ls
   
   IF (k_d /= 2) THEN
   
      WRITE (*,*) ' qc_ty0_sp_s implemented only  in' 
      WRITE (*,*) ' the two-dimensional case: STOP.'
      STOP
       
   ENDIF
   
   DO mm = 1, SIZE(ms0);  ms = ms0(mm)

      DO ls = 1, l_Gs
         
         gzls = SUM(gzs(iis(:,ms)) * wws(:,ls))
             
         x = gzls * JACs(ms) * pp_ws(ls) * yy_Gs(ls,ms) 
             
         v0(1, jjs(:,ms)) = v0(1, jjs(:,ms)) - (normals(2, ls,ms)*x) * wws(:,ls) 
     
         v0(2, jjs(:,ms)) = v0(2, jjs(:,ms)) + (normals(1, ls,ms)*x) * wws(:,ls) 
     
      ENDDO

   ENDDO


END SUBROUTINE qc_ty0_sp_s

!******************************************************************************
!******************************************************************************
! 3D subroutines follow
!******************************************************************************
!******************************************************************************

SUBROUTINE ComputeJacobianMatrix_3d(np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, beta, CC,  Re, UU)
!===============================================================================

   IMPLICIT NONE
   
   !-----------------------------------------------------------------------!
   INTEGER,                            INTENT(IN) :: np 
   INTEGER,            DIMENSION(:),   INTENT(IN) :: mm
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: jj_L
   INTEGER,   POINTER, DIMENSION(:)               :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:),   INTENT(IN) :: js_D 
   LOGICAL,                            INTENT(IN) :: DESINGULARIZE
   INTEGER,                            INTENT(IN) :: beta
   
   TYPE(CSR_MUMPS_Complex_Matrix)                 :: CC

   REAL(KIND=8),                           INTENT(IN) :: Re
   REAL(KIND=8), OPTIONAL, DIMENSION(:,:), INTENT(IN) :: UU
   !-----------------------------------------------------------------------!

   REAL(KIND=8), PARAMETER :: one = 1d0
   INTEGER :: Nx, p
   

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  Re [(U.G)_ + (_.G)U)]  +  K_  +  G_ (weak)  
   !             CC  <---  - V._   
   !-------------ONLY THE CONSTANT CONTRIBUTION---------------------- 

   Nx = SIZE(CC%i) - 1

   CC%e = 0

   CALL qc_1y1_sp_gg_3d_M (mm, jj,  1d0/Re, beta,  CC) ! + stiffness (GRAD:GRAD)

   CALL qc_1y0_sp_3d_M (mm, jj, jj_L, -one, beta,  CC) ! + pressure gradient (ibp)  

   CALL qc_0y1_sp_3d_M (mm, jj, jj_L, -one, beta,  CC) ! - velocity divergence
      

   !-----------------------------------------------------------------
   !-------------GENERATION OF THE JACOBIAN MATRIX-------------------
   !-------------OF THE COUPLED EQUATION SYSTEM----------------------            
   !             CC  <---  [(U.G)_ + (_.G)U)]  + 1/Re *  K_  +  G_ (ibp)  
   !-------------ADD THE ITERATION DEPENDENDT PART------------------- 


   IF (PRESENT(UU)) THEN 
      CALL qc_oseen2y_sp_3d_M (mm, jj, UU, beta,  CC) ! linearized terms
   END IF
   
   CALL Dirichlet_rc_3d_M (np, js_Axis, js_D, 1d0,  CC)

   IF (DESINGULARIZE) THEN
      ! reduction of the row of the last equation to 
      ! the diagonal element alone, set equal to 1   
      DO p = CC%i(Nx), CC%i(Nx + 1) - 1
         CC%e(p) = 0
         IF (CC%j(p) == Nx) CC%e(p) = 1
      ENDDO
   ENDIF
   
END SUBROUTINE ComputeJacobianMatrix_3d

!------------------------------------------------------------------------------

SUBROUTINE qc_1y1_sp_gg_3d_M (m0, jj, alpha, beta,  CC) 
!==========================================

!  +  alpha [<< (Dw), y D_ >> ]   ===>   CC
!
!  beta is the azimutal wave number
!
!  ===>   CC   cumulative

!  TWO-DIMENSIONAL VERSION ONLY

   USE Gauss_points

   IMPLICIT NONE

   INTEGER, DIMENSION(:),         INTENT(IN) :: m0    ! (me)      elements indices
   INTEGER, DIMENSION(:,:),       INTENT(IN) :: jj    ! (n_w, me) nodes of the parabolic element
   REAL(KIND=8),                  INTENT(IN) :: alpha
   INTEGER,                       INTENT(IN) :: beta
   TYPE(CSR_MUMPS_Complex_Matrix)            :: CC  

   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl
   
   REAL(KIND=8)    :: alpha_pyJ
   COMPLEX(KIND=8) :: x
   INTEGER         :: np, mm, k, l, m, n, ni, nj, i, j, p, i_, j_

   IF (k_d /= 2) THEN
   
      WRITE (*,*) 'qc_1y1_sp_gg_3d_M  is implemented only for 2d meshes'
      WRITE (*,*) 'STOP.'
      STOP
   
   ENDIF

   np = MAXVAL(jj) ! number of parabolic nodes
   
	! outer cycle on elements.  m = current element index
	!
   DO mm = 1, SIZE(m0);  m = m0(mm)

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l)) ! dwl (2,n_w) = derivatives of the parabolic weight functions
            ENDDO
         ENDDO

         ! coefficient
         !
         alpha_pyJ = alpha * pp_w(l) * yy_G(l,m) * JAC(m)

			! cycle on parabolic nodes of the current element.  i = current node index
			!
         DO ni = 1, n_w;  i = jj(ni, m)

			   ! cycle on parabolic nodes of the current element.  j = current node index
			   !
            DO nj = 1, n_w;  j = jj(nj, m)

               ! FIRST BLOCK-ROW
               !================
                   
               ! diagonal block (1,1)
               
               ! y * ( dwx/dx . dvx/dx  +  dwx/dy . dvx/dy )  + 1/R * ( beta**2 wx . vx )
               ! y * ( dwx/dx . d_x/dx  +  dwx/dy . d_x/dy )  + 1/R * ( beta**2 wx . _x )
                          
               !x = alpha_pyJ * SUM(dwl(:,ni) * dwl(:,nj)) / JAC(m)**2
               x = alpha_pyJ *  (            dwl(1,ni) * dwl(1,nj) / JAC(m)**2  &
                                +            dwl(2,ni) * dwl(2,nj) / JAC(m)**2  &
                                +  beta**2 * ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2)

               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
               
               ! off-diagonal block (1,2)  ZERO
                     
               
               ! off-diagonal block (1,3)  ZERO
             
              
               ! SECOND BLOCK-ROW
               !=================
              
               i_ = i + np;
               
               ! off-diagonal block (2,1) ZERO  

               ! diagonal block (2,2)
              
               j_ = j + np 
              
               ! y * ( dwy/dx . dvy/dx  +  dwy/dy . dvy/dy )  + 1/R * ( wy . vy  +  beta**2 wy . vy )
               ! y * ( dwy/dx . d_y/dx  +  dwy/dy . d_y/dy )  + 1/R * ( wy . _y  +  beta**2 wy . _y )
             
               x = alpha_pyJ *  (            dwl(1,ni) * dwl(1,nj) / JAC(m)**2    &
                                +            dwl(2,ni) * dwl(2,nj) / JAC(m)**2    &
                                +            ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2 &
                                +  beta**2 * ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2)
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
              
               ! off-diagonal block (2,3)

               j_ = j + 2*np

               ! + 1/R * 2 i beta wy . vt )
               ! + 1/R * 2 i beta wy . _t )

               x = alpha_pyJ * 2 * CMPLX(0d0,1d0,KIND=8) * beta * ww(ni,l) * ww(nj,l) / yy_G(l,m)**2

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

                  
               ! THIRD BLOCK-ROW
               !================
              
               i_ = i + 2*np;
               
               ! off-diagonal block (3,1) ZERO
               
               ! off-diagonal block (3,2)

               j_ = j + np

               ! - 1/R * 2 i beta wt . vy )
               ! - 1/R * 2 i beta wt . _y )

               x = - alpha_pyJ * 2 * CMPLX(0d0,1d0,KIND=8) * beta * ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
            
               ! diagonal block (3,3)
               
               j_ = j + 2*np

               ! y * ( dwy/dx . dvy/dx  +  dwy/dy . dvy/dy )  + 1/R * ( wy . vy  +  beta**2 wy . vy )
               ! y * ( dwy/dx . d_y/dx  +  dwy/dy . d_y/dy )  + 1/R * ( wy . _y  +  beta**2 wy . _y )
             
               x = alpha_pyJ *  (            dwl(1,ni) * dwl(1,nj) / JAC(m)**2    &
                                +            dwl(2,ni) * dwl(2,nj) / JAC(m)**2    &
                                +            ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2 &
                                +  beta**2 * ww(ni,l)  * ww(nj,l)  / yy_G(l,m)**2)
              
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_1y1_sp_gg_3d_M

!------------------------------------------------------------------------------

SUBROUTINE qc_oseen2y_sp_3d_M (m0, jj, gg, beta,  CC) 
!===========================================

!  +  << w, y (g.D)_ >>  +  << w, y (_.D)g >>   ===>   CC    /ALL VECTORS/
!
!  beta is the azimutal wave number
!    
!  ACCUMULATES CONTRIBUTIONS TO ALL NINE VELOCITY BLOCKS
!                      
!  ===>   CC   cumulative      

   USE Gauss_points

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),    INTENT(IN) :: m0   ! (me)      elements indices
   INTEGER,      DIMENSION(:,:),  INTENT(IN) :: jj   ! (n_w, me) nodes of the parabolic element
   REAL(KIND=8), DIMENSION(:,:),  INTENT(IN) :: gg   ! (3,np)    vector function known on the parabolic nodes
   INTEGER,                       INTENT(IN) :: beta
   TYPE(CSR_MUMPS_Complex_Matrix)            :: CC  
   
   REAL(KIND=8), DIMENSION(k_d, n_w) :: dwl   
   REAL(KIND=8), DIMENSION(3, n_w)   :: ggm
   REAL(KIND=8), DIMENSION(k_d, 3)   :: dglo_py
   REAL(KIND=8), DIMENSION(3)        :: gl
   REAL(KIND=8), DIMENSION(n_w)      :: dgl_py
   INTEGER,      DIMENSION(n_w)      :: jjm

   INTEGER      :: np, mm, k3, k, l, m, n, ni, nj, i, j, p, i_, j_
   REAL(KIND=8) :: wij, xd, xo  ! xd --> diagonal contribution
                                ! xo --> off diagonal contribution
   COMPLEX(KIND=8) :: x3 ! azimutal contribution

   WRITE(*,*) '    qc_oseen2y_sp_3d_M'
   WRITE(*,*) '*** assuming NON SWIRLING and AXISYMMETRIC base flow ***'

   np = SIZE(gg, 2) ! number of parabolic nodes


	! outer cycle on elements.  m = current element index
	!
   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm = jj(:,m)    ! jjm (n_w)   = indices of the parabolic nodes of the current element
      ggm = gg(:, jjm) ! ggm (3,n_w) = vector function known on the parabolic nodes of the current element

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G
         
         DO k3 = 1, 3
            gl(k3) = SUM(ggm(k3,:) * ww(:,l)) ! gl (3) = g integrated over the current element [ ww (n_w,l_G) = weight functions ]
         ENDDO

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l)) ! dwl (2,n_w) = derivatives of the parabolic weight functions
            ENDDO
         ENDDO
    
         ! nell'operatore d'advezione intervengono solo le prime
         ! due componenti g(1) = gx = gz  e  g(2) = gy = gR
         !
			! y * ( g(1) * d_ / dx  +  g(2) * d_ / dy )
			!
         DO n = 1, n_w
            dgl_py(n) = SUM(gl(1:2) * dwl(:,n)) * pp_w(l) * yy_G(l,m)
         ENDDO

			! dglo_py(:,1) = y  * dgx / dx
			! dglo_py(:,2) = y  * dgy / dy
			! dglo_py(:,3) = y  * dgz / dz
			!
         DO k = 1, k_d
            DO k3 = 1, 3
               dglo_py(k, k3) = SUM(dwl(k,:) * ggm(k3,:)) * pp_w(l) * yy_G(l,m)
            ENDDO
         ENDDO
         
			! cycle on parabolic nodes of the current element.  i = current node index
			!
         DO ni = 1, n_w;  i = jjm(ni)
                           
			   ! cycle on parabolic nodes of the current element.  j = current node index
				! 
            DO nj = 1, n_w;  j = jjm(nj)
               
					! xd = y * ( g(1) * w_ * d_ / dx  +  g(2) * w_ * d_ / dy )
					!
               xd = ww(ni,l) * dgl_py(nj)

               wij = ww(ni,l) * ww(nj,l)
               
               ! FIRST BLOCK ROW
               !================

               ! block (1,1)
                                           ! *** ONLY IF SWIRLING BASE FLOW ***
               ! xd  +  y * vx * ( dgx/dx  +  i beta/R * gt )
                
               xo = wij * dglo_py(1,1)

               x3 = 0d0 !+ CMPLX(0d0,1d0,KIND=8) * beta * wij * gl(3) * pp_w(l) * JAC(m)
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xd + xo + x3;  EXIT;  ENDIF
               ENDDO

               ! block (1,2)

               j_ = j + np

               ! y * vy * dgx/dy

               xo = wij * dglo_py(2,1)
               
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO
                              
               ! *** ONLY IF NON-AXISYMMETRIC BASE FLOW ***
               ! block (1,3)
               !
               !j_ = j + 2*np
             
                  
               ! SECOND BLOCK ROW
               !=================
               i_ = i + np 
         
               ! block (2,1)

               ! y * vx * dgy/dx 
       
               xo = wij * dglo_py(1,2)     
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               ENDDO

               ! block (2,2)

               j_ = j + np
                                           ! *** ONLY IF SWIRLING BASE FLOW ***
               ! xd  +  y * vy * ( dgy/dy  +  i beta/R * gt )
            
               xo = wij * dglo_py(2,2)

               x3 = 0d0 !+ CMPLX(0d0,1d0,KIND=8) * beta * wij * gl(3) * pp_w(l) * JAC(m)
             
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd + xo + x3;  EXIT;  ENDIF
               ENDDO

               ! *** ONLY IF SWIRLING (and NON-AXISYMMETRIC) BASE FLOW ***
               ! block (2,3)
               !
               !j_ = j + 2*np
               !
               !! y * ( - 2 * vt * gt/R )
               !
               !xo = - 2 * wij * gl(3) * pp_w(l) * JAC(m)
               !
               !DO p = CC%i(i_),  CC%i(i_+1) - 1
               !   IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               !ENDDO


               ! THIRD BLOCK ROW
               !================
               i_ = i + 2*np 
            
               ! *** ONLY IF SWIRLING BASE FLOW ***
               ! block (3,1)
               !
               !! y * vx * dgt/dx
               !
               !xo = wij * dglo_py(1,3) 
               !
               !DO p = CC%i(i_),  CC%i(i_+1) - 1
               !   IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               !ENDDO
               
               ! *** ONLY IF SWIRLING BASE FLOW ***
               ! block (3,2)
               !
               !j_ = j + np
               !
               !! y * vy * ( gt/R  +  dgt/dy )
               !
               !xo = wij * ( gl(3) * pp_w(l) * JAC(m)  +  dglo_py(2,3)) 
               !
               !DO p = CC%i(i_),  CC%i(i_+1) - 1
               !   IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xo;  EXIT;  ENDIF
               !ENDDO

            
               ! block (3,3)

               j_ = j + 2*np
                                         ! *** ONLY IF SWIRLING BASE FLOW ***
               ! xd  +  y * vt * ( gy/R  +  i beta/R * gy )
            
               xo = wij * gl(2) * pp_w(l) * JAC(m)

               x3 = 0d0 !+ CMPLX(0d0,1d0,KIND=8) * beta * wij * gl(3) * pp_w(l) * JAC(m)

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + xd + xo + x3;  EXIT;  ENDIF
               ENDDO

            
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_oseen2y_sp_3d_M

!------------------------------------------------------------------------------

SUBROUTINE qc_1y0_sp_3d_M (m0, jj, jj_L, alpha, beta,  CC)
!================================================

!  +  alpha ( < D.w, y _L >  +  < w(2), _L > )   ===>   CC
!
!  beta is the azimutal wave number
!
!  ===>   CC   cumulative

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),    INTENT(IN) :: m0       ! (me)      elements indices
   INTEGER,      DIMENSION(:,:),  INTENT(IN) :: jj, jj_L
   REAL(KIND=8),                  INTENT(IN) :: alpha
   INTEGER,                       INTENT(IN) :: beta
   TYPE(CSR_MUMPS_Complex_Matrix)            :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER         :: np, mm, m, l, n, k, ni, nj, i, j, p, j_ 
   COMPLEX(KIND=8) :: x
 
 
   np = MAXVAL(jj)  


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


	! outer cycle on elements.  m = current element index
	!
   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)    ! jjm   (n_w)     = indices of the parabolic nodes of the current element
      jjm_L = jj_L(:,m)    ! jjm_L (n_w_L)   = indices of the linear    nodes of the current element

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l)) ! dwl (2,n_w) = derivatives of the parabolic weight functions
            ENDDO
         ENDDO

			! cycle on parabolic nodes of the current element.  i = current node index
			!
         DO ni = 1, n_w;   i = jjm(ni)

			   ! cycle on linear nodes of the current element.  j = current node index
				! 
            DO nj = 1, n_w_L;   j = jjm_L(nj)

               ! last block-column
               !
               j_ = j + 3*np

               ! first rectangular off-diagonal block of the last block-column  

               ! x = alpha * dwl(1,ni) * w_L(nj,l) * pp_w(l) * yy_G(l,m)
                 
               x = alpha * dwl(1,ni) * w_L(nj,l) * pp_w(l) * yy_G(l,m)
                 
               DO p = CC%i(i),  CC%i(i+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
 
               ! second rectangular off-diagonal block of the last block-column  
   
               ! x1 = alpha * dwl(2,ni) * w_L(nj,l) * pp_w(l) * yy_G(l,m)
               ! x2 = alpha *  ww(ni,l) * w_L(nj,l) * pp_w(l) * JAC(m)
              
               x = alpha * w_L(nj,l) * pp_w(l) &
                         * (dwl(2,ni) * yy_G(l,m)  +  ww(ni,l) * JAC(m))

               DO p = CC%i(i+np),  CC%i(i+np+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
             
               ! third rectangular off-diagonal block of the last block-column  

               ! x = alpha * i beta * ww(ni,l) * w_L(nj,l) * pp_w(l) * JAC(m)
              
               x = alpha * CMPLX(0d0,1d0,KIND=8) * beta * ww(ni,l) * w_L(nj,l) * pp_w(l) * JAC(m) 
               
               DO p = CC%i(i+2*np),  CC%i(i+2*np+1) - 1
                  IF (CC%j(p) == j_) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
   
            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_1y0_sp_3d_M

!------------------------------------------------------------------------------

SUBROUTINE qc_0y1_sp_3d_M (m0, jj, jj_L, alpha, beta,  CC)
!================================================

!  +  alpha  ( < w_L, y D._ >  +  < w_L, __(2) > )  ===>   CC
!
!  beta is the azimutal wave number
!
!  ===>   CC   cumulative

   USE Gauss_points
   
   USE Gauss_points_L

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),    INTENT(IN) :: m0      ! (me)      elements indices
   INTEGER,      DIMENSION(:,:),  INTENT(IN) :: jj, jj_L
   REAL(KIND=8),                  INTENT(IN) :: alpha
   INTEGER,                       INTENT(IN) :: beta
   TYPE(CSR_MUMPS_Complex_Matrix)            :: CC  

   REAL(KIND=8), DIMENSION(SIZE(jj_L,1), l_G) :: w_L

   REAL(KIND=8), DIMENSION(k_d, n_w)     :: dwl
   INTEGER,      DIMENSION(n_w)          :: jjm
   INTEGER,      DIMENSION(SIZE(jj_L,1)) :: jjm_L
  
   INTEGER         :: np, mm, m, l, n, k, ni, nj, i, j, p, i_
   COMPLEX(KIND=8) :: x  


   np = MAXVAL(jj)  
  

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


	! outer cycle on elements.  m = current element index
	!
   DO mm = 1, SIZE(m0);  m = m0(mm)

      jjm   = jj  (:,m)    ! jjm   (n_w)     = indices of the parabolic nodes of the current element
      jjm_L = jj_L(:,m)    ! jjm_L (n_w_L)   = indices of the linear    nodes of the current element

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G

         DO k = 1, k_d
            DO n = 1, n_w
               dwl(k,n) = SUM(MNR(k,:,m) * Dw_re(:,n,l)) ! dwl (2,n_w) = derivatives of the parabolic weight functions
            ENDDO
         ENDDO
 
			! cycle on linear nodes of the current element.  i = current node index
			!
         DO ni = 1, n_w_L;   i = jjm_L(ni)

            ! bottom block-row
            !
            i_ = i + 3*np

			   ! cycle on parabolic nodes of the current element.  j = current node index
				! 
            DO nj = 1, n_w;   j = jjm(nj)
               
               ! first rectangular off-diagonal block of the bottom block-row  
              
               ! x = alpha * w_L(ni,l) * pp_w(l) * yy_G(l,m) * dwl(1,nj)

               x = alpha * w_L(ni,l) * pp_w(l) * yy_G(l,m) * dwl(1,nj)
                 
               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO
  
               ! second rectangular off-diagonal block of the bottom block-row

               ! x1 = alpha * w_L(ni,l) * pp_w(l) * yy_G(l,m) * dwl(2,nj)
               ! x2 = alpha * w_L(ni,l) * pp_w(l) *  ww(nj,l) * JAC(m)

               x = alpha * w_L(ni,l) * pp_w(l)  &
                         * (dwl(2,nj) * yy_G(l,m)  +  ww(nj,l) * JAC(m))

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j+np) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

               ! third rectangular off-diagonal block of the bottom block-row

               ! x = alpha * i beta * w_L(ni,l) * pp_w(l) * ww(nj,l) * JAC(m)

               x = alpha * CMPLX(0d0,1d0,KIND=8) * beta * w_L(ni,l) * pp_w(l) * ww(nj,l) * JAC(m)

               DO p = CC%i(i_),  CC%i(i_+1) - 1
                  IF (CC%j(p) == j+2*np) THEN;  CC%e(p) = CC%e(p) + x;  EXIT;  ENDIF
               ENDDO

            ENDDO

         ENDDO

      ENDDO

   ENDDO

END SUBROUTINE qc_0y1_sp_3d_M

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_3d (np, js_Axis, js_D, us_D,  xx)
!====================================================

   IMPLICIT NONE
   
   INTEGER,                           INTENT(IN)  :: np
   INTEGER,             DIMENSION(:), INTENT(IN)  :: js_Axis
   TYPE(dyn_int_line),  DIMENSION(:), INTENT(IN)  :: js_D  
   TYPE(dyn_real_line), DIMENSION(:), INTENT(IN)  :: us_D  
   COMPLEX(KIND=8),     DIMENSION(:)              :: xx

   INTEGER :: k

   DO k = 1, 3
  
      IF (k == 2  .OR.  k == 3) xx(js_Axis + (k-1)*np) = CMPLX(0d0,0d0,KIND=8)
  
      xx(js_D(k)%DIL + (k-1)*np) = us_D(k)%DRL
  
   ENDDO 

END SUBROUTINE Dirichlet_c_3d

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_c_3d_M (np, js_Axis, js_D,  CC)
!================================================

!  Modification of elements of selected rows of
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes js_Axis, jsx_D, jsy_D and jst_D

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                          INTENT(IN)    :: np 
   INTEGER,            DIMENSION(:), INTENT(IN)    :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN)    :: js_D   
   TYPE(CSR_MUMPS_Complex_Matrix)                  :: CC  

   INTEGER :: k, n, i_, p

   
   DO k = 2, 3

      DO n = 1, SIZE(js_Axis)
    
         i_ = js_Axis(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0d0
         
            IF (CC%j(p) == i_) CC%e(p) = 1d0
                
         ENDDO  
         
      ENDDO
   
   ENDDO
   
   
   DO k = 1, 3

      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n) + (k-1)*np
    
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0d0
         
            IF (CC%j(p) == i_) CC%e(p) = 1d0
         
         ENDDO
         
      ENDDO
   
   ENDDO
   
   
END SUBROUTINE Dirichlet_c_3d_M

!------------------------------------------------------------------------------

SUBROUTINE Dirichlet_rc_3d_M (np, js_Axis, js_D, diagE,  CC)
!========================================================

!  Modification of elements of selected rows AND columns
!  matrix CC to impose Dirichlet boundary conditions 
!  at the nodes  js_D(1)%DIL  and  js_D(2)%DIL  

   USE Gauss_points

   IMPLICIT NONE
  
   INTEGER,                          INTENT(IN)    :: np
   INTEGER,            DIMENSION(:), INTENT(IN)    :: js_Axis
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN)    :: js_D
   REAL(KIND=8),                     INTENT(IN)    :: diagE
   TYPE(CSR_MUMPS_Complex_Matrix)                  :: CC 

   INTEGER :: k, n, i_, p

   DO k = 2, 3

      DO n = 1, SIZE(js_Axis)
    
         i_ = js_Axis(n) + (k-1)*np
    
         ! column
         WHERE ( CC%j == i_ )

            CC%e = 0d0

         ENDWHERE
    
         ! row
         DO p = CC%i(i_), CC%i(i_+1) - 1
        
            CC%e(p) = 0d0
         
            IF (CC%j(p) == i_) CC%e(p) = 1d0
                
         ENDDO  
         
      ENDDO
   
   ENDDO
   
  
   DO k = 1, 3
  
      DO n = 1, SIZE(js_D(k)%DIL)
    
         i_ = js_D(k)%DIL(n)  +  (k-1) * np

         ! column
         WHERE ( CC%j == i_ )

            CC%e = 0d0

         ENDWHERE
    
         ! row
         DO p = CC%i(i_), CC%i(i_+1) - 1
         
            CC%e(p) = 0d0
         
            IF (CC%j(p) == i_) CC%e(p) = CMPLX(diagE, 0d0, KIND=8)
         
         ENDDO
    
      ENDDO
  
   ENDDO 
   
  
END SUBROUTINE Dirichlet_rc_3d_M

!==============================================================================
!==============================================================================
!==============================================================================

END MODULE qc_sp_M
