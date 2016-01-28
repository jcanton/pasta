MODULE newton

   USE par_solve_mumps
   USE qv_sp
   USE qc_sp_M
   USE qs_L_sp
   USE dynamic_structures
   USE sparse_matrix_profiles
   USE sparse_matrix_operations
   USE global_variables
   USE miscellaneous_subroutines
   USE prep_mesh_p1p2_sp
   USE start_sparse_kit
   USE Dirichlet_Neumann
   USE axisym_boundary_values
   USE case_dependent
   !
   USE loca_types
   USE loca_pd
   USE loca_bord


CONTAINS
!=======

FUNCTION nonlinear_solver_conwrap (x_vec, step_num, lambda, delta_s, con_ptr) &
   RESULT(num_newt_its)
!
! Put the call to your nonlinear solver here.
! Input:
!    x_vec     solution vector
!    con_ptr   pointer to continuation structure, cast to (void *)
!              must be passed to nonlinear solver and then passed
!              to bordering algorithms.
!    step_num  Continuation step number
! 
! Output:
!    x_vec     solution vector
! 
! Return Value:
!    num_newt_its  Number of Newton iterations needed for
!                  convergence, used to pick next step size.
!                  Negative value means nonlinear solver didn't converge.
! 

   IMPLICIT NONE

   ! input variables
   REAL(KIND=8), DIMENSION(0:Nx-1) :: x_vec ! FixMe check (or change) the indices of the vector
   INTEGER                         :: step_num
   REAL(KIND=8)                    :: lambda
   REAL(KIND=8)                    :: delta_s
   TYPE(con_struct), OPTIONAL      :: con_ptr
   ! output variables
   INTEGER                             :: num_newt_its

   ! common variables used
   ! Jacobian
   ! p_in
   ! velCmpnnts, np, np_L, Nx
   ! x0, u0, p0
   ! mm, jj, jj_L, js_D, zero_bvs_D
   ! DESINGULARIZE
   ! Re

   ! local variables
   INTEGER                                   :: n
   REAL(KIND=8)                              :: residual
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::     vv
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dx, ww, rhs, dx0


   INTEGER                         :: continuation_converged=0
   REAL(KIND=8), DIMENSION(0:Nx-1) :: x_hook       ! FixMe check (or change) the indices of the vector
   REAL(KIND=8), DIMENSION(0:Nx-1) :: delta_x_hook ! FixMe check (or change) the indices of the vector
   REAL(KIND=8)                    :: Reltol = 1d-3
   REAL(KIND=8)                    :: Abstol = 1d-8
   ! (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)


   ! executable statements

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to nonlinear_solver_conwrap'
   WRITE(*,*)


   ALLOCATE ( vv (velCmpnnts, np) )
   ALLOCATE ( ww (np_L) )
   ALLOCATE ( dx (Nx) )
   ALLOCATE ( rhs (Nx) )
   ALLOCATE ( dx0 (Nx) )
   dx0 = 1.29


   x0 = x_vec
   xx = x_vec

!************************************************************************************
!*** Newton't method in bi-incremental form is reported at the end of this module ***
!************************************************************************************

   !=============================================================
   ! start of NEWTON'S ITERATIONS IN NON-INCREMENTAL FORM
   !
   WRITE(*,*) '    Start of Newton''s iterations'
   WRITE(*,*) '    in NON-INCREMENTAL form'
   !
   DO n = 1, p_in%nwtn_maxite
      
      WRITE(*,*)
      WRITE(*,*) '    n = ', n

      IF (n == p_in%nwtn_maxite) THEN
       
         WRITE(*,*) '   ************************************************'
         WRITE(*,*) '   *Maximum number of Newton''s iterations reached:'
         WRITE(*,*) '   *ite_max = ', p_in%nwtn_maxite
         WRITE(*,*) '   ************************************************'
         WRITE(*,*)
      
      ENDIF
      
      CALL extract (x0,  u0)

!      CALL extract (x0,  u0, p0)
!      CALL vtk_plot_P2 (rr, jj, jj_L, u0, p0, trim(p_in%plot_directory) // 'iteSol.vtk')

      !------------------------------------------------------------------
      ! call case dependent subroutine
      !
      CALL case_newton_iteprocess(n, continuation_converged)

      !------------------------------------------------------------------
      !-------------GENERATION OF THE RIGHT-HAND SIDE--------------------
      !
      ! NON-INCREMENTAL FORM
      ! rhs <---  (u0 \dot \nabla)u0 + f
      !           0
      vv = 0
!      CALL extract (x0,  u0)
      CALL qv_0y01_sp (mm, jj, u0,  vv)
!      CALL qc_ty0_sp_s (ms_2, jjs, iis,  c_2,  vv)  !  cumulative
!      CALL qc_ny0_sp_s (ms_3, jjs, iis, -q_3,  vv)  !  cumulative

      u0(1,:) = volumeForcing(1,1)
      u0(2,:) = volumeForcing(2,1)
      u0(3,:) = volumeForcing(3,1)
      CALL qv_0y0_sp   (mm, jj, u0, 1d0, vv)

      u0(1,:) = volumeForcing(1,2)
      u0(2,:) = volumeForcing(2,2)
      u0(3,:) = volumeForcing(3,2)
      CALL qv_0y0_dR_sp(mm, jj, u0, 1d0, vv)
      
      ww = 0
      CALL collect (vv, ww,  dx) ! here dx is the RHS

      !------------------------------------------------------------------
      !-------------ENFORCING DIRICHLET BOUNDARY CONDITIONS ON THE RHS---
      !
      ! NON-INCREMENTAL FORM
      ! non-homogeneous boundary conditions
      !
      CALL Dirichlet_c (np, js_Axis, js_D, bvs_D,  dx)
      IF (DESINGULARIZE) dx(Nx) = 0d0

      !------------------------------------------------------------------
      !-------------COMPUTE RESIDUAL-------------------------------------

      ! rhs <-- (du0 \dot \nabla)du0
      !         div{u0}
      vv = 0
      CALL extract (dx0,  u0)
      CALL qv_0y01_sp   (mm, jj, u0, vv)
      ww = 0
      CALL collect (vv, ww,  rhs)
      CALL Dirichlet_c (np, js_Axis, js_D, zero_bvs_D,  rhs)
      IF (DESINGULARIZE) rhs(Nx) = 0d0

      residual = MAXVAL(ABS(rhs))
      
      WRITE(*,*) '    |res|_L-infty = ', residual

      !===================================================================
      !===================================================================
      IF (residual < p_in%nwtn_tol .AND. continuation_converged == 1) THEN
         x_vec = x0
         WRITE (*,*)
         WRITE (*,*) '    Residual on the converged solution:'
         WRITE (*,*) '    |res|_L-infty = ', residual
         WRITE (*,*) '    End of Newton''s iterations.'
         WRITE (*,*)
         EXIT
      ENDIF
      !===================================================================
      !===================================================================

      !------------------------------------------------------------------
      !-------------GENERATION OF THE JACOBIAN MATRIX--------------------
      !-------------OF THE COUPLED EQUATION SYSTEM-----------------------
      !             Jacobian  <--- [(u0.V)_ + (_.V)u0)]  +  K_  +  V_ (ibp)
      !

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    Re = ', Re
#endif
      CALL extract (x0,  u0)
      CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
      CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)

      !------------------------------------------------------------------
      !-------------DIRECT SOLUTION OF THE COUPLED SYSTEM----------------

      CALL par_mumps_master (DIRECT_SOLUTION, 1, Jacobian, 0, dx)
      ! as Newton's method is implemented in NON-INCREMENTAL FORM,
      ! we actually compute x_n even though we call it dx,
      ! dx has then to be computed as dx = x_n - x_n-1
      !
      dx = dx - x0

      !------------------------------------------------------------------
      !-------------LOCA'S STUFF-----------------------------------------
      
      IF ( PRESENT(con_ptr) ) THEN

         WRITE(*,*)
         WRITE(*,*) '    --> CALL to continuation_hook'
#if DEBUG > 2
         WRITE(*,*) '    |x0|_L-infty = ', MAXVAL(ABS(x0))
         WRITE(*,*) '    |dx|_L-infty = ', MAXVAL(ABS(dx))
#endif
         WRITE(*,*)
         !-------------------------------------
         ! WARNING: continuation_hook expects the solution of
         ! J(-x) = + R   and not
         ! J( x) = - R
         !-------------------------------------

         x_hook       = x0
         delta_x_hook = - dx

         continuation_converged = continuation_hook(x_hook, delta_x_hook, &
                                                      con_ptr, Reltol, Abstol);
         dx = - delta_x_hook

         WRITE(*,*) '    done.'
         WRITE(*,*)

      ELSE
         continuation_converged = 1
      ENDIF

      !------------------------------------------------------------------
      !-------------UPDATE SOLUTION VECTOR-------------------------------

      x_vec = x0 + dx

      x0  = x_vec

      dx0 = dx

   ENDDO
   !
   ! end of NEWTON'S ITERATIONS IN NON-INCREMENTAL FORM
   !=============================================================


   !------------------------------------------------------------------
   !-------------UPDATE EVERYTHING------------------------------------

   x0   = x_vec
   xx   = x_vec
   pd%x => xx
   CALL extract (x0,  u0, p0)
   CALL extract (xx,  uu, pp)


   IF ( n <= p_in%nwtn_maxite ) THEN
      num_newt_its = n
   ELSE
      num_newt_its = -1
   ENDIF

   !------------------------------------------------------------------
   ! call case dependent subroutine
   !
   CALL case_newton_postprocess()



   DEALLOCATE( vv, ww, dx, dx0, rhs )


END FUNCTION nonlinear_solver_conwrap

!==============================================================================

END MODULE newton

!------------------------------------------------------------------------------
! qui sotto e` riportato qualcosa che non deve essre perso

!!   FUNCTION nonlinear_solver_conwrap (x_vec, con_ptr, step_num, lambda, delta_s) &
!!      RESULT(num_newt_its) BIND(C, NAME='nonlinear_solver_conwrap')
!!   !
!!   ! Put the call to your nonlinear solver here.
!!   ! Input:
!!   !    x_vec     solution vector
!!   !    con_ptr   pointer to continuation structure, cast to (void *)
!!   !              must be passed to nonlinear solver and then passed
!!   !              to bordering algorithms.
!!   !    step_num  Continuation step number
!!   ! 
!!   ! Output:
!!   !    x_vec     solution vector
!!   ! 
!!   ! Return Value:
!!   !    num_newt_its  Number of Newton iterations needed for
!!   !                  convergence, used to pick next step size.
!!   !                  Negative value means nonlinear solver didn't converge.
!!   ! 
!!   
!!      USE ISO_C_BINDING
!!   
!!      IMPLICIT NONE
!!   
!!      INTERFACE
!!         FUNCTION continuation_hook(x_hook, delta_x_hook, con_ptr, Reltol, Abstol) &
!!            RESULT(continuation_converged) BIND(C, NAME='continuation_hook')
!!            USE ISO_C_BINDING
!!            REAL(KIND=C_DOUBLE)        :: x_hook
!!            REAL(KIND=C_DOUBLE)        :: delta_x_hook
!!            TYPE(C_PTR),         VALUE :: con_ptr
!!            REAL(KIND=C_DOUBLE), VALUE :: Reltol
!!            REAL(KIND=C_DOUBLE), VALUE :: Abstol
!!            INTEGER(KIND=C_INT)        :: continuation_converged
!!         END FUNCTION continuation_hook
!!      END INTERFACE
!!   
!!      ! input variables
!!      REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
!!      TYPE(C_PTR), VALUE                 :: con_ptr
!!      INTEGER(KIND=C_INT), VALUE         :: step_num
!!      REAL(KIND=C_DOUBLE), VALUE         :: lambda
!!      REAL(KIND=C_DOUBLE), VALUE         :: delta_s
!!      ! output variables
!!      INTEGER(KIND=C_INT)                :: num_newt_its
!!   
!!      ! common variables used
!!      ! Jacobian
!!      ! p_in
!!      ! velCmpnnts, np, np_L, Nx
!!      ! x0, u0, p0
!!      ! mm, jj, jj_L, js_D, zero_bvs_D
!!      ! DESINGULARIZE
!!      ! Re
!!   
!!      ! local variables
!!      INTEGER                                   :: n
!!      REAL(KIND=8)                              :: residual
!!      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::     vv, du0
!!      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dx, ww, dx0
!!   
!!   
!!      INTEGER(KIND=C_INT)                       :: continuation_converged=0
!!      REAL(KIND=C_DOUBLE), DIMENSION(Nx)        :: x_hook
!!      REAL(KIND=C_DOUBLE), DIMENSION(Nx)        :: delta_x_hook
!!      REAL(KIND=C_DOUBLE)                       :: Reltol = 1d-3
!!      REAL(KIND=C_DOUBLE)                       :: Abstol = 1d-8
!!      ! (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)
!!   
!!   
!!      ! executable statements
!!   
!!      WRITE(*,*)
!!      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!      WRITE(*,*) '--> CALL to nonlinear_solver_conwrap'
!!      WRITE(*,*)
!!   
!!   
!!      ALLOCATE ( vv (velCmpnnts, np) )
!!      ALLOCATE ( ww (np_L) )
!!      ALLOCATE ( dx (Nx) )
!!      ALLOCATE ( dx0(Nx) )
!!      ALLOCATE ( du0(velCmpnnts, np) )
!!   
!!   
!!      x0 = x_vec
!!      xx = x_vec
!!   
!!   
!!      !============================================
!!      ! start of FIRST STEP IN NON-INCREMENTAL FORM
!!      !
!!      WRITE(*,*) '    First step in NON-INCREMENTAL form'
!!      !
!!      !------------------------------------------------------------------
!!      !-------------GENERATION OF THE RIGHT-HAND SIDE--------------------
!!      !
!!      ! NON-INCREMENTAL FORM
!!      ! rhs <---  (u0 \dot \nabla)u0
!!      !           0
!!      vv = 0
!!   !   CALL extract (x0,  u0, p0)
!!      CALL extract (x0,  u0)
!!      CALL qv_0y01_sp (mm, jj, u0,  vv)
!!   !   CALL qc_ty0_sp_s (ms_2, jjs, iis,  c_2,  vv)  !  cumulative
!!   !   CALL qc_ny0_sp_s (ms_3, jjs, iis, -q_3,  vv)  !  cumulative
!!   
!!      u0(1,:) = volumeForcing(1,1)
!!      u0(2,:) = volumeForcing(2,1)
!!      u0(3,:) = volumeForcing(3,1)
!!      CALL qv_0y0_sp   (mm, jj, u0, 1d0, vv)
!!   
!!      u0(1,:) = volumeForcing(1,2)
!!      u0(2,:) = volumeForcing(2,2)
!!      u0(3,:) = volumeForcing(3,2)
!!      CALL qv_0y0_dR_sp(mm, jj, u0, 1d0, vv)
!!   
!!      ww = 0
!!      CALL collect (vv, ww,  dx) ! here dx is the RHS
!!      !------------------------------------------------------------------
!!      !-------------ENFORCING DIRICHLET BOUNDARY CONDITIONS ON THE RHS---
!!      !
!!      ! NON-INCREMENTAL FORM
!!      ! non-homogeneous boundary conditions
!!      !
!!      CALL Dirichlet_c (np, js_Axis, js_D, bvs_D,  dx)
!!      IF (DESINGULARIZE) dx(Nx) = 0d0
!!      !------------------------------------------------------------------
!!      !-------------GENERATION OF THE JACOBIAN MATRIX--------------------
!!      !-------------OF THE COUPLED EQUATION SYSTEM-----------------------
!!      !             Jacobian  <--- [(u0.V)_ + (_.V)u0)]  +  K_  +  V_ (ibp)
!!      !
!!   !write(*,*) '*check*'
!!   !write(*,*) '    Re = ', Re
!!      CALL extract (x0,  u0)
!!      CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
!!      CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)
!!      !------------------------------------------------------------------
!!      !-------------DIRECT SOLUTION OF THE COUPLED SYSTEM----------------
!!      CALL par_mumps_master (DIRECT_SOLUTION, 1, Jacobian, 0, dx)
!!      ! as the first step is done in NON-INCREMENTAL FORM,
!!      ! we actually compute x_1 even though we call it dx,
!!      ! dx has then to be computed as dx = x_1 - x_0
!!      !
!!      dx = dx - x0
!!      !------------------------------------------------------------------
!!      !-------------UPDATE SOLUTION VECTOR-------------------------------
!!      x_vec = x0 + dx 
!!      x0  = x_vec
!!      dx0 = dx
!!      !
!!      ! end of FIRST STEP IN NON-INCREMENTAL FORM
!!      !============================================
!!   
!!   
!!      WRITE(*,*)
!!   
!!   
!!      !=============================================================
!!      ! start of NEWTON'S ITERATIONS IN BI-INCREMENTAL FORM
!!      !
!!      WRITE(*,*) '    Start of Newton''s iterations'
!!      WRITE(*,*) '    in BI-INCREMENTAL form'
!!      !
!!      DO n = 1, p_in%nwtn_maxite
!!         
!!         WRITE(*,*)
!!         WRITE(*,*) '    n = ', n
!!   
!!         IF (n == p_in%nwtn_maxite) THEN
!!          
!!            WRITE(*,*) '   ************************************************'
!!            WRITE(*,*) '   *Maximum number of Newton''s iterations reached:'
!!            WRITE(*,*) '   *ite_max = ', p_in%nwtn_maxite
!!            WRITE(*,*) '   ************************************************'
!!            WRITE(*,*)
!!         
!!         ENDIF
!!         
!!         CALL extract (x0,  u0)
!!   
!!   !      CALL extract (x0,  u0, p0)
!!   !      CALL vtk_plot_P2 (rr, jj, jj_L, u0, p0, trim(p_in%plot_directory) // 'iteSol.vtk')
!!   
!!         !------------------------------------------------------------------
!!         ! call case dependent subroutine
!!         !
!!         CALL case_newton_iteprocess(n, continuation_converged)
!!   
!!         !------------------------------------------------------------------
!!         !-------------GENERATION OF THE RIGHT-HAND SIDE--------------------
!!         !
!!         ! BI-INCREMENTAL FORM
!!         ! rhs  <---  - (du0 \dot \nabla)du0
!!         !            0
!!         !
!!         vv = 0
!!         CALL extract (dx0,  du0)
!!         CALL qv_0y01_sp (mm, jj, du0,  vv)
!!   
!!         ww = 0
!!   
!!         CALL collect (-vv, ww,  dx) ! here dx is the RHS
!!   
!!         !------------------------------------------------------------------
!!         !-------------ENFORCING DIRICHLET BOUNDARY CONDITIONS ON THE RHS---
!!         !
!!         ! BI-INCREMENTAL FORM
!!         ! differential type boundary conditions
!!         ! (homogeneous if not calling LOCA)
!!         !
!!         CALL extract_Dirichlet_c (np, js_Axis, js_D, x0,  old_bvs_D)
!!         CALL Dirichlet_c_DIFF (np, js_Axis, js_D, bvs_D, old_bvs_D,  dx)
!!   
!!         IF (DESINGULARIZE) dx(Nx) = 0d0
!!   
!!   !write(*,*) '*check*'
!!   !do ii = 1, SIZE(du0, 1)
!!   !   write(*,*) 'MAXdelta_bvs_D = ', MAXVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
!!   !   write(*,*) 'MINdelta_bvs_D = ', MINVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
!!   !enddo
!!   
!!         !------------------------------------------------------------------
!!         !-------------COMPUTE RESIDUAL-------------------------------------
!!   
!!         residual = MAXVAL(ABS(dx))
!!         
!!         WRITE(*,*) '    |res|_L-infty = ', residual
!!   
!!         !===================================================================
!!         !===================================================================
!!         IF (residual < p_in%nwtn_tol .AND. continuation_converged == 1) THEN
!!            x_vec = x0
!!            WRITE (*,*)
!!            WRITE (*,*) '    Residual on the converged solution:'
!!            WRITE (*,*) '    |res|_L-infty = ', residual
!!            WRITE (*,*) '    End of Newton''s iterations.'
!!            WRITE (*,*)
!!            EXIT
!!         ENDIF
!!         !===================================================================
!!         !===================================================================
!!   
!!         !------------------------------------------------------------------
!!         !-------------GENERATION OF THE JACOBIAN MATRIX--------------------
!!         !-------------OF THE COUPLED EQUATION SYSTEM-----------------------
!!         !             Jacobian  <--- [(u0.V)_ + (_.V)u0)]  +  K_  +  V_ (ibp)
!!         !
!!   
!!   !write(*,*) '*check*'
!!   !write(*,*) '    Re = ', Re
!!         CALL extract (x0,  u0)
!!         CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
!!         CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)
!!   
!!         !------------------------------------------------------------------
!!         !-------------DIRECT SOLUTION OF THE COUPLED SYSTEM----------------
!!   
!!         CALL par_mumps_master (DIRECT_SOLUTION, 1, Jacobian, 0, dx)
!!   
!!         !------------------------------------------------------------------
!!         !-------------LOCA'S STUFF-----------------------------------------
!!         
!!         IF ( C_ASSOCIATED(con_ptr) ) THEN
!!   
!!            WRITE(*,*)
!!            WRITE(*,*) '    --> CALL to continuation_hook'
!!            !WRITE(*,*) '    |x0|_L-infty = ', MAXVAL(ABS(x0))
!!            !WRITE(*,*) '    |dx|_L-infty = ', MAXVAL(ABS(dx))
!!            WRITE(*,*)
!!            !-------------------------------------
!!            ! WARNING: continuation_hook expects the solution of
!!            ! J(-x) = + R   and not
!!            ! J( x) = - R
!!            !-------------------------------------
!!   
!!            x_hook       = x0
!!            delta_x_hook = - dx
!!   
!!            continuation_converged = continuation_hook(x_hook(1), delta_x_hook(1), &
!!                                                         con_ptr, Reltol, Abstol);
!!            dx = - delta_x_hook
!!   
!!            WRITE(*,*) '    done.'
!!            WRITE(*,*)
!!   
!!         ELSE
!!            continuation_converged = 1
!!         ENDIF
!!   
!!         !------------------------------------------------------------------
!!         !-------------UPDATE SOLUTION VECTOR-------------------------------
!!   
!!         x_vec = x0 + dx
!!   
!!         x0  = x_vec
!!         dx0 = dx
!!   
!!      ENDDO
!!      !
!!      ! end of NEWTON'S ITERATIONS IN BI-INCREMENTAL FORM
!!      !=============================================================
!!   
!!   
!!      !------------------------------------------------------------------
!!      !-------------UPDATE EVERYTHING------------------------------------
!!   
!!      x0   = x_vec
!!      xx   = x_vec
!!      pd%x = C_LOC(xx)
!!      CALL extract (x0,  u0, p0)
!!      CALL extract (xx,  uu, pp)
!!   
!!   
!!      IF ( n <= p_in%nwtn_maxite ) THEN
!!         num_newt_its = n
!!      ELSE
!!         num_newt_its = -1
!!      ENDIF
!!   
!!      !------------------------------------------------------------------
!!      ! call case dependent subroutine
!!      !
!!      CALL case_newton_postprocess()
!!   
!!   
!!   
!!      DEALLOCATE( vv, ww, dx, dx0, du0 )
!!   
!!   
!!   END FUNCTION nonlinear_solver_conwrap
!!   
!!   !------------------------------------------------------------------------------
