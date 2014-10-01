MODULE loca_wrappers
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 26/8/2013
!
   USE ISO_C_BINDING
   USE loca_parameters
   USE dynamic_structures
   USE sparse_matrix_profiles
   USE sparse_matrix_operations
   USE global_variables
   USE miscellaneous_subroutines
   USE prep_mesh_p1p2_sp
   USE start_sparse_kit
   USE Dirichlet_Neumann
   USE read_input_files
   USE qv_sp
   USE qc_sp_M
   USE qs_L_sp
   USE par_solve_mumps
   USE EigenSolve
   USE axisym_boundary_values
   USE vtk_plot
   USE case_dependent


!------------------------------------------------------------------------------


   IMPLICIT NONE

   !---------------------------------------------------------------------------
   ! LOCA's interfaces and variables
   
   !**********************************
   !* WARNING:
   !* passdown_struct is also defined
   !* in loca_interface.c
   !* Update both files if changes are
   !* needed
   !**********************************
   TYPE, BIND(C) :: passdown_struct
      TYPE(C_PTR)         :: x ! pointer to xx set in main.f90 at line \approx 136
      REAL(KIND=C_DOUBLE) :: reynolds, oscar, romeo, whisky, h, tol
      INTEGER(KIND=C_INT) :: bif_param, param, maxiter, ldz, &
                             num_linear_its, debug
   END TYPE passdown_struct

   INTERFACE
      SUBROUTINE do_loca(pd) BIND(C, NAME='do_loca')
         USE ISO_C_BINDING
         TYPE(C_PTR), VALUE :: pd
      END SUBROUTINE do_loca
   END INTERFACE

   TYPE(C_PTR) :: my_null_ptr
   INTEGER     :: ite_num

   TYPE(passdown_struct), TARGET, BIND(C) :: pd

   !---------------------------------------------------------------------------


   
   TYPE(CSR_MUMPS_Complex_Matrix) :: JmoM ! shifted matrix [J-i*omega*M]
   !TYPE(CSR_MUMPS_Matrix)         :: JmsM ! shifted matrix [J-sigma_loca*M]

   LOGICAL :: Mass_init=.FALSE.
   LOGICAL :: JmsM_init=.FALSE.


   INTEGER :: ii

   

CONTAINS



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LOCA'S WRAPPERS


FUNCTION nonlinear_solver_conwrap (x_vec, con_ptr, step_num, lambda, delta_s) &
   RESULT(num_newt_its) BIND(C, NAME='nonlinear_solver_conwrap')
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

   USE ISO_C_BINDING

   IMPLICIT NONE

   INTERFACE
      FUNCTION continuation_hook(x_hook, delta_x_hook, con_ptr, Reltol, Abstol) &
         RESULT(continuation_converged) BIND(C, NAME='continuation_hook')
         USE ISO_C_BINDING
         REAL(KIND=C_DOUBLE)        :: x_hook
         REAL(KIND=C_DOUBLE)        :: delta_x_hook
         TYPE(C_PTR),         VALUE :: con_ptr
         REAL(KIND=C_DOUBLE), VALUE :: Reltol
         REAL(KIND=C_DOUBLE), VALUE :: Abstol
         INTEGER(KIND=C_INT)        :: continuation_converged
      END FUNCTION continuation_hook
   END INTERFACE

   ! input variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
   TYPE(C_PTR), VALUE                 :: con_ptr
   INTEGER(KIND=C_INT), VALUE         :: step_num
   REAL(KIND=C_DOUBLE), VALUE         :: lambda
   REAL(KIND=C_DOUBLE), VALUE         :: delta_s
   ! output variables
   INTEGER(KIND=C_INT)                :: num_newt_its

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


   INTEGER(KIND=C_INT)                       :: continuation_converged=0
   REAL(KIND=C_DOUBLE), DIMENSION(Nx)        :: x_hook
   REAL(KIND=C_DOUBLE), DIMENSION(Nx)        :: delta_x_hook
   REAL(KIND=C_DOUBLE)                       :: Reltol = 1d-3
   REAL(KIND=C_DOUBLE)                       :: Abstol = 1d-8
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

!write(*,*) '*check*'
!write(*,*) '    Re = ', Re
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
      
      IF ( C_ASSOCIATED(con_ptr) ) THEN

         WRITE(*,*)
         WRITE(*,*) '    --> CALL to continuation_hook'
!write(*,*) '    |x0|_L-infty = ', MAXVAL(ABS(x0))
!write(*,*) '    |dx|_L-infty = ', MAXVAL(ABS(dx))
         WRITE(*,*)
         !-------------------------------------
         ! WARNING: continuation_hook expects the solution of
         ! J(-x) = + R   and not
         ! J( x) = - R
         !-------------------------------------

         x_hook       = x0
         delta_x_hook = - dx

         continuation_converged = continuation_hook(x_hook(1), delta_x_hook(1), &
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
   pd%x = C_LOC(xx)
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

!------------------------------------------------------------------------------

FUNCTION linear_solver_conwrap(xsol, jac_flag, tmp) &
   RESULT(ires) BIND(C, NAME='linear_solver_conwrap')
!
! Put the call to your linear solver here. There are three options
! about the reuse of the preconditioner. It is always safe to
! just solve the matrix from scratch.
! Input:
!    xsol   Right hand side
!    jac_flag  Flag indicating the status of the Jacobian so that
!     preconditioners can be used:
!     NEW_JACOBIAN:  recalculate preconditioner
!     OLD_JACOBIAN:  reuse preconditioner
!     SAME_BUT_UNSCALED_JACOBIAN: Must rescale the matrix and
!           can reuse preconditioner. This happens
!           when the matrix has been recalculated
!           at the same conditions as before.
!    tmp Work space array same length as x, only used for
!     the SAME_BUT_UNSCALED_JACOBIAN option.
!
! Output:
!    xsol   Solution vector
!
! Return Value:
!    Negative value means there was an error.

   USE ISO_C_BINDING

   IMPLICIT NONE
   
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: xsol, tmp
   INTEGER(KIND=C_INT), VALUE         :: jac_flag

   INTEGER(KIND=C_INT)                :: ires
   
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to linear_solver_conwrap'

!write(*,*) '*check*'
!write(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(xsol))
   
   pd%num_linear_its = pd%num_linear_its + 1

   IF (jac_flag == NEW_JACOBIAN .OR. jac_flag == SAME_BUT_UNSCALED_JACOBIAN) THEN
    
      WRITE(*,*) '   *NOT updating Jacobian matrix*'
      ! WRITE(*,*) '    updating Jacobian matrix'

      ! fixme
      !      I think that there is no need to update the Jacobian matrix as it
      !      gets updated in matrix_residual_fill_conwrap!
      !      The NEW_JACOBIAN flag only means that the preconditioner has to be
      !      changed since the jacobian has changed.
      ! fixme

      ! Build the Jacobian matrix and factorize it
      ! WRITE(*,*) '*check*'
      ! WRITE(*,*) '    Re = ', Re
      ! CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
      ! CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)

   ENDIF
   
   CALL par_mumps_master (DIRECT_SOLUTION, 1, Jacobian, 0, xsol)

!write(*,*) '*check*'
!write(*,*) '    |sol|_L-infty = ', MAXVAL(ABS(xsol))
   

   ires = 1

END FUNCTION linear_solver_conwrap

!------------------------------------------------------------------------------

FUNCTION komplex_linear_solver_conwrap(c, d, jac_flag, omega, tmp) & 
         RESULT(ires) BIND(C, NAME='komplex_linear_solver_conwrap')
!
! Put the call to your linear solver here. There are three options
! about the reuse of the preconditioner. It is always safe to
! just solve the matrix from scratch.
! Input:
!    c   Right hand side of real part
!    d   Right hand side of imaginary part
!    jac_flag  Flag indicating the status of the Jacobian so that
!     preconditioners can be used:
!     NEW_JACOBIAN:  recalculate preconditioner
!     OLD_JACOBIAN:  reuse preconditioner
!     OLD_JACOBIAN_DESTROY:   reuse preconditioner,
!              then destroy the preconditioner
!    omega : shift
!
! Output:
!    c, d   Solution vectors
!
! Return Value:
!    Negative value means linear solver didn't converge.
!
   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input/output variables 
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: c, d, tmp
   INTEGER(KIND=C_INT), VALUE         :: jac_flag
   REAL(KIND=C_DOUBLE), INTENT(IN)    :: omega
   INTEGER(KIND=C_INT)                :: ires
   ! local variables
   COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs
   INTEGER       :: i
   LOGICAL, SAVE :: JmoM_init=.FALSE.

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to komplex_linear_solver_conwrap'

   ALLOCATE(rhs(Nx))
   
   rhs = CMPLX(c, d)

!write(*,*) '*check*'
!write(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(DBLE(rhs))), MAXVAL(ABS(AIMAG(rhs)))
   

   IF (jac_flag == NEW_JACOBIAN) THEN

      IF ( .NOT.Mass_init ) THEN
         WRITE(*,*) '    creating Mass matrix'

         ALLOCATE( Mass%i      (SIZE(Jacobian%i))       ); Mass%i       = Jacobian%i
         ALLOCATE( Mass%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass%i_mumps = Jacobian%i_mumps
         ALLOCATE( Mass%j      (SIZE(Jacobian%j))       ); Mass%j       = Jacobian%j
         ALLOCATE( Mass%e      (SIZE(Jacobian%e))       ); Mass%e       = 0d0
         CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Mass)
         ! impose boundary conditions on the Mass Matrix
         CALL Dirichlet_c_M_MASS (np, js_Axis, js_D,  Mass)
         Mass_init = .TRUE.
!write(*,*) '*check*'
!write(*,*) '    |Mass%e|_L-infty = ', MAXVAL(ABS(Mass%e))
      ENDIF


      IF ( .NOT.JmoM_init ) THEN
         WRITE(*,*) '    creating [J-i*omega*M] matrix'

         ! Create the matrix [J-i*omega*M] and store it in position 2 of the MUMPS
         ! array id
         ! JmoM <-- [J-i*omega*M]
         ALLOCATE( JmoM%i      (SIZE(Jacobian%i))       ); JmoM%i       = Jacobian%i
         ALLOCATE( JmoM%i_mumps(SIZE(Jacobian%i_mumps)) ); JmoM%i_mumps = Jacobian%i_mumps
         ALLOCATE( JmoM%j      (SIZE(Jacobian%j))       ); JmoM%j       = Jacobian%j
         ALLOCATE( JmoM%e      (SIZE(Jacobian%e))       ); JmoM%e       = CMPLX(0d0, 0d0, KIND=8)
         CALL par_mumps_master (INITIALIZATION, 2, JmoM, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   2, JmoM, 0)
         JmoM_init=.TRUE.
      ENDIF

write(*,*) '*check*'
write(*,*) '    Re    = ', Re
write(*,*) '    beta  = ', beta
write(*,*) '    omega = ', omega

      IF ( beta /= 0 ) THEN
         WRITE(*,*) '*******************************************************************'
         WRITE(*,*) ''
         WRITE(*,*) 'gran cazzata!'
         WRITE(*,*) 'la continuazione della biforcazione di Hopf non puo` essere fatta'
         WRITE(*,*) 'per beta /= 0 per via di come e` scritta LOCA'
         WRITE(*,*) ''
         WRITE(*,*) '*******************************************************************'
         WRITE(*,*) 'STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
      ENDIF

      WRITE(*,*) '    filling [J-i*omega*M] matrix'
      DO i = 1, SIZE(Jacobian%e)
         ! we can do this as J and M have the same sparsity pattern 
         ! because they have been lazily built
         JmoM%e(i) = CMPLX( Jacobian%e(i), - omega*Mass%e(i) )
      ENDDO

!write(*,*) '*check*'
!write(*,*) '    |[J-i*omega*M]%e|_L-infty = ', MAXVAL(ABS(DBLE(JmoM%e))), MAXVAL(ABS(AIMAG(JmoM%e)))

      CALL par_mumps_master (NUMER_FACTOR, 2, JmoM, 0)

   ENDIF

   CALL par_mumps_master (DIRECT_SOLUTION, 2, JmoM, 0, rhs)
!write(*,*) '*check*'
!write(*,*) '    |sol|_L-infty = ', MAXVAL(ABS(DBLE(rhs))), MAXVAL(ABS(AIMAG(rhs)))


   IF (jac_flag == OLD_JACOBIAN_DESTROY) THEN
    
      IF ( JmoM_init ) THEN
         WRITE(*,*) '    destroying [J-i*omega*M] matrix'

         CALL par_mumps_master (DEALLOCATION, 2, JmoM, 0)
         DEALLOCATE(JmoM%i, JmoM%i_mumps, JmoM%j, JmoM%e)
         JmoM_init=.FALSE.

      ELSE

         WRITE(*,*) '****STRANGE THING HAPPENING: ************'
         WRITE(*,*) '    [J-i*omega*M] was already unallocated'
         WRITE(*,*) '*****************************************'

      ENDIF

   ENDIF

   c = DBLE (rhs)
   d = AIMAG(rhs)

   DEALLOCATE(rhs)

   ires = 0


END FUNCTION komplex_linear_solver_conwrap

!------------------------------------------------------------------------------

SUBROUTINE matrix_residual_fill_conwrap(xsol, rhs, matflag) &
   BIND(C, NAME='matrix_residual_fill_conwrap')
!
! Put the call to your matrix/residual fill routine here.
! Input:
!    xsol      Solution vector
!    matflag   Flag indicating residual (RHS_ONLY), matrix (MATRIX_ONLY),
!              or both (RHS_MATRIX) are requested.
! 
! Output:
!    rhs       Right hand side
! 
! Return Value:
! 
   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: xsol
   INTEGER(KIND=C_INT), VALUE         :: matflag
   ! output variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: rhs
   ! local variables
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vv, ff
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: ww

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to matrix_residual_fill_conwrap'

!write(*,*) '*check*'
!write(*,*) '    |xsol|_L-infty = ', MAXVAL(ABS(xsol))

   IF (matflag == RHS_ONLY .OR. matflag == RHS_MATRIX .OR. matflag == RHS_MATRIX_SAVE) THEN

      WRITE(*,*) '    generating new RHS'

      ALLOCATE ( vv(SIZE(u0,1), SIZE(u0,2)) )
      ALLOCATE ( ff(SIZE(u0,1), SIZE(u0,2)) )
      ALLOCATE ( ww(SIZE(p0)) )
      !------------------------------------------------------------------
      !-------------GENERATION OF THE RIGHT-HAND SIDE--------------------
      !
      ! INCREMENTAL FORM
      ! rhs <-- - (u0 \dot \nabla)u0 + 1/Re lapl{u0} - grad{p0} + f
      !      ~~ - div{u0} ~~
! WRITE(*,*) '*check*'
! WRITE(*,*) '    Re = ', Re
      vv = 0
      CALL extract (xsol,  u0, p0)

      CALL qv_0y01_sp   (mm, jj, u0,              vv) ! (u0 \dot \nabla)u0
      CALL qc_1y1_sp_gg (mm, jj, u0, 1d0/Re,      vv) ! - 1/Re lapl{u0} !!!!!! ibp
      CALL qv_y_10_hybrid_sp (mm, jj, jj_L, -p0,  vv) ! grad{p0} !!!!!!!!!!!!! ibp
      ff(1,:) = volumeForcing(1,1)
      ff(2,:) = volumeForcing(2,1)
      ff(3,:) = volumeForcing(3,1)
      CALL qv_0y0_sp   (mm, jj, ff, 1d0,  vv)
      ff(1,:) = volumeForcing(1,2)
      ff(2,:) = volumeForcing(2,2)
      ff(3,:) = volumeForcing(3,2)
      CALL qv_0y0_dR_sp(mm, jj, ff, 1d0,  vv)
      !do ii = 1, size(vv,2)
      !   write(*,*) vv(:,ii)
      !enddo
      !stop

      ww = 0
      ! no need to use this (furthermore the subroutine commented here
      ! is for 2d cartesian problems)
      !CALL qs_01_hybrid_L_sp (mm, jj, jj_L, u0,   ww) ! div{u0}

      CALL collect (-vv, -ww,  rhs)
      !------------------------------------------------------------------
      !-------------ENFORCING DIRICHLET BOUNDARY CONDITIONS ON THE RHS---
      !
      ! INCREMENTAL FORM
      ! differential type boundary conditions
      !
      CALL extract_Dirichlet_c (np, js_Axis, js_D, xsol,  old_bvs_D)
      CALL Dirichlet_c_DIFF (np, js_Axis, js_D, bvs_D, old_bvs_D, rhs)

      IF (DESINGULARIZE) rhs(Nx) = 0d0

!write(*,*) '*check*'
!do ii = 1, SIZE(u0, 1)
!   write(*,*) 'MAXdelta_bvs_D = ', MAXVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
!   write(*,*) 'MINdelta_bvs_D = ', MINVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
!enddo

      DEALLOCATE( ff, vv, ww )

      !-------------------------------------
      ! WARNING: continuation_hook expects the solution of
      ! J(dx) = - R   and not
      ! J(dx) = + R
      !-------------------------------------
      rhs = -rhs

!write(*,*) '*check*'
!write(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(rhs))
   ENDIF  

   IF (matflag == MATRIX_ONLY .OR. matflag == RHS_MATRIX .OR. matflag == RECOVER_MATRIX ) THEN

      WRITE(*,*) '    filling new Jacobian matrix'

      CALL extract (xsol,  u0) 
      !------------------------------------------------------------------
      !-------------GENERATION OF THE JACOBIAN MATRIX--------------------
      !-------------OF THE COUPLED EQUATION SYSTEM-----------------------
      !             Jacobian  <--- [(u0.V)_ + (_.V)u0)]  +  K_  +  V_ (ibp)

! WRITE(*,*) '*check*'
! WRITE(*,*) '    Re = ', Re
      CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
      CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)

!write(*,*) '*check*'
!write(*,*) '    |Jacobian%e|_L-infty = ', MAXVAL(ABS(Jacobian%e))
   ENDIF

   WRITE(*,*)

END SUBROUTINE matrix_residual_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE mass_matrix_fill_conwrap(xsol, rhs) &
   BIND(C, NAME='mass_matrix_fill_conwrap')
!
! Put the call to your matrix/residual fill routine here.
! Input:
!    xsol      Solution vector (dummy variable)
!    rhs       Right hand side (dummy variable)
!
! Output:
!    Creates mass matrix
!
! Return Value:
!
   USE ISO_C_BINDING

   IMPLICIT NONE

   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: xsol, rhs


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to mass_matrix_fill_conwrap'


   IF ( .NOT.Mass_init ) THEN
      WRITE(*,*) '    creating Mass matrix'

      ALLOCATE( Mass%i      (SIZE(Jacobian%i))       ); Mass%i       = Jacobian%i
      ALLOCATE( Mass%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass%i_mumps = Jacobian%i_mumps
      ALLOCATE( Mass%j      (SIZE(Jacobian%j))       ); Mass%j       = Jacobian%j
      ALLOCATE( Mass%e      (SIZE(Jacobian%e))       ); Mass%e       = 0d0
      CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Mass)
      ! impose boundary conditions on the Mass Matrix
      CALL Dirichlet_c_M_MASS (np, js_Axis, js_D,  Mass)
      Mass_init = .TRUE.

!write(*,*) '*check*'
!write(*,*) '    |Mass%e|_L-infty = ', MAXVAL(ABS(Mass%e))
   ENDIF

END SUBROUTINE mass_matrix_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE matvec_mult_conwrap(xxx, yyy) &
   BIND(C, NAME='matvec_mult_conwrap')
!
! Put the call to your matrix-vector multiply here.
! Input:
!    xxx         Vector of length number of unknowns
!
! Output:
!    yyy         Jacobian matrix times xxx.
!
! Return Value:
!

   USE ISO_C_BINDING

   IMPLICIT NONE
 
   REAL(KIND=C_DOUBLE), DIMENSION(Nx)         :: xxx

   REAL(KIND=C_DOUBLE), DIMENSION(Nx)         :: yyy

   ! local variables
   REAL(KIND=8) :: yi
   INTEGER      :: number_of_rows, i, p

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to matvec_mult_conwrap'

!write(*,*) '*check*'
!write(*,*) '    |x|_L-infty = ', MAXVAL(ABS(xxx))

   CALL dAtimx(yyy, Jacobian%e, Jacobian%j, Jacobian%i, xxx)

!write(*,*) '*check*'
!write(*,*) '    |y|_L-infty = ', MAXVAL(ABS(yyy))

END SUBROUTINE matvec_mult_conwrap

!------------------------------------------------------------------------------

SUBROUTINE mass_matvec_mult_conwrap(xxx, yyy) &
   BIND(C, NAME='mass_matvec_mult_conwrap')
!
! Put the call to your matrix-vector multiply here.
! Input:
!    xxx         Vector of length number of unknowns
!
! Output:
!    yyy         Mass matrix times xxx.
!
! Return Value:
!

   USE ISO_C_BINDING

   IMPLICIT NONE
 
   REAL(KIND=C_DOUBLE), DIMENSION(Nx)         :: xxx

   REAL(KIND=C_DOUBLE), DIMENSION(Nx)         :: yyy

   ! local variables
   REAL(KIND=8) :: yi
   INTEGER      :: number_of_rows, i, p

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to mass_matvec_mult_conwrap'

   IF ( .NOT.Mass_init ) THEN
      WRITE(*,*) '***************************'
      WRITE(*,*) '********* STRANGE *********'
      WRITE(*,*) '    nonexistent Mass matrix'
      WRITE(*,*) '********* STRANGE *********'
      WRITE(*,*) '***************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

!write(*,*) '*check*'
!write(*,*) '    |x|_L-infty = ', MAXVAL(ABS(xxx))

   CALL dAtimx(yyy, Mass%e, Mass%j, Mass%i, xxx)

!write(*,*) '*check*'
!write(*,*) '    |y|_L-infty = ', MAXVAL(ABS(yyy))

END SUBROUTINE mass_matvec_mult_conwrap

!------------------------------------------------------------------------------

SUBROUTINE create_shifted_matrix_conwrap() &
   BIND(C, NAME='create_shifted_matrix_conwrap')
!
! Allocates a sets sparsity pointers for shifted matrix.
! Only used by eigensolver
! 

   USE ISO_C_BINDING

   IMPLICIT NONE
 
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to create_shifted_matrix_conwrap'
   WRITE(*,*)
! 
!   IF ( .NOT.JmsM_init ) THEN
!      ! Create the matrix [J-sigma_loca*M] and store it in position 3 of the MUMPS
!      ! array id
!      ! JmsM <-- [J-sigma_loca*M]
!      ALLOCATE( JmsM%i      (SIZE(Jacobian%i))       ); JmsM%i       = Jacobian%i
!      ALLOCATE( JmsM%i_mumps(SIZE(Jacobian%i_mumps)) ); JmsM%i_mumps = Jacobian%i_mumps
!      ALLOCATE( JmsM%j      (SIZE(Jacobian%j))       ); JmsM%j       = Jacobian%j
!      ALLOCATE( JmsM%e      (SIZE(Jacobian%e))       ); JmsM%e       = 0d0
!      CALL par_mumps_master (INITIALIZATION, 3, JmsM, 0)
!      CALL par_mumps_master (SYMBO_FACTOR,   3, JmsM, 0)
!      JmsM_init=.TRUE.
!   ENDIF
   
END SUBROUTINE create_shifted_matrix_conwrap

!------------------------------------------------------------------------------

SUBROUTINE shifted_matrix_fill_conwrap (sigma_loca) &
   BIND(C, NAME='shifted_matrix_fill_conwrap')
!
! Routine to create shifted matrix  J-sigma_loca*M
! Only called by eigensolver
! 

   USE ISO_C_BINDING

   IMPLICIT NONE
 
   REAL(KIND=C_DOUBLE), INTENT(IN) :: sigma_loca

!   INTEGER :: i
!
!   WRITE(*,*)
!   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!   WRITE(*,*) '--> CALL to shifted_matrix_fill_conwrap'
!   WRITE(*,*)
! 
!   DO i = 1, SIZE(Jacobian%e)
!         ! we can do this as J and M have the same sparsity pattern 
!         ! because they have been lazily built
!         JmsM%e(i) = Jacobian%e(i) - sigma_loca*Mass%e(i)
!   END DO

END SUBROUTINE shifted_matrix_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE shifted_linear_solver_conwrap(xvec, yvec, jac_flag, tol) &
   BIND(C, NAME='shifted_linear_solver_conwrap')
!
! Put the call to your linear solver here. There are three options
! about the reuse of the preconditioner. It is always safe to
! just solve the matrix from scratch.
! Input:
!    xvec       Right hand side
!    jac_flag   Flag indicating the status of the Jacobian so that
!               preconditioners can be used:
!               NEW_JACOBIAN:   recalculate preconditioner
!               OLD_JACOBIAN:   reuse preconditioner
!    tol        linear solver tolerance
!
! Output:
!    yvec       Solution vector
!
! Return Value:
!    Negative value means linear solver didn't converge.
!

   USE ISO_C_BINDING

   IMPLICIT NONE
 
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: xvec
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: yvec
   INTEGER(KIND=C_INT), VALUE         :: jac_flag
   REAL(KIND=C_DOUBLE), VALUE         :: tol

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to shifted_linear_solver_conwrap'
   WRITE(*,*)


!   IF ( jac_flag == NEW_JACOBIAN ) THEN
! 
!      CALL par_mumps_master (NUMER_FACTOR, 3, JmsM, 0)
! 
!   END IF
! 
!   CALL par_mumps_master (DIRECT_SOLUTION, 3, JmsM, 0, xvec)
! 
!   yvec = xvec

END SUBROUTINE shifted_linear_solver_conwrap

!------------------------------------------------------------------------------

SUBROUTINE destroy_shifted_matrix_conwrap() &
   BIND(C, NAME='destroy_shifted_matrix_conwrap')
!
! Deletes memory for shifted matrix.
! Only used by eigensolver
!

   USE ISO_C_BINDING

   IMPLICIT NONE

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to destroy_shifted_matrix_conwrap'
   WRITE(*,*)

!   IF ( JmsM_init ) THEN
! 
!      DEALLOCATE(JmsM%i, JmsM%i_mumps, JmsM%j, JmsM%e)
!      JmsM_init=.FALSE.
!      ! the MUMPS module will take care of deallocating isave(3) the next time
!      ! it will be used
! 
!   ELSE
! 
!      WRITE(*,*) 'JmsM was already unallocated'
! 
!   ENDIF

END SUBROUTINE destroy_shifted_matrix_conwrap

!------------------------------------------------------------------------------

SUBROUTINE assign_parameter_conwrap(param) &
   BIND(C, NAME='assign_parameter_conwrap')
! 
! Put the call to a routine to assign the continuation parameter here.
! Input:
!    param     New value of continuation parameter.
! 
! Output:
! 
! Return Value:
! 

   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   REAL(KIND=C_DOUBLE), VALUE :: param

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to assign_parameter_conwrap'
   WRITE(*,*)

   WRITE(*,*) '    assigning parameter'

   SELECT CASE( pd%param )

      CASE( REYNOLDS )
         WRITE(*,*) '    Reynolds = ', param
         WRITE(*,*) '    oldRe    = ', Re
         Re          = param
         pd%reynolds = param

      CASE( OSCAR )
         WRITE(*,*) '    oscar    = ', param
         WRITE(*,*) '    oldOscar = ', flow_parameters(1)

         CALL case_loca_changeOscar(param)
         flow_parameters(1)  = param
         pd%oscar            = param

      CASE( ROMEO )
         WRITE(*,*) '    romeo    = ', param
         WRITE(*,*) '    oldRomeo = ', flow_parameters(2)

         CALL case_loca_changeRomeo(param)
         flow_parameters(2) = param
         pd%romeo           = param

      CASE( WHISKY )
         WRITE(*,*) '    whisky    = ', param
         WRITE(*,*) '    oldWhisky = ', flow_parameters(3)

         CALL case_loca_changeWhisky(param)
         flow_parameters(3) = param
         pd%whisky          = param

   END SELECT

   WRITE(*,*) ''

END SUBROUTINE assign_parameter_conwrap

!------------------------------------------------------------------------------

SUBROUTINE assign_bif_parameter_conwrap(bif_param) &
   BIND(C, NAME='assign_bif_parameter_conwrap')
! 
! Put the call to a routine to assign the continuation parameter here.
! Input:
!    param     New value of continuation parameter.
! 
! Output:
! 
! Return Value:
! 

   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   REAL(KIND=C_DOUBLE), VALUE :: bif_param


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to assign_bif_parameter_conwrap'
   WRITE(*,*)

   WRITE(*,*) '    assigning bifurcation parameter'

   SELECT CASE( pd%bif_param )

      CASE( REYNOLDS )
         WRITE(*,*) '    Reynolds = ', bif_param
         WRITE(*,*) '    oldRe    = ', Re
         Re          = bif_param
         pd%reynolds = bif_param

      CASE( OSCAR )
         WRITE(*,*) '    oscar    = ', bif_param
         WRITE(*,*) '    oldOscar = ', flow_parameters(1)
       
         CALL case_loca_changeOscar(bif_param)
         flow_parameters(1)  = bif_param
         pd%oscar            = bif_param

      CASE( ROMEO )
         WRITE(*,*) '    romeo    = ', bif_param
         WRITE(*,*) '    oldRomeo = ', flow_parameters(2)

         CALL case_loca_changeRomeo(bif_param)
         flow_parameters(2) = bif_param
         pd%romeo           = bif_param

      CASE( WHISKY )
         WRITE(*,*) '    whisky    = ', bif_param
         WRITE(*,*) '    oldWhisky = ', flow_parameters(3)

         CALL case_loca_changeWhisky(bif_param)
         flow_parameters(3) = bif_param
         pd%whisky          = bif_param

   END SELECT

   WRITE(*,*) ''

END SUBROUTINE assign_bif_parameter_conwrap



! END OF LOCA'S WRAPPERS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





!------------------------------------------------------------------------------
! other LOCA's subroutines which are not wrappers

SUBROUTINE param_output(param) &
   BIND(C, NAME='param_output')
! 
! Print parameter(s) to file
! 
   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   REAL(KIND=C_DOUBLE), VALUE :: param
   ! local variables
   INTEGER           :: fid = 22
   CHARACTER(LEN=50) :: filenm


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to param_output'
   WRITE(*,*)

   SELECT CASE( pd%param )
      CASE( REYNOLDS )
         filenm = './locaOut/reynolds.dat'
         WRITE(*,*) 'writing file: ', trim(filenm)
      CASE( OSCAR )
         filenm = './locaOut/oscar.dat'
         WRITE(*,*) 'writing file: ', trim(filenm)
      CASE( ROMEO )
         filenm = './locaOut/romeo.dat'
         WRITE(*,*) 'writing file: ', trim(filenm)
      CASE( WHISKY )
         filenm = './locaOut/whisky.dat'
         WRITE(*,*) 'writing file: ', trim(filenm)
   END SELECT

   OPEN(UNIT= fid, FILE= trim(filenm), ACCESS= 'APPEND')
   WRITE(fid,*) param
   CLOSE(fid)

   ! print all parameters to file
   filenm = './locaOut/all.dat'
   OPEN(UNIT= fid, FILE= trim(filenm), ACCESS= 'APPEND')
   WRITE(*,*) 'writing file: ', trim(filenm)
   WRITE(fid,*) pd%reynolds, pd%oscar, pd%romeo, pd%whisky

   CLOSE(fid)

   CALL case_loca_paramout()

END SUBROUTINE param_output

!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_loca(x_vec, filenm, filenmLen) &
   BIND(C, NAME='vtk_plot_loca')

   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
   CHARACTER(KIND=C_CHAR)             :: filenm
   INTEGER(KIND=C_INT), VALUE         :: filenmLen
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_plot
   REAL(KIND=8), DIMENSION(np_L)          :: p_plot


   IF ( p_in%write_plots_flag ) THEN

      CALL extract(x_vec, u_plot, p_plot)

      CALL vtk_plot_P2 (rr, jj, jj_L, u_plot, p_plot, &
        trim(p_in%plot_directory)//'locaContSolution'//filenm(1:filenmLen)//'.vtk' )

   ENDIF


END SUBROUTINE vtk_plot_loca

!------------------------------------------------------------------------------

SUBROUTINE compute_eigen(x_vec, filenm, filenmLen, shiftIm) &
   BIND(C, NAME='compute_eigen')
!
! Compute eigenvalues and eigenvectors and save them to file
!
   USE ISO_C_BINDING

   IMPLICIT NONE

   ! input variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
   CHARACTER(KIND=C_CHAR)             :: filenm
   INTEGER(KIND=C_INT), VALUE         :: filenmLen
   REAL(KIND=C_DOUBLE), VALUE         :: shiftIm

   ! local variables
   TYPE(CSR_MUMPS_Complex_Matrix)     :: Lns_cmplx, Mass_cmplx

   LOGICAL, DIMENSION(velCmpnnts, number_of_sides) :: Dir_eigen
   TYPE(dyn_int_line), DIMENSION(velCmpnnts)       :: js_D_eigen
   LOGICAL                                         :: DESINGULARIZE_eigen

   REAL(KIND=8), DIMENSION(SIZE(Jacobian%e)) :: Jacobian_save

   INTEGER :: i, k, eigen_plotNum

   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: directEigenvalues,  adjointEigenvalues
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: directEigenvectors, adjointEigenvectors
   REAL(KIND=8),    DIMENSION(:),   ALLOCATABLE :: structuralsens

   INTEGER :: statusMsg

   CHARACTER(LEN=128) :: shiftName, & ! used to insert the shift in file names
                         shiftNameRe, shiftNameIm
   INTEGER            :: shiftNameTemp

   REAL(KIND=8) :: dataIdentifier = 313d0 ! Donald Duck's plate number!
   INTEGER      :: shiftsNumber = 1
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: shifts

   COMPLEX(KIND=8), DIMENSION(Nx) :: tmpEigen1, tmpEigen2 ! used to check residuals

   !-------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to compute_eigen'
   WRITE(*,*) ''
   !-------------------------------------------------


   ! if "present" 'shiftIm' override read imaginary part of complex shift
   IF ( ABS(shiftIm) > 1d-8 ) THEN
      WRITE(*,*) '    overriding read complex shift'
      WRITE(*,*) '    LOCA suggests: ', shiftIM
      p_in%eigen_sigma = CMPLX(DBLE(p_in%eigen_sigma), shiftIm, KIND=8)

   ELSEIF (       DBLE(p_in%eigen_sigma) == dataIdentifier &
           .AND. AIMAG(p_in%eigen_sigma) == dataIdentifier) THEN
      ! read shifts from file
      WRITE(*,*) '--> Reading shifts from file shifts.in'

      OPEN( UNIT = 20, FILE = 'shifts.in' )
      READ(20,*) shiftsNumber
      ALLOCATE(shifts(shiftsNumber))
      DO i = 1, shiftsNumber
         READ(20,*) shifts(i)
      ENDDO
      CLOSE(20)
#if DEBUG > 0
      WRITE(*,*) '    read ', shiftsNumber, ' shifts'
      WRITE(*,*) '    Done.'
#endif
   ENDIF

   !----------------------------------------------------
   ! PREPARE COMPLEX MATRICES FOR THE EIGENVALUE PROBLEM
   !
   Jacobian_save = Jacobian%e

   ! (1)
   ! prepare boundary conditions
   !
#if DEBUG > 2
   WRITE(*,*) '*check*'
   WRITE(*,*) '    number_of_sides = ', number_of_sides
#endif

   IF ( p_in%eigen_BC == 1 ) THEN
      ! homogeneous Dirichlet on every border
#if DEBUG > 0
      WRITE(*,*) '    boundary conditions: zero velocity on every border'
#endif
      Dir_eigen = .TRUE.
      DESINGULARIZE_eigen = .TRUE.
   ELSEIF ( p_in%eigen_BC == 2 ) THEN
      ! same BCs as for the base flow but homogeneous
      ! Dirichlet where the base flow has nonhomogeneous Dirichlet
#if DEBUG > 0
      WRITE(*,*) '    boundary conditions: same as base flow'
#endif
      Dir_eigen = Dir
      DESINGULARIZE_eigen = DESINGULARIZE
   ELSE
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** p_in % eigen_BC               ***'
      WRITE(*,*) '*** set to: ', p_in%eigen_BC
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

   DO k = 1, velCmpnnts
      CALL Dirichlet_nodes_gen (jjs, sides, Dir_eigen(k,:),  js_D_eigen(k)%DIL)
   ENDDO

   ! (2)
   ! allocate the complex matrices
   !
#if DEBUG > 0
   WRITE(*,*) '--> Creating complex matrices'
#endif

   ALLOCATE( Lns_cmplx%i      (SIZE(Jacobian%i))       ); Lns_cmplx%i       = Jacobian%i
   ALLOCATE( Lns_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) ); Lns_cmplx%i_mumps = Jacobian%i_mumps
   ALLOCATE( Lns_cmplx%j      (SIZE(Jacobian%j))       ); Lns_cmplx%j       = Jacobian%j
   ALLOCATE( Lns_cmplx%e      (SIZE(Jacobian%e))       ); Lns_cmplx%e       = CMPLX(0d0, 0d0,KIND=8)

   ALLOCATE( Mass_cmplx%i      (SIZE(Jacobian%i))       ); Mass_cmplx%i       = Jacobian%i
   ALLOCATE( Mass_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass_cmplx%i_mumps = Jacobian%i_mumps
   ALLOCATE( Mass_cmplx%j      (SIZE(Jacobian%j))       ); Mass_cmplx%j       = Jacobian%j
   ALLOCATE( Mass_cmplx%e      (SIZE(Jacobian%e))       ); Mass_cmplx%e       = CMPLX(0d0, 0d0, KIND=8)

   ! (3)
   ! fill the Lns matrix
   ! NOTE: the sign of the Lns matrix has to be changed
   ! since the problem we are solving is:
   ! lambda*Mass*x = -Lns*x
   !
#if DEBUG > 2
   WRITE(*,*) '*check*'
   WRITE(*,*) '    Re   = ', Re
   WRITE(*,*) '    beta = ', beta
#endif
   CALL extract(x_vec, u0)
   CALL qc_1y1_sp_gg_3d_M  (mm, jj,                  1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
   CALL qc_oseen2y_sp_3d_M (mm, jj,                   u0,    beta,  Lns_cmplx) ! + linearized terms
   CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,           -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
   CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,           -1d0,    beta,  Lns_cmplx) ! - velocity divergence
   CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_eigen, 1d0,           Lns_cmplx) ! Dirichlet BCs

   Lns_cmplx%e = - Lns_cmplx%e

   IF (DESINGULARIZE_eigen) THEN
      ! column
      WHERE (Lns_cmplx%j == Nx)
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
      ENDWHERE
      ! row
      DO i = Lns_cmplx%i(Nx), Lns_cmplx%i(Nx + 1) - 1
         Lns_cmplx%e(i) = CMPLX(0d0,0d0,KIND=8)
         IF (Lns_cmplx%j(i) == Nx) Lns_cmplx%e(i) = CMPLX(1d0,0d0,KIND=8)
      ENDDO
   ENDIF

   ! (4)
   ! create the real Mass matrix with qc_0y0_zero_sp_M and
   ! convert the elements of the real Mass matrix (Mass) to Complex type
   ! copying them into the new CRS_MUMPS_Complex_Matrix Mass_cmplx
   !
   ! EXPLAINATION: we use the Jacobian matrix to save memory in case the Mass
   !               matrix hasn't been created yet
   !
   Jacobian%e = 0d0
   CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Jacobian)

   ! impose boundary conditions on the Mass Matrix
   CALL Dirichlet_rc_M (np, js_Axis, js_D_eigen, 0d0,  Jacobian)

   Mass_cmplx%e = CMPLX(Jacobian%e, 0d0, KIND=8)


!+++
!write(*,*) 'exporting matrices'
!open(unit=20,file='Jacobian.dat',form='formatted')
!do i = 1,size(Lns_cmplx%e)-1
!   write(20,*) Lns_cmplx%i_mumps(i), Lns_cmplx%j(i), dble(Lns_cmplx%e(i))
!enddo
!close(20)
!open(unit=20,file='Mass.dat',form='formatted')
!do i = 1,size(Mass_cmplx%e)-1
!   write(20,*) Mass_cmplx%i_mumps(i), Mass_cmplx%j(i), dble(Mass_cmplx%e(i))
!enddo
!close(20)
!write(*,*) 'done'
!return
!+++



   DO i = 1, shiftsNumber

      IF ( shiftsNumber /= 1 ) THEN
         p_in%eigen_sigma = shifts(i)
         WRITE(*,*) '    Shift number: ', i
      ENDIF

      ! Create very elaborate and stupid file name
      shiftNameTemp = NINT(ABS(DBLE(p_in%eigen_sigma))*1e3)
      CALL intToChar6 (shiftNameTemp,  shiftNameRe)
      IF ( DBLE(p_in%eigen_sigma) >= 0d0 ) THEN
         WRITE(shiftNameRe,'(a7)') '+'//trim(shiftNameRe)
      ELSE
         WRITE(shiftNameRe,'(a7)') '-'//trim(shiftNameRe)
      ENDIF
      shiftNameTemp = NINT(ABS(AIMAG(p_in%eigen_sigma))*1e3)
      CALL intToChar6 (shiftNameTemp,  shiftNameIm)
      IF ( AIMAG(p_in%eigen_sigma) >= 0d0 ) THEN
         WRITE(shiftNameIm,'(a7)') '+'//trim(shiftNameIm)
      ELSE
         WRITE(shiftNameIm,'(a7)') '-'//trim(shiftNameIm)
      ENDIF

      WRITE(shiftName, '(a26)') '-Shift_Re'//trim(shiftNameRe)//'_Im'//trim(shiftNameIm)




      !=====================================
      ! COMPUTE EIGENVALUES AND EIGENVECTORS
      !
      IF (      p_in%eigen_directAdjoint_flag == 1 &
           .OR. p_in%eigen_directAdjoint_flag == 3 &
           .OR. p_in%eigen_compute_structSens_flag ) THEN
         !
         ! direct eigenvalues and eigenvectors
         !
         WRITE(*,*)
         WRITE(*,*) '    Computing eigensolutions for the DIRECT problem'
         WRITE(*,*)

         CALL eigensComplexShiftInvert(p_in%eigen_nev,   &
                                       p_in%eigen_maxit, &
                                       p_in%eigen_tol,   &
                                       p_in%eigen_sigma, &
                          Lns_cmplx, Mass_cmplx, 1,  statusMsg, directEigenvalues, directEigenvectors)

         IF ( statusMsg .EQ. 0 ) THEN
            !----------------
            ! CHECK RESIDUALS
            !
            WRITE(*,*)
            WRITE(*,*) '    Check residuals on the direct problem'
            DO k = 1, SIZE(directEigenvalues)
            
               CALL zAtimx (tmpEigen1, Mass_cmplx%e, Mass_cmplx%j, Mass_cmplx%i, directEigenvectors(:,k))
               CALL zAtimx (tmpEigen2, Lns_cmplx%e,  Lns_cmplx%j,  Lns_cmplx%i,  directEigenvectors(:,k))
            
               WRITE(*,*) '    eig', k, MAXVAL(ABS( tmpEigen1*directEigenvalues(k) - tmpEigen2 ))
            
            ENDDO
            WRITE(*,*)

            !-------------
            ! SAVE RESULTS
            !
            CALL Save_eigenvalues  (directEigenvalues,  &
                               trim(p_in%eigen_output_directory)// &
                                'directEigenvalues'//filenm(1:filenmLen)//trim(shiftName)//'.dat')

!            CALL Save_eigenvectors (directEigenvectors, &
!                               trim(p_in%eigen_output_directory)// &
!                               'directEigenvectors'//filenm(1:filenmLen)//trim(shiftName)//'.dat')

            !------------------
            ! PLOT EIGENVECTORS
            !
            IF ( p_in%write_plots_flag ) THEN

               IF ( SIZE(directEigenvectors,2) .LT. p_in%eigen_plotNumber ) THEN
                  eigen_plotNum = SIZE(directEigenvectors,2)
               ELSE
                  eigen_plotNum = p_in%eigen_plotNumber
               ENDIF

               CALL vtk_plot_eigenvectors (rr, jj,  DBLE(directEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'directEigenvectorsRe'//filenm(1:filenmLen)//trim(shiftName)//'.vtk')
               CALL vtk_plot_eigenvectors (rr, jj, AIMAG(directEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'directEigenvectorsIm'//filenm(1:filenmLen)//trim(shiftName)//'.vtk')
            ENDIF

         ENDIF

      ENDIF

      IF (      p_in%eigen_directAdjoint_flag == 2 &
           .OR. p_in%eigen_directAdjoint_flag == 3 &
           .OR. p_in%eigen_compute_structSens_flag ) THEN
         !
         ! adjoint eigenvalues and eigenvectors
         !
         WRITE(*,*)
         WRITE(*,*) '    Computing eigensolutions for the ADJOINT problem'
         WRITE(*,*)

         CALL eigensComplexShiftInvert(p_in%eigen_nev,   &
                                       p_in%eigen_maxit, &
                                       p_in%eigen_tol,   &
                                       p_in%eigen_sigma, &
                          Lns_cmplx, Mass_cmplx, 2,  statusMsg, adjointEigenvalues, adjointEigenvectors)

         IF ( statusMsg .EQ. 0 ) THEN
            !----------------
            ! CHECK RESIDUALS
            !
            WRITE(*,*)
            WRITE(*,*) '    Check residuals on the adjoint problem'
            DO k = 1, SIZE(adjointEigenvalues)
            
               CALL zAtimx_T (tmpEigen1, Mass_cmplx%e, Mass_cmplx%j, Mass_cmplx%i, adjointEigenvectors(:,k))
               CALL zAtimx_T (tmpEigen2, Lns_cmplx%e,  Lns_cmplx%j,  Lns_cmplx%i,  adjointEigenvectors(:,k))
            
               WRITE(*,*) '    eig', k, MAXVAL(ABS( tmpEigen1*adjointEigenvalues(k) - tmpEigen2 ))
            
            ENDDO
            WRITE(*,*)

            !-------------
            ! SAVE RESULTS
            !
            CALL Save_eigenvalues  (adjointEigenvalues,  &
                                trim(p_in%eigen_output_directory)// &
                                 'adjointEigenvalues'//filenm(1:filenmLen)//trim(shiftName)//'.dat')

!            CALL Save_eigenvectors (adjointEigenvectors, &
!                                trim(p_in%eigen_output_directory)// &
!                                'adjointEigenvectors'//filenm(1:filenmLen)//trim(shiftName)//'.dat')

            !------------------
            ! PLOT EIGENVECTORS
            !
            IF ( p_in%write_plots_flag ) THEN

               IF ( SIZE(adjointEigenvectors,2) .LT. p_in%eigen_plotNumber ) THEN
                  eigen_plotNum = SIZE(adjointEigenvectors,2)
               ELSE
                  eigen_plotNum = p_in%eigen_plotNumber
               ENDIF

               CALL vtk_plot_eigenvectors (rr, jj,  DBLE(adjointEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'adjointEigenvectorsRe'//filenm(1:filenmLen)//trim(shiftName)//'.vtk')
               CALL vtk_plot_eigenvectors (rr, jj, AIMAG(adjointEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'adjointEigenvectorsIm'//filenm(1:filenmLen)//trim(shiftName)//'.vtk')
            ENDIF

         ENDIF

      ENDIF

      IF ( p_in%eigen_compute_structSens_flag ) THEN
         !
         ! structural sensitivity
         !
         CALL structuralSensitivity(Mass_cmplx, adjointEigenvectors(:,p_in%structSens_eigenNumber), &
                                                 directEigenvectors(:,p_in%structSens_eigenNumber), &
                                                 velCmpnnts, np,  structuralsens)
         CALL vtk_plot_scalar_P2 (rr, jj, structuralsens, &
                           './structSensOut/'//'structuralSensitivity'//filenm(1:filenmLen)//'.vtk')

      ENDIF

      !----------------------------------------
      ! DEALLOCATE EIGENVALUES AND EIGENVECTORS
      !
      IF ( p_in%eigen_directAdjoint_flag == 1 .AND. .NOT.p_in%eigen_compute_structSens_flag .AND. statusMsg == 0 ) THEN
         DEALLOCATE(  directEigenvalues,  directEigenvectors )
      ENDIF
      IF ( p_in%eigen_directAdjoint_flag == 2 .AND. .NOT.p_in%eigen_compute_structSens_flag .AND. statusMsg == 0 ) THEN
         DEALLOCATE( adjointEigenvalues, adjointEigenvectors )
      ENDIF
      IF ( p_in%eigen_directAdjoint_flag == 3 .OR.  p_in%eigen_compute_structSens_flag .AND. statusMsg == 0 ) THEN
         DEALLOCATE(  directEigenvalues,  directEigenvectors )
         DEALLOCATE( adjointEigenvalues, adjointEigenvectors )
      ENDIF
      IF ( p_in%eigen_compute_structSens_flag ) THEN
         DEALLOCATE( structuralsens )
      ENDIF
   
   ENDDO ! shiftsNumber



   !------------------------------------------
   ! RESTORE Jacobian MATRIX FOR THE BASE FLOW
   !
   Jacobian%e = Jacobian_save


   !--------------------------------
   ! DEALLOCATE THE COMPLEX MATRICES
   !
   DEALLOCATE( Lns_cmplx%i,  Lns_cmplx%i_mumps,  Lns_cmplx%j,  Lns_cmplx%e )
   DEALLOCATE( Mass_cmplx%i, Mass_cmplx%i_mumps, Mass_cmplx%j, Mass_cmplx%e )


END SUBROUTINE compute_eigen

!-----------------------------------------------------------------------------

SUBROUTINE compute_structSens (Jacobian)

   IMPLICIT NONE
   ! input variables
   TYPE(CSR_MUMPS_Matrix) :: Jacobian ! only used for its sparsity pattern to create the Mass_cmplx matrix

   ! local variables
   TYPE(CSR_MUMPS_Complex_Matrix) :: Mass_cmplx
   COMPLEX(KIND=8), DIMENSION(Nx) :: directEigenvector, adjointEigenvector
   REAL(KIND=8),    DIMENSION(:), ALLOCATABLE :: structuralsens
   REAL(KIND=8),    DIMENSION(Nx) :: tempEigenvectorRe, tempEigenvectorIm ! needed because the subroutine
                                                                          ! read_eigenvector needs them

   LOGICAL, DIMENSION(velCmpnnts, number_of_sides) :: Dir_eigen
   TYPE(dyn_int_line), DIMENSION(velCmpnnts)       :: js_D_eigen
   LOGICAL                                         :: DESINGULARIZE_eigen

   INTEGER :: k

!----------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Start of compute_structSens'
   WRITE(*,*) ''
!----------------------------------------------------


!-------------------------------
! CREATE THE complex MASS MATRIX

   ! (1)
   ! prepare HOMOGENEOUS boundary conditions on ALL BORDERS!
   !
   Dir_eigen = .TRUE.
   DESINGULARIZE_eigen = .TRUE.
   DO k = 1, velCmpnnts
      CALL Dirichlet_nodes_gen (jjs, sides, Dir_eigen(k,:),  js_D_eigen(k)%DIL)
   ENDDO

   ! (2)
   ! allocate Mass_cmplx
   !
   ALLOCATE( Mass_cmplx%i      (SIZE(Jacobian%i))       )
   ALLOCATE( Mass_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) )
   ALLOCATE( Mass_cmplx%j      (SIZE(Jacobian%j))       )
   ALLOCATE( Mass_cmplx%e      (SIZE(Jacobian%e))       )

   ! (3)
   ! create the real Mass matrix with qc_0y0_zero_sp_M and
   ! convert the elements of the real Mass matrix (Mass) to Complex type
   ! copying them into the new CRS_MUMPS_Complex_Matrix Mass_cmplx
   !
   ! EXPLAINATION: we use the Jacobian matrix to save memory in case the Mass
   !               matrix hasn't been created (yet)
   !
   Jacobian%e = 0
   CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Jacobian)

   ! impose HOMOGENEOUS boundary conditions on the Mass Matrix
   CALL Dirichlet_c_M_MASS (np, js_Axis, js_D_eigen,  Jacobian)

   Mass_cmplx%i       = Jacobian%i
   Mass_cmplx%i_mumps = Jacobian%i_mumps
   Mass_cmplx%j       = Jacobian%j
   Mass_cmplx%e       = CMPLX(Jacobian%e, 0d0, KIND=8)

!----------------------
! READ THE EIGENVECTORS

   CALL read_eigenvector (Nx, p_in%structSens_eigenNumber,             &
                          trim(p_in%structSens_directEigen_name),      &
                          LEN(trim(p_in%structSens_directEigen_name)), &
                          tempEigenvectorRe, tempEigenvectorIm)

   directEigenvector  = CMPLX(tempEigenvectorRe, tempEigenvectorIm, KIND=8)



   CALL read_eigenvector (Nx, p_in%structSens_eigenNumber,              &
                          trim(p_in%structSens_adjointEigen_name),      &
                          LEN(trim(p_in%structSens_adjointEigen_name)), &
                          tempEigenvectorRe, tempEigenvectorIm)

   adjointEigenvector = CMPLX(tempEigenvectorRe, tempEigenvectorIm, KIND=8)

!-------------------------------
! COMPUTE STRUCTURAL SENSITIVITY

   CALL structuralSensitivity(Mass_cmplx, adjointEigenvector, &
                                           directEigenvector, &
                                           velCmpnnts, np,  structuralsens)

!-------------
! PLOT RESULTS

   CALL vtk_plot_scalar_P2 (rr, jj, structuralsens, &
                     './structSensOut/'//'structuralSensitivity'//'SteadyState'//'.vtk')

!---------------------
! DEALLOCATE VARIABLES

   DEALLOCATE( structuralsens )
   DEALLOCATE( Mass_cmplx%i, Mass_cmplx%i_mumps, Mass_cmplx%j, Mass_cmplx%e )

END SUBROUTINE compute_structSens

!==============================================================================
!==============================================================================
!==============================================================================

END MODULE loca_wrappers

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
