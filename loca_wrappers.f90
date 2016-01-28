MODULE loca_wrappers
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/10/2014
!
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
   USE restart_io
   !
   USE loca_types
   USE loca_pd


!------------------------------------------------------------------------------


   IMPLICIT NONE

   TYPE(CSR_MUMPS_Complex_Matrix) :: JmoM ! shifted matrix [J-i*omega*M]
   !TYPE(CSR_MUMPS_Matrix)         :: JmsM ! shifted matrix [J-sigma_loca*M] ! USELESS

   !LOGICAL :: Mass_init=.FALSE. ! moved to global_variables
   !LOGICAL :: JmsM_init=.FALSE. ! USELESS

   INTEGER :: ii

   

CONTAINS



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! LOCA'S WRAPPERS

FUNCTION linear_solver_conwrap(xsol, jac_flag, tmp) &
   RESULT(ires)
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

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(Nx) :: xsol, tmp
   INTEGER                     :: jac_flag

   INTEGER                     :: ires
   
#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to linear_solver_conwrap'
#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(xsol))
#endif
   
   pd%num_linear_its = pd%num_linear_its + 1

   IF (jac_flag == NEW_JACOBIAN .OR. jac_flag == SAME_BUT_UNSCALED_JACOBIAN) THEN
    
#if DEBUG > 0
      WRITE(*,*) '   *NOT updating Jacobian matrix*'
      ! WRITE(*,*) '    updating Jacobian matrix'
#endif
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

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |sol|_L-infty = ', MAXVAL(ABS(xsol))
#endif
   

   ires = 1

END FUNCTION linear_solver_conwrap

!------------------------------------------------------------------------------

FUNCTION komplex_linear_solver_conwrap(c, d, jac_flag, omega, tmp) & 
         RESULT(ires)
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
   IMPLICIT NONE
   ! input/output variables 
   REAL(KIND=8), DIMENSION(Nx) :: c, d, tmp
   INTEGER,      INTENT(IN)    :: jac_flag
   REAL(KIND=8), INTENT(IN)    :: omega
   INTEGER                     :: ires
   ! local variables
   COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs
   INTEGER       :: i
   LOGICAL, SAVE :: JmoM_init=.FALSE.

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to komplex_linear_solver_conwrap'
#endif

   ALLOCATE(rhs(Nx))
   
   rhs = CMPLX(c, d, KIND=8)

!#if DEBUG > 3
!   CALL save_vector(c,'rhs_r.vector')
!   CALL save_vector(d,'rhs_i.vector')
!#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(DBLE(rhs))), MAXVAL(ABS(AIMAG(rhs)))
   WRITE(*,*) '    rhs(1)        = ',  DBLE(rhs(1)),  AIMAG(rhs(1))
   WRITE(*,*) '    rhs(end)      = ',  DBLE(rhs(Nx)), AIMAG(rhs(Nx))
#endif
   

   IF (jac_flag == NEW_JACOBIAN) THEN

      IF ( .NOT.Mass_init ) THEN
         WRITE(*,*) 'should not get here'
         STOP
      ENDIF


      IF ( .NOT.JmoM_init ) THEN
#if DEBUG > 0
         WRITE(*,*) '    creating [J-i*omega*M] matrix'
#endif

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

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    Re    = ', Re
      WRITE(*,*) '    beta  = ', beta
      WRITE(*,*) '    omega = ', omega
#endif

#if DEBUG > 0
      WRITE(*,*) '    filling [J-i*omega*M] matrix'
#endif

!#if DEBUG > 3
!   CALL save_vector( Jacobian%e,'Jacobian.vector')
!#endif
      DO i = 1, SIZE(Jacobian%e)
         ! we can do this as J and M have the same sparsity pattern 
         ! because they have been lazily built
         JmoM%e(i) = CMPLX( Jacobian%e(i), - omega*Mass%e(i), KIND=8 )
      ENDDO

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    |[J-i*omega*M]%e|_L-infty = ', MAXVAL(ABS(DBLE(JmoM%e))), MAXVAL(ABS(AIMAG(JmoM%e)))
#endif

      CALL par_mumps_master (NUMER_FACTOR, 2, JmoM, 0)

   ENDIF

   IF (jac_flag == NEW_LNS_BETA) THEN

      IF ( .NOT.Mass_init ) THEN
         WRITE(*,*) 'should not get here 1'
         STOP
      ENDIF
      IF ( .NOT.Lns_cmplx_init ) THEN
         WRITE(*,*) 'should not get here 2'
         STOP
      ENDIF

      IF ( .NOT.JmoM_init ) THEN
#if DEBUG > 0
         WRITE(*,*) '    allocating [Lns-i*omega*M] matrix'
#endif

         ! Create the matrix [Lns-i*omega*M] and store it in position 2 of the MUMPS
         ! array id
         ! JmoM <-- [Lns-i*omega*M]
         ALLOCATE( JmoM%i      (SIZE(Jacobian%i))       ); JmoM%i       = Jacobian%i
         ALLOCATE( JmoM%i_mumps(SIZE(Jacobian%i_mumps)) ); JmoM%i_mumps = Jacobian%i_mumps
         ALLOCATE( JmoM%j      (SIZE(Jacobian%j))       ); JmoM%j       = Jacobian%j
         ALLOCATE( JmoM%e      (SIZE(Jacobian%e))       ); JmoM%e       = CMPLX(0d0, 0d0, KIND=8)
         CALL par_mumps_master (INITIALIZATION, 2, JmoM, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   2, JmoM, 0)
         JmoM_init=.TRUE.
      ENDIF

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    omega = ', omega
#endif

#if DEBUG > 0
      WRITE(*,*) '    filling [Lns-i*omega*M] matrix'
#endif

      !JmoM%e = Lns_cmplx%e - CMPLX( 0d0, omega*Mass%e )
      DO i = 1, SIZE(Lns_cmplx%e)
         ! we can do this as Lns and M have the same sparsity pattern 
         ! because they have been lazily built
         JmoM%e(i) = Lns_cmplx%e(i) - CMPLX( 0d0, omega, KIND=8 )*Mass%e(i)
      ENDDO

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    |[Lns-i*omega*M]%e|_L-infty = ', MAXVAL(ABS(DBLE(JmoM%e))), MAXVAL(ABS(AIMAG(JmoM%e)))
#endif

      CALL par_mumps_master (NUMER_FACTOR, 2, JmoM, 0)

   ENDIF

!#if DEBUG > 3
!   CALL save_vector( DBLE(JmoM%e),'JmoM_r.vector')
!   CALL save_vector(AIMAG(JmoM%e),'JmoM_i.vector')
!#endif
!#if DEBUG > 3
!   CALL save_vector( Mass%e,'Mass.vector')
!#endif

   CALL par_mumps_master (DIRECT_SOLUTION, 2, JmoM, 0, rhs)


#if DEBUG > 0
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |sol|_L-infty = ', MAXVAL(ABS(DBLE(rhs))), MAXVAL(ABS(AIMAG(rhs)))
   WRITE(*,*) '    sol(1)        = ',  DBLE(rhs(1)),  AIMAG(rhs(1))
   WRITE(*,*) '    sol(end)      = ',  DBLE(rhs(Nx)), AIMAG(rhs(Nx))
#endif


   IF (jac_flag == OLD_JACOBIAN_DESTROY) THEN
    
      IF ( JmoM_init ) THEN
#if DEBUG > 0
         WRITE(*,*) '    destroying [J-i*omega*M] matrix'
#endif

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

!#if DEBUG > 3
!   CALL save_vector(c,'sol_r.vector')
!   CALL save_vector(d,'sol_i.vector')
!#endif

   DEALLOCATE(rhs)

   ires = 0

END FUNCTION komplex_linear_solver_conwrap

!------------------------------------------------------------------------------

SUBROUTINE matrix_residual_fill_conwrap(xsol, rhs, matflag)
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
   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: xsol
   INTEGER                     :: matflag
   ! output variables
   REAL(KIND=8), DIMENSION(Nx) :: rhs
   ! local variables
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vv, ff
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: ww

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to matrix_residual_fill_conwrap'
#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |xsol|_L-infty = ', MAXVAL(ABS(xsol))
#endif

   IF (matflag == RHS_ONLY .OR. matflag == RHS_MATRIX .OR. matflag == RHS_MATRIX_SAVE) THEN

#if DEBUG > 0
      WRITE(*,*) '    generating new RHS'
#endif

      ALLOCATE ( vv(SIZE(u0,1), SIZE(u0,2)) )
      ALLOCATE ( ff(SIZE(u0,1), SIZE(u0,2)) )
      ALLOCATE ( ww(SIZE(p0)) )
      !------------------------------------------------------------------
      !-------------GENERATION OF THE RIGHT-HAND SIDE--------------------
      !
      ! INCREMENTAL FORM
      ! rhs <-- - (u0 \dot \nabla)u0 + 1/Re lapl{u0} - grad{p0} + f
      !      ~~ - div{u0} ~~
#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    Re = ', Re
#endif

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

      ww = 0d0
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

#if DEBUG > 2
      WRITE(*,*) '*check*'
      DO ii = 1, SIZE(u0, 1)
         WRITE(*,*) 'MAXdelta_bvs_D = ', MAXVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
         WRITE(*,*) 'MINdelta_bvs_D = ', MINVAL(bvs_D(ii)%DRL-old_bvs_D(ii)%DRL)
      ENDDO
#endif

      DEALLOCATE( ff, vv, ww )

      !-------------------------------------
      ! WARNING: continuation_hook expects the solution of
      ! J(dx) = - R   and not
      ! J(dx) = + R
      !-------------------------------------
      rhs = -rhs

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    |rhs|_L-infty = ', MAXVAL(ABS(rhs))
#endif

   ENDIF  

   IF (matflag == MATRIX_ONLY .OR. matflag == RHS_MATRIX .OR. matflag == RECOVER_MATRIX ) THEN

#if DEBUG > 0
      WRITE(*,*) '    filling new Jacobian matrix'
#endif

      CALL extract (xsol,  u0) 
      !------------------------------------------------------------------
      !-------------GENERATION OF THE JACOBIAN MATRIX--------------------
      !-------------OF THE COUPLED EQUATION SYSTEM-----------------------
      !             Jacobian  <--- [(u0.V)_ + (_.V)u0)]  +  K_  +  V_ (ibp)

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    Re = ', Re
#endif

      CALL ComputeJacobianMatrix (np, mm, jj, jj_L, js_Axis, js_D, DESINGULARIZE, Jacobian, Re, u0)
      CALL par_mumps_master (NUMER_FACTOR, 1, Jacobian, 0)

#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    |Jacobian%e|_L-infty = ', MAXVAL(ABS(Jacobian%e))
#endif

   ENDIF

END SUBROUTINE matrix_residual_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE mass_matrix_fill_conwrap(xsol, rhs)
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
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(Nx) :: xsol, rhs

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to mass_matrix_fill_conwrap'
#endif


   IF ( .NOT.Mass_init ) THEN
#if DEBUG > 0
      WRITE(*,*) '    allocating Mass matrix'
#endif

      ALLOCATE( Mass%i      (SIZE(Jacobian%i))       ); Mass%i       = Jacobian%i
      ALLOCATE( Mass%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass%i_mumps = Jacobian%i_mumps
      ALLOCATE( Mass%j      (SIZE(Jacobian%j))       ); Mass%j       = Jacobian%j
      ALLOCATE( Mass%e      (SIZE(Jacobian%e))       ); Mass%e       = 0d0
      Mass_init = .TRUE.

   ENDIF

   Mass%e = 0d0
   CALL qc_0y0_zero_sp_M (mm, jj, 1d0,  Mass)
   ! impose boundary conditions on the Mass Matrix
   CALL Dirichlet_c_M_MASS (np, js_Axis, js_D,  Mass)

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |Mass%e|_L-infty = ', MAXVAL(ABS(Mass%e))
#endif

END SUBROUTINE mass_matrix_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE Lns_matrix_fill_conwrap(xsol, matflag)
!
! Input:
!    xsol      Solution vector
!    matflag
! 
! Return Value:
! 
   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: xsol
   INTEGER, INTENT(IN)         :: matflag
   ! local variables
   INTEGER       :: i
   LOGICAL, SAVE :: base_alloc = .FALSE.
   COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: Lns_save

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to Lns_matrix_fill_conwrap'
   WRITE(*,*) '    matflag = ', matflag
#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |xsol|_L-infty = ', MAXVAL(ABS(xsol))
#endif

   IF (matflag == NEW_BASE_LNS .OR. matflag == PERT_LNS) THEN

      IF (.NOT. Lns_cmplx_init) THEN
         ALLOCATE( Lns_cmplx%i      (SIZE(Jacobian%i))       ); Lns_cmplx%i       = Jacobian%i
         ALLOCATE( Lns_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) ); Lns_cmplx%i_mumps = Jacobian%i_mumps
         ALLOCATE( Lns_cmplx%j      (SIZE(Jacobian%j))       ); Lns_cmplx%j       = Jacobian%j
         ALLOCATE( Lns_cmplx%e      (SIZE(Jacobian%e))       ); Lns_cmplx%e       = CMPLX(0d0, 0d0,KIND=8)
         Lns_cmplx_init = .TRUE.
      ENDIF

      ! fill the Lns matrix
      !
#if DEBUG > 1
      WRITE(*,*) '*check*'
      WRITE(*,*) '    Re   = ', Re
      WRITE(*,*) '    beta = ', beta
#endif

      CALL extract(xsol, u0)
      Lns_cmplx%e = CMPLX(0d0, 0d0,KIND=8)
      CALL qc_1y1_sp_gg_3d_M  (mm, jj,            1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
      CALL qc_oseen2y_sp_3d_M (mm, jj,             u0,    beta,  Lns_cmplx) ! + linearized terms
      CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,     -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
      CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,     -1d0,    beta,  Lns_cmplx) ! - velocity divergence
      CALL Dirichlet_c_3d_M   (np, js_Axis, js_D,                Lns_cmplx) ! Dirichlet BCs
   
      IF (DESINGULARIZE) THEN
         !! column
         !WHERE (Lns_cmplx%j == Nx)
         !   Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
         !ENDWHERE
         ! row
         DO i = Lns_cmplx%i(Nx), Lns_cmplx%i(Nx + 1) - 1
            Lns_cmplx%e(i) = CMPLX(0d0,0d0,KIND=8)
            IF (Lns_cmplx%j(i) == Nx) Lns_cmplx%e(i) = CMPLX(1d0,0d0,KIND=8)
         ENDDO
      ENDIF

      IF (matflag == NEW_BASE_LNS) THEN

         IF (.NOT.base_alloc) THEN
            ALLOCATE(Lns_save(SIZE(Lns_cmplx%e)))
            base_alloc = .TRUE.
         ENDIF

         Lns_save = Lns_cmplx%e

      ENDIF

   ELSEIF (matflag == BASE_LNS) THEN

      Lns_cmplx%e = Lns_save

   ENDIF

!#if DEBUG > 3
!   CALL save_vector( DBLE(Lns_cmplx%e),'Lns_r.vector')
!   CALL save_vector(AIMAG(Lns_cmplx%e),'Lns_i.vector')
!#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |Lns_cmplx%e|_L-infty = ', MAXVAL(ABS(DBLE(Lns_cmplx%e))), MAXVAL(ABS(AIMAG(Lns_cmplx%e)))
#endif

END SUBROUTINE Lns_matrix_fill_conwrap

!------------------------------------------------------------------------------

SUBROUTINE matvec_mult_conwrap(xxx, yyy)
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
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(Nx)         :: xxx

   REAL(KIND=8), DIMENSION(Nx)         :: yyy

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to matvec_mult_conwrap'
#endif

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |x|_L-infty = ', MAXVAL(ABS(xxx))
#endif

   CALL dAtimx(yyy, Jacobian%e, Jacobian%j, Jacobian%i, xxx)

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |y|_L-infty = ', MAXVAL(ABS(yyy))
#endif

END SUBROUTINE matvec_mult_conwrap

!------------------------------------------------------------------------------

SUBROUTINE mass_matvec_mult_conwrap(xxx, yyy)
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
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(Nx)         :: xxx

   REAL(KIND=8), DIMENSION(Nx)         :: yyy

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to mass_matvec_mult_conwrap'
#endif

   IF ( .NOT.Mass_init ) THEN
      WRITE(*,*) '***************************'
      WRITE(*,*) '********* STRANGE *********'
      WRITE(*,*) '    nonexistent Mass matrix'
      WRITE(*,*) '********* STRANGE *********'
      WRITE(*,*) '***************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |x|_L-infty = ', MAXVAL(ABS(xxx))
#endif

   CALL dAtimx(yyy, Mass%e, Mass%j, Mass%i, xxx)

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |y|_L-infty = ', MAXVAL(ABS(yyy))
#endif

END SUBROUTINE mass_matvec_mult_conwrap

!------------------------------------------------------------------------------

SUBROUTINE Lns_matvec_mult_conwrap(xxx_r, xxx_i, yyy_r, yyy_i)
!
! Put the call to your matrix-vector multiply here.
! Input:
!    xxx         Vector of length number of unknowns
!
! Output:
!    yyy         Lns matrix times xxx.
!
! Return Value:
!
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(Nx) :: xxx_r, xxx_i

   REAL(KIND=8), DIMENSION(Nx) :: yyy_r, yyy_i

   COMPLEX(KIND=8), DIMENSION(Nx) :: xxx, yyy

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to Lns_matvec_mult_conwrap'
#endif

   xxx = CMPLX(xxx_r, xxx_i, KIND=8)
   yyy = 0d0

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |x|_L-infty = ', MAXVAL(ABS(DBLE(xxx))), MAXVAL(ABS(AIMAG(xxx)))
#endif

   CALL zAtimx(yyy, Lns_cmplx%e, Lns_cmplx%j, Lns_cmplx%i, xxx)

#if DEBUG > 1
   WRITE(*,*) '*check*'
   WRITE(*,*) '    |y|_L-infty = ', MAXVAL(ABS(DBLE(yyy))), MAXVAL(ABS(AIMAG(yyy)))
#endif

   yyy_r = DBLE(yyy)
   yyy_i = AIMAG(yyy)

END SUBROUTINE Lns_matvec_mult_conwrap

!------------------------------------------------------------------------------

SUBROUTINE assign_parameter_conwrap(param)
! 
! Put the call to a routine to assign the continuation parameter here.
! Input:
!    param     New value of continuation parameter.
! 
! Output:
! 
! Return Value:
! 

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8) :: param

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to assign_parameter_conwrap'
   WRITE(*,*)
   WRITE(*,*) '    assigning parameter'
#endif


   SELECT CASE( pd%param )

      CASE( REYNOLDS )
#if DEBUG > 0
         WRITE(*,*) '    Reynolds = ', param
         WRITE(*,*) '    oldRe    = ', Re
#endif

         Re          = param
         pd%reynolds = param

      CASE( OSCAR )
#if DEBUG > 0
         WRITE(*,*) '    oscar    = ', param
         WRITE(*,*) '    oldOscar = ', flow_parameters(1)
#endif

         CALL case_loca_changeOscar(param)
         flow_parameters(1)  = param
         pd%oscar            = param

      CASE( ROMEO )
#if DEBUG > 0
         WRITE(*,*) '    romeo    = ', param
         WRITE(*,*) '    oldRomeo = ', flow_parameters(2)
#endif

         CALL case_loca_changeRomeo(param)
         flow_parameters(2) = param
         pd%romeo           = param

      CASE( WHISKY )
#if DEBUG > 0
         WRITE(*,*) '    whisky    = ', param
         WRITE(*,*) '    oldWhisky = ', flow_parameters(3)
#endif

         CALL case_loca_changeWhisky(param)
         flow_parameters(3) = param
         pd%whisky          = param

   END SELECT


END SUBROUTINE assign_parameter_conwrap

!------------------------------------------------------------------------------

SUBROUTINE assign_bif_parameter_conwrap(bif_param)
! 
! Put the call to a routine to assign the continuation parameter here.
! Input:
!    param     New value of continuation parameter.
! 
! Output:
! 
! Return Value:
! 

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8) :: bif_param

#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to assign_bif_parameter_conwrap'
   WRITE(*,*)
   WRITE(*,*) '    assigning bifurcation parameter'
#endif

   SELECT CASE( pd%bif_param )

      CASE( REYNOLDS )
#if DEBUG > 0
         WRITE(*,*) '    Reynolds = ', bif_param
         WRITE(*,*) '    oldRe    = ', Re
#endif

         Re          = bif_param
         pd%reynolds = bif_param

      CASE( OSCAR )
#if DEBUG > 0
         WRITE(*,*) '    oscar    = ', bif_param
         WRITE(*,*) '    oldOscar = ', flow_parameters(1)
#endif
       
         CALL case_loca_changeOscar(bif_param)
         flow_parameters(1)  = bif_param
         pd%oscar            = bif_param

      CASE( ROMEO )
#if DEBUG > 0
         WRITE(*,*) '    romeo    = ', bif_param
         WRITE(*,*) '    oldRomeo = ', flow_parameters(2)
#endif

         CALL case_loca_changeRomeo(bif_param)
         flow_parameters(2) = bif_param
         pd%romeo           = bif_param

      CASE( WHISKY )
#if DEBUG > 0
         WRITE(*,*) '    whisky    = ', bif_param
         WRITE(*,*) '    oldWhisky = ', flow_parameters(3)
#endif

         CALL case_loca_changeWhisky(bif_param)
         flow_parameters(3) = bif_param
         pd%whisky          = bif_param

   END SELECT


END SUBROUTINE assign_bif_parameter_conwrap

!******************************************************************************
!******************************************************************************
!**** what follows is a collection of conwraps that have not been modified ****
!******************************************************************************
!******************************************************************************

SUBROUTINE assign_multi_parameter_conwrap(param_vec)
! Put the call to a routine to assign the continuation parameters here.
! Input:
!    param_vec     New values of continuation parameters.
!

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:), INTENT(IN) ::  param_vec  
   
   WRITE(*,*) 'multi_param not implemented'
   STOP
  
END SUBROUTINE assign_multi_parameter_conwrap

!------------------------------------------------------------------------------

SUBROUTINE calc_scale_vec_conwrap(x, scale_vec, numUnks)
! Put the call to a routine to calculate a scaling vector here.
! Input:
!    x          New value of continuation parameter.
!    numUnks    Number of unknowns on this proc, the length of x
!               and scale_vec.
!
! Output:
!    scale_vec  Vector of length number of unknowns used to scale
!               variables so that one type of unknown (e.g. pressure)
!               doesn't dominate over others. Used to balance the
!               variables and the arc-length variable in arc-length
!               continuation, and for scaling the NULL vector in
!               turning point tracking. Using reciprocal of the average
!               value of that variable type is a good choice. Vector
!               of all ones should suffice for most problems.
!

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
   REAL(KIND=8), DIMENSION(:)             :: scale_vec
   INTEGER,                    INTENT(IN) :: numUnks


   scale_vec = 1d0 

END SUBROUTINE calc_scale_vec_conwrap

!------------------------------------------------------------------------------

FUNCTION gsum_double_conwrap(summ) RESULT(output)
! Put the call to a routine to calculate a global sum.
! Just return summ for single processor jobs.
! Input:
!    summ    Value of double on this processor to be summed on all procs.
!
! Return Value:
!    The global sum is returned on all processors.

   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: summ
   REAL(KIND=8)             :: output
 
   output = summ
  
END FUNCTION gsum_double_conwrap

!------------------------------------------------------------------------------

FUNCTION  gmax_int_conwrap(maxx) RESULT(output)
! Put the call to a routine to calculate a global max.
! Just return maxx for single processor jobs.
! Input:
!    maxx    Value of integer on this processor to be maxed on all procs.
!
! Return Value:
!    The global max is returned on all processors.
!
! Only used by Eigensolver

   IMPLICIT NONE
 
   INTEGER, INTENT(IN) :: maxx
   INTEGER :: output
 
   output = maxx
  
END FUNCTION gmax_int_conwrap

!------------------------------------------------------------------------------

SUBROUTINE random_vector_conwrap(x, numOwnedUnks)
! Put a routine to calculate a random vector.
! Input:
!    numOwnedUnks  Length of owned nodes part of x.
!
! Output:
!    x             Random vector.
!
! Used by eigensolver only

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:) :: x
   INTEGER :: numOwnedUnks

   CALL RANDOM_SEED()    ! initialize seed based on date and time (at least with Intel compilers)
   CALL RANDOM_NUMBER(x) ! generate random number
  
END SUBROUTINE random_vector_conwrap

!------------------------------------------------------------------------------

SUBROUTINE perturb_solution_conwrap(x, x_old, scale_vec, numOwnedUnks)
! Put a routine to perturb the solution vector a little bit here.
! This is to move a converged solution at a singularity off
! of the singularity before doing continuation. This ain't pretty
! but has helped convergence on some turning point tracking problems.
! Input:
!    x_old         Current solution vector.
!    scale_vec     Work space for a vector to scale x.
!    numOwnedUnks  Length of owned nodes part of x, x_old, scale_vec
!
! Output:
!    x             Solution vector perturbed a bit.

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: x_old
   REAL(KIND=8), DIMENSION(0:) :: scale_vec
   INTEGER :: numOwnedUnks   
                  
   CALL random_vector_conwrap(x, numOwnedUnks)

   x = x_old * (1d0 + 1d-4 * (x - 0.5d0))

END SUBROUTINE perturb_solution_conwrap

!------------------------------------------------------------------------------

FUNCTION free_energy_diff_conwrap(x, x2) RESULT(output)
! Call to return the free energy difference betwen two solutions
! Input:
!    x    One solution vector
!    x2   Second solution vector
!
! Output:
!
! Return Value:
!    The difference in the free energy beween the two solutions

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: x2   

   REAL(KIND=8) :: output
 
   output = 1.0d0
  
END FUNCTION free_energy_diff_conwrap

!------------------------------------------------------------------------------

SUBROUTINE eigenvector_output_conwrap(j, num_soln_flag, xr, evr, xi, evi, step_num)
! Call to write out eigenvectors
! Input:
!    j    Eigenvector index/mode
!    num_soln_flag  =1 for real vector, real eigenvalue
!                    =2 for complex (imaginary part has info)
!    xr   Real part of eigenvector
!    evr  Real part of eigenvalue
!    xi   Imaginary part of eigenvector (NULL if num_soln_flag==1)
!    evi  Imaginary part of eigenvalue
!    step_num  integer step number for use in output

   IMPLICIT NONE

   INTEGER :: j
   INTEGER :: num_soln_flag
   REAL(KIND=8), DIMENSION(0:) :: xr
   REAL(KIND=8) :: evr
   REAL(KIND=8), DIMENSION(0:) :: xi
   REAL(KIND=8) :: evi
   INTEGER :: step_num

   ! dummy because we use our own

END SUBROUTINE eigenvector_output_conwrap

!------------------------------------------------------------------------------

SUBROUTINE create_shifted_matrix_conwrap()
!
! Allocates a sets sparsity pointers for shifted matrix.
! Only used by eigensolver
! 
   IMPLICIT NONE
 
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

SUBROUTINE shifted_matrix_fill_conwrap (sigma_loca)
!
! Routine to create shifted matrix  J-sigma_loca*M
! Only called by eigensolver
! 
   IMPLICIT NONE
 
   REAL(KIND=8), INTENT(IN) :: sigma_loca

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

SUBROUTINE shifted_linear_solver_conwrap(xvec, yvec, jac_flag, tol)
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
   IMPLICIT NONE
 
   REAL(KIND=8), DIMENSION(Nx) :: xvec
   REAL(KIND=8), DIMENSION(Nx) :: yvec
   INTEGER                     :: jac_flag
   REAL(KIND=8)                :: tol

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

SUBROUTINE destroy_shifted_matrix_conwrap()
!
! Deletes memory for shifted matrix.
! Only used by eigensolver
!
   IMPLICIT NONE

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

! END OF LOCA'S WRAPPERS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





!------------------------------------------------------------------------------
! other LOCA's subroutines which are not wrappers

SUBROUTINE param_output(param, omega)
! 
! Print parameter(s) to file
! 
   IMPLICIT NONE
   ! input variables
   REAL(KIND=8)           :: param
   REAL(KIND=8), OPTIONAL :: omega ! for Hopf tracking
   ! local variables
   INTEGER           :: fid = 22
   CHARACTER(LEN=50) :: filenm


#if DEBUG > 0
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to param_output'
   WRITE(*,*)
#endif

   SELECT CASE( pd%param )
      CASE( REYNOLDS )
         filenm = './locaOut/reynolds.dat'
      CASE( OSCAR )
         filenm = './locaOut/oscar.dat'
      CASE( ROMEO )
         filenm = './locaOut/romeo.dat'
      CASE( WHISKY )
         filenm = './locaOut/whisky.dat'
   END SELECT

#if DEBUG > 0
   WRITE(*,*) 'writing file: ', trim(filenm)
#endif

   OPEN(UNIT= fid, FILE= trim(filenm), POSITION= 'APPEND')
   WRITE(fid,*) param
   CLOSE(fid)

   ! print all parameters to file
   filenm = './locaOut/all.dat'
   OPEN(UNIT= fid, FILE= trim(filenm), POSITION= 'APPEND')
   WRITE(*,*) 'writing file: ', trim(filenm)

   IF (PRESENT(omega)) THEN
      WRITE(fid,*) pd%reynolds, pd%oscar, pd%romeo, pd%whisky, omega
   ELSE
      WRITE(fid,*) pd%reynolds, pd%oscar, pd%romeo, pd%whisky
   ENDIF

   CLOSE(fid)

   CALL case_loca_paramout()

END SUBROUTINE param_output
!
! C compatible version follows
!
!! SUBROUTINE param_output(param) &
!!    BIND(C, NAME='param_output')
!! ! 
!! ! Print parameter(s) to file
!! ! 
!!    USE ISO_C_BINDING
!! 
!!    IMPLICIT NONE
!!    ! input variables
!!    REAL(KIND=C_DOUBLE), VALUE :: param
!!    ! local variables
!!    INTEGER           :: fid = 22
!!    CHARACTER(LEN=50) :: filenm
!! 
!! 
!!    WRITE(*,*)
!!    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!    WRITE(*,*) '--> CALL to param_output'
!!    WRITE(*,*)
!! 
!!    SELECT CASE( pd%param )
!!       CASE( REYNOLDS )
!!          filenm = './locaOut/reynolds.dat'
!!          WRITE(*,*) 'writing file: ', trim(filenm)
!!       CASE( OSCAR )
!!          filenm = './locaOut/oscar.dat'
!!          WRITE(*,*) 'writing file: ', trim(filenm)
!!       CASE( ROMEO )
!!          filenm = './locaOut/romeo.dat'
!!          WRITE(*,*) 'writing file: ', trim(filenm)
!!       CASE( WHISKY )
!!          filenm = './locaOut/whisky.dat'
!!          WRITE(*,*) 'writing file: ', trim(filenm)
!!    END SELECT
!! 
!!    OPEN(UNIT= fid, FILE= trim(filenm), POSITION= 'APPEND')
!!    WRITE(fid,*) param
!!    CLOSE(fid)
!! 
!!    ! print all parameters to file
!!    filenm = './locaOut/all.dat'
!!    OPEN(UNIT= fid, FILE= trim(filenm), POSITION= 'APPEND')
!!    WRITE(*,*) 'writing file: ', trim(filenm)
!!    WRITE(fid,*) pd%reynolds, pd%oscar, pd%romeo, pd%whisky
!! 
!!    CLOSE(fid)
!! 
!!    CALL case_loca_paramout()
!! 
!! END SUBROUTINE param_output

!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_loca(x_vec, filenm)

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: x_vec
   CHARACTER(*)                :: filenm
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_plot
   REAL(KIND=8), DIMENSION(np_L)          :: p_plot


   IF ( p_in%write_plots_flag ) THEN

      CALL extract(x_vec, u_plot, p_plot)

      CALL vtk_plot_P2 (rr, jj, jj_L, u_plot, p_plot, &
        trim(p_in%plot_directory)//'locaContSolution'//trim(filenm)//'.vtk' )

   ENDIF


END SUBROUTINE vtk_plot_loca
!
! C compatible version follows
!
!! SUBROUTINE vtk_plot_loca(x_vec, filenm, filenmLen) &
!!    BIND(C, NAME='vtk_plot_loca')
!! 
!!    USE ISO_C_BINDING
!! 
!!    IMPLICIT NONE
!!    ! input variables
!!    REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
!!    CHARACTER(KIND=C_CHAR)             :: filenm
!!    INTEGER(KIND=C_INT), VALUE         :: filenmLen
!!    ! local variables
!!    REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_plot
!!    REAL(KIND=8), DIMENSION(np_L)          :: p_plot
!! 
!! 
!!    IF ( p_in%write_plots_flag ) THEN
!! 
!!       CALL extract(x_vec, u_plot, p_plot)
!! 
!!       CALL vtk_plot_P2 (rr, jj, jj_L, u_plot, p_plot, &
!!         trim(p_in%plot_directory)//'locaContSolution'//filenm(1:filenmLen)//'.vtk' )
!! 
!!    ENDIF
!! 
!! 
!! END SUBROUTINE vtk_plot_loca

!------------------------------------------------------------------------------

SUBROUTINE compute_eigen(x_vec, filenm, shiftIm)
!
! Compute eigenvalues and eigenvectors and save them to file
!

   IMPLICIT NONE

   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: x_vec
   CHARACTER(*)                :: filenm
   REAL(KIND=8)                :: shiftIm

   ! local variables
   !TYPE(CSR_MUMPS_Complex_Matrix) :: Lns_cmplx, Mass_cmplx ! moved to global_variables

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
   !
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
#if DEBUG > 1
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

   IF (.NOT. Lns_cmplx_init) THEN
      ALLOCATE( Lns_cmplx%i      (SIZE(Jacobian%i))       ); Lns_cmplx%i       = Jacobian%i
      ALLOCATE( Lns_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) ); Lns_cmplx%i_mumps = Jacobian%i_mumps
      ALLOCATE( Lns_cmplx%j      (SIZE(Jacobian%j))       ); Lns_cmplx%j       = Jacobian%j
      ALLOCATE( Lns_cmplx%e      (SIZE(Jacobian%e))       ); Lns_cmplx%e       = CMPLX(0d0, 0d0,KIND=8)
      Lns_cmplx_init = .TRUE.
   ENDIF

   IF (.NOT. Mass_cmplx_init) THEN
      ALLOCATE( Mass_cmplx%i      (SIZE(Jacobian%i))       ); Mass_cmplx%i       = Jacobian%i
      ALLOCATE( Mass_cmplx%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass_cmplx%i_mumps = Jacobian%i_mumps
      ALLOCATE( Mass_cmplx%j      (SIZE(Jacobian%j))       ); Mass_cmplx%j       = Jacobian%j
      ALLOCATE( Mass_cmplx%e      (SIZE(Jacobian%e))       ); Mass_cmplx%e       = CMPLX(0d0, 0d0, KIND=8)
      Mass_cmplx_init = .TRUE.
   ENDIF

   ! (3)
   ! fill the Lns matrix
   ! NOTE: the sign of the Lns matrix has to be changed
   ! since the problem we are solving is:
   ! lambda*Mass*x = -Lns*x
   !
#if DEBUG > 1
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
                                'directEigenvalues'//trim(filenm)//trim(shiftName)//'.dat')

#ifdef SAVEEIGENVECTOR
            CALL Save_eigenvectors (directEigenvectors, &
                               trim(p_in%eigen_output_directory)// &
                               'directEigenvectors'//trim(filenm)//trim(shiftName)//'.dat')
#endif

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
                                           'directEigenvectorsRe'//trim(filenm)//trim(shiftName)//'.vtk')
               CALL vtk_plot_eigenvectors (rr, jj, AIMAG(directEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'directEigenvectorsIm'//trim(filenm)//trim(shiftName)//'.vtk')
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
                                 'adjointEigenvalues'//trim(filenm)//trim(shiftName)//'.dat')

#ifdef SAVEEIGENVECTOR
            CALL Save_eigenvectors (adjointEigenvectors, &
                                trim(p_in%eigen_output_directory)// &
                                'adjointEigenvectors'//trim(filenm)//trim(shiftName)//'.dat')
#endif

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
                                           'adjointEigenvectorsRe'//trim(filenm)//trim(shiftName)//'.vtk')
               CALL vtk_plot_eigenvectors (rr, jj, AIMAG(adjointEigenvectors(:,1:eigen_plotNum)), &
                                           trim(p_in%plot_directory)// &
                                           'adjointEigenvectorsIm'//trim(filenm)//trim(shiftName)//'.vtk')
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
                           './structSensOut/'//'structuralSensitivity'//trim(filenm)//'.vtk')

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
   Lns_cmplx_init  = .FALSE.
   Mass_cmplx_init = .FALSE.


END SUBROUTINE compute_eigen
!
! part of C version follows
!
!!SUBROUTINE compute_eigen(x_vec, filenm, filenmLen, shiftIm) &
!!   BIND(C, NAME='compute_eigen')
!!!
!!! Compute eigenvalues and eigenvectors and save them to file
!!!
!!   USE ISO_C_BINDING
!!
!!   IMPLICIT NONE
!!
!!   ! input variables
!!   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x_vec
!!   CHARACTER(KIND=C_CHAR)             :: filenm
!!   INTEGER(KIND=C_INT), VALUE         :: filenmLen
!!   REAL(KIND=C_DOUBLE), VALUE         :: shiftIm
!!
!!   ! local variables
!!   TYPE(CSR_MUMPS_Complex_Matrix)     :: Lns_cmplx, Mass_cmplx
!!
!!   LOGICAL, DIMENSION(velCmpnnts, number_of_sides) :: Dir_eigen
!!   TYPE(dyn_int_line), DIMENSION(velCmpnnts)       :: js_D_eigen
!!   LOGICAL                                         :: DESINGULARIZE_eigen
!!
!!   REAL(KIND=8), DIMENSION(SIZE(Jacobian%e)) :: Jacobian_save
!!
!!   INTEGER :: i, k, eigen_plotNum
!!
!!   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: directEigenvalues,  adjointEigenvalues
!!   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: directEigenvectors, adjointEigenvectors
!!   REAL(KIND=8),    DIMENSION(:),   ALLOCATABLE :: structuralsens
!!
!!   INTEGER :: statusMsg
!!
!!   CHARACTER(LEN=128) :: shiftName, & ! used to insert the shift in file names
!!                         shiftNameRe, shiftNameIm
!!   INTEGER            :: shiftNameTemp
!!
!!   REAL(KIND=8) :: dataIdentifier = 313d0 ! Donald Duck's plate number!
!!   INTEGER      :: shiftsNumber = 1
!!   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: shifts
!!
!!   COMPLEX(KIND=8), DIMENSION(Nx) :: tmpEigen1, tmpEigen2 ! used to check residuals
!!END SUBROUTINE compute_eigen

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
                          tempEigenvectorRe, tempEigenvectorIm)

   directEigenvector  = CMPLX(tempEigenvectorRe, tempEigenvectorIm, KIND=8)



   CALL read_eigenvector (Nx, p_in%structSens_eigenNumber,              &
                          trim(p_in%structSens_adjointEigen_name),      &
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
