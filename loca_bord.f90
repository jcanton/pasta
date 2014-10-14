! -----------------------------------------------------------------------------
!   LOCA 1.0: Library of Continuation Algorithms
!   Copyright (C) 2001, Sandia National Laboratories
! -----------------------------------------------------------------------------

MODULE Loca_bord

   USE Loca_types
   USE Loca_util     ! scaled_dot_prod, dp
   USE loca_wrappers ! linear_solver_conwrap, assign_parameter_conwrap, assign_bif_parameter_conwrap, calc_scale_vec_conwrap, ...

   IMPLICIT NONE

   INTEGER :: AGS_option = 0


CONTAINS
!=======

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION continuation_hook(x, delta_x, con, reltol, abstol) RESULT (output)
! return TRUE for non-continuation problems, as well as zero and
! first order continuation runs

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct)            :: con
   REAL(KIND=8)                :: reltol
   REAL(KIND=8)                :: abstol

   INTEGER :: output

   INTEGER :: converged = TRUE

   ! If arc_length continuation is being done, then branch here
   ! to do the rest of the bordering algorithm. This involves a
   ! second matrix fill and linear solve. (See JNS thesis
   ! eqns 7.18-7.19 for equations, 7.21-7.23 for method.
   IF (con%general_info%method == ARC_LENGTH_CONTINUATION .AND. con%private_info%step_num > 0) THEN

      converged = arc_length_bordering_alg(x, delta_x, con, reltol, abstol)

   ! If turning point calculations are being done, then branch here
   ! to do the rest of the algorithm. This involves 6 more matrix fills
   ! and 3 more linear solves in the current implementation.
   ELSE IF (con%general_info%method == TURNING_POINT_CONTINUATION) THEN

#if DEBUG > 2
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1  ~~ max delta_x =', MAXVAL(ABS(delta_x))
#endif

      converged = turning_point_alg(x, delta_x, con, reltol, abstol)

#if DEBUG > 2
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.2  ~~ converged =', converged
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.3  ~~ max delta_x =', MAXVAL(ABS(delta_x))
#endif

   ! If pitchfork bifurcation calculations are being done, then branch
   ! here to do the rest of the algorithm.
   ELSE IF (con%general_info%method == PITCHFORK_CONTINUATION) THEN

      converged = pitchfork_alg(x, delta_x, con, reltol, abstol)

   ! If Hopf bifurcation calculations are being done, then branch
   ! here to do the rest of the algorithm.
   ELSE IF (con%general_info%method == HOPF_CONTINUATION) THEN

#if DEBUG > 2
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.4  ~~ max delta_x =', MAXVAL(ABS(delta_x))
#endif

      converged = hopf_alg(x, delta_x, con, reltol, abstol)

#if DEBUG > 2
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.5  ~~ converged =', converged
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.6  ~~ max delta_x =', MAXVAL(ABS(delta_x))
#endif

   ! If phase_transition tracking is being done, then branch here
   ! to do the rest of the algorithm. This involves
   ! 3 more linear solves in the current implementation.
   ELSE IF (con%general_info%method == PHASE_TRANSITION_CONTINUATION) THEN

      converged = phase_transition_alg(x, delta_x, con, reltol, abstol)

   ! No bordering alg needed for zero and first order continuations
   ELSE IF (     con%general_info%method == ZERO_ORDER_CONTINUATION &
            .OR. con%general_info%method == FIRST_ORDER_CONTINUATION &
            .OR. con%general_info%method == LOCA_LSA_ONLY &
	         .OR. (con%general_info%method == ARC_LENGTH_CONTINUATION .AND. con%private_info%step_num == 0) &
            .OR. (con%general_info%method == MANIFOLD_CONTINUATION .AND. con%private_info%step_num == 0)) THEN

      converged = TRUE

   ELSE IF (con%general_info%method == MANIFOLD_CONTINUATION) THEN

      converged = manifold_alg(x, delta_x, con, reltol, abstol)

   ! perform error check
   ELSE
      IF (con%general_info%printproc > 1) &
         WRITE(*,*) '\nERROR continuation_hook: Unknown method ', con%general_info%method
      output = -1
      RETURN
   END IF

   ! Turn off flag that is true for first Newton iteration of each solve */

   con%private_info%first_iter = FALSE

   ! Return flag indicating convergence of the Newton iter */

   output = converged

END FUNCTION continuation_hook

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION arc_length_bordering_alg(x, delta_x, con, reltol, abstol) RESULT (output)
! Bordering algorithm for arc-length continuation is performed.
! Notation is from page 181 of JNS thesis. The vector "w" is
! already calculated (eq 7-21b), and in array delta_x.
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct) :: con
   REAL(KIND=8) :: reltol
   REAL(KIND=8) :: abstol
   INTEGER :: output

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dx2
   REAL(KIND=8) :: y
   REAL(KIND=8) ::arc_length_equation = 0, param_update, scaled_resid
   REAL(KIND=8) ::al_eq_soln_contrib, al_eq_param_contrib
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL
   INTEGER :: i

   !************************** BEGIN EXECUTION *******************************


   ! Next, "u" is calculated and stored as dx2. (eq. 7-21a) This is exactly
   ! the same linear solve as a continuation predictor step.

   ALLOCATE(dx2(0:con%general_info%numUnks-1))
   ALLOCATE(NULL(0:con%general_info%numUnks-1))

   CALL calc_rhs_continuation(ARC_CONT_SOL2, x, dx2, NULL, NULL, NULL, con%general_info%param, &
                              con%general_info%perturb, NULL, con%general_info%numUnks, con%general_info%numOwnedUnks)
   i = linear_solver_conwrap(dx2, CHECK_JACOBIAN, NULL)

   ! Calculate the arc-length equation ("g" in bordering algorithm notation)

   al_eq_soln_contrib = &
        scaled_dot_prod(x, con%private_info%x_tang, con%private_info%scale_vec, con%general_info%numOwnedUnks) &
      - scaled_dot_prod(con%private_info%x_old,con%private_info%x_tang,con%private_info%scale_vec, con%general_info%numOwnedUnks)

   al_eq_param_contrib = con%private_info%dp_ds * (con%general_info%param - con%private_info%param_old)

   arc_length_equation = al_eq_soln_contrib + al_eq_param_contrib - con%private_info%arc_step

   ! Calculate "y", the -update to the continuation parameter (JNS eq. 7-22) */
   y = (arc_length_equation -  &
        scaled_dot_prod(con%private_info%x_tang, delta_x, con%private_info%scale_vec, con%general_info%numOwnedUnks)) &
      /(con%private_info%dp_ds -  &
        scaled_dot_prod(con%private_info%x_tang, dx2, con%private_info%scale_vec, con%general_info%numOwnedUnks))

   ! y is subtracted because the minus sign in Newton's method has been
   ! switched to the update term -- the same reason delta_x is subtracted
   ! in update_soln
   con%general_info%param = con%general_info%param - y

   CALL assign_parameter_conwrap(con%general_info%param)

   ! Modify delta_x, which will be returned to the nonlinear solver (eq.7-23)
   delta_x = delta_x - y * dx2

   DEALLOCATE(dx2)
   DEALLOCATE(NULL)

   ! Check whether or not the arc-length equation and continuation param
   ! are converged

   param_update = ABS(y) / (reltol * ABS(con%general_info%param) + abstol)
   scaled_resid = ABS(arc_length_equation) &
                / (reltol * ABS(con%private_info%arc_step) + abstol)

   IF (con%general_info%printproc > 7) THEN
      WRITE(*,*) 'Fraction of arc-step due to soln, param = ', &
                  al_eq_soln_contrib/con%private_info%arc_step, &
                  al_eq_param_contrib/con%private_info%arc_step
   END IF

   IF (con%general_info%printproc > 4) THEN
      WRITE(*,*) 'Arc Length Continuation: Convergence Criteria'
      WRITE(*,*) 'Variable     Scaled Update (<1)  Unscaled Update  New Value'
      WRITE(*,*) '***********************************************************'
      WRITE(*,*) 'parameter    ', param_update, -y, con%general_info%param
      WRITE(*,*) 'arc length   ', scaled_resid, arc_length_equation, &
                                  con%private_info%arc_step + arc_length_equation
      WRITE(*,*) '***********************************************************'
   END IF

   IF ((param_update < 1.0) .AND. (scaled_resid < 1.0)) THEN
      output = TRUE
   ELSE
      output = FALSE
   END IF

END FUNCTION arc_length_bordering_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION turning_point_alg(x, delta_x, con, reltol, abstol) RESULT (output)
! Algorithm for locking on to a turning point.
! Theory currently comes from a TEX document of Louis Romero. (AGS 1/98)
! Lines labeled  SCALED  are additions for new scaling (section 4 of same
! write-up). The SCALED lines can be commented out to recover unscaled version.

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct) :: con
   REAL(KIND=8) :: reltol
   REAL(KIND=8) :: abstol
   INTEGER :: output

!#define SCALE_TP

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a, b, c, d, x_tmp
   REAL(KIND=8) :: dt_p
!#ifdef SCALE_TP
   REAL(KIND=8) ::a_big, b_big, c_big, d_big
!#endif
   REAL(KIND=8) ::param_update, r_update, vecnorm, gnum_unks, RayQ
   INTEGER :: first
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL
   INTEGER :: i

   first = TRUE

! This flag can be 0 or 1 and effects scaling
! Seems to work better as 0 for TP calcs
AGS_option=0

   !************************** BEGIN EXECUTION *******************************/

   ! Allocate arrays for turning point bordering algorithm

   ALLOCATE(    a(0:con%general_info%numUnks-1))
   ALLOCATE(    b(0:con%general_info%numUnks-1))
   ALLOCATE(    c(0:con%general_info%numUnks-1))
   ALLOCATE(    d(0:con%general_info%numUnks-1))
   ALLOCATE(x_tmp(0:con%general_info%numUnks-1))
   ALLOCATE( NULL(0:con%general_info%numUnks-1))

   ! construct "a" vector from delta_x
   a = - delta_x

   ! Next, "b" is calculated. This is exactly
   ! the same linear solve as a continuation predictor step except
   ! the bif_param is perturbed, not the continuation param.
   CALL calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, &
                              con%turning_point_info%bif_param, con%general_info%perturb, &
                              NULL, con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL)

   ! First time through this routine, initialize null vector y and
   ! norm vector phi with scaled b.  First Newton iter of subsequent
   ! continuation steps, set guess for y and phi equal to converged
   ! null vector of previous converged step.
   ! If using a previous null vector, it has already been read in.
   IF (first == TRUE) THEN

      IF (con%general_info%nv_restart == FALSE) THEN

         CALL calc_scale_vec_conwrap(x, con%private_info%scale_vec, con%general_info%numUnks)

         con%private_info%x_tang = con%private_info%scale_vec

         vecnorm = ltransnorm(b, con%private_info%scale_vec)

         con%private_info%x_tang = b/vecnorm

      ! Load in old null vector if specified */
      ELSE

         IF (con%general_info%printproc > 7) WRITE(*,*) 'Loading in previous null vector ...'

         con%private_info%x_tang = con%turning_point_info%nv

         CALL calc_scale_vec_conwrap(x, con%private_info%scale_vec, con%general_info%numUnks)

         vecnorm = ltransnorm(con%private_info%x_tang, con%private_info%scale_vec)

         con%private_info%x_tang = con%private_info%x_tang / vecnorm

      END IF

      first = FALSE

   ELSE

      ! This section is optional, changing scale vector between Newton iters
      CALL calc_scale_vec_conwrap(x, con%private_info%scale_vec, con%general_info%numUnks)

      vecnorm = ltransnorm(con%private_info%x_tang, con%private_info%scale_vec)

      IF (con%general_info%printproc > 7) WRITE(*,*) 'Rescaling r_vec by ',vecnorm,' to make its length = 1'

      con%private_info%x_tang = con%private_info%x_tang /vecnorm

   END IF


   ! Rescale a and b vectors as in Louie's write-up, section 4

!#ifdef SCALE_TP
  b_big = 1.d0 / ltransnorm(b, con%private_info%scale_vec)   ! SCALED
  a_big = -ltransnorm(a, con%private_info%scale_vec) * b_big ! SCALED
  a = a + a_big * b                   ! SCALED
  b = b * b_big                       ! SCALED
!#endif

   ! Next, "c" is calculated as a function of a and y.
   CALL calc_rhs_continuation(TP_CONT_SOL3, x, c, a, con%private_info%scale_vec, x_tmp, &
                              con%turning_point_info%bif_param, con%general_info%perturb, &
                              con%private_info%x_tang, con%general_info%numUnks, con%general_info%numOwnedUnks)

   ! Get null vector residual now.
   RayQ = null_vector_resid(0.d0, 0.d0, con%private_info%x_tang, NULL, FALSE)

   i = linear_solver_conwrap(c, SAME_BUT_UNSCALED_JACOBIAN, x_tmp)

   ! Next, "d" is calculated as a function of b and y.
   CALL calc_rhs_continuation(TP_CONT_SOL4, x, d, b, con%private_info%scale_vec, x_tmp,  &
                              con%turning_point_info%bif_param, con%general_info%perturb,  &
                              con%private_info%x_tang, con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp)

   ! Calculate the updates to bif_param (stored in dt_p),
   ! y (stored in c), and x (stored in -delta_x).

   ! Rescale c and d vectors as in Louie's write-up, section 4

!#ifdef SCALE_TP
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1.1  ~~ d =', MAXVAL(ABS(d))
WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1.2  ~~ scale_vec =', MAXVAL(ABS(con%private_info%scale_vec))
WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1.3  ~~ ltransnorm(d, con%private_info%scale_vec) =', &
           ltransnorm(d, con%private_info%scale_vec)
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
   d_big = 1.d0 / ltransnorm(d, con%private_info%scale_vec)   ! SCALED
   c_big = -ltransnorm(c, con%private_info%scale_vec) * d_big ! SCALED
   c = c + c_big * d                   ! SCALED
   d = d * d_big                       ! SCALED
!#endif
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1.4  ~~ d_big =', d_big
!WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~  Passato 2.2.1.2  ~~ max b =', MAXVAL(ABS(b))
! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG

   dt_p = ( (1.d0 - AGS_option * ltransnorm(con%private_info%x_tang, con%private_info%scale_vec)) &
            - ltransnorm(c, con%private_info%scale_vec) ) &
        / ltransnorm(d, con%private_info%scale_vec)

   c = c + dt_p * d + (AGS_option-1) * con%private_info%x_tang

   ! dt_p change meaning from Beta to alpha here */

!#ifdef SCALE_TP
   dt_p = dt_p * d_big + c_big ! SCALED
!#endif

   delta_x = -a

   delta_x = delta_x - dt_p * b

   ! dt_p change meaning from alpha to dt_p here

!#ifdef SCALE_TP
   dt_p = a_big + dt_p * b_big ! SCALED
!#endif

   ! Update the bif_param and the null vector y

   con%turning_point_info%bif_param = con%turning_point_info%bif_param + dt_p
   CALL assign_bif_parameter_conwrap(con%turning_point_info%bif_param)

   con%private_info%x_tang = con%private_info%x_tang + c

   ! Check whether or not the arc-length equation and continuation param
   ! are converged

   param_update = ABS(dt_p) / (reltol * ABS(con%turning_point_info%bif_param) + abstol)

   ! get update norm of y, to check for convergence */

   r_update = SUM((c / (ABS(con%private_info%x_tang)*reltol + abstol))**2)

   gnum_unks = gsum_double_conwrap(1.d0*con%general_info%numOwnedUnks)

   r_update = sqrt( gsum_double_conwrap(r_update)/gnum_unks )

   IF (RayQ == -1.0) THEN
      IF (con%general_info%printproc > 1) WRITE(*,*) 'Rayleigh Quotient Error: zero denom'
   ELSE
      IF (con%general_info%printproc > 4) WRITE(*,*) 'Rayleigh Quotient before updates: ', RayQ
   END IF

   IF (con%general_info%printproc > 4) THEN
      WRITE(*,*) 'Turning Point Continuation: Convergence Criteria'
      WRITE(*,*) 'Variable     Scaled Update (<1)  Unscaled Update  New Value'
      WRITE(*,*) '***********************************************************'
      WRITE(*,*) 'parameter    ', param_update, dt_p, con%turning_point_info%bif_param
      WRITE(*,*) 'Null vector  ', r_update
      WRITE(*,*) '***********************************************************'
   END IF

   DEALLOCATE(a)
   DEALLOCATE(b)
   DEALLOCATE(c)
   DEALLOCATE(d)
   DEALLOCATE(x_tmp)
   DEALLOCATE(NULL)

   ! return convergence status of the parameter and null vector */

   IF ((param_update < 1.0) .AND. (r_update < 10.0)) THEN
      output = TRUE
   ELSE
      output = FALSE
   END IF

END FUNCTION turning_point_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION pitchfork_alg(x, delta_x, con, reltol, abstol) RESULT (output)

! Algorithm for locking on to a Pitchfork bifurcation
! This assumes that con%private_info%x_tang contains an asymmetric vector
! the first time through, and an initial guess for the
! null vector every time (with both meanings the first time).
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct) :: con
   REAL(KIND=8) :: reltol
   REAL(KIND=8) :: abstol
   INTEGER :: output

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a,  b, c, d, e, f, x_tmp
   REAL(KIND=8) :: dt_p
   REAL(KIND=8) :: param_update, r_update, vecnorm, gnum_unks, tmp, RayQ
   REAL(KIND=8) :: eps_update, d_eps, ipa, ipb, ipc, ipx, ltd, lte, ltf
   REAL(KIND=8) :: t1, t2, t3, t4
   INTEGER :: first
   REAL(KIND=8) :: eps
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL
   INTEGER :: i

   first = TRUE

! This flag can be 0 or 1 and effects scaling
! Seems to work better as 1 for PF calcs
AGS_option=1

   !*************************** BEGIN EXECUTION ********************************

   ! first time through this routine, set psi asymmetric vector to con-tang
   !   and set epsilon "asymmetry variable" to 0.0 */

   IF (first == TRUE) THEN

      IF (con%general_info%printproc > 7) WRITE(*,*) 'In pitchfork alg, AGS_option= ', AGS_option

      eps = 0.d0

      tmp = SQRT(ip(con%pitchfork_info%psi, con%pitchfork_info%psi))

      IF (tmp == 0.d0) THEN
        IF (con%general_info%printproc > 1) THEN
          WRITE(*,*) 'ERROR in pitchfork alg: Asymmetric vector must be supplied'
        END IF
        output = -1
        RETURN
      END IF

      con%pitchfork_info%psi = con%pitchfork_info%psi / tmp

      ! con%pitchfork_info%psi is used as phi vector as well
      con%private_info%x_tang = con%pitchfork_info%psi

      first = FALSE

      t1 = ip(x, con%pitchfork_info%psi)
      t2 = ip(x, x)
      t3 = ip(con%pitchfork_info%psi, con%pitchfork_info%psi)
      t4 = t1 / SQRT(t2*t3)
      IF (con%general_info%printproc > 7) THEN
        WRITE(*,*) 'Pitchfork Alg: On input <x,psi> = ', t4,' Should be 0.0.'
      END IF
   END IF

   ! calculate variable averages, to be used in ltransnorm calls below
   CALL calc_scale_vec_conwrap(x, con%private_info%scale_vec, con%general_info%numOwnedUnks)

   ! Make sure the null vector phi is length 1 to begin with
   vecnorm = ltransnorm(con%private_info%x_tang, con%private_info%scale_vec)

   IF (con%general_info%printproc > 7) WRITE(*,*) 'Rescaling phi by ',vecnorm,' to make it length 1'

   con%private_info%x_tang = con%private_info%x_tang / vecnorm

   ! Allocate arrays for pitchfork bordering algorithm */
   ALLOCATE(    a(0:con%general_info%numUnks-1))
   ALLOCATE(    b(0:con%general_info%numUnks-1))
   ALLOCATE(    c(0:con%general_info%numUnks-1))
   ALLOCATE(    d(0:con%general_info%numUnks-1))
   ALLOCATE(    e(0:con%general_info%numUnks-1))
   ALLOCATE(    f(0:con%general_info%numUnks-1))
   ALLOCATE(x_tmp(0:con%general_info%numUnks-1))

   ALLOCATE(NULL(0:con%general_info%numUnks-1))

   ! Begin real stuff
   ! construct "a" vector from delta_x
   a = - delta_x

   ! Next, "b" is calculated. This is exactly
   ! the same linear solve as a continuation predictor step except
   ! the bif_param is perturbed, not the continuation param.
   CALL calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, &
                              con%pitchfork_info%bif_param, con%general_info%perturb, &
                              NULL, con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL)

   ! Next, "c" is calculated using just con%pitchfork_info%psi as rhs

   c = - con%pitchfork_info%psi

   i = linear_solver_conwrap(c, OLD_JACOBIAN, NULL)

   ! Next, "d" is calculated as a function of a and phi.

   CALL calc_rhs_continuation(TP_CONT_SOL3, x, d, a, con%private_info%scale_vec, x_tmp, &
                              con%pitchfork_info%bif_param, con%general_info%perturb,   &
                              con%private_info%x_tang, con%general_info%numUnks,        &
                              con%general_info%numOwnedUnks)

   ! Get null vector residual now.
   RayQ = null_vector_resid(0.d0, 0.d0, con%private_info%x_tang, NULL, FALSE)

   i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp)

   IF (AGS_option==1) d = d + con%private_info%x_tang

   ! Next, "e" is calculated as a function of b and phi.
   CALL calc_rhs_continuation(TP_CONT_SOL4, x, e, b, con%private_info%scale_vec, x_tmp, &
                              con%pitchfork_info%bif_param, con%general_info%perturb,   &
                              con%private_info%x_tang, con%general_info%numUnks,        &
                              con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(e, SAME_BUT_UNSCALED_JACOBIAN, x_tmp)

   ! Next, "f" is calculated as a function of c and phi.
   CALL calc_rhs_continuation(TP_CONT_SOL3, x, f, c, con%private_info%scale_vec, x_tmp, &
                              con%pitchfork_info%bif_param, con%general_info%perturb,   &
                              con%private_info%x_tang, con%general_info%numUnks,        &
                              con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(f, SAME_BUT_UNSCALED_JACOBIAN, x_tmp)

   IF (AGS_option==1)  f = f + con%private_info%x_tang

   ! Calculate the updates to bif_param (stored in dt_p),
   ! phi (stored in d), and x (stored in -delta_x).
   ipx = ip(x, con%pitchfork_info%psi)
   ipa = ip(a, con%pitchfork_info%psi)
   ipb = ip(b, con%pitchfork_info%psi)
   ipc = ip(c, con%pitchfork_info%psi)
   ltd = ltransnorm(d, con%private_info%scale_vec)
   lte = ltransnorm(e, con%private_info%scale_vec)
   ltf = ltransnorm(f, con%private_info%scale_vec)

   d_eps = - eps + ( (ipx+ipa)*lte + ipb*(1.D0 - ltd) ) / (ipb*ltf - ipc*lte)

   dt_p = ( 1.d0 - ltd - ltf*(d_eps + eps)) / lte

   ! Negative of update vector here
   delta_x = - a - c*(eps + d_eps) - b*dt_p

   ! use c space for delta phi
   c = - con%private_info%x_tang + d + e*dt_p + f*(eps + d_eps)

   ! Update the bif_param, eps,  and the null vector phi
   eps = eps + d_eps

   con%pitchfork_info%bif_param = con%pitchfork_info%bif_param + dt_p

   CALL assign_bif_parameter_conwrap(con%pitchfork_info%bif_param)

   con%private_info%x_tang = con%private_info%x_tang + c

   ! Check whether or not the continuation param, null vector, and eps
   ! are converged
   param_update = ABS(dt_p)  / (reltol * ABS(con%pitchfork_info%bif_param) + abstol)

   eps_update   = ABS(d_eps) / (reltol * ABS(eps) + abstol)

   ! get update norm of phi, to check for convergence */
   r_update = SUM( ( c / (ABS(con%private_info%x_tang)*reltol + abstol) )**2 )

   gnum_unks = gsum_double_conwrap(1.d0*con%general_info%numOwnedUnks)

   r_update = SQRT(gsum_double_conwrap(r_update)/gnum_unks)

   IF (RayQ == -1.d0) THEN
      IF (con%general_info%printproc > 1)  WRITE(*,*) 'Rayleigh Quotient Error: zero denom'
   ELSE
      IF (con%general_info%printproc > 4)  WRITE(*,*) 'Rayleigh Quotient before updates: ', RayQ
   END IF

   IF (con%general_info%printproc > 4) THEN
      WRITE(*,*) 'Pitchfork Continuation: Convergence Criteria'
      WRITE(*,*) 'Variable     Scaled Update (<1)  Unscaled Update  New Value'
      WRITE(*,*) '***********************************************************'
      WRITE(*,*) 'parameter    ', param_update, dt_p, con%pitchfork_info%bif_param
      WRITE(*,*) 'eps          ', eps_update, d_eps, eps
      WRITE(*,*) 'Null vector  ', r_update
      WRITE(*,*) '***********************************************************'
   END IF

   DEALLOCATE(a)
   DEALLOCATE(b)
   DEALLOCATE(c)
   DEALLOCATE(d)
   DEALLOCATE(e)
   DEALLOCATE(f)
   DEALLOCATE(x_tmp)
   DEALLOCATE(NULL)

   ! return convergence status of the parameter, eps and null vector */

   if (param_update < 1.d0 .AND. eps_update < 1.d0 .AND. r_update < 10.d0) THEN
      output = TRUE
   else
      output = FALSE
   END IF

END FUNCTION pitchfork_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION hopf_alg(x, delta_x, con, reltol, abstol) RESULT (output)

! Algorithm for locking on to and tracking a Hopf point.
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct)            :: con
   REAL(KIND=8)                :: reltol
   REAL(KIND=8)                :: abstol

   INTEGER :: output

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a, b, c, d, e, f, g, h, x_tmp, M_1, M_2
   REAL(KIND=8) :: dt_p
   INTEGER      :: i
   REAL(KIND=8) :: x_update, param_update, gnum_unks, tmp, RayQ, tmp2, alpha, beta
   REAL(KIND=8) :: omega_update, y_update, z_update, y_tmp, z_tmp
   INTEGER      :: first=TRUE, first2=TRUE
   REAL(KIND=8) :: test, count

   REAL(KIND=8) :: ltc, ltd, lte, ltf, ltg, lth, delta_omega, delta_y, delta_z
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL

   !************************** BEGIN EXECUTION *******************************
   first=TRUE; first2=TRUE

   ! Allocate arrays for hopf bordering algorithm */
   ALLOCATE(    a(0:con%general_info%numUnks-1))
   ALLOCATE(    b(0:con%general_info%numUnks-1))
   ALLOCATE(    c(0:con%general_info%numUnks-1))
   ALLOCATE(    d(0:con%general_info%numUnks-1))
   ALLOCATE(    e(0:con%general_info%numUnks-1))
   ALLOCATE(    f(0:con%general_info%numUnks-1))
   ALLOCATE(    g(0:con%general_info%numUnks-1))
   ALLOCATE(    h(0:con%general_info%numUnks-1))
   ALLOCATE(x_tmp(0:con%general_info%numUnks-1))
   ALLOCATE( NULL(0:con%general_info%numUnks-1))

   ! First time through this routine, check (1) that eigenvectors are
   ! non-zero and (2) the functional dependence of the Mass Matrix.
   IF (first == TRUE) THEN

      IF (con%general_info%printproc > 4) &
         WRITE(*,*) '   Entering Hopf Algorithm'

      tmp  = SQRT(dp(con%hopf_info%y_vec, con%hopf_info%y_vec))

      tmp2 = SQRT(dp(con%hopf_info%z_vec, con%hopf_info%z_vec))

      IF ((tmp == 0d0) .OR. (tmp2 == 0d0)) THEN

         IF (con%general_info%printproc > 1) &
            WRITE(*,*) 'ERROR in Hopf alg: complex vector pair must be supplied'

         output = -1
         RETURN
      END IF

      ! Test to see if the Mass Matrix is a function of the
      ! parameter or of the solution vector. Set mass_param and mass_x
      ! accordingly - 1=true, 0=false
      IF (con%hopf_info%mass_flag == 1) THEN

         con%private_info%mass_param = 1
         con%private_info%mass_x     = 1

      ELSE IF (con%hopf_info%mass_flag == 2) THEN

         con%private_info%mass_param = 0
         con%private_info%mass_x     = 1

      ELSE IF (con%hopf_info%mass_flag == 3) THEN

         con%private_info%mass_param = 1
         con%private_info%mass_x     = 0

      ELSE IF (con%hopf_info%mass_flag == 4) THEN

         con%private_info%mass_param = 0
         con%private_info%mass_x     = 0

      ELSE IF (con%hopf_info%mass_flag == 5) THEN

         ! Heuristic determination
         con%private_info%mass_param = 0
         con%private_info%mass_x     = 0

         IF (con%general_info%printproc > 7) &
            WRITE(*,*) 'Hopf Continuation: Determining  Mass Matrix Dependencies'

         CALL mass_matrix_fill_conwrap(x, x_tmp)
         CALL mass_matvec_mult_conwrap(con%hopf_info%y_vec, a)
         CALL matvec_mult_conwrap(con%hopf_info%z_vec, c)

         DO i = 0, con%general_info%numUnks-1
            x_tmp(i) = x(i) + scalar_perturbation(x(i), con%general_info%perturb)
         END DO

         CALL mass_matrix_fill_conwrap(x_tmp, f)
         CALL mass_matvec_mult_conwrap(con%hopf_info%y_vec, b)
         CALL matvec_mult_conwrap(con%hopf_info%z_vec, d)

         count = 0d0

         DO i = 0, con%general_info%numOwnedUnks-1
            IF (d(i)-c(i) /= 0d0) THEN
               count = count + 1d0
               e(i) = con%hopf_info%omega * (b(i)-a(i)) / (d(i)-c(i))
            ELSE
               e(i) = 0d0
            END IF
         END DO

         test = SQRT(dp(e, e)) / count

#if DEBUG > 2
         WRITE(*,*) 'heur: M/x = ', test
#endif
         IF (test > 0.01d0) &
            con%private_info%mass_x = 1

         dt_p = scalar_perturbation(con%hopf_info%bif_param, con%general_info%perturb)

         CALL assign_bif_parameter_conwrap(con%hopf_info%bif_param + dt_p)

         CALL mass_matrix_fill_conwrap(x, x_tmp)

         CALL mass_matvec_mult_conwrap(con%hopf_info%y_vec, b)

         CALL matvec_mult_conwrap(con%hopf_info%z_vec, d)

         CALL assign_bif_parameter_conwrap(con%hopf_info%bif_param)

         count = 0d0

         DO i = 0, con%general_info%numOwnedUnks - 1
            IF (d(i)-c(i) /= 0d0) THEN
               count = count + 1d0
               e(i)  = con%hopf_info%omega * (b(i)-a(i)) / (d(i)-c(i))
            ELSE
               e(i) = 0d0
            END IF
         END DO

         count = gsum_double_conwrap(count)
         test = SQRT(dp(e, e)) / count

#if DEBUG > 2
         WRITE(*,*) 'heur: M/param = ', test
#endif
         IF (test > 0.01d0) &
            con%private_info%mass_param = 1

         IF (con%general_info%printproc > 7) THEN
            IF(con%private_info%mass_x==1) THEN
               WRITE(*,*) 'Mass Matrix is a Function of Solution Vector!'
            ELSE
               WRITE(*,*) 'Mass Matrix is Independent of Solution Vector!'
            END IF
            IF(con%private_info%mass_param == 1) THEN
               WRITE(*,*) 'Mass Matrix is a Function of Bifurcation Parameter!'
            ELSE
               WRITE(*,*) 'Mass Matrix is Independent of Bifurcation Parameter!'
            END IF
         END IF

      ELSE

         IF (con%general_info%printproc > 1) &
               WRITE(*,*) 'ERROR: Mass Matrix Derivatives Flag not set in solve_continuation()'

         output = -1
         RETURN ! FixMe I added this, should this be here?

      END IF

      first = FALSE

   END IF ! IF (first)

   IF (con%general_info%printproc > 7) THEN
      IF(con%private_info%mass_x == 1) THEN
         WRITE(*,*) 'dM/dx is included in Komplex solves'
      ELSE
         WRITE(*,*) 'dM/dx is NOT included in Komplex solves'
      END IF
      IF(con%private_info%mass_param == 1) THEN
         WRITE(*,*) 'dM/d(param) is included in Komplex solves'
      ELSE
         WRITE(*,*) 'dM/d(param) is NOT included in Komplex solves'
      END IF
   END IF

   ! If Mass Matrix is not constant, allocate extra work arrays
   IF (con%private_info%mass_param == TRUE .OR. con%private_info%mass_x == TRUE) THEN
      ALLOCATE(M_1(0:con%general_info%numUnks-1))
      ALLOCATE(M_2(0:con%general_info%numUnks-1))
   END IF

   IF (first2 == TRUE) THEN
      ! calculate variable averages, to be used in ltransnorm calls below
      con%private_info%scale_vec = con%hopf_info%y_vec

      ! Make sure the eigenvectors y_vec and z_vec are orthogonal
      alpha = ltransnorm(con%hopf_info%y_vec, con%private_info%scale_vec)
      beta  = ltransnorm(con%hopf_info%z_vec, con%private_info%scale_vec)

      DO i=0, con%general_info%numOwnedUnks-1
         y_tmp = con%hopf_info%y_vec(i)
         z_tmp = con%hopf_info%z_vec(i)
         con%hopf_info%y_vec(i) =  (alpha * y_tmp + beta  * z_tmp)/(alpha*alpha + beta*beta)
         con%hopf_info%z_vec(i) = -(beta  * y_tmp - alpha * z_tmp)/(alpha*alpha + beta*beta)
      END DO

      first2 = FALSE
   END IF

   ! construct "a" vector from delta_x
   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Constructed *a* vector!'

   a = - delta_x

   ! Next, "b" is calculated. This is exactly
   ! the same linear solve as a continuation predictor step except
   ! the bif_param is perturbed, not the continuation param.
   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Calculating *b* vector!'

   CALL calc_rhs_continuation(TP_CONT_SOL2, x, b,                                      &
                              NULL, NULL, NULL,                                        &
                              con%hopf_info%bif_param, con%general_info%perturb, NULL, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(b, CHECK_JACOBIAN, NULL)

   ! Fill the Mass Matrix for RHS and Komplex Solves
   CALL mass_matrix_fill_conwrap(x, x_tmp)

   ! Next, "c" and "d" are calculated as a function of y_vec, z_vec,
   ! and the Mass Matrix M. c and d hold the rhs vector on input to
   ! solver and have solution on exit.
   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Calculating *c* and *d* vectors (Komplex solve)!'

   CALL calc_rhs_continuation(HP_CONT_SOL3, x, c,                                   &
                              con%hopf_info%z_vec, con%hopf_info%y_vec, x_tmp,      &
                              con%hopf_info%bif_param, con%general_info%perturb, d, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = komplex_linear_solver_conwrap(c, d, NEW_JACOBIAN, con%hopf_info%omega, x_tmp)

   ! Next, "e" and "f" are calculated as a function of a, y_vec, and z_vec.
   ! e and f hold the rhs vector on input to  solver and have solution
   ! on exit.
   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Calculating *e* and *f* vectors (Komplex solve)!'

   CALL calc_rhs_continuation(TP_CONT_SOL3, x, e,                                                     &
                              a, con%private_info%scale_vec, x_tmp,                                   &
                              con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%y_vec, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   CALL calc_rhs_continuation(TP_CONT_SOL3, x, f,                                                     &
                              a, con%private_info%scale_vec, x_tmp,                                   &
                              con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%z_vec, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   ! Get null vector residual now.
   RayQ = null_vector_resid(0d0, con%hopf_info%omega, con%hopf_info%y_vec, con%hopf_info%z_vec, TRUE)

   IF (con%private_info%mass_x == TRUE) THEN
      CALL calc_rhs_continuation(HP_CONT_DMDX, x, M_1,                                                   &
                                 a, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%z_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      CALL calc_rhs_continuation(HP_CONT_DMDX, x, M_2,                                                   &
                                 a, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%y_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      ! FixMe I edited this, check it
      !DO i = 0, con%general_info%numUnks-1
      !   e(i) = e(i) + M_1(i) * con%hopf_info%omega
      !   f(i) = f(i) - M_2(i) * con%hopf_info%omega
      !END DO
      e = e + M_1 * con%hopf_info%omega
      f = f - M_2 * con%hopf_info%omega

   END IF

   i = komplex_linear_solver_conwrap(e, f, OLD_JACOBIAN, con%hopf_info%omega, x_tmp)

   ! Next, "g" and "h" are calculated as a function of a, y_vec, and z_vec.
   ! g and h hold the rhs vector on input to  solver and have solution
   ! on exit.
   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Calculating *g* and *h* vectors  (Komplex solve)!'

   CALL calc_rhs_continuation(TP_CONT_SOL4, x, g,                                                     &
                              b, con%private_info%scale_vec, x_tmp,                                   &
                              con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%y_vec, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   CALL calc_rhs_continuation(TP_CONT_SOL4, x, h,                                                     &
                              b, con%private_info%scale_vec, x_tmp,                                   &
                              con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%z_vec, &
                              con%general_info%numUnks, con%general_info%numOwnedUnks)

   IF (con%private_info%mass_x == TRUE) THEN
      CALL calc_rhs_continuation(HP_CONT_DMDX, x, M_1,                                                   &
                                 b, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%z_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      CALL calc_rhs_continuation(HP_CONT_DMDX, x, M_2,                                                   &
                                 b, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%y_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      ! FixMe I edited this, check it
      !DO i = 0, con%general_info%numUnks-1
      !   g(i) = g(i) + M_1(i) * con%hopf_info%omega
      !   h(i) = h(i) - M_2(i) * con%hopf_info%omega
      !END DO
      g = g + M_1 * con%hopf_info%omega
      h = h - M_2 * con%hopf_info%omega

   END IF

   IF (con%private_info%mass_param == TRUE) THEN
      M_1 = x
      M_2 = x

      CALL calc_rhs_continuation(HP_CONT_DMDPARAM, x, M_1,                                               &
                                 a, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%z_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      CALL calc_rhs_continuation(HP_CONT_DMDPARAM, x, M_2,                                               &
                                 a, con%private_info%scale_vec, x_tmp,                                   &
                                 con%hopf_info%bif_param, con%general_info%perturb, con%hopf_info%y_vec, &
                                 con%general_info%numUnks, con%general_info%numOwnedUnks)

      ! FixMe I edited this, check it
      !DO i = 0, con%general_info%numUnks-1
      !  g(i) = g(i) + M_1(i) * con%hopf_info%omega
      !  h(i) = h(i) - M_2(i) * con%hopf_info%omega
      !END DO
      g = g + M_1 * con%hopf_info%omega
      h = h - M_2 * con%hopf_info%omega

   END IF

   i = komplex_linear_solver_conwrap(g, h, OLD_JACOBIAN_DESTROY, con%hopf_info%omega, x_tmp)

   IF (con%general_info%printproc > 4) &
      WRITE(*,*) '   Hopf Continuation: Finished Komplex solves,  Updating solution'

   ! Calculate the updates and scaled updates to bif_param,
   ! omega, y_vec, z_vec, and x.

   ! Calculate Update parameters
   ltc = ltransnorm(c, con%private_info%scale_vec)
   ltd = ltransnorm(d, con%private_info%scale_vec)
   lte = ltransnorm(e, con%private_info%scale_vec)
   ltf = ltransnorm(f, con%private_info%scale_vec)
   ltg = ltransnorm(g, con%private_info%scale_vec)
   lth = ltransnorm(h, con%private_info%scale_vec)

   ! Update turning point param
   dt_p = (ltc*ltf - lte*ltd + ltd) / (ltd*ltg - ltc*lth)
   con%hopf_info%bif_param = con%hopf_info%bif_param + dt_p
   CALL assign_bif_parameter_conwrap(con%hopf_info%bif_param)
   param_update = ABS(dt_p) / (reltol * ABS(con%hopf_info%bif_param) + abstol)

   ! Update imaginary eigenvalue omega
   delta_omega = (lth*dt_p + ltf) / ltd
   con%hopf_info%omega = con%hopf_info%omega + delta_omega
   omega_update = ABS(delta_omega) / (reltol * ABS(con%hopf_info%omega) + abstol)

   ! Update delta_x - state variables: Negative of update vector here */
   delta_x = - (a + dt_p*b)

   x_update = SUM(  ABS(delta_x) / (ABS(x)*reltol + abstol) &
                  * ABS(delta_x) / (ABS(x)*reltol + abstol) )

   ! Update eigenvectors for real and imaginary parts */
   y_update = 0d0
   z_update = 0d0

   DO i = 0, con%general_info%numOwnedUnks-1

      e(i) = - con%hopf_info%y_vec(i) + e(i) + g(i)*dt_p - c(i)*delta_omega
      f(i) = - con%hopf_info%z_vec(i) + f(i) + h(i)*dt_p - d(i)*delta_omega

      delta_y = e(i)
      delta_z = f(i)

      y_update = y_update + ABS(delta_y)/(ABS(con%hopf_info%y_vec(i))*reltol + abstol) &
                          * ABS(delta_y)/(ABS(con%hopf_info%y_vec(i))*reltol + abstol)
      z_update = z_update + ABS(delta_z)/(ABS(con%hopf_info%z_vec(i))*reltol + abstol) &
                          * ABS(delta_z)/(ABS(con%hopf_info%z_vec(i))*reltol + abstol)

      con%hopf_info%y_vec(i) = con%hopf_info%y_vec(i) + delta_y
      con%hopf_info%z_vec(i) = con%hopf_info%z_vec(i) + delta_z

   END DO

   gnum_unks = gsum_double_conwrap(1d0*con%general_info%numOwnedUnks)

   y_update  = SQRT(gsum_double_conwrap(y_update) / gnum_unks)
   z_update  = SQRT(gsum_double_conwrap(z_update) / gnum_unks)
   x_update  = SQRT(gsum_double_conwrap(x_update) / gnum_unks)

   ! Check whether or not continuation param, omega and eigenvectors
   ! are converged
   IF (RayQ == -1d0) THEN
      IF (con%general_info%printproc > 1) WRITE(*,*) '   Rayleigh Quotient Error: zero denom'
   ELSE
      IF (con%general_info%printproc > 4) WRITE(*,*) '   Rayleigh Quotient before updates: ', RayQ
   END IF

   IF (con%general_info%printproc > 4) THEN
      WRITE(*,*) '   Hopf Continuation: Convergence Criteria'
      WRITE(*,*) '   Variable     Scaled Update (<1)  Unscaled Update  New Value'
      WRITE(*,*) '   ***********************************************************'
      WRITE(*,*) '   X vector     ', x_update
      WRITE(*,*) '   parameter    ', param_update, dt_p, con%hopf_info%bif_param
      WRITE(*,*) '   omega        ', omega_update, delta_omega, con%hopf_info%omega
      WRITE(*,*) '   Y vector     ', y_update
      WRITE(*,*) '   Z vector     ', z_update
      WRITE(*,*) '   ***********************************************************'
   END IF

   DEALLOCATE(a)
   DEALLOCATE(b)
   DEALLOCATE(c)
   DEALLOCATE(d)
   DEALLOCATE(e)
   DEALLOCATE(f)
   DEALLOCATE(g)
   DEALLOCATE(h)
   DEALLOCATE(x_tmp)
   DEALLOCATE(NULL)
   IF (con%private_info%mass_param == TRUE .OR. con%private_info%mass_x == TRUE) THEN
      DEALLOCATE(M_1)
      DEALLOCATE(M_2)
   END IF

   ! return convergence status of the parameters and vectors
   IF ((x_update < 1.0) .AND. (param_update < 1.0) .AND. (omega_update < 1.0) &
                        .AND. (y_update < 100.0)   .AND. (z_update < 100.0)) THEN
      first2 = TRUE
      output = TRUE
   ELSE
      output = FALSE
   END IF

END FUNCTION hopf_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION phase_transition_alg(x, delta_x, con, reltol, abstol) RESULT (output)

! Algorithm for locking on to a phase transition.
! This involves finding two solutions that exist at the
! same parameters and at the same energy.
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(0:) ::x
   REAL(KIND=8), DIMENSION(0:) ::delta_x
   TYPE(con_struct) :: con
   REAL(KIND=8) ::reltol
   REAL(KIND=8) ::abstol
   INTEGER :: output


   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a, b, c, d
   REAL(KIND=8) :: g, dg_dtp, dt_p
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_tmp, x2_tmp
   REAL(KIND=8) :: dgac, dgbd, eps
   REAL(KIND=8) :: gnum_unks, x2_update, param_update
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL
   INTEGER :: i

   !************************** BEGIN EXECUTION *******************************
   ALLOCATE(     a(0:con%general_info%numUnks-1))
   ALLOCATE(     b(0:con%general_info%numUnks-1))
   ALLOCATE(     c(0:con%general_info%numUnks-1))
   ALLOCATE(     d(0:con%general_info%numUnks-1))
   ALLOCATE( x_tmp(0:con%general_info%numUnks-1))
   ALLOCATE(x2_tmp(0:con%general_info%numUnks-1))

   a = - delta_x

   ! Construct b vector using tangent calculation */
   CALL calc_rhs_continuation(TP_CONT_SOL2, x, b, NULL, NULL, NULL, &
                              con%phase_transition_info%bif_param, con%general_info%perturb, &
                              NULL, con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(b, OLD_JACOBIAN, NULL)

   ! Now do Newton iteration on second vector for c
   CALL matrix_residual_fill_conwrap(con%phase_transition_info%x2, c, RHS_MATRIX)

   i = linear_solver_conwrap(c, NEW_JACOBIAN, NULL)

   c = - c

   ! Construct d vector using tangent calculation
   CALL calc_rhs_continuation(TP_CONT_SOL2, con%phase_transition_info%x2, d, NULL, NULL, NULL, &
                              con%phase_transition_info%bif_param, con%general_info%perturb, &
                              NULL, con%general_info%numUnks, con%general_info%numOwnedUnks)

   i = linear_solver_conwrap(d, CHECK_JACOBIAN, NULL)

   ! Now start bordering alg for equal-energy constraint
   ! g is energy difference to be driven to zero,
   ! dg_dtp is the derivative of the energy w.r.t. the parameter
   ! dgac is the directional derivative in of g the a:c direction
   ! dgbd is the directional derivative of g in the b:d direction
   g = free_energy_diff_conwrap(x, con%phase_transition_info%x2)

   dt_p = scalar_perturbation(con%phase_transition_info%bif_param, con%general_info%perturb)

   CALL assign_bif_parameter_conwrap(con%phase_transition_info%bif_param + dt_p)

   dg_dtp = (free_energy_diff_conwrap(x, con%phase_transition_info%x2) - g ) / dt_p

   CALL assign_bif_parameter_conwrap(con%phase_transition_info%bif_param)

   eps = con%general_info%perturb * SQRT( dp(x,x) / (dp(a,a) + con%general_info%perturb) &
       + dp(con%phase_transition_info%x2,con%phase_transition_info%x2) &
         / (dp(c,c) + con%general_info%perturb))

   x_tmp  = x  + eps * a
   x2_tmp = con%phase_transition_info%x2 + eps * c

   dgac = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g ) / eps

   eps = con%general_info%perturb * SQRT(dp(x,x) / (dp(b,b) + con%general_info%perturb) &
       + dp(con%phase_transition_info%x2,con%phase_transition_info%x2) &
         / (dp(d,d) + con%general_info%perturb))

   x_tmp  = x  + eps * b
   x2_tmp = con%phase_transition_info%x2 + eps * d

   dgbd = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g ) / eps

   dt_p =  -(g  + dgac) / (dgbd + dg_dtp)

   ! update continuation parameter
   con%phase_transition_info%bif_param = con%phase_transition_info%bif_param + dt_p
   CALL assign_bif_parameter_conwrap(con%phase_transition_info%bif_param)

   ! update con%phase_transition_info%x2, checking for convergence */
   delta_x = c + dt_p * d

   ! get update norm to check for convergence

   ! FA: original code probably wrong
   ! x2_update = 0.d0
   ! for (i=0; i < con%general_info%numUnks; i++) {
   !  tmp      += delta_x[i]/(ABS(con%phase_transition_info%x2[i])*reltol + abstol);
   !           ^ this plus is very strange
   !  x2_update += tmp*tmp;
   ! }
   ! FA
   x2_update = SUM((delta_x/(ABS(con%phase_transition_info%x2)*reltol + abstol))**2)

   gnum_unks = gsum_double_conwrap(1.d0 * con%general_info%numOwnedUnks)

   x2_update  = SQRT( gsum_double_conwrap(x2_update/gnum_unks) )

   param_update = ABS(dt_p) &
                / (reltol * ABS(con%phase_transition_info%bif_param) + abstol)

   ! update con%phase_transition_info%x2 and modify delta_x for update in usual Newton routine */

   con%phase_transition_info%x2 = con%phase_transition_info%x2 + delta_x
   delta_x = - a - dt_p * b

   DEALLOCATE(a)
   DEALLOCATE(b)
   DEALLOCATE(c)
   DEALLOCATE(d)
   DEALLOCATE(x_tmp)
   DEALLOCATE(x2_tmp)

   IF (con%general_info%printproc > 4) THEN
      WRITE(*,*) 'Phase transition algorithm scaled parameter update  : ', param_update
      WRITE(*,*) 'Phase transition algorithm scaled solution #2 update: ', x2_update
   END IF

   ! x2_update and param_update must be < 1 for convergence */
   IF (x2_update < 1.d0 .AND. param_update < 1.d0) THEN
      output = TRUE
   ELSE
      output = FALSE
   END IF

END FUNCTION phase_transition_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION manifold_alg(x, delta_x, con, reltol, abstol) RESULT (converged)
! Bordering algorithm for multi-parameter arc-lenght continuation

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: delta_x
   TYPE(con_struct) :: con
   REAL(KIND=8) :: reltol
   REAL(KIND=8) :: abstol

   INTEGER :: converged

   converged = FALSE

!#ifdef LOCA_MF
!
!   extern int MFLOCANVNX(MFNVector)
!   extern int MFLOCANVNP(MFNVector)
!   extern double *MFLOCANVX(MFNVector)
!   extern double *MFLOCANVP(MFNVector)
!
!
!   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: a, rhs
!   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: b
!   INTEGER :: i, j
!   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: ux, u0x, up, u0p
!   REAL(KIND=8) :: param_update = 0
!   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: kbyk, rhsp, col_ix, col_ip
!   INTEGER,      DIMENSION(:),   ALLOCATABLE :: ipiv
!   INTEGER :: info = 0, one = 1
!   CHARACTER, PARAMETER :: cflag = 'N'
!
!   TYPE(MFNVector col)
!
!   ALLOCATE(  a(0:con%general_info%numUnks-1))
!   ALLOCATE(rhs(0:con%general_info%numUnks-1))
!
!   ALLOCATE(b(0:con%general_info%numUnks-1,0:con%manifold_info%k-1))
!
!   a = - delta_x
!
!   ux  = MFLOCANVX(con%manifold_info%u)
!   u0x = MFLOCANVX(con%manifold_info%u0)
!   up  = MFLOCANVP(con%manifold_info%u)
!   u0p = MFLOCANVP(con%manifold_info%u0)
!
!   ! Calculate dx/dp  by solving inv(J)*dR/dp for each parameter
!   CALL matrix_residual_fill_conwrap(x, rhs, RHS_ONLY)
!
!   DO i = 0, con%manifold_info%k-1
!
!      CALL calc_rhs_multi(up, i, con%manifold_info%k, b(:,i), rhs, x, &
!                          con%general_info%perturb, con%general_info%numUnks)
!
!      ! Solver for dx/dp tangent vector
!      i = linear_solver_conwrap(b(:,i), OLD_JACOBIAN, NULL)
!
!   END DO
!
!   ! Set up k-by-k matric problem for delta_p vector */
!
!   ALLOCATE(kbyk(0:(con%manifold_info%k * con%manifold_info%k)-1))
!   ALLOCATE(ipiv(0:con%manifold_info%k-1))
!   ALLOCATE(rhsp(0:con%manifold_info%k-1))
!
!   DO i = 0, con%manifold_info%k-1
!      col     = MFMColumn(con%manifold_info%phi, i)
!      col_ix  = MFLOCANVX(col)
!      col_ip  = MFLOCANVP(col)
!      ipiv(i) = 0.d0
!
!      ! scaled_dp(x,y) is used in MF code distance alg: needs to be the same */
!      rhsp(i) = - (scaled_dp(col_ix, ux) - scaled_dp(col_ix, u0x))
!
!      rhsp(i) = rhsp(i) - SUM(col_ip(j) * (up - u0p))
!
!      rhsp(i) = rhsp(i) - scaled_dp(col_ix, a)
!
!      DO j = 0, con%manifold_info%k - 1
!         kbyk(i + j*con%manifold_info%k) = scaled_dp(col_ix, b(:,j)) + col_ip(j)
!      END DO
!
!      CALL MFFreeNVector(col)
!
!   END DO
!
!   CALL  dgetrf(con%manifold_info%k, con%manifold_info%k, kbyk, con%manifold_info%k, ipiv, info)
!
!   IF (info==0) &
!      CALL dgetrs(cflag, con%manifold_info%k, one, kbyk, con%manifold_info%k, ipiv, rhsp, con%manifold_info%k, info)
!
!   ! Newton update for parameter vector p, and set -update in delta_x
!   param_update = SUM(ABS(rhsp) / (reltol * ABS(up) + abstol))
!   up = up + rhsp
!   delta_x = delta_x - MATMUL(b, rhsp)
!
!   CALL assign_multi_parameter_conwrap(up)
!
!   IF (param_update < 1.d0) THEN
!      converged = TRUE
!   ELSE
!      converged = FALSE
!   END IF
!
!   IF (con%general_info%printproc) THEN
!
!      WRITE(*,*) 'Manifold Continuation: Convergence Criteria'
!      WRITE(*,*) 'Variable     Scaled Update (<1)  Unscaled Update  New Value'
!      WRITE(*,*) '***********************************************************'
!      DO i = 0, con%manifold_info%k-1
!         WRITE(*,*) 'parameter(%d)    ', i, param_update, rhsp(i), up(i)
!      END DO
!      WRITE(*,*) '***********************************************************'
!   END IF
!
!   DEALLOCATE(kbyk)
!   DEALLOCATE(ipiv)
!   DEALLOCATE(rhsp)
!
!   DEALLOCATE(a)
!   DEALLOCATE(rhs)
!   DEALLOCATE(b)
!
!#else
  converged =  0
!#endif

END FUNCTION manifold_alg

!******************************************************************************
!******************************************************************************
!******************************************************************************

FUNCTION scalar_perturbation(param, perturb) RESULT (output)
! Routine for picking a perturbation to a scalar for simple forward
! difference derivatives.

   IMPLICIT NONE

   REAL(KIND=8) :: param
   REAL(KIND=8) :: perturb
   REAL(KIND=8) :: output

   output = perturb * (perturb + ABS(param))

END FUNCTION scalar_perturbation

!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE calc_rhs_continuation(rhs_type, x, resid_vector, ab_vec, scale_vec, x_tmp, &
			                           param, perturb, r_vec, numUnks, numOwnedUnks)

! routine to pre-calculate the non-standard right-hand-sides for solves
! of continuation runs. This routine returns the matrix fill time.
! The matrix at the current solution is also calculated.

   IMPLICIT NONE

   INTEGER                     :: rhs_type
   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8), DIMENSION(0:) :: resid_vector
   REAL(KIND=8), DIMENSION(0:) :: ab_vec
   REAL(KIND=8), DIMENSION(0:) :: scale_vec
   REAL(KIND=8), DIMENSION(0:) :: x_tmp
   REAL(KIND=8)                :: param
   REAL(KIND=8)                :: perturb
   REAL(KIND=8), DIMENSION(0:) :: r_vec
   INTEGER                     :: numUnks
   INTEGER                     :: numOwnedUnks

   REAL(KIND=8) :: abdp, dc_p, dc_p1
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: resid_delta

   ! Allocate and initialize resid_delta
   !CALL vec_init(resid_vector)
   resid_vector = 0d0
   ALLOCATE(resid_delta(0:numUnks-1))

   ! For the first 2 options, resid_delta is the right hand side at the
   ! perturbed value of the continuation parameter

   IF (rhs_type == CONT_TANGENT) THEN

      dc_p = scalar_perturbation(param, perturb)

      CALL assign_parameter_conwrap(param + dc_p)

      CALL matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY)

      CALL assign_parameter_conwrap(param)

      CALL matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX)

      resid_vector = (resid_delta - resid_vector)/dc_p

   ELSE IF (rhs_type == ARC_CONT_SOL2) THEN

      dc_p = scalar_perturbation(param, perturb)

      CALL assign_parameter_conwrap(param + dc_p)

      CALL matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY)

      CALL assign_parameter_conwrap(param)

      CALL matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX_SAVE)

      resid_vector = (resid_delta - resid_vector)/dc_p

   ! For TP_CONT_SOL2, resid_delta is the right hand side at the
   ! perturbed value of the turning point parameter
   ELSE IF (rhs_type == TP_CONT_SOL2 ) THEN

      dc_p = scalar_perturbation(param, perturb)

      CALL assign_bif_parameter_conwrap(param + dc_p)

      CALL matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY)

      CALL assign_bif_parameter_conwrap(param)

      CALL matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX_SAVE)

      resid_vector = -(resid_delta - resid_vector)/dc_p

   ! For TP_CONT_SOL3, resid_delta is the matrix at perturbed value of
   ! solution in the direction of ab_vec multiplied by the vector r_vec
   ELSE IF (rhs_type == TP_CONT_SOL3 ) THEN

      abdp = dp(ab_vec, ab_vec)

      IF (ABS(abdp) < 1d-99) THEN
        WRITE(*,*) 'Fatal error. ab_vec too small: ',abdp,', try a larger perturbation!'
        STOP ! FixMe
      ELSE
        dc_p1 = perturb * (perturb + SQRT(dp(x,x)/abdp))
      END IF

      x_tmp = x + dc_p1 * ab_vec

      CALL matrix_residual_fill_conwrap(x_tmp, resid_delta, RHS_MATRIX)

      CALL matvec_mult_conwrap(r_vec, resid_delta)

      CALL matrix_residual_fill_conwrap(x, resid_vector, RECOVER_MATRIX)

      CALL matvec_mult_conwrap(r_vec, resid_vector)

      resid_vector = - AGS_option * resid_vector - (resid_delta - resid_vector)/dc_p1

   ! For TP_CONT_SOL4, resid_delta is the matrix at perturbed value of
   ! solution in the direction of ab_vec multiplied by the vector r_vec
   ! plus the matrix perturbed in the direction of the bif_parameter
   ! multiplied by the vector r_vec.
   ELSE IF (rhs_type == TP_CONT_SOL4 ) THEN

      dc_p = scalar_perturbation(param, perturb)

      abdp = dp(ab_vec, ab_vec)

      IF (ABS(abdp) < 1d-99) THEN
        WRITE(*,*) 'Fatal error. ab_vec too small: ',abdp,', try a larger perturbation!'
        STOP ! FixMe
      ELSE
        dc_p1 = perturb * (perturb + sqrt(dp(x,x)/abdp))
      END IF

      x_tmp = x + dc_p1 * ab_vec

      CALL matrix_residual_fill_conwrap(x_tmp, resid_delta, MATRIX_ONLY)

      CALL matvec_mult_conwrap(r_vec, x_tmp)

      CALL assign_bif_parameter_conwrap(param + dc_p)

      CALL matrix_residual_fill_conwrap(x, resid_delta, MATRIX_ONLY)

      CALL assign_bif_parameter_conwrap(param)

      CALL matvec_mult_conwrap(r_vec, resid_delta)

      ! Two numerical differences with different perturbations are summed
      ! together
      resid_delta = resid_delta/dc_p + x_tmp/dc_p1

      ! also sum together the perturbations
      dc_p = 1d0 / (1d0/dc_p1 + 1d0/dc_p)

      CALL matrix_residual_fill_conwrap(x, resid_vector, RECOVER_MATRIX)

      CALL matvec_mult_conwrap(r_vec, resid_vector)

      resid_vector = - resid_delta + resid_vector/dc_p

   ELSE IF (rhs_type == HP_CONT_SOL3 ) THEN

      CALL mass_matvec_mult_conwrap(ab_vec, resid_vector)
      CALL mass_matvec_mult_conwrap(scale_vec, resid_delta)
      r_vec = - resid_delta

   ELSE IF (rhs_type == HP_CONT_DMDX ) THEN

      abdp = dp(ab_vec, ab_vec)

      IF (ABS(abdp) < 1d-99) THEN
        WRITE(*,*) 'Fatal error. ab_vec too small: ',abdp,', try a larger perturbation!'
        STOP ! FixMe
      ELSE
        dc_p1 = perturb * (perturb + SQRT(dp(x,x)/abdp))
      END IF

      x_tmp = x + dc_p1 * ab_vec

      CALL mass_matrix_fill_conwrap(x_tmp, resid_delta)
      CALL mass_matvec_mult_conwrap(r_vec, resid_delta)
      CALL mass_matrix_fill_conwrap(x, resid_vector)
      CALL mass_matvec_mult_conwrap(r_vec, resid_vector)

      resid_vector = -(resid_delta - resid_vector)/dc_p1

   ELSE IF (rhs_type == HP_CONT_DMDPARAM ) THEN

      dc_p = scalar_perturbation(param, perturb)

      CALL assign_bif_parameter_conwrap(param + dc_p)

      CALL mass_matrix_fill_conwrap(x, resid_vector)
      CALL mass_matvec_mult_conwrap(r_vec, resid_delta)

      CALL assign_bif_parameter_conwrap(param)

      CALL mass_matrix_fill_conwrap(x, resid_vector)
      CALL mass_matvec_mult_conwrap(r_vec, resid_vector)

      resid_vector = -(resid_delta - resid_vector)/dc_p

   END IF

   DEALLOCATE(resid_delta)

END SUBROUTINE calc_rhs_continuation

!******************************************************************************
!******************************************************************************
!******************************************************************************

!#ifdef LOCA_MF
!SUBROUTINE calc_rhs_multi(param_vec, i, k, phix, rhs, x, perturb, numUnks)
!
!   IMPLICIT NONE
!
!   REAL(KIND=8), DIMENSION(0:) ::param_vec
!   INTEGER :: i
!   INTEGER :: k
!   REAL(KIND=8), DIMENSION(0:) ::phix
!   REAL(KIND=8), DIMENSION(0:) ::rhs
!   REAL(KIND=8), DIMENSION(0:) ::x
!   REAL(KIND=8) ::perturb
!   INTEGER :: numUnks
!
!   REAL(KIND=8) ::dc_p, p_init
!   INTEGER :: j
!
!   p_init = param_vec(i)
!
!   dc_p = scalar_perturbation(param_vec(i), perturb)
!
!   param_vec(i) = param_vec(i) + dc_p
!   CALL assign_multi_parameter_conwrap(param_vec)
!
!!   CALL vec_init(phix)
!   phix = 0d0
!   CALL matrix_residual_fill_conwrap(x, phix, RHS_ONLY)
!
!   param_vec(i) = p_init
!   CALL assign_multi_parameter_conwrap(param_vec)
!
!   phix = - (phix - rhs) / dc_p
!
!END SUBROUTINE calc_rhs_multi
!#endif

!******************************************************************************
!******************************************************************************
!******************************************************************************

END MODULE Loca_bord
