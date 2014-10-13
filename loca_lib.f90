! -----------------------------------------------------------------------------
!   LOCA 1.0: Library of Continuation Algorithms

MODULE Loca_lib

   USE Loca_types       ! definition of all the types of variables needed
   USE Loca_util        ! initialize_util_routines
   USE Loca_bord        ! calc_rhs_continuation
   USE Loca_interface_o ! solution_output_conwrap

   USE Loca_wrappers    ! assign_parameter_conwrap, assign_bif_parameter_conwrap, perturb_solution_conwrap, linear_solver_conwrap
   USE newton


!   USE Loca_eigenvalue
!   USE Loca_mf ! loca_mf not implemented yet

   IMPLICIT NONE

!****************************************************************************
!****************************************************************************
!****************************************************************************

CONTAINS
!=======

FUNCTION con_lib(con) RESULT(output)

   IMPLICIT NONE
!****************************************************************************
!*
!*  Input Variables:
!*
!****************************************************************************
   TYPE(con_struct), POINTER :: con
   INTEGER                   :: output

   ! Local Variables
   !
   INTEGER :: n     ! Loop index
   INTEGER :: order ! Continuation order flag:
				  	     ! 0 - zero-order continuation
					     ! 1 - first-order continuation
					     ! 2 - arc-length continuation
				        ! This flag is always 0 on the first
					     ! solve, and 0 for turning point or any
					     ! other special continuation

   INTEGER :: num_newt_conv = 0 ! Number of newton iterations to reach
                                ! convergence for last nonlinear solve
                                ! -- used to pick next step size
				                    ! ALSO error flag, when < 0

   CHARACTER(LEN=7), PARAMETER :: yo = 'con_lib'

   INTEGER :: sn_old, sn_new ! Sign of con%private_info%dp_ds, to check for a
                             ! turning point

   INTEGER :: tan_flag ! Set to zero when tang_factor is smaller
                       ! than step limit specified: step will be
                       ! halved and repeated even though converged

   ! These quantities are used for arc length step control:
   REAL(KIND=8) :: delta_s, end_passed, max_step, temp_step
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_tang_old, x2_old, y_vec_old, z_vec_old
   REAL(KIND=8) :: bif_param_old
   REAL(KIND=8) :: ds_ratio = 1d0, first_arc_step = 0d0, arc_step_adj = 0d0
   REAL(KIND=8) :: tang_factor = 1d0
   REAL(KIND=8) :: step                ! step size of continuation parameter
   REAL(KIND=8) :: step_old            ! old step size of continuation parameter
   REAL(KIND=8) :: omega_old           ! old Hopf tracking frequency
   TYPE(arc_scale_struct) :: arc_scale ! allocate scaling params in struct
   REAL(KIND=8), DIMENSION(:), POINTER :: NULL
   INTEGER :: i


   !******************************* First Executable Statment *****************
   !
   WRITE(*,*) 'LOCA90 v0.1, suck it Sandia Corporation'

   ! Send vector length information to utilities file, so routines such
   ! as vector allocs and dot products won't need length information
   !
   CALL initialize_util_routines(con%general_info%numOwnedUnks, con%general_info%numUnks)

   ! Fork to MF library if doing multiparameter continuation
   IF (con%general_info%method == MANIFOLD_CONTINUATION) THEN
      WRITE(*,*) 'FATAL ERROR: mf library not implemented yet'
      STOP
      ! CALL loca_mf_driver(con)
      output = 1
   END IF

   ! Initialize arrays to store predicted and old solutions. Save the current
   ! solution, x, into con%private_info%x_old
   !
   ALLOCATE(con%private_info%x_old(0:con%general_info%numUnks-1))
   ALLOCATE(con%private_info%scale_vec(0:con%general_info%numUnks-1))
   ALLOCATE(NULL(0:con%general_info%numUnks-1))

   con%private_info%param_old = con%general_info%param
   step     = 0d0
   step_old = 0d0

   con%private_info%x_old = con%general_info%x

   SELECT CASE (con%general_info%method)

      CASE(ZERO_ORDER_CONTINUATION)
         ! FixMe this is not allocated in the C version, I don't know why Franco allocated it
         ALLOCATE(con%private_info%x_tang(0:con%general_info%numUnks-1))

      CASE(FIRST_ORDER_CONTINUATION)
         ALLOCATE(con%private_info%x_tang(0:con%general_info%numUnks-1))

      CASE(ARC_LENGTH_CONTINUATION)
         con%private_info%arc_step = 0d0
         ALLOCATE(con%private_info%x_tang(0:con%general_info%numUnks-1))
         ALLOCATE(x_tang_old(0:con%general_info%numUnks-1))

      CASE(TURNING_POINT_CONTINUATION)
         ! FixMe these are not allocated in the C version, I don't know why Franco allocated them
         ALLOCATE(con%private_info%x_tang(0:con%general_info%numUnks-1))
         ALLOCATE(x_tang_old(0:con%general_info%numUnks-1))

      CASE(PITCHFORK_CONTINUATION)
         ALLOCATE(con%private_info%x_tang(0:con%general_info%numUnks-1))
         ALLOCATE(x_tang_old(0:con%general_info%numUnks-1))

      CASE(HOPF_CONTINUATION)
         omega_old = con%hopf_info%omega
         ALLOCATE(y_vec_old(0:con%general_info%numUnks-1))
         ALLOCATE(z_vec_old(0:con%general_info%numUnks-1))

      CASE(PHASE_TRANSITION_CONTINUATION)
         ALLOCATE(x2_old(0:con%general_info%numUnks-1))

   END SELECT

   ! Initialize variables used in arc length step control

   IF (con%general_info%method == ARC_LENGTH_CONTINUATION) THEN

      arc_scale%dx_fac = 1d0

      IF (con%arclength_info%dp_ds2_goal < 1d-6) &
        arc_scale%dx_fac = 1d2

      arc_scale%dx0         = arc_scale%dx_fac
      arc_scale%dx_fac_max  = 1d8
      arc_scale%dp_ds_goal  = SQRT(con%arclength_info%dp_ds2_goal)
      arc_scale%dp_ds_limit = arc_scale%dx_fac_max * arc_scale%dp_ds_goal &
                              / SQRT(1d0 + con%arclength_info%dp_ds2_goal &
                              * (arc_scale%dx_fac_max * arc_scale%dx_fac_max - 1d0))

      IF (con%arclength_info%dp_ds_max < arc_scale%dp_ds_goal) &
         con%arclength_info%dp_ds_max = arc_scale%dp_ds_goal

   END IF

   ! Adjust the BCs/Properties/whatever that the con param really represents
   !
   CALL assign_parameter_conwrap(con%general_info%param)

   IF      (con%general_info%method == TURNING_POINT_CONTINUATION) THEN

      CALL assign_bif_parameter_conwrap(con%turning_point_info%bif_param)

   ELSE IF (con%general_info%method == PITCHFORK_CONTINUATION) THEN

      CALL assign_bif_parameter_conwrap(con%pitchfork_info%bif_param)

   ELSE IF (con%general_info%method == HOPF_CONTINUATION) THEN

      CALL assign_bif_parameter_conwrap(con%hopf_info%bif_param)

   ELSE IF (con%general_info%method == PHASE_TRANSITION_CONTINUATION) THEN

      CALL assign_bif_parameter_conwrap(con%phase_transition_info%bif_param)

   END IF

   ! In tp_continuation, perturb initial guess off of potential singularity
   !
   IF (con%general_info%method == TURNING_POINT_CONTINUATION .OR. &
       con%general_info%method == PITCHFORK_CONTINUATION) THEN

      IF (con%general_info%printproc > 4) &
         WRITE(*,*) yo, ': Adding random perturbation for continuation'

      CALL perturb_solution_conwrap(con%general_info%x, con%private_info%x_old, &
                                    con%private_info%scale_vec, con%general_info%numOwnedUnks)
   END IF

   ! Print out general time integration information

   IF (con%general_info%printproc > 1) THEN

      WRITE(*,*)
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(*,*) yo, ': Start Continuation'
      WRITE(*,*) 'Initial step size = ', con%stepping_info%first_step
      WRITE(*,*) 'Max number of continuation steps = ', con%stepping_info%max_steps
      WRITE(*,*) 'Max parameter value = ', con%stepping_info%max_param
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(*,*)

   END IF

   !***************************** CONTINUATION LOOP **************************/

   ! Initialize the time step counter to 0. Set order flag to zero-order
   ! continuation through first solution.

   order = 0
   con%private_info%step_num   = 0
   con%private_info%nstep      = con%stepping_info%base_step
   con%stepping_info%last_step = FALSE

   ! Loop through a number of continuation steps - note the loop index may not
   ! represent the actual step counter, due to failed steps.  The convention
   ! here is that solution 0 is not a step, so there will be
   ! con%stepping_info%max_steps+1 solutions, numbered 0 through
   ! con->stepping_info.max_steps.

   n = -1

   DO WHILE(n < con%stepping_info%max_steps)
      ! FixMe in the C version this is <= con%stepping_info%max_steps

      n = n + 1

      ! Print out an initial statement about the step.
      IF (con%general_info%printproc > 1) &
         CALL print_cont_step1(order, step, step_old, con)

      !  Set flag for new solve to detect first Newton iter

      con%private_info%first_iter = TRUE

      ! Solve the system of equations at the current step.
      ! Note - con%general_info%x is considered to be updated, on return from this
      ! solution.
      num_newt_conv = nonlinear_solver_conwrap(con%general_info%x, con%private_info%step_num, &
                                               con%general_info%param, step, con)

      ! If tan_factor changes too much, tan_flag tells to halve step & reset.

      tan_flag = TRUE

      ! Check for convergence

      IF (num_newt_conv < 0) THEN

         ! Convergence Failure!
         !
         ! If initial guess did not converge, abort entire continuation run

         IF (con%private_info%step_num == 0) THEN

            n = con%stepping_info%max_steps ! Force IO and exit
            con%private_info%step_num = - 1 ! Set failure flag

            IF (con%general_info%printproc > 1) THEN
               WRITE(*,*) yo, ': INITIAL GUESS DID NOT CONVERGE'
               WRITE(*,*) 'ABORTING CONTINUATION RUN'
            END IF

         ! If this convergence failure wasn't the first or last step, cut step
         ! size in half and calculate a new initial guess.  New guess is
         ! the old solution plus the tangent to the previous prediction
         ! times the halved step size.
         !
         ELSE

            IF (n < con%stepping_info%max_steps) THEN ! This is not the last step

               step = step * 0.5d0

               ! Check for step size too small, abort if below min_delta_p
               IF (ABS(step) < con%stepping_info%min_delta_p) THEN

                  n = con%stepping_info%max_steps ! Force IO and exit
                  con%private_info%step_num = - 1 ! Set failure flag

                  IF (con%general_info%printproc > 1) THEN
                     WRITE(*,*) yo, ': CONTINUATION STEP SIZE TOO SMALL'
                     WRITE(*,*) 'ABORTING CONTINUATION RUN'
                  END IF
               END IF

               con%general_info%param = con%private_info%param_old + step

               CALL assign_parameter_conwrap(con%general_info%param)

               SELECT CASE (order)

                  CASE(0)
                     con%general_info%x = con%private_info%x_old

                     SELECT CASE (con%general_info%method)

                        CASE(TURNING_POINT_CONTINUATION)
                           con%private_info%x_tang = x_tang_old
                           CALL assign_bif_parameter_conwrap(bif_param_old)
                           con%turning_point_info%bif_param = bif_param_old

                        CASE(PITCHFORK_CONTINUATION)
                           con%private_info%x_tang = x_tang_old
                           CALL assign_bif_parameter_conwrap(bif_param_old)
                           con%pitchfork_info%bif_param = bif_param_old

                        CASE(PHASE_TRANSITION_CONTINUATION)
                           con%phase_transition_info%x2 = x2_old
                           CALL assign_bif_parameter_conwrap(bif_param_old)
                           con%phase_transition_info%bif_param = bif_param_old

                        CASE(HOPF_CONTINUATION)
                           con%hopf_info%y_vec = y_vec_old
                           con%hopf_info%z_vec = z_vec_old
                           CALL assign_bif_parameter_conwrap(bif_param_old)
                           con%hopf_info%bif_param = bif_param_old
                           con%hopf_info%omega = omega_old

                       END SELECT

                  CASE(1)
                     con%general_info%x =  con%private_info%x_old - step * con%private_info%x_tang

                  CASE(2)
                     con%private_info%arc_step = 0.5d0 * con%private_info%arc_step
                     con%general_info%x =  con%private_info%x_old + con%private_info%arc_step * con%private_info%x_tang

               END SELECT

               ! Also, reset last_step flag to continue to final parameter value.
               con%stepping_info%last_step = FALSE

            ELSE ! If it was the last step, however, reset to previous solution.
               con%general_info%param = con%private_info%param_old
               CALL assign_parameter_conwrap(con%general_info%param)
               step = 0d0
               con%private_info%x_old = con%general_info%x
               con%stepping_info%last_step = TRUE
            END IF ! (n < con%stepping_info%max_steps)

            ! Print out failure message
            IF (con%general_info%printproc > 1) &
               CALL print_cont_step_fail(order, step, con)

         END IF ! (con%private_info%step_num == 0)

      ELSE ! Solver did Converge!!

         ! Check to see if parameter value passed end value (either direction)
         end_passed = (con%general_info%param - con%stepping_info%max_param) &
                    * (con%private_info%param_old - con%stepping_info%max_param)

         IF (con%private_info%step_num /= 0 .AND. end_passed <= 0d0) &
            con%stepping_info%last_step = TRUE

         ! For arc-length continuation using the tangent factor,
         ! calculate con%private_info%x_tang and new tang_factor in advance,
         ! then compare to value for previous step.  If this
         ! difference exceeds max_tang_step, treat this as a failed
         ! step and halve arc_step as above.

         IF (con%stepping_info%last_step == FALSE .AND. con%general_info%method == ARC_LENGTH_CONTINUATION) THEN

            IF (con%general_info%printproc > 4) THEN
               WRITE(*,*) '   Doing Pseudo Arc-length continuation --'
               WRITE(*,*) '   Calculating tangent vector by one linear solve'
            END IF

            IF (con%private_info%step_num == 0) &
               step = con%stepping_info%first_step

            CALL calc_rhs_continuation(CONT_TANGENT, con%general_info%x, con%private_info%x_tang, &
                                       NULL, NULL, NULL,                                          &
                                       con%general_info%param, con%general_info%perturb, NULL,    &
                                       con%general_info%numUnks, con%general_info%numOwnedUnks)

            i = linear_solver_conwrap(con%private_info%x_tang, NEW_JACOBIAN, NULL)

            IF (con%private_info%step_num > 0) THEN
               tang_factor =  ABS(scaled_dot_prod(x_tang_old, con%private_info%x_tang, &
                                                  con%private_info%scale_vec, con%general_info%numOwnedUnks)) &
                           / SQRT(scaled_dot_prod(x_tang_old, x_tang_old, &
                                                  con%private_info%scale_vec, con%general_info%numOwnedUnks) &
                                * scaled_dot_prod(con%private_info%x_tang, con%private_info%x_tang, &
                                                  con%private_info%scale_vec, con%general_info%numOwnedUnks))
            END IF

            tang_factor = EXP(con%arclength_info%tang_exp * LOG(tang_factor))

            IF (con%general_info%printproc > 7) &
               WRITE(*,*) ' Tangent factor is ', tang_factor

            ! Repeat step if tang_factor is too small

            IF (con%private_info%step_num > 1 .AND. tang_factor < con%arclength_info%tang_step_limit) THEN

               IF (con%general_info%printproc > 7) &
                  WRITE(*,*) ' Step limit exceeded: Retrying ...'
               
               tan_flag = FALSE
               step = 0.5d0 * step
               con%general_info%param = con%private_info%param_old + step
               
               CALL assign_parameter_conwrap(con%general_info%param)
               
               con%private_info%arc_step = 0.5d0 * con%private_info%arc_step
               
               con%general_info%x =  con%private_info%x_old + con%private_info%arc_step * x_tang_old
               
               con%private_info%x_tang = x_tang_old

            END IF

         END IF

         ! If tan_flag has not been set to zero, proceed with continuation

         IF (tan_flag == TRUE) THEN

            ! Print out final results of a successful time step

            IF (order == 2) &
               step = con%general_info%param - con%private_info%param_old

            IF (con%general_info%printproc > 4) &
               CALL print_cont_step2(order, step, con)

            ! If first continuation step, set to value from input file.  If
            ! controlled steps, use # of Newton iters to pick next step size
            ! For arc-length continuation, impose maximum step size as
            ! approximated by arc_step and dp_ds.
            ! Note:  without time step control, step size can never increase.

            step_old = step

            IF (con%private_info%step_num == 0) THEN

               step = con%stepping_info%first_step

            ELSE

               ! normal step control
               IF (con%stepping_info%last_step == FALSE .AND. con%stepping_info%step_ctrl > 0d0) THEN

                  IF (order == 2) THEN
                     con%private_info%arc_step = con%private_info%arc_step &
                                               * simple_step_control(num_newt_conv,                    &
                                                                     con%stepping_info%max_newton_its, &
                                                                     con%stepping_info%step_ctrl)
                  ELSE
                     step = step * simple_step_control(num_newt_conv,                    &
                                                       con%stepping_info%max_newton_its, &
                                                       con%stepping_info%step_ctrl)

                     IF (ABS(step) > con%stepping_info%max_delta_p) &
                        step = SIGN(1d0,step) * con%stepping_info%max_delta_p

                  END IF

               ! for constant step runs where the step has been cut, let it
               ! increase again with step control of 0.5
               ELSE IF (order < 2 .AND. ABS(step) < ABS(con%stepping_info%first_step)) THEN

                  step = step * simple_step_control(num_newt_conv,                    &
                                                    con%stepping_info%max_newton_its, &
                                                    0.5d0)

               ELSE IF (order == 2 .AND. ABS(con%private_info%arc_step) < arc_step_adj) THEN

                  con%private_info%arc_step = con%private_info%arc_step &
                                            * simple_step_control(num_newt_conv,                    &
                                                                  con%stepping_info%max_newton_its, &
                                                                  0.5d0)

                  IF (con%private_info%arc_step > arc_step_adj) &
                     con%private_info%arc_step = arc_step_adj

               END IF

            END IF


            IF (con%general_info%method == ARC_LENGTH_CONTINUATION) THEN
               delta_s = con%private_info%arc_step
            ELSE
               delta_s = step
            END IF

            ! Output information at the end of every successful time step
            ! Depending on the solution method, there can be up to three
            ! solution vectors and parameter values to write out at a solution
            !
            SELECT CASE (con%general_info%method)

               CASE(ZERO_ORDER_CONTINUATION)
                  CALL solution_output_conwrap(1, con%general_info%x, con%general_info%param, NULL, &
                                               0d0, NULL, 0d0,                                      &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(FIRST_ORDER_CONTINUATION)
                  CALL solution_output_conwrap(1, con%general_info%x, con%general_info%param, NULL, &
                                               0d0, NULL, 0d0,                                      &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(ARC_LENGTH_CONTINUATION)
                  CALL solution_output_conwrap(1, con%general_info%x, con%general_info%param, NULL, &
                                               0d0, NULL, 0d0,                                      &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(TURNING_POINT_CONTINUATION)
                  CALL solution_output_conwrap(2, con%general_info%x, con%general_info%param, con%private_info%x_tang, &
                                               con%turning_point_info%bif_param, NULL, 0d0,                            &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(PITCHFORK_CONTINUATION)
                  CALL solution_output_conwrap(2, con%general_info%x, con%general_info%param, con%private_info%x_tang, &
                                               con%pitchfork_info%bif_param, NULL, 0d0,                                &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(HOPF_CONTINUATION)
                  CALL solution_output_conwrap(3, con%general_info%x, con%general_info%param, con%hopf_info%y_vec, &
                                               con%hopf_info%bif_param, con%hopf_info%z_vec, con%hopf_info%omega,  &
                                               con%private_info%step_num, num_newt_conv, con)

               CASE(PHASE_TRANSITION_CONTINUATION)
                  CALL solution_output_conwrap(2, con%general_info%x, con%general_info%param, con%phase_transition_info%x2, &
                                               con%phase_transition_info%bif_param, NULL, 0d0,                              &
                                               con%private_info%step_num, num_newt_conv, con)

            END SELECT

            ! Check current parameter value against the maximum.
            !
            IF (con%stepping_info%last_step == TRUE) THEN

               n = con%stepping_info%max_steps ! Force IO and exit

               IF (con%general_info%printproc > 1) THEN
                  WRITE(*,*) 'Simulation completed continuation in ', con%private_info%nstep,' steps'
                  WRITE(*,*) 'Final Parameter Value: ', con%general_info%param
                  WRITE(*,*)
               END IF

            END IF

            IF (n < con%stepping_info%max_steps) THEN

               ! Finally, its time to do some continuation, since the previous step
               ! converged and it wasn't the last step.

               IF (con%general_info%printproc > 4) THEN
                  WRITE(*,*)
                  WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                  WRITE(*,*) 'Calculating initial guess for next continuation step'
               END IF

               ! Set the continuation order 0, 1 (Euler-Newton), or 2 (Arc-length)
               ! Always do only zero order continuation of a turning point
               ! or any other special continuation.

               SELECT CASE (con%general_info%method)
                 CASE(FIRST_ORDER_CONTINUATION)
                    order = 1
                 CASE(ARC_LENGTH_CONTINUATION)
                    order = 2
                 CASE DEFAULT
                    order = 0
               END SELECT

               ! Possibly adjust the step value for this step so it hits maximum
               ! value exactly (zero or first order continuation).
               ! For arc length continuation, this is done after solution scaling.

               IF (order /= 2 .AND. con%private_info%step_num /= 0) THEN

                  temp_step = delta_s
                  end_passed = (con%general_info%param + temp_step - con%stepping_info%max_param) &
                             * (con%general_info%param - con%stepping_info%max_param)

                  ! If end_passed <= 0, next step would take param past end value

                  IF (end_passed <= 0) THEN

                     step = con%stepping_info%max_param - con%general_info%param
                     con%stepping_info%last_step = TRUE

                     IF (con%general_info%printproc > 7) &
                        WRITE(*,*) '******** LAST PATH STEP!'

                  END IF

               END IF

               ! Calculate the tangent to the solution, con%private_info%x_tang, for the
               ! current step. This is trivial for 0-order continuation, requires
               ! 1 linear solve for 1st order continuation and for arc-length
               ! continuation. Use tangent to predict new solution in con%general_info%x.

               SELECT CASE (order)

                  CASE(0)
                     IF (con%general_info%printproc > 4) THEN
                        WRITE(*,*) '   Doing Zeroth-order continuation --'
                        WRITE(*,*) '   Previous solution used as initial guess'
                     END IF

                     ! NO definition of con%private_info%x_tang needed for zero order.
                     ! Don't set it to zero because that will mess up the
                     ! turning point and phase transition tracking algorithms,
                     ! which use that space for other purposes.

                     ! Save the old solution, before overwriting with the new solution
                     ! for use in restarting after failed steps

                     con%private_info%x_old = con%general_info%x

                     SELECT CASE (con%general_info%method)

                        CASE(TURNING_POINT_CONTINUATION)
                           x_tang_old    = con%private_info%x_tang
                           bif_param_old = con%turning_point_info%bif_param

                        CASE(PITCHFORK_CONTINUATION)
                           x_tang_old    = con%private_info%x_tang
                           bif_param_old = con%pitchfork_info%bif_param

                        CASE(PHASE_TRANSITION_CONTINUATION)
                           x2_old        = con%phase_transition_info%x2
                           bif_param_old = con%phase_transition_info%bif_param

                        CASE(HOPF_CONTINUATION)
                           y_vec_old     = con%hopf_info%y_vec
                           z_vec_old     = con%hopf_info%z_vec
                           bif_param_old = con%hopf_info%bif_param
                           omega_old     = con%hopf_info%omega

                     END SELECT

                     ! perturb guess off of singularity in tp_continuation */
                     !
                     IF (con%general_info%method == TURNING_POINT_CONTINUATION .OR. &
                         con%general_info%method == PITCHFORK_CONTINUATION) THEN

                        IF (con%general_info%printproc > 4) &
                           WRITE(*,*) 'con_lib: Adding random perturbation for continuation'

                        CALL perturb_solution_conwrap(con%general_info%x, con%private_info%x_old, &
                                                      con%private_info%scale_vec, con%general_info%numOwnedUnks)
                     END IF

                  CASE(1)
                     IF (con%general_info%printproc > 4) THEN
                       WRITE(*,*) '   Doing First-order continuation --'
                       WRITE(*,*) '   Calculating tangent vector by one linear solve'
                     END IF

                     ! Choose perturbation for numerical derivative of Residuals w.r.t
                     ! continuation parameter, and solve for the tangent vector as in
                     ! eq. 7.13 in John's thesis.  The continuation parameter and
                     ! perturbation, con%general_info%param and delta_param, are passed to the
                     ! linear solver in the spots for the time and CJ, which aren't'
                     ! needed for steady problems.

                     CALL calc_rhs_continuation(CONT_TANGENT, con%general_info%x, con%private_info%x_tang, &
                                                NULL, NULL, NULL,                                          &
                                                con%general_info%param, con%general_info%perturb, NULL,    &
                                                con%general_info%numUnks, con%general_info%numOwnedUnks)

                     i = linear_solver_conwrap(con%private_info%x_tang, NEW_JACOBIAN, NULL)

                     ! Save the old solution, before overwriting with the new solution

                     con%private_info%x_old = con%general_info%x

                     ! Multiply the tangent vector, con%private_info%x_tang initially, by the step
                     ! length, and add to con%general_info%x, to obtain an initial guess at
                     ! the next parameter value.

                     con%general_info%x = con%general_info%x - con%private_info%x_tang * step

                  CASE(2) ! Arclength ccontinuation predictor and step control

                     ! con%private_info%x_tang vector found above.

                     ! Function "solution_scale" rescales solution as necessary.
                     !
                     ds_ratio = solution_scale(con, arc_scale)

                     IF (con%general_info%printproc > 7) &
                        WRITE(*,*) 'Solution scale factor is ', arc_scale%dx_fac

                     ! Adjust arc_step for current scale factor.
                     con%private_info%arc_step = con%private_info%arc_step / ds_ratio
                     arc_step_adj = ABS(first_arc_step * arc_scale%ds_fac)

                     ! Adjust arc_step for current tangent factor.
                     con%private_info%arc_step = con%private_info%arc_step * tang_factor

                     ! Also reduce arc_step if too large.
                     max_step = ABS(con%stepping_info%max_delta_p / con%private_info%dp_ds)

                     IF (con%private_info%arc_step > max_step) &
                        con%private_info%arc_step = max_step

                     ! Readjust the step value for this step so it hits maximum
                     ! or end value (approximately) for arc length continuation.

                     IF (con%private_info%step_num /= 0) THEN

                        temp_step = delta_s * con%private_info%dp_ds
                        IF (step < 0) temp_step = -1d0 * temp_step
                        end_passed = (con%general_info%param + temp_step - con%stepping_info%max_param) &
                                   * (con%general_info%param - con%stepping_info%max_param)

                        ! If end_passed < 0, next step would take param past end value

                        IF (end_passed < 0) THEN
                           temp_step = con%stepping_info%max_param - con%general_info%param
                           IF (step < 0) temp_step = -1d0 * temp_step
                           con%private_info%arc_step = ABS(temp_step / con%private_info%dp_ds)
                           con%stepping_info%last_step = TRUE
                           IF (con%general_info%printproc > 7) WRITE(*,*) ' ******** LAST PATH STEP!'
                        END IF

                     END IF

                     ! If this is the first step, pick con%private_info%arc_step so that this step
                     ! will progress the parameter by approximately con%stepping_info%step
                     IF (con%private_info%step_num == 0) THEN

                        IF (step < 0) THEN
                           con%private_info%dp_ds = - 1d0 * con%private_info%dp_ds
                           sn_old = -1
                        ELSE
                           sn_old = 1
                        END IF

                        con%private_info%arc_step = step / con%private_info%dp_ds
                        first_arc_step = con%private_info%arc_step

                     ELSE

                        ! Pick sign of con%private_info%dp_ds according to eq. 7.14b in JNS thesis
                        ! NOTE: -1.0 factor multiplying solution terms is because
                        ! con%private_info%x_tang is currently the negative of the tangent vector.
                        !
                        ! and check if a turning point was passed --
                        IF (-1d0 * (scaled_dot_prod(con%private_info%x_tang, con%general_info%x, &
                                                    con%private_info%scale_vec, con%general_info%numOwnedUnks) &
                                  - scaled_dot_prod(con%private_info%x_tang, con%private_info%x_old, &
                                                    con%private_info%scale_vec, con%general_info%numOwnedUnks)) &
                            + con%general_info%param - con%private_info%param_old < 0d0) THEN

                           con%private_info%dp_ds = -1d0 * con%private_info%dp_ds
                           sn_new = -1
                        ELSE
                           sn_new = 1
                        END IF

                        IF ((con%general_info%printproc > 1) .AND. sn_old /= sn_new) &
                           WRITE(*,*) 'A turning point was passed !!!!!!!'

                        sn_old = sn_new

                     END IF

                     ! Save the old solution, before overwriting with the new solution
                     con%private_info%x_old = con%general_info%x

                     ! Calculate prediction for next step from Eqs. 7.15&7.16 in JNS
                     ! thesis (leaving con%private_info%x_tang = u_dot).
                     con%private_info%x_tang  = - con%private_info%x_tang * con%private_info%dp_ds
                     con%general_info%x = con%general_info%x + con%private_info%x_tang * con%private_info%arc_step

                     step       = con%private_info%dp_ds * con%private_info%arc_step
                     x_tang_old = con%private_info%x_tang

               END SELECT ! SELECT CASE (order)

               ! Increment the continuation parameter.  Update the
               ! BCs/Properties/whatever that the continuation parameter really
               ! represents.
               con%private_info%param_old = con%general_info%param
               con%general_info%param = con%general_info%param + step
               CALL assign_parameter_conwrap(con%general_info%param)

               ! Increment the step counter. Print final message.
               con%private_info%step_num = con%private_info%step_num + 1
               con%private_info%nstep    = con%private_info%nstep + 1

            END IF ! of:  if (n < con%stepping_info%max_steps)

         END IF ! of:  if (tan_flag == TRUE)

      END IF ! of else section for converged solves

   END DO ! loop over continuation step attempts --- for (n = 0; ... ---

   !*********************CLEAN-UP AREA*****************************************/

   ! Free auxiliary vectors no matter what happened

   DEALLOCATE (con%private_info%x_old)
   DEALLOCATE (con%private_info%scale_vec)
   DEALLOCATE (NULL)

   SELECT CASE (con%general_info%method)

      CASE(ZERO_ORDER_CONTINUATION)
          DEALLOCATE (con%private_info%x_tang)

      CASE(FIRST_ORDER_CONTINUATION)
          DEALLOCATE (con%private_info%x_tang)

      CASE(ARC_LENGTH_CONTINUATION)
          IF (con%general_info%nv_save == FALSE) DEALLOCATE (con%private_info%x_tang)
          DEALLOCATE (x_tang_old)

      CASE(TURNING_POINT_CONTINUATION)
          IF (con%general_info%nv_save == FALSE) DEALLOCATE (con%private_info%x_tang)
          DEALLOCATE (x_tang_old)

      CASE(PITCHFORK_CONTINUATION)
          IF (con%general_info%nv_save == FALSE) DEALLOCATE (con%private_info%x_tang)
          DEALLOCATE (x_tang_old)

      CASE(PHASE_TRANSITION_CONTINUATION)
          DEALLOCATE (x2_old)

      CASE(HOPF_CONTINUATION)
          DEALLOCATE (y_vec_old)
          DEALLOCATE (z_vec_old)

   END SELECT

   ! Send back the overall result of the time step

   output = con%private_info%step_num

END FUNCTION con_lib

!***************************************************************************
!***************************************************************************
!***************************************************************************

FUNCTION solution_scale(con, arc) RESULT(output)

   IMPLICIT NONE

   TYPE(con_struct)       :: con
   TYPE(arc_scale_struct) :: arc
   REAL(KIND=8) :: output, ds_ratio = 1d0, umag

   ! Calculate average of each variable for scaling of arc-length eq
   
   CALL calc_scale_vec_conwrap(con%general_info%x, con%private_info%scale_vec, con%general_info%numUnks)
   
   ! Get dot product before solution scaling.
   
   umag = arc%dx0 * arc%dx0 &
        * scaled_dot_prod(con%private_info%x_tang, con%private_info%x_tang, &
                          con%private_info%scale_vec, con%general_info%numOwnedUnks)
   
   ! Adjust con%private_info%scale_vec by current scale factor arc%dx_fac.
   
   con%private_info%scale_vec = con%private_info%scale_vec * arc%dx_fac
   
   arc%umag2 = scaled_dot_prod (con%private_info%x_tang, con%private_info%x_tang, &
                                con%private_info%scale_vec, con%general_info%numOwnedUnks)
   
   ! Calculate deriv of parameter w.r.t arc length, Eq.7.14a in JNS
   ! thesis.  Save actual value in arc%dp_ds_act.
   
   con%private_info%dp_ds = 1d0 / SQRT(1d0 + arc%umag2)
   
   arc%dp_ds_old = con%private_info%dp_ds
   
   ! On the first step, set arc%dx_fac to the value which will correspond
   ! to the desired dp_ds value.
   
   IF (con%private_info%step_num == 0 .AND. con%arclength_info%dp_ds2_goal > 1d-6) THEN
   
      arc%dx_fac = con%private_info%dp_ds                                              &
                   * SQRT(   (1d0 - con%arclength_info%dp_ds2_goal)                    &
                           / (1d0 - con%private_info%dp_ds * con%private_info%dp_ds) ) &
                   / arc%dp_ds_goal
      
      arc%dx0 = arc%dx_fac
      
      con%private_info%scale_vec = con%private_info%scale_vec * arc%dx_fac
      
      arc%umag2 = scaled_dot_prod(con%private_info%x_tang, con%private_info%x_tang, &
                                  con%private_info%scale_vec, con%general_info%numOwnedUnks)
      
      con%private_info%dp_ds = 1d0 / SQRT(1d0 + arc%umag2)
      
      ds_ratio = con%private_info%dp_ds / arc%dp_ds_old
   
   END IF
   
   ! If dp_ds is too large, calculate factor to adjust arc_step
   ! such that (dp)^2 / (ds)^2 doesn't exceed dp_ds2_goal.
   
   IF (con%private_info%step_num > 0 .AND. con%arclength_info%dp_ds2_goal > 1d-6 &
        .AND. con%private_info%dp_ds > con%arclength_info%dp_ds_max) THEN
   
   
      IF(con%general_info%printproc > 7) &
         WRITE(*,*)  'dp_ds out of limits at',  con%private_info%dp_ds, ' Rescaling...'
   
      ! Calculate scale factor for con%private_info%scale_vec (arc%dx_fac).
      ! Limit arc%dx_fac to arc%dx_fac_max to avoid division by zero.
      arc%dx_fac_old = arc%dx_fac
   
      IF(con%private_info%dp_ds > arc%dp_ds_limit) THEN
   
         arc%dx_fac = arc%dx_fac_max
   
      ELSE

         arc%dx_fac = arc%dx_fac * con%private_info%dp_ds                                 &
                      * SQRT(   (1d0 - con%arclength_info%dp_ds2_goal)                    &
                              / (1d0 - con%private_info%dp_ds * con%private_info%dp_ds) ) &
                      / arc%dp_ds_goal
      END IF
   
      ! Multiply con%private_info%scale_vec through by arc%dx_fac.
      con%private_info%scale_vec = con%private_info%scale_vec * arc%dx_fac / arc%dx_fac_old
   
      ! Recalculate unknown dot product (arc%umag2) and dp_ds.
      arc%umag2 = scaled_dot_prod(con%private_info%x_tang, con%private_info%x_tang, &
                                  con%private_info%scale_vec, con%general_info%numOwnedUnks)
   
      con%private_info%dp_ds = 1d0 / SQRT(1d0 + arc%umag2)
   
      ds_ratio = con%private_info%dp_ds / arc%dp_ds_old
   
   END IF
   
   ! Get arc%ds_fac
   
   arc%ds_fac = 1d0 / ( con%private_info%dp_ds * SQRT(1d0 + umag) )
   
   ! Return ds_ratio
   
   output = ds_ratio

END FUNCTION solution_scale

!***************************************************************************
!***************************************************************************
!***************************************************************************

FUNCTION simple_step_control(num_newt_conv, max_Newton_steps, step_ctrl) RESULT(output)
! Function to calculate the increase in time step for the pseudo time-
! integration option based on:
!
!    num_newt_conv,      the number of Newton iterations the last step
!                        required to reach convergence, and
!    max_Newton_steps,   the maximum number of Newton steps allowed.
!    step_ctrl           aggressiveness of step size routine,
!                        0.0 for constant step, 2.0 is very big
!
! This simplistic function will increase the time step if the last step
! converged. It is a quadratic function in the ratio of the number of
! Newton steps taken compared to the maximum number of steps allowed
! up to a maximum of 1+aggressiveness times the previous step.

   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: num_newt_conv, max_Newton_steps
   REAL(KIND=8), INTENT(IN) :: step_ctrl
   REAL(KIND=8)             :: factor, output

   factor = (max_Newton_steps - num_newt_conv) / (max_Newton_steps - 1d0)

   output = 1d0 + step_ctrl * factor * factor

END FUNCTION simple_step_control

!***************************************************************************
!***************************************************************************
!***************************************************************************

SUBROUTINE print_cont_step1(order, step, step_old, con)
! Print out for relevant time step information

   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: order
   REAL(KIND=8), INTENT(IN) :: step, step_old
   TYPE(con_struct)         :: con

   CHARACTER(LEN=70)        :: string


   IF (con%general_info%method == TURNING_POINT_CONTINUATION) THEN

      string = 'Zero-order Turning Point Continuation'

   ELSE IF (con%general_info%method == PHASE_TRANSITION_CONTINUATION) THEN

      string = 'Zero-order Phase Transition Continuation'

   ELSE IF (con%general_info%method == PITCHFORK_CONTINUATION) THEN

      string = 'Zero-order Pitchfork Continuation'

   ELSE IF (con%general_info%method == HOPF_CONTINUATION) THEN

      string = 'Zero-order Hopf Continuation'

   ELSE IF (order == 0)  THEN

      string = 'Zero-order Continuation'

   ELSE IF (order == 1)  THEN

      string = 'First-order Continuation'

   ELSE IF (order == 2)  THEN

      string = 'Pseudo Arc-length Continuation'

   END IF

   WRITE(*,*)
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*) 'Start of Step: ', con%private_info%nstep, ' Continuation Param = ', con%general_info%param, &
              ' from ', con%general_info%param - step
   WRITE(*,*) 'Continuation method = ', string
   WRITE(*,*) '   delta_c_p        = ', step
   WRITE(*,*) '   delta_c_p_old    = ', step_old

   IF (order == 2) &
      WRITE(*,*) '   delta_s = ', con%private_info%arc_step

   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)

END SUBROUTINE print_cont_step1

!***************************************************************************
!***************************************************************************
!***************************************************************************

SUBROUTINE print_cont_step2(order, step, con)
! Print out for relevant continuation step information

   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: order
   REAL(KIND=8), INTENT(IN) :: step
   TYPE(con_struct)         :: con

   WRITE(*,*)
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*) 'Continuation Step Number ', con%private_info%nstep, &
              ' was a success: Param = ', con%general_info%param

   IF (con%general_info%method == TURNING_POINT_CONTINUATION) THEN

      WRITE(*,*) '   Turning Point located at: ', con%general_info%param, con%turning_point_info%bif_param

   ELSE IF (con%general_info%method == PHASE_TRANSITION_CONTINUATION) THEN

      WRITE(*,*) '   Phase Transition located at: ', con%general_info%param, con%phase_transition_info%bif_param

   ELSE IF (con%general_info%method == PITCHFORK_CONTINUATION) THEN

      WRITE(*,*) '   Pitchfork Bifurcation located at: ', con%general_info%param, con%pitchfork_info%bif_param

   ELSE IF (con%general_info%method == HOPF_CONTINUATION) THEN

      WRITE(*,*) '   Hopf Bifurcation located at: ', con%general_info%param, con%hopf_info%bif_param, con%hopf_info%omega

   ELSE

      WRITE(*,*) '   Order   = ', order

   END IF

   IF (order == 2) &
      WRITE(*,*) '   delta_s = ', con%private_info%arc_step

   WRITE(*,*) '   delta_c_p = ', step
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)

END SUBROUTINE print_cont_step2

!***************************************************************************
!***************************************************************************
!***************************************************************************

SUBROUTINE print_cont_step_fail(order, step, con)
! Print Out descriptive information on why the current step failed

   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: order
   REAL(KIND=8), INTENT(IN) :: step
   TYPE(con_struct)         :: con

   WRITE(*,*)
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*) 'Continuation Step Number ', con%private_info%nstep, &
              ' experienced a convergence failure in the nonlinear or linear solver'

   WRITE(*,*) '   Value of parameter at failed step =', con%general_info%param

   IF (order < 2) THEN
      WRITE(*,*) '   delta_c_p of the failed step      = ', 2*step
      WRITE(*,*) '   Next value of delta_c_p           = ', step

   ELSE
      WRITE(*,*) '   delta_s of the failed step        = ', 2.0*con%private_info%arc_step
      WRITE(*,*) '   Next value of delta_s             = ', con%private_info%arc_step
   END IF
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)

END SUBROUTINE print_cont_step_fail

!***************************************************************************
!***************************************************************************
!***************************************************************************

END MODULE Loca_lib
