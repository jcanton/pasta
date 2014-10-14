! -----------------------------------------------------------------------------
!   LOCA 1.0: Library of Continuation Algorithms

MODULE Loca_interface_i

!******************************************************************************
!******************************************************************************
!******************************************************************************

   USE Loca_types       ! definition of all the types of variables needed
   USE Loca_pd          ! declaration of the passdown structure 'pd'
   USE Loca_lib         ! con_lib
   USE Loca_interface_o ! cvarsparser, print_con_struct, print_final

   IMPLICIT NONE


CONTAINS
!=======
!
!  SUBROUTINE do_loca(pd)
!  FUNCTION continuation_hook(x, delta_x, con, reltol, abstol) RESULT (output)


!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE do_loca(pd)
   !
   !*************************************
   ! edit: removed con_struct from input,
   ! use the passdown_struct in its
   ! place and read con in cvarsparser
   !*************************************

   IMPLICIT NONE
   ! input variable
   TYPE(passdown_struct) :: pd
   ! local variables
   TYPE(con_struct) :: con
   INTEGER                   :: nstep

   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) ' --> CALL to do_loca'
   WRITE(*,*) ''

   ! read input file
   !
   CALL cvarsparser(con, pd, "loca_data.in")

   ! Assign problem information coming through argument list to elements
   ! of the passdown structure. This makes them global to the file and
   ! can be accessed below in the conwrap routines.
   !
   !pd%matrix_save    = FALSE

   ! Fill all the necessary information into the 'con' structure
   ! for running the continuation routines.

   ! First load the general info structure
   !
!   con%general_info%method          = con_input%general_info%method
!   con%general_info%param           = con_input%general_info%param
!   con%general_info%x               = con_input%general_info%x
!   con%general_info%numUnks         = con_input%general_info% numUnks
!   con%general_info%numOwnedUnks    = con_input%general_info%numOwnedUnks
!   con%general_info%printproc       = con_input%general_info%printproc
   con%general_info%nv_restart      = FALSE
   con%general_info%nv_save         = FALSE
   con%general_info%perturb         = 1.0d-6

   ! Then load the stepping info structure
   !
!   con%stepping_info%first_step     = con_input%stepping_info%first_step
   con%stepping_info%base_step      = 0
!   con%stepping_info%max_steps      = con_input%stepping_info%max_steps
!   con%stepping_info%max_param      = con_input%stepping_info%max_param
!   con%stepping_info%max_delta_p    = con_input%stepping_info%max_delta_p
!   con%stepping_info%step_ctrl      = con_input%stepping_info%step_ctrl
!   con%stepping_info%max_newton_its = con_input%stepping_info%max_newton_its

   ! Then load one of the method dependent structures
   !
   SELECT CASE (con%general_info%method)

      CASE(ZERO_ORDER_CONTINUATION)

      CASE (FIRST_ORDER_CONTINUATION)

      CASE (ARC_LENGTH_CONTINUATION)
!         pd%matrix_save = FALSE
!         con%arclength_info%dp_ds2_goal     = con_input%arclength_info%dp_ds2_goal
!         con%arclength_info%dp_ds_max       = con_input%arclength_info%dp_ds_max
!         con%arclength_info%tang_exp        = con_input%arclength_info%tang_exp
!         con%arclength_info%tang_step_limit = con_input%arclength_info%tang_step_limit

      CASE (TURNING_POINT_CONTINUATION)
!         con%turning_point_info%bif_param = con_input%turning_point_info%bif_param
         NULLIFY(con%turning_point_info%nv)

      CASE (PITCHFORK_CONTINUATION)
!         con%pitchfork_info%bif_param  = con_input%pitchfork_info%bif_param
!         con%pitchfork_info%psi        = con_input%pitchfork_info%psi

      CASE (HOPF_CONTINUATION)
!         pd%matrix_save = FALSE
!         con%hopf_info%bif_param   = con_input%hopf_info%bif_param
!         con%hopf_info%omega       = con_input%hopf_info%omega
!         con%hopf_info%y_vec       = con_input%hopf_info%y_vec
!         con%hopf_info%z_vec       = con_input%hopf_info%z_vec
!         con%hopf_info%mass_flag   = con_input%hopf_info%mass_flag

      CASE (HOPF_BETA_CONTINUATION)
!         pd%matrix_save = FALSE
!         con%hopf_info%bif_param   = con_input%hopf_info%bif_param
!         con%hopf_info%omega       = con_input%hopf_info%omega
!         con%hopf_info%y_vec       = con_input%hopf_info%y_vec
!         con%hopf_info%z_vec       = con_input%hopf_info%z_vec
!         con%hopf_info%mass_flag   = con_input%hopf_info%mass_flag

      CASE (MANIFOLD_CONTINUATION)
         WRITE(*,*) 'manifold continuation not yet implemented'
         con%manifold_info%k =  2                      ! number of parameters
         NULLIFY(con%manifold_info%param_vec)          ! vec of length k
         NULLIFY(con%manifold_info%param_lower_bounds) ! vec of length k
         NULLIFY(con%manifold_info%param_upper_bounds) ! vec of length k
         NULLIFY(con%manifold_info%param_steps)        ! vec of length k

      CASE (PHASE_TRANSITION_CONTINUATION)
         WRITE(*,*) 'phase transition continuation not yet implemented'
!         con%phase_transition_info%bif_param = con_input%phase_transition_info%bif_param
!         con%phase_transition_info%x2        = con_input%phase_transition_info%x2

      CASE DEFAULT
         WRITE(*,*) 'ERROR: Unknown LOCA input method: ', con%general_info%method
         STOP ! FixMe change with MPI_ABORT?

   END SELECT

   ! Finally, load the eigensolver structures
   !
!   con%eigen_info%Num_Eigenvalues   = con_input%eigen_info%Num_Eigenvalues
   ! This is a temporary default:
   con%eigen_info%Num_Eigenvectors  = con%eigen_info%Num_Eigenvalues
   con%eigen_info%sort              = TRUE
!   con%eigen_info%Shift_Point(0)    = con_input%eigen_info%Shift_Point(0)
!   con%eigen_info%Shift_Point(1)    = con_input%eigen_info%Shift_Point(1)
!   con%eigen_info%Shift_Point(2)    = con_input%eigen_info%Shift_Point(2)
!   con%eigen_info%Arnoldi           = con_input%eigen_info%Arnoldi
!   con%eigen_info%Residual_Tol(0)   = con_input%eigen_info%Residual_Tol(0)
!   con%eigen_info%Residual_Tol(1)   = con_input%eigen_info%Residual_Tol(1)
!   con%eigen_info%Max_Iter          = con_input%eigen_info%Max_Iter
!   con%eigen_info%Every_n_Steps     = con_input%eigen_info%Every_n_Steps

   ! print out continuation structure
   !
   CALL print_con_struct(con)


   ! CALL LOCA, now that all information is set
   !
   nstep = con_lib(con)


   ! Final printing and cleanup here
   !
   IF (con%general_info%printproc == TRUE) THEN
      CALL print_final(con%general_info%param, nstep)
!      CALL print_final(con%general_info%param, nstep, pd%num_mat_fills, &
!                       pd%num_res_fills, pd%num_linear_its)
   END IF

   SELECT CASE(con%general_info%method)
      CASE(PITCHFORK_CONTINUATION)
         DEALLOCATE(con%pitchfork_info%psi)

      CASE(HOPF_CONTINUATION)
         DEALLOCATE(con%hopf_info%y_vec)
         DEALLOCATE(con%hopf_info%z_vec)

      CASE(HOPF_BETA_CONTINUATION)
         DEALLOCATE(con%hopf_info%y_vec)
         DEALLOCATE(con%hopf_info%z_vec)

   END SELECT

!   DEALLOCATE(con)

END SUBROUTINE do_loca

!******************************************************************************
!******************************************************************************
!******************************************************************************

END MODULE Loca_interface_i
