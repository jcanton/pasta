! -----------------------------------------------------------------------------
!   LOCA 1.0: Library of Continuation Algorithms

MODULE Loca_interface_o

!******************************************************************************
!******************************************************************************
!******************************************************************************

   USE loca_types    ! definition of all the types of variables needed
   USE loca_pd       ! declaration of the passdown structure 'pd'
   USE eigensolve    ! read_eigenvector
   USE restart_io    ! write_restart, write_QP_restart
   USE loca_wrappers ! param_out, vtk_plot_loca, compute_eigen

   IMPLICIT NONE


CONTAINS
!=======
!
!  SUBROUTINE cvarsparser(con, pd, filename)
!  SUBROUTINE print_con_struct(con)
!  SUBROUTINE print_final(param, step_num)
!  SUBROUTINE solution_output_conwrap(num_soln_flag, x,  param, x2, param2, x3, param3, step_num, num_its, con)


!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE cvarsparser(con, pd, filename)
! Input file reader code
! The top sections are problem specific information
! The bottom sections can be used by any LOCA application
! since they read directly into the LOCA con structure

   IMPLICIT NONE
   ! input variables
   TYPE(passdown_struct) :: pd
   CHARACTER(*)          :: filename
   ! output variable
   !TYPE(con_struct) :: con
   TYPE(con_struct), POINTER :: con
   ! local variables
   INTEGER           :: i
   INTEGER           :: hopf_nev,  hopf_ind
   INTEGER           :: pitch_nev, pitch_ind
   REAL(KIND=8)      :: bif_param_init = 0
   CHARACTER(LEN=90) :: tmp, temp1, temp2
   CHARACTER(LEN=90) :: temp ! FixMe this string may need to be trimmed to fit exactly the input file format

   OPEN(UNIT=11, FILE=filename,  FORM='formatted', STATUS='unknown')

   !... Problem specific Info ...
   !*************************************
   ! <!--Passdown-->
   READ(11,*) tmp
   !*************************************
!   READ(11,*) temp, pd%lambda
!   READ(11,*) temp, pd%gamma
!   READ(11,*) temp, pd%mu
!   READ(11,*) temp, pd%alpha

!   READ(11,*) temp, pd%problem
!   READ(11,*) temp, pd%n
!   READ(11,*) temp, pd%maxiter
!   READ(11,*) temp, pd%tol
   READ(11,*) temp, pd%debug

   ! For serial runs, this processor is always the printing processsor
   ! con%general_info%printproc = 1;
   !
   ! temporary fix for turning off printing for debug=0
   ! if (pd%debug==0) con%general_info%printproc=0;
   ! Doing it this way will allow the print level to be set (serial only)
   ! FixMe in the old version here I had written con%general_info%printproc = 8, but I don't know why...
   con%general_info%printproc = pd%debug

   !pd%num_res_fills  = 0
   !pd%num_mat_fills  = 0
   pd%num_linear_its = 0

   !... General Info ...*/
   con%general_info%numUnks      = pd%ldz
   con%general_info%numOwnedUnks = pd%ldz
   con%general_info%x            => pd%x

#if DEBUG > 2
   ! FixMe need to check if the indices of this vector go from 0 to ldz-1 or from 1 to ldz
   WRITE(*,*) 'Check to have correctly received the initial solution'
   WRITE(*,*) '   pd->x[0]   = ', con%general_info%x(0)
   WRITE(*,*) '   pd->x[end] = ', con%general_info%x(pd%ldz-1)
#endif

   !... General Info ...
   !*************************************
   ! <!--General_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, con%general_info%method
   READ(11,*)  temp, pd%bif_param

   ! Assign initial value of bifurcation param based on flag value
   !
   WRITE(*,*) 'Bifurcation parameter:'
   !
   SELECT CASE (pd%bif_param)

      CASE (REYNOLDS)
         bif_param_init = pd%reynolds
         WRITE(*,*) '    Reynolds = ', bif_param_init

      CASE (OSCAR)
         bif_param_init = pd%oscar
         WRITE(*,*) '    oscar = ', bif_param_init

      CASE (ROMEO)
         bif_param_init = pd%romeo
         WRITE(*,*) '    romeo = ', bif_param_init

      CASE (WHISKY)
         bif_param_init = pd%whisky
         WRITE(*,*) '    whisky = ', bif_param_init

   END SELECT

   READ(11,*) temp, pd%param

   WRITE(*,*) 'Continuation parameter:'
   !
   SELECT CASE (pd%param)

      CASE (REYNOLDS)
         con%general_info%param = pd%reynolds
         WRITE(*,*) '    Reynolds = ', con%general_info%param

      CASE (OSCAR)
         con%general_info%param = pd%oscar
         WRITE(*,*) '    oscar = ', con%general_info%param

      CASE (ROMEO)
         con%general_info%param = pd%romeo
         WRITE(*,*) '    romeo = ', con%general_info%param

      CASE (WHISKY)
         con%general_info%param = pd%whisky
         WRITE(*,*) '    whisky = ', con%general_info%param

   END SELECT

   WRITE(*,*) 'Method:'
   !
   SELECT CASE (con%general_info%method)

      CASE (ZERO_ORDER_CONTINUATION)
         WRITE(*,*) '    Zero-order continuation'

      CASE (FIRST_ORDER_CONTINUATION)
         WRITE(*,*) '    First-order continuation'

      CASE (ARC_LENGTH_CONTINUATION)
         WRITE(*,*) '    Arc length continuation'

      CASE (TURNING_POINT_CONTINUATION)
         WRITE(*,*) '    Turning point continuation'
         con%turning_point_info%bif_param = bif_param_init

      CASE (PITCHFORK_CONTINUATION)
         WRITE(*,*) '    Pitckfork continuation'
         con%pitchfork_info%bif_param = bif_param_init

      CASE (HOPF_CONTINUATION)
         WRITE(*,*) '    Hopf continuation'
         con%hopf_info%bif_param = bif_param_init

      CASE (PHASE_TRANSITION_CONTINUATION)
         WRITE(*,*) '    Phase transition continuation'
         WRITE(*,*) '*** WARNING ***'

      CASE (AUGMENTING_CONDITION)
         WRITE(*,*) '    Augmenting condition'
         WRITE(*,*) '*** WARNING ***'

      CASE (MANIFOLD_CONTINUATION)
         WRITE(*,*) '    Manifold continuation'
         WRITE(*,*) '*** WARNING ***'

      CASE DEFAULT
         WRITE(*,*) 'Unknown method in input'
         WRITE(*,*) 'STOP'
         STOP ! FixMe change with MPI_ABORT?

   END SELECT

   !... Stepping Info ...
   !*************************************
   ! <!--Stepping_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, con%stepping_info%max_steps
   READ(11,*)  temp, con%stepping_info%max_param
   READ(11,*)  temp, con%stepping_info%first_step
   READ(11,*)  temp, con%stepping_info%step_ctrl
   READ(11,*)  temp, con%stepping_info%max_delta_p

   WRITE(*,*) 'Continuation parameter goal value = ', con%stepping_info%max_param

   con%stepping_info%max_newton_its = pd%maxiter;

   !... Arclength Info ...
   !*************************************
   ! <!--Arclength_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, con%arclength_info%dp_ds2_goal
   READ(11,*)  temp, con%arclength_info%dp_ds_max
   READ(11,*)  temp, con%arclength_info%tang_exp
   READ(11,*)  temp, con%arclength_info%tang_step_limit

   !... Pitchfork Info ...
   !*************************************
   ! <!--Pitchfork_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, pitch_nev
   READ(11,*)  temp, pitch_ind
   READ(11,*)  temp, temp1 ! load directory FixMe change these in the 'loca_data.in' files
   READ(11,*)  temp, temp2 ! file name      FixMe change these in the 'loca_data.in' files
   !
   IF (con%general_info%method == PITCHFORK_CONTINUATION) THEN
      IF(ADJUSTL(temp1) /= "none") THEN
         ALLOCATE(con%pitchfork_info%psi(0:pd%ldz-1))
         IF (pitch_nev == 1) THEN
            ! FixMe need to understand what needs to be loaded and use an appropriate routine
            !CALL load(con%pitchfork_info%psi, pd%ldz,1, trim(temp1)//trim(temp2))
            WRITE(*,*) 'loaded psi vector (for pitchfork) from filename: ', trim(temp1)//trim(temp2)
         ELSE
            ! FixMe need to understand what's needed here
            ! in the old version I had a commented call to 'pitch_input_arpack'
            ! for pitch_nev > 1, I would not know what to do for pitch_nev = 0...
            !WRITE(*,*) 'loaded psi vector (for pitchfork) from arpack sol. in filename: ', trim(temp1)//trim(temp2)
         ENDIF
      ELSE
         ! FixMe why would I end up here?
         WRITE(*,*) '*** Not loading psi vector?? ***'
      END IF
   END IF

   !... Hopf Info ...
   !*************************************
   ! <!--Hopf_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, hopf_nev
   READ(11,*)  temp, temp1

   IF (con%general_info%method == HOPF_CONTINUATION) THEN
      IF (ADJUSTL(temp1) /= "none") THEN
         ALLOCATE(con%hopf_info%y_vec(0:pd%ldz-1))
         ALLOCATE(con%hopf_info%z_vec(0:pd%ldz-1))
         IF (hopf_nev == 1) THEN
            CALL read_eigenvector(pd%ldz, hopf_nev, trim(temp1), con%hopf_info%y_vec, con%hopf_info%z_vec)
            WRITE(*,*) 'loaded y and z vectors (for hopf) from filename: ', trim(temp1)
#if DEBUG > 2
            WRITE(*,*) 'Check to have correctly received the eigenvector for hopf'
            WRITE(*,*) '    y_vec[0] = ', con%hopf_info%y_vec(0), ' y_vec[end] = ', con%hopf_info%y_vec(pd%ldz-1)
            WRITE(*,*) '    z_vec[0] = ', con%hopf_info%z_vec(0), ' z_vec[end] = ', con%hopf_info%z_vec(pd%ldz-1)
#endif
         ELSE
            ! FixMe need to understand what's needed here
         ENDIF
      ELSE
         ! FixMe why would I end up here?
         WRITE(*,*) '*** Not loading hopf eigenvector?? ***'
      END IF
   END IF

   READ(11,*)  temp, con%hopf_info%omega
   READ(11,*)  temp, con%hopf_info%mass_flag

   !... Eigen Info ...
   !*************************************
   ! <!--Eigen_Info-->
   READ(11,*)  tmp
   !*************************************
   READ(11,*)  temp, con%eigen_info%Num_Eigenvalues ! used as on/off flag for eigenvalues computation
   READ(11,*)  temp, con%eigen_info%Every_n_Steps


   CLOSE(11)

END SUBROUTINE cvarsparser

!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE print_con_struct(con)
! Routine for printing out the con structure to the screen

   IMPLICIT NONE

   TYPE(con_struct), POINTER :: con

   WRITE(*,*)
   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*) 'Continuation structure:'
   WRITE(*,*)

   WRITE(*,*) '    con%general_info%param =             ', con%general_info%param
   WRITE(*,*) '    con%general_info%numUnks =           ', con%general_info%numUnks
   WRITE(*,*) '    con%general_info%numOwnedUnks =      ', con%general_info%numOwnedUnks
   WRITE(*,*) '    con%general_info%printproc =         ', con%general_info%printproc

   WRITE(*,*) '    con%stepping_info%first_step =       ', con%stepping_info%first_step
   WRITE(*,*) '    con%stepping_info%max_steps =        ', con%stepping_info%max_steps
   WRITE(*,*) '    con%stepping_info%max_param =        ', con%stepping_info%max_param
   WRITE(*,*) '    con%stepping_info%max_delta_p =      ', con%stepping_info%max_delta_p
   WRITE(*,*) '    con%stepping_info%step_ctrl =        ', con%stepping_info%step_ctrl
   WRITE(*,*) '    con%stepping_info%max_newton_its =   ', con%stepping_info%max_newton_its

   IF(con%general_info%method==ARC_LENGTH_CONTINUATION) THEN
      WRITE(*,*) '    con%arclength_info%dp_ds2_goal =     ', con%arclength_info%dp_ds2_goal
      WRITE(*,*) '    con%arclength_info%dp_ds_max =       ', con%arclength_info%dp_ds_max
      WRITE(*,*) '    con%arclength_info%tang_exp =        ', con%arclength_info%tang_exp
      WRITE(*,*) '    con%arclength_info%tang_step_limit = ', con%arclength_info%tang_step_limit
   END IF

   IF(con%general_info%method==TURNING_POINT_CONTINUATION) &
      WRITE(*,*) '    con%turning_point_info%bif_param =   ', con%turning_point_info%bif_param

   IF(con%general_info%method==PITCHFORK_CONTINUATION) &
      WRITE(*,*) '    con%pitchfork_info%bif_param =       ', con%pitchfork_info%bif_param

   IF(con%general_info%method==HOPF_CONTINUATION) &
      WRITE(*,*) '    con%hopf_info%bif_param =            ', con%hopf_info%bif_param

   IF(con%eigen_info%Num_Eigenvalues > 0) THEN
!      WRITE(*,*) 'con%eigen_info%Shift_Point(0) =      ', con%eigen_info%Shift_Point(0)
!      WRITE(*,*) 'con%eigen_info%Shift_Point(1) =      ', con%eigen_info%Shift_Point(1)
!      WRITE(*,*) 'con%eigen_info%Arnoldi =             ', con%eigen_info%Arnoldi
!      WRITE(*,*) 'con%eigen_info%Residual_Tol(0) =     ', con%eigen_info%Residual_Tol(0)
      WRITE(*,*) '    con%eigen_info%Every_n_Steps =       ', con%eigen_info%Every_n_Steps
   END IF

   WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)

END SUBROUTINE print_con_struct

!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE print_final(param, step_num) !, num_mat_fills, num_res_fills, total_linear_its)
! Print out the final results and counters

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: param
   INTEGER,      INTENT(IN) :: step_num
!   INTEGER,      INTENT(IN) :: num_mat_fills
!   INTEGER,      INTENT(IN) :: num_res_fills
!   INTEGER,      INTENT(IN) :: total_linear_its

   WRITE(*,*)
   WRITE(*,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*) 'CONTINUATION ROUTINE HAS FINISHED: '
   WRITE(*,*) 'Ending Parameter value     = ', param
   WRITE(*,*) 'Number of steps            = ', step_num
!   WRITE(*,*) 'Number of Matrix fills     = ', num_mat_fills
!   WRITE(*,*) 'Number of Residual fills   = ', num_res_fills
!   WRITE(*,*) 'Number of linear solve its = ', total_linear_its
   WRITE(*,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(*,*)

END SUBROUTINE print_final

!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE solution_output_conwrap(num_soln_flag, x,  param,  &
                                                  x2, param2, &
                                                  x3, param3, &
                                                  step_num, num_its, con)

! Put the call to your solution output (both file and screen) routines here.
! Input:
!    num_soln_flag  Flag for how many solution vectors are being passed for
!                   output. For parameter continuation and turning point
!                   tracking there is just 1, for pitchfork and phase
!                   transitions there are 2, and for Hopfs there are 3.
!    x            First solution vector for output.
!    param        Continuation parameter value.
!    x2           Second solution vector for output (y_vec or x2)
!    param2       Bifurcation parameter value.
!    x3           Third solution vector for output (z_vec for Hopfs)
!    param3       Third Parameter value (frequency Hopfs)
!    step_num+1   Time index to output to (step_num is 0 based).
!    num_its      Number of Newton iterations used for for convergence
!    con          pointer to continuation structure, for passing to
!                 the eigensolver.

   IMPLICIT NONE

   INTEGER                     :: num_soln_flag
   REAL(KIND=8), DIMENSION(0:) :: x
   REAL(KIND=8)                :: param
   REAL(KIND=8), DIMENSION(0:) :: x2
   REAL(KIND=8)                :: param2
   REAL(KIND=8), DIMENSION(0:) :: x3
   REAL(KIND=8)                :: param3
   INTEGER                     :: step_num
   INTEGER                     :: num_its
   TYPE(con_struct), POINTER   :: con

   CHARACTER(LEN=64) :: filenm
   CHARACTER(LEN=6)  :: str1, str2
   CHARACTER(LEN=14) :: str3
   REAL(KIND=8)      :: shiftIm=0d0


   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) ' --> CALL to solution_output_conwrap'
   WRITE(*,*) ''

   ! prepare part of the file name
   !
   IF ( step_num < 10 ) THEN
      WRITE(str1, '(i1)') step_num
   ELSEIF (step_num < 100 ) THEN
      WRITE(str1, '(i2)') step_num
   ELSEIF (step_num < 1000 ) THEN
      WRITE(str1, '(i3)') step_num
   ELSEIF (step_num < 10000 ) THEN
      WRITE(str1, '(i4)') step_num
   ELSEIF (step_num < 100000 ) THEN
      WRITE(str1, '(i5)') step_num
   ELSEIF (step_num < 1000000 ) THEN
      WRITE(str1, '(i6)') step_num
   ELSE
      WRITE(*,*) '*************************************************'
      WRITE(*,*) '*** Error:                                    ***'
      WRITE(*,*) '*** Integer too large.                        ***'
      WRITE(*,*) '*** Change subroutine solution_output_conwrap ***'
      WRITE(*,*) '*** in module loca_interface                  ***'
      WRITE(*,*) '*************************************************'
      WRITE(*,*) 'STOP.'
      STOP ! FixMe
      !CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   IF ( con%stepping_info%max_steps < 10 ) THEN
      WRITE(str2, '(i1)') con%stepping_info%max_steps
   ELSEIF (con%stepping_info%max_steps < 100 ) THEN
      WRITE(str2, '(i2)') con%stepping_info%max_steps
   ELSEIF (con%stepping_info%max_steps < 1000 ) THEN
      WRITE(str2, '(i3)') con%stepping_info%max_steps
   ELSEIF (con%stepping_info%max_steps < 10000 ) THEN
      WRITE(str2, '(i4)') con%stepping_info%max_steps
   ELSEIF (con%stepping_info%max_steps < 100000 ) THEN
      WRITE(str2, '(i5)') con%stepping_info%max_steps
   ELSEIF (con%stepping_info%max_steps < 1000000 ) THEN
      WRITE(str2, '(i6)') con%stepping_info%max_steps
   ELSE
      WRITE(*,*) '*************************************************'
      WRITE(*,*) '*** Error:                                    ***'
      WRITE(*,*) '*** Integer too large.                        ***'
      WRITE(*,*) '*** Change subroutine solution_output_conwrap ***'
      WRITE(*,*) '*** in module loca_interface                  ***'
      WRITE(*,*) '*************************************************'
      WRITE(*,*) 'STOP.'
      STOP ! FixMe
      !CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   WRITE (str3,*) '_'//trim(str1)//'_'//trim(str2)

   ! Any Continuation run
   !
   IF(num_soln_flag==1 .OR. num_soln_flag==2 .OR. num_soln_flag==3) THEN

      WRITE (filenm,*) 'locaCont'//trim(str3)//'.dat'
      CALL write_restart(x, param, step_num, con%stepping_info%max_steps, trim(filenm))

      WRITE (filenm,*) 'suite'//trim(str3)//'.QPrestart'
      CALL write_QP_restart(x, trim(filenm))

      ! Print out parameter to file
      !
      CALL param_output(param)

   END IF

   ! 2-parameter bifurcation tracking runs
   !
   IF(num_soln_flag==2 .OR. num_soln_flag==3) THEN

      WRITE (filenm,*) 'locaBifTrack'//trim(str3)//'.dat'
      CALL write_restart(x2, param2, step_num, con%stepping_info%max_steps, trim(filenm))

   END IF

   ! For Hopf tracking, a third solution and parameter to write out
   !
   IF(num_soln_flag==3) THEN

      WRITE (filenm,*) 'locaHopf'//trim(str3)//'.dat'
      CALL write_restart(x3, param3, step_num, con%stepping_info%max_steps, trim(filenm))

      shiftIm = param3

   END IF

   WRITE(*,*)

   ! plot the solution
   !
   CALL vtk_plot_loca(x, trim(str3))

   ! compute eigenvalues
   !
   IF (con%eigen_info%Num_Eigenvalues > 0) THEN
      IF (step_num == 0 .OR. MODULO(step_num,con%eigen_info%Every_n_Steps) == 0) THEN
         CALL compute_eigen(x, trim(str3), shiftIm)
      END IF
   END IF

END SUBROUTINE solution_output_conwrap

!******************************************************************************
!******************************************************************************
!******************************************************************************

END MODULE Loca_interface_o
