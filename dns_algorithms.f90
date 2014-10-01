MODULE dns_algorithms
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 2/9/2013
!
   USE dynamic_structures
   USE sparse_matrix_profiles
   USE global_variables
   USE miscellaneous_subroutines
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE start_sparse_kit  ! for collect and extract subroutines
   USE qc_sp_M
   USE qv_sp
   USE par_solve_mumps
   USE restart_io
   USE vtk_plot

   IMPLICIT NONE
   ! variables "global" to this module
   TYPE(CSR_MUMPS_Matrix)                    :: Matr
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: xOld2
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u02

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE dns (x0, x00)
!
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | x0 is an array with (\ub, p) used both as initial    |
!     |    condition and as solution at the previous time    |
!     |    step.                                             |
!     |                                                      |
!     | Global variables used:                               |
!     |                                                      |
!     | u0, p0 solution at previous time step                |
!     |                                                      |
!     | xx, uu, pp solution at current time step             |
!     |                                                      |
!     %------------------------------------------------------%
!
   IMPLICIT NONE

   ! input variables
   REAL(KIND=8), DIMENSION(:) :: x0
   REAL(KIND=8), DIMENSION(:), OPTIONAL :: x00
   ! local variables
   REAL(KIND=8) :: tInit
   REAL(KIND=8) :: tEnd
   REAL(KIND=8) :: dt
   REAL(KIND=8) :: dtPlot
   INTEGER      :: i, nSteps, itPlot
   CHARACTER(LEN=128) :: restart_name

!----------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to dns'
   WRITE(*,*) ''
!----------------------------------------------------

!--------------------------------------
! CHANGE SOME NAMES TO EASE THE READING

   tInit  = p_in % dns_tInit
   tEnd   = p_in % dns_tEnd
   dt     = p_in % dns_dt
   dtPlot = p_in % dns_dtPlot

   WRITE(restart_name, '(A)') trim(p_in%dns_output_directory)//'restartDNS.bin'

!-----------
! PARAMETERS

   nSteps = NINT( (tEnd - tInit) / dt )
   itPlot = NINT( dtPlot / dt )

!----------------
! INITIALIZE Matr

   ALLOCATE( Matr%i      (SIZE(Jacobian%i))       ); Matr%i       = Jacobian%i
   ALLOCATE( Matr%i_mumps(SIZE(Jacobian%i_mumps)) ); Matr%i_mumps = Jacobian%i_mumps
   ALLOCATE( Matr%j      (SIZE(Jacobian%j))       ); Matr%j       = Jacobian%j
   ALLOCATE( Matr%e      (SIZE(Jacobian%e))       ); Matr%e       = 0d0

   SELECT CASE (p_in%dns_method)
      CASE (1)
         CALL par_mumps_master (INITIALIZATION, 7, Matr, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   7, Matr, 0)

      CASE (2)
         CALL par_mumps_master (INITIALIZATION, 7, Matr, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   7, Matr, 0)
         ! assemble matrix
         !
         CALL qc_0y0_zero_sp_M (mm, jj,        1.5d0/dt, Matr) !   Mass
         CALL qc_1y1_sp_gg_M   (mm, jj,        1d0/Re,   Matr) ! + stifness (GRAD:GRAD)
         CALL qc_1y0_sp_M      (mm, jj, jj_L, -1d0,      Matr) ! + pressure gradient (ibp)
         CALL qc_0y1_sp_M      (mm, jj, jj_L, -1d0,      Matr) ! - velocity divergence
         CALL Dirichlet_c_M    (np, js_Axis, js_D,  Matr)
         CALL par_mumps_master (NUMER_FACTOR, 7, Matr, 0)

      CASE DEFAULT
         WRITE(*,*) '*************************************'
         WRITE(*,*) '*** Wrong parameter:              ***'
         WRITE(*,*) '*** p_in % dns_method             ***'
         WRITE(*,*) '*** set to: ', p_in%dns_method
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END SELECT
         
!-------------------------------
! INITIALIZE previous time steps

   ALLOCATE( xOld2(Nx)          )
   ALLOCATE( u02(velCmpnnts,np) )

   IF (PRESENT(x00)) THEN
      ! using the restart
      xOld2 = x00
   ELSE
      ! starting from a steady solution
      xOld2 = x0
   ENDIF

!----------
! TIME LOOP

   WRITE(*,*) '    Start of time loop'

   DO i = 1, nSteps

      WRITE(*,*) ''
      WRITE(*,*) '    time = ', tInit + i*dt

      ! compute solution at new time
      CALL advanceInTime (p_in%dns_method, dt, x0,  xx)

      ! plot and save restart
      IF ( MODULO(i, itPlot) == 0 ) THEN

         CALL plot_frame(xx, tInit + i*dt)

         CALL write_restart_bin(x0, tInit + (i-1)*dt, i, nSteps, trim(restart_name))

      ENDIF

      ! update solution
      xOld2 = x0
      x0    = xx

   ENDDO

   WRITE(*,*) '    End of time loop'


END SUBROUTINE dns

!------------------------------------------------------------------------------

SUBROUTINE advanceInTime (method, dt, xOld, x)

   IMPLICIT NONE
   ! input variables
   INTEGER                    :: method
   REAL(KIND=8)               :: dt
   REAL(KIND=8), DIMENSION(:) :: xOld
   ! output variable
   REAL(KIND=8), DIMENSION(:) :: x
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: uRHS
   REAL(KIND=8), DIMENSION(np_L)          :: pRHS


   SELECT CASE (method)

      CASE (1)
      !------------------------
      ! BDF2 with semi-explicit
      ! convective term
      !------------------------
         ! (1)
         ! extract previous velocity fields
         !
         CALL extract (xOld2,  u02)
         CALL extract (xOld,   u0)
         
         ! (2)
         ! assemble matrix
         !
         Matr%e = 0d0
         CALL qc_0y0_zero_sp_M (mm, jj,        1.5d0/dt, Matr) !   Mass
         CALL qc_advecy_sp_M   (mm, jj,        2*u0-u02, Matr) ! + convective term
         CALL qc_1y1_sp_gg_M   (mm, jj,        1d0/Re,   Matr) ! + stifness (GRAD:GRAD)
         CALL qc_1y0_sp_M      (mm, jj, jj_L, -1d0,      Matr) ! + pressure gradient (ibp)
         CALL qc_0y1_sp_M      (mm, jj, jj_L, -1d0,      Matr) ! - velocity divergence
         
         CALL Dirichlet_c_M (np, js_Axis, js_D,  Matr)
         CALL par_mumps_master (NUMER_FACTOR, 7, Matr, 0)
         
         ! (3)
         ! assemble RHS
         !
         uRHS = 0d0
         CALL qv_0y0_sp (mm, jj, u0,   2/dt,    uRHS)
         CALL qv_0y0_sp (mm, jj, u02, -0.5/dt,  uRHS)
         pRHS = 0d0
         CALL collect (uRHS, pRHS,  xx)
         
         ! (4)
         ! enforce Dirichlet boundary conditions
         !
         CALL Dirichlet_c (np, js_Axis, js_D, bvs_D,  xx)
         IF (DESINGULARIZE) xx(Nx) = 0d0
         
         ! (5)
         ! solve linear system
         CALL par_mumps_master (DIRECT_SOLUTION, 7, Matr, 0, xx)

      CASE (2)
      !------------------------
      ! BDF2 with explicit
      ! convective term
      !------------------------
         ! (1)
         ! extract previous velocity fields
         !
         CALL extract (xOld2,  u02)
         CALL extract (xOld,   u0)
         
         ! (3)
         ! assemble RHS
         !
         uRHS = 0d0
         CALL qv_0y0_sp  (mm, jj, u0,       -2/dt,   uRHS)
         CALL qv_0y0_sp  (mm, jj, u02,       0.5/dt, uRHS)
         CALL qv_0y01_sp (mm, jj, 2*u0-u02,          uRHS)
         pRHS = 0d0
         CALL collect (-uRHS, pRHS,  xx)
         
         ! (4)
         ! enforce Dirichlet boundary conditions
         !
         CALL Dirichlet_c (np, js_Axis, js_D, bvs_D,  xx)
         IF (DESINGULARIZE) xx(Nx) = 0d0
         
         ! (5)
         ! solve linear system
         CALL par_mumps_master (DIRECT_SOLUTION, 7, Matr, 0, xx)

   END SELECT


END SUBROUTINE advanceInTime

!------------------------------------------------------------------------------

SUBROUTINE plot_frame(xx, time)

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:) :: xx
   REAL(KIND=8)               :: time
   ! local variables
   INTEGER           :: time_plot
   CHARACTER(LEN=12) :: counter


   ! generate file name
   time_plot = NINT( time*1d3 )
   CALL intToChar6 (time_plot,  counter)
   ! plot
   CALL extract (xx, uu, pp)
   CALL vtk_plot_P2 (rr, jj, jj_L, uu, pp, trim(p_in%dns_output_directory)//'sol'//trim(counter)//'.vtk')

END SUBROUTINE plot_frame

!==============================================================================

END MODULE dns_algorithms
