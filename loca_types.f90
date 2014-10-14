MODULE Loca_types

   IMPLICIT NONE

!****************************************************************************
!                       PARAMETER STATEMENTS                                *
!****************************************************************************
   INTEGER, PARAMETER :: TRUE = 1
   INTEGER, PARAMETER :: FALSE = 0

   ! Choices for continuation method
   INTEGER, PARAMETER :: ZERO_ORDER_CONTINUATION       = 0
   INTEGER, PARAMETER :: FIRST_ORDER_CONTINUATION      = 1
   INTEGER, PARAMETER :: ARC_LENGTH_CONTINUATION       = 2
   INTEGER, PARAMETER :: TURNING_POINT_CONTINUATION    = 3
   INTEGER, PARAMETER :: PITCHFORK_CONTINUATION        = 4
   INTEGER, PARAMETER :: HOPF_CONTINUATION             = 5
   INTEGER, PARAMETER :: HOPF_BETA_CONTINUATION        = 51
   INTEGER, PARAMETER :: PHASE_TRANSITION_CONTINUATION = 6
   INTEGER, PARAMETER :: AUGMENTING_CONDITION          = 7
   INTEGER, PARAMETER :: MANIFOLD_CONTINUATION         = 8
   INTEGER, PARAMETER :: LOCA_LSA_ONLY                 = 9

   ! Choices for matrix fill -- what quantities to calculate
   INTEGER, PARAMETER ::  RHS_ONLY        = 100
   INTEGER, PARAMETER ::  MATRIX_ONLY     = 101
   INTEGER, PARAMETER ::  RHS_MATRIX      = 102
   INTEGER, PARAMETER ::  RHS_MATRIX_SAVE = 103
   INTEGER, PARAMETER ::  RECOVER_MATRIX  = 104
   ! FixMe add new parameters for complex Hopf

   ! Choices for linear solve about the state of the Jacobian matrix
   INTEGER, PARAMETER ::  NEW_JACOBIAN               = 200
   INTEGER, PARAMETER ::  OLD_JACOBIAN               = 201
   INTEGER, PARAMETER ::  SAME_BUT_UNSCALED_JACOBIAN = 202
   INTEGER, PARAMETER ::  OLD_JACOBIAN_DESTROY       = 203
   INTEGER, PARAMETER ::  CHECK_JACOBIAN             = 204
   ! FixMe add new parameters for complex Hopf

   ! Internal flags for rhs calculation for continuation linear solves
   INTEGER, PARAMETER ::  CONT_TANGENT     = 300
   INTEGER, PARAMETER ::  ARC_CONT_SOL2    = 301
   INTEGER, PARAMETER ::  TP_CONT_SOL2     = 302
   INTEGER, PARAMETER ::  TP_CONT_SOL3     = 303
   INTEGER, PARAMETER ::  TP_CONT_SOL4     = 304
   INTEGER, PARAMETER ::  HP_CONT_SOL3     = 305
   INTEGER, PARAMETER ::  HP_CONT_DMDX     = 306
   INTEGER, PARAMETER ::  HP_CONT_DMDPARAM = 307
   ! FixMe add new parameters for complex Hopf

   INTEGER, PARAMETER :: REYNOLDS = 0
   INTEGER, PARAMETER :: OSCAR    = 1
   INTEGER, PARAMETER :: ROMEO    = 2
   INTEGER, PARAMETER :: WHISKY   = 3

!***************************************************************************
!                       STRUCTURE DEFINITIONS
!***************************************************************************

!
! The first structures are sub-structures of con_struct, the
! structure containing all the user-supplied info.


   TYPE general_info_struct
      INTEGER :: method                        ! Solution strategy: see method choices above
      REAL(KIND=8) :: param                    ! value of continuation parameter
      REAL(KIND=8), DIMENSION(:), POINTER :: x ! current solution vector FixMe remove pointer?
      REAL(KIND=8) :: perturb                  ! Perturbation magnitude
      INTEGER :: numUnks                       ! Num unknowns on this Proc (including externals)
      INTEGER :: numOwnedUnks                  ! Num unknowns on this Proc (NOT incl. externals)
      INTEGER :: printproc                     ! Logical indicating if this Proc prints to stdout
      INTEGER :: nv_restart                    ! Restarted null vector flag
      INTEGER :: nv_save                       ! Null vector save flag
   END TYPE general_info_struct

   TYPE stepping_info_struct
      REAL(KIND=8) :: first_step  ! Initial step size
      INTEGER :: base_step        ! Number of the first step (for output)
      INTEGER :: max_steps        ! Maximum # of continuation steps
      INTEGER :: last_step        ! Last step flag
      REAL(KIND=8) :: max_param   ! parameter value to stop at
      REAL(KIND=8) :: max_delta_p ! continuation parameter step limit
      REAL(KIND=8) :: min_delta_p ! continuation parameter step limit
      REAL(KIND=8) :: step_ctrl   ! step aggressiveness -- 0.0 for constant step
      INTEGER :: max_newton_its   ! Max # Newton steps, used only for step control
   END TYPE stepping_info_struct

   TYPE arclength_info_struct
      REAL(KIND=8) :: dp_ds2_goal     ! Square of target dp_ds value for rescaling (desired solution contribution to arc length)
      REAL(KIND=8) :: dp_ds_max       ! High dp_ds value at which to rescale
      REAL(KIND=8) :: tang_exp        ! Power to which tang_factor is raised
      REAL(KIND=8) :: tang_step_limit ! Minimum value of tang_factor between steps
   END TYPE arclength_info_struct

   TYPE turning_point_info_struct
      REAL(KIND=8) :: bif_param                 ! Initial guess of bifurcation parameter
      REAL(KIND=8), DIMENSION(:), POINTER :: nv ! Restarted null vector (read in from file)
   END TYPE turning_point_info_struct

   TYPE pitchfork_info_struct
      REAL(KIND=8) :: bif_param                  ! Initial guess of bifurcation parameter
      REAL(KIND=8), DIMENSION(:), POINTER :: psi ! Antisymmetry vector (also init guess for the null vector, y_vec)
   END TYPE pitchfork_info_struct

   TYPE hopf_info_struct
      REAL(KIND=8) :: bif_param                    ! Initial guess of bifurcation parameter
      REAL(KIND=8) :: omega                        ! Initial guess of Hopf frequency
      REAL(KIND=8), DIMENSION(:), POINTER :: y_vec ! Initial guess of null vector (real)
      REAL(KIND=8), DIMENSION(:), POINTER :: z_vec ! Initial guess of null vector (imag)
      INTEGER :: mass_flag
   END TYPE hopf_info_struct

   TYPE phase_transition_info_struct
      REAL(KIND=8) :: bif_param                 ! Initial guess of bifurcation parameter
      REAL(KIND=8), DIMENSION(:), POINTER :: x2 ! Initial guess of second_solution_vector
   END TYPE phase_transition_info_struct

   TYPE manifold_info_struct
      INTEGER ::  k                                             ! Manifold Dimension
      REAL(KIND=8), DIMENSION(:), POINTER :: param_vec          ! Parameter Vector
      REAL(KIND=8), DIMENSION(:), POINTER :: param_lower_bounds ! Parameter lower bounds
      REAL(KIND=8), DIMENSION(:), POINTER :: param_upper_bounds ! Parameter upper bounds
      REAL(KIND=8), DIMENSION(:), POINTER :: param_steps        ! Parameter step guess
      ! FixMe commented because the link to the MF library is missing
      !MFNKMatrix phi            ! Matrix of null vectors
      !MFNVector u0              ! Previous solution
      !MFNVector u               ! Current solution
      !MFNSpace space            ! Pointer to MF's LOCA NSpace
   END TYPE manifold_info_struct

   TYPE eigen_info_struct
      INTEGER :: Num_Eigenvalues                   ! Number of Eigenvalues to Calculate
      INTEGER :: Num_Eigenvectors                  ! Number of Eigenvectors to Write
      INTEGER :: sort                              ! Flag to sort eigenvalues by real part
      REAL(KIND=8), DIMENSION(0:2) :: Shift_Point  ! Point for shift and invert (Real, Imaginary) and for cayley delta
      INTEGER :: Arnoldi                           ! Arnoldi Space Size
      REAL(KIND=8), DIMENSION(0:1) :: Residual_Tol ! Convergence Tolerance for the Eigenvalue Residual Equation, and linear solver tol
      INTEGER :: Max_Iter                          ! Maximum number of iterations of eigensolver
      INTEGER :: Every_n_Steps                     ! Allow for eigenvalue calc every n steps along a continuation run
   END TYPE eigen_info_struct

   TYPE private_info_struct
      INTEGER :: mass_x                                ! flag that turns on dM/dx in komplex solves
      INTEGER :: mass_param                            ! flag that turns on dM/d(param) in komplex solves
      INTEGER :: first_iter                            ! flag for first Newton iter of each solve
      INTEGER :: step_num                              ! Current continuation step number
      INTEGER :: nstep                                 ! Current step number (for output)
      REAL(KIND=8) :: param_old                        ! old value of continuation parameter
      REAL(KIND=8) :: arc_step                         ! step size of arc length variable
      REAL(KIND=8) :: dp_ds                            ! derivative of parameter w.r.t. arc length
      REAL(KIND=8), DIMENSION(:), POINTER :: x_old     ! previous solution vector
      REAL(KIND=8), DIMENSION(:), POINTER :: x_tang    ! tangent to the solution vector w.r.t. parameter
      REAL(KIND=8), DIMENSION(:), POINTER :: scale_vec ! scaling vector for better arclength conditioning
   END TYPE private_info_struct


!
!  con_struct: structure of parameters for controlling the continuation routines
!       Most of these are info to be set by the user.
!       The private_info stuct contains internal info not
!       intended to be accessed or set by the user.
!
   TYPE con_struct
      TYPE(general_info_struct)          :: general_info
      TYPE(stepping_info_struct)         :: stepping_info
      TYPE(arclength_info_struct)        :: arclength_info
      TYPE(turning_point_info_struct)    :: turning_point_info
      TYPE(pitchfork_info_struct)        :: pitchfork_info
      TYPE(hopf_info_struct)             :: hopf_info
      TYPE(phase_transition_info_struct) :: phase_transition_info
      TYPE(manifold_info_struct)         :: manifold_info
      TYPE(eigen_info_struct)            :: eigen_info

      TYPE(private_info_struct)          :: private_info
   END TYPE con_struct


   TYPE arc_scale_struct
      REAL(KIND=8) :: dx0         ! Initial scale factor for solution
      REAL(KIND=8) :: dx_fac      ! scale factor for solution
      REAL(KIND=8) :: dx_fac_max  ! maximum allowable value for arc_scale%dx_fac
      REAL(KIND=8) :: umag2       ! square of scaled solution vector magnitude
      REAL(KIND=8) :: dx_fac_old  ! scale factor from previous step
      REAL(KIND=8) :: dp_ds_limit ! dp_ds value at which to default to maximum value of arc_scale%dx_fac = arc_scale%dx_fac_max
      REAL(KIND=8) :: dp_ds_goal  ! dp_ds value to reset to by rescaling
      REAL(KIND=8) :: dp_ds_old   ! dp_ds value from previous step
      REAL(KIND=8) :: ds_fac      ! arc length scale factor
   END TYPE arc_scale_struct


   TYPE passdown_struct
      REAL(KIND=8), DIMENSION(:), POINTER :: x ! FixMe remove pointer and transform in allocatable?
      INTEGER                             :: beta
      REAL(KIND=8)                        :: reynolds, oscar, romeo, whisky, h, tol
      INTEGER                             :: bif_param, param, maxiter, ldz, num_linear_its, debug
   END TYPE passdown_struct

!****************************************************************************
!****************************************************************************
!****************************************************************************
END MODULE Loca_types
