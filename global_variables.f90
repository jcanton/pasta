MODULE global_variables
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 26/8/2013
!
   USE read_input_files       ! for the definition of the type 'program_input'
   USE dynamic_structures     ! for the definition of various vector types
   USE sparse_matrix_profiles ! for the definition of real and complex CSR matrices


   PUBLIC

   TYPE(program_input) :: p_in ! parameters read from file 'program_data.in'

   INTEGER, PARAMETER :: velCmpnnts=3 ! number of velocity components
   INTEGER            :: number_of_sides ! = MAXVAL(sides)

   LOGICAL, ALLOCATABLE,  DIMENSION(:) :: Axi, Axis
   INTEGER,      POINTER, DIMENSION(:) :: js_Axis ! For homogeneous bcs on

   TYPE(dyn_int_line),   DIMENSION(velCmpnnts) :: js_D ! List of Dirichlet nodes

   LOGICAL, ALLOCATABLE, DIMENSION(:,:)        :: Dir
   TYPE(dyn_real_line),  DIMENSION(velCmpnnts) :: bvs_D, & ! Dirichlet boundary values
                                             zero_bvs_D, &
                                              old_bvs_D

   INTEGER,      POINTER, DIMENSION(:) :: ms_2, ms_3 ! List of boundary elements
   REAL(KIND=8), POINTER, DIMENSION(:) :: c_2,  q_3  ! boundary values

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: in_bvs_D

   TYPE(CSR_MUMPS_Matrix) :: Jacobian, Mass

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: uu, u0
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: pp, p0

   REAL(KIND=8), DIMENSION(:), POINTER             :: x0 ! used for passing from C to F
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:), TARGET :: xx ! used for passing from F to C
                                                         !
   REAL(KIND=8) :: Re                                    ! when used inside F they
   INTEGER      :: beta, Nx                              ! behave in the same way

   REAL(KIND=8), DIMENSION(3) :: flow_parameters

   REAL(KIND=8), DIMENSION(3,2) :: volumeForcing

   LOGICAL :: DESINGULARIZE

   LOGICAL, ALLOCATABLE, DIMENSION(:)        :: Dir_psi
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: psi
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: zz

   INTEGER :: parMumpsJob
   INTEGER :: matrID

   INCLUDE 'mpif.h'
   INTEGER :: myRank, nProc, mpiErrC, mpiIerr

END MODULE global_variables
