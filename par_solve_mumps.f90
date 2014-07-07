!==============================================================================
!
! MODULE par_solve_mumps
! ==================
!
! Interface for MUMPS (Multifrontal Massively Parallel Solver),
! a package for solving linear systems with unsymmetric, symmetric
! positive definite and general symmetric matrices. The solution
! is based on a multifrontal technique which is a direct method
! using the LU or the LDL' factorization of the matrix. For further
! informations about MUMPS refer to the project web sites:
!
! http://www.enseeiht.fr/apo/MUMPS/
! http://www.ens-lyon.fr/~jylexcel/MUMPS
!
! The interface allows to store up to max_save factorized linear
! systems for subsequent solutions with different rhs.
!
!
! AUTHOR: Jacopo Canton (jcanton@mech.kth.se) - 12/09/2013
!
! LAST REVISION: 12/09/2013
!
!==============================================================================
MODULE par_solve_mumps

   USE global_variables
   USE sparse_matrix_profiles 
 
   IMPLICIT NONE

   INCLUDE 'dmumps_struc.h'
   INCLUDE 'zmumps_struc.h'

   ! mumps settings
   INTEGER, PARAMETER :: hostWorking = 1

   
   ! public methods
   PUBLIC :: par_mumps_master, par_mumps

   ! stored mumps instances
   INTEGER, PARAMETER :: max_save = 20
   TYPE(DMUMPS_STRUC), DIMENSION(max_save), SAVE :: id_real
   TYPE(ZMUMPS_STRUC), DIMENSION(max_save), SAVE :: id_cmpl
   INTEGER, SAVE :: realCmpl

   ! parameters
   INTEGER, PARAMETER, PUBLIC :: INITIALIZATION  =  0
   INTEGER, PARAMETER, PUBLIC :: SYMBO_FACTOR    =  1
   INTEGER, PARAMETER, PUBLIC :: NUMER_FACTOR    =  2
   INTEGER, PARAMETER, PUBLIC :: DIRECT_SOLUTION =  3
   INTEGER, PARAMETER, PUBLIC :: TRANSP_SOLUTION =  4
   INTEGER, PARAMETER, PUBLIC :: DEALLOCATION    = -1
   INTEGER, PARAMETER, PUBLIC :: FINALIZATION    =  99
   INTEGER, PARAMETER :: REAL_FLAG =  1
   INTEGER, PARAMETER :: CMPL_FLAG =  2

   ! memory parameters
   INTEGER, PARAMETER :: memRelax = 5


   ! interface for real/complex matrix and rhs
   !
   ! this is the subroutine that needs to be called by the master
   ! process
   !
   INTERFACE par_mumps_master

      MODULE PROCEDURE par_mumps_master_real,  par_mumps_master_cmpl

   END INTERFACE par_mumps_master

 CONTAINS
!==============================================================================


!------------------------------------------------------------------------------
! SUBROUTINE par_mumps_master_real
! ---------------------------------
!
! This is the subroutine that is executed only by the master
! process
!
!------------------------------------------------------------------------------
!
!   parMumpsJob ----> Parameter indicating the type of computation
!                     to be performed
!
!   matrID      ----> Parameter identifying the matrix 
!
!   A%i         ----> Pointers to the beginning of each row of A
!   A%j         ----> Column indices
!   A%e         ----> Non zero entries of the matrix REAL NUMBERS
!
!   symFlag     ----> 0 -> unsymmetric matrix
!                     1 -> symmetric definite positive matrix
!                     2 -> general symmetric matrix
!
!   rhsSol      ----> OPTIONAL vector containing the rhs in input
!                     and the solution in output, only needed when
!                     solving the linear system
!
!------------------------------------------------------------------------------
SUBROUTINE par_mumps_master_real (parMumpsJob, matrID, A, symFlag, rhsSol)

   IMPLICIT NONE
   ! input variables
   INTEGER,                INTENT(IN) :: parMumpsJob
   INTEGER,                INTENT(IN) :: matrID
   TYPE(CSR_MUMPS_Matrix), INTENT(IN) :: A
   INTEGER,                INTENT(IN) :: symFlag
   REAL(KIND=8), DIMENSION(:), TARGET, OPTIONAL :: rhsSol
   ! local variables
   INTEGER :: mpiIerr

   realCmpl = REAL_FLAG

   ! Check supplied parameters
   IF (matrID < 1 .OR. matrID > max_save) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
      WRITE(*,*) '* par_mumps_master_real           ***'
      WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** matrID out of range           ***'
      WRITE(*,*) '*** set to : ', matrID
      WRITE(*,*) '*** range is: 0 to ', max_save
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

   
   SELECT CASE (parMumpsJob)

!---->
      CASE (INITIALIZATION)
      ! initialize
      !
      WRITE(*,*) ''
      WRITE(*,*) 'PAR_SOLVE_MUMPS_real:'
      IF (hostWorking == 1) THEN
      WRITE(*,*) '             host participating to computations'
      ELSE
      WRITE(*,*) '             host NOT participating to computations'
      ENDIF
      WRITE(*,*) '             using a centralized assembled matrix'
      WRITE(*,*) '             using a centralized dense rhs_sol'
      WRITE(*,*) '             not using out-of-core'
      WRITE(*,*) '             on : error messages'
      WRITE(*,*) '             on : diagnostic, statistics and warnings'
      WRITE(*,*) '             on : global information'
      id_real(matrID)%SYM = symFlag
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)

		! diagnostic, statistics and warnings
      id_real(matrID)%ICNTL(2)  = 0
      id_real(matrID)%ICNTL(4)  = 0
      ! centralized assembled matrix
      id_real(matrID)%ICNTL(5)  = 0
      id_real(matrID)%ICNTL(18) = 0
      ! centralized dense rhs_sol
      id_real(matrID)%ICNTL(20) = 0
      id_real(matrID)%ICNTL(21) = 0
      ! deactivating out-of-core capability
      id_real(matrID)%ICNTL(22) = 0

      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Initialization successful.'

!---->
      CASE (SYMBO_FACTOR)
      ! symbolic factorization
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: starting symbolic factorization ...'
      IF (ASSOCIATED(id_real(matrID)%IRN)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** id_real already in use        ***'
         WRITE(*,*) '*** id_real: ', matrID
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_real(matrID)%N   = SIZE(A%i) - 1
      id_real(matrID)%NZ  = SIZE(A%j)
      id_real(matrID)%IRN => A%i_mumps
      id_real(matrID)%JCN => A%j
      
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_real(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Symbolic factorization computed successfully.'
      ELSE
         CALL mumps_status(id_real(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (NUMER_FACTOR)
      ! numerical factorization
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: starting numerical factorization ...'
      IF (SIZE(A%i) - 1 /= id_real(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** linear system                 ***'
         WRITE(*,*) '*** n_syst           : ', SIZE(A%i)-1
         WRITE(*,*) '*** id_real(matrID)%N: ', id_real(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      id_real(matrID)%A => A%e

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_real(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been symbolically ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_real(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Numerical factorization computed successfully.'
      ELSE
         CALL mumps_status(id_real(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (DIRECT_SOLUTION)
      ! direct solution
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: starting direct solution ...'
      IF (.NOT.PRESENT(rhsSol)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to solve a      ***'
         WRITE(*,*) '*** system without supplying the  ***'
         WRITE(*,*) '*** right had side                ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      IF (SIZE(rhsSol) /= id_real(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** right hand side               ***'
         WRITE(*,*) '*** n_rhs       : ', SIZE(rhsSol)
         WRITE(*,*) '*** id_real(matrID)%N: ', id_real(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_real(matrID)%RHS => rhsSol
      id_real(matrID)%ICNTL(9) = 1

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_real(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been numerically  ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_real(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Solution computed successfully.'
      ELSE
         CALL mumps_status(id_real(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (TRANSP_SOLUTION)
      ! solution with transposed matrix
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: starting transposed solution ...'
      IF (.NOT.PRESENT(rhsSol)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to solve a      ***'
         WRITE(*,*) '*** system without supplying the  ***'
         WRITE(*,*) '*** right had side                ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      IF (SIZE(rhsSol) /= id_real(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** right hand side               ***'
         WRITE(*,*) '*** n_rhs            : ', SIZE(rhsSol)
         WRITE(*,*) '*** id_real(matrID)%N: ', id_real(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_real(matrID)%RHS => rhsSol
      id_real(matrID)%ICNTL(9) = 0

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)

      ! Print operation result
      IF ( id_real(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been numerically  ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_real(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Transposed solution computed successfully.'
      ELSE
         CALL mumps_status(id_real(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (DEALLOCATION)
      ! deallocate space
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: starting deallocation ...'
      IF (ASSOCIATED(id_real(matrID)%IRN)) THEN
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
         CALL par_mumps (parMumpsJob, matrID)
      ELSE
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_real           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to deallocate   ***'
         WRITE(*,*) '*** the space for a matrix which  ***'
         WRITE(*,*) '*** has never been allocated      ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      WRITE (*,*) 'PAR_SOLVE_MUMPS_real: Matrix deallocated successfully.'

!---->
      CASE (FINALIZATION)
      ! finalization
      !
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)


!---->
      CASE DEFAULT
      ! error
      !
      WRITE(*,*) '*************************************'
      WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
      WRITE(*,*) '* par_mumps_master_real           ***'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** parMumpsJob                   ***'
      WRITE(*,*) '*** set to: ', parMumpsJob
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP

   END SELECT

END SUBROUTINE par_mumps_master_real


!------------------------------------------------------------------------------
! SUBROUTINE par_mumps_master_cmpl
! ---------------------------------
!
! This is the subroutine that is executed only by the master
! process
!
!------------------------------------------------------------------------------
!
!   parMumpsJob ----> Parameter indicating the type of computation
!                     to be performed
!
!   matrID      ----> Parameter identifying the matrix 
!
!   A%i         ----> Pointers to the beginning of each row of A
!   A%j         ----> Column indices
!   A%e         ----> Non zero entries of the matrix COMPLEX NUMBERS
!
!   symFlag     ----> 0 -> unsymmetric matrix
!                     1 -> symmetric definite positive matrix
!                     2 -> general symmetric matrix
!
!   rhsSol      ----> OPTIONAL vector containing the rhs in input
!                     and the solution in output, only needed when
!                     solving the linear system
!
!------------------------------------------------------------------------------
SUBROUTINE par_mumps_master_cmpl (parMumpsJob, matrID, A, symFlag, rhsSol)

   IMPLICIT NONE
   ! input variables
   INTEGER,                        INTENT(IN) :: parMumpsJob
   INTEGER,                        INTENT(IN) :: matrID
   TYPE(CSR_MUMPS_Complex_Matrix), INTENT(IN) :: A
   INTEGER,                        INTENT(IN) :: symFlag
   COMPLEX(KIND=8), DIMENSION(:), TARGET, OPTIONAL :: rhsSol
   ! local variables
   INTEGER :: mpiIerr

   realCmpl = CMPL_FLAG

   ! Check supplied parameters
   IF (matrID < 1 .OR. matrID > max_save) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
      WRITE(*,*) '* par_mumps_master_cmpl           ***'
      WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** matrID out of range           ***'
      WRITE(*,*) '*** set to : ', matrID
      WRITE(*,*) '*** range is: 0 to ', max_save
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

   
   SELECT CASE (parMumpsJob)

!---->
      CASE (INITIALIZATION)
      ! initialize
      !
      WRITE(*,*) ''
      WRITE(*,*) 'PAR_SOLVE_MUMPS_cmpl:'
      IF (hostWorking == 1) THEN
      WRITE(*,*) '             host participating to computations'
      ELSE
      WRITE(*,*) '             host NOT participating to computations'
      ENDIF
      WRITE(*,*) '             using a centralized assembled matrix'
      WRITE(*,*) '             using a centralized dense rhs_sol'
      WRITE(*,*) '             not using out-of-core'
      WRITE(*,*) '             on : error messages'
      WRITE(*,*) '             on : diagnostic, statistics and warnings'
      WRITE(*,*) '             on : global information'
      id_cmpl(matrID)%SYM = symFlag
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)

		! diagnostic, statistics and warnings
      id_cmpl(matrID)%ICNTL(2)  = 0
      id_cmpl(matrID)%ICNTL(4)  = 0
      ! centralized assembled matrix
      id_cmpl(matrID)%ICNTL(5)  = 0
      id_cmpl(matrID)%ICNTL(18) = 0
      ! centralized dense rhs_sol
      id_cmpl(matrID)%ICNTL(20) = 0
      id_cmpl(matrID)%ICNTL(21) = 0
      ! deactivating out-of-core capability
      id_cmpl(matrID)%ICNTL(22) = 0

      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: Initialization successful.'

!---->
      CASE (SYMBO_FACTOR)
      ! symbolic factorization
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: starting symbolic factorization ...'
      IF (ASSOCIATED(id_cmpl(matrID)%IRN)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** id_cmpl already in use        ***'
         WRITE(*,*) '*** id_cmpl: ', matrID
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_cmpl(matrID)%N   = SIZE(A%i) - 1
      id_cmpl(matrID)%NZ  = SIZE(A%j)
      id_cmpl(matrID)%IRN => A%i_mumps
      id_cmpl(matrID)%JCN => A%j
      
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_cmpl(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: Symbolic factorization computed successfully.'
      ELSE
         CALL mumps_status(id_cmpl(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (NUMER_FACTOR)
      ! numerical factorization
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: starting numerical factorization ...'
      IF (SIZE(A%i) - 1 /= id_cmpl(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** linear system                 ***'
         WRITE(*,*) '*** n_syst           : ', SIZE(A%i)-1
         WRITE(*,*) '*** id_cmpl(matrID)%N: ', id_cmpl(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      id_cmpl(matrID)%A => A%e

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_cmpl(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been symbolically ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_cmpl(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: Numerical factorization computed successfully.'
      ELSE
         CALL mumps_status(id_cmpl(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (DIRECT_SOLUTION)
      ! direct solution
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: starting direct solution ...'
      IF (.NOT.PRESENT(rhsSol)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to solve a      ***'
         WRITE(*,*) '*** system without supplying the  ***'
         WRITE(*,*) '*** right had side                ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      IF (SIZE(rhsSol) /= id_cmpl(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** right hand side               ***'
         WRITE(*,*) '*** n_rhs       : ', SIZE(rhsSol)
         WRITE(*,*) '*** id_cmpl(matrID)%N: ', id_cmpl(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_cmpl(matrID)%RHS => rhsSol
      id_cmpl(matrID)%ICNTL(9) = 1

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)
      ! Print operation result
      IF ( id_cmpl(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been numerically  ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_cmpl(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: Solution computed successfully.'
      ELSE
         CALL mumps_status(id_cmpl(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (TRANSP_SOLUTION)
      ! solution with transposed matrix
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: starting transposed solution ...'
      IF (.NOT.PRESENT(rhsSol)) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to solve a      ***'
         WRITE(*,*) '*** system without supplying the  ***'
         WRITE(*,*) '*** right had side                ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      IF (SIZE(rhsSol) /= id_cmpl(matrID)%N) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** incoherent size of supplied   ***'
         WRITE(*,*) '*** right hand side               ***'
         WRITE(*,*) '*** n_rhs            : ', SIZE(rhsSol)
         WRITE(*,*) '*** id_cmpl(matrID)%N: ', id_cmpl(matrID)%N
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF

      id_cmpl(matrID)%RHS => rhsSol
      id_cmpl(matrID)%ICNTL(9) = 0

      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)

      ! Print operation result
      IF ( id_cmpl(matrID)%INFOG(1) .EQ. -3 ) THEN
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** the problem referred to by    ***'
         WRITE(*,*) '*** matrID: ', matrID
         WRITE(*,*) '*** has not yet been numerically  ***'
         WRITE(*,*) '*** factorized                    ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ELSEIF ( id_cmpl(matrID)%INFOG(1) .EQ. 0 ) THEN
         WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: Transposed solution computed successfully.'
      ELSE
         CALL mumps_status(id_cmpl(matrID)%INFOG(1))
      ENDIF

!---->
      CASE (DEALLOCATION)
      ! deallocate space
      !
      WRITE (*,*) 'PAR_SOLVE_MUMPS_cmpl: starting deallocation ...'
      IF (ASSOCIATED(id_cmpl(matrID)%IRN)) THEN
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
         CALL par_mumps (parMumpsJob, matrID)
      ELSE
         WRITE(*,*) '*************************************'
         WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
         WRITE(*,*) '* par_mumps_master_cmpl           ***'
         WRITE(*,*) '*** ERROR:                        ***'
         WRITE(*,*) '*** you''re trying to deallocate   ***'
         WRITE(*,*) '*** the space for a matrix which  ***'
         WRITE(*,*) '*** has never been allocated      ***'
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         STOP
      ENDIF
      WRITE (*,*) 'PAR_SOLVE_MUMPS: Matrix deallocated successfully.'

!---->
      CASE (FINALIZATION)
      ! finalization
      !
      CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)
      CALL par_mumps (parMumpsJob, matrID)


!---->
      CASE DEFAULT
      ! error
      !
      WRITE(*,*) '*************************************'
      WRITE(*,*) '* PAR SOLVE MUMPS                 ***'
      WRITE(*,*) '* par_mumps_master_cmpl           ***'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** parMumpsJob                   ***'
      WRITE(*,*) '*** set to: ', parMumpsJob
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP

   END SELECT

END SUBROUTINE par_mumps_master_cmpl

!------------------------------------------------------------------------------
! SUBROUTINE par_mumps
! ---------------------------------
!
! This is the subroutine that is called by both the master process
! and the slaves
!
!------------------------------------------------------------------------------
!
!   parMumpsJob ----> Parameter indicating the type of computation
!                     to be performed
!
!   matrID      ----> Parameter identifying the matrix 
!
!------------------------------------------------------------------------------
SUBROUTINE par_mumps (parMumpsJob, matrID)

   IMPLICIT NONE
   ! input variables
   INTEGER, INTENT(IN) :: parMumpsJob
   INTEGER, INTENT(IN) :: matrID
   ! local variables
   INTEGER :: myRank, mpiIerr
    
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpiIerr)

   CALL MPI_BCAST (parMumpsJob, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
   CALL MPI_BCAST (matrID,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
   CALL MPI_BCAST (realCmpl,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)

   SELECT CASE (parMumpsJob)

!---->
      CASE (INITIALIZATION)
      ! initialize
      !
      IF (realCmpl == REAL_FLAG) THEN
         id_real(matrID)%COMM = MPI_COMM_WORLD
         id_real(matrID)%PAR  = hostWorking
         CALL MPI_BCAST (id_real(matrID)%SYM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)

         id_real(matrID)%JOB = -1
         CALL DMUMPS(id_real(matrID))

         ! output stream for error messages
         id_real(matrID)%ICNTL(1) = 6
         ! output stream for diagnostic, statistics and warnings
         id_real(matrID)%ICNTL(2) = 0
         ! output stream for global information
         id_real(matrID)%ICNTL(3) = 0
         ! level of printing
         id_real(matrID)%ICNTL(4) = 1
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         id_cmpl(matrID)%COMM = MPI_COMM_WORLD
         id_cmpl(matrID)%PAR  = hostWorking
         CALL MPI_BCAST (id_cmpl(matrID)%SYM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)

         id_cmpl(matrID)%JOB = -1
         CALL ZMUMPS(id_cmpl(matrID))

         ! output stream for error messages
         id_cmpl(matrID)%ICNTL(1) = 6
         ! output stream for diagnostic, statistics and warnings
         id_cmpl(matrID)%ICNTL(2) = 0
         ! output stream for global information
         id_cmpl(matrID)%ICNTL(3) = 0
         ! level of printing
         id_cmpl(matrID)%ICNTL(4) = 1
      ENDIF

!---->
      CASE (SYMBO_FACTOR)
      ! symbolic factorization
      !
      IF (realCmpl == REAL_FLAG) THEN
         CALL MPI_BCAST (id_real(matrID)%N,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
         CALL MPI_BCAST (id_real(matrID)%NZ,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)

         id_real(matrID)%JOB = 1
         CALL DMUMPS(id_real(matrID))
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         CALL MPI_BCAST (id_cmpl(matrID)%N,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
         CALL MPI_BCAST (id_cmpl(matrID)%NZ,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiIerr)
         CALL MPI_BARRIER (MPI_COMM_WORLD, mpiIerr)

         id_cmpl(matrID)%JOB = 1
         CALL ZMUMPS(id_cmpl(matrID))
      ENDIF
      
!---->
      CASE (NUMER_FACTOR)
      ! numerical factorization
      !
      IF (realCmpl == REAL_FLAG) THEN
         id_real(matrID)%ICNTL(23) = id_real(matrID)%INFOG(16) * memRelax
         !write(*,*) 'infog(16) = ', id_real(matrID)%INFOG(16)
         !write(*,*) 'icntl(23) = ', id_real(matrID)%ICNTL(23)
         id_real(matrID)%JOB  = 2
         CALL DMUMPS(id_real(matrID))
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         id_cmpl(matrID)%ICNTL(23) = id_cmpl(matrID)%INFOG(16) * memRelax
         id_cmpl(matrID)%JOB  = 2
         CALL ZMUMPS(id_cmpl(matrID))
      ENDIF

!---->
      CASE (DIRECT_SOLUTION)
      ! direct solution
      !
      IF (realCmpl == REAL_FLAG) THEN
         id_real(matrID)%JOB  = 3
         CALL DMUMPS(id_real(matrID))
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         id_cmpl(matrID)%JOB  = 3
         CALL ZMUMPS(id_cmpl(matrID))
      ENDIF

!---->
      CASE (TRANSP_SOLUTION)
      ! solution with transposed matrix
      !
      IF (realCmpl == REAL_FLAG) THEN
         id_real(matrID)%JOB  = 3
         CALL DMUMPS(id_real(matrID))
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         id_cmpl(matrID)%JOB  = 3
         CALL ZMUMPS(id_cmpl(matrID))
      ENDIF

!---->
      CASE (DEALLOCATION)
      ! deallocate space
      !
      IF (realCmpl == REAL_FLAG) THEN
         id_real(matrID)%JOB  = -2
         CALL DMUMPS(id_real(matrID))
      ELSEIF (realCmpl == CMPL_FLAG) THEN
         id_cmpl(matrID)%JOB  = -2
         CALL ZMUMPS(id_cmpl(matrID))
      ENDIF

!---->
      CASE (FINALIZATION)
      ! finalization
      !
      CALL MPI_FINALIZE(mpiIerr)
      WRITE(*,*) '============================================================'
      WRITE(*,*) '             END of program for process.', myRank
      WRITE(*,*) '============================================================'
      STOP

   END SELECT

END SUBROUTINE par_mumps

!------------------------------------------------------------------------------
! SUBROUTINE mumps_status
! -----------------------
!
! This subroutine print the description of the MUMPS status
! after the last operation given by the integer 'ierr'.
!
!------------------------------------------------------------------------------
!
!   ierr        ----> MUMPS status returned in INFOG(1)
!
!------------------------------------------------------------------------------

SUBROUTINE  mumps_status (ierr)

   INTEGER,  INTENT(IN) :: ierr

   ! Print result of the factorization/solution operation
   SELECT CASE(ierr)

      CASE(-2)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): NZ is out of range.'
 
      CASE(-5,-7,-13)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): Out of memory.'
 
      CASE(-6)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): Matrix singular in structure.'
 
      CASE(-8,-14,-15,-18)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): MUMPS parameter MAXIS '// &
                     'is too small.'
 
      CASE(-9,-11,-12,-19)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): MUMPS parameter MAXS '// &
                     'is too small.'
 
      CASE(-10)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): Matrix numerically singular.'
 
      CASE(-16)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): N is out of range.'
 
      CASE(-22)
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): A supplied pointer array '// &
                     'is either not associated or has insufficient size.'
 
      CASE DEFAULT
      WRITE (*,*) 'PAR_SOLVE_MUMPS(',ierr,'): Unknown error.'
 
   END SELECT 

   STOP

END SUBROUTINE  mumps_status

!==============================================================================

END MODULE  par_solve_mumps
