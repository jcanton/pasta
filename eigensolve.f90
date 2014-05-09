MODULE EigenSolve
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 23/8/2013
!
   USE sparse_matrix_profiles
   USE sparse_matrix_operations
   USE par_solve_mumps
   USE ISO_C_BINDING
   USE miscellaneous_subroutines

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE eigensComplexShiftInvert(nev, maxit, tol, sigma, A, M, directAdjoint, eigenvalues, eigenvectors)

!     Simple program to illustrate the idea of reverse communication
!     in shift and invert mode for a generalized complex nonsymmetric 
!     eigenvalue problem.
!
!     ... Suppose we want to solve A*x = lambda*M*x in shift-invert mode,
!
!     ... where the shift sigma is a complex number.
!
!     ... OP = inv[A-SIGMA*M]*M
!
!     ... Use mode 3 of ZNAUPD .
!
!\Routines called:
!     znaupd   ARPACK reverse communication interface routine.
!     zneupd   ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dlapy2   LAPACK routine to compute sqrt(x**2 + y**2) carefully (AVOID UNECESSARY OVERFLOW).
!     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
!
!-----------------------------------------------------------------------

   IMPLICIT NONE

!  %-----------------%
!  | Input variables |
!  %-----------------%
   INTEGER,         INTENT(IN) :: nev, maxit
   REAL(KIND=8),    INTENT(IN) :: tol
   COMPLEX(KIND=8), INTENT(IN) :: sigma

   TYPE(CSR_MUMPS_Complex_Matrix), INTENT(IN) :: A, M

   INTEGER, INTENT(IN) :: directAdjoint ! 1 = direct problem, 2 = adjoint problem

!  %------------------%
!  | Output variables |
!  %------------------%
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: eigenvalues
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: eigenvectors

!  %--------------%
!  | Local Arrays |
!  %--------------%
   INTEGER, DIMENSION(11) :: iparam
   INTEGER, DIMENSION(14) :: ipntr
   
   LOGICAL, DIMENSION(:), ALLOCATABLE :: lselect
   
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: Ax, Mx, d, workd, resid, workev, workl
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: v
   
   REAL(KIND=8),    DIMENSION(:),   ALLOCATABLE :: rwork, tmp_r, tmp_i
   REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: rd

   TYPE(CSR_MUMPS_Complex_Matrix) :: AmsM ! AmsM = A - shift*M

!  %---------------%
!  | Local Scalars |
!  %---------------%
   COMPLEX(KIND=8)  :: shift
   CHARACTER(LEN=1) :: bmat
   CHARACTER(LEN=2) :: which

   INTEGER          :: ido, n, ncv, lworkl, info, i, j, nconv

   LOGICAL          :: rvec 

   !LOGICAL, SAVE    :: symbolic_init=.FALSE.

!  %-----------------------------%
!  | BLAS & LAPACK routines used |
!  %-----------------------------%
   REAL(KIND=8), EXTERNAL ::  dznrm2 , dlapy2 

!------------------------------------------------------------------------------
!
!  %-----------------------%
!  | Executable statements |
!  %-----------------------%
!
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Start of eigensComplexShiftInvert'
   WRITE(*,*) ''

   IF ( directAdjoint == 1 ) THEN
      ! direct problem
      shift = sigma
   ELSEIF ( directAdjoint == 2 ) THEN
      ! adjoint problem
      shift = DCONJG(sigma)
   ELSE
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** ERROR                         ***'
      WRITE(*,*) '*** wrong directAdjoint parameter ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF


   IF ( sigma /= CMPLX(0d0,0d0,KIND=8) ) THEN

      ! Create the matrix [A-shift*M] and store it in position 4 of the MUMPS
      ! array id, then deallocate it to save RAM:

      ALLOCATE( AmsM%i      (SIZE(A%i))       ); AmsM%i       = A%i
      ALLOCATE( AmsM%i_mumps(SIZE(A%i_mumps)) ); AmsM%i_mumps = A%i_mumps
      ALLOCATE( AmsM%j      (SIZE(A%j))       ); AmsM%j       = A%j
      ALLOCATE( AmsM%e      (SIZE(A%e))       ); AmsM%e       = A%e

!      IF ( .NOT.symbolic_init ) THEN
         CALL par_mumps_master (INITIALIZATION, 4, AmsM, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   4, AmsM, 0)
!         symbolic_init=.TRUE.
!      ENDIF

      ! AmsM = [A-shift*M]
      CALL zAlinB_s ( CMPLX(1d0,0d0,KIND=8), A%e, -shift, M%e,  AmsM%e )

      CALL par_mumps_master (NUMER_FACTOR, 4, AmsM, 0)

      DEALLOCATE( AmsM%i, AmsM%i_mumps, AmsM%j, AmsM%e )

      WRITE(*,*) ''
      WRITE(*,*) '    Shifted matrix assembled correctly.'
      WRITE(*,*) ''

   ELSE

      ! [A-shift*M] = A
      ! store it into position 4 of the MUMPS array id
!      IF ( .NOT.symbolic_init ) THEN
         CALL par_mumps_master (INITIALIZATION, 4, A, 0)
         CALL par_mumps_master (SYMBO_FACTOR,   4, A, 0)
!         symbolic_init=.TRUE.
!      ENDIF
      CALL par_mumps_master (NUMER_FACTOR, 4, A, 0)
   
      WRITE(*,*) ''
      WRITE(*,*) '    No shift requested.'
      WRITE(*,*) ''

   END IF
!
!  %----------------------------------------------------%
!  | The number N is the dimension of the matrix.  A    |
!  | generalized eigenvalue problem is solved (BMAT =   |
!  | 'G').  NEV is the number of eigenvalues (closest   |
!  | to SIGMAR) to be approximated.  Since the          |
!  | shift-invert mode is used,  WHICH is set to 'LM'.  |
!  | The user can modify NEV, NCV, SIGMA to solve       |
!  | problems of different sizes, and to get different  |
!  | parts of the spectrum.  However, The following     |
!  | condition must be satisfied:                       |
!  |               NEV + 2 <= NCV                       |
!  %----------------------------------------------------%
!
   n   = SIZE( A%i ) - 1
   
   ncv   = 2 * nev + 2
!   ncv   = 10 * nev + 2
   
   bmat  = 'G'
   which = 'LM'
!
!  %-----------------------------------------------------%
!  | The work array WORKL is used in ZNAUPD  as          |
!  | workspace.  Its dimension LWORKL is set as          |
!  | illustrated below.  The parameter TOL determines    |
!  | the stopping criterion. If TOL<=0, machine          |
!  | precision is used.  The variable IDO is used for    |
!  | reverse communication, and is initially set to 0.   |
!  | Setting INFO=0 indicates that a random vector is    |
!  | generated in ZNAUPD  to start the Arnoldi iteration.|
!  %-----------------------------------------------------%
!
   lworkl = 3*ncv**2+5*ncv 
   ido    = 0
   info   = 0
!
!  %---------------------------------------------------%
!  | This program uses exact shifts with respect to    |
!  | the current Hessenberg matrix (IPARAM(1) = 1).    |
!  | IPARAM(3) specifies the maximum number of Arnoldi |
!  | iterations allowed. Mode 3 of ZNAUPD  is used     |
!  | (IPARAM(7) = 3).  All these options can be        |
!  | changed by the user. For details see the          |
!  | documentation in ZNAUPD .                         |
!  %---------------------------------------------------%
!
   iparam(1) = 1     ! ishfts = use exact shifts
   iparam(3) = maxit ! maxit  = max number of Arnoldi iterations allowed
   iparam(7) = 3     ! mode   = ZNAUPD mode
   
   ALLOCATE(lselect(ncv))
   ALLOCATE(Ax(n), Mx(n), workd(3*n), resid(n), tmp_r(n), tmp_i(n)) 
   ALLOCATE(d(ncv), v(n,ncv), workev(2*ncv), workl(lworkl))
   ALLOCATE(rwork(ncv), rd(ncv,3))
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
   i = 0
   DO
      i = i+1
   
      WRITE(*,*) 'Eigensolve, iteration n. ', i
!
!     %----------------------------------------------%
!     | Repeatedly call the routine ZNAUPD  and take |
!     | actions indicated by parameter IDO until     |
!     | either convergence is indicated or maxit     |
!     | has been exceeded.                           |
!     %----------------------------------------------%
!
      CALL znaupd  ( ido  , bmat , n     , which, &
                     nev  , tol  , resid , ncv,   &
                     v    , n    , iparam, ipntr, &
                     workd, workl, lworkl, rwork, info )

      ! WRITE(*,*) 'Eigensolve, requested action ', ido

      IF (ido .eq. -1) THEN
!
!        %-------------------------------------------%
!        | Perform  y <--- OP*x = inv[A-sigma*M]*M*x |
!        | to force starting vector into the range   |
!        | of OP.   The user should supply his/her   |
!        | own matrix vector multiplication routine  |
!        | and a linear system solver.  The matrix   |
!        | vector multiplication routine should take |
!        | workd(ipntr(1)) as the input. The final   |
!        | result should be returned to              |
!        | workd(ipntr(2)).                          |
!        %-------------------------------------------%
!
         IF ( directAdjoint == 1 ) THEN
            ! direct problem
            CALL zAtimx (Mx, M%e, M%j, M%i, workd(ipntr(1):ipntr(1)+n-1))
            CALL par_mumps_master (DIRECT_SOLUTION, 4, A, 0, Mx)
         ELSE
            ! adjoint problem
            CALL zAtimx_T (Mx, M%e, M%j, M%i, workd(ipntr(1):ipntr(1)+n-1))
            CALL par_mumps_master (TRANSP_SOLUTION, 4, A, 0, Mx)
         ENDIF

         workd(ipntr(2):ipntr(2)+n-1) = Mx
!
!        %------------------------------------------%
!        | L O O P   B A C K to call ZNAUPD  again. |
!        %------------------------------------------%
!
         CYCLE

      ELSE IF ( ido .eq. 1) THEN
!
!        %-----------------------------------------%
!        | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
!        | M*x has been saved in workd(ipntr(3)).  |
!        | The user only need the linear system    |
!        | solver here that takes workd(ipntr(3))  |
!        | as input, and returns the result to     |
!        | workd(ipntr(2)).                        |
!        %-----------------------------------------%
!
         workd(ipntr(2):ipntr(2)+n-1) = workd(ipntr(3):ipntr(3)+n-1)

         IF ( directAdjoint == 1 ) THEN
            ! direct problem
            CALL par_mumps_master (DIRECT_SOLUTION, 4, A, 0, workd(ipntr(2):ipntr(2)+n-1))
         ELSE
            ! adjoint problem
            CALL par_mumps_master (TRANSP_SOLUTION, 4, A, 0, workd(ipntr(2):ipntr(2)+n-1))
         ENDIF
!
!        %------------------------------------------%
!        | L O O P   B A C K to call ZNAUPD  again. |
!        %------------------------------------------%
!
         CYCLE

      ELSE IF ( ido .eq. 2) THEN
!
!        %---------------------------------------------%
!        |          Perform  y <--- M*x                |
!        | Need matrix vector multiplication routine   |
!        | here that takes workd(ipntr(1)) as input    |
!        | and returns the result to workd(ipntr(2)).  |
!        %---------------------------------------------%
!
         IF ( directAdjoint == 1 ) THEN
            ! direct problem
            CALL zAtimx (workd(ipntr(2):ipntr(2)+n-1), M%e, M%j, M%i, workd(ipntr(1):ipntr(1)+n-1))
         ELSE
            ! adjoint problem
            CALL zAtimx_T (workd(ipntr(2):ipntr(2)+n-1), M%e, M%j, M%i, workd(ipntr(1):ipntr(1)+n-1))
         ENDIF
!
!        %------------------------------------------%
!        | L O O P   B A C K to call ZNAUPD  again. |
!        %------------------------------------------%
!
         CYCLE
         
      ELSE IF ( ido .eq. 99) THEN
      
         WRITE(*,*) 'Done eigensolution, exit info: ', info
         EXIT

      END IF 
         
   END DO  ! End main loop
! 
!  %-----------------------------------------%
!  | Either we have convergence, or there is |
!  | an error.                               |
!  %-----------------------------------------%
!  
   IF ( info .lt. 0 ) THEN ! There is an error.
!
!     %----------------------------%
!     |  Error message, check the  |
!     |  documentation in ZNAUPD   |
!     %----------------------------%
!
      WRITE(*,*) ' Error with ZNAUPD, info = ', info
      STOP       ' Check the documentation of ZNAUPD.'
!
   ELSE ! We made it.
!
!     %-------------------------------------------%
!     | No fatal errors occurred.                 |
!     | Post-Process using ZNEUPD .               |
!     |                                           |
!     | Computed eigenvalues may be extracted.    |
!     |                                           |
!     | Eigenvectors may also be computed now if  |
!     | desired.  (indicated by rvec = .true.)    |
!     %-------------------------------------------%
!
      rvec = .true.

      CALL zneupd (rvec, 'A', lselect, d, v, n, shift, &
                   workev, bmat, n, which, nev, tol, resid, ncv, v, &
                   n, iparam, ipntr, workd, workl, lworkl, rwork, info)
!
!     %----------------------------------------------%
!     | Eigenvalues are returned in the one          |
!     | dimensional array D.  The corresponding      |
!     | eigenvectors are returned in the first NCONV |
!     | (=IPARAM(5)) columns of the two dimensional  |
!     | array V if requested.  Otherwise, an         |
!     | orthogonal basis for the invariant subspace  |
!     | corresponding to the eigenvalues in D is     |
!     | returned in V.                               |
!     %----------------------------------------------%
!
      IF ( info .ne. 0) THEN
! 
!        %------------------------------------%
!        | Error condition:                   |
!        | Check the documentation of ZNEUPD. |
!        %------------------------------------%
!
         WRITE(*,*) ' Error with ZNEUPD, info = ', info
         STOP       ' Check the documentation of ZNEUPD.'

      ELSE

         ! Number of converged eigenvalues
         nconv = iparam(5)
          
         ALLOCATE(eigenvalues(nconv))
         ALLOCATE(eigenvectors(n,nconv))
                           
         eigenvalues  = d(1:nconv)
         eigenvectors = v(:,1:nconv)
          
             
!        ! (a) Metodo di calcolo con inversione della matrice NON IMPLEMENTATO
!        DO j = 1, nconv
!
!!!
! todo :: add the matrix vector multiplication
!        v_temp = M * v(1:n,j), M is the mass matrix of the full problem
!
!           CALL zAtimx (y,      a,   ja,  ia,  x)
!           CALL zAtimx (v_temp, M%e, M%j, M%i, v(1:n,j))
!
!
!!!
! todo :: call the MUMPS complex solver
!        Mx = J^-1 v_temp
! se e` necessario fare il loop con l'inversione della matrice bisogna
! fattorizzare numericamente la matrice A e salvarla in un'istanza dell'array id
! di MUMPS
!        
!           Ax = v(1:n,j)
!
!           CALL zaxpy (n, -d(j), Mx, 1, Ax, 1) ! Ax = A v_j - d(j) M v_j 
!            
!           rd(j,1) = DBLE (d(j))
!           rd(j,2) = AIMAG (d(j))
!           rd(j,3) = dznrm2 (n, Ax, 1)
!           rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
!            
!        END DO
!
!        %-----------------------------%
!        | Display computed residuals. |
!        %-----------------------------%
!
!        CALL dmout (6, nconv, 3, rd, ncv, -6, 'Ritz values (Real, Imag) and direct residuals')
!            
 
         ! (b)
         ! Metodo di calcolo senza inversione della matrice
         !
         DO j = 1, nconv

            ! Mx = M * v(1:n,j)
            IF ( directAdjoint == 1 ) THEN
               ! direct problem
               CALL zAtimx (Mx, M%e, M%j, M%i, v(1:n,j))
            ELSE
               ! adjoint problem
               CALL zAtimx_T (Mx, M%e, M%j, M%i, v(1:n,j))
            ENDIF

            ! Ax = A v(1:n,j)
            IF ( directAdjoint == 1 ) THEN
               ! direct problem
               CALL zAtimx (Ax, A%e, A%j, A%i, v(1:n,j))
            ELSE
               ! adjoint problem
               CALL zAtimx_T (Ax, A%e, A%j, A%i, v(1:n,j))
            ENDIF

            CALL zaxpy (n, -d(j), Mx, 1, Ax, 1) ! Ax = A v_j - d(j) M v_j 
             
            rd(j,1) = DBLE (d(j))
            rd(j,2) = AIMAG (d(j))
            rd(j,3) = dznrm2 (n, Ax, 1)
            rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
             
         END DO
!
!        %-----------------------------%
!        | Display computed residuals. |
!        %-----------------------------%
!
         CALL dmout (6, nconv, 3, rd, ncv, -6, 'Ritz values (Real, Imag) and direct residuals')

      END IF
!
!     %-------------------------------------------%
!     | Print additional convergence information. |
!     %-------------------------------------------%
!
      IF ( info .eq. 1) THEN
          WRITE(*,*) ''
          WRITE(*,*) ' Maximum number of iterations reached.'
          WRITE(*,*) ''
      ELSE IF ( info .eq. 3) THEN
          WRITE(*,*) ''
          WRITE(*,*) ' No shifts could be applied during implicit',&
                     ' Arnoldi update, try increasing NCV.'
          WRITE(*,*) ''
      END IF

      WRITE(*,*) ''
      WRITE(*,*) '====== '
      WRITE(*,*) ' '
      WRITE(*,*) ' Size of the matrix is ', n
      WRITE(*,*) ' The number of Ritz values requested is ', nev
      WRITE(*,*) ' The number of Arnoldi vectors generated (NCV) is ', ncv
      WRITE(*,*) ' What portion of the spectrum: ', which
      WRITE(*,*) ' The number of converged Ritz values is ', nconv
      WRITE(*,*) ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
      WRITE(*,*) ' The number of OP*x is ', iparam(9)
      WRITE(*,*) ' The convergence criterion is ', tol
      WRITE(*,*) ' The number of main loop iterations is ', i
      WRITE(*,*) ' Shift used: ', shift
      WRITE(*,*) ' '

   END IF

   CALL par_mumps_master (DEALLOCATION, 4, A, 0)
   DEALLOCATE(lselect)
   DEALLOCATE(Ax, Mx, workd, resid, tmp_r, tmp_i)
   DEALLOCATE(d, v, workev, workl)
   DEALLOCATE(rwork, rd)

END SUBROUTINE eigensComplexShiftInvert

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

SUBROUTINE Save_eigenvalues (eigenvalues, file_name)

!-----------------------------------------------------------
   IMPLICIT NONE

   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: eigenvalues
   CHARACTER(*), OPTIONAL :: file_name
   
   ! local variables
   INTEGER :: i, nev
!-----------------------------------------------------------

   IF (PRESENT(file_name)) THEN
      OPEN (UNIT = 19, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')
      WRITE(*,*) '--> Writing eigenvalues file: ' // trim(file_name) // ' ...'
   ELSE
      OPEN (UNIT = 19, FILE = 'eigenvectors', FORM = 'formatted', STATUS = 'unknown')
      WRITE(*,*) '--> Writing eigenvalues file: eigenvalues ...'
   END IF

   nev = SIZE(eigenvalues)

   DO i = 1, nev
      ! good format for gnuplot and octave
      WRITE(19,'(g17.10,1x,g17.10)') DBLE(eigenvalues(i)), AIMAG(eigenvalues(i))
   END DO
   CLOSE(19)

   WRITE (*,*) '    Done.'

END SUBROUTINE Save_eigenvalues
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

SUBROUTINE Save_eigenvectors (eigenvectors, file_name)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/5/2013
!
!-----------------------------------------------------------
   IMPLICIT NONE

   COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: eigenvectors
   CHARACTER(*), OPTIONAL :: file_name
   
   ! local variables
   INTEGER :: i, nev, Nx
!-----------------------------------------------------------
   

   IF (PRESENT(file_name)) THEN
      OPEN (UNIT = 19, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')
      WRITE(*,*) '--> Writing eigenvector file: ' // trim(file_name) // ' ...'
   ELSE
      OPEN (UNIT = 19, FILE = 'eigenvectors', FORM = 'formatted', STATUS = 'unknown')
      WRITE(*,*) '--> Writing eigenvector file: eigenvectors ...'
   END IF

   Nx  = SIZE(eigenvectors,1)
   nev = SIZE(eigenvectors,2)

   WRITE (19,*) nev, Nx
   
   DO i = 1, Nx
      WRITE (19,*) eigenvectors(i,:)
   END DO
   
   CLOSE (19)

   WRITE (*,*) '    Done.'

END SUBROUTINE Save_eigenvectors

!---------------------------------------------------------------------------

SUBROUTINE read_eigenvector (Nx, nevRead, filenm, filenmLen, &
                                                eigenvectorRe, eigenvectorIm) &
   BIND(C, NAME='read_eigenvector')
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/5/2013
!
!-----------------------------------------------------------
   USE ISO_C_BINDING

   IMPLICIT NONE
   ! input variables
   INTEGER(KIND=C_INT), VALUE :: Nx
   INTEGER(KIND=C_INT), VALUE :: nevRead
   CHARACTER(KIND=C_CHAR)     :: filenm
   INTEGER(KIND=C_INT), VALUE :: filenmLen
   ! output variables
   REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: eigenvectorRe, eigenvectorIm
   
   ! local variables
   INTEGER :: i, nev, NxRead
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: eigenvectors
!-----------------------------------------------------------

   WRITE(*,*) '--> Reading eigenvector file: ' // filenm(1:filenmLen) // ' ...'
   WRITE(*,*) '    eigenvector number: ', nevRead

   OPEN (UNIT = 19, FILE = filenm(1:filenmLen))

   READ (19,*) nev, NxRead
   
   IF ( Nx /= NxRead ) THEN
      WRITE(*,*) '    inconsistent dimensions:'
      WRITE(*,*) '    Nx = ', Nx, ' Nx read from this file: ', NxRead
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

   WRITE(*,*) '    Nx read from this file: ', NxRead

   ALLOCATE( eigenvectors(Nx,nev) )

   DO i = 1, Nx
      READ (19,*) eigenvectors(i,:)
   ENDDO

   CLOSE(19)

   eigenvectorRe = DBLE ( eigenvectors(:,nevRead) )
   eigenvectorIm = AIMAG( eigenvectors(:,nevRead) )

   DEALLOCATE( eigenvectors )

   WRITE(*,*) '    Done.'

END SUBROUTINE read_eigenvector


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


SUBROUTINE structuralSensitivity(M, left_eigenvector, right_eigenvector, &
                                                velCmpnnts, np,  structuralsens)
!
!     Compute the structural sensitivity |v| |u| / (v,u).
!
!     ... M: mass matrix
!
!-----------------------------------------------------------------------

   IMPLICIT NONE

!  %-----------------%
!  | Input variables |
!  %-----------------%

   TYPE(CSR_MUMPS_Complex_Matrix),  INTENT(IN) :: M
   COMPLEX(KIND=8), DIMENSION(:),   INTENT(IN) :: left_eigenvector, right_eigenvector
   INTEGER :: velCmpnnts, np

!  %------------------%
!  | Output variables |
!  %------------------%
   REAL(KIND=8),    DIMENSION(:),  ALLOCATABLE :: structuralsens   

!  %--------------%
!  | Local Arrays |
!  %--------------%

   COMPLEX(KIND=8), DIMENSION(:),    ALLOCATABLE :: temp
   COMPLEX(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: u_lvect, u_rvect

!  %---------------%
!  | Local Scalars |
!  %---------------%

   INTEGER :: Nx
   REAL(KIND=8) :: normalization

!  %--------------------------%
!  | Executable instructions  |
!  %--------------------------%

   Nx   = SIZE( M%i ) - 1
   
   ALLOCATE(temp(Nx))
   ALLOCATE(u_lvect(velCmpnnts,np), u_rvect(velCmpnnts,np))
   ALLOCATE(structuralsens(np))


   ! Compute the denominator
   CALL zAtimx (temp, M%e, M%j, M%i, right_eigenvector)
   
   normalization = ABS(SUM(CONJG(left_eigenvector) * temp))
   
   ! Compute the structural sensitivity
   CALL extract_cmplx (left_eigenvector,   u_lvect)
   CALL extract_cmplx (right_eigenvector,  u_rvect)
   
   structuralsens = SQRT(  ABS(u_lvect(1,:)*u_rvect(1,:))**2 &
                         + ABS(u_lvect(1,:)*u_rvect(2,:))**2 &
                         + ABS(u_lvect(1,:)*u_rvect(3,:))**2 &
                         + ABS(u_lvect(2,:)*u_rvect(1,:))**2 &
                         + ABS(u_lvect(2,:)*u_rvect(2,:))**2 &
                         + ABS(u_lvect(2,:)*u_rvect(3,:))**2 &
                         + ABS(u_lvect(3,:)*u_rvect(1,:))**2 &
                         + ABS(u_lvect(3,:)*u_rvect(2,:))**2 &
                         + ABS(u_lvect(3,:)*u_rvect(3,:))**2)/ normalization



   DEALLOCATE(temp)
   DEALLOCATE(u_lvect, u_rvect)

END SUBROUTINE structuralSensitivity

!==============================================================================

END MODULE EigenSolve  
