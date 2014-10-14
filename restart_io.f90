MODULE restart_io
!
! Contains two subroutine for reading and writing velocity
! and pressure fields to file.
!
! - write_restart
! - read_restart
!
! plus one subroutine to write a restart file for the projection
! program by Quartapelle and Guermond
!
! - write_QP_restart
!
!==============================================================================

USE global_variables
USE prep_mesh_p1p2_sp ! for some global variables as np
USE miscellaneous_subroutines  ! for collect and extract subroutines

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE write_restart(x, param, step_num, max_steps, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 2/10/2014
!
! - x         :: solution vector
! - param     :: value of important parameter
! - step_num  :: number of iterations taken to compute the solution
! - max_steps :: maximum number of iterations allowed
!
   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: x
   REAL(KIND=8)                :: param
   INTEGER                     :: step_num, max_steps
   CHARACTER(*)                :: filenm

   LOGICAL :: existFlag

   INQUIRE( FILE = trim(p_in%restart_directory)//trim(filenm), EXIST = existFlag )
   IF (.NOT.existFlag) THEN

      WRITE(*,*)
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) '--> Writing restart file: '//trim(p_in%restart_directory)//trim(filenm)//' ...'

   ELSE

      ! may not be very portable
      !CALL RENAME(trim(p_in%restart_directory)//filenm(1:filenmLen), trim(p_in%restart_directory)//filenm(1:filenmLen)//'.bak')
      WRITE(*,*)
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) '--> Restart file '//trim(p_in%restart_directory)//trim(filenm)//' exists, OVERWRITING'

   ENDIF

   OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(filenm), FORM = 'UNFORMATTED' )


   WRITE(20) param
   WRITE(20) step_num, max_steps
   WRITE(20) velCmpnnts
   WRITE(20) np
   WRITE(20) np_L

!   WRITE(*,*) '    param      = ', param
!   WRITE(*,*) '    ite        = ', step_num, '/', max_steps
!   WRITE(*,*) '    velCmpnnts = ', velCmpnnts
!   WRITE(*,*) '    np         = ', np
!   WRITE(*,*) '    np_L       = ', np_L

!   CALL extract(x, u_save, p_save)
!   ! write fields ASCII
!   DO i = 1, np
!      WRITE(20, *) u_save(:,i)
!   END DO
!   DO i = 1, np_L
!      WRITE(20, *) p_save(i)
!   END DO
   !
   ! write fields BINARY
   WRITE(20) x

   CLOSE(20)


   WRITE(*,*) '    Done.'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! TEMPORARY
   !WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   !WRITE(*,*) '--> Writing matrix file: '//trim(p_in%restart_directory)//'matrix_'//filenm(1:filenmLen)//' ...'
   !OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//'matrix_'//filenm(1:filenmLen) )
   !do i = 1, size(Jacobian%e)
   !   write(20,*) Jacobian%i_mumps(i), Jacobian%j(i), Jacobian%e(i)
   !enddo
   !close(20)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE write_restart
!
! C compatible version follows
!
!! SUBROUTINE write_restart(x, param, step_num, max_steps, filenm, filenmLen) &
!!    BIND(C, NAME='write_restart')
!! !
!! ! Author: Jacopo Canton
!! ! E-mail: jcanton@mech.kth.se
!! ! Last revision: 9/9/2014
!! !
!! ! - x         :: solution vector
!! ! - param     :: value of important parameter
!! ! - step_num  :: number of iterations taken to compute the solution
!! ! - max_steps :: maximum number of iterations allowed
!! !
!!    USE ISO_C_BINDING
!! 
!!    IMPLICIT NONE
!!    ! input variables
!!    REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x
!!    REAL(KIND=C_DOUBLE), VALUE         :: param
!!    INTEGER(KIND=C_INT), VALUE         :: step_num, max_steps
!!    CHARACTER(KIND=C_CHAR)             :: filenm
!!    INTEGER(KIND=C_INT), VALUE         :: filenmLen
!!    !CHARACTER(KIND=C_CHAR),DIMENSION(filenmLen) :: filenm
!!    ! local variables
!! !   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_save
!! !   REAL(KIND=8), DIMENSION(np_L)          :: p_save
!! !   INTEGER :: i
!!    LOGICAL :: existFlag
!!    !CHARACTER(LEN=128) :: Ffilenm
!! 
!!    !Ffilenm = transfer(filenm(1:filenmLen), Ffilenm)
!! 
!!    INQUIRE( FILE = trim(p_in%restart_directory)//filenm(1:filenmLen), EXIST = existFlag )
!!    IF (.NOT.existFlag) THEN
!! 
!!       WRITE(*,*)
!!       WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!       WRITE(*,*) '--> Writing restart file: '//trim(p_in%restart_directory)//filenm(1:filenmLen)//' ...'
!!       !WRITE(*,*) '--> Writing restart file: '//trim(p_in%restart_directory)//trim(Ffilenm)//' ...'
!! 
!!    ELSE
!! 
!!       ! may not be very portable
!!       !CALL RENAME(trim(p_in%restart_directory)//filenm(1:filenmLen), trim(p_in%restart_directory)//filenm(1:filenmLen)//'.bak')
!!       WRITE(*,*)
!!       WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!       WRITE(*,*) '--> Restart file '//trim(p_in%restart_directory)//filenm(1:filenmLen)//' exists, overwriting'
!! 
!!    ENDIF
!! 
!!    OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//filenm(1:filenmLen), FORM = 'UNFORMATTED' )
!!    !OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(Ffilenm) )
!! 
!! 
!!    WRITE(20) param
!!    WRITE(20) step_num, max_steps
!!    WRITE(20) velCmpnnts
!!    WRITE(20) np
!!    WRITE(20) np_L
!! 
!! !   WRITE(*,*) '    param      = ', param
!! !   WRITE(*,*) '    ite        = ', step_num, '/', max_steps
!! !   WRITE(*,*) '    velCmpnnts = ', velCmpnnts
!! !   WRITE(*,*) '    np         = ', np
!! !   WRITE(*,*) '    np_L       = ', np_L
!! 
!! !   CALL extract(x, u_save, p_save)
!! !   ! write fields ASCII
!! !   DO i = 1, np
!! !      WRITE(20, *) u_save(:,i)
!! !   END DO
!! !   DO i = 1, np_L
!! !      WRITE(20, *) p_save(i)
!! !   END DO
!!    !
!!    ! write fields BINARY
!!    WRITE(20) x
!! 
!!    CLOSE(20)
!! 
!! 
!!    WRITE(*,*) '    Done.'
!! 
!!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    ! TEMPORARY
!!    !WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!    !WRITE(*,*) '--> Writing matrix file: '//trim(p_in%restart_directory)//'matrix_'//filenm(1:filenmLen)//' ...'
!!    !OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//'matrix_'//filenm(1:filenmLen) )
!!    !do i = 1, size(Jacobian%e)
!!    !   write(20,*) Jacobian%i_mumps(i), Jacobian%j(i), Jacobian%e(i)
!!    !enddo
!!    !close(20)
!!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! END SUBROUTINE write_restart

!------------------------------------------------------------------------------

SUBROUTINE read_restart(x, param, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/10/2014
!
! - x      :: solution vector
! - param  :: value of important parameter
!
   IMPLICIT NONE
   ! input variables
   CHARACTER(*) :: filenm
   ! output variables
   REAL(KIND=8), DIMENSION(Nx) :: x
   REAL(KIND=8)                :: param
   !CHARACTER(KIND=C_CHAR),DIMENSION(filenmLen) :: filenm
   ! local variables
   INTEGER :: step_num, max_steps
   INTEGER :: i, uc_in, np_in, np_L_in
#ifdef ASCIIRESTART
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_save
   REAL(KIND=8), DIMENSION(np_L)          :: p_save
#endif

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Reading restart file: '//trim(p_in%restart_directory)//trim(filenm)//' ...'

#ifdef ASCIIRESTART
   OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(filenm), FORM = 'FORMATTED' )

   READ(20, *) param
   READ(20, *) step_num, max_steps
   READ(20, *) uc_in   ! number of velocity components
   READ(20, *) np_in   ! number of P2 nodes
   READ(20, *) np_L_in ! number of P1 nodes
#else
   OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(filenm), FORM = 'UNFORMATTED' )

   READ(20) param
   READ(20) step_num, max_steps
   READ(20) uc_in   ! number of velocity components
   READ(20) np_in   ! number of P2 nodes
   READ(20) np_L_in ! number of P1 nodes
#endif

#if DEBUG > 0
   WRITE(*,*) '    param      = ', param
   WRITE(*,*) '    velCmpnnts = ', uc_in
   WRITE(*,*) '    np         = ', np_in
   WRITE(*,*) '    np_L       = ', np_L_in
#endif

   ! check dimension consistency
   ! uu
   IF ( velCmpnnts /= uc_in .OR. np /= np_in ) THEN
      WRITE(*,*) '    Error: wrong uu size'
      WRITE(*,*) '    SIZE(uu,1) = ', velCmpnnts, ', saved in this file: ', uc_in
      WRITE(*,*) '    SIZE(uu,2) = ', np,         ', saved in this file: ', np_in
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF
   ! pp
   IF ( np_L /= np_L_in ) THEN
      WRITE(*,*) '    Error: wrong pp size'
      WRITE(*,*) '    SIZE(pp) = ', np_L, ', saved in this file: ', np_L_in
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF

   ! read fields
#ifdef ASCIIRESTART
   DO i = 1, np
      READ(20, *) u_save(:,i)
   END DO
   DO i = 1, np_L
      READ(20, *) p_save(i)
   END DO
   CALL collect(u_save, p_save, x)
#else
   READ(20) x
#endif

   CLOSE(20)

   WRITE(*,*) '    Done.'

END SUBROUTINE read_restart
!
! C compatible version follows
!
!! SUBROUTINE read_restart(x, param, filenm, filenmLen) &
!!    BIND(C, NAME='read_restart')
!! !
!! ! Author: Jacopo Canton
!! ! E-mail: jcanton@mech.kth.se
!! ! Last revision: 9/9/2014
!! !
!! ! - x      :: solution vector
!! ! - param  :: value of important parameter
!! !
!!    USE ISO_C_BINDING
!! 
!!    IMPLICIT NONE
!!    ! input variables
!!    CHARACTER(KIND=C_CHAR)             :: filenm
!!    INTEGER(KIND=C_INT), VALUE         :: filenmLen
!!    ! output variables
!!    REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x
!!    REAL(KIND=C_DOUBLE)                :: param
!!    !CHARACTER(KIND=C_CHAR),DIMENSION(filenmLen) :: filenm
!!    ! local variables
!!    INTEGER                            :: step_num, max_steps
!!    INTEGER :: i, uc_in, np_in, np_L_in
!! !   CHARACTER(LEN=128) :: Ffilenm
!! #ifdef ASCIIRESTART
!!    REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_save
!!    REAL(KIND=8), DIMENSION(np_L)          :: p_save
!! #endif
!! 
!!    !Ffilenm = transfer(filenm(1:filenmLen), Ffilenm)
!! 
!!    WRITE(*,*)
!!    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!    WRITE(*,*) '--> Reading restart file: '//trim(p_in%restart_directory)//filenm(1:filenmLen)//' ...'
!!    !WRITE(*,*) '--> Reading restart file: '//trim(p_in%restart_directory)//trim(Ffilenm)//' ...'
!! 
!! #ifdef ASCIIRESTART
!!    OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//filenm(1:filenmLen), FORM = 'FORMATTED' )
!!    !OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(Ffilenm) )
!! 
!!    READ(20, *) param
!!    READ(20, *) step_num, max_steps
!!    READ(20, *) uc_in   ! number of velocity components
!!    READ(20, *) np_in   ! number of P2 nodes
!!    READ(20, *) np_L_in ! number of P1 nodes
!! #else
!!    OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//filenm(1:filenmLen), FORM = 'UNFORMATTED' )
!!    !OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(Ffilenm) )
!! 
!!    READ(20) param
!!    READ(20) step_num, max_steps
!!    READ(20) uc_in   ! number of velocity components
!!    READ(20) np_in   ! number of P2 nodes
!!    READ(20) np_L_in ! number of P1 nodes
!! #endif
!! 
!! #if DEBUG > 0
!!    WRITE(*,*) '    param      = ', param
!!    WRITE(*,*) '    velCmpnnts = ', uc_in
!!    WRITE(*,*) '    np         = ', np_in
!!    WRITE(*,*) '    np_L       = ', np_L_in
!! #endif
!! 
!!    ! check dimension consistency
!!    ! uu
!!    IF ( velCmpnnts /= uc_in .OR. np /= np_in ) THEN
!!       WRITE(*,*) '    Error: wrong uu size'
!!       WRITE(*,*) '    SIZE(uu,1) = ', velCmpnnts, ', saved in this file: ', uc_in
!!       WRITE(*,*) '    SIZE(uu,2) = ', np,         ', saved in this file: ', np_in
!!       WRITE(*,*) '    STOP.'
!!       CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!!    END IF
!!    ! pp
!!    IF ( np_L /= np_L_in ) THEN
!!       WRITE(*,*) '    Error: wrong pp size'
!!       WRITE(*,*) '    SIZE(pp) = ', np_L, ', saved in this file: ', np_L_in
!!       WRITE(*,*) '    STOP.'
!!       CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!!    END IF
!! 
!!    ! read fields
!! #ifdef ASCIIRESTART
!!    DO i = 1, np
!!       READ(20, *) u_save(:,i)
!!    END DO
!!    DO i = 1, np_L
!!       READ(20, *) p_save(i)
!!    END DO
!!    CALL collect(u_save, p_save, x)
!! #else
!!    READ(20) x
!! #endif
!! 
!!    CLOSE(20)
!! 
!!    WRITE(*,*) '    Done.'
!! 
!! END SUBROUTINE read_restart

!------------------------------------------------------------------------------

SUBROUTINE write_restart_bin(x, param, step_num, max_steps, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 15/1/2014
!
! - x         :: solution vector
! - param     :: value of important parameter
! - step_num  :: number of iterations taken to compute the solution
! - max_steps :: maximum number of iterations allowed
!

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:) :: x
   REAL(KIND=8)               :: param
   INTEGER                    :: step_num, max_steps
   CHARACTER(*)               :: filenm


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Writing binary restart file: '//trim(filenm)//' ...'


   OPEN( UNIT = 20, FILE = trim(filenm), FORM='UNFORMATTED' )

   WRITE(20) param
   WRITE(20) step_num, max_steps
   WRITE(20) SIZE(x)

   WRITE(*,*) '    param   = ', param
   WRITE(*,*) '    ite     = ', step_num, '/', max_steps
   WRITE(*,*) '    size(x) = ', SIZE(x)

   ! write field
   WRITE(20) x

   CLOSE(20)


   WRITE(*,*) '    Done.'

END SUBROUTINE write_restart_bin

!------------------------------------------------------------------------------

SUBROUTINE read_restart_bin(x, param, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/1/2014
!
! - x      :: solution vector
! - param  :: value of important parameter
!
   ! input variables
   CHARACTER(*)               :: filenm
   ! output variables
   REAL(KIND=8), DIMENSION(:) :: x
   REAL(KIND=8)               :: param
   ! local variables
   INTEGER :: sizeX


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Reading binary restart file: '//trim(filenm)//' ...'

   OPEN( UNIT = 20, FILE = trim(filenm), FORM='UNFORMATTED' )

   READ(20) param
   READ(20) ! jump this line
   READ(20) sizeX ! vector length

   WRITE(*,*) '    param      = ', param

   ! check dimension consistency
   !
   IF ( SIZE(x) /= sizeX ) THEN
      WRITE(*,*) '    Error: wrong x size'
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ', saved in this file: ', sizeX
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF

   ! read field
   READ(20) x

   CLOSE(20)

   WRITE(*,*) '    Done.'

END SUBROUTINE read_restart_bin

!------------------------------------------------------------------------------

SUBROUTINE write_cmplx_restart_bin(x, param, step_num, max_steps, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/03/2014
!
! - x         :: solution vector
! - param     :: value of important parameter
! - step_num  :: number of iterations taken to compute the solution
! - max_steps :: maximum number of iterations allowed
!

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:) :: x
   REAL(KIND=8)                  :: param
   INTEGER                       :: step_num, max_steps
   CHARACTER(*)                  :: filenm


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Writing binary restart file: '//trim(filenm)//' ...'


   OPEN( UNIT = 20, FILE = trim(filenm), FORM='UNFORMATTED' )

   WRITE(20) param
   WRITE(20) step_num, max_steps
   WRITE(20) SIZE(x)

   WRITE(*,*) '    param   = ', param
   WRITE(*,*) '    ite     = ', step_num, '/', max_steps
   WRITE(*,*) '    size(x) = ', SIZE(x)

   ! write field
   WRITE(20) x

   CLOSE(20)


   WRITE(*,*) '    Done.'

END SUBROUTINE write_cmplx_restart_bin

!------------------------------------------------------------------------------

SUBROUTINE read_cmplx_restart_bin(x, param, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/03/2014
!
! - x      :: solution vector
! - param  :: value of important parameter
!
   ! input variables
   CHARACTER(*)                  :: filenm
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:) :: x
   REAL(KIND=8)                  :: param
   ! local variables
   INTEGER :: sizeX


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Reading binary restart file: '//trim(filenm)//' ...'

   OPEN( UNIT = 20, FILE = trim(filenm), FORM='UNFORMATTED' )

   READ(20) param
   READ(20) ! jump this line
   READ(20) sizeX ! vector length

   WRITE(*,*) '    param      = ', param

   ! check dimension consistency
   !
   IF ( SIZE(x) /= sizeX ) THEN
      WRITE(*,*) '    Error: wrong x size'
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ', saved in this file: ', sizeX
      WRITE(*,*) '    STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END IF

   ! read field
   READ(20) x

   CLOSE(20)

   WRITE(*,*) '    Done.'

END SUBROUTINE read_cmplx_restart_bin

!------------------------------------------------------------------------------

SUBROUTINE write_QP_restart(x, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 2/10/2014
!
! - x         :: solution vector
!

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: x
   CHARACTER(*)                :: filenm
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_save
   REAL(KIND=8), DIMENSION(np_L)          :: p_save


   IF ( p_in%write_QP_restart_flag ) THEN
      ! WRITE QP RESTART
      WRITE(*,*)
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) '--> Writing QP restart file: '//trim(p_in%restart_directory)//trim(filenm)//' ...'
      
      WRITE(*,*) '    velCmpnnts = ', velCmpnnts
      WRITE(*,*) '    np         = ', np
      WRITE(*,*) '    np_L       = ', np_L
      

      OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//trim(filenm), STATUS='unknown', FORM='UNFORMATTED')
      
      CALL extract(x, u_save, p_save)

      WRITE(20) 0.d0, 1d-1, SIZE(uu,2), SIZE(pp)

      WRITE(20) u_save;  WRITE(20) u_save;  WRITE(20) u_save
      WRITE(20) p_save;  WRITE(20) p_save;  WRITE(20) p_save

      CLOSE(20)

      WRITE(*,*) '    Done.'

   ENDIF

   IF ( p_in%write_BVS_flag ) THEN
      ! WRITE BOUNDARY VALUES TO FILE
      !CALL write_BVS (8, u_save, rr, jjs, sides, filenm(1:filenmLen))
      !CALL write_BVS (9, u_save, rr, jjs, sides, filenm(1:filenmLen))
   END IF

END SUBROUTINE write_QP_restart
!
! C compatible version follows
!
!! SUBROUTINE write_QP_restart(x, filenm, filenmLen) &
!!    BIND(C, NAME='write_QP_restart')
!! !
!! ! Author: Jacopo Canton
!! ! E-mail: jcanton@mech.kth.se
!! ! Last revision: 15/7/2013
!! !
!! ! - x         :: solution vector
!! !
!!    USE ISO_C_BINDING
!! 
!!    IMPLICIT NONE
!!    ! input variables
!!    REAL(KIND=C_DOUBLE), DIMENSION(Nx) :: x
!!    CHARACTER(KIND=C_CHAR)             :: filenm
!!    INTEGER(KIND=C_INT), VALUE         :: filenmLen
!!    ! local variables
!!    REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_save
!!    REAL(KIND=8), DIMENSION(np_L)          :: p_save
!! 
!! 
!!    IF ( p_in%write_QP_restart_flag ) THEN
!!       ! WRITE QP RESTART
!!       WRITE(*,*)
!!       WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!!       WRITE(*,*) '--> Writing QP restart file: '//trim(p_in%restart_directory)//filenm(1:filenmLen)//' ...'
!!       
!!       WRITE(*,*) '    velCmpnnts = ', velCmpnnts
!!       WRITE(*,*) '    np         = ', np
!!       WRITE(*,*) '    np_L       = ', np_L
!!       
!! 
!!       OPEN( UNIT = 20, FILE = trim(p_in%restart_directory)//filenm(1:filenmLen), STATUS='unknown', FORM='UNFORMATTED')
!!       
!!       CALL extract(x, u_save, p_save)
!! 
!!       WRITE(20) 0.d0, 1d-1, SIZE(uu,2), SIZE(pp)
!! 
!!       WRITE(20) u_save;  WRITE(20) u_save;  WRITE(20) u_save
!!       WRITE(20) p_save;  WRITE(20) p_save;  WRITE(20) p_save
!! 
!!       CLOSE(20)
!! 
!!       WRITE(*,*) '    Done.'
!! 
!!    ENDIF
!! 
!!    IF ( p_in%write_BVS_flag ) THEN
!!       ! WRITE BOUNDARY VALUES TO FILE
!!       !CALL write_BVS (8, u_save, rr, jjs, sides, filenm(1:filenmLen))
!!       !CALL write_BVS (9, u_save, rr, jjs, sides, filenm(1:filenmLen))
!!    END IF
!! 
!! END SUBROUTINE write_QP_restart

!==============================================================================

END MODULE restart_io
