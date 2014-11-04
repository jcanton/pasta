MODULE vorticity_stream

   USE global_variables
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE Dirichlet_Neumann ! for Dirichlet_nodes_gen subroutine
   USE start_sparse_kit  ! for start_matrix_2d_p2
   USE qs_sp
   USE qs_sp_M
   USE qv_sp
   USE par_solve_mumps

   IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE  compute_vorticity_stream (jj, jjs, js, uu, rr, sides, Axis, Dir_psi,  zz, psi)

!  Compute the vorticity field  zz  and Stokes stream function  psi
!  corresponding to the 2D solenoidal velocity field  uu

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jjs
   INTEGER,      DIMENSION(:),   INTENT(IN) :: js
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:),   INTENT(IN) :: sides
   LOGICAL,      DIMENSION(:),   INTENT(IN) :: Axis
   LOGICAL,      DIMENSION(:),   INTENT(IN) :: Dir_psi

   REAL(KIND=8), DIMENSION(:)               :: zz, psi

   LOGICAL, SAVE :: first_time = .TRUE.

   TYPE(CSR_MUMPS_matrix), SAVE :: vortMatr, psiMatr

   INTEGER,      DIMENSION(:), POINTER,     SAVE :: js_psi_D, js_Axis
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: as_psi_D


   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to compute_vorticity_stream'

!------------------------------------------------------------------------------
!-------------MATRICES ALLOCATION AND SYMBOLIC FACTORIZATION-------------------

   IF (first_time) THEN

      CALL Dirichlet_nodes_gen (jjs, sides, Axis,  js_Axis)

      CALL Dirichlet_nodes_gen (jjs, sides, Dir_psi .AND. .NOT.Axis,  js_psi_D)
      ALLOCATE (as_psi_D(SIZE(js_psi_D)))


      WRITE(*,*)
      WRITE(*,*) '    Structuring of the matrix for the vorticity problem'

      CALL start_matrix_2d_p2 (SIZE(uu,2), jj, js,  vortMatr)

!      CALL symbolic_factorization (vortMatr, 1, 5)
      CALL par_mumps_master (INITIALIZATION, 10, vortMatr, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   10, vortMatr, 0)

!      WRITE (*,*) '    Symbolic factorization of vortMatr matrix for vorticity computation'

      ALLOCATE (vortMatr%e(SIZE(vortMatr%j)))

      CALL qs_0y0_sp_M (1.d0,     vortMatr)
      CALL Dirichlet_M (js_Axis,  vortMatr)

!      CALL numerical_factorization (vortMatr, 5)
      CALL par_mumps_master (NUMER_FACTOR,   10, vortMatr, 0)

!      WRITE (*,*) '    Numerical factorization of problem  vortMatr zz = k.rot u '


      WRITE(*,*)
      WRITE(*,*) '    Structuring of the matrix for the Stokes stream function problem'

!      psiMatr%i => vortMatr%i
!      psiMatr%j => vortMatr%j
      ALLOCATE (psiMatr%i      (SIZE(vortMatr%i))      ); psiMatr%i       = vortMatr%i
      ALLOCATE (psiMatr%i_mumps(SIZE(vortMatr%i_mumps))); psiMatr%i_mumps = vortMatr%i_mumps
      ALLOCATE (psiMatr%j      (SIZE(vortMatr%j))      ); psiMatr%j       = vortMatr%j
      ALLOCATE (psiMatr%e      (SIZE(vortMatr%e))      ); psiMatr%e       = 0d0

!      CALL symbolic_factorization (psiMatr, 1, 6)
      CALL par_mumps_master (INITIALIZATION, 11, psiMatr, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   11, psiMatr, 0)

!      WRITE (*,*) ' Symbolic factorization of matrix of psiMatr = (Dw).R D + 1/R '

      CALL qs_1y1_sp_M (1.d0,  psiMatr, 1d0) ! + SINGULAR TERM

      CALL Dirichlet_M (js_Axis,   psiMatr)
      CALL Dirichlet_M (js_psi_D,  psiMatr)

!      CALL numerical_factorization (psiMatr, 6)
      CALL par_mumps_master (NUMER_FACTOR,   11, psiMatr, 0)

!      WRITE (*,*) ' Numerical factorization of matrix '
!      WRITE (*,*) ' of problem  [(Dw).R D + 1/R] psi = R k.Rot u '


      first_time = .FALSE.

   ENDIF


!------------------------------------------------------------------------------
!-------------VORTICITY COMPUTATION--------------------------------------------

   ! right hand side for the vorticity equation

   CALL qs_0y1_sp_c (uu,  zz)  !  zz <--- (w, k.Rot u)

   CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  zz)

!   CALL direct_solution (zz, 5)
   CALL par_mumps_master (DIRECT_SOLUTION, 10, vortMatr, 0, zz)


!   WRITE (*,*) ' Solution of problem  vortMatr zz = k.Rot u '
   WRITE (*,*) '   Vorticity field computed'


!------------------------------------------------------------------------------
!-------------STOKES STREAM FUNCTION COMPUTATION-------------------------------

   ! right hand side for the Stokes stream function elliptic equation

   CALL qs_0y1_sp_c (uu,  psi)  !  psi <--- (w, y k.Rot u)

   as_psi_D = stream_boundary_values (js_psi_D, rr)

   CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  psi, .true.)
   CALL Dirichlet (js_psi_D,                    as_psi_D,  psi, .true.)

!   CALL direct_solution (psi, 6)
   CALL par_mumps_master (DIRECT_SOLUTION, 11, psiMatr, 0, psi)

!   WRITE (*,*) ' Solution of problem  [(Dw).R D + 1/R] psi = R k.Rot u '
   WRITE (*,*) '    Stokes stream function computed'


!   CALL par_mumps_master (DEALLOCATION, 10, vortMatr, 0)
!   CALL par_mumps_master (DEALLOCATION, 11, psiMatr, 0)


END SUBROUTINE compute_vorticity_stream

!-----------------------------------------------------------------------------

SUBROUTINE  compute_axial_plane_vorticity (jj, jjs, js, ww, Axis,  zz_R, zz_z)

!  Compute the vorticity components  zz_R  and  zz_z  of an
!  axisymmetric swrirling velocity component ww

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jjs
   INTEGER,      DIMENSION(:),   INTENT(IN) :: js
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww
   LOGICAL,      DIMENSION(:),   INTENT(IN) :: Axis

   REAL(KIND=8), DIMENSION(:)               :: zz_R, zz_z

   REAL(KIND=8), DIMENSION(2, SIZE(ww)) :: uu

   LOGICAL, SAVE :: first_time = .TRUE.

   TYPE(CSR_MUMPS_matrix), SAVE :: vortMatr_R, vortMatr_z

   INTEGER, DIMENSION(:), POINTER, SAVE :: js_Axis

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to compute_axial_plane_vorticity'

!------------------------------------------------------------------------------
!-------------MATRICES ALLOCATION AND SYMBOLIC FACTORIZATION-------------------

   IF (first_time) THEN

      CALL Dirichlet_nodes_gen (jjs, sides, Axis,  js_Axis)

      WRITE(*,*)
      WRITE(*,*) '    Structuring of the vortMatr_R matrix for zz_R'

      CALL start_matrix_2d_p2 (SIZE(ww), jj, js,  vortMatr_R)

!      CALL symbolic_factorization (vortMatr_R, 1, 7)
      CALL par_mumps_master (INITIALIZATION, 12, vortMatr_R, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   12, vortMatr_R, 0)

!      WRITE (*,*) ' Symbolic factorization of vortMatr_R matrix for zz_R'

      ALLOCATE (vortMatr_R%e(SIZE(vortMatr_R%j)))

      CALL qs_00_sp_M  (1.d0,     vortMatr_R, .true.)
      CALL Dirichlet_M (js_Axis,  vortMatr_R, .true.)

!      CALL numerical_factorization (vortMatr_R, 7)
      CALL par_mumps_master (NUMER_FACTOR,   12, vortMatr_R, 0)

!      WRITE (*,*) ' Numerical factorization for problem  vortMatr_R zz_R = - dw/dz '


 !     WRITE (*,*) ' Structuring of vortMatr_z matrix for zz_z'

!      vortMatr_z%i => vortMatr_R%i
!      vortMatr_z%j => vortMatr_R%j
      ALLOCATE (vortMatr_z%i      (SIZE(vortMatr_R%i))      ); vortMatr_z%i       = vortMatr_R%i
      ALLOCATE (vortMatr_z%i_mumps(SIZE(vortMatr_R%i_mumps))); vortMatr_z%i_mumps = vortMatr_R%i_mumps
      ALLOCATE (vortMatr_z%j      (SIZE(vortMatr_R%j))      ); vortMatr_z%j       = vortMatr_R%j
      ALLOCATE (vortMatr_z%e      (SIZE(vortMatr_R%e))      ); vortMatr_z%e       = 0d0

!      CALL symbolic_factorization (vortMatr_z, 1, 8)
      CALL par_mumps_master (INITIALIZATION, 13, vortMatr_z, 0)
      CALL par_mumps_master (SYMBO_FACTOR,   13, vortMatr_z, 0)

!      WRITE (*,*) ' Symbolic factorization of vortMatr_z matrix for zz_z'

      CALL qs_00_sp_M  (1d0,  vortMatr_z, .true.)

!      CALL numerical_factorization (vortMatr_z, 8)
      CALL par_mumps_master (NUMER_FACTOR,   13, vortMatr_z, 0)

!      WRITE (*,*) ' Numerical factorization for problem  vortMatr_z zz_z = dw/dR '


      first_time = .FALSE.

   ENDIF


   ! right hand sides of the equations for the
   ! vorticity components in the axial plane

   CALL qv_01_sp (mm, jj, ww,  uu)

   zz_R = - uu(1,:)
   zz_z =   uu(2,:)


!------------------------------------------------------------------------------
!-------------RADIAL VORTICITY COMPUTATION-------------------------------------

   CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  zz_R, .true.)

!   CALL direct_solution (zz_R, 7)
   CALL par_mumps_master (DIRECT_SOLUTION, 12, vortMatr_R, 0, zz_R)

!   WRITE (*,*) ' Solution of problem  vortMatr_R zz_R = - dw/dz '
   WRITE(*,*) '    Radial vorticity component computed'


!------------------------------------------------------------------------------
!-------------AXIAL VORTICITY COMPUTATION--------------------------------------

!   CALL direct_solution (zz_z, 8)
   CALL par_mumps_master (DIRECT_SOLUTION, 13, vortMatr_z, 0, zz_z)

!   WRITE (*,*) ' Solution of problem  vortMatr_z zz_z = dw/dR '
   WRITE(*,*) '    Axial vorticity component computed'


END SUBROUTINE  compute_axial_plane_vorticity

!------------------------------------------------------------------------------

FUNCTION  stream_boundary_values (js, rr) RESULT(psis)

!  This program defines boundary values for Stokes stream function

!  Here, the velocity boundary values are assumed not to depend on time

   IMPLICIT NONE

   INTEGER,      DIMENSION(:),   INTENT(IN) :: js
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(SIZE(js)) :: psis


   psis = 0d0


   ! Inflow boundary

!!   WHERE (rr(2,js)  >  4*a - eps)
!!
!!      ! per il segno meno, non dimenticare che l'integrazione
!!      ! lungo il contorno deve essere fatta con la variabile
!!      ! che aumenta percorrendo il contorno nel verso positivo
!!
!!      psis = - uR_IN_max * 4*a  &
!!             * (   (1/3.0d0) * (rr(1,js)/a)**3  &
!!                 + (3/2.0d0) * (rr(1,js)/a)**2  &
!!                 +     2     *  rr(1,js)/a      &
!!                 +  2/3.0d0  )
!!
!!   END WHERE
!!
!!
!!   ! Superior wall of the inlet section and of the tube
!!
!!   WHERE ( (rr(1,js) > -a - eps  .AND.  rr(2,js) > a - eps)  .OR.  &
!!
!!           (rr(1,js) >  H - eps  .AND.  rr(2,js) > b - eps) )
!!
!!      psis = uR_IN_max  * 4 * (2/3.0d0) * a**2 / rr(2,js)
!!
!!    ! psis = uz_OUT_max * b**2 / rr(2,js)  !  alternativa equivalente
!!
!!   END WHERE
!!
!!
!!   ! Outflow boundary
!!
!!   WHERE ( rr(1,js) > H + a - eps  .AND.  &
!!           rr(2,js) >= 0  .AND.  rr(2,js) <= b )
!!
!!      psis = uz_OUT_max * (rr(2,js)/2) * (1 - (rr(2,js)/b)**2/2)
!!
!!
!!   END WHERE


END FUNCTION  stream_boundary_values

!=============================================================================

END MODULE vorticity_stream
