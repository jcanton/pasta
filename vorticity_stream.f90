MODULE vorticity_stream

   USE global_variables
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE Dirichlet_Neumann ! for Dirichlet_nodes_gen subroutine
   USE start_sparse_kit  ! for start_matrix_2d_p2
   USE qs_sp
   USE qs_sp_M
   USE par_solve_mumps

   IMPLICIT NONE

CONTAINS

! added to avoid errors when compiling with intel compilers
SUBROUTINE dummySubroutine()
   IMPLICIT NONE
END SUBROUTINE dummySubroutine

!-----------------------------------------------------------------------------

!!  SUBROUTINE  compute_vorticity_stream (jj, jjs, js, uu, rr, t, sides, Axis, Dir_psi,  zz, psi)
!!  
!!  !  Compute the vorticity field  zz  and Stokes stream function  psi
!!  !  corresponding to the 2D solenoidal velocity field  uu
!!  
!!     IMPLICIT NONE
!!  
!!     INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jjs
!!     INTEGER,      DIMENSION(:),   INTENT(IN) :: js
!!     REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu
!!     REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
!!     REAL(KIND=8),                 INTENT(IN) :: t ! time
!!     INTEGER,      DIMENSION(:),   INTENT(IN) :: sides
!!     LOGICAL,      DIMENSION(:),   INTENT(IN) :: Axis
!!     LOGICAL,      DIMENSION(:),   INTENT(IN) :: Dir_psi
!!  
!!     REAL(KIND=8), DIMENSION(:),   INTENT(OUT):: zz, psi
!!  
!!  
!!     LOGICAL, SAVE :: first_time = .TRUE.
!!  
!!     TYPE(CSR_matrix), SAVE :: MM, KK
!!  
!!     INTEGER,      DIMENSION(:), POINTER,     SAVE :: js_psi_D, js_Axis
!!     REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: as_psi_D
!!  
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------MATRICES ALLOCATION AND SYMBOLIC FACTORIZATION-------------------
!!  
!!     IF (first_time) THEN
!!     
!!        CALL Dirichlet_nodes_gen (jjs, sides, Axis,  js_Axis)
!!  
!!        CALL Dirichlet_nodes_gen (jjs, sides, Dir_psi .AND. .NOT.Axis,  js_psi_D)
!!        ALLOCATE (as_psi_D(SIZE(js_psi_D)))
!!    
!!  
!!        WRITE (*,*) ' Structuring of the matrix for the vorticity problem'
!!  
!!        CALL start_matrix_2d_p2 (SIZE(uu,2), jj, js,  MM)
!!  
!!        CALL symbolic_factorization (MM, 1, 5)
!!  
!!        WRITE (*,*) ' Symbolic factorization of MM matrix for vorticity computation'
!!  
!!        ALLOCATE (MM%e(SIZE(MM%j)))
!!  
!!        CALL qs_00_sp_M  (1.d0,     MM, .true.) 
!!        CALL Dirichlet_M (js_Axis,  MM, .true.) 
!!  
!!        CALL numerical_factorization (MM, 5)
!!  
!!        WRITE (*,*) ' Numerical factorization of problem  MM zz = k.rot u '
!!  
!!  
!!        WRITE (*,*) ' Structuring of the matrix for the Stokes stream function problem'
!!        
!!        KK%i => MM%i;  KK%j => MM%j;  ALLOCATE (KK%e(SIZE(MM%e)))
!!  
!!        CALL symbolic_factorization (KK, 1, 6)
!!  
!!        WRITE (*,*) ' Symbolic factorization of matrix of KK = (Dw).R D + 1/R '
!!  
!!        CALL qs_1y1_sp_M (1.d0,  KK, 1.0d0, .true.) ! + SINGULAR TERM
!!  
!!        CALL Dirichlet_M (js_Axis,   KK, .true.)
!!        CALL Dirichlet_M (js_psi_D,  KK, .true.)
!!  
!!        CALL numerical_factorization (KK, 6) 
!!  
!!        WRITE (*,*) ' Numerical factorization of matrix '
!!        WRITE (*,*) ' of problem  [(Dw).R D + 1/R] psi = R k.Rot u '
!!  
!!  
!!        first_time = .FALSE.
!!  
!!     ENDIF
!!  
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------VORTICITY COMPUTATION--------------------------------------------
!!  
!!     ! right hand side for the vorticity equation 
!!     
!!     CALL qs_01_sp_c (uu,  zz)  !  zz <--- (w, k.Rot u)
!!  
!!     CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  zz, .true.)
!!    
!!     CALL direct_solution (zz, 5)  
!!     
!!     WRITE (*,*) ' Solution of problem  MM zz = k.Rot u '  
!!     WRITE (*,*) 'Vorticity field computed'
!!  
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------STOKES STREAM FUNCTION COMPUTATION-------------------------------
!!  
!!     ! right hand side for the Stokes stream function elliptic equation
!!  
!!     CALL qs_0y1_sp_c (uu,  psi)  !  psi <--- (w, y k.Rot u)
!!  
!!     as_psi_D = Stokes_stream_boundary_values (js_psi_D, rr, t)
!!  
!!     CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  psi, .true.)
!!     CALL Dirichlet (js_psi_D,                    as_psi_D,  psi, .true.)
!!  
!!     CALL direct_solution (psi, 6)
!!     
!!     WRITE (*,*) ' Solution of problem  [(Dw).R D + 1/R] psi = R k.Rot u ' 
!!     WRITE (*,*) 'Stokes stream function computed'
!!  
!!  
!!  END SUBROUTINE compute_vorticity_stream
!!  
!!  !-----------------------------------------------------------------------------
!!  
!!  SUBROUTINE  axial_plane_vorticity (jj, jjs, js, ww, Axis,  zz_R, zz_z)
!!  
!!  !  Compute the vorticity components  zz_R  and  zz_z  of an
!!  !  axisymmetric swrirling velocity component ww 
!!  
!!     IMPLICIT NONE
!!  
!!     INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jjs 
!!     INTEGER,      DIMENSION(:),   INTENT(IN) :: js
!!     REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww
!!     LOGICAL,      DIMENSION(:),   INTENT(IN) :: Axis
!!    
!!     REAL(KIND=8), DIMENSION(:),   INTENT(OUT):: zz_R, zz_z
!!     
!!     REAL(KIND=8), DIMENSION(2, SIZE(ww)) :: uu
!!  
!!     LOGICAL, SAVE :: first_time = .TRUE.
!!  
!!     TYPE(CSR_matrix), SAVE :: MM, NN
!!  
!!     INTEGER, DIMENSION(:), POINTER, SAVE :: js_Axis
!!   
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------MATRICES ALLOCATION AND SYMBOLIC FACTORIZATION-------------------
!!  
!!     IF (first_time) THEN
!!     
!!        CALL Dirichlet_nodes_gen (jjs, sides, Axis,  js_Axis)
!!  
!!        WRITE (*,*) ' Structuring of the MM matrix for zz_R'
!!  
!!        CALL start_matrix_2d_p2 (SIZE(ww), jj, js,  MM)
!!  
!!        CALL symbolic_factorization (MM, 1, 7)
!!  
!!        WRITE (*,*) ' Symbolic factorization of MM matrix for zz_R'
!!  
!!        ALLOCATE (MM%e(SIZE(MM%j)))
!!  
!!        CALL qs_00_sp_M  (1.d0,     MM, .true.)                                                
!!        CALL Dirichlet_M (js_Axis,  MM, .true.)
!!  
!!        CALL numerical_factorization (MM, 7)
!!  
!!        WRITE (*,*) ' Numerical factorization for problem  MM zz_R = - dw/dz '
!!  
!!  
!!        WRITE (*,*) ' Structuring of NN matrix for zz_z'
!!        
!!        NN%i => MM%i;  NN%j => MM%j;  ALLOCATE (NN%e(SIZE(MM%e)))
!!  
!!        CALL symbolic_factorization (NN, 1, 8)
!!  
!!        WRITE (*,*) ' Symbolic factorization of NN matrix for zz_z'
!!  
!!        CALL qs_00_sp_M  (1.d0,  NN, .true.) 
!!  
!!        CALL numerical_factorization (NN, 8) 
!!  
!!        WRITE (*,*) ' Numerical factorization for problem  NN zz_z = dw/dR '
!!  
!!  
!!        first_time = .FALSE.
!!  
!!     ENDIF
!!  
!!  
!!     ! right hand sides of the equations for the 
!!     ! vorticity components in the axial plane  
!!  
!!     CALL qv_01_sp (ww,  uu)  
!!  
!!     zz_R = - uu(1,:) 
!!     zz_z =   uu(2,:) 
!!  
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------RADIAL VORTICITY COMPUTATION-------------------------------------
!!  
!!     CALL Dirichlet (js_Axis, SPREAD(0.d0,1,SIZE(js_Axis)),  zz_R, .true.)
!!  
!!     CALL direct_solution (zz_R, 7)  
!!     
!!     WRITE (*,*) ' Solution of problem  MM zz_R = - dw/dz '  
!!     WRITE (*,*) 'Radial vorticity component computed'
!!  
!!  
!!  !------------------------------------------------------------------------------
!!  !-------------AXIAL VORTICITY COMPUTATION--------------------------------------
!!    
!!     CALL direct_solution (zz_z, 8)  
!!     
!!     WRITE (*,*) ' Solution of problem  NN zz_z = dw/dR '
!!     WRITE (*,*) 'Axial vorticity component computed'
!!  
!!  
!!  END SUBROUTINE  axial_plane_vorticity
!!  
!!  !------------------------------------------------------------------------------
!!  
!!  FUNCTION  Stokes_stream_boundary_values (js, rr, t) RESULT(psis)
!!  
!!  !  This program defines boundary values for Stokes stream function 
!!  !  at time  t (optional)
!!  
!!  !  Here, the velocity boundary values are assumed not to depend on time  
!!  
!!     IMPLICIT NONE
!!    
!!     INTEGER,      DIMENSION(:),   INTENT(IN) :: js
!!     REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
!!     REAL(KIND=8), OPTIONAL,       INTENT(IN) :: t
!!     
!!     REAL(KIND=8), DIMENSION(SIZE(js)) :: psis
!!  
!!     
!!     psis = 0      
!!  
!!  
!!     ! Inflow boundary 
!!  
!!     WHERE (rr(2,js)  >  4*a - eps)
!!     
!!        ! per il segno meno, non dimenticare che l'integrazione
!!        ! lungo il contorno deve essere fatta con la variabile
!!        ! che aumenta percorrendo il contorno nel verso positivo 
!!      
!!        psis = - uR_IN_max * 4*a  & 
!!               * (   (1/3.0d0) * (rr(1,js)/a)**3  &
!!                   + (3/2.0d0) * (rr(1,js)/a)**2  &
!!                   +     2     *  rr(1,js)/a      &
!!                   +  2/3.0d0  ) 
!!     
!!     END WHERE
!!  
!!  
!!     ! Superior wall of the inlet section and of the tube 
!!  
!!     WHERE ( (rr(1,js) > -a - eps  .AND.  rr(2,js) > a - eps)  .OR.  &
!!             
!!             (rr(1,js) >  H - eps  .AND.  rr(2,js) > b - eps) )
!!      
!!        psis = uR_IN_max  * 4 * (2/3.0d0) * a**2 / rr(2,js)
!!        
!!      ! psis = uz_OUT_max * b**2 / rr(2,js)  !  alternativa equivalente  
!!  
!!     END WHERE  
!!  
!!  
!!     ! Outflow boundary 
!!  
!!     WHERE ( rr(1,js) > H + a - eps  .AND.  &
!!             rr(2,js) >= 0  .AND.  rr(2,js) <= b )
!!  
!!        psis = uz_OUT_max * (rr(2,js)/2) * (1 - (rr(2,js)/b)**2/2)           
!!              
!!  
!!     END WHERE
!!  
!!     
!!     IF (PRESENT(t)) psis = psis * time_dep(t)
!!  
!!     
!!  END FUNCTION  Stokes_stream_boundary_values

!=============================================================================

END MODULE vorticity_stream
