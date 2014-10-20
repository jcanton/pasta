MODULE transient_growth
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 4/5/2014
!
   USE dynamic_structures
   USE sparse_matrix_profiles
   USE sparse_matrix_operations
   USE global_variables
   USE miscellaneous_subroutines
   USE prep_mesh_p1p2_sp ! for some global variables as jj
   USE Dirichlet_Neumann ! for Dirichlet_nodes_gen subroutine
   USE qc_sp_M
   USE qv_sp
   USE qs_l_m
   USE qs_l_sp
   USE par_solve_mumps
   USE vtk_plot
   USE eigensolve ! needed if one wants to use an eigenvector as initial guess
   USE dns_algorithms

   IMPLICIT NONE
   ! variables "global" to this module
   TYPE(CSR_MUMPS_Complex_Matrix)               :: Wd, Zd, Wa, Za, W0, Z0, MassV
   COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: uu_tg, u0_tg
   COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:)   :: pp_tg, p0_tg
   COMPLEX(KIND=8), DIMENSION(:), POINTER       :: xx_tg, x0_tg
   LOGICAL, ALLOCATABLE, DIMENSION(:,:)         :: Dir_tg
   TYPE(dyn_int_line), DIMENSION(velCmpnnts)    :: js_D_tg
   LOGICAL                                      :: DESINGULARIZE_tg
   TYPE(dyn_real_line), DIMENSION(velCmpnnts)   :: zero_bvs_D_tg
   REAL(KIND=8) :: dtE, dtH

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE compute_transientGrowth(x_vec, Lns, filenm)

   IMPLICIT NONE

   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x_vec  ! computed base flow
   TYPE(CSR_MUMPS_Matrix)                 :: Lns    ! Linearized Navier--Stokes operator (Jacobian matrix)
   CHARACTER(*),               INTENT(IN) :: filenm
   ! "output" variables
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: singularValue  ! these two have a dimension more than necessary
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: singularVector ! because the algorithm for the SVD is "general"
   COMPLEX(KIND=8), DIMENSION(np)               :: Ek0, Ek        ! kinetic energy fields
   COMPLEX(KIND=8)                              :: Ek0_s, Ek_s    ! kinetic energy
   ! local variables
   INTEGER      :: nsv
   COMPLEX(KIND=8) :: sigma
   INTEGER      :: nSteps
   INTEGER      :: i, k, ii
   LOGICAL      :: existFlag
   COMPLEX(KIND=8) :: resLinfty
   INTEGER      :: iteNum
   COMPLEX(KIND=8), DIMENSION(velCmpnnts,np) :: uInitGuess

   CHARACTER(LEN=128) :: restart_name
   REAL(KIND=8) :: dummy

   REAL(KIND=8) :: dataIdentifier = 313d0 ! Donald Duck's plate number!
   INTEGER      :: tausNumber = 1
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: taus
   CHARACTER(LEN=128) :: tauName ! used to insert tau in file names
   INTEGER            :: tauNameTemp

   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmpVector

!----------------------------------------------------
   WRITE(*,*) ''
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> CALL to compute_transientGrowth'
   WRITE(*,*) ''
!----------------------------------------------------

!!!!----------------------------------------------------
!!!!----------------------------------------------------
!!!!----------------------------------------------------
!!!   Ek0_s = 0*1
!!!   write(*,*) Ek0_s
!!!   Ek0_s = CMPLX(0d0,0d0,KIND=8)
!!!   write(*,*) Ek0_s
!!!   Ek0_s = Ek0_s + 1
!!!   write(*,*) Ek0_s
!!!   Ek0_s = CMPLX(0d0,0d0,KIND=8) + CMPLX(0d0,1d0,KIND=8) * 0.2d0
!!!   write(*,*) Ek0_s
!!!   stop
!!!   write(*,*) Ek0_s
!!!   write(*,*) Ek0_s
!!!
!!!!----------------------------------------------------
!!!!----------------------------------------------------
!!!!----------------------------------------------------
   ALLOCATE (uu_tg(velCmpnnts, np), u0_tg(velCmpnnts, np))
   ALLOCATE (pp_tg(np_L),           p0_tg(np_L))
   ALLOCATE (xx_tg(Nx),             x0_tg(Nx))

!--------------------------------------------------------
! DEFINE HERE SOME PARAMETERS WHICH SHOULD NOT BE CHANGED

   nsv   = 1   ! we only want the largest singular value
   sigma = CMPLX(0d0,0d0,KIND=8) ! and we don't use any shift

!--------------------------------------------------
! PREPARE MATRICES FOR THE TRANSIENT GROWTH PROBLEM

   ! (1)
   ! prepare boundary conditions
   !
!   WRITE(*,*) '*check*'
!   WRITE(*,*) '    beta = ', beta
!   WRITE(*,*) '    number_of_sides = ', number_of_sides

   ALLOCATE ( Dir_tg(velCmpnnts, number_of_sides) )
   IF ( p_in%tranGrowth_BC == 1 ) THEN
      ! homogeneous Dirichlet on every border
      WRITE(*,*) '    boundary conditions: zero velocity on every border'
      Dir_tg = .TRUE.
      DESINGULARIZE_tg = .TRUE.
   ELSEIF ( p_in%tranGrowth_BC == 2 ) THEN
      ! same BCs as for the base flow but homogeneous
      ! Dirichlet where the base flow has nonhomogeneous Dirichlet
      WRITE(*,*) '    boundary conditions: same as base flow'
      Dir_tg = Dir
      DESINGULARIZE_tg = DESINGULARIZE
   ELSE
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong parameter:              ***'
      WRITE(*,*) '*** p_in % tranGrowth_BC          ***'
      WRITE(*,*) '*** set to: ', p_in%tranGrowth_BC
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

   DO k = 1, velCmpnnts
      CALL Dirichlet_nodes_gen (jjs, sides, Dir_tg(k,:) .AND. .NOT.Axis,  js_D_tg(k)%DIL)
      ALLOCATE (zero_bvs_D_tg(k)%DRL(SIZE(js_D_tg(k)%DIL)))
      zero_bvs_D_tg(k)%DRL  = 0d0
   ENDDO

   ! (2)
   ! create the Mass matrix
   !
   ALLOCATE( Mass%i      (SIZE(Lns%i))       ); Mass%i       = Lns%i
   ALLOCATE( Mass%i_mumps(SIZE(Lns%i_mumps)) ); Mass%i_mumps = Lns%i_mumps
   ALLOCATE( Mass%j      (SIZE(Lns%j))       ); Mass%j       = Lns%j
   ALLOCATE( Mass%e      (SIZE(Lns%e))       ); Mass%e       = 0d0

   CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Mass)

   ! (2b)
   ! create the MassV matrix (non-singolar mass only for the velocity)
   !
   CALL zEssM (CMPLX(Mass%e,0d0,KIND=8), Mass%j,  Mass%i,  velCmpnnts*np, &
                                MassV%e, MassV%j, MassV%i, MassV%i_mumps)

   CALL par_mumps_master (INITIALIZATION, 9, MassV, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   9, MassV, 0)
   CALL par_mumps_master (NUMER_FACTOR,   9, MassV, 0)

   ! (2c)
   ! enforce boundary conditions on Mass
   !
   CALL Dirichlet_rc_M (np, js_Axis, js_D_tg, 0d0,  Mass)

   ALLOCATE( Mass_cmplx%i      (SIZE(Lns%i))       ); Mass_cmplx%i       = Lns%i
   ALLOCATE( Mass_cmplx%i_mumps(SIZE(Lns%i_mumps)) ); Mass_cmplx%i_mumps = Lns%i_mumps
   ALLOCATE( Mass_cmplx%j      (SIZE(Lns%j))       ); Mass_cmplx%j       = Lns%j
   ALLOCATE( Mass_cmplx%e      (SIZE(Lns%e))       ); Mass_cmplx%e       = CMPLX(Mass%e,0d0,KIND=8)

   DEALLOCATE( Mass%i, Mass%i_mumps, Mass%j, Mass%e )

   ! (3)
   ! create the matrices for the time stepper. Store Wd in
   ! position 5 of the MUMPS array id and Wa in position 8
   ! Prepare also Zd and Za, matrices for the RHSs
   !
   ALLOCATE( Lns_cmplx%i      (SIZE(Lns%i))       ); Lns_cmplx%i       = Lns%i
   ALLOCATE( Lns_cmplx%i_mumps(SIZE(Lns%i_mumps)) ); Lns_cmplx%i_mumps = Lns%i_mumps
   ALLOCATE( Lns_cmplx%j      (SIZE(Lns%j))       ); Lns_cmplx%j       = Lns%j
   ALLOCATE( Lns_cmplx%e      (SIZE(Lns%e))       ); Lns_cmplx%e       = CMPLX(Lns%e,0d0,KIND=8)

   ALLOCATE( Wd%i      (SIZE(Lns%i))       ); Wd%i       = Lns%i
   ALLOCATE( Wd%i_mumps(SIZE(Lns%i_mumps)) ); Wd%i_mumps = Lns%i_mumps
   ALLOCATE( Wd%j      (SIZE(Lns%j))       ); Wd%j       = Lns%j
   ALLOCATE( Wd%e      (SIZE(Lns%e))       ); Wd%e       = CMPLX(0d0,0d0,KIND=8)

   ALLOCATE( Zd%i      (SIZE(Lns%i))       ); Zd%i       = Lns%i
   ALLOCATE( Zd%i_mumps(SIZE(Lns%i_mumps)) ); Zd%i_mumps = Lns%i_mumps
   ALLOCATE( Zd%j      (SIZE(Lns%j))       ); Zd%j       = Lns%j
   ALLOCATE( Zd%e      (SIZE(Lns%e))       ); Zd%e       = CMPLX(0d0,0d0,KIND=8)

   CALL par_mumps_master (INITIALIZATION, 5, Wd, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   5, Wd, 0)

   ALLOCATE( Wa%i      (SIZE(Lns%i))       ); Wa%i       = Lns%i
   ALLOCATE( Wa%i_mumps(SIZE(Lns%i_mumps)) ); Wa%i_mumps = Lns%i_mumps
   ALLOCATE( Wa%j      (SIZE(Lns%j))       ); Wa%j       = Lns%j
   ALLOCATE( Wa%e      (SIZE(Lns%e))       ); Wa%e       = CMPLX(0d0,0d0,KIND=8)

   ALLOCATE( Za%i      (SIZE(Lns%i))       ); Za%i       = Lns%i
   ALLOCATE( Za%i_mumps(SIZE(Lns%i_mumps)) ); Za%i_mumps = Lns%i_mumps
   ALLOCATE( Za%j      (SIZE(Lns%j))       ); Za%j       = Lns%j
   ALLOCATE( Za%e      (SIZE(Lns%e))       ); Za%e       = CMPLX(0d0,0d0,KIND=8)

   CALL par_mumps_master (INITIALIZATION, 8, Wa, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   8, Wa, 0)

   ! (4)
   ! prepare matrices for first step of the
   ! time stepper with implicit Euler
   !
   ! NOTE: Z0 = Z0'
   !
   ALLOCATE( W0%i      (SIZE(Wd%i))       ); W0%i       = Wd%i
   ALLOCATE( W0%i_mumps(SIZE(Wd%i_mumps)) ); W0%i_mumps = Wd%i_mumps
   ALLOCATE( W0%j      (SIZE(Wd%j))       ); W0%j       = Wd%j
   ALLOCATE( W0%e      (SIZE(Wd%e))       ); W0%e       = CMPLX(0d0,0d0,KIND=8)

   ALLOCATE( Z0%i      (SIZE(Zd%i))       ); Z0%i       = Zd%i
   ALLOCATE( Z0%i_mumps(SIZE(Zd%i_mumps)) ); Z0%i_mumps = Zd%i_mumps
   ALLOCATE( Z0%j      (SIZE(Zd%j))       ); Z0%j       = Zd%j
   ALLOCATE( Z0%e      (SIZE(Zd%e))       ); Z0%e       = CMPLX(0d0,0d0,KIND=8)

   CALL par_mumps_master (INITIALIZATION, 6, W0, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   6, W0, 0)

!-------------------------------
! IF CHOSEN, READ TAUs FROM FILE

   IF ( p_in%tranGrowth_tau == dataIdentifier ) THEN

      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** WARNING:                      ***'
      WRITE(*,*) '*** this feature has a bug that   ***'
      WRITE(*,*) '*** prevents it from working.     ***'
      WRITE(*,*) '*** The first tau is correctly    ***'
      WRITE(*,*) '*** computed but not the          ***'
      WRITE(*,*) '*** following ones.               ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)

      WRITE(*,*)
      WRITE(*,*) '--> Reading taus from file taus.in'

      OPEN( UNIT = 20, FILE = 'taus.in' )
      READ(20,*) tausNumber
      READ(20,*)
      ALLOCATE(taus(tausNumber))
      DO ii = 1, tausNumber
         READ(20,*) taus(ii)
      ENDDO
      CLOSE(20)
      WRITE(*,*) '    read ', tausNumber, ' taus'
      WRITE(*,*) '    Done.'

   ENDIF

!--------------
! INITIAL GUESS

   WRITE(*,*)

   ! look for restart file (this feauture is not implemented at the moment due to work in
   !                        progress on the 3d formulation, restart files are never saved)
   !
   tauNameTemp = NINT(p_in%tranGrowth_tau*1e3)
   CALL intToChar6 (tauNameTemp,  tauName)
   WRITE(restart_name, '(A)') './tranGrowthOut/tranGrowthRestart-tau'//trim(tauName)//'.bin'
   INQUIRE (FILE=restart_name , EXIST=existFlag)

   IF ( existFlag ) THEN
      !
      ! read restart file and use it as initial guess
      !
      IF ( p_in%tranGrowth_initGuess == 1 ) THEN
         WRITE(*,*) '************************************************'
         WRITE(*,*) '*** Warning:                                 ***'
         WRITE(*,*) '*** I found a restart but you told me to use ***'
         WRITE(*,*) '*** ARPACK''s random guess.                   ***'
         WRITE(*,*) '*** These two choiches are not compatible.   ***'
         WRITE(*,*) '************************************************'
         WRITE(*,*) 'STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
      ENDIF
      CALL read_cmplx_restart_bin(xx_tg(1:velCmpnnts*np), dummy, trim(restart_name))
      CALL extract_cmplx (xx_tg,  uInitGuess)
      WRITE(*,*)
      WRITE(*,*) '    initial guess: ****** READ FROM RESTART FILE ******'
      WRITE(*,*)

   ELSE
      !
      ! create initial guess
      !
      SELECT CASE ( p_in%tranGrowth_initGuess )
         CASE (1)
            ! (a)
            ! random guess created by ARPACK
            ! WARNING: the velocity is NOT zero on borders
            WRITE(*,*) '    initial guess: ARPACK''s random'
            WRITE(*,*)
         CASE (2)
            ! (b)
            ! pseudo-random guess
            !
            WRITE(*,*) '    initial guess: pseudo-random'
            WRITE(*,*)
            CALL init_random_seed()
            CALL RANDOM_NUMBER(xx)
            CALL extract (xx,  uu)
            CALL projectDiv (uu,  u0)
            uInitGuess = CMPLX(u0,0d0,KIND=8)
         CASE (3)
            ! (c)
            ! base flow
            WRITE(*,*) '    initial guess: base flow'
            WRITE(*,*)
            uInitGuess = CMPLX(u0,0d0,KIND=8)
         CASE (4)
            ! (d)
            ! already computed eigenvector
            WRITE(*,*) '    initial guess: eigenvector'
            WRITE(*,*)
            CALL read_eigenvector (Nx, 1, './tranGrowthOut/eigenTranGrowth.dat',  xx, x0)
            xx_tg = CMPLX(xx,x0,KIND=8)
            CALL extract_cmplx (xx_tg,  uInitGuess)
         CASE DEFAULT
            WRITE(*,*) '*************************************'
            WRITE(*,*) '*** Wrong parameter:              ***'
            WRITE(*,*) '*** p_in % tranGrowth_initGuess   ***'
            WRITE(*,*) '*** set to: ', p_in%tranGrowth_initGuess
            WRITE(*,*) '*************************************'
            WRITE(*,*) 'STOP.'
            CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
      END SELECT

   ENDIF
   p0_tg = CMPLX(0d0,0d0,KIND=8)
   !CALL vtk_plot_P2 (rr, jj, jj_L,  DBLE(uInitGuess), 0*p0, './tranGrowthOut/tranGrowthInitGuess_Re.vtk')
   !CALL vtk_plot_P2 (rr, jj, jj_L, AIMAG(uInitGuess), 0*p0, './tranGrowthOut/tranGrowthInitGuess_Im.vtk')


!--------------------
! OUTER CYCLE ON TAUs
! note: at the moment the matrices for the time stepper are filled
!       each time tau is changed, this is not very convenient but is done
!       this way to end at exactly t = tau

DO ii = 1, tausNumber

   IF ( tausNumber /= 1 ) THEN
      p_in%tranGrowth_tau = taus(ii)
      WRITE(*,*) '    tau number: ', ii
   ENDIF

   ! Create very elaborate and stupid file name
   tauNameTemp = NINT(p_in%tranGrowth_tau*1e3)
   CALL intToChar6 (tauNameTemp,  tauName)

!---------------
! SET TIME STEPS

   nSteps = NINT(p_in%tranGrowth_tau / p_in%tranGrowth_dt)
   dtE    = p_in%tranGrowth_dt / 100
   dtH    = (p_in%tranGrowth_tau - dtE)/nSteps

!-----------------------------------
! FILL MATRICES FOR THE TIME STEPPER

   WRITE(*,*)
   SELECT CASE ( p_in%tranGrowth_method )
      CASE (1)
         ! explicit Euler
         !
         WRITE(*,*) '    method: explicit Euler'
         WRITE(*,*) '            NOT IMPLEMENTED'
         WRITE(*,*) '            AND PROBABLY UNUSABLE'
         WRITE(*,*) '            STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
         !
         ! Wd  = [ 1/dtH Mass  + div(u) ]
         !
         !
         ! Zd  = [ 1/dtH Mass  -   Lns (EXCEPT div(u)) ]
         !

      CASE (2)
         ! Crank Nicolson
         !
         WRITE(*,*) '    method: Crank Nicolson'
         WRITE(*,*)
         !
         ! Wd  = [ 1/dtH Mass  +  1/2 Lns ]
         !
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
         CALL extract(x_vec, u0)
         CALL qc_1y1_sp_gg_3d_M  (mm, jj,               1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
         CALL qc_oseen2y_sp_3d_M (mm, jj,                u0,    beta,  Lns_cmplx) ! + linearized terms
         CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,        -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
         CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,        -2d0,    beta,  Lns_cmplx) ! - velocity divergence
         CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg, 1d0,           Lns_cmplx) ! Dirichlet BCs
         IF (DESINGULARIZE_tg) THEN
            ! row
            DO i = Lns_cmplx%i(Nx), Lns_cmplx%i(Nx + 1) - 1
               Lns_cmplx%e(i) = CMPLX(0d0,0d0,KIND=8)
               IF (Lns_cmplx%j(i) == Nx) Lns_cmplx%e(i) = CMPLX(1d0,0d0,KIND=8)
            ENDDO
         ENDIF
!write(*,*) 'maxvals0d ', MAXVAL(DBLE(Lns_cmplx%e)), MAXVAL(AIMAG(Lns_cmplx%e))
!write(*,*) 'minvals0d ', MINVAL(DBLE(Lns_cmplx%e)), MINVAL(AIMAG(Lns_cmplx%e))
         CALL zAlinB_s (CMPLX(1d0/dtH,0d0,KIND=8), Mass_cmplx%e, &
                        CMPLX(0.5d0,0d0,KIND=8),   Lns_cmplx%e,  Wd%e)
         CALL par_mumps_master (NUMER_FACTOR, 5, Wd, 0)
         !
         ! Zd  = [ 1/dtH Mass  -  1/2 Lns (EXCEPT div(u)) ]
         !
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
         CALL extract(x_vec, u0)
         CALL qc_1y1_sp_gg_3d_M  (mm, jj,               1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
         CALL qc_oseen2y_sp_3d_M (mm, jj,                u0,    beta,  Lns_cmplx) ! + linearized terms
         CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,        -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
         CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg, 1d0,           Lns_cmplx) ! Dirichlet BCs
         CALL zAlinB_s (CMPLX(1d0/dtH,0d0,KIND=8), Mass_cmplx%e, &
                        CMPLX(-0.5d0,0d0,KIND=8),  Lns_cmplx%e,  Zd%e)
         !
         ! Wa  = [ 1/dtH Mass'  +  1/2 Lns' ]
         !
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
         CALL extract(x_vec, u0)
         CALL qc_1y1_sp_gg_3d_M  (mm, jj,                 1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
         CALL qc_oseen2y_sp_3d_M (mm, jj,                  u0,    beta,  Lns_cmplx) ! + linearized terms
         CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,          -2d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
         CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! - velocity divergence
         CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg,   1d0,           Lns_cmplx) ! Dirichlet BCs
         IF (DESINGULARIZE_tg) THEN
            ! column
            WHERE (Lns_cmplx%j == Nx)
               Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
            ENDWHERE
            IF (Lns_cmplx%j(SIZE(Lns_cmplx%j,1)) == Nx) Lns_cmplx%e(SIZE(Lns%j,1)) = CMPLX(1d0,0d0,KIND=8)
         ENDIF
         CALL zAlinB_s (CMPLX(1d0/dtH,0d0,KIND=8), Mass_cmplx%e, &
                        CMPLX(0.5d0,0d0,KIND=8),   Lns_cmplx%e,  Wa%e)
         CALL par_mumps_master (NUMER_FACTOR, 8, Wa, 0)
         !
         ! Za  = [ 1/dtH Mass'  -  1/2 Lns' (EXCEPT div(u)) ]
         !
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
         CALL extract(x_vec, u0)
         CALL qc_1y1_sp_gg_3d_M  (mm, jj,                 1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
         CALL qc_oseen2y_sp_3d_M (mm, jj,                  u0,    beta,  Lns_cmplx) ! + linearized terms
         CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,          -1d0,    beta,  Lns_cmplx) ! - velocity divergence
         CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg,   1d0,           Lns_cmplx) ! Dirichlet BCs
         CALL zAlinB_s (CMPLX(1d0/dtH,0d0,KIND=8), Mass_cmplx%e, &
                        CMPLX(-0.5d0,0d0,KIND=8),  Lns_cmplx%e,  Za%e)

      CASE (3)
         ! implicit Euler
         !
         WRITE(*,*) '    method: implicit Euler'
         WRITE(*,*) '            NOT IMPLEMENTED'
         WRITE(*,*) '            STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
         !
         ! Wd  = [ 1/dtH Mass  +  Lns ]
         !
         !
         ! Zd  = [ 1/dtH Mass ]
         !

      CASE DEFAULT
         WRITE(*,*) '*************************************'
         WRITE(*,*) '*** Wrong parameter:              ***'
         WRITE(*,*) '*** p_in % tranGrowth_method      ***'
         WRITE(*,*) '*** set to: ', p_in%tranGrowth_method
         WRITE(*,*) '*************************************'
         WRITE(*,*) 'STOP.'
         CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   END SELECT
   WRITE(*,*)

   !
   ! W0  = [ 1/dtE Mass  +  Lns ]
   !
   Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
   CALL extract(x_vec, u0)
   CALL qc_1y1_sp_gg_3d_M  (mm, jj,               1d0/Re, beta,  Lns_cmplx) ! stifness (GRAD:GRAD)
   CALL qc_oseen2y_sp_3d_M (mm, jj,                u0,    beta,  Lns_cmplx) ! + linearized terms
   CALL qc_1y0_sp_3d_M     (mm, jj, jj_L,        -1d0,    beta,  Lns_cmplx) ! + pressure gradient (ibp)
   CALL qc_0y1_sp_3d_M     (mm, jj, jj_L,        -1d0,    beta,  Lns_cmplx) ! - velocity divergence
   CALL Dirichlet_rc_3d_M  (np, js_Axis, js_D_tg, 1d0,           Lns_cmplx) ! Dirichlet BCs
   IF (DESINGULARIZE_tg) THEN
      ! column
      WHERE (Lns_cmplx%j == Nx)
         Lns_cmplx%e = CMPLX(0d0,0d0,KIND=8)
      ENDWHERE
      ! row
      DO i = Lns_cmplx%i(Nx), Lns_cmplx%i(Nx + 1) - 1
         Lns_cmplx%e(i) = CMPLX(0d0,0d0,KIND=8)
         IF (Lns_cmplx%j(i) == Nx) Lns_cmplx%e(i) = CMPLX(1d0,0d0,KIND=8)
      ENDDO
   ENDIF
!write(*,*) 'maxvals00 ', MAXVAL(DBLE(Lns_cmplx%e)), MAXVAL(AIMAG(Lns_cmplx%e))
!write(*,*) 'minvals00 ', MINVAL(DBLE(Lns_cmplx%e)), MINVAL(AIMAG(Lns_cmplx%e))
   CALL zAlinB_s (CMPLX(1d0/dtE,0d0,KIND=8), Mass_cmplx%e, &
                  CMPLX(1d0,0d0,KIND=8),     Lns_cmplx%e,  W0%e)
   CALL par_mumps_master (NUMER_FACTOR, 6, W0, 0)
   !
   ! Z0  = [ 1/dtE Mass ]
   !
   Z0%e =  CMPLX(1d0/dtE,0d0,KIND=8) * Mass_cmplx%e

!do i = 1, size( js_D_tg(2)%DIL )
!   do k = Wd%i(js_D_tg(2)%DIL(i)+np), Wd%i(js_D_tg(2)%DIL(i)+np+1)-1
!      if ( Wd%e(k) /= 0d0 ) then
!         write(*,*) i, js_D_tg(2)%DIL(i)+np, Wd%i_mumps(k), Wd%j(k), Wd%e(k)
!      endif
!   enddo
!enddo
!CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)

   WRITE(*,*) '    Matrices assembled correctly.'
   WRITE(*,*)

!-------------------------
! COMPUTE TRANSIENT GROWTH

   CALL singularValueDecomposition(nsv,   uInitGuess, p_in%tranGrowth_maxit, p_in%tranGrowth_tol, &
                                   sigma, nSteps,     singularValue, singularVector, iteNum)

!--------------------------------
! NORMALIZE RIGHT-SINGULAR VECTOR

   ! there is no need as ARPACK already normalizes it so that
   ! v'*v = 1

!-------------------
! ALLOCATE tmpVector

   ALLOCATE( tmpVector(SIZE(singularVector,1), 2*SIZE(singularVector,2)) )

!----------------------------------
! REBUILD THE USUAL VELOCITY VECTOR

   DO k = 1, velCmpnnts
      u0_tg(k,:) = singularVector( (k-1)*np+1 : k*np, 1 )
   ENDDO

!---------------------------------------------------
! COMPUTE KINETIC ENERGY OF THE OPTIMAL PERTURBATION

   Ek0 = SQRT(  (u0_tg(1,:)*u0_tg(1,:))**2 &
              + (u0_tg(1,:)*u0_tg(2,:))**2 &
              + (u0_tg(1,:)*u0_tg(3,:))**2 &
              + (u0_tg(2,:)*u0_tg(1,:))**2 &
              + (u0_tg(2,:)*u0_tg(2,:))**2 &
              + (u0_tg(2,:)*u0_tg(3,:))**2 &
              + (u0_tg(3,:)*u0_tg(1,:))**2 &
              + (u0_tg(3,:)*u0_tg(2,:))**2 &
              + (u0_tg(3,:)*u0_tg(3,:))**2 ) / 2

   CALL zAtimx (tmpVector(:,1), MassV%e, MassV%j, MassV%i, singularVector(:,1))
   Ek0_s = SUM( singularVector(:,1) * tmpVector(:,1) ) / 2

!--------------------------
! PLOT OPTIMAL PERTURBATION

   IF ( p_in%write_plots_flag ) THEN
      ! PLOT OUTPUT IN VTK FORMAT
      CALL vtk_plot_P2        (rr, jj, jj_L, DBLE(u0_tg), DBLE(p0_tg), &
           './tranGrowthOut/'//'tranGrowthShape0-tau_'//trim(tauName)//trim(filenm)//'_Re.vtk')
      CALL vtk_plot_scalar_P2 (rr, jj,       DBLE(Ek0),    &
           './tranGrowthOut/'//'tranGrowthEk0-tau_'//trim(tauName)//trim(filenm)//'_Re.vtk')

      CALL vtk_plot_P2        (rr, jj, jj_L, AIMAG(u0_tg), AIMAG(p0_tg), &
           './tranGrowthOut/'//'tranGrowthShape0-tau_'//trim(tauName)//trim(filenm)//'_Im.vtk')
      CALL vtk_plot_scalar_P2 (rr, jj,       AIMAG(Ek0),    &
           './tranGrowthOut/'//'tranGrowthEk0-tau_'//trim(tauName)//trim(filenm)//'_Im.vtk')
   END IF

!-------------------
! CHECK THE RESIDUAL

   WRITE(*,*)
   WRITE(*,*) '    Check residual ...'

   tmpVector(:,1) = singularVector(:,1)
   CALL timeStepper(tmpVector(:,1), nSteps, 1,  tmpVector(:,2))

      ! rebuild the usual velocity vector
      DO k = 1, velCmpnnts
         uu_tg(k,:) = tmpVector( (k-1)*np+1 : k*np, 2 )
      ENDDO
      ! compute kinetic energy of the evolution of the optimal perturbation
      Ek = SQRT(  (uu_tg(1,:)*uu_tg(1,:))**2 &
                + (uu_tg(1,:)*uu_tg(2,:))**2 &
                + (uu_tg(1,:)*uu_tg(3,:))**2 &
                + (uu_tg(2,:)*uu_tg(1,:))**2 &
                + (uu_tg(2,:)*uu_tg(2,:))**2 &
                + (uu_tg(2,:)*uu_tg(3,:))**2 &
                + (uu_tg(3,:)*uu_tg(1,:))**2 &
                + (uu_tg(3,:)*uu_tg(2,:))**2 &
                + (uu_tg(3,:)*uu_tg(3,:))**2 ) / 2
      CALL zAtimx (tmpVector(:,1), MassV%e, MassV%j, MassV%i, tmpVector(:,2))
      Ek_s = SUM( tmpVector(:,2) * tmpVector(:,1) ) / 2
      ! plot evolution of the optimal perturbation
      IF ( p_in%write_plots_flag ) THEN
         ! PLOT OUTPUT IN VTK FORMAT
         CALL vtk_plot_P2        (rr, jj, jj_L, DBLE(uu_tg), DBLE(pp_tg), &
              './tranGrowthOut/'//'tranGrowthShape1-tau_'//trim(tauName)//trim(filenm)//'_Re.vtk')
         CALL vtk_plot_scalar_P2 (rr, jj,       DBLE(Ek),       &
              './tranGrowthOut/'//'tranGrowthEk1-tau_'//trim(tauName)//trim(filenm)//'_Re.vtk')

         CALL vtk_plot_P2        (rr, jj, jj_L, AIMAG(uu_tg), AIMAG(pp_tg), &
              './tranGrowthOut/'//'tranGrowthShape1-tau_'//trim(tauName)//trim(filenm)//'_Im.vtk')
         CALL vtk_plot_scalar_P2 (rr, jj,       AIMAG(Ek),       &
              './tranGrowthOut/'//'tranGrowthEk1-tau_'//trim(tauName)//trim(filenm)//'_Im.vtk')
      END IF

   CALL timeStepper(tmpVector(:,2), nSteps, 2,  tmpVector(:,1))

   resLinfty = MAXVAL(ABS( tmpVector(:,1) - singularValue(1)*singularVector(:,1) ))

   WRITE(*,*) '    |res|_L-infty = ', resLinfty
   WRITE(*,*)
   WRITE(*,*) '    Ek_0        = u0''*M*u0 = ', Ek0_s
   WRITE(*,*) '    Ek_1        = uu''*M*uu = ', Ek_s
   WRITE(*,*) '    Ek_1 / Ek_0 = (uu''*M*uu) / ( u0''*M*u0 )  = ', Ek_s / Ek0_s

!-------------
! SAVE RESULTS

   WRITE(*,*)
   WRITE(*,*) '    Saving results'

   INQUIRE (FILE='./tranGrowthOut/tranGrowth.dat', EXIST=existFlag)
   IF (.NOT.existFlag) THEN
      OPEN(UNIT=20, FILE='./tranGrowthOut/tranGrowth.dat', STATUS='new', ACTION='write')
      WRITE(20,*)   '         tau               method                dt' &
                  //'                      G                      Ek/Ek0                      resLinfty' &
                  //'                ite'
      WRITE(20,*)
   ELSE
      OPEN(UNIT=20, FILE='./tranGrowthOut/tranGrowth.dat', STATUS='old', POSITION='append', ACTION='write')
   ENDIF

   WRITE(20,*) p_in%tranGrowth_tau, p_in%tranGrowth_method, p_in%tranGrowth_dt, singularValue(1), Ek_s / Ek0_s, resLinfty, iteNum
   CLOSE(20)

!---------------------------------------------
! DEALLOCATE VARIABLES WHICH ARE tau DEPENDENT

   DEALLOCATE( singularValue, singularVector, tmpVector )

!---------------------------
! END OF OUTER CYCLE ON TAUs
! note: at the moment the matrices for the time stepper are filled
!       each time tau is changed, this is not very convenient but is done
!       this way to end at exactly t = tau

ENDDO

!-----------------------------------------------
! DEALLOCATE VARIABLES WHICH ARE tau INDEPENDENT

   IF ( p_in%tranGrowth_tau == dataIdentifier ) THEN
      DEALLOCATE(taus)
   ENDIF

   DEALLOCATE( Dir_tg )
   DO k = 1, velCmpnnts
      DEALLOCATE ( zero_bvs_D_tg(k)%DRL )
   ENDDO

   DEALLOCATE( Mass_cmplx%i, Mass_cmplx%i_mumps, Mass_cmplx%j, Mass_cmplx%e )

   CALL par_mumps_master (DEALLOCATION, 6, W0, 0)
   DEALLOCATE( W0%i, W0%i_mumps, W0%j, W0%e )
   DEALLOCATE( Z0%i, Z0%i_mumps, Z0%j, Z0%e )
   CALL par_mumps_master (DEALLOCATION, 5, Wd, 0)
   DEALLOCATE( Wd%i, Wd%i_mumps, Wd%j, Wd%e )
   DEALLOCATE( Zd%i, Zd%i_mumps, Zd%j, Zd%e )
   CALL par_mumps_master (DEALLOCATION, 8, Wa, 0)
   DEALLOCATE( Wa%i, Wa%i_mumps, Wa%j, Wa%e )
   DEALLOCATE( Za%i, Za%i_mumps, Za%j, Za%e )

   CALL par_mumps_master (DEALLOCATION, 9, MassV, 0)
   DEALLOCATE( MassV%i, MassV%i_mumps, MassV%j, MassV%e )

   DEALLOCATE( uu_tg, u0_tg, pp_tg, p0_tg, xx_tg, x0_tg )


END SUBROUTINE compute_transientGrowth

!------------------------------------------------------------------------------

SUBROUTINE singularValueDecomposition(nsv, uInitGuess, maxit, tol, sigma, nSteps,  singularValues, singularVectors, iteNum)
!
!     Call ARPACK to find a few of the
!     largest singular values(sigma) and corresponding right singular 
!     vectors (v) for the the matrix A by solving the symmetric problem:
!          
!                        (A'*A)*v = sigma*v
! 
!     where A is an m by n COMPLEX matrix.
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The NEV right singular vectors will be computed in   |
!     | the N by NCV array V.                                |
!     |                                                      |
!     | The NEV left singular vectors will be computed in    |
!     | the M by NEV array U.                                |
!     |                                                      |
!     | NSV is the number of singular values requested.      |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     %------------------------------------------------------%
!
   IMPLICIT NONE

!  %-----------------%
!  | Input variables |
!  %-----------------%
   COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN) :: uInitGuess
   INTEGER,                         INTENT(IN) :: nsv, maxit
   REAL(KIND=8),                    INTENT(IN) :: tol
   COMPLEX(KIND=8),                 INTENT(IN) :: sigma
   INTEGER,                         INTENT(IN) :: nSteps

!  %------------------%
!  | Output variables |
!  %------------------%
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: singularValues
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: singularVectors
   INTEGER                                      :: iteNum

!  %--------------%
!  | Local Arrays |
!  %--------------%
   INTEGER, DIMENSION(11) :: iparam
   INTEGER, DIMENSION(14) :: ipntr
   
   LOGICAL, DIMENSION(:), ALLOCATABLE :: lselect
   
   COMPLEX(KIND=8), DIMENSION(:),   ALLOCATABLE :: Ax, d, workd, workev, resid, workl
   COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: v
   REAL(KIND=8),    DIMENSION(:),   ALLOCATABLE :: rwork
   REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: rd

!  %---------------%
!  | Local Scalars |
!  %---------------%
   CHARACTER(LEN=1) :: bmat
   CHARACTER(LEN=2) :: which

   INTEGER          :: ido, m, n, ncv, lworkl, info, ierr, &
                       i, j, nconv

   LOGICAL          :: rvec

   INTEGER            :: tauNameTemp
   CHARACTER(LEN=6)   :: tauName
   CHARACTER(LEN=128) :: restart_name
 
!  %-----------------------------%
!  | BLAS & LAPACK routines used |
!  %-----------------------------%
   REAL(KIND=8), EXTERNAL :: dznrm2, zaxpy, dlapy2

!-----------------------------------------------------------------------
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Start of singularValueDecomposition'
   WRITE(*,*) ''

   

!
!  %-------------------------------------------------%
!  | The following sets dimensions for this problem. |
!  %-------------------------------------------------%
!
   !m = SIZE( Wd%i ) - 1
   !n = SIZE( Wd%i ) - 1 ! WARNING: we suppose Wd is a square matrix
   ! NO: the above would mean considering the pressure in the perturbation,
   !     we don't want this because we consider a "pseudo"-operator for the
   !     LNSE where we have  u_dot = A u  and no pressure
   !     thus:
   m = velCmpnnts*np
   n = velCmpnnts*np
!
!  %------------------------------------------------%
!  | Specifications for ARPACK usage are set        | 
!  | below:                                         |
!  |                                                |
!  |    1) NCV sets the length of the Arnoldi       |
!  |       factorization                            |
!  |                                                |
!  |    2) This is a standard problem               |
!  |         (indicated by bmat  = 'I')             |
!  |                                                |
!  |    4) Ask for the NSV singular values of       |
!  |       largest magnitude                        |
!  |         (indicated by which = 'LM')            |
!  |       See documentation in DSAUPD for the      |
!  |       other options SM, BE.                    | 
!  |                                                |
!  %------------------------------------------------%
!
   ncv   = 2 * nsv + 2
   bmat  = 'I'
   which = 'LM'
!
!  %-----------------------------------------------------%
!  | Specification of stopping rules and initial         |
!  | conditions before calling DSAUPD                    |
!  |                                                     |
!  |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
!  |               computed   true                       |
!  |                                                     |
!  |      If TOL .le. 0,  then TOL <- macheps            |
!  |              (machine precision) is used.           |
!  |                                                     |
!  | IDO  is the REVERSE COMMUNICATION parameter         |
!  |      used to specify actions to be taken on return  |
!  |      from DSAUPD. (See usage below.)                |
!  |                                                     |
!  |      It MUST initially be set to 0 before the first |
!  |      call to DSAUPD.                                |
!  |                                                     |
!  | INFO on entry specifies starting vector information |
!  |      and on return indicates error codes            |
!  |                                                     |
!  |      Initially, setting INFO=0 indicates that a     |
!  |      random starting vector is requested to         |
!  |      start the ARNOLDI iteration.  Setting INFO to  |
!  |      a nonzero value on the initial call is used    |
!  |      if you want to specify your own starting       |
!  |      vector (This vector must be placed in RESID.)  |
!  |                                                     |
!  | The work array WORKL is used in DSAUPD as           |
!  | workspace.  Its dimension LWORKL is set as          |
!  | illustrated below.                                  |
!  %-----------------------------------------------------%
!
   !lworkl = ncv*(ncv+8)    ! SYMMETRIC     -> DSAUPD
   !lworkl = 3*ncv**2+6*ncv ! NON SYMMETRIC -> DNAUPD
   lworkl = 3*ncv**2+5*ncv  ! SYMMETRIC -> ZNAUPD
   ido    = 0
   SELECT CASE ( p_in%tranGrowth_initGuess )
      CASE (1)
         ! (a)
         ! random guess created by ARPACK
         !
         info = 0
      CASE (2:)
         ! (b)
         ! user provided initial guess
         !
         info = 1
   END SELECT

   ALLOCATE(lselect(ncv))
   ALLOCATE(Ax(m), d(ncv), workd(3*n), workev(2*ncv), resid(n), workl(lworkl))
   ALLOCATE(v(n,ncv))
   ALLOCATE(rwork(ncv), rd(ncv,3))
!
!  %---------------------------------------------------%
!  | Specification of the initial guess                |
!  | the choice is based on the parameter              |
!  | p_in % tranGrowth_initGuess                       |
!  | read from input file 'program_data.in'            |
!  %---------------------------------------------------%
!
   SELECT CASE ( p_in%tranGrowth_initGuess )
      CASE (1)
         ! (a)
         ! random guess created by ARPACK
         !
      CASE (2:)
         ! (b)
         ! user provided initial guess
         !
         DO i = 1, velCmpnnts
            resid( (i-1)*np+1 : i*np ) = uInitGuess(i,:)
         ENDDO
   END SELECT
!
!  %---------------------------------------------------%
!  | Specification of Algorithm Mode:                  |
!  |                                                   |
!  | This program uses the exact shift strategy        |
!  | (indicated by setting IPARAM(1) = 1.)             |
!  | IPARAM(3) specifies the maximum number of Arnoldi |
!  | iterations allowed.  Mode 1 of DSAUPD is used     |
!  | (IPARAM(7) = 1). All these options can be changed |
!  | by the user. For details see the documentation in |
!  | DSAUPD.                                           |
!  %---------------------------------------------------%
!
   iparam(1) = 1     ! ishifts = use exact shifts
   iparam(3) = maxit ! maxit   = max number of Arnoldi iterations allowed
   iparam(7) = 1     ! mode    = DSAUPD mode
!
!
!
!  %------------------------------------------------%
!  | M A I N   L O O P (Reverse communication loop) |
!  %------------------------------------------------%
!
!
!
   ! tentative restart, does not work
   !
   !tauNameTemp = NINT(p_in%tranGrowth_tau*1e3)
   !CALL intToChar6 (tauNameTemp,  tauName)
   !WRITE(restart_name, '(A)') './tranGrowthOut/tranGrowthRestart-tau'//trim(tauName)//'.bin'

   i = 0
   DO
      i = i+1
   
      WRITE(*,*) 'SVD, iteration n. ', i
!
!     %---------------------------------------------%
!     | Repeatedly call the routine ZNAUPD and take | 
!     | actions indicated by parameter IDO until    |
!     | either convergence is indicated or maxit    |
!     | has been exceeded.                          |
!     %---------------------------------------------%
!
      CALL znaupd ( ido  , bmat , n     , which, &
                    nsv  , tol  , resid , ncv,   &
                    v    , n    , iparam, ipntr, &
                    workd, workl, lworkl, rwork, info )

      IF (ido .eq. -1 .OR. ido .eq. 1) THEN
!
!        %---------------------------------------%
!        | Perform matrix vector multiplications |
!        |              w <--- A*x       (av())  |
!        |              y <--- A'*w      (atv()) |
!        | The user should supply his/her own    |
!        | matrix vector multiplication routines |
!        | here that takes workd(ipntr(1)) as    |
!        | the input, and returns the result in  |
!        | workd(ipntr(2)).                      |
!        %---------------------------------------%
!
         ! call av (m, n, workd(ipntr(1)), ax) 
         CALL timeStepper(workd(ipntr(1):ipntr(1)+n-1), nSteps, 1,  Ax)
         ! call atv (m, n, ax, workd(ipntr(2)))
         CALL timeStepper(Ax, nSteps, 2,  workd(ipntr(2):ipntr(2)+n-1))
!
!        %-----------------------------------------%
!        | L O O P   B A C K to call DSAUPD again. |
!        %-----------------------------------------%

         ! save "restart", does not work
         !CALL write_restart_bin (workd(ipntr(2):ipntr(2)+n-1), p_in%tranGrowth_tau, i, maxit, restart_name)
!
         CYCLE

      ELSE IF ( ido .eq. 99) THEN
      
         WRITE(*,*) 'Done SVD, exit info: ', info
         EXIT

      END IF 

   END DO ! End main loop
!
!  %----------------------------------------%
!  | Either we have convergence or there is |
!  | an error.                              |
!  %----------------------------------------%
!
   IF ( info .lt. 0 ) THEN ! There is an error.
!
!     %--------------------------%
!     | Error message. Check the |
!     | documentation in DSAUPD. |
!     %--------------------------%
!
      WRITE(*,*) ' Error with DSAUPD, info = ', info
      STOP       ' Check the documentation of ZNAUPD.'
!
   ELSE ! We made it.
!
!     %--------------------------------------------%
!     | No fatal errors occurred.                  |
!     | Post-Process using ZNEUPD.                 |
!     |                                            |
!     | Computed singular values may be extracted. |  
!     |                                            |
!     | Singular vectors may also be computed now  |
!     | if desired.  (indicated by rvec = .true.)  | 
!     |                                            |
!     | The routine DSEUPD now called to do this   |
!     | post processing                            | 
!     %--------------------------------------------%
!        
      rvec = .true.

      CALL zneupd ( rvec , 'A'   , lselect, d     , v,      &
                    n    , sigma , workev , bmat  , n,      &
                    which, nsv   , tol    , resid , ncv,    &
                    v    , n     , iparam , ipntr , workd,  &
                    workl, lworkl, rwork  , ierr )
!
!     %-----------------------------------------------%
!     | Singular values are returned in the first     |
!     | column of the two dimensional array S         |
!     | and the corresponding right singular vectors  | 
!     | are returned in the first NEV columns of the  |
!     | two dimensional array V as requested here.    |
!     %-----------------------------------------------%
!
      IF ( ierr .ne. 0) THEN
!
!        %------------------------------------%
!        | Error condition:                   |
!        | Check the documentation of DSEUPD. |
!        %------------------------------------%
!
         WRITE(*,*) ' Error with ZNEUPD, info = ', ierr
         STOP       ' Check the documentation of DSEUPD.'

      ELSE

         ! Number of converged singular values
         nconv =  iparam(5)

         ALLOCATE(singularValues(nconv))
         ALLOCATE(singularVectors(n,nconv))

         singularValues  = d(1:nconv)
         singularVectors = v(:,1:nconv)

! compute the residual norm, not used for now
!!!           do 20 j=1, nconv
!!!
!!!                %---------------------------%
!!!                | Compute the residual norm |
!!!                |                           |
!!!                |   ||  A*x - lambda*x ||   |
!!!                |                           |
!!!                | for the NCONV accurately  |
!!!                | computed eigenvalues and  |
!!!                | eigenvectors.  (iparam(5) |
!!!                | indicates how many are    |
!!!                | accurate to the requested |
!!!                | tolerance)                |
!!!                %---------------------------%
!!!
!!!                call av(nx, v(1,j), ax)
!!!                call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
!!!                rd(j,1) = dble (d(j))
!!!                rd(j,2) = dimag (d(j))
!!!                rd(j,3) = dznrm2 (n, ax, 1)
!!!                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
!!!20          continue


!
!        %-------------------------------%
!        | Display computed residuals    |
!        %-------------------------------%
!
         CALL dmout(6, nconv, 3, rd, ncv, -6, &
                   'Ritz values (Real, Imag) and realtive residuals')
      END IF
!
!     %------------------------------------------%
!     | Print additional convergence information |
!     %------------------------------------------%
!
      IF ( info .eq. 1) THEN
         WRITE(*,*) ' '
         WRITE(*,*) ' Maximum number of iterations reached.'
         WRITE(*,*) ' '
      ELSE IF ( info .eq. 3) THEN
         WRITE(*,*) ' ' 
         WRITE(*,*) ' No shifts could be applied during implicit', &
                    ' Arnoldi update, try increasing NCV.'
         WRITE(*,*) ' '
      END IF

      WRITE(*,*) ' '
      WRITE(*,*) ' SVD '
      WRITE(*,*) ' ==== '
      WRITE(*,*) ' '
      WRITE(*,*) ' Size of the matrix is ', n
      WRITE(*,*) ' The number of Ritz values requested is ', nsv
      WRITE(*,*) ' The number of Arnoldi vectors generated (NCV) is ', ncv
      WRITE(*,*) ' What portion of the spectrum: ', which
      WRITE(*,*) ' The number of converged Ritz values is ', nconv 
      WRITE(*,*) ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
      WRITE(*,*) ' The number of OP*x is ', iparam(9)
      WRITE(*,*) ' The convergence criterion is ', tol
      WRITE(*,*) ' Shift used: ', sigma
      WRITE(*,*) ' '

      iteNum = iparam(9)

   END IF

   DEALLOCATE(rwork, rd)
   DEALLOCATE(v)
   DEALLOCATE(Ax, d, workd, workev, resid, workl)
   DEALLOCATE(lselect)


END SUBROUTINE singularValueDecomposition

!------------------------------------------------------------------------------

SUBROUTINE timeStepper (u0v, nSteps, dirAdj,  u1v)
!
! u0v: input  variable containing only the velocity components, NO PRESSURE
!      organized as a single column vector
!
! u1v: output variable containing only the velocity components, NO PRESSURE
!      organized as a single column vector
!
!++++++++++++++++++++++++CRANK--NICOLSON METHOD++++++++++++++++++++++++++++++++
!
   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: u0v
   INTEGER,                       INTENT(IN) :: nSteps
   INTEGER,                       INTENT(IN) :: dirAdj ! 1=direct, 2=adjoint
   ! output variable
   COMPLEX(KIND=8), DIMENSION(:) :: u1v
   ! local variables
   INTEGER           :: i, k
   ! temp
   CHARACTER(LEN=12) :: counter

   ! executable statements
!+++
if ( testForNaN(DBLE(u0v)).OR.testForNaN(AIMAG(u0v)) ) then
   write(*,*) 'NaN detected'
   write(*,*) 'in u0v'
   stop
endif
!+++

   ! (0)
   ! rebuild the usual velocity vector
   !
   DO k = 1, velCmpnnts
      u0_tg(k,:) = u0v( (k-1)*np+1 : k*np )
   ENDDO


   ! (1)
   ! build a full unknowns vector x0_tg = (u0_tg, p0_tg)
   !
   CALL collect_cmplx (u0_tg, 0d0*p0_tg,  x0_tg) !!! p0_tg ??? At the moment this is ok
                                                 ! because the first time step is computed
                                                 ! with implicit Euler which doesn't
                                                 ! need the pressure

!!+++
!if ( testForNaN(DBLE(x0_tg)).OR.testForNaN(AIMAG(x0_tg)) ) then
!   write(*,*) 'NaN detected'
!   write(*,*) 'in x0_tg collect'
!   stop
!endif
!!+++
   ! (2)
   ! time steps
   !
   IF ( dirAdj == 1 ) THEN
      !-------------------
      ! direct problem
      !-------------------
      !
      ! initial REALLY SMALL step with implicit Euler
      !
      WRITE(*,*) '    DIR - time = ', dtE
      ! (a)
      ! assemble RHS
!write(*,*) 'maxvals1 ', MAXVAL(DBLE(x0_tg)), MAXVAL(AIMAG(x0_tg))
      CALL zAtimx (xx_tg, Z0%e, Z0%j, Z0%i, x0_tg)
!write(*,*) 'maxvals2 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
      ! (b)
      ! impose homogeneous Dirichlet BC on RHS
      CALL Dirichlet_c_3d (np, js_Axis, js_D_tg, zero_bvs_D_tg,  xx_tg)
!write(*,*) 'maxvals3 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
      IF (DESINGULARIZE_tg) xx_tg(Nx) = 0d0
      ! (c)
      ! solve system
      CALL par_mumps_master (DIRECT_SOLUTION, 6, W0, 0, xx_tg)
!write(*,*) 'maxvals4 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
      ! (d)
      ! update solution
      x0_tg = xx_tg
      !
      DO i = 1, nSteps
         !
         WRITE(*,*) '    DIR - time = ', dtE + i*dtH
         ! (a)
         ! assemble RHS
!write(*,*) 'maxvals5 ', MAXVAL(DBLE(x0_tg)), MAXVAL(AIMAG(x0_tg))
         CALL zAtimx (xx_tg, Zd%e, Zd%j, Zd%i, x0_tg)
!write(*,*) 'maxvals6 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
         ! (b)
         ! impose homogeneous Dirichlet BC on RHS
         CALL Dirichlet_c_3d (np, js_Axis, js_D_tg, zero_bvs_D_tg,  xx_tg)
!write(*,*) 'maxvals7 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
         IF (DESINGULARIZE_tg) xx_tg(Nx) = 0d0
         ! (c)
         ! solve system
         CALL par_mumps_master (DIRECT_SOLUTION, 5, Wd, 0, xx_tg)
!write(*,*) 'maxvals8 ', MAXVAL(DBLE(xx_tg)), MAXVAL(AIMAG(xx_tg))
         ! (d)
         ! update solution
         x0_tg = xx_tg
         ! temporary
         ! plot
         !CALL intToChar6 (i,  counter)
         !CALL extract (xx_tg,  uu_tg, p0_tg)
         !CALL vtk_plot_P2 (rr, jj, jj_L, uu_tg, p0_tg, './tranGrowthOut/tranGrowthDir'//trim(counter)//'.vtk')
      ENDDO

   ELSEIF ( dirAdj == 2 ) THEN
      !-------------------
      ! adjoint problem
      !-------------------
      !
      ! initial REALLY SMALL step with implicit Euler
      !
      WRITE(*,*) '    DIR - time = ', dtE
      ! (a)
      ! assemble RHS
      CALL zAtimx_T (xx_tg, Z0%e, Z0%j, Z0%i, x0_tg)
      ! (b)
      ! impose homogeneous Dirichlet BC on RHS
      CALL Dirichlet_c_3d (np, js_Axis, js_D_tg, zero_bvs_D_tg,  xx_tg)
      IF (DESINGULARIZE_tg) xx_tg(Nx) = 0d0
      ! (c)
      ! solve system
      CALL par_mumps_master (TRANSP_SOLUTION, 6, W0, 0, xx_tg)
      ! (d)
      ! update solution
      x0_tg = xx_tg
      !
      DO i = 1, nSteps
         !
         WRITE(*,*) '    ADJ - time = ', dtE + i*dtH
         ! (a)
         ! assemble RHS
         CALL zAtimx_T (xx_tg, Za%e, Za%j, Za%i, x0_tg)
         ! (b)
         ! impose homogeneous Dirichlet BC on RHS
         CALL Dirichlet_c_3d (np, js_Axis, js_D_tg, zero_bvs_D_tg,  xx_tg)
         IF (DESINGULARIZE_tg) xx_tg(Nx) = 0d0
         ! (c)
         ! solve system
         CALL par_mumps_master (TRANSP_SOLUTION, 8, Wa, 0, xx_tg)
         ! (d)
         ! update solution
         x0_tg = xx_tg
         ! temporary
         ! plot
         !CALL intToChar6 (i,  counter)
         !CALL extract (xx_tg,  uu_tg, p0_tg)
         !CALL vtk_plot_P2 (rr, jj, jj_L, uu_tg, p0_tg, './tranGrowthOut/tranGrowthAdj'//trim(counter)//'.vtk')
      ENDDO

   ENDIF

   ! (3)
   ! output only the velocity field in single column format
   ! and extract the pressure field so that it can be used for the next steps
   !
   CALL extract_cmplx (xx_tg,  uu_tg, p0_tg)
   DO k = 1, velCmpnnts
      u1v( (k-1)*np+1 : k*np ) = uu_tg(k,:)
   ENDDO


   ! temporary
   ! plot results to check for mistakes
   !IF ( dirAdj == 1 ) THEN
   !   CALL vtk_plot_P2 (rr, jj, jj_L,  DBLE(uu_tg),  DBLE(p0_tg), &
   !        './tranGrowthOut/'//'tranGrowthDir_Re.vtk')
   !   CALL vtk_plot_P2 (rr, jj, jj_L, AIMAG(uu_tg), AIMAG(p0_tg), &
   !        './tranGrowthOut/'//'tranGrowthDir_Im.vtk')
   !CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   !ELSEIF ( dirAdj == 2 ) THEN
   !   CALL vtk_plot_P2 (rr, jj, jj_L,  DBLE(uu_tg),  DBLE(p0_tg), &
   !        './tranGrowthOut/'//'tranGrowthAdj_Re.vtk')
   !   CALL vtk_plot_P2 (rr, jj, jj_L, AIMAG(uu_tg), AIMAG(p0_tg), &
   !        './tranGrowthOut/'//'tranGrowthAdj_Im.vtk')
   !ENDIF

   WRITE(*,*)

END SUBROUTINE timeStepper

!------------------------------------------------------------------------------

SUBROUTINE init_random_seed()

   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid, t(2), s
   integer(8) :: count, tms
   
   write(*,*) '    Generating random seed ...'

   call random_seed(size = n)
   allocate(seed(n))
!  ! First try if the OS provides a random number generator
!  open(newunit=un, file="/dev/urandom", access="stream", &
!       form="unformatted", action="read", status="old", iostat=istat)
!  if (istat == 0) then
!     write(*,*) '    random number generator provided by the OS'
!     read(un) seed
!     close(un)
!  else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      write(*,*) '    using current time and pID'
      call system_clock(count)
      if (count /= 0) then
         t = transfer(count, t)
      else
         call date_and_time(values=dt)
         tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24 * 60 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
         t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (n >= 3) then
         seed(1) = t(1) + 36269
         seed(2) = t(2) + 72551
         seed(3) = pid
         if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
         end if
      else
         seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
      end if
!  end if
   call random_seed(put=seed)

   write(*,*) '    Done.'
   write(*,*)

END SUBROUTINE init_random_seed

!------------------------------------------------------------------------------

SUBROUTINE projectDiv (u_approx,  uDiv)

   IMPLICIT NONE
   ! input variable
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u_approx
   ! output variable
   REAL(KIND=8), DIMENSION(:,:) :: uDiv
   ! local variables
   TYPE(CSR_MUMPS_Matrix)                 :: prMatr  ! projection matrix
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: tmpVel
   REAL(KIND=8), DIMENSION(Nx)            :: rhsSol
   INTEGER :: i


   WRITE(*,*) '    Projecting random guess in divergence free space ...'

   ! (1)
   ! create projection matrix
   !
   ALLOCATE( prMatr%i      (SIZE(Wd%i))       ); prMatr%i       = Wd%i
   ALLOCATE( prMatr%i_mumps(SIZE(Wd%i_mumps)) ); prMatr%i_mumps = Wd%i_mumps
   ALLOCATE( prMatr%j      (SIZE(Wd%j))       ); prMatr%j       = Wd%j
   ALLOCATE( prMatr%e      (SIZE(Wd%e))       ); prMatr%e       = 0d0

   CALL par_mumps_master (INITIALIZATION, 7, prMatr, 0)
   CALL par_mumps_master (SYMBO_FACTOR,   7, prMatr, 0)

   CALL qc_0y0_zero_sp_M (mm, jj,        1d0,  prMatr) ! mass matrix
   CALL qc_1y0_sp_M      (mm, jj, jj_L, -1d0,  prMatr) ! + lambda gradient (ibp)
   CALL qc_0y1_sp_M      (mm, jj, jj_L, -1d0,  prMatr) ! - velocity divergence

   CALL Dirichlet_rc_M (np, js_Axis, js_D_tg, 1d0,  prMatr)
   IF (DESINGULARIZE_tg) THEN
      DO i = prMatr%i(Nx), prMatr%i(Nx + 1) - 1
         prMatr%e(i) = 0
         IF (prMatr%j(i) == Nx) prMatr%e(i) = 1
      ENDDO
   ENDIF

   CALL par_mumps_master (NUMER_FACTOR, 7, prMatr, 0)

   ! (2)
   ! create rhs
   !
   tmpVel = 0d0
   CALL qv_0y0_sp (mm, jj, u_approx, 1d0,  tmpVel) ! alternatively CALL qv_mass_grad_sp (mm, jj, u_approx, 0*uu_tg,  tmpVel)
   CALL collect (tmpVel, 0*p0,  rhsSol)
   CALL Dirichlet_c (np, js_Axis, js_D_tg, zero_bvs_D_tg,  rhsSol)

   ! (3)
   ! solve linear system
   !
   CALL par_mumps_master (DIRECT_SOLUTION, 7, prMatr, 0, rhsSol)

   ! (4)
   ! output solution
   !
   CALL extract (rhsSol,  uDiv)

!    ! (4b) TEMP toDo remove
!    ! div(u) control 
!    !
!    CALL par_mumps_master (DEALLOCATION, 7, prMatr, 0)
!    DEALLOCATE (prMatr%i, prMatr%i_mumps, prMatr%j, prMatr%e)
!    CALL start_matrix_2d_p1(np_L,  jj_L , js_L , prMatr)
!    CALL par_mumps_master (INITIALIZATION, 7, prMatr, 0)
!    CALL par_mumps_master (SYMBO_FACTOR,   7, prMatr, 0)
! 
!    prMatr%e = 0d0
!    CALL qs_00_L_M (mm, jj_L, 1.0d0,  prMatr)
!    CALL par_mumps_master (NUMER_FACTOR, 7, prMatr, 0)
!    CALL qs_01_hybrid_L_sp (mm, jj, jj_L, uDiv,  pp)
!    CALL par_mumps_master (DIRECT_SOLUTION, 7, prMatr, 0, pp)
!    CALL vtk_plot_P2 (rr, jj, jj_L, uDiv, pp, './tranGrowthOut/tranGrowthDiv.vtk')


   ! (5)
   ! deallocate matrix
   !
   CALL par_mumps_master (DEALLOCATION, 7, prMatr, 0)
   DEALLOCATE (prMatr%i, prMatr%i_mumps, prMatr%j, prMatr%e)

   WRITE(*,*) '    Done.'

END SUBROUTINE projectDiv

!------------------------------------------------------------------------------

! not changed to complex variables yet
!
! SUBROUTINE evolve_transientGrowth(x_vec)
! 
!    IMPLICIT NONE
! 
!    ! input variables
!    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x_vec   ! computed base flow
!    ! "output" variables
!    REAL(KIND=8), DIMENSION(np) :: Ek0, Ek     ! kinetic energy fields
!    REAL(KIND=8)                :: Ek0_s, Ek_s ! kinetic energy
!    REAL(KIND=8), DIMENSION(Nx) :: x_opt       ! computed optimal perturbation
!    ! local variables
!    REAL(KIND=8), DIMENSION(Nx) :: x_bfl   ! computed base flow SAVE
!    REAL(KIND=8), DIMENSION(Nx) :: x_optE  ! already evolved optimal perturbation (may not be used)
!    REAL(KIND=8), DIMENSION(velCmpnnts,np) :: uE
!    REAL(KIND=8), DIMENSION(np_L)          :: pE
!    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmpVector
!    LOGICAL :: existFlag
!    INTEGER :: i
!    INTEGER           :: time_plot=0
!    CHARACTER(LEN=12) :: counter
!    REAL(KIND=8) :: dummy
!    CHARACTER(LEN=128) :: restart_name
! 
! !----------------------------------------------------
!    WRITE(*,*) ''
!    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
!    WRITE(*,*) '--> CALL to evolve_transientGrowth'
!    WRITE(*,*) ''
! !----------------------------------------------------
! 
! !---------------
! ! SAVE BASE FLOW
! 
!    x_bfl = x_vec
! 
! !--------------------------
! ! READ OPTIMAL PERTURBATION
! 
!    CALL vtk_read_P2 (p_in%etg_opt_perturb, rr, jj, jj_L,  u0_tg)
!    CALL collect (u0_tg, p0_tg,  x_opt)
! 
! !-------------------------------------------
! ! CHECK IF THE EVOLUTION HAS ALREADY STARTED
! 
!    time_plot = NINT( p_in%dns_dtPlot*1d3 )
!    CALL intToChar6 (time_plot,  counter)
!    INQUIRE (FILE=trim(p_in%dns_output_directory)//'sol'//trim(counter)//'.vtk' , EXIST=existFlag)
! 
!    IF ( existFlag ) THEN
! 
!       DO WHILE ( existFlag )
!          time_plot = time_plot + NINT( p_in%dns_dtPlot*1d3 )
!          CALL intToChar6 (time_plot,  counter)
!          INQUIRE (FILE=trim(p_in%dns_output_directory)//'sol'//trim(counter)//'.vtk' , EXIST=existFlag)
!       ENDDO
! 
!       p_in%dns_tInit = time_plot/1d3 - p_in%dns_dtPlot
! 
!       time_plot = time_plot - NINT( p_in%dns_dtPlot*1d3 )
!       CALL intToChar6 (time_plot,  counter)
! 
!       WRITE(*,*)
!       WRITE(*,*) '    Using restart:'
!       WRITE(*,*) '    starting evolution from time: ', p_in%dns_tInit
!       WRITE(*,*)
! 
!       pE  = 0 ! unused
! 
!       CALL vtk_read_P2 (trim(p_in%dns_output_directory)//'sol'//trim(counter)//'.vtk', rr, jj, jj_L,  uE)
!       CALL collect (uE, pE,  x_optE)
! 
!       !-----------------------------------------------
!       ! x0_tg IS THE ALREADY EVOLVED OPTIMAL PERTURBATION
!       
!       x0_tg = x_optE
! 
!       ! read previous time step
!       WRITE(restart_name, '(A)') trim(p_in%dns_output_directory)//'restartDNS.bin'
!       INQUIRE (FILE=restart_name , EXIST=existFlag)
!       IF ( existFlag ) THEN
!          CALL read_restart_bin(x_optE, dummy, trim(restart_name))
!       ELSE
!          WRITE(*,*) '***************************************'
!          WRITE(*,*) '*** Error:                          ***'
!          WRITE(*,*) '*** could not find DNS restart file ***'
!          WRITE(*,*) '*** ', trim(restart_name)
!          WRITE(*,*) '***************************************'
!          WRITE(*,*) 'STOP.'
!          CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
!       ENDIF
! 
!       !---------------------------------------------------
!       ! EVOLVE OPTIMAL PERTURBATION FROM WHERE WE LEFT OFF
!       
!       CALL dns(x0_tg, x_optE)
! 
!    ELSE
! 
!       !---------------------------------------
!       ! SUM BASE FLOW AND OPTIMAL PERTURBATION
!       
!          x0_tg = x_bfl + x_opt
! 
!       !----------------------------
!       ! EVOLVE OPTIMAL PERTURBATION
!       
!          CALL dns(x0_tg)
! 
!    ENDIF
! 
!    existFlag = .FALSE.
! 
! 
! !-----------------------
! ! CREATE THE MASS MATRIX
! 
!    ALLOCATE( Mass%i      (SIZE(Jacobian%i))       ); Mass%i       = Jacobian%i
!    ALLOCATE( Mass%i_mumps(SIZE(Jacobian%i_mumps)) ); Mass%i_mumps = Jacobian%i_mumps
!    ALLOCATE( Mass%j      (SIZE(Jacobian%j))       ); Mass%j       = Jacobian%j
!    ALLOCATE( Mass%e      (SIZE(Jacobian%e))       ); Mass%e       = 0d0
!    CALL qc_0y0_zero_sp_M (mm, jj, 1d0, Mass)
!    ! create the MassV matrix (non-singolar mass only for the velocity)
!    !
!    CALL dEssM (Mass%e, Mass%j, Mass%i, velCmpnnts*np, MassV%e, MassV%j, MassV%i, MassV%i_mumps)
!    !
!    DEALLOCATE( Mass%i, Mass%i_mumps, Mass%j, Mass%e )
! 
! !-------------------
! ! ALLOCATE tmpVector
! 
!    ALLOCATE( tmpVector(velCmpnnts*np, 2) )
! 
! !-----------------------
! ! COMPUTE KINETIC ENERGY
! 
!    ! compute kinetic energy of the optimal perturbation
!    CALL extract (x_opt,  u0_tg)
!    DO i = 1, velCmpnnts
!       tmpVector( (i-1)*np+1 : i*np, 1 ) = u0_tg(i,:)
!    ENDDO
! 
!    Ek0 = SQRT(  (u0_tg(1,:)*u0_tg(1,:))**2 &
!               + (u0_tg(1,:)*u0_tg(2,:))**2 &
!               + (u0_tg(1,:)*u0_tg(3,:))**2 &
!               + (u0_tg(2,:)*u0_tg(1,:))**2 &
!               + (u0_tg(2,:)*u0_tg(2,:))**2 &
!               + (u0_tg(2,:)*u0_tg(3,:))**2 &
!               + (u0_tg(3,:)*u0_tg(1,:))**2 &
!               + (u0_tg(3,:)*u0_tg(2,:))**2 &
!               + (u0_tg(3,:)*u0_tg(3,:))**2 ) / 2
! 
!    CALL dAtimx (tmpVector(:,2), MassV%e, MassV%j, MassV%i, tmpVector(:,1))
!    Ek0_s = SUM( tmpVector(:,1) * tmpVector(:,2) ) / 2
! 
! 
!    ! compute kinetic energy of the evolution of the optimal perturbation
!    xx_tg = x0_tg - x_bfl
!    CALL extract (xx_tg,  uu_tg)
!    DO i = 1, velCmpnnts
!       tmpVector( (i-1)*np+1 : i*np, 1 ) = uu_tg(i,:)
!    ENDDO
! 
!    Ek = SQRT(  (uu_tg(1,:)*uu_tg(1,:))**2 &
!              + (uu_tg(1,:)*uu_tg(2,:))**2 &
!              + (uu_tg(1,:)*uu_tg(3,:))**2 &
!              + (uu_tg(2,:)*uu_tg(1,:))**2 &
!              + (uu_tg(2,:)*uu_tg(2,:))**2 &
!              + (uu_tg(2,:)*uu_tg(3,:))**2 &
!              + (uu_tg(3,:)*uu_tg(1,:))**2 &
!              + (uu_tg(3,:)*uu_tg(2,:))**2 &
!              + (uu_tg(3,:)*uu_tg(3,:))**2 ) / 2
!    CALL dAtimx (tmpVector(:,2), MassV%e, MassV%j, MassV%i, tmpVector(:,1))
!    Ek_s = SUM( tmpVector(:,1) * tmpVector(:,2) ) / 2
! 
!    WRITE(*,*)
!    WRITE(*,*) '    u0''*u0   = ', SUM(Ek0)
!    WRITE(*,*) '    u0''*M*u0 = ', Ek0_s
!    WRITE(*,*) '    uu''*uu   = ', SUM(Ek)
!    WRITE(*,*) '    uu''*M*uu = ', Ek_s
!    WRITE(*,*) '    (uu''*uu) / ( u0''*u0 )  = ', SUM(Ek) / SUM(Ek0)
!    WRITE(*,*) '    (uu''*M*uu) / ( u0''*M*u0 )  = ', Ek_s / Ek0_s
! 
! !-------------
! ! SAVE RESULTS
! 
!    WRITE(*,*)
!    WRITE(*,*) '    Saving results'
! 
!    INQUIRE (FILE='./tranGrowthOut/tranGrowthEvolution.dat', EXIST=existFlag)
!    IF (.NOT.existFlag) THEN
!       OPEN(UNIT=20, FILE='./tranGrowthOut/tranGrowthEvolution.dat', STATUS='new', ACTION='write')
!       WRITE(20,*)   '         tau             Ek/Ek0'
!       WRITE(20,*)
!    ELSE
!       OPEN(UNIT=20, FILE='./tranGrowthOut/tranGrowthEvolution.dat', STATUS='old', POSITION='append', ACTION='write')
!    ENDIF
! 
!    WRITE(20,*) p_in%dns_tEnd, Ek_s / Ek0_s
!    CLOSE(20)
! 
! !---------------------
! ! DEALLOCATE VARIABLES
! 
!    DEALLOCATE( tmpVector )
!    DEALLOCATE( MassV%i, MassV%i_mumps, MassV%j, MassV%e )
! 
! 
! END SUBROUTINE evolve_transientGrowth

!==============================================================================

END MODULE transient_growth
