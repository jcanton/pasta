MODULE axisym_boundary_values

   USE global_variables

   USE dynamic_structures

   USE Dirichlet_Neumann

   USE Gauss_points

   IMPLICIT NONE

CONTAINS

SUBROUTINE gen_dirichlet_boundary_values (rr, sides, Dir, jjs, js_D, in_bvs_D, bvs_D)

   ! Hypotesis:
   ! - Dir is false on the axis

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:),     INTENT(IN) :: rr
   INTEGER,      DIMENSION(:),       INTENT(IN) :: sides
   LOGICAL,      DIMENSION(:,:),     INTENT(IN) :: Dir
   INTEGER,      DIMENSION(:,:),     INTENT(IN) :: jjs
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: js_D ! List of Dirichlet nodes
   REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN) :: in_bvs_D

   TYPE(dyn_real_line), DIMENSION(:) :: bvs_D ! Dirichlet boundary values

   ! local variables
   INTEGER :: i, k, j
   LOGICAL, DIMENSION(number_of_sides) :: this_side
   INTEGER, DIMENSION(:), POINTER      :: nodes_on_this_side

!+++
   INTEGER                                   :: dataIdentifier = 313 ! Donald Duck's plate number!
   LOGICAL, DIMENSION(:), ALLOCATABLE        :: dataSides
!+++

   ! executable statements


   DO k = 1, velCmpnnts
      ! initialization
      bvs_D(k)%DRL = 0
   END DO

!+++
   ALLOCATE( dataSides(number_of_sides) )
   dataSides = .false.
   DO i = 1, number_of_sides
      IF ( in_bvs_D(1,i,1) == dataIdentifier ) THEN
         ! For this side read boundary values from file
         dataSides(i) = .true.
      ENDIF
   ENDDO
!+++

   DO i = 1, number_of_sides ! cycle on side number

      ! get the indexes of the nodes on the i-th side
      this_side = .false.;  this_side(i) = .true.
      CALL Dirichlet_nodes_gen (jjs, sides, this_side, nodes_on_this_side)

      IF ( .NOT.dataSides(i) ) THEN
         ! assign velocity profile with polinomial coefficients
         !
         DO k = 1, velCmpnnts ! cycle on velocity components
            IF ( Dir(k,i) ) THEN

               DO j = 1, SIZE(nodes_on_this_side) ! cycle on nodes on this side
                  
                  WHERE ( js_D(k)%DIL == nodes_on_this_side(j) )
                     bvs_D(k)%DRL = & 
                       in_bvs_D(k,i,1)                                &
                     + in_bvs_D(k,i,2)*rr(1,nodes_on_this_side(j))    &
                     + in_bvs_D(k,i,3)*rr(1,nodes_on_this_side(j))**2 &
                     + in_bvs_D(k,i,4)*rr(2,nodes_on_this_side(j))    &
                     + in_bvs_D(k,i,5)*rr(2,nodes_on_this_side(j))**2
                  END WHERE

               ENDDO

            ENDIF
         ENDDO
!+++
      ELSE
         ! assign velocity profile reading from file 
         !
         CALL read_BVS (i, in_bvs_D(1,i,2), nodes_on_this_side, velCmpnnts, rr, Dir, js_D,  bvs_D)

      ENDIF
!+++
      DEALLOCATE(nodes_on_this_side)
   ENDDO

   

END SUBROUTINE gen_dirichlet_boundary_values

!------------------------------------------------------------------------------

SUBROUTINE read_BVS (sideNumber, scaleFactor, nodes_on_this_side, velCmpnnts, rr, Dir, js_D,  bvs_D)

   IMPLICIT NONE
   ! input variables
   INTEGER,                          INTENT(IN) :: sideNumber
   REAL(KIND=8),                     INTENT(IN) :: scaleFactor
   INTEGER, DIMENSION(:), POINTER,   INTENT(IN) :: nodes_on_this_side
   INTEGER,                          INTENT(IN) :: velCmpnnts
   REAL(KIND=8), DIMENSION(:,:),     INTENT(IN) :: rr
   LOGICAL, DIMENSION(:,:),          INTENT(IN) :: Dir
   TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: js_D ! List of Dirichlet nodes
   ! output variables
   TYPE(dyn_real_line), DIMENSION(:) :: bvs_D ! Dirichlet boundary values
   ! local variables
   CHARACTER(LEN=64)                         :: dataFileName
   INTEGER                                   :: dataNumberOfPoints
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dataSideData
   INTEGER                                   :: j, k
   INTEGER                                   :: locXa, locXb, &
                                                locYa, locYb, &
                                                locA,  locB
   ! executable statements

   ! (1)
   ! create file name
   !
   IF ( sideNumber<10 ) THEN
      WRITE(dataFileName,'(a19,i1,a4)') './boundaryData/side', sideNumber, '.dat'
   ELSEIF ( sideNumber<100 ) THEN
      WRITE(dataFileName,'(a19,i2,a4)') './boundaryData/side', sideNumber, '.dat'
   ELSE
      WRITE(*,*) '*** Increase possible side number in read_BVS ***'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

   ! (2)
   ! read point values
   !
   WRITE(*,*) '--> Reading boundary values from file ', trim(dataFileName)

   OPEN( UNIT = 20, FILE = trim(dataFileName) )

   READ(20, *) dataNumberOfPoints

   ALLOCATE ( dataSideData(velCmpnnts + k_d, dataNumberOfPoints) )

   READ(20, *) ! jump one line

   DO j = 1, dataNumberOfPoints
      READ(20, *) dataSideData(:,j)
   ENDDO

   CLOSE(20)
   WRITE(*,*) '    Done.'

   dataSideData(1:velCmpnnts,:) = dataSideData(1:velCmpnnts,:) * scaleFactor

   ! (3)
   ! check that the end points of the read side correspond to the end
   ! point of the side I'm using
   !
   ! xmin
   IF ( dataSideData(velCmpnnts+1,1) /= MINVAL(rr(1, nodes_on_this_side)) ) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong xmin                    ***'
      WRITE(*,*) '*** read xmin = ', dataSideData(velCmpnnts+1,1)
      WRITE(*,*) '*** mesh xmin = ', MINVAL(rr(1, nodes_on_this_side))
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   ! ymin
   IF ( dataSideData(velCmpnnts+2,1) /= MINVAL(rr(2, nodes_on_this_side)) ) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong ymin                    ***'
      WRITE(*,*) '*** read ymin = ', dataSideData(velCmpnnts+2,1)
      WRITE(*,*) '*** mesh ymin = ', MINVAL(rr(2, nodes_on_this_side))
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   ! xmax
   IF ( dataSideData(velCmpnnts+1,dataNumberOfPoints) /= MAXVAL(rr(1, nodes_on_this_side)) ) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong xmax                    ***'
      WRITE(*,*) '*** read xmax = ', dataSideData(velCmpnnts+1,dataNumberOfPoints)
      WRITE(*,*) '*** mesh xmax = ', MAXVAL(rr(1, nodes_on_this_side))
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF
   ! ymax
   IF ( dataSideData(velCmpnnts+2,dataNumberOfPoints) /= MAXVAL(rr(2, nodes_on_this_side)) ) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** Wrong ymax                    ***'
      WRITE(*,*) '*** read ymax = ', dataSideData(velCmpnnts+2,dataNumberOfPoints)
      WRITE(*,*) '*** mesh ymax = ', MAXVAL(rr(2, nodes_on_this_side))
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

   ! (4)
   ! interpolate read values
   !
   WRITE(*,*) '    Interpolating values ...'

   DO k = 1, velCmpnnts
      IF ( Dir(k,sideNumber) ) THEN
         DO j = 1, SIZE(nodes_on_this_side)
            locXa = MINLOC(ABS( rr(1,nodes_on_this_side(j)) - dataSideData(velCmpnnts+1,:) ),1)
            locYa = MINLOC(ABS( rr(2,nodes_on_this_side(j)) - dataSideData(velCmpnnts+2,:) ),1)


            IF ( locXa /=1 .AND. locXa /= SIZE(dataSideData(velCmpnnts+1,:)) ) THEN
               IF ( MINLOC(ABS( rr(1,nodes_on_this_side(j)) &
                                - (/dataSideData(velCmpnnts+1,locXa-1), dataSideData(velCmpnnts+1,locXa+1)/) ),1) == 1 ) THEN
                  locXb = locXa-1
               ELSE
                  locXb = locXa+1
               ENDIF
            ELSE
               locXb = locXa
            ENDIF
            IF ( locYa /=1 .AND. locYa /= SIZE(dataSideData(velCmpnnts+2,:)) ) THEN
               IF ( MINLOC(ABS( rr(2,nodes_on_this_side(j)) &
                                - (/dataSideData(velCmpnnts+2,locYa-1), dataSideData(velCmpnnts+2,locYa+1)/) ),1) == 1 ) THEN
                  locYb = locYa-1
               ELSE
                  locYb = locYa+1
               ENDIF
            ELSE
               locYb = locYa
            ENDIF
            
            locA = MAXVAL( (/locXa,locYa/) )
            locB = MAXVAL( (/locXb,locYb/) )

            IF ( locA /= locB ) THEN
               where( js_D(k)%DIL == nodes_on_this_side(j) )
               bvs_D(k)%DRL = dataSideData(k,locA) + (dataSideData(k,locB) - dataSideData(k,locA)) &
                  / SQRT( (dataSideData(velCmpnnts+1,locB)-dataSideData(velCmpnnts+1,locA))**2 &
                        + (dataSideData(velCmpnnts+2,locB)-dataSideData(velCmpnnts+2,locA))**2 ) &
                  * SQRT( (rr(1,nodes_on_this_side(j))-dataSideData(velCmpnnts+1,locA))**2 &
                        + (rr(2,nodes_on_this_side(j))-dataSideData(velCmpnnts+2,locA))**2 )
               endwhere
            ELSE
               where( js_D(k)%DIL == nodes_on_this_side(j) )
                  bvs_D(k)%DRL = dataSideData(k,locA)
               endwhere
            ENDIF

         ENDDO
      ENDIF
   ENDDO
   DEALLOCATE( dataSideData )

   WRITE(*,*) '    Done.'


END SUBROUTINE read_BVS

!------------------------------------------------------------------------------

SUBROUTINE write_BVS (sideNumber, uu, rr, jjs, sides, file_name)

   IMPLICIT NONE
   ! input variables
   INTEGER,                          INTENT(IN) :: sideNumber
   REAL(KIND=8), DIMENSION(:,:),     INTENT(IN) :: uu
   REAL(KIND=8), DIMENSION(:,:),     INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:),     INTENT(IN) :: jjs
   INTEGER,      DIMENSION(:),       INTENT(IN) :: sides
   CHARACTER(*),                     INTENT(IN) :: file_name
   ! local variables
   LOGICAL, DIMENSION(:), ALLOCATABLE        :: this_side
   INTEGER, DIMENSION(:), POINTER            :: nodes_on_this_side
   CHARACTER(LEN=64)                         :: dataFolder='./boundaryData/'
   CHARACTER(LEN=64)                         :: dataFileName
   INTEGER                                   :: dataNumberOfNodes
   INTEGER                                   :: j

   ! executable statements

   ALLOCATE( this_side(MAXVAL(sides)) )

   ! get the indexes of the nodes on this side
   this_side = .false.;  this_side(sideNumber) = .true.
   CALL Dirichlet_nodes_gen (jjs, sides, this_side, nodes_on_this_side)

   dataNumberOfNodes = SIZE(nodes_on_this_side)

   IF ( sideNumber<10 ) THEN
      WRITE(dataFileName,'(a5,i1,a4)') '_side', sideNumber, '.dat'
   ELSEIF ( sideNumber<100 ) THEN
      WRITE(dataFileName,'(a5,i2,a4)') '_side', sideNumber, '.dat'
   ELSE
      WRITE(*,*) '*** Increase possible side number in read_BVS ***'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF


   WRITE(*,*) '--> Writing boundary values to file ', trim(dataFolder)//trim(file_name)//trim(dataFileName)

   OPEN( UNIT = 20, FILE = trim(dataFolder)//trim(file_name)//trim(dataFileName) )

   WRITE(20, *) dataNumberOfNodes

   WRITE(20, *) ! jump one line

   DO j = 1, dataNumberOfNodes
      WRITE(20, *) uu(:,nodes_on_this_side(j)), rr(:,nodes_on_this_side(j))
   ENDDO

   CLOSE(20)

   WRITE(*,*) '    Done.'

END SUBROUTINE write_BVS

!------------------------------------------------------------------------------

SUBROUTINE inner_vel_profile(U_med, aa,  coeffs)

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), INTENT(IN) :: U_med ! mean velocity in the tube
   REAL(KIND=8), INTENT(IN) :: aa    ! radius of the tube
   ! output variables
   REAL(KIND=8), DIMENSION(6) :: coeffs

   coeffs(1) = U_med*3/2        ! cost
   coeffs(2) = 0                ! z
   coeffs(3) = 0                ! z**2
   coeffs(4) = 0                ! R
   coeffs(5) = -U_med*3/2/aa**2 ! R**2
   coeffs(6) = 0                ! LOG(R)

END SUBROUTINE inner_vel_profile

!------------------------------------------------------------------------------

SUBROUTINE outer_vel_profile(U_med, aa, bb,  coeffs)

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), INTENT(IN) :: U_med ! mean velocity in the tube
   REAL(KIND=8), INTENT(IN) :: aa    ! outer radius
   REAL(KIND=8), INTENT(IN) :: bb    ! inner radius
   ! output variables
   REAL(KIND=8), DIMENSION(6) :: coeffs
   ! local variables
   REAL(KIND=8) :: constPart

   constPart = (aa-bb) / ( aa**2*bb - 2/3*aa**3 -bb**3/3 &
                         - 1/LOG(bb/aa)*(aa**2-bb**2)*(aa-bb-bb*LOG(aa)+bb*LOG(bb)) )

   coeffs(1) = U_med*constPart*( -aa**2 - (aa**2-bb**2)/LOG(bb/aa)*LOG(aa) ) ! cost
   coeffs(2) = 0                                                             ! z
   coeffs(3) = 0                                                             ! z**2
   coeffs(4) = 0                                                             ! R
   coeffs(5) = U_med*constPart                                               ! R**2
   coeffs(6) = U_med*constPart*(aa**2-bb**2)/LOG(bb/aa)                      ! LOG(R)

END SUBROUTINE outer_vel_profile

!==============================================================================

END MODULE axisym_boundary_values
