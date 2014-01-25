MODULE vtk_plot
!
! Contains three subroutines to generate .vtk output files.
!
USE start_sparse_kit ! needed for the extract subroutine
!
!==============================================================================

CONTAINS
! vtk_plot_P2         :: scalar + vector fields on P2 grid
!                        plots the velocity field on the P2 grid stored in 'rr'
!
! vtk_plot_scalar_P2  :: scalar field  on P2 grid
!                        plots a scalar field on the P2 grid stored in 'rr'
!
! vtk_plot_P1         :: scalar field  on P1 grid
!                        plots the pressure field on the P1 grid stored in 'rr_L'


!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_P2 (rr, jj, jj_L, uu, pp, file_name)
!
! Author: Jacopo Canton
! E-mail: jacopo.canton@mail.polimi.it
! Last revision: 6/10/2013
!
!++++++++ 2D/3D version +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! - rr   :: coordinates of the nodes of the P2 grid (triangle vertices plus
!           midpoints)
! - jj   :: connectivity matrix of the P2 grid
! - uu   :: a vector field whose values are defined on the P2 nodes
! - pp   :: a scalar field whose values are defined on the P1 nodes
!
! - file_name :: a string for the file name
!

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jj_L
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: pp

   CHARACTER(*),                 INTENT(IN) :: file_name

   ! internal variables
   INTEGER :: m
   INTEGER :: me, np
   INTEGER :: uc
   REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: pp_P2
   CHARACTER(LEN=128) :: s_buffer
   character(1), parameter :: end_rec = char(10) ! end-character for binary-record finalize

   me = SIZE(jj, 2) ! number of volume elements
   np = SIZE(rr, 2) ! number of P2 nodes
   uc = SIZE(uu, 1) ! number of velocity components

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Plotting ' // trim(file_name) // ' ...'

   OPEN( UNIT       = 20,            &
         FILE       = file_name,     &
         FORM       = 'unformatted', &
         STATUS     = 'unknown',     &
         access     = 'SEQUENTIAL',  &
         action     = 'WRITE',       &
         convert    = 'BIG_ENDIAN',  &
         recordtype = 'STREAM',      &
         buffered   = 'YES')


   WRITE(20) '# vtk DataFile Version 3.0'//end_rec
   WRITE(20) trim(file_name)//end_rec
   WRITE(20) 'BINARY'//end_rec
   WRITE(20) 'DATASET UNSTRUCTURED_GRID'//end_rec

   WRITE(20) ! jump one line

   ! nodes coordinates
   WRITE(s_buffer, '(a7,i7,a7)') 'POINTS ', np, ' double'   
   write(20) trim(s_buffer)//end_rec
   DO m = 1, np
      WRITE(20) rr(1, m), rr(2, m), 0.0
   END DO
   
   WRITE(20) ! jump one line

   ! connectivity matrix
   WRITE(s_buffer, '(a6,i7,a1,i7)') 'CELLS ', me, ' ', me*7 ! 7 = number of numbers per row
   write(20) trim(s_buffer)//end_rec
   DO m = 1, me
      WRITE(20) 6, jj(1, m)-1, jj(2, m)-1, jj(3, m)-1, &
                   jj(6, m)-1, jj(4, m)-1, jj(5, m)-1
   END DO

   WRITE(20) ! jump one line

   ! type of cells
   WRITE(s_buffer, '(a11,i7)') 'CELL_TYPES ', me
   write(20) trim(s_buffer)//end_rec
   DO m = 1, me
      WRITE(20) '22' ! 22 = P2 triangles
   END DO

   WRITE(20) ! jump one line

   ! point data
   WRITE(s_buffer, '(a11,i7)') 'POINT_DATA ', np
   write(20) trim(s_buffer)//end_rec


   ! scalar field
   CALL plot_inject_p1_p2 (jj_L, jj, pp,  pp_P2)
   WRITE(s_buffer, '(a)') 'SCALARS p double'
   write(20) trim(s_buffer)//end_rec
   WRITE(s_buffer, '(a)') 'LOOKUP_TABLE default'
   write(20) trim(s_buffer)//end_rec

   DO m = 1, np
      WRITE(20) pp_P2(m)
   END DO
   WRITE(20) ! jump one line

   ! vector field
   WRITE(s_buffer, '(a)') 'VECTORS u double'
   write(20) trim(s_buffer)//end_rec

   IF ( uc == 3 ) THEN
      DO m = 1, np
         WRITE(20) uu(1, m), uu(2, m), uu(3,m)
      END DO
   ELSEIF ( uc == 2 ) THEN
      DO m = 1, np
         WRITE(20) uu(1, m), uu(2, m), '0.0'
      END DO
   ELSEIF ( uc == 1 ) THEN
      DO m = 1, np
         WRITE(20) uu(1, m), '0.0 0.0'
      END DO
   ELSE
      ! what kind of vectorial field would you like to plot?!?
      WRITE(*,*) '    wrong vectorial field dimensions'
      WRITE(*,*) '    STOP.'
      STOP
   END IF


   CLOSE(20)

   
   WRITE(*,*) '    Done.'

END SUBROUTINE vtk_plot_P2

!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_scalar_P2 (rr, jj, psi, file_name)
!
! Author: Jacopo Canton
! E-mail: jacopo.canton@mail.polimi.it
! Last revision: 14/3/2013
!
! - rr   :: coordinates of the nodes of the P2 grid (triangle vertices plus
!           midpoints)
! - jj   :: connectivity matrix of the P2 grid
! - psi  :: a scalar field whose values are defined on the P2 nodes
!
! - file_name :: a string for the file name

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: psi

   CHARACTER(*),                 INTENT(IN) :: file_name

   INTEGER :: m
   INTEGER :: me, np


   me = SIZE(jj, 2) ! number of volume elements
   np = SIZE(rr, 2) ! number of P2 nodes

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Plotting ' // trim(file_name) // ' ...'

   OPEN( UNIT = 20, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')

   WRITE(20, '(a)') '# vtk DataFile Version 3.0'
   WRITE(20, '(a)') file_name
   WRITE(20, '(a)') 'ASCII'
   WRITE(20, '(a)') 'DATASET UNSTRUCTURED_GRID'

   WRITE(20, *) ! jump one line

   ! nodes coordinates
   WRITE(20, '(a7,i7,a7)') 'POINTS ', np, ' double'
   DO m = 1, np
      WRITE(20, '(2(g17.10,1x),g17.10)') rr(1, m), rr(2, m), 0.0
   END DO
   
   WRITE(20, *) ! jump one line

   ! connectivity matrix
   WRITE(20, '(a6,i7,a1,i7)') 'CELLS ', me, ' ', me*7 ! 7 = number of numbers per row
   DO m = 1, me
      WRITE(20, '(i1,1x,5(i7,1x),i7)') 6, jj(1, m)-1, jj(2, m)-1, jj(3, m)-1, &
                                          jj(6, m)-1, jj(4, m)-1, jj(5, m)-1
   END DO

   WRITE(20, *) ! jump one line

   ! type of cells
   WRITE(20, '(a11,i7)') 'CELL_TYPES ', me
   DO m = 1, me
      WRITE(20, '(a)') '22' ! 22 = P2 triangles
   END DO

   WRITE(20, *) ! jump one line

   ! point data
   WRITE(20, '(a11,i7)') 'POINT_DATA ', np


   ! scalar field
   WRITE(20, '(a)') 'SCALARS psi double'
   WRITE(20, '(a)') 'LOOKUP_TABLE default'
   DO m = 1, np
      WRITE(20, '(g17.10)') psi(m)
   END DO
   WRITE(20, *) ! jump one line


   CLOSE(20)


   WRITE(*,*) '    Done.'

END SUBROUTINE vtk_plot_scalar_P2

!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_P1 (rr_L, jj_L, pp, file_name)
!
! Author: Jacopo Canton
! Last revision: 13/3/2013
!
! - rr_L :: coordinates of the nodes of the P1 grid (triangle vertices)
! - jj_L :: connectivity matrix of the P1 grid
! - pp   :: a scalar field whose values are defined on the P1 nodes
!
! - file_name :: a string for the file name

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr_L
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj_L
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: pp

   CHARACTER(*),                 INTENT(IN) :: file_name

   INTEGER :: m
   INTEGER :: me, np_L


   me   = SIZE(jj_L, 2) ! number of volume elements
   np_L = SIZE(rr_L, 2) ! number of P1 nodes

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Plotting ' // trim(file_name) // ' ...'

   OPEN( UNIT = 20, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')

   WRITE(20, '(a)') '# vtk DataFile Version 3.0'
   WRITE(20, '(a)') file_name
   WRITE(20, '(a)') 'ASCII'
   WRITE(20, '(a)') 'DATASET UNSTRUCTURED_GRID'

   WRITE(20, *) ! jump one line

   ! nodes coordinates
   WRITE(20, '(a7,i7,a7)')     'POINTS ', np_L, ' double'
   DO m = 1, np_L
      WRITE(20, '(2(g17.10,1x),g17.10)') rr_L(1, m), rr_L(2, m), 0.0
   END DO
   
   WRITE(20, *) ! jump one line

   ! connectivity matrix
   WRITE(20, '(a6,i7,a1,i7)') 'CELLS ', me, ' ', me*4 ! 4 = number of numbers per row
   DO m = 1, me
      WRITE(20, '(i1,1x,2(i7,1x),i7)') 3, jj_L(1, m)-1, jj_L(2, m)-1, jj_L(3, m)-1
   END DO

   WRITE(20, *) ! jump one line

   ! type of cells
   WRITE(20, '(a11,i7)') 'CELL_TYPES ', me
   DO m = 1, me
      WRITE(20, '(a)') '5' ! 5 = P1 triangles
   END DO

   WRITE(20, *) ! jump one line

   ! point data
   WRITE(20, '(a11,i7)') 'POINT_DATA ', np_L
   ! scalar field
   WRITE(20, '(a)') 'SCALARS p double'
   WRITE(20, '(a)') 'LOOKUP_TABLE default'
   DO m = 1, np_L
      WRITE(20, '(g17.10)') pp(m)
   END DO


   CLOSE(20)


   WRITE(*,*) '    Done.'

END SUBROUTINE vtk_plot_P1

!------------------------------------------------------------------------------

SUBROUTINE vtk_plot_eigenvectors (rr, jj, eigvecs, file_name)
!
! Author: Jacopo Canton
! E-mail: jacopo.canton@mail.polimi.it
! Last revision: 21/5/2013
!
!******************************************************************************
!******** 3D version ONLY: the velocity field has 3 components.  **************
!******** The declaration of variable u_plot and the actual      **************
!******** writing, both marked with asterisks, need to be        **************
!******** changed depending on the number of velocity components **************
!******************************************************************************
!
! - rr      :: coordinates of the nodes of the P2 grid (triangle vertices plus
!              midpoints)
! - jj      :: connectivity matrix of the P2 grid
! - eigvecs :: a matrix containing 'nev' eigenvectors, defined on the P2 nodes,
!              organized by columns
!
! - file_name :: a string for the file name
!
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: eigvecs

   CHARACTER(*),                 INTENT(IN) :: file_name

   INTEGER :: m, i
   INTEGER :: me, np, nev
   REAL(KIND=8), DIMENSION(3,SIZE(rr,2)) :: u_plot
   !                       ^
   !                       |
   !                      This value has to be changed between 2 and 3
   !                      depending on the number of components of the
   !                      velocity field
   !******************************************************************

   me  = SIZE(jj, 2)     ! number of volume elements
   np  = SIZE(rr, 2)     ! number of P2 nodes
   nev = SIZE(eigvecs,2) ! number of eigenvectors

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Plotting ' // trim(file_name) // ' ...'

   OPEN( UNIT = 20, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')

   WRITE(20, '(a)') '# vtk DataFile Version 3.0'
   WRITE(20, '(a)') file_name
   WRITE(20, '(a)') 'ASCII'
   WRITE(20, '(a)') 'DATASET UNSTRUCTURED_GRID'

   WRITE(20, *) ! jump one line

   ! nodes coordinates
   WRITE(20, '(a7,i7,a7)') 'POINTS ', np, ' double'
   DO m = 1, np
      WRITE(20, '(2(g17.10,1x),g17.10)') rr(1, m), rr(2, m), 0.0
   END DO
   
   WRITE(20, *) ! jump one line

   ! connectivity matrix
   WRITE(20, '(a6,i7,a1,i7)') 'CELLS ', me, ' ', me*7 ! 7 = number of numbers per row
   DO m = 1, me
      WRITE(20, '(i1,1x,5(i7,1x),i7)') 6, jj(1, m)-1, jj(2, m)-1, jj(3, m)-1, &
                                          jj(6, m)-1, jj(4, m)-1, jj(5, m)-1
   END DO

   WRITE(20, *) ! jump one line

   ! type of cells
   WRITE(20, '(a11,i7)') 'CELL_TYPES ', me
   DO m = 1, me
      WRITE(20, '(a)') '22' ! 22 = P2 triangles
   END DO

   WRITE(20, *) ! jump one line

   ! point data
   WRITE(20, '(a11,i7)') 'POINT_DATA ', np

   ! vector fields
   DO i = 1, nev
      u_plot = 0
      CALL extract(eigvecs(:,i), u_plot)
      IF ( i < 10 ) THEN
         WRITE(20, '(a11,i1,a7)') 'VECTORS eig', i, ' double'
      ELSEIF ( i < 100 ) THEN
         WRITE(20, '(a11,i2,a7)') 'VECTORS eig', i, ' double'
      ELSEIF ( i < 1000 ) THEN
         WRITE(20, '(a11,i3,a7)') 'VECTORS eig', i, ' double'
      ENDIF
      DO m = 1, np
         !*********************************************************************
         ! WRITE(20, '(2(g17.10,1x),a3)') u_plot(1, m), u_plot(2, m), '0.0'          ! 2 velocity components
         WRITE(20, '(2(g17.10,1x),g17.10)') u_plot(1, m), u_plot(2, m), u_plot(3, m) ! 3 velocity components
         !*********************************************************************
      END DO
      WRITE(20, *) ! jump one line
   END DO


   CLOSE(20)

   
   WRITE(*,*) '    Done.'

END SUBROUTINE vtk_plot_eigenvectors

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE plot_inject_p1_p2 (jj_L, jj_P, pp_L,  pp_P)

   ! Injection from Linear interpolation to Parabolic interpolation
   ! Vertex nodes of the parabolic triangle are the first three
   ! nodes of the element 

!****************** 2D grid only **********************************************

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_L, jj_P
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_L
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_P

   pp_P(1:SIZE(pp_L)) = pp_L

   pp_P(jj_P(4,:)) = (pp_L(jj_L(2,:)) + pp_L(jj_L(3,:)))*0.5
   pp_P(jj_P(5,:)) = (pp_L(jj_L(3,:)) + pp_L(jj_L(1,:)))*0.5
   pp_P(jj_P(6,:)) = (pp_L(jj_L(1,:)) + pp_L(jj_L(2,:)))*0.5

END SUBROUTINE plot_inject_p1_p2

!==============================================================================

SUBROUTINE vtk_read_P2 (file_name, rr, jj, jj_L, uu)
!
! Author: Jacopo Canton
! E-mail: jacopo.canton@mail.polimi.it
! Last revision: 5/12/2013
!
!++++++++ 2D/3D version +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! - rr   :: coordinates of the nodes of the P2 grid (triangle vertices plus
!           midpoints)
! - jj   :: connectivity matrix of the P2 grid
! - uu   :: a vector field whose values are defined on the P2 nodes
! - pp   :: a scalar field whose values are defined on the P1 nodes
!
! - file_name :: a string for the file name
!

   IMPLICIT NONE

   CHARACTER(*),                 INTENT(IN) :: file_name

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jj_L
   REAL(KIND=8), DIMENSION(:,:)             :: uu
   !REAL(KIND=8), DIMENSION(:)               :: pp

   ! internal variables
   INTEGER :: m
   INTEGER :: me, np
   INTEGER :: uc
   !REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: pp_P2


   me = SIZE(jj, 2) ! number of volume elements
   np = SIZE(rr, 2) ! number of P2 nodes
   uc = SIZE(uu, 1) ! number of velocity components

   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Reading ' // trim(file_name) // ' ...'

   OPEN( UNIT = 20, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')

   READ(20, *)
   READ(20, *)
   READ(20, *)
   READ(20, *)
   READ(20, *)

   ! nodes coordinates
   READ(20, *) ! 'POINTS ', np, ' double'
   DO m = 1, np
      READ(20, *)
   END DO
   
   READ(20, *) ! jump one line

   ! connectivity matrix
   READ(20, *) ! 'CELLS ', me, ' ', me*7 ! 7 = number of numbers per row
   DO m = 1, me
      READ(20, *)
   END DO

   READ(20, *) ! jump one line

   ! type of cells
   READ(20, *) ! 'CELL_TYPES ', me
   DO m = 1, me
      READ(20, *) ! '22' ! 22 = P2 triangles
   END DO

   READ(20, *) ! jump one line

   ! point data
   READ(20, *) ! 'POINT_DATA ', np

   ! scalar field
   READ(20, *) ! 'SCALARS p double'
   READ(20, *) ! 'LOOKUP_TABLE default'
   DO m = 1, np
      READ(20, *) ! pp_P2(m)
   END DO

   READ(20, *) ! jump one line

   ! vector field
   READ(20, *) ! 'VECTORS u double'

   DO m = 1, np
      READ(20, *) uu(1, m), uu(2, m), uu(3,m)
   END DO

   CLOSE(20)

   
   WRITE(*,*) '    Done.'

END SUBROUTINE vtk_read_P2

!------------------------------------------------------------------------------

!==============================================================================

END MODULE vtk_plot
