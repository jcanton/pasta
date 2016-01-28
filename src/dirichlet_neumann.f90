MODULE  Dirichlet_Neumann

!!!!!!!!!  USE sparse_matrix_profiles  ! eliminato


CONTAINS


!------------------------------------------------------------------------------


SUBROUTINE Dir_nodes_size (jjs, sides, Dir,  js_D_size)

   IMPLICIT NONE

   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjs
   INTEGER, DIMENSION(:),   INTENT(IN)  :: sides
   LOGICAL, DIMENSION(:),   INTENT(IN)  :: Dir
   INTEGER,                 INTENT(OUT) :: js_D_size

   INTEGER, DIMENSION(:), ALLOCATABLE :: j_all, j_dif
   INTEGER :: nj, ms, ns

!  count all nodes of the boundary elements which belong to
!  the sides where a Dirichlet boundary condition is prescribed

!  If the same node belongs to different sides of the boundary
!  the algorithm takes it into account only once.

   nj = 0
   DO ms = 1, SIZE(jjs, 2)
      IF ( Dir(sides(ms)) )  nj = nj + SIZE(jjs, 1)
   ENDDO

   ALLOCATE (j_all(nj), j_dif(nj))

!  collect all nodes of the boundary elements which belong to
!  the sides where a Dirichlet boundary condition is prescribed

   nj = 0

   DO ms = 1, SIZE(jjs, 2)

      IF ( Dir(sides(ms)) ) THEN

         DO ns = 1, SIZE(jjs, 1)
            nj = nj + 1
            j_all(nj) = jjs(ns, ms)
         ENDDO

      ENDIF

   ENDDO


   CALL sort_diff (j_all,  j_dif, js_D_size)

   DEALLOCATE (j_all, j_dif)

CONTAINS

SUBROUTINE sort_diff (a,  a_d, n_a_d)

!  sort in ascending order of the integer array  a  and generation
!  of the integer array  a_d  whose first  n_a_d  leading entries
!  contain different values in ascending order, while all the
!  remaining entries are set to zero

!  sorting by Shell's method.

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(INOUT) :: a
   INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
   INTEGER,               INTENT(OUT)   :: n_a_d

   INTEGER :: n, na, inc, i, j, k, ia

   na = SIZE(a)

!  sort phase

   IF (na == 0) THEN
      n_a_d = 0
      RETURN
   ENDIF

   inc = 1
   DO WHILE (inc <= na)
      inc = inc * 3
      inc = inc + 1
   ENDDO

   DO WHILE (inc > 1)
      inc = inc/3
      DO i = inc + 1, na
         ia = a(i)
         j = i
         DO WHILE (a(j-inc) > ia)
            a(j) = a(j-inc)
            j = j - inc
            IF (j <= inc) EXIT
         ENDDO
         a(j) = ia
      ENDDO
   ENDDO

!  compression phase

   n = 1
   a_d(n) = a(1)
   DO k = 2, na
      IF (a(k) > a(k-1)) THEN
         n = n + 1
         a_d(n) = a(k)
      ENDIF
   ENDDO

   n_a_d = n

   a_d(n_a_d + 1 : na) = 0

END SUBROUTINE sort_diff

END SUBROUTINE Dir_nodes_size


!------------------------------------------------------------------------------


SUBROUTINE Dir_nodes_gen (jjs, sides, Dir,  js_D)

   IMPLICIT NONE

   INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs
   INTEGER, DIMENSION(:),   INTENT(IN) :: sides
   LOGICAL, DIMENSION(:),   INTENT(IN) :: Dir
   INTEGER, DIMENSION(:),  INTENT(OUT) :: js_D

   INTEGER, DIMENSION(:), ALLOCATABLE :: j_all, j_dif
   INTEGER :: nj, ms, ns, nj_dif

!  count all nodes of the boundary elements which belong to
!  the sides where a Dirichlet boundary condition is prescribed

!  If the same node belongs to different sides of the boundary
!  the algorithm takes it into account only once.

   nj = 0
   DO ms = 1, SIZE(jjs, 2)
      IF ( Dir(sides(ms)) )  nj = nj + SIZE(jjs, 1)
   ENDDO

   ALLOCATE (j_all(nj), j_dif(nj))

!  collect all nodes of the boundary elements which belong to
!  the sides where a Dirichlet boundary condition is prescribed

   nj = 0
   DO ms = 1, SIZE(jjs, 2)

      IF ( Dir(sides(ms)) ) THEN

         DO ns = 1, SIZE(jjs, 1)
            nj = nj + 1
            j_all(nj) = jjs(ns, ms)  ! js(ns, ms)  errore ? 
         ENDDO

      ENDIF

   ENDDO

   CALL sort_diff (j_all,  j_dif, nj_dif)

   js_D = j_dif(1:nj_dif)

   DEALLOCATE (j_all, j_dif)


CONTAINS

SUBROUTINE sort_diff (a,  a_d, n_a_d)

!  sort in ascending order of the integer array  a  and generation
!  of the integer array  a_d  whose first  n_a_d  leading entries
!  contain different values in ascending order, while all the
!  remaining entries are set to zero

!  sorting by Shell's method.

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(INOUT) :: a
   INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
   INTEGER,               INTENT(OUT)   :: n_a_d

   INTEGER :: n, na, inc, i, j, k, ia

   na = SIZE(a)

!  sort phase

   IF (na == 0) THEN
      n_a_d = 0
      RETURN
   ENDIF

   inc = 1
   DO WHILE (inc <= na)
      inc = inc * 3
      inc = inc + 1
   ENDDO

   DO WHILE (inc > 1)
      inc = inc/3
      DO i = inc + 1, na
         ia = a(i)
         j = i
         DO WHILE (a(j-inc) > ia)
            a(j) = a(j-inc)
            j = j - inc
            IF (j <= inc) EXIT
         ENDDO
         a(j) = ia
      ENDDO
   ENDDO

!  compression phase

   n = 1
   a_d(n) = a(1)
   DO k = 2, na
      IF (a(k) > a(k-1)) THEN
         n = n + 1
         a_d(n) = a(k)
      ENDIF
   ENDDO

   n_a_d = n

   a_d(n_a_d + 1 : na) = 0

END SUBROUTINE sort_diff

END SUBROUTINE Dir_nodes_gen


!-------------------------------------------------------------------------------


SUBROUTINE Dirichlet_nodes_gen (jjs, sides, Dir,  js_D)
  ! qui NEWTON non sta passando l'intera variabile Dir
  ! ma solamente Dir(k,:) dove k = 1, 2
  ! cioe` prima (k=1) passa la riga contenente l'informazione
  ! sulla componente x della velocita` per tutti i segmenti di
  ! bordo (le colonne di Dir), poi (k=2) passa la riga contenente
  ! l'informazione sulla componente y della velocita` per tutti
  ! i segmenti di bordo (le colonne di Dir)

!  Generation of the list of surface nodes js_D for
!  enforcing Dirichlet boundary conditions 

   INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs
   INTEGER, DIMENSION(:),   INTENT(IN) :: sides
   LOGICAL, DIMENSION(:),   INTENT(IN) :: Dir
   INTEGER, DIMENSION(:),   POINTER    :: js_D
 
   INTEGER :: j_size

   CALL Dir_nodes_size (jjs, sides, Dir,  j_size)

!  IF (ASSOCIATED(js_D)) DEALLOCATE(js_D)

   ALLOCATE (js_D(j_size))

   CALL Dir_nodes_gen (jjs, sides, Dir,  js_D)

END SUBROUTINE Dirichlet_nodes_gen 


!-------------------------------------------------------------------------------


SUBROUTINE Neumann_elements_gen (mms, sides, Neu,  ms_N)

!  Generation of the list of surface elements ms_N for
!  enforcing NONhomogeneous Neumann conditions in weak  
!  form through the suorface integral   

!  mms is any (sub)set of the boundary elements

   INTEGER, DIMENSION(:), INTENT(IN) :: mms
   INTEGER, DIMENSION(:), INTENT(IN) :: sides
   LOGICAL, DIMENSION(:), INTENT(IN) :: Neu
   INTEGER, DIMENSION(:), POINTER    :: ms_N
 
   INTEGER :: m_size, m, ms

   m_size = COUNT( Neu(sides(mms)) )

   ALLOCATE ( ms_N(m_size) )
   
   m = 0
   
   DO ms = 1, SIZE(mms)
     
      IF ( Neu(sides(mms(ms))) )  THEN 
        
         m = m + 1
      
         ms_N(m) = mms(ms) 
      
      ENDIF
   
   ENDDO
   
   
END SUBROUTINE Neumann_elements_gen 


END MODULE  Dirichlet_Neumann
