     !> \brief modulo che usa le funzioni di fem_2d_grid.f90
     !!     per leggere la mesh; poi crea i nodi per l'interpolazione
     !!     P2.
     !!
     !!  jj (n_w,  me)   nodes of the  volume_elements               
     !!                                            
     !!  jjs(n_ws, mes)  nodes of the surface_elements in volume numbering
     !!                                                           
     !!  iis(n_ws, mes)  nodes of the surface_elements in surface numbering         
     !!               anche nota come: matrice di connettivita` degli         
     !!               elementi di superficie, partendo dalla matrice di       
     !!               connettivita` jjs degli elementi di superficie          
     !!               'basata' sulla numerazione dei nodi di volume, e        
     !!               dall'array js dei nodi di bordo according alla          
     !!               numerazione dei nodi di volume.                   
     !!                                                           
     !!   js(nps)        nodes on boundary (surface)                          
     !!                                                           
     !!   mm(me)         (possibly sub) set of elements for quadrature              
     !!                                                           
     !!  mms(mes)        (possibly sub) set of surface_elements for surf_quadrature 
     !!                                                           
      !! INTEGER,      ALLOCATABLE, DIMENSION(:,:) :: jj, jjs, iis                     
      !! INTEGER,      ALLOCATABLE, DIMENSION(:)   :: js                         
      !! REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr                               
      !!                                                                               
      !! INTEGER,      ALLOCATABLE, DIMENSION(:,:) :: jj_L, jjs_L, iis_L, neigh     
      !! INTEGER,      ALLOCATABLE, DIMENSION(:)   :: js_L                             
      !! REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_L                             
      !!                                                                               
      !! INTEGER,      ALLOCATABLE, DIMENSION(:)   :: mm, mms, sides, neighs           
      !! INTEGER                                   :: me, mes, np, nps                 
      !! INTEGER                                   :: np_L, nps_L                      
MODULE  prep_mesh_p1p2_sp

   USE Gauss_points_L;      USE Gauss_points
   ! che ovviamente si trovano dentro a:
   !
   !  gauss_points_2d_l.f90  e  gauss_points_2d_p_sp.f90
   !
   ! e contengono:
   !  Gauss_gen_L          Gauss_gen

   
   USE Dirichlet_Neumann    
   
   
   USE fem_2d_grid  !  fem_2d_openGrid
                    !  fem_2d_closeGrid
                    !  fem_2d_readGrid
                    !  fem_2d_getElemType
                    !  fem_2d_getGridSizes
                    

   IMPLICIT NONE

!------------------------------------------------------------------------------
!
!  jj (n_w,  me)   nodes of the  volume_elements
!
!  jjs(n_ws, mes)  nodes of the surface_elements in volume numbering
!
!  iis(n_ws, mes)  nodes of the surface_elements in surface numbering
!
!   js(nps)        nodes on surface
!
!   mm(me)         (possibly sub) set of elements for quadrature
!
!  mms(mes)        (possibly sub) set of surface_elements for surf_quadrature
!
!------------------------------------------------------------------------------

! Autore:  Andrea Marra 
!        a quanto pare col cazzo che e` stato Marra,
!        questo modulo e` opera di Quartapelle


PRIVATE

INTEGER,      ALLOCATABLE, DIMENSION(:,:), PUBLIC :: jj
!< *jj (n_w, me)* : indici dei nodi degli elementi nel dominio (volume).
INTEGER,      ALLOCATABLE, DIMENSION(:,:), PUBLIC :: jjs
!< *jjs (n_ws, mes)* : indici degli elementi di bordo, numerati accordingly a quelli di volume.
INTEGER,      ALLOCATABLE, DIMENSION(:,:), PUBLIC :: iis
!< *iis (n_ws, mes)* : indici degli elementi di bordo (surface) in surface numbering. Nota anche come matrice di connettivita` degli elementi di superficie, partendo dalla matrice jjs degli elementi di superficie numerati accordingly a quelli di volume, e dall'array js dei nodi di bordo accordingly alla numerazione dei nodi di volume.
INTEGER,      POINTER,     DIMENSION(:),   PUBLIC :: js
!< *js (nps)* : indici dei nodi di bordo in volume numbering.  
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: rr

INTEGER,      ALLOCATABLE, DIMENSION(:,:), PUBLIC :: jj_L, jjs_L, iis_L, neigh
INTEGER,      POINTER,     DIMENSION(:),   PUBLIC :: js_L ! boundary nodes --> volume numbering
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: rr_L

INTEGER,      ALLOCATABLE, DIMENSION(:),   PUBLIC :: sides, neighs
INTEGER,      ALLOCATABLE, DIMENSION(:),   PUBLIC :: mm, mms

INTEGER, PUBLIC :: me, mes, np, nps, np_L, nps_L


PUBLIC :: read_p1_gen_p2_sp



CONTAINS


!------------------------------------------------------------------------------
     !> \brief Usa le funzioni di fem_2d_grid.f90
     !!     per leggere la mesh; poi crea i nodi per l'interpolazione
     !!     P2.
SUBROUTINE  read_p1_gen_p2_sp (dir, fil)
  
   IMPLICIT NONE
   
   CHARACTER(LEN=64), INTENT(IN) :: dir, fil
   
   INTEGER :: nw_L, nws_L,  m, d_end, f_end, loops
  
   LOGICAL :: reverse_boundary = .TRUE. 
            !!< se gli elementi del contorno hanno i nodi 
            !!  ordinati in senso opposto al verso antiorario
            !!  del contorno esterno, che e' il verso 
            !!  positivo naturale per avere le normali 
            !!  uscenti dal dominio e rispettare quindi
            !!  il teorema di Gauss, occorre invertire l'ordine
            !!  dei nodi estremi di ciascun elemento.

   d_end = last_c_leng (64, dir)
   f_end = last_c_leng (64, fil)


   WRITE(*,*)   
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Loading mesh-file ', trim(fil), ' ...'


!  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

   CALL fem_2d_openGrid (dir(1:d_end), fil(1:f_end))
   
   CALL fem_2d_getElemType (nw_L, nws_L)
   
   CALL fem_2d_getGridSizes (me, mes, np_L, nps_L)
   
   IF ( nw_L /= n_w_L  .OR.  nws_L /= n_ws_L )  THEN
!      write(*,*) '    MOD : prep_mesh_p1p2_sp'
!      write(*,*) '    SUB : read_p1_gen_p2_sp'
!      write(*,*) '    CAL : newton'
      WRITE(*,*) '    grid data incompatible with finite-element data'
!      write(*,*)
      STOP
   ENDIF
   
   ALLOCATE (jj_L (n_w_L,  me),   neigh(n_w_L,  me))
   ALLOCATE (jjs_L(n_ws_L, mes),  iis_L(n_ws_L, mes)) 
   ALLOCATE (sides(mes), neighs(mes))
   ALLOCATE (js_L(nps_L))
   ALLOCATE (rr_L(k_d, np_L))
   
   CALL fem_2d_readGrid (rr_L, jj_L, jjs_L, iis_L, js_L, sides, neigh, neighs)
   ! sides : side index of the surface element. array di interi di dimensione
   !         (mes) che contiene per ogni elemento di superficie l'indice che
   !         indica a quale segmento di bordo appartiene quell'elemento
  
   CALL fem_2d_closeGrid
!   write(*,*) '    MOD : prep_mesh_p1p2_sp'
!   write(*,*) '    SUB : read_p1_gen_p2_sp'
!   write(*,*) '    CAL : newton'
   WRITE(*,*) '    End of p1 grid reading '
!   write(*,*)
   
   
   
   !  END OF GRID READING -----------------------------------------------------
   
   
   !------ Reverse the P1 Grid if necessary -----------------------------------
                          
   IF (reverse_boundary)  & 
   CALL invert_nodes_of_bound_elements (jjs_L, iis_L)
   
   !------ Creating P2 Grid starting from P1 Grid -----------------------------
   
   
   np = np_L + number_of_edges(neigh) ! Number of parabolic nodes
   
   
   !  Use of Euler relation for counting the number of internal loops
   
   loops = np - np_L - (me + np_L - 1)
!   write(*,*) '    MOD : prep_mesh_p1p2_sp'
!   write(*,*) '    SUB : read_p1_gen_p2_sp' 
!   write(*,*) '    CAL : newton'
   WRITE(*,*) '    number of connected internal boundaries:' 
   WRITE(*,*) '    loops = ', loops
!   write(*,*)
   
   
   ALLOCATE (jj(n_w, me),  jjs(n_ws, mes), iis(n_ws, mes)) 
   
   ALLOCATE (rr(k_d, np))
   
   CALL gen_p2_from_p1 (jj_L, jjs_L, rr_L, neigh, neighs,  jj, jjs, rr)   
  
   CALL Dirichlet_nodes_gen (jjs, SPREAD(1,1,mes), SPREAD(.TRUE.,1,1),  js)
  
   CALL surf_nodes_i (jjs, js,  iis)
   
!   write(*,*) '    MOD : prep_mesh_p1p2_sp'
!   write(*,*) '    SUB : read_p1_gen_p2_sp'
!   write(*,*) '    CAL : newton'
   WRITE(*,*) '    End of p2 grid creating'
!   write(*,*)
  
   
   nps = SIZE(js)
   
   nps_L = SIZE(js_L)

    
   ALLOCATE (mm(me), mms(mes)) 

   DO m = 1, me;    mm(m) = m;  ENDDO  ! elements for volume integration
   DO m = 1, mes;  mms(m) = m;  ENDDO  ! elements for surface integration

   CALL Gauss_gen_L(np_L, me, nps_L, mes, jj_L, jjs_L, rr_L)
   ! in file : gauss_points_2d_l.f90

   CALL Gauss_gen  (np,   me, nps,   mes, jj,   jjs,   rr)
   ! in file : gauss_points_2d_p_sp.f90
! Parabolic-Linear element with SubParametric transformation
! Storage of the derivative at the reference element


   CLOSE(30)

   WRITE(*,*)
   WRITE(*,*) '--> mesh data:'
   WRITE(*,*) '    me   = ', me,        ', Number of domain   (volume)  elements'
   WRITE(*,*) '    mes  = ', mes,       ', Number of boundary (surface) elements'
   WRITE(*,*) '    np_L = ', np_L,      ', Number of nodes in the domain'
   WRITE(*,*) '    nps  = ', nps,       ', Number of nodes on the boundary =SIZE(js)' 
   WRITE(*,*) '    np   = ', np,        ', Number of parabolic nodes'
   WRITE(*,*) '    Nx   = ', 3*np+np_L, ', Number of unknowns'

END SUBROUTINE  read_p1_gen_p2_sp


!------------------------------------------------------------------------------


FUNCTION  last_c_leng (len_str, string) RESULT(leng)

   IMPLICIT NONE

   INTEGER,                 INTENT(IN) :: len_str
   CHARACTER (LEN=len_str), INTENT(IN) :: string
   INTEGER :: leng

   INTEGER :: i

   leng = len_str

   DO i = 1, len_str
      
      IF (string(i:i) == ' ') THEN
         leng = i - 1;  EXIT
      ENDIF
      
   ENDDO

END FUNCTION  last_c_leng

   
!------------------------------------------------------------------------------


SUBROUTINE  surf_nodes_i (jjs, js,  iis)

!  generation of the surface element connectivity matrix  iis
!  based on the surface node numbering, starting from the
!  connectivity matrix  jjs  of the surface elements based on
!  the volume node numbering, and from the array  js  of the
!  boundary nodes according to the volume node numbering

   IMPLICIT NONE

   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjs
   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
   INTEGER, DIMENSION(:,:), INTENT(OUT) :: iis

   INTEGER :: ms, ls, j, i

   DO ms = 1, SIZE(jjs,2)
     
      DO ls = 1, SIZE(jjs,1)
   
         j = jjs(ls,ms)
   
         DO i = 1, SIZE(js)
            IF ( j == js(i) )  iis(ls,ms) = i
         ENDDO
   
      ENDDO
   
   ENDDO

END SUBROUTINE  surf_nodes_i


!-------------------------------------------------------------------------------


SUBROUTINE  gen_p2_from_p1 (jj_L, jjs_L, rr_L, neigh, neighs, jj, jjs, rr)
!^^^^^^^^^^^^^^^^^^^^^^~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!  jj_L(:, :)  nodes of the  volume_elements of the input p1 grid        
! jjs_L(:, :)  nodes of the surface_elements of the input p1 grid 
!  rr_L(:, :)  Cartesian coordinates of the nodes of the input p1 grid  

!  jj(:, :)  nodes of the  volume_elements of the output p2 grid   
! jjs(:, :)  nodes of the surface_elements of the output p2 grid   
!  rr(:, :)  Cartesian coordinates of the nodes of the output p2 grid 

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_L, jjs_L, neigh
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: rr_L
   INTEGER,      DIMENSION(:),   INTENT(IN)  :: neighs
   INTEGER,      DIMENSION(:,:), INTENT(OUT) :: jj, jjs
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: rr

   LOGICAL, DIMENSION(SIZE(jj_L,2)) :: virgin
 
   INTEGER, DIMENSION(SIZE(jj_L,1), SIZE(jj_L,2)) :: m_op
   INTEGER, DIMENSION(SIZE(jj_L,1), SIZE(jj_L,2)) :: j_mid
   INTEGER, DIMENSION(SIZE(jjs_L,2))              :: js_mid
  
   INTEGER :: np_L, me, nw, nws, kd,         & 
              n, n1, n2, ms, k, k1, k2, kk,  &
              m_op_k, i, m, mm 
  
   REAL(KIND=8), DIMENSION(SIZE(rr_L,1)) :: r_mid

   nw  = SIZE( jj_L, 1)  ! nodes in each volume element
   nws = SIZE(jjs_L, 1)  ! nodes in each surface element
   kd  = SIZE( rr_L, 1)  ! spatial dimensions

   m_op = neigh


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!  GENERATION OF THE P2 GRID

   np_L = SIZE(rr_L, 2)
   me   = SIZE(jj_L, 2)

   rr(:, 1:np_L) = rr_L

   jj(1:nw, :) = jj_L

   virgin = .true.

   n = np_L
  
   DO m = 1, me ! loop on the elements
     
      DO k = 1, nw ! loop on the nodes (sides) of the element

         m_op_k = m_op(k, m)
         k1 = MODULO(k,   nw) + 1;  n1 = jj_L(k1, m)
         k2 = MODULO(k+1, nw) + 1;  n2 = jj_L(k2, m)
         r_mid = (rr_L(:, n1) + rr_L(:, n2))/2

         IF (m_op_k == 0) THEN  !  the side is on the boundary

            n = n + 1
            j_mid(k, m) = n
            rr(:, n) = r_mid

            ! surface elements of the p2 grid are defined later

         ELSE  !  the side is internal
                  
            IF ( virgin(m_op_k) ) THEN  !  the side is new

               n = n + 1
               j_mid(k, m) = n
               rr(:, n) = r_mid

            ELSE  !  the side has been already considered

               mm = m_op_k
               DO i = 1, nw
                  IF (m_op(i, mm) == m) kk = i
               ENDDO
               j_mid(k, m) = j_mid(kk, mm)

            ENDIF

         ENDIF

      ENDDO

      virgin(m) = .FALSE.

   ENDDO

!  connectivity matrix of the p2 grid

   jj(nw+1 : SIZE(jj,1), :) = j_mid


!  connectivity matrix of the surface elements of the p2 grid

   DO ms = 1, SIZE(jjs_L, 2);  mm = neighs(ms)
      DO i = 1,nw
         IF (m_op(i, mm) == 0) kk = i
      ENDDO
      js_mid(ms) = j_mid(kk, mm)
   ENDDO

   jjs(1:SIZE(jjs_L,1), :) = jjs_L

   jjs(3, :) = js_mid


END SUBROUTINE  gen_p2_from_p1

!------------------------------------------------------------------------

SUBROUTINE  invert_nodes_of_bound_elements (jjs_L, iis_L)

   IMPLICIT NONE

   INTEGER, DIMENSION(:,:), INTENT(INOUT) :: jjs_L, iis_L

   INTEGER :: ms, j1, i1 
   
!  reorder the two extreme nodes of the boundary elements for 2D 
!  problems to assure a positive orientation.
!  WARNING: internal nodes are not reordered

   DO ms = 1, SIZE(jjs_L, 2)

               j1 = jjs_L(1,ms)
      jjs_L(1,ms) = jjs_L(2,ms)
      jjs_L(2,ms) = j1        
      
               i1 = iis_L(1,ms)
      iis_L(1,ms) = iis_L(2,ms)
      iis_L(2,ms) = i1
      
   ENDDO


END SUBROUTINE  invert_nodes_of_bound_elements


!-----------------------------------------------------------------------


FUNCTION  number_of_edges(neigh) RESULT(n_e)
   
   ! determina il numero di mid-side nodes di un reticolo 
   ! di triangoli fornito tramite la matrice di connettivita' 
   ! neigh degli elementi opposti ai nodi 
   ! se non esiste elemento opposto al nodo  neigh(k,m) = 0
  
   INTEGER, DIMENSION(:,:), INTENT(IN) :: neigh
   INTEGER :: n_e
    
   LOGICAL, DIMENSION(0 : SIZE(neigh,2)) :: virgin
   INTEGER :: m
   
   virgin = .TRUE.
   virgin(0) = .TRUE.
   
   n_e = 0 
   
   DO m = 1, SIZE(neigh,2)
   
      n_e = n_e + COUNT(virgin(neigh(:,m)))
     
      virgin(m) = .FALSE. 
   
   ENDDO

   
END FUNCTION  number_of_edges


!------- Eulero's relations for counting the number of internal loops
!------- in 2 and 3 dimensions
!
!   SELECT CASE (k_d)  !  Viva Eulero
!
!
!      CASE (2)
!       
!         np     =  np_L  +  (me + np_L - 1)  +  loops  !  in 2D
!        
!         nps_L  =  mes
!        
!         nps    =  2*mes
!
!
!      CASE (3)
!
!         k_faces  =  2*me  +  mes/2
!        
!         l_edges  =  k_faces  +  np_L - 1  -  me  ! =  me  +  mes/2  +  np - 1
!        
!         np       =  np_L  +  l_edges             !  in 3D
!
!
!         nps_L  =  mes/2  +  2
!        
!         nps    =  nps_L  +  mes  +  nps_L  -  2
!
!
!   END SELECT



END MODULE  prep_mesh_p1p2_sp
