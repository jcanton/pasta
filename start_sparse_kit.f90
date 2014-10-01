MODULE  start_sparse_kit
   
   USE sparse_matrix_profiles
   
   IMPLICIT NONE
   
   SAVE
   
   PRIVATE :: sort_pick


CONTAINS


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE  start_matrix_2d_p1 (np, jj, js,  AA)
!+++++++++++++++++++++++++++++++++++++++++++++++

!  Construction of the CSR data structure for matrix  AA
!  Matrix AA can be nonsymmetric.  The column indices of
!  the nonzero matrix elements are ordered increasingly

!  LINEAR TRIANGULAR ELEMENTS 

!  The diagonal element in each row is assumed to be different from zero

!   np  --->  number of nodes of the grid  ==  order of matrix AA
!
!   jj  --->  standard element-to-node topology of the fem grid
!
!   js  --->  array of the node numbers of all boundary nodes (global numbering)
!
!   AA  --->  matrix in CSR format (DERIVED TYPE)

   IMPLICIT NONE

   INTEGER,                 INTENT(IN)  :: np
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
   TYPE(CSR_MUMPS_Matrix),  INTENT(OUT) :: AA

   INTEGER, DIMENSION(np) :: nz_el_in_each_row
   INTEGER :: i, j, m, NA_nz, v, w, k, l

   ALLOCATE (AA % i(np+1))


   nz_el_in_each_row = 1

   DO m = 1, SIZE(jj,2)

      nz_el_in_each_row(jj(:,m))  =  nz_el_in_each_row(jj(:,m))  +  1

   ENDDO

   ! correction for the boundary nodes

   nz_el_in_each_row(js)  =  nz_el_in_each_row(js)  +  1


   ! definition of the first (nz) element of each row

   AA % i(1) = 1

   DO i = 1, np

      AA % i(i+1)  =  AA % i(i)  +  nz_el_in_each_row(i)

   ENDDO

   NA_nz  =  AA % i(np+1)  -  1

   ALLOCATE (AA % j(NA_nz), AA % e(NA_nz))


   !  column indices of all nonzero elements

   AA % j = 0

   DO m = 1, SIZE(jj,2)

      DO v = 1, SIZE(jj,1);   i = jj(v,m)

         k  =  AA % i(i)           ! first element of the row
         l  =  AA % i(i+1)  -  1   ! last  element of the row

         DO w = 1, SIZE(jj,1);   j = jj(w,m)

            IF (ANY(AA % j(k:l) == j))  CYCLE

            ! forward shift and stack of the new element
            AA % j(l : k+1 : -1)  =  AA % j(l-1 : k : -1)
            AA % j(k)  =  j

         ENDDO

      ENDDO

   ENDDO


   ! sort in scending order of the column indices of each row

   DO i = 1, np

      k  =  AA % i(i)           ! first element of the row
      l  =  AA % i(i+1)  -  1   ! last  element of the row

      CALL  sort_pick ( AA % j(k : l) )

   ENDDO

   !======================================================================
   !!!!!!!!!!!! MODIFIED !!!!!!!!!!!!!!!!!
   
   ALLOCATE ( AA % i_mumps (SIZE(AA%j)) )
   
   DO i = 1, SIZE(AA%i)-1

      DO j = AA%i(i), AA%i(i+1)-1

         AA%i_mumps(j) = i

      ENDDO

   ENDDO
   
   !======================================================================

   WRITE (*,*) '    CSR matrix structure for the P1 grid completed'
   

END SUBROUTINE  start_matrix_2d_p1


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE  start_matrix_2d_p2 (np, jj, js,  AA)
!+++++++++++++++++++++++++++++++++++++++++++++++

!  Construction of the CSR data structure for matrix  AA
!  Matrix AA can be nonsymmetric.  The column indices of
!  the nonzero matrix elements are ordered increasingly

!  PARABOLIC TRIANGULAR ELEMENTS 

!  The diagonal element in each row is assumed to be different from zero

!   np  --->  number of nodes of the grid  ==  order of matrix AA
!
!   jj  --->  standard element-to-node topology of the fem grid
!
!   js  --->  array of the node numbers of all boundary nodes (global numbering)
!
!   AA  --->  matrix in CSR format (DERIVED TYPE)
   
   IMPLICIT NONE
  
   INTEGER,                 INTENT(IN)  :: np
   INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
   INTEGER, DIMENSION(:),   INTENT(IN)  :: js
   TYPE(CSR_MUMPS_Matrix),  INTENT(OUT) :: AA
   
   INTEGER, DIMENSION(np) :: nz_el_in_each_row
   INTEGER :: i, j, m, NA_nz, v, w, k, l, ii, c

   ALLOCATE (AA % i(np+1))

 
   ! preliminary determination of the boundary nodes that
   ! are mid-side nodes to have the right correction for
   ! all boundary nodes, including the mid-side ones. 
 
   nz_el_in_each_row = 0

   DO m = 1, SIZE(jj,2)

      nz_el_in_each_row(jj(:,m))  =  nz_el_in_each_row(jj(:,m))  +  1

   ENDDO

   WHERE (nz_el_in_each_row == 1)
      nz_el_in_each_row = - 1
   ELSEWHERE
      nz_el_in_each_row = 0
   END WHERE

   ! end of the preliminary countercorrection phase


   nz_el_in_each_row = nz_el_in_each_row + 1

   DO m = 1, SIZE(jj,2)

      nz_el_in_each_row(jj(1:3,m))  =  nz_el_in_each_row(jj(1:3,m))  +  3
      
      nz_el_in_each_row(jj(4:6,m))  =  nz_el_in_each_row(jj(4:6,m))  +  4
   
   ENDDO

   ! correction for all boundary nodes

   nz_el_in_each_row(js)  =  nz_el_in_each_row(js)  +  2


   ! definition of the first (nz) element of each row

   AA % i(1) = 1

   DO i = 1, np

      AA % i(i+1)  =  AA % i(i)  +  nz_el_in_each_row(i)

   ENDDO

   NA_nz  =  AA % i(np+1)  -  1

   ALLOCATE (AA % j(NA_nz), AA % e(NA_nz))

!+++ Edited by Pier   
   ALLOCATE ( AA % i_mumps (NA_nz) )
   
   c = 1
   
   DO i = 1, np
   
      DO ii = AA%i(i), AA%i(i+1) -1
      
         AA%i_mumps(c) = i
         
         c = c + 1
      
      END DO
      
   END DO
!+++
 
   !  column indices of all nonzero elements

   AA % j = 0

   DO m = 1, SIZE(jj,2)

      DO v = 1, SIZE(jj,1);   i = jj(v,m)

         k  =  AA % i(i)           ! first element of the row
         l  =  AA % i(i+1)  -  1   ! last  element of the row

         DO w = 1, SIZE(jj,1);   j = jj(w,m)

            IF (ANY(AA % j(k:l) == j))  CYCLE

            ! forward shift and stack of the new element
            AA % j(l : k+1 : -1)  =  AA % j(l-1 : k : -1)
            AA % j(k)  =  j

         ENDDO

      ENDDO

   ENDDO


   ! sort in scending order of the column indices of each row

   DO i = 1, np

      k  =  AA % i(i)           ! first element of the row
      l  =  AA % i(i+1)  -  1   ! last  element of the row

      CALL  sort_pick ( AA % j(k : l) )

   ENDDO

   
   WRITE (*,*) '    CSR matrix structure for the P2 grid completed'


END SUBROUTINE  start_matrix_2d_p2


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE  start_coupled_system_axisym (np, np_L, jj, js,  CC)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Cylindrical coordinate for Axisymmetric problems with swirl
 
! generation of the CSR structure of the matrix CC
! of the coupled velocity-pressure system with 
! interpolations of different order (p2-p1) from that
! of the matrix AA of the uncoupled system with the
! parabolic interpolation (p2) of the velocity component. 
 
! The CSR structure of the coupling rectangular part
! is based on the assumption that all the linear nodes
! of the grid precede all the remaining midside-nodes
! of the parabolic interpolation
    
   IMPLICIT NONE
   ! input variables
   INTEGER,                 INTENT(IN) :: np, np_L
   INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
   INTEGER, DIMENSION(:),   INTENT(IN) :: js
   ! output variables
   TYPE(CSR_MUMPS_Matrix) :: CC
  
   ! local variables
   TYPE(CSR_MUMPS_Matrix)       :: AA

   INTEGER, DIMENSION(np)       :: nz_RR   
   INTEGER, DIMENSION(np + 1)   :: RR_i 
   
   INTEGER, DIMENSION(np_L)     :: nz_RT 
   INTEGER, DIMENSION(np_L + 1) :: RT_i
   
   INTEGER, DIMENSION(:),   ALLOCATABLE :: RR_j, RT_j        
    
!   INTEGER, DIMENSION(:,:), ALLOCATABLE :: RR
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 2 modifiche fatte da Auteri

!  INTEGER, DIMENSION(np) :: column
   INTEGER, DIMENSION(np_L) :: col_num

! e` anche stato eliminato l'IF (SAVE_MEMORY)
! perche' porlo uguale a .FALSE. comporta un uso
! di memoria eccessivo.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER :: NR, size_CC, k, i, j, iA, ii, p, n, c


   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++'
   WRITE(*,*) '--> Creating Jacobian matrix sparse structure'


!+++
   CALL start_matrix_2d_p2 (np, jj, js,  AA)
!+++


   ! preliminar counting of the coupling rectangular part  

   DO i = 1, np
   
      k = 0
  
      DO p = AA%i(i),  AA%i(i+1) - 1
  
         IF (AA%j(p) <= np_L)  k = k + 1
  
      ENDDO
  
      nz_RR(i) = k
  
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! modifica 1 di Auteri:
    
   nz_RT = 0

   DO n = 1, SIZE(AA%j)
      
      j = AA%j(n)
      IF (j <= np_L) nz_RT(j) = nz_RT(j) + 1

   END DO
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   ! Consistency check
  
   NR = SUM(nz_RR)
  
   ALLOCATE (RR_j(NR),  RT_j(NR)) 
  
   WRITE (*,*) '--> Consistency check'
   WRITE (*,*) '    SUM(nz_RR) = ', SUM(nz_RR)
   WRITE (*,*) '    SUM(nz_RT) = ', SUM(nz_RT)
   WRITE (*,*) '    tot from the first np_L rows of AA = ', AA%i(np_L+1) - 1
   WRITE (*,*)

   IF (NR /= AA%i(np_L+1) - 1) THEN
   
      WRITE (*,*) 'incoherence between NR and AA%i(np_L+1) - 1'
      WRITE (*,*) 'in start_coupled_system.  STOP'
      
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr) 
   
   ENDIF

   ! ALLOCATION OF THE STRUCTURE FOR CC OF THE COUPLED SYSTEM   

   size_CC = 9 * SIZE(AA%j)  +  6 * NR  
  
   ! bottom-right last diagonal element    
   size_CC = size_CC  +  1 
   
   WRITE (*,*) '    size_AA = ', SIZE(AA%j)
   WRITE (*,*) '    size_CC = ', size_CC
  
   ALLOCATE( CC%i(3*np + np_L + 1) )
   ALLOCATE( CC%i_mumps(size_CC), CC%j(size_CC), CC%e(size_CC) )

  
   ! DEFINE THE COLUMN INDICES OF THE RECTANGULAR ARRAY RR
  
   ! first, the pointer vector RR_i 
  
   RR_i(1) = 1
   
   DO i = 1, np
  
      RR_i(i+1) = RR_i(i)  +  nz_RR(i)
   
   ENDDO
      
   ! then, the column indices of RR   
      
   k = 0   
      
   DO i = 1, np
      
      DO p = AA%i(i),  AA%i(i+1) - 1
      
         IF (AA%j(p) <= np_L) THEN 
      
            k = k + 1 
      
            RR_j(k) = AA%j(p)  
      
         ENDIF 
      
      ENDDO
      
   ENDDO   
     
   ! 1 <= RR_j <= np_L  


   ! CONSTRUCTION OF THE TRANSPOSED RECTANGULAR ARRAY RT 
 
   ! first, the pointer vector RT_i  
 
   RT_i(1) = 1
  
   DO i = 1, np_L
   
      RT_i(i+1) = RT_i(i)  +  nz_RT(i)  
   
   ENDDO

   ! then, the column index vector RT_j  ! More difficult
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! modifica 2 di Auteri:

   col_num = 0
   
   DO ii = 1, np

      DO p = RR_i(ii), RR_i(ii+1)-1

         j = RR_j(p)

         col_num(j) = col_num(j) + 1

         RT_j(RT_i(j) + col_num(j) - 1) = ii

      END DO

   END DO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
   ! SET THE NUMBER OF nz-ELEMENTS OF THE ROWS 
   ! OF THE MATRIX CC OF THE COUPLED SYSTEM 
 
   CC%i(1) = 1  
  
   ! first block-row 
    
   DO i = 1, np 
 
      CC%i(i+1) = CC%i(i)  +  3 * (AA%i(i+1) - AA%i(i))  +  nz_RR(i)

   ENDDO

  
   ! second block-row 
  
   DO i = np + 1,  2*np
   
      iA = i - np 

      CC%i(i+1) = CC%i(i)  +  3 * (AA%i(iA+1) - AA%i(iA))  +  nz_RR(iA)
   
   ENDDO
   
   
   ! third block-row 
  
   DO i = 2*np + 1,  3*np
   
      iA = i - 2*np 

      CC%i(i+1) = CC%i(i)  +  3 * (AA%i(iA+1) - AA%i(iA))  +  nz_RR(iA)
   
   ENDDO


   ! bottom rectangular block-row: USE TRANSPOSED ARRAY RT 

   DO i = 3*np + 1,  3*np + np_L
    
      iA = i - 3*np

      CC%i(i+1) = CC%i(i)  +  3 * nz_RT(iA) ! == (RT_i(iA+1) - RT_i(iA))   

   ENDDO
   
   ! bottom-right last diagonal element 
   CC%i(3*np + np_L + 1)  =  CC%i(3*np + np_L + 1)  +  1 


   ! DEFINITION OF THE COLUMN INDICES
   !=================================
   
   ! first block-row 

   k = 0

   DO i = 1, np 
    
      DO p = AA%i(i),  AA%i(i+1) - 1          
         k = k + 1
         CC%j(k) = AA%j(p)
      ENDDO
      
      DO p = AA%i(i),  AA%i(i+1) - 1    
         k = k + 1
         CC%j(k) = AA%j(p) + np
      ENDDO
     
      DO p = AA%i(i),  AA%i(i+1) - 1    
         k = k + 1
         CC%j(k) = AA%j(p) + 2*np
      ENDDO
     
      DO p = AA%i(i),  AA%i(i+1) - 1 
     
         IF (AA%j(p) <= np_L) THEN
            k = k + 1
            CC%j(k) = AA%j(p) + 3*np
         ENDIF 
      
      ENDDO
 
   ENDDO
 
    
   ! second block-row 
   
   DO i = np + 1,  2*np
   
      iA = i - np
    
      DO p = AA%i(iA),  AA%i(iA+1) - 1              
         k = k + 1 
         CC%j(k) = AA%j(p)
      ENDDO
      
      DO p = AA%i(iA),  AA%i(iA+1) - 1    
         k = k + 1 
         CC%j(k) = AA%j(p) + np
      ENDDO
     
      DO p = AA%i(iA),  AA%i(iA+1) - 1    
         k = k + 1 
         CC%j(k) = AA%j(p) + 2*np
      ENDDO
 
      DO p = AA%i(iA),  AA%i(iA+1) - 1  
      
         IF (AA%j(p) <= np_L) THEN
            k = k + 1 
            CC%j(k) = AA%j(p) + 3*np
         ENDIF
      
      ENDDO
 
   ENDDO
 
  
   ! third block-row 
   
   DO i = 2*np + 1,  3*np
   
      iA = i - 2*np
    
      DO p = AA%i(iA),  AA%i(iA+1) - 1              
         k = k + 1 
         CC%j(k) = AA%j(p)
      ENDDO
      
      DO p = AA%i(iA),  AA%i(iA+1) - 1    
         k = k + 1 
         CC%j(k) = AA%j(p) + np
      ENDDO
  
      DO p = AA%i(iA),  AA%i(iA+1) - 1    
         k = k + 1 
         CC%j(k) = AA%j(p) + 2*np
      ENDDO
 
      DO p = AA%i(iA),  AA%i(iA+1) - 1  
      
         IF (AA%j(p) <= np_L) THEN
            k = k + 1 
            CC%j(k) = AA%j(p) + 3*np
         ENDIF
      
      ENDDO
 
   ENDDO
 

   ! bottom rectangular block-row

   DO i = 3*np + 1,  3*np + np_L
   
      iA = i - 3*np
    
      DO p = RT_i(iA),  RT_i(iA+1) - 1    
         k = k + 1
         CC%j(k) = RT_j(p) 
      ENDDO
      
      DO p = RT_i(iA),  RT_i(iA+1) - 1    
         k = k + 1   
         CC%j(k) = RT_j(p) + np
      ENDDO
  
      DO p = RT_i(iA),  RT_i(iA+1) - 1    
         k = k + 1   
         CC%j(k) = RT_j(p) + 2*np
      ENDDO
 
   ENDDO
  
   ! bottom-right last diagonal element
   k = k + 1
   CC%j(k) = 3*np + np_L


!+++ Edited by Pier
   c = 1   ! COUNTER !
   
   DO i = 1, 3*np + np_L
   
      DO ii = CC%i(i), CC%i(i+1) - 1
      
         CC%i_mumps(c) = i
         
         c = c + 1
         
      ENDDO
      
   END DO   
!+++

   WRITE (*,*) '    CSR  start_coupled_system  completed.'

END SUBROUTINE  start_coupled_system_axisym

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE  sort_pick (arr)

!  sorts an integer array  arr  into ascending numerical order, by stright
!  insertion.   arr  is replaced on output by its sorted rearrangement.
!  Not parallelizable.  Very inefficient for  SIZE(arr) > 20.

!  Adapted from:  Numerical Recipes in Fortran 90,  p. 1167

   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(INOUT) :: arr

   INTEGER :: aj, i, j

   DO j = 2, SIZE(arr);   aj = arr(j)

      DO i = j - 1, 1, -1
         IF (arr(i) <= aj) EXIT
         arr(i+1) = arr(i)
      ENDDO

      arr(i+1) = aj

   ENDDO

END SUBROUTINE  sort_pick



END MODULE  start_sparse_kit
