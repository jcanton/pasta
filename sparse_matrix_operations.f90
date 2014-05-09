MODULE sparse_matrix_operations
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
   IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
! subroutines for real matrices
!------------------------------------------------------------------------------
!
!
SUBROUTINE dAlinB_s (alpha, a, beta, b,  c)
!
!-----------------------------------------------------------------------
!         C  =  alpha A  +  beta B
!----------------------------------------------------------------------- 
! linearly combine two matrices
! Matrices are stored in compressed sparse row storage.
! The three matrices need to have the same sparsity pattern
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, b  = input real matrices in compressed sparse row format
!
! alpha,
! beta  = real coefficients of the linear combination
!
!
! on return:
!-----------
! c     = real matrix, containing the sum C  =  alpha A  +  beta B
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8),               INTENT(IN) :: alpha, beta
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a,     b
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: c
   ! local variables
   INTEGER :: i

!$OMP PARALLEL DO
   DO i = 1, SIZE(a)

      c(i) = alpha*a(i) + beta*b(i)

   END DO
!$OMP END PARALLEL DO

END SUBROUTINE dAlinB_s

!------------------------------------------------------------------------------

SUBROUTINE dAtimx (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array, containing the product y = A*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: y
   ! local variables
   REAL(KIND=8) :: yi
   INTEGER      :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! not necessary to check x-y consistency because A has to be square


   ! compute y <-- A*x
!$OMP PARALLEL DO
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x

      yi = 0d0

      DO p = ia(i), ia(i+1)-1

         yi = yi + a(p)*x(ja(p))

      END DO

      ! store result in y(i)

      y(i) = yi

   END DO
!$OMP END PARALLEL DO


END SUBROUTINE dAtimx

!------------------------------------------------------------------------------

SUBROUTINE dAtimx_T (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         transpose(A) times a vector
!----------------------------------------------------------------------- 
! multiplies the transpose of a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array, containing the product y = A'*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   REAL(KIND=8), DIMENSION(:)             :: y
   ! local variables
   INTEGER :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx_T Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> dAtimx_T Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! not necessary to check x-y consistency because A has to be square

   ! compute y <-- A'*x
   y = 0

!  $OMP PARALLEL DO &    ! not working yet with big matrices,
!  $OMP PRIVATE(p)  &    ! don't know why
!  $OMP REDUCTION (+: y)
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x
      DO p = ia(i), ia(i+1)-1

         y(ja(p)) = y(ja(p)) + a(p)*x(i)

      END DO

   END DO
!  $OMP END PARALLEL DO


END SUBROUTINE dAtimx_T

!------------------------------------------------------------------------------

SUBROUTINE dEssM (a, ja, ia, n, b, jb, ib, i_mumpsb)
!
!-----------------------------------------------------------------------
!         Extract square sub-Matrix  nxn
!----------------------------------------------------------------------- 
! Extract a square submatrix B from matrix A
! Matrix A and B are stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! b, jb,
!    ib = output matrix in compressed sparse row format.
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,      DIMENSION(:), INTENT(IN) :: ja, ia
   INTEGER,                    INTENT(IN) :: n
   ! output variables
   REAL(KIND=8), DIMENSION(:), POINTER           :: b
   INTEGER,      DIMENSION(:), POINTER           :: jb, ib
   INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb
   ! local variables
   INTEGER      :: number_of_rows, i, p, c


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   !
   IF ( n > number_of_rows ) THEN
      WRITE(*,*) '--> dEssM Error: n > rows(A)' 
      WRITE(*,*) '    n = ', n, ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF


   ! count elements
   !
   c = 0
   DO i = 1, n

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

         ENDIF

      ENDDO

   ENDDO

   ALLOCATE( ib(n+1) )
   ALLOCATE( jb(c), b(c) )

   ! fill matrix
   !
   c = 0
   DO i = 1, n

      ib(i) = c + 1

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

            jb(c) = ja(p)
             b(c) =  a(p)

         ENDIF

      ENDDO

   ENDDO
   ib(n+1) = c + 1


   IF (PRESENT(i_mumpsb)) THEN

      ALLOCATE( i_mumpsb(SIZE(jb)) )

      DO i = 1, SIZE(ib) - 1
      
         DO p = ib(i), ib(i+1) - 1

            i_mumpsb(p) = i

         ENDDO

      END DO

   ENDIF

END SUBROUTINE dEssM
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutines for complex matrices
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
SUBROUTINE zAlinB_s (alpha, a, beta, b,  c)
!
!-----------------------------------------------------------------------
!         C  =  alpha A  +  beta B
!----------------------------------------------------------------------- 
! linearly combine two matrices
! Matrices are stored in compressed sparse row storage.
! The three matrices need to have the same sparsity pattern
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 25/8/2013
!
! on entry:
!----------
! a, b  = input complex matrices in compressed sparse row format
!
! alpha,
! beta  = complex coefficients of the linear combination
!
!
! on return:
!-----------
! c     = complex matrix, containing the sum C  =  alpha A  +  beta B
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8),               INTENT(IN) :: alpha, beta
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a,     b
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: c
   ! local variables
   INTEGER :: i

!$OMP PARALLEL DO
   DO i = 1, SIZE(a)

      c(i) = alpha*a(i) + beta*b(i)

   END DO
!$OMP END PARALLEL DO

END SUBROUTINE zAlinB_s

!------------------------------------------------------------------------------

SUBROUTINE zAtimx (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 22/3/2013
!
! on entry:
!----------
! x     = complex array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = complex array, containing the product y = A*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: y
   ! local variables
   COMPLEX(KIND=8) :: yi
   INTEGER         :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! not necessary to check x-y consistency because A has to be square

   ! compute y <-- A*x
!$OMP PARALLEL DO
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x

      yi = CMPLX(0d0, 0d0, KIND=8)

      DO p = ia(i), ia(i+1)-1

         yi = yi + a(p)*x(ja(p))

      END DO

      ! store result in y(i)

      y(i) = yi

   END DO
!$OMP END PARALLEL DO


END SUBROUTINE zAtimx

!------------------------------------------------------------------------------

SUBROUTINE zAtimx_T (y, a, ja, ia, x)
!
!-----------------------------------------------------------------------
!         transpose(A) times a vector
!----------------------------------------------------------------------- 
! multiplies the transpose of a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 22/3/2013
!
! on entry:
!----------
! x     = complex array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = complex array, containing the product y = A'*x
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: x
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:)             :: y
   ! local variables
   INTEGER :: number_of_rows, i, p


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   ! A - x
   IF ( SIZE(x) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx_T Error: wrong A - x size' 
      WRITE(*,*) '    SIZE(x) = ', SIZE(x), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! A - y
   IF ( SIZE(y) /= number_of_rows ) THEN
      WRITE(*,*) '--> zAtimx_T Error: wrong A - y size' 
      WRITE(*,*) '    SIZE(y) = ', SIZE(y), ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF
   ! not necessary to check x-y consistency because A has to be square

   ! compute y <-- A'*x
   y = CMPLX(0d0, 0d0, KIND=8)
   
!  $OMP PARALLEL DO &    ! not working yet with big matrices,
!  $OMP PRIVATE(p)  &    ! don't know why
!  $OMP REDUCTION (+: y)
   DO i = 1, number_of_rows

   ! compute the inner product of row i with vector x
      DO p = ia(i), ia(i+1)-1

         y(ja(p)) = y(ja(p)) + a(p)*x(i)

      END DO

   END DO
!  $OMP END PARALLEL DO


END SUBROUTINE zAtimx_T

!------------------------------------------------------------------------------

SUBROUTINE zEssM (a, ja, ia, n, b, jb, ib, i_mumpsb)
!
!-----------------------------------------------------------------------
!         Extract square sub-Matrix  nxn
!----------------------------------------------------------------------- 
! Extract a square submatrix B from matrix A
! Matrix A and B are stored in compressed sparse row storage.
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 16/03/2014
!
! on entry:
!----------
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! b, jb,
!    ib = output matrix in compressed sparse row format.
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   ! input variables
   COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: a
   INTEGER,         DIMENSION(:), INTENT(IN) :: ja, ia
   INTEGER,                       INTENT(IN) :: n
   ! output variables
   COMPLEX(KIND=8), DIMENSION(:), POINTER        :: b
   INTEGER,      DIMENSION(:), POINTER           :: jb, ib
   INTEGER,      DIMENSION(:), POINTER, OPTIONAL :: i_mumpsb
   ! local variables
   INTEGER      :: number_of_rows, i, p, c


   number_of_rows = SIZE(ia)-1

   ! check dimension consistency
   !
   IF ( n > number_of_rows ) THEN
      WRITE(*,*) '--> zEssM Error: n > rows(A)' 
      WRITE(*,*) '    n = ', n, ' SIZE(A,1) = ', number_of_rows
      WRITE(*,*) '    STOP.'
      STOP
   END IF


   ! count elements
   !
   c = 0
   DO i = 1, n

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

         ENDIF

      ENDDO

   ENDDO

   ALLOCATE( ib(n+1) )
   ALLOCATE( jb(c), b(c) )

   ! fill matrix
   !
   c = 0
   DO i = 1, n

      ib(i) = c + 1

      DO p = ia(i), ia(i+1) - 1

         IF ( ja(p) <= n ) THEN

            c = c + 1

            jb(c) = ja(p)
             b(c) =  a(p)

         ENDIF

      ENDDO

   ENDDO
   ib(n+1) = c + 1


   IF (PRESENT(i_mumpsb)) THEN

      ALLOCATE( i_mumpsb(SIZE(jb)) )

      DO i = 1, SIZE(ib) - 1
      
         DO p = ib(i), ib(i+1) - 1

            i_mumpsb(p) = i

         ENDDO

      END DO

   ENDIF

END SUBROUTINE zEssM
!
!

!==============================================================================

END MODULE sparse_matrix_operations
