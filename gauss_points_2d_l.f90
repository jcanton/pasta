MODULE Gauss_points_L


!  AXISYMMETRIC


   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d_L = 2,  n_w_L  = 3,  l_G_L  = 3,   &
                                             n_ws_L = 2,  l_Gs_L = 2
 
   REAL(KIND=8), DIMENSION(n_w_L,  l_G_L),           PUBLIC :: ww_L
   REAL(KIND=8), DIMENSION(n_ws_L, l_Gs_L),          PUBLIC :: wws_L

   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dw_L
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: jac_py_L
   
   REAL(KIND=8), DIMENSION(:, :),       ALLOCATABLE, PUBLIC :: jac_psy_L
   REAL(KIND=8), DIMENSION(:,    :, :), ALLOCATABLE, PUBLIC :: normals_L
   REAL(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE, PUBLIC :: dws_L, dw_psy_L

  !REAL(KIND=8), DIMENSION(k_d_L,  n_w_L,  l_G_L,   me),  PUBLIC :: dw_L 
  !REAL(KIND=8), DIMENSION(l_G_L,   me),                  PUBLIC :: jac_py_L
  
  !REAL(KIND=8), DIMENSION(k_d_L,          l_Gs_L,  mes), PUBLIC :: normals_L
  !REAL(KIND=8), DIMENSION(l_Gs_L,  mes),                 PUBLIC :: jac_psy_L
  !REAL(KIND=8), DIMENSION(k_d_L-1,n_ws_L, l_Gs_L,  mes), PUBLIC :: dws_L, dw_psy_L

   PUBLIC Gauss_gen_L


CONTAINS


SUBROUTINE Gauss_gen_L (np, me, nps, mes, jj, jjs, rr)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: np, me, nps, mes

   INTEGER,      DIMENSION(n_w_L,  me),  INTENT(IN) :: jj
   INTEGER,      DIMENSION(n_ws_L, mes), INTENT(IN) :: jjs
   REAL(KIND=8), DIMENSION(k_d_L,  np),  INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(k_d_L, n_w_L,  l_G_L) :: dd
   REAL(KIND=8), DIMENSION( 1 ,  n_ws_L, l_Gs_L) :: dds
   REAL(KIND=8), DIMENSION(l_G_L)                :: pp
   REAL(KIND=8), DIMENSION(l_Gs_L)               :: pps

   REAL(KIND=8), DIMENSION(n_w_L)        :: r
   REAL(KIND=8), DIMENSION(n_ws_L)       :: rs
   REAL(KIND=8), DIMENSION(k_d_L, k_d_L) :: dr
   REAL(KIND=8), DIMENSION( 1 , k_d_L)   :: drs

   REAL(KIND=8) :: jac, jacs
   INTEGER      :: m, l, k, h, n,  ms, ls

   m = nps ! otherwise nps is not used 

   IF (ALLOCATED(dw_L)) THEN
    
      DEALLOCATE(dw_L, jac_py_L,  jac_psy_L, normals_L, dws_L, dw_psy_L)
  
   END IF

   ALLOCATE(dw_L(k_d_L, n_w_L, l_G_L,  me))
   ALLOCATE(jac_py_L(l_G_L,  me))
   
   ALLOCATE(jac_psy_L(l_Gs_L, mes))
   ALLOCATE(normals_L(k_d_L,  l_Gs_L, mes))
   ALLOCATE(dws_L  (1,  n_ws_L,  l_Gs_L,  mes),  &  ! 1 = k_d_L - 1
            dw_psy_L(1,  n_ws_L,  l_Gs_L,  mes))     ! 1 =  2  -  1
  
!  evaluate and store the values of derivatives and of the
!  jacobian determinant at Gauss points of all volume elements

!  volume elements

   CALL element_2d (ww_L, dd, pp)

   DO m = 1, me

      DO l = 1, l_G_L

         DO k = 1, k_d_L
            r = rr(k, jj(:,m))
            DO h = 1, k_d_L
               dr(k, h) = SUM(r * dd(h,:,l))
            ENDDO
         ENDDO

         jac = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)

         DO n = 1, n_w_L
            dw_L(1, n, l, m)   &
               = (+ dd(1,n,l)*dr(2,2) - dd(2,n,l)*dr(2,1))/jac
            dw_L(2, n, l, m)   &
               = (- dd(1,n,l)*dr(1,2) + dd(2,n,l)*dr(1,1))/jac
         ENDDO

         jac_py_L(l, m) = jac * pp(l)
         
         ! modification for axisymmetric equations 

         jac_py_L(l, m) = SUM(rr(2, jj(:,m)) * ww_L(:,l)) * jac * pp(l)
                        !yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
     
      ENDDO

   ENDDO


!  surface elements

   CALL element_1d (wws_L, dds, pps)

   DO ms = 1, mes

      DO ls = 1, l_Gs_L

         DO k = 1, k_d_L
            rs = rr(k, jjs(:,ms))
            drs(1, k) = SUM(rs * dds(1,:,ls))
         ENDDO

         jacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )

         normals_L(1, ls, ms) = + drs(1,2)/jacs   ! outward normal
         normals_L(2, ls, ms) = - drs(1,1)/jacs   ! outward normal

         jac_psy_L(ls, ms) = jacs * pps(ls)
         
         ! modification for axisymmetric equations 
         
         jac_psy_L(ls, ms) = SUM(rr(2, jjs(:,ms)) * wws_L(:,ls)) * jacs * pps(ls)
                            !yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
      ENDDO

   ENDDO


   DO ms = 1, mes  ! necessary only for evaluating gradient
                   ! tangential to the surface (ex COMMON Gauss_tan)

      DO ls = 1, l_Gs_L
      
           dws_L(1, :, ls, ms) = dds(1, :, ls)
         dw_psy_L(1, :, ls, ms) = dds(1, :, ls) * pps(ls)
      
      ENDDO

   ENDDO


   PRINT*, 'end of gen_Gauss_L'


   CONTAINS
   !=======

   SUBROUTINE element_2d (w, d, p)

!     triangular element with linear interpolation

!     Degree 2, 3 Points formula for a triangular domain
!     Hammer and Stroud
!
!     A. H. Stroud, 
!     Approximate Calculation of Multiple Integrals   
!     Prentice Hall, Englewood Cliffs, 1971
!     page 307 and 368  

!        w(n_w_L, l_G_L) : values of shape functions at Gauss points
!     d(2, n_w_L, l_G_L) : derivatives values of shape functions at Gauss points
!             p(l_G_L) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w_L, l_G_L), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_w_L, l_G_L), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G_L),           INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G_L) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6

      REAL(KIND=8) :: f1, f2, f3, x, y
      f1(x, y) = one - x - y
      f2(x, y) = x
      f3(x, y) = y

      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

      DO j = 1, l_G_L

            w(1, j) = f1(xx(j), yy(j))
         d(1, 1, j) = - one
         d(2, 1, j) = - one

            w(2, j) = f2(xx(j), yy(j))
         d(1, 2, j) = one
         d(2, 2, j) = zero

            w(3, j) = f3(xx(j), yy(j))
         d(1, 3, j) = zero
         d(2, 3, j) = one

               p(j) = one/six

      ENDDO

   END SUBROUTINE element_2d

!------------------------------------------------------------------------------

   SUBROUTINE element_1d (w, d, p)

!     one-dimensional element with linear interpolation

!     Degree 3, 2 Points Gauss integration formula

!        w(n_w_L, l_G_L) : values of shape functions at Gauss points
!     d(1, n_w_L, l_G_L) : derivatives values of shape functions at Gauss points
!             p(l_G_L) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws_L, l_Gs_L), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(1, n_ws_L, l_Gs_L), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs_L),            INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs_L) :: xx
      INTEGER :: j

      REAL(KIND=8) :: one = 1,  two = 2,  three = 3

      REAL(KIND=8) :: f1, f2, x
      f1(x) = (one - x)/two
      f2(x) = (x + one)/two

      xx(1) = - one/SQRT(three)
      xx(2) = + one/SQRT(three)

      DO j = 1, l_Gs_L

            w(1, j) = f1(xx(j))
         d(1, 1, j) = - one/two

            w(2, j) = f2(xx(j))
         d(1, 2, j) = + one/two

               p(j) = one

      ENDDO

   END SUBROUTINE element_1d


END SUBROUTINE Gauss_gen_L

END MODULE Gauss_points_L
