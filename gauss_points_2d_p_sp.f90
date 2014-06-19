MODULE Gauss_points


!  AXISYMMETRIC


! Parabolic-Linear element with SubParametric transformation
! Storage of the derivative at the reference element
! The Derivative at the reference element are denoted by UPPER CASE


   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: k_d = 2,  n_w  = 6,  l_G  = 7,   &
                                           n_ws = 3,  l_Gs = 3

   INTEGER, PARAMETER, PUBLIC :: n_w_v = 3   ! number of vertex nodes
                              !! =========   ! the vertex nodes are assumed to
                                             ! be the first three nodes of the
                                             ! parabolic triangular element

   INTEGER, PARAMETER, PUBLIC :: n_ws_vs = 2 ! similarly for the surface edge
                              !! ===========

   REAL(KIND=8), DIMENSION(n_w,   l_G),           PUBLIC :: ww
   REAL(KIND=8), DIMENSION(n_ws,  l_Gs),          PUBLIC :: wws
   REAL(KIND=8), DIMENSION(k_d,   n_w,  l_G),     PUBLIC :: Dw_re    ! <=== NEW
   REAL(KIND=8), DIMENSION(1,     n_ws, l_Gs),    PUBLIC :: Dws_res  ! <=== NEW
   REAL(KIND=8), DIMENSION(l_G),                  PUBLIC :: pp_w
   REAL(KIND=8), DIMENSION(l_Gs),                 PUBLIC :: pp_ws

   REAL(KIND=8), DIMENSION(:),       ALLOCATABLE, PUBLIC :: JAC
   REAL(KIND=8), DIMENSION(:),       ALLOCATABLE, PUBLIC :: JACs
   REAL(KIND=8), DIMENSION(:, :, :), ALLOCATABLE, PUBLIC :: MNR
   REAL(KIND=8), DIMENSION(:, :),    ALLOCATABLE, PUBLIC :: MNRs
   REAL(KIND=8), DIMENSION(:, :, :), ALLOCATABLE, PUBLIC :: normals
   
   REAL(KIND=8), DIMENSION(:, :),    ALLOCATABLE, PUBLIC :: yy_G
   REAL(KIND=8), DIMENSION(:, :),    ALLOCATABLE, PUBLIC :: yy_Gs


!  REAL(KIND=8), DIMENSION(me),                   PUBLIC :: JAC
!  REAL(KIND=8), DIMENSION(mes),                  PUBLIC :: JACs
!  REAL(KIND=8), DIMENSION(k_d,  k_d,  me),       PUBLIC :: MNR
!  REAL(KIND=8), DIMENSION(k_d,        mes),      PUBLIC :: MNRs
!  REAL(KIND=8), DIMENSION(k_d,  l_Gs, mes),      PUBLIC :: normals

!  REAL(KIND=8), DIMENSION(l_G,  me),             PUBLIC :: yy_G
!  REAL(KIND=8), DIMENSION(l_Gs, mes),            PUBLIC :: yy_Gs




   PUBLIC Gauss_gen

CONTAINS

SUBROUTINE Gauss_gen (np, me, nps, mes, jj, jjs, rr)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: np, me, nps, mes

   INTEGER,      DIMENSION(n_w,  me),  INTENT(IN) :: jj
   INTEGER,      DIMENSION(n_ws, mes), INTENT(IN) :: jjs
   REAL(KIND=8), DIMENSION(k_d,  np),  INTENT(IN) :: rr

   REAL(KIND=8), DIMENSION(k_d,  n_w_v)   :: dw_v
   REAL(KIND=8), DIMENSION(1,    n_ws_vs) :: dws_vs
   REAL(KIND=8), DIMENSION(k_d,      k_d) :: dr
   REAL(KIND=8), DIMENSION(1,        k_d) :: drs

   INTEGER :: m, k, h, l, ms, ls

   m = nps ! otherwise nps is not used

   IF (ALLOCATED(JAC)) THEN

      DEALLOCATE (JAC, JACs, MNR, MNRs, normals, yy_G, yy_Gs)

   ENDIF

   ALLOCATE (JAC (me),  MNR (k_d, k_d, me))
   ALLOCATE (JACs(mes), MNRs(k_d, mes))
   ALLOCATE (normals(k_d, l_Gs, mes))

   ALLOCATE (yy_G(l_G, me), yy_Gs(l_Gs, mes))

!  volume elements

!  evaluate the values and derivatives of all weigthing functions
!  at all Gauss points of the parabolic reference element

   CALL element_2d_p2 (ww, Dw_re, pp_w)

   CALL element_2d_p1_v (dw_v) ! for linear elements the derivatives 
                               ! are constant

!  evaluate and store the constant values of the Jacobian determinant
!  and of the minors for each tetrahedral element
!  SUB-parametric transformation is performed

   DO m = 1, me

      ! dr(k,h) = dr_k^phys / dr_h^re 
     
      DO k = 1, k_d
         DO h = 1, k_d
            dr(k, h) = SUM( rr(k, jj(1:n_w_v, m)) * dw_v(h, :) )
         ENDDO
      ENDDO

      MNR(1,1, m) =   dr(2,2);   MNR(1,2, m) = - dr(2,1)
      MNR(2,1, m) = - dr(1,2);   MNR(2,2, m) =   dr(1,1)
      
      JAC(m) = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)


      ! Modifications for axisymmetric problems 

      DO l = 1, l_G 

         yy_G(l, m) = SUM(rr(2, jj(:, m)) * ww(:,l))
      
      ENDDO
   

   ENDDO


!  boundary elements

   CALL element_1d_p2 (wws, Dws_res, pp_ws)

   CALL element_1d_p1_vs (dws_vs) ! for linear elements the derivatives 
                                  ! are constant

   DO ms = 1, mes

      DO k = 1, k_d
         drs(1, k) = SUM( rr(k, jjs(1:n_ws_vs, ms)) * dws_vs(1, :) )
      ENDDO

      MNRs(1, ms) = + drs(1,2)   
      MNRs(2, ms) = - drs(1,1)
     
      JACs(ms) = SQRT( drs(1,1)**2 + drs(1,2)**2 )

      DO ls = 1, l_Gs
         normals(:, ls, ms) = MNRs(:, ms)/JACs(ms)  ! outward normal
      ENDDO


      ! Modifications for axisymmetric problems 

      DO ls = 1, l_Gs 
       
         yy_Gs(ls, ms) = SUM(rr(2, jjs(:, ms)) * wws(:,ls))
      
      ENDDO
  

   ENDDO


!   PRINT*, 'end of gen_Gauss_2d_p2_sp (sub-parametric)'


   CONTAINS
   !=======


   SUBROUTINE element_2d_p2 (w, d, p)

!     triangular element with quadratic interpolation
!     and seven Gauss integration points

!     Degree 5, 7 Points formula for a triangular domain
!     Radon, Hammer, Marlowe and Stroud
!
!     A. H. Stroud, 
!     Approximate Calculation of Multiple Integrals   
!     Prentice Hall, Englewood Cliffs, 1971
!     page 314   

!        w(n_w, l_G) : values of shape functions at Gauss points
!     d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
!             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
      INTEGER :: j

      REAL(KIND=8) :: zero = 0,  half  = 0.5,  one  = 1,  &
                      two  = 2,  three = 3,    four = 4,  &
                      five = 5,   nine = 9

      REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
                      df1x, df2x, df3x, df4x, df5x, df6x, &
                      df1y, df2y, df3y, df4y, df5y, df6y, &
                      x, y, r, a,  s1, s2, t1, t2, b1, b2,  area, sq


      f1(x, y) = (half - x - y) * (one - x - y) * two
      f2(x, y) = x * (x - half) * two
      f3(x, y) = y * (y - half) * two
      f4(x, y) = x * y * four
      f5(x, y) = y * (one - x - y) * four
      f6(x, y) = x * (one - x - y) * four

      df1x(x, y) = -three + four * (x + y)
      df2x(x, y) = (two*x - half) * two
      df3x(x, y) = zero
      df4x(x, y) =  y * four
      df5x(x, y) = -y * four
      df6x(x, y) = (one - two*x - y) * four

      df1y(x, y) = -three + four * (x + y)
      df2y(x, y) = zero
      df3y(x, y) = (two*y - half) * two
      df4y(x, y) =  x * four
      df5y(x, y) = (one - x - two*y) * four
      df6y(x, y) = -x * four

!     Degree 5, 7 Points;  Stroud: p. 314, Approximate calculation of
!                          Multiple integrals (Prentice--Hall), 1971.

      area = one/two

      sq = SQRT(three*five)

      r  = one/three;                          a = area * nine/40

      s1 = (6 - sq)/21;  t1 = (9 + 2*sq)/21;  b1 = area * (155 - sq)/1200

      s2 = (6 + sq)/21;  t2 = (9 - 2*sq)/21;  b2 = area * (155 + sq)/1200

      xx(1) = r;    yy(1) = r;    p(1) = a

      xx(2) = s1;   yy(2) = s1;   p(2) = b1
      xx(3) = s1;   yy(3) = t1;   p(3) = b1
      xx(4) = t1;   yy(4) = s1;   p(4) = b1

      xx(5) = s2;   yy(5) = s2;   p(5) = b2
      xx(6) = s2;   yy(6) = t2;   p(6) = b2
      xx(7) = t2;   yy(7) = s2;   p(7) = b2


      DO j = 1, l_G

            w(1, j) =  f1 (xx(j), yy(j))
         d(1, 1, j) = df1x(xx(j), yy(j))
         d(2, 1, j) = df1y(xx(j), yy(j))

            w(2, j) =  f2 (xx(j), yy(j))
         d(1, 2, j) = df2x(xx(j), yy(j))
         d(2, 2, j) = df2y(xx(j), yy(j))

            w(3, j) =  f3 (xx(j), yy(j))
         d(1, 3, j) = df3x(xx(j), yy(j))
         d(2, 3, j) = df3y(xx(j), yy(j))

            w(4, j) =  f4 (xx(j), yy(j))
         d(1, 4, j) = df4x(xx(j), yy(j))
         d(2, 4, j) = df4y(xx(j), yy(j))

            w(5, j) =  f5 (xx(j), yy(j))
         d(1, 5, j) = df5x(xx(j), yy(j))
         d(2, 5, j) = df5y(xx(j), yy(j))

            w(6, j) =  f6 (xx(j), yy(j))
         d(1, 6, j) = df6x(xx(j), yy(j))
         d(2, 6, j) = df6y(xx(j), yy(j))

      ENDDO

   END SUBROUTINE element_2d_p2


!------------------------------------------------------------------------------


   SUBROUTINE element_1d_p2(w, d, p)

!     one-dimensional element with quadratic interpolation
!     and three Gauss integration points

!     Degree 5, 3 Points Gauss integration formula

!        w(n_w, l_G) : values of shape functions at Gauss points
!     d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
!             p(l_G) : weight for Gaussian quadrature at Gauss points

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
      REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p

      REAL(KIND=8), DIMENSION(l_Gs) :: xx
      INTEGER :: j
 
      REAL(KIND=8) :: zero = 0,   one = 1,    two = 2,  three = 3,  &
                      five = 5, eight = 8,   nine = 9

      REAL(KIND=8) :: f1, f2, f3, df1, df2, df3, x

      f1(x) = (x - one)*x/two
      f2(x) = (x + one)*x/two
      f3(x) = (x + one)*(one - x)

      df1(x) = (two*x - one)/two
      df2(x) = (two*x + one)/two
      df3(x) = -two*x

      xx(1) = -SQRT(three/five)
      xx(2) =  zero
      xx(3) =  SQRT(three/five) 

      p(1)  =  five/nine
      p(2)  =  eight/nine 
      p(3)  =  five/nine 


      DO j = 1, l_Gs

            w(1, j) =  f1(xx(j))
         d(1, 1, j) = df1(xx(j))

            w(2, j) =  f2(xx(j))
         d(1, 2, j) = df2(xx(j))

            w(3, j) =  f3(xx(j))
         d(1, 3, j) = df3(xx(j))


      ENDDO

   END SUBROUTINE element_1d_p2


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


   SUBROUTINE element_2d_p1_v (d)

!     triangular element with linear interpolation

!     Degree 2, 3 Points formula for a triangular domain
!     Hammer and Stroud
!
!     A. H. Stroud, 
!     Approximate Calculation of Multiple Integrals   
!     Prentice Hall, Englewood Cliffs, 1971
!     page 307 and 368  

!     d(1:2, n_w_v) : CONSTANT derivatives of the shape functions

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(2, n_w_v), INTENT(OUT) :: d

      REAL(KIND=8) :: zero = 0,  one = 1

!     REAL(KIND=8) :: f1, f2, f3, x, y
!     f1(x, y) = one - x - y 
!     f2(x, y) = x
!     f3(x, y) = y

     
      d(1,1) = - one;   d(1,2) = one;    d(1,3) = zero
      d(2,1) = - one;   d(2,2) = zero;   d(2,3) = one
      

   END SUBROUTINE element_2d_p1_v

!------------------------------------------------------------------------------

   SUBROUTINE element_1d_p1_vs (d)

!     edge element with linear interpolation

!     Degree 3, 2 Points Gauss integration formula

!     d(1, n_ws_vs) : CONSTANT derivatives of the shape functions

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(1, n_ws_vs), INTENT(OUT) :: d

      REAL(KIND=8) :: two = 2

!     REAL(KIND=8) :: f1, f2, x, y
!     f1(x) = (one - x)/two
!     f2(x) = (x + one)/two
      
     
      d(1,1) = - 1/two;   d(1,2) = 1/two
     

   END SUBROUTINE element_1d_p1_vs


END SUBROUTINE Gauss_gen


END MODULE Gauss_points
