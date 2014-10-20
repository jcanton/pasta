! -----------------------------------------------------------------------------
!   LOCA 1.0: Library of Continuation Algorithms
!   Copyright (C) 2001, Sandia National Laboratories
! -----------------------------------------------------------------------------

MODULE Loca_util

   USE Loca_types
   USE loca_wrappers

   IMPLICIT NONE

   INTEGER      :: nn_o = -1  ! Length of vector that is acted on (owned unknowns)
   INTEGER      :: nn_t = -1  ! Length of vector that is allcoated (total unknowns)
   REAL(KIND=8) :: nn_g = -1  ! Total number of unknowns on all procs, cast to a double, which is the sum of nn_o's
   REAL(KIND=8), DIMENSION(:), POINTER :: SclVec ! Pointer to the scaling vector, so that routines that use it won't have to pass it in


CONTAINS
!=======

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
SUBROUTINE initialize_util_routines(n_o, n_t) ! ok
! This routine sets up information that is used by the rouitnes in this file
! and must be called before any of the routines within it are called.
! The current input are "int n_o" the number of owned unknowns for this
! processors, which is the length of a vector that is acted upon, and
! "int n_t" which is the length that the vector must be allocated

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: n_o, n_t

   nn_o = n_o
   nn_t = n_t
   nn_g = gsum_double_conwrap(DBLE(nn_o))

END SUBROUTINE initialize_util_routines

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!SUBROUTINE init_scale_utils(scale_vec)
!! Make scale_vec accessible within this file to avoid passing it through MF
!! Used by 'loca_mf'
!
!   IMPLICIT NONE
!
!   REAL(KIND=8), DIMENSION(:), TARGET :: scale_vec
!
!   SclVec => scale_vec
!
!END SUBROUTINE init_scale_utils

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!SUBROUTINE vec_init(u) USELESS
!! Initialize the vector to zero.
!!
!
!   IMPLICIT NONE
!
!   REAL(KIND=8), DIMENSION(:) :: u
!
!   u = 0.d0
!
!END SUBROUTINE vec_init


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!SUBROUTINE vec_copy(dx, dy) USELESS
!! This function copies the value of the vector of type double, dx, into the
!! vector of type double, dy.  No error checking is done. If N is negative or
!! zero, dy is not changed.  A stride of 1 is assumed, unlike the linpack
!! routine called DCOPY (which is why the name is dcopy1).
!
!   IMPLICIT NONE
!
!   REAL(KIND=8), DIMENSION(:) :: dx, dy
!
!   dy = dx
!
!END SUBROUTINE vec_copy

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
FUNCTION dp(x, y) RESULT(output)
! simple dot product

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x, y

   REAL(KIND=8) :: output

   output = SUM(x * y)

   output = gsum_double_conwrap(output)

END FUNCTION dp

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
FUNCTION scaled_dp(x, y) RESULT(output)
! simple dot product scaled by scaling vector and length

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:) :: x, y

   REAL(KIND=8) :: output

   output = SUM(x*SclVec * y*SclVec)

   output = gsum_double_conwrap(output/nn_g)

END FUNCTION scaled_dp

!***************************************************************************
!***************************************************************************
!***************************************************************************

FUNCTION scaled_dot_prod(x, y, scale_vec, ownedUnknowns) RESULT(output) ! ok

   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: x, y, scale_vec
   INTEGER, OPTIONAL :: ownedUnknowns

   REAL(KIND=8) :: output

   output = SUM(x * y * scale_vec * scale_vec)

   output = gsum_double_conwrap(output)

END FUNCTION scaled_dot_prod

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
FUNCTION ip(x, y) RESULT(output)
! inner product for pitchfork tracking. This inner product must
! have the same symmetry as the pitchfork, so a simple dot product
! only works if the mesh is symmetric. Otherwise, this product
! must be weighted by the lumped mass at each node */

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:) :: x, y

   REAL(KIND=8) :: output

   output = dp(x,y)

END FUNCTION ip

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
FUNCTION ltransnorm(x, scale_vec) RESULT(output) ! ok
! rescale by factor to make averages array element size one

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(0:) :: x, scale_vec

   REAL(KIND=8) :: output

   output = dp(x, scale_vec) / nn_g

END FUNCTION ltransnorm

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
FUNCTION null_vector_resid(r_val, i_val, r_vec, i_vec, mm_flag, b_flag) RESULT(output)
! Calculate bifurcation equation residual norm (as Rayleigh quotient of null vector).

   IMPLICIT NONE
   REAL(KIND=8)                :: r_val, i_val
   REAL(KIND=8), DIMENSION(0:) :: r_vec, i_vec
   INTEGER                     :: mm_flag
   INTEGER, OPTIONAL           :: b_flag

   REAL(KIND=8) :: output

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::  vr, vi, v1
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::  Jy, Jz, By, Bz
   REAL(KIND=8) :: yJy, yJz, zJy, zJz, yBy, yBz, zBy, zBz
   REAL(KIND=8) :: num_r, num_i, denom, norm_r, norm_i, norm

   ! Proceed according to bifurcation method
   IF (i_val == 0d0) THEN

      ! For these cases, the null vector is real:
      !   When a mass matrix is available:
      !       R.Q. = (n.J.n) / (n.B.n)
      !
      !   When a mass matrix is NOT available:
      !       R.Q. = (n.J.n) / (n.n)

      ! r_vec needs to be copied to vector with space for externals
      ALLOCATE(vr(0:nn_o-1))
      vr = r_vec

      ALLOCATE(v1(0:nn_o-1))
      CALL matvec_mult_conwrap(vr, v1)
      num_r = dp(vr, v1)

      ! v1 is either n or B.n, depending on mm_flag

      IF (mm_flag == TRUE) THEN
         CALL mass_matvec_mult_conwrap(vr, v1)
      ELSE
         v1 = vr
      END IF

      denom = dp(vr, v1)

      IF(DABS(denom) > 1d-20) THEN
         norm_r = r_val - num_r / denom
         norm = DABS(norm_r)
      ELSE
         norm = -1d0
      END IF

      DEALLOCATE(vr)
      DEALLOCATE(v1)

   ELSE
      !
      ! At a Hopf bifurcation, the eigenvalue is imaginary and the null vector is
      ! complex, so the following formula is used for the complex Rayleigh quotient:
      !
      !                           h          h
      !       R.Q. = (omega)i - (n .J.n) / (n .B.n)
      !
      !        h
      ! where n  is the conjugate transpose of null vector n, J is the LSA
      ! Jacobian matrix, and B is the mass matrix.
      !
      ! Note: In this routine, y and z refer to the real and imaginary parts
      ! of n, respectively.

      ! r_vec needs to be copied to vector with space for externals
      ALLOCATE(vr(0:nn_o-1))

      vr = r_vec

      ALLOCATE(vi(0:nn_o-1))

      vi = i_vec

      ! Allocate work vectors
      ALLOCATE(Jy(0:nn_o-1))
      ALLOCATE(Jz(0:nn_o-1))
      ALLOCATE(By(0:nn_o-1))
      ALLOCATE(Bz(0:nn_o-1))

      ! Use work vectors for matrix-vector products
      IF (PRESENT(b_flag)) THEN
         CALL Lns_matvec_mult_conwrap(vr, vi, Jy, Jz)
      ELSE
         CALL matvec_mult_conwrap(vr, Jy)
         CALL matvec_mult_conwrap(vi, Jz)
      ENDIF

      CALL mass_matvec_mult_conwrap(vr, By)
      CALL mass_matvec_mult_conwrap(vi, Bz)

      ! Intermediate terms are then obtained by dot products
      yJy = dp(vr, Jy)
      yJz = dp(vr, Jz)
      zJy = dp(vi, Jy)
      zJz = dp(vi, Jz)
      yBy = dp(vr, By)
      yBz = dp(vr, Bz)
      zBy = dp(vi, By)
      zBz = dp(vi, Bz)

      ! First, get denominator and check if it is too small to divide by
      denom = (yBy + zBz) * (yBy + zBz) + (yBz - zBy) * (yBz - zBy)

      IF (DABS(denom) == 0d0) THEN

         norm = -1d0

      ELSE

         ! Get real and imaginary parts of numerator
         num_r = (yJy + zJz) * (yBy + zBz) + (yJz - zJy) * (yBz - zBy)
         num_i = (yJz - zJy) * (yBy + zBz) - (yJy + zJz) * (yBz - zBy)

         ! Construct the complex norm: (norm_r) + (norm_i)i
         norm_r = r_val - num_r / denom
         norm_i = i_val - num_i / denom

         norm = SQRT(norm_r*norm_r + norm_i*norm_i)

      END IF

      ! Clean up and return */
      DEALLOCATE(vr)
      DEALLOCATE(vi)
      DEALLOCATE(Jy)
      DEALLOCATE(Jz)
      DEALLOCATE(By)
      DEALLOCATE(Bz)

   END IF

   ! Flag value of -1.0 means denominator of Rayleigh Quotient is zero,
   ! which is a zero mass matrix or zero eigenvectors
   output = norm

END FUNCTION null_vector_resid

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

!SUBROUTINE sort_by_real(nconv, ncv, ldv, d, v)
!! This routine takes the converged eigenvalue array d and assigns
!! index array j_index for sorting from highest to lowest real part.
!! Complex eigenvalue pairs will be identified with equal index values.
!
!   IMPLICIT NONE
!
!   INTEGER :: nconv, ncv, ldv
!   REAL(KIND=8), DIMENSION(:) :: d, v
!
!   INTEGER :: i, j, k, l, m, skip, skip_next
!   INTEGER :: ntd, ntv
!   INTEGER,      DIMENSION(:), ALLOCATABLE :: count
!   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: td, tv
!
!   ! First check if nconv is smaller than two.
!   IF (nconv < 2) RETURN
!
!   ! Array count[i] will record how many eigenvalues are more than the i'th */
!   ntd = 2 * ncv
!
!   ntv = (1 + nconv) * ldv
!
!   ALLOCATE(count(0:ncv))
!
!   count = 0
!
!   ! Mark final element of count array with -1
!   count(ncv) = -1
!
!   ! Arrays td and tv will hold the sorted eigenvalues and eigenvectors
!   ALLOCATE(td(0:ntd-1))
!
!   td = 0.d0
!
!   ALLOCATE(tv(0:ntv-1))
!
!   tv = 0.d0
!
!   ! Compare real parts for all i,j eigenvalue pairs, where j has the larger
!   ! initial index. Increment count for the smaller real part.
!   ! Here, when eigenvalue j is the first of a complex eigenpair, "skip_next"
!   ! will be TRUE; if it is the second of the pair, "skip" will be TRUE.
!   skip = FALSE
!   skip_next = FALSE
!   DO i=0, nconv-1
!
!      ! Determine if this is the first of a complex eigenpair
!      ! If this is the second of a complex eigenpair, reset the skip_next flag
!      IF (skip == TRUE) THEN
!         skip_next = FALSE
!      ELSE
!         IF (d(i+ncv) == 0.d0) THEN
!            skip_next = FALSE
!         ELSE
!            skip_next = TRUE
!         END IF
!      END IF
!
!      DO j=0, i-1
!
!         ! Do not compare complex conjugates - this ensures
!         ! that both will have the same value of count
!         IF (skip == FALSE .OR. i /= j-1) THEN
!
!           ! If d values are different, increment count for the smaller one. */
!           IF (d(j) < d(i)) THEN
!             count(j) = count(j) + 1
!           ELSE IF (d(j) > d(i)) THEN
!             count(i) = count(i) + 1
!           ! If d values are the same but not a complex eigenpair,
!           ! increment count for the larger index (always j)
!           ELSE
!             count(i) = count(i) + 1
!           END IF
!
!         END IF
!
!      END DO
!
!      ! set skip for next pass of j
!      skip = skip_next
!   END DO
!
!   ! Now copy eigenvalues and eigenvectors into temporary arrays td and tv
!   ! in their sorted positions as determined by the count array */
!   skip = FALSE
!   skip_next = FALSE
!   DO j=0, nconv-1
!
!      ! Initial position: j --> Sorted position: i = count[j]
!      i = count(j)
!
!      ! Determine if this is the second of a complex eigenpair;
!      ! if so, copying was done on previous pass */
!      IF (skip == FALSE) THEN
!
!         ! Determine if this is the first of a complex eigenpair, set skip_next
!         IF (count(j+1) == i) THEN
!            skip_next = TRUE
!         ELSE
!            skip_next = FALSE
!         END IF
!
!         ! Copy eigenvalue into td
!         td(i)     = d(j)
!         td(i+ncv) = d(j+ncv)
!
!         ! Copy complex conjugate eigenvalue into next td position if applicable
!         IF (skip_next == TRUE) THEN
!            td(i+1)     = d(j)
!            td(i+ncv+1) = d(j+ncv+1)
!         END IF
!
!         ! Copy eigenvector into tv
!         DO k = 0, ldv-1
!
!           ! Assign indices into v and tv
!           l = i * ldv + k
!           m = j * ldv + k
!           tv(l) = v(m)
!
!           ! Copy corresponding element of complex conjugate eigenvector
!           ! into next tv position if applicable
!           IF (skip_next == TRUE)  tv(l+ldv) = v(m+ldv)
!
!         END DO
!
!      ! If this is the second of a complex eigenpair, just reset skip_next
!      ELSE
!        skip_next = FALSE
!      END IF
!
!      ! set skip_next for next pass
!      skip = skip_next
!
!   END DO
!
!   ! Now rewrite d and v arrays in sorted order
!   DO i = 0, nconv-1
!      d(i) = td(i)
!      d(i+ncv) = td(i+ncv)
!      DO j = 0, ldv-1
!         k = i * ldv + j
!         v(k) = tv(k)
!      END DO
!   END DO
!
!   ! Free temporary arrays.
!   DEALLOCATE(count)
!   DEALLOCATE(td)
!   DEALLOCATE(tv)
!
!END SUBROUTINE sort_by_real

END MODULE Loca_util

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
