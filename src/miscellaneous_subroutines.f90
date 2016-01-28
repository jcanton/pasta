MODULE miscellaneous_subroutines
!
! Contains subroutines which did not fit anywhere else
!
   USE prep_mesh_p1p2_sp ! for some global variables as mm, jj
   USE global_variables
   USE Gauss_points ! needed for the field average subroutine
!
!==============================================================================

   IMPLICIT NONE

   INTERFACE computeFieldAverage
      MODULE PROCEDURE computeFieldAverage_rv
   END INTERFACE computeFieldAverage

   INTERFACE computeFieldIntegral
      MODULE PROCEDURE computeFieldIntegral_r
   END INTERFACE computeFieldIntegral

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE collect (uu, ppL,  xx)

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ppL
   
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: xx

   INTEGER :: np, np_L
 
   np   = SIZE(uu, 2)
   np_L = SIZE(ppL)
   
   xx(1:np) = uu(1,:)
  
   xx(np+1 : 2*np) = uu(2,:)
   
   xx(2*np + 1 : 3*np) = uu(3,:)
  
   xx(3*np + 1 : 3*np + np_L) = ppL
 
END SUBROUTINE collect

!------------------------------------------------------------------------------

SUBROUTINE extract (xx,  uu,  ppL)

   IMPLICIT NONE
    
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: xx
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: uu
   REAL(KIND=8), DIMENSION(:),   INTENT(OUT), OPTIONAL :: ppL
   
   INTEGER :: np, np_L
 
   np = SIZE(uu, 2)
   
   uu(1,:) = xx(1:np) 
  
   uu(2,:) = xx(np+1 : 2*np)
  
   uu(3,:) = xx(2*np+1 : 3*np)
 
   IF (PRESENT(ppL)) THEN 
     
      np_L = SIZE(ppL)
     
      ppL = xx(3*np+1 : 3*np + np_L) 
 
   ENDIF
 
END SUBROUTINE extract

!------------------------------------------------------------------------------

SUBROUTINE collect_cmplx (uu, ppL,  xx)

   IMPLICIT NONE
   
   COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN)  :: uu
   COMPLEX(KIND=8), DIMENSION(:),   INTENT(IN)  :: ppL
   
   COMPLEX(KIND=8), DIMENSION(:)                :: xx

   INTEGER :: np, np_L
 
   np   = SIZE(uu, 2)
   np_L = SIZE(ppL)

   xx(1:np) = uu(1,:)
  
   xx(np+1 : 2*np) = uu(2,:)
   
   xx(2*np + 1 : 3*np) = uu(3,:)
  
   xx(3*np + 1 : 3*np + np_L) = ppL
 
!   write(*,*) 'collect_cmplx'
!   write(*,*) 'np   = ', np
!   write(*,*) 'np_L = ', np_L
!   write(*,*)
!   write(*,*) uu(1,:)
!   write(*,*)
!   write(*,*) uu(2,:)
!   write(*,*)
!   write(*,*) uu(3,:)
!   write(*,*)
!   write(*,*) ppL
   
END SUBROUTINE collect_cmplx

!------------------------------------------------------------------------------

SUBROUTINE extract_cmplx (xx,  uu,  ppL)

   IMPLICIT NONE
    
   COMPLEX(KIND=8), DIMENSION(:),   INTENT(IN)  :: xx
   COMPLEX(KIND=8), DIMENSION(:,:)              :: uu
   COMPLEX(KIND=8), DIMENSION(:),   OPTIONAL    :: ppL
   
   INTEGER :: np, np_L
 
   np = SIZE(uu, 2)
   
   uu(1,:) = xx(1:np) 
  
   uu(2,:) = xx(np+1 : 2*np)
  
   uu(3,:) = xx(2*np+1 : 3*np)
 
   IF (PRESENT(ppL)) THEN 
     
      np_L = SIZE(ppL)
     
      ppL = xx(3*np+1 : 3*np + np_L) 
 
   ENDIF
 
END SUBROUTINE extract_cmplx


!------------------------------------------------------------------------------

SUBROUTINE intToChar6 (intNumber, sixCharString)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 26/1/2014
!
! - intNumber     :: integer number
! - sixCharString :: six characters string with zeros
!
   IMPLICIT NONE
   ! input variables
   INTEGER, INTENT(IN) :: intNumber
   ! output variables
   CHARACTER(*) :: sixCharString

   IF ( intNumber < 10 ) THEN
      WRITE(sixCharString, '(a5,i1)') '00000', intNumber
   ELSEIF (intNumber < 100 ) THEN
      WRITE(sixCharString, '(a4,i2)') '0000',  intNumber
   ELSEIF (intNumber < 1000 ) THEN
      WRITE(sixCharString, '(a3,i3)') '000',   intNumber
   ELSEIF (intNumber < 10000 ) THEN
      WRITE(sixCharString, '(a2,i4)') '00',    intNumber
   ELSEIF (intNumber < 100000 ) THEN
      WRITE(sixCharString, '(a1,i5)') '0',     intNumber
   ELSEIF (intNumber < 1000000 ) THEN
      WRITE(sixCharString, '(i6)')             intNumber
   ELSE
      WRITE(*,*) '*******************************************'
      WRITE(*,*) '*** Error:                              ***'
      WRITE(*,*) '*** Integer too large.                  ***'
      WRITE(*,*) '*** Change subroutine intToChar6        ***'
      WRITE(*,*) '*** in module miscellaneous_subroutines ***'
      WRITE(*,*) '*******************************************'
      WRITE(*,*) 'STOP.'
      CALL MPI_ABORT(MPI_COMM_WORLD, mpiErrC, mpiIerr)
   ENDIF

END SUBROUTINE intToChar6

!------------------------------------------------------------------------------

FUNCTION testForNaN (variable) RESULT (flag)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 4/5/2014
!
! - variable: real type variable to be tested

   IMPLICIT NONE
   !
   REAL(KIND=8), DIMENSION(:) :: variable
   !
   LOGICAL :: flag

   IF ( .NOT. ALL(variable.EQ.variable) ) THEN
      flag = .TRUE.
   ELSE
      flag = .FALSE.
   ENDIF

END FUNCTION

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE computeFieldIntegral_r (fld,  intFld)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 29/9/2015
!
! - fld        : field to be integrated
!
! - intFld     : integrated values of the field

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(np), INTENT(IN) :: fld
   ! output variable
   REAL(KIND=8) :: intFld
   ! local variables
   REAL(KIND=8), DIMENSION(np) :: intvFld
   INTEGER :: mi, m, l

   intFld = 0
	! outer cycle on elements.  m = current element index
	!
   DO mi = 1, SIZE(mm);  m = mm(mi)

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G

         intFld = intFld + SUM(fld(jj(:,m)) * ww(:,l)) * JAC(m) * pp_w(l)

      ENDDO

   ENDDO

END SUBROUTINE computeFieldIntegral_r

!------------------------------------------------------------------------------

SUBROUTINE computeFieldAverage_rv (fld,  avgFld_v)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 29/9/2015
!
! - fld          : field to be averaged
!
! - avgFld_v     : averaged values of the field

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np), INTENT(IN) :: fld
   ! output variable
   REAL(KIND=8), DIMENSION(velCmpnnts) :: avgFld_v
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: ones
   REAL(KIND=8), DIMENSION(velCmpnnts)    :: avgFld, area
   INTEGER :: mi, m, l, k

   ones   = 1d0
   avgFld = 0d0

	! outer cycle on elements.  m = current element index
	!
   DO mi = 1, SIZE(mm);  m = mm(mi)

		! cycle on parabolic Gauss points [ l_G = 7 ]
		!
      DO l = 1, l_G

         ! cycle on velocity components
         !
         DO k = 1, velCmpnnts

            avgFld(k) = avgFld(k) + SUM(fld(k, jj(:,m))  * ww(:,l)) * JAC(m) * pp_w(l)
            area(k)   = area(k)   + SUM(ones(k, jj(:,m)) * ww(:,l)) * JAC(m) * pp_w(l)

         ENDDO

      ENDDO

   ENDDO

   DO k = 1, velCmpnnts

      avgFld_v(k) = avgFld(k) / area(k)

   ENDDO


END SUBROUTINE computeFieldAverage_rv

!==============================================================================

END MODULE miscellaneous_subroutines
