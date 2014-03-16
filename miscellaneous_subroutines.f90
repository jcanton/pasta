MODULE miscellaneous_subroutines
!
! Contains subroutines which did not fit anywhere else
!
!==============================================================================

   IMPLICIT NONE

CONTAINS

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
      STOP
   ENDIF

END SUBROUTINE intToChar6

!==============================================================================

END MODULE miscellaneous_subroutines
