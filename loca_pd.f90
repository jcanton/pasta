MODULE Loca_pd
!******************************************************************************
!******************************************************************************
!******************************************************************************

! Define passdown structure: this structure is global to this file, and
! provides a way to pass variables from the top solve_continuation routine
! down to the various wrapper routines below, without passing through the
! continuation library. For instance, the Jacobian matrix is never seen
! in the continuation routines, but is needed by several of the wrapper
! routines in this file. The passdown structure is a way to make the
! matrix global to this file.
   
   USE Loca_types

   IMPLICIT NONE
   
   TYPE(passdown_struct) :: pd

END MODULE Loca_pd

