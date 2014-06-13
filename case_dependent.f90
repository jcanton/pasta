MODULE case_dependent
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This module contains routines which are case dependent and can be adapted
! to suit the needs of a particular flow analysis.
!
   USE global_variables
   USE miscellaneous_subroutines

!------------------------------------------------------------------------------

   IMPLICIT NONE


CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE case_problemset()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after reading the 'problem_data.in' file

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

   !***coaxial jets
   IF (flow_parameters(1) /= 0d0) THEN
      WRITE(*,*) '    assigning velocity ratio based on header input'
      in_bvs_D(1,1,2) = 1d0
      in_bvs_D(1,5,2) = flow_parameters(1)
   ENDIF

END SUBROUTINE case_problemset

!------------------------------------------------------------------------------

SUBROUTINE case_preprocess()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed before the analysis starts but after the program
! initialization

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables
   ! executable statements

END SUBROUTINE case_preprocess

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis1()
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This routine is executed after analysis 1:
! Steady state computation

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts) :: u_avg

   ! executable statements
   CALL computeFieldAverage(uu,  u_avg)
   WRITE(*,*)
   WRITE(*,*) '--> Average velocity field'
   WRITE(*,*) '    avg(u_z) = ', u_avg(1)
   WRITE(*,*) '    avg(u_r) = ', u_avg(2)
   WRITE(*,*) '    avg(u_t) = ', u_avg(3)

END SUBROUTINE case_postprocess_analysis1

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis3()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after analysis 3:
! Eigenvalue computation on an already computed base flow

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

END SUBROUTINE case_postprocess_analysis3

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis4()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after analysis 4:
! Structural sensitivity analysis on an already computed base flow, and already
! computed both direct and adjoint eigenvectors

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

END SUBROUTINE case_postprocess_analysis4

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis5()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after analysis 5:
! Transient growth computation on an already computed base flow

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

END SUBROUTINE case_postprocess_analysis5

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis6()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after analysis 6
! DNS

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

END SUBROUTINE case_postprocess_analysis6

!------------------------------------------------------------------------------

SUBROUTINE case_postprocess_analysis7()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed after analysis 7
! Evolution of optimal linear perturbations an already computed base flow

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables

   ! executable statements

END SUBROUTINE case_postprocess_analysis7

!------------------------------------------------------------------------------


!==============================================================================

END MODULE case_dependent
