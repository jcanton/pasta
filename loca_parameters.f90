MODULE Loca_parameters

   IMPLICIT NONE

   !**************************************************************************!
   !*                       DEFINE STATEMENTS                                *!
   !**************************************************************************!

   !**********************************
   !* WARNING:
   !* All of this definitions are also
   !* present in loca_const.h
   !* Update both files if changes are
   !* needed
   !**********************************


   INTEGER, PARAMETER ::  TRUE  = 1
   INTEGER, PARAMETER ::  FALSE = 0

   ! Choices for continuation method
   INTEGER, PARAMETER ::  ZERO_ORDER_CONTINUATION       = 0
   INTEGER, PARAMETER ::  FIRST_ORDER_CONTINUATION      = 1
   INTEGER, PARAMETER ::  ARC_LENGTH_CONTINUATION       = 2
   INTEGER, PARAMETER ::  TURNING_POINT_CONTINUATION    = 3
   INTEGER, PARAMETER ::  PITCHFORK_CONTINUATION        = 4
   INTEGER, PARAMETER ::  HOPF_CONTINUATION             = 5
   INTEGER, PARAMETER ::  PHASE_TRANSITION_CONTINUATION = 6
   INTEGER, PARAMETER ::  AUGMENTING_CONDITION          = 7
   INTEGER, PARAMETER ::  MANIFOLD_CONTINUATION         = 8
   INTEGER, PARAMETER ::  LOCA_LSA_ONLY                 = 9

   ! Choices for matrix fill -- what quantities to calculate
   INTEGER, PARAMETER ::  RHS_ONLY        = 100
   INTEGER, PARAMETER ::  MATRIX_ONLY     = 101
   INTEGER, PARAMETER ::  RHS_MATRIX      = 102
   INTEGER, PARAMETER ::  RHS_MATRIX_SAVE = 103
   INTEGER, PARAMETER ::  RECOVER_MATRIX  = 104

   ! Choices for linear solve about the state of the Jacobian matrix
   INTEGER, PARAMETER ::  NEW_JACOBIAN               = 200
   INTEGER, PARAMETER ::  OLD_JACOBIAN               = 201
   INTEGER, PARAMETER ::  SAME_BUT_UNSCALED_JACOBIAN = 202
   INTEGER, PARAMETER ::  OLD_JACOBIAN_DESTROY       = 203
   INTEGER, PARAMETER ::  CHECK_JACOBIAN             = 204

   ! Internal flags for rhs calculation for continuation linear solves
   INTEGER, PARAMETER ::  CONT_TANGENT     = 300
   INTEGER, PARAMETER ::  ARC_CONT_SOL2    = 301
   INTEGER, PARAMETER ::  TP_CONT_SOL2     = 302
   INTEGER, PARAMETER ::  TP_CONT_SOL3     = 303
   INTEGER, PARAMETER ::  TP_CONT_SOL4     = 304
   INTEGER, PARAMETER ::  HP_CONT_SOL3     = 305
   INTEGER, PARAMETER ::  HP_CONT_DMDX     = 306
   INTEGER, PARAMETER ::  HP_CONT_DMDPARAM = 307

   ! Problem parameters
   INTEGER, PARAMETER ::  REYNOLDS = 0
   INTEGER, PARAMETER ::  OSCAR    = 1
   INTEGER, PARAMETER ::  ROMEO    = 2
   INTEGER, PARAMETER ::  WHISKY   = 3


END MODULE Loca_parameters
