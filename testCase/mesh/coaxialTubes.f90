program coaxialTubes

!==============================================================================
! geometry definitions
!
! (X) : border number
!                     _____________________L3_(8)______________________
!                     |                                                |
!                     L2 (7)                                           |
!                     |                                                |
!      _______L1_(6)__|                                                |
!      |                                                               |
!      H (5)                                                           |
!      |_______(4)____                                                 |
!      _______________|S (3)                                           L4 (9)
!      |       (2)                                                     |
!      |                                                               |
!      R (1)                                                           |
!      |                                                               |
!  __ .| __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __| . __ 
!                     (0,0)        axis of symmetry (10)
!
!==============================================================================


!==============================================================================
! variables definition

  IMPLICIT NONE
   ! geometry variables
   REAL(KIND=8)   :: R          ! central jet radius
   REAL(KIND=8)   :: S          ! distance between jets
   REAL(KIND=8)   :: H          ! lateral jet height
   REAL(KIND=8)   :: L1         ! tubes length
   REAL(KIND=8)   :: L2         ! mesh border height = L4 - R - S - H
   REAL(KIND=8)   :: L3         ! mesh length
   REAL(KIND=8)   :: L4         ! mesh height
   REAL(KIND=8)   :: Roi        ! outer jet inner radius
   REAL(KIND=8)   :: Roo        ! outer jet outer radius
   ! mesh variables
   REAL(KIND=8)   :: hTb        ! mesh characteristic dimension on the borders of the tubes
   REAL(KIND=8)   :: hTc        ! mesh characteristic dimension on the centerline of the tubes
   REAL(KIND=8)   :: BL         ! boundary layer in the tubes
   REAL(KIND=8)   :: BLo        ! "outer" boundary layer in the tubes
   REAL(KIND=8)   :: hA         ! mesh characteristic dimension on the axis
   REAL(KIND=8)   :: hB         ! mesh characteristic dimension on the borders
   REAL(KIND=8)   :: hS         ! mesh characteristic dimension on the wall between the tubes
   INTEGER        :: nS         ! number of elements on the wall between the tubes
   REAL(KIND=8)   :: Lf         ! refinement segment
   REAL(KIND=8)   :: alphaWakeU ! angle of the wake behind the wall, upper segment (-90 turns it off) [deg]
   REAL(KIND=8)   :: alphaWakeL ! angle of the wake behind the wall, lower segment, assigned
   REAL(KIND=8)   :: expWake    ! exponent for the characteristic dimension of the mesh in the wake
   ! miscellaneous
   INTEGER        :: laplSmooth  ! laplacian smoothing Y/N?
   REAL(KIND=8)   :: x, dxW, dxT ! used to build the refinement areas
   INTEGER        :: i, n_pointsW, n_pointsT ! idem
   INTEGER        :: fid1 = 100, fid2 = 101, fid3 = 102


!==============================================================================
! read input file 'INPUTgeometry'

   OPEN(UNIT=fid1, FILE='./INPUTgeometry', FORM='formatted', STATUS='old')

   READ(fid1,*) R
   READ(fid1,*) S
   READ(fid1,*) H
   READ(fid1,*) L1
   READ(fid1,*) L4
   READ(fid1,*) L3
   READ(fid1,*) ! jump one line
   READ(fid1,*) hTb
   READ(fid1,*) hTc
   READ(fid1,*) BL
   READ(fid1,*) hA
   READ(fid1,*) hB
   READ(fid1,*) nS
   READ(fid1,*) Lf
   READ(fid1,*) alphaWakeU
   READ(fid1,*) expWake
   READ(fid1,*) ! jump one line
   READ(fid1,*) laplSmooth

   WRITE(*,*) 'R = ', R
   WRITE(*,*) 'S = ', S
   WRITE(*,*) 'H = ', H
   WRITE(*,*)

   CLOSE(fid1)

   ! short calculations
   L2  = L4 - R - S - H
   Roi = R + S
   Roo = R + S + H
   hS  = S/nS

   BLo       = 0.5
   n_pointsT = (L1-Lf)/hTb
   dxT       = (L1-Lf)/n_pointsT

   alphaWakeL = (R-S)/L3
   n_pointsW  = L3/0.5
   dxW        = L3/n_pointsW

!==============================================================================
! write on file 'boundary.coaxialTubes'

   OPEN(UNIT=fid2, FILE='./boundary.coaxialTubes')

   WRITE(fid2,*) 'NV      NE'
   WRITE(fid2,*) '10      10'
   WRITE(fid2,*) '#'

! geometry
   WRITE(fid2,*) '######## GEOMETRY###########'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 1'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', -L1, ' 0.0 ', -L1, R
   WRITE(fid2,*) 'END 1'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 2'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', -L1, R, ' 0.0 ', R
   WRITE(fid2,*) 'END 2'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 3'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 0.0 ', R, ' 0.0 ', Roi
   WRITE(fid2,*) 'END 3'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 4'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 0.0 ', Roi, -L1, Roi
   WRITE(fid2,*) 'END 4'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 5'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', -L1, Roi, -L1, Roo
   WRITE(fid2,*) 'END 5'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 6'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', -L1, Roo, ' 0.0 ', Roo
   WRITE(fid2,*) 'END 6'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 7'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 0.0 ', Roo, ' 0.0 ', L4
   WRITE(fid2,*) 'END 7'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 8'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 0.0 ', L4, L3, L4
   WRITE(fid2,*) 'END 8'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 9'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', L3, L4, L3, ' 0.0'
   WRITE(fid2,*) 'END 9'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 10'
   WRITE(fid2,*) 'line'
   WRITE(fid2,*) '1000 ', L3, ' 0.0 ', -L1, ' 0.0'
   WRITE(fid2,*) 'END 10'
   WRITE(fid2,*) '#'

! topology
   WRITE(fid2,*) '######## TOPOLOGY ###########'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 1'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       1       4       1       2'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hTc
   WRITE(fid2,*) '       1 ',  (R-BL)*(1-BLo)/R, hTc
   WRITE(fid2,*) '       1 ',  (R-BL)/R, hTb
   WRITE(fid2,*) '       1     1.0 ',    hTb
   WRITE(fid2,*) 'END 1'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 2'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       2       3       2       3'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hTb
   WRITE(fid2,*) '       1 ',  1-Lf/L1,  hS
   WRITE(fid2,*) '       1     1.0 ',    hS
   WRITE(fid2,*) 'END 2'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 3'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       3       2       3       4'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hS
   WRITE(fid2,*) '       1     1.0 ',    hS
   WRITE(fid2,*) 'END 3'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 4'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       4       3       4       5'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hS
   WRITE(fid2,*) '       1 ',  Lf/L1,    hS
   WRITE(fid2,*) '       1     1.0 ',    hTb
   WRITE(fid2,*) 'END 4'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 5'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       5       6       5       6'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hTb
   WRITE(fid2,*) '       1 ',  BL/H,     hTb
   WRITE(fid2,*) '       1 ',  (BL+(H-2*BL)/2*BLo)/H,   hTc
   WRITE(fid2,*) '       1 ',  1-(BL+(H-2*BL)/2*BLo)/H, hTc
   WRITE(fid2,*) '       1 ',  1-BL/H,   hTb
   WRITE(fid2,*) '       1     1.0 ',    hTb
   WRITE(fid2,*) 'END 5'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 6'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       6       3       6       7'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hTb
   WRITE(fid2,*) '       1 ',  1-Lf/L1,  hS
   WRITE(fid2,*) '       1     1.0 ',    hS
   WRITE(fid2,*) 'END 6'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 7'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       7       3       7       8'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hS
   WRITE(fid2,*) '       1 ',  Lf/L2,    hS
   WRITE(fid2,*) '       1     1.0 ',    hB
   WRITE(fid2,*) 'END 7'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 8'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       8       2       8       9'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hB
   WRITE(fid2,*) '       1     1.0 ',    hB
   WRITE(fid2,*) 'END 8'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 9'
   WRITE(fid2,*) '       e       n      bv      ev'
   WRITE(fid2,*) '       9       3       9      10'
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',    hB
   WRITE(fid2,*) '       1 ',  1-R/L4,   hA
   WRITE(fid2,*) '       1     1.0 ',    hA
   WRITE(fid2,*) 'END 9'
   WRITE(fid2,*) '#'
    WRITE(fid2,*) 'BEGIN 10'
   WRITE(fid2,*) '       e       n      bv      ev'
    IF ( alphaWakeU /= -90 ) THEN
       WRITE(fid2,*) '       10 ',n_pointsW-1+3,'      10       1' 
    ELSE
       WRITE(fid2,*) '       10      3      10       1' 
    ENDIF
   WRITE(fid2,*) '       i       s       h'
   WRITE(fid2,*) '       1     0.0 ',           hA
    IF ( alphaWakeU /= -90 ) THEN
      DO i = 1, n_pointsW-1
         x = i*dxW
         WRITE(fid2,*) '       1 ', x/(L1+L3), hA - (hA-hS)*(x/L3)**(1/expWake)
      ENDDO
    ENDIF
   WRITE(fid2,*) '       1 ',  L3/(L1+L3), hS
   WRITE(fid2,*) '       1     1.0 ',      hTc
   WRITE(fid2,*) 'END 10'
   WRITE(fid2,*) '#'

   CLOSE(fid2)

   WRITE(*,*) 'Done writing boundary.coaxialTubes'

!==============================================================================
! write on file 'domain.coaxialTubes'

   OPEN(UNIT=fid3, FILE='./domain.coaxialTubes')

   WRITE(fid3,*) '   nstep   lsmooth(1/0)'
   WRITE(fid3,*) '   10000   ', laplSmooth
   WRITE(fid3,*) '#'

!------------------------------------------
! Backgrid points

   WRITE(fid3,*) '   n   (backgrid points)'

   IF ( alphaWakeU == -90 ) THEN
      WRITE(fid3,*) 6*(n_pointsT-1)
   ELSE
      WRITE(fid3,*) 6*(n_pointsT-1) + 3*(n_pointsW-1)
   ENDIF

   WRITE(fid3,*) ' x y data'
   
   ! tubes points
   DO i = 1, n_pointsT-1
      x = i*dxT
      ! inner tube
      WRITE(fid3,*) -L1 + x, (R - BL) * (1-BLo), hTc + (hS-hTc)*(x/(L1-Lf))
      WRITE(fid3,*) -L1 + x,  R - BL,            hTb + (hS-hTb)*(x/(L1-Lf))
      ! outer tube
      WRITE(fid3,*) -L1 + x, Roi + BL,                  hTb + (hS-hTb)*(x/(L1-Lf))
      WRITE(fid3,*) -L1 + x, Roi + BL + (H-2*BL)/2*BLo, hTc + (hS-hTc)*(x/(L1-Lf))
      WRITE(fid3,*) -L1 + x, Roo - BL - (H-2*BL)/2*BLo, hTc + (hS-hTc)*(x/(L1-Lf))
      WRITE(fid3,*) -L1 + x, Roo - BL,                  hTb + (hS-hTb)*(x/(L1-Lf))
   ENDDO

   ! wake points
   IF ( alphaWakeU /= -90 ) THEN
      alphaWakeU = TAN(alphaWakeU * ACOS(-1.0) / 180)
      IF (alphaWakeU > (L4-(R+1.5*S))/L3 ) THEN
         alphaWakeU = (L4-(R+1.5*S))/L3
      ENDIF

      DO i = 1, n_pointsW-1
         x = i*dxW
         !WRITE(fid3,*) x, R   - Lf - x*alphaWakeL, hS + (hB-hS)*(x/L3)**expWake ! lower segment
         WRITE(fid3,*) x, R   - Lf,                hS + (hA-hS)*(x/L3)**expWake ! lower segment
         WRITE(fid3,*) x, Roi + Lf + x*alphaWakeU, hS + (hA-hS)*(x/L3)**expWake ! center lower segment
         !WRITE(fid3,*) x, Roo - Lf - x*alphaWakeL, hS + (hB-hS)*(x/L3)**expWake ! center upper segment
         WRITE(fid3,*) x, Roo + Lf + x*alphaWakeU, hS + (hA-hS)*(x/L3)**expWake ! upper segment
      ENDDO
   ENDIF

   WRITE(fid3,*) '#'
!------------------------------------------
! Steiner points
   WRITE(fid3,*) '   n   (steiner points)'
   WRITE(fid3,*) '   0'
   WRITE(fid3,*) ' x y'

   CLOSE(fid3)


!==============================================================================
! end of program

   WRITE(*,*) 'ALL DONE.'

end program coaxialTubes
