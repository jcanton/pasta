MODULE read_input_files

   IMPLICIT NONE

   TYPE program_input

      ! program data
      INTEGER            :: method
      INTEGER            :: nwtn_maxite
      REAL(KIND=8)       :: nwtn_tol
      CHARACTER(LEN=128) :: mesh_directory
      CHARACTER(LEN=128) :: mesh_name
      CHARACTER(LEN=128) :: plot_directory
      CHARACTER(LEN=128) :: restart_directory
      CHARACTER(LEN=128) :: input_restart_file
      CHARACTER(LEN=128) :: output_restart_file
      LOGICAL            :: read_restart_flag
      LOGICAL            :: write_restart_flag
      LOGICAL            :: write_QP_restart_flag
      LOGICAL            :: write_BVS_flag
      LOGICAL            :: write_plots_flag
      ! eigenvalue data
      INTEGER            :: eigen_BC
      INTEGER            :: eigen_nev
      INTEGER            :: eigen_maxit
      REAL(KIND=8)       :: eigen_tol
      COMPLEX(KIND=8)    :: eigen_sigma
      CHARACTER(LEN=128) :: eigen_output_directory
      INTEGER            :: eigen_plotNumber
      INTEGER            :: eigen_directAdjoint_flag
      LOGICAL            :: eigen_compute_structSens_flag
      ! structural sensitivity data
      INTEGER            :: structSens_eigenNumber
      CHARACTER(LEN=128) :: structSens_directEigen_name
      CHARACTER(LEN=128) :: structSens_adjointEigen_name
      ! transient growth data
      INTEGER            :: tranGrowth_method
      INTEGER            :: tranGrowth_initGuess
      INTEGER            :: tranGrowth_BC
      REAL(KIND=8)       :: tranGrowth_tau
      REAL(KIND=8)       :: tranGrowth_dt
      INTEGER            :: tranGrowth_maxit
      REAL(KIND=8)       :: tranGrowth_tol
      ! dns data
      INTEGER            :: dns_method
      REAL(KIND=8)       :: dns_tInit
      REAL(KIND=8)       :: dns_tEnd
      REAL(KIND=8)       :: dns_dt
      REAL(KIND=8)       :: dns_dtPlot
      CHARACTER(LEN=128) :: dns_output_directory
      ! evolve transient growth data
      CHARACTER(LEN=128) :: etg_opt_perturb

   END TYPE program_input

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE read_program_data(file_name,  prog_input)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 24/8/2013
!
! - file_name  :: file to be read
! - prog_input :: structured defined at the beginning of this module

   IMPLICIT NONE
   ! input variables
   CHARACTER(*), INTENT(IN) :: file_name
   ! output variables
   TYPE(program_input) :: prog_input
   ! local variables
   INTEGER :: fid = 22

   ! executable statements
   OPEN (UNIT = fid, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! program data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % method
   READ  (fid,*) prog_input % nwtn_maxite
   READ  (fid,*) prog_input % nwtn_tol
   READ  (fid,*) prog_input % mesh_directory,     prog_input % mesh_name
   READ  (fid,*) prog_input % plot_directory
   READ  (fid,*) prog_input % restart_directory
   READ  (fid,*) prog_input % input_restart_file, prog_input % output_restart_file
   READ  (fid,*) prog_input % read_restart_flag,  prog_input % write_restart_flag
   READ  (fid,*) prog_input % write_QP_restart_flag
   READ  (fid,*) prog_input % write_BVS_flag
   READ  (fid,*) prog_input % write_plots_flag
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! eigenvalue data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % eigen_BC
   READ  (fid,*) prog_input % eigen_nev
   READ  (fid,*) prog_input % eigen_maxit
   READ  (fid,*) prog_input % eigen_tol
   READ  (fid,*) prog_input % eigen_sigma
   READ  (fid,*) prog_input % eigen_output_directory
   READ  (fid,*) prog_input % eigen_plotNumber
   READ  (fid,*) prog_input % eigen_directAdjoint_flag
   READ  (fid,*) prog_input % eigen_compute_structSens_flag
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! structural sensitivity data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % structSens_eigenNumber
   READ  (fid,*) prog_input % structSens_directEigen_name
   READ  (fid,*) prog_input % structSens_adjointEigen_name
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! transient growth data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % tranGrowth_method
   READ  (fid,*) prog_input % tranGrowth_initGuess
   READ  (fid,*) prog_input % tranGrowth_BC
   READ  (fid,*) prog_input % tranGrowth_tau
   READ  (fid,*) prog_input % tranGrowth_dt
   READ  (fid,*) prog_input % tranGrowth_maxit
   READ  (fid,*) prog_input % tranGrowth_tol
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! dns data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % dns_method
   READ  (fid,*) prog_input % dns_tInit
   READ  (fid,*) prog_input % dns_tEnd
   READ  (fid,*) prog_input % dns_dt
   READ  (fid,*) prog_input % dns_dtPlot
   READ  (fid,*) prog_input % dns_output_directory
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) ! evolve transient growth data
   READ  (fid,*) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   READ  (fid,*) prog_input % etg_opt_perturb
   CLOSE(fid)


END SUBROUTINE read_program_data

!==============================================================================

END MODULE read_input_files
