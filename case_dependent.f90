MODULE case_dependent
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This module contains routines which are case dependent and can be adapted
! to suit the needs of a particular flow analysis.
! The idea is to have here ALL what is case dependent so that it's more
! difficult to make stupid memory mistakes to the core of the program.
!
   USE global_variables
   USE miscellaneous_subroutines
   USE axisym_boundary_values
   USE Gauss_points
   USE Gauss_points_L
   USE vorticity_stream
   USE vtk_plot


!------------------------------------------------------------------------------

   IMPLICIT NONE


CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE case_problemset()
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This routine is executed after reading the 'problem_data.in' file

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables
   REAL(KIND=8) :: Rt, Rp=0.5

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_problemset'

   !***torus
   WRITE(*,*) '    curvature = ', flow_parameters(1)
   WRITE(*,*) '    R_t       = ', Rp/flow_parameters(1)

   ! the mesh has R_t = 0.5
   rr_L(2,:) = rr_L(2,:) - 0.5d0 + Rp/flow_parameters(1)
   rr(2,:)   = rr(2,:)   - 0.5d0 + Rp/flow_parameters(1)
   CALL Gauss_gen_L(np_L, me, nps_L, mes, jj_L, jjs_L, rr_L)
   CALL Gauss_gen  (np,   me, nps,   mes, jj,   jjs,   rr)

!write(*,*) '*** check'
!write(*,*) '    maxR = ', MAXVAL(rr(2,:))
!write(*,*) '    minR = ', MINVAL(rr(2,:))


   WRITE(*,*) '    done: case_problemset'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_preprocess'
   WRITE(*,*) '    done: case_preprocess'
   WRITE(*,*)
END SUBROUTINE case_preprocess

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE case_newton_iteprocess(ite, continuation_converged)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This routine is called inside the nonlinear solver (Newton's method) before
! the generation of the right-hand side and Jacobian matrix

   IMPLICIT NONE
   ! input variables
   INTEGER, INTENT(IN) :: ite
   ! output variables
   INTEGER :: continuation_converged
   ! local variables
   REAL(KIND=8) :: Ub, Ub_tol
   REAL(KIND=8), DIMENSION(velCmpnnts) :: u_avg
   REAL(KIND=8), SAVE :: fn1, Ubn1, fn2, Ubn2

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_newton_iteprocess'

   !***torus
   Ub = 1d0
   Ub_tol = 1e-3
   !
   CALL computeFieldAverage(u0,  u_avg)
   fn1  = volumeForcing(3,2)
   Ubn1 = u_avg(3)

   IF ( ABS(Ub-Ubn1) > Ub_tol ) THEN
      continuation_converged=0
      WRITE(*,*) '    tol NOT satisfied'
   ELSE
      continuation_converged=1
      WRITE(*,*) '    tol satisfied'
   ENDIF
!   WRITE(*,*) '    Average velocity field'
!   WRITE(*,*) '    avg(u_z) = ', u_avg(1)
!   WRITE(*,*) '    avg(u_r) = ', u_avg(2)
   WRITE(*,*) '    avg(u_t) = ', u_avg(3)

   IF (ite>1) THEN
      ! secant method
      volumeForcing(3,2) = fn1 - (Ubn1-Ub)*(fn1-fn2)/(Ubn1-Ubn2)
!      write(*,*) 'secant', fn1, Ubn1, fn2, Ubn2
   ELSE
      volumeForcing(3,2) = fn1 + (Ub - Ubn1) * 0.2
!      write(*,*) 'stupid', fn1, Ubn1, fn2, Ubn2
   ENDIF

   fn2  = fn1
   Ubn2 = Ubn1

   WRITE(*,*) '    force    = ', volumeForcing(3,2)

   WRITE(*,*) '    done: case_newton_iteprocess'
   WRITE(*,*)
END SUBROUTINE case_newton_iteprocess

!------------------------------------------------------------------------------

SUBROUTINE case_newton_postprocess()
!
! Author:
! E-mail:
! Last revision:
!
! This routine is called after Newton's method has reached convergence

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables


   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_newton_postprocess'
   WRITE(*,*) '    done: case_newton_postprocess'
   WRITE(*,*)
END SUBROUTINE case_newton_postprocess

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
   INTEGER           :: fid = 23
   REAL(KIND=8), DIMENSION(velCmpnnts) :: u_avg

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis1'

   !***torus
   CALL computeFieldAverage(uu,  u_avg)
   WRITE(*,*)
   WRITE(*,*) '--> Average velocity field'
   WRITE(*,*) '    avg(u_z) = ', u_avg(1)
   WRITE(*,*) '    avg(u_r) = ', u_avg(2)
   WRITE(*,*) '    avg(u_t) = ', u_avg(3)

   OPEN(UNIT= fid, FILE='./locaOut/paramsUb.dat', ACCESS= 'APPEND')
   WRITE(fid,*) Re, flow_parameters(1), flow_parameters(2), flow_parameters(3), &
                u_avg(1), u_avg(2), u_avg(3), volumeForcing(3,2)
   CLOSE(fid)

   ! computation of vorticity and stream function
   ! as the boundary conditions are not imposed in a general form (yet),
   ! this part needs to be modified according to the geometry and
   ! BCs of the problem being solved
   ALLOCATE (Dir_psi(number_of_sides))
   ALLOCATE (zz(np), psi(np))
   Dir_psi = (/.TRUE./)
   CALL compute_vorticity_stream (jj, jjs, js, uu, rr, sides, Axis, Dir_psi,  zz, psi)
   CALL vtk_plot_scalar_P2 (rr, jj,  zz, trim(p_in%plot_directory) // 'steadyStateVorticity.vtk')
   CALL vtk_plot_scalar_P2 (rr, jj, psi, trim(p_in%plot_directory) // 'steadyStateStream.vtk')
   DEALLOCATE(Dir_psi)
   DEALLOCATE(zz, psi)



   WRITE(*,*) '    done: case_postprocess_analysis1'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis3'
   WRITE(*,*) '    done: case_postprocess_analysis3'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis4'
   WRITE(*,*) '    done: case_postprocess_analysis4'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis5'
   WRITE(*,*) '    done: case_postprocess_analysis5'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis6'
   WRITE(*,*) '    done: case_postprocess_analysis6'
   WRITE(*,*)
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
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_postprocess_analysis7'
   WRITE(*,*) '    done: case_postprocess_analysis7'
   WRITE(*,*)
END SUBROUTINE case_postprocess_analysis7

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE case_loca_solution_output(x_vec, filenm)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 24/11/2014
!
! This routine is executed when LOCA has converged on a steady solution
! and saves the results

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), DIMENSION(Nx) :: x_vec
   CHARACTER(*)                :: filenm
   ! output variables
   ! local variables
   REAL(KIND=8), DIMENSION(velCmpnnts,np) :: u_zzpsi

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_loca_solution_output'

   !***torus
!   ! computation of vorticity and stream function
!   ! as the boundary conditions are not imposed in a general form (yet),
!   ! this part needs to be modified according to the geometry and
!   ! BCs of the problem being solved
!   CALL extract(x_vec, u_zzpsi)
!   ALLOCATE (Dir_psi(number_of_sides))
!   ALLOCATE (zz(np), psi(np))
!   Dir_psi = (/.TRUE./)
!   CALL compute_vorticity_stream (jj, jjs, js, u_zzpsi, rr, sides, Axis, Dir_psi,  zz, psi)
!   CALL vtk_plot_scalar_P2 (rr, jj,  zz, trim(p_in%plot_directory)//'locaContVorticity'//trim(filenm)//'.vtk')
!   CALL vtk_plot_scalar_P2 (rr, jj, psi, trim(p_in%plot_directory)//'locaContStream'//trim(filenm)//'.vtk')
!   DEALLOCATE(Dir_psi)
!   DEALLOCATE(zz, psi)


   WRITE(*,*) '    done: case_loca_solution_output'
   WRITE(*,*)
END SUBROUTINE case_loca_solution_output

!------------------------------------------------------------------------------

SUBROUTINE case_loca_paramout()
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This routine is executed when LOCA has converged on a steady solution
! and saves the results

   IMPLICIT NONE
   ! input variables
   ! output variables
   ! local variables
   INTEGER           :: fid = 23
   REAL(KIND=8), DIMENSION(velCmpnnts) :: u_avg

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_loca_paramout'

   !***torus
   CALL extract (xx,  uu)
   CALL computeFieldAverage(uu,  u_avg)
!   WRITE(*,*) '    Average velocity field'
!   WRITE(*,*) '    avg(u_z) = ', u_avg(1)
!   WRITE(*,*) '    avg(u_r) = ', u_avg(2)
   WRITE(*,*) '    avg(u_t) = ', u_avg(3)

   OPEN(UNIT= fid, FILE='./locaOut/paramsUb.dat', ACCESS= 'APPEND')
   WRITE(fid,*) Re, flow_parameters(1), flow_parameters(2), flow_parameters(3), &
                u_avg(1), u_avg(2), u_avg(3), volumeForcing(3,2)
   CLOSE(fid)


   WRITE(*,*) '    done: case_loca_paramout'
   WRITE(*,*)
END SUBROUTINE case_loca_paramout

!------------------------------------------------------------------------------

SUBROUTINE case_loca_changeOscar(oscar)
!
! Author: Jacopo Canton
! E-mail: jcanton@mech.kth.se
! Last revision: 13/6/2014
!
! This routine is executed when LOCA needs to change parameter Oscar

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), INTENT(IN) :: oscar
   ! output variables
   ! local variables
   REAL(KIND=8) :: Rt, Rp=0.5

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_loca_changeOscar'

   !***torus
   WRITE(*,*) '    old curvature = ', flow_parameters(1)
   WRITE(*,*) '    old R_t       = ', Rp/flow_parameters(1)
   WRITE(*,*) '    new curvature = ', oscar
   WRITE(*,*) '    new R_t       = ', Rp/oscar

   rr_L(2,:) = rr_L(2,:) + Rp/oscar-Rp/flow_parameters(1)
   rr(2,:) = rr(2,:) + Rp/oscar-Rp/flow_parameters(1)
   CALL Gauss_gen_L(np_L, me, nps_L, mes, jj_L, jjs_L, rr_L)
   CALL Gauss_gen  (np,   me, nps,   mes, jj,   jjs,   rr)

write(*,*) '*** check'
write(*,*) '    dR   = ', Rp/oscar-Rp/flow_parameters(1)
write(*,*) '    maxR = ', MAXVAL(rr(2,:))
write(*,*) '    minR = ', MINVAL(rr(2,:))


   WRITE(*,*) '    done: case_loca_changeOscar'
   WRITE(*,*)
END SUBROUTINE case_loca_changeOscar

!------------------------------------------------------------------------------

SUBROUTINE case_loca_changeRomeo(romeo)
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed when LOCA needs to change parameter Romeo

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), INTENT(IN) :: romeo
   ! output variables
   ! local variables

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_loca_changeRomeo'
   WRITE(*,*) '    done: case_loca_changeRomeo'
   WRITE(*,*)
END SUBROUTINE case_loca_changeRomeo

!------------------------------------------------------------------------------

SUBROUTINE case_loca_changeWhisky(whisky)
!
! Author:
! E-mail:
! Last revision:
!
! This routine is executed when LOCA needs to change parameter Whisky

   IMPLICIT NONE
   ! input variables
   REAL(KIND=8), INTENT(IN) :: whisky
   ! output variables
   ! local variables

   ! executable statements
   WRITE(*,*)
   WRITE(*,*) '--> call to: case_loca_changeWhisky'
   WRITE(*,*) '    done: case_loca_changeWhisky'
   WRITE(*,*)
END SUBROUTINE case_loca_changeWhisky


!==============================================================================

END MODULE case_dependent
