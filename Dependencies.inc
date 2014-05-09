###############################################################################
# Objects
#
OBJ    = dirichlet_neumann.o         \
         dynamic_structures.o        \
         gauss_points_2d_l.o         \
         gauss_points_2d_p_sp.o      \
         prep_mesh_p1p2_sp.o         \
         qv_sp.o                     \
         qc_sp_m.o                   \
         qs_l_m.o                    \
         qs_l_sp.o                   \
         qs_sp.o                     \
         qs_sp_m.o                   \
         sparse_matrix_profiles.o    \
         sparse_matrix_operations.o  \
         global_variables.o          \
         miscellaneous_subroutines.o \
         start_sparse_kit.o          \
         fem_2d_grid.o               \
         par_solve_mumps.o           \
         axisym_boundary_values.o    \
         vtk_plot.o                  \
         restart_io.o                \
         read_input_files.o          \
         loca_wrappers.o             \
         loca_interface.o            \
         loca_parameters.o           \
         eigensolve.o                \
         transient_growth.o          \
         dns_algorithms.o            \
         vorticity_stream.o
           
###############################################################################
# Dependencies

# Fortran
main.o: $(OBJ)
   
read_input_files.o:
dirichlet_neumann.o:
dynamic_structures.o:
gauss_points_2d_l.o:
gauss_points_2d_p_sp.o:
fem_2d_grid.o:
sparse_matrix_profiles.o:
sparse_matrix_operations.o:
loca_parameters.o:
global_variables.o: read_input_files.o   \
                    dynamic_structures.o \
                    sparse_matrix_profiles.o

miscellaneous_subroutines.o:

par_solve_mumps.o: global_variables.o \
                   sparse_matrix_profiles.o

start_sparse_kit.o: sparse_matrix_profiles.o

qs_l_m.o: sparse_matrix_profiles.o \
          gauss_points_2d_p_sp.o   \
          gauss_points_2d_l.o

qs_l_sp.o: gauss_points_2d_p_sp.o

qs_sp.o: gauss_points_2d_p_sp.o

qs_sp_m.o: sparse_matrix_profiles.o \
           gauss_points_2d_p_sp.o

qv_sp.o: gauss_points_2d_p_sp.o

qc_sp_m.o: sparse_matrix_profiles.o \
           gauss_points_2d_l.o      \
           gauss_points_2d_p_sp.o

prep_mesh_p1p2_sp.o: gauss_points_2d_l.o    \
                     gauss_points_2d_p_sp.o \
                     dirichlet_neumann.o    \
                     fem_2d_grid.o

vtk_plot.o: miscellaneous_subroutines.o

restart_io.o: sparse_matrix_profiles.o \
              global_variables.o       \
              prep_mesh_p1p2_sp.o      \
              miscellaneous_subroutines.o

axisym_boundary_values.o: global_variables.o   \
                          dynamic_structures.o \
                          dirichlet_neumann.o  \
                          gauss_points_2d_p_sp.o

loca_wrappers.o: loca_parameters.o           \
                 par_solve_mumps.o           \
                 qv_sp.o                     \
                 qc_sp_m.o                   \
                 qs_l_sp.o                   \
                 dynamic_structures.o        \
                 sparse_matrix_profiles.o    \
                 global_variables.o          \
                 miscellaneous_subroutines.o \
                 prep_mesh_p1p2_sp.o         \
                 start_sparse_kit.o          \
                 dirichlet_neumann.o         \
                 read_input_files.o          \
                 eigensolve.o                \
                 axisym_boundary_values.o    \
                 vtk_plot.o

eigensolve.o: sparse_matrix_profiles.o   \
              sparse_matrix_operations.o \
              miscellaneous_subroutines.o

transient_growth.o: dynamic_structures.o        \
                    sparse_matrix_profiles.o    \
                    sparse_matrix_operations.o  \
                    global_variables.o          \
                    miscellaneous_subroutines.o \
                    prep_mesh_p1p2_sp.o         \
                    dirichlet_neumann.o         \
                    qc_sp_m.o                   \
                    qv_sp.o                     \
                    qs_l_m.o                    \
                    qs_l_sp.o                   \
                    dns_algorithms.o            \
                    par_solve_mumps.o           \
                    vtk_plot.o

dns_algorithms.o: dynamic_structures.o        \
                  sparse_matrix_profiles.o    \
                  global_variables.o          \
                  miscellaneous_subroutines.o \
                  prep_mesh_p1p2_sp.o         \
                  start_sparse_kit.o          \
                  qc_sp_m.o                   \
                  qv_sp.o                     \
                  par_solve_mumps.o           \
                  restart_io.o                \
                  vtk_plot.o

vorticity_stream.o: global_variables.o  \
                    prep_mesh_p1p2_sp.o \
                    dirichlet_neumann.o \
                    start_sparse_kit.o  \
                    qs_sp.o             \
                    qs_sp_m.o           \
                    par_solve_mumps.o

# C
loca_interface.o: loca_wrappers.o
