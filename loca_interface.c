/*
 -----------------------------------------------------------------------------
 Interface between the FEM f90 solver and loca continuation algorithms.
 -----------------------------------------------------------------------------
 *
 * The continuation library currently consists of 9 files.
 * (1) [loca_interface.c]  The beam version of the interface
 *     to the continuation library. This file (and only this file)
 *     must be extensively modified for every new code that
 *     uses the library. It consists of do_loca, the
 *     top level routine called by the application code, and ~10
 *     wrapper routines, which must be filled for your specific
 *     code. The most noteworthy of these include
 *     nonlinear_solver_conwrap, linear_solver_conwrap,
 *     matrix_residual_fill_conwrap, and assign_parameter_conwrap.
 *     Each wrapper has its own comments.
 * (2) [loca_lib.c]  This is the stepping algorithm for zeroth-order,
 *     first-order, and arc-length continuation. The bifurcation
 *     tracking algorithms are only zeroth-order. This routine includes
 *     such things as step-size control and predictor steps. Ideally,
 *     this file does not need to be modified when linking to a new
 *     application code.
 * (3) [loca_bord.c] This file contains the bordering algorithms,
 *     currently for arc-length continuation, turning point tracking,
 *     pitchfork tracking, hopf tracking, and phase transition tracking.
 *     These routines are accessed from within the Newton iteration,
 *     and make frequent calls to the wrapper routines (matrix fill
 *     and solve). Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (4) [loca_const.h] This header file includes definitions of the
 *     continuation structures, define statements for flags used
 *     by the continuation library, and prototypes for the wrapper
 *     routines. Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (5) [loca_util.c] This file contains various utility routines,
 *     particularly vector operations. The vector operations in loca
 *     are all being migrated to this file, so that the rest of loca
 *     will become independent of the data structure for storing
 *     vectors. 
 * (6) [loca_util.h] Extern statements for loca_util.c functions.
 *(789)[loca_eigenvalue.c, loca_eigen_c2f.F, loca_eigen_cayley.F]
 *     Files containing routines for calculating leading eigenalues,
 *     including the Cayley transformation, interface to ARPACK and
 *     PARPACK, and post processing routines.
 *
 *
 * How to interface to this library:
 * (0) Have a steady-state code that uses Newton's method
 * (1) Call do_loca (top routine of this file) from your code, in
 *     the place where you normally call the steady-state or transient
 *     drivers.
 * (2) In your current nonlinear solver, add "(void *) con_ptr"
 *     to the argument list. Pass NULL in that place for all
 *     calls to the nonlinear solver besides the one from
 *     nonlinear_solver_conwrap below. In your Newton loop,
 *     after the matrix equation has been solved for the
 *     update vector, delta_x, but before the solution vector, x,
 *     has been updated by delta_x, put in the following call:
 *      if (con_ptr != NULL) continuation_converged =
 *         continuation_hook(x, delta_x, con_ptr, Reltol, Abstol);
 *     (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)
 *     When checking convergence of your Newton's method, also
 *     check that (continuation_converged == TRUE).
 * (3) Change the contents of all routines in this file. In
 *     solve_continuation, many fields of the con and con structures
 *     must be set. Follow the template and the comments.
 *     Also, the passdown structure can be defined and set with
 *     any information that needs to be passed from the
 *     solve_continuation routine to the wrapper routines below
 *     that isn't needed by the continuation library. All
 *     of the wrapper routines (all routines in this file besides
 *     solve_continuation, which all have names ending in _conwrap)
 *     must be filled in with the corresponding call to your
 *     application code. Each has its own comments.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/**********************************
 * WARNING:
 * passdown_struct is also defined
 * in loca_wrappers.f90
 * Update both files if changes are
 * needed
 **********************************/
struct passdown_struct {
   double *x;
   double reynolds, vratio, mu, alpha, h, tol; 
   int bif_param, param, maxiter, ldz, num_linear_its, debug;
   };

// extern pd;

/* This include file is for the continuation structures, defines, prototypes */
/* con_struct is defined in the library , so it cannot be modified.*/

#include "loca_const.h"

static void print_con_struct(const struct con_struct *con);
void cvarsparser(struct con_struct *con, struct passdown_struct *pd, const char* filename);
void print_c_line(const char *charstr, int ntimes);
static void print_final(const double param, const int step_num);




/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

//void do_loca(struct con_struct *con_input, struct problem_info *pi)
//
// ******************************************
// edit 1: remove the problem_info structure
// ******************************************
//
// void do_loca(struct con_struct *con_input)
//
// ******************************************
// edit 2: use the passdown_struc as input
//         and read 'con' in cvarsparser
// ******************************************
//
void do_loca(struct passdown_struct *pd)

{
  int nstep;
  struct con_struct con;


  printf("\n +++++++++++++++++++++++++++++++++++++\n");
  printf(" --> CALL to do_loca\n\n");



  /* read input file */
  cvarsparser(&con, pd, "loca_data.in");



  /* First load the general info structure */
  con.general_info.nv_restart      = FALSE;
  con.general_info.nv_save         = FALSE;
  con.general_info.perturb         = 1.0e-6;

  /* Then load the stepping info structure */
  con.stepping_info.base_step      = 0;

  /* Then load one of the method dependent structures */
  switch (con.general_info.method) {
    case ZERO_ORDER_CONTINUATION:
    case FIRST_ORDER_CONTINUATION:
      break;
    case ARC_LENGTH_CONTINUATION:
      // there is nothing here because con_input was identical to con
      break;
    case TURNING_POINT_CONTINUATION: 
      con.turning_point_info.nv = NULL;
      break;
    case PITCHFORK_CONTINUATION: 
      // there is nothing here because con_input was identical to con
      break; 
    case HOPF_CONTINUATION: 
      // there is nothing here because con_input was identical to con
      break;
    case MANIFOLD_CONTINUATION:
      printf("\nmanifold continuation not yet implemented\n");
      con.manifold_info.k =  2; /* number of parameters */
      con.manifold_info.param_vec          = NULL; /* vec of length k */
      con.manifold_info.param_lower_bounds = NULL; /* vec of length k */
      con.manifold_info.param_upper_bounds = NULL; /* vec of length k */
      con.manifold_info.param_steps        = NULL; /* vec of length k */
      break;
    case PHASE_TRANSITION_CONTINUATION:
      printf("\nphase transition continuation not yet implemented\n");
      break;
    default:
      printf("ERROR: Unknown LOCA input method: %d\n",con.general_info.method);
      exit(-1);
  }

  /* Finally, load the eigensolver structures */
  /* This is a temporary default: */
  con.eigen_info.Num_Eigenvectors  = con.eigen_info.Num_Eigenvalues;
  con.eigen_info.sort              = TRUE;

  /* print out continuation structure */
  print_con_struct(&con);

  /* CALL LOCA, now that all information is set */
  nstep = con_lib(&con);

  /* Final printing and cleanup here */

  if (con.general_info.printproc)
    print_final(con.general_info.param, nstep); //, pd.num_mat_fills, pd.num_res_fills, pd.num_linear_its);

  switch (con.general_info.method) {
    case PITCHFORK_CONTINUATION: 
      free((void *)con.pitchfork_info.psi);
      break;
    case HOPF_CONTINUATION: 
      free((void *)con.hopf_info.y_vec);
      free((void *)con.hopf_info.z_vec);
      break;
  }

  return;
  
} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void cvarsparser(struct con_struct *con, struct passdown_struct *pd,
                 const char* filename)
/* Input file reader code
 * The top sections are problem specific information
 * The bottom sections can be used by any LOCA application 
 *   since they read directly into the LOCA con structure
 */
{

  int i;
  int hopf_nev;
  int hopf_ind;
  int pitch_nev;
  int pitch_ind;
  double bif_param_init=0;
  char *temp;
  char *temp2;
  char *tmp;
  temp = (char *) malloc(100 * sizeof(char));
  temp2 = (char *) malloc(100 * sizeof(char));
  tmp = (char *) malloc(100 * sizeof(char));  
  FILE* data;
  int dummy; // just to avoid boring warnings in compilation


  extern void read_eigenvector (int, int, char*, int, double*, double*);

//  extern void pitch_input_arpack(double*,double*,char*,char*,int,int);

//  extern void pitch_input_loca(double*,double*,char*,char*);

//  extern void hopf_input_arpack(double*,double*, double*, char*,char*, int, int);

//  extern void hopf_input_loca(double*,double*, double*, char*,char*);

  data = fopen(filename, "r");
  if(data==NULL) { printf("\tCan't open file: %s\n", filename); exit(1); }

  /*... Problem specific Info ...*/
  
  //*******************************************
  // <!--Passdown-->
  dummy=fscanf(data, "%s\n", tmp);
  //*******************************************
//  dummy=fscanf(data, "%s\t%d\n",  temp, &pd->maxiter); // removed from here because
//  dummy=fscanf(data, "%s\t%le\n", temp, &pd->tol);     // they are read from program_data.in
  dummy=fscanf(data, "%s\t%d\n",  temp, &pd->debug);

  /*... General Info ...*/

  con->general_info.numUnks = pd->ldz;
  con->general_info.numOwnedUnks = pd->ldz;
  // pd->x = (double *) malloc(pd->ldz * sizeof(double)) ;
  con->general_info.x = pd->x;

  printf("Check to have correctly received the initial solution:\n");
  printf("\tpd->x[0]   = %20.10f \n",   pd->x[0]);
  printf("\tpd->x[end] = %20.10f \n\n", pd->x[pd->ldz-1]);
 
  //*******************************************
  // <!--General_Info_H5_P4_T3_A2_F1_Z0-->
  dummy=fscanf(data, "%s\n", tmp);
  //*******************************************
  //
  dummy=fscanf(data, "%s\t%d\n", temp, &con->general_info.method);
  //
  dummy=fscanf(data, "%s\t%d\n", temp, &pd->bif_param);
  //
  /* Assign initial value of bifurcation param based on flag value */
  //
  printf("Bifurcation parameter:\t");
  switch (pd->bif_param) {
    case REYNOLDS : bif_param_init = pd->reynolds; printf("Reynolds, "); break;
    case VRATIO   : bif_param_init = pd->vratio;   printf("vRatio, ");   break;
    case MU       : bif_param_init = pd->mu;       printf("mu, ");       break;
    case ALPHA    : bif_param_init = pd->alpha;    printf("alpha, ");    break;
  }
  printf("initial value = %5.4f\n", bif_param_init);
  //
  dummy=fscanf(data, "%s\t%d\n", temp, &pd->param);
  //
  printf("Continuation parameter:\t");
  switch (pd->param) {
    case REYNOLDS : con->general_info.param = pd->reynolds;  printf("Reynolds, ");break;
    case VRATIO   : con->general_info.param = pd->vratio;    printf("vRatio, ");  break;
    case MU       : con->general_info.param = pd->mu;        printf("mu, ");      break;
    case ALPHA    : con->general_info.param = pd->alpha;     printf("alpha, ");   break;
  }
  printf("initial value = %5.4f\n", con->general_info.param);
  //
  printf("Method:\t\t\t");
  switch (con->general_info.method){
     case ZERO_ORDER_CONTINUATION:
        printf("Zero-order continuation\n");
        break;
     case FIRST_ORDER_CONTINUATION:
        printf("First-order continuation\n");
        break;
     case ARC_LENGTH_CONTINUATION:
        printf("Arc length continuation\n");
        break;
     case TURNING_POINT_CONTINUATION:
        printf("Turning point continuation\n");
        con->turning_point_info.bif_param = bif_param_init;
        break;
     case PITCHFORK_CONTINUATION:
        printf("Pitchfork continuation\n");
        con->pitchfork_info.bif_param = bif_param_init;
        break;
     case HOPF_CONTINUATION:
        printf("Hopf continuation\n");
        con->hopf_info.bif_param = bif_param_init;
        break;
     case PHASE_TRANSITION_CONTINUATION:
        printf("Phase transition continuation\n");
        printf("\n+++\nWARNING\n+++\n");
        break;
     case AUGMENTING_CONDITION:
        printf("Augmenting condition\n");
        printf("\n+++\nWARNING\n+++\n");
        break;
     case MANIFOLD_CONTINUATION:
        printf("Manifold continuation\n");
        printf("\n+++\nWARNING\n+++\n");
        break;
     default:
        printf("Unknown method in input\nSTOP.\n\n");
        exit(-1);
  }


  /* For serial runs, this processor is always the printing processsor */
  /* con->general_info.printproc = 1; */

  /* temporary fix for turning off printing for debug=0 */
  /* if (pd->debug==0) con->general_info.printproc=0; */
  /* Doing it this way will allow the print level to be set (serial only) */
  pd->debug = 8;
  con->general_info.printproc = pd->debug;


  /*... Stepping Info ...*/

  //*******************************************
  // <!--Stepping_Info-->
  dummy=fscanf(data,"%s\n", tmp);
  //*******************************************
  dummy=fscanf(data,"%s\t%d\n", temp, &con->stepping_info.max_steps);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->stepping_info.max_param);
  printf("Continuation parameter goal value = %5.4f\n", con->stepping_info.max_param);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->stepping_info.first_step);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->stepping_info.step_ctrl);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->stepping_info.max_delta_p);
  con->stepping_info.max_newton_its = pd->maxiter;

  /*... Arclength Info ...*/

  //*******************************************
  // <!--Arclength_Info-->
  dummy=fscanf(data,"%s\n", tmp);
  //*******************************************
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->arclength_info.dp_ds2_goal);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->arclength_info.dp_ds_max);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->arclength_info.tang_exp);
  dummy=fscanf(data,"%s\t%lf\n", temp, &con->arclength_info.tang_step_limit);

  /*... Pitchfork Info ...*/

  //*******************************************
  // <!--Pitchfork_Info-->
  dummy=fscanf(data,"%s\n", tmp);
  //*******************************************
  dummy=fscanf(data,"%s\t%d\n", tmp, &pitch_nev);
  dummy=fscanf(data,"%s\t%d\n", tmp, &pitch_ind);
  dummy=fscanf(data,"%s\t%s\n", tmp, temp);
  dummy=fscanf(data,"%s\t%s\n", tmp, temp2);
  if (con->general_info.method==PITCHFORK_CONTINUATION) {
   con->pitchfork_info.psi = (double *) malloc(pd->ldz * sizeof(double));
   if (pitch_nev > 1) {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // pitch_input_arpack(pd->x,con->pitchfork_info.psi,temp,temp2,pitch_nev,pitch_ind);
    //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    printf("\tloaded psi vector (for pitchfork) from arpack sol. in directory:%s\n", temp2);
   }
   else {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // pitch_input_loca(pd->x,con->pitchfork_info.psi,temp,temp2);
    //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    printf("\tloaded psi vector (for pitchfork) from directory:%s\n", temp2);
    };
  };

  
  /*... Hopf Info ...*/

  //*******************************************
  // <!--Hopf_Info-->
  dummy=fscanf(data,"%s\n", tmp);
  //*******************************************
  dummy=fscanf(data,"%s\t%d\n", tmp, &hopf_nev);
  //dummy=fscanf(data,"%s\t%d\n", tmp, &hopf_ind);
  dummy=fscanf(data,"%s\t%s\n", tmp, temp);
  //dummy=fscanf(data,"%s\t%s\n", tmp, temp2);
  if (con->general_info.method==HOPF_CONTINUATION) {
   con->hopf_info.y_vec = (double *) malloc(pd->ldz * sizeof(double));
   con->hopf_info.z_vec = (double *) malloc(pd->ldz * sizeof(double)); 
   if (hopf_nev >= 1) {
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    read_eigenvector (pd->ldz, hopf_nev, temp, strlen(temp), con->hopf_info.y_vec, con->hopf_info.z_vec);
    printf("Check to have correctly received the eigenvector:\n");
    printf("\ty_vec[0] = %20.10f y_vec[end] = %20.10f \n",   con->hopf_info.y_vec[0], con->hopf_info.y_vec[pd->ldz-1]);
    printf("\tz_vec[0] = %20.10f z_vec[end] = %20.10f \n\n", con->hopf_info.z_vec[0], con->hopf_info.z_vec[pd->ldz-1]);
    //
    // hopf_input_arpack(pd->x,con->hopf_info.y_vec,con->hopf_info.z_vec,temp,temp2,hopf_nev,hopf_ind);
    // printf("\tloaded eigenvector (for hopf) from arpack sol. in file%s\n", temp);
    //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
   }
   else {
    printf("\n****************************************\n");
    printf("****************************************\n");
    printf("DISASTER jack didn't know what to do here");
    printf("\n****************************************\n");
    printf("****************************************\n");
   };
    
  };

  dummy=fscanf(data,"%s\t%lf\n", tmp, &con->hopf_info.omega);
  dummy=fscanf(data,"%s\t%d\n", tmp, &con->hopf_info.mass_flag);


  /*... Eigen Info ...*/

  //*******************************************
  // <!--Eigen_Info-->
  dummy=fscanf(data,"%s\n", tmp);
  //*******************************************
  dummy=fscanf(data,"%s\t%d\n", tmp, &con->eigen_info.Num_Eigenvalues); // Used as on/off flag for eigenvalues comp.
  dummy=fscanf(data,"%s\t%d\n", tmp, &con->eigen_info.Every_n_Steps);
  fclose(data);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void solution_output_conwrap(int num_soln_flag, double *x, double param,
                             double *x2, double param2,
                             double *x3, double param3,
                             int step_num, int num_its, struct con_struct *con)

/* Put the call to your solution output (both file and screen) routines here.
 * Input:
 *    num_soln_flag  Flag for how many solution vectors are being passed for
 *                   output. For parameter continuation and turning point
 *                   tracking there is just 1, for pitchfork and phase
 *                   transitions there are 2, and for Hopfs there are 3.
 *    x            First solution vector for output.
 *    param        Continuation parameter value.
 *    x2           Second solution vector for output (y_vec or x2)
 *    param2       Bifurcation parameter value.
 *    x3           Third solution vector for output (z_vec for Hopfs)
 *    param3       Third Parameter value (frequency Hopfs)
 *    step_num+1   Time index to output to (step_num is 0 based).
 *    num_its      Number of Newton iterations used for for convergence
 *    con          pointer to continuation structure, for passing to
 *                 the eigensolver.
 *
 * Output:
 *
 * Return Value:
 */
{
 
   extern void write_restart(double*, double, int, int, char*, int);

   extern void write_QP_restart(double*, char*, int);
  
   extern void param_output(double);

   extern void compute_eigen(double*, char*, int, double);
   
   extern void vtk_plot_loca(double*, char*, int);


   char filenm[64];
   double shiftIm=0.;


   printf("\n+++++++++++++++++++++++++++++++++++++\n");
   printf("--> CALL to solution_output_conwrap\n");

   /* Any Continuation run */

   if(num_soln_flag==1 || num_soln_flag==2 || num_soln_flag==3) {

      sprintf(filenm, "locaCont_%d_%d.dat", step_num, con->stepping_info.max_steps);
   
      write_restart(x, param, step_num, con->stepping_info.max_steps, filenm, strlen(filenm));

      sprintf(filenm, "suite_%d_%d.QPrestart", step_num, con->stepping_info.max_steps);

      write_QP_restart(x, filenm, strlen(filenm));
   
      /* Print out parameter to file */
   
      param_output(param);

    }

    /* 2-parameter bifurcation tracking runs */

    if(num_soln_flag==2 || num_soln_flag==3) {

      sprintf(filenm, "locaBifTrack_%d_%d.dat", step_num, con->stepping_info.max_steps);
      write_restart(x2, param2, step_num, con->stepping_info.max_steps, filenm, strlen(filenm));

    }

    /* For Hopf tracking, a third solution and parameter to write out */

    if(num_soln_flag==3) {

      sprintf(filenm, "locaHopf_%d_%d.dat", step_num, con->stepping_info.max_steps);
      write_restart(x3, param3, step_num, con->stepping_info.max_steps, filenm, strlen(filenm));

      shiftIm = param3;

    }


    printf("\n");

    // plot the solution
    sprintf(filenm, "_%d_%d", step_num, con->stepping_info.max_steps);
    vtk_plot_loca(x, filenm, strlen(filenm));

    if (con->eigen_info.Num_Eigenvalues > 0) {
      if (step_num == 0 || step_num%con->eigen_info.Every_n_Steps == 0) {
        sprintf(filenm, "_%d_%d", step_num, con->stepping_info.max_steps);
        compute_eigen(x, filenm, strlen(filenm), shiftIm);
      }
    }

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void eigenvector_output_conwrap(int j, int num_soln_flag, double *xr, double evr,
                                double *xi, double evi, int step_num)
/* Call to write out eigenvectors
 * Input:
 *    j    Eigenvector index/mode
 *    num_soln_flag  =1 for real vector, real eigenvalue
                     =2 for complex (imaginary part has info)
 *    xr   Real part of eigenvector
 *    evr  Real part of eigenvalue
 *    xi   Imaginary part of eigenvector (NULL if num_soln_flag==1)
 *    evi  Imaginary part of eigenvalue
 *    step_num  integer step number for use in output
 *
 * Output:
 *
 * Return Value:
 */
{
   /*
    ***** WARNING ***************
    * dummy subroutine because we
    * are using our own
    *****************************
    */
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_final(const double param, const int step_num) //, const int num_mat_fills,
              // const int num_res_fills, const int total_linear_its)

/*
 * Print out the final results and counters
 */

{
  printf("\n"); print_c_line("~", 80); 
  printf("\n"); print_c_line("~", 80); 
  printf("CONTINUATION ROUTINE HAS FINISHED: \n");
  printf("\tEnding Parameter value     = %g\n", param);
  printf("\tNumber of steps            = %d\n", step_num);
  //printf("\tNumber of Matrix fills     = %d\n", num_mat_fills);
  //printf("\tNumber of Residual fills   = %d\n", num_res_fills);
  //printf("\tNumber of linear solve its = %d\n", total_linear_its);
  printf("\n"); print_c_line("~", 80); 
  printf("\n"); print_c_line("~", 80); 
  printf("\n"); 

} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void writetofile(const double *x, const int n, const int m,
                 const char *filename)
{
  FILE* data; 
  int i,j,k=0;
  data = fopen(filename, "w"); 
  if(data==NULL) { printf("\tCan't open file: %s\n", filename); exit(1); }  
  for (i=0;i<n;i++) for(j=0;j<m;j++) { 
   fprintf(data, "%30.12lf\n", x[k]);
   k++;
  } 
  fclose(data);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void writetwofile(const double *x1, const double *x2,
                  const int n, const int m, const char *filename)
{
  FILE* data;
  int i,j,k=0;
  data = fopen(filename, "w");
  if(data==NULL) { printf("\tCan't open file: %s\n", filename); exit(1); }
  for (i=0;i<n;i++) for(j=0;j<m;j++) {
   fprintf(data, "%30.12g%30.12g\n", x1[k], x2[k]);
   k++;
  }
  fclose(data);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_con_struct(const struct con_struct* con)
/* Routine for printing out the con structure to the screen */

{
  printf("\n"); print_c_line("~", 80); 
  printf("Continuation structure:\n\n");

  printf("\tcon->general_info.param=             %10.4g\n",
    con->general_info.param);  
  printf("\tcon->general_info.numUnks=           %10d\n",
    con->general_info.numUnks);  
  printf("\tcon->general_info.numOwnedUnks=      %10d\n",
    con->general_info.numOwnedUnks);  
  printf("\tcon->general_info.printproc=         %10d\n",
    con->general_info.printproc);  

  printf("\tcon->stepping_info.first_step=       %10.4g\n",
    con->stepping_info.first_step);  
  printf("\tcon->stepping_info.max_steps=        %10d\n",
    con->stepping_info.max_steps);  
  printf("\tcon->stepping_info.max_param=        %10.4g\n",
    con->stepping_info.max_param);  
  printf("\tcon->stepping_info.max_delta_p=      %10.4g\n",
    con->stepping_info.max_delta_p);  
  printf("\tcon->stepping_info.step_ctrl=        %10.4g\n",
    con->stepping_info.step_ctrl);  
  printf("\tcon->stepping_info.max_newton_its=   %10d\n",
    con->stepping_info.max_newton_its);  

  if(con->general_info.method==ARC_LENGTH_CONTINUATION) {
    printf("\tcon->arclength_info.dp_ds2_goal=     %10.4g\n",
      con->arclength_info.dp_ds2_goal);
    printf("\tcon->arclength_info.dp_ds_max=       %10.4g\n",
      con->arclength_info.dp_ds_max);
    printf("\tcon->arclength_info.tang_exp=        %10.4g\n",
      con->arclength_info.tang_exp);
    printf("\tcon->arclength_info.tang_step_limit= %10.4g\n",
      con->arclength_info.tang_step_limit);
  }

  if(con->general_info.method==TURNING_POINT_CONTINUATION) 
    printf("\tcon->turning_point_info.bif_param=   %10.4g\n",
      con->turning_point_info.bif_param);
 
  if(con->general_info.method==PITCHFORK_CONTINUATION) 
    printf("\tcon->pitchfork_info.bif_param=       %10.4g\n",
      con->pitchfork_info.bif_param);
 
  if(con->general_info.method==HOPF_CONTINUATION) 
    printf("\tcon->hopf_info.bif_param=            %10.4g\n",
      con->hopf_info.bif_param);
 
  if(con->eigen_info.Num_Eigenvalues > 0) {
    printf("\tcon->eigen_info.Every_n_Steps=       %10d\n",
      con->eigen_info.Every_n_Steps);
  }
  print_c_line("~", 80); printf("\n"); 
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void print_c_line (const char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++) printf("%c", *charstr);
  printf("\n");
}

/**********************END OF INPUT/OUTPUT WRAPPERS***************************/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_multi_parameter_conwrap(double *param_vec)

   // OK

/* Put the call to a routine to assign the continuation parameters here.
 * Input:
 *    param_vec     New values of continuation parameters.
 *
 * Output:
 *
 * Return Value:
 */
{
  printf("multi_param not implemented in test_driver\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void calc_scale_vec_conwrap(double *x, double *scale_vec, int numUnks)

   // OK

/* Put the call to a routine to calculate a scaling vector here.
 * Input:
 *    x          New value of continuation parameter.
 *    numUnks    Number of unknowns on this proc, the length of x
 *               and scale_vec.
 *
 * Output:
 *    scale_vec  Vector of length number of unknowns used to scale
 *               variables so that one type of unknown (e.g. pressure)
 *               doesn't dominate over others. Used to balance the
 *               variables and the arc-length variable in arc-length
 *               continuation, and for scaling the NULL vector in
 *               turning point tracking. Using reciprocal of the average
 *               value of that variable type is a good choice. Vector
 *               of all ones should suffice for most problems.
 *
 * Return Value:
 */
{
  int i;
  for (i=0;i<numUnks;i++) scale_vec[i] = 1.0; 
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double gsum_double_conwrap(double sum)

   // OK

/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    sum     Value of double on this processor to be summed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global sum is returned on all processors.
 */
{
  return sum;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int gmax_int_conwrap(int max)

   // OK

/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    max     Value of integer on this processor to be maxed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global max is returned on all processors.
 *
 * Only used by Eigensolver
 */
{
  return max;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void random_vector_conwrap(double *x, int numOwnedUnks)

   // OK

/* Put a routine to calculate a random vector.
 * Input:
 *    numOwnedUnks  Length of owned nodes part of x.
 *
 * Output:
 *    x             Random vector.
 *
 * Used by eigensolver only
 */
{
  static int iseed[4]={1422,0432,3432,2131};

  dlaruv_(iseed, &numOwnedUnks, x);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void perturb_solution_conwrap(double *x, double *x_old,
                    double *scale_vec, int numOwnedUnks)

   // OK

/* Put a routine to perturb the solution vector a little bit here.
 * This is to move a converged solution at a singularity off
 * of the singularity before doing continuation. This ain't pretty
 * but has helped convergence on some turning point tracking problems.
 * Input:
 *    x_old         Current solution vector.
 *    scale_vec     Work space for a vector to scale x.
 *    numOwnedUnks  Length of owned nodes part of x, x_old, scale_vec
 *
 * Output:
 *    x             Solution vector perturbed a bit.
 *
 * Return Value:
 */
{
  int i;

  random_vector_conwrap(x, numOwnedUnks);

  for(i=0;i<numOwnedUnks;i++) {
    x[i] = x_old[i] * (1.0 + 1.0e-4*(x[i]-0.5));
  } 
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double free_energy_diff_conwrap(double *x, double *x2)

   // OK

/* Call to return the free energy difference betwen two solutions
 * Input:
 *    x    One solution vector
 *    x2   Second solution vector
 *
 * Output:
 *
 * Return Value:
 *    The difference in the free energy beween the two solutions
 */
{
  return 1.0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
