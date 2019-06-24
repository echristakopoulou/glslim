/**************************************************************/
/*! \file main_glslimr0.c
    \brief This is the file containing the main entry for GLSLIMr0.

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2016
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*! \brief The main entry for learning in gu-mode. Gu-mode is 
  an iterative mode, in which the clustering assignment of the
  users remains fixed, but the gu vector is updated in every 
  iteration, along with the model learnt.
 */
/**************************************************************/
int main(int argc, char *argv[])
{

  ctrl_t *ctrl;
  double *g;
  int *participation;
  gk_csr_t *train, *test;
  char model_bin_string[150] = "";
  char gu_string[150] = "";
  char gu_readable_string[150] = "";
  char initial_participation_string[150] = "";
  size_t nrows;
  int my_id, my_num_procs;
  gk_csr_t *model = NULL;
  int iter;
   gk_csr_t *prev_model = NULL;
  FILE *fpout = NULL;
  double error = 0;
  double objective = 0;
  double *norms = NULL;
  double objective_diff = 1;
  double old_objective = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &my_num_procs);

  srand(0);

  /* parse command line */
  ctrl = create_ctrl();
  parse_cmdline(ctrl, argc, argv);

  ctrl->id = my_id;
  ctrl->num_procs = my_num_procs;

  /* Reading train and test */
  train = gk_csr_Read(ctrl->train_file, GK_CSR_FMT_CSR, 1, 1);
  test = gk_csr_Read(ctrl->test_file, GK_CSR_FMT_CSR, 1, 1);

  /* read participation file */
  nrows = (size_t) train->nrows + 1;
  sprintf(initial_participation_string, "%s_v0", ctrl->participation_file);
  participation = gk_i32readfile(initial_participation_string, &nrows);

  /* read gu file */
  sprintf(gu_string, "%s_v%d", ctrl->gu_file, ctrl->start_iteri);
  g = gk_dreadfilebin(gu_string, &nrows);

  /* If it is not the 1st iteration, read the previous model 
     If it is, use the mode learnt with different regularization. */
  if (ctrl->start_iteri != 0) {
    sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file,
	    ctrl->start_iteri - 1);
    prev_model = gk_csr_Read(model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
    if ((ctrl->size) * (ctrl->num_clusters + 1) > prev_model->ncols)
      prev_model->ncols = (ctrl->size) * (ctrl->num_clusters + 1);
    gk_csr_CreateIndex(prev_model, GK_CSR_COL);
  }
  else {  
    if (ctrl->prev_model_file) {
      sprintf(model_bin_string, "%s_bin_v0", ctrl->prev_model_file);
      prev_model = gk_csr_Read(model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
      if ((ctrl->size) * (ctrl->num_clusters + 1) > prev_model->ncols)
	prev_model->ncols = (ctrl->size) * (ctrl->num_clusters + 1); 
      gk_csr_CreateIndex(prev_model, GK_CSR_COL);
    }
  }

  iter = ctrl->start_iteri;

  while ((iter < ctrl->end_iteri) && (objective_diff > ctrl->threshold)) {  
    /* Learn the model for this iteration and then save it as 
       prev_model for the next iteration */
    
    model = learn(ctrl, train, participation, g, prev_model);
    
    if (prev_model != NULL) {
      gk_csr_Free(&prev_model);
    }
    prev_model = model;
    
    /* Column indexing the model */
    if ((ctrl->size) * (ctrl->num_clusters + 1) > model->ncols)
      model->ncols = (ctrl->size) * (ctrl->num_clusters + 1);
    gk_csr_CreateIndex(model, GK_CSR_COL);
    
    /* Compute training error */
    error = train_predict(ctrl, train, model, g, participation);

    /* Compute frobenius and l1 norm of the model 
       and value of objective function */
    if (ctrl->id == 0) {
      norms = norm(ctrl, model);
      fpout = gk_fopen(ctrl->stats_file, "a", "stats open");
      objective = norms[0] + norms[1] + error;
      fprintf(fpout, "error %f frob %f l1 %f obj %f\n",
	      error, norms[0], norms[1], objective);
      gk_fclose(fpout);
      old_objective = objective;
    }
    
    /* Compute HR and ARHR */
    slim_test(ctrl, model, g, train, test, participation);

    /* Learn gu for this iteration */
    gk_free((void **) &g, LTERM);
    g = learn_gu_all(ctrl, train, model, participation);
    
    /* Compute training error, objective value and HR/ARHR 
       for the same model with the new gu */
      error = train_predict(ctrl, train, model, g, participation);

    if (ctrl->id == 0) {
      fpout = gk_fopen(ctrl->stats_file, "a", "stats open");
      objective = norms[0] + norms[1] + error;
      fprintf(fpout, "error %f frob %f l1 %f obj %f\n",
	      error, norms[0], norms[1], objective);
      gk_fclose(fpout);
      objective_diff = (old_objective - objective) / objective;
    }

    MPI_Bcast(&objective_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    slim_test(ctrl, model, g, train, test, participation);

    if (ctrl->id == 0){
      sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file, iter);
      gk_csr_Write(model, model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
      sprintf(gu_string, "%s_v%d", ctrl->gu_file, iter+1);
      sprintf(gu_readable_string, "%s_readable_v%d", ctrl->gu_file, iter+1);
      gk_dwritefilebin(gu_string, train->nrows, g);
      gk_dwritefile(g, gu_readable_string, train->nrows);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    iter++;
  }
  
  /* Write */
    if (ctrl->id == 0) {
    sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file, iter - 1);
    sprintf(gu_string, "%s_v%d", ctrl->gu_file, iter);
    sprintf(gu_readable_string, "%s_readable_v%d", ctrl->gu_file, iter);
    gk_csr_Write(model, model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
    gk_dwritefilebin(gu_string, train->nrows, g);
    gk_dwritefile(g, gu_readable_string, train->nrows);
    }

  /* Freeing */
   gk_csr_Free(&model);
  gk_free((void **) &g, &participation, LTERM);
  if(ctrl->id==0){
    gk_free((void **) &norms, LTERM);
  }
  gk_csr_Free(&train);
  gk_csr_Free(&test);
  free_ctrl(ctrl);
  
  MPI_Finalize();
}
