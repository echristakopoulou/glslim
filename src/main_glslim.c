/**************************************************************/
/*! \file main_glslim.c
    \brief This is the file containing the main entry for GLSLIM.

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2016
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*! \brief The main entry for learning in participation mode.
  Participation-mode is an iterative mode, in which the clustering
  assignment as well as the gu vector is updated in every
  iteration, along with the model learnt.
 */
/**************************************************************/
int main(int argc, char *argv[])
{

  ctrl_t *ctrl;
  gk_csr_t *train, *test;
  char model_bin_string[150] = "";
  //  char model_string[150]="";
  char participation_string[150] = "";
  char gu_string[150] = "";
  char gu_readable_string[150] = "";
  char indifference_string[150]="";
  int *participation;
  int * participation2;
  double *g;
  int * indifference;
  int i;
  int iter, my_id, my_num_procs,  participation_diff;
  size_t nrows;
  gk_csr_t *model = NULL;
  gk_csr_t *prev_model = NULL;
  double *norms = NULL;
  double error = 0;
  double objective = 0;
  FILE *fpout = NULL;
  double threshold = 0;
  
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
  nrows = (size_t) train->nrows;
  test = gk_csr_Read(ctrl->test_file, GK_CSR_FMT_CSR, 1, 1);

  iter = ctrl->start_iteri;

  /* Reading gu file */
  sprintf(gu_string, "%s_v%d", ctrl->gu_file, iter);
  g = gk_dreadfilebin(gu_string, &nrows);
  //  g = gk_dsmalloc(train->nrows, 0.5, "malloc g");

  indifference = gk_imalloc(nrows, "malloc indifference");
  /* Reading participation file */
  sprintf(participation_string, "%s_v%d", ctrl->participation_file, iter);
  participation = gk_i32readfile(participation_string, &nrows);

  participation_diff = train->nrows;

  /* If it is not the 1st iteration, read the previous model. If it is, 
     use the model learnt with different regularization as previous */

  if (iter != 0) {
    sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file, iter - 1);
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
  
    threshold = 0.01 * train->nrows;
  
    while ((iter < ctrl->end_iteri) && (participation_diff > threshold)) {
  
    /* Learn the model for this iteration and then save it as
       prev model for the next iteration */
      model = learn(ctrl, train, participation, g, prev_model);
    if (prev_model != NULL) {
      gk_csr_Free(&prev_model);
    }
    prev_model = model;
    
    MPI_Barrier(MPI_COMM_WORLD);
  
    /* Model column indexing */
      if ((ctrl->size) * (ctrl->num_clusters + 1) > model->ncols)
      model->ncols = (ctrl->size) * (ctrl->num_clusters + 1);
    gk_csr_CreateIndex(model, GK_CSR_COL);
  
    /* Compute the frobenius and l1 norm of the model */
      if (ctrl->id == 0) {
      norms = norm(ctrl, model);
    }
  
    /* Compute training error */
      error = train_predict(ctrl, train, model, g, participation);

    /* Compute the value of the objective function */
   if (ctrl->id == 0) {
      fpout = gk_fopen(ctrl->stats_file, "a", "stats open");
      objective = norms[0] + norms[1] + error;
      fprintf(fpout, "error %f frob %f l1 %f obj %f\n",
	      error, norms[0], norms[1], objective);
      gk_fclose(fpout);
      }

    /* Compute the HR and ARHR */
   slim_test(ctrl, model, g, train, test, participation);

   gk_free((void **) &g, LTERM);
   participation2 = participation;
   
   gk_iset(train->nrows, 0, indifference);
    
    /* Learn the clustering assignment for this iteration
       (The gu vector gets updated at the same time) */
   participation = learn_pu_all(ctrl, train, model, &g, participation2, 
				&indifference);

    /* Computes how many users changed clustering assignment */
   participation_diff = 0;
   for (i = 0; i < train->nrows; i++) {
     if (participation[i] != participation2[i]) {
       participation_diff++;
     }
   }
   gk_free((void **) &participation2, LTERM);
  
    /* Compute training error, objective and HR/ARHR 
       for the same model , with the new gu and participation. */
      error = train_predict(ctrl, train, model, g, participation);

    if (ctrl->id == 0) {
      fpout = gk_fopen(ctrl->stats_file, "a", "stats open");
      objective = norms[0] + norms[1] + error;
      fprintf(fpout, "error %f frob %f l1 %f obj %f user_diff %d\n",
	      error, norms[0], norms[1], objective, participation_diff);
      gk_fclose(fpout);
    }
    
    
    
    slim_test(ctrl, model, g, train, test, participation);
    
    if (ctrl->id == 0){
     //      sprintf(model_string, "%s_v%d", ctrl->model_file, iter);
      //      gk_csr_Write(model, model_string, GK_CSR_FMT_CSR, 1, 1);
      sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file, iter);
      gk_csr_Write(model, model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
      sprintf(gu_string, "%s_v%d", ctrl->gu_file, iter+1);
      sprintf(gu_readable_string, "%s_readable_v%d", ctrl->gu_file, iter+1);
      sprintf(participation_string, "%s_v%d", ctrl->participation_file,
	      iter+1);
      sprintf(indifference_string, "indiff_%s_v%d", ctrl->participation_file,
              iter+1);
      gk_dwritefilebin(gu_string, train->nrows, g);
      gk_dwritefile(g, gu_readable_string, train->nrows);
      gk_i32writefile(participation, participation_string, train->nrows);
      gk_i32writefile(indifference, indifference_string, train->nrows);
    }
    
    iter++;
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
    
    /* Write */
    if (ctrl->id == 0) {
      sprintf(model_bin_string, "%s_bin_v%d", ctrl->model_file, iter - 1);
      sprintf(gu_string, "%s_v%d", ctrl->gu_file, iter);
      sprintf(gu_readable_string, "%s_readable_v%d", ctrl->gu_file, iter);
      sprintf(participation_string, "%s_v%d", ctrl->participation_file,
	      iter);
      sprintf(indifference_string, "indiff_%s_v%d", ctrl->participation_file,
	      iter+1);
      gk_csr_Write(model, model_bin_string, GK_CSR_FMT_BINROW, 1, 1);
      gk_dwritefilebin(gu_string, train->nrows, g);
      gk_dwritefile(g, gu_readable_string, train->nrows);
      gk_i32writefile(participation, participation_string, train->nrows);
      gk_i32writefile(indifference, indifference_string, train->nrows);
    }
    
  /* Freeing */
    gk_csr_Free(&model);
  gk_csr_Free(&train);
  gk_csr_Free(&test);
  gk_free((void **) &g, &participation, &indifference, LTERM);
  if(ctrl->id==0){
    gk_free((void **)&norms, LTERM);
  }
  free_ctrl(ctrl);
  
  MPI_Finalize();
}
