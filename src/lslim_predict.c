/**************************************************l************/
/*! \file lslim_predict.c
  
    \brief This file contains all the routines for testing and
    computing training error, used by LSLIM.  

*/
/**************************************************************/


#include<slim.h>

/**************************************************************/
/*! \brief Top-N recommendations and evaluations
    
    \param[in] ctrl          A ctrl structure
    \param[in] model         A model
    \param[in] train         The training data from which the model is learned
    \param[in] test          The testing data
    \param[in] participation The vector of the clustering assignment.
 */
/**************************************************************/
void lslim_test(ctrl_t * ctrl, gk_csr_t * model,
	       gk_csr_t * train, gk_csr_t * test, int *participation)
{

  int nu, ni, nhits, n, u, jj, kk, i;
  int datasize, step, startu, endu, nrcmd;
  double hr = 0;
  double arhr = 0;
  double precision = 0;
  double recall = 0;
  double local_precision, local_recall;
  double *eval, *overall_eval;
  int *iidx = NULL;
  gk_dkv_t *rcmd = NULL;
  FILE *fpout = NULL;

  /* evaluation results for return */
  eval = gk_dsmalloc(4, 0, "malloc eval");

  /* set up MPI - deciding which subset of users each processor will predict */

  nu = test->nrows;
  datasize = nu;
  step = (datasize / ctrl->num_procs) +
      (ctrl->id < (datasize % ctrl->num_procs) ? 1 : 0);
  startu =
      ((datasize / ctrl->num_procs) * ctrl->id) + gk_min(ctrl->id,
							 datasize %
							 ctrl->num_procs);
  endu = startu + step;

  if ((endu < datasize) && (ctrl->id == ctrl->num_procs - 1)) {
    endu = datasize;
    step = datasize - startu;
  }

  /* Computing the recommendations for the subset of users decided before */

  n = 0;
  ni = ctrl->size;
  iidx = gk_imalloc(ni, "malloc *iidx");
  gk_iset(ni, -1, iidx);
  rcmd = gk_dkvmalloc(ni, "malloc rcmd");

  for (u = startu; u < endu; u++) {
    nhits = 0;

    /* no testing instances for this user */
    if (test->rowptr[u + 1] - test->rowptr[u] == 0) {
      continue;
    }
    n++;

    /* top-n recommendation */
    nrcmd = lslim_recommend(ctrl, model, iidx, train, u, &rcmd,
		      participation[u]);

    /* stats for the recommendation */
    local_precision = 0;
    local_recall = 0;
    /* evaluations */
    for (jj = 0; jj < nrcmd; jj++) {
      for (kk = test->rowptr[u]; kk < test->rowptr[u + 1]; kk++) {
	/* hit hit */
	if (rcmd[jj].val == test->rowind[kk]) {
	  nhits++;
	  eval[0] += 1.0;	// hit rate
	  eval[1] += 1.0 / (double) (jj + 1);	// arhr
	  break;
	}
      }
    }
    local_precision = nhits / (double) nrcmd;
    local_recall =
	nhits / (double) (test->rowptr[u + 1] - test->rowptr[u]);


    /* clean up rcmd for this user */
    eval[2] += local_precision;
    eval[3] += local_recall;
  }

  /* Combining stats of the different subsets of users to the main ones */
  if (ctrl->id == 0) {
    overall_eval = gk_dsmalloc((4 * ctrl->num_procs), 0, "malloc sizes");
  }

  MPI_Gather(eval, 4, MPI_DOUBLE, overall_eval, 4, MPI_DOUBLE, 0,
	     MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  if (ctrl->id == 0) {
    hr = 0;
    arhr = 0;
    precision = 0;
    recall = 0;

    for (i = 0; i < (ctrl->num_procs * 4); i += 4) {
      hr += overall_eval[i];
      arhr += overall_eval[i + 1];
      precision += overall_eval[i + 2];
      recall += overall_eval[i + 3];
    }

    hr /= (double) nu;
    arhr /= (double) nu;
    precision /= (double) nu;
    recall /= (double) nu;

    fpout = gk_fopen(ctrl->hr_file, "a", "gk_csr_Write: fpout");
    fprintf(fpout,
	    "HR = %.5f ARHR = %.5f Precision = %.5f Recall = %.5f\n",
	    hr, arhr, precision, recall);
    fclose(fpout);
  }

  /* finish up */
  gk_free((void **) &rcmd, &iidx, &eval, LTERM);

  if (ctrl->id == 0) {
    gk_free((void **) &overall_eval, LTERM);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}



/**************************************************************/
/*! \brief Top-N recommendation for a user
  
    \param[in] ctrl            A ctrl structure
    \param[in] model           A model
    \param[in] iidx            An auxiliary array for efficient recommendations
    \param[in] train           Training data from which the model is learned
    \param[in] u               The index of the user for which the top-n 
                               recommendations are generated
    \param[out] rcmd           The list of recommendations, in which the 
                               keys are the recommendation scores and the 
		               values are the item indices
   \param[in] cluster          The clustering assignment of the user.
			       
    \return int                The actual number of recommendations
 */
/**************************************************************/
int lslim_recommend(ctrl_t * ctrl, gk_csr_t * model,
	      int *iidx,
	      gk_csr_t * train, int u, gk_dkv_t ** rcmd, int cluster)
{
  int ni, ii, nuitrn, i, nrcmd, nrcmd2;
  double *global = NULL;

  ni = ctrl->size;

  /* Finding out which of the items the user has not rated/bought */
  for (ii = train->rowptr[u]; ii < train->rowptr[u + 1]; ii++) {
    if (train->rowind[ii] < ni)
      iidx[train->rowind[ii]] -= 1;
    else
      break;
  }

  /* When a user has no rated items, he has no recommnedations either */
  nuitrn = train->rowptr[u + 1] - train->rowptr[u];
  if (nuitrn == 0) {
    *rcmd = NULL;
    return 0;
  }


  /* computing the predictions for all items for user u */
  global = lslim_all_predict(ctrl, model, train, u, cluster);

  /* creating the recommendation list for the user 
     (for only the items he has not rated) */

  nrcmd = 0;

  for (i = 0; i < ni; i++) {
    if (iidx[i] >= -1) {
      (*rcmd)[nrcmd].key = global[i];
      (*rcmd)[nrcmd].val = i;
      nrcmd++;
    }
  }

  /* Populating again iidx with -1 to be ready for the next user */
  for (ii = train->rowptr[u]; ii < train->rowptr[u + 1]; ii++) {
    if (train->rowind[ii] < ni)
      iidx[train->rowind[ii]] += 1;
    else
      break;
  }

  /* sorting */
  gk_dkvsortd(nrcmd, *rcmd);
  nrcmd2 = gk_min(nrcmd, ctrl->topn);

  /* finishing up */
  gk_free((void **) &global, LTERM);

  return nrcmd2;

}


/**************************************************************/
/*! \brief Compute all possible predictions for user u for all items.

  \param[in] ctrl             Ctrl structure.
  \param[in] model            The model.
  \param[in] train            Training data.
  \param[in] u                The user in question.
  \param[in] cluster_id       The cluster assignment of user u.
  
  \return    rcmd2            The vector of size ni containing all the 
                              prediction values.
  
 */
/****************************************************************/
double *lslim_all_predict(ctrl_t * ctrl, gk_csr_t * model, gk_csr_t * train,
		    int u, int cluster_id)
{

  int ni, i, ii, j, jj, index;
  double *rcmd2;
  double fixed_rcmd_part;

  ni = ctrl->size;

  rcmd2 = gk_dsmalloc(ni, 0, "malloc rcmd2");

  for (i = train->rowptr[u]; i < train->rowptr[u + 1]; i++) {
    ii = train->rowind[i];
    fixed_rcmd_part = train->rowval[i];
    index = (ctrl->size) * (cluster_id) + ii;
    for (j = model->colptr[index]; j < model->colptr[index + 1]; j++) {
      jj = model->colind[j];
      rcmd2[jj] += fixed_rcmd_part * model->colval[j];
    }
  }


  return rcmd2;
}

/******************************************************************************/
/*! \brief Compute the training error for a specific user.

  \param[in] ctrl            Control structure.
  \param[in] model           The model.
  \param[in] train           The training data.
  \param[in] u               The user in question.
  \param[in] cluster         The clustering assignment of the user in question.
  
  \return    error           The training error for this user.
 */
/******************************************************************************/
double lslim_training_error(ctrl_t * ctrl, gk_csr_t * model,
		      gk_csr_t * train, int u,
		      int cluster)
{
  int ni, i;
  double *global = NULL;
  double *training_row;
  double error;

  ni = ctrl->size;

  training_row = gk_dsmalloc(train->ncols, 0, "malloc training_row");
  get_row(train, u, training_row);

  global = lslim_all_predict(ctrl, model, train, u, cluster);

  error = 0;

  for (i = 0; i < ni; i++) {
    error += (global[i] - training_row[i]) *(global[i] - training_row[i]);
  }


  /* Freeing up */

  gk_free((void **) &global, &training_row, LTERM);

  return error;
}


/******************************************************************************/
/*! \brief Compute the training error for all users.

  \param[in] ctrl           The control structure. 
  \param[in] train          The training data.
  \param[in] model          The model.
  \param[in] participation  The clustering assignment vector of all users.
  
  \return    error          The overall training error.
*/
/******************************************************************************/
double lslim_train_predict(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		     int *participation)
{
  int u, datasize, step, startu, endu, i; 
  double total_error = 0;
  double overall_error = 0;
  double error = 0;
  double *overall_errors = NULL;

  /* Setting up MPI */

  datasize = train->nrows;
  step = (datasize / ctrl->num_procs) +
      (ctrl->id < (datasize % ctrl->num_procs) ? 1 : 0);
  startu = ((datasize / ctrl->num_procs) * ctrl->id) +
      gk_min(ctrl->id, datasize % ctrl->num_procs);
  endu = startu + step;

  if ((endu < datasize) && (ctrl->id == ctrl->num_procs - 1)) {
    endu = datasize;
    step = datasize - startu;
  }

  /* Each processor computes the training error for a subset of the users. */

  for (u = startu; u < endu; u++) {
    error = lslim_training_error(ctrl, model, train, u, participation[u]);
    total_error += error;
  }

  /* Combining the training errors from the different processors to the main one
     and printing the error */
  if (ctrl->id == 0) {
    overall_errors = gk_dsmalloc(ctrl->num_procs, 0, "malloc sizes");
  }

  MPI_Gather(&total_error, 1, MPI_DOUBLE, overall_errors, 1, MPI_DOUBLE, 0,
	     MPI_COMM_WORLD);

  if (ctrl->id == 0) {
    overall_error = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      overall_error += overall_errors[i];
    }

    gk_free((void **) &overall_errors, LTERM);
  }

  return overall_error;
}
