/**************************************************************/
/*! \file util.c
  
    \brief This file contains all the utility routines. 
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*! 
  \brief Create a ctrl structure wich contains all the default
         parameters.
 
  \return ctrl_t* A pointer to a created ctrl structure
*/
/**************************************************************/
ctrl_t *create_ctrl()
{

  ctrl_t *ctrl = gk_malloc(sizeof(ctrl_t), "malloc ctrl");

  ctrl->train_file = NULL;
  ctrl->test_file = NULL;
  ctrl->participation_file = NULL;

  ctrl->gu_file = NULL;

  ctrl->model_file = NULL;

  ctrl->prev_model_file = NULL;

  ctrl->hr_file = NULL;

  ctrl->stats_file = NULL;

  ctrl->dbglvl = 0;

  ctrl->beta = 0.0;
  ctrl->lambda = 0.0;

  ctrl->local_beta = 0.0;
  ctrl->local_lambda = 0.0;

  ctrl->starti = -1;
  ctrl->endi = -1;

  ctrl->start_iteri = 0;
  ctrl->end_iteri = 60;

  ctrl->num_clusters = 5;

  ctrl->optTol = 1e-5;
  
  ctrl->threshold = 0.0001;
  
  ctrl->max_bcls_niters = 10000;

  ctrl->bl = 0;
  ctrl->bu = 1e20;

  ctrl->size = 1682;

  ctrl->topn = 10;

  ctrl->num_threads = 1;

  ctrl->num_procs = 1;

  ctrl->id = 0;

  return ctrl;

}


/**************************************************************/
/*! 
  \brief Free a ctrl structure
  
  \param[in] ctrl A pointer to a ctrl structure to be freed
*/
/**************************************************************/
void free_ctrl(ctrl_t * ctrl)
{

  gk_free((void **) &ctrl->model_file, LTERM);
  gk_free((void **) &ctrl->prev_model_file, LTERM);
  gk_free((void **) &ctrl->train_file, LTERM);
  gk_free((void **) &ctrl->test_file, LTERM);
  gk_free((void **) &ctrl->gu_file, LTERM);
  gk_free((void **) &ctrl->participation_file, LTERM);
  gk_free((void **) &ctrl->hr_file, LTERM);
  gk_free((void **) &ctrl->stats_file, LTERM);

  gk_free((void **) &ctrl, LTERM);

}


/**************************************************************/
/*! 
  \brief Get a column from a csr matrix
  
  \param[in]  constraint A matrix from which one column is to 
                         be retrievd
  \param[in]  i          The index of the column to be retrieved
  \param[out] w          The output vector which saves the 
                         retrieved column
*/
/**************************************************************/
void get_column(gk_csr_t * constraint, int i, double *w)
{
  int nnz, j, k;

  if (i >= constraint->ncols) {
    gk_dset(constraint->nrows, 0, w);
  }
  else {
    nnz = constraint->colptr[i + 1] - constraint->colptr[i];

    for (j = 0; j < nnz; j++) {
      k = *(constraint->colptr[i] + j + constraint->colind);
      w[k] = *(constraint->colptr[i] + j + constraint->colval);
    }
  }
}


/**************************************************************/
/*!
  \brief Get a row from a csr matrix

  \param[in]  constraint A matrix from which one row is to
                         be retrieved
  \param[in]   i         The index of the row to be retrieved
  \param[out]  w         The output vector which saves the
                         retrieved row
*/
/****************************************************************/
void get_row(gk_csr_t * constraint, int i, double *w)
{
  int nnz, j, k;

  if (i >= constraint->nrows) {
    gk_dset(constraint->nrows, 0, w);
  }
  else {
    nnz = constraint->rowptr[i + 1] - constraint->rowptr[i];

    for (j = 0; j < nnz; j++) {
      k = *(constraint->rowptr[i] + j + constraint->rowind);
      w[k] = *(constraint->rowptr[i] + j + constraint->rowval);
    }
  }
}


/*****************************************************************************/
/*!

\brief Write an array of integers into a file, in human readable format.

\param[in] array      The array of integers we want to write.
\param[in] filename   The name of the file we write this array to.
\param[in] size       The size of the array we write.

 */
/*****************************************************************************/
void gk_i32writefile(int *array, char *filename, int size)
{

  FILE *fp;
  int u;

  fp = fopen(filename, "w");
  for (u = 0; u < size; u++) {
    fprintf(fp, "%d\n", array[u]);
  }
  fclose(fp);

}

/*****************************************************************************/
/*!

\brief Write an array of doubles into a file, in human readable format.

\param[in] array      The array of doubles we want to write.
\param[in] filename   The name of the file we write this array to.
\param[in] size       The size of the array we write.

*/
/*****************************************************************************/

void gk_dwritefile(double *array, char *filename, int size)
{

  FILE *fp;
  int u;

  fp = fopen(filename, "w");
  for (u = 0; u < size; u++) {
    fprintf(fp, "%f\n", array[u]);
  }
  fclose(fp);

}

/**************************************************************************/
