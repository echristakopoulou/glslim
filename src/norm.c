/**************************************************l************/
/*! \file norm.c
    \brief This file computes the Frobenius and l1 norm.
*/
/**************************************************************/
#include<slim.h>

/**************************************************************/
/*! \brief The main entry for computing the norms.
  
  \param_in  ctrl   The control structure
  \param_in  model  The model 
 
  \return    norms  The matrix containing the Frobenius and l1 norm.
  
 */
/**************************************************************/
double *norm(ctrl_t * ctrl, gk_csr_t * model)
{
  int i;
  double frobenius_norm, l1_norm;
  double *norms;

  norms = gk_dmalloc(2, "malloc norms");
  gk_csr_ComputeSums(model, GK_CSR_COL);
  gk_csr_ComputeSquaredNorms(model, GK_CSR_COL);
  frobenius_norm = 0;
  l1_norm = 0;

  for (i = 0; i < ctrl->size; i++) {
    frobenius_norm += (ctrl->beta) * (model->cnorms[i]);
    l1_norm += (ctrl->lambda) * (model->csums[i]);
  }

  for (i = ctrl->size; i < model->ncols; i++) {
    frobenius_norm += (ctrl->local_beta) * (model->cnorms[i]);
    l1_norm += (ctrl->local_lambda) * (model->csums[i]);
  }

  norms[0] = frobenius_norm;
  norms[1] = l1_norm;

  return norms;
}
