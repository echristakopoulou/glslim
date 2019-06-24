/******************************************************************************/
/*! \file gu_learn.c

\brief This file contains the functions for 
       finding the weight of every user (g[u])
       used by GLSLIM and GLSLIMr0.
*//****************************************************************************/

#include <slim.h>

/*****************************************************************************/
/*! \brief A function which learns the weights of the users.
  
  \param[in] ctrl            Control structure
  \param[in] train           Training data
  \param[in] model           Model
  \param[in] participation   The cluster assignment
  
  \return     g              The vector of weights
*/
/*****************************************************************************/
double *learn_gu_all(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		     int *participation)
{
  int u, datasize, step, startu, endu;
  double *new_g;
  int i;
  int *sizes = NULL;
  int *displays = NULL;
  double *g = NULL;

  /* Set up MPI */
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

  new_g = gk_dmalloc(step, "malloc new g");

  /* Learn the weight gu for the set of users corresponding to this node. */
  for (u = startu; u < endu; u++) {
    new_g[u - startu] = learn_gu(ctrl, train, model, u, participation[u]);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* Gather the gus of the subset of the users of the different nodes
     to one main gu vector of node 0, which then is broadcasted to all nodes. */
  if (ctrl->id == 0) {
    sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
  }

  MPI_Gather(&step, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (ctrl->id == 0) {
    displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
	displays[i + 1] = displays[i] + sizes[i];
    }
  }

  g = gk_dmalloc(train->nrows, "malloc g");

  MPI_Gatherv(new_g, step, MPI_DOUBLE, g,
	      sizes, displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(g, train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Freeing */
  gk_free((void **) &new_g, LTERM);

  if (ctrl->id == 0) {
    gk_free((void **) &sizes, &displays, LTERM);
  }

  return g;

}


/*****************************************************************************/
/*! \brief A function which learns the gu of a specific user u.

  \param[in] ctrl         The ctrl structure
  \param[in] train        The training data
  \param[in] model        The model
  \param[in] u            The user u
  \param[in] cluster_id   The cluster assignment of the user u

  \return g               The personalized weight of a user u.
*/
/****************************************************************************/
double learn_gu(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		int u, int cluster_id)
{

  double *b;
  int ncols = ctrl->size;
  int i;
  double nom = 0;
  double denom = 0;
  double result = 0;
  double *global = NULL;
  double *local = NULL;
  int flag = 0;

  /* Compute training vector */
  b = gk_dsmalloc(ctrl->size, 0, "malloc b");

  for (i = train->rowptr[u]; i < train->rowptr[u + 1]; i++) {
    if (train->rowind[i] < ncols) {
      b[train->rowind[i]] = train->rowval[i];
    }
  }

  /* Compute vector of predictions based on the global and the local 
     portions of the model. */

  global = gk_dsmalloc(ncols, 0, "malloc global");
  local = gk_dsmalloc(ncols, 0, "malloc global");

  clean_predict(ctrl, model, train, u, &global, &local, cluster_id);

  /* Compute gu */
  for (i = 0; i < ncols; i++) {
    nom += (global[i] - local[i]) *(b[i] - local[i]);
    denom += (global[i] - local[i]) *(global[i] - local[i]);
  }
  if(denom!=0){
    result = nom / denom;
  }
  else{
    flag = 1;
    result = 2;
  }
  
  /* Scale gu to be in the [0-1] range */
  if ((result > 1)&&(flag == 0))
    result = 1;
  if (result < 0)
    result = 0;

  /* Freeing */
  gk_free((void **) &b, (void **) &global, (void **) &local, LTERM);

  return result;
}
