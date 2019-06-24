/**************************************************************/
/*! \File learn.c
    \brief This is the file containing the fundamental functions 
    to learn the model, used by GLSLIM and GLSLIMr0.

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2016
*/
/**************************************************************/

#include<slim.h>

/*****************************************************************************/
/*! \brief This function creates a new training matrix A'' and then 
           calls the main function local_learn, for computing the model.
	   Original training matrix A is of size n * m (rows * cols)

	   An intermediate training matrix A' is created which has:
           nrows = rows of the original training matrix A
           ncols = (number_of_clusters + 1) * cols of the original matrix A
	   Every user has then exactly the nnzs of A
	   plus the same nnzs copied for the cluster in which he belongs.
	   For all the other clusters he has 0 entries.
	   
	   In order to be able to regularize properly, we add underneath A'
	   a diagonal matrix which is of size: ncols * ncols.
	   This diagonal matrix contains the regularization params across the
	   diagonal. For the submatrix m * m (cols * cols) it has the global 
	   regularization (lg) across the diagonal and for the rest, 
	   it has the local regularization (ll).
	   The final matrix is A''.

           Example: for 3 clusters:
                m     m     m      m                                 
	   n  [                         ]       [     ]           [  b  ] n
           m  [ lg                      ]       [  w  ]           [  0  ] m
           m  [       ll                ]   *   [     ]    =      [  0  ] m
           m  [             ll          ]       [     ]           [  0  ] m
           m  [                    ll   ]                         [  0  ] m
                                                w is of size
                           A''                   4m * 1.


    \param[in] ctrl           The ctrl structure.
    \param[in] train          The training data A.
    \param[in] participation  The assignment of users to clusters.
    \param[in] g              The vector with the weights of the users
    \param[in] prev_model     The previous model (used to expedite the learning)

    \return    model          The model learnt.

 */
/*****************************************************************************/
gk_csr_t *learn(ctrl_t * ctrl, gk_csr_t * train, int *participation,
		double *g, gk_csr_t * prev_model)
{

  int i, pos, j, offset, global_nnz;
  gk_csr_t *mat, *model = NULL;

  global_nnz = train->rowptr[train->nrows] - train->rowptr[0];
  mat = gk_csr_Create();
  mat->ncols = (ctrl->num_clusters + 1) * train->ncols;
  mat->nrows = train->nrows + mat->ncols;
  mat->rowptr = gk_zmalloc(mat->nrows + 1, "gk_csr_Dep: local_rowptr");
  mat->rowind =
      gk_imalloc(2 * global_nnz + mat->ncols, "gk_csr_Dep: local_rowind");
  mat->rowval =
      gk_fmalloc(2 * global_nnz + mat->ncols, "gk_csr_Dep: local_rowval");
  mat->rowptr[0] = 0;

  for (i = 0, pos = 0; i < train->nrows; i++) {
    /* Copying the original training matrix A */
    for (j = train->rowptr[i]; j < train->rowptr[i + 1]; j++) {
      mat->rowind[pos] = train->rowind[j];
      mat->rowval[pos] = train->rowval[j] * g[i];
      pos++;
    }
    /* Copying the nnzs of the user to the cluster to which he belongs */
    offset = (participation[i] + 1) * train->ncols;

    for (j = train->rowptr[i]; j < train->rowptr[i + 1]; j++) {
      mat->rowind[pos] = train->rowind[j] + offset;
      mat->rowval[pos] = train->rowval[j] * (1 - g[i]);
      pos++;
    }
    mat->rowptr[i + 1] = pos;
  }

  /* Adding the diagonal matrix underneath A' */
  /* Global regularization */
  for (i = train->nrows; i < train->nrows + train->ncols; i++, pos++) {
    mat->rowind[pos] = i - train->nrows;
    mat->rowval[pos] = ctrl->beta;
    mat->rowptr[i] = pos;
  }
  /* Local regularization */
  for (i = train->nrows + train->ncols; i < mat->nrows; i++, pos++) {
    mat->rowind[pos] = i - train->nrows;
    mat->rowval[pos] = ctrl->local_beta;
    mat->rowptr[i] = pos;
  }
  mat->rowptr[mat->nrows] = pos;

  gk_csr_CreateIndex(mat, GK_CSR_COL);

  /* Learning the model */

  model = local_learn(ctrl, mat, prev_model);

  gk_csr_Free(&mat);

  return model;
}

/**************************************************************/
/*! \brief Learning
  
    \details This routine contains the learning algorithm used by GLSLIM and 
             GLSLIMr0.
    
    \param[in] ctrl       A ctrl structure which contains all the parameters
    \param[in] train      The training data A''		    
    \param[in] prev_model The previous model- used to speed up 
                          learning. It is optional argument.
   
    \return model         The model returned.

 */
/**************************************************************/
gk_csr_t *local_learn(ctrl_t * ctrl, gk_csr_t * train,
		      gk_csr_t * prev_model)
{
  int i, nr, nc, ni, pos, j, ii;
  int global_nnz;
  int basestart, baseend, datasize, step, starti, endi;
  gk_csr_t *mat = NULL;
  double *bl, *bu, *c;
  int *iinds, *jinds;
  float *vals;
  int *nnzs = NULL;
  int *rinds1 = NULL;
  int *rinds = NULL;
  int *rjinds = NULL;
  float *rvals = NULL;
  int rank = 0;
  int max_nnzs = 0;
  double tmr;
  int *rrowcnt = NULL;

  /* set up timers */
  gk_clearwctimer(tmr);
  gk_startwctimer(tmr);

  /* constants used across all problems */
  nr = train->nrows;
  nc = train->ncols;
  ni = train->ncols / (ctrl->num_clusters + 1);

  /* mallocing */
  bl = gk_dsmalloc(nc, ctrl->bl, "malloc bl");	/*lower bound */
  bu = gk_dsmalloc(nc, ctrl->bu, "malloc bu");	/*upper bound */
  c = gk_dmalloc(nc, "malloc c");	/*linear vector */

  gk_dset(ni, ctrl->lambda, c);
  gk_dset(nc - ni, ctrl->local_lambda, c + ni);

  /*starting and ending columns */
  basestart = (ctrl->starti >= 0) ? ctrl->starti : 0;
  baseend = (ctrl->endi >= 0) ? ctrl->endi : ni;
  datasize = baseend - basestart;
  step = (datasize / ctrl->num_procs) +
      (ctrl->id < (datasize % ctrl->num_procs) ? 1 : 0);
  starti = ((datasize / ctrl->num_procs) * ctrl->id) +
      gk_min(ctrl->id, datasize % ctrl->num_procs);
  endi = starti + step;

  if ((endi < datasize) && (ctrl->id == ctrl->num_procs - 1)) {
    endi = datasize;
    step = datasize - starti;
  }

  pos = 0;

  iinds = gk_ismalloc(step * nc, 0, "malloc iinds");
  jinds = gk_ismalloc(step * nc, 0, "malloc jinds");
  vals = gk_fsmalloc(step * nc, 0, "malloc vals");

  /* go through all columns  */
#pragma omp parallel num_threads(ctrl->num_threads)
  {
    int mypos;
    double *w, *b;
    wspace_t *myws;
    BCLS *ls;

    myws = (wspace_t *) gk_malloc(sizeof(wspace_t), "myws");
    myws->mat = train;
    myws->ncols = ni;

    ls = bcls_create_prob(nr, nc);

    w = gk_dsmalloc(nc, 0, "malloc w");
    b = gk_dsmalloc(nr, 0, "malloc b");
#pragma omp for private(i, j) schedule(dynamic)
    for (i = starti; i < endi; i++) {

      // this column is totally empty 
      if (train->colptr[i + 1] - train->colptr[i] == 0) {
	continue;
      }

      /**********************************************************/
      /*                BCLS learning                           */
      /**********************************************************/
      /* get the i-th column from A */

      for (j = train->colptr[i]; j < train->colptr[i + 1] - 1; j++) {
	ii = train->colind[j];
	b[ii] = 1;
      }

      myws->max_bcls_niters = gk_min(ctrl->max_bcls_niters,
				     50 * (train->colptr[i + 1] -
					   train->colptr[i]));

      gk_dset(nc, 0, w);

      // disable  

      myws->acol = i;
      if (prev_model != NULL) {
	get_row(prev_model, i, w);
      }

      bcsol(ctrl, b, w, myws, bl, bu, 0, c, ls);

      for (j = train->colptr[i]; j < train->colptr[i + 1] - 1; j++) {
	ii = train->colind[j];
	b[ii] = 0;
      }

      /**********************************************************/
      /*              dump the data                             */
      /**********************************************************/

      /* compute the triplets */
#pragma omp critical
      {
	for (j = 0; j < nc; j++) {
	  if (w[j] > EPSILON) {
	    mypos = pos++;
	    iinds[mypos] = i;
	    jinds[mypos] = j;
	    vals[mypos] = w[j];
	  }
	}
      }
    }				// end of starti - endi 
    bcls_free_prob(ls);
    gk_free((void **) &myws, (void **) &b, (void **) &w, LTERM);
  }
  gk_stopwctimer(tmr);
  
  /*  if(ctrl->id == 0){
    printf("time passed is %f\n", gk_getwctimer(tmr));
    }*/
  /**********************************************************/
  /*         Combine all the mat of the different processes
     to the total_mat with MPI                  */
  /**********************************************************/

  if (ctrl->id == 0) {
    nnzs = gk_imalloc(ctrl->num_procs, "malloc nnzs");
    gk_iset(ctrl->num_procs, 0, nnzs);
  }

  MPI_Gather(&pos, 1, MPI_INT, nnzs, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (ctrl->id == 0) {
    global_nnz = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      global_nnz += nnzs[i];
    }
  }

  MPI_Bcast(&global_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Finding the max nnzs between all nodes. This covers the case
     when a node which is not 0 has a bigger iinds/jinds/vals matrix 
     than the one in node 0 */

  if (ctrl->id == 0) {

    max_nnzs = nnzs[0];

    for (rank = 1; rank < ctrl->num_procs; rank++)
      if (nnzs[rank] > max_nnzs)
	max_nnzs = nnzs[rank];
  }

  /* Each node creates its own row count. This gets sent to node 0, in order
     to create the total rowptr. */

  if (global_nnz / ctrl->num_procs >= ni) {
    rrowcnt = gk_ismalloc(ni, 0, "rrowcnt");
    for (i = 0; i < pos; i++) {
      rrowcnt[iinds[i]]++;
    }
  }

  /* Every node sends its own iinds, jinds and vals. */
  if (ctrl->id != 0) {
    if (global_nnz / ctrl->num_procs < ni) {
      MPI_Send(iinds, pos, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else {
      MPI_Send(rrowcnt, ni, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Send(iinds, pos, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(jinds, pos, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(vals, pos, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
  }

  if (ctrl->id == 0) {
    if (global_nnz / ctrl->num_procs < ni) {
      rinds1 = gk_icopy(nnzs[0], iinds, gk_imalloc(max_nnzs, "rinds1"));
    }
    rinds = gk_icopy(nnzs[0], iinds, gk_imalloc(max_nnzs, "rinds"));
    rjinds = gk_icopy(nnzs[0], jinds, gk_imalloc(max_nnzs, "rjinds"));
    rvals = gk_fcopy(nnzs[0], vals, gk_fmalloc(max_nnzs, "rvals"));
  }

  gk_free((void **) &iinds, &jinds, &vals, &bl, &bu, &c, LTERM);

  /* Allocate and populate matrix */
  mat = gk_csr_Create();
  mat->nrows = ni;
  mat->ncols = nc;
  mat->rowptr = gk_zsmalloc(ni + 1, 0, "rowptr");
  mat->rowind = gk_imalloc(global_nnz, "rowind");
  mat->rowval = gk_fmalloc(global_nnz, "rowval");

  if (ctrl->id == 0) {
    if (global_nnz / ctrl->num_procs < ni) {
      for (rank = 0; rank < ctrl->num_procs; rank++) {
	if (rank != 0) {
	  MPI_Recv(rinds1, nnzs[rank], MPI_INT, rank, 0,
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	for (i = 0; i < nnzs[rank]; i++) {
	  mat->rowptr[rinds1[i]]++;
	}
      }
    }
    else {
      for (rank = 0; rank < ctrl->num_procs; rank++) {
	if (rank != 0) {
	  MPI_Recv(rrowcnt, ni, MPI_INT, rank, 0, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	}
	for (i = 0; i < ni; i++) {
	  mat->rowptr[i] += rrowcnt[i];
	}
      }
    }
    MAKECSR(i, mat->nrows, mat->rowptr);

    for (rank = 0; rank < ctrl->num_procs; rank++) {
      if (rank != 0) {
	MPI_Recv(rinds, nnzs[rank], MPI_INT, rank, 0,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(rjinds, nnzs[rank], MPI_INT, rank, 0,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(rvals, nnzs[rank], MPI_FLOAT, rank, 0,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      for (i = 0; i < nnzs[rank]; i++) {
	mat->rowind[mat->rowptr[rinds[i]]] = rjinds[i];
	mat->rowval[mat->rowptr[rinds[i]]] = rvals[i];
	mat->rowptr[rinds[i]]++;
      }
    }

    SHIFTCSR(i, mat->nrows, mat->rowptr);

    gk_free((void **) &rinds, &rjinds, &rvals,
	    &nnzs, LTERM);

    if (global_nnz / ctrl->num_procs < ni) {
      gk_free((void **) &rinds1, LTERM);
    }
  }

  if (rrowcnt != NULL) {
    gk_free((void **) &rrowcnt, LTERM);
  }

  /* Broadcast the matrix to the other nodes */
  MPI_Bcast(mat->rowptr, ni + 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(mat->rowind, global_nnz, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(mat->rowval, global_nnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

  return mat;
}
