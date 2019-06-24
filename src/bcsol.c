/**************************************************************/
/*! \file bcsol.c
  
    \brief This file contains all the routines needed for 
           BCLS optimization. 
*/
/**************************************************************/



#include<slim.h>

/*! \brief Number of iterations */
int bcls_niters = 0;

// Type of projected search.
static int proj_search = BCLS_PROJ_SEARCH_EXACT;
/* static int proj_search = BCLS_PROJ_SEARCH_APPROX; */
// Method for Newton step computation.                           
static int newton_step = BCLS_NEWTON_STEP_LSQR;



/**************************************************************/
/*! \brief call_back function, periodically called by BCLS to 
   test if the user wants to exit. This is from BCLS.
*/
/**************************************************************/
int call_back(BCLS * ls, void *UsrWrk)
{
  int err;
  err = 0;			/* No error. */
  return err;
}

/**************************************************************/
/*! \brief call_back function, immediately terminate BCLS 
           iterations based on how many iterations it runs
*/
/**************************************************************/
int call_back_it(BCLS * ls, void *UsrWrk)
{
  wspace_t *wspace = (wspace_t *)UsrWrk;
  wspace->niters++;

  if (wspace->niters == ((wspace_t *) UsrWrk)->max_bcls_niters){
    printf("max niters reached\n");
    return 1;
  }
  else
    return 0;
}



/**************************************************************/
/*! \brief Pretty_printer, this is the print-routine that will 
           be used by BCLS for its output. This is from BCLS.
*/
/**************************************************************/
int pretty_printer(void *io_file, char *msg)
{
  fprintf(io_file, "%s", msg);
  return 0;
}


/**************************************************************/
/*! \brief Aprod, matrix-vector products. This is from BCLS.

  \details
  If     mode == BCLS_PROD_A  (0), compute y <- A *x, with x untouched;
  and if mode == BCLS_PROD_At (1), compute x <- A'*y, with y untouched.

*/
/**************************************************************/
int Aprod(const int mode, const int m, const int n, const int nix,
	  int ix[], double x[], double y[], void *UsrWrk)
{

  int j, k, l, ncols, acol;
  double xj, sum;
  wspace_t *Wrk;
  int *colind;
  ssize_t *colptr;
  float *colval;

  Wrk = (wspace_t *) UsrWrk;
  ncols = Wrk->ncols;
  colind = Wrk->mat->colind;
  colptr = Wrk->mat->colptr;
  colval = Wrk->mat->colval;
  acol = Wrk->acol;

  if (mode == BCLS_PROD_A) {

    memset(y, 0, sizeof(double) * m);

    for (l = 0; l < nix; l++) {
      j = ix[l];

      /* skip the inactive column */
      if (j % ncols == acol)
	continue;

      if ((xj = x[j]) != 0.0) {
	for (k = colptr[j]; k < colptr[j + 1]; k++) {
	  y[colind[k]] += xj * colval[k];
	}
      }
    }
  }

  else if (mode == BCLS_PROD_At) {
    for (l = 0; l < nix; l++) {
      j = ix[l];

      /* skip the inactive column */
      if (j % ncols == acol) {
	x[j] = 0;
	continue;
      }

      for (sum = 0.0, k = colptr[j]; k < colptr[j + 1]; k++) {
	sum += colval[k] * y[colind[k]];
      }
      x[j] = sum;
    }
  }


  return 0;
}


/**************************************************************/
/*! \brief BCLS learning. This is from BCLS

    \details This is to solve the problem
    \f[
    \begin{array}{ll} 
    \displaystyle\mathop{\hbox{minimize}}_x & \frac12\|Ax-a_i\|_2^2 + 
    \frac{1}{2}\beta\|x\|_2^2 + \lambda \|x\|_1 \\ 
    \hbox{subject to} & 0 \le x \\
                      & x_i = 0 \\
    \end{array} 
    \f]
 */
/**************************************************************/
void bcsol(ctrl_t * ctrl, double *bb, double *x, wspace_t * Wrk,
	   double *bl, double *bu, double beta, double *c, BCLS * ls)
{

  Wrk->niters = 0;
  /*  Problem dimensions. */
  int m = Wrk->mat->nrows;
  int n = Wrk->mat->ncols;

  /* init a bcls problem */

  bcls_init_prob(ls);
  bcls_set_problem_data(ls, m, n, Aprod, Wrk, beta, x, bb, c, bl, bu);

  /* set up tolerance */
  ls->optTol = ctrl->optTol;

  bcls_set_print_hook(ls, stdout, pretty_printer);
  ls->proj_search = proj_search;
  ls->newton_step = newton_step;

  if (ctrl->max_bcls_niters > 0)
    ls->CallBack = call_back_it;
  else
    ls->CallBack = call_back;

  /* call the solver */
  bcls_solve_prob(ls);
  /* solution */
  if (ctrl->dbglvl > 1) {
    int nnzx = 0;
    printf("\n Solution\n --------\n");
    printf("%4s  %18s %1s %18s %1s %18s  %18s\n",
	   "Var", "Lower", "", "Value", "", "Upper", "Gradient");
    for (int j = 0; j < n; j++) {
      if (x[j] > 1e-10) {
	nnzx++;
	char *blActiv = "";
	char *buActiv = "";
	if (x[j] - bl[j] < ls->epsx)
	  blActiv = "=";
	if (bu[j] - x[j] < ls->epsx)
	  buActiv = "=";
	printf("%4d  %18.11e %1s %18.11e %1s %18.11e  %18.11e\n",
	       j + 1, bl[j], blActiv, x[j], buActiv, bu[j], (ls->g)[j]);
      }
    }
    printf("%d nnz solution values\n", nnzx);
  }

}
