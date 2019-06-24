/**************************************************************/
/*! \file
  
    \brief This file contains all prototypes.
*/
/**************************************************************/


#ifndef __PROTO_H__
#define __PROTO_H__

/* cmd.c */
void parse_cmdline(ctrl_t * ctrl, int argc, char *argv[]);

/* util.c */
ctrl_t *create_ctrl();
void free_ctrl(ctrl_t * ctrl);
void get_column(gk_csr_t * constraint, int i, double *w);
void get_row(gk_csr_t * constraint, int i, double *w);
void gk_i32writefile(int *array, char *filename, int size);
void gk_dwritefile(double *array, char *filename, int size);


/* bcsol.c */
int call_back(BCLS * ls, void *UsrWrk);
int call_back_it(BCLS * ls, void *UsrWrk);
int pretty_printer(void *io_file, char *msg);
int Aprod(int mode, int m, int n, int nix, int ix[],
	  double x[], double y[], void *UsrWrk);
void bcsol(ctrl_t * ctrl, double *bb, double *x, wspace_t * Wrk,
	   double *bl, double *bu, double beta, double *c, BCLS * ls);


/* predict.c */
void slim_test(ctrl_t * ctrl, gk_csr_t * model, double *g,
	       gk_csr_t * train, gk_csr_t * test, int *participation);
int recommend(ctrl_t * ctrl, gk_csr_t * model,
	      double g_participation, int *iidx,
	      gk_csr_t * train, int u, gk_dkv_t ** rcmd, int cluster);
double *all_predict(ctrl_t * ctrl, gk_csr_t * model, gk_csr_t * train,
		    int u, double g_participation, int cluster);
void clean_predict(ctrl_t * ctrl, gk_csr_t * model, gk_csr_t * train,
		   int u, double **global_rcmd, double **local_rcmd,
		   int cluster);
double training_error(ctrl_t * ctrl, gk_csr_t * model,
		      double g_participation, gk_csr_t * train, int u,
		      int cluster);
double train_predict(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		     double *g, int *participation);

/* lslim_predict.c*/ 
void lslim_test(ctrl_t * ctrl, gk_csr_t * model, 
               gk_csr_t * train, gk_csr_t * test, int *participation);
int lslim_recommend(ctrl_t * ctrl, gk_csr_t * model,
              int *iidx,
              gk_csr_t * train, int u, gk_dkv_t ** rcmd, int cluster);
double *lslim_all_predict(ctrl_t * ctrl, gk_csr_t * model, gk_csr_t * train,
                    int u, int cluster);
double lslim_training_error(ctrl_t * ctrl, gk_csr_t * model,
                      gk_csr_t * train, int u,
                      int cluster);
double lslim_train_predict(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
                     int *participation);


/* learn.c */
gk_csr_t *learn(ctrl_t * ctrl, gk_csr_t * train, int *participation,
		double *g, gk_csr_t * prev_model);
gk_csr_t *local_learn(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * prev_model);


/* lslim_learn.c */
gk_csr_t *lslim_learn(ctrl_t * ctrl, gk_csr_t * train, int *participation,
                 gk_csr_t * prev_model);
gk_csr_t *lslim_local_learn(ctrl_t * ctrl, gk_csr_t * orig_train, 
			    gk_csr_t * train, gk_csr_t * prev_model);


/* participation_learn.c */
int *learn_pu_all(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		  double **g, int *participation, int **indifference);
int learn_pu(ctrl_t * ctrl, gk_csr_t * train, int u, gk_csr_t * model, 
	     double *g_participation, int participation, int *indiff);

/* lslim_participation_learn.c */
int *lslim_learn_pu_all(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
			  int * participation, int **indifference);
int lslim_learn_pu(ctrl_t * ctrl, gk_csr_t * train, int u, gk_csr_t * model,
             int participation, int *indiff);

/* gu_learn.c */
double *learn_gu_all(ctrl_t * ctrl, gk_csr_t * train,
		     gk_csr_t * model, int *participation);
double learn_gu(ctrl_t * ctrl, gk_csr_t * train, gk_csr_t * model,
		int u, int cluster_id);

/* norm.c */
double *norm(ctrl_t * ctrl, gk_csr_t * model);

#endif
