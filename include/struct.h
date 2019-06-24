/**************************************************************/
/*! \file
  
    \brief This file contains all the necessary data structures
*/
/**************************************************************/


#ifndef __STRUCT_H__
#define __STRUCT_H__



/**************************************************************/
/*!   
    \brief A data structure for ctrl parameters
 */
/**************************************************************/
typedef struct {

  /*! a file name that contains the training data in csr format */
  char *train_file;
  /*! a file name that contains the testing data in csr format */
  char *test_file;
  /*! a file name into which the model in csr format will be output */
  char *model_file;
  /*! a file name which has the model with the previous regularization */
  /* that will be input */
  char *prev_model_file;
  /*! a file name containing the clustering assignment */
  char *participation_file;
  /* a file name containing the vector of weights */
  char *gu_file;
  /*! a file name into which the training error & the norms will be output */
  char *stats_file;
  /*! a file into which the hr and arhr will be output */
  char *hr_file;
  /*! debug level, default 0 */
  int dbglvl;
  /*! the regularization parameter for L-1 norm for the global component */
  double lambda;
  /*! the regularization parameter for L-2 norm for the global component */
  double beta;
  /*! the regularization parameter for L-1 norm for the local component */
  double local_lambda;
  /*! the regularization parameter for L-2 norm for the local component */
  double local_beta;
  /*!the starting column index from which the coefficient matrix is calculated*/
  int starti;
  /*! the ending column index from which the coefficient matrix is calculated */
  int endi;
  /* the iteration id at which we will start (in that way the different
     iterations dont have to be run consecutively - but we can start 
     running the program, pause it and continue from where we stopped) */
  int start_iteri;
  /* the iteration id at which we will stop */
  int end_iteri;
  /*! optimality tolerance */
  float optTol;
  /* threshold on the difference between the objective value between two subsequent 
     iterations. It is used in glslimr0. */
  float threshold;
  /*! max number of iterations allowed in BCLS solver */
  int max_bcls_niters;
  /*! lower bound for BCLS */
  double bl;
  /*! upper bound for BCLS */
  double bu;
  /* the number of columns of the original train
     It needs to be specified since we extend the train and the model, so 
     we need to keep track of the original number of columns */
  int size;
  /*! the number of recommendations to be recommended */
  int topn;
  /*! the number of clusters */
  int num_clusters;
  /*! the number of threads */
  int num_threads;
  /*! the number of processors */
  int num_procs;
  /*! the id of the processor */
  int id;

} ctrl_t;

/**************************************************************/
/*!                                                                    
  \brief A workspace structure used for BCLS. This is adopated          
         from BCLS.                                            
*/
/**************************************************************/
typedef struct {
  int ncols;
  int acol;
  int max_bcls_niters;
  int niters;
  gk_csr_t *mat;
} wspace_t;

#endif
