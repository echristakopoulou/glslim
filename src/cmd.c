/**************************************************************/
/*! \file cmd.c
  
    \brief This file contains all the routines for parameter
           setup from the user
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*!
  \brief A structure for command-line options
 */
/**************************************************************/
static struct gk_option slim_options[] = {
  {"train_file", 1, 0, CMD_TRAIN_FILE},
  {"test_file", 1, 0, CMD_TEST_FILE},
  {"model_file", 1, 0, CMD_MODEL_FILE},
  {"prev_model_file", 1, 0, CMD_PREV_MODEL_FILE},
  {"participation_file", 1, 0, CMD_PARTICIPATION_FILE},
  {"gu_file", 1, 0, CMD_GU_FILE},
  {"stats_file", 1, 0, CMD_STATS_FILE},
  {"hr_file", 1, 0, CMD_HR_FILE},
  {"dbglvl", 1, 0, CMD_DBGLVL},
  {"lambda", 1, 0, CMD_LAMBDA},
  {"beta", 1, 0, CMD_BETA},
  {"local_beta", 1, 0, CMD_LOCAL_BETA},
  {"local_lambda", 1, 0, CMD_LOCAL_LAMBDA},
  {"start_iteri", 1, 0, CMD_START_ITERI},
  {"end_iteri", 1, 0, CMD_END_ITERI},
  {"starti", 1, 0, CMD_STARTI},
  {"endi", 1, 0, CMD_ENDI},
  {"optTol", 1, 0, CMD_OPTTOL},
  {"threshold", 1, 0, CMD_THRESHOLD},
  {"max_bcls_niters", 1, 0, CMD_MAX_BCLS_NITERS},
  {"bu", 1, 0, CMD_BU},
  {"bl", 1, 0, CMD_BL},
  {"size", 1, 0, CMD_SIZE},
  {"topn", 1, 0, CMD_TOPN},
  {"num_clusters", 1, 0, CMD_NUM_CLUSTERS},
  {"num_threads", 1, 0, CMD_NUM_THREADS},
  {"num_procs", 1, 0, CMD_NUM_PROCS},
  {"id", 1, 0, CMD_ID},
  {"help", 0, 0, CMD_HELP},
  {0, 0, 0, 0}
};



/**************************************************************/
/*! \brief Mini help
 */
/**************************************************************/
static char helpstr[][512] = {
  " ",
  " Usage",
  " 	executable [options]",
  " ",
  " 	 -train_file=string",
  " 		Specifies the input file which contains the training data. ",
  "             This file should be in .csr format. ",
  " 		",
  "      -test_file=string",
  "             Specifies the input file which contains the testing data. ",
  "             This file should be in .csr format.",
  " 	",
  "      -model_file=string",
  "             Specifies the output file which will contains a model matrix. ",
  "             The output file will be in .csr format. ",
  " ",
  "      -prev_model_file=string",
  "             Specifies the model file of the previous regularization ",
  "             which will be used as input in order to speed up learning. ",
  " ",
  "      -participation_file=string",
  "             Specifies the file which contrains the clustering assignment.",
  "             The output file will be an array of integers. ",
  " ",
  "      -gu_file=string",
  "             Specifies the file which contrains the user weights.",
  "             The output file will be an array of doubles. ",
  " ",
  "      -stats_file=string",
  "             Specifies the statistics file which will contain the training",
  "             error, the frobenius norm, the l1 norm and the objective value",
  " ",
  "      -hr_file=string",
  "             Specifies the output file which will contain the hr and arhr ",
  "             on the test set.",
  " ",
  "      -lambda=float",
  "             Specifies the regularization parameter for the $\ell_1$ norm",
  "             for the global component. ",
  "             The default value is 0.0.",
  " 		",
  "      -beta=float",
  "             Specifies the regularizationi parameter for the $\ell_2$ norm",
  "             for the global component. ",
  "             The default value is 0.0.",
  " ",
  "      -local_lambda=float",
  "             Specifies the regularization parameter for the $\ell_1$ norm",
  "             for the local component. ",
  "             The default value is 0.0.",
  "             ",
  "      -local_beta=float",
  "             Specifies the regularization parameter for the $\ell_2$ norm",
  "             for the local component. ",
  "             The default value is 0.0.",
  "             ",
  "      -starti=int",
  "             Specifies the index of the first column (C-style indexing) ",
  "             from which the sparse coefficient matrix will be calculated. ",
  "             The default value is 0.",
  " ",
  "      -endi=int",
  "             Specifies the index of the last column (exclusively) up to ",
  "             which the sparse coefficient matrix will be calculated. The ",
  "             default value is the number of total columns. ",
  " ",
  "      -dbglvl=int",
  "             Specifies the debug level. The default value is 0.",
  " ",
  "      -optTol=float",
  "             Specifies the threshold which control the optimization. Once ",
  "             the error from two optimization iterations is smaller than ",
  "             this value, the optimization process will be terminated. ",
  "             The default value is 1e-5.",
  " ",
  "      -max_bcls_niters=int",
  "             Specifies the maximum number of iterations that is allowed ",
  "             for optimization. Once the number of iterations reaches this ",
  "             value, the optimization process will be terminated. ",
  "             The default value is 1e4.",
  " ",
  "      -threshold= float",
  "             Specifies the ratio of the objective value between two subsequent ",
  "             iterations. Based on this value, it will specify when ",
  "             GLSLIMr0 will converge. The default value is 0.0001.",
  "             ",
  "      -bu = float",
  "             Specifies the upper bound for BCLS.",
  "             The default value is 1e20.",
  "             ",
  "      -bl = float",
  "             Specifies the lower bound for BCLS.",
  "             The default value is 0.",
  "             ",
  "      -topn=int",
  "             Specifies the number of recommendations to be produced for ",
  "             each user. The default value is 10.",
  " ",
  "      -start_iteri=int",
  "             Specifies the starting iteration number. ",
  "             The default value is 0.",
  " ",
  "      -end_iteri=int",
  "             Specifies the ending iteration number. ",
  "             The default value is 60.",
  " ",
  "      -size=int",
  "             Specifies the number of columns of the train file",
  "             in other words this is the column space of recommendations.",
  "             The default value is 1682.",
  " ",
  "      -num_clusters=int",
  "             Specifies the number of clusters. ",
  "             The default value is 5.",
  " ",
  "      -num_threads=int",
  "             Specifies the number of threads. ",
  "             The default value is 1.",
  " ",
  "      -num_procs=int",
  "             Specifies the number of processors. ",
  "             The default value is 1.",
  "             ",
  "      -id = int",
  "             Specifies the id of the processor. ",
  "             The default value is 0.",
  "             ",
  "      -help",
  "             Print this message.",
  " 			",
  ""
};


/**************************************************************/
/*! \brief A short help 
*/
/**************************************************************/
static char shorthelpstr[][100] = {
  " ",
  "   Usage: executable [options] ",
  "          use 'executable -help' for a summary of the options.",
  "          e.g. './glslim -help'",
  ""
};


/**************************************************************/
/*! \brief Entry point of the command-line argument parsing
  
    \param[out] ctrl  A ctrl structure to be filled out
    \param[in]  argc  Number of arguments
    \param[in]  argv  A list of arguments
*/
/**************************************************************/
void parse_cmdline(ctrl_t * ctrl, int argc, char *argv[])
{

  int c = -1, option_index = -1;

  if (ctrl == NULL)
    ctrl = create_ctrl();

  while ((c =
	  gk_getopt_long_only(argc, argv, "", slim_options,
			      &option_index)) != -1) {
    switch (c) {

    case CMD_TRAIN_FILE:
      ctrl->train_file = gk_strdup(gk_optarg);
      break;

    case CMD_TEST_FILE:
      ctrl->test_file = gk_strdup(gk_optarg);
      break;

    case CMD_PARTICIPATION_FILE:
      ctrl->participation_file = gk_strdup(gk_optarg);
      break;

    case CMD_GU_FILE:
      ctrl->gu_file = gk_strdup(gk_optarg);
      break;

    case CMD_MODEL_FILE:
      ctrl->model_file = gk_strdup(gk_optarg);
      break;

    case CMD_PREV_MODEL_FILE:
      ctrl->prev_model_file = gk_strdup(gk_optarg);
      break;

    case CMD_STATS_FILE:
      ctrl->stats_file = gk_strdup(gk_optarg);
      break;

    case CMD_HR_FILE:
      ctrl->hr_file = gk_strdup(gk_optarg);
      break;

    case CMD_DBGLVL:
      ctrl->dbglvl = atoi(gk_optarg);
      break;

    case CMD_LAMBDA:
      ctrl->lambda = atof(gk_optarg);
      break;

    case CMD_BETA:
      ctrl->beta = atof(gk_optarg);
      break;

    case CMD_LOCAL_BETA:
      ctrl->local_beta = atof(gk_optarg);
      break;

    case CMD_LOCAL_LAMBDA:
      ctrl->local_lambda = atof(gk_optarg);
      break;

    case CMD_STARTI:
      ctrl->starti = atoi(gk_optarg);
      break;

    case CMD_ENDI:
      ctrl->endi = atoi(gk_optarg);
      break;

    case CMD_START_ITERI:
      ctrl->start_iteri = atoi(gk_optarg);
      break;

    case CMD_END_ITERI:
      ctrl->end_iteri = atoi(gk_optarg);
      break;

    case CMD_OPTTOL:
      ctrl->optTol = atof(gk_optarg);
      break;
      
    case CMD_THRESHOLD:
      ctrl->threshold = atof(gk_optarg);
      break;

    case CMD_MAX_BCLS_NITERS:
      ctrl->max_bcls_niters = atoi(gk_optarg);
      break;
      
    case CMD_BU:
      ctrl->bu = atof(gk_optarg);
      break;
      
    case CMD_BL:
      ctrl->bl = atof(gk_optarg);
      break;
      
    case CMD_SIZE:
      ctrl->size = atoi(gk_optarg);
      break;

    case CMD_TOPN:
      ctrl->topn = atoi(gk_optarg);
      break;

    case CMD_NUM_CLUSTERS:
      ctrl->num_clusters = atoi(gk_optarg);
      break;

    case CMD_NUM_THREADS:
      ctrl->num_threads = atoi(gk_optarg);
      break;

    case CMD_NUM_PROCS:
      ctrl->num_procs = atoi(gk_optarg);
      break;

    case CMD_ID:
      ctrl->id = atoi(gk_optarg);
      break;

    case CMD_HELP:
      for (int i = 0; strlen(helpstr[i]) > 0; i++)
	printf("%s\n", helpstr[i]);
      exit(0);

    case '?':
    default:
      printf("Illegal command-line option(s) %s\n", gk_optarg);
      exit(0);

    }
  }

  if (argc - gk_optind != 0 || argc == 1) {
    for (int i = 0; strlen(shorthelpstr[i]) > 0; i++)
      printf("%s\n", shorthelpstr[i]);
    exit(0);
  }


}
