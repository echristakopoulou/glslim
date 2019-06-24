/**************************************************************/
/*! \File main_create_gu.c
    \brief This is the file containing the main entry for 
    creating the binary gu file. 

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2016
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*! \brief The main entry for creating the initial binary gu file.
  For every user, the initial gu = 0.5. 
  The gu file has as many rows as the number of users.
 */
/**************************************************************/
int main(int argc, char *argv[])
{
  srand(0);
  int my_id, my_num_procs;
  ctrl_t *ctrl;
  gk_csr_t *train = NULL;
  char gu_string[50] = "";
  double *g;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &my_num_procs);


  /* parse command line */
  ctrl = create_ctrl();
  parse_cmdline(ctrl, argc, argv);

  ctrl->id = my_id;
  ctrl->num_procs = my_num_procs;

  /* read train and test file */
  train = gk_csr_Read(ctrl->train_file, GK_CSR_FMT_CSR, 1, 1);

  /* read gu file */
  sprintf(gu_string, "%s_v0", ctrl->gu_file);
  //g = gk_dreadfilebin(gu_string, &nrows);

  g = gk_dsmalloc(train->nrows, 0.5, "mcalloc g");
  gk_dwritefilebin(gu_string, train->nrows, g);

  /* freeing */
  gk_free((void **) &g, LTERM);
  gk_csr_Free(&train);
  free_ctrl(ctrl);

  MPI_Finalize();
}
