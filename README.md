LSLIM: Global and Local Sparse LInear Methods for Top-N Recommendation
________________________________________________________________________
________________________________________________________________________

GLSLIM is a toolkit of MPI software for top-N recommendation which implements 
the methods: GLSLIM, GLSLIMr0, LSLIM presented in the paper: 
Local Item-Item Models for Top-N Recommendation.

Building & Installing
----------------------
The code has been tested on Linux (x86_64) and it is MPI-dependent.
In order to compile the code, gcc (https://gcc.gnu.org/), 
CMake (https://cmake.org/), gsl (https://www.gnu.org/software/gsl/)
and MPI (https://www.open-mpi.org/) are required.

The code can run on MacOS as well, if the compiler used is gcc.
Also, the Command Line Tools will be of use: xcode-select --install

Once GLSLIM is downloaded, it is uncompressed with:

     $ tar -xzf glslim.tar.gz

which will create the directory glslim. 

First, the the BCLS library needs to be compiled and then we copy the necessary 
files to the folder BCLSlib:

In order to install GLSLIM, please run the command:

       $ ./build.sh 

Then, the executables create_gu, glslim, glslimr0 and lslim are created in 
glslim/build/src.

How to Run GLSLIM Methods
-------------------------

LSLIM
------

	$ mpirun -np 24 ./lslim -train_file=train_1 -test_file=test_1 -size=150
	  -participation_file=participation -model_file=model -hr_file=pred 
	  -stats_file=stats [-num_clusters=10 -local_beta=1 -local_lambda=0.1 
	  -prev_model_file=other_model -start_iteri=0 -end_iteri=10 -topn=20]

* This is an example command on how to run LSLIM with 24 cores for the 
  train_file 'train_1' and the test_file 'test_1', that have 150 items/columns 
  (specified by size) with the initial assignment of users to clusters in 
  participation_file 'participation_v0'. The files 'train_1', 'test_1' and 
  'participation_v0' need to be provided by the user. 
 
* The new assignment of users in iterations 1, 2 ,.. will be output in files 
  'participation_v1', 'participation_v2',... Also, files 
  'indiff_participation_v1', 'indiff_participation_v2',.. will be created that 
  will show for the users who remained in the same cluster, 
  whether this happened because there is no difference in the training error, 
  or because there is no cluster for which the training error is smaller.

GLSLIMr0
--------

	$ mpirun -np 24 ./glslimr0 -train_file=train_1 -test_file=test_1 
	   -size=150 -participation_file=participation -gu_file=gu 
	   -model_file=model -hr_file=pred -stats_file=stats [-num_clusters=10 
	   -beta=10 -lambda=0.01 -local_beta=1 -local_lambda=0.1 
	   -prev_model_file=other_model -start_iteri=0 -end_iteri=10 -topn=20]

*  This is an example command on how to run GLSLIMr0 with 24 cores for the 
   train_file 'train_1' and the test_file 'test_1', with 150 items/columns 
   (specified by size) with the initial assignment of users to clusters in 
   participation_file 'participation_v0' and the initial personalized weights of
   the users in gu_file 'gu_v0'. The files 'train_1', 'test_1', 
   'participation_v0' and 'gu_v0' need to be provided by the user.

*  The updated personalized weights of iterations 1, 2, e.t.c. will be output
   in files 'gu_v1', 'gu_v2', e.t.c. (in binary format). The updated 
   personalized weights in human readable format are output in files
   'gu_readable_v1', 'gu_readable_v2',e.t.c. 

GLSLIM
-------

	$ mpirun -np 24 ./glslim -train_file=train_1 -test_file=test_1 -size=150
	   -participation_file=participation -gu_file=gu -model_file=model 
	   -hr_file=pred -stats_file=stats [-num_clusters=10 -beta=10 
	   -lambda=0.01 -local_beta=1 -local_lambda=0.1 
	   -prev_model_file=other_model -start_iteri=0 -end_iteri=10 -topn=20]

*  This is an example command on how to run GLSLIM with 24 cores for 'train_1' 
   and 'test_1', that have 150 items/columns with the initial assignment of users
   to clusters in participation_file 'participation_v0' and the initial 
   personalized weights of the users in gu_file 'gu_v0'. The files 'train_1', 
   'test_1', 'participation_v0' and 'gu_v0' are provided by the user.

*  The updated personalized weights of iterations 1, 2, e.t.c. will be output in
   files 'gu_v1', 'gu_v2', e.t.c. These files are in binary format. The updated 
   personalized weights in human readable format are output in files 
   'gu_readable_v1', 'gu_readable_v2',e.t.c.

*  The new assignment of users in iterations 1,2,.. will be output in files 
   'participation_v1', 'participation_v2',... Also, files 
   'indiff_participation_v1','indiff_participation_v2',.. will be created that 
   will show for the users who remained in the same cluster, 
   whether this happened because there is no difference in the training error, 
   or because there is no cluster for which the training error is smaller.

The parameters described above for all the methods (train_file, test_file, 
gu_file, participation_file) were the input parameters. Below we will
present the rest of the parameters. With the exception of the parameter 
prev_model_file they are output parameters.
The parameters presented in brackets are optional. The rest are required. 

*  The option model_file=model specifies that the model of the first iteration 
   will be output in binary format in file 'model_bin_v0'. The model of the 
   second iteration will be output in file 'model_bin_v1', and so on. 

*  The option hr_file=pred means that the evaluation measures (HR, ARHR, 
   Precision and Recall) of every iteration will be output in the file 'pred'. 

*  The option stats_file=stats specifies that the different statistics per 
   iteration will be output in file 'stats'. The statistics are: training error,
   frobenius norm of the the models, l1 norm of the  models, value of objective 
   function and the number of users who switched clusters. 

*  The parameter num_clusters indicates the number of clusters we would like to 
   have. The default value is 5. 

*  The parameter local_beta corresponds to the value of the l2 regularization 
   parameter for the local component. The default value is 0.

*  The parameter local_lambda corresponds to the value of the l1 regularization 
   parameter for the local component. The default value is 0.

*  The parameter global_beta corresponds to the value of the l2 regularization 
   parameter for the global component. The default value is 0. 
   It is used only in GLSLIMr0 and GLSLIM.

*  The parameter global_lambda corresponds to the value of the l1 regularization
   parameter for the global component. The default value is 0. 
   It is used only in GLSLIMr0 and GLSLIM.

*  The option prev_model_file=other_model allows us to use another model, stored
   in the file 'other_model_bin_v0'. This model should be trained 
   on the same training file (usually with different regularization parameters)
   and it will be used as a warm-start for the first iteration, 
   to speed up the learning of the model. 

*  The option start_iteri specifies from which iteration the algorithm will start.
   This parameter is used when for some reason, our algorithm stopped, 
   before achieving convergence. Then, for example if our algorithm stopped at 
   iteration 2, we can select start_iteri=2, and then the algorithm will read 
   the model file 'model_bin_v1' and go on from there. The default value is 0.

*  The option end_iteri specifies the maximum number of iterations before the 
   algorithm exits. The default value is 60.

*  The option topn shows the size of the top-N list. The default value is 10.

To see a complete list of the available parameters, along with a small 
description, please open the manpage for GLSLIM - GLSLIMr0 - LSLIM.

         $ ./lslim -help
	 $ ./glslimr0 -help
	 $ ./glslim -help

File Formats
-------------
The train, test and model files are stored in CSR
(https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29) format.

We present the file formats for a toy example of 5 users and 3 items, using 2 
clusters, with the algorithm converging in 3 iterations.

Also please look at the example files provided in the folder 'example_files'. 
The example given is using 10 clusters.

* Train File
  
  The train file has a separate row per user. Every row has #item_id #rating 
  #item_id #rating e.t.c. Since we are using implicit data, all our ratings are 1.
  The example train file below shows that the first user purchased items 2 and 3,
  the second user purchased item 3, the third user purchased item 1, the fourth 
  user puchased items 1 and 2, and finally the fifth user purchased items 1 and 3.
  
  2 1 3 1
  3 1
  1 1
  1 1 2 1
  1 1 3 1

* Test File
  
  The test file has a separate row per user. Every row has #item_id #rating, 
  specifying which items have been hidden per user.
  The ratings are '1' here, too. 
  The example test file shown below shows that for the first user, the first 
  item was put in the test set, for the second, the third and the fifth users,
  the second item was put in the test set, while for the fourth user the third 
  item was hidden.

  1 1
  2 1
  2 1
  3 1
  2 1

* Model File

  The model file, learned by GLSLIM methods, has a separate row per item.
  Every row has #item_id #value #item_id #value e.t.c., showing the similarities
  of the item depicted in the row with the various #item_ids.
  In the example model file shown below, item 1 has similarity of value 1.20 
  with item 1 and of value 0.10 with item 2. The second item has similarity 0.01
  with the first item, similarity 0.95 with the second item and similarity 0.75 
  with the third item. Finally, the third item has similarity 0.90 with the 
  third item.

  1 1.20 2 0.10 
  1 0.01 2 0.95 3 0.75
  3 0.90

* Participation File
  
  The participation file has a separate row per user.
  Every row contains the id of the cluster the user is assigned to: 
  {0,1,...,#num_clusters-1}. The number of clusters presented in this file should 
  be the same as the number of clusters specified by the parameter num_clusters.
  The example participation file is shown below:

  0
  0
  1
  1
  0

* Gu File
  
  The gu file has a separate row per user, showing the personalized weight of the 
  user, which takes values between 0 and 1. The example gu file is shown below:
  
  0.3
  0.7
  0
  1
  0

  In order to create the initial binary gu_file 'gu_v0', the executable 
  create_gu is provided (after compiling). The command
  
	$ ./create_gu -train_file=train_1 -gu_file=gu

  will assign to every user an initial personalized weight of 0.5.

* Hr File
  
  The hr file has a new row per iteration, showing the evaluation measures. 
  The example hr file is shown below:

  HR = 0.91429 ARHR = 0.78140 Precision = 0.09143 Recall = 0.91429
  HR = 0.91584 ARHR = 0.78707 Precision = 0.09158 Recall = 0.91584
  HR = 0.91719 ARHR = 0.78584 Precision = 0.09172 Recall = 0.91719

* Stats File
  
  The stats file has a new row per iteration. The example stats file 
  is shown below. The second column indicates the training error, the fourth 
  column indicates the frobenius norm of the model, the sixth column indicates 
  the l1 norm of the model, the eighth column indicates the value of the 
  objective function and the tenth column shows the number of users who switched
  clusters.
  
```
  error 331045.401105 frob 559.244558 l1 2002.114266 obj 333606.759929
  error 309886.794087 frob 559.244558 l1 2002.114266 obj 312448.152911 user_diff 17699
  error 294140.738849 frob 859.860888 l1 2535.343362 obj 297535.943099
  error 287945.493466 frob 859.860888 l1 2535.343362 obj 291340.697715 user_diff 7399
  error 281689.162127 frob 996.559348 l1 2809.266735 obj 285494.988210
  error 278022.718709 frob 996.559348 l1 2809.266735 obj 281828.544792 user_diff 5559
```

Acknowledgements
-----------------
We would like to thank Professor Friedlander (http://www.cs.ubc.ca/~mpf/) for 
providing the BCLS library.

We would also like to thank Professor George Karypis (http://glaros.dtc.umn.edu/)
for providing the GKlib library.

Finally, we would like to thank Professor Xia Ning (http://cs.iupui.edu/~xning/)
for providing the SLIM code (http://www-users.cs.umn.edu/~xning/slim/html/) 
presented in the paper "SLIM: Sparse Linear Methods for Top-N Recommender Systems".

License
-------
Please look at the file LICENSE. The BCLS and GKLIB libraries are published under their own licenses.

Citation
---------
If you use our code in your research, please cite us:

```
inproceedings{christakopoulou2016local,
title={Local item-item models for top-n recommendation},
author={Christakopoulou, Evangelia and Karypis, George},
booktitle={Proceedings of the 10th ACM Conference on Recommender Systems},
pages={67--74},
year={2016},
organization={ACM}
}
```

Contact Information
-------------------
GLSLIM was written by Evangelia Christakopoulou and George Karypis.
If you have any issues, please send us an email to: evangel@cs.umn.edu, 
karypis@cs.umn.edu. 
