
##### The brief explanation of files in the dropbox #####

##### 1. Original format #####

 @@@ Functions_lapply_unit.R @@@

   : This R file contains individual functions which are used in the simulation code.
     There are 11 functions in this file and there are three groups according to the role of the functions.

     The first group is simulation-related functions and there are 3 functions.
     The second group is random data generating functions and there are 2 functions.
     The third group is test-statistics and their p-value deriving functions and there are 6 functions.

  --------------------------------------------------------------------------------------------------------------
  1) Simulation procedure

    * simulation2 : This function is the Monte-Carlo simulation code. Set inputs of the function and 
                   then according to the "mod.type"(size or power), it derives the results of size or power
                   for the significant level 0.05 and 0.1 respectively in the paper. That is, it calculates 
                   size or power for each significant level of RBB sup-wald, ADF tests by using simu.null and simu.alt functions below.
   
    * simu.null : This function is the Monte-Carlo simuation code for fixed r1(number of simulation), r2(number of bootstrap),
                  rho, theta and etc. It generates result of size test of table1 in the paper. This function is included in 
                  simulation function above. The difference with simu.alt function is it uses dgp_lin function below when generates random data.
                  There are two results according to the significant level 0.05 and 0.1 each. 

    * simu.alt : This function is the Monte-Carlo simuation code for fixed r1(number of simulation), r2(number of bootstrap),
                  rho, theta, addtinoally, gamma and alpha and etc. It generates result of power test of table1 in the paper.
                  This function is included in simulation function above. The difference with simu.null function is it uses 
                 dgp_thr function below when generates random data. There are two results according to the significant level 0.05 and 0.1 each. 


  2) Random data generating procedure for Monte-Carlo simulation

    * dgp_lin  : This function generates data x according to the random univariate standard normal data.
                 Data generating model is ARMA model. It is used in the simu.null function.

    * dgp_thr  : This function generates data x according to random univariate standard normal data.
                 Data generating model is band-TVECM model without lag-terms and with ARMA residual. 
                 It is used in the simu.alt function.

  3) Test-statistic and its p-value calculating procedure
    
    * block.boot : This function is rbb based sup-wald or adf test excuting function. There are 2 distinguished inputs
                   in this function. First, "mod.b" variable which has 2 options : "fix" or "opt". This variable determines
                   the type of deciding the number of block. If one chooses "fix", block size is determined from "bl" variable 
                   which is input of the function. On the other hand, if one chooses "opt", block size is determined not by "bl"
                   but by the function "b.star". This function is in the package "np" and it finds the optimal number of block.
                   Second, "mod.boot" variable which also has 2 options : "bnd" or "adf". This variable determines the type of test.
                   If one chooses "bnd", it executes rbb sup-wald test and "adf", it executes rbb adf test.
                   The final result is the same as rbb.bnd or rbb.adf function below.

    * rbb.bnd : This function is a rbb sup-wald test code. It calculates p-value of rbb sup-wald test.
                It uses bndadf function below.

    * rbb.adf : This function is a rbb adf test code. It calculates p-value of adf tset.
                This function uses adf function below.

    * bndadf  : It is a function which finds the sup-wald statistic. 
               It finds the argmax(for finding sup-wald) gamma1 and gamma2 by iterating according to each pre-set grid. 
               It uses find.w2 function on its execution.
               It also finds SETAR residual data which is used in the rbb.bnd function for generating y star.

    * find.w2 : This function is the code which calculates estimation of ECT terms(alpha)
                and wald test statistic for a given gamma1, gamma2 ,data and kn.

    * adf   : It is a function which derives ADF based residual and ADF test statistic for a given data.

 =================================================================================================================
 =================================================================================================================
 
 @@@ Simulation_unit.R @@@

   : This R file which executes Monte-Carlo simulation. In this file, one can set some main input parameters
     for the simulation function and then get the result of table1 of the paper which is explained above Functions_lapply_unit.R file.
     

 =================================================================================================================
 =================================================================================================================
 
 @@@ Parameter_unit.R @@@

   : This R file has 6 individual functions which are "block.boot.detail", "rbbbnd.detail","rbbadf.detail", "bndadf.det", "adf.det", "find.w2.det" codes.
     These functions are similar to the functions in "3) Test-statistic and its p-value calculating procedure". 
     The differeces are only on the result of them ( they give more detail and specified information )
     See r file if you want to know more details of the results. 
     

 =================================================================================================================
 =================================================================================================================
 










##### 2.Simple format #####

   : The difference is at the method of testing. The original format executes 200 number of simulations and there are 200 number of bootstrap procedure
     for each simulation. That is, the original method derives the empirical distribution of RBB sup-wald / RBB ADF test statistics
     by bootstrapping and calculates p-value of initial sup-wald and ADF statistics(non-bootstrapped statistics). And then, 
     by simulation, we can find the ratio of rejection.

     However, this approach takes very long time and so, we can use the alternative method. The simple method executes 200 number of 
     simulations but just 1 number of bootstrap procedure. For 200 number of simulations, we can derive 200 number of bootstrapped 
     RBB sup-wald and RBB ADF statistics from each of simulation. By using these statistics, we can find 0.9 and 0.95 quantile of statistics.
     And we regard these values as critical values and the test statistics are non-bootstrapped sup-wald and ADF test statistics for 
     each simulation. That is, we can find the alternative empirical distribution and critical values of sup-wald and ADF statistics 
     from the above 200 number of bootstarpped statistics, and then can calculate the ratio of rejection(we can judge whether the non-bootstrapped 
     initial test statistics is larger or smaller than this critical value and the rejection of null hypothesis).

     So the difference is in these 4 functions below : rbb.bnd, rbb,adf
     The other functions are the same as the functions in the original method.


  --------------------------------------------------------------------------------------------------------------

    * simu.null.b1 : This function is the same as "simu.null" function in the original format except some algoritm in the function.
                     It takes bootstrapped sup-wald / adf statistics and non-bootstrapped sup-wald / adf statistics from rbb.bnd.b1
                     and rbb.adf.b1 functions below. Then, it calculates the empirical critical values for significant level 0.05, 0.1 from
                     200 bootstrapped sup-wald and adf test statistics. Finally, finds the ratio of rejection. The only difference with 
                    "simu.alt.b1" is it uses dgp_lin function when generates random data.


    * simu.alt.b1 : This function is the same as "simu.alt" function in the original format except some algoritm in the function.
                    It takes bootstrapped sup-wald / adf statistics and non-bootstrapped sup-wald / adf statistics from rbb.bnd.b1
                    and rbb.adf.b1 functions below. Then, it calculates the empirical critical values for significant level 0.05, 0.1 from
                    200 bootstrapped sup-wald and adf test statistics. Finally, finds the ratio of rejection. The only difference with 
                    "simu.null.b1" is it uses dgp_thr function when generates random data.

    * rbb.bnd.b1 : This function is a rbb sup-wald test code. It derives non-bootstrapped initial sup-wald test statistics
                   and one bootstrapped(RBB) sup-wald statistics. This function uses bndadf function which is the same as in the original format.  

    * rbb.adf.b1 : This function is a rbb adf test code. It derives non-bootstrapped initial adf test statistics
                   and one bootstrapped(RBB) adf statistics. This function uses adf function which is the same as in the original format. 
 