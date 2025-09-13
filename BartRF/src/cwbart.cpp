/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "Rcpp.h"
#include "randomeffect.h"

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

//[[Rcpp::export]]
RcppExport SEXP cwbart(
    SEXP sexp_in,
    SEXP sexp_ip,
    SEXP sexp_inp,
    SEXP sexp_ix,
    SEXP sexp_iy,
    SEXP sexp_ixp,
    SEXP sexp_im,
    SEXP sexp_inc,
    SEXP sexp_ind,
    SEXP sexp_iburn,
    SEXP sexp_ipower,
    SEXP sexp_ibase,
    SEXP sexp_itau,
    SEXP sexp_inu,
    SEXP sexp_ilambda,
    SEXP sexp_isigest,
    SEXP sexp_iw,
    SEXP sexp_idart,
    SEXP sexp_itheta,
    SEXP sexp_iomega,
    SEXP sexp_igrp,
    SEXP sexp_ia,
    SEXP sexp_ib,
    SEXP sexp_irho,
    SEXP sexp_iaug,
    SEXP sexp_inkeeptrain,
    SEXP sexp_inkeeptest,
    SEXP sexp_inkeeptestme,
    SEXP sexp_inkeeptreedraws,
    SEXP sexp_inprintevery,
    SEXP sexp_Xinfo,
    SEXP sexp_Z,
    SEXP sexp_q,
    SEXP sexp_group_num,
    SEXP sexp_id,
    SEXP sexp_z_c0,
    SEXP sexp_z_d0,
    SEXP sexp_z_b,
    SEXP sexp_z_gamma,
    SEXP sexp_z_gamma_mean,
    SEXP sexp_z_gamma_cov,
    SEXP sexp_z_lambda,
    SEXP sexp_z_lambda_mean,
    SEXP sexp_z_lambda_cov,
    SEXP sexp_z_alpha,
    SEXP sexp_z_beta,
    SEXP sexp_seed
)
{
  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(sexp_in);
  size_t p = Rcpp::as<int>(sexp_ip);
  size_t np = Rcpp::as<int>(sexp_inp);
  Rcpp::NumericVector  xv(sexp_ix);
  double *ix = &xv[0];
  Rcpp::NumericVector  yv(sexp_iy); 
  double *iy = &yv[0];
  Rcpp::NumericVector  xpv(sexp_ixp);
  double *ixp = &xpv[0];
  size_t m = Rcpp::as<int>(sexp_im);
  Rcpp::IntegerVector _nc(sexp_inc);
  int *numcut = &_nc[0];
  size_t nd = Rcpp::as<int>(sexp_ind);
  size_t burn = Rcpp::as<int>(sexp_iburn);
  double mybeta = Rcpp::as<double>(sexp_ipower);
  double alpha = Rcpp::as<double>(sexp_ibase);
  double tau = Rcpp::as<double>(sexp_itau);
  double nu = Rcpp::as<double>(sexp_inu);
  double lambda = Rcpp::as<double>(sexp_ilambda);
  double sigma=Rcpp::as<double>(sexp_isigest);
  Rcpp::NumericVector  wv(sexp_iw); 
  double *iw = &wv[0];
  bool dart;
  if(Rcpp::as<int>(sexp_idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(sexp_ia);
  double b = Rcpp::as<double>(sexp_ib);
  double rho = Rcpp::as<double>(sexp_irho);
  bool aug;
  if(Rcpp::as<int>(sexp_iaug)==1) aug=true;
  else aug=false;
  double theta = Rcpp::as<double>(sexp_itheta);
  double omega = Rcpp::as<double>(sexp_iomega);
  Rcpp::IntegerVector _grp(sexp_igrp);
  int *grp = &_grp[0];
  size_t nkeeptrain = Rcpp::as<int>(sexp_inkeeptrain);
  size_t nkeeptest = Rcpp::as<int>(sexp_inkeeptest);
  size_t nkeeptestme = Rcpp::as<int>(sexp_inkeeptestme);
  size_t nkeeptreedraws = Rcpp::as<int>(sexp_inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(sexp_inprintevery);
  Rcpp::NumericMatrix Xinfo(sexp_Xinfo);
  
  //random effect
  int group_num = Rcpp::as<int>(sexp_group_num);
  int q = Rcpp::as<int>(sexp_q);
  
  Rcpp::IntegerVector idv(sexp_id);
  int *id = &idv[0];
  
  Rcpp::NumericVector Zv(sexp_Z);
  std::vector<std::vector<double>> Z_2d(n, std::vector<double>(q));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < q; j++) {
      Z_2d[i][j] = Zv[i * q + j];
    }
  }
  
  //max second dimension
  int max_second_dim = 0;
  for (int i = 0; i < group_num; i++) {
    if (id[i+1] - id[i] > max_second_dim) {
      max_second_dim = id[i+1] - id[i];
    }
  }
  std::vector<std::vector<std::vector<double>>> Z_3d(group_num, std::vector<std::vector<double>>(max_second_dim, std::vector<double>(q)));
  for (int i = 0; i < group_num; i++) {
    for (int j = 0; j < id[i+1] - id[i]; j++) {
      for (int l = 0; l < q; l++) {
        Z_3d[i][j][l] = Z_2d[j+id[i]][l];
      }
    }
  }
  
  double z_c0 = Rcpp::as<double>(sexp_z_c0);
  double z_d0 = Rcpp::as<double>(sexp_z_d0);
  
  Rcpp::NumericVector z_bv(sexp_z_b);
  std::vector<std::vector<double>> z_b(group_num, std::vector<double>(q));
  for (int i = 0; i < group_num; i++) {
    for (int j = 0; j < q; j++) {
      z_b[i][j] = z_bv[i * q + j];
      //z_b[i][j] = 1;
    }
  }
  
  Rcpp::NumericVector z_gamma_v(sexp_z_gamma);
  std::vector<double> z_gamma(q * (q - 1) / 2);
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    z_gamma[i] = z_gamma_v[i];
    //z_gamma[i] = 1;
  }
  std::vector<std::vector<double>> z_Gamma(q, std::vector<double>(q));
  z_Gamma = low_matrix(z_gamma, q);
  
  Rcpp::NumericVector z_gamma_meanv(sexp_z_gamma_mean);
  double *z_gamma_mean = &z_gamma_meanv[0];
  
  Rcpp::NumericMatrix z_gamma_covv(sexp_z_gamma_cov);
  double **z_gamma_cov = new double*[q * (q - 1) / 2];
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    z_gamma_cov[i] = new double[q * (q - 1) / 2];
    for (int j = 0; j < q * (q - 1) / 2; j++) {
      z_gamma_cov[i][j] = z_gamma_covv[i * q * (q - 1) / 2 + j];
    }
  }
  
  Rcpp::NumericVector z_lambda_v(sexp_z_lambda);
  std::vector<double> z_lambda(q);
  for (int i = 0; i < q; i ++) {
    z_lambda[i] = z_lambda_v[i];
    //z_lambda[i] = 1;
  }
  
  Rcpp::NumericVector z_lambda_meanv(sexp_z_lambda_mean);
  double *z_lambda_mean = &z_lambda_meanv[0];
  
  Rcpp::NumericVector z_lambda_covv(sexp_z_lambda_cov);
  double *z_lambda_cov = &z_lambda_covv[0];
  
  Rcpp::NumericVector z_alpha_v(sexp_z_alpha);
  double *z_alpha = &z_alpha_v[0];
  Rcpp::NumericVector z_beta_v(sexp_z_beta);
  double *z_beta = &z_beta_v[0];
  
  double *z_p_posterior = new double[q];
  for (int i = 0; i < q; i++) {
    if(z_lambda[i] == 0)
      z_p_posterior[i] = 0;
    else
      z_p_posterior[i] = 1;
  }
  
  double *z_indicator = new double[q];
  for (int i = 0; i < q; i++) {
    if(z_lambda[i] == 0)
      z_indicator[i] = 0;
    else
      z_indicator[i] = 1;
  }
  
  //sample from beta distribution beta(z_alpha, z_beta)
  std::vector<double> z_p(q);
  for (int i = 0; i < q; i++) {
    z_p[i] = R::rbeta(z_alpha[i] + z_indicator[i], z_beta[i] + 1 - z_indicator[i]);
  }
  
  double seed = Rcpp::as<double>(sexp_seed);
  
  // update z_Gamma
  for (int i = 0; i < q; i++) {
    if (z_lambda[i] == 0) {
      for (int j = 0; j < i; j++) {
        z_Gamma[i][j] = 0;
      }
      for (int j = i; j < q; j++) {
        z_Gamma[j][i] = 0;
      }
      z_Gamma[i][i] = 1;
    }
  }
  int gamma_index = 0;
  for (int i = 1; i < q; i++) {
    for (int j = 0; j < i; j++) {
      z_gamma[gamma_index] = z_Gamma[i][j];
      gamma_index ++;
    }
  }

  // Rcpp::NumericVector true_fix_v(sexp_true_fix);
  // std::vector<double> true_fix(n);
  // for (int i = 0; i < n; i++) {
  //   true_fix[i] = true_fix_v[i];
  // }

  // printf("n: %d\n", n);
  // printf("p: %d\n", p);
  // 
  // //print id
  // printf("id: \n");
  // for (int i = 0; i < group_num + 1; i++) {
  //   printf("%d ", id[i]);
  // }
  // printf("\n");
  // 
  // printf("y: \n");
  // for (int i = 0; i < n; i++) {
  //   printf("%lf ", iy[i]);
  // }
  // printf("\n");
  // 
  // printf("x: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < p; j++) {
  //     printf("%lf ", ix[i * p + j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("Z_2d: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < q; j++) {
  //     printf("%lf ", Z_2d[i][j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("Z_3d: \n");
  // for (int i = 0; i < group_num; i++) {
  //   for (int j = 0; j < id[i+1] - id[i]; j++) {
  //     for (int l = 0; l < q; l++) {
  //       printf("%lf ", Z_3d[i][j][l]);
  //     }
  //     printf("\n");
  //   }
  //   printf("\n");
  // }
  // 
  // printf("q: %d\n", q);
  // printf("group_num: %d\n", group_num);
  // printf("id: \n");
  // for (int i = 0; i < group_num + 1; i++) {
  //   printf("%d ", id[i]);
  // }
  // printf("\n");
  // 
  // printf("z c0: %lf\n", z_c0);
  // printf("z d0: %lf\n", z_d0);
  // 
  // printf("z b: \n");
  // for (int i = 0; i < group_num; i++) {
  //   for (int j = 0; j < q; j++) {
  //     printf("%lf ", z_b[i][j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("z gamma: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   printf("%lf ", z_gamma[i]);
  // }
  // printf("\n");
  // 
  // printf("z Gamma: \n");
  // for (int i = 0; i < q; i++) {
  //   for (int j = 0; j < q; j++) {
  //     printf("%lf ", z_Gamma[i][j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("z gamma mean: \n");
  // for (int i = 0; i < q; i++) {
  //   printf("%lf ", z_gamma_mean[i]);
  // }
  // printf("\n");
  // 
  // printf("z gamma cov: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%lf ", z_gamma_cov[i][j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("z lambda: \n");
  // for (int i = 0; i < q; i ++) {
  //   printf("%lf ", z_lambda[i]);
  // }
  // printf("\n");
  // //
  // printf("z lambda mean: \n");
  // for (int i = 0; i < q; i++) {
  //   printf("%lf ", z_lambda_mean[i]);
  // }
  // printf("\n");
  // 
  // printf("z lambda cov: \n");
  // for (int i = 0; i < q; i++) {
  //   printf("%lf ", z_lambda_cov[i]);
  // }
  // printf("\n");
  // 
  // printf("z alpha: %lf\n", z_alpha);
  // printf("z beta: %lf\n", z_beta);

   //return data structures (using Rcpp)
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector randmean(n); //random effect
   Rcpp::NumericVector fixmean(n); //fixed effect
   Rcpp::NumericVector temean(np);
   Rcpp::NumericVector sdraw(nd);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix fixdraw(nkeeptrain,n);
   Rcpp::NumericMatrix randdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::List mr_vecs(nkeeptreedraws);

   //random number generation
   arn gen;

   bart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<(size_t)numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }


   for(size_t i=0;i<n;i++) fixmean[i]=0.0;
   for(size_t i=0;i<n;i++) randmean[i]=0.0;
   for(size_t i=0;i<n;i++) trmean[i]=0.0;
   for(size_t i=0;i<np;i++) temean[i]=0.0;

   printf("*****Into main of wbart\n");
   //-----------------------------------------------------------

   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:beta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
                   mybeta,alpha,tau,nu,lambda);
   printf("*****sigma: %lf\n",sigma);
   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
   printf("*****printevery: %zu\n",printevery);
   printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);

   //--------------------------------------------------
   //heterbart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,iy,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);

   //--------------------------------------------------
   //sigma
   //gen.set_df(n+nu);
   double *svec = new double[n];
   for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;

   //--------------------------------------------------

   std::stringstream treess;  //string stream to write trees to  
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; //posterior mean for prediction
   if(np) { fhattest = new double[np]; }

   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   size_t temecnt=0; //count test draws into posterior mean
   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,keeptestme,keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();
   size_t total=nd+burn;
   
   // Rcpp::NumericMatrix b_save(nd, group_num * q);
   Rcpp::NumericMatrix gamma_save(nd, q * (q - 1) / 2);
   Rcpp::NumericMatrix lambda_save(nd, q);
   // Rcpp::NumericMatrix p_save(nd, q);
   
   double *bart_predictions = new double[n];
   for(size_t i=0;i<n;i++) bart_predictions[i] = bm.f(i);
   //center the bart predictions
   // double bart_mean = 0;
   // for (size_t i = 0; i < n; i++) {
   //   bart_mean += bart_predictions[i];
   // }
   // bart_mean = bart_mean / n;
   // for (size_t i = 0; i < n; i++) {
   //   bart_predictions[i] = bart_predictions[i] - bart_mean;
   //   //bart_predictions[i] = 0;
   // }

   for(size_t i=0;i<total;i++) {
     
      // printf("iteration %zu\n",i);
      if(i%printevery==0) printf("done %zu (out of %zu)\n",i,total);
      if(i==(burn/2)&&dart) bm.startdart();
      
      // printf("draw sigma\n");
      //random effect result
      std::vector<double> rf_predictions = rf_predict(n, group_num, id, Z_3d, z_lambda, z_Gamma, z_b, q);
      // printf("rf prediction: \n");
      // for (int j = 0; j < n; j++) {
      //   printf("%lf ", rf_predictions[j]);
      // }
      
      //draw bart (run once every 10 iterations)
      // printf("draw bart\n");
      double *bart_predictions_temp = new double[n];
      double bart_mean = 0;
      for (size_t j = 0; j < n; j++) {
        bart_predictions_temp[j] = 0;
      }
      if(i%1==0) {
        for (size_t k = 0; k < 1; k++) {
          bm.draw(sigma,rf_predictions,gen);
          for (size_t j = 0; j < n; j++) {
            bart_predictions_temp[j] += bm.f(j);
          }
        }
        for (size_t j = 0; j < n; j++) {
          bart_predictions[j] = bart_predictions_temp[j] / 1;
          //bart_predictions[j] = 0;
        }
  
        //center the bart predictions
        for (size_t j = 0; j < n; j++) {
          bart_mean += bart_predictions[j];
        }
        bart_mean = bart_mean / n;
        for (size_t j = 0; j < n; j++) {
          bart_predictions[j] = bart_predictions[j] - bart_mean;
        }
      }
      delete [] bart_predictions_temp;
      
      // for (size_t j = 0; j < n; j++) {
      //   bart_predictions[j] = 0;
      // }
      
      // 
      // printf("rf prediction: \n");
      // for (int j = 0; j < n; j++) {
      //   printf("%lf ", rf_predictions[j]);
      // }
      // printf("\n");
      //draw b
      // printf("draw b\n");
      z_b = sample_b(n, group_num, id, iy, bart_predictions, Z_3d, z_Gamma, z_lambda, pow(sigma,2), q, seed);
      
      
      //for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k)); rss += restemp*restemp;}
      //sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      // printf("sigma: %lf\n", sigma);
      
      // printf("z gamma: \n");
      // for (int i = 0; i < q * (q - 1) / 2; i++) {
      //   printf("%lf ", z_gamma[i]);
      // }
      // printf("\n");
      
      //draw lambda (lambda should behind gamma for the update of indicator)
      // printf("draw lambda\n");
      std::vector<std::vector<std::vector<double>>> z_T = calculate_T(group_num, id, iy, bart_predictions, Z_3d, z_Gamma, z_lambda, z_b, pow(sigma,2), q);
      z_lambda = sample_lambda(z_lambda_mean, z_lambda_cov, z_indicator, z_p, z_p_posterior, group_num, id, iy, bart_predictions, Z_3d, z_T, z_Gamma, z_lambda, z_b, pow(sigma,2), q, seed);
      
      if (i % 100 == 0) {
        printf("z lambda: \n");
        for (int i = 0; i < q; i ++) {
          printf("%lf ", z_lambda[i]);
        }
        printf("\n");
        printf("bart mean: ");
        printf("%lf\n", bart_mean);
      }
      
      //draw gamma
      // printf("draw gamma\n");
      // printf("calculate U\n");
      std::vector<std::vector<std::vector<double>>> z_U = calculate_U(group_num, id, iy, bart_predictions, Z_3d, z_Gamma, z_lambda, z_b, pow(sigma,2), q);
      z_gamma = sample_gamma(z_gamma_mean, z_gamma_cov, group_num, id, iy, bart_predictions, Z_3d, z_lambda, z_indicator, z_b, pow(sigma,2), z_U, q, seed);
      z_Gamma = low_matrix(z_gamma, q);
      
      //sigma
      sigma = sample_sigma(lambda,nu, n, group_num, id, iy, bart_predictions, rf_predictions,q,seed);
      
      //draw p
      z_p = sample_p(z_indicator, z_alpha, z_beta, q, seed);
      
      if(i>=burn) {
        
        sdraw[i-burn]=sigma;
        //save b, gamma, lambda
        // for (int j = 0; j < group_num; j++) {
        //   for (int k = 0; k < q; k++) {
        //     b_save(i-burn, j * q + k) = z_b[j][k];
        //   }
        // }
        int gamma_index = 0;
        for (int j = 1; j < q; j++) {
          for (int k = 0; k < j; k++) {
            gamma_save(i-burn, gamma_index) = z_Gamma[j][k];
            gamma_index ++;
          }
        }
        for (int j = 0; j < q; j++) {
          lambda_save(i-burn, j) = z_lambda[j];
          // p_save(i-burn, j) = z_p_posterior[j];
        }
        
        
         for(size_t k=0;k<n;k++) trmean[k] += (bart_predictions[k] + rf_predictions[k]);
         for(size_t k=0;k<n;k++) fixmean[k] += bart_predictions[k];
         for(size_t k=0;k<n;k++) randmean[k] += rf_predictions[k];
    
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k) + rf_predictions[k];
            for(size_t k=0;k<n;k++) fixdraw(trcnt,k)=bm.f(k);
            for(size_t k=0;k<n;k++) randdraw(trcnt,k)=rf_predictions[k];
            
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
            tecnt+=1;
         }
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
//	   #ifndef NoRcpp
//	   Rcpp::List lists(m*treesaslists);
//	   #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
/*      
	      #ifndef NoRcpp
	      varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
	      #endif
*/
	    }
      #ifndef NoRcpp
//	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    mr_vecs[k] = bm.getmrvec();
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      //varcnt(i-burn,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	      //varprb(i-burn,j)=ivarprb[j];
	    }
      #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif

            treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   for(size_t k=0;k<n;k++) fixmean[k]/=nd;
   for(size_t k=0;k<n;k++) randmean[k]/=nd;
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
   printf("check counts\n");
   printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   if(svec) delete [] svec;
   if(bart_predictions) delete [] bart_predictions;
   if(z_indicator) delete [] z_indicator;
   for (int i = 0; i < q * (q - 1) / 2; i++) {
     delete [] z_gamma_cov[i];
   }
   delete [] z_gamma_cov;
   delete [] z_p_posterior;
   
   

   //--------------------------------------------------
   //return
#ifndef NoRcpp
   Rcpp::List ret;
   ret["sigma"]=sdraw;
   ret["fix.train"]=fixdraw;
   ret["fix.train.mean"]=fixmean;
   ret["random.train"]=randdraw;
   ret["random.train.mean"]=randmean;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;
   ret["mr_vecs"]=mr_vecs;
   // ret["rf_b"] = b_save;
   ret["rf_gamma"] = gamma_save;
   ret["rf_lambda"] = lambda_save;
   // ret["rf_p"] = p_save;

   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
//   if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}
