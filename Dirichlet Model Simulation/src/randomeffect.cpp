#include <Rcpp.h>
#include <eigen3/Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
using namespace Rcpp;

/*** R
library(MASS)
library(invgamma)
library(truncnorm)
library(stats)
*/


// transfer a vector to a low matrix (checked 06.25)
std::vector<std::vector<double>> low_matrix(std::vector<double> gamma, int p) {
  std::vector<std::vector<double>> low(p, std::vector<double>(p));
  int temp = 0;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      if(i == j) {
        low[i][j] = 1;
      } 
      else if(i < j) { 
        low[i][j] = 0;
      }
      else {
        low[i][j] = gamma[temp + j];
      }
    }
    temp = temp + i;
  }
  return low;
}

// random effect predict (checked 06.25 checked by code)
std::vector<double> rf_predict(int N, int n, int *id, std::vector<std::vector<std::vector<double>>> z,
                               std::vector<double> lambda, std::vector<std::vector<double>> Gamma, std::vector<std::vector<double>> b, int q) {
  Eigen::MatrixXd Lambda = Eigen::MatrixXd::Zero(q, q);
  for (int i = 0; i < q; i++) {
    Lambda(i, i) = lambda[i];
  }
  Eigen::MatrixXd Gamma_matrix = Eigen::MatrixXd::Zero(q, q);
  for (int i = 0; i < q; i++) {
    for (int j = 0; j < q; j++) {
      Gamma_matrix(i, j) = Gamma[i][j];
    }
  }
  Eigen::MatrixXd b_matrix = Eigen::MatrixXd::Zero(n, q);
  for (int i = 0; i < n; i++) {
    b_matrix.row(i) = Eigen::Map<Eigen::VectorXd>(&b[i][0], q);
  }
  
  Eigen::MatrixXd beta_matrix = Eigen::MatrixXd::Zero(n, q);
  beta_matrix = (Lambda * Gamma_matrix * b_matrix.transpose()).transpose();
  
  std::vector<double> result(N);
  for(int i = 0; i < n; i++) {
    for(int j = id[i]; j < id[i+1]; j++) {
      Eigen::VectorXd z_vector = Eigen::VectorXd::Zero(q);
      z_vector = Eigen::Map<Eigen::VectorXd>(&z[i][j-id[i]][0], q);
      result[j] =  (beta_matrix.row(i)) * z_vector;
    }
  }
  return result;
}

// sample sigma (be careful for z) (checked 06.26)
double sample_sigma(double lambda, double nu, int N, int n, int *id, double *y, double *bart_predictions,
                    std::vector<double> rf_predictions, int q, int seed) {
  Eigen::VectorXd y_vector = Eigen::VectorXd::Zero(N);
  y_vector = Eigen::Map<Eigen::VectorXd>(y, N);

  Eigen::VectorXd bart_predictions_vector = Eigen::VectorXd::Zero(N);
  bart_predictions_vector = Eigen::Map<Eigen::VectorXd>(bart_predictions, N);
  
  Eigen::VectorXd rf_predictions_vector = Eigen::VectorXd::Zero(N);
  rf_predictions_vector = Eigen::Map<Eigen::VectorXd>(&rf_predictions[0], N);
  
  Eigen::VectorXd residuals = y_vector - bart_predictions_vector - rf_predictions_vector;
  
  // printf("c0: %f\n", c_0);
  // printf("d0: %f\n", d_0);
  
  // double c = c_0 + N / 2.0;
  // double d = d_0 + 0.5 * (residuals.transpose() * residuals).value();
  // printf("residuals.transpose() * residuals: %f\n", (residuals.transpose() * residuals).value());
  
  double sigma = sqrt((nu * lambda + (residuals.transpose() * residuals).value()) / R::rchisq(nu + N));
  
  // double c = c_0 + N / 2.0;
  // double d = d_0 + 0.5 * (residuals.transpose() * residuals).value();
  // 
  // //inverse gamma sampler, c is shape, d is scale
  // Rcpp::Environment mvn("package:invgamma");
  // Rcpp::Function rinvgamma = mvn["rinvgamma"];
  // Rcpp::NumericVector sigma_vector = rinvgamma(1, c, d);
  // double sigma = sqrt(sigma_vector[0]);
  
  return sigma;
}

//sample b (n * q) prior: N(0,I) (checked)
std::vector<std::vector<double>> sample_b(int N, int n, int *id, double *y, double *bart_predictions,
                                          std::vector<std::vector<std::vector<double>>> z, std::vector<std::vector<double>> Gamma, std::vector<double> lambda,
                                          double sigma, int q, int seed) {
  Eigen::VectorXd y_vector = Eigen::VectorXd::Zero(N);
  y_vector = Eigen::Map<Eigen::VectorXd>(y, N);
  
  Eigen::VectorXd bart_predictions_vector = Eigen::VectorXd::Zero(N);
  bart_predictions_vector = Eigen::Map<Eigen::VectorXd>(bart_predictions, N);
  
  Eigen::MatrixXd Gamma_matrix = Eigen::MatrixXd::Zero(q, q);
  for (int i = 0; i < q; i++) {
    Gamma_matrix.row(i) = Eigen::Map<Eigen::VectorXd>(&Gamma[i][0], q);
  }
  
  //diagonal matrix for lambda
  Eigen::MatrixXd Lambda = Eigen::MatrixXd::Zero(q, q);
  for (int i = 0; i < q; i++) {
    Lambda(i, i) = lambda[i];
  }
  
  Eigen::MatrixXd b_matrix = Eigen::MatrixXd::Zero(n, q);
  for (int i = 0; i < n; i++) {
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(id[i+1] - id[i], q);
    for (int j = 0; j < id[i+1] - id[i]; j++) {
      Z.row(j) = Eigen::Map<Eigen::VectorXd>(&z[i][j][0], q);
    }
    Eigen::MatrixXd H_i = Eigen::MatrixXd::Zero(q, q);
    Eigen::RowVectorXd h_i = Eigen::VectorXd::Zero(q);
    Eigen::MatrixXd v_i = Eigen::MatrixXd::Zero(q, q);
    
    for (int j = id[i]; j < id[i+1]; j++) {
      Eigen::RowVectorXd v_ij = Eigen::VectorXd::Zero(q);
      v_ij = (Z.row(j-id[i]) * Lambda * Gamma_matrix);
      v_i = v_i + v_ij.transpose() * v_ij;
      h_i = h_i + v_ij * (y[j] - bart_predictions[j]);
    }
    
    // printf("v_i is: \n");
    // for (int j = 0; j < q; j++) {
    //   for (int k = 0; k < q; k++) {
    //      printf("%f ", v_i(j, k));
    //   }
    //    printf("\n");
    // }
    // arma::mat v_i_arma(q, q);
    // for (int j = 0; j < q; j++) {
    //   for (int k = 0; k < q; k++) {
    //     v_i_arma(j, k) = v_i(j, k) / sigma;
    //   }
    // }
    // arma::mat H_i_arma(q, q);
    // H_i_arma = inv_sympd(v_i_arma + arma::eye(q, q));
    // 
    // for (int j = 0; j < q; j++) {
    //   for (int k = 0; k < q; k++) {
    //     H_i(j, k) = H_i_arma(j, k);
    //   }
    // }
    Eigen::MatrixXd H_i_copy = Eigen::MatrixXd::Zero(q, q);
    H_i_copy = (1 / sigma) * v_i + Eigen::MatrixXd::Identity(q, q);
    
    // printf("H_i_copy is: \n");
    // for (int j = 0; j < q; j++) {
    //   for (int k = 0; k < q; k++) {
    //     printf("%f ", H_i_copy(j, k));
    //   }
    //   printf("\n");
    // }
    
    H_i = ((1 / sigma) * v_i + Eigen::MatrixXd::Identity(q, q)).llt().solve(Eigen::MatrixXd::Identity(q,q));
    
    h_i = (1 / sigma) * h_i * H_i;
    
    // printf("H_i inverse is: \n");
    // for (int j = 0; j < q; j++) {
    //   for (int k = 0; k < q; k++) {
    //     printf("%f ", H_i(j, k));
    //   }
    //   printf("\n");
    // }
    // printf("h_i is: \n");
    // for (int j = 0; j < q; j++) {
    //   printf("%f ", h_i(j));
    // }
    // printf("\n");
    
    Rcpp::NumericVector h_i_vector(q);
    for (int j = 0; j < q; j++) {
      h_i_vector[j] = h_i(j);
    }
    Rcpp::NumericMatrix H_i_matrix(q, q);
    for (int j = 0; j < q; j++) {
      for (int k = 0; k < q; k++) {
        H_i_matrix(j, k) = H_i(j, k);
      }
    }
    
    //mvn sampler
    Rcpp::Environment mvn("package:MASS");
    Rcpp::Function mvrnorm = mvn["mvrnorm"];
    Rcpp::NumericVector b_vector = mvrnorm(1, h_i_vector, H_i_matrix);
    for (int j = 0; j < q; j++) {
      b_matrix(i, j) = b_vector[j];
    }
  }
  std::vector<std::vector<double>> b(n, std::vector<double>(q));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < q; j++) {
      b[i][j] = b_matrix(i, j);
    }
  }
  
  //printf
  // printf("b is: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < q; j++) {
  //     printf("%f ", b[i][j]);
  //   }
  //   printf("\n");
  // }
  
  return b;
}

//calculate U (a three-dimensional array) (checked 06.25)
std::vector<std::vector<std::vector<double>>> calculate_U(int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                      std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>>b, double sigma, int q) {
  //max second dimension of U:n*id[n]
  int max_id = 0;
  for (int i = 0; i < n; i++) {
    if (id[i+1] - id[i] > max_id) {
      max_id = id[i+1] - id[i];
    }
  }
  //U:n*id[n]*(q*(q-1)/2)
  std::vector<std::vector<std::vector<double>>> U(n, std::vector<std::vector<double>>(max_id, std::vector<double>(q * (q - 1) / 2)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < id[i+1] - id[i]; j++) {
      int temp = 0;
      for (int h = 1; h < q; h++) {
        for (int l = 0; l < h; l++){
          U[i][j][temp] = b[i][l] * lambda[h] * z[i][j][h];
          temp ++;
        }
      }
    }
  }
  
  // printf("U is: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < id[i+1] - id[i]; j++) {
  //     for (int l = 0; l < q * (q - 1) / 2; l++) {
  //       printf("%f ", U[i][j][l]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  
  return U;
}

//calculate T (a three-dimensional array) (checked 06.25)
std::vector<std::vector<std::vector<double>>> calculate_T(int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                      std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>> b, double sigma, int q) {
  int max_id = 0;
  for (int i = 0; i < n; i++) {
    if (id[i+1] - id[i] > max_id) {
      max_id = id[i+1] - id[i];
    }
  }
  std::vector<std::vector<std::vector<double>>> T(n, std::vector<std::vector<double>>(max_id, std::vector<double>(q)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < id[i+1] - id[i]; j++) {
      for (int l = 0; l < q; l++) {
        double med_val = b[i][l];
        for (int h = 0; h < l; h++){
          med_val = med_val + b[i][h] * Gamma[l][h];
        }
        T[i][j][l] = med_val * z[i][j][l];
      }
    }
  }
  return T;
}

//calculate Eta (a three-dimensional array) (checked 06.25)
std::vector<std::vector<std::vector<double>>> calculate_Eta(int n, int *id, double *y, double *bart_predictions,
                        std::vector<double> lambda, std::vector<std::vector<std::vector<double>>> T, int q) {
  int max_id = 0;
  for (int i = 0; i < n; i++) {
    if (id[i+1] - id[i] > max_id) {
      max_id = id[i+1] - id[i];
    }
  }
  
  // //print bart_predictions
  // printf("bart_predictions is: \n");
  // for (int i = 0; i < id[n]; i++) {
  //   printf("%f ", bart_predictions[i]);
  // }
  // printf("\n");
  // 
  // //print y
  // printf("y is: \n");
  // for (int i = 0; i < id[n]; i++) {
  //   printf("%f ", y[i]);
  // }
  // printf("\n");

  // //print T
  // printf("T is: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < id[i+1] - id[i]; j++) {
  //     for (int l = 0; l < q; l++) {
  //       printf("%f ", T[i][j][l]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");

  std::vector<std::vector<std::vector<double>>> Eta(n, std::vector<std::vector<double>>(max_id, std::vector<double>(q)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < id[i+1] - id[i]; j++) {
      for (int l = 0; l < q; l++) {
        Eta[i][j][l] = y[j + id[i]] - bart_predictions[j + id[i]];
        for (int h = 0; h < q; h++){
          if (h != l) {
            Eta[i][j][l] = Eta[i][j][l] - T[i][j][h] * lambda[h];
          }
        }
      }
    }
  }
  
  // printf("Eta is: \n");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < id[i+1] - id[i]; j++) {
  //     for (int l = 0; l < q; l++) {
  //       printf("%f ", Eta[i][j][l]);
  //     }
  //     printf("\n");
  //   }
  // }
  // printf("\n");
  
  return Eta;
}

//sample lambda (prior: N(mean,sigma2))
std::vector<double> sample_lambda(double *mean_prior,double *sigma_prior,double *indicator, std::vector<double> p, double *posterior_p,
                      int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                      std::vector<std::vector<std::vector<double>>> T,
                      std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>> b, double sigma, int q, int seed) {
  
  for(int l = 0; l < q; l++) {
    // printf("l is: %d\n", l);
    Eigen::MatrixXd T_matrix = Eigen::MatrixXd::Zero(id[n], q);
    for (int i = 0; i < n; i++) {
      for (int j = id[i]; j < id[i+1]; j++) {
        T_matrix.row(j) = Eigen::Map<Eigen::VectorXd>(&T[i][j-id[i]][0], q);
      }
    }
    
    std::vector<std::vector<std::vector<double>>> Eta = calculate_Eta(n, id, y, bart_predictions, lambda, T, q);
    
    Eigen::MatrixXd Eta_matrix = Eigen::MatrixXd::Zero(id[n], q);
    for (int i = 0; i < n; i++) {
      for (int j = id[i]; j < id[i+1]; j++) {
        Eta_matrix.row(j) = Eigen::Map<Eigen::VectorXd>(&Eta[i][j-id[i]][0], q);
      }
    }
    
    double med_val1 = T_matrix.col(l).transpose() * T_matrix.col(l);
    double med_val2 = T_matrix.col(l).transpose() * Eta_matrix.col(l);
    
    double inversew2 = med_val1 / sigma;
    
    double lambda_tilde = med_val2 / med_val1;
    
    double sigma2_hat = 1 / (inversew2 + 1 / sigma_prior[l]);
    double lambda_hat = sigma2_hat * (inversew2 * lambda_tilde + (1 / sigma_prior[l]) * mean_prior[l]);
    
    // printf("lambda_hat is: %f\n", lambda_hat);
    // printf("sigma2_hat is: %f\n", sigma2_hat);
    // printf("p[l] is: %f\n", p[l]);
    double p_w0 = 0;
    Rcpp::Environment base("package:stats");
    Rcpp::Function pnorm = base["pnorm"];
    Rcpp::NumericVector pnorm_result = pnorm(-lambda_hat / sqrt(sigma2_hat));
    double cdf1 = pnorm_result[0];
    Rcpp::NumericVector pnorm_result2 = pnorm(-mean_prior[l] / sqrt(sigma_prior[l]));
    double cdf2 = pnorm_result2[0];
    // printf("cdf is: %f\n", log(1 - cdf));
    double p_w1 = log(sqrt(sigma2_hat)) - log(sqrt(sigma_prior[l])) + log(1 - cdf1) - log(1 - cdf2) - (pow(mean_prior[l],2) / sigma_prior[l]) / 2 + pow(lambda_hat,2) / (2 * sigma2_hat);
    double log_p_hat = - log(1 + exp(-log(p[l])- p_w1 +log(1-p[l]) + p_w0));
    double p_hat = exp(log_p_hat);
    
    //print l and p_hat
    // printf("l is: %d\n", l);
    // printf("p_hat is: %f\n", p_hat);
    
    posterior_p[l] = p_hat;
    
    if (p_hat > 1) p_hat = 1;
    if (p_hat < 0.0001) p_hat = 0;
    
    //bernoulli sampler
    indicator[l] = R::rbinom(1, p_hat);
  // }
  // 
  // for(int l = 0; l < q; l++) {
  //   // printf("l is: %d\n", l);
  //   
  //   Eigen::MatrixXd Eta_matrix = Eigen::MatrixXd::Zero(id[n], q);
  //   for (int i = 0; i < n; i++) {
  //     for (int j = id[i]; j < id[i+1]; j++) {
  //       Eta_matrix.row(j) = Eigen::Map<Eigen::VectorXd>(&Eta[i][j-id[i]][0], q);
  //     }
  //   }
    // 
    // double med_val1 = T_matrix.col(l).transpose() * T_matrix.col(l);
    // double med_val2 = T_matrix.col(l).transpose() * Eta_matrix.col(l);
    // 
    // double inversew2 = med_val1 / sigma;
    // 
    // double lambda_tilde = med_val2 / med_val1;
    // 
    // double sigma2_hat = 1 / (inversew2 + 1 / sigma_prior[l]);
    // double lambda_hat = sigma2_hat * (inversew2 * lambda_tilde + (1 / sigma_prior[l]) * mean_prior[l]);
    
    if (indicator[l] == 1) {
      //truncated normal distribution (0, infinity)
      Rcpp::Environment truncnorm("package:truncnorm");
      Rcpp::Function rtruncnorm = truncnorm["rtruncnorm"];
      Rcpp::NumericVector lambda_temp = rtruncnorm(1, 0, R_PosInf, lambda_hat, sqrt(sigma2_hat));
      lambda[l] = lambda_temp[0];
      //lambda[l] = r_lefttruncnorm(0, lambda_hat, sqrt(sigma2_hat), seed);
    } else {
      lambda[l] = 0;
    }
  }
  return lambda;
}

//sample p

std::vector<double> sample_p(double *indicator, double *alpha, double *beta, int q, int seed) {  

  std::vector<double> p(q);
  for (int i = 0; i < q; i++) {
    p[i] = R::rbeta(alpha[i] + indicator[i], beta[i] + 1 - indicator[i]);
  }
  return p;
}

//sample Gamma (prior: N(mean,sigma2)) (checked 06.25)
std::vector<double> sample_gamma(double *mean_prior,double **sigma_prior, int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                                 std::vector<double> lambda,double *indicator, std::vector<std::vector<double>> b, double sigma, std::vector<std::vector<std::vector<double>>> U, int q, int seed) {
  //q(q-1)/2
  Eigen::VectorXd lambda_vector = Eigen::VectorXd::Zero(q);
  lambda_vector = Eigen::Map<Eigen::VectorXd>(&lambda[0], q);

  Eigen::MatrixXd med_mat = Eigen::MatrixXd::Zero(q * (q - 1) / 2, q * (q - 1) / 2);
  Eigen::VectorXd med_vec = Eigen::VectorXd::Zero(q * (q - 1) / 2);
  Eigen::MatrixXd b_matrix(n,q);
  for (int i = 0; i < n; i ++){
    b_matrix.row(i) = Eigen::Map<Eigen::VectorXd>(&b[i][0], q);
  }
  
  for(int i = 0; i < n; i++) {
    Eigen::MatrixXd z_matrix(id[i+1] - id[i],q);
    for (int z_i = 0; z_i < id[i+1] - id[i]; z_i ++) {
      z_matrix.row(z_i) = Eigen::Map<Eigen::VectorXd>(&z[i][z_i][0], q);
    }

    for(int j = 0; j < id[i+1] - id[i]; j++) {
      //every row of b and z become row vector
      Eigen::VectorXd U_vector = Eigen::VectorXd::Zero(q * (q - 1) / 2);
      U_vector = Eigen::Map<Eigen::VectorXd>(&U[i][j][0], q * (q - 1) / 2);
      // printf("U_vector is: \n");
      // for (int l = 0; l < q * (q - 1) / 2; l++) {
      //   printf("%f ", U_vector(l));
      // }
      // printf("\n");
      med_mat = med_mat + U_vector * U_vector.transpose();
      med_vec = med_vec + U_vector * (y[j + id[i]] - bart_predictions[j + id[i]] -
        (((b_matrix.row(i).array()) * (z_matrix.row(j).array()) * (lambda_vector.transpose().array())).matrix().sum()));
    }
  }
  
  Eigen::VectorXd mean_prior_vec = Eigen::Map<Eigen::VectorXd>(mean_prior, q * (q - 1) / 2);
  
  Eigen::MatrixXd sigma_prior_matrix = Eigen::MatrixXd::Zero(q * (q - 1) / 2, q * (q - 1) / 2);
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    for (int j = 0; j < q * (q - 1) / 2; j++) {
      sigma_prior_matrix(i,j) = sigma_prior[i][j];
    }
  }
  
  // //print med_mat
  // printf("med_mat is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%f ", med_mat(i, j));
  //   }
  //   printf("\n");
  // }
  
  // printf("med_vec is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   printf("%f ", med_vec(i));
  // }
  // printf("\n");
  // printf("med_mat is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%f ", med_mat(i, j));
  //   }
  //   printf("\n");
  // }
  // 
  Eigen::MatrixXd med_mat_copy = Eigen::MatrixXd::Zero(q * (q - 1) / 2, q * (q - 1) / 2);
  med_mat_copy = (1 / sigma * med_mat + sigma_prior_matrix.inverse());
  // printf("med_mat_copy is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%f ", med_mat_copy(i, j));
  //   }
  //   printf("\n");
  // }
  // 
  
  Eigen::MatrixXd covariance_matrix = (1 / sigma * med_mat + sigma_prior_matrix.inverse()).llt().solve(Eigen::MatrixXd::Identity(q * (q - 1) / 2, q * (q - 1) / 2));
  
  // printf("covariance_matrix is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%f ", covariance_matrix(i, j));
  //   }
  //   printf("\n");
  // }
  // 
  Eigen::VectorXd mean_vector = covariance_matrix * (1 / sigma * med_vec + sigma_prior_matrix.inverse() * mean_prior_vec);
  
  
  //sample gamma from the multivariate normal distribution
  // printf("sample gamma using MASS");
  Rcpp::Environment mvn("package:MASS");
  Rcpp::Function mvrnorm = mvn["mvrnorm"];
  Rcpp::NumericVector mean_vector_R(q * (q - 1) / 2);
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    mean_vector_R[i] = mean_vector(i);
  }
  Rcpp::NumericMatrix covariance_matrix_R(q * (q - 1) / 2, q * (q - 1) / 2);
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    for (int j = 0; j < q * (q - 1) / 2; j++) {
      covariance_matrix_R(i, j) = covariance_matrix(i, j);
    }
  }
  
  // //printf the mean and covariance matrix
  // printf("mean_vector is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   printf("%f ", mean_vector(i));
  // }
  // printf("\n");
  // printf("covariance_matrix is: \n");
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   for (int j = 0; j < q * (q - 1) / 2; j++) {
  //     printf("%f ", covariance_matrix(i, j));
  //   }
  //   printf("\n");
  // }
  
  Rcpp::NumericVector gamma_vector = mvrnorm(1, mean_vector_R, covariance_matrix_R);
  std::vector<double> gamma(q * (q - 1) / 2);
  for (int i = 0; i < q * (q - 1) / 2; i++) {
    gamma[i] = gamma_vector[i];
  }
  
  // printf("Original gamma is: \n"); 
  // for (int i = 0; i < q * (q - 1) / 2; i++) {
  //   printf("%f ", gamma[i]);
  // }
  // printf("\n");
  
  std::vector<std::vector<double>> Gamma = low_matrix(gamma, q);
  //if lambda[l] is 0, then Gamma[l,] and Gamma[,l] should be 0, diagnoal elements should be 1
  for (int i = 0; i < q; i++) {
    if (lambda[i] < 1e-16) {
      for (int j = 0; j < i; j++) {
        Gamma[i][j] = 0;
      }
      for (int j = i; j < q; j++) {
        Gamma[j][i] = 0;
      }
      Gamma[i][i] = 1;
    }
  }
  
  //get new gamma from the low triangle Gamma
  int temp = 0;
  for (int i = 1; i < q; i++) {
    for (int j = 0; j < i; j++) {
      gamma[temp] = Gamma[i][j];
      temp ++;
    }
  }
  
  return gamma;
}

// //[[Rcpp::export]]
// SEXP random_effect(SEXP sexp_n, SEXP sexp_id, SEXP sexp_p, SEXP sexp_q, SEXP sexp_y, SEXP sexp_bart_prediction, SEXP sexp_Z,
//                    SEXP sexp_sigma, SEXP sexp_c0, SEXP sexp_d0,
//                    SEXP sexp_b,
//                    SEXP sexp_gamma, SEXP sexp_gamma_mean, SEXP sexp_gamma_cov,
//                    SEXP sexp_lambda, SEXP sexp_lambda_mean, SEXP sexp_lambda_cov, SEXP sexp_p,
//                    SEXP iters, SEXP burn,
//                    SEXP sexp_seed) {
//   int n = as<int>(sexp_n);
//   
//   Rcpp::IntegerVector idv(sexp_id);
//   int *id = &idv[0];
// 
//   int p = as<int>(sexp_p);
//   int q = as<int>(sexp_q);
// 
//   Rcpp::NumericVector yv(sexp_y);
//   double *y = &yv[0];
// 
//   Rcpp::NumericVector bart_predictionv(sexp_bart_prediction);
//   double *bart_prediction = &bart_predictionv[0];
// 
//   Rcpp::NumericVector Zv(sexp_Z);
//   double **Z_2d = new double*[id[n]];
//   for (int i = 0; i < id[n]; i++) {
//     Z_2d[i] = new double[q];
//     for (int j = 0; j < q; j++) {
//       Z_2d[i][j] = Zv[i * q + j];
//     }
//   }
//   double ***Z_3d = new double**[n];
//   for (int i = 0; i < n; i++) {
//     Z_3d[i] = new double*[id[i+1] - id[i]];
//     for (int j = 0; j < id[i+1] - id[i]; j++) {
//       Z_3d[i][j] = new double[q];
//       for (int l = 0; l < q; l++) {
//         Z_3d[i][j][l] = Z_2d[j+id[i]][l];
//       }
//     }
//   }
// 
//   double c0 = as<double>(sexp_c0);
//   double d0 = as<double>(sexp_d0);
// 
//   Rcpp::NumericVector gamma_meanv(sexp_gamma_mean);
//   double *gamma_mean = &gamma_meanv[0];
// 
//   Rcpp::NumericVector gamma_covv(sexp_gamma_cov);
//   double **gamma_cov = new double*[q * (q - 1) / 2];
//   for (int i = 0; i < q * (q - 1) / 2; i++) {
//     gamma_cov[i] = new double[q * (q - 1) / 2];
//     for (int j = 0; j < q * (q - 1) / 2; j++) {
//       gamma_cov[i][j] = gamma_covv[i * q * (q - 1) / 2 + j];
//     }
//   }
//   
//   double *indicator = new double[q];
//   for (int i = 0; i < q; i++) {
//     indicator[i] = 1;
//   }
// 
//   Rcpp::NumericVector lambda_meanv(sexp_lambda_mean);
//   double *lambda_mean = &lambda_meanv[0];
// 
//   Rcpp::NumericVector lambda_covv(sexp_lambda_cov);
//   double *lambda_cov = &lambda_covv[0];
// 
//   Rcpp::NumericVector pv(sexp_p);
//   double *p = &pv[0];
//   
//   int niters = as<int>(iters);
//   int burnin = as<int>(burn);
// 
//   int seed = as<int>(sexp_seed);
//   
//   //initialize the random effect
//   double sigma = as<double>(sexp_sigma);
//   
//   Rcpp::NumericVector b_v(sexp_b);
//   double **b = new double*[n];
//   for (int i = 0; i < n; i++) {
//     b[i] = new double[q];
//     for (int j = 0; j < q; j++) {
//       b[i][j] = b_v[i * q + j];
//     }
//   }
//   
//   Rcpp::NumericVector gamma_v(sexp_gamma);
//   double *gamma = new double[q * (q - 1) / 2];
//   for (int i = 0; i < q * (q - 1) / 2; i++) {
//     gamma[i] = gamma_v[i];
//   }
//   double **Gamma = low_matrix(gamma, q);
//   
//   Rcpp::NumericVector lambda_v(sexp_lambda);
//   double *lambda = new double[q];
//   for (int i = 0; i < q; i++) {
//     lambda[i] = lambda_v[i];
//   }
//   
//   for (int i = 0; i < q; i++) {
//     if (lambda[i] == 0) {
//       for (int j = 0; j < i; j++) {
//         Gamma[i][j] = 0;
//       }
//       for (int j = i; j < q; j++) {
//         Gamma[j][i] = 0;
//       }
//       Gamma[i][i] = 1;
//     }
//   }
//   
//   int temp = 0;
//   for (int i = 1; i < q; i++) {
//     for (int j = 0; j < i; j++) {
//       gamma[temp] = Gamma[i][j];
//       temp ++;
//     }
//   }
// 
//   //save result after burnin th iterations
//   Rcpp::NumericMatrix sigma_save(niters - burnin, 1);
//   Rcpp::NumericMatrix b_save(niters - burnin, n * q);
//   Rcpp::NumericMatrix gamma_save(niters - burnin, q * (q - 1) / 2);
//   Rcpp::NumericMatrix lambda_save(niters - burnin, q);
//   Rcpp::NumericMatrix D_save(niters - burnin, q);
//   
//   double *random_effect_prediction =  rf_predict(id[n], n, id, Z_3d, lambda, Gamma, b, q);
//     
//   // 
//   for (int iter = 0; iter < niters; iter++) {
//     printf("iter: %d\n", iter);
//     //sample sigma
//     sigma = sample_sigma(c0, d0, id[n], n, id, y, bart_prediction, Z_2d, b, Gamma, lambda, q, seed);
//     printf("sigma: %f\n", sigma);
//     //
//     // sample b
//     b = sample_b(id[n], n, id, y, bart_prediction, Z_3d, Gamma, lambda, sigma, q, seed);
//      
//     //calculate U
//     double ***U = calculate_U(n, id, y, bart_prediction, Z_3d, Gamma, lambda, b, sigma, q);
// 
//     //sample gamma
//     gamma = sample_gamma(gamma_mean, gamma_cov, n, id, y, bart_prediction, Z_3d, lambda, b, sigma, U, q, seed);
//     printf("gamma: ");
//     printf("\n");
//     for (int i = 0; i < q * (q - 1) / 2; i++) {
//       printf("%f ", gamma[i]);
//     }
//     printf("\n");
// 
//     Gamma = low_matrix(gamma, q);
// 
//     //sample lambda
//     lambda = sample_lambda(lambda_mean, lambda_cov, indicator, p, n, id, y, bart_prediction, Z_3d, Gamma, lambda, b, sigma, q, seed);
//     printf("lambda: ");
//     printf("\n");
//     for (int i = 0; i < q; i++) {
//       printf("%f ", lambda[i]);
//     }
//     printf("\n");
//     
//     //C++: diag_D <- (diag(lambda) %*% Gamma) %*% t(Gamma) %*% diag(lambda)
//     Eigen::MatrixXd Lambda = Eigen::MatrixXd::Zero(q, q);
//     for (int i = 0; i < q; i++) {
//       Lambda(i, i) = lambda[i];
//     }
//     Eigen::MatrixXd Gamma_matrix = Eigen::MatrixXd::Zero(q, q);
//     for (int i = 0; i < q; i++) {
//       for (int j = 0; j < q; j++) {
//         Gamma_matrix(i, j) = Gamma[i][j];
//       }
//     }
//     Eigen::MatrixXd D = (Lambda * Gamma_matrix) * Gamma_matrix.transpose() * Lambda;
//     //get diagonal of D
//     double *D_diag = new double[q];
//     for (int i = 0; i < q; i++) {
//       D_diag[i] = D(i, i);
//     }
// 
//     //when iter > 5000, save the result for the following iterations
//     if (iter >= burnin) {
//       sigma_save(iter - burnin, 0) = sigma;
//       for (int i = 0; i < n; i++) {
//         for (int j = 0; j < q; j++) {
//           b_save(iter - burnin, i * q + j) = b[i][j];
//         }
//       }
//       for (int i = 0; i < q * (q - 1) / 2; i++) {
//         gamma_save(iter - burnin, i) = gamma[i];
//       }
//       for (int i = 0; i < q; i++) {
//         lambda_save(iter - burnin, i) = lambda[i];
//       }
//       for (int i = 0; i < q; i++) {
//         D_save(iter - burnin, i) = D_diag[i];
//       }
//     }
//     
//   }
//   
//   return Rcpp::List::create(Rcpp::Named("sigma") = sigma_save,
//                             Rcpp::Named("b") = b_save,
//                             Rcpp::Named("gamma") = gamma_save,
//                             Rcpp::Named("lambda") = lambda_save,
//                             Rcpp::Named("D") = D_save);
//   return R_NilValue;
// }
// 
// 
