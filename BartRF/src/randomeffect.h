#include <Rcpp.h>
#include <eigen3/Eigen/Dense>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
using namespace Rcpp;

std::vector<std::vector<double>> low_matrix(std::vector<double> gamma, int p);

std::vector<double> rf_predict(int N, int n, int *id, std::vector<std::vector<std::vector<double>>> z,
                               std::vector<double> lambda, std::vector<std::vector<double>> Gamma, std::vector<std::vector<double>> b, int q);

double sample_sigma(double lambda, double nu, int N, int n, int *id, double *y, double *bart_predictions,
                    std::vector<double> rf_predictions, int q, int seed);

std::vector<std::vector<double>> sample_b(int N, int n, int *id, double *y, double *bart_predictions,
                                          std::vector<std::vector<std::vector<double>>> z, std::vector<std::vector<double>> Gamma, std::vector<double> lambda,
                                          double sigma, int q, int seed);

std::vector<std::vector<std::vector<double>>> calculate_U(int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                                                          std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>>b, double sigma, int q);

std::vector<std::vector<std::vector<double>>> calculate_T(int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                                                          std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>> b, double sigma, int q);

std::vector<std::vector<std::vector<double>>> calculate_Eta(int n, int *id, double *y, double *bart_predictions,
                                                            std::vector<double> lambda, std::vector<std::vector<std::vector<double>>> T, int q);

std::vector<double> sample_lambda(double *mean_prior,double *sigma_prior,double *indicator, std::vector<double> p, double *posterior_p,
                                  int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                                  std::vector<std::vector<std::vector<double>>> T,
                                  std::vector<std::vector<double>> Gamma, std::vector<double> lambda, std::vector<std::vector<double>> b, double sigma, int q, int seed);

std::vector<double> sample_p(double *indicator, double *alpha, double *beta, int q, int seed);

std::vector<double> sample_gamma(double *mean_prior,double **sigma_prior, int n, int *id, double *y, double *bart_predictions, std::vector<std::vector<std::vector<double>>> z,
                                 std::vector<double> lambda,double *indicator, std::vector<std::vector<double>> b, double sigma, std::vector<std::vector<std::vector<double>>> U, int q, int seed);
