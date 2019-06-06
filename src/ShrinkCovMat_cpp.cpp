#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat centerdata(arma::mat X) {
  int N = X.n_cols;
  arma::mat Ans = X;
  arma::mat X_mean = mean(X, 1);
  for(int i=1; i<N+1; i++){
    Ans.col(i-1) -= X_mean;
  }
  return Ans;
}




#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec trace_stats_uncentered(arma::mat X) {
  int p = X.n_rows;
  int N = X.n_cols;
  arma::mat X_i = arma::zeros(p, 1);
  arma::mat X_i_square = arma::zeros(p, 1);
  arma::mat X_i_cube = arma::zeros(p, 1);
  arma::mat X_square = square(X);
  arma::mat X_cube = pow(X, 3);
  arma::mat X_sum = sum(X, 1);
  arma::mat X_mean = X_sum/N;
  arma::mat X_sum_minus_i = X_sum;
  arma::mat X_square_sum = sum(X_square, 1);
  arma::mat X_square_sum_minus_i = X_square_sum;
  arma::mat X_cube_sum_minus_i = sum(X_cube, 1);
  double Y_1N = accu(X_square_sum_minus_i)/N;
  double Y_4N = 0;
  double Y_2N = 0;
  double trace_S2_1 = accu(pow(sum(X_square, 0),2));
  double trace_S2_3 = 0;
  double trace_S2_4 = pow(N, 2) * pow(accu(square(X_mean)),2);
  double Q = accu(pow(sum(square(centerdata(X))), 2)) / (N - 1);
  double Y_3N = 0;
  arma::mat Sum1 = arma::zeros(p, 1);
  arma::mat Sum21 = arma::zeros(p, 1);
  arma::mat Sum22 = arma::zeros(p, 1);
  for(int i=1; i<N+1; i++){
    X_i = X.col(i-1);
    X_i_square = X_square.col(i-1);
    X_i_cube = X_cube.col(i-1);
    X_sum -= X_i;
    X_sum_minus_i -= X_i;
    Y_4N += accu(X_i % X_sum_minus_i);
    trace_S2_3 += pow(accu(X_i % X_mean), 2) ; 
    Sum1 += sum(X_i % X_sum, 1);
    Sum21 += sum(X_i_cube % X_sum, 1);
    X_cube_sum_minus_i -= X_i_cube;
    Sum22 += sum(X_cube_sum_minus_i % X_i, 1);
    X_square_sum_minus_i -= X_i_square;
    Y_3N += accu(X_square_sum_minus_i % X_i_square);
    for(int j=i+1; j<N+1; j++){
      Y_2N += pow(accu(X_i % X.col(j-1)), 2);
    }
  }
  Y_4N = 2 * Y_4N/N/(N-1);
  double T_1N = Y_1N - Y_4N;
  double trace_S2_2 = 2 * Y_2N;
  trace_S2_3 = -2 * N * trace_S2_3;
  double trace_S2 = (trace_S2_1 + trace_S2_2 + trace_S2_3 + trace_S2_4)/pow(N-1, 2);
  double T_2N = (N-1)*((N-1)*(N-2)*trace_S2 + pow(T_1N, 2) - N*Q)/(N*(N-2)*(N-3));
  double Y_7N = 2 * (accu(Sum1 % X_square_sum) - accu(Sum21 + Sum22));
  double Y_8N = 4 * (accu(square(Sum1)) - Y_3N - Y_7N);
  Y_3N = 2 * Y_3N/N/(N - 1);
  Y_7N = Y_7N/N/(N - 1)/(N - 2);
  Y_8N = Y_8N/N/(N - 1)/(N - 2)/(N - 3);
  double T_3N = Y_3N - 2 * Y_7N + Y_8N;
  arma::vec ans = arma::vec(3);
  ans[0] = T_1N;
  ans[1] = T_2N;
  ans[2] = T_3N;
  return ans;
}




#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec trace_stats_centered(arma::mat X) {
  int p = X.n_rows;
  int N = X.n_cols;
  arma::mat X_i = arma::zeros(p, 1);
  arma::mat X_i_square = arma::zeros(p, 1);
  arma::mat X_square_sum = sum(square(X), 1);
  double Y_1N = accu(X_square_sum)/N;
  double Y_2N = 0;
  double Y_3N = 0;
  for(int i=1; i<N; i++){
    X_i = X.col(i-1);
    X_i_square = square(X_i);
    X_square_sum -= X_i_square;
    Y_3N += accu(X_square_sum % X_i_square);
    for(int j=i+1; j<N+1; j++){
      Y_2N += pow(accu(X_i % X.col(j - 1)), 2);
    }
  }
  Y_2N = 2 * Y_2N/N/(N - 1);
  Y_3N = 2 * Y_3N/N/(N - 1);
  arma::vec ans = arma::vec(3);
  ans[0] = Y_1N;
  ans[1] = Y_2N;
  ans[2] = Y_3N;
  return ans;
}
