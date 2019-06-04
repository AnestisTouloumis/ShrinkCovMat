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
arma::vec optimal_intensities_uncentered(arma::mat X, arma::mat S_matrix) {
  int p = X.n_rows;
  int N = X.n_cols;
  double T_1N = trace(S_matrix);
  double trace_S2 = accu(square(S_matrix));
  double Q = accu(pow(sum(square(centerdata(X))), 2)) / (N - 1);
  double T_2N = (N-1)*((N-1)*(N-2)*trace_S2 + pow(T_1N, 2) - N*Q)/(N*(N-2)*(N-3));
  arma::mat Sum1 = arma::zeros(p, 1);
  arma::mat Sum21 = arma::zeros(p, 1);
  arma::mat Sum22 = arma::zeros(p, 1);
  arma::mat Sum3 = arma::zeros(p, 1);
  arma::mat X_i = arma::zeros(p, 1);
  arma::mat X_i_square = arma::zeros(p, 1);
  arma::mat X_i_cube = arma::zeros(p, 1);
  arma::mat X_square = square(X);
  arma::mat X_cube = pow(X, 3);
  arma::mat X_sum = sum(X, 1);
  arma::mat X_mean = X_sum/N;
  arma::mat X_square_sum = sum(X_square, 1);
  arma::mat X_cube_sum = sum(X_cube, 1);
  for(int i=1; i<N+1; i++){
    X_i = X.col(i-1);
    X_i_square = X_square.col(i-1);
    X_i_cube = X_cube.col(i-1);
    X_sum -= X_i;
    Sum1 += sum(X_i % X_sum, 1);
    Sum21 += sum(X_i_cube % X_sum, 1);
    X_cube_sum -= X_i_cube;
    Sum22 += sum(X_cube_sum % X_i, 1);
    X_square_sum -= X_i_square;
    Sum3 += sum(X_square_sum % X_i_square, 1);
  }
  double Term1 = 2 * accu(Sum3)/N/(N - 1);
  double Term2 = 2 * (accu(Sum1 % sum(X_square, 1)) - accu(Sum21 + Sum22));
  double Term3 = 4 * (accu(square(Sum1)) - accu(Sum3) - Term2);
  Term2 = Term2/N/(N - 1)/(N - 2);
  Term3 = Term3/N/(N - 1)/(N - 2)/(N - 3);
  double T_3N = Term1 - 2 * Term2 + Term3;
  arma::vec ans = arma::vec(3);
  ans[0] = T_1N;
  ans[1] = T_2N;
  ans[2] = T_3N;
  return ans;
}




#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec optimal_intensities_centered(arma::mat X) {
  int p = X.n_rows;
  int N = X.n_cols;
  arma::mat X_i = arma::zeros(p, 1);
  arma::mat X_i_square = arma::zeros(p, 1);
  arma::mat X_square = square(X);
  arma::mat X_square_sum = sum(X_square, 1);
  double Y_1N = accu(X_square_sum)/N;
  double Sum1=0;
  arma::mat Sum2 = arma::zeros(p, 1);
  for(int i=1; i<N; i++){
    X_i = X.col(i-1);
    X_i_square = X_square.col(i-1);
    X_square_sum -= X_i_square;
    Sum2 += sum(X_square_sum % X_i_square, 1);
    for(int j=i+1; j<N+1; j++){
    Sum1 += pow(accu(X_i % X.col(j-1)), 2);
      }
    }
  double Y_2N = 2*Sum1/N/(N-1);
  double Y_3N = 2 * accu(Sum2)/N/(N - 1);
  arma::vec ans = arma::vec(3);
  ans[0] = Y_1N;
  ans[1] = Y_2N;
  ans[2] = Y_3N;
  return ans;
}
