/*
 *
 */

#ifndef COMPUTE_ANOMALY_H
#define COMPUTE_ANOMALY_H

#endif //COMPUTE_ANOMALY_H

#include <chrono>
#include <armadillo>

using namespace arma;
void compute_anomaly_score(const Mat<double> &X,  // (N x M)
                           const Cube<double> &A, // (M x M x K)
                           const Mat<double> &m, // (K x M)
                        //    const Mat<double> &g_mat, // (N x K)
                           Mat<double> &theta, // (M x K)
                           // Return arguments
                           Mat<double> &anomaly_score, // (N, M)
                           bool verbose = false){

    int N = X.n_rows, M = X.n_cols, K = theta.n_cols;

    Cube<double> U(N, M, K, fill::zeros),
                 W(N, M, K, fill::zeros);

    std::chrono::steady_clock sc;

    if(verbose){
        std::cout<<termcolor::blue<<"============ Computing anomaly score ============="<<endl;
    }

    auto start = sc.now();     // start timer
    for(int k=0; k < K; ++k){
        W.slice(k).each_row() = 1/A.slice(k).diag().as_row();
    }

    Cube<double> tmp(N, M, K, fill::zeros);
    for (int  k=0; k < K; ++k) {
        tmp.slice(k) = (X.each_row() - m.row(k)) * A.slice(k); // NxM

        for (int i = 0; i < M; ++i) {
            U.slice(k).col(i) = X.col(i) - (tmp.slice(k).col(i) / A(i, i, k));
        }
    }

    // Compute the Gating function (Posterior Distribution)
    // Equation (17) in section 3.2
    Mat<double> g_unnorm(N, K, fill::zeros);
    Mat<double> g_mat(N, K, fill::zeros);
    // assert(g_mat.is_finite() && "GMRF GMat has NANs.");

    // Equation (22) in section 3.3
    anomaly_score.resize(N, M);
    vec score(N, fill::zeros);
    for(int i=0; i < M; ++i){
        score.fill(0.0);
        // for (int i=0; i< M; ++i){
        for(int k =0; k < K; ++k) {
            g_unnorm.col(k) = theta(k) * normpdf(X.col(i),
                                                            U.slice(k).col(i),
                                                            W.slice(k).col(i));
            g_mat.col(k) = g_unnorm.col(k);
        }
        // }
        vec denom = sum(g_mat, 1); // Sum each row
        g_mat.each_col() /= (denom + 1e-8);

        for(int k=0; k < K; ++k){
            score += g_mat.col(k) % normpdf(X.col(i),
                                            U.slice(k).col(i),
                                            W.slice(k).col(i));
        }
        anomaly_score.col(i) = -log(score);
    }
    auto end = sc.now();
    auto time_span = static_cast<std::chrono::duration<double>>(end - start);
    // std::cout <<termcolor::blue<<"Operation took: " << time_span.count() << " seconds"<<endl;
    if(verbose){
        std::cout<<termcolor::blue<<"=================================================="<<endl;
    }
}