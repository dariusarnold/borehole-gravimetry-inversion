#include "Norms.h"
#include "utils.h"



LinearInterpolationNorm::LinearInterpolationNorm(double _a, double _b, const std::vector<double>& x_values, const std::vector<double>& y_values) :
        Norm(x_values, y_values),
        a(_a),
        b(_b){}

void LinearInterpolationNorm::gram_matrix_analytical() {
    // get top left gram matrix, containing Gamma_jk
    Norm::gram_matrix_analytical();
    // calculate the additional data given as (gj, 1)
    Eigen::VectorXd additional(gram_matrix.cols()+1);
    for (size_t i = 0; i < measurement_depths.size(); ++i){
        additional(i) = 1.;
    }
    additional(additional.size()-1) = 0.;
    // add one more row and column to the matrix
    gram_matrix.conservativeResize(gram_matrix.rows()+1, gram_matrix.cols()+1);
    // fill additional space, Gram matrix is symmetrical so the same vector is inserted twice
    gram_matrix.rightCols(1) = additional;
    gram_matrix.bottomRows(1) = additional.transpose();
}

double LinearInterpolationNorm::gram_entry_analytical(double xj, double xk) {
    return 1. + std::min(xj, xk) - a;
}

double LinearInterpolationNorm::representant_function(double xj, double x) {
    return -ramp(xj-x) + 1. + xj - a;
}

void LinearInterpolationNorm::solve_for_alpha() {
    // extend data by the additional constant 0
    measurement_data.emplace_back(0.);
    Norm::solve_for_alpha();
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> LinearInterpolationNorm::calculate_density_distribution(uint64_t num_steps) {
    // fill depth vector with ascending values
    std::vector<double>  interpolated_x_values;
    interpolated_x_values.reserve(num_steps);
    double stepsize = (b-a)  / num_steps;
    // num_steps + 2 to get one step past the end to see the L2 norm falling to zero
    for (size_t i = 0; i != num_steps+2; ++i){
        interpolated_x_values.emplace_back(stepsize * i);
    }
    // discretize density distribution by evaluating the following formula
    // rho(z) = sum_k alpha_k g_k(z)
    std::vector<double> x;
    std::vector<double> y;
    for (auto x_value : interpolated_x_values){
        // get beta out of alpha
        double f_of_x = alpha(alpha.rows()-1);

        for (long int j = 0; j != alpha.size()-1; ++j){
            f_of_x += alpha[j] * representant_function(measurement_depths[j], x_value);
        }
        x.push_back(x_value);
        y.push_back(f_of_x);
    }
    Eigen::VectorXd x_eigen = std_to_eigen(x);
    Eigen::VectorXd y_eigen = std_to_eigen(y);
    return std::make_pair(x_eigen, y_eigen);
}