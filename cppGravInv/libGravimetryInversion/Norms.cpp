#include <iostream>
#include "Norms.h"
#include "MeasurementData.h"
#include "utils.h"
#include "ModelParameters.h"


Norm::Norm(const std::vector<double>& _measurement_depths, const std::vector<double>& _measurement_data) :
    measurement_depths(_measurement_depths),
    measurement_data(_measurement_data),
    gram_matrix(),
    alpha() {}

void Norm::gram_matrix_analytical() {
    gram_matrix.resize(measurement_depths.size(), measurement_depths.size());
    for (size_t column_index = 0; column_index < measurement_depths.size(); ++column_index){
        for (size_t row_index = 0; row_index < measurement_depths.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(measurement_depths[row_index], measurement_depths[column_index]);
        }
    }
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> Norm::calculate_density_distribution(uint64_t num_steps) {
    // fill depth vector with ascending values
    std::vector<double>  depth_meters;
    double stepsize = measurement_depths.back() / num_steps;
    // num_steps + 2 to get one step past the end to see the L2 norm falling to zero
    for (size_t i = 0; i != num_steps+2; ++i){
        depth_meters.push_back(stepsize * i);
    }
    // discretize density distribution by evaluating the following formula
    // rho(z) = sum_k alpha_k g_k(z)
    std::vector<double> depth;
    std::vector<double> density;
    for (auto discretization_depth : depth_meters){
        double dens = 0;
        for (long int j = 0; j != alpha.size(); ++j){
            dens += alpha[j] * representant_function(measurement_depths[j], discretization_depth);
        }
        depth.push_back(discretization_depth);
        density.push_back(dens);
    }
    //TODO create these as Eigen types instead of converting them
    Eigen::VectorXd depth_eigen = std_to_eigen(depth);
    Eigen::VectorXd density_eigen = std_to_eigen(density);
    return std::make_pair(depth_eigen, density_eigen);
}

void Norm::solve_for_alpha() {
    // initialize eigen::vector from std::vectors data containing measurement results corrected for free air gradient
    Eigen::VectorXd data_vec = std_to_eigen(measurement_data);
    // use Eigen to solve the matrix equation
    alpha = gram_matrix.colPivHouseholderQr().solve(data_vec);
}

std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> Norm::do_work(uint64_t discretization_steps) {
    gram_matrix_analytical();
    solve_for_alpha();
    auto result = calculate_density_distribution(discretization_steps);
    ModelParameters params = {-1, -1, -1};
    return std::make_tuple(result, params);

}



ErrorNorm::ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors) :
    Norm(depth, data),
    measurement_errors(errors),
    sigma_matrix(Eigen::MatrixXd::Constant(measurement_errors.size(), measurement_errors.size(), 0.)){
    // create sigma matrix from the measurement errors
    // fill sigma diagonal with values
    Eigen::VectorXd sigma_vec = std_to_eigen(measurement_errors);
    sigma_matrix.diagonal() = sigma_vec;
}

std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> ErrorNorm::do_work(uint64_t discretization_steps, double nu){
    // either nu was given from user or calculated using bisection
    solve_for_alpha(nu);
    auto density = calculate_density_distribution(discretization_steps);
    double misfit_squared_over_n = calculate_misfit(nu)/measurement_data.size();
    double norm = ErrorNorm::calculate_norm();
    ModelParameters params = {nu, misfit_squared_over_n, norm};
    return std::make_tuple(density, params);
}

std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> ErrorNorm::do_work(uint64_t discretization_steps) {
    // calculate the lagrange multiplicator using bysection
    double nu_left_start = 0.00001;
    double nu_right_start = 10000;
    // set threshold to number of data points
    double threshold = measurement_data.size();
    double nu = calc_nu_bysection(nu_left_start, nu_right_start, threshold);
    return ErrorNorm::do_work(discretization_steps, nu);
}

void ErrorNorm::solve_for_alpha(double nu) {
    auto sigma_squared = sigma_matrix * sigma_matrix;
    ErrorNorm::gram_matrix_analytical();
    // now set up the equation to be solved to calculate alpha
    Eigen::MatrixXd term = 1/nu * sigma_squared + gram_matrix;
    Eigen::VectorXd data_vec = std_to_eigen(measurement_data);
    alpha = term.colPivHouseholderQr().solve(data_vec);

}

double ErrorNorm::calculate_norm() {
    return alpha.transpose() * gram_matrix * alpha;
}

double ErrorNorm::calculate_misfit(double nu) {
    double norm = (sigma_matrix*alpha).norm();
    double misfit = (1./(nu*nu)) * norm*norm;
    return misfit;
}

double ErrorNorm::calc_nu_bysection(double nu_left, double nu_right, double desired_misfit, double accuracy) {
    // calc misfits for left and right end of interval to check if the middle is within this interval
    solve_for_alpha(nu_left);
    double misfit_left = calculate_misfit(nu_left);
    solve_for_alpha(nu_right);
    double misfit_right = calculate_misfit(nu_right);
    // if the misfit mid is outside of the interval spanned by the two start values, error
    if (misfit_left < desired_misfit){
        throw std::range_error("Start value nu left to big");
    }
    if (misfit_right > desired_misfit){
        throw std::range_error("Start value nu right to small");
    }
    // calc misfit for center of interval
    nu_left = log(nu_left);
    nu_right = log(nu_right);
    double misfit_mid, nu_mid;
    // else use bisection to search the optimal nu
    do {
        // calc misfit for center of interval
        nu_mid = (nu_right + nu_left)/2.;
        solve_for_alpha(exp(nu_mid));
        misfit_mid = calculate_misfit(exp(nu_mid));
        // compare with desired misfit and half interval size by taking the left or the right part
        misfit_mid > desired_misfit ? nu_left = nu_mid : nu_right = nu_mid;
    }while (std::abs(misfit_mid - desired_misfit) > accuracy);
    return exp(nu_mid);
}