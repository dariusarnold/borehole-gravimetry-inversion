#include <iostream>
#include "Norms.h"
#include "MeasurementData.h"
#include "utils.h"


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

std::vector<Result> Norm::calculate_density_distribution(uint64_t num_steps) {
    // fill depth vector with ascending values
    std::vector<double>  depth_meters;
    depth_meters.reserve(num_steps);
    double stepsize = measurement_depths.back() / num_steps;
    // num_steps + 2 to get one step past the end to see the L2 norm falling to zero
    for (size_t i = 0; i != num_steps+2; ++i){
        depth_meters.emplace_back(stepsize * i);
    }
    // discretize density distribution by evaluating the following formula
    // rho(z) = sum_k alpha_k g_k(z)
    std::vector<Result> density;
    density.reserve(num_steps);
    for (auto discretization_depth : depth_meters){
        double dens = 0;
        for (long int j = 0; j != alpha.size(); ++j){
            dens += alpha[j] * representant_function(measurement_depths[j], discretization_depth);
        }
        density.emplace_back(Result{discretization_depth, dens});
    }
    return density;
}

void Norm::solve_for_alpha() {
    // initialize eigen::vector from std::vectors data containing measurement results corrected for free air gradient
    Eigen::Map<const Eigen::VectorXd> data_vec(measurement_data.data(), measurement_data.size());
    // use Eigen to solve the matrix equation
    alpha = gram_matrix.colPivHouseholderQr().solve(data_vec);
}

std::vector<Result> Norm::do_work(uint64_t discretization_steps) {
    gram_matrix_analytical();
    solve_for_alpha();
    auto density = calculate_density_distribution(discretization_steps);
    return density;
}



ErrorNorm::ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors) :
    Norm(depth, data),
    measurement_errors(errors),
    sigma_matrix(){}

std::vector<Result> ErrorNorm::do_work(double nu, uint64_t discretization_steps){
    // create sigma² matrix from the measurement errors
    sigma_matrix.resize(measurement_errors.size(), measurement_errors.size());
    // fill sigma squared with values
    Eigen::VectorXd sigma_vec = Eigen::Map<const Eigen::VectorXd>(measurement_errors.data(), measurement_errors.size());
    sigma_matrix.diagonal() = sigma_vec;
    solve_for_alpha(nu);
    auto density = calculate_density_distribution(discretization_steps);
    std::cout <<  "Misfit squared/N: " << ErrorNorm::calculate_misfit(nu)/measurement_data.size() << std::endl;
    std::cout << "Norm: " << ErrorNorm::calculate_norm() << std::endl;
    return density;
}

std::vector<Result> ErrorNorm::do_work(uint64_t discretization_steps) {
    // create sigma² matrix from the measurement errors
    sigma_matrix.resize(measurement_errors.size(), measurement_errors.size());
    // fill sigma squared with values
    Eigen::VectorXd sigma_vec = std_to_eigen(measurement_errors);
    sigma_matrix.diagonal() = sigma_vec;
    // calculate the lagrange multiplicator using bysection
    double nu_left_start = 0.01;
    double nu_right_start = 10000;
    double threshold = measurement_data.size();
    double nu = calc_nu_bysection(nu_left_start, nu_right_start, threshold);
    solve_for_alpha(nu);
    auto density = calculate_density_distribution(discretization_steps);
    std::cout << "Nu: " << nu << std::endl;
    std::cout <<  "Misfit squared/N: " << ErrorNorm::calculate_misfit(nu)/measurement_data.size() << std::endl;
    std::cout << "Norm: " << ErrorNorm::calculate_norm() << std::endl;
    return density;
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

double ErrorNorm::calc_nu_bysection(double nu_left, double nu_right, double desired_misfit) {
    nu_left = log(nu_left);
    nu_right = log(nu_right);
    double accuracy = 0.01;
    //TODO check whether desired misfit is within the range of left/right
    /*
    // calc misfits for left and right end of interval
    solve_for_alpha(nu_left);
    double misfit_left = calculate_misfit(nu_left);
    solve_for_alpha(nu_right);
    double misfit_right = calculate_misfit(nu_right);
     */
    double nu_mid, misfit_mid;
    do {
        // calc misfit for center of interval
        nu_mid = (nu_right + nu_left)/2.;
        solve_for_alpha(exp(nu_mid));
        misfit_mid = calculate_misfit(exp(nu_mid));
        // compare with desired misfit and half interval size by taking the left or the right part
        misfit_mid > desired_misfit ? nu_left = nu_mid : nu_right = nu_mid;
    }while ((misfit_mid - desired_misfit) > accuracy);
    return exp(nu_mid);
}