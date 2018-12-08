//
// Created by darius on 11/11/18.
//

#include <Eigen/Dense>
#include <vector>
#include "Norms.h"
#include "utils.h"
#include "MeasurementData.h"



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
        for (size_t j = 0; j != alpha.size(); ++j){
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
    Eigen::VectorXd alpha_eigen = gram_matrix.colPivHouseholderQr().solve(data_vec);
    //convert result from Eigen type to std::vector and return it
    alpha = std::vector<double>(&alpha_eigen[0], alpha_eigen.data()+alpha_eigen.cols()*alpha_eigen.rows());
}

std::vector<Result> Norm::do_work(uint64_t discretization_steps) {
    gram_matrix_analytical();
    solve_for_alpha();
    auto density = calculate_density_distribution(discretization_steps);
    return density;
}


L2_Norm::L2_Norm(const std::vector<double> &depth, const std::vector<double> &data) :
    Norm(depth, data){}

double L2_Norm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}


double L2_Norm::representant_function(double zj, double z){
    return -gamma*heaviside(zj-z);
}


W12_Norm::W12_Norm(const std::vector<double> &depth, const std::vector<double> &data) :
        Norm(depth, data){}

double W12_Norm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}


double W12_Norm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}


double Seminorm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}


Seminorm::Seminorm(const std::vector<double>& depth, const std::vector<double>& data) :
    Norm(depth, data){}


double Seminorm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}


void Seminorm::gram_matrix_analytical() {
    // get top left gram matrix, containing Gamma_jk
    Norm::gram_matrix_analytical();
    // calculate the additional data given as (gj, 1)
    Eigen::VectorXd additional(gram_matrix.cols()+1);
    for (size_t i = 0; i < measurement_depths.size(); ++i){
        additional(i) = -gamma*measurement_depths[i];
    }
    additional(additional.size()-1) = 0.;
    // add one more row and column to the matrix
    gram_matrix.conservativeResize(gram_matrix.rows()+1, gram_matrix.cols()+1);
    // fill additional space, Gram matrix is symmetrical so the same vector is inserted twice
    gram_matrix.rightCols(1) = additional;
    gram_matrix.bottomRows(1) = additional.transpose();
}


void Seminorm::solve_for_alpha() {
    // extend data by the additional constant 0
    measurement_data.emplace_back(0.);
    Norm::solve_for_alpha();

}


std::vector<Result> Seminorm::calculate_density_distribution(uint64_t num_steps) {
    std::vector<Result> density_variable = Norm::calculate_density_distribution(num_steps);
    auto density_constant = alpha.back();
    std::for_each(density_variable.begin(), density_variable.end(), [density_constant](Result& x){x.density +=density_constant; });
    return density_variable;
}



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

std::vector<Result> LinearInterpolationNorm::calculate_density_distribution(uint64_t num_steps) {
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
    std::vector<Result> discretized_interpolated_function;
    discretized_interpolated_function.reserve(num_steps);
    for (auto x_value : interpolated_x_values){
        // get beta out of alpha
        double f_of_x = alpha.back();

        for (size_t j = 0; j != alpha.size()-1; ++j){
            f_of_x += alpha[j] * representant_function(measurement_depths[j], x_value);
        }
        discretized_interpolated_function.emplace_back(Result{x_value, f_of_x});
    }
    return discretized_interpolated_function;
}


ErrorNorm::ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors) :
    measurement_depths(depth),
    measurement_data(data),
    measurement_errors(errors),
    gram_matrix(),
    sigma_matrix(),
    alpha(){}

ErrorNorm::~ErrorNorm() {}

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
    Eigen::VectorXd sigma_vec = Eigen::Map<const Eigen::VectorXd>(measurement_errors.data(), measurement_errors.size());
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
    Eigen::VectorXd data_vec = Eigen::Map<const Eigen::VectorXd>(measurement_data.data(), measurement_data.size());
    alpha = term.colPivHouseholderQr().solve(data_vec);

}

std::vector<Result> ErrorNorm::calculate_density_distribution(uint64_t num_steps) {
    // same as Norm::calculate_density_distribution
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

void ErrorNorm::gram_matrix_analytical() {
    // same as Norm::gram_matrix_analytical
    gram_matrix.resize(measurement_depths.size(), measurement_depths.size());
    for (size_t column_index = 0; column_index < measurement_depths.size(); ++column_index){
        for (size_t row_index = 0; row_index < measurement_depths.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(measurement_depths[row_index], measurement_depths[column_index]);
        }
    }
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
    //TODO check whether desired misfit is within th range of left/right
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
        // compare with desired misfit and hhalf interval size by taking the left or the right part
        misfit_mid > desired_misfit ? nu_left = nu_mid : nu_right = nu_mid;
    }while ((misfit_mid - desired_misfit) > accuracy);
    return exp(nu_mid);
}


L2ErrorNorm::L2ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors) :
    ErrorNorm(depth, data, errors){}

L2ErrorNorm::~L2ErrorNorm() {}

double L2ErrorNorm::representant_function(double zj, double z) {
    // same as L2Norm
    return -gamma*heaviside(zj-z);
}

double L2ErrorNorm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}
