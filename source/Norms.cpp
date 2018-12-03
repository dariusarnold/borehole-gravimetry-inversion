//
// Created by darius on 11/11/18.
//

#include <Eigen/Dense>
#include <vector>
#include "Norms.h"
#include "utils.h"
#include "MeasurementData.h"


ErrorNorm::ErrorNorm() :
    gram_matrix(),
    sigma_matrix(),
    alpha(){}

ErrorNorm::~ErrorNorm() {}

void ErrorNorm::do_work(const std::vector<double>& depth, const std::vector<double>& data, double nu, const std::vector<double>& sigma, double threshold_squared){
    // create sigma² matrix from the measurement errors
    sigma_matrix.resize(sigma.size(), sigma.size());
    // fill sigma squared with values
    Eigen::VectorXd sigma_vec = Eigen::Map<const Eigen::VectorXd>(sigma.data(), sigma.size());
    sigma_matrix.diagonal() = sigma_vec;
    auto sigma_squared = sigma_matrix * sigma_matrix;

    ErrorNorm::gram_matrix_analytical(depth);
    // now set up the equation to be solved to calculate alpha
    Eigen::MatrixXd term = 1/nu * sigma_squared + gram_matrix;
    Eigen::VectorXd data_vec = Eigen::Map<const Eigen::VectorXd>(data.data(), data.size());
    Eigen::VectorXd alpha_eigen = term.colPivHouseholderQr().solve(data_vec);
    //convert result from Eigen type to std::vector and return it
    alpha = std::vector<double>(&alpha_eigen[0], alpha_eigen.data()+alpha_eigen.cols()*alpha_eigen.rows());
    std::cout <<  "Misfit squared/N: " << ErrorNorm::calculate_misfit(nu)/data.size() << std::endl;
    std::cout << "Norm: " << ErrorNorm::calculate_norm() << std::endl;
}

std::vector<Result> ErrorNorm::calculate_density_distribution(const std::vector<double> &depth, uint64_t num_steps) {
    // same as Norm::calculate_density_distribution
    // fill depth vector with ascending values
    std::vector<double>  depth_meters;
    depth_meters.reserve(num_steps);
    double stepsize = depth.back() / num_steps;
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
            dens += alpha[j] * representant_function(depth[j], discretization_depth);
        }
        density.emplace_back(Result{discretization_depth, dens});
    }
    return density;
}

void ErrorNorm::gram_matrix_analytical(const std::vector<double> &depth) {
    // same as Norm::gram_matrix_analytical
    gram_matrix.resize(depth.size(), depth.size());
    for (size_t column_index = 0; column_index < depth.size(); ++column_index){
        for (size_t row_index = 0; row_index < depth.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(depth[row_index], depth[column_index]);
        }
    }
}

double ErrorNorm::calculate_norm() {
    //TODO Save alpha as eigen type and only transform to vector later
    return std_to_eigen(alpha).transpose() * gram_matrix * std_to_eigen(alpha);
}

double ErrorNorm::calculate_misfit(double nu) {
    double norm = (sigma_matrix*std_to_eigen(alpha)).norm();
    return (1./(nu*nu)) * norm*norm;
}

L2ErrorNorm::L2ErrorNorm() {}

L2ErrorNorm::~L2ErrorNorm() {}

double L2ErrorNorm::representant_function(double zj, double z) {
    // same as L2Norm
    return -gamma*heaviside(zj-z);
}

double L2ErrorNorm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}

void L2ErrorNorm::do_work(const std::vector<double>& depth, const std::vector<double>& data, double nu, const std::vector<double>& sigma, double threshold_squared) {
    ErrorNorm::do_work(depth, data, nu, sigma, threshold_squared);
}

