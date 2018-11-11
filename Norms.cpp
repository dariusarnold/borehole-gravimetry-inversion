//
// Created by darius on 11/11/18.
//

#include <Eigen/Dense>
#include <vector>
#include "Norms.h"
#include "utils.h"
#include "MeasurementData.h"



Eigen::MatrixXd Norm::gram_matrix_analytical(const std::vector<double>& depth) {
    Eigen::MatrixXd gram_matrix(depth.size(), depth.size());
    for (size_t column_index = 0; column_index < depth.size(); ++column_index){
        for (size_t row_index = 0; row_index < depth.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(depth[row_index], depth[column_index]);
        }
    }
    return gram_matrix;
}


std::vector<Result> Norm::calculate_density_distribution(const std::vector<double> &alpha, const std::vector<double> &depth, uint64_t num_steps) {
    // fill depth vector with ascending values
    std::vector<double>  depth_meters;
    depth_meters.reserve(num_steps);
    double stepsize = depth.back() / num_steps;
    for (size_t i = 0; i != num_steps; ++i){
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


double L2_Norm::gram_entry_analytical(double zj, double zk) {
    return gamma*gamma*std::min(zj, zk);
}


double L2_Norm::representant_function(double zj, double z){
    return -gamma*heaviside(zj-z);
}


double W12_Norm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}


double W12_Norm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}
