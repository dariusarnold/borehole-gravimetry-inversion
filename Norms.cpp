//
// Created by darius on 11/11/18.
//

#include <Eigen/Dense>
#include <vector>
#include "Norms.h"
#include "utils.h"
#include "MeasurementData.h"


Eigen::MatrixXd L2_Norm::gram_matrix_analytical(const std::vector<double> &depth) {
    std::vector<double> gram_matrix_diag_elements;
    gram_matrix_diag_elements.reserve(depth.size());
    // only unique elements in the gram matrix are the diagonal elements because of
    // g_ij = min(g_i, g_j)
    // calculate these unique elements first and fill gram matrix with them later
    // use the analytical expression for the gram matrix Gamma_jk = gamma²*min(z_j, z_k)
    auto gram_analytical = [this](double zj, double zk){ return gamma*gamma*std::min(zj, zk); };
    for (auto el : depth){
        gram_matrix_diag_elements.emplace_back(gram_analytical(el, el));
    }
    // create matrix with as many columns/rows as data entries read from file
    // and fill it with the values from the diagonals in this pattern:
    // g11 g11 g11
    // g11 g22 g22
    // g11 g22 g33
    Eigen::MatrixXd gram_matrix(depth.size(), depth.size());
    for (size_t column_index = 0; column_index < depth.size(); ++column_index){
        for (size_t row_index = 0; row_index < depth.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_matrix_diag_elements[std::min(row_index, column_index)];
        }
    }
    return gram_matrix;
}

std::vector<Result> L2_Norm::calculate_density_distribution(const std::vector<double> &alpha, const std::vector<double >& depth, uint64_t num_steps) {
    // representant function of L2 norm
    auto representant_function = [this](double zj, double z){
        return -gamma*heaviside(zj-z);
    };
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
    for (size_t i = 0; i != num_steps; ++i){
        double dens = 0;
        for (size_t j = 0; j != alpha.size(); ++j){
            dens += alpha[j] * representant_function(depth[j], depth_meters[i]);
        }
        density.emplace_back(Result{depth_meters[i], dens});
    }
    return density;
}

/*
Eigen::MatrixXd W12_Norm::gram_matrix_analytical(const std::vector<double>& depth) {
    Eigen::MatrixXd gram_matrix(depth.size(), depth.size());
    for (size_t column_index = 0; column_index < depth.size(); ++column_index){
        for (size_t row_index = 0; row_index < depth.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(depth[row_index], depth[column_index]);
        }
    }
    return gram_matrix;
}

double W12_Norm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}

std::vector<std::pair<double, double>>  W12_Norm::calculate_density_distribution(const std::vector<double>& alpha, const std::vector<double>& representants, uint64_t num_steps) {
    return std::vector<std::pair<double, double>>();
}
*/