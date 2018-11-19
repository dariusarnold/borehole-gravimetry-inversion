//
// Created by darius on 11/11/18.
//

#include <Eigen/Dense>
#include <vector>
#include "Norms.h"
#include "utils.h"
#include "MeasurementData.h"



void Norm::gram_matrix_analytical(const std::vector<double>& depth) {
    gram_matrix.resize(depth.size(), depth.size());
    for (size_t column_index = 0; column_index < depth.size(); ++column_index){
        for (size_t row_index = 0; row_index < depth.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_entry_analytical(depth[row_index], depth[column_index]);
        }
    }
}


std::vector<Result> Norm::calculate_density_distribution(const std::vector<double> &depth, uint64_t num_steps) {
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


void Norm::solve_for_alpha(const std::vector<double> &data) {
    // initialize eigen::vector from std::vectors data containing measurement results corrected for free air gradient
    Eigen::Map<const Eigen::VectorXd> data_vec(data.data(), data.size());
    // use Eigen to solve the matrix equation
    Eigen::VectorXd alpha_eigen = gram_matrix.colPivHouseholderQr().solve(data_vec);
    //convert result from Eigen type to std::vector and return it
    alpha = std::vector<double>(&alpha_eigen[0], alpha_eigen.data()+alpha_eigen.cols()*alpha_eigen.rows());
}

void Norm::do_work(const std::vector<double>& depth, const std::vector<double>& data) {
    gram_matrix_analytical(depth);
    solve_for_alpha(data);
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


double Seminorm::gram_entry_analytical(double zj, double zk) {
    double zmin = std::min(zj, zk);
    double gamma_square = gamma*gamma;
    return gamma_square*zj*zk + gamma_square*(1./3. * std::pow(zmin, 3) - 1./2. * std::pow(zmin, 2) * (zj+zk) + zmin*zj*zk);
}



double Seminorm::representant_function(double zj, double z) {
    return gamma/2. * (zj - z) * (zj - z) * heaviside(zj-z) - gamma * (zj + 1./2. * zj*zj);
}


void Seminorm::gram_matrix_analytical(const std::vector<double> &depth) {
    // get top left gram matrix, containing Gamma_jk
    Norm::gram_matrix_analytical(depth);
    // calculate the additional data given as (gj, 1)
    Eigen::VectorXd additional(gram_matrix.cols()+1);
    for (size_t i = 0; i < depth.size(); ++i){
        additional(i) = -gamma*depth[i];
    }
    additional(additional.size()-1) = 0.;
    // add one more row and column to the matrix
    gram_matrix.conservativeResize(gram_matrix.rows()+1, gram_matrix.cols()+1);
    // fill additional space, Gram matrix is symmetrical so the same vector is inserted twice
    gram_matrix.rightCols(1) = additional;
    gram_matrix.bottomRows(1) = additional.transpose();
}


void Seminorm::solve_for_alpha(const std::vector<double> &data) {
    // extend data by the additional constant 0
    // TODO: data is extended in this function but stays the same in GravimetriyInversion
    auto data_extended = data;
    data_extended.emplace_back(0.);
    Norm::solve_for_alpha(data_extended);;

}


std::vector<Result>
Seminorm::calculate_density_distribution(const std::vector<double> &depth, uint64_t num_steps) {
    std::vector<Result> density_variable = Norm::calculate_density_distribution(depth, num_steps);
    auto density_constant = alpha.back();
    std::for_each(density_variable.begin(), density_variable.end(), [density_constant](Result& x){x.density +=density_constant; });
    return density_variable;
}


