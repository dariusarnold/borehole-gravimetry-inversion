//
// Class to perform a gravimetry inversion
//


// standard library includes
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <experimental/filesystem>
// additional includes from installed libraries
#include <Eigen/Dense>
#include <Eigen/LU>
// includes from own code
#include "MeasurementData.h"
#include "GravimetryInversion.h"
#include "utils.h"
#include "FileWriter.h"


GravimetryInversion::GravimetryInversion(uint64_t _discretization_steps) :
        discretization_steps(_discretization_steps),
        data(),
        gram_matrix(),
        alpha(),
        density(),
        depth_meters(){}


void GravimetryInversion::invert_data_from_file_L2_norm(fs::path &filepath, uint64_t steps) {
    GravimetryInversion mr(steps);
    mr.read_measurements_file(filepath);
    mr.calculate_gram_matrix_L2_norm();
    mr.solve_alpha();
    mr.calculate_density_distribution();
    filepath.replace_extension({".dens"});
    mr.write_density_distribution_to_file(filepath);
}


void GravimetryInversion::read_measurements_file(const fs::path& filepath) {
    std::ifstream matrixFile(filepath);
    std::vector<MeasurementData> output;
    if (matrixFile.is_open()){
        MeasurementData row;
        while (matrixFile >> row) {
            output.push_back(row);
        }
    }
    matrixFile.close();
    data = output;
}


void GravimetryInversion::print_data(){
    std::cout << data;
}


void GravimetryInversion::print_gram() {
    std::cout << gram_matrix << std::endl;
}


void GravimetryInversion::print_alpha() {
    std::cout << alpha << std::endl;
}


void GravimetryInversion::calculate_gram_matrix_L2_norm() {
    std::vector<double> gram_matrix_diag_elements;
    gram_matrix_diag_elements.reserve(data.size());
    // only unique elements in the gram matrix are the diagonal elements because of
    // g_ij = min(g_i, g_j)
    // calculate these unique elements first and fill gram matrix with them later
    // use the analytical expression for the gram matrix Gamma_jk = gamma²*min(z_j, z_k)
    auto gram_analytical = [this](double zj, double zk){ return gamma*gamma*std::min(zj, zk); };
    for (auto el : data){
        gram_matrix_diag_elements.emplace_back(gram_analytical(el.depth, el.depth));
    }
    // create matrix with as many columns/rows as data entries read from file
    // and fill it with the values from the diagonals in this pattern:
    // g11 g11 g11
    // g11 g22 g22
    // g11 g22 g33
    gram_matrix.resize(data.size(), data.size());
    for (vec_size_t column_index = 0; column_index < data.size(); ++column_index){
        for (vec_size_t row_index = 0; row_index < data.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_matrix_diag_elements[std::min(row_index, column_index)];
        }
    }
}


void GravimetryInversion::solve_alpha(){
    // create eigen::vector and copy gravity measurements into it
    // TODO find better way to create Eigen vector
    Eigen::VectorXd data_vec(data.size());
    for (vec_size_t i = 0; i < data.size(); ++i){
        data_vec(i) = data[i].grav;
    }
    // use Eigen to solve the matrix equation
    Eigen::VectorXd alpha_eigen = gram_matrix.colPivHouseholderQr().solve(data_vec);
    //convert result from Eigen type to std::vector
    alpha = std::vector<double>(&alpha_eigen[0], alpha_eigen.data()+alpha_eigen.cols()*alpha_eigen.rows());
}


void GravimetryInversion::write_density_distribution_to_file(const fs::path& filepath) {
    FileWriter fw;
    fw.writeData(depth_meters, density, filepath);
}

void GravimetryInversion::calculate_density_distribution() {
    // fill depth vector with ascending values
    depth_meters.reserve(discretization_steps);
    double stepsize = data.back().depth / discretization_steps;
    for (vec_size_t i = 0; i != discretization_steps; ++i){
        depth_meters.emplace_back(stepsize * i);
    }
    // discretize density distribution by evaluating the following formula
    // rho(z) = sum_k alpha_k g_k(z)
    density.reserve(discretization_steps);
    for (vec_size_t i = 0; i != discretization_steps; ++i){
        double dens = 0;
        for (vec_size_t j = 0; j != alpha.size(); ++j){
            dens += alpha[j] * heaviside(data[j].depth - depth_meters[i]);
        }
        density.emplace_back(-gamma*dens);
    }
}