//
// Class to perform a gravimetry inversion
//

#define DEBUG 1

// standard library includes
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
// additional includes from installed libraries
#include <Eigen/Dense>
#include <Eigen/LU>
// includes from own code
#include "MeasurementData.h"
#include "GravimetryInversion.h"
#include "utils.h"
#include "FileWriter.h"



void GravimetryInversion::invert_data_from_file(const std::string &filepath) {
    GravimetryInversion mr;
    mr.read_measurements_file(filepath);
    #ifdef DEBUG
        mr.print_data();
    #endif
    mr.calculate_gram_matrix();
    #ifdef DEBUG
        mr.print_gram();
    #endif
    mr.solve_alpha();
    #ifdef DEBUG
        mr.print_alpha();
    #endif
    mr.write_density_distribution_to_file(filepath);
}


void GravimetryInversion::read_measurements_file(const std::string& filepath) {
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


void GravimetryInversion::calculate_gram_matrix() {
    // create Representants for all data values
    representant_functions.reserve(data.size());
    for (auto el : data){
        representant_functions.emplace_back(el.depth);
    }

    // create integrator and bind common arguments (limits, steps) to it,
    // creating a new function that only takes one argument, the function to integrate
    // set maximum integration limit to the maximum measured depth
    Integrator integrate;
    auto gram_integrate = std::bind(integrate, std::placeholders::_1, LOWER_LIMIT, data.back().depth, INTEGRAL_STEPS);
    // only unique elements in the gram matrix are the diagonal elements because of
    // g_ij = min(g_i, g_j)
    // calculate these unique elements first and fill gram matrix with them later
    std::vector<double> gram_matrix_diag_elements;
    gram_matrix_diag_elements.reserve(data.size());

    std::transform(representant_functions.begin(), representant_functions.end(), std::back_inserter(gram_matrix_diag_elements), [&gram_integrate](Representant r){return gram_integrate(r*r);});

    // create matrix with as many columns/rows as data entries read from file
    // and fill it with the values from the diagonals in this pattern:
    // g11 g11 g11
    // g11 g22 g22
    // g11 g22 g33
    gram_matrix.resize(data.size(), data.size());
    for (int column_index = 0; column_index < data.size(); ++column_index){
        for (int row_index = 0; row_index < data.size(); ++row_index){
            gram_matrix(row_index, column_index) = gram_matrix_diag_elements[std::min(row_index, column_index)];
        }
    }
}


void GravimetryInversion::solve_alpha(){
    // create eigen::vector and copy gravity measurements into it
    Eigen::VectorXd data_vec(data.size());
    for (int i = 0; i < data.size(); ++i){
        data_vec(i) = data[i].grav;
    }
    // use Eigen to solve the matrix equation
    Eigen::VectorXd alpha_eigen = gram_matrix.colPivHouseholderQr().solve(data_vec);
    //convert result from Eigen type to std::vector
    alpha = std::vector<double>(&alpha_eigen[0], alpha_eigen.data()+alpha_eigen.cols()*alpha_eigen.rows());
}


void GravimetryInversion::write_density_distribution_to_file(const std::string& filepath) {
    // fill depth vector with ascending values
    depth_meters.reserve(GravimetryInversion::INTEGRAL_STEPS);
    double stepsize = data.back().depth / INTEGRAL_STEPS;
    for (std::vector<double>::size_type i = 0; i != INTEGRAL_STEPS; ++i){
        depth_meters.emplace_back(stepsize * i);
    }
    density.reserve(GravimetryInversion::INTEGRAL_STEPS);
    for (int i = 0; i != INTEGRAL_STEPS; ++i){
        double dens = 0;
        for (int j = 0; j != representant_functions.size(); ++j){
            dens += alpha[j] * representant_functions[j](depth_meters[i]);
        }
        density.emplace_back(dens);
    }
    FileWriter fw;
    fw.writeData(depth_meters, density, filepath);
}
}


Representant::Representant(double zj) : zj(zj) {}


double Representant::operator()(double z) const{
    return -gamma * heaviside(zj -z);
}


