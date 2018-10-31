//
// Class to perform a gravimetry inversion
//

// standard library includes
#include <fstream>
#include <iostream>
#include <algorithm>
// additional includes from installed libraries
#include <Eigen/Dense>
#include <Eigen/LU>
// includes from own code
#include "MeasurementData.h"
#include "GravimetryInversion.h"
#include "utils.h"



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


std::function<double(double)> operator*(const Representant& a, const Representant& b){
    return [&a, &b](double arg) { return a(arg) * b(arg);};
}

void GravimetryInversion::print_data(){
    std::cout << data;
}


void GravimetryInversion::calculate_gram_matrix() {
    // create Representants for all data values
    std::vector<Representant> repr;
    repr.reserve(data.size());
    for (auto el : data){
        repr.push_back(Representant(el.depth));
    }
    // create integrator and bind common arguments (limits, steps) to it
    // set maximum integration limit to the maximum measured depth
    Integrator integrate;
    auto gram_integrate = std::bind(integrate, std::placeholders::_1, LOWER_LIMIT, data.back().depth, INTEGRAL_STEPS);
    // only unique elements in the gram matrix are the diagonal elements because of
    // g_ij = min(g_i, g_j)
    std::vector<double> gram_matrix_diag_elements;
    gram_matrix_diag_elements.reserve(data.size());
    for (auto el : repr){
        gram_matrix_diag_elements.push_back(gram_integrate(el*el));
    }
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


Representant::Representant(double zj) : zj(zj) {}


double Representant::operator()(double z) const{
    return -gamma * heaviside(zj -z);
}


