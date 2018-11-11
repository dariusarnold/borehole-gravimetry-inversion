//
// Class to perform a gravimetry inversion
//


// standard library includes
#include <fstream>
// additional includes from installed libraries
#include <Eigen/Dense>
#include <Eigen/LU>
// includes from own code
#include "MeasurementData.h"
#include "GravimetryInversion.h"
#include "utils.h"
#include "FileWriter.h"
#include "Norms.h"


GravimetryInversion::GravimetryInversion(std::unique_ptr<Norm> _norm, uint64_t _discretization_steps) :
        norm(std::move(_norm)),
        discretization_steps(_discretization_steps),
        data(),
        gram_matrix(),
        alpha(),
        result(){}


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


void GravimetryInversion::calculate_gram_matrix() {
    std::vector<double> depth;
    for (auto el : data){
        depth.emplace_back(el.depth);
    }
    gram_matrix = norm->gram_matrix_analytical(depth);
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


void GravimetryInversion::calculate_density_distribution() {
    std::vector<double> depth;
    for (auto el : data){
        depth.emplace_back(el.depth);
    }
    result = norm->calculate_density_distribution(alpha, depth, discretization_steps);
}


void GravimetryInversion::write_density_distribution_to_file(const fs::path& filepath) {
    FileWriter fw;
    fw.writeData(result, filepath);
}