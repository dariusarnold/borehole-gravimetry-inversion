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
        measurement_depths(),
        measurement_data(),
        gram_matrix(),
        alpha(),
        result(){}


void GravimetryInversion::read_measurements_file(const fs::path& filepath) {
    std::ifstream matrixFile(filepath);
    if (matrixFile.is_open()){
        MeasurementData row;
        while (matrixFile >> row) {
            measurement_depths.push_back(row.depth);
            measurement_data.push_back(row.grav);
        }
    }
    matrixFile.close();
}


void GravimetryInversion::calculate_gram_matrix() {
    gram_matrix = norm->gram_matrix_analytical(measurement_depths);
}


void GravimetryInversion::solve_alpha(){
    alpha = norm->solve_for_alpha(measurement_data, gram_matrix);
}


void GravimetryInversion::calculate_density_distribution() {
    result = norm->calculate_density_distribution(alpha, measurement_depths, discretization_steps);
}


void GravimetryInversion::write_density_distribution_to_file(const fs::path& filepath) {
    FileWriter fw;
    fw.writeData(result, filepath);
}