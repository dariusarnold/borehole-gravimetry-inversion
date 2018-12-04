//
// Class to perform a gravimetry inversion
//


// standard library includes
#include <fstream>
// additional includes from installed libraries
#include <Eigen/Dense>
// includes from own code
#include "GravimetryInversion.h"
#include "utils.h"
#include "FileIO.h"
#include "Norms.h"


GravimetryInversion::GravimetryInversion(std::unique_ptr<ErrorNorm> _norm, uint64_t _discretization_steps) :
        norm(std::move(_norm)),
        discretization_steps(_discretization_steps),
        result(){}


void GravimetryInversion::calculate_density_distribution() {
    result = norm->calculate_density_distribution(discretization_steps);
}


void GravimetryInversion::write_density_distribution_to_file(const fs::path& filepath) {
    FileIO fw;
    fw.writeData(result, filepath);
}