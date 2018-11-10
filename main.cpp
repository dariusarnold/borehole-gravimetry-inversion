#include <iostream>

#include "MeasurementData.h"
#include "GravimetryInversion.h"



void gravimetry_inversion(const std::string& filepath, uint64_t steps){
    GravimetryInversion mr(steps);
    mr.read_measurements_file(filepath);
    mr.calculate_gram_matrix_L2_norm();
    mr.solve_alpha();
    // find file ending
    auto s = filepath.rfind(".dat");
    // create new string since filepath is const, change file ending to .dens and save results
    std::string out_fname = filepath;
    out_fname.replace(s, 4, std::string(".dens"));
    mr.write_density_distribution_to_file(out_fname);
}


/**
 * Call programm with positional arguments:
 * 1. input measurement data filepath
 * 2. number of discretization steps
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    if (argc > 1){
        std::string filepath = argv[1];
        uint64_t steps = atoi(argv[2]);
        // input checking of steps, if less than 0 use default value
        steps > 0 ? gravimetry_inversion(filepath, steps) : gravimetry_inversion(filepath, 10000);
        return 0;
    }else {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps" << std::endl;
        return -1;
    }
}