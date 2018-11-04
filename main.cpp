#include <iostream>

#include "MeasurementData.h"
#include "GravimetryInversion.h"



void gravimetry_inversion(const std::string& filepath, uint64_t steps){
    GravimetryInversion mr(steps);
    mr.read_measurements_file(filepath);
    //mr.print_data();
    mr.calculate_gram_matrix();
    //mr.print_gram();
    mr.solve_alpha();
    //mr.print_alpha();
    // find file ending
    auto s = filepath.rfind(".dat");
    // create new string since filepath is const, append _density to file ending
    std::string out_fname = filepath;
    out_fname.replace(s, 4, std::string("_density.dat"));
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
        std::cout << steps << std::endl;
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