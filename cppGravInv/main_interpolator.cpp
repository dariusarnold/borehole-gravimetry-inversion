#include <iostream>
#include <experimental/filesystem>

#include "GravimetryInversion.h"
#include "Norms.h"


namespace fs = std::experimental::filesystem;

//TODO refactor to use only the Library.h functions

/**
 * Call programm with positional arguments:
 * 1. input measurement data filepath
 * 2. number of discretization steps
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input data filepath\n"
                  << "2. number of discretization steps\n"
                  << "3. parameters a b (lower and upper bound of evaluation interval of interpolated function)"
                  << std::endl;
        return -1;
    }
    else{
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        double a = std::stod(argv[3]);
        double b = std::stod(argv[4]);
        GravimetryInversion<LinearInterpolationNorm>::interpolate_data_from_file(filepath, steps, a, b);
        return 0;
    }
}