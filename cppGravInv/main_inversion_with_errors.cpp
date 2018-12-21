#include <iostream>
#include <experimental/filesystem>

#include "FileIO.h"
#include "Library.h"


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
    if ((argc != 4) and (argc != 5)) {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps\n"
                  << "3. Id of norm: L2 (0), W12Semi (1)\n"
                  << "4. Optional: Lagrange multiplicator nu, if not given optimum is determined"
                  << std::endl;
        return -1;
    }
    else{
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        // read lagrange multiplicator
        double nu = -1;
        if (argc == 5){
            nu = std::stod(argv[4]);
        }
        std::pair<Eigen::VectorXd, Eigen::VectorXd> depth_dens;
        ModelParameters params;
        switch (atoi(argv[3])){
            case 0:
                //std::tie(depth_dens, params) = GravimetryInversion<L2ErrorNorm>::invert_data_from_file_with_errors(filepath, steps, nu);
                std::tie(depth_dens, params) = inversion_error(filepath, ErrorNorms::L2ErrorNorm, steps, nu);
                break;
            case 1:
                //std::tie(depth_dens, params) = GravimetryInversion<SemiErrorNorm>::invert_data_from_file_with_errors(filepath, steps, nu);
                std::tie(depth_dens, params) = inversion_error(filepath, ErrorNorms::SemiErrorNorm, steps, nu);
                break;
            default:
                std::cout << "Enter valid norm id" << std::endl;
                return -1;
        }
        // print the Inversion parameters
        std::cout << params;
        // save Inversion results in file
        FileIO fw;
        filepath.replace_extension(".dens");
        fw.writeData(depth_dens, filepath);
        return 0;
    }
}