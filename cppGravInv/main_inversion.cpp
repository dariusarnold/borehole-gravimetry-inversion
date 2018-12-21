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
    if (argc != 4) {
        std::cout << "Call programm with positional arguments:\n"
                  << "1. input measurement data filepath\n"
                  << "2. number of discretization steps\n"
                  << "3. Id of norm: L2 (0), W12 (1), Semi (2)"
                  << std::endl;
        return -1;
    }
    else{
        // read command line parameters
        fs::path filepath{argv[1]};
        uint64_t steps = std::stoul(argv[2]);
        std::pair<Eigen::VectorXd, Eigen::VectorXd> depth_dens;
        ModelParameters params;
        switch (atoi(argv[3])){
            case 0:
                std::tie(depth_dens, params) = GravimetryInversion<L2_Norm>::invert_data_from_file(filepath, steps);
                break;
            case 1:
                std::tie(depth_dens, params) = GravimetryInversion<W12_Norm>::invert_data_from_file(filepath, steps);
                break;
            case 2:
                std::tie(depth_dens, params) = GravimetryInversion<Seminorm>::invert_data_from_file(filepath, steps);
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