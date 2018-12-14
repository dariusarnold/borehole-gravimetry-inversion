#include "Library.h"
#include "GravimetryInversion.h"
#include <experimental/filesystem>



namespace fs = std::experimental::filesystem;

std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>
        inversion_error(const std::string& filepath, Norms norm_id, uint64_t discretization_steps, double nu){
    auto fpath = fs::path(filepath);
    switch(norm_id){
        case Norms::L2ErrorNorm:
            return GravimetryInversion<L2ErrorNorm>::invert_data_from_file_with_errors(fpath, discretization_steps, nu);
            break;
        case Norms::SemiErrorNorm:
            return GravimetryInversion<SemiErrorNorm>::invert_data_from_file_with_errors(fpath, discretization_steps, nu);
            break;
    }
}


int f(std::string s){
    return s.length();
}

