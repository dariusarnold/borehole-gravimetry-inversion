#include <experimental/filesystem>
#include "Library.h"
#include "GravimetryInversion.h"
#include "Norms.h"


namespace fs = std::experimental::filesystem;

std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>
        inversion(const std::string& filepath, Norms norm_id, uint64_t discretization_steps){
    auto fpath = fs::path(filepath);
    switch (norm_id){
        case Norms::L2Norm:
            return GravimetryInversion<L2_Norm>::invert_data_from_file(fpath, discretization_steps);
        case Norms::W12Norm:
            return GravimetryInversion<W12_Norm>::invert_data_from_file(fpath, discretization_steps);
        case Norms::SemiNorm:
            return GravimetryInversion<Seminorm>::invert_data_from_file(fpath, discretization_steps);
        default:
            return GravimetryInversion<L2_Norm>::invert_data_from_file(fpath, discretization_steps);
    }
}


std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>
        inversion_error(const std::string& filepath, ErrorNorms norm_id, uint64_t discretization_steps, double nu){
    auto fpath = fs::path(filepath);
    switch(norm_id){
        case ErrorNorms::L2ErrorNorm:
            return GravimetryInversion<L2ErrorNorm>::invert_data_from_file_with_errors(fpath, discretization_steps, nu);
        case ErrorNorms::SemiErrorNorm:
            return GravimetryInversion<SemiErrorNorm>::invert_data_from_file_with_errors(fpath, discretization_steps, nu);
        default:
            return GravimetryInversion<L2ErrorNorm>::invert_data_from_file_with_errors(fpath, discretization_steps, nu);
    }
}


std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>
        interpolation(const std::string& filepath, InterpolationNorms norm_id, double a, double b, uint64_t discretization_steps){
    auto fpath = fs::path(filepath);
    switch(norm_id){
        case InterpolationNorms::LinearInterpolationNorm:
            return GravimetryInversion<LinearInterpolationNorm>::interpolate_data_from_file(fpath, discretization_steps, a, b);
        default:
            return GravimetryInversion<LinearInterpolationNorm>::interpolate_data_from_file(fpath, discretization_steps, a, b);
    }
}
