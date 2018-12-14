#include <iostream>
#include <vector>

#include "utils.h"

Eigen::VectorXd std_to_eigen(const std::vector<double> &vec) {
    Eigen::VectorXd eigen_vec = Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
    return eigen_vec;
}

std::vector<double> eigen_to_std(Eigen::VectorXd vec) {
    auto std_vec = std::vector<double>(&vec[0], vec.data()+vec.cols()*vec.rows());
    return std_vec;
}
