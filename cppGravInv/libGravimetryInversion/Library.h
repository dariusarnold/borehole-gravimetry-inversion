#ifndef GRAVITYINVERSION_LIBRARY_H
#define GRAVITYINVERSION_LIBRARY_H


#include <Eigen/Dense>
#include "ModelParameters.h"
#include <string>

enum class Norms{L2ErrorNorm, SemiErrorNorm};


std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>
        inversion_error(const std::string& filepath, Norms norm_id, uint64_t discretization_steps=10000, double nu=-1);


int f(std::string s);



#endif //GRAVITYINVERSION_LIBRARY_H
