#ifndef GRAVITYINVERSION_LIBRARY_H
#define GRAVITYINVERSION_LIBRARY_H


#include <Eigen/Dense>
#include "ModelParameters.h"


enum class Norms{L2Norm, W12Norm, SemiNorm};
enum class ErrorNorms{L2ErrorNorm, SemiErrorNorm};
enum class InterpolationNorms{LinearInterpolationNorm};

using InvResult = std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters>;



/**
 * Do an inversion with measurement error free data and return density model and inversion model parameters.
 * @param filepath Depth, data pairs are read from this file.
 * @param norm_id Gives the norm to be used
 * @param discretization_steps Number of points to use for the discretization of the resulting density model.
 * @return
 */
InvResult inversion(const std::string& filepath, Norms norm_id, uint64_t discretization_steps=10000);


/**
 * Do an inversion with measurement data with errors and return density model and inversion model parameters.
 * @param filepath Depth, data, measurement error triples are read from this file.
 * @param norm_id Gives the norm to be used
 * @param discretization_steps Number of points to use for the discretization of the resulting density model.
 * @param nu Lagrange multiplicator nu. If not given the default value -1 means that an optimal value should be
 * found
 * @return Two vectors containing the depth/density values, a ModelParameters instance that holds other results of the
 * inversion.
 */
InvResult inversion_error(const std::string& filepath, ErrorNorms norm_id, uint64_t discretization_steps=10000, double nu=-1);


/**
 * Interpolate between x, f(x) pairs from a file and return the discretized interpolation result.
 * @param filepath
 * @param norm_id
 * @param a lower boundary of evaluation interval
 * @param b upper boundary of evaluation interval
 * @param discretization_steps
 * @return
 */
InvResult interpolation(const std::string& filepath, InterpolationNorms norm_id, double a, double b, uint64_t discretization_steps=10000);



#endif //GRAVITYINVERSION_LIBRARY_H
