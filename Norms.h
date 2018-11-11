//
// Created by darius on 11/11/18.
//

#ifndef GRAVITYINVERSION_NORM_H
#define GRAVITYINVERSION_NORM_H

#include <vector>
#include <Eigen/Dense>
#include "Result.h"

struct Norm {
    Norm() = default;
    virtual ~Norm() = default;

    /**
     * Calculate gram matrix analytically
     * @param depth Vector containing measuremnt depth values in m
     * @return Gram matrix
     */
    virtual Eigen::MatrixXd gram_matrix_analytical(const std::vector<double>& depth) = 0;

    /**
     * Evaluate alpha and representants to discretize a density distribution
     * @param alpha Vector of coefficients alpha
     * @param representants vector of representant functions g_j
     * @param num_steps Number of steps to use for discretizing density distribution
     * @return
     */
    virtual std::vector<Result>  calculate_density_distribution(const std::vector<double>& alpha, const std::vector<double>& depth, uint64_t num_steps) = 0;

    const double gamma = 0.08382;
};


struct L2_Norm : Norm{
    L2_Norm() = default;
    ~L2_Norm() = default;
    Eigen::MatrixXd gram_matrix_analytical(const std::vector<double>& depth);

    std::vector<Result> calculate_density_distribution(const std::vector<double>& alpha, const std::vector<double>& depth, uint64_t num_steps);
};

/*
struct W12_Norm : Norm{
    W12_Norm() = default;
    ~W12_Norm() = default;
    Eigen::MatrixXd gram_matrix_analytical(const std::vector<double>& depth);

    std::vector<std::pair<double, double>> calculate_density_distribution(const std::vector<double>& alpha, const std::vector<double>& depth, uint64_t num_steps);
private:
    double gram_entry_analytical(double zj, double zk);
};
*/

#endif //GRAVITYINVERSION_NORM_H
