//
// Created by darius on 11/11/18.
//

#ifndef GRAVITYINVERSION_NORM_H
#define GRAVITYINVERSION_NORM_H

#include <vector>
#include <Eigen/Dense>
#include "Result.h"



struct ErrorNorm{
    ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    virtual ~ErrorNorm() = 0;
    /**
     * Do the inversion for a constant nu and a given threshold
     * @param nu Lagrange multiplicator, larger nu weights the misfit over the norm,
     * smaller nu weights the norm over the misfit
     * @param threshold_squared Threshold of the misfit.
     */
    virtual void do_work(double nu);
    /**
     * Do the inversion for a misfit threshold of T² = N, where N is the number of measurement points.
     * Optimal lagrange multiplicator is determined by bisection search, so that the misfit fully uses the threshold.
     */
    virtual void do_work();
    virtual std::vector<Result> calculate_density_distribution(uint64_t num_steps);
    /**
     * Calculate the misfit X² using the formula nu^-2 | sigma * alpha |²
     * @param nu
     * @return
     */
    virtual double calculate_misfit(double nu);
protected:
    virtual double calculate_norm();
    virtual double representant_function(double zj, double z) = 0;
    virtual double gram_entry_analytical(double zj, double zk) = 0;
    virtual void gram_matrix_analytical();
    virtual void solve_for_alpha(double nu);
    double calc_nu_bysection(double nu_left, double nu_right, double desired_misfit);

    const std::vector<double>& measurement_depths;
    const std::vector<double>& measurement_data;
    const std::vector<double>& measurement_errors;
    const double gamma = 0.08382;   // mGal m^−1 g^−1 cm^3
    // Will hold gram matrix and if required be extended to include additional columns/rows
    Eigen::MatrixXd gram_matrix;
    // will hold error matrix sigma. Sigma is diagonal with error values sigma_1 sigma_2 etc on the diagonal
    Eigen::MatrixXd sigma_matrix;
    // will hold coefficients alpha
    Eigen::VectorXd alpha;
};


struct L2ErrorNorm : public ErrorNorm{
    L2ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    ~L2ErrorNorm() override;
    double representant_function(double zj, double z) override;
    double gram_entry_analytical(double zj, double zk) override;
};


#endif //GRAVITYINVERSION_NORM_H
