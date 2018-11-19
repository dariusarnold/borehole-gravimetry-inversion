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
     * Do all the work required by this norm to set up and solve the equation system
     * @param depth vector of depths or x values
     * @param data vector of measurement data points or function values associated with the x values
     */
    virtual void do_work(const std::vector<double>& depth, const std::vector<double>& data);
    /**
     * Calculate gram matrix analytically. Derived classes have to implement gram_entry_analytical
     * @param depth Vector containing measurement depth values in m
     * @return Gram matrix
     */
    virtual std::vector<Result>  calculate_density_distribution(const std::vector<double>& depth, uint64_t num_steps);
    /**
     * Calculate the coefficients alpha by solving the system of equations given as:
     * \vec{d} = \Gamma \vec{\alpha}}
     * This solves it appropiately in the default case, if a norm modifies the data vector,
     * it will have to override this function to implement it
     * @param data
     * @param gram_matrix
     * @return
     */
protected:
    virtual void solve_for_alpha(const std::vector<double>& data);

    virtual void gram_matrix_analytical(const std::vector<double>& depth);
    /**
     * Evaluate alpha and representants to discretize a density distribution
     * @param alpha Vector of coefficients alpha
     * @param representants vector of representant functions g_j
     * @param num_steps Number of steps to use for discretizing density distribution
     * @return
     */
    /**
     * Calculate a single entry (row index j, column index k) of the Gram matrix
     * @param zj depth in meters
     * @param zk depth in meters
     * @return
     */
    virtual double gram_entry_analytical(double zj, double zk) = 0;
    /**
     * Evaluate representant g_j(z) at depth z
     * @param zj constant depth used to build representant function
     * @param z depth in meters
     * @return
     */
    virtual double representant_function(double zj, double z) = 0;
    /**
     * Constant for 4 * pi * Gravity constant
     */

    Eigen::MatrixXd gram_matrix;
    std::vector<double> alpha;

    const double gamma = 0.08382;   // mGal m^−1 g^−1 cm^3

};


struct L2_Norm : public Norm{
    L2_Norm() = default;
    ~L2_Norm() override = default;

    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;
};


struct W12_Norm : public Norm{
    W12_Norm() = default;
    ~W12_Norm() override = default;

    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;
};


struct Seminorm : public Norm{
    Seminorm() = default;
    ~Seminorm() override = default;

    /**
     * Override base class since gram matrix is build differently
     * @param depth
     * @return
     */
    void gram_matrix_analytical(const std::vector<double>& depth) override;
    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;
    /**
     * Override solving the linear equation system since data vector has to be modified before solving
     * @param data
     * @param gram_matrix
     * @return
     */
    void solve_for_alpha(const std::vector<double> &data) override;

    std::vector<Result>  calculate_density_distribution(const std::vector<double>& depth, uint64_t num_steps) override;
};


#endif //GRAVITYINVERSION_NORM_H
