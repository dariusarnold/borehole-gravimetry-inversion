#ifndef GRAVITYINVERSION_NORM_H
#define GRAVITYINVERSION_NORM_H

#include <vector>
#include <Eigen/Dense>
#include "Result.h"


struct Norm {
    Norm(const std::vector<double>& depth, const std::vector<double>& data);
    virtual ~Norm() = default;
    /**
     * Do all the work required by this norm to set up and solve the equation system
     * @param depth vector of depths or x values
     * @param data vector of measurement data points or function values associated with the x values
     */
    virtual std::vector<Result> do_work(uint64_t discretization_steps);
protected:
    /**
    * Evaluate alpha and representants to discretize a density distribution
    * @param alpha Vector of coefficients alpha
    * @param num_steps Number of steps to use for discretizing density distribution
    * @return
    */
    virtual std::vector<Result>  calculate_density_distribution(uint64_t num_steps);
    /**
     * Calculate the coefficients alpha by solving the system of equations given as:
     * \vec{d} = \Gamma \vec{\alpha}}
     * This solves it appropiately in the default case, if a norm modifies the data vector,
     * it will have to override this function to implement it
     * @param data measurement or function values
     * @return
     */
    virtual void solve_for_alpha();
    /**
     * Calculate gram matrix analytically. Derived classes have to implement gram_entry_analytical
     * @param depth Vector containing measurement depth values in m
     * @return Gram matrix
     */
    virtual void gram_matrix_analytical();
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
     * Constant for 4 * pi * Gravity constant, used for inverting borehole gravimetry data
     */
    const double gamma = 0.08382;   // mGal m^−1 g^−1 cm^3
    std::vector<double> measurement_depths;
    std::vector<double> measurement_data;
    // Will hold gram matrix and if required be extended to include additional columns/rows
    Eigen::MatrixXd gram_matrix;
    // will hold coefficients alpha
    // TODO make this an eigen vector
    std::vector<double> alpha;
};


struct L2_Norm : public Norm{
    L2_Norm(const std::vector<double>& depth, const std::vector<double>& data);
    ~L2_Norm() override = default;
    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;
};


struct W12_Norm : public Norm{
    W12_Norm(const std::vector<double>& depth, const std::vector<double>& data);
    ~W12_Norm() override = default;
    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;
};


struct Seminorm : public Norm{
    Seminorm(const std::vector<double>& depth, const std::vector<double>& data);
    ~Seminorm() override = default;
    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;

    /**
 * Override base class since gram matrix is build differently
 * @param depth
 * @return
 */
    void gram_matrix_analytical() override;
    /**
     * Override solving the linear equation system since data vector has to be modified before solving
     * @param data
     * @param gram_matrix
     * @return
     */
    void solve_for_alpha() override;

    std::vector<Result>  calculate_density_distribution(uint64_t num_steps) override;
};


struct LinearInterpolationNorm : public Norm{
    /**
     * Constructor specifying bounds of interval in which interpolated function is evaluated
     * @param _a Lower bound
     * @param _b Upper bound
     */

    LinearInterpolationNorm(double _a, double _b, const std::vector<double>& x_values, const std::vector<double>& y_values);
    ~LinearInterpolationNorm() override = default;
    double gram_entry_analytical(double zj, double zk) override;
    double representant_function(double zj, double z) override;

    void gram_matrix_analytical() override;
    void solve_for_alpha() override;
    std::vector<Result> calculate_density_distribution(uint64_t num_steps) override;
private:
    double a;
    double b;
};


struct ErrorNorm{
    ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    virtual ~ErrorNorm() = 0;
    /**
     * Do the inversion for a constant nu and a given threshold
     * @param nu Lagrange multiplicator, larger nu weights the misfit over the norm,
     * smaller nu weights the norm over the misfit
     * @param discretization_steps number of steps used to discretize the inverted density model
     */
    virtual std::vector<Result> do_work(double nu, uint64_t discretization_steps);
    /**
     * Do the inversion for a misfit threshold of T² = N, where N is the number of measurement points.
     * Optimal lagrange multiplicator is determined by bisection search, so that the misfit fully uses the threshold.
     * @param discretization_steps number of steps used to discretize the inverted density model
     */
    virtual std::vector<Result> do_work(uint64_t discretization_steps);
protected:
    virtual std::vector<Result> calculate_density_distribution(uint64_t num_steps);
    /**
     * Calculate the misfit X² using the formula nu^-2 | sigma * alpha |²
     * @param nu
     * @return
     */
    virtual double calculate_misfit(double nu);
    virtual double calculate_norm();
    virtual double representant_function(double zj, double z) = 0;
    virtual double gram_entry_analytical(double zj, double zk) = 0;
    virtual void gram_matrix_analytical();
    virtual void solve_for_alpha(double nu);
    double calc_nu_bysection(double nu_left, double nu_right, double desired_misfit);

    std::vector<double> measurement_depths;
    std::vector<double> measurement_data;
    std::vector<double> measurement_errors;
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
