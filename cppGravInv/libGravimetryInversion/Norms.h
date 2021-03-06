#ifndef GRAVITYINVERSION_NORM_H
#define GRAVITYINVERSION_NORM_H

#include <vector>
#include <Eigen/Dense>


// forward declaration
struct ModelParameters;


struct Norm {
    Norm(const std::vector<double>& depth, const std::vector<double>& data);
    virtual ~Norm() = default;
    /**
     * Do all the work required by this norm to set up and solve the equation system
     * @param depth vector of depths or x values
     * @param data vector of measurement data points or function values associated with the x values
     */
    virtual std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> do_work(uint64_t discretization_steps);
protected:
    /**
    * Evaluate alpha and representants to discretize a density distribution
    * @param alpha Vector of coefficients alpha
    * @param num_steps Number of steps to use for discretizing density distribution
    * @return
    */
    virtual std::pair<Eigen::VectorXd, Eigen::VectorXd> calculate_density_distribution(uint64_t num_steps);
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
    Eigen::VectorXd alpha;
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

    std::pair<Eigen::VectorXd, Eigen::VectorXd> calculate_density_distribution(uint64_t num_steps) override;
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
    std::pair<Eigen::VectorXd, Eigen::VectorXd> calculate_density_distribution(uint64_t num_steps) override;
private:
    double a;
    double b;
};


// abstract base class for all norms that handle data with measurement errors
struct ErrorNorm : public Norm{
    ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    ~ErrorNorm() override = default;
    /**
     * Do the inversion for a constant nu and a given threshold
     * @param nu Lagrange multiplicator, larger nu weights the misfit over the norm,
     * smaller nu weights the norm over the misfit
     * @param discretization_steps number of steps used to discretize the inverted density model
     */
    std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> do_work(uint64_t discretization_steps, double nu);
    /**
     * Do the inversion for a misfit threshold of T² = N, where N is the number of measurement points.
     * Optimal lagrange multiplicator is determined by bisection search, so that the misfit fully uses the threshold.
     * @param discretization_steps number of steps used to discretize the inverted density model
     */
    std::tuple<std::pair<Eigen::VectorXd, Eigen::VectorXd>, ModelParameters> do_work(uint64_t discretization_steps) override;
protected:
    /**
     * Calculate the misfit X² using the formula nu^-2 | sigma * alpha |²
     * @param nu
     * @return
     */
    virtual double calculate_misfit(double nu);
    virtual double calculate_norm();
    virtual void solve_for_alpha(double nu);
    /**
     * The Lagrange parameter nu weighs between the norm and the misfit of a model.
     * Misfit is a strictly monotonically falling function with regard to the lagrange parameter nu.
     * Therefore a bisection search can be used to find the ideal Lagrange parameter where the misfit of the inversion 
     * model is equal to the desired_misfit. 
     * The function will throw a std::range_error when the desired misfit can't be reached with the values of nu in the
     * interval between nu_left and nu_right
     * @param nu_left Start value of the left border of the interval
     * @param nu_right Start value of the right border of the interval
     * @param desired_misfit Target misfit value to reach
     * @return Optimal Lagrange parameter nu. 
     */
    double calc_nu_bysection(double nu_left, double nu_right, double desired_misfit, double accuracy = 0.01);

    std::vector<double> measurement_errors;
    // will hold error matrix sigma. Sigma is diagonal with error values sigma_1 sigma_2 etc on the diagonal
    Eigen::MatrixXd sigma_matrix;
};


struct L2ErrorNorm : public ErrorNorm{
    L2ErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    ~L2ErrorNorm() override;
    double representant_function(double zj, double z) override;
    double gram_entry_analytical(double zj, double zk) override;
};


struct SemiErrorNorm : public ErrorNorm{
    SemiErrorNorm(const std::vector<double>& depth, const std::vector<double>& data, const std::vector<double>& errors);
    ~SemiErrorNorm() override = default;
    double representant_function(double zj, double z) override;
    double gram_entry_analytical(double zj, double zk) override;
    void gram_matrix_analytical(double nu);
    double calculate_misfit(double nu) override;
    void solve_for_alpha(double nu) override;
    std::pair<Eigen::VectorXd, Eigen::VectorXd> calculate_density_distribution(uint64_t num_steps) override;
};

#endif //GRAVITYINVERSION_NORM_H
