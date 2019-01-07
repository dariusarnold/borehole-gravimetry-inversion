import numpy as np
import math
import pyGravInv
import matplotlib.pyplot as plt



def matrix_S_inverted(delta_x, gamma, shape):
    """
    Create matrix S⁻1 that represents the derivative
    :param delta_x: Distance between discretization steps on depth axis
    :param gamma: 0 < gamma << 1, otherwise we would have a Seminorm and S would not be invertible
    :param shape: tuple of number of rows/columns to construct S with
    :return:
    """
    # Fill matrix with ones, then zero out top triangle
    S_inverted = np.tril(np.ones(shape))
    # fill left column with 1/gamma
    S_inverted[:, 0] = 1/gamma
    return delta_x * S_inverted


def matrix_Sigma_inverted(errors):
    """
    Create inverted error matrix Sigma, which holds 1/error in the diagonals and is zero everywhere else
    :param errors: Numpy array of measurement errors
    :return:
    """
    shape = (len(errors), len(errors))
    Sigma_inverted = np.zeros(shape)
    # get indices of diagonal and reassign values
    Sigma_inverted[np.diag_indices_from(Sigma_inverted)] = 1/errors
    return Sigma_inverted


def matrix_B(measurement_depths, inversion_depths):
    """
    Calculate B matrix from B rows.
    :return: B matrix
    """

    def calculate_B_row(depth_j):
        """
        Calculate a row of the B matrix
        :param depth_distribution: ndarray containing inversion depths
        :param depth_j: measurement depth for which to calculate row
        :return: row of B
        """
        gamma = 0.08382  # (mGal cm^3)/(m g)
        z = inversion_depths
        alpha = np.searchsorted(z, depth_j) - 1
        del_j = (depth_j - z[alpha])/(z[alpha+1] - z[alpha])

        B_j = np.zeros(z.size)
        if alpha == 0:
            B_j[alpha] = (2-del_j)*del_j*(z[alpha+1] - z[alpha])
        else:
            B_j[0] = (z[1]-z[0])
            for i in range(1, alpha):
                B_j[i] = (z[i]-z[i-1]) + (z[i+1]-z[i])
            B_j[alpha] = (z[alpha]-z[alpha-1]) + (2 - del_j)*del_j*(z[alpha+1] - z[alpha])
        B_j[alpha+1] = del_j**2*(z[alpha+1] - z[alpha])
        B_j *= -0.5*gamma
        return B_j

    B = np.zeros((measurement_depths.size, inversion_depths.size))
    for i, zj in enumerate(measurement_depths):
        B[i] = calculate_B_row(zj)
    return B


def matrix_A(measurement_depths, measurement_errors, inversion_depths):
    """
    Create NxL matrix A
    :param measurement_depths: np array of depths where values were measured
    :param measurement_errors: np array of measurement errors, same size as measurement_depths
    :param inversion_depths: np array of depths where the inversion result was discretized
    :return:
    """
    Sigma_inverted = matrix_Sigma_inverted(measurement_errors)
    B = matrix_B(measurement_depths, inversion_depths)
    S_shape = (len(inversion_depths), len(inversion_depths))
    dx = inversion_depths[1] - inversion_depths[0]
    S_inverted = matrix_S_inverted(gamma=1E-3, shape=S_shape, delta_x=dx)
    A = Sigma_inverted @ B @ S_inverted
    return A


def do_SVD(fname):
    """

    :param fname:
    :return:
    """
    # read/calculate required values
    measurement_depths, _, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths)

    # create Matrix A
    A = matrix_A(measurement_depths, measurement_errors, inversion_depths)
    # do svd
    V, Lambda, U_transposed = np.linalg.svd(A, full_matrices=False)
    return V, Lambda, U_transposed


def get_NL(fname):
    measurement_depths, _, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths)
    # N: number of measurement points
    N = len(measurement_depths)
    # L : number of inversion depths
    L = len(inversion_depths)
    return N, L


def print_check_result(name, check):
    print_helper = lambda condition: "OK" if condition else "FAIL"
    s = f"{name:<14} {print_helper(check):.>6}"
    print(s)


def ex11_1(fname):
    V, Lambda, U_transposed = do_SVD(fname)
    N, L = get_NL(fname)

    # check if the conditions are fullfilled
    V_transposed = np.transpose(V)
    VTV = V_transposed @ V
    VVT = V @ V_transposed
    I_N = np.identity(N)
    I_L = np.identity(L)
    print_check_result(name="V^T*V == V*V^T", check=np.allclose(VTV, VVT))
    print_check_result(name="V^T*V == I_N", check=np.allclose(VTV, I_N))
    #next condition
    U = np.transpose(U_transposed)
    UUT = U @ U_transposed
    UTU = U_transposed @ U
    print_check_result(name="U^T*U == I_N", check=np.allclose(UTU, I_N))
    print_check_result(name="U*U^T != I_L", check=not np.allclose(UUT, I_L))


def get_inversion_depth_steps(measurement_depths, num_steps=10000):
    """
    Get the discretization depths of the inversion result
    :param measurement_depths: np array containing measurement depths
    :param num_steps: Number of steps used for discretozation
    :return:
    """
    # add 2 to get the same result as the C++ code, where two additional steps are added so the effect for some norms
    # that occurs when going past the last measurement point is seen
    num_steps += 2
    return np.linspace(0, measurement_depths[-1], num_steps)


def new_data(old_data, V_transposed, measurement_errors):
    return (1/measurement_errors) * V_transposed @ old_data


def misfit_squared(new_data, nu, Lambda):
    """
    Calculate new misfit
    :param new_data: d''
    :param nu: Lagrange parameter
    :param Lambda: Vector containing the Eigenvalues of the SVD
    :return: Chi²
    """
    zwischen = (new_data / (1 + nu*Lambda**2))**2
    return np.sum(zwischen)


def ex11_2(fname):
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    V, Lambda, U_transposed = do_SVD(fname)
    # calculate new data
    d_double_prime = new_data(measurement_data, np.transpose(V),  measurement_errors)
    # calculate misfit for multiple values of nu
    range_of_nus = np.logspace(-2, 2, 1000)
    range_of_chis = np.array([misfit_squared(d_double_prime, nu, Lambda) for nu in range_of_nus])
    plt.plot(range_of_nus, range_of_chis, ".")
    plt.title(r"$\nu/\chi^2$")
    plt.xlabel(r"Lagrangeparameter $\nu$")
    plt.ylabel(r"Misfit $\chi^2$")
    plt.show()


def calculate_coeffs(new_data, nu, Lambda):
    return new_data / ( (1/nu * 1/Lambda) + Lambda )


def calculate_model(S_inverted, U, new_coefficients):
    return S_inverted @ U @ new_coefficients


def optimal_nu_bysection(target_misfit, new_data, Lambda, accuracy=0.01):
    """
    Calculate the optimal Lagrange parameter using a bisection method.
    The optimum nu maximizes the misfit.
    :param target_misfit:
    :return:
    """
    import sys
    # set borders of interval in which to search
    nu_left = sys.float_info.epsilon
    nu_right = 1E6
    chi_squared_left = misfit_squared(new_data, nu_left, Lambda)
    chi_squared_right = misfit_squared(new_data, nu_right, Lambda)
    # misfit over nu is a monotonically falling function, check if target_misfit can be reached
    assert chi_squared_left > target_misfit > chi_squared_right, "Target misfit was set outside of possible range"

    nu_mid = (nu_left + nu_right) / 2
    chi_squared = misfit_squared(new_data, nu_mid, Lambda)
    while abs(chi_squared - target_misfit) > accuracy:
        nu_mid = (nu_left + nu_right) / 2
        chi_squared = misfit_squared(new_data, nu_mid, Lambda)
        if chi_squared > target_misfit:
            nu_left = nu_mid
        else:
            nu_right = nu_mid
    print(f"Optimal Lagrange parameter: {nu_mid}")
    return nu_mid


def correct_free_air_gradient(measurement_depths, measurement_data):
    f = 0.308 # mGal / m
    return measurement_data - f * measurement_depths


def read_data(fname):
    """
    Read measurement data with errors from fname and correct the data for the free air gradient
    :param fname:
    :return:
    """
    measurement_depths, measurement_data, measurement_errors = np.loadtxt(fname, unpack=True)
    measurement_data = correct_free_air_gradient(measurement_depths, measurement_data)
    return measurement_depths, measurement_data, measurement_errors


def ex11_3(fname):
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths)
    V, Lambda, U_transposed = do_SVD(fname)
    # calculate new data
    d_double_prime = new_data(measurement_data, np.transpose(V),  measurement_errors)
    # calculate desired misfit using bisection
    desired_misfit = len(measurement_depths)
    optimal_nu = optimal_nu_bysection(desired_misfit, d_double_prime, Lambda)
    # calculate new coefficients
    alpha_double_prime = calculate_coeffs(d_double_prime, optimal_nu, Lambda)
    S_shape = (len(inversion_depths), len(inversion_depths))
    dx = inversion_depths[1] - inversion_depths[0]
    S_inverted = matrix_S_inverted(gamma=1E-3, shape=S_shape, delta_x=dx)
    dens_model = calculate_model(S_inverted, np.transpose(U_transposed), new_coefficients=alpha_double_prime)
    plt.plot(inversion_depths, dens_model)
    plt.title("Density Model")
    plt.xlabel(r"Depth $z$ (m)")
    plt.ylabel(r"Density $\rho$ (g/cm³)")
    plt.show()


def ex11_4(fname):
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths)
    V, Lambda, U_transposed = do_SVD(fname)
    S_shape = (len(inversion_depths), len(inversion_depths))
    dx = inversion_depths[1] - inversion_depths[0]
    S_inverted = matrix_S_inverted(gamma=1E-3, shape=S_shape, delta_x=dx)
    U = np.transpose(U_transposed)
    representants = S_inverted @ U
    # representants now holds the representants in its columns, but it is set up so that indexing it accesses its rows
    # transpose it to access the representants by index
    representants = np.transpose(representants)
    f, ax_arr = plt.subplots(len(representants), sharex=True)
    for i, repr in enumerate(representants):
        ax_arr[i].plot(repr)
    plt.show()


def main():
    selection = input("Which task? (1, 2, 3, 4)\n")
    selection = list(selection.strip())
    for s in selection:
        try:
            function_name = getattr(__import__(__name__), f"ex11_{s}")
        except AttributeError:
            print(f"Task {s} is not available")
            return
        function_name("data/grav15.dat")


if __name__ == '__main__':
    main()